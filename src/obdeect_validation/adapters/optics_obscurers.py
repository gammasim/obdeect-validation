"""Effective-area loss from the obscurers of an optical model (secondary baffle and structure).

Two parts, selected by the first argument.

``analytic`` (obdeect only, runs on the host). For each field angle the effective area of the
optical model (all obscurers on) over the effective area of its obscurer-free copy is traced with
the same photons and seed, geometry only. The obdeect census (``OpticalModel.obscurer_hits``) counts
the light each obscurer absorbs as the first obscurer on the path; the mass balance
``detected = unobscured - first`` says the ratio equals ``1 - sum(first) / unobscured``.
Gate: exact identity, tolerance 1e-6 absolute on the ratio. Both sides sum the same float32 weights
of the same photons in float64, so only round-off can separate them; the tolerance is that
round-off bound, chosen before the run.

``compare`` (container: simtools-validate-optics with both backends). Three variants per seed:
as modelled (``base``), secondary baffle removed (``secondary_mirror_baffle`` set to r1 = 0,
the sim_telarray keyword ``secondary_baffle``) and parameterised structure shadowing removed
(``telescope_transmission`` set to 1, i.e. no shadowing). The modified parameters are passed with
simtools ``overwrite_model_parameters``; the checkout is not changed. Per offset and backend:

    rho_baffle    = A(base) / A(no baffle)
    rho_structure = A(base) / A(no structure)

averaged over the seeds. Gate per offset and variant: |rho_obdeect - rho_sim_telarray| <=
limit, with limit 0.02, the optics requirement set by the project (goal 1 %, acceptable at most
2 %; offsets above the goal are reported). sigma is the pooled single-run seed-to-seed spread of rho
over the seeds of both backends (sqrt((s_obd^2 + s_sim^2) / 2)); the limit is not widened by it: an
offset with 3 sigma above the limit is reported as underpowered and fails. The
limit is applied to the ratio itself, not to the lost light: a 1 % relative gate on the loss would
be stricter (about 0.02 % absolute for the SCT baffle at 1 deg); that choice is left open.

The sim_telarray side of ``compare`` has not been run in this workspace (no sim_telarray build on
the host); see the case documentation.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

from obdeect_validation.adapters import common
from obdeect_validation.adapters.simtools_validate_optics import (
    AdapterError,
    _table,
    run_simtools,
)

LIMIT = 0.02  # optics requirement: acceptable at most 2 % (relative)
GOAL = 0.01  # optics requirement: goal 1 % (reported, not gated)
SIGMA_FACTOR = 3.0
ANALYTIC_TOLERANCE = 1e-6
BACKENDS = ("sim_telarray", "obdeect")
VARIANTS = ("base", "no_baffle", "no_structure")
SEEDS_REQUIRED = 3
BAFFLE = "secondary_mirror_baffle"
TRANSMISSION = "telescope_transmission"


StudyError = AdapterError


# --- configuration --------------------------------------------------------------------------


def _config(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise StudyError(f"cannot read config {path}: {exc}") from exc
    required = {"simtools_config", "analytic", "compare"}
    if not isinstance(value, dict) or set(value) != required:
        raise StudyError(f"config must contain exactly {sorted(required)}")
    base = value["simtools_config"]
    if (
        not isinstance(base, dict)
        or base.get("application") != "simtools-validate-optics"
    ):
        raise StudyError("simtools_config.application must be simtools-validate-optics")
    options = base.get("configuration")
    if not isinstance(options, dict) or "ray_tracing_backend" in options:
        raise StudyError("simtools_config.configuration must omit ray_tracing_backend")
    analytic = value["analytic"]
    if not isinstance(analytic, dict) or not {
        "model_version",
        "photons",
        "seed",
    } <= set(analytic):
        raise StudyError("analytic needs model_version, photons and seed")
    compare = value["compare"]
    if not isinstance(compare, dict) or not {"seeds", "offsets_deg"} <= set(compare):
        raise StudyError("compare needs seeds and offsets_deg")
    if len(compare["seeds"]) < SEEDS_REQUIRED:
        raise StudyError(f"compare needs at least {SEEDS_REQUIRED} seeds")
    return value


def _offsets(values: list[list[float]]) -> list[tuple[float, float]]:
    """Field offsets (x, y) in degrees; the star is on the axis for (0, 0)."""
    offsets = []
    for pair in values:
        if len(pair) != 2 or any(
            isinstance(v, bool) or not math.isfinite(v) for v in pair
        ):
            raise StudyError("offsets_deg entries must be [x, y] finite numbers")
        offsets.append((float(pair[0]), float(pair[1])))
    return offsets


# --- analytic part (obdeect only) -----------------------------------------------------------


def _towards(offset: tuple[float, float]) -> tuple[float, float, float]:
    """Unit vector to a star at the field offset (same convention as ``obdeect ray-tracing``)."""
    fx, fy = (math.radians(a) for a in offset)
    x, y, z = math.tan(fx), math.tan(fy), 1.0
    norm = math.sqrt(x * x + y * y + z * z)
    return (x / norm, y / norm, z / norm)


def _area(arrivals: dict[str, Any]) -> float:
    """Summed input weights of the detected photons (geometry only)."""
    import numpy as np
    from obdeect.ray_tracing.model import DETECTED

    detected = np.asarray(arrivals["status"]) == DETECTED
    return float(np.asarray(arrivals["weight"], dtype=float)[detected].sum())


def _without_obscurers(model) -> Any:
    """Obscurer-free copy of an optical model: no obscurers and no parameterised transmission."""
    from obdeect.ray_tracing.hashing import model_hash
    from obdeect.ray_tracing.model import OpticalModel

    document = model.document
    document["optics"]["obscurers"] = []
    document["optics"]["telescope_transmission"] = None
    document["sha256"] = model_hash(document)
    return OpticalModel(document)


def analytic_rows(model, offsets, photons: int, seed: int) -> list[dict[str, float]]:
    """Per offset: area ratio with the obscurers over without, and the census identity residual."""
    free = _without_obscurers(model)
    rows = []
    for offset in offsets:
        star = model.star_photons(photons, seed=seed, towards_source=_towards(offset))
        full = model.trace(star, seed=seed, apply_responses=False)
        obscured = _area(full)
        unobscured_trace = free.trace(star, seed=seed, apply_responses=False)
        unobscured = _area(unobscured_trace)
        if not unobscured > 0.0:
            raise StudyError(f"no light reaches the focal surface at offset {offset}")
        hits = model.obscurer_hits(star, seed=seed)
        first_share = float(sum(hits["first"])) / float(hits["unobscured"])
        ratio = obscured / unobscured
        rows.append(
            {
                "off_x_deg": offset[0],
                "off_y_deg": offset[1],
                "angle_deg": math.hypot(*offset),
                "area_ratio": ratio,
                "first_share": first_share,
                "residual": ratio - (1.0 - first_share),
            }
        )
    return rows


def execute_analytic(config: dict[str, Any], output: Path) -> dict[str, float | int]:
    from obdeect.ray_tracing.model import OpticalModel

    options = config["simtools_config"]["configuration"]
    analytic = config["analytic"]
    offsets = _offsets(config["compare"]["offsets_deg"])
    model = OpticalModel.from_simulation_models(
        str(common.models_root()),
        options["telescope"],
        analytic["model_version"],
        configuration_seed=int(analytic["seed"]),
        zenith_deg=float(options["zenith_angle"]),
    )
    rows = analytic_rows(
        model, offsets, int(analytic["photons"]), int(analytic["seed"])
    )
    residual = max(abs(row["residual"]) for row in rows)
    metrics: dict[str, float | int] = {
        "analytic_n_offsets": len(rows),
        "analytic_max_abs_residual": residual,
        "analytic_limit_exceeded": int(residual > ANALYTIC_TOLERANCE),
    }
    details = {
        "telescope": options["telescope"],
        "model_version": analytic["model_version"],
        "model_sha256": model.sha256,
        "photons_per_offset": int(analytic["photons"]),
        "seed": int(analytic["seed"]),
        "tolerance": ANALYTIC_TOLERANCE,
        "rows": rows,
        "definition": (
            "area_ratio = A(all obscurers) / A(obscurer-free copy), geometry only; "
            "first_share = sum of the first-obscurer weights / unobscured weight "
            "(OpticalModel.obscurer_hits)"
        ),
    }
    output.mkdir(parents=True, exist_ok=True)
    common.write_results(metrics, details, root=output)
    return metrics


# --- compare part (container: simtools with both backends) ----------------------------------


def variant_overrides(telescope: str) -> dict[str, dict[str, Any] | None]:
    """simtools overwrite_model_parameters per variant (None: as modelled).

    Baffle r1 = 0 means no baffle (schema of secondary_mirror_baffle). A transmission of 1 with
    method 0 and no amplitude is no shadowing.
    """
    return {
        "base": None,
        "no_baffle": {telescope: {BAFFLE: {"value": [0.0, 0.0, 0.0, 0.0, 0.0]}}},
        "no_structure": {
            telescope: {TRANSMISSION: {"value": [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]}}
        },
    }


def _run_simtools(
    options, overrides, backend: str, output: Path
) -> dict[tuple[float, float], float]:
    """Run simtools-validate-optics once; effective area [m^2] keyed by field offset (x, y) [deg]."""
    run_simtools(
        "simtools-validate-optics", options, backend, output, overrides=overrides
    )
    table = _table(output)
    areas: dict[tuple[float, float], float] = {}
    for x, y, value in zip(
        table["off_x"], table["off_y"], table["eff_area"], strict=True
    ):
        if value is None:
            raise StudyError(f"undefined eff_area at offset ({x}, {y}) in {output}")
        areas[(round(x, 6), round(y, 6))] = float(value)
    return areas


def _ratios(areas: dict[str, dict]) -> dict[str, dict[tuple[float, float], float]]:
    """rho_baffle and rho_structure per field offset of one seed and backend."""
    base = areas["base"]
    return {
        "rho_baffle": {k: base[k] / areas["no_baffle"][k] for k in base},
        "rho_structure": {k: base[k] / areas["no_structure"][k] for k in base},
    }


def compare_rows(
    per_backend: dict[str, list[dict[str, dict[tuple[float, float], float]]]],
    offsets: list[tuple[float, float]],
) -> list[dict[str, Any]]:
    """Per field offset and ratio: seed means of both backends, the gate and its inputs."""
    rows = []
    for offset in offsets:
        row: dict[str, Any] = {
            "off_x_deg": offset[0],
            "off_y_deg": offset[1],
            "angle_deg": math.hypot(*offset),
        }
        for ratio in ("rho_baffle", "rho_structure"):
            stats = {}
            for backend, runs in per_backend.items():
                values = [run[ratio][offset] for run in runs]
                mean, std = common.mean_std(values)
                stats[backend] = {"mean": mean, "std": std}
            obd, sim = stats["obdeect"], stats["sim_telarray"]
            sigma = math.sqrt((obd["std"] ** 2 + sim["std"] ** 2) / 2.0)
            delta = obd["mean"] - sim["mean"]
            row[ratio] = {
                "sim_telarray": sim["mean"],
                "obdeect": obd["mean"],
                "delta": delta,
                "sigma_seed": sigma,
                "threshold": LIMIT,
                "exceeded": abs(delta) > LIMIT,
                "underpowered": SIGMA_FACTOR * sigma > LIMIT,
                "goal_exceeded": abs(delta) > GOAL,
            }
        rows.append(row)
    return rows


def execute_compare(config: dict[str, Any], output: Path) -> dict[str, float | int]:
    options = dict(config["simtools_config"]["configuration"])
    seeds = [int(s) for s in config["compare"]["seeds"]]
    overrides = variant_overrides(options["telescope"])
    per_backend: dict[str, list[dict[str, dict[tuple[float, float], float]]]] = {
        b: [] for b in BACKENDS
    }
    for backend in BACKENDS:
        for seed in seeds:
            areas = {}
            for variant in VARIANTS:
                run_options = {
                    **options,
                    "instrument_seed": seed,
                    "simulation_seed": [seed],
                }
                run_dir = output / backend / f"seed{seed}" / variant
                areas[variant] = _run_simtools(
                    run_options, overrides[variant], backend, run_dir
                )
            per_backend[backend].append(_ratios(areas))
    offsets = sorted(per_backend[BACKENDS[0]][0]["rho_baffle"])
    rows = compare_rows(per_backend, offsets)
    metrics: dict[str, float | int] = {"n_offsets": len(rows)}
    for ratio, label in (("rho_baffle", "baffle"), ("rho_structure", "structure")):
        metrics[f"{label}_max_abs_delta"] = max(
            abs(row[ratio]["delta"]) for row in rows
        )
        metrics[f"{label}_limit_exceeded"] = sum(
            int(row[ratio]["exceeded"]) for row in rows
        )
        metrics[f"{label}_underpowered"] = sum(
            int(row[ratio]["underpowered"]) for row in rows
        )
        metrics[f"{label}_goal_exceeded_offsets"] = sum(
            int(row[ratio]["goal_exceeded"]) for row in rows
        )
    metrics["exceeded_limits"] = (
        metrics["baffle_limit_exceeded"] + metrics["structure_limit_exceeded"]
    )
    details = {
        "telescope": options["telescope"],
        "limit": LIMIT,
        "sigma_factor": SIGMA_FACTOR,
        "seeds": seeds,
        "variants": overrides,
        "rule": "|rho_obdeect - rho_sim_telarray| <= limit per offset; underpowered if 3 * pooled seed sigma > limit",
        "goal": GOAL,
        "rows": rows,
    }
    output.mkdir(parents=True, exist_ok=True)
    common.write_results(metrics, details, root=output)
    return metrics


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("part", choices=("analytic", "compare"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path.cwd())
    args = parser.parse_args(argv)
    try:
        config = _config(args.config)
        run = execute_analytic if args.part == "analytic" else execute_compare
        print(json.dumps(run(config, args.output), indent=2, sort_keys=True))
        return 0
    except (StudyError, OSError, ValueError) as exc:
        print(f"optics obscurer study: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
