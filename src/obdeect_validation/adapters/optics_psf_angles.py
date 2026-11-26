"""On-axis cumulative PSF (LSTN-01) and mean incidence angles (SSTS-01): simtools with both backends.

Both parts run the same simtools applications with ``ray_tracing_backend`` set to sim_telarray and
to obdeect, on three seeds (instrument and simulation seed equal), and compare the seed means.

``psf``: ``simtools-validate-optics --type cumulative_psf`` with the measured LST PSF curve
(``PSFcurve_data_v2.ecsv``, simtools integration resource; zenith 20 deg). Quantity: on-axis D80
from the ray-tracing results table. The simtools application writes the cumulative curve only as a
figure, so the curve itself is not compared here (the D80 is its 80 % point).

``angles``: ``simtools-derive-incident-angle`` for SSTS-01 at 1 deg off axis, 10^6 photons per run.
Quantities: mean incidence angle on the focal surface, on the primary and on the secondary mirror
(photons that complete the optical path, unweighted; mean of the histogram that simtools writes:
1000 bins for the focal surface, 100 bins for the mirrors).

Gates. Let m_b be the mean over the seeds of backend b and s_b the sample standard deviation over
the seeds (the seed-to-seed spread of one run). Pooled spread sigma = sqrt((s_obd^2 + s_sim^2)/2).
Gate: |m_obd - m_sim| <= requirement; the requirement is the optics requirement set by the project
for the PSF (D80 goal 1 %, acceptable at most 2 %, relative) and 0.05 deg for the angles (the
agreement stated in docs/ray-tracing/VALIDATION.md). The limit is not widened by the seeds: a case
with 3 sigma above the requirement is reported as underpowered and fails (more seeds needed).

The sim_telarray side needs the container; the obdeect side runs wherever simtools runs.
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

BACKENDS = ("sim_telarray", "obdeect")
PSF_LIMIT = 0.02  # optics requirement: acceptable at most 2 % (relative D80)
PSF_GOAL = 0.01  # optics requirement: goal 1 % (reported, not gated)
ANGLE_LIMIT_DEG = 0.05
SIGMA_FACTOR = 3.0
MIRRORS = ("focal", "primary", "secondary")


StudyError = AdapterError


def _config(path: Path, part: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise StudyError(f"cannot read config {path}: {exc}") from exc
    if not isinstance(value, dict) or part not in value:
        raise StudyError(f"config must contain a '{part}' section")
    section = value[part]
    base = section.get("simtools_config")
    expected = {
        "psf": "simtools-validate-optics",
        "angles": "simtools-derive-incident-angle",
    }[part]
    if not isinstance(base, dict) or base.get("application") != expected:
        raise StudyError(f"{part}.simtools_config.application must be {expected}")
    options = base.get("configuration")
    if not isinstance(options, dict) or "ray_tracing_backend" in options:
        raise StudyError(
            f"{part}.simtools_config.configuration must omit ray_tracing_backend"
        )
    seeds = section.get("seeds")
    if not isinstance(seeds, list) or len(seeds) < 3:
        raise StudyError(f"{part}.seeds needs at least three seeds")
    return section


def _seed_options(options: dict[str, Any], seed: int) -> dict[str, Any]:
    return {**options, "instrument_seed": seed, "simulation_seed": [seed]}


def _gate(stats: dict[str, dict[str, float]], requirement: float) -> dict[str, Any]:
    obd, sim = stats["obdeect"], stats["sim_telarray"]
    sigma = math.sqrt((obd["std"] ** 2 + sim["std"] ** 2) / 2.0)
    delta = obd["mean"] - sim["mean"]
    return {
        "sim_telarray": sim["mean"],
        "obdeect": obd["mean"],
        "delta": delta,
        "sigma_seed": sigma,
        "threshold": requirement,
        "exceeded": abs(delta) > requirement,
        "underpowered": SIGMA_FACTOR * sigma > requirement,
    }


# --- psf -------------------------------------------------------------------------------------


def _on_axis_d80_cm(output: Path) -> float:
    table = _table(output)
    matches = [
        value
        for x, y, value in zip(
            table["off_x"], table["off_y"], table["psf_cm"], strict=True
        )
        if abs(x) < 1e-10 and abs(y) < 1e-10
    ]
    if len(matches) != 1 or matches[0] is None:
        raise StudyError(f"expected exactly one finite on-axis psf_cm in {output}")
    return float(matches[0])


def execute_psf(section: dict[str, Any], output: Path) -> dict[str, float | int]:
    options = section["simtools_config"]["configuration"]
    values: dict[str, list[float]] = {}
    for backend in BACKENDS:
        values[backend] = []
        for seed in section["seeds"]:
            run = output / backend / f"seed{seed}"
            run_simtools(
                "simtools-validate-optics",
                _seed_options(options, int(seed)),
                backend,
                run,
            )
            values[backend].append(_on_axis_d80_cm(run))
    stats = {}
    for backend in BACKENDS:
        mean, std = common.mean_std(values[backend])
        stats[backend] = {"mean": mean, "std": std}
    # Relative gate: the spreads are taken on D80 relative to the seed means.
    relative = {
        b: {"mean": 1.0, "std": stats[b]["std"] / stats[b]["mean"]} for b in BACKENDS
    }
    gate = _gate(relative, PSF_LIMIT)
    gate["d80_cm"] = {b: stats[b]["mean"] for b in BACKENDS}
    gate["d80_cm_seeds"] = values
    gate["delta"] = stats["obdeect"]["mean"] / stats["sim_telarray"]["mean"] - 1.0
    metrics = {
        "psf_max_abs_rel_delta": abs(gate["delta"]),
        "psf_limit_exceeded": int(gate["exceeded"]),
        "psf_underpowered": int(gate["underpowered"]),
        # Reported only: the D80 goal of the optics requirement is 1 %.
        "psf_goal_exceeded": int(abs(gate["delta"]) > PSF_GOAL),
    }
    details = {"requirement": PSF_LIMIT, "sigma_factor": SIGMA_FACTOR, "gate": gate}
    output.mkdir(parents=True, exist_ok=True)
    common.write_results(metrics, details, root=output)
    return metrics


# --- angles ----------------------------------------------------------------------------------


def _histogram_mean(path: Path) -> float:
    from astropy.table import QTable

    table = QTable.read(path, format="ascii.ecsv")
    angle = table["incidence_angle"].to_value("deg")
    weight = [float(w) for w in table["fraction"]]
    total = sum(weight)
    if not total > 0:
        raise StudyError(f"empty incidence-angle histogram {path}")
    return sum(a * w for a, w in zip(angle, weight, strict=True)) / total


def _angle_means(run: Path, telescope: str) -> dict[str, float]:
    focal = sorted(run.glob("incident_angles/*.ecsv"))
    if len(focal) != 1:
        raise StudyError(f"expected one incident_angles/*.ecsv under {run}")
    means = {"focal": _histogram_mean(focal[0])}
    for mirror in ("primary", "secondary"):
        path = run / telescope / f"{mirror}_mirror_incidence_angle.ecsv"
        if not path.is_file():
            raise StudyError(f"missing {path}")
        means[mirror] = _histogram_mean(path)
    return means


def execute_angles(section: dict[str, Any], output: Path) -> dict[str, float | int]:
    options = section["simtools_config"]["configuration"]
    telescope = options["telescope"]
    values: dict[str, dict[str, list[float]]] = {
        b: {m: [] for m in MIRRORS} for b in BACKENDS
    }
    for backend in BACKENDS:
        for seed in section["seeds"]:
            run = output / backend / f"seed{seed}"
            run_simtools(
                "simtools-derive-incident-angle",
                _seed_options(options, int(seed)),
                backend,
                run,
            )
            for mirror, mean in _angle_means(run, telescope).items():
                values[backend][mirror].append(mean)
    gates = {}
    for mirror in MIRRORS:
        stats = {}
        for backend in BACKENDS:
            mean, std = common.mean_std(values[backend][mirror])
            stats[backend] = {"mean": mean, "std": std}
        gates[mirror] = _gate(stats, ANGLE_LIMIT_DEG)
    metrics: dict[str, float | int] = {}
    for mirror in MIRRORS:
        metrics[f"angle_{mirror}_max_abs_delta_deg"] = abs(gates[mirror]["delta"])
        metrics[f"angle_{mirror}_limit_exceeded"] = int(gates[mirror]["exceeded"])
        metrics[f"angle_{mirror}_underpowered"] = int(gates[mirror]["underpowered"])
    metrics["angle_exceeded_limits"] = sum(
        metrics[f"angle_{m}_limit_exceeded"] for m in MIRRORS
    )
    details = {
        "requirement_deg": ANGLE_LIMIT_DEG,
        "sigma_factor": SIGMA_FACTOR,
        "seeds": [int(s) for s in section["seeds"]],
        "values_deg": values,
        "gates": gates,
    }
    output.mkdir(parents=True, exist_ok=True)
    common.write_results(metrics, details, root=output)
    return metrics


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("part", choices=("psf", "angles"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path.cwd())
    args = parser.parse_args(argv)
    try:
        section = _config(args.config, args.part)
        run = execute_psf if args.part == "psf" else execute_angles
        print(json.dumps(run(section, args.output), indent=2, sort_keys=True))
        return 0
    except (StudyError, OSError, ValueError) as exc:
        print(f"optics PSF and angle study: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
