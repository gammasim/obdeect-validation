"""Run simtools optics validation with both ray-tracing backends and compare per offset."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

METRICS = ("psf_cm", "eff_area", "eff_flen")
LIMIT_STATUS = (
    "requirement set by the project (2026-11-24): goal 1 %, acceptable at most 2 %"
)
LIMIT_BASIS = (
    "D80 and effective area: goal 1 %, acceptable at most 2 % per offset (limits = acceptance, "
    "goals reported, not gated). Effective focal length: 0.2 % (stricter than the requirement; "
    "measured agreement 0.1 %, phase1.md). The limit is not widened by the statistical error: an "
    "offset whose 3-sigma error exceeds the limit is reported as underpowered and fails"
)
STAT_BASIS = {
    "psf_cm": "not available: simtools output has no D80 uncertainty",
    "eff_area": "not available: simtools output has no effective-area uncertainty",
    "eff_flen": "eff_flen_err (centroid-propagated) per backend, combined in quadrature",
}
SIGMA_FACTOR = 3.0


class AdapterError(ValueError):
    """Invalid optics adapter configuration or result."""


def _config(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AdapterError(f"cannot read config {path}: {exc}") from exc
    if not isinstance(value, dict) or not (
        {"simtools_config", "limits"}
        <= set(value)
        <= {"simtools_config", "limits", "goals"}
    ):
        raise AdapterError(
            "config must contain simtools_config, limits and optionally goals"
        )
    app = value["simtools_config"]
    if (
        not isinstance(app, dict)
        or app.get("application") != "simtools-validate-optics"
    ):
        raise AdapterError(
            "simtools_config.application must be simtools-validate-optics"
        )
    options = app.get("configuration")
    if not isinstance(options, dict) or any(not isinstance(k, str) for k in options):
        raise AdapterError("simtools_config.configuration must be an object")
    if "ray_tracing_backend" in options or "output_path" in options:
        raise AdapterError(
            "configuration must omit ray_tracing_backend and output_path; adapter sets them per backend"
        )
    limits = value["limits"]
    if not isinstance(limits, dict) or set(limits) != set(METRICS):
        raise AdapterError(
            "limits must define psf_cm, eff_area, and eff_flen relative limits"
        )
    for key, limit in limits.items():
        if (
            isinstance(limit, bool)
            or not isinstance(limit, (int, float))
            or not math.isfinite(limit)
            or limit < 0
        ):
            raise AdapterError(f"limits.{key} must be a finite non-negative number")
    return value


def _column(
    table, name: str, unit, path: Path, offsets: list[float]
) -> list[float | None]:
    """Return column values in ``unit``; NaN is allowed only on axis (eff_flen and its error)."""
    import astropy.units as u

    column = table[name]
    if getattr(column, "unit", None) is None:
        if name not in ("eff_flen", "eff_flen_err"):
            raise AdapterError(f"{path} column {name} has no unit")
        values = [
            float(v) for v in column
        ]  # simtools stores eff_flen as plain cm (ray_tracing.py)
    else:
        try:
            values = [float(v) for v in column.to_value(unit)]
        except u.UnitConversionError as exc:
            raise AdapterError(
                f"{path} column {name} cannot be converted to {unit}"
            ) from exc
    normalized: list[float | None] = []
    for index, value in enumerate(values):
        if math.isfinite(value):
            normalized.append(value)
        elif (
            name in ("eff_flen", "eff_flen_err")
            and math.isnan(value)
            and offsets[index] <= 1e-12
        ):
            normalized.append(None)
        else:
            raise AdapterError(
                f"{path} has undefined {name} at row {index} (offset {offsets[index]:.3g} deg)"
            )
    return normalized


def _table(path: Path) -> dict[str, list[float | None]]:
    try:
        import astropy.units as u
        from astropy.table import QTable
    except ImportError as exc:
        raise AdapterError(
            "astropy is required to compare simtools ECSV output"
        ) from exc
    files = sorted(path.glob("results/*.ecsv"))
    if len(files) != 1:
        raise AdapterError(
            f"expected exactly one results/*.ecsv under {path}; found {len(files)}"
        )
    table = QTable.read(files[0], format="ascii.ecsv")
    required = {"off_x", "off_y", "psf_cm", "eff_area", "eff_flen"}
    missing = required - set(table.colnames)
    if missing:
        raise AdapterError(
            f"{files[0]} is missing columns: {', '.join(sorted(missing))}"
        )
    if len(table) == 0:
        raise AdapterError(f"{files[0]} has no rows")
    off_x = _column(table, "off_x", u.deg, files[0], [0.0] * len(table))
    off_y = _column(table, "off_y", u.deg, files[0], [0.0] * len(table))
    offsets = [math.hypot(x, y) for x, y in zip(off_x, off_y, strict=True)]
    columns: dict[str, list[float | None]] = {"off_x": off_x, "off_y": off_y}
    units = {"psf_cm": u.cm, "eff_area": u.m**2, "eff_flen": u.cm}
    for name, unit in units.items():
        columns[name] = _column(table, name, unit, files[0], offsets)
    columns["eff_flen_err"] = (
        _column(table, "eff_flen_err", u.cm, files[0], offsets)
        if "eff_flen_err" in table.colnames
        else [None] * len(table)
    )
    return columns


def _sorted(table: dict[str, list[float | None]]) -> dict[str, list[float | None]]:
    order = sorted(
        range(len(table["off_x"])), key=lambda i: (table["off_x"][i], table["off_y"][i])
    )
    return {key: [values[i] for i in order] for key, values in table.items()}


def _rows(sim: dict, obd: dict, limits: dict[str, float]) -> list[dict[str, Any]]:
    """Per aligned offset: signed relative difference (obd/sim - 1), statistical error and flags."""
    rows = []
    for index in range(len(sim["off_x"])):
        off_x, off_y = sim["off_x"][index], sim["off_y"][index]
        row: dict[str, Any] = {
            "off_x_deg": off_x,
            "off_y_deg": off_y,
            "angle_deg": math.hypot(off_x, off_y),
        }
        for metric in METRICS:
            left, right = sim[metric][index], obd[metric][index]
            if left is None and right is None:
                row[metric] = None
                continue
            if left is None or right is None:
                raise AdapterError(
                    f"{metric} undefined for only one backend at aligned row {index}"
                )
            if left == 0:
                raise AdapterError(
                    f"reference {metric} is zero at aligned row {index}; relative difference is undefined"
                )
            delta = right / left - 1.0
            se = None
            if metric == "eff_flen":
                err_sim, err_obd = (
                    sim["eff_flen_err"][index],
                    obd["eff_flen_err"][index],
                )
                se = math.hypot(err_obd, err_sim * right / left) / abs(left)
            threshold = limits[metric]
            row[metric] = {
                "sim": left,
                "obd": right,
                "delta": delta,
                "rel_se": se,
                "threshold": threshold,
                "exceeded": abs(delta) > threshold,
                # The comparison cannot resolve the limit: its 3-sigma error is larger.
                "underpowered": SIGMA_FACTOR * (se or 0.0) > threshold,
            }
        rows.append(row)
    return rows


def _summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for metric in METRICS:
        values = [row[metric] for row in rows if row[metric] is not None]
        worst = max(values, key=lambda v: abs(v["delta"]))
        summary[metric] = {
            "max_abs_delta": abs(worst["delta"]),
            "signed_delta_at_max": worst["delta"],
            "exceeded_offsets": [
                {
                    "angle_deg": row["angle_deg"],
                    "off_x_deg": row["off_x_deg"],
                    "off_y_deg": row["off_y_deg"],
                    "delta": row[metric]["delta"],
                    "threshold": row[metric]["threshold"],
                }
                for row in rows
                if row[metric] is not None and row[metric]["exceeded"]
            ],
        }
    return summary


VALUE_LABELS = {
    "psf_cm": "PSF D80 [cm]",
    "eff_area": "effective mirror area [m²]",
    "eff_flen": "effective focal length [cm]",
}
# Offset directions of the star: along x (E/W, filled markers) or along y (N/S, open markers).
DIRECTIONS = (
    (True, "o", "-", "offset along x (E/W)"),
    (False, "s", "--", "offset along y (N/S)"),
)


def plot_values(rows: list[dict[str, Any]], path: Path, label: str) -> None:
    """obdeect and sim_telarray values per metric against off-axis angle, one panel per metric.

    Colour is the program, marker and line the offset direction of the star: the optics are
    not rotationally symmetric, so the two directions give different values off axis.
    """
    from obdeect_validation.adapters import plotting

    fig, axes = plotting.figure(len(METRICS), 1, (6, 8))
    for axis, metric in zip(axes[:, 0], METRICS, strict=True):
        points = [row for row in rows if row[metric] is not None]
        for program, key, color in (
            ("obdeect", "obd", plotting.NATIVE),
            ("sim_telarray", "sim", plotting.REFERENCE),
        ):
            for east_west, marker, style, direction in DIRECTIONS:
                selected = sorted(
                    (r for r in points if (abs(r["off_y_deg"]) < 1e-9) == east_west),
                    key=lambda row: row["angle_deg"],
                )
                axis.plot(
                    [row["angle_deg"] for row in selected],
                    [row[metric][key] for row in selected],
                    marker=marker,
                    linestyle=style,
                    color=color,
                    markerfacecolor=color if east_west else "white",
                    label=f"{program}, {direction}",
                )
        axis.set_ylabel(VALUE_LABELS[metric])
    axes[0, 0].legend(fontsize="small")
    axes[-1, 0].set_xlabel("off-axis angle [deg]")
    fig.suptitle(f"{label}: obdeect and sim_telarray", fontsize="small")
    plotting.save(fig, path)


def _plot(
    rows: list[dict[str, Any]], limits: dict[str, float], path: Path, label: str
) -> None:
    """Relative obd/sim difference per metric against off-axis angle.

    Colour is the telescope (same colour in every figure), marker and line the offset
    direction; the shaded band is the limit.
    """
    from obdeect_validation.adapters import plotting

    color = plotting.case_color(label)
    fig, axes = plotting.figure(len(METRICS), 1, (6, 8))
    for axis, metric in zip(axes[:, 0], METRICS, strict=True):
        axis.axhspan(
            -100 * limits[metric],
            100 * limits[metric],
            color=plotting.LIMIT,
            alpha=0.18,
            label=f"limit ±{100 * limits[metric]:g} %",
        )
        axis.axhline(0, color=plotting.ANALYTIC, linewidth=0.6)
        points = [row for row in rows if row[metric] is not None]
        for east_west, marker, style, name in DIRECTIONS:
            selected = sorted(
                (row for row in points if (abs(row["off_y_deg"]) < 1e-9) == east_west),
                key=lambda row: row["angle_deg"],
            )
            axis.plot(
                [row["angle_deg"] for row in selected],
                [100 * row[metric]["delta"] for row in selected],
                marker=marker,
                linestyle=style,
                color=color,
                markerfacecolor=color if east_west else "white",
                label=name,
            )
        axis.set_ylabel(f"{metric} obd/sim - 1 [%]")
        axis.legend(ncols=3)
    axes[-1, 0].set_xlabel("off-axis angle [deg]")
    fig.suptitle(f"{label}: shaded band = limit", fontsize="small")
    plotting.save(fig, path)


def _versions() -> dict[str, str | None]:
    from importlib import metadata

    out: dict[str, str | None] = {}
    for name in ("gammasimtools", "obdeect"):
        try:
            out[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            out[name] = None
    return out


def execute(config_path: Path, output: Path) -> dict[str, float | int]:
    config = _config(config_path)
    limits = config["limits"]
    output.mkdir(parents=True, exist_ok=True)
    result_tables = {}
    backends: dict[str, Any] = {}
    for backend in ("sim_telarray", "obdeect"):
        backend_dir = output / backend
        backend_dir.mkdir(parents=True, exist_ok=True)
        app_config = {
            "application": "simtools-validate-optics",
            "configuration": {
                **config["simtools_config"]["configuration"],
                "ray_tracing_backend": backend,
                "output_path": str(backend_dir),
            },
        }
        config_file = backend_dir / "simtools-config.json"
        config_file.write_text(
            json.dumps({"applications": [app_config]}, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        command = ["simtools-validate-optics", "--config", str(config_file)]
        try:
            with (
                (backend_dir / "stdout.log").open("wb") as stdout,
                (backend_dir / "stderr.log").open("wb") as stderr,
            ):
                proc = subprocess.run(
                    command, stdout=stdout, stderr=stderr, check=False
                )
        except OSError as exc:
            raise AdapterError(f"cannot run {command[0]}: {exc}") from exc
        if proc.returncode:
            raise AdapterError(
                f"{command[0]} ({backend}) exited {proc.returncode}; see {backend_dir}"
            )
        backends[backend] = {"dir": str(backend_dir), "returncode": proc.returncode}
        result_tables[backend] = _table(backend_dir)

    sim = _sorted(result_tables["sim_telarray"])
    obd = _sorted(result_tables["obdeect"])
    if len(sim["off_x"]) != len(obd["off_x"]):
        raise AdapterError("backend offset rows do not align in length")
    for key in ("off_x", "off_y"):
        if any(
            not math.isclose(a, b, rel_tol=0, abs_tol=1e-10)
            for a, b in zip(sim[key], obd[key], strict=True)
        ):
            raise AdapterError("backend offset rows do not align")

    rows = _rows(sim, obd, limits)
    summary = _summary(rows)
    metrics: dict[str, float | int] = {"n_offsets": len(rows)}
    for metric in METRICS:
        metrics[f"max_rel_{metric}"] = summary[metric]["max_abs_delta"]
        metrics[f"{metric}_limit_exceeded"] = int(
            bool(summary[metric]["exceeded_offsets"])
        )
    metrics["exceeded_limits"] = sum(
        metrics[f"{metric}_limit_exceeded"] for metric in METRICS
    )
    goals = config.get("goals", {})
    for metric in METRICS:
        values = [row[metric] for row in rows if row[metric] is not None]
        metrics[f"{metric}_underpowered"] = sum(int(v["underpowered"]) for v in values)
        if metric in goals:
            # Reported only: offsets above the goal (below the limit they still pass).
            metrics[f"{metric}_goal_exceeded_offsets"] = sum(
                int(abs(v["delta"]) > goals[metric]) for v in values
            )

    details = {
        "status": LIMIT_STATUS,
        "limits_basis": LIMIT_BASIS,
        "limits": limits,
        "goals": config.get("goals", {}),
        "statistical_uncertainty": STAT_BASIS,
        "exceeded_rule": "abs(obd/sim - 1) > limit per offset",
        "underpowered_rule": f"{SIGMA_FACTOR:g} * rel_se > limit per offset",
        "configuration": config["simtools_config"]["configuration"],
        "versions": _versions(),
        "backends": backends,
        "summary": summary,
        "rows": rows,
    }
    (output / "details.json").write_text(
        json.dumps(details, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    telescope = config["simtools_config"]["configuration"].get("telescope", "optics")
    _plot(rows, limits, output / "optics_relative_difference.png", telescope)
    plot_values(rows, output / "optics_values.png", telescope)
    (output / "metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metrics


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path.cwd())
    args = parser.parse_args(argv)
    try:
        print(json.dumps(execute(args.config, args.output), indent=2, sort_keys=True))
        return 0
    except (AdapterError, OSError, ValueError) as exc:
        print(f"simtools optics adapter: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
