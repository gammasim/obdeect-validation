"""Measure the LST camera-window effect with obdeect only.

The window parameters come from the simulation-models checkout (LSTN-design).
Model versions without them get them in a temporary copy of the production
table; the checkout is not changed.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

from obdeect_validation.adapters.simtools_validate_optics import _table

WINDOW_PARAMETERS = (
    "camera_window_parameters",
    "camera_window_refractive_index",
    "camera_window_internal_transmittance",
)
WINDOW_VERSION = "1.0.0"
METRICS = ("psf_cm", "eff_area", "eff_flen")


class StudyError(ValueError):
    """Invalid study configuration or optics result."""


def _config(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise StudyError(f"cannot read study config {path}: {exc}") from exc
    if not isinstance(value, dict) or set(value) != {"simtools_config", "study"}:
        raise StudyError("config must contain simtools_config and study")
    base = value["simtools_config"]
    if (
        not isinstance(base, dict)
        or base.get("application") != "simtools-validate-optics"
    ):
        raise StudyError("simtools_config.application must be simtools-validate-optics")
    options = base.get("configuration")
    if not isinstance(options, dict):
        raise StudyError("simtools_config.configuration must be an object")
    if options.get("telescope") != "LSTN-01":
        raise StudyError("the window study is defined for telescope LSTN-01")
    study = value["study"]
    kinds = {"reference", "fixed-window", "focus-scan"}
    if not isinstance(study, dict) or study.get("kind") not in kinds:
        raise StudyError("study.kind must be reference, fixed-window, or focus-scan")
    return value


def _focus_scan(study: dict[str, Any]) -> list[float]:
    scan = study.get("scan", {})
    start, stop, step = (scan.get(key) for key in ("start_cm", "stop_cm", "step_cm"))
    if any(
        isinstance(v, bool) or not isinstance(v, (int, float))
        for v in (start, stop, step)
    ):
        raise StudyError("scan start_cm, stop_cm, and step_cm must be numbers")
    if step <= 0 or start < stop:
        raise StudyError("scan requires start_cm >= stop_cm and step_cm > 0")
    count = round((start - stop) / step)
    if count > 100:
        raise StudyError("focus scan is limited to 101 points")
    values = [round(start - i * step, 8) for i in range(count + 1)]
    if not math.isclose(values[-1], stop, abs_tol=1e-7):
        raise StudyError("scan range must be divisible by step_cm")
    return values


def _focus_override(focus_cm: float) -> dict[str, Any]:
    """simtools overwrite_model_parameters entry for the focus offset."""
    return {"LSTN-design": {"focus_offset": {"value": [focus_cm, 0.0, 0.0, 0.0]}}}


def _model_overlay(models: Path, destination: Path, version: str) -> Path:
    """Copy the models and add the window parameters to the LSTN-design production table."""
    source = (
        models / "simulation-models"
        if (models / "simulation-models").is_dir()
        else models
    )
    manifest_path = source / "productions" / version / "LSTN-design.json"
    if not manifest_path.is_file():
        raise StudyError(
            f"model production {version}/LSTN-design.json is missing under {models}"
        )
    for name in WINDOW_PARAMETERS:
        record = (
            source
            / "model_parameters"
            / "LSTN-design"
            / name
            / f"{name}-{WINDOW_VERSION}.json"
        )
        if not record.is_file():
            raise StudyError(f"window parameter {record} is missing")
    base = destination / "simulation-models" if source != models else destination
    shutil.copytree(source, base, symlinks=True)
    manifest_path = base / "productions" / version / "LSTN-design.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["parameters"]["LSTN-design"].update(
        dict.fromkeys(WINDOW_PARAMETERS, WINDOW_VERSION)
    )
    manifest_path.write_text(json.dumps(manifest, indent=4) + "\n", encoding="utf-8")
    return destination


def _run_one(
    options: dict[str, Any], output: Path, models: Path, focus_cm: float
) -> dict[str, list[float | None]]:
    output.mkdir(parents=True, exist_ok=True)
    overrides = json.dumps(_focus_override(focus_cm), separators=(",", ":"))
    config = {
        "applications": [
            {
                "application": "simtools-validate-optics",
                "configuration": {
                    **options,
                    "ray_tracing_backend": "obdeect",
                    "simulation_models_path": str(models),
                    "output_path": str(output),
                    "overwrite_model_parameters": overrides,
                },
            }
        ]
    }
    config_path = output / "simtools-config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    with (
        (output / "stdout.log").open("wb") as stdout,
        (output / "stderr.log").open("wb") as stderr,
    ):
        try:
            process = subprocess.run(
                ["simtools-validate-optics", "--config", str(config_path)],
                stdout=stdout,
                stderr=stderr,
                check=False,
            )
        except OSError as exc:
            raise StudyError(f"cannot run simtools-validate-optics: {exc}") from exc
    if process.returncode:
        raise StudyError(
            f"simtools-validate-optics exited {process.returncode}; see {output}"
        )
    return _table(output)


def _on_axis(table: dict[str, list[float | None]]) -> float:
    matches = [
        value
        for x, y, value in zip(
            table["off_x"], table["off_y"], table["psf_cm"], strict=True
        )
        if abs(x) < 1e-10 and abs(y) < 1e-10
    ]
    if len(matches) != 1 or matches[0] is None:
        raise StudyError("expected exactly one finite on-axis D80 value")
    return float(matches[0])


def _relative_delta(
    metric: str, baseline: float | None, window: float | None
) -> float | None:
    if baseline is None and window is None:
        return None
    if baseline is None or window is None or baseline == 0:
        raise StudyError(f"{metric} rows do not align or contain undefined values")
    return window / baseline - 1


def _paired_metrics(reference, candidate) -> list[dict[str, float | None]]:
    """Window-on vs window-off metrics, row by row at equal field offsets."""
    by_offset = {
        (round(x, 10), round(y, 10)): index
        for index, (x, y) in enumerate(
            zip(reference["off_x"], reference["off_y"], strict=True)
        )
    }
    if len(by_offset) != len(candidate["off_x"]):
        raise StudyError("optics rows do not align in length")
    rows = []
    for index, (x, y) in enumerate(
        zip(candidate["off_x"], candidate["off_y"], strict=True)
    ):
        baseline_index = by_offset.get((round(x, 10), round(y, 10)))
        if baseline_index is None:
            raise StudyError("optics rows do not align in field offset")
        row: dict[str, float | None] = {
            "off_x_deg": x,
            "off_y_deg": y,
            "angle_deg": math.hypot(x, y),
        }
        for metric in METRICS:
            baseline, window = (
                reference[metric][baseline_index],
                candidate[metric][index],
            )
            row[f"reference_{metric}"] = baseline
            row[f"window_{metric}"] = window
            row[f"{metric}_relative_delta"] = _relative_delta(metric, baseline, window)
        rows.append(row)
    return rows


def _max_abs(rows: list[dict[str, Any]], key: str) -> float:
    return max((abs(row[key]) for row in rows if row[key] is not None), default=0.0)


def _reference(options, output, models, focus):
    table = _run_one(options, output / "reference", models, focus)
    metrics = {"n_offsets": len(table["off_x"]), "focus_offset_cm": focus}
    metrics["on_axis_d80_cm"] = _on_axis(table)
    return metrics, {"study": "reference", "focus_offset_cm": focus, "table": table}


def _fixed_window(options, output, models, window_models, focus):
    reference = _run_one(options, output / "window-off", models, focus)
    window = _run_one(options, output / "window-on", window_models, focus)
    rows = _paired_metrics(reference, window)
    metrics = {"n_offsets": len(window["off_x"]), "focus_offset_cm": focus}
    metrics["on_axis_d80_cm"] = _on_axis(window)
    for metric in METRICS:
        metrics[f"max_abs_rel_{metric}"] = _max_abs(rows, f"{metric}_relative_delta")
    details = {
        "study": "fixed-window",
        "comparison": "window on vs off at the same focus offset",
        "focus_offset_cm": focus,
        "off_axis_metrics": rows,
    }
    return metrics, details


def _focus_scan_study(options, output, models, window_models, focus, study):
    reference = _run_one(options, output / "window-off", models, focus)
    scan = []
    for candidate in _focus_scan(study):
        table = _run_one(
            options,
            output / f"window-on-focus-{candidate:.2f}",
            window_models,
            candidate,
        )
        scan.append(
            {
                "focus_offset_cm": candidate,
                "on_axis_d80_cm": _on_axis(table),
                "table": table,
            }
        )
    best = min(
        scan,
        key=lambda row: (row["on_axis_d80_cm"], abs(row["focus_offset_cm"] - focus)),
    )
    off_axis = _paired_metrics(reference, best["table"])
    at_1_to_4 = [row for row in off_axis if 1 - 1e-9 <= row["angle_deg"] <= 4 + 1e-9]
    metrics = {
        "n_focus_points": len(scan),
        "best_focus_offset_cm": best["focus_offset_cm"],
        "best_on_axis_d80_cm": best["on_axis_d80_cm"],
        "n_off_axis_1_to_4_deg": len(at_1_to_4),
        "max_abs_rel_d80_1_to_4_deg": _max_abs(at_1_to_4, "psf_cm_relative_delta"),
        "n_offsets_over_2pct": sum(
            abs(row["psf_cm_relative_delta"]) > 0.02 for row in at_1_to_4
        ),
    }
    details = {
        "study": "focus-scan",
        "comparison": f"best window-on focus vs window-off reference at {focus} cm",
        "focus_scan": [{k: v for k, v in row.items() if k != "table"} for row in scan],
        "best_focus_offset_cm": best["focus_offset_cm"],
        "best_scan_table": best["table"],
        "best_focus_off_axis_d80": at_1_to_4,
    }
    return metrics, details


def execute(config_path: Path, output: Path) -> dict[str, float | int]:
    config = _config(config_path)
    options = config["simtools_config"]["configuration"]
    study = config["study"]
    output.mkdir(parents=True, exist_ok=True)
    models = Path(options["simulation_models_path"])
    focus = float(study.get("focus_offset_cm", 6.55))
    with tempfile.TemporaryDirectory(prefix="lst-window-models-") as temp:
        if study["kind"] == "reference":
            metrics, details = _reference(options, output, models, focus)
        else:
            window_models = _model_overlay(
                models, Path(temp) / "models", str(options["model_version"])
            )
            if study["kind"] == "fixed-window":
                metrics, details = _fixed_window(
                    options, output, models, window_models, focus
                )
            else:
                metrics, details = _focus_scan_study(
                    options, output, models, window_models, focus, study
                )
    details["backend"] = "obdeect"
    (output / "metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "details.json").write_text(
        json.dumps(details, indent=2, sort_keys=True) + "\n", encoding="utf-8"
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
    except (StudyError, OSError, ValueError) as exc:
        print(f"LST window optics study: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
