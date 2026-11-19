import json
from pathlib import Path

import pytest

from obdeect_validation.adapters import optics_window_study as study


def _models(root: Path, *, with_window: bool = True) -> Path:
    base = root / "simulation-models"
    production = base / "productions" / "7.0.0"
    production.mkdir(parents=True)
    (production / "LSTN-design.json").write_text(
        json.dumps({"parameters": {"LSTN-design": {"focus_offset": "1.0.0"}}})
    )
    parameters = base / "model_parameters" / "LSTN-design"
    names = study.WINDOW_PARAMETERS if with_window else study.WINDOW_PARAMETERS[:1]
    for name in names:
        (parameters / name).mkdir(parents=True)
        (parameters / name / f"{name}-1.0.0.json").write_text("{}")
    return root


def test_focus_scan_includes_both_ends_and_rejects_bad_step():
    assert study._focus_scan(
        {"scan": {"start_cm": 6.95, "stop_cm": 6.15, "step_cm": 0.05}}
    ) == pytest.approx([6.95 - 0.05 * i for i in range(17)])
    with pytest.raises(study.StudyError, match="divisible"):
        study._focus_scan({"scan": {"start_cm": 6.5, "stop_cm": 6.0, "step_cm": 0.3}})


def test_focus_override_uses_simtools_overwrite_format():
    assert study._focus_override(6.3) == {
        "LSTN-design": {"focus_offset": {"value": [6.3, 0.0, 0.0, 0.0]}}
    }


def test_model_overlay_adds_window_to_production_without_touching_source(
    tmp_path: Path,
):
    models = _models(tmp_path / "models")
    overlay = study._model_overlay(models, tmp_path / "overlay", "7.0.0")

    manifest = json.loads(
        (overlay / "simulation-models/productions/7.0.0/LSTN-design.json").read_text()
    )
    for name in study.WINDOW_PARAMETERS:
        assert manifest["parameters"]["LSTN-design"][name] == "1.0.0"
    source = models / "simulation-models/productions/7.0.0/LSTN-design.json"
    assert "camera_window" not in source.read_text()


def test_model_overlay_requires_window_parameter_files(tmp_path: Path):
    models = _models(tmp_path / "models", with_window=False)
    with pytest.raises(study.StudyError, match="camera_window_refractive_index"):
        study._model_overlay(models, tmp_path / "overlay", "7.0.0")


def test_paired_metrics_reports_relative_change():
    reference = {
        "off_x": [0.0, 1.0],
        "off_y": [0.0, 0.0],
        "psf_cm": [3.0, 4.0],
        "eff_area": [300.0, 300.0],
        "eff_flen": [2800.0, 2800.0],
    }
    window = {**reference, "psf_cm": [3.3, 4.0]}
    rows = study._paired_metrics(reference, window)
    assert rows[0]["psf_cm_relative_delta"] == pytest.approx(0.1)
    assert rows[1]["psf_cm_relative_delta"] == pytest.approx(0.0)
    assert rows[1]["angle_deg"] == pytest.approx(1.0)


def test_study_config_rejects_other_telescopes(tmp_path: Path):
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "simtools_config": {
                    "application": "simtools-validate-optics",
                    "configuration": {"telescope": "LSTS-01"},
                },
                "study": {"kind": "reference"},
            }
        )
    )
    with pytest.raises(study.StudyError, match="LSTN-01"):
        study._config(config)
