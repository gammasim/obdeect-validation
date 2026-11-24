import json
from pathlib import Path

import pytest

from obdeect_validation.adapters import common
from obdeect_validation.adapters import optics_psf_angles as adapter


def test_gate_threshold_is_the_requirement_and_noise_is_underpowered():
    stats = {
        "obdeect": {"mean": 10.0, "std": 0.0},
        "sim_telarray": {"mean": 10.0, "std": 0.0},
    }
    gate = adapter._gate(stats, 0.05)
    assert gate["threshold"] == pytest.approx(0.05)
    assert gate["exceeded"] is False

    stats["obdeect"]["mean"] = 10.2
    assert adapter._gate(stats, 0.05)["exceeded"] is True

    noisy = {
        "obdeect": {"mean": 10.0, "std": 0.1},
        "sim_telarray": {"mean": 10.0, "std": 0.1},
    }
    gate = adapter._gate(noisy, 0.05)
    assert gate["sigma_seed"] == pytest.approx(0.1)
    assert gate["threshold"] == pytest.approx(0.05)
    assert gate["underpowered"] is True


def test_mean_std_and_seed_options():
    assert common.mean_std([2.0, 4.0]) == pytest.approx((3.0, 2**0.5))
    options = adapter._seed_options({"telescope": "LSTN-01"}, 13)
    assert options == {
        "telescope": "LSTN-01",
        "instrument_seed": 13,
        "simulation_seed": [13],
    }


def test_config_checks_application_and_backend(tmp_path: Path):
    section = {
        "seeds": [11, 12, 13],
        "simtools_config": {
            "application": "simtools-derive-incident-angle",
            "configuration": {"ray_tracing_backend": "obdeect"},
        },
    }
    path = tmp_path / "cfg.json"
    path.write_text(json.dumps({"angles": section}))
    with pytest.raises(adapter.StudyError, match="omit ray_tracing_backend"):
        adapter._config(path, "angles")

    section["simtools_config"]["configuration"] = {}
    path.write_text(json.dumps({"psf": section}))
    with pytest.raises(adapter.StudyError, match="must be simtools-validate-optics"):
        adapter._config(path, "psf")


def test_config_needs_three_seeds(tmp_path: Path):
    section = {
        "seeds": [11, 12],
        "simtools_config": {
            "application": "simtools-validate-optics",
            "configuration": {},
        },
    }
    path = tmp_path / "cfg.json"
    path.write_text(json.dumps({"psf": section}))
    with pytest.raises(adapter.StudyError, match="at least three seeds"):
        adapter._config(path, "psf")


def test_histogram_mean_is_weighted_and_rejects_empty(tmp_path: Path):
    import astropy.units as u
    from astropy.table import QTable

    table = QTable()
    table["incidence_angle"] = [10.0, 20.0] * u.deg
    table["fraction"] = [1.0, 3.0]
    path = tmp_path / "h.ecsv"
    table.write(path, format="ascii.ecsv")
    assert adapter._histogram_mean(path) == pytest.approx(17.5)

    table["fraction"] = [0.0, 0.0]
    table.write(path, format="ascii.ecsv", overwrite=True)
    with pytest.raises(adapter.StudyError, match="empty"):
        adapter._histogram_mean(path)
