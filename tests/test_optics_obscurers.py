import json
from pathlib import Path

import pytest

from obdeect_validation.adapters import common
from obdeect_validation.adapters import optics_obscurers as adapter

OFFSETS = [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0)]


def runs(values_per_seed):
    """Per seed: rho dicts keyed by offset, from {offset: value} for rho_baffle and rho_structure."""
    return [
        {"rho_baffle": dict(b), "rho_structure": dict(s)} for b, s in values_per_seed
    ]


def test_variant_overrides_set_baffle_to_no_baffle_and_transmission_to_identity():
    overrides = adapter.variant_overrides("SCTS-design")
    assert overrides["base"] is None
    assert (
        overrides["no_baffle"]["SCTS-design"]["secondary_mirror_baffle"]["value"]
        == [0.0] * 5
    )
    transmission = overrides["no_structure"]["SCTS-design"]["telescope_transmission"][
        "value"
    ]
    assert transmission == [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]


def test_ratios_divide_base_by_each_variant_per_offset():
    areas = {
        "base": {(0.0, 0.0): 40.0, (0.0, 1.0): 30.0},
        "no_baffle": {(0.0, 0.0): 50.0, (0.0, 1.0): 35.0},
        "no_structure": {(0.0, 0.0): 80.0, (0.0, 1.0): 60.0},
    }
    ratios = adapter._ratios(areas)
    assert ratios["rho_baffle"] == pytest.approx({(0.0, 0.0): 0.8, (0.0, 1.0): 30 / 35})
    assert ratios["rho_structure"] == pytest.approx({(0.0, 0.0): 0.5, (0.0, 1.0): 0.5})


def test_mean_std_uses_sample_standard_deviation():
    mean, std = common.mean_std([1.0, 2.0, 3.0])
    assert mean == pytest.approx(2.0)
    assert std == pytest.approx(1.0)
    assert common.mean_std([5.0]) == (5.0, 0.0)


def test_compare_gate_uses_limit_or_three_pooled_seed_sigmas():
    # Identical seeds: no spread; the 2 % acceptance limit is the threshold.
    same = [runs([({o: 0.95 for o in OFFSETS}, {o: 0.9 for o in OFFSETS})] * 3)] * 2
    rows = adapter.compare_rows({"sim_telarray": same[0], "obdeect": same[1]}, OFFSETS)
    assert all(not row["rho_baffle"]["exceeded"] for row in rows)
    assert rows[0]["rho_baffle"]["threshold"] == pytest.approx(adapter.LIMIT)

    # A 2 % offset between the backend means exceeds the 1 % limit at every offset.
    sim = [runs([({o: 0.95 for o in OFFSETS}, {o: 0.9 for o in OFFSETS})] * 3)]
    obd = [runs([({o: 0.97 for o in OFFSETS}, {o: 0.9 for o in OFFSETS})] * 3)]
    rows = adapter.compare_rows({"sim_telarray": sim[0], "obdeect": obd[0]}, OFFSETS)
    assert all(row["rho_baffle"]["exceeded"] for row in rows)
    assert all(not row["rho_structure"]["exceeded"] for row in rows)


def test_compare_limit_is_fixed_and_noisy_seeds_are_underpowered():
    noisy_sim = runs(
        [
            ({o: 0.95 + d for o in OFFSETS}, {o: 0.9 for o in OFFSETS})
            for d in (-0.02, 0.0, 0.02)
        ]
    )
    obd = runs([({o: 0.97 for o in OFFSETS}, {o: 0.9 for o in OFFSETS})] * 3)
    rows = adapter.compare_rows({"sim_telarray": noisy_sim, "obdeect": obd}, OFFSETS)
    gate = rows[0]["rho_baffle"]
    # The limit is not widened by the seed spread; a spread whose 3 sigma exceeds the limit makes
    # the comparison underpowered (it fails, more seeds are needed).
    assert gate["sigma_seed"] > 0
    assert gate["threshold"] == pytest.approx(adapter.LIMIT)
    assert gate["exceeded"] is (abs(gate["delta"]) > adapter.LIMIT)
    assert gate["underpowered"] is (3 * gate["sigma_seed"] > adapter.LIMIT)


def test_config_rejects_missing_sections_and_backend_in_options(tmp_path: Path):
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"simtools_config": {}, "analytic": {}}))
    with pytest.raises(adapter.StudyError, match="exactly"):
        adapter._config(bad)

    base = {
        "simtools_config": {
            "application": "simtools-validate-optics",
            "configuration": {"ray_tracing_backend": "obdeect"},
        },
        "analytic": {"model_version": "7.0.0", "photons": 10, "seed": 1},
        "compare": {"seeds": [1, 2, 3], "offsets_deg": [[0, 0]]},
    }
    path = tmp_path / "cfg.json"
    path.write_text(json.dumps(base))
    with pytest.raises(adapter.StudyError, match="omit ray_tracing_backend"):
        adapter._config(path)


def test_config_needs_three_seeds(tmp_path: Path):
    value = {
        "simtools_config": {
            "application": "simtools-validate-optics",
            "configuration": {"site": "South"},
        },
        "analytic": {"model_version": "7.0.0", "photons": 10, "seed": 1},
        "compare": {"seeds": [1, 2], "offsets_deg": [[0, 0]]},
    }
    path = tmp_path / "cfg.json"
    path.write_text(json.dumps(value))
    with pytest.raises(adapter.StudyError, match="at least 3 seeds"):
        adapter._config(path)


def test_towards_is_unit_vector_on_axis_and_tilted():
    x, y, z = adapter._towards((0.0, 0.0))
    assert (x, y, z) == pytest.approx((0.0, 0.0, 1.0))
    x, y, z = adapter._towards((1.0, 0.0))
    assert x * x + y * y + z * z == pytest.approx(1.0)
    assert x > 0
