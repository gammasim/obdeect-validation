"""Pure helpers of the performance regression case (no obdeect run)."""

import math

import pytest

from obdeect_validation.adapters import performance_common as pcm
from obdeect_validation.adapters import performance_regression as pr


def test_host_key_is_platform_machine_and_brand_sanitised() -> None:
    assert (
        pcm.host_key("Darwin", "arm64", "Apple M3 Max") == "darwin-arm64-apple-m3-max"
    )
    assert (
        pcm.host_key("Linux", "x86_64", "Intel(R) Xeon(R) CPU @ 2.40GHz")
        == "linux-x86-64-intel-r-xeon-r-cpu-2-40ghz"
    )


def test_host_key_rejects_an_empty_identity() -> None:
    with pytest.raises(ValueError):
        pcm.host_key("", "", "")


def test_summarise_takes_the_best_repetition_and_the_relative_spread() -> None:
    record = pr.summarise([0.5, 0.6, 0.55], photons=1_000_000)
    assert record["ns_per_photon"] == pytest.approx([500.0, 600.0, 550.0])
    assert record["best_ns_per_photon"] == pytest.approx(500.0)
    assert record["spread"] == pytest.approx(0.2)


def test_summarise_rejects_a_non_positive_photon_count() -> None:
    with pytest.raises(ValueError):
        pr.summarise([1.0], photons=0)


def _baseline(best: float, spread: float) -> dict:
    return {"best_ns_per_photon": best, "spread": spread}


def test_classify_allows_the_floor_of_five_percent_without_spread() -> None:
    result = pr.classify(_baseline(104.0, 0.01), _baseline(100.0, 0.01))
    assert result["threshold"] == pytest.approx(1.05)
    assert result["regressed"] == 0


def test_classify_flags_a_slowdown_above_the_floor() -> None:
    result = pr.classify(_baseline(106.0, 0.01), _baseline(100.0, 0.01))
    assert result["ratio"] == pytest.approx(1.06)
    assert result["regressed"] == 1


def test_classify_widens_the_tolerance_to_three_spreads_of_the_noisier_side() -> None:
    # spread 0.05 on the baseline side: tolerance 1 + 3 x 0.05 = 1.15
    result = pr.classify(_baseline(110.0, 0.01), _baseline(100.0, 0.05))
    assert result["threshold"] == pytest.approx(1.15)
    assert result["spread_used"] == pytest.approx(0.05)
    assert result["regressed"] == 0


def test_speedups_are_not_regressions() -> None:
    result = pr.classify(_baseline(50.0, 0.01), _baseline(100.0, 0.01))
    assert result["ratio"] == pytest.approx(0.5)
    assert result["regressed"] == 0


def test_compare_counts_regressions_and_missing_baselines() -> None:
    measured = {
        "a": pr.summarise([1.0, 1.0], photons=1_000_000),
        "b": pr.summarise([2.0, 2.0], photons=1_000_000),
        "c": pr.summarise([1.0, 1.0], photons=1_000_000),
    }
    baseline = {
        "a": {"best_ns_per_photon": 1000.0, "spread": 0.0},
        "b": {"best_ns_per_photon": 1000.0, "spread": 0.0},
    }
    metrics, details = pr.compare(measured, baseline)
    assert metrics["regressions"] == 1  # b: 2000 ns against 1000 ns
    assert metrics["baseline_missing"] == 1  # c
    assert details["c"]["baseline"] is None
    assert metrics["a_ratio"] == pytest.approx(1.0)


def test_dish_radius_of_facets_adds_the_circumradius_to_the_centre() -> None:
    document = {
        "optics": {
            "primary": {
                "facets": [
                    {"centre_m": [3.0, 4.0, 0.0], "diameter_m": 2.0, "shape": "round"},
                    {
                        "centre_m": [0.0, 0.0, 0.0],
                        "diameter_m": 2.0 * math.sqrt(3.0),
                        "shape": "hexagon",
                    },
                ]
            }
        }
    }
    assert pr.dish_radius_m(document) == pytest.approx(6.0)


def test_dish_radius_of_an_asphere_is_its_outer_radius() -> None:
    document = {"optics": {"primary": {"asphere": {"outer_radius_m": 2.12}}}}
    assert pr.dish_radius_m(document) == pytest.approx(2.12)


def test_dish_radius_rejects_a_primary_without_geometry() -> None:
    with pytest.raises(ValueError):
        pr.dish_radius_m({"optics": {"primary": {}}})
