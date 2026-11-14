import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from obdeect_validation.adapters import source_comparison as sc  # noqa: E402


def test_ks_statistic_identical_and_disjoint_samples():
    a = np.linspace(0.0, 1.0, 101)
    assert sc.ks_statistic(a, a) == 0.0
    assert sc.ks_statistic(a, a + 2.0) == pytest.approx(1.0)


def test_ks_one_sample_uniform_grid():
    u = (np.arange(1000) + 0.5) / 1000.0
    assert sc.ks_one_sample(u, lambda x: x) == pytest.approx(5e-4, abs=1e-12)


def test_ks_critical_alpha_1e3():
    assert sc.ks_critical(10000) == pytest.approx(1.9486 / 100.0, rel=1e-3)
    assert sc.ks_critical(100, 100) == pytest.approx(1.9486 * math.sqrt(2.0 / 100.0), rel=1e-3)


def test_ratio_limit_floor_and_statistical_part():
    assert sc.ratio_limit(0.0) == 0.01
    assert sc.ratio_limit(0.005) == pytest.approx(0.025)


def test_frank_tamm_matches_closed_form_and_magnitude():
    n_minus_one = 2.9e-4
    value = sc.frank_tamm_per_metre(n_minus_one, 300.0, 600.0)
    cos_c = 1.0 / (1.0 + n_minus_one)
    expected = 2 * math.pi * sc.ALPHA * (1 - cos_c**2) * (1 / 300 - 1 / 600) * 1e9
    assert value == pytest.approx(expected, rel=1e-12)
    assert 40.0 < value < 50.0  # about 44.6 per metre for a sea-level-like index


def test_cherenkov_cos_threshold_and_beta():
    assert sc.cherenkov_cos(0.0) == 1.0
    assert sc.cherenkov_cos(2.9e-4, beta=0.5) == pytest.approx(2.0 / (1.00029), rel=1e-12)


def test_night_sky_photons_reproduce_configured_rate():
    photons = sc.night_sky_photons(0.24, 10.0, 15.0, 1.0)
    omega = 2 * math.pi * (1 - math.cos(math.radians(10.0)))
    area_cm2 = math.pi * 1500.0**2
    assert photons == pytest.approx(0.24 * omega * area_cm2, rel=1e-12)
    assert abs(photons - 161937) < 1.0  # configured photons per event used in the case


def test_observation_crossing_and_impact_parameter():
    p = np.array([[0.0, 0.0, 10.0], [20.0, 0.0, 5.0]])
    d = np.array([[0.0, 0.0, -1.0], [0.0, 0.0, -1.0]])
    crossing = sc.observation_crossing(p, d)
    assert np.allclose(crossing[:, 2], 0.0)
    assert np.allclose(crossing[:, :2], [[0.0, 0.0], [20.0, 0.0]])
    assert np.allclose(sc.impact_parameter(p, d), [0.0, 20.0])


def test_camera_comparison_flags_bias():
    rng = np.random.default_rng(1)
    a = rng.normal(2000.0, 50.0, 200)
    same = sc.camera_comparison(a, a.copy())
    assert same["camera_pe_ratio_limit_exceeded"] == 0
    assert same["camera_pe_ratio"] == pytest.approx(1.0)
    biased = sc.camera_comparison(a * 1.1, a)
    assert biased["camera_pe_ratio_limit_exceeded"] == 1


def test_camera_comparison_requires_two_events():
    with pytest.raises(ValueError):
        sc.camera_comparison(np.array([1.0]), np.array([1.0, 2.0]))


def test_muon_expected_hits_without_sphere_limit_equals_track_yield(monkeypatch):
    # With a huge sphere every Cherenkov photon hits it: expected = yield per metre x 3000 m.
    monkeypatch.setitem(sc.MUON, "radius_m", 1.0e5)
    heights = np.arange(0.0, 3001.0, 100.0)
    table = np.column_stack((heights, np.full(heights.size, 2.9e-4)))
    expected = sc.muon_expected_hits(table, step_m=5.0, azimuths=720)
    per_metre = sc.frank_tamm_per_metre(2.9e-4, 300.0, 600.0)
    assert expected == pytest.approx(per_metre * 3000.0, rel=2e-3)
