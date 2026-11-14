"""Pure helpers of the calibration-source adapter (source_camera)."""

import numpy as np
import pytest

from obdeect_validation.adapters import source_camera as sc


def _stream(n, seed, shift=0.0, tshift=0.0):
    rng = np.random.default_rng(seed)
    return {
        "cx": rng.normal(shift, 0.01, n),
        "cy": rng.normal(0.0, 0.01, n),
        "cz": -np.ones(n),
        "t": rng.normal(53.4 + tshift, 0.04, n),
        "z": np.full(n, 16.0),
        "w": np.ones(n),
        "wavelength": np.full(n, 400.0),
        "r": np.full(n, 10.0),
    }


def test_ks_snaps_round_off_zeros():
    d, p = sc._ks(np.zeros(100), np.full(100, 6.1e-17))
    assert d == 0.0
    assert p == pytest.approx(1.0)


def test_ks_detects_shape_difference():
    rng = np.random.default_rng(1)
    d, p = sc._ks(rng.normal(0, 1, 20000), rng.normal(0.2, 1, 20000))
    assert d > 0.05 and p < 1e-3


def test_emission_identical_streams_have_no_flags():
    a, b = _stream(5000, 1), _stream(5000, 2)
    m = sc.emission_metrics([2500.0, 2500.0], a, [2500.0, 2500.0], b, "none")
    assert m["photon_ratio_limit_exceeded"] == 0
    assert m["t_mean_limit_exceeded"] == 0
    assert m["event_count_delta"] == 0


def test_emission_flags_normalisation_and_timing():
    a, b = _stream(5000, 1), _stream(5000, 2, tshift=0.5)
    m = sc.emission_metrics([2600.0, 2600.0], a, [2500.0, 2500.0], b, "none")
    assert m["photon_ratio_limit_exceeded"] == 1
    assert m["t_mean_limit_exceeded"] == 1


def test_emission_single_event_uses_poisson_error():
    a, b = _stream(5000, 1), _stream(5000, 2)
    m = sc.emission_metrics([400000.0], a, [399959.0], b, "none")
    assert np.isfinite(m["photon_ratio_se"])
    assert m["photon_ratio_limit_exceeded"] == 0


def test_emission_excludes_reference_bunches_beyond_fiducial_radius():
    a, b = _stream(5000, 1), _stream(5000, 2)
    b["t"] = b["t"].copy()
    b["t"][:50] = -2.2e8
    b["r"] = b["r"].copy()
    b["r"][:50] = 2.0e4
    m = sc.emission_metrics([2500.0, 2500.0], a, [2500.0, 2500.0], b, "none")
    assert m["t_reference_outside_radius_fraction"] == pytest.approx(50 / 5000)
    assert abs(m["t_mean_reference"]) < 100
    assert m["t_mean_limit_exceeded"] == 0


def test_camera_identical_samples_pass():
    rng = np.random.default_rng(3)
    pe = {1: (rng.poisson(100, 500).astype(float), rng.poisson(5, 1855).astype(float))}
    m = sc.camera_metrics(pe, pe)
    assert m["pe_ratio"] == pytest.approx(1.0)
    assert m["pe_ratio_limit_exceeded"] == 0
    assert m["image_reduced_chi2"] == pytest.approx(0.0)
    assert m["image_chi2_limit_exceeded"] == 0


def test_camera_flags_gain_difference():
    rng = np.random.default_rng(4)
    a = {1: (rng.poisson(100, 800).astype(float), rng.poisson(5, 1855).astype(float))}
    b = {1: (rng.poisson(120, 800).astype(float), rng.poisson(5, 1855).astype(float))}
    m = sc.camera_metrics(a, b)
    assert m["pe_ratio"] == pytest.approx(1.2, rel=0.03)
    assert m["pe_ratio_limit_exceeded"] == 1


def test_camera_requires_common_telescope():
    with pytest.raises(ValueError):
        sc.camera_metrics({1: (np.ones(3), np.ones(4))}, {2: (np.ones(3), np.ones(4))})


def test_accepted_ratio_ignores_photons_outside_the_radius():
    # Native photons from z = 1000 m on axis land at r = 1000 m * |cx|; half land beyond 14 m.
    cx = np.r_[np.full(500, 0.010), np.full(500, 0.0150)]
    native = {"x": np.zeros(1000), "y": np.zeros(1000), "z": np.full(1000, 1000.0), "cx": cx,
              "cy": np.zeros(1000), "cz": -np.sqrt(1 - cx**2), "w": np.ones(1000)}
    reference = {"r": np.full(500, 1000.0), "w": np.ones(500)}
    ratio, se = sc.accepted_ratio(native, 1, reference, 1, 1400.0)
    assert ratio == pytest.approx(1.0)
    assert se > 0
