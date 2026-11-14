"""Time model of the flasher-secondary comparison with ff-gct (source_camera.ff_gct_time_metrics)."""

import numpy as np
import pytest

from obdeect_validation.adapters import source_camera as sc


def _samples(cz_values, extra_ref_cz=()):
    """Native photons at z = 0.5 m, zero emission time; reference bunches for the same photons."""
    cz = np.asarray(cz_values, dtype=np.float64)
    cx = np.sqrt(1.0 - cz**2)
    z = np.full(cz.size, 0.5)
    path = -z / cz  # m, signed
    native = {
        "x": np.zeros(cz.size),
        "y": np.zeros(cz.size),
        "z": z,
        "cx": cx,
        "cy": np.zeros(cz.size),
        "cz": cz,
        "t": np.zeros(cz.size),
        "w": np.ones(cz.size),
        "wavelength": np.full(cz.size, 400.0),
    }
    reference = {
        "x": path * cx * 100.0,
        "y": np.zeros(cz.size),
        "r": np.abs(path * cx) * 100.0,
        "cx": cx,
        "cy": np.zeros(cz.size),
        "cz": -np.sqrt(1.0 - cx**2),
        "t": np.abs(path) / sc.C_M_NS,
        "w": np.ones(cz.size),
        "wavelength": np.full(cz.size, 400.0),
    }
    return native, reference


def test_time_block_identical_samples_have_no_flags():
    x = np.linspace(1.0, 5.0, 200)
    block = sc._time_block(x, np.ones(200), x, np.ones(200))
    assert block["t_mean_delta"] == pytest.approx(0.0)
    assert block["t_rms_delta"] == pytest.approx(0.0)
    assert block["t_mean_limit_exceeded"] == 0
    assert block["t_rms_limit_exceeded"] == 0


def test_ff_gct_time_keeps_photons_heading_away_and_applies_same_selection():
    # Upward photon (cz > 0) is kept with its backward-extended impact and time |path| / c.
    # The grazing photon (|cz| = 0.1 < 0.3) is excluded on the native side and on the reference side.
    native, reference = _samples([-0.5, 0.5, -0.1, -0.9])
    reference["t"] = np.append(reference["t"], 100.0)
    reference["r"] = np.append(reference["r"], 10.0)
    reference["cz"] = np.append(reference["cz"], -0.1)
    reference["cx"] = np.append(reference["cx"], 0.0)
    reference["cy"] = np.append(reference["cy"], 0.0)
    reference["x"] = np.append(reference["x"], 0.0)
    reference["y"] = np.append(reference["y"], 0.0)
    reference["w"] = np.append(reference["w"], 1.0)
    reference["wavelength"] = np.append(reference["wavelength"], 400.0)
    m = sc.ff_gct_time_metrics(native, reference)
    assert m["native_photons_moving_away_fraction"] == pytest.approx(0.25)
    assert m["t_native_excluded_fraction"] == pytest.approx(0.25)
    assert m["t_reference_excluded_fraction"] == pytest.approx(0.4)
    assert m["time_correction_ns_mean"] == pytest.approx(
        np.mean([1.0, 1.0, 0.5 / 0.9]) / sc.C_M_NS  # selected photons only
    )
    # Identical arrival times on both sides (the excluded reference bunch at 100 ns is not used).
    assert m["t_mean_delta"] == pytest.approx(0.0, abs=1e-12)
    assert m["t_rms_delta"] == pytest.approx(0.0, abs=1e-12)
    assert m["t_mean_limit_exceeded"] == 0
    assert m["t_rms_limit_exceeded"] == 0


def test_ff_gct_time_flags_a_shifted_reference():
    rng = np.random.default_rng(3)
    cz = rng.uniform(0.3, 1.0, 4000) * rng.choice([-1.0, 1.0], 4000)
    native, reference = _samples(cz)
    reference["t"] = reference["t"] + 0.5
    m = sc.ff_gct_time_metrics(native, reference)
    assert m["t_mean_delta"] == pytest.approx(-0.5)
    assert m["t_mean_limit_exceeded"] == 1


def test_ff_gct_time_rejects_photons_parallel_to_plane():
    native, reference = _samples([-0.5, 0.5])
    native["cz"] = np.array([-0.5, 0.0])
    with pytest.raises(ValueError):
        sc.ff_gct_time_metrics(native, reference)


def test_camera_single_event_uses_poisson_error():
    reference = {0: (np.array([10000.0]), np.full(100, 100.0))}
    candidate = {0: (np.array([10100.0]), np.full(100, 101.0))}
    m = sc.camera_metrics(reference, candidate)
    expected = (10100.0 / 10000.0) * np.hypot(
        1 / np.sqrt(10000.0), 1 / np.sqrt(10100.0)
    )
    assert m["pe_ratio"] == pytest.approx(1.01)
    assert m["pe_ratio_se"] == pytest.approx(expected)
    assert np.isfinite(m["image_reduced_chi2"])
