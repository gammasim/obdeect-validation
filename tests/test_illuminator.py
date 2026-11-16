"""Geometry and per-pixel map helpers of the illuminator adapter (illuminator.py)."""

import math

import numpy as np
import pytest

from obdeect_validation.adapters import illuminator as ill


def test_pointing_is_unit_vector_from_illuminator_to_telescope():
    unit, theta, phi = ill.pointing([100.0, 0.0, 50.0], [0.0, 0.0, 0.0])
    assert np.allclose(unit, [-0.894427, 0.0, -0.447214], atol=1e-6)
    # Telescope zenith: 180 - arccos(dz/|d|); here dz/|d| = -0.447 -> 180 - 116.57 = 63.43 deg
    assert theta == pytest.approx(
        180.0 - math.degrees(math.acos(-50 / math.sqrt(12500)))
    )
    # atan2(0, -100) = 180 deg, and simtools uses 180 - atan2 -> 0 deg.
    assert phi == pytest.approx(0.0)


def test_scaled_illuminator_stays_on_the_line_of_sight():
    telescope = np.array([1.0, 2.0, 3.0])
    source = np.array([101.0, 2.0, 103.0])
    for scale in (0.5, 1.0, 2.0):
        moved = ill.scaled_illuminator(source, telescope, scale)
        assert np.allclose(moved - telescope, scale * (source - telescope))
        unit_before, _, _ = ill.pointing(source, telescope)
        unit_after, _, _ = ill.pointing(moved, telescope)
        assert np.allclose(unit_before, unit_after)


def test_pixel_map_consistent_with_identical_samples():
    rng = np.random.default_rng(1)
    mean = np.full(50, 30.0)
    ref = rng.poisson(mean, size=(40, 50)).astype(float)
    cand = rng.poisson(mean, size=(40, 50)).astype(float)
    result = ill.pixel_map_metrics(ref, cand)
    assert result["pixel_map_pixels_used"] == 50
    assert 0.3 < result["pixel_map_reduced_chi2"] < 2.0
    assert result["pixel_map_reduced_chi2_limit_exceeded"] == 0


def test_pixel_map_flags_a_5_percent_scale_difference_with_enough_events():
    rng = np.random.default_rng(2)
    mean = np.full(200, 200.0)
    ref = rng.poisson(mean, size=(400, 200)).astype(float)
    cand = rng.poisson(1.05 * mean, size=(400, 200)).astype(float)
    result = ill.pixel_map_metrics(ref, cand)
    assert result["pixel_map_reduced_chi2_limit_exceeded"] == 1


def test_pixel_map_needs_two_events():
    with pytest.raises(ValueError):
        ill.pixel_map_metrics(np.ones((1, 5)) * 10, np.ones((3, 5)) * 10)


def test_bent_time_mode_uses_arrival_times_and_landing_radius():
    from obdeect_validation.adapters import source_camera as sc

    native = {
        "t": np.array([10.0, 20.0, 30.0]),
        "x": np.array([0.0, 0.1, 0.0]),  # landing [m]: inside 15 m, inside, outside
        "y": np.array([0.0, 0.0, 16.0]),
    }
    ref = {"t": np.array([1.0, 2.0]), "r": np.array([0.0, 2000.0])}
    t_native, native_ok, ref_ok, toward = sc.emission_times(native, ref, "bent")
    assert np.array_equal(t_native, native["t"])
    assert native_ok.tolist() == [True, True, False]
    assert ref_ok.tolist() == [True, False]
    assert toward is None


def _exponential_profile():
    heights = np.linspace(0.0, 30000.0, 301)
    index_minus_one = 2.9e-4 * np.exp(-heights / 8000.0)
    return heights, index_minus_one


def test_bent_excess_vertical_ray_is_straight_air_delay():
    from obdeect_validation.adapters import source_camera as sc

    heights, nm1 = _exponential_profile()
    z = 63.0
    exact = ill.bent_excess_ns(z, 1.0, heights, nm1, obs_m=2150.0)
    h = np.linspace(2150.0, 2150.0 + z, 20001)
    air = np.trapezoid(np.exp(np.interp(h, heights, np.log(nm1))), h) / sc.C_M_NS
    assert exact == pytest.approx(air, rel=1e-6, abs=1e-9)


def test_bent_excess_grazing_differs_from_straight_air_delay():
    from obdeect_validation.adapters import source_camera as sc

    heights, nm1 = _exponential_profile()
    z, cz = 63.0, 0.0562
    exact = ill.bent_excess_ns(z, cz, heights, nm1, obs_m=2150.0)
    h = np.linspace(2150.0, 2150.0 + z, 20001)
    straight = np.trapezoid(np.exp(np.interp(h, heights, np.log(nm1))), h) / cz
    straight /= sc.C_M_NS
    assert abs(exact - straight) > 0.1  # ray bending matters at grazing emission


def test_cone_ratio_expected_is_close_to_margin_squared_at_large_distance():
    assert ill.cone_ratio_expected(1120.0) == pytest.approx(1.02**2, rel=2e-4)
    assert ill.cone_ratio_expected(200.0) > 1.0404
