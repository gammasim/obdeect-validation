"""Unit tests of the independent curved-ray tracer (numpy only, no obdeect).

The tracer is the reference for propagation.refraction-curved-ray. These tests check it against
closed forms and against independently written quadratures, so that a failing comparison with
obdeect can be attributed to one side or the other.
"""

import ast
import math
from pathlib import Path

import numpy as np
import pytest

from obdeect_validation.adapters import refraction_curved_reference as ref

SCALE_HEIGHT_M = 8000.0
NU0 = 3.0e-4
EARTH = ref.EARTH_RADIUS_M
FLAT = 1.0e12


def exponential(h):
    return NU0 * np.exp(-np.asarray(h, dtype=np.float64) / SCALE_HEIGHT_M)


def impact(radius, h_t, zenith_deg):
    return (radius + h_t) * math.sin(math.radians(zenith_deg))


def plane_parallel_shift_trapezoid(theta_deg, h_t, h_e, nu, points=2_000_001):
    """Plane-parallel shift by trapezoid in h: integral of tan(theta(h)) - tan(theta_s).

    Independent of the spherical-shell code: a constant vacuum angle, Snell's law in layers
    sin(theta(h)) = sin(theta_s) n_e / n(h), and a plain trapezoid rule.
    """
    sin_s = math.sin(math.radians(theta_deg))
    n_e = 1.0 + nu(h_e)
    heights = np.linspace(h_t, h_e, points)
    sin_h = sin_s * n_e / (1.0 + nu(heights))
    tan_h = sin_h / np.sqrt(1.0 - sin_h**2)
    tan_s = sin_s / math.sqrt(1.0 - sin_s**2)
    values = tan_h - tan_s
    return float(np.trapezoid(values, heights))


def test_flat_limit_shift_matches_first_order_plane_parallel():
    # First order in nu: shift = sin(theta)/cos^3(theta) * int (nu_e - nu) dh with the
    # exponential integral int nu dh = H (nu_t - nu_e). The next order is ~ nu / cos^2(theta).
    h_t, h_e = 2000.0, 20000.0
    for zenith in (30.0, 60.0, 80.0):
        theta = math.radians(zenith)
        b = impact(FLAT, h_t, zenith)
        got = ref.ray_integrals(exponential, FLAT, h_t, h_e, b)["shift_m"]
        nu_e, nu_t = (
            NU0 * math.exp(-h_e / SCALE_HEIGHT_M),
            NU0 * math.exp(-h_t / SCALE_HEIGHT_M),
        )
        first = (
            math.sin(theta)
            / math.cos(theta) ** 3
            * ((nu_e * (h_e - h_t)) - SCALE_HEIGHT_M * (nu_t - nu_e))
        )
        assert got == pytest.approx(first, rel=2.0 * NU0 / math.cos(theta) ** 2)


def test_flat_limit_shift_matches_independent_plane_parallel_integral():
    # The spherical shift differs from the plane-parallel one by a term of order dh/R: measured
    # relative differences are 1.2e-7 (70 deg) at R = 1e12, 1.2e-8 at 1e13 and 1.2e-9 at 1e14.
    h_t, h_e = 2000.0, 20000.0
    for zenith in (20.0, 45.0, 70.0, 80.0):
        b = impact(FLAT, h_t, zenith)
        got = ref.ray_integrals(exponential, FLAT, h_t, h_e, b)["shift_m"]
        independent = plane_parallel_shift_trapezoid(zenith, h_t, h_e, exponential)
        assert got == pytest.approx(independent, rel=1e-6)


def test_flat_limit_path_is_the_straight_distance():
    # The O(dh/R) correction of the path is 2.7e-8 relative at R = 1e12 and 2.7e-10 at 1e14.
    h_t, h_e, zenith = 2000.0, 20000.0, 60.0
    b = impact(1.0e14, h_t, zenith)
    got = ref.ray_integrals(exponential, 1.0e14, h_t, h_e, b)["path_m"]
    expected = (h_e - h_t) / math.cos(math.radians(zenith))
    assert got == pytest.approx(expected, rel=1e-7)


def test_spherical_path_is_closed_form():
    # int dr / cos(theta_vac) = int r dr / sqrt(r^2 - b^2) = sqrt(r_e^2 - b^2) - sqrt(r_t^2 - b^2)
    h_t, h_e, zenith = 2000.0, 60000.0, 85.0
    b = impact(EARTH, h_t, zenith)
    r_t, r_e = EARTH + h_t, EARTH + h_e
    expected = math.sqrt(r_e * r_e - b * b) - math.sqrt(r_t * r_t - b * b)
    got = ref.ray_integrals(exponential, EARTH, h_t, h_e, b)["path_m"]
    assert got == pytest.approx(expected, rel=1e-9)


def test_vertical_ray_has_no_shift_and_excess_is_the_air_delay():
    h_t, h_e = 0.0, 20000.0
    result = ref.ray_integrals(exponential, EARTH, h_t, h_e, 0.0)
    assert result["shift_m"] == 0.0
    # theta = 0: n / cos(theta_ref) - 1 / cos(theta_vac) = nu, so excess equals the air delay.
    assert result["excess_ns"] == pytest.approx(result["straight_ns"], rel=1e-12)


def test_straight_delay_vertical_matches_exponential_closed_form():
    # (1/c) int_0^dh nu0 exp(-h/H) dh = nu0 H (1 - exp(-dh/H)) / c
    dh = 20000.0
    result = ref.ray_integrals(exponential, EARTH, 0.0, dh, 0.0)
    closed = (
        NU0 * SCALE_HEIGHT_M * (1.0 - math.exp(-dh / SCALE_HEIGHT_M)) / ref.C_M_PER_NS
    )
    assert result["straight_ns"] == pytest.approx(closed, rel=1e-10)


def test_straight_delay_slanted_flat_limit_matches_closed_form():
    # (1/c) int nu / cos(theta) dh = nu0 H (1 - exp(-dh/H)) / (c cos theta) in the flat limit
    h_t, h_e, zenith = 0.0, 20000.0, 50.0
    b = impact(FLAT, h_t, zenith)
    result = ref.ray_integrals(exponential, FLAT, h_t, h_e, b)
    closed = (
        NU0
        * SCALE_HEIGHT_M
        * (1.0 - math.exp(-(h_e - h_t) / SCALE_HEIGHT_M))
        / (ref.C_M_PER_NS * math.cos(math.radians(zenith)))
    )
    assert result["straight_ns"] == pytest.approx(closed, rel=1e-8)


def test_shift_is_negative_for_positive_nu_and_bounded_by_plane_parallel():
    h_t, h_e = 2000.0, 20000.0
    for zenith in (20.0, 60.0, 80.0):
        b = impact(EARTH, h_t, zenith)
        shift = ref.ray_integrals(exponential, EARTH, h_t, h_e, b)["shift_m"]
        assert shift < 0.0
        flat = ref.ray_integrals(exponential, FLAT, h_t, h_e, impact(FLAT, h_t, zenith))
        # The curved shift is smaller in magnitude than the plane-parallel one (docs: 7 % at 80 deg).
        assert abs(shift) < abs(flat["shift_m"])


def test_quadrature_converges_under_refinement():
    h_t, h_e, zenith = 2000.0, 60000.0, 85.0
    b = impact(EARTH, h_t, zenith)
    values = {
        n: ref.ray_integrals(exponential, EARTH, h_t, h_e, b, intervals=n)
        for n in (2**10, 2**12, 2**14)
    }
    for key in ("shift_m", "excess_ns", "straight_ns"):
        coarse = abs(values[2**10][key] - values[2**14][key])
        middle = abs(values[2**12][key] - values[2**14][key])
        assert middle < coarse  # the error decreases as the step is refined
        finer = ref.ray_integrals(exponential, EARTH, h_t, h_e, b, intervals=2**15)[key]
        change = abs(finer - values[2**14][key]) / max(abs(values[2**14][key]), 1e-9)
        assert change < 1e-8  # halving the step at the production resolution


def test_simpson_weights_integrate_polynomials_exactly():
    # Weights scaled by the interval length 2 integrate u^3 over [0, 2] exactly: int = 4.
    intervals = 64
    u = np.linspace(0.0, 2.0, intervals + 1)
    integral = float(np.sum(ref.simpson_weights(intervals) * 2.0 * u**3))
    assert integral == pytest.approx(4.0, rel=1e-12)


def test_simpson_weights_reject_odd_intervals():
    with pytest.raises(ValueError):
        ref.simpson_weights(3)


def test_log_spline_is_exact_at_knots_and_for_log_linear_profiles():
    x = np.linspace(0.0, 10000.0, 11)
    nu = np.exp(-x / 3000.0) * 1e-4
    spline = ref.LogSpline(x, nu)
    np.testing.assert_allclose(np.exp(spline(x)), nu, rtol=1e-12)
    between = np.linspace(0.0, 10000.0, 1001)
    np.testing.assert_allclose(
        np.exp(spline(between)), np.exp(-between / 3000.0) * 1e-4, rtol=1e-12
    )
    assert spline.m[0] == pytest.approx(0.0, abs=1e-12)  # natural boundary conditions
    assert spline.m[-1] == pytest.approx(0.0, abs=1e-12)


def test_log_spline_is_natural_for_nonlinear_log():
    x = np.array([0.0, 1000.0, 2000.0, 3000.0, 4000.0])
    y = np.array([0.0, 1.0, 0.0, 1.0, 0.0])
    spline = ref.LogSpline(x, np.exp(y))
    np.testing.assert_allclose(spline(x), y, atol=1e-12)
    assert spline.m[0] == pytest.approx(0.0, abs=1e-12)
    assert spline.m[-1] == pytest.approx(0.0, abs=1e-12)


def test_emission_point_has_requested_altitude_and_zenith():
    for radius in (EARTH, FLAT):
        for height in (3000.0, 60000.0):
            for zenith in (0.0, 40.0, 85.0):
                point, direction = ref.emission_point(height, zenith, radius)
                centre = np.array([0.0, 0.0, -(radius + ref.SITE_ALTITUDE_M)])
                altitude = float(np.linalg.norm(point - centre)) - radius
                # float64 resolution at |P - C| ~ 1e12 m is about 1e-4 m.
                tolerance = 1e-5 if radius == EARTH else 1e-3
                assert altitude == pytest.approx(height, abs=tolerance)
                assert np.linalg.norm(direction) == pytest.approx(1.0, abs=1e-14)
                b = float(np.linalg.norm(np.cross(point - centre, direction)))
                sin_s = b / (radius + ref.SITE_ALTITUDE_M)
                assert math.degrees(math.asin(sin_s)) == pytest.approx(zenith, abs=1e-6)


def test_limits_document_the_rules():
    assert set(ref.LIMITS) == {
        "position",
        "time",
        "direction",
        "straight_delay",
        "plane_limit",
        "convergence",
    }
    assert ref.position_limit(80.0, -1000.0) == pytest.approx(1.0)
    assert ref.position_limit(85.0, -1000.0) == 0.05
    assert ref.position_limit(80.0, -1.0) == 0.005


def test_tracer_has_no_module_level_obdeect_import():
    # obdeect is imported inside the candidate-run functions only, never by the tracer.
    tree = ast.parse(Path(ref.__file__).read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Import):
            assert all(alias.name.split(".")[0] != "obdeect" for alias in node.names)
        if isinstance(node, ast.ImportFrom):
            assert node.module is None or node.module.split(".")[0] != "obdeect"
