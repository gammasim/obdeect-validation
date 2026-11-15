"""Reference functions of the curvature cases (no obdeect needed).

The Chapman integral against its closed forms (zenith, horizon, flat limit), the exponential
atmosphere against the integral of its own density, the table writer and reader, and the
geometry helpers of the sky-frame case.
"""

import math

import numpy as np
import pytest
from obdeect_validation.adapters import curvature_chapman as chapman_case
from obdeect_validation.adapters import curvature_geometry as geometry_case


def _bessel_k1(x):
    """K_1(x) = integral_0^inf exp(-x cosh t) cosh t dt (integral representation)."""
    t = np.linspace(0.0, 12.0, 240001)
    f = np.exp(-x * np.cosh(t)) * np.cosh(t)
    return _simpson_grid(f, t)


def _simpson_grid(values, grid):
    intervals = len(grid) - 1
    weights = np.ones(len(grid))
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    return (grid[-1] - grid[0]) / (3.0 * intervals) * float(np.sum(weights * values))


def test_chapman_at_zenith_is_one():
    # Ch(x, 1) = 1: the column straight up from r is H exactly.
    assert chapman_case.chapman(800.0, 1.0, intervals=200000) == pytest.approx(1.0, rel=1e-7)


@pytest.mark.parametrize("x", [5.0, 20.0])
def test_chapman_at_horizon_matches_bessel_closed_form(x):
    # Ch(x, 0) = x e^x K_1(x).
    expected = x * math.exp(x) * _bessel_k1(x)
    actual = chapman_case.chapman(x, 0.0, intervals=400000, upper_y=40.0)
    assert actual == pytest.approx(expected, rel=1e-7)


def test_chapman_horizon_limit_sqrt_pi_x_over_2():
    # Asymptotically Ch(x, 0) = sqrt(pi x / 2) (1 + 3/(8x) + ...).
    x = 2000.0
    actual = chapman_case.chapman(x, 0.0, intervals=400000, upper_y=2.0)
    assert actual == pytest.approx(math.sqrt(math.pi * x / 2.0) * (1.0 + 3.0 / (8.0 * x)), rel=1e-4)


def test_chapman_flat_limit_is_secant_for_large_x():
    # Ch(x, c) -> 1/c for x -> infinity (plane-parallel airmass), correction of order 1/x.
    c = 0.5
    value = chapman_case.chapman(1e5, c, intervals=2000000, upper_y=2.0)
    assert value == pytest.approx(1.0 / c, rel=1e-4)


def test_slant_depth_chapman_matches_exact_exponential_integral():
    """The closed-form Chapman slant depth against a direct integral of the exponential."""
    h_t = chapman_case.SITE_M
    for zenith in (0.0, 45.0, 80.0):
        cos_t = math.cos(math.radians(zenith))
        z = 15000.0

        def density(h):
            return (chapman_case.VERTICAL_DEPTH / chapman_case.SCALE_HEIGHT_M) * np.exp(
                -(h - h_t) / chapman_case.SCALE_HEIGHT_M
            )

        def integrand(u):
            h = h_t + u * u
            return density(h) * 2.0 * u / np.sqrt(chapman_case.cos_theta_squared(h, h_t, cos_t))

        direct = chapman_case.simpson(integrand, 0.0, math.sqrt(z - h_t), 400000)
        closed = chapman_case.slant_depth_chapman(h_t, z, cos_t)
        assert closed == pytest.approx(direct, rel=1e-7)


def test_cos_theta_squared_stable_form_matches_naive_form():
    h_t = 2000.0
    for zenith in (30.0, 60.0):
        cos_t = math.cos(math.radians(zenith))
        for h in (2100.0, 5000.0, 40000.0):
            naive = (
                1.0
                - (1.0 - cos_t * cos_t)
                * ((chapman_case.EARTH_RADIUS_M + h_t) / (chapman_case.EARTH_RADIUS_M + h)) ** 2
            )
            stable = chapman_case.cos_theta_squared(h, h_t, cos_t)
            assert stable == pytest.approx(naive, rel=1e-12)


def test_cos_theta_squared_keeps_grazing_rays_small_corrections():
    # At 89.5 deg, 1 mm above the telescope, the curvature correction is sin^2 * 2 dh / R, about
    # 3e-10 on cos^2 = 7.6e-5. The stable form resolves it; the naive form in double would not.
    cos_t = math.cos(math.radians(89.5))
    sin_t = math.sin(math.radians(89.5))
    dh = 1e-3
    correction = chapman_case.cos_theta_squared(2000.0 + dh, 2000.0, cos_t) - cos_t * cos_t
    assert correction == pytest.approx(
        sin_t * sin_t * 2.0 * dh / chapman_case.EARTH_RADIUS_M, rel=1e-4
    )


def test_plane_parallel_limit_of_table_slant_depth(monkeypatch):
    """R -> infinity: the table form gives (tau(z) - tau(h_t)) / cos(theta_t)."""
    monkeypatch.setattr(chapman_case, "EARTH_RADIUS_M", 1e12)
    heights = np.array([2000.0, 2500.0, 4000.0, 8000.0, 16000.0, 40000.0])
    tau = chapman_case.exact_depth(heights)
    table = chapman_case.PiecewiseTable(heights, tau)
    for zenith in (0.0, 30.0, 60.0):
        cos_t = math.cos(math.radians(zenith))
        table_form = chapman_case.slant_depth_table(table, 2000.0, 40000.0, cos_t, intervals=200000)
        plane = chapman_case.plane_parallel_depth(table, 2000.0, 40000.0, cos_t)
        assert table_form == pytest.approx(plane, rel=1e-6)


def test_exact_depth_is_the_exponential_with_unit_vertical_depth():
    assert chapman_case.exact_depth(2000.0) == pytest.approx(0.0)
    assert chapman_case.exact_depth(1e9) == pytest.approx(1.0)


def test_table_file_round_trip_and_format(tmp_path):
    path = tmp_path / "table.dat"
    chapman_case.write_table(path)
    lines = path.read_text().splitlines()
    header = next(line for line in lines if line.startswith("# H2="))
    assert header.startswith("# H2= 2.000, H1=")
    heights, tau = chapman_case.read_table(path)
    assert heights.size == chapman_case.N_HEIGHTS + 1 == 100
    assert heights[0] == 2000.0 and heights[-1] == pytest.approx(100000.0)
    assert np.all(np.diff(heights) > 0.0)
    row = [line for line in lines if line.split()[0] == "400"][0].split()
    assert len(row) == chapman_case.N_HEIGHTS + 1
    assert np.allclose(tau, chapman_case.exact_depth(heights), atol=1e-10)


def test_table_interpolation_error_is_measured_and_small():
    """Log-linear interpolation of the exact depth: the table-only error in details.json."""
    heights = chapman_case.table_heights_km() * 1000.0
    table = chapman_case.PiecewiseTable(
        np.concatenate(([2000.0], heights)),
        chapman_case.exact_depth(np.concatenate(([2000.0], heights))),
    )
    probe = np.linspace(2000.0, 40000.0, 40001)
    error = np.abs(table.depth(probe) - chapman_case.exact_depth(probe))
    assert error.max() < 1e-4
    nodes = np.concatenate(([2000.0], heights))
    assert np.allclose(table.depth(nodes), chapman_case.exact_depth(nodes), atol=1e-12)


def test_table_derivative_is_zero_above_top_and_segmentwise():
    heights = np.array([2000.0, 4000.0, 8000.0])
    tau = np.array([0.0, 0.5, 0.75])
    table = chapman_case.PiecewiseTable(heights, tau)
    assert table.derivative(9000.0) == 0.0
    h = 3000.0
    expected = (0.5 / (math.log10(4000.0) - math.log10(2000.0))) / (h * math.log(10.0))
    assert table.derivative(h) == pytest.approx(expected)


def test_stored_direction_is_downward_and_stream_definition_is_consistent():
    cx, cy = chapman_case.stored_direction(60.0)
    cos_t = chapman_case.zenith_cosine(cx, cy)
    assert float(cos_t) == pytest.approx(0.5, rel=1e-6)
    cx0, cy0 = chapman_case.stored_direction(0.0)
    assert float(chapman_case.zenith_cosine(cx0, cy0)) == pytest.approx(1.0)


def test_sample_geometry_has_thirty_samples_in_order():
    samples = chapman_case.sample_geometry()
    assert len(samples) == 30
    assert samples[0] == (0, 5000.0) and samples[-1] == (89.5, 40000.0)


def test_telescope_table_and_photon_table_counts():
    telescopes = geometry_case.telescope_table()
    assert len(telescopes) == 12
    assert [t[0] for t in telescopes] == list(range(1, 13))
    photons = geometry_case.photon_table()
    assert len(photons) == 120
    assert sum(p["upward"] for p in photons) == 12
    assert all(float(p["cz"]) > 0 for p in photons if p["upward"])
    assert all(float(p["cz"]) < 0 for p in photons if not p["upward"])


def test_emission_altitude_on_axis_and_lateral_shift():
    # On the axis the sphere altitude is the flat one: zem = z + H.
    assert geometry_case.emission_altitude_sphere(0.0, 0.0, 4000.0) == pytest.approx(
        6000.0, abs=1e-6
    )
    # Lateral: the shift is d^2 / (2R) to (z + H) / R.
    radius = geometry_case.EARTH_RADIUS_M
    shift = float(geometry_case.emission_altitude_sphere(5000.0, 0.0, 3000.0)) - 5000.0
    expected = 5000.0**2 / (2.0 * radius)
    assert shift == pytest.approx(expected, rel=1e-2)


def test_lateral_shift_closed_form_against_d2_over_2r():
    radius = geometry_case.EARTH_RADIUS_M
    closed = geometry_case.lateral_shift_closed_form(5000.0, 1000.0)
    assert closed == pytest.approx(5000.0**2 / (2.0 * radius), rel=1e-2)
    # The closed form is the sphere result, exactly.
    direct = float(geometry_case.emission_altitude_sphere(5000.0, 0.0, 1000.0)) - 3000.0
    assert closed == pytest.approx(direct, rel=1e-9)


def test_straight_arrival_vertical_and_inclined():
    # Vertical photon: arrival below the emission, time = path / c.
    x, y, t = geometry_case.straight_arrival(100.0, -50.0, 1000.0, 0.0, 0.0, 0.0)
    assert x == pytest.approx(100.0) and y == pytest.approx(-50.0)
    assert t == pytest.approx(1000.0 / geometry_case.LIGHT_SPEED_M_PER_S * 1e9)
    # 45 deg zenith: path = dz / cos(45 deg), lateral shift = path sin(45 deg) = dz.
    cx = math.sin(math.radians(45.0))
    cz = math.cos(math.radians(45.0))
    x, y, t = geometry_case.straight_arrival(0.0, 0.0, 1000.0, cx, 0.0, 0.0)
    assert x == pytest.approx(1000.0, rel=1e-7)
    path = 1000.0 / cz
    assert t == pytest.approx(path / geometry_case.LIGHT_SPEED_M_PER_S * 1e9, rel=1e-7)


def test_cz_magnitude_matches_unit_vector():
    assert float(geometry_case.cz_magnitude(0.6, 0.0)) == pytest.approx(0.8)
