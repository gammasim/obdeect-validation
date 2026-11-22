"""Unit tests of the upward-photon propagation references (independent of obdeect), then the obdeect checks.

The references are checked first: profile interpolation, transmission table, Snell integrals against a
closed form, critical angle, slant depth against the plane-parallel limit. Adapter runs with obdeect are
skipped without obdeect.
"""

import math

import numpy as np
import pytest

from obdeect_validation.adapters import direction_losses as dl
from obdeect_validation.adapters import upward_common as uc
from obdeect_validation.adapters import upward_extinction as uext
from obdeect_validation.adapters import upward_geometry as ug


@pytest.fixture(scope="module")
def prof():
    return uc.Profile()


@pytest.fixture(scope="module")
def tab():
    return uc.Transmission()


def test_profile_reproduces_the_levels_and_nodes_are_at_most_3_m_apart(prof):
    level_values = np.log(prof.nu(prof.altitude[1:-1]))
    table = np.loadtxt(uc.PROFILE, comments="#", usecols=(0, 3))
    np.testing.assert_allclose(
        np.log(prof.nu(table[1:-1, 0] * 1000.0)), np.log(table[1:-1, 1]), rtol=1e-10
    )
    assert np.max(np.diff(prof.grid)) <= 3.0 + 1e-9
    assert level_values.size > 0


def test_profile_is_zero_above_the_table(prof):
    assert float(prof.nu(prof.top + 1.0)) == 0.0


def test_transmission_is_zero_at_h2_and_matches_the_table_at_h1(tab):
    assert float(tab.tau(500, tab.h2)) == 0.0
    row = tab.rows[500]
    assert float(tab.tau(500, tab.h1[3])) == pytest.approx(row[3], rel=1e-12)


def test_transmission_is_log10_linear_between_table_heights(tab):
    h_a, h_b = tab.h1[0], tab.h1[1]
    mid = 10 ** (0.5 * (math.log10(h_a) + math.log10(h_b)))
    expected = 0.5 * (tab.rows[500][0] + tab.rows[500][1])
    assert float(tab.tau(500, mid)) == pytest.approx(expected, rel=1e-12)


def test_wavelength_is_linear_between_rows(tab):
    a, b = float(tab.tau(400, 3000.0)), float(tab.tau(401, 3000.0))
    assert float(tab.tau(400.25, 3000.0)) == pytest.approx(
        0.75 * a + 0.25 * b, rel=1e-12
    )


def test_dtau_dh_matches_a_numerical_derivative(tab):
    h = 3300.0
    eps = 1e-3
    numeric = (float(tab.tau(500, h + eps)) - float(tab.tau(500, h - eps))) / (2 * eps)
    assert tab.dtau_dh(500, h) == pytest.approx(numeric, rel=1e-6)


def test_snell_without_refraction_is_the_straight_line(prof):
    # A 1 m path: the ray is straight to within the index gradient over 1 m.
    h_a, h_b = 2000.0, 2001.0
    sin_a = 0.5
    res = uc.snell_integrals(prof, h_a, sin_a, h_b)
    tan_a = sin_a / math.sqrt(1 - sin_a**2)
    assert res["horizontal"] == pytest.approx(tan_a * (h_b - h_a), rel=1e-6)
    # The time integral carries the index n = 1 + nu, evaluated at the mid height over the 1 m interval.
    n_mid = 1.0 + float(prof.nu(0.5 * (h_a + h_b)))
    expected = n_mid * (h_b - h_a) / math.sqrt(1 - sin_a**2) / uc.C_M_NS
    assert res["time_ns"] == pytest.approx(expected, rel=1e-7)


def test_critical_sin_gives_sin_one_at_the_plane(prof):
    h_e, h_t = uc.SITE_M + 100.0, uc.SITE_M + 500.0
    s = uc.critical_sin(prof, h_e, h_t)
    res = uc.snell_integrals(prof, h_e, s, h_t)
    assert res["sin_b"] == pytest.approx(1.0, abs=1e-12)


def test_turning_ray_is_reported_as_turning(prof):
    h_e, h_t = uc.SITE_M + 100.0, uc.SITE_M + 500.0
    s = uc.critical_sin(prof, h_e, h_t) + 0.5 * (1 - uc.critical_sin(prof, h_e, h_t))
    res = uc.snell_integrals(prof, h_e, s, h_t)
    assert res["turns"] is True


def test_stored_direction_is_unit_to_float32_rounding():
    cx, cy, cz = uc.stored_direction(87.0, 130.0)
    norm = float(cx) ** 2 + float(cy) ** 2 + float(cz) ** 2
    assert norm == pytest.approx(1.0, abs=2e-7)


def test_straight_line_geometry_reference_for_a_vertical_photon():
    data = {
        "x": np.array([2.0], dtype=np.float32),
        "y": np.array([-1.0], dtype=np.float32),
        "z": np.array([100.0], dtype=np.float32),
        "cx": np.array([0.0], dtype=np.float32),
        "cy": np.array([0.0], dtype=np.float32),
        "cz": np.array([1.0], dtype=np.float32),
        "t": np.array([3.0], dtype=np.float32),
    }
    ref = ug.reference(data)
    assert ref["x"][0] == pytest.approx(2.0)
    assert ref["t"][0] == pytest.approx(3.0 + 400.0 / uc.C_M_NS)


def test_spherical_slant_depth_tends_to_the_plane_value_for_a_large_radius(tab):
    rng = np.random.default_rng(3)
    cols = uext.upward_sample(rng, 6)
    plane = uext.plane_reference(tab, cols)
    flat = uext.spherical_reference(tab, cols, uext.FLAT_RADIUS_M)
    np.testing.assert_allclose(flat, plane, atol=1e-6)


def test_direction_population_has_the_expected_categories():
    _cols, names, expected, turn_sin = dl.population()
    assert len(names) == 8 * dl.N
    assert 0.0 < turn_sin < 1.0
    reasons = [e for e in expected if e is not None]
    assert len(reasons) == 6 * dl.N
    assert reasons.count("moving_away") == 2 * dl.N
    assert reasons.count("below_ground") == 2 * dl.N


def test_obdeect_geometry_adapter_runs(tmp_path):
    pytest.importorskip("obdeect.propagation")
    metrics, _ = ug.analyse(tmp_path)
    assert metrics["sample_count"] == 140
    assert metrics["dropped"] == 0
    assert metrics["max_direction_error"] == 0.0
    assert metrics["max_position_error_ulp"] <= 1.0
