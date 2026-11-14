import numpy as np
import pytest

from obdeect_validation.adapters import source_comparison as sc


def _constant_table(n_minus_one, top_m=3000.0):
    heights = np.arange(0.0, top_m + 100.0, 100.0)
    return np.column_stack((heights, np.full(heights.size, n_minus_one)))


def test_air_integral_constant_index_is_index_times_height():
    table = _constant_table(2.0e-4)
    integral = sc.muon_air_integral_m(table, [0.0, 500.0, 3000.0])
    assert integral == pytest.approx([0.0, 0.1, 0.6], rel=1e-12, abs=1e-12)


def test_arrival_vacuum_vertical_photon_is_emission_time_plus_zero():
    # Vacuum and cz = 1: the photon reaches z = 0 at its emission time -h/c plus h/c.
    table = _constant_table(0.0)
    arrival = sc.muon_arrival_expected_ns([300.0], [1.0], table)
    assert arrival == pytest.approx([0.0], abs=1e-12)


def test_arrival_vacuum_slanted_photon_is_geometric_delay():
    # Vacuum: arrival = h (1/|cz| - 1) / c, the slant path minus the vertical path.
    table = _constant_table(0.0)
    cz = 0.99
    arrival = sc.muon_arrival_expected_ns([1000.0], [cz], table)
    expected = 1000.0 * (1.0 / cz - 1.0) / sc.C_M_NS
    assert arrival == pytest.approx([expected], rel=1e-12)


def test_arrival_constant_index_vertical_is_air_delay():
    # Vertical, constant n-1 = a: the extra time over vacuum is a h / c.
    a = 2.0e-4
    table = _constant_table(a)
    arrival = sc.muon_arrival_expected_ns([2000.0], [1.0], table)
    assert arrival == pytest.approx([a * 2000.0 / sc.C_M_NS], rel=1e-9)


def test_arrival_matches_straight_path_over_numerical_integral():
    # Independent check: integrate n(h)/|cz| along the path with a fine grid.
    heights = np.arange(0.0, 3001.0, 100.0)
    n_minus_one = 2.2e-4 - 1.9e-8 * heights
    table = np.column_stack((heights, n_minus_one))
    h, cz = 1500.0, 0.98
    z = np.linspace(0.0, h, 200001)
    path_time = (
        np.trapezoid((1.0 + np.interp(z, heights, n_minus_one)) / cz, z) / sc.C_M_NS
    )
    emission = -h / sc.C_M_NS
    arrival = sc.muon_arrival_expected_ns([h], [cz], table)
    assert arrival == pytest.approx([emission + path_time], rel=1e-7)
