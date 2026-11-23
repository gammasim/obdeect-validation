"""Unit tests of the sources-with-directions validation (isotropic to plane, upward beams).

The analytic references are checked first; the obdeect runs need obdeect.sources and the obdeect
command, and are skipped otherwise.
"""

import itertools
import math

import numpy as np
import pytest

from obdeect_validation.adapters import sources_directions

try:
    import obdeect.sources  # noqa: F401
except ImportError:
    HAVE_OBDEECT = False
else:
    HAVE_OBDEECT = True
needs_obdeect = pytest.mark.skipif(
    not HAVE_OBDEECT, reason="obdeect.sources is not installed"
)
D = sources_directions


def test_equal_probability_edges_hold_the_analytic_cdf():
    height = 100.0
    edges = D.radial_edges(height)
    assert edges[0] == 0.0 and math.isinf(edges[-1])
    cdf = 1.0 - height / np.sqrt(edges[1:-1] ** 2 + height**2)
    np.testing.assert_allclose(
        cdf, np.arange(1, D.RADIAL_BINS) / D.RADIAL_BINS, rtol=1e-12
    )


def test_annulus_fractions_sum_to_the_hemisphere():
    height = 100.0
    edges = D.radial_edges(height)
    fractions = [
        D.hemisphere_fraction(height, a, b) for a, b in itertools.pairwise(edges)
    ]
    assert sum(fractions) == pytest.approx(0.5, rel=1e-12)
    assert fractions == pytest.approx([0.5 / D.RADIAL_BINS] * D.RADIAL_BINS, rel=1e-12)


def test_radial_density_integrates_to_the_hemisphere():
    height = 100.0
    r = np.linspace(0.0, 2e5, 2_000_001)
    integrand = D.radial_density(height, 1, r) * 2.0 * np.pi * r
    integral = float(np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r)))
    # The grid stops at 2e5 m: the analytic fraction beyond it is missing from the sum.
    beyond = D.hemisphere_fraction(height, 2e5, math.inf)
    assert integral == pytest.approx(0.5 - beyond, rel=1e-6)


def test_radial_chi2_accepts_samples_from_the_analytic_distribution():
    rng = np.random.default_rng(1)
    height = 100.0
    # Inverse CDF of P(r < R | kept) = 1 - H / sqrt(R^2 + H^2).
    u = rng.uniform(size=200_000)
    r = height * np.sqrt(1.0 / (1.0 - u) ** 2 - 1.0)
    result = D.radial_test(r, height)
    assert result["chi2_ndf"] == pytest.approx(1.0, abs=0.5)
    assert sum(result["counts"]) == len(r)


def test_radial_chi2_rejects_a_wrong_height():
    rng = np.random.default_rng(2)
    u = rng.uniform(size=200_000)
    r = 100.0 * np.sqrt(1.0 / (1.0 - u) ** 2 - 1.0)
    assert D.radial_test(r, 150.0)["chi2_ndf"] > 2.0


def test_time_reference_and_ulp_of_exact_times():
    x = np.array([0.0, 30.0, 300.0], dtype=np.float32)
    y = np.zeros(3, dtype=np.float32)
    t = D.arrival_time_reference(x, y, 100.0).astype(np.float32)
    data = {"t": t, "x": x, "y": y}
    error_ns, error_ulp = D.time_test(data, 100.0)
    assert error_ulp <= 0.5 + 1e-12
    assert error_ns < 1e-4


def test_beam_direction_is_unit_and_follows_zenith():
    d = D.beam_direction(20.0, 30.0)
    assert np.linalg.norm(d) == pytest.approx(1.0, rel=1e-12)
    assert d[2] == pytest.approx(math.cos(math.radians(20.0)), rel=1e-12)
    assert D.beam_direction(0.0, 0.0)[2] == pytest.approx(1.0)


@needs_obdeect
def test_isotropic_part_keeps_half_and_balances_photons(tmp_path):
    cases = D.isotropic_part(tmp_path, photons=20_000)
    for case in cases.values():
        assert case["balance"] == 0
        assert case["run_end_delta"] == 0
        assert case["rows_delta"] == 0
        assert abs(case["kept_fraction_pull"]) < 4.0
        assert case["on_plane_z_max_abs_m"] == 0.0
        assert case["kept_rows"] + case["summary"]["not_on_plane"] == 20_000


@needs_obdeect
def test_plane_wave_rejections_name_the_problem():
    result = D.plane_wave_rejections()
    assert result["horizontal_towards_source"]["message_found"]
    assert result["start_plane_wrong_side"]["message_found"]


@needs_obdeect
def test_laser_upward_beam_matches_the_reference(tmp_path):
    result = D.laser_part(tmp_path)
    assert result["direction_error_deg"] < 1e-4
    assert result["start_error_m"] < 0.01
    assert result["start_height_error_m"] < 0.01
    assert result["time_error_ulp"] < 2
    assert result["tilt_at_limit_accepted"]
    assert result["tilt_over_rejected"]
    assert result["horizontal_rejected"]
    assert result["start_below_rejected"]
