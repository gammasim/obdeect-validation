"""Analytic and brute-force checks of the curved and plane-parallel optical depths (no obdeect)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from obdeect_validation.adapters import curved_extinction as ce

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
R_EARTH = ce.EARTH_RADIUS_M


def _fixture_table() -> ce.Table:
    return ce.parse_table((FIXTURES / "atm_trans_2000_1_10_0_0_2000.dat").read_text())


def _exponential_table(
    tau0=1.0, scale_m=8000.0, h0=2156.0, top_m=5.0e5, n=4000
) -> ce.Table:
    """tau(h) = tau0 (1 - exp(-(h - h0)/H)) on a fine log-spaced table (all wavelengths equal)."""
    # Nodes geometrically spaced in height above H2, so that the first segment is 1 m wide
    # (a wide first segment gives a log-linear error near the telescope of about 1e-5).
    heights = h0 + np.geomspace(1.0, top_m - h0, n)
    tau = tau0 * (1.0 - np.exp(-(heights - h0) / scale_m))
    return ce.Table(
        h2_m=h0,
        heights_m=heights,
        wavelengths=np.array([300]),
        tau=tau[None, :],
    )


def _exponential_tau(h, tau0=1.0, scale_m=8000.0, h0=2156.0):
    return tau0 * (1.0 - np.exp(-(np.asarray(h) - h0) / scale_m))


def _chapman_reference(h_t, zem, c, tau0=1.0, scale_m=8000.0, h0=2156.0, nodes=4000):
    """Slant depth for the exponential tau: integral of (dtau/dh) / cos theta(h) over [h_t, zem].

    Gauss-Legendre in y = (h - h_t) / H on [0, (zem - h_t) / H]; independent of the table code.
    The integrand is the Chapman-type integral of an exponential atmosphere.
    """
    x, w = np.polynomial.legendre.leggauss(nodes)
    y_max = (zem - h_t) / scale_m
    y = 0.5 * y_max * (x + 1.0)
    weights = 0.5 * y_max * w
    h = h_t + scale_m * y
    p = (R_EARTH + h_t) * math.sqrt(1.0 - c * c)
    cos_h = np.sqrt(1.0 - (p / (R_EARTH + h)) ** 2)
    density = tau0 / scale_m * np.exp(-(h - h0) / scale_m)
    return float(np.sum(weights * scale_m * density / cos_h))


def test_parse_fixture_table_structure():
    table = _fixture_table()
    text = (FIXTURES / "atm_trans_2000_1_10_0_0_2000.dat").read_text()
    header = next(line for line in text.splitlines() if "H1=" in line)
    assert table.h2_m == pytest.approx(2000.0)
    assert table.heights_m.size == len(header.split("H1=")[1].split())
    assert table.wavelengths[0] == 200 and table.wavelengths[-1] == 999
    assert table.tau.shape == (800, table.heights_m.size)
    assert np.all(np.diff(table.heights_m) > 0.0)


def test_parse_rejects_bad_tables():
    with pytest.raises(ValueError):
        ce.parse_table("# no header\n200 1.0\n")
    text = "# H2= 2.000, H1=    2.050     2.100\n200 1.0 2.0\n202 1.0 2.0\n"
    with pytest.raises(ValueError):
        ce.parse_table(text)


def test_vertical_depth_zero_at_and_below_h2_and_constant_above_top():
    table = _fixture_table()
    below = ce.vertical_depth(table, [0.0, 1999.0, 2000.0])
    assert np.all(below == 0.0)
    top = ce.vertical_depth(table, [table.heights_m[-1], 2.0 * table.heights_m[-1]])
    assert np.allclose(top[:, 0], table.tau[:, -1])
    assert np.allclose(top[:, 1], table.tau[:, -1])


def test_vertical_depth_log_linear_between_heights():
    table = _fixture_table()
    h1, h2 = table.heights_m[3], table.heights_m[4]
    mid = math.sqrt(h1 * h2)  # midpoint in log h
    value = ce.vertical_depth(table, [mid])[:, 0]
    assert np.allclose(value, 0.5 * (table.tau[:, 3] + table.tau[:, 4]), rtol=1e-12)


def test_radius_to_infinity_equals_plane_parallel():
    table = _fixture_table()
    zem = np.array([5000.0, 9000.0, 12000.0])
    c = np.array([0.5, 0.2, 0.9])
    h_t = 2167.0
    curved = ce.curved_optical_depth(table, zem, c, h_t, radius_m=1.0e12)
    plane = ce.plane_optical_depth(table, zem, c, h_t)
    assert np.allclose(curved, plane, rtol=1e-6, atol=0.0)


def test_vertical_incidence_equals_plane_parallel():
    table = _fixture_table()
    zem = np.array([3000.0, 7000.0, 15000.0])
    h_t = 2167.0
    curved = ce.curved_optical_depth(table, zem, np.ones(3), h_t)
    plane = ce.plane_optical_depth(table, zem, np.ones(3), h_t)
    assert np.allclose(curved, plane, rtol=1e-6, atol=1e-12)


def test_exponential_table_flat_limit_closed_form():
    table = _exponential_table()
    h_t, zem = 2167.0, 12000.0
    for c in (0.6, 0.3):
        slant = ce.curved_optical_depth(table, [zem], [c], h_t, radius_m=1.0e12)[0, 0]
        closed = (_exponential_tau(zem) - _exponential_tau(h_t)) / c
        assert slant == pytest.approx(closed, rel=1e-6)


def test_chapman_numeric_for_real_radius():
    """Real R: the slant depth matches the Chapman-type integral of the exponential atmosphere."""
    table = _exponential_table(n=6000)
    h_t = 2167.0
    zem = (
        h_t + 40.0 * 8000.0
    )  # e^-40 beyond the emission range: effectively the full ray
    for c in (0.5, 0.2, 0.1):
        slant = ce.curved_optical_depth(table, [zem], [c], h_t)[0, 0]
        reference = _chapman_reference(h_t, zem, c)
        assert slant == pytest.approx(reference, rel=1e-4)


def test_slant_depth_converged_in_quadrature():
    table = _exponential_table()
    h_t, zem = 2167.0, 20000.0
    coarse = ce.curved_optical_depth(table, [zem], [0.15], h_t, intervals=64)
    fine = ce.curved_optical_depth(table, [zem], [0.15], h_t, intervals=1024)
    assert np.allclose(coarse, fine, rtol=1e-9)


def test_curved_is_bounded_by_plane_and_finite_at_horizon():
    table = _fixture_table()
    zem = np.full(6, 9000.0)
    c = np.array([0.9, 0.5, 0.2, 0.1, 0.01, 0.0])
    h_t = 2167.0
    curved = ce.curved_optical_depth(table, zem, c, h_t)
    assert curved.shape == (6, 800)
    assert np.all(np.isfinite(curved))
    # cz = 0 is finite (grazing photons); the plane-parallel airmass diverges there.
    plane = ce.plane_optical_depth(table, zem[:-1], c[:-1], h_t)
    assert np.all(curved[:-1] <= plane * (1.0 + 1e-12))


def test_slant_depth_monotonic_in_emission_height_and_zenith():
    table = _exponential_table()
    h_t = 2167.0
    zem = np.linspace(3000.0, 60000.0, 25)
    values = ce.curved_optical_depth(table, zem, np.full(zem.size, 0.3), h_t)[:, 0]
    assert np.all(np.diff(values) >= 0.0)
    c = np.linspace(0.05, 1.0, 20)
    by_zenith = ce.curved_optical_depth(table, np.full(c.size, 20000.0), c, h_t)[:, 0]
    assert np.all(np.diff(by_zenith) <= 0.0)


def test_transmission_pair_ordering():
    table = _fixture_table()
    zem = np.array([8000.0, 12000.0])
    c = np.array([0.25, 0.1])
    curved, plane = ce.transmission_pair(table, zem, c, 2167.0, pathlen_m=28.0)
    assert np.all(curved >= plane)
    assert np.all((curved >= 0.0) & (curved <= 1.0))


def test_ratio_helper():
    ratio, error = ce.ratio_with_error(2.0, 0.02, 1.0, 0.0)
    assert ratio == pytest.approx(2.0)
    assert error == pytest.approx(0.02)  # 1 % relative error of the numerator
    ratio, error = ce.ratio_with_error(3.0, 0.3, 3.0, 0.3)
    assert ratio == pytest.approx(1.0)
    assert error == pytest.approx(math.sqrt(2) * 0.1)
    with pytest.raises(ValueError):
        ce.ratio_with_error(0.0, 0.1, 1.0, 0.1)


def test_band_bins_probabilities_normalised_and_bounded():
    n, p = ce.band_bins((240.5, 1000.0))
    assert n[0] == 240 and n[-1] == 999
    assert p.sum() == pytest.approx(1.0, rel=1e-12)
    assert np.all(p > 0.0)
    n, p = ce.band_bins((300.0, 300.0 + 1.0))
    assert list(n) == [300] and p[0] == pytest.approx(1.0)


def test_band_expectation_matches_brute_force_sampling():
    """Exact bin probabilities against sampled lambda (1/lambda^2, floor to nm) per bunch."""
    table = _fixture_table()
    band = (240.0, 900.0)
    zem = np.array([6000.0, 11000.0, 17500.0])
    c = np.array([0.45, 0.2, 0.7])
    h_t = 2167.0
    expected = ce.band_expectation(table, band, zem, c, h_t, pathlen_m=28.0)

    # Transmission of every table wavelength (rows 200..999); the draws index into it.
    curved_all, plane_all = ce.transmission_pair(table, zem, c, h_t, 28.0)
    rng = np.random.default_rng(7)
    draws = 200000
    u = rng.random(draws)
    wavelength = np.floor(1.0 / (1.0 / band[0] - u * (1.0 / band[0] - 1.0 / band[1])))
    column = (wavelength - table.wavelengths[0]).astype(int)
    for b in range(zem.size):
        values_c = curved_all[b, column]
        values_p = plane_all[b, column]
        se_c = values_c.std(ddof=1) / math.sqrt(draws)
        se_p = values_p.std(ddof=1) / math.sqrt(draws)
        assert abs(values_c.mean() - expected["curved"][b]) < 5.0 * se_c
        assert abs(values_p.mean() - expected["plane"][b]) < 5.0 * se_p


def test_weighted_band_expectation_matches_sampling():
    """Detection-weighted mean sum(P w T) / sum(P w) against sampled photons with acceptance w(lambda)."""
    table = _fixture_table()
    band = (240.0, 900.0)
    zem = np.array([6000.0, 13000.0])
    c = np.array([0.5, 0.3])
    h_t = 2167.0

    def weight(nm):
        # A smooth, strongly wavelength-dependent acceptance (zero above 800 nm).
        nm = np.asarray(nm, dtype=float)
        return np.clip((800.0 - nm) / 500.0, 0.0, 1.0) * (
            0.2 + 0.8 * np.exp(-(((nm - 380.0) / 60.0) ** 2))
        )

    expected = ce.band_expectation(
        table, band, zem, c, h_t, pathlen_m=28.0, weight=weight
    )
    curved_all, _plane_all = ce.transmission_pair(table, zem, c, h_t, 28.0)
    rng = np.random.default_rng(11)
    draws = 400000
    u = rng.random(draws)
    wavelength = np.floor(1.0 / (1.0 / band[0] - u * (1.0 / band[0] - 1.0 / band[1])))
    column = (wavelength - table.wavelengths[0]).astype(int)
    acceptance = weight(wavelength)
    for b in range(zem.size):
        sampled = curved_all[b, column] * acceptance
        mean = sampled.mean() / acceptance.mean()
        # Ratio estimator: standard error from the sample of w*T and w (first order).
        se = sampled.std(ddof=1) / math.sqrt(draws) / acceptance.mean()
        assert abs(mean - expected["curved"][b]) < 5.0 * se + 1e-6
    assert expected["live_fraction"] == pytest.approx(1.0)
    unit = ce.unit_variance(band, weight)
    sampled_w = acceptance
    assert unit == pytest.approx(
        sampled_w.var(ddof=0) / sampled_w.mean() ** 2, rel=2e-2
    )


def test_band_expectation_chunking_invariant_and_variance():
    table = _fixture_table()
    zem = np.linspace(6000.0, 15000.0, 9)
    c = np.linspace(0.2, 0.9, 9)
    whole = ce.band_expectation(table, (250.0, 800.0), zem, c, 2167.0, 28.0)
    pieces = ce.band_expectation(table, (250.0, 800.0), zem, c, 2167.0, 28.0, chunk=2)
    for key in whole:
        assert np.allclose(whole[key], pieces[key], rtol=1e-12, atol=0.0)
    assert np.all(whole["curved_variance"] >= 0.0)
    assert np.all(whole["curved"] >= whole["plane"])


def test_band_above_table_limit_has_zero_transmission():
    table = _fixture_table()
    out = ce.band_expectation(table, (1000.0, 1200.0), [9000.0], [0.4], 2167.0)
    assert out["curved"][0] == 0.0
    assert out["plane"][0] == 0.0


def test_weighted_helpers():
    photons = np.array([10.0, 30.0])
    assert ce.weighted_mean([0.5, 0.25], photons) == pytest.approx(0.3125)
    se = ce.weighted_standard_error([0.04, 0.09], photons)
    assert se == pytest.approx(math.sqrt(100 * 0.04 + 900 * 0.09) / 40.0)


def test_emission_at_or_below_telescope_rejected():
    table = _fixture_table()
    with pytest.raises(ValueError):
        ce.curved_optical_depth(table, [2100.0], [0.5], 2167.0)
    with pytest.raises(ValueError):
        ce.plane_optical_depth(table, [9000.0], [0.0], 2167.0)
