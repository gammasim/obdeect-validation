"""Reference distributions of the emission-profile adapter (analytic CDFs against their densities)."""

import math

import numpy as np
import pytest
from scipy import stats

from obdeect_validation.adapters import emission_profiles as ep


def integrate_cdf(pdf, lo, hi, v):
    """CDF from a density by trapezoid integration on a fine grid (independent of the adapter)."""
    x = np.linspace(lo, hi, 400_001)
    f = pdf(x)
    cumulative = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))))
    return np.interp(v, x, cumulative / cumulative[-1])


def test_lambertian_closed_form_matches_density():
    ref = ep.ref_lambertian()
    v = np.linspace(0.0, 1.0, 101)
    np.testing.assert_allclose(
        ref.cdf(v), integrate_cdf(lambda u: 1.0 - u, 0.0, 1.0, v), atol=1e-8
    )
    mean, var = ep.grid_moments(0.0, 1.0, lambda u: 1.0 - u)
    assert mean == pytest.approx(ref.mean, abs=1e-8)
    assert var == pytest.approx(ref.var, abs=1e-8)


def test_rayleigh_closed_form_matches_density():
    ref = ep.ref_rayleigh()
    v = np.linspace(0.0, 2.0, 101)
    np.testing.assert_allclose(
        ref.cdf(v),
        integrate_cdf(lambda u: 1.0 + (1.0 - u) ** 2, 0.0, 2.0, v),
        atol=1e-8,
    )
    mean, var = ep.grid_moments(0.0, 2.0, lambda u: 1.0 + (1.0 - u) ** 2)
    assert mean == pytest.approx(ref.mean, abs=1e-8)
    assert var == pytest.approx(ref.var, abs=1e-8)


def test_inverse_square_closed_form_matches_density():
    ref = ep.ref_inverse_square(300.0, 600.0)
    v = np.linspace(300.0, 600.0, 101)
    np.testing.assert_allclose(
        ref.cdf(v), integrate_cdf(lambda lam: lam**-2.0, 300.0, 600.0, v), atol=1e-8
    )
    mean, var = ep.grid_moments(300.0, 600.0, lambda lam: lam**-2.0)
    assert mean == pytest.approx(ref.mean, rel=1e-8)
    assert var == pytest.approx(ref.var, rel=1e-6)


def test_exponential_pulse_matches_density():
    ref = ep.ref_exponential_pulse(2.0, 1.5)
    v = np.linspace(0.5, 20.0, 101)
    # Integrate to 60 ns: the truncated tail beyond it is exp(-59.5 / 1.5) ~ 1e-18.
    np.testing.assert_allclose(
        ref.cdf(v),
        integrate_cdf(lambda t: np.exp(-(t - 0.5) / 1.5) / 1.5, 0.5, 60.0, v),
        atol=1e-7,
    )
    assert ref.cdf(np.array([0.4]))[0] == 0.0


def test_piecewise_linear_cdf_matches_trapezoid():
    x = [0.0, 1.0, 3.0, 4.0]
    y = [0.5, 2.0, 0.0, 1.0]
    cdf, pdf = ep.piecewise_linear_cdf(x, y)
    v = np.linspace(-1.0, 5.0, 301)
    np.testing.assert_allclose(
        cdf(v), integrate_cdf(lambda t: pdf(t), 0.0, 4.0, v), atol=1e-6
    )
    assert cdf(np.array([-1.0]))[0] == 0.0
    assert cdf(np.array([5.0]))[0] == pytest.approx(1.0)


def test_gauss_angular_reference_is_a_truncated_gaussian_in_angle():
    sigma = math.radians(5.0)
    ref = ep.ref_gauss(5.0)
    # The density in u corresponds to exp(-theta^2 / 2 sigma^2) per solid angle with
    # u = 1 - cos(theta): check the CDF against the angular integral directly.
    theta = np.linspace(0.0, min(5 * sigma, math.pi), 200_001)
    weight = np.exp(-(theta**2) / (2 * sigma**2)) * np.sin(theta)
    cumulative = np.concatenate(
        ([0.0], np.cumsum(0.5 * (weight[1:] + weight[:-1]) * np.diff(theta)))
    )
    cumulative /= cumulative[-1]
    u_at = 1.0 - np.cos(theta)
    v = np.array([0.2, 0.5, 0.8]) * ref.hi
    expected = np.interp(v, u_at, cumulative)
    np.testing.assert_allclose(ref.cdf(v), expected, atol=1e-6)


def test_complex_pulse_cdf_matches_monte_carlo():
    rng = np.random.default_rng(12345)
    n = 2_000_000
    flat, rms, decay, delay = 4.0, 1.0, 2.0, 1.0
    sample = (
        delay
        + flat * (rng.random(n) - 0.5)
        + rms * rng.standard_normal(n)
        + decay * (rng.exponential(size=n) - 1.0)
    )
    ref = ep.ref_complex_pulse(flat, rms, decay, delay)
    statistic, _ = stats.kstest(sample, ref.cdf)
    # The Monte Carlo sample has 2e6 draws: the KS statistic of a correct CDF is ~1e-3.
    assert statistic < 3e-3
    assert ref.mean == pytest.approx(float(np.mean(sample)), abs=0.01)
    assert ref.var == pytest.approx(float(np.var(sample)), rel=0.01)


def test_ks_control_rejects_isotropic_against_lambertian():
    rng = np.random.default_rng(7)
    u = 2.0 * rng.random(50_000)
    _, p_value = ep.ks_p_value(u, ep.ref_lambertian().cdf)
    assert p_value < ep.P_THRESHOLD


def test_ks_accepts_correct_distribution():
    rng = np.random.default_rng(8)
    cz = rng.uniform(-1.0, 1.0, 50_000)
    _, p_value = ep.ks_p_value(1.0 - cz, ep.ref_isotropic().cdf)
    assert p_value > 0.001


def test_azimuth_sample_is_uniform_in_turns():
    angles = np.linspace(-math.pi, math.pi, 1001)[:-1]
    cx, cy = np.cos(angles), np.sin(angles)
    phase = ep.azimuth_sample(cx, cy)
    assert np.all((phase >= 0.0) & (phase < 1.0))
    np.testing.assert_allclose(
        np.sort(phase), np.sort(np.mod(angles / (2 * math.pi), 1.0)), atol=1e-9
    )


def test_ground_expected_entry_of_a_radial_ray():
    # A photon from P = (0, 0, 40) straight down meets the sphere (R = 12.5) at z = 12.5 after
    # path 40 - 12.5 = 27.5 m; a photon pointing away (up) misses.
    sky = {
        "cx": np.array([0.0, 0.0], dtype=np.float32),
        "cy": np.array([0.0, 0.0], dtype=np.float32),
        "cz": np.array([-1.0, 1.0], dtype=np.float32),
        "t": np.array([0.0, 0.0], dtype=np.float32),
    }
    expect = ep.ground_expected(sky)
    assert expect["hit"].tolist() == [True, False]
    assert expect["s"][0] == pytest.approx(27.5, abs=1e-12)
    np.testing.assert_allclose(expect["entry"][0], [0.0, 0.0, 12.5], atol=1e-12)
