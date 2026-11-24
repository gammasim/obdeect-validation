"""Unit tests of the sources-targets validation: reference fractions and the statistics.

The closed forms are checked against direct quadrature first; the obdeect runs need obdeect.sources
and are skipped otherwise.
"""

import math

import pytest
from scipy import integrate

from obdeect_validation.adapters import sources_targets as T

try:
    import obdeect.sources  # noqa: F401
except ImportError:
    HAVE_OBDEECT = False
else:
    HAVE_OBDEECT = True
needs_obdeect = pytest.mark.skipif(
    not HAVE_OBDEECT, reason="obdeect.sources is not installed"
)


def quad_cap(intensity, alpha, theta_max=math.pi):
    """Direct quadrature of the emitted light inside a cone of half-angle alpha (on axis)."""
    total = integrate.quad(lambda t: intensity(t) * math.sin(t), 0.0, theta_max)[0]
    part = integrate.quad(
        lambda t: intensity(t) * math.sin(t), 0.0, min(alpha, theta_max)
    )[0]
    return part / total


@pytest.mark.parametrize("ratio", T.DISC_RATIOS)
def test_closed_forms_match_quadrature(ratio):
    alpha = math.atan(ratio)
    assert T.cap_fraction("isotropic", alpha) == pytest.approx(
        quad_cap(lambda t: 1.0, alpha), rel=1e-9
    )
    assert T.cap_fraction("lambertian", alpha) == pytest.approx(
        quad_cap(lambda t: math.cos(t), alpha, 0.5 * math.pi), rel=1e-9
    )
    for name, deg in T.CONE_DEG.items():
        theta_c = math.radians(deg)
        # Step profile: integrate over the cone itself (the step is at theta_c).
        part = integrate.quad(math.sin, 0.0, min(alpha, theta_c))[0]
        total = integrate.quad(math.sin, 0.0, theta_c)[0]
        assert T.cap_fraction(name, alpha) == pytest.approx(part / total, rel=1e-12)


def test_lambertian_fraction_is_radius_ratio():
    for ratio in T.DISC_RATIOS:
        alpha = math.atan(ratio)
        assert T.cap_fraction("lambertian", alpha) == pytest.approx(
            ratio**2 / (1.0 + ratio**2), rel=1e-12
        )


def test_gauss_fraction_is_a_quadrature_of_the_truncated_profile():
    sigma = math.radians(T.GAUSS_SIGMA_DEG)
    theta_max = 5.0 * sigma
    alpha = math.atan(1.0)
    expected = quad_cap(lambda t: math.exp(-0.5 * (t / sigma) ** 2), alpha, theta_max)
    assert T.cap_fraction("gauss_5deg", alpha) == pytest.approx(expected, rel=1e-9)
    assert T.cap_fraction("gauss_5deg", math.pi) == pytest.approx(1.0, rel=1e-12)


def test_sphere_fraction_matches_isotropic_cone_quadrature():
    for factor in T.SPHERE_DISTANCES:
        cos_max = math.sqrt(1.0 - factor**-2)
        expected = 0.5 * (1.0 - cos_max)
        assert T.sphere_fraction(factor) == pytest.approx(expected, rel=1e-12)


def test_tilted_cones_miss_every_disc():
    for ratio in T.DISC_RATIOS:
        alpha = math.degrees(math.atan(ratio))
        for name in T.TILTED_PROFILES:
            assert alpha + T.CONE_DEG[name] < T.TILT_DEG


def test_pull_and_sigma_conventions():
    assert T.pull(10.0, 12.0, 1.0) == pytest.approx(2.0)
    assert T.pull(3.0, 3.0, 0.0) == 0.0
    assert math.isinf(T.pull(3.0, 4.0, 0.0))
    assert T.count_sigma(5.0) == 0.5
    assert T.count_sigma(5.25) == 0.5  # rounding bound, independent of the fraction


def test_weight_sigma_counts_rounding_for_unit_weights():
    rows = 1000
    t = {"rows": rows, "w_sum": float(rows), "w2_sum": float(rows)}
    assert T.weight_sigma(t) == pytest.approx(
        0.5
    )  # constant weights: rounding of one row


def test_threshold_pair_straddles_the_automatic_target():
    lo, hi = T.THRESHOLD_PAIR
    assert lo < 1.02 < hi


@needs_obdeect
def test_threshold_check_is_per_side_against_own_light():
    # Each side agrees with its own analytic sphere fraction (the light is continuous across 1.02 R);
    # the two sides do not agree with each other (their true light differs: the old check was wrong).
    res = T.threshold_jump(T.SEED + 1000)
    assert res["violations"] == 0
    assert abs(res["pull_below"]) <= T.GATE_PULL
    assert abs(res["pull_above"]) <= T.GATE_PULL
    light_difference = res["sphere_light_below"] - res["sphere_light_above"]
    assert light_difference > 100.0  # F(1.019 R) - F(1.021 R) = 0.0047 of the photons
    assert res["difference_diagnostic"] == pytest.approx(light_difference, abs=5.0)


def test_disc_configs_cover_tilted_and_on_axis_cases():
    configs = T.disc_configs()
    tilted = [c for c in configs if c["id"].endswith("-tilted")]
    assert len(tilted) == len(T.TILTED_PROFILES) * len(T.DISC_RATIOS)
    assert all(c["fraction"] == 0.0 for c in tilted)
    on_axis = [c for c in configs if not c["id"].endswith("-tilted")]
    assert all(c["fraction"] > 0.0 for c in on_axis)


@needs_obdeect
def test_small_disc_emission_has_exact_expected_light():
    # R/h = 0.2 with an isotropic point: the weight sum must be N F within 5 sigma.
    config = next(c for c in T.disc_configs() if c["id"] == "disc-isotropic-R_h0.2")
    result = T.evaluate_config(config, seed=7)
    assert abs(result["weight_pull"]) < T.GATE_PULL
    assert result["zero_weight_rows"] == 0
    assert result["expected_rows_api"] == pytest.approx(result["rows"], abs=1.0)


@needs_obdeect
def test_sphere_emission_writes_one_row_per_hit_photon():
    config = T.sphere_configs()[0]
    result = T.evaluate_config(config, seed=8)
    assert result["rows"] == pytest.approx(result["expected_light"], abs=1.0)
    assert result["rows_per_hit"] == pytest.approx(1.0, abs=1e-5)


@needs_obdeect
def test_isotropic_disc_at_large_ratio_writes_one_row_per_hit():
    # Rejection sampling of the profile (docs/sources/spec/source-core.md): a disc target writes
    # one row of weight 1 per photon that reaches the disc, at any R/h. (Uniform-area sampling,
    # which wrote rows ~ tan^2(alpha) / (2 (1 - cos alpha)) per hit, is no longer used.)
    config = next(c for c in T.disc_configs() if c["id"] == "disc-isotropic-R_h1")
    result = T.evaluate_config(config, seed=9)
    # rows / (N F) fluctuates with the count of hits: 5 sigma of a Poisson count with mean N F.
    assert result["rows_per_hit"] == pytest.approx(
        1.0, abs=5.0 / math.sqrt(result["expected_light"])
    )


@needs_obdeect
def test_nonisotropic_sphere_is_rejected():
    assert T.nonisotropic_sphere_accepted() == 0
