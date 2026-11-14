import numpy as np
import pytest

from obdeect_validation.adapters import source_comparison as sc


def _unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v, axis=1)[:, None]


def test_closest_approach_is_the_same_for_every_point_on_the_ray():
    rng = np.random.default_rng(1)
    d = _unit(rng.normal(size=(50, 3)))
    p = rng.normal(size=(50, 3)) * 10.0
    steps = rng.uniform(-100.0, 100.0, 50)
    q = p + steps[:, None] * d
    assert sc.closest_approach(q, d) == pytest.approx(
        sc.closest_approach(p, d), abs=1e-9
    )


def test_closest_approach_equals_distance_from_origin_to_line():
    # Ray through (3, 4, 0) along z: closest approach 5 m, at (3, 4, z) for any z.
    d = np.array([[0.0, 0.0, 1.0]])
    for z in (-7.0, 0.0, 12.0):
        assert sc.closest_approach(np.array([[3.0, 4.0, z]]), d) == pytest.approx([5.0])


def test_closest_approach_matches_impact_parameter_helper():
    rng = np.random.default_rng(2)
    d = _unit(rng.normal(size=(20, 3)))
    p = rng.normal(size=(20, 3))
    assert sc.closest_approach(p, d) == pytest.approx(
        sc.impact_parameter(p, d), abs=1e-12
    )


def test_centre_plane_time_of_sphere_point_is_flat_for_declared_geometry():
    # Native emission: entry point P on the hemisphere facing -d, time -zi R / c. Then the centre-plane
    # time is the emission time for every ray (the pulse is flat there).
    rng = np.random.default_rng(3)
    radius = 15.0
    d = _unit(
        np.column_stack(
            [rng.normal(0, 0.1, 100), rng.normal(0, 0.1, 100), -np.ones(100)]
        )
    )
    # Points on the sphere with P.d < 0 (entry side): P = R * unit(-d * zi + perpendicular part).
    perp = rng.normal(size=(100, 3))
    perp -= np.sum(perp * d, axis=1)[:, None] * d
    perp = _unit(perp)
    zi = rng.uniform(0.0, 1.0, 100)
    p_unit = zi[:, None] * (-d) + np.sqrt(1 - zi**2)[:, None] * perp
    p = radius * p_unit
    t = -zi * radius / sc.C_M_NS
    assert sc.centre_plane_time(t, p, d) == pytest.approx(np.zeros(100), abs=1e-12)


def test_centre_plane_time_is_the_same_on_the_observation_plane_crossing():
    # The arrival time at the observation plane, transformed with a different point of the same ray,
    # must give the same centre-plane time (reference and native use one definition).
    d = np.array([[0.1, 0.0, -np.sqrt(1 - 0.01)]])
    p = np.array([[10.0, 0.5, 2.0]])
    speed = 0.3
    t_p = np.array([1.0])
    s = -p[:, 2] / d[:, 2]  # distance along d to z = 0
    q = p + s[:, None] * d
    t_q = t_p + s / speed
    centre_p = sc.centre_plane_time(t_p, p, d, speed_m_per_ns=speed)
    centre_q = sc.centre_plane_time(t_q, q, d, speed_m_per_ns=speed)
    assert centre_q == pytest.approx(centre_p, abs=1e-12)


def test_projected_area_expectation_for_uniform_disc_is_b_squared_uniform():
    # Isotropic flux inside a cone: for a fixed direction the flux is uniform over the disc
    # perpendicular to d, so (b/R)^2 is uniform on [0, 1].
    rng = np.random.default_rng(4)
    n = 20000
    b = np.sqrt(rng.random(n))
    stat = sc.ks_one_sample(b**2, lambda x: x)
    assert stat < sc.ks_critical(n)


def test_impact_parameter_cdf_of_area_uniform_sampler_is_r_squared():
    # Sampler ri = sqrt(U) (obdeect sources spec): b/R = ri has CDF r^2 (the analytic check of
    # impact_parameter_native_vs_sampler). The old ri = 1 - sqrt(U) would give 2r - r^2.
    rng = np.random.default_rng(5)
    n = 20000
    ri = np.sqrt(rng.random(n))
    assert sc.ks_one_sample(ri, lambda r: r * r) < sc.ks_critical(n)
    assert sc.ks_one_sample(ri, lambda r: 2 * r - r * r) > sc.ks_critical(n)
