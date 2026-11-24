"""Unit tests of the high-zenith refraction reference and gates (no obdeect run needed)."""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from obdeect_validation.adapters import refraction_high_zenith as rhz
from obdeect_validation.adapters import upward_common as uc


@pytest.fixture(scope="module")
def prof():
    return uc.profile()


def test_gates_follow_the_docstring():
    assert rhz.gates(60.0) == (1e-6, 1e-6)
    assert rhz.gates(82.0) == (1e-6, 1e-6)
    assert rhz.gates(85.0) == (4e-6, 1.4e-5)
    assert rhz.gates(87.0) == (4e-6, 1.4e-5)


def test_vertical_ray_has_no_displacement(prof):
    dx, excess_ns, theta_obs = rhz.ray_reference(prof, 10000.0, 0.0)
    assert dx == 0.0
    assert theta_obs == 0.0
    # vertical: excess = integral of nu dh / c over the path
    nu_integral = quad(lambda h: float(prof.nu(h)), uc.SITE_M, 10000.0, limit=800)[0]
    assert excess_ns == pytest.approx(nu_integral / uc.C_M_NS, rel=1e-9)


@pytest.mark.parametrize(
    ("height", "zenith"), [(5000.0, 60.0), (20000.0, 80.0), (10000.0, 87.0)]
)
def test_ray_reference_matches_quadrature(prof, height, zenith):
    """The solve_ivp reference agrees with scipy quad of the Snell integrals (both use one profile)."""
    s = math.sin(math.radians(zenith))
    dx, excess, _ = rhz.ray_reference(prof, height, s)
    quad_ref = uc.snell_integrals(prof, height, s, uc.SITE_M)
    cos_e = math.sqrt(1 - s * s)
    dx_quad = quad_ref["horizontal"] - (height - uc.SITE_M) * s / cos_e
    t_quad = quad_ref["time_ns"] - (height - uc.SITE_M) / (uc.C_M_NS * cos_e)
    assert dx == pytest.approx(dx_quad, abs=1e-5)
    assert excess == pytest.approx(t_quad, abs=1e-5)


def test_bending_is_towards_the_vertical(prof):
    """Snell: the exit direction is steeper than the emission direction, and the displacement negative."""
    s = math.sin(math.radians(70.0))
    dx, _, theta_obs = rhz.ray_reference(prof, 10000.0, s)
    assert dx < 0.0
    assert math.sin(theta_obs) < s


def test_total_time_adds_vacuum_path():
    s = math.sin(math.radians(60.0))
    cos_e = math.sqrt(1 - s * s)
    total = rhz.total_time_ns(10000.0, s, 2.0)
    assert total == pytest.approx(
        (10000.0 - uc.SITE_M) / (uc.C_M_NS * cos_e) + 2.0, rel=1e-14
    )


def test_build_rays_stores_float32_directions():
    rays = rhz.build_rays()
    assert len(rays) == len(rhz.HEIGHTS_M) * len(rhz.ZENITHS_DEG)
    for ray in rays:
        assert ray["cx"] == float(np.float32(ray["cx"]))
        assert ray["sin_e"] == pytest.approx(
            math.sin(math.radians(ray["zenith"])), abs=1e-6
        )
        assert ray["cz_doc"] == pytest.approx(math.sqrt(1 - ray["sin_e"] ** 2))


def test_profile_grid_is_the_obdeect_model(prof):
    """Inside a grid cell ln(n - 1) is linear, the model of obdeect's exponential fine cells."""
    g = prof.grid
    k = len(g) // 2
    mid = 0.5 * (g[k] + g[k + 1])
    slope = (prof.lnnu[k + 1] - prof.lnnu[k]) / (g[k + 1] - g[k])
    assert float(prof.nu(mid)) == pytest.approx(
        math.exp(prof.lnnu[k] + slope * (mid - g[k])), rel=1e-12
    )
