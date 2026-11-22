"""Unit tests of the upward-lidar reference (independent of obdeect), then the obdeect adapter run.

The reference is checked first: table interpolation and its slope, disc quadrature, phase
normalisation, the trapezoid grid against scipy quad, and the time bins. The obdeect run is skipped
without obdeect.
"""

import math
import numpy as np
import pytest
from scipy.integrate import quad

from obdeect_validation.adapters import upward_lidar as ul


@pytest.fixture(scope="module")
def scat():
    return ul.Scattering()


@pytest.fixture(scope="module")
def tel():
    return ul.Transmission()


def test_table_is_the_exponential_example_at_the_nodes(scat):
    # Node at 2 km a.s.l. of the 355 nm rows: tau_mol = tau_R0 exp(-2/8), tau_aer = tau_A0 exp(-2/1.5).
    z = 2000.0 - ul.SITE_M  # altitude 2 km; the reference uses altitudes a.s.l.
    tau_r0 = 0.00865 * (355 / 1000.0) ** -4.05
    tau_a0 = 0.1 * (355 / 550.0) ** -1.3
    h = np.array([2000.0])
    assert scat._value(scat.tau_mol, h)[0] == pytest.approx(
        tau_r0 * math.exp(-2.0 / 8.0), rel=1e-6
    )
    assert scat._value(scat.tau_aer, h)[0] == pytest.approx(
        tau_a0 * math.exp(-2.0 / 1.5), rel=1e-6
    )
    assert scat.ssa[0] == pytest.approx(0.9) and scat.g[0] == pytest.approx(0.65)
    assert z == 0.0


def test_extinction_per_metre_is_the_slope_of_the_interpolant(scat):
    h = np.array([2300.0, 5400.0, 17000.0])
    eps = 1e-3
    numeric = (
        scat._value(scat.tau_mol, h + eps) - scat._value(scat.tau_mol, h - eps)
    ) / (2 * eps)
    np.testing.assert_allclose(-scat.state(h)["beta_mol"], numeric, rtol=1e-6)
    assert np.all(scat.state(h)["beta_mol"] > 0.0)


def test_no_extinction_above_the_top_node(scat):
    h = np.array([scat.top_m + 100.0])
    assert scat.state(h)["beta_mol"][0] == 0.0
    assert scat.tau(h)[0] == pytest.approx(scat.tau(np.array([scat.top_m]))[0])


def test_disc_quadrature_is_the_uniform_disc():
    x, y, w = ul.disc_nodes()
    assert w.sum() == pytest.approx(1.0, abs=1e-14)
    assert np.sum(w * (x**2 + y**2)) == pytest.approx(
        0.5, abs=1e-10
    )  # mean r^2 of a unit disc
    assert np.sum(w * x) == pytest.approx(0.0, abs=1e-14)


def test_phase_functions_are_normalised():
    integral_r = quad(lambda m: 2 * math.pi * ul.rayleigh_phase(m), -1, 1)[0]
    integral_hg = quad(
        lambda m: 2 * math.pi * ul.henyey_greenstein_phase(m, 0.65), -1, 1
    )[0]
    assert integral_r == pytest.approx(1.0, rel=1e-10)
    assert integral_hg == pytest.approx(1.0, rel=1e-10)


def test_disc_average_of_the_phase_tends_to_the_centre_value_for_a_small_disc(
    scat, tel, monkeypatch
):
    cfg = ul.CONFIGS[0]
    s = np.linspace(100.0, 400.0, 7)
    monkeypatch.setattr(
        ul, "RADIUS_M", 1e-4
    )  # same radius for both: only the aim averaging differs
    small, _, _ = ul._geometry(cfg, s, scat, tel, disc=True)
    centre, _, _ = ul._geometry(cfg, s, scat, tel, disc=False)
    np.testing.assert_allclose(small, centre, rtol=1e-7)


def test_reference_grid_agrees_with_scipy_quad(scat, tel):
    """Trapezoid on the path against adaptive quad of the same integrand (upward, 100 m)."""
    cfg = ul.CONFIGS[0]
    s_end = ul.path_end(cfg)

    def integrand(s):
        value, _, _ = ul._geometry(cfg, np.array([s]), scat, tel, disc=True)
        return float(value[0])

    # Break the path at the table nodes and the transmission nodes, as the reference does.
    z_nodes = np.concatenate([scat.h_m, tel.table.h1]) - ul.SITE_M
    vz = float(cfg.direction()[2])
    nodes = (z_nodes - cfg.start_z_m) / vz
    nodes = np.sort(nodes[(nodes > 0.0) & (nodes < s_end)])
    edges = np.concatenate([[0.0], nodes, [s_end]])
    total_quad = sum(
        quad(integrand, a, b, limit=400, epsabs=0.0, epsrel=1e-9)[0]
        for a, b in zip(edges[:-1], edges[1:])
        if b > a
    )
    total_grid = ul.reference(cfg, scat, tel, grid=200000)["total"]
    assert total_grid == pytest.approx(total_quad, rel=2e-5)


def test_time_bins_partition_the_reference(scat, tel):
    cfg = ul.CONFIGS[0]
    ref = ul.reference(cfg, scat, tel, grid=50000)
    assert np.all(ref["bins"] >= 0.0)
    assert float(np.sum(ref["bins"])) == pytest.approx(ref["total"], rel=1e-6)


def test_arrival_time_is_monotonic_along_the_path(scat, tel):
    for cfg in ul.CONFIGS:
        ref = ul.reference(cfg, scat, tel, grid=20000)
        assert np.all(np.diff(ref["arrival_m"]) >= -1e-9)


def test_chi2_counts_only_bins_with_reference_and_error():
    obs = np.array([1.0, 2.0, 3.0, 0.0])
    err = np.array([0.5, 0.5, 0.0, 0.5])
    ref = np.array([1.0, 0.0, 3.0, 0.0])
    chi2, used, pulls = ul.chi2_bins(obs, err, ref)
    assert list(used) == [0]
    assert chi2 == pytest.approx(0.0)
    assert pulls.size == 1


def test_obdeect_adapter_run_is_consistent(tmp_path, monkeypatch):
    pytest.importorskip("obdeect.propagation")
    monkeypatch.chdir(tmp_path)
    metrics, details = ul.analyse(tmp_path, n_override=40, grid=20000)
    assert metrics["sample_count"] == 3 * 40
    for key in ("up_100m", "up_1km", "down_100m"):
        assert metrics[f"{key}_plane_up_count_error"] == 0
        assert metrics[f"{key}_plane_down_count_error"] == 0
        assert metrics[f"{key}_moving_away_error"] == 0
        assert math.isfinite(metrics[f"{key}_total_yield_pull"])
    assert metrics["empty_bin_content"] == 0.0
    assert details["configurations"]["up-100m"]["bins_ndf"] >= 0
