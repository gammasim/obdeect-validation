"""Pure helpers of the scattered-light and night-sky camera adapters (analytic expectations)."""

import math
from pathlib import Path

import numpy as np
import pytest

from obdeect_validation.adapters import night_sky_camera as nsc
from obdeect_validation.adapters import scattered_light as sl

ROOT = Path(__file__).resolve().parents[1]


def test_table_value_interpolates_one_column(tmp_path):
    table = tmp_path / "t.dat"
    table.write_text("# comment\n200 1.0 5.0\n400 3.0 7.0\n", encoding="utf-8")
    assert nsc._table_value(table, 1, 300.0) == pytest.approx(2.0)
    assert nsc._table_value(table, 2, 300.0) == pytest.approx(6.0)
    assert nsc._table_value(table, 2, 400.0) == pytest.approx(7.0)


def test_weighted_moments_of_constant_values_have_zero_error():
    w = np.array([1.0, 2.0, 3.0])
    cos = np.full(3, 0.5)
    z = np.full(3, 100.0)
    moments = sl.weighted_moments(w, cos, z)
    assert moments["cos_theta"][0] == pytest.approx(0.5)
    assert moments["cos2_theta"][0] == pytest.approx(0.25)
    assert moments["z"][0] == pytest.approx(100.0)
    assert moments["z"][1] == pytest.approx(0.0)


def test_scattered_expectation_scales_with_primary_weight_and_needs_downward_beam():
    direction = np.array([0.95, 0.0, -0.3156])
    direction /= np.linalg.norm(direction)
    origin = np.array([-1000.0, 0.0, 3000.0])
    s = np.linspace(9000.0, 11000.0, 201)
    w1 = np.ones(s.size)
    base = sl.analytic_expectation(s, w1, direction, origin)
    doubled = sl.analytic_expectation(s, 2.0 * w1, direction, origin)
    assert doubled["yield"] == pytest.approx(2.0 * base["yield"])
    assert base["yield"] > 0.0
    with pytest.raises(ValueError):
        sl.analytic_expectation(s, w1, -direction, origin)


def test_rayleigh_table_is_exponential_in_altitude():
    # Analytic Rayleigh depth tau(h) = tau0 exp(-h/8 km): beta = tau/8000 m (see module docstring).
    h_km = 2.0
    tau = sl.TAU_0 * math.exp(-h_km / 8.0)
    beta = tau / sl.SCALE_HEIGHT_M
    assert beta == pytest.approx(sl.TAU_0 * math.exp(-0.25) / 8000.0)


def test_depolarized_rayleigh_phase_is_used_by_scattering(tmp_path):
    propagation = pytest.importorskip("obdeect.propagation")
    rho = sl.DEPOLARIZATION
    gamma = rho / (2.0 - rho)
    rows = []
    for wavelength in (354, 355, 356):
        for height in (1, 60):
            tau_mol = 0.1 * math.exp(-height / 8.0)
            rows.append(f"{wavelength} {height} {tau_mol:.12g} 0 1 0")
    tables = []
    for name, depolarization in (("unpolarized", 0.0), ("depolarized", rho)):
        path = tmp_path / f"{name}.dat"
        path.write_text(
            "\n".join([*rows, f"depolarization {depolarization}"]) + "\n",
            encoding="utf-8",
        )
        tables.append(propagation.ScatteringTable(str(path)))
    transmission = propagation.TransmissionTable(
        str(ROOT / "fixtures/atm_trans_2000_1_10_0_0_2000.dat"), 2000.0
    )
    primaries = {
        "x": np.zeros(16, dtype=np.float32),
        "y": np.zeros(16, dtype=np.float32),
        "z": np.full(16, 10000.0, dtype=np.float32),
        "cx": np.zeros(16, dtype=np.float32),
        "cy": np.zeros(16, dtype=np.float32),
        "t": np.zeros(16, dtype=np.float32),
        "w": np.ones(16, dtype=np.float32),
        "wavelength": np.full(16, 355.0, dtype=np.float32),
    }
    telescope = (1, 3000.0, 0.0, 0.0, 12.5)
    outputs = [
        propagation.scatter(primaries, table, transmission, [telescope], 7, 1, 0, 16)[1]
        for table in tables
    ]
    plain, depolarized = outputs
    assert len(plain["w"]) > 0
    for key in plain:
        if key != "w":
            np.testing.assert_array_equal(plain[key], depolarized[key])
    mu = plain["cz"].astype(np.float64)
    phase_rho = (
        3.0 / (8.0 * math.pi * (2.0 + gamma)) * ((1.0 + gamma) + (1.0 - gamma) * mu**2)
    )
    phase_zero = 3.0 / (16.0 * math.pi) * (1.0 + mu**2)
    np.testing.assert_allclose(
        depolarized["w"] / plain["w"], phase_rho / phase_zero, rtol=2e-6
    )


def test_zero_extinction_table_moves_bottom_to_observation_level(tmp_path, monkeypatch):
    fixture = tmp_path / "fixture.dat"
    fixture.write_text(
        "# H2= 2.156, H1=    2.206     2.256\n200 0 0\n", encoding="utf-8"
    )
    monkeypatch.setattr(nsc, "UNIT_TRANSMISSION", fixture)
    path = nsc.zero_extinction_table(tmp_path)
    text = path.read_text(encoding="utf-8")
    assert "H2= 2.140" in text
    assert "H1=    2.206     2.256" in text
    assert "200 0 0" in text


def test_funnel_average_of_constant_table_is_the_constant(tmp_path, monkeypatch):
    mirrors = tmp_path / "mirrors.dat"
    mirrors.write_text("100 0 151 2900 3 0\n-300 200 151 2900 3 0\n", encoding="utf-8")
    table = tmp_path / "lg.dat"
    table.write_text("0 0.8\n10 0.8\n40 0.8\n", encoding="utf-8")
    monkeypatch.setattr(nsc, "MIRROR_LIST", mirrors)
    monkeypatch.setattr(nsc, "LG_TABLE", table)
    average, normal = nsc.funnel_efficiency_average()
    assert average == pytest.approx(0.8)
    assert normal == pytest.approx(0.8)


def test_funnel_average_uses_facet_geometry_and_increases_with_angle(
    tmp_path, monkeypatch
):
    if not nsc.MIRROR_LIST.exists():
        pytest.skip("LST mirror list not available")
    facets = nsc.mirror_facets()
    assert facets.shape == (198, 3)
    # Hexagon areas of the mirror list: 390.98 m^2 (list header).
    assert facets[:, 1].sum() * 1e-4 == pytest.approx(390.98, abs=0.05)
    # A table rising with angle gives an average above the normal-incidence value.
    table = tmp_path / "lg.dat"
    table.write_text("0 0.75\n10 0.80\n20 0.85\n40 0.90\n", encoding="utf-8")
    monkeypatch.setattr(nsc, "LG_TABLE", table)
    average, normal = nsc.funnel_efficiency_average()
    assert normal == pytest.approx(0.75)
    assert average > normal


def test_configuration_values_used_by_the_expectation():
    if not nsc.LST_CFG.exists():
        pytest.skip("LST configuration not available")
    assert nsc._cfg_number("mirror_degraded_reflection") == pytest.approx(0.80)
    assert nsc._cfg_number("telescope_transmission") == pytest.approx(0.9690)
    assert nsc.traced_area_cm2() == pytest.approx(377.71e4)
    assert nsc.projected_area_cm2() == pytest.approx(388.38e4)
