"""Independent checks for wavelength interpolation and dispersive air delay."""

from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCATTERING = ROOT / "fixtures/scattering_exponential.dat"


def _ciddor_refractivity(wavelength_nm):
    """Ciddor standard dry-air refractivity, independent of the obdeect implementation."""
    s2 = (1000.0 / wavelength_nm) ** 2
    return 1e-8 * (5792105.0 / (238.0185 - s2) + 167917.0 / (57.362 - s2))


def _ciddor_group_refractivity(wavelength_nm):
    s2 = (1000.0 / wavelength_nm) ** 2
    phase = _ciddor_refractivity(wavelength_nm)
    derivative_term = (
        2.0 * s2 * (5792105.0 / (238.0185 - s2) ** 2 + 167917.0 / (57.362 - s2) ** 2)
    )
    return phase + 1e-8 * derivative_term


def test_fractional_scattering_wavelength_interpolates_adjacent_rows():
    ScatteringTable = pytest.importorskip("obdeect.propagation").ScatteringTable

    table = ScatteringTable(str(SCATTERING))
    altitude = np.array([5000.0, 5000.0, 5000.0])
    wavelengths = np.array([354.0, 354.5, 355.0])
    tau_mol, tau_aer, ssa, asymmetry = table.profile(wavelengths, altitude)
    lower = table.profile(np.array([354.0]), np.array([5000.0]))
    upper = table.profile(np.array([355.0]), np.array([5000.0]))
    for values, lo, hi in zip((tau_mol, tau_aer, ssa, asymmetry), lower, upper):
        assert values[0] == pytest.approx(lo[0])
        assert values[1] == pytest.approx((lo[0] + hi[0]) / 2.0, rel=2e-7)
        assert values[2] == pytest.approx(hi[0])


def test_dispersion_reference_scales_production_refraction_and_group_delay(tmp_path):
    propagation = pytest.importorskip("obdeect.propagation")
    from obdeect.stream import Reader, Writer

    from obdeect_validation.adapters.refraction_curved_reference import (
        C_M_PER_NS,
        NuProfile,
        profile_path,
    )

    site_m, height_m, wavelength_nm, zenith_deg = 2000.0, 10000.0, 250.0, 70.0
    zenith_rad = np.radians(zenith_deg)
    input_cx = np.float32(np.sin(zenith_rad))
    source = tmp_path / "dispersion-input.pst"
    with Writer(str(source), "dispersion validation") as writer:
        writer.write(
            "run",
            {
                "run_number": 1,
                "observation_level_m": site_m,
                "telescope_x_m": [0.0],
                "telescope_y_m": [0.0],
                "telescope_z_m": [0.0],
                "telescope_r_m": [12.0],
                "refraction_applied": 0,
                "wavelength_min_nm": wavelength_nm,
                "wavelength_max_nm": wavelength_nm,
            },
        )
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1, "frame": 0})
        writer.write_photons(
            1,
            0,
            0,
            {
                "x": np.array([0.0], dtype=np.float32),
                "y": np.array([0.0], dtype=np.float32),
                "z": np.array([height_m], dtype=np.float32),
                "cx": np.array([input_cx], dtype=np.float32),
                "cy": np.array([0.0], dtype=np.float32),
                "t": np.array([0.0], dtype=np.float32),
                "w": np.array([1.0], dtype=np.float32),
                "wavelength": np.array([wavelength_nm], dtype=np.float32),
            },
        )
        writer.write("telescope_end", {"telescope_id": 1, "photons": 1.0})
        writer.write("event_end", {})
        writer.write("run_end", {})

    profile = profile_path()
    outputs = []
    for reference_nm in (500.0, wavelength_nm):
        output = tmp_path / f"dispersion-{reference_nm:g}.pst"
        propagation.propagate(
            source,
            output,
            site_altitude_m=site_m,
            seed=1,
            profile=profile,
            extinction=False,
            dispersion_reference_nm=reference_nm,
        )
        photon_chunks = [
            record["columns"]
            for record in Reader(str(output))
            if record["type"] == "photons"
        ]
        assert photon_chunks
        outputs.append(
            {
                key: np.concatenate([chunk[key] for chunk in photon_chunks])
                for key in photon_chunks[0]
            }
        )

    heights = np.linspace(site_m, site_m + height_m, 200001)
    nu = NuProfile(profile)(heights)
    theta_em = np.arcsin(float(input_cx))
    expected = []
    for reference_nm in (500.0, wavelength_nm):
        k_phase = _ciddor_refractivity(wavelength_nm) / _ciddor_refractivity(
            reference_nm
        )
        k_group = _ciddor_group_refractivity(wavelength_nm) / _ciddor_refractivity(
            reference_nm
        )
        sin_theta = np.sin(theta_em) * (1.0 + k_phase * nu[-1]) / (1.0 + k_phase * nu)
        cos_theta = np.sqrt(1.0 - sin_theta**2)
        expected.append(
            (
                np.trapezoid(sin_theta / cos_theta, heights),
                np.trapezoid((1.0 + k_group * nu) / cos_theta, heights) / C_M_PER_NS,
            )
        )
    for photons, (x_ref, t_ref) in zip(outputs, expected):
        assert photons["x"][0] == pytest.approx(x_ref, abs=0.05)
        assert photons["t"][0] == pytest.approx(t_ref, abs=0.1)
    assert outputs[0]["x"][0] != pytest.approx(outputs[1]["x"][0], abs=0.01)
    assert outputs[0]["t"][0] != pytest.approx(outputs[1]["t"][0], abs=0.01)
