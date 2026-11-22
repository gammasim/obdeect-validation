"""Independent checks for wavelength interpolation and dispersive air delay."""

from pathlib import Path

import numpy as np
import pytest

WORKSPACE = Path(__file__).resolve().parents[2]
SCATTERING = WORKSPACE / "obdeect-validation/fixtures/scattering_exponential.dat"


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


def test_reference_wavelength_group_factor_is_independently_defined():
    reference_nm = 500.0
    k_group = _ciddor_group_refractivity(reference_nm) / _ciddor_refractivity(
        reference_nm
    )
    assert k_group == pytest.approx(1.04724, rel=2e-4)
