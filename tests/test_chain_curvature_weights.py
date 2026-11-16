"""Detection weights of the curved chain case on synthetic tables (no sim_telarray reference).

A CORSIKA Cherenkov band of 100 to 2000 nm is wider than the mirror table (200 to 1099.5 nm).
The weight is quantum efficiency times reflectivity, so it is exactly zero where the quantum
efficiency is zero; the table only has to cover the non-zero part of the band.
"""

import numpy as np
import pytest

from obdeect_validation.adapters import chain_curvature


def _write(path, rows):
    path.write_text("".join(f"{nm} {value}\n" for nm, value in rows), encoding="utf-8")


@pytest.fixture
def weights(tmp_path):
    (tmp_path / "cfg" / "CTA").mkdir(parents=True)
    _write(
        tmp_path / chain_curvature.QE_TABLE,
        [(200.0, 0.0), (220.0, 0.2), (300.0, 0.3), (800.0, 0.0)],
    )
    _write(
        tmp_path / chain_curvature.MIRROR_TABLE,
        [(200.0, 0.8), (1100.0, 0.5)],
    )
    _write(tmp_path / chain_curvature.WINDOW_TABLE, [(100.0, 1.0), (2000.0, 1.0)])
    return chain_curvature.DetectionWeights(reference=tmp_path)


def test_band_wider_than_mirror_table_is_accepted_where_quantum_is_zero(weights):
    nm = np.array([100.0, 150.0, 250.0, 2000.0])
    values = weights.detection(nm)
    assert values[0] == 0.0 and values[1] == 0.0 and values[3] == 0.0
    reflectivity_250 = np.interp(250.0, [200.0, 1100.0], [0.8, 0.5])
    assert values[2] == pytest.approx(float(weights.quantum(250.0)) * reflectivity_250)


def test_camera_window_multiplies_the_weight(tmp_path):
    (tmp_path / "cfg" / "CTA").mkdir(parents=True)
    _write(tmp_path / chain_curvature.QE_TABLE, [(200.0, 0.5), (800.0, 0.5)])
    _write(tmp_path / chain_curvature.MIRROR_TABLE, [(200.0, 1.0), (1100.0, 1.0)])
    _write(
        tmp_path / chain_curvature.WINDOW_TABLE,
        [(259.0, 0.0), (300.0, 0.5), (1000.0, 0.5)],
    )
    weights = chain_curvature.DetectionWeights(reference=tmp_path)
    # the window is evaluated at lambda + 0.5 (sim_config.c): a linear ramp from 259 to 300 nm
    expected = 0.5 * np.interp(260.5, [259.0, 300.0], [0.0, 0.5])
    assert weights.detection(np.array([260.0]))[0] == pytest.approx(expected)


def test_detection_equals_quantum_times_reflectivity(weights):
    nm = np.array([220.0, 300.0, 500.0])
    expected = weights.quantum(nm) * np.interp(nm, [200.0, 1100.0], [0.8, 0.5])
    assert weights.detection(nm) == pytest.approx(expected)


def test_non_zero_quantum_outside_mirror_table_is_rejected(tmp_path):
    (tmp_path / "cfg" / "CTA").mkdir(parents=True)
    _write(
        tmp_path / chain_curvature.QE_TABLE,
        [(200.0, 0.1), (800.0, 0.1)],
    )
    _write(tmp_path / chain_curvature.MIRROR_TABLE, [(300.0, 0.8), (1100.0, 0.5)])
    _write(tmp_path / chain_curvature.WINDOW_TABLE, [(100.0, 1.0), (2000.0, 1.0)])
    weights = chain_curvature.DetectionWeights(reference=tmp_path)
    with pytest.raises(ValueError, match="does not cover"):
        weights.detection(np.array([250.0]))
