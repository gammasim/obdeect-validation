"""Unit tests of the ground-frame sources case (sources.ground-muon)."""

import math

import numpy as np
import pytest

from obdeect_validation.adapters import common
from obdeect_validation.adapters import sphere_sources as ss

needs_sphere = pytest.mark.skipif(
    not common.obdeect_supports_sphere(),
    reason="needs obdeect with the sphere representation (stream format 1.2)",
)


def test_chord_fraction_is_the_track_inside_the_sphere():
    impact = math.hypot(ss.LINE_START[1], ss.LINE_START[2])
    chord = 2.0 * math.sqrt(ss.RADIUS_M**2 - impact**2)
    assert ss.chord_fraction() == pytest.approx(chord / ss.LINE_LENGTH_M, rel=1e-9)


def test_pull_is_zero_at_the_expectation_and_signed():
    assert ss.pull(50, 100, 0.5) == pytest.approx(0.0)
    assert ss.pull(60, 100, 0.5) == pytest.approx(2.0)


def _data(rows):
    keys = ("x", "y", "z", "cx", "cy", "cz", "t", "zem", "w", "path", "wavelength")
    columns = {k: np.array([r[k] for r in rows], dtype=np.float32) for k in keys}
    return {
        "columns": columns,
        "run_end": {"misses_sphere": 5},
        "run": {"ground_representation": "sphere"},
    }


def test_geometry_checks_accept_consistent_inside_and_entry_rows():
    u = ss.LINE_DIRECTION
    along_in = 50.0
    inside_pos = ss.LINE_START + along_in * u  # (-10, 2, 3) m: |r| = 10.6 < R
    along_out = 20.0
    emission = ss.LINE_START + along_out * u  # (-40, 2, 3) m, outside
    entry_x = -math.sqrt(ss.RADIUS_M**2 - 2.0**2 - 3.0**2)
    entry = np.array([entry_x, 2.0, 3.0])
    path_out = float(np.linalg.norm(entry - emission))
    rows = [
        {
            "x": inside_pos[0],
            "y": inside_pos[1],
            "z": inside_pos[2],
            "cx": 0.0,
            "cy": 0.0,
            "cz": 1.0,
            "t": along_in / ss.C_M_NS,
            "zem": ss.OBSERVATION_M + 3.0,
            "w": 1.0,
            "path": 0.0,
            "wavelength": 400.0,
        },
        {
            "x": entry[0],
            "y": entry[1],
            "z": entry[2],
            "cx": 1.0,
            "cy": 0.0,
            "cz": 0.0,
            "t": (along_out + path_out) / ss.C_M_NS,
            "zem": ss.OBSERVATION_M + 3.0,
            "w": 1.0,
            "path": path_out,
            "wavelength": 400.0,
        },
    ]
    metrics, _ = ss.check(_data(rows), {"misses_sphere": 5})
    assert metrics["line_on_track_max_m"] < 1e-5
    assert metrics["inside_path_nonzero"] == 0 and metrics["outside_path_zero"] == 0
    assert metrics["entry_radius_max_error_m"] < 1e-5
    assert metrics["time_budget_ratio_max"] < 1.0
    assert metrics["zem_budget_ratio_max"] < 1.0


def test_a_row_with_a_wrong_time_breaks_the_time_gate():
    u = ss.LINE_DIRECTION
    pos = ss.LINE_START + 50.0 * u
    rows = [
        {
            "x": pos[0],
            "y": pos[1],
            "z": pos[2],
            "cx": 0.0,
            "cy": 0.0,
            "cz": 1.0,
            "t": 50.0 / ss.C_M_NS + 1.0,
            "zem": ss.OBSERVATION_M + 3.0,
            "w": 1.0,
            "path": 0.0,
            "wavelength": 400.0,
        }
    ]
    metrics, _ = ss.check(_data(rows), {"misses_sphere": 0})
    assert metrics["time_budget_ratio_max"] > 1.0


@needs_sphere
def test_small_run_balances_rows_and_misses(tmp_path):
    from obdeect.sources import run

    description = ss.description()
    for emission in description["events"][0]["emissions"]:
        emission["photons"] = 2000
    summary = run.write(description, tmp_path / "small.pst")
    data = ss.read_rows(tmp_path / "small.pst")
    rows = len(data["columns"]["x"])
    assert (
        rows + data["run_end"]["misses_sphere"] == 4000
    )  # every emitted photon is a row or a miss
    metrics, _ = ss.check(data, summary)
    assert metrics["misses_sphere_summary_mismatch"] == 0
    assert metrics["ground_representation_sphere"] == 1
    assert metrics["inside_path_nonzero"] == 0
