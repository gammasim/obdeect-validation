"""Joint array chain helpers: union efficiency bound, stereo trigger file."""

import numpy as np
import pytest

from obdeect_validation.adapters import performance_array_joint as joint
from obdeect_validation.adapters import performance_chain as pc


def _elements(tmp_path):
    low = tmp_path / "low.dat"
    high = tmp_path / "high.dat"
    low.write_text("300 0.10\n600 0.20\n")
    high.write_text("300 0.30\n600 0.05\n")
    elements = {}
    for _, element in joint.ARRAY:
        elements[element] = {"qe": low, "headroom": 0.7}
    elements["MSTN-03"] = {"qe": high, "headroom": 0.8}
    elements["MSTN-15"] = {"qe": high, "headroom": 0.8}
    return elements


def test_union_bound_covers_every_qe(tmp_path):
    elements = _elements(tmp_path)
    files = joint.union_files(tmp_path / "w", elements)
    grid = pc.WAVELENGTH_GRID
    bound = np.maximum(
        pc.sim_telarray_qe(elements["LSTN-01"]["qe"]) / 0.7,
        pc.sim_telarray_qe(elements["MSTN-03"]["qe"]) / 0.8,
    )
    for name, headroom in (("low.dat", 0.7), ("high.dat", 0.8)):
        scaled = np.loadtxt(files["table_dir"] / name)[:, 1]
        qe = pc.sim_telarray_qe(tmp_path / name)
        assert np.allclose(scaled * bound, qe)
        assert scaled.max() <= headroom + 1e-9
    assert files["bound_max"] == pytest.approx(bound.max())
    assert len(grid) == len(bound)


def test_stereo_trigger_changes_only_the_lst_line(tmp_path):
    (tmp_path / "array_triggers.dat").write_text(
        "# x\nTrigger 1 of 1, 2, 3, 4 width 120.0\nTrigger 1 of 5, 6 width 200.0 minsep 40.0\n"
    )
    out = joint.stereo_triggers(tmp_path, tmp_path / "stereo.dat")
    lines = out.read_text().splitlines()
    assert lines[1].startswith("Trigger 2 of 1, 2, 3, 4")
    assert lines[2].startswith("Trigger 1 of 5, 6")


def test_stereo_trigger_needs_the_lst_line(tmp_path):
    (tmp_path / "array_triggers.dat").write_text("Trigger 1 of 5, 6 width 200.0\n")
    with pytest.raises(RuntimeError):
        joint.stereo_triggers(tmp_path, tmp_path / "stereo.dat")
