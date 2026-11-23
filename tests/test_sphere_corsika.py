"""Unit tests of the CORSIKA 1217 sphere-row case (stream.corsika-photons3d).

The EventIO walk, the 1217 decoding and the sphere rules are checked on analytic bunches first; the
obdeect comparisons need the sphere representation (stream format 1.2) and skip otherwise.
"""

import numpy as np
import pytest

from obdeect_validation.adapters import common
from obdeect_validation.adapters import sphere_corsika as sc

needs_sphere = pytest.mark.skipif(
    not common.obdeect_supports_sphere(),
    reason="needs obdeect with the sphere representation (stream format 1.2)",
)


def test_eventio_walk_finds_the_run_telescopes_and_bunches():
    run = sc.read_run()
    assert run["observation_level_m"] == pytest.approx(2200.0)
    assert [round(t["r_cm"]) for t in run["telescopes"]] == [1200, 800]
    assert [round(t["z_cm"]) for t in run["telescopes"]] == [1800, 2000]
    assert [len(b["values"]) for b in run["bunches"]] == [8, 3]
    assert run["blocks"] == [1200, 1201, 1202, 1203, 1204, 1209, 1210]


def test_pyeventio_agrees_with_the_walk_on_top_level_blocks():
    cross = sc.pyeventio_cross_check()
    assert cross["types"] == sc.read_run()["blocks"]
    assert cross["n_telescopes"] == 2


def _bunch(values, telescope=0):
    return {
        "array": 0,
        "telescope": telescope,
        "photons": 1.0,
        "values": np.array([values], dtype=np.float32),
    }


def _run(bunches):
    telescopes = [{"x_cm": 0.0, "y_cm": 0.0, "z_cm": 1800.0, "r_cm": 1200.0}]
    return {"observation_level_m": 2200.0, "telescopes": telescopes, "bunches": bunches}


def test_downward_bunch_emitted_above_the_sphere_enters_at_its_top():
    # A = (0, 0, 1200) cm, d = (0, 0, -1), dist = 3000 cm: E = (0, 0, 42) m, entry (0, 0, 12) m, path 30 m
    row = sc.reference_rows(_run([_bunch([0, 0, 1200, 0, 0, -1, 100, 3000, 2, 0])]))[0]
    assert not row["miss"] and not row["inside"]
    assert (row["x"], row["y"], row["z"]) == pytest.approx((0.0, 0.0, 12.0), abs=1e-6)
    assert row["path"] == pytest.approx(30.0)
    assert row["t"] == pytest.approx(100.0, abs=1e-9)
    assert row["zem"] == pytest.approx(2200.0 + 18.0 + 42.0)


def test_bunch_emitted_inside_is_a_row_at_its_emission_point():
    # A = (300, 0, -400) cm, d = (0, 0, -1), dist = 800 cm: E = (3, 0, 4) m, inside, path 0, t = ctime - dist/c
    row = sc.reference_rows(_run([_bunch([300, 0, -400, 0, 0, -1, 50, 800, 5, 350])]))[
        0
    ]
    assert row["inside"] and row["path"] == 0.0
    assert (row["x"], row["y"], row["z"]) == pytest.approx((3.0, 0.0, 4.0), abs=1e-6)
    assert row["t"] == pytest.approx(50.0 - 8.0 / sc.C_M_NS)


def test_bunch_passing_the_sphere_at_a_distance_is_a_miss():
    # horizontal line at y = 15 m passes the 12 m sphere: misses_sphere
    row = sc.reference_rows(_run([_bunch([0, 1500, 0, 1, 0, 0, 500, 3000, 9, 0])]))[0]
    assert row["miss"] is True


def test_bunch_pointing_away_from_the_sphere_is_a_miss():
    # A = (0, 0, 30) m above the centre, d = (0, 0, 1), dist = 10 m: E = 20 m above, moving up (away)
    row = sc.reference_rows(_run([_bunch([0, 0, 3000, 0, 0, 1, 600, 1000, 7, 0])]))[0]
    assert row["miss"] is True


def test_ulp_error_counts_float32_spacings():
    assert float(sc.ulp_error(np.float32(1.0), 1.0)) == 0.0
    assert float(
        sc.ulp_error(np.nextafter(np.float32(1.0), np.float32(2.0)), 1.0)
    ) == pytest.approx(1.0)


def test_plane_reference_moves_a_downward_row_to_z_zero():
    rows = [
        {
            "telescope": 1,
            "x": 0.0,
            "y": 0.0,
            "z": 12.0,
            "cx": 0.0,
            "cy": 0.0,
            "cz": -1.0,
            "t": 100.0,
            "zem": 2260.0,
            "w": 2.0,
            "wavelength": 0.0,
        }
    ]
    out = sc.plane_reference(rows)[0]
    assert out["x"] == pytest.approx(0.0) and out["y"] == pytest.approx(0.0)
    assert out["t"] == pytest.approx(100.0 + 12.0 / sc.C_M_NS)
    assert out["zem"] == 2260.0


@needs_sphere
def test_from_corsika_rows_match_the_reference(tmp_path):
    run = sc.read_run()
    reference = sc.reference_rows(run)
    pst = tmp_path / "photons3d.pst"
    common.command(
        [common.tool("obdeect"), "stream", "from-corsika", sc.FIXTURE],
        tmp_path / "log",
        stdout_file=pst,
    )
    stream = sc.read_stream(pst)
    compared = sc.compare_rows(reference, stream["rows"])
    assert compared["row_count_obdeect"] == compared["row_count_reference"] == 9
    assert max(compared["max_ulp"].values()) <= 1.0
    assert stream["run"]["ground_representation"] == "sphere"
    assert (
        stream["run_end"]["misses_sphere"]
        == sum(1 for r in reference if r["miss"])
        == 2
    )
