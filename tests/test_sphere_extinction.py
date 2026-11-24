"""Unit tests of the sphere-extinction references (propagation.sphere-extinction): analytic limits."""

import numpy as np
import pytest

from obdeect_validation.adapters import common
from obdeect_validation.adapters import sphere_extinction as se

needs_sphere = pytest.mark.skipif(
    not common.obdeect_supports_sphere(),
    reason="needs obdeect with the sphere representation (stream format 1.2)",
)


def test_flat_chord_at_constant_height_is_k_times_length():
    wl, h, length = 410.0, 2100.0, 300.0
    assert se.flat_chord_depth(wl, h, h, length) == pytest.approx(
        abs(se.table().dtau_dh(wl, h)) * length, rel=1e-9
    )


def test_plane_formula_and_flat_chord_agree_for_a_sloped_chord():
    wl = 410.0
    h_e, h_q = 2130.0, 2030.5
    u = 0.1
    length = (h_e - h_q) / u
    assert se.flat_chord_depth(wl, h_e, h_q, length) == pytest.approx(
        se.plane_depth(wl, h_e, h_q, u), abs=1e-9
    )


def test_depth_is_zero_for_zero_length_and_chord_is_linear_in_length_at_constant_height():
    assert se.flat_row_depth(410.0, 2100.0, 0.0, -0.05, 0.0) == 0.0
    # emission and row at the same altitude (z_row = 70 m above the centre, centre at 30 m): a horizontal chord
    short = se.flat_row_depth(410.0, 2100.0, 70.0, -0.05, 10.0)
    long = se.flat_row_depth(410.0, 2100.0, 70.0, -0.05, 20.0)
    assert long == pytest.approx(2.0 * short, rel=1e-9)


def test_spherical_chord_reduces_to_flat_for_a_short_horizontal_chord():
    # a 200 m horizontal chord at the telescope altitude: the sphere dip is L^2 / 2R ~ 1.6 mm
    wl = 410.0
    q = np.array([0.0, 0.0, 0.0])
    d = np.array([1.0, 0.0, 0.0])
    _, spherical = se.spherical_chord_depth(wl, q, d, 200.0)
    flat = se.flat_chord_depth(wl, 2030.0 + 0.0, 2030.0, 200.0)
    assert spherical == pytest.approx(flat, rel=1e-4)


def test_photon_groups_cover_every_class():
    sets = se.photon_sets()
    groups = set(sets["groups"].tolist())
    assert groups == {
        "a_steep",
        "a_upward",
        "b_near_horizontal",
        "b_horizontal",
        "c_continuity",
        "d_equator",
        "e_inside",
    }
    cz = sets["columns"]["cz"]
    assert np.all(cz[sets["groups"] == "b_horizontal"] == 0.0)


@needs_sphere
def test_plane_and_sphere_rows_of_one_photon_agree(tmp_path):
    from obdeect.propagation import propagate

    sets = se.photon_sets()
    cols = sets["columns"]
    mask = sets["groups"] == "d_equator"
    sub = {k: v[mask] for k, v in cols.items()}
    inp = se.write_input(tmp_path / "in.pst", sub)
    propagate(
        inp,
        tmp_path / "sph.pst",
        transmission=str(se.uc.TRANSMISSION),
        site_altitude_m=se.OBS_M,
        seed=1,
        ground_representation="sphere",
        mode="weighted",
        store_p=True,
        pathlen_m=0.0,
    )
    propagate(
        inp,
        tmp_path / "pl.pst",
        transmission=str(se.uc.TRANSMISSION),
        site_altitude_m=se.OBS_M,
        seed=1,
        mode="weighted",
        store_p=True,
        pathlen_m=0.0,
    )
    sphere = se.read_rows(tmp_path / "sph.pst")
    plane = se.read_rows(tmp_path / "pl.pst")
    assert len(sphere["p"]) == len(plane["p"]) == int(mask.sum())
    assert np.max(np.abs(sphere["p"] - plane["p"]) / plane["p"]) < 1e-6
