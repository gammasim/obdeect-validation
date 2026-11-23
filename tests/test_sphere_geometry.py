"""Unit tests of the sphere-geometry case (propagation.sphere-geometry): references and helpers.

The ray-integration reference is checked against the closed-form straight line (no refraction: the
profile is zero above its top). The obdeect comparison needs the sphere representation and skips otherwise.
"""

import math

import numpy as np
import pytest

from obdeect_validation.adapters import common
from obdeect_validation.adapters import sphere_geometry as sg
from obdeect_validation.adapters import upward_common as uc

needs_sphere = pytest.mark.skipif(
    not common.obdeect_supports_sphere(),
    reason="needs obdeect with the sphere representation (stream format 1.2)",
)


def closed_form_entry(origin, direction, centre, radius):
    rel = np.asarray(origin, float) - np.asarray(centre, float)
    d = sg.unit(direction)
    b = rel @ d
    disc = b * b - (rel @ rel - radius * radius)
    s = -b - math.sqrt(disc)
    return s, rel + s * d


def test_stored_cz_convention_gives_the_documented_sign_and_magnitude():
    assert float(sg.stored_cz(0.6, 0.0, True)) == pytest.approx(-0.8, abs=1e-7)
    assert float(sg.stored_cz(0.6, 0.0, False)) == pytest.approx(0.8, abs=1e-7)
    assert sg.stored_cz(0.0, 0.0, True) == np.float32(-1.0)


def test_ray_integration_reproduces_the_straight_line_above_the_profile_top():
    ref = sg.RayReference(uc.profile(), uc.SITE_M)
    centre = np.array([0.0, 0.0, 20.0])
    origin = np.array([30.0, -20.0, 130000.0])  # above the top (120 km): no refraction
    direction = (
        np.array([5.0, -3.0, 20.0]) - origin
    )  # aimed at the sphere, 5 m off its centre
    entered, s, rel, _ = ref.trace(
        origin, direction, 0.0, centre, sg.RADIUS_M, 300000.0
    )
    s_exact, pos_exact = closed_form_entry(origin, direction, centre, sg.RADIUS_M)
    assert entered
    # relative 1e-8 over a 130 km path: the coarse phase locates the shell crossing as an event
    assert s == pytest.approx(s_exact, rel=1e-8)
    assert rel == pytest.approx(pos_exact, abs=2e-3)


def test_no_refraction_reference_classifies_drops_and_inside_rows():
    centre = np.asarray(sg.TELESCOPE_M)
    cols = {
        "x": np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32),
        "y": np.array([0.0, 0.0, 0.0, 40.0], dtype=np.float32),
        "z": np.array([-5.0, 20.0, 300.0, 300.0], dtype=np.float32),
        "cx": np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32),
        "cy": np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32),
        "cz": np.array([1.0, -1.0, -1.0, -1.0], dtype=np.float32),
        "t": np.zeros(4, dtype=np.float32),
    }
    ref = sg.reference_no_refraction(cols)
    assert list(ref["reason"]) == ["below_ground", "ok", "ok", "misses_sphere"]
    assert ref["values"]["path"][1] == 0.0
    assert ref["values"]["x"][1] == pytest.approx(centre[0])
    assert ref["values"]["path"][2] == pytest.approx(300.0 - 20.0 - sg.RADIUS_M)


def test_match_by_index_places_survivors_and_rejects_unknown_identities():
    rows = {
        "wavelength": np.array(
            [sg.ID_BASE_NM, sg.ID_BASE_NM + 2 * sg.ID_STEP_NM], dtype=np.float32
        ),
        "x": np.array([1.0, 2.0], dtype=np.float32),
        "y": np.zeros(2, dtype=np.float32),
        "z": np.zeros(2, dtype=np.float32),
        "cx": np.zeros(2, dtype=np.float32),
        "cy": np.zeros(2, dtype=np.float32),
        "cz": -np.ones(2, dtype=np.float32),
        "t": np.zeros(2, dtype=np.float32),
        "zem": np.zeros(2, dtype=np.float32),
        "w": np.ones(2, dtype=np.float32),
        "path": np.zeros(2, dtype=np.float32),
    }
    out = sg.match_by_index(rows, 3)
    assert list(out["present"]) == [True, False, True]
    assert out["x"][2] == pytest.approx(2.0)
    rows["wavelength"] = np.array(
        [sg.ID_BASE_NM + 5 * sg.ID_STEP_NM, sg.ID_BASE_NM], dtype=np.float32
    )
    with pytest.raises(ValueError):
        sg.match_by_index(rows, 3)


def test_float32_budget_is_one_spacing_plus_the_direction_budget():
    one = float(sg.float32_budget(1.0, 0.0))
    assert one == pytest.approx(np.spacing(np.float32(1.0)))
    assert float(sg.float32_budget(1.0, 3000.0)) == pytest.approx(
        one + 2 * sg.EPS32 * 3000.0
    )
    assert float(sg.float32_budget(1.0, 3000.0, time=True)) == pytest.approx(
        one + 2 * sg.EPS32 * 3000.0 / sg.C_M_NS
    )


def test_refraction_reference_drops_points_below_ground():
    cols = {
        "x": np.array([0.0], dtype=np.float32),
        "y": np.array([0.0], dtype=np.float32),
        "z": np.array([-3.0], dtype=np.float32),
        "cx": np.array([0.0], dtype=np.float32),
        "cy": np.array([0.0], dtype=np.float32),
        "cz": np.array([1.0], dtype=np.float32),
        "t": np.zeros(1, dtype=np.float32),
    }
    out = sg.reference_refraction(cols, uc.profile(), sg.TELESCOPE_M, s_max=6000.0)
    assert out["reason"][0] == "below_ground"


@needs_sphere
def test_obdeect_entry_of_one_straight_photon_is_on_the_sphere(tmp_path):
    from obdeect_validation.adapters import sphere_geometry as adapter

    cols = {
        "x": np.array([0.0], dtype=np.float32),
        "y": np.array([0.0], dtype=np.float32),
        "z": np.array([300.0], dtype=np.float32),
        "cx": np.array([0.0], dtype=np.float32),
        "cy": np.array([0.0], dtype=np.float32),
        "cz": np.array([-1.0], dtype=np.float32),
        "t": np.array([5.0], dtype=np.float32),
        "w": np.ones(1, dtype=np.float32),
        "wavelength": np.array([adapter.ID_BASE_NM], dtype=np.float32),
    }
    _summary, raw, _ = adapter.run_obdeect(tmp_path, "single", cols)
    rows = adapter.match_by_index(raw, 1)
    assert rows["present"][0]
    assert np.linalg.norm([rows["x"][0], rows["y"][0], rows["z"][0]]) == pytest.approx(
        sg.RADIUS_M, abs=1e-5
    )
    assert rows["path"][0] == pytest.approx(
        300.0 - sg.TELESCOPE_M[2] - sg.RADIUS_M, abs=1e-4
    )
    assert rows["t"][0] == pytest.approx(5.0 + rows["path"][0] / sg.C_M_NS, abs=1e-3)
