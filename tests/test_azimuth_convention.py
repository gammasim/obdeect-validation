"""Unit tests of the azimuth-convention validation (analytic directions, counting rules).

Only numpy is needed; the obdeect runs are in the adapter and the cases.
"""

import math

import numpy as np
import pytest

from obdeect_validation.adapters import azimuth_convention as A


def test_corsika_direction_axes():
    assert np.allclose(A.corsika_direction(0.0, 37.0), [0.0, 0.0, -1.0])
    assert np.allclose(A.corsika_direction(90.0, 0.0), [1.0, 0.0, 0.0])  # north
    assert np.allclose(A.corsika_direction(90.0, 90.0), [0.0, 1.0, 0.0])  # west
    assert np.allclose(A.corsika_direction(30.0, 180.0), [-0.5, 0.0, -math.sqrt(3) / 2])


def test_astronomical_azimuth_is_clockwise_from_north_to_east():
    assert np.allclose(A.astro_unit(90.0, 0.0), [1.0, 0.0, 0.0])  # north
    assert np.allclose(
        A.astro_unit(90.0, 90.0), [0.0, -1.0, 0.0]
    )  # east is -y (y west)
    assert np.allclose(A.astro_unit(0.0, 123.0), [0.0, 0.0, 1.0])  # up


def test_expected_roles():
    source = next(p for p in A.PRESETS if p.name == "muon")
    beam = next(p for p in A.PRESETS if p.name == "laser-beam")
    night = next(p for p in A.PRESETS if p.name == "night-sky")
    assert beam.zenith_input == pytest.approx(150.0)  # 30 deg from down
    # sky sources: the light comes from the east (az 90) and travels west and down (+y, -z)
    assert np.allclose(source.expected(90.0), [0.0, 0.5, -math.sqrt(3) / 2])
    assert np.allclose(night.expected(90.0), source.expected(90.0))
    # beams: the pointing is the travel; az 90 points east (-y) and down
    assert np.allclose(beam.expected(90.0), [0.0, -0.5, -math.sqrt(3) / 2])
    assert np.allclose(beam.expected(0.0), [0.5, 0.0, -math.sqrt(3) / 2])


def test_named_direction_is_source_for_sky_and_pointing_for_beams():
    source = next(p for p in A.PRESETS if p.name == "muon")
    beam = next(p for p in A.PRESETS if p.name == "laser-beam")
    travel = np.array([0.0, 0.5, -math.sqrt(3) / 2])
    assert np.allclose(source.named(travel), -travel)
    assert np.allclose(beam.named(travel), travel)


def test_mean_direction_of_a_single_direction_has_no_spread():
    u, se, rows = A.mean_direction(
        [0.5, 0.5], [0.0, 0.0], [-math.sqrt(0.75)] * 2, [1.0, 3.0]
    )
    assert np.allclose(u, [0.5, 0.0, -math.sqrt(0.75)])
    assert se == pytest.approx(0.0, abs=1e-12)
    assert rows == 2


def test_mean_direction_recovers_cone_axis():
    rng = np.random.default_rng(1)
    theta = 0.024
    phi = rng.uniform(0, 2 * math.pi, 200000)
    axis = A.corsika_direction(30.0, 90.0)
    # cone of half angle theta around the axis
    e1 = np.array([math.cos(math.radians(30)) * 0 + 1.0, 0.0, 0.0])
    e1 = e1 - np.dot(e1, axis) * axis
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    v = math.cos(theta) * axis + math.sin(theta) * (
        np.outer(np.cos(phi), e1) + np.outer(np.sin(phi), e2)
    )
    u, se, _ = A.mean_direction(v[:, 0], v[:, 1], v[:, 2], np.ones(len(phi)))
    assert A.angle_rad(u, axis) < 5 * se + 1e-9
    assert se < 1e-3


def test_violation_and_disagreement_rules():
    expected = A.corsika_direction(30.0, 90.0)
    wrong = A.corsika_direction(30.0, 270.0)
    assert A.angle_rad(expected, wrong) == pytest.approx(
        math.pi / 3, abs=1e-9
    )  # cos = 0.5
    base = {"se_rad": 0.0, "named": expected.tolist()}
    other_same = {"se_rad": 0.0, "named": expected.tolist()}
    other_wrong = {"se_rad": 0.0, "named": wrong.tolist()}
    assert not A.pair_disagrees(base, other_same)
    assert A.pair_disagrees(base, other_wrong)
    # an azimuth offset of 1e-3 rad at zenith 30 deg is an angle of 0.5e-3 rad: below the requirement
    near = A.corsika_direction(30.0, 90.0 + math.degrees(1e-3))
    assert A.angle_rad(expected, near) < A.ANGLE_TOL_RAD
    assert not A.pair_disagrees(base, {"se_rad": 0.0, "named": near.tolist()})


def test_tolerance_takes_the_larger_of_requirement_and_k_se():
    assert A.tolerance(0.0) == A.ANGLE_TOL_RAD
    assert A.tolerance(1e-3) == pytest.approx(A.K_SE * 1e-3)


def test_corsika_angles_from_log(tmp_path):
    log = tmp_path / "x.corsika.log"
    log.write_text("some header\nTHETAP 20 20\nPHIP 180 180\n", encoding="utf-8")
    assert A.parse_corsika_angles(log) == (20.0, 180.0)
    bad = tmp_path / "bad.corsika.log"
    bad.write_text("nothing here\n", encoding="utf-8")
    with pytest.raises(ValueError):
        A.parse_corsika_angles(bad)


def test_azimuth_wrap():
    assert A._wrap_deg(359.0) == pytest.approx(-1.0)
    assert A._wrap_deg(-179.0 + 360.0) == pytest.approx(-179.0)
    assert A._wrap_deg(0.0) == 0.0
