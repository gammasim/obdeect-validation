"""Unit tests of the optics sphere-vs-plane case (optics.sphere-vs-plane): beams, comparison, enclosure."""

import numpy as np
import pytest

from obdeect_validation.adapters import common
from obdeect_validation.adapters import sphere_optics as so

needs_sphere = pytest.mark.skipif(
    not (
        common.obdeect_supports_sphere() and common.obdeect_has("obdeect.ray_tracing")
    ),
    reason="needs obdeect with the sphere representation (stream format 1.2) and ray_tracing",
)


def test_parallel_beam_is_downward_along_the_star_axis_and_on_the_disc():
    cols = so.parallel_beam(500, 2.0, 1)
    d = np.stack([cols["cx"], cols["cy"], cols["cz"]], axis=1).astype(np.float64)
    assert np.allclose(np.linalg.norm(d, axis=1), 1.0, atol=1e-6)
    assert np.allclose(d, -so.star_axis(), atol=1e-6)
    assert float(cols["cz"][0]) < 0.0
    pos = np.stack([cols["x"], cols["y"], cols["z"]], axis=1).astype(np.float64)
    lateral = pos - so.START_M * so.star_axis()
    assert np.max(np.linalg.norm(lateral, axis=1)) <= 2.0 + 1e-4
    assert np.all(np.abs(lateral @ so.star_axis()) < 1e-3)


def test_horizontal_rows_have_cz_zero_and_unique_identities():
    cols = so.horizontal_rows(300, 2)
    assert np.all(cols["cz"] == 0.0)
    assert len(np.unique(cols["wavelength"])) == 300


def test_budget_is_four_spacings_of_the_value():
    assert float(so.budget(1.0)) == pytest.approx(4.0 * np.spacing(np.float32(1.0)))


def _camera(ids, x, t):
    n = len(ids)
    return {
        "wavelength": so.WAVELENGTH_BASE_NM
        + so.WAVELENGTH_STEP_NM * np.asarray(ids, dtype=np.float64),
        "x": np.asarray(x, dtype=np.float32),
        "y": np.zeros(n, dtype=np.float32),
        "t": np.asarray(t, dtype=np.float32),
        "zem": np.zeros(n, dtype=np.float32),
        "cx": np.zeros(n, dtype=np.float32),
        "cy": np.zeros(n, dtype=np.float32),
        "cz": -np.ones(n, dtype=np.float32),
    }


def test_compare_camera_identical_sets_have_no_difference():
    plane = _camera([1, 4, 7], [0.5, 1.5, 2.5], [100.0, 200.0, 300.0])
    result = so.compare_camera(
        plane, _camera([7, 1, 4], [2.5, 0.5, 1.5], [300.0, 100.0, 200.0])
    )
    assert result["set_mismatch"] == 0
    assert result["max_position_difference_m"] == 0.0
    assert result["max_time_difference_ns"] == 0.0


def test_compare_camera_counts_a_photon_detected_on_one_side_only():
    plane = _camera([1, 4, 7], [0.5, 1.5, 2.5], [100.0, 200.0, 300.0])
    sphere = _camera([1, 4], [0.5, 1.5], [100.0, 200.0])
    assert so.compare_camera(plane, sphere)["set_mismatch"] == 1


def test_compare_camera_reports_a_time_offset_against_the_budget():
    plane = _camera([1], [0.5], [6768.0])
    sphere = _camera([1], [0.5], [6768.0 + 2 * float(np.spacing(np.float32(6768.0)))])
    result = so.compare_camera(plane, sphere)
    assert result["time_budget_ratio"] == pytest.approx(0.5, rel=1e-6)


@needs_sphere
def test_small_sphere_is_refused_and_the_fixture_sphere_is_accepted(tmp_path):
    from obdeect.propagation import propagate
    from obdeect.ray_tracing import OpticalModel
    from obdeect.ray_tracing.stream import trace_stream

    model = OpticalModel.load(str(so.FIXTURE))
    small = so._write(tmp_path / "in.pst", so.parallel_beam(200, 2.1, 3), 0.5)
    propagate(
        small,
        tmp_path / "small.pst",
        extinction=False,
        site_altitude_m=0.0,
        seed=1,
        ground_representation="sphere",
    )
    with pytest.raises(Exception, match="must enclose the optical model"):
        trace_stream(
            model,
            str(tmp_path / "small.pst"),
            str(tmp_path / "small-cam.pst"),
            seed=5,
            pointing=(so.ZENITH_DEG, so.AZIMUTH_DEG),
        )
    ok = so._write(
        tmp_path / "ok.pst", so.parallel_beam(200, 2.1, 3), so.FIXTURE_RADIUS_M
    )
    propagate(
        ok,
        tmp_path / "ok-sphere.pst",
        extinction=False,
        site_altitude_m=0.0,
        seed=1,
        ground_representation="sphere",
    )
    summary = trace_stream(
        model,
        str(tmp_path / "ok-sphere.pst"),
        str(tmp_path / "ok-cam.pst"),
        seed=5,
        pointing=(so.ZENITH_DEG, so.AZIMUTH_DEG),
    )
    assert summary["photons_in"] > 0
