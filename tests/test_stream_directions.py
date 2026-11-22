"""Unit tests of the direction round trip.

Needs obdeect.stream; the CLI test also needs the obdeect command.
"""

import shutil

import numpy as np
import pytest

pytest.importorskip("obdeect.stream")
stream_directions = pytest.importorskip("obdeect_validation.adapters.stream_directions")


def test_grid_has_exact_vertical_and_horizontal_directions():
    grid = stream_directions.directions()
    cz = grid["cz"]
    assert cz.dtype == np.float32
    assert cz[stream_directions.UP] == 1.0
    assert cz[stream_directions.DOWN] == -1.0
    horizontal = slice(
        stream_directions.FIRST_HORIZONTAL,
        stream_directions.FIRST_HORIZONTAL + stream_directions.HORIZONTAL,
    )
    assert np.all(cz[horizontal] == 0.0)
    norm = grid["cx"].astype(np.float64) ** 2 + grid["cy"] ** 2 + grid["cz"] ** 2
    np.testing.assert_allclose(norm, 1.0, atol=1e-6)
    counts = stream_directions.sign_counts(cz)
    assert counts["horizontal"] == stream_directions.HORIZONTAL
    assert counts["upward"] > 0 and counts["downward"] > 0
    assert sum(counts.values()) == stream_directions.GRID_SIZE


def test_round_trip_is_bit_exact(tmp_path):
    grid = stream_directions.directions()
    stream = tmp_path / "directions.pst"
    written = stream_directions.write_round_trip(stream, grid)
    header, chunks = stream_directions.read_back(stream)
    assert header["format_minor"] == 1
    result = stream_directions.compare(written, chunks)
    assert result["photon_count"] == 4 * stream_directions.GRID_SIZE
    assert result["bit_mismatch_count"] == 0
    assert result["max_abs_difference"] == 0.0


def test_compare_counts_a_single_changed_bit(tmp_path):
    grid = stream_directions.directions()
    stream = tmp_path / "directions.pst"
    written = stream_directions.write_round_trip(stream, grid)
    _, chunks = stream_directions.read_back(stream)
    telescope = 3  # ids 1..4: sky, ground, telescope, camera
    column = chunks[telescope][0]["columns"]["cz"].copy()
    column[5] = np.nextafter(column[5], np.float32(2.0))
    chunks[telescope][0]["columns"]["cz"] = column
    result = stream_directions.compare(written, chunks)
    assert result["bit_mismatch_count"] == 1
    assert result["per_frame"]["telescope"]["bit_mismatches"] == 1


def test_sign_counts_split_upward_downward_horizontal():
    counts = stream_directions.sign_counts(np.array([0.5, -0.5, 0.0, 0.0, 1.0]))
    assert counts == {"upward": 2, "downward": 1, "horizontal": 2}


@pytest.mark.skipif(shutil.which("obdeect") is None, reason="obdeect CLI not installed")
def test_to_simtel_rejects_upward_and_horizontal_outside_camera(tmp_path):
    rejected, accepted, camera_mixed = stream_directions.rejection_checks(tmp_path)
    assert set(rejected) == set(stream_directions.REJECTING_FRAMES)
    for frame in stream_directions.REJECTING_FRAMES:
        assert not rejected[frame]["converted"], frame
        assert rejected[frame]["message_found"], frame
    assert set(camera_mixed) == {stream_directions.ANY_SIGN_FRAME}
    assert camera_mixed["camera"]["converted"], "camera frame must accept any cz sign"
    for frame in stream_directions.SIMTEL_FRAMES:
        assert accepted[frame]["converted"], frame
        assert accepted[frame]["output_bytes"] > 0, frame
