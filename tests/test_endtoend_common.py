"""Unit tests of the end-to-end helpers: verdicts, pipelines with timeouts, exact photon comparison,
efficiency status edits and the byte-level stream damage of the robustness case."""

import shutil
import struct

import numpy as np
import pytest

from obdeect_validation.adapters import endtoend_common as e2e


def stage(returncode=0, stderr="", timed_out=False):
    return e2e.Stage(
        argv=["x"], returncode=returncode, stderr=stderr, timed_out=timed_out
    )


@pytest.mark.parametrize(
    ("kwargs", "verdict"),
    [
        ({"returncode": 0}, "ok"),
        ({"returncode": 1, "stderr": "obdeect: error: bad input"}, "error"),
        ({"returncode": 2}, "error"),
        ({"returncode": -11}, "crash"),
        ({"returncode": 139}, "crash"),
        (
            {"returncode": 1, "stderr": "Traceback (most recent call last):\n  File"},
            "crash",
        ),
        ({"returncode": None}, "crash"),
        ({"returncode": -9, "timed_out": True}, "hang"),
    ],
)
def test_classify(kwargs, verdict):
    assert e2e.classify(stage(**kwargs)) == verdict


def test_pipeline_ok_and_output(tmp_path):
    if shutil.which("printf") is None:
        pytest.skip("printf not available")
    stages = e2e.pipeline(
        [["printf", "abc"], ["cat"]],
        output=tmp_path / "out",
        log_dir=tmp_path / "logs",
        timeout=30,
    )
    assert [e2e.classify(s) for s in stages] == ["ok", "ok"]
    assert (tmp_path / "out").read_bytes() == b"abc"


def test_pipeline_hang_is_killed(tmp_path):
    stages = e2e.pipeline([["sleep", "30"]], log_dir=tmp_path / "logs", timeout=0.5)
    assert e2e.classify(stages[0]) == "hang"


def test_pipeline_early_consumer_exit_does_not_hang(tmp_path):
    if shutil.which("yes") is None:
        pytest.skip("yes not available")
    stages = e2e.pipeline(
        [["yes"], ["head", "-c", "10"]],
        output=tmp_path / "out",
        log_dir=tmp_path / "logs",
        timeout=30,
    )
    assert not any(s.timed_out for s in stages)
    content = (tmp_path / "out").read_bytes()
    assert content == b"y\ny\ny\ny\ny\n"


def photons(x, cx=None, event=1, telescope=1):
    n = len(x)
    return {
        "event": np.full(n, event, dtype=np.int64),
        "telescope": np.full(n, telescope, dtype=np.int64),
        "index": np.arange(n, dtype=np.int64),
        "x": np.asarray(x, dtype=np.float32),
        "cx": np.asarray(cx if cx is not None else np.zeros(n), dtype=np.float32),
    }


def test_compare_rows_identical_and_one_bit_difference():
    left = photons([1.0, 2.0, 3.0])
    assert e2e.compare_rows(left, photons([1.0, 2.0, 3.0]))["photon_mismatches"] == 0
    one_ulp = np.nextafter(np.float32(2.0), np.float32(3.0))
    result = e2e.compare_rows(left, photons([1.0, one_ulp, 3.0]))
    assert result["photon_mismatches"] == 1
    assert result["photons_compared"] == 3


def test_compare_rows_sign_of_zero_and_nan_are_bit_differences():
    assert e2e.compare_rows(photons([0.0]), photons([-0.0]))["photon_mismatches"] == 1
    nan_a = np.float32(np.nan)
    nan_b = np.frombuffer(struct.pack("<I", 0x7FC00001), dtype=np.float32)[0]
    assert (
        e2e.compare_rows(photons([nan_a]), photons([nan_b]))["photon_mismatches"] == 1
    )


def test_compare_rows_counts_photons_on_one_side_only():
    result = e2e.compare_rows(photons([1.0, 2.0, 3.0]), photons([1.0]))
    assert result["photon_mismatches"] == 2
    assert result["photons_compared"] == 3


def test_compare_rows_detects_other_event_and_telescope():
    assert (
        e2e.compare_rows(photons([1.0], event=1), photons([1.0], event=2))[
            "photon_mismatches"
        ]
        == 1
    )
    assert (
        e2e.compare_rows(photons([1.0], telescope=1), photons([1.0], telescope=2))[
            "photon_mismatches"
        ]
        == 1
    )


def test_compare_rows_column_sets():
    result = e2e.compare_rows(
        photons([1.0]), {**photons([1.0]), "cy": np.zeros(1, dtype=np.float32)}
    )
    assert (
        result["column_sets_equal"] == 0
    )  # different column sets are reported, not silently compared


def test_drop_events_removes_only_those_events():
    rows = {**photons([1.0, 2.0], event=1), "_chunks": 2}
    rows2 = {**photons([3.0], event=2)}
    merged = {
        k: np.concatenate([rows[k], rows2[k]]) for k in rows2 if not k.startswith("_")
    }
    merged["_chunks"] = 2
    kept = e2e.drop_events(merged, [2])
    assert len(kept["x"]) == 2
    assert set(kept["event"].tolist()) == {1}
    assert kept["_chunks"] == 2


def test_status_helpers():
    meta = {"factors": "atmosphere", "atmosphere.state": "applied"}
    added = e2e.with_status_factor(meta, "qe", "envelope", "src")
    assert added["factors"] == "atmosphere,qe"
    assert added["qe.state"] == "envelope"
    assert "qe.state" not in meta  # input not modified
    removed = e2e.without_factor(added, "atmosphere")
    assert removed["factors"] == "qe"
    assert "atmosphere.state" not in removed
    assert removed["qe.state"] == "envelope"


def test_summary_figure_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    target = tmp_path / "summary.png"
    e2e.summary_figure(target, {"a_photon_mismatch": 0, "b_accepted": 1}, "test case")
    assert target.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_robustness_byte_patches_on_synthetic_stream():
    pytest.importorskip("eventio")
    from obdeect_validation.adapters import robustness_stream as rs

    header = rs.HEADER
    stream = bytearray()
    stream += header.pack(rs.MAGIC, 1, 0, 0)  # empty stream record, length 0
    payload = struct.pack("<IIII", 1, 3, 0b1, 4) + struct.pack("<Q", 0) + bytes(16)
    stream += header.pack(rs.MAGIC, rs.PHOTONS, 0, len(payload)) + payload
    assert [kind for _, kind, _, _ in rs.records(stream)] == [1, rs.PHOTONS]
    _, first = rs.first_photons(stream)
    assert (
        first == header.size + header.size
    )  # payload of the photons record (after the 16-byte header)

    rows = bytearray(stream)
    rs.rows_too_large(rows)
    assert struct.unpack_from("<I", rows, first + 12)[0] == 70000

    magic = bytearray(stream)
    rs.record_magic(magic)
    assert bytes(magic[header.size : header.size + 4]) == b"XXXX"
    assert [kind for _, kind, _, _ in rs.records(magic)] == [
        1
    ]  # reading stops at the bad magic

    cut = bytearray(stream)
    rs.cut_header(cut)
    assert len(cut) == 20


def test_robustness_format_major_patch():
    pytest.importorskip("eventio")
    from obdeect_validation.adapters import robustness_stream as rs

    data = bytearray(b"xx" + b"format_major" + b"\x01" + struct.pack("<q", 1) + b"yy")
    rs.format_major(data)
    assert struct.unpack_from("<q", data, 2 + len(b"format_major") + 1)[0] == 2
