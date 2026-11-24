"""Unit tests of the rejection adapter: input builders, byte edits and the result summary."""

import json
import struct
from pathlib import Path

import numpy as np
import pytest

from obdeect_validation.adapters import robustness_rejections as rr

CASE = Path(__file__).resolve().parents[1] / "cases/robustness/stream-rejections.json"


def _reader_available():
    try:
        import obdeect.stream  # noqa: F401
    except ImportError:
        return False
    return True


needs_obdeect = pytest.mark.skipif(
    not _reader_available(), reason="obdeect.stream not installed"
)


def test_eventio_block_layout():
    block = rr.eventio_block(1202, 7, b"abcde")
    sync, type_version, ident, length = struct.unpack("<IIII", block[:16])
    assert sync == rr.EVENTIO_SYNC
    assert type_version == 1202
    assert ident == 7
    assert length == 8  # padded to 4 bytes
    assert block[16:] == b"abcde\0\0\0"


def test_nested_block_has_no_sync():
    block = rr.eventio_block(1210, 0, b"", top_level=False)
    assert len(block) == 12


def test_float_block_count_and_values():
    body = rr.float_block([1.0, 2.5])
    assert struct.unpack("<i", body[:4])[0] == 2
    assert struct.unpack("<2f", body[4:]) == (1.0, 2.5)


def test_corsika_run_blocks():
    open_run = rr.corsika_run(1, run_end=False)
    closed = rr.corsika_run(1)
    types = [struct.unpack("<I", b[4:8])[0] & 0xFFFF for b in closed]
    assert types == [1200, 1202, 1209, 1210]
    assert [struct.unpack("<I", b[4:8])[0] & 0xFFFF for b in open_run] == [
        1200,
        1202,
        1209,
    ]


def test_format_minor_edit_keeps_length(tmp_path):
    if not _reader_available():
        pytest.skip("obdeect.stream not installed")
    path = rr.write_stream(
        tmp_path / "s.pst", rr.valid_ops(rr.chunk_columns(rr.CAMERA_NAMES))
    )
    size = path.stat().st_size
    rr.set_int_field(path, "format_minor", 0)
    assert path.stat().st_size == size
    from obdeect.stream import Reader

    first = next(iter(Reader(str(path))))
    assert first["metadata"]["format_minor"] == 0


def test_chunk_columns_unit_direction():
    cols = rr.chunk_columns(rr.CAMERA_NAMES)
    norm = cols["cx"] ** 2 + cols["cy"] ** 2
    assert np.all(norm < 1.0)
    assert "cz" not in cols
    upward = rr.chunk_columns(rr.CAMERA_NAMES + ("cz",), upward=True)
    assert np.all(upward["cz"] > 0)
    np.testing.assert_allclose(upward["cz"] ** 2 + norm, 1.0, rtol=1e-6)


@needs_obdeect
def test_nonunit_direction_is_rejected_on_read(tmp_path):
    # The writer writes the value as given; the reader rejects a direction that is not a unit vector.
    names = rr.CAMERA_NAMES + ("cz",)
    path = rr.write_stream(
        tmp_path / "n.pst", rr.valid_ops(rr.chunk_columns(names, cz_value=-0.5))
    )
    from obdeect.stream import Reader

    with pytest.raises(RuntimeError, match="unit vector"):
        list(Reader(str(path)))


@needs_obdeect
@pytest.mark.parametrize("column", ["x", "y", "t", "w"])
def test_nan_variant_is_rejected_on_read(tmp_path, column):
    # The writer writes NaN as given; the reader rejects non-finite column values.
    cols = rr.chunk_columns(rr.CAMERA_NAMES)
    cols[column][2] = np.float32("nan")
    path = rr.write_stream(tmp_path / "nan.pst", rr.valid_ops(cols))
    from obdeect.stream import Reader

    with pytest.raises(RuntimeError, match="not finite"):
        list(Reader(str(path)))


@needs_obdeect
def test_raw_stream_reads_like_the_writer(tmp_path):
    from obdeect.stream import Reader

    ops = rr.valid_ops(rr.chunk_columns(rr.CAMERA_NAMES))
    written = list(Reader(str(rr.write_stream(tmp_path / "w.pst", ops))))
    raw = list(Reader(str(rr.write_raw_stream(tmp_path / "r.pst", ops))))
    assert [r["type"] for r in raw] == [r["type"] for r in written]
    assert raw[1]["metadata"]["run_number"] == written[1]["metadata"]["run_number"]
    assert set(raw[1]["metadata"]) == set(written[1]["metadata"])


@needs_obdeect
def test_raw_nesting_violation_is_rejected_by_the_reader(tmp_path):
    from obdeect.stream import Reader

    ops = [
        ("run", rr.RUN),
        ("event", rr.event_meta(1)),
        ("event", rr.event_meta(2)),  # event inside an event: the Writer refuses it
        ("event_end", {}),
        ("event_end", {}),
        ("run_end", {}),
    ]
    path = rr.write_raw_stream(tmp_path / "nested.pst", ops)
    with pytest.raises(RuntimeError, match="malformed stream"):
        list(Reader(str(path)))


def test_summarize_counts_variants_and_controls():
    results = [
        {
            "scenario": "nesting",
            "name": "a",
            "consumer": "cat",
            "role": "variant",
            "accepted": True,
        },
        {
            "scenario": "nesting",
            "name": "a",
            "consumer": "reader",
            "role": "variant",
            "accepted": False,
        },
        {
            "scenario": "missing_wavelength",
            "name": "c",
            "consumer": "to-simtel",
            "role": "control",
            "accepted": False,
        },
        {
            "scenario": "unknown_major",
            "name": "major-2",
            "consumer": "reader",
            "role": "variant",
            "accepted": False,
        },
        {
            "scenario": "unknown_major",
            "name": "major-2-later-calls",
            "consumer": "reader",
            "role": "variant",
            "accepted": True,
        },
    ]
    metrics = rr.summarize(results)
    assert metrics["nesting_violations_accepted"] == 1
    assert metrics["control_rejected"] == 1
    assert metrics["reader_first_call_accepted"] == 0
    assert metrics["reader_continues_after_rejection"] == 1
    assert metrics["write_failure_exit_zero"] == 0


def test_case_metrics_match_adapter_output():
    case = json.loads(CASE.read_text())
    expected = set(rr.summarize([]).keys())
    listed = {m["metric"] for m in case["metrics"]}
    assert listed == expected
    assert all(m["expected"] == 0 and m["atol"] == 0 for m in case["metrics"])
