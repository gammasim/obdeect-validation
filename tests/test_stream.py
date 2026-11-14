"""Unit tests for the stream adapters' pure helpers (needs obdeect and numpy)."""

import numpy as np
from obdeect.stream import Writer

from obdeect_validation.adapters import stream_chain

COLUMNS = ["x", "y", "cx", "cy", "t", "zem", "w", "wavelength"]


def _write(path, order, shift=0.0):
    values = np.arange(12, dtype=np.float32) + shift
    columns = {name: (values + index)[order] for index, name in enumerate(COLUMNS)}
    with Writer(str(path), "test") as writer:
        writer.write("run", {"run_number": 1})
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1, "frame": "camera"})
        writer.write_photons(1, "camera", 0, columns)
        writer.write("telescope_end", {"telescope_id": 1, "photons": len(order)})
        writer.write("event_end", {})
        writer.write("run_end", {})
    return path


def test_canonical_ignores_photon_order(tmp_path):
    ordered = _write(tmp_path / "ordered.gst", np.arange(12))
    shuffled = _write(
        tmp_path / "shuffled.gst", np.random.default_rng(3).permutation(12)
    )
    digests = stream_chain.canonical(ordered)
    assert list(digests) == ["0:1"]
    assert stream_chain.differing(digests, stream_chain.canonical(shuffled)) == 0


def test_canonical_detects_changed_photon(tmp_path):
    reference = stream_chain.canonical(_write(tmp_path / "a.gst", np.arange(12)))
    shifted = stream_chain.canonical(
        _write(tmp_path / "b.gst", np.arange(12), shift=1.0)
    )
    assert (
        stream_chain.differing(reference, shifted) == 2
    )  # changed key counts on both sides


def test_differing_counts_missing_keys():
    assert stream_chain.differing({"0:1": "a", "0:2": "b"}, {"0:1": "a"}) == 1
    assert stream_chain.differing({"0:1": "a"}, {"0:1": "b"}) == 2
