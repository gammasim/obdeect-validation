import importlib

import numpy as np
import pytest

pytest.importorskip("eventio")
chain_optics = importlib.import_module("obdeect_validation.adapters.chain_optics")


def test_centroid_offsets_zero_for_identical_images():
    event = {"photoelectrons": {"1": {"photoelectrons": [0, 2, 4, 0]}}}
    result = chain_optics.centroid_offsets({"100": event}, {"100": event})
    assert result["centroid_matched_samples"] == 1
    assert result["centroid_max_abs_index_offset"] == 0.0


def test_centroid_offset_sign():
    ref = {"100": {"photoelectrons": {"1": {"photoelectrons": [0, 1, 0, 0]}}}}
    cand = {"100": {"photoelectrons": {"1": {"photoelectrons": [0, 0, 1, 0]}}}}
    result = chain_optics.centroid_offsets(ref, cand)
    assert result["centroid_mean_index_offset"] == pytest.approx(1.0)


def test_trigger_difference_paired(monkeypatch):
    monkeypatch.setattr(chain_optics, "TELESCOPES", (1, 2))
    reference = {
        "100": {"telescope_events": {"1": {}}},
        "101": {"telescope_events": {}},
    }
    candidate = {
        "100": {"telescope_events": {"1": {}, "2": {}}},
        "101": {"telescope_events": {"1": {}}},
    }
    rows = chain_optics.trigger_difference(reference, candidate, ["100", "101"])
    assert rows[1]["difference"] == pytest.approx(0.5)
    assert rows[2]["difference"] == pytest.approx(0.5)
    assert rows[1]["se"] == pytest.approx(np.std([0, 1], ddof=1) / np.sqrt(2))


def test_trigger_difference_identical_samples_not_flagged(monkeypatch):
    monkeypatch.setattr(chain_optics, "TELESCOPES", (1,))
    events = {"100": {"telescope_events": {"1": {}}}, "101": {"telescope_events": {}}}
    rows = chain_optics.trigger_difference(events, events, ["100", "101"])
    assert rows[1] == {"difference": 0.0, "se": 0.0, "flag": 0}
