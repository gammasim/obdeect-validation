"""Unit tests for the IACT summary helper (reader replaced by synthetic bunches)."""

import numpy as np
import pytest

pytest.importorskip("obdeect.stream")
summarize_iact = pytest.importorskip("obdeect_validation.adapters.summarize_iact")


def test_summary_over_all_bunches(monkeypatch):
    data = {
        "x": np.array([0.0, 2.0]),
        "y": np.array([1.0, 1.0]),
        "wavelength": np.array([300.0, 500.0]),
    }
    monkeypatch.setattr(
        summarize_iact.readers, "read_iact", lambda path: ([3.0, 1.0], data)
    )
    summary = summarize_iact.summarise("unused.iact")
    assert summary["events"] == 2
    assert summary["photons"] == pytest.approx(4.0)
    assert summary["bunches"] == 2
    assert summary["wavelength_mean_nm"] == pytest.approx(400.0)
    assert summary["x_mean_m"] == pytest.approx(1.0)
    assert summary["x_rms_m"] == pytest.approx(1.0)
    assert summary["y_rms_m"] == pytest.approx(0.0)


def test_empty_file_has_no_bunch_statistics(monkeypatch):
    empty = {"x": np.empty(0), "y": np.empty(0), "wavelength": np.empty(0)}
    monkeypatch.setattr(
        summarize_iact.readers, "read_iact", lambda path: ([0.0], empty)
    )
    summary = summarize_iact.summarise("unused.iact")
    assert summary == {"events": 1, "photons": 0.0, "bunches": 0}
    assert "bunches 0" in summarize_iact.format_summary(summary)
