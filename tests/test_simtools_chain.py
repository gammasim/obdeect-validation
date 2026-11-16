"""Tests of the simtools chain comparison helpers (no simulation needed)."""

import gzip

import pytest

from obdeect_validation.adapters import chain_optics, simtools_chain


def test_triggered_accepts_integer_telescope_keys():
    # eventio keys telescope_events by int
    events = {"100": {"telescope_events": {1: object()}}}

    assert chain_optics.triggered(events, "100", 1)
    assert not chain_optics.triggered(events, "101", 1)


def test_triggered_accepts_string_telescope_keys():
    events = {"100": {"telescope_events": {"1": object()}}}

    assert chain_optics.triggered(events, "100", 1)


def test_triggered_without_telescope_events_is_false():
    assert not chain_optics.triggered({"100": {}}, "100", 1)


def test_trigger_difference_sees_one_sided_trigger():
    reference = {"100": {"telescope_events": {1: object()}}, "101": {}}
    candidate = {"100": {}, "101": {"telescope_events": {1: object()}}}

    rows = chain_optics.trigger_difference(reference, candidate, ["100", "101"])

    assert rows[1]["difference"] == 0.0
    assert rows[1]["se"] == pytest.approx(0.5**0.5)  # paired: sqrt(2) / 2


def test_atmosphere_check_accepts_expected_table_only(tmp_path):
    good = tmp_path / "good.err.gz"
    with gzip.open(good, "wt") as handle:
        handle.write(
            "Read atmospheric transmission data from file /r/atmosphere-zero-transmission.dat\n"
        )
    bad = tmp_path / "bad.err.gz"
    with gzip.open(bad, "wt") as handle:
        handle.write("Read atmospheric transmission data from file /r/other.dat\n")

    assert (
        simtools_chain.atmosphere_check(good, "atmosphere-zero-transmission.dat") == 0
    )
    assert simtools_chain.atmosphere_check(bad, "atmosphere-zero-transmission.dat") == 1


def test_atmosphere_check_without_read_line_is_a_failure(tmp_path):
    empty = tmp_path / "empty.err.gz"
    with gzip.open(empty, "wt") as handle:
        handle.write("no table here\n")

    assert simtools_chain.atmosphere_check(empty, "x.dat") == 1
