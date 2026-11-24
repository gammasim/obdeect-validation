"""Unit tests of the event-number robustness adapter: the source geometry, the binomial pull and
the identity and missing-output rules. No obdeect run is needed."""

import math

import numpy as np
import pytest

from obdeect_validation.adapters import robustness_event_numbers as adapter


def columns(values, shift=0.0):
    x = np.asarray(values, dtype=np.float32)
    return {
        "x": x,
        "y": x * 0.5,
        "w": np.ones(len(x), dtype=np.float32),
        "t": (x + np.float32(shift)).astype(np.float32),
    }


def test_event_list_covers_the_required_numbers_and_stops_at_the_allowed_maximum():
    required = {1, 2, 255, 256, 300, 65535, 65536, 100000, 2**31 - 1, 2**32 - 1}
    assert set(adapter.EVENTS) == required
    assert max(adapter.EVENTS) == 2**32 - 1
    assert adapter.BASELINE == 1


def test_source_lies_on_the_pointing_axis():
    x, y, z = adapter.source_position_m()
    distance = math.sqrt(x * x + y * y + z * z)
    assert distance == pytest.approx(adapter.SOURCE_DISTANCE_M)
    # pointing zenith 20 deg, azimuth 180 deg: the source is on the optical axis, south of the telescope
    assert x == pytest.approx(-1000.0 * math.sin(math.radians(20.0)), abs=1e-9)
    assert abs(y) < 1e-9
    assert z == pytest.approx(1000.0 * math.cos(math.radians(20.0)))


def test_binomial_fraction_and_variance():
    fraction, variance = adapter.binomial_fraction(60.0, 100.0)
    assert fraction == pytest.approx(0.6)
    assert variance == pytest.approx(0.6 * 0.4 / 100.0)
    with pytest.raises(ValueError):
        adapter.binomial_fraction(0.0, 0.0)


def test_pull_is_zero_for_equal_fractions_and_matches_the_formula():
    assert adapter.pull(0.5, 0.01, 0.5, 0.01) == 0.0
    assert adapter.pull(0.6, 0.0004, 0.5, 0.0001) == pytest.approx(
        0.1 / math.sqrt(0.0005)
    )
    assert adapter.pull(0.0, 0.0, 0.0, 0.0) == 0.0
    assert adapter.pull(0.2, 0.0, 0.1, 0.0) == adapter.PULL_MISSING


def test_pull_maximum_is_over_the_events_against_the_baseline():
    fractions = {1: (0.5, 0.0001), 2: (0.51, 0.0001), 3: (0.55, 0.0001)}
    value, pulls = adapter.detected_fraction_pull_max(fractions)
    assert set(pulls) == {2, 3}
    assert value == pytest.approx(pulls[3])
    assert value == pytest.approx(0.05 / math.sqrt(0.0002))


def test_missing_output_fails_the_pull_gate():
    value, pulls = adapter.detected_fraction_pull_max({1: (0.5, 0.0001), 256: None})
    assert value == adapter.PULL_MISSING
    assert pulls == {256: None}
    value, _ = adapter.detected_fraction_pull_max({1: None, 2: (0.5, 0.0001)})
    assert value == adapter.PULL_MISSING
    assert value > adapter.PULL_GATE


def test_identical_pairs_need_bit_identical_photons_and_photons():
    same = columns([1.0, 2.0, 3.0])
    other = columns([1.0, 2.0, 3.5])
    rows = {
        1: _rows(same),
        2: _rows(dict(same)),
        3: _rows(other),
        4: None,
        5: _rows(columns([])),
    }
    assert adapter.identical_pairs(rows) == [(1, 2)]


def test_identical_pairs_see_a_one_ulp_difference():
    base = columns([1.0, 2.0, 3.0])
    shifted = dict(base)
    shifted["t"] = np.nextafter(base["t"], np.float32(np.inf)).astype(np.float32)
    assert adapter.identical_pairs({1: _rows(base), 2: _rows(shifted)}) == []


def test_count_and_centroid_of_camera_photons():
    rows = _rows(columns([1.0, 3.0]))
    assert adapter.count_photons(rows) == 2
    assert adapter.count_photons(None) == 0
    assert adapter.centroid(rows) == pytest.approx([2.0, 1.0])
    assert adapter.centroid(None) is None


def _rows(cols):
    """photon_rows-like dict: row keys plus the columns."""
    count = len(cols["x"])
    rows = {
        "event": np.ones(count, dtype=np.int64),
        "telescope": np.ones(count, dtype=np.int64),
        "index": np.arange(count, dtype=np.int64),
        "_chunks": 1,
    }
    rows.update(
        {name: np.asarray(value, dtype=np.float32) for name, value in cols.items()}
    )
    return rows
