"""Pure helpers of the curved-shower extinction and chain adapters (no obdeect, CORSIKA or sim_telarray)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("eventio")

from obdeect_validation.adapters import chain_curvature as chain
from obdeect_validation.adapters import curvature_corsika_extinction as adapter
from obdeect_validation.adapters import curved_extinction as ce
from obdeect_validation.adapters import curved_showers

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"


def _table():
    return ce.parse_table((FIXTURES / "atm_trans_2000_1_10_0_0_2000.dat").read_text())


def test_curved_cases_match_fixture_cards():
    cases = adapter.curved_cases()
    assert [case[0] for case in cases] == [65, 75]
    assert cases[0][1:] == (300, 5000)
    assert all(low == 1000 and high == 10000 for _, low, high in cases[1:])


def test_telescope_altitude_is_site_plus_card_telescope_height():
    # The legacy card places the telescopes at z = 0 (observation level), so h_t = 2156 m.
    assert adapter.telescope_altitude_m() == pytest.approx(2156.0)
    assert adapter.telescope_altitude_m() == pytest.approx(
        2156.0 + curved_showers.layout()["telescopes"][0][2] / 100.0
    )


def test_check_inputs_names_every_missing_file(tmp_path):
    with pytest.raises(adapter.MissingInputError) as info:
        adapter.check_inputs(tmp_path)
    message = str(info.value)
    for zenith, low, high in adapter.curved_cases():
        name = curved_showers.shower_name(
            zenith, "curved", energy_min=low, energy_max=high
        )
        assert name in message
    # Present showers but no CORSIKA logs: still missing.
    for zenith, low, high in adapter.curved_cases():
        (
            tmp_path
            / curved_showers.shower_name(
                zenith, "curved", energy_min=low, energy_max=high
            )
        ).touch()
    with pytest.raises(adapter.MissingInputError, match="corsika.log"):
        adapter.check_inputs(tmp_path)
    for zenith, low, high in adapter.curved_cases():
        log = adapter.shower_file(tmp_path, "curved", zenith, low, high).with_suffix(
            ".corsika.log"
        )
        log.touch()
    adapter.check_inputs(tmp_path)


def test_ratio_violations_counts_below_one_and_falls():
    zeniths = [75, 65]
    assert adapter.ratio_violations(zeniths, [1.1, 1.001]) == 0
    assert adapter.ratio_violations([65, 75], [1.01, 1.00]) == 1  # fall with zenith
    assert adapter.ratio_violations([65, 75], [0.9, 1.01]) == 1  # ratio below one
    assert adapter.ratio_violations([65, 75], [1.0, 1.0]) == 0


def test_tolerance_floor_and_five_standard_errors():
    assert adapter.tolerance(0.0) == pytest.approx(1.0e-3)
    assert adapter.tolerance(1.0e-4) == pytest.approx(1.0e-3)
    assert adapter.tolerance(1.0e-3) == pytest.approx(5.0e-3)


def test_expectation_is_photon_weighted_and_ordered():
    table = _table()
    band = (250.0, 800.0)
    zem = np.array([6000.0, 6000.0, 12000.0])
    cz = np.array([0.4, 0.4, 0.3])
    photons = np.array([3.0, 7.0, 0.0])
    h_t = 2167.05
    single = adapter.expectation(zem[:1], cz[:1], photons[:1], table, band, h_t)
    doubled = adapter.expectation(zem[:2], cz[:2], photons[:2], table, band, h_t)
    assert doubled["curved"] == pytest.approx(single["curved"])
    assert doubled["plane"] == pytest.approx(single["plane"])
    assert doubled["photon_total"] == pytest.approx(10.0)
    assert doubled["predicted_ratio"] >= 1.0
    assert doubled["curved"] >= doubled["plane"]
    assert doubled["curved_standard_error"] >= 0.0


def test_stage_total_requires_the_stage():
    stats = {"Detected": [3, 4], "Not absorbed/max. Q.E.": [10, 12]}
    assert chain.stage_total(stats, "Detected") == 7.0
    with pytest.raises(ValueError, match="Missing"):
        chain.stage_total(stats, "Missing")


def test_event_and_trigger_comparison():
    left = {
        "1": {"trigger_information": {"trigger_pattern": 3}},
        "2": {"trigger_information": {"trigger_pattern": 1}},
    }
    right = {
        "2": {"trigger_information": {"trigger_pattern": 2}},
        "3": {"trigger_information": {"trigger_pattern": 1}},
    }
    assert chain.event_id_delta(left, right) == 2
    assert chain.trigger_mismatches(left, right) == 1
    assert chain.trigger_mismatches(left, left) == 0


def test_ratio_check_accepts_agreement_and_flags_offsets():
    agree = chain.ratio_check(1.2, 1.2, 0.0, 1.0e6, 1.0e6)
    assert agree["limit_exceeded"] == 0
    assert agree["limit"] == pytest.approx(
        0.01
    )  # the 1 % floor dominates at these counts
    off = chain.ratio_check(1.2 * 1.05, 1.2, 0.0, 1.0e6, 1.0e6)
    assert off["limit_exceeded"] == 1


def test_ratio_check_limit_grows_with_few_counts_and_draws():
    result = chain.ratio_check(1.03, 1.0, 0.0, 400.0, 400.0)
    assert result["limit"] > 0.01
    assert result["limit_exceeded"] == 0
    drawn = chain.ratio_check(1.0, 1.0, 0.02, 1.0e9, 1.0e9)
    assert drawn["limit"] == pytest.approx(0.1, rel=1e-4)  # 5 x 2 % draw error


def test_plane_check_would_flag_a_curved_sim_telarray():
    # A curved sim_telarray lying 10 % from the plane prediction must be flagged.
    plane_count = 0.45 * 1.0e6
    agree = chain.ratio_check(1.001, 1.0, 1e-5, plane_count, 1.0e6)
    assert agree["limit_exceeded"] == 0
    curved_like = chain.ratio_check(1.10, 1.0, 1e-5, plane_count, 1.0e6)
    assert curved_like["limit_exceeded"] == 1
    assert curved_like["relative_difference"] == pytest.approx(0.10, rel=1e-6)


def test_read_two_column_ignores_comments_and_extra_columns(tmp_path):
    table = tmp_path / "table.dat"
    table.write_text("# header\n300 0.5 0.1  # note\n200 0.2\nbad line\n400 0.6\n")
    data = chain.read_two_column(table)
    assert data.tolist() == [[200.0, 0.2], [300.0, 0.5], [400.0, 0.6]]
    unsorted = tmp_path / "unsorted.dat"
    unsorted.write_text("300 0.1\n300 0.2\n")
    with pytest.raises(ValueError):
        chain.read_two_column(unsorted)


def test_detection_weights_are_quantum_efficiency_times_reflectivity(tmp_path):
    (tmp_path / chain.QE_TABLE).parent.mkdir(parents=True)
    (tmp_path / chain.QE_TABLE).write_text("200 0.1\n400 0.4\n800 0.0\n")
    (tmp_path / chain.MIRROR_TABLE).parent.mkdir(parents=True, exist_ok=True)
    (tmp_path / chain.MIRROR_TABLE).write_text("200 0.5 0.0\n1000 0.9 0.0\n")
    (tmp_path / chain.WINDOW_TABLE).parent.mkdir(parents=True, exist_ok=True)
    (tmp_path / chain.WINDOW_TABLE).write_text("200 1.0\n1000 1.0\n")  # open window
    weights = chain.DetectionWeights(tmp_path)
    assert weights.quantum([300.0, 900.0]).tolist() == pytest.approx([0.25, 0.0])
    # 400 nm: QE 0.4, reflectivity 0.5 + (200 / 800) * 0.4 = 0.6 (linear between 200 and 1000 nm)
    detection = weights.detection([400.0])
    assert detection[0] == pytest.approx(0.4 * 0.6)
    # Outside the quantum table the weight is zero by definition, even without mirror data there.
    assert weights.detection([1200.0])[0] == 0.0
    # Non-zero quantum efficiency beyond the mirror table is an error.
    short = tmp_path / "short"
    for name, text in (
        (chain.QE_TABLE, "200 0.1\n400 0.4\n800 0.0\n"),
        (chain.MIRROR_TABLE, "200 0.5 0.0\n300 0.9 0.0\n"),
        (chain.WINDOW_TABLE, "200 1.0\n1000 1.0\n"),
    ):
        (short / name).parent.mkdir(parents=True, exist_ok=True)
        (short / name).write_text(text)
    with pytest.raises(ValueError, match="does not cover"):
        chain.DetectionWeights(short).detection([400.0])


def test_ratio_violations_and_helpers_are_pure():
    assert adapter.ratio_violations([], []) == 0
    assert adapter.LIMITS["predicted_ratio"][0].startswith("curved / plane")
    assert chain.LIMITS["simtel_plane_parallel"][0].startswith("|N_ref")
