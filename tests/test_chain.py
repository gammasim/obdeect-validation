import importlib
import itertools
import math

import numpy as np
import pytest

pytest.importorskip("eventio")
chain_optics = importlib.import_module("obdeect_validation.adapters.chain_optics")
atmosphere_chain = importlib.import_module(
    "obdeect_validation.adapters.atmosphere_chain"
)
curved_extinction = importlib.import_module(
    "obdeect_validation.adapters.curved_extinction"
)
extinction_reference = importlib.import_module(
    "obdeect_validation.adapters.extinction_reference"
)

H2_M = 2156.0
SCALE_HEIGHT_M = 8000.0
TAU_INF = 0.3  # saturating vertical depth above H2 (tau increases with height, as in the table)
MONO_BAND = (300.0, 301.0)  # narrow band: every draw gives wavelength 300 nm


def exponential_depth(wavelength, height_m):
    """Synthetic table: tau = 0 at and below H2, TAU_INF (1 - exp(-(h - H2) / H_s)) above."""
    if height_m <= H2_M:
        return 0.0
    return TAU_INF * (1.0 - math.exp(-(height_m - H2_M) / SCALE_HEIGHT_M))


def test_expected_transmission_matches_analytic_telescope_and_focus_terms():
    zem_m, cz, pathlen = 15000.0, math.cos(math.radians(20.0)), 28.0
    telescope = H2_M + 44.0
    mean, _ = atmosphere_chain.expected_transmission(
        [(zem_m, cz, 1.0)],
        exponential_depth,
        MONO_BAND,
        [0.5],
        telescope_altitude_m=telescope,
        pathlen_m=pathlen,
    )

    def tau(height):
        return TAU_INF * (1.0 - math.exp(-(height - H2_M) / SCALE_HEIGHT_M))

    od_tel = tau(telescope)
    od_focus = tau(telescope + pathlen) - od_tel
    expected = math.exp(-(abs(tau(zem_m) - od_tel) / cz + od_focus))
    assert mean == pytest.approx(expected, rel=1e-12)
    previous = math.exp(-tau(zem_m) / cz)
    assert (
        abs(mean - previous) > 1e-4
    )  # the telescope and focus terms shift it by about 5e-4


def test_expected_transmission_reduces_to_airmass_formula_at_h2_without_path():
    samples = [(3000.0, 1.0, 2.0), (9000.0, 0.5, 3.0), (25000.0, 0.8, 1.5)]
    draws = [0.2, 0.5, 0.9]
    mean, _ = atmosphere_chain.expected_transmission(
        samples,
        exponential_depth,
        (250.0, 700.0),
        draws,
        telescope_altitude_m=H2_M,
        pathlen_m=0.0,
    )
    # Previous formula: exp(-tau(zem) / |cz|), with the same wavelength draws.
    weights = []
    total = 0.0
    for (zem_m, cz, photons), u in zip(samples, draws):
        wavelength = math.floor(1.0 / (1.0 / 250.0 - u * (1.0 / 250.0 - 1.0 / 700.0)))
        probability = math.exp(-exponential_depth(wavelength, zem_m) / cz)
        weights.append(photons * probability)
        total += photons
    assert mean == pytest.approx(sum(weights) / total, rel=1e-12)


def test_expected_transmission_default_is_site_level_with_pathlen():
    samples = [(15000.0, 0.9, 1.0), (4000.0, 0.6, 2.0)]
    draws = [0.1, 0.8]
    default, _ = atmosphere_chain.expected_transmission(
        samples, exponential_depth, (250.0, 700.0), draws
    )
    explicit, _ = atmosphere_chain.expected_transmission(
        samples,
        exponential_depth,
        (250.0, 700.0),
        draws,
        telescope_altitude_m=atmosphere_chain.SITE_ALTITUDE_M
        + atmosphere_chain.TELESCOPE_Z_M,
        pathlen_m=atmosphere_chain.PATHLEN_M,
    )
    assert default == pytest.approx(explicit, rel=1e-15)


def test_expected_transmission_decreases_with_path_and_emission_height():
    cz = math.cos(math.radians(20.0))
    telescope = H2_M + 10.0

    def mean_for(zem_m, pathlen):
        return atmosphere_chain.expected_transmission(
            [(zem_m, cz, 1.0)],
            exponential_depth,
            MONO_BAND,
            [0.5],
            telescope_altitude_m=telescope,
            pathlen_m=pathlen,
        )[0]

    by_path = [mean_for(15000.0, path) for path in [0.0, 10.0, 28.0, 60.0]]
    assert all(later < earlier for earlier, later in itertools.pairwise(by_path))
    by_height = [
        mean_for(height, 28.0) for height in [2300.0, 5000.0, 15000.0, 40000.0]
    ]
    assert all(later < earlier for earlier, later in itertools.pairwise(by_height))


def test_table_reference_follows_curved_extinction_convention(tmp_path):
    # Convention shared by sim_telarray and the chain: tau = 0 at and below H2, linear in
    # log10(h) between table heights, constant above the last height.
    text = (
        "# H2= 2.156, H1=    2.206     2.256     2.356\n"
        "300 0.10 0.20 0.35\n"
        "301 0.05 0.12 0.30\n"
    )
    path = tmp_path / "table.dat"
    path.write_text(text, encoding="utf-8")
    depth = extinction_reference.table_reference(path)[1]
    expected = curved_extinction.vertical_depth(
        curved_extinction.parse_table(text),
        [1000.0, 2156.0, 2180.0, 2206.0, 2230.0, 2300.0, 2356.0, 5000.0, 20000.0],
    )
    heights = [1000.0, 2156.0, 2180.0, 2206.0, 2230.0, 2300.0, 2356.0, 5000.0, 20000.0]
    for column, height in enumerate(heights):
        assert depth(300, height) == pytest.approx(expected[0, column], abs=1e-12)
        assert depth(301, height) == pytest.approx(expected[1, column], abs=1e-12)
    assert depth(300, 1000.0) == 0.0
    assert depth(300, 20000.0) == pytest.approx(0.35)


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
