import math

import pytest

from obdeect_validation.adapters import atmosphere_chain as chain

TABLE = """# H2= 2.156, H1=    2.206     2.256
# Atmospheric model: 1 (Tropical atmosphere)
200 0.1 0.2
999 0.3 0.4
1000 0.5 0.6
"""


def test_cherenkov_band_reads_corsika_log_line():
    text = " CHERENKOV WAVELENGTH BAND FROM    100.000000      TO    2000.00000      NANOMETER\n"
    assert chain.cherenkov_band(text) == (100.0, 2000.0)


def test_cherenkov_band_requires_the_line():
    with pytest.raises(ValueError):
        chain.cherenkov_band("no band here")


def test_cut_table_removes_rows_at_or_above_limit_and_keeps_header():
    cut = chain.cut_table(TABLE, 1000)
    assert cut.splitlines()[:2] == TABLE.splitlines()[:2]
    assert [line.split()[0] for line in cut.splitlines()[2:]] == ["200", "999"]


def test_zero_table_keeps_structure_with_unit_transmission():
    zero = chain.zero_table(TABLE).splitlines()
    assert zero[0] == TABLE.splitlines()[0]
    assert zero[2:] == ["200 0 0", "999 0 0", "1000 0 0"]


def test_photon_statistics_parses_block_up_to_trigger_section():
    log = (
        "Photon statistics:\n"
        "  All photons:             18605134  17093240  18511833\n"
        "  Not absorbed/max. Q.E.:    703464    656510    678727\n"
        "  Detected:                  325062    300780    328071\n"
        "Trigger statistics:\n"
        "  Tel. triggered:    1 2 3\n"
    )
    stats = chain.photon_statistics(log)
    assert stats["All photons"] == [18605134, 17093240, 18511833]
    assert stats["Not absorbed/max. Q.E."] == [703464, 656510, 678727]
    assert stats["Detected"] == [325062, 300780, 328071]
    assert "Tel. triggered" not in stats


def test_propagation_summary_parses_obdeect_message():
    log = "obdeect propagation run: photons in 10888516, out 5391964 (weight 1.32392e+07), dropped 0\n"
    assert chain.propagation_summary(log) == {
        "photons_in": 10888516, "photons_out": 5391964, "weight_out": 1.32392e7, "dropped": 0}


def test_limits_use_five_sigma_or_floor():
    assert chain.ratio_limit(1e6, 1e6) == pytest.approx(max(0.01, 5 * math.sqrt(2e-6)))
    assert chain.ratio_limit(1e12, 1e12) == 0.01
    assert chain.limit_exceeded(0.011, 0.01) == 1
    assert chain.limit_exceeded(0.009, 0.01) == 0


def test_expected_transmission_with_constant_depth_is_exponential_in_airmass():
    samples = [(0.0, 1.0, 2.0), (0.0, 0.5, 3.0)]  # zem_m, |cz|, photons
    draws = [0.3, 0.7]
    mean, se = chain.expected_transmission(samples, lambda wavelength, height: 0.2, (250.0, 700.0), draws)
    expected = (2.0 * math.exp(-0.2) + 3.0 * math.exp(-0.4)) / 5.0
    assert mean == pytest.approx(expected)
    assert se >= 0.0


def test_expected_transmission_zero_above_max_lambda():
    samples = [(0.0, 1.0, 1.0)]
    mean, _ = chain.expected_transmission(samples, lambda wavelength, height: 0.0, (1000.0, 2000.0), [0.0])
    assert mean == 0.0
