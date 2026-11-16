"""Unit tests of the pass/fail arithmetic of the mode-3 envelope chain (no sim_telarray needed)."""

from pathlib import Path

import numpy as np
import pytest

from obdeect_validation.adapters import chain_mode3 as c3
from obdeect_validation.adapters import performance_chain as pc


def test_clustered_pe_bias_of_a_constant_ratio_has_no_scatter() -> None:
    # Candidate is 2 % above the reference in every shower: bias 0.02, residuals zero.
    rows = [(0.02 * ref, ref, None) for ref in (10.0, 30.0, 60.0)]
    bias, se = c3.clustered_pe_bias(rows)
    assert bias == pytest.approx(0.02)
    assert se == pytest.approx(0.0, abs=1e-12)


def test_clustered_pe_bias_standard_error_uses_showers_as_clusters() -> None:
    rows = [(1.0, 10.0, None), (-1.0, 10.0, None)]
    bias, se = c3.clustered_pe_bias(rows)
    assert bias == pytest.approx(0.0)
    # residuals +1 and -1, K/(K-1) = 2: sqrt(2 * 2) / 20
    assert se == pytest.approx(np.sqrt(4.0) / 20.0)


def test_clustered_pe_bias_needs_two_showers() -> None:
    with pytest.raises(ValueError):
        c3.clustered_pe_bias([(1.0, 10.0, None)])


@pytest.mark.parametrize(
    ("value", "se", "floor", "expected"),
    [
        (0.009, 0.0, 0.01, 0),  # inside the 1 % floor
        (0.011, 0.001, 0.01, 1),  # above the floor, 11 se
        (0.011, 0.005, 0.01, 0),  # inside 3 se (0.015)
        (0.016, 0.005, 0.01, 1),  # above 3 se
    ],
)
def test_bias_limit_is_the_larger_of_floor_and_three_standard_errors(
    value, se, floor, expected
) -> None:
    assert c3.bias_limit_exceeded(value, se, floor) == expected


def test_clustered_time_bias_is_the_mean_over_showers() -> None:
    rows = [(0.0, 1.0, 0.1), (0.0, 1.0, None), (0.0, 1.0, -0.1), (0.0, 1.0, 0.1)]
    mean, se = c3.clustered_time_bias(rows)
    assert mean == pytest.approx(1.0 / 30.0)
    assert se == pytest.approx(np.std([0.1, -0.1, 0.1], ddof=1) / np.sqrt(3))


def test_thinned_expectation_uses_the_envelope_at_the_integer_nanometre() -> None:
    # Envelope 0.5 on the 300 nm bin and 0.25 on the 301 nm bin, constant inside each bin.
    envelope = np.array(
        [[300, 0.5], [300.999999, 0.5], [301, 0.25], [301.999999, 0.25]]
    )
    weights = np.array([10.0, 4.0])
    wavelengths = np.array([300.7, 301.2])
    expected, sigma = c3.thinned_expectation(weights, wavelengths, envelope)
    assert expected == pytest.approx(10.0 * 0.5 + 4.0 * 0.25)
    assert sigma == pytest.approx(np.sqrt(6.0))


def _configuration(tmp_path: Path, *, table_value, order_table_first=True):
    qe_name = "quantum_efficiency-test.dat"
    config_dir = tmp_path / "config"
    table_dir = tmp_path / "work" / "qe-table"
    config_dir.mkdir(parents=True)
    table_dir.mkdir(parents=True)
    grid = np.arange(250, 260)
    np.savetxt(
        table_dir / qe_name, np.column_stack([grid, np.full(grid.size, table_value)])
    )
    inputs = {"qe": config_dir / qe_name, "config_dir": config_dir}
    files = {"table_dir": table_dir, "headroom": 0.7}
    first, second = f"-I{table_dir}", f"-I{config_dir}"
    argv = [first, second] if order_table_first else [second, first]
    argv = ["sim_telarray", *argv, "-C", "Bypass_Optics=3", "-"]
    return inputs, files, argv


def _texts(qe_name="quantum_efficiency-test.dat", *, atmosphere=True, bypass=True):
    lines = [f"Reading quantum efficiency from file {qe_name}."]
    if atmosphere:
        lines.append(
            f"Read atmospheric transmission data from file /x/{pc.ATMOSPHERE_FILE}"
        )
    candidate = "\n".join(lines) + "\n"
    if bypass:
        candidate += f"{c3.BYPASS3_LINE} for telescope 1\n"
    reference = "\n".join(lines) + "\n"
    return {"reference": reference, "candidate": candidate}


def test_configuration_checks_pass_when_all_evidence_is_present(tmp_path: Path) -> None:
    inputs, files, argv = _configuration(tmp_path, table_value=0.7)
    failures = c3.configuration_failures(inputs, files, argv, _texts())
    assert failures == {"atmosphere": 0, "qe_table": 0, "bypass_mode": 0}


def test_configuration_checks_flag_wrong_atmosphere_and_missing_bypass_line(
    tmp_path: Path,
) -> None:
    inputs, files, argv = _configuration(tmp_path, table_value=0.7)
    texts = _texts(atmosphere=False, bypass=False)
    failures = c3.configuration_failures(inputs, files, argv, texts)
    assert failures["atmosphere"] == 2
    assert failures["bypass_mode"] == 1


def test_configuration_checks_flag_search_order_and_unscaled_table(
    tmp_path: Path,
) -> None:
    inputs, files, argv = _configuration(
        tmp_path, table_value=0.8, order_table_first=False
    )
    failures = c3.configuration_failures(inputs, files, argv, _texts())
    assert (
        failures["qe_table"] == 2
    )  # config directory first, and the maximum is not the headroom


def test_candidate_and_reference_run_without_the_photon_threshold():
    extra = pc.telescope_arguments(Path("/cfg/model"), "/cfg/atmosphere.dat")
    assert "min_photons=0" in extra
    assert extra[extra.index("min_photons=0") - 1] == "-C"
