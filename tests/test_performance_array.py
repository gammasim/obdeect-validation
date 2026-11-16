"""Pure helpers of the multi-telescope air-shower timing (no sim_telarray, no obdeect run)."""

import pytest

from obdeect_validation.adapters import performance_array as pa


def _stages(walls, cpus, compile_wall=5.0):
    """Stage records of two repetitions: one compile and one pipe per telescope."""
    stages = {}
    for k, (wall, cpu) in enumerate(zip(walls, cpus), start=1):
        stages[f"compile_tel{k}"] = [
            {"wall_s": compile_wall, "cpu_s": compile_wall}
        ] * 2
        stages[f"pipe_tel{k}"] = [
            {"wall_s": wall, "cpu_s": cpu},
            {"wall_s": wall * 2, "cpu_s": cpu * 2},
        ]
    return stages


def test_chain_medians_sum_pipes_per_repetition_and_skip_compiles():
    stages = _stages(walls=[1.0, 2.0], cpus=[0.5, 1.5])
    # repetition 0: 1 + 2 = 3; repetition 1: 2 + 4 = 6; median of two values is 4.5
    assert pa.chain_medians(stages, "wall_s") == pytest.approx(4.5)
    assert pa.chain_medians(stages, "cpu_s") == pytest.approx((2.0 + 4.0) / 2)


def _metrics(chain_wall, chain_cpu, bias_flags):
    metrics = {
        "reference_wall_s": 10.0,
        "reference_cpu_s": 9.0,
        "photon_checks_failed": 0,
    }
    for t, flag in zip(range(1, 7), bias_flags):
        metrics[f"physics_tel{t}_photoelectron_bias_limit_exceeded"] = flag
    return metrics, {"wall_s": chain_wall, "cpu_s": chain_cpu}


def test_budgets_pass_when_chain_within_reference_and_all_bias_within_limits():
    metrics, chain = _metrics(9.0, 8.0, [0] * 6)
    out, reasons = pa.array_budgets(metrics, chain)
    assert out["array_chain_wall_budget_exceeded"] == 0
    assert out["array_chain_cpu_budget_exceeded"] == 0
    assert out["equal_statistics_budget_exceeded"] == 0
    assert out["array_chain_wall_ratio"] == pytest.approx(0.9)
    assert reasons == {}


def test_budget_exceeded_when_chain_slower_or_one_telescope_biased():
    metrics, chain = _metrics(11.0, 8.0, [0, 0, 0, 0, 1, 0])
    out, _ = pa.array_budgets(metrics, chain)
    assert out["array_chain_wall_budget_exceeded"] == 1
    assert out["array_chain_cpu_budget_exceeded"] == 0
    assert out["equal_statistics_budget_exceeded"] == 1


def test_unmeasured_budget_is_exceeded_with_reason():
    metrics, _ = _metrics(9.0, 8.0, [0] * 6)
    out, reasons = pa.array_budgets(metrics, {})
    assert out["array_chain_wall_budget_exceeded"] == 1
    assert "wall" in reasons and "cpu" in reasons


def test_camera_map_names_one_telescope(tmp_path):
    path = pa.camera_map(tmp_path, 5, 10.893)
    lines = [line for line in path.read_text().splitlines() if not line.startswith("#")]
    assert lines == [f"5 {pa.MIRROR_CLASS} 10.893"]


STDOUT = """Photon statistics:
  All photons:                    0         0         0         0      2089         0
  Used photons:                   0         0         0         0         0         0
  Not absorbed/max. Q.E.:         0         0         0         0         0         0
"""


def test_photon_statistics_reads_the_column_of_the_telescope():
    assert pa.photon_statistics(STDOUT, 5) == (2089, 0)
    assert pa.photon_statistics(STDOUT, 1) == (0, 0)


def test_photon_statistics_rejects_missing_or_short_lines():
    with pytest.raises(ValueError, match="found 0"):
        pa.photon_statistics("Used photons: 1 2\n", 1)
    with pytest.raises(ValueError, match="no column for telescope 6"):
        pa.photon_statistics("All photons: 1 2\nUsed photons: 1 2\n", 6)


def test_trace_out_rows_is_the_output_of_the_summary():
    text = "trace-stream: 2 events, 6604 photons in, 2089 out\n"
    assert pa.trace_out_rows(text) == 2089
    with pytest.raises(ValueError, match="no trace-stream summary"):
        pa.trace_out_rows("nothing\n")


def _photon_outputs(tmp_path, all_values, trace_rows):
    outputs = {}
    for telescope, _ in pa.ARRAY:
        work = tmp_path / f"candidate-tel{telescope}"
        work.mkdir()
        all_photons = all_values[telescope]
        (work / "candidate.stdout").write_text(
            "  All photons: "
            + " ".join(str(all_photons if t == telescope else 0) for t in range(1, 7))
            + "\n"
            "  Used photons: "
            + " ".join(str(all_photons if t == telescope else 0) for t in range(1, 7))
            + "\n"
        )
        (work / "trace_stream.log").write_text(
            f"trace-stream: 2 events, 9 photons in, {trace_rows[telescope]} out\n"
        )
        outputs[telescope] = (work, None)
    return outputs


def test_photon_checks_pass_when_nothing_is_dropped_and_bunches_match_trace(tmp_path):
    counts = {t: 100 * t for t, _ in pa.ARRAY}
    metrics, details = {}, {}
    pa.photon_checks(_photon_outputs(tmp_path, counts, counts), metrics, details)
    assert metrics["photon_checks_failed"] == 0
    assert metrics["photon_used_ne_all_tel5"] == 0
    assert metrics["to_simtel_bunch_ne_trace_rows_tel5"] == 0
    assert metrics["to_simtel_bunches_tel5"] == 500
    assert details["photon_checks"]["tel5"] == {
        "all": 500,
        "used": 500,
        "trace_rows": 500,
    }


def test_photon_checks_flag_a_gated_or_mismatched_telescope(tmp_path):
    counts = {t: 100 * t for t, _ in pa.ARRAY}
    rows = dict(counts)
    rows[6] = 3  # trace rows differ from the bunches sim_telarray read
    metrics, _ = {}, {}
    pa.photon_checks(_photon_outputs(tmp_path, counts, rows), metrics, {})
    assert metrics["to_simtel_bunch_ne_trace_rows_tel6"] == 1
    assert metrics["photon_checks_failed"] == 1


def test_sim_options_disable_the_photon_threshold_on_both_sides():
    reference = pa.reference_options("/cfg/atmosphere.dat")
    candidate = pa.candidate_options("/cfg/atmosphere.dat")
    assert reference["min_photons"] == 0 == candidate["min_photons"]
    assert reference["Bypass_Optics"] == 0 and candidate["Bypass_Optics"] == 3


def test_measured_flags_fewer_than_twenty_matched_samples():
    details = {"failures": {}}
    few = pa.measured(
        {"matched_event_telescope_samples": 19, "photoelectron_relative_bias": 0.0},
        "tel5",
        details,
    )
    assert few == {
        "matched_event_telescope_samples": 19,
        "photoelectron_bias_limit_exceeded": 1,
    }
    assert (
        "not measured: 19 matched samples < 20" in details["failures"]["physics_tel5"]
    )
    enough = {"matched_event_telescope_samples": 20, "photoelectron_relative_bias": 0.0}
    assert pa.measured(enough, "tel5", details) is enough
