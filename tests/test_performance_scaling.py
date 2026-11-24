"""Pure parts of the scaling measurement: fits, fixed share, thread metrics, metric assembly, chain flow."""

import pytest

from obdeect_validation.adapters import performance_chain
from obdeect_validation.adapters import performance_scaling as ps

INPUTS = {"rotation": 0.0, "atmosphere": "atm.dat", "models": "/models"}


def record(wall, cpu=None):
    return {"wall_s": wall, "cpu_s": wall if cpu is None else cpu}


def summary(**stages):
    """Point summary from {stage: (wall_s, cpu_s)}."""
    return {
        stage: {"wall_s": wall, "cpu_s": cpu, "reps": 1}
        for stage, (wall, cpu) in stages.items()
    }


def test_fit_line_recovers_three_point_line():
    fixed, per_unit = ps.fit_line([0, 10, 20], [2.0, 7.0, 12.0])
    assert fixed == pytest.approx(2.0)
    assert per_unit == pytest.approx(0.5)


def test_fit_line_two_points_is_exact():
    fixed, per_unit = ps.fit_line([10, 20], [7.0, 12.0])
    assert fixed == pytest.approx(2.0)
    assert per_unit == pytest.approx(0.5)


def test_fit_line_needs_two_distinct_counts():
    with pytest.raises(ValueError):
        ps.fit_line([10, 10], [1.0, 2.0])


def test_fit_stage_fixed_only_uses_mean_and_zero_slope():
    fit = ps.fit_stage([(0, 1.0, 0.5), (10, 3.0, 1.5), (20, 2.0, 1.0)], fixed_only=True)
    assert fit == {
        "fixed_wall_s": pytest.approx(2.0),
        "per_shower_wall_s": 0.0,
        "fixed_cpu_s": pytest.approx(1.0),
        "per_shower_cpu_s": 0.0,
    }


def test_fixed_share_and_zero_total():
    assert ps.fixed_share(2.0, 0.5, 4) == pytest.approx(0.5)
    assert ps.fixed_share(0.0, 0.0, 20) is None


def test_stage_metrics_names_and_share_at_n():
    fit = {
        "fixed_wall_s": 2.0,
        "per_shower_wall_s": 0.5,
        "fixed_cpu_s": 1.0,
        "per_shower_cpu_s": 0.25,
    }
    metrics = ps.stage_metrics("propagation", fit, 20)
    assert metrics["propagation_fixed_wall_s"] == 2.0
    assert metrics["propagation_per_shower_wall_s"] == 0.5
    assert metrics["propagation_fixed_cpu_s"] == 1.0
    assert metrics["propagation_per_shower_cpu_s"] == 0.25
    assert metrics["propagation_fixed_share_at_N"] == pytest.approx(
        1 / 6
    )  # 2 / (2 + 0.5 x 20)


def test_thread_counts_capped_at_nproc():
    assert ps.thread_counts(1) == (1,)
    assert ps.thread_counts(0) == (1,)
    assert ps.thread_counts(3) == (1, 2)
    assert ps.thread_counts(4) == (1, 2, 4)
    assert ps.thread_counts(64) == (1, 2, 4, 8)


def test_thread_metrics_speedup_and_efficiency():
    metrics = ps.thread_metrics({1: (8.0, 8.0), 2: (4.0, 8.1), 4: (2.5, 9.0)}, 4)
    assert metrics["trace_stream_threads1_speedup"] == pytest.approx(1.0)
    assert metrics["trace_stream_threads2_speedup"] == pytest.approx(2.0)
    assert metrics["trace_stream_threads2_efficiency"] == pytest.approx(1.0)
    assert metrics["trace_stream_threads4_speedup"] == pytest.approx(3.2)
    assert metrics["trace_stream_threads4_efficiency"] == pytest.approx(0.8)
    assert metrics["trace_stream_threads4_cpu_s"] == 9.0
    assert metrics["host_nproc"] == 4


def test_median_of_records():
    assert ps.median_of([record(1.0, 2.0), record(3.0, 4.0), record(2.0, 3.0)]) == (
        2.0,
        3.0,
    )


def full_points(events=10, zero=True):
    """Points with a known line: wall = 1 + 0.1 x showers, cpu = 2 + 0.2 x showers."""

    def point(showers):
        return (
            showers,
            summary(
                **{
                    stage: (1.0 + 0.1 * showers, 2.0 + 0.2 * showers)
                    for stage in (
                        "from_corsika",
                        "propagation",
                        "trace_stream",
                        "to_simtel",
                    )
                },
                compile_model=(0.5, 0.4),
            ),
        )

    points = {events: point(events), 2 * events: point(2 * events)}
    if zero:
        points[0] = point(0)
    return points


def test_assemble_fits_all_stages_and_reports_fixed_compile():
    metrics, fits = ps.assemble(full_points(), 10, {1: (8.0, 8.0)}, 1)
    assert metrics["propagation_fixed_wall_s"] == pytest.approx(1.0)
    assert metrics["propagation_per_shower_wall_s"] == pytest.approx(
        0.1
    )  # 1 + 0.1 x showers
    assert metrics["propagation_fixed_cpu_s"] == pytest.approx(2.0)
    assert metrics["propagation_per_shower_cpu_s"] == pytest.approx(0.2)
    assert metrics["propagation_fixed_share_at_N"] == pytest.approx(
        0.5
    )  # 1.0 / (1.0 + 0.1 x 10)
    assert metrics["compile_model_per_shower_wall_s"] == 0.0
    assert metrics["compile_model_fixed_share_at_N"] == pytest.approx(1.0)
    assert metrics["compile_model_fixed_wall_s"] == pytest.approx(0.5)
    for stage in (
        "from_corsika",
        "propagation",
        "compile_model",
        "trace_stream",
        "to_simtel",
    ):
        assert metrics[f"{stage}_scaling_completed"] == 1
    assert len(fits["propagation"]["points"]) == 3
    assert metrics["thread_scaling_completed"] == 1  # nproc 1: only threads=1


def test_assemble_two_points_when_zero_shower_point_failed():
    points = full_points(zero=False)
    metrics, fits = ps.assemble(points, 10, {}, 1)
    assert metrics["propagation_scaling_completed"] == 1
    assert len(fits["propagation"]["points"]) == 2
    assert metrics["propagation_fixed_wall_s"] == pytest.approx(1.0)


def test_assemble_stage_not_completed_without_n_or_2n():
    points = {10: (10, summary(from_corsika=(1.0, 1.0)))}  # only N: 2N missing
    metrics, fits = ps.assemble(points, 10, {}, 1)
    assert metrics["from_corsika_scaling_completed"] == 0
    assert metrics["trace_stream_scaling_completed"] == 0
    assert "trace_stream_fixed_wall_s" not in metrics
    assert "trace_stream" not in fits


def test_assemble_threads_completed_only_with_all_counts():
    rows = {k: (8.0 / k, 8.0) for k in (1, 2, 4)}
    metrics, _ = ps.assemble(full_points(), 10, rows, 4)
    assert metrics["thread_scaling_completed"] == 1
    assert metrics["trace_stream_threads4_speedup"] == pytest.approx(4.0)
    assert metrics["host_nproc"] == 4

    missing, _ = ps.assemble(full_points(), 10, {1: (8.0, 8.0)}, 4)
    assert missing["thread_scaling_completed"] == 0
    assert "trace_stream_threads2_wall_s" not in missing


def test_chain_commands_order_and_native_outputs(tmp_path):
    commands = ps.chain_commands(tmp_path, tmp_path / "in.iact", INPUTS)
    assert [stage for stage, *_ in commands] == list(ps.STAGES)
    stage, argv, _, stdout = commands[0]
    assert (
        argv[:3] == ["obdeect", "stream", "from-corsika"]
        and stdout == tmp_path / "showers.gst"
    )
    trace = dict((s, a) for s, a, *_ in commands)["trace_stream"]
    assert trace[trace.index("--threads") + 1] == 1


def test_trace_argv_threads(tmp_path):
    argv = ps.trace_argv("m", "p.gst", "t.gst", 8)
    assert argv[argv.index("--threads") + 1] == 8
    assert argv[argv.index("-i") + 1] == "p.gst"
    assert argv[argv.index("-o") + 1] == "t.gst"


def test_run_chain_once_stops_at_failure_and_reports_skipped(tmp_path, monkeypatch):
    calls = []

    def fake_run(argv, log, stdout_file=None):
        calls.append(argv[:3])
        if argv[1] == "propagation":
            raise RuntimeError("propagation exited 1")
        return record(0.1)

    monkeypatch.setattr(performance_chain, "run_timed", fake_run)
    records = {stage: [] for stage in ps.STAGES}
    failures = {}
    skipped = ps.run_chain_once(
        tmp_path, tmp_path / "in.iact", INPUTS, records, failures
    )
    assert list(failures) == ["propagation"]
    assert skipped == ["compile_model", "trace_stream", "to_simtel"]
    assert len(records["from_corsika"]) == 1
    assert records["compile_model"] == []
    assert calls == [
        ["obdeect", "stream", "from-corsika"],
        ["obdeect", "propagation", "run"],
    ]


def test_measure_points_breaks_repetitions_after_failure(tmp_path, monkeypatch):
    def fake_run(argv, log, stdout_file=None):
        if argv[1] == "ray-tracing" and argv[2] == "trace-stream":
            raise RuntimeError("trace-stream failed")
        return record(1.0, 2.0)

    monkeypatch.setattr(performance_chain, "run_timed", fake_run)
    state = ps.measure_points({"only": (tmp_path, tmp_path / "in.iact")}, INPUTS, 3)
    records, failures, skipped = state["only"]
    assert set(failures) == {"trace_stream"}
    assert skipped == ["to_simtel"]
    assert len(records["propagation"]) == 1  # no repetition after the failure
    summary_ = ps.point_summary(records, failures)
    assert "trace_stream" not in summary_
    assert summary_["propagation"] == {
        "wall_s": 1.0,
        "cpu_s": 2.0,
        "wall_s_min": 1.0,
        "cpu_s_min": 2.0,
        "reps": 1,
    }


def test_measure_points_alternates_repetitions(tmp_path, monkeypatch):
    passes = []

    def fake_run(argv, log, stdout_file=None):
        if log.name == "from-corsika.log":  # first stage: one entry per chain pass
            passes.append(log.parent.name)
        return record(1.0)

    monkeypatch.setattr(performance_chain, "run_timed", fake_run)
    runs = {
        name: (tmp_path / name, tmp_path / name / "in.iact") for name in ("a", "b", "c")
    }
    for work, _ in runs.values():
        work.mkdir()
    ps.measure_points(runs, INPUTS, 2)
    assert passes == ["a", "b", "c", "a", "b", "c"]


def test_measure_points_failing_point_does_not_stop_others(tmp_path, monkeypatch):
    passes = []

    def fake_run(argv, log, stdout_file=None):
        if log.name == "from-corsika.log":
            passes.append(log.parent.name)
        if log.parent.name == "bad" and argv[1] == "propagation":
            raise RuntimeError("propagation exited 1")
        return record(1.0)

    monkeypatch.setattr(performance_chain, "run_timed", fake_run)
    runs = {
        name: (tmp_path / name, tmp_path / name / "in.iact") for name in ("bad", "good")
    }
    for work, _ in runs.values():
        work.mkdir()
    state = ps.measure_points(runs, INPUTS, 3)
    assert passes == ["bad", "good", "good", "good"]  # bad stops after its first pass
    bad_records, bad_failures, bad_skipped = state["bad"]
    assert list(bad_failures) == ["propagation"]
    assert bad_skipped == ["compile_model", "trace_stream", "to_simtel"]
    assert len(bad_records["from_corsika"]) == 1
    good_records, good_failures, good_skipped = state["good"]
    assert good_failures == {} and good_skipped == []
    assert len(good_records["to_simtel"]) == 3


def test_measure_threads_alternates_and_keeps_completed_counts(tmp_path, monkeypatch):
    order = []

    def fake_run(argv, log, stdout_file=None):
        threads = argv[argv.index("--threads") + 1]
        order.append(threads)
        if threads == 2 and len(order) > 3:  # threads=2 fails on its second repetition
            raise RuntimeError("trace-stream exited 1")
        return record(8.0 / threads)

    monkeypatch.setattr(performance_chain, "run_timed", fake_run)
    rows, failures = ps.measure_threads(tmp_path, 4, 3)
    assert order == [1, 2, 4, 1, 2, 4, 1, 4]
    assert set(failures) == {2}
    assert set(rows) == {1, 4}  # counts that did not fail keep all their repetitions
    assert rows[4] == (2.0, 2.0)
    assert ps.assemble(full_points(), 10, rows, 4)[0]["thread_scaling_completed"] == 0
