import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from obdeect_validation.adapters import chain_optics as co
from obdeect_validation.adapters import performance, performance_chain


def option(argv, *names):
    """Value following the first of names in an argument array."""
    for name in names:
        if name in argv:
            return argv[argv.index(name) + 1]
    raise KeyError(names)


def test_ratio_is_none_without_a_positive_reference() -> None:
    assert performance_chain.ratio(2.0, 4.0) == 0.5
    assert performance_chain.ratio(1.0, 0.0) is None
    assert performance_chain.ratio(1.0, None) is None


def test_summarise_sums_wall_and_cpu_and_takes_peak_rss() -> None:
    children = [
        {"wall_s": 1.0, "cpu_s": 0.5, "peak_rss_mb": 10.0, "argv": ["a"]},
        {"wall_s": 2.0, "cpu_s": 1.5, "peak_rss_mb": 30.0, "argv": ["b"]},
    ]
    stage = performance_chain.summarise(children)
    assert (stage["wall_s"], stage["cpu_s"], stage["peak_rss_mb"]) == (3.0, 2.0, 30.0)
    assert stage["argv"] == [["a"], ["b"]]
    with pytest.raises(RuntimeError):
        performance_chain.summarise([])


def test_total_is_sequential_sum_and_maximum_peak() -> None:
    stages = [
        {"wall_s": 1.0, "cpu_s": 1.0, "peak_rss_mb": 50.0},
        {"wall_s": 4.0, "cpu_s": 3.0, "peak_rss_mb": 20.0},
    ]
    assert performance_chain.total(stages) == {
        "wall_s": 5.0,
        "cpu_s": 4.0,
        "peak_rss_mb": 50.0,
    }


def test_measure_runs_steps_in_order_and_records_each_repetition(monkeypatch) -> None:
    monkeypatch.setattr(performance, "RUSAGE", [])
    order = []

    def step(name, wall):
        def run():
            order.append(name)
            performance.RUSAGE.append(
                {
                    "argv": [name],
                    "wall_s": wall,
                    "cpu_s": wall,
                    "peak_rss_mb": wall * 10,
                }
            )

        return name, run

    per_stage, totals = performance_chain.measure(
        [step("a", 1.0), step("b", 2.0)], reps=2
    )
    assert order == ["a", "b", "a", "b"]
    assert [len(per_stage[name]) for name in ("a", "b")] == [2, 2]
    assert totals == [
        {"wall_s": 3.0, "cpu_s": 3.0, "peak_rss_mb": 20.0},
        {"wall_s": 3.0, "cpu_s": 3.0, "peak_rss_mb": 20.0},
    ]


def test_pipeline_metrics_are_medians_and_ratios_are_candidate_over_reference() -> None:
    totals = [
        {"wall_s": 1.0, "cpu_s": 1.0, "peak_rss_mb": 10.0},
        {"wall_s": 3.0, "cpu_s": 2.0, "peak_rss_mb": 30.0},
        {"wall_s": 2.0, "cpu_s": 4.0, "peak_rss_mb": 20.0},
    ]
    per_stage = {"sim": [{"wall_s": 2.0, "cpu_s": 2.0, "peak_rss_mb": 5.0}] * 3}
    metrics = performance_chain.pipeline_metrics("reference", per_stage, totals)
    assert metrics["reference_completed"] == 1
    assert metrics["reference_reps"] == 3
    assert metrics["reference_wall_s"] == 2.0
    assert metrics["reference_peak_rss_mb"] == 20.0
    assert metrics["reference_sim_wall_s"] == 2.0

    candidate = {"wall_s": 4.0, "cpu_s": 2.0, "peak_rss_mb": 10.0}
    reference = {"wall_s": 2.0, "cpu_s": 4.0, "peak_rss_mb": 20.0}
    assert performance_chain.ratio_metrics("threads1", candidate, reference) == {
        "ratio_wall_threads1": 2.0,
        "ratio_cpu_threads1": 0.5,
        "ratio_peak_rss_threads1": 0.5,
    }


def test_failed_pipeline_is_recorded_and_gives_no_ratio(monkeypatch) -> None:
    def boom(steps, reps):
        raise RuntimeError("sim_telarray exited 1")

    monkeypatch.setattr(performance_chain, "measure", boom)
    metrics, details = {}, {"failures": {}, "stages": {}, "pipelines": {}}
    ok = performance_chain._run("candidate_threads1", [], 1, metrics, details)
    assert ok is False
    assert metrics == {"candidate_threads1_completed": 0}
    assert "sim_telarray exited 1" in details["failures"]["candidate_threads1"]["error"]
    metrics.update({f"reference_{key}": 1.0 for key in performance_chain.PHOTON_KEYS})
    performance_chain._ratios(metrics, True, ["candidate_threads1"])
    assert not any(key.startswith("ratio_") for key in metrics)


def _timed_step(calls, label, wall=1.0, cpu=1.0, rss=5.0):
    """A (name, step) pair that records its label in calls and one timed child in RUSAGE."""

    def run():
        calls.append(label)
        performance.RUSAGE.append(
            {"argv": [label], "wall_s": wall, "cpu_s": cpu, "peak_rss_mb": rss}
        )

    return "sim", run


def test_run_interleaved_alternates_pipelines_and_keeps_the_repetition_count(
    monkeypatch,
) -> None:
    monkeypatch.setattr(performance, "RUSAGE", [])
    calls = []
    metrics, details = {}, {"failures": {}, "stages": {}, "pipelines": {}}
    ok = performance_chain.run_interleaved(
        [
            ("reference", [_timed_step(calls, "reference")]),
            ("candidate_threads1", [_timed_step(calls, "candidate_threads1")]),
        ],
        reps=3,
        metrics=metrics,
        details=details,
    )
    assert ok == {"reference": True, "candidate_threads1": True}
    assert calls == ["reference", "candidate_threads1"] * 3
    assert metrics["reference_reps"] == metrics["candidate_threads1_reps"] == 3
    assert len(details["pipelines"]["reference"]) == 3
    assert len(details["stages"]["candidate_threads1"]["sim"]) == 3


def test_run_interleaved_failure_skips_only_the_failing_label(monkeypatch) -> None:
    monkeypatch.setattr(performance, "RUSAGE", [])
    calls = []

    def missing():
        calls.append("candidate_bad")
        raise OSError("missing binary")

    metrics, details = {}, {"failures": {}, "stages": {}, "pipelines": {}}
    ok = performance_chain.run_interleaved(
        [
            ("candidate_bad", [("sim", missing)]),
            ("reference", [_timed_step(calls, "reference")]),
        ],
        reps=3,
        metrics=metrics,
        details=details,
    )
    assert ok == {"candidate_bad": False, "reference": True}
    assert calls == ["candidate_bad", "reference", "reference", "reference"]
    assert metrics["candidate_bad_completed"] == 0
    assert metrics["reference_completed"] == 1
    assert metrics["reference_reps"] == 3
    assert "missing binary" in details["failures"]["candidate_bad"]["error"]
    assert "candidate_bad" not in details["pipelines"]


def test_run_single_label_records_the_same_metrics_as_pipeline_metrics(
    monkeypatch,
) -> None:
    monkeypatch.setattr(performance, "RUSAGE", [])
    calls = []
    step = _timed_step(calls, "reference", wall=2.0, cpu=1.0, rss=5.0)
    metrics, details = {}, {"failures": {}, "stages": {}, "pipelines": {}}
    assert performance_chain._run("reference", [step], 2, metrics, details) is True
    per_stage, totals = performance_chain.measure([step], 2)
    assert metrics == performance_chain.pipeline_metrics("reference", per_stage, totals)
    assert details["pipelines"]["reference"] == totals


def test_ratios_include_ratio_of_minima_over_repetitions() -> None:
    metrics = {
        "reference_wall_s": 2.0,
        "reference_cpu_s": 4.0,
        "reference_peak_rss_mb": 20.0,
        "reference_wall_s_min": 1.0,
        "reference_cpu_s_min": 3.0,
        "candidate_threads1_completed": 1,
        "candidate_threads1_wall_s": 4.0,
        "candidate_threads1_cpu_s": 2.0,
        "candidate_threads1_peak_rss_mb": 10.0,
        "candidate_threads1_wall_s_min": 2.0,
        "candidate_threads1_cpu_s_min": 1.5,
    }
    performance_chain._ratios(metrics, True, ["candidate_threads1"])
    assert metrics["ratio_wall_threads1"] == 2.0
    assert metrics["ratio_wall_threads1_min"] == 2.0
    assert metrics["ratio_cpu_threads1_min"] == 0.5
    assert "ratio_peak_rss_threads1_min" not in metrics
    assert performance_chain.ratio_min_metrics("candidate_x", {}) == {}


def test_telescope_arguments_use_one_telescope_and_optional_atmosphere() -> None:
    plain = performance_chain.telescope_arguments(Path("/cfg"))
    assert plain[:1] == ["-I/cfg"] and "maximum_telescopes=1" in plain
    assert not any("atmospheric_transmission" in arg for arg in plain)
    with_unit = performance_chain.telescope_arguments(
        Path("/cfg"), Path("/work/unit.dat")
    )
    assert "atmospheric_transmission=/work/unit.dat" in with_unit


@pytest.mark.parametrize("name", sorted(performance_chain.SOURCES))
def test_source_pipelines_use_the_same_photon_counts(name) -> None:
    spec = performance_chain.SOURCES[name]
    candidate = spec["candidate"](spec["events"], spec["photons"])
    reference = spec["reference"](spec["events"], spec["photons"])
    for names in (("--events", "-N"), ("--photons", "-n")):
        assert option(candidate, *names) == option(reference, *names)
    assert option(candidate, "--photons", "-n") == spec["photons"]
    assert option(candidate, "--events", "-N") == spec["events"]


def test_threads_variants_cover_one_thread_and_the_obdeect_default() -> None:
    assert dict(performance_chain.THREAD_VARIANTS) == {
        "threads1": 1,
        "threadsdefault": 0,
    }


def test_shower_candidate_passes_the_thread_setting_to_trace_stream(
    monkeypatch, tmp_path
) -> None:
    calls = []

    def fake_run(argv, log, cwd=None, stdout_file=None):
        calls.append(list(map(str, argv)))
        return {
            "argv": argv,
            "wall_s": 0.0,
            "cpu_s": 0.0,
            "peak_rss_mb": 0.0,
            "returncode": 0,
        }

    monkeypatch.setattr(performance_chain, "run_timed", fake_run)
    monkeypatch.setattr(performance_chain.cs, "simulate", lambda *a, **k: None)
    inputs = {
        "rotation": 0.0,
        "atmosphere": tmp_path / "atmosphere.dat",
        "models": "/models",
        "config": tmp_path / "CTAO.cfg",
        "config_dir": tmp_path,
        "unit": tmp_path / "unit.dat",
    }
    steps = performance_chain.shower_candidate(
        tmp_path, tmp_path / "in.iact", inputs, threads=1
    )
    for _, step in steps:
        step()
    trace = next(argv for argv in calls if "trace-stream" in argv)
    assert option(trace, "--threads") == "1"
    assert option(trace, "--telescopes") == "1"
    assert "--telescope-transmission" in trace
    assert [name for name, _ in steps] == [
        "from_corsika",
        "propagation",
        "compile_model",
        "trace_stream",
        "to_simtel",
        "sim_telarray",
    ]
    assert co.TELESCOPES == (1,)


def test_bin_envelope_bounds_the_qe_and_the_table_restores_it() -> None:
    qe = np.array([0.0, 0.5, 1.0])
    headroom = performance_chain.headroom_for(0.063)
    bound, table = performance_chain.bin_envelope(qe, headroom)
    assert np.all(bound >= qe)
    assert np.all(table <= 1.0)
    assert bound[0] == performance_chain.ENVELOPE_FLOOR and table[0] == 0.0
    assert np.allclose(bound * table, qe)
    assert table[1] == table[2] == pytest.approx(headroom)


def test_envelope_is_constant_on_each_nm_bin(tmp_path: Path) -> None:
    grid = np.array([250, 251])
    bound = np.array([0.1, 0.2])
    path = tmp_path / "envelope.dat"
    performance_chain.write_envelope(path, grid, bound)
    nodes = np.loadtxt(path)
    assert np.all(np.diff(nodes[:, 0]) > 0)
    # sim_telarray bins a photon at floor(lambda): any wavelength inside bin n gets bound n.
    for n, value in zip(grid, bound):
        for wavelength in (n, n + 0.25, n + 0.999):
            assert np.interp(wavelength, nodes[:, 0], nodes[:, 1]) == pytest.approx(
                value
            )


def test_sim_telarray_qe_reads_the_first_two_columns_at_integer_nm(
    tmp_path: Path,
) -> None:
    source = tmp_path / "qe.dat"
    source.write_text("250 0.1 0.01\n260 0.3 0.02\n")
    values = performance_chain.sim_telarray_qe(source, grid=np.array([250, 255, 260]))
    assert values == pytest.approx([0.1, 0.2, 0.3])


def test_envelope_files_write_the_table_under_the_config_file_name(
    tmp_path: Path,
) -> None:
    qe = tmp_path / "config" / "quantum_efficiency-test.dat"
    qe.parent.mkdir()
    qe.write_text(
        "\n".join(f"{n} {0.2 if n >= 250 else 0.0}" for n in range(200, 801, 50))
    )
    headroom = performance_chain.headroom_for(0.063)
    files = performance_chain.envelope_files(tmp_path / "work", qe, headroom)
    assert files["table_dir"] == tmp_path / "work" / "qe-table"
    table = np.loadtxt(files["table_dir"] / qe.name)
    assert table.shape == (len(performance_chain.WAVELENGTH_GRID), 2)
    assert set(np.round(table[table[:, 0] >= 250, 1], 6)) == {round(headroom, 6)}
    assert files["bound_max"] == pytest.approx(0.2 / headroom)
    assert files["headroom"] == headroom


SIM_ARGV = ["sim_telarray", "-", "-C", "Bypass_Optics=3"]


def test_envelope_pipeline_thins_with_the_qe_bound_on_one_telescope() -> None:
    stages = performance_chain.envelope_pipeline(
        "in.iact",
        "atm.dat",
        "envelope.dat",
        "lst.model",
        "camera-map.txt",
        SIM_ARGV,
        threads=0,
    )
    assert [name for name, _ in stages] == [
        "from_corsika",
        "propagation",
        "trace_stream",
        "to_simtel",
        "sim_telarray",
    ]
    stage = dict(stages)
    assert option(stage["propagation"], "--envelope") == "qe=envelope.dat"
    assert option(stage["propagation"], "--telescopes") == "1"
    assert (
        "--mode" not in stage["propagation"]
        and "--wavelength-range" not in stage["propagation"]
    )
    assert stage["propagation"][-4:] == ["-i", "-", "-o", "-"]
    assert option(stage["trace_stream"], "--threads") == 0
    assert option(stage["trace_stream"], "--telescopes") == "1"
    assert option(stage["to_simtel"], "--scaled-factor") == "qe"
    assert stage["to_simtel"][-2:] == ["-o", "-"]
    assert stage["sim_telarray"] == SIM_ARGV


def test_unthinned_pipeline_is_weighted_and_has_no_scaled_table() -> None:
    stage = dict(
        performance_chain.envelope_pipeline(
            "in.iact",
            "atm.dat",
            None,
            "lst.model",
            "camera-map.txt",
            SIM_ARGV,
            threads=0,
        )
    )
    assert "--envelope" not in stage["propagation"]
    assert option(stage["propagation"], "--mode") == "weighted"
    assert "--scaled-factor" not in stage["to_simtel"]


def test_headroom_keeps_the_first_cut_below_one_for_the_model_qe_variation() -> None:
    headroom = performance_chain.headroom_for(0.063)
    assert headroom == pytest.approx(1 / 1.378)
    assert headroom * (1 + 6 * 0.063) == pytest.approx(1.0)


def test_candidate_steps_and_camera_stage_precedence(tmp_path: Path) -> None:
    files = {"envelope": tmp_path / "e.dat", "table_dir": tmp_path / "qe-table"}
    inputs = {
        "atmosphere": tmp_path / "atm.dat",
        "reference_atmosphere": tmp_path / "OBS.dat",
        "rotation": 0.0,
        "config": tmp_path / "CTAO.cfg",
        "config_dir": tmp_path,
        "models": "/models",
    }
    steps = performance_chain.envelope_candidate(
        tmp_path, tmp_path / "in.iact", inputs, files, threads=0
    )
    assert [name for name, _ in steps] == ["compile_model", "pipeline"]
    argv, output = performance_chain.simtel_candidate(tmp_path, inputs, files)
    assert argv[-1] == "-"  # standard input from to-simtel
    table = argv.index(f"-I{files['table_dir']}")
    assert table < argv.index(f"-I{tmp_path}")  # the table shadows the model's QE file
    assert "Bypass_Optics=3" in argv
    assert output.name == "candidate.simtel"


def test_run_pipeline_pipes_stages_and_records_each_member(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(performance, "RUSAGE", [])
    output = tmp_path / "out.bin"
    logs = {"producer": tmp_path / "p.log", "consumer": tmp_path / "c.log"}
    record = performance_chain.run_pipeline(
        [("producer", ["printf", "abc"]), ("consumer", ["cat"])], output, logs
    )
    assert output.read_bytes() == b"abc"
    assert [child["name"] for child in record["children"]] == ["producer", "consumer"]
    assert (
        record["wall_s"] >= max(child["wall_s"] for child in record["children"]) - 1e-9
    )
    assert performance.RUSAGE[-1] is record


def test_run_pipeline_raises_on_a_failed_member(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(performance, "RUSAGE", [])
    logs = {"bad": tmp_path / "b.log", "sink": tmp_path / "s.log"}
    with pytest.raises(RuntimeError, match="bad exited 1"):
        performance_chain.run_pipeline(
            [("bad", ["false"]), ("sink", ["cat"])], tmp_path / "out.bin", logs
        )


def test_run_pipeline_kills_started_members_when_setup_fails(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(performance, "RUSAGE", [])
    started = []
    real_popen = subprocess.Popen

    def popen(argv, **kwargs):
        if argv[0] == "no-such-stage":
            raise FileNotFoundError(argv[0])
        proc = real_popen(argv, **kwargs)
        started.append(proc)
        return proc

    monkeypatch.setattr(performance_chain.subprocess, "Popen", popen)
    logs = {"sleeper": tmp_path / "s.log", "missing": tmp_path / "m.log"}
    stages = [
        ("sleeper", [sys.executable, "-c", "import time; time.sleep(60)"]),
        ("missing", ["no-such-stage"]),
    ]
    with pytest.raises(FileNotFoundError):
        performance_chain.run_pipeline(stages, tmp_path / "out.bin", logs)
    assert started[0].returncode < 0  # killed and reaped, not left running
    assert performance.RUSAGE == []


def test_summarise_keeps_pipeline_members() -> None:
    members = [{"name": "propagation", "wall_s": 1.0, "cpu_s": 0.5, "peak_rss_mb": 3.0}]
    stage = performance_chain.summarise(
        [
            {
                "wall_s": 2.0,
                "cpu_s": 1.0,
                "peak_rss_mb": 3.0,
                "argv": [],
                "children": members,
            }
        ]
    )
    assert stage["children"] == members


def test_child_metrics_are_medians_per_member() -> None:
    summaries = [
        {"children": [{"name": "trace_stream", "wall_s": w, "cpu_s": 2 * w}]}
        for w in (1.0, 3.0, 2.0)
    ]
    assert performance_chain.child_metrics("candidate_envelope", summaries) == {
        "candidate_envelope_trace_stream_exit_s": 2.0,
        "candidate_envelope_trace_stream_cpu_s": 4.0,
    }


def test_aggregate_gives_range_and_relative_spread() -> None:
    records = [
        {"wall_s": w, "cpu_s": 1.0, "peak_rss_mb": 1.0} for w in (9.0, 10.0, 11.0)
    ]
    metrics = performance_chain._aggregate("reference", records)
    assert metrics["reference_wall_s_min"] == 9.0
    assert metrics["reference_wall_s_max"] == 11.0
    assert metrics["reference_cpu_s_min"] == metrics["reference_cpu_s_max"] == 1.0
    assert metrics["reference_wall_s_spread"] == pytest.approx(2.0 / 10.0)
    zero = [{"wall_s": 0.0, "cpu_s": 0.0, "peak_rss_mb": 0.0}] * 2
    assert "x_wall_s_spread" not in performance_chain._aggregate("x", zero)


def test_timing_reliability_flags_spread_and_load() -> None:
    quiet = {"reference_wall_s_spread": 0.01, "host_loadavg_before": 0.4}
    assert performance_chain.timing_reliability(quiet, ["reference"])[0] == 0
    noisy = {
        "reference_wall_s_spread": 0.08,
        "candidate_threads1_wall_s_spread": 0.01,
        "host_loadavg_before": 1.5,
    }
    flag, details = performance_chain.timing_reliability(
        noisy, ["reference", "candidate_threads1"]
    )
    assert flag == 1
    assert len(details["reasons"]) == 2
    assert "reference wall spread" in details["reasons"][0]
    assert "host_loadavg_before" in details["reasons"][1]
    assert performance_chain.timing_reliability({}, ["reference"]) == (
        0,
        {
            "spread_limit": performance_chain.SPREAD_LIMIT,
            "load_limit": performance_chain.LOAD_LIMIT,
            "reasons": [],
        },
    )


def test_loadavg_is_none_where_the_host_does_not_report_it(monkeypatch) -> None:
    def unavailable():
        raise OSError("no load average")

    monkeypatch.setattr(performance_chain.os, "getloadavg", unavailable)
    assert performance_chain.loadavg() is None
    monkeypatch.delattr(performance_chain.os, "getloadavg")
    assert performance_chain.loadavg() is None


def test_normalised_metrics_per_shower_per_photon_and_cpu_shares() -> None:
    metrics = {
        "reference_wall_s": 100.0,
        "reference_cpu_s": 50.0,
        "candidate_threads1_wall_s": 40.0,
        "candidate_threads1_cpu_s": 20.0,
        "candidate_threads1_compile_model_cpu_s": 4.0,
        "candidate_threads1_propagation_cpu_s": 4.0,
        "candidate_threads1_trace_stream_cpu_s": 8.0,
        "candidate_threads1_to_simtel_cpu_s": 2.0,
        "candidate_threads1_sim_telarray_cpu_s": 2.0,
        "candidate_threads1_from_corsika_cpu_s": 0.0,
        "candidate_envelope_threads1_wall_s": 10.0,
        "candidate_envelope_threads1_cpu_s": 12.0,
        "candidate_envelope_threads1_pipeline_cpu_s": 9.0,
        "candidate_envelope_threads1_compile_model_cpu_s": 3.0,
        "candidate_envelope_threads1_propagation_cpu_s": 3.0,
        "candidate_envelope_threads1_sim_telarray_cpu_s": 6.0,
        "telescope_photons": 1e6,
    }
    details = {
        "settings": {"showers": 10},
        "pipelines": {"reference": [], "candidate_threads1": []},
        "stages": {
            "candidate_threads1": {
                "compile_model": [],
                "propagation": [],
                "trace_stream": [],
                "to_simtel": [],
                "sim_telarray": [],
                "from_corsika": [],
            },
            "candidate_envelope_threads1": {
                "compile_model": [],
                "pipeline": [
                    {
                        "children": [
                            {"name": "propagation"},
                            {"name": "sim_telarray"},
                        ]
                    }
                ],
            },
        },
    }
    out = performance_chain.normalised_metrics(metrics, details)
    assert out["reference_s_per_shower"] == 10.0
    assert out["candidate_threads1_cpu_s_per_shower"] == 2.0
    assert out["reference_cpu_ns_per_photon"] == pytest.approx(50.0 * 1e9 / 1e6)
    # candidate_threads1: denominator is its CPU without compile_model (16 s).
    assert out["candidate_threads1_propagation_cpu_share"] == pytest.approx(4 / 16)
    assert out["candidate_threads1_trace_stream_cpu_share"] == pytest.approx(0.5)
    assert "candidate_threads1_compile_model_cpu_share" not in out
    # envelope: pipe members over the pipe stage (9 s).
    assert out["candidate_envelope_threads1_propagation_cpu_share"] == pytest.approx(
        3 / 9
    )
    assert out["candidate_envelope_threads1_sim_telarray_cpu_share"] == pytest.approx(
        6 / 9
    )
    assert "candidate_envelope_threads1_compile_model_cpu_share" not in out


def test_normalised_metrics_leave_out_what_is_not_measured() -> None:
    details = {
        "settings": {},
        "pipelines": {"reference": []},
        "stages": {},
    }
    metrics = {"reference_wall_s": 5.0, "reference_cpu_s": 3.0}
    out = performance_chain.normalised_metrics(metrics, details)
    assert out == {}


def test_reference_optics_is_the_difference_of_the_two_sim_telarray_runs() -> None:
    metrics = {
        "reference_sim_telarray_cpu_s": 10.0,
        "reference_sim_telarray_wall_s": 12.0,
        "candidate_threads1_sim_telarray_cpu_s": 2.0,
        "candidate_threads1_sim_telarray_wall_s": 3.0,
        "candidate_threads1_trace_stream_cpu_s": 4.0,
    }
    details = {}
    out = performance_chain.reference_optics_metrics(metrics, details)
    assert out == {
        "reference_optics_cpu_s": 8.0,
        "reference_optics_wall_s": 9.0,
        "reference_camera_cpu_s": 2.0,
        "trace_stream_vs_reference_optics_cpu_ratio": 0.5,
    }
    assert details == {}


def test_reference_optics_ratio_is_omitted_and_explained_when_not_positive() -> None:
    metrics = {
        "reference_sim_telarray_cpu_s": 2.0,
        "candidate_threads1_sim_telarray_cpu_s": 2.5,
        "candidate_threads1_trace_stream_cpu_s": 4.0,
    }
    details = {}
    out = performance_chain.reference_optics_metrics(metrics, details)
    assert out["reference_optics_cpu_s"] == pytest.approx(-0.5)
    assert "trace_stream_vs_reference_optics_cpu_ratio" not in out
    assert details["reference_optics"]["trace_stream_ratio"].startswith("omitted")


def test_gain_metrics_thinning_and_piping_without_compile_model() -> None:
    metrics = {
        "candidate_unthinned_threads1_cpu_s": 30.0,
        "candidate_unthinned_threads1_wall_s": 15.0,
        "candidate_unthinned_threads1_compile_model_cpu_s": 2.0,
        "candidate_unthinned_threads1_compile_model_wall_s": 1.0,
        "candidate_envelope_threads1_cpu_s": 10.0,
        "candidate_envelope_threads1_wall_s": 5.0,
        "candidate_envelope_threads1_compile_model_cpu_s": 2.0,
        "candidate_envelope_threads1_compile_model_wall_s": 1.0,
        "candidate_threads1_cpu_s": 50.0,
        "candidate_threads1_wall_s": 40.0,
        "candidate_threads1_compile_model_cpu_s": 4.0,
        "candidate_threads1_compile_model_wall_s": 2.0,
    }
    out = performance_chain.gain_metrics(metrics)
    assert out["thinning_gain_cpu"] == pytest.approx(28.0 / 8.0)
    assert out["thinning_gain_wall"] == pytest.approx(14.0 / 4.0)
    assert out["piping_gain_cpu"] == pytest.approx(46.0 / 28.0)
    assert out["piping_gain_wall"] == pytest.approx(38.0 / 14.0)
    assert performance_chain.gain_metrics({}) == {}


def test_unthinned_candidate_is_one_pipe_without_envelope_or_scaled_table(
    tmp_path: Path, monkeypatch
) -> None:
    captured = []
    monkeypatch.setattr(
        performance_chain,
        "_checked_pipeline",
        lambda stages, output, logs: captured.append(dict(stages)),
    )
    inputs = {
        "rotation": 0.0,
        "atmosphere": tmp_path / "atmosphere.dat",
        "models": "/models",
        "config": tmp_path / "CTAO.cfg",
        "config_dir": tmp_path,
        "reference_atmosphere": tmp_path / "reference.dat",
    }
    steps = performance_chain.envelope_candidate(
        tmp_path, tmp_path / "in.iact", inputs, None, threads=1
    )
    assert [name for name, _ in steps] == ["compile_model", "pipeline"]
    steps[1][1]()
    stage = captured[0]
    assert "--envelope" not in stage["propagation"]
    assert option(stage["propagation"], "--mode") == "weighted"
    assert option(stage["trace_stream"], "--threads") == 1
    assert "--scaled-factor" not in stage["to_simtel"]
    assert not any("qe-table" in str(arg) for arg in stage["sim_telarray"])


def test_variant_counts_parse_the_trace_summary(tmp_path: Path, monkeypatch) -> None:
    (tmp_path / "trace_stream.log").write_text(
        "trace-stream: 20 events, 1399647 photons in, 780587 out\n"
    )
    monkeypatch.setattr(performance_chain, "iact_weight", lambda path: 17.0)
    assert performance_chain.variant_counts(tmp_path, tmp_path / "c.iact") == {
        "into_simtel_weight": 17.0,
        "into_trace_photons": 1399647,
        "out_of_trace_photons": 780587,
    }
