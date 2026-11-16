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
    assert "--telescopes" not in stage["propagation"]
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
        "candidate_envelope_trace_stream_wall_s": 2.0,
        "candidate_envelope_trace_stream_cpu_s": 4.0,
    }


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
