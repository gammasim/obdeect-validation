import sys
from pathlib import Path

from obdeect_validation.adapters import performance, performance_profile

REPORT = """# Samples: 1K of event 'cycles:P'
# Event count (approx.): 1234567
#
# Overhead  Command  Shared Object      Symbol
# ........  .......  .................  ..............................
#
    41.25%  obdeect  library.so          [.] trace_photon
            |
            |--40.10%--trace_photon_loop
            |          main
            ---trace_photon_loop
    12.50%  obdeect  libc.so.6          [.] __memmove_avx_unaligned
     3.00%  obdeect  [kernel.kallsyms]  [k] clear_page_erms
"""

SCRIPT = "obdeect 1234 1.0: 1 cycles:\n\tfeed0 trace_photon\n\n" * 3


def test_parse_report_keeps_top_level_symbols_only() -> None:
    symbols = performance_profile.parse_report(REPORT)
    assert [s["symbol"] for s in symbols] == [
        "trace_photon",
        "__memmove_avx_unaligned",
        "clear_page_erms",
    ]
    assert symbols[0] == {
        "overhead_percent": 41.25,
        "command": "obdeect",
        "shared_object": "library.so",
        "symbol": "trace_photon",
    }
    assert symbols[2]["shared_object"] == "[kernel.kallsyms]"


def test_count_samples_counts_unindented_header_lines() -> None:
    assert performance_profile.count_samples(SCRIPT) == 3
    assert performance_profile.count_samples("") == 0


def test_stage_argv_trace_stream_is_single_thread_file_mode(tmp_path) -> None:
    inputs = {"atmosphere": tmp_path / "atm.dat"}
    argv, stdout_file = performance_profile.stage_argv(
        "trace_stream", tmp_path, tmp_path / "in.iact", inputs, {}
    )
    argv = list(map(str, argv))
    assert argv[:3] == ["obdeect", "ray-tracing", "trace-stream"]
    assert argv[argv.index("--threads") + 1] == "1"
    assert argv[argv.index("-i") + 1] == str(tmp_path / "propagated.gst")
    assert argv[argv.index("-o") + 1] == str(tmp_path / "traced.gst")
    assert stdout_file is None


def test_stage_argv_envelope_propagation_uses_the_bound_table(tmp_path) -> None:
    inputs = {"atmosphere": tmp_path / "atm.dat"}
    envelope = {"envelope": tmp_path / "qe-envelope.dat"}
    argv, _ = performance_profile.stage_argv(
        "propagation_envelope", tmp_path, tmp_path / "in.iact", inputs, envelope
    )
    argv = list(map(str, argv))
    assert argv[argv.index("--envelope") + 1] == f"qe={tmp_path / 'qe-envelope.dat'}"
    assert "--mode" not in argv
    assert argv[argv.index("-o") + 1] == str(tmp_path / "propagated-envelope.gst")


def test_stage_argv_from_corsika_writes_stdout_to_native_stream(tmp_path) -> None:
    argv, stdout_file = performance_profile.stage_argv(
        "from_corsika", tmp_path, tmp_path / "in.iact", {}, {}
    )
    argv = list(map(str, argv))
    assert argv == ["obdeect", "stream", "from-corsika", str(tmp_path / "in.iact")]
    assert stdout_file == tmp_path / "showers.gst"


def test_stage_argv_rejects_unknown_stage(tmp_path) -> None:
    try:
        performance_profile.stage_argv("nope", tmp_path, tmp_path, {}, {})
    except ValueError:
        return
    raise AssertionError("unknown stage accepted")


def test_execute_sends_stdout_to_file_and_returns_exit_code(tmp_path) -> None:
    code = performance_profile._execute(
        [sys.executable, "-c", "print('payload')"],
        tmp_path / "log.txt",
        tmp_path / "out.txt",
    )
    assert code == 0
    assert (tmp_path / "out.txt").read_text().strip() == "payload"
    failing = performance_profile._execute(
        [sys.executable, "-c", "import sys; sys.exit(3)"], tmp_path / "log2.txt"
    )
    assert failing == 3


def test_flamegraph_reports_missing_tools(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(performance_profile.shutil, "which", lambda name: None)
    reason = performance_profile.flamegraph(SCRIPT, tmp_path / "x-flame.svg")
    assert "not on PATH" in reason
    assert not (tmp_path / "x-flame.svg").exists()


def test_profile_stage_falls_back_to_timing_without_perf(monkeypatch, tmp_path) -> None:
    timed = []

    def fake_run_timed(argv, log, stdout_file=None, **kwargs):
        timed.append((argv, stdout_file))
        return {"wall_s": 1.5}

    monkeypatch.setattr(performance, "run_timed", fake_run_timed)
    metrics, details = {}, {"profile": {}}
    performance_profile.profile_stage(
        "propagation",
        ["true-tool", "--x"],
        None,
        None,
        tmp_path,
        tmp_path,
        metrics,
        details,
    )
    entry = details["profile"]["propagation"]
    assert entry["mode"] == "timing"
    assert entry["reason"] == "perf not found on PATH"
    assert entry["wall_s"] == 1.5
    assert timed == [(["true-tool", "--x"], None)]
    assert metrics == {
        "propagation_profiled": 0,
        "propagation_perf_samples": 0,
        "propagation_profile_completed": 1,
    }


def test_profile_stage_records_blocked_perf_and_uses_perf_stat(
    monkeypatch, tmp_path
) -> None:
    calls = []

    def fake_execute(argv, log, stdout_file=None, prefix=()):
        calls.append(list(prefix))
        if "record" in prefix:
            Path(log).write_text(
                "Access to performance monitoring is limited. "
                "Consider /proc/sys/kernel/perf_event_paranoid\n"
            )
            return 255
        return 0

    monkeypatch.setattr(performance_profile, "_execute", fake_execute)
    metrics, details = {}, {"profile": {}}
    performance_profile.profile_stage(
        "trace_stream",
        ["true-tool"],
        None,
        "/usr/bin/perf",
        tmp_path,
        tmp_path,
        metrics,
        details,
    )
    entry = details["profile"]["trace_stream"]
    assert entry["mode"] == "perf-stat"
    assert "blocked" in entry["reason"]
    assert "perf_event_paranoid" in entry["reason"]
    assert calls[1][:2] == ["/usr/bin/perf", "stat"]
    assert performance_profile.PERF_STAT_EVENTS in calls[1]
    assert metrics["trace_stream_profiled"] == 0
    assert metrics["trace_stream_profile_completed"] == 1


def test_profile_stage_record_writes_report_metrics_and_top_symbols(
    monkeypatch, tmp_path
) -> None:
    def fake_execute(argv, log, stdout_file=None, prefix=()):
        data = prefix[prefix.index("-o") + 1]
        Path(data).write_bytes(b"perf")
        return 0

    class Done:
        def __init__(self, stdout, returncode=0):
            self.stdout = stdout
            self.returncode = returncode

    def fake_run(command, **kwargs):
        if "report" in command:
            return Done(REPORT)
        return Done(SCRIPT)

    monkeypatch.setattr(performance_profile, "_execute", fake_execute)
    monkeypatch.setattr(performance_profile.subprocess, "run", fake_run)
    monkeypatch.setattr(performance_profile.shutil, "which", lambda name: None)
    metrics, details = {}, {"profile": {}}
    performance_profile.profile_stage(
        "propagation",
        ["true-tool"],
        None,
        "/usr/bin/perf",
        tmp_path,
        tmp_path,
        metrics,
        details,
    )
    entry = details["profile"]["propagation"]
    assert entry["mode"] == "perf-record"
    assert entry["report"] == "propagation-perf-report.txt"
    assert (tmp_path / "propagation-perf-report.txt").read_text() == REPORT
    assert len(entry["top"]) == 3
    assert metrics == {
        "propagation_profiled": 1,
        "propagation_perf_samples": 3,
        "propagation_top1_percent": 41.25,
        "propagation_profile_completed": 1,
    }


def test_profile_stage_top_symbols_are_capped(monkeypatch, tmp_path) -> None:
    many = "".join(
        f"    {10 - i % 9}.00%  obdeect  libx.so  [.] sym{i}\n" for i in range(30)
    )

    def fake_execute(argv, log, stdout_file=None, prefix=()):
        Path(prefix[prefix.index("-o") + 1]).write_bytes(b"perf")
        return 0

    class Done:
        def __init__(self, stdout):
            self.stdout = stdout
            self.returncode = 0

    monkeypatch.setattr(performance_profile, "_execute", fake_execute)
    monkeypatch.setattr(
        performance_profile.subprocess,
        "run",
        lambda command, **kwargs: Done(many if "report" in command else ""),
    )
    metrics, details = {}, {"profile": {}}
    performance_profile.profile_stage(
        "to_simtel",
        ["true-tool"],
        None,
        "/usr/bin/perf",
        tmp_path,
        tmp_path,
        metrics,
        details,
    )
    assert (
        len(details["profile"]["to_simtel"]["top"]) == performance_profile.TOP_SYMBOLS
    )
