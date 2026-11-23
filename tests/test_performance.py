import json
from pathlib import Path

import pytest

from obdeect_validation.adapters import common, performance


def test_restrict_table_keeps_comments_and_rows_below_limit(tmp_path: Path) -> None:
    source = tmp_path / "source.dat"
    source.write_text("# header\n\n 999 0.1 0.2\n1000 0.3 0.4\n 1200 0.5\n")
    target = tmp_path / "target.dat"
    assert performance.restrict_table(source, target) == 3
    assert target.read_text() == "# header\n\n 999 0.1 0.2\n"


def test_restrict_table_keeps_non_numeric_lines(tmp_path: Path) -> None:
    source = tmp_path / "source.dat"
    source.write_text("lambda 1 2\n500 0.1\n")
    target = tmp_path / "target.dat"
    performance.restrict_table(source, target)
    assert target.read_text() == "lambda 1 2\n500 0.1\n"


def test_stream_weight_parses_obdeect_summary(tmp_path: Path, monkeypatch) -> None:
    class Result:
        stdout = (
            "records:\n  photons 3\ntotal: 10888516 photons, weight sum 54210207.61\n"
        )

    monkeypatch.setattr(performance.subprocess, "run", lambda *args, **kwargs: Result())
    assert performance.stream_weight(tmp_path / "x.gst") == (10888516, 54210207.61)


def test_run_timed_records_child_usage_and_raises_on_failure(tmp_path: Path) -> None:
    record = performance.run_timed(["true"], tmp_path / "log.txt")
    assert record["returncode"] == 0 and record["wall_s"] >= 0 and record["cpu_s"] >= 0
    try:
        performance.run_timed(["false"], tmp_path / "log2.txt")
    except RuntimeError as error:
        assert "exited 1" in str(error)
    else:
        raise AssertionError("failing child must raise")


def test_parse_vmhwm_reads_peak_resident_set_in_kb() -> None:
    status = "Name:\tsleep\nVmPeak:\t  9000 kB\nVmHWM:\t    1536 kB\nVmRSS:\t 1200 kB\n"
    assert performance.parse_vmhwm_kb(status) == 1536
    assert performance.parse_vmhwm_kb("Name:\tzombie\n") is None


def test_rss_sampler_keeps_the_largest_sample() -> None:
    values = iter([100, 300, 200])
    seen = []

    def read():
        seen.append(None)
        return next(values, None)

    sampler = performance.RssSampler(pid=0, read=read, interval=0)
    sampler.start()
    while len(seen) < 5:
        pass
    assert sampler.stop() == 300


def test_rss_sampler_without_samples_returns_none() -> None:
    sampler = performance.RssSampler(pid=0, read=lambda: None, interval=0)
    sampler.start()
    assert sampler.stop() is None


def test_run_timed_small_child_does_not_report_parent_peak(tmp_path: Path) -> None:
    # A large parent: on Linux a forked child inherits the parent's high-water mark through exec.
    parent_ballast = bytearray(200 * 1024 * 1024)
    for offset in range(0, len(parent_ballast), 4096):
        parent_ballast[offset] = 1
    record = performance.run_timed(["sleep", "0.3"], tmp_path / "log.txt")
    assert record["argv"] == ["sleep", "0.3"]
    assert record["peak_rss_mb"] < 100, record
    assert record["peak_rss_source"] in {"proc-VmHWM-sampled", "wait4-ru_maxrss"}


OBS = "atmospheric_transmission-2.0.0-OBS-North.dat"


def test_check_atmosphere_accepts_the_expected_table_by_file_name() -> None:
    text = (
        f"Read atmospheric transmission data from file /results/cfg/model/7.0.0/{OBS}\n"
    )
    assert common.check_atmosphere(text, OBS) == [OBS]


def test_check_atmosphere_rejects_the_fallback_table(tmp_path: Path) -> None:
    text = "Read atmospheric transmission data from file atm_trans_2147_1_10_2_0_2147.dat\n"
    with pytest.raises(RuntimeError, match="atm_trans_2147"):
        common.check_atmosphere(text, OBS)


def test_check_atmosphere_requires_the_read_line() -> None:
    with pytest.raises(RuntimeError, match="no 'Read atmospheric transmission' line"):
        common.check_atmosphere("Starting sim_telarray\n", OBS)


def test_check_atmosphere_log_reads_the_file(tmp_path: Path) -> None:
    log = tmp_path / "reference.log"
    log.write_text(f"Read atmospheric transmission data from file {OBS}\n")
    assert common.check_atmosphere_log(log, OBS) == [OBS]


def fake_runner(walls, cpu=1.0):
    """run_once that appends one child record per call, with the given wall times."""
    walls = iter(walls)

    def run_once():
        performance.RUSAGE.append(
            {
                "argv": ["fake"],
                "wall_s": next(walls),
                "cpu_s": cpu,
                "peak_rss_mb": 10.0,
                "loadavg_before": None,
            }
        )

    return run_once


def test_measure_reports_wall_spread_and_rate_only_with_photons(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(performance, "RUSAGE", [])
    stage_in = tmp_path / "in.gst"
    stage_in.write_bytes(b"12")
    stage_out = tmp_path / "out.gst"
    stage_out.write_bytes(b"1234")
    metrics, details = {}, {"failures": {}, "stages": {}}
    assert performance.measure(
        "stage",
        fake_runner([2.0, 1.0, 4.0], cpu=3.0),
        100.0,
        stage_in,
        stage_out,
        metrics,
        details,
        3,
    )
    assert metrics["stage_wall_s_median"] == 2.0
    assert metrics["stage_wall_s_min"] == 1.0
    assert metrics["stage_wall_s_max"] == 4.0
    assert metrics["stage_wall_s_spread"] == pytest.approx(1.5)  # (4 - 1) / 2
    assert metrics["stage_cpu_s_median"] == 3.0
    assert metrics["stage_photons_per_s"] == pytest.approx(50.0)
    assert metrics["stage_bytes_in"] == 2 and metrics["stage_bytes_out"] == 4
    assert details["stages"]["stage"]["records"][0]["wall_s"] == 2.0

    no_photons = {}
    assert performance.measure(
        "compile",
        fake_runner([1.0]),
        None,
        stage_in,
        stage_out,
        no_photons,
        details,
        1,
    )
    assert no_photons["compile_wall_s_spread"] == 0.0
    assert "compile_photons_per_s" not in no_photons
    assert "compile_photons_in" not in no_photons


def test_measure_records_a_failed_stage(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(performance, "RUSAGE", [])

    def fail():
        raise RuntimeError("bad input")

    metrics, details = {}, {"failures": {}, "stages": {}}
    assert not performance.measure(
        "stage", fail, 1.0, tmp_path, tmp_path, metrics, details, 2
    )
    assert details["failures"] == {"stage": "bad input"}
    assert metrics == {"stage_completed": 0}


def test_run_timed_records_load_average_before(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(performance, "loadavg", lambda: 0.25)
    record = performance.run_timed(["true"], tmp_path / "log.txt")
    assert record["loadavg_before"] == 0.25


def test_loadavg_is_none_without_os_support(monkeypatch) -> None:
    def unavailable():
        raise OSError("no load average")

    monkeypatch.setattr(performance.os, "getloadavg", unavailable)
    assert performance.loadavg() is None


def test_cpu_metrics_share_of_envelope_chain_and_ns_per_photon() -> None:
    metrics = {
        "from_corsika_cpu_s_median": 1.0,
        "propagation_envelope_cpu_s_median": 2.0,
        "trace_stream_cpu_s_median": 4.0,
        "to_simtel_cpu_s_median": 1.0,
        "propagation_cpu_s_median": 8.0,
        "propagation_photons_in": 1e6,
        "propagation_envelope_photons_in": 1e6,
    }
    out = performance.cpu_metrics(
        metrics, ["from_corsika", "propagation_envelope", "propagation"]
    )
    assert out["from_corsika_cpu_share"] == pytest.approx(1 / 8)
    assert out["propagation_envelope_cpu_share"] == pytest.approx(2 / 8)
    assert out["propagation_cpu_share"] == pytest.approx(1.0)  # not in the chain: 8 / 8
    assert out["propagation_ns_per_photon_cpu"] == pytest.approx(8000.0)
    assert "from_corsika_ns_per_photon_cpu" not in out  # no photon count


def test_timing_unreliable_flags_spread_or_busy_host() -> None:
    stages = ["a", "b"]
    calm = {
        "a_wall_s_spread": 0.01,
        "b_wall_s_spread": 0.02,
        "host_loadavg_before": 0.3,
    }
    assert performance.timing_unreliable(calm, stages) == 0
    assert performance.timing_unreliable({**calm, "b_wall_s_spread": 0.06}, stages) == 1
    assert (
        performance.timing_unreliable({**calm, "host_loadavg_before": 1.5}, stages) == 1
    )
    assert performance.timing_unreliable({"a_wall_s_spread": 0.01}, stages) == 0


def test_stages_case_uses_the_air_shower_setup() -> None:
    case = json.loads(
        (Path(__file__).parents[1] / "cases/performance/stages.json").read_text()
    )
    measure = next(c for c in case["commands"] if c["name"] == "measure-stages")
    assert measure["argv"][-2:] == ["--config", "cfg/model/7.0.0/CTAO.cfg"]
    completed = next(m for m in case["metrics"] if m["metric"] == "stages_completed")
    assert completed["expected"] == 6
