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
