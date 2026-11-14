from pathlib import Path

from obdeect_validation.adapters import performance


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
