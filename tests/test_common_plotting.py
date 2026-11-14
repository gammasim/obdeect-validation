"""Shared adapter helpers, the plotting module and the runner fixes made in review."""

import hashlib
import json
import signal
import subprocess
import sys
import time

import numpy as np
import pytest

from obdeect_validation import runner
from obdeect_validation.adapters import common


def test_event_deltas_counts_unmatched_events_in_full():
    assert common.event_deltas([10.0, 5.0], [9.0]) == [1.0, 5.0]
    assert common.event_deltas([1.0], [1.0, 7.0]) == [0.0, 7.0]


def test_models_root_prefers_override_then_mounted_checkout(monkeypatch, tmp_path):
    monkeypatch.setenv("VALIDATION_MODELS", str(tmp_path))
    assert common.models_root() == tmp_path
    monkeypatch.delenv("VALIDATION_MODELS")
    assert common.models_root().name in ("models", "simulation-models")


def test_write_json_rejects_nan_and_converts_numpy(tmp_path):
    common.write_json(tmp_path / "a.json", {"x": np.float64(1.5), "p": tmp_path})
    assert json.loads((tmp_path / "a.json").read_text())["x"] == 1.5
    with pytest.raises(ValueError):
        common.write_json(tmp_path / "b.json", {"x": float("nan")})


def test_command_error_carries_log_tail(tmp_path):
    script = "import sys; print('boom detail', file=sys.stderr); sys.exit(3)"
    with pytest.raises(RuntimeError, match="exited 3.*boom detail"):
        common.command([sys.executable, "-c", script], tmp_path / "log")


def test_sha256_matches_hashlib(tmp_path):
    path = tmp_path / "f"
    path.write_bytes(b"x" * 3_000_000)
    assert common.sha256(path) == hashlib.sha256(b"x" * 3_000_000).hexdigest()


def test_histograms_skip_empty_samples_and_handle_constant_values(tmp_path):
    plotting = pytest.importorskip("obdeect_validation.adapters.plotting")
    fig, axis = plotting.figure()
    plotting.histograms(
        axis,
        [("a", np.zeros(10), plotting.NATIVE), ("b", np.array([]), plotting.REFERENCE)],
        "x",
    )
    plotting.save(fig, tmp_path / "h.png")
    assert (tmp_path / "h.png").stat().st_size > 0


def test_camera_bars_state_missing_samples(tmp_path):
    plotting = pytest.importorskip("obdeect_validation.adapters.plotting")
    fig, axis = plotting.figure()
    plotting.camera_bars(axis, np.array([]), np.array([3.0, 5.0]), "ref")
    plotting.save(fig, tmp_path / "bars.png")
    assert (tmp_path / "bars.png").is_file()


def test_case_colors_are_stable_per_family():
    plotting = pytest.importorskip("obdeect_validation.adapters.plotting")
    assert plotting.case_color("LSTN-01") == plotting.case_color("lst")
    assert plotting.case_color("LSTN-01") != plotting.case_color("MSTN-04")
    assert plotting.case_color("unknown-case") == plotting.case_color("unknown-case")


def test_git_provenance_hashes_binary_diffs(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    repo.mkdir()
    git = ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t"]
    subprocess.run([*git, "init", "-q"], check=True)
    (repo / "blob.bin").write_bytes(b"\x00\x01")
    subprocess.run([*git, "add", "."], check=True)
    subprocess.run([*git, "commit", "-qm", "x"], check=True)
    (repo / "blob.bin").write_bytes(b"\xb0\xff\x00\x02")  # not UTF-8
    monkeypatch.chdir(repo)
    entry = next(
        e for e in runner._git_provenance() if e["path"] == str(repo.resolve())
    )
    assert entry["dirty"] is True
    assert len(entry["diff_sha256"]) == 64


def test_stop_process_group_survives_permission_error(monkeypatch):
    """macOS reports PermissionError for a group that only holds zombies."""
    process = subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(30)"], start_new_session=True
    )
    real_killpg = runner.os.killpg

    def killpg(pid, sig):
        if sig == 0:
            raise PermissionError
        real_killpg(pid, sig)

    monkeypatch.setattr(runner.os, "killpg", killpg)
    started = time.monotonic()
    runner._stop_process_group(process)
    assert process.returncode is not None
    assert time.monotonic() - started < 10
    assert process.returncode in (-signal.SIGTERM, -signal.SIGKILL)


def test_report_html_shows_status_reason_and_figures(tmp_path):
    case_dir = tmp_path / "case-a"
    case_dir.mkdir()
    (case_dir / "plot.png").write_bytes(b"png")
    summary = {
        "suite": "s",
        "status": "BLOCKED",
        "cases": [
            {
                "id": "case-a",
                "name": "Case A",
                "status": "BLOCKED",
                "reason": "run interrupted",
                "missing": ["nsbls"],
                "steps": [],
                "metrics": [],
            }
        ],
    }
    runner._write_reports(tmp_path, summary)
    page = (tmp_path / "report.html").read_text()
    assert "class='s BLOCKED'" in page
    assert "run interrupted" in page and "nsbls" in page
    assert "case-a/plot.png" in page


def _fake_stream(records):
    return lambda _path: iter(records)


def _photons(n=2):
    columns = {k: np.zeros(n) for k in ("x", "y", "cx", "cy", "t", "wavelength")}
    return {"type": "photons", "columns": {**columns, "w": np.ones(n)}}


def test_native_reader_rejects_photons_before_first_event(monkeypatch):
    from obdeect_validation.adapters import readers

    monkeypatch.setattr(readers, "Reader", _fake_stream([_photons()]))
    with pytest.raises(ValueError, match="before its first event"):
        readers.read_native("x")
    counts, columns = readers.read_native("x", require_events=False)
    assert counts == [2.0] and columns["w"].size == 2


def test_native_reader_counts_weight_per_event(monkeypatch):
    from obdeect_validation.adapters import readers

    records = [{"type": "event"}, _photons(3), {"type": "event"}, _photons(1)]
    monkeypatch.setattr(readers, "Reader", _fake_stream(records))
    assert readers.read_native("x")[0] == [3.0, 1.0]


def test_camera_metrics_need_two_events_per_sample():
    from obdeect_validation.adapters import source_camera as sc

    one = {1: (np.array([5.0]), np.ones(4))}
    two = {1: (np.array([5.0, 6.0]), np.ones(4))}
    with pytest.raises(ValueError, match="at least 2"):
        sc.camera_metrics(one, two)


def test_extinction_figure_colours_follow_the_photon_wavelengths(tmp_path):
    plotting = pytest.importorskip("obdeect_validation.adapters.plotting")
    from obdeect_validation.adapters import extinction_reference as er

    expected = np.linspace(0.1, 0.9, 6)
    er._plot(
        expected,
        expected,
        np.array([300, 300, 400, 400, 500, 500.0]),
        tmp_path / "e.png",
    )
    assert (tmp_path / "e.png").is_file() and plotting.CYCLE
