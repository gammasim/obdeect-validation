import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON = sys.executable


def write_json(path: Path, value: object) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def invoke(*args: str, cwd: Path) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    return subprocess.run(
        [PYTHON, "-m", "obdeect_validation", *args],
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def metric_case(value: float = 2.0, *, case_id: str = "sample") -> dict:
    code = (
        "import json,pathlib; "
        f"pathlib.Path('metrics.json').write_text(json.dumps({{'energy': {value}}}))"
    )
    return {
        "id": case_id,
        "commands": [{"name": "produce", "argv": [PYTHON, "-c", code]}],
        "metrics": [{"metric": "energy", "expected": value, "atol": 0.01}],
    }


def suite_file(tmp_path: Path, case: dict | str) -> Path:
    cases_dir = tmp_path / "cases"
    case_file = (
        write_json(cases_dir / "sample.json", case) if isinstance(case, dict) else case
    )
    return write_json(
        tmp_path / "suite.json", {"id": "unit", "cases": [str(case_file)]}
    )


def test_list_and_doctor_json(tmp_path: Path) -> None:
    write_json(tmp_path / "case.json", metric_case())
    suite = write_json(tmp_path / "suite.json", {"id": "unit", "cases": ["case.json"]})
    listed = invoke("list", "--suite", str(suite), "--json", cwd=tmp_path)
    assert listed.returncode == 0, listed.stderr
    assert "sample" in listed.stdout

    environment = write_json(
        tmp_path / "environment.json", {"env": {}, "requires": [PYTHON]}
    )
    diagnosed = invoke(
        "doctor", "--environment", str(environment), "--json", cwd=tmp_path
    )
    assert diagnosed.returncode == 0, diagnosed.stderr
    assert PYTHON in diagnosed.stdout


def test_run_pass_and_report_html_escapes_case_text(tmp_path: Path) -> None:
    suite = suite_file(
        tmp_path, metric_case(case_id="x</title><script>alert(1)</script>")
    )
    output = tmp_path / "run"
    result = invoke("run", "--suite", str(suite), "--output", str(output), cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert "PASS" in json.dumps(summary)
    html = (output / "report.html").read_text(encoding="utf-8")
    assert "&lt;script&gt;" in html
    assert "<script>alert(1)</script>" not in html
    compared = invoke("report", "--run", str(output), cwd=tmp_path)
    assert compared.returncode == 0, compared.stderr


def test_metric_mismatch_returns_one(tmp_path: Path) -> None:
    suite = suite_file(tmp_path, metric_case(value=3.0))
    result = invoke(
        "run", "--suite", str(suite), "--output", str(tmp_path / "run"), cwd=tmp_path
    )
    assert result.returncode == 0, result.stderr

    # Change the expected value while keeping the produced measurement fixed.
    case_path = tmp_path / "cases" / "sample.json"
    case = metric_case(value=2.0)
    case["metrics"][0]["expected"] = 3.0
    write_json(case_path, case)
    result = invoke(
        "run",
        "--case",
        str(case_path),
        "--output",
        str(tmp_path / "failed"),
        cwd=tmp_path,
    )
    assert result.returncode == 1


def test_missing_prerequisite_returns_three(tmp_path: Path) -> None:
    case = metric_case()
    case["requires"] = ["no-such-obdeect-validation-tool"]
    suite = suite_file(tmp_path, case)
    result = invoke(
        "run", "--suite", str(suite), "--output", str(tmp_path / "run"), cwd=tmp_path
    )
    assert result.returncode == 3
    summary = json.loads(
        (tmp_path / "run" / "summary.json").read_text(encoding="utf-8")
    )
    assert "BLOCKED" in json.dumps(summary)


def test_missing_environment_input_and_module_block_before_command(
    tmp_path: Path,
) -> None:
    marker = tmp_path / "launched"
    case = {
        "id": "environment-check",
        "commands": [
            {
                "name": "mark-launch",
                "argv": [PYTHON, "-c", f"open({str(marker)!r}, 'w').write('started')"],
            }
        ],
    }
    suite = suite_file(tmp_path, case)
    environment = write_json(
        tmp_path / "environment.json",
        {
            "inputs": ["required-input.dat"],
            "modules": ["obdeect_validation_module_missing_for_test"],
        },
    )
    output = tmp_path / "run"
    result = invoke(
        "run",
        "--suite",
        str(suite),
        "--environment",
        str(environment),
        "--output",
        str(output),
        cwd=tmp_path,
    )
    assert result.returncode == 3
    assert not marker.exists()
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    blocked = summary["cases"][0]
    assert blocked["status"] == "BLOCKED"
    assert blocked["steps"] == []
    assert "obdeect_validation_module_missing_for_test" in blocked["missing_modules"]
    assert any(
        path.endswith("required-input.dat") for path in blocked["missing_inputs"]
    )


def test_missing_case_input_blocks_only_that_case(tmp_path: Path) -> None:
    launched = tmp_path / "launched"
    shower = tmp_path / "shower.iact"
    command = {
        "name": "mark-launch",
        "argv": [PYTHON, "-c", f"open({str(launched)!r}, 'w').write('started')"],
    }
    blocked_case = {
        "id": "needs-shower",
        "inputs": [str(shower)],
        "commands": [command],
    }
    other_marker = tmp_path / "other-launched"
    available = tmp_path / "available.iact"
    available.write_bytes(b"ready")
    other_case = {
        "id": "available",
        "inputs": [str(available)],
        "commands": [
            {
                "name": "mark",
                "argv": [
                    PYTHON,
                    "-c",
                    f"open({str(other_marker)!r}, 'w').write('started')",
                ],
            }
        ],
    }
    first_case = write_json(tmp_path / "blocked.json", blocked_case)
    second_case = write_json(tmp_path / "available.json", other_case)
    suite = write_json(
        tmp_path / "suite.json",
        {"id": "unit", "cases": [str(first_case), str(second_case)]},
    )
    output = tmp_path / "run"
    result = invoke("run", "--suite", str(suite), "--output", str(output), cwd=tmp_path)
    assert result.returncode == 3
    assert not launched.exists()
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    case = summary["cases"][0]
    assert case["status"] == "BLOCKED"
    assert case["steps"] == []
    assert case["missing_inputs"] == [str(shower)]
    assert other_marker.exists()
    assert summary["cases"][1]["status"] == "PASS"

    shower.write_bytes(b"x")  # present: the case runs
    output = tmp_path / "run-present"
    result = invoke("run", "--suite", str(suite), "--output", str(output), cwd=tmp_path)
    assert result.returncode == 0
    assert launched.exists()
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert summary["cases"][0]["missing_inputs"] == []


def test_case_inputs_must_be_a_list_of_paths(tmp_path: Path) -> None:
    command = {"name": "noop", "argv": [PYTHON, "-c", "pass"]}
    valid = suite_file(tmp_path, {"id": "ok", "inputs": [], "commands": [command]})
    assert invoke("list", "--suite", str(valid), cwd=tmp_path).returncode == 0
    for bad in ("shower.iact", [""], [1]):
        case = {"id": "bad-inputs", "inputs": bad, "commands": [command]}
        suite = suite_file(tmp_path, case)
        result = invoke("list", "--suite", str(suite), cwd=tmp_path)
        assert result.returncode == 2, bad


def test_timeout_returns_execution_error(tmp_path: Path) -> None:
    case = {
        "id": "timeout",
        "timeout": 0.05,
        "commands": [
            {"name": "wait", "argv": [PYTHON, "-c", "import time; time.sleep(5)"]}
        ],
    }
    suite = suite_file(tmp_path, case)
    result = invoke(
        "run", "--suite", str(suite), "--output", str(tmp_path / "run"), cwd=tmp_path
    )
    assert result.returncode == 2


def test_timeout_terminates_command_descendants(tmp_path: Path) -> None:
    marker = tmp_path / "child-finished"
    child = f"import time,pathlib; time.sleep(0.8); pathlib.Path({str(marker)!r}).write_text('alive')"
    parent = f"import subprocess,sys,time; subprocess.Popen([sys.executable, '-c', {child!r}]); time.sleep(5)"
    case = {
        "id": "timeout-child",
        "timeout": 0.1,
        "commands": [{"name": "spawn", "argv": [PYTHON, "-c", parent]}],
    }
    suite = suite_file(tmp_path, case)
    result = invoke(
        "run", "--suite", str(suite), "--output", str(tmp_path / "run"), cwd=tmp_path
    )
    assert result.returncode == 2
    subprocess.run([PYTHON, "-c", "import time; time.sleep(1.0)"], check=True)
    assert not marker.exists()


def test_sigint_writes_partial_reports_and_stops_process_group(tmp_path: Path) -> None:
    child_pid_file = tmp_path / "child.pid"
    child_marker = tmp_path / "child-survived"
    child_code = (
        "import os,pathlib,time; "
        f"pathlib.Path({str(child_pid_file)!r}).write_text(str(os.getpid())); "
        "time.sleep(1.5); "
        f"pathlib.Path({str(child_marker)!r}).write_text('survived')"
    )
    parent_code = (
        "import pathlib,subprocess,sys,time; "
        f"subprocess.Popen([sys.executable, '-c', {child_code!r}]); "
        "time.sleep(60)"
    )
    case = {
        "id": "interrupt-me",
        "commands": [{"name": "spawn-child", "argv": [PYTHON, "-c", parent_code]}],
    }
    second = {
        "id": "must-not-start",
        "commands": [{"name": "noop", "argv": [PYTHON, "-c", "pass"]}],
    }
    suite = write_json(
        tmp_path / "suite.json", {"id": "interrupt", "cases": [case, second]}
    )
    output = tmp_path / "run"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    runner = subprocess.Popen(
        [
            PYTHON,
            "-m",
            "obdeect_validation",
            "run",
            "--suite",
            str(suite),
            "--output",
            str(output),
        ],
        cwd=tmp_path,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    deadline = time.monotonic() + 10
    try:
        while (
            not child_pid_file.exists()
            and runner.poll() is None
            and time.monotonic() < deadline
        ):
            time.sleep(0.05)
        assert child_pid_file.exists(), (
            "runner never started the command's child process"
        )
        time.sleep(0.2)
        runner.send_signal(signal.SIGINT)
        stdout, stderr = runner.communicate(timeout=10)
        assert runner.returncode == 130, f"stdout={stdout!r} stderr={stderr!r}"
    finally:
        if runner.poll() is None:
            runner.kill()
            runner.communicate(timeout=5)

    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert summary["interrupted"] is True
    assert [case["status"] for case in summary["cases"]] == ["ERROR", "BLOCKED"]
    assert summary["cases"][0]["steps"][0]["error"] == "interrupted"
    assert (output / "junit.xml").is_file()
    assert (output / "report.html").is_file()
    assert "interrupt-me" in (output / "report.html").read_text(encoding="utf-8")
    time.sleep(1.7)
    assert not child_marker.exists()
    child_pid = int(child_pid_file.read_text(encoding="utf-8"))
    proc_stat = Path(f"/proc/{child_pid}/stat")
    if proc_stat.exists():
        state = proc_stat.read_text(encoding="utf-8").rsplit(")", 1)[1].split()[0]
        assert state == "Z", (
            f"child process {child_pid} is still running (state {state})"
        )
    else:
        with pytest.raises(ProcessLookupError):
            os.kill(child_pid, 0)


def test_suite_loads_sibling_cases_through_parent_path(tmp_path: Path) -> None:
    write_json(tmp_path / "cases" / "sample.json", metric_case())
    suite = write_json(
        tmp_path / "suites" / "suite.json",
        {"id": "siblings", "cases": ["../cases/sample.json"]},
    )
    result = invoke("list", "--suite", str(suite), "--json", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert "sample" in result.stdout


def test_command_only_infrastructure_case_passes_without_scientific_metrics(
    tmp_path: Path,
) -> None:
    case = {
        "id": "tool-probe",
        "commands": [{"name": "version", "argv": [PYTHON, "--version"]}],
    }
    suite = suite_file(tmp_path, case)
    result = invoke(
        "run", "--suite", str(suite), "--output", str(tmp_path / "run"), cwd=tmp_path
    )
    assert result.returncode == 0, result.stderr
    summary = json.loads(
        (tmp_path / "run" / "summary.json").read_text(encoding="utf-8")
    )
    assert summary["cases"][0]["status"] == "PASS"
    assert summary["cases"][0]["metrics"] == []


def test_resume_reuses_pass_but_invalidates_changed_case(tmp_path: Path) -> None:
    marker = tmp_path / "runs.txt"
    case = metric_case()
    case["commands"][0]["argv"][-1] += f"; open({str(marker)!r}, 'a').write('run\\n')"
    suite = suite_file(tmp_path, case)
    output = tmp_path / "run"
    first = invoke("run", "--suite", str(suite), "--output", str(output), cwd=tmp_path)
    assert first.returncode == 0, first.stderr
    second = invoke(
        "run", "--suite", str(suite), "--output", str(output), "--resume", cwd=tmp_path
    )
    assert second.returncode == 0, second.stderr
    assert marker.read_text(encoding="utf-8").splitlines() == ["run"]

    case_path = tmp_path / "cases" / "sample.json"
    changed = metric_case()
    changed["name"] = "changed fingerprint"
    changed["commands"][0]["argv"][-1] += (
        f"; open({str(marker)!r}, 'a').write('run\\n')"
    )
    write_json(case_path, changed)
    suite = write_json(
        tmp_path / "suite.json", {"id": "unit", "cases": ["cases/sample.json"]}
    )
    third = invoke(
        "run", "--suite", str(suite), "--output", str(output), "--resume", cwd=tmp_path
    )
    assert third.returncode == 0, third.stderr
    assert marker.read_text(encoding="utf-8").splitlines() == ["run", "run"]
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert "PASS" in json.dumps(summary)


def test_resume_invalidates_changed_and_missing_case_input(tmp_path: Path) -> None:
    marker = tmp_path / "runs.txt"
    input_file = tmp_path / "input.dat"
    input_file.write_text("first", encoding="utf-8")
    case = metric_case()
    case["inputs"] = [str(input_file)]
    case["commands"][0]["argv"][-1] += f"; open({str(marker)!r}, 'a').write('run\\n')"
    suite = suite_file(tmp_path, case)
    output = tmp_path / "run"
    args = ("run", "--suite", str(suite), "--output", str(output))
    assert invoke(*args, cwd=tmp_path).returncode == 0
    assert invoke(*args, "--resume", cwd=tmp_path).returncode == 0
    assert marker.read_text().splitlines() == ["run"]
    input_file.write_text("second", encoding="utf-8")
    assert invoke(*args, "--resume", cwd=tmp_path).returncode == 0
    assert marker.read_text().splitlines() == ["run", "run"]
    input_file.unlink()
    assert invoke(*args, "--resume", cwd=tmp_path).returncode == 3
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert summary["cases"][0]["status"] == "BLOCKED"
    assert summary["cases"][0]["missing_inputs"] == [str(input_file)]


def test_resume_invalidates_when_tracked_source_changes(tmp_path: Path) -> None:
    def git(*args: str) -> None:
        subprocess.run(["git", *args], cwd=tmp_path, check=True, capture_output=True)

    git("init", "-q")
    git("config", "user.email", "validation-tests@example.invalid")
    git("config", "user.name", "Validation tests")
    source = tmp_path / "source.txt"
    source.write_text("before\n", encoding="utf-8")
    git("add", "source.txt")
    git("commit", "-qm", "initial source")

    marker = tmp_path / "source-runs.txt"
    case = metric_case()
    case["commands"][0]["argv"][-1] += f"; open({str(marker)!r}, 'a').write('run\\n')"
    suite = suite_file(tmp_path, case)
    output = tmp_path / "run"
    first = invoke("run", "--suite", str(suite), "--output", str(output), cwd=tmp_path)
    assert first.returncode == 0, first.stderr

    source.write_text("after\n", encoding="utf-8")
    resumed = invoke(
        "run", "--suite", str(suite), "--output", str(output), "--resume", cwd=tmp_path
    )
    assert resumed.returncode == 0, resumed.stderr
    assert marker.read_text(encoding="utf-8").splitlines() == ["run", "run"]


def test_rerun_does_not_accept_stale_metrics(tmp_path: Path) -> None:
    case = metric_case()
    suite = suite_file(tmp_path, case)
    output = tmp_path / "run"
    first = invoke("run", "--suite", str(suite), "--output", str(output), cwd=tmp_path)
    assert first.returncode == 0, first.stderr

    case["commands"][0]["argv"] = [PYTHON, "-c", "pass"]
    write_json(tmp_path / "cases" / "sample.json", case)
    rerun = invoke("run", "--suite", str(suite), "--output", str(output), cwd=tmp_path)
    assert rerun.returncode == 2
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert summary["cases"][0]["status"] == "ERROR"


def test_invalid_config_returns_two(tmp_path: Path) -> None:
    suite = write_json(
        tmp_path / "suite.json", {"id": "bad", "cases": [{"id": "broken"}]}
    )
    result = invoke(
        "run", "--suite", str(suite), "--output", str(tmp_path / "run"), cwd=tmp_path
    )
    assert result.returncode == 2
