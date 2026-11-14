"""Configuration validation, command execution, and report writing."""

from __future__ import annotations

import hashlib
import copy
import importlib.metadata
import html
import json
import math
import os
import platform
import signal
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

from . import __version__


class ConfigError(ValueError):
    """Raised when a suite or case violates the documented JSON contract."""


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigError(f"cannot read JSON config {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ConfigError(f"{path}: top level must be a JSON object")
    return value


def _keys(
    value: dict[str, Any], required: set[str], allowed: set[str], where: str
) -> None:
    missing, extra = required - value.keys(), value.keys() - allowed
    if missing or extra:
        details = []
        if missing:
            details.append("missing " + ", ".join(sorted(missing)))
        if extra:
            details.append("unknown " + ", ".join(sorted(extra)))
        raise ConfigError(f"{where}: " + "; ".join(details))


def _string(value: Any, where: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"{where}: expected a non-empty string")
    return value


def validate_case(case: dict[str, Any], where: str = "case") -> dict[str, Any]:
    _keys(
        case,
        {"id", "commands"},
        {"id", "name", "requires", "timeout", "commands", "metrics"},
        where,
    )
    _string(case["id"], f"{where}.id")
    if not isinstance(case["commands"], list) or not case["commands"]:
        raise ConfigError(f"{where}.commands: expected a non-empty list")
    if "name" in case:
        _string(case["name"], f"{where}.name")
    if "requires" in case:
        _require_list(case["requires"], f"{where}.requires")
    _timeout(case.get("timeout", 3600), f"{where}.timeout")
    for i, command in enumerate(case["commands"]):
        label = f"{where}.commands[{i}]"
        if not isinstance(command, dict):
            raise ConfigError(f"{label}: expected an object")
        _keys(command, {"name", "argv"}, {"name", "argv", "timeout", "cwd"}, label)
        _string(command["name"], f"{label}.name")
        argv = command["argv"]
        if (
            not isinstance(argv, list)
            or not argv
            or any(not isinstance(arg, str) for arg in argv)
        ):
            raise ConfigError(f"{label}.argv: expected a non-empty list of strings")
        if not argv[0]:
            raise ConfigError(f"{label}.argv[0]: executable cannot be empty")
        if "timeout" in command:
            _timeout(command["timeout"], f"{label}.timeout")
        if "cwd" in command:
            rel = Path(_string(command["cwd"], f"{label}.cwd"))
            if rel.is_absolute() or ".." in rel.parts:
                raise ConfigError(
                    f"{label}.cwd: must stay within the case output directory"
                )
    metrics = case.get("metrics", [])
    if not isinstance(metrics, list):
        raise ConfigError(f"{where}.metrics: expected a list")
    for i, check in enumerate(metrics):
        label = f"{where}.metrics[{i}]"
        if not isinstance(check, dict):
            raise ConfigError(f"{label}: expected an object")
        _keys(
            check, {"metric", "expected", "atol"}, {"metric", "expected", "atol"}, label
        )
        _string(check["metric"], f"{label}.metric")
        for key in ("expected", "atol"):
            if (
                isinstance(check[key], bool)
                or not isinstance(check[key], (int, float))
                or not math.isfinite(check[key])
            ):
                raise ConfigError(f"{label}.{key}: expected a finite number")
        if check["atol"] < 0:
            raise ConfigError(f"{label}.atol: must be non-negative")
    return case


def _timeout(value: Any, where: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise ConfigError(f"{where}: expected a positive finite number")
    return float(value)


def _require_list(value: Any, where: str) -> list[str]:
    if not isinstance(value, list) or any(
        not isinstance(item, str) or not item for item in value
    ):
        raise ConfigError(f"{where}: expected a list of executable names")
    return value


def load_suite(path: Path) -> dict[str, Any]:
    suite = _read_json(path)
    _keys(suite, {"id", "cases"}, {"id", "cases"}, str(path))
    _string(suite["id"], f"{path}.id")
    if not isinstance(suite["cases"], list) or not suite["cases"]:
        raise ConfigError(f"{path}.cases: expected a non-empty list")
    cases = []
    for i, entry in enumerate(suite["cases"]):
        if isinstance(entry, str):
            case_path = (path.parent / entry).resolve()
            cases.append(validate_case(_read_json(case_path), str(case_path)))
        elif isinstance(entry, dict):
            cases.append(validate_case(entry, f"{path}.cases[{i}]"))
        else:
            raise ConfigError(
                f"{path}.cases[{i}]: expected a case object or relative path"
            )
    ids = [case["id"] for case in cases]
    if len(set(ids)) != len(ids):
        raise ConfigError(f"{path}: case ids must be unique")
    if len({_safe_id(case_id) for case_id in ids}) != len(ids):
        raise ConfigError(
            f"{path}: case ids collide after conversion to output directory names"
        )
    return {"id": suite["id"], "cases": cases}


def _environment(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {"env": {}, "requires": [], "modules": [], "inputs": []}
    value = _read_json(path)
    _keys(value, set(), {"env", "requires", "modules", "inputs"}, str(path))
    env = value.get("env", {})
    if not isinstance(env, dict) or any(
        not isinstance(k, str) or not k or not isinstance(v, str)
        for k, v in env.items()
    ):
        raise ConfigError(
            f"{path}.env: expected a mapping of variable names to strings"
        )
    requires = _require_list(value.get("requires", []), f"{path}.requires")
    modules = _require_list(value.get("modules", []), f"{path}.modules")
    inputs = value.get("inputs", [])
    if not isinstance(inputs, list) or any(
        not isinstance(item, str) or not item for item in inputs
    ):
        raise ConfigError(f"{path}.inputs: expected a list of file paths")
    return {
        "env": env,
        "requires": requires,
        "modules": modules,
        "inputs": [((path.parent / item).resolve()) for item in inputs],
    }


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_provenance() -> list[dict[str, Any]]:
    roots: set[Path] = set()
    here = Path.cwd().resolve()
    candidates = [here, *here.parents]
    candidates.extend(child for child in here.parent.iterdir() if child.is_dir())
    for candidate in candidates:
        try:
            check = subprocess.run(
                ["git", "-C", str(candidate), "rev-parse", "--show-toplevel"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if check.returncode == 0:
                roots.add(Path(check.stdout.strip()).resolve())
        except (OSError, subprocess.TimeoutExpired):
            continue
    result = []
    for root in sorted(roots):

        def git(*args: str) -> str | None:
            # Large checkouts (model repositories) can be slow; record unknown instead of failing.
            try:
                proc = subprocess.run(
                    ["git", "-C", str(root), *args],
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
            except subprocess.TimeoutExpired:
                return None
            return proc.stdout.strip() if proc.returncode == 0 else ""

        head = git("rev-parse", "HEAD")
        status = git("status", "--porcelain")
        diff = git("diff", "HEAD", "--binary") if status else ""
        result.append(
            {
                "path": str(root),
                "head": head or None,
                "dirty": None if status is None else bool(status),
                "diff_sha256": None
                if not status or diff is None
                else hashlib.sha256(diff.encode()).hexdigest(),
            }
        )
    return result


def _tool_fingerprints(
    case: dict[str, Any], environment: dict[str, Any]
) -> tuple[dict[str, str | None], list[str]]:
    names = set(case.get("requires", [])) | set(environment["requires"])
    names.update(command["argv"][0] for command in case["commands"])
    found: dict[str, str | None] = {}
    missing = []
    for name in sorted(names):
        executable = shutil.which(
            name, path=environment["env"].get("PATH", os.environ.get("PATH"))
        )
        if executable is None:
            found[name] = None
            missing.append(name)
        else:
            found[name] = _hash_file(Path(executable))
    return found, missing


def _missing_modules(environment: dict[str, Any]) -> list[str]:
    missing = []
    env = {**os.environ, **environment["env"]}
    for module in environment["modules"]:
        try:
            proc = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    "import importlib,sys; importlib.import_module(sys.argv[1])",
                    module,
                ],
                env=env,
                capture_output=True,
                timeout=30,
                check=False,
            )
            available = proc.returncode == 0
        except subprocess.TimeoutExpired:
            available = False
        if not available:
            missing.append(module)
    return missing


def _stop_process_group(process: subprocess.Popen[bytes]) -> None:
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    deadline = time.monotonic() + 2
    while time.monotonic() < deadline:
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            break
        time.sleep(0.05)
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def _provenance(environment: dict[str, Any]) -> dict[str, Any]:
    source_dir = Path(__file__).resolve().parent
    sources = {
        str(path.relative_to(source_dir)): _hash_file(path)
        for path in sorted(source_dir.rglob("*"))
        if path.is_file() and "__pycache__" not in path.parts
    }
    inputs = {
        str(path): _hash_file(path) if path.is_file() else None
        for path in environment["inputs"]
    }
    effective_env = os.environ.copy()
    effective_env.update(environment["env"])
    env_fingerprint = hashlib.sha256(
        json.dumps(dict(sorted(effective_env.items())), separators=(",", ":")).encode()
    ).hexdigest()
    distributions = {}
    for name in (
        "obdeect",
        "obdeect-stream",
        "obdeect-propagation",
        "obdeect-sources",
        "obdeect-ray-tracing",
        "gammasimtools",
        "numpy",
        "astropy",
        "eventio",
    ):
        try:
            distribution = importlib.metadata.distribution(name)
            native = {}
            for entry in distribution.files or ():
                if (
                    str(entry).endswith((".so", ".dylib", ".dll", ".pyd"))
                    or "bin" in Path(str(entry)).parts
                ):
                    installed = Path(distribution.locate_file(entry))
                    if installed.is_file():
                        native[str(entry)] = _hash_file(installed)
            distributions[name] = {
                "version": distribution.version,
                "native_sha256": native,
            }
        except importlib.metadata.PackageNotFoundError:
            distributions[name] = None
    return {
        "runner_version": __version__,
        "python": sys.version,
        "platform": platform.platform(),
        "environment_sha256": env_fingerprint,
        "inputs": inputs,
        "distributions": distributions,
        "sources": sources,
        "git": _git_provenance(),
    }


def _fingerprint(
    case: dict[str, Any], provenance: dict[str, Any], tools: dict[str, str | None]
) -> str:
    payload = {"case": case, "provenance": provenance, "executables": tools}
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def doctor(path: Path) -> dict[str, Any]:
    environment = _environment(path)
    checks = [{"name": "python", "status": "PASS", "detail": sys.executable}]
    for name in sorted(set(environment["requires"])):
        resolved = shutil.which(
            name, path=environment["env"].get("PATH", os.environ.get("PATH"))
        )
        checks.append(
            {
                "name": name,
                "status": "PASS" if resolved else "BLOCKED",
                "detail": resolved or "executable not found",
            }
        )
    for module in sorted(set(environment["modules"])):
        missing = _missing_modules({**environment, "modules": [module]})
        checks.append(
            {
                "name": module,
                "status": "BLOCKED" if missing else "PASS",
                "detail": "module not importable" if missing else "import available",
            }
        )
    for input_path in environment["inputs"]:
        checks.append(
            {
                "name": str(input_path),
                "status": "PASS" if input_path.is_file() else "BLOCKED",
                "detail": "file available"
                if input_path.is_file()
                else "file not found",
            }
        )
    return {
        "status": "BLOCKED"
        if any(c["status"] == "BLOCKED" for c in checks)
        else "PASS",
        "checks": checks,
    }


def _safe_id(value: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value)
    if cleaned in ("", ".", ".."):
        raise ConfigError(f"unsafe case id: {value!r}")
    return cleaned


def _expand_case(case: dict[str, Any]) -> dict[str, Any]:
    expanded = copy.deepcopy(case)
    for command in expanded["commands"]:
        command["argv"] = [
            arg.replace("{python}", sys.executable) for arg in command["argv"]
        ]
    return expanded


def _run_case(
    case: dict[str, Any],
    output: Path,
    environment: dict[str, Any],
    provenance: dict[str, Any],
) -> dict[str, Any]:
    tools, missing = _tool_fingerprints(case, environment)
    fingerprint = _fingerprint(case, provenance, tools)
    missing_modules = _missing_modules(environment)
    missing_inputs = [str(path) for path in environment["inputs"] if not path.is_file()]
    result: dict[str, Any] = {
        "id": case["id"],
        "name": case.get("name", case["id"]),
        "fingerprint": fingerprint,
        "status": "BLOCKED" if missing or missing_inputs or missing_modules else "PASS",
        "missing": missing,
        "missing_modules": missing_modules,
        "missing_inputs": missing_inputs,
        "steps": [],
        "metrics": [],
    }
    if missing or missing_inputs or missing_modules:
        return result
    case_dir = output / _safe_id(case["id"])
    case_dir.mkdir(parents=True, exist_ok=True)
    base_env = os.environ.copy()
    base_env.update(environment["env"])
    for command in case["commands"]:
        step_dir = case_dir / command.get("cwd", ".")
        step_dir.mkdir(parents=True, exist_ok=True)
        timeout = command.get("timeout", case.get("timeout", 3600))
        started = time.monotonic()
        step_number = len(result["steps"]) + 1
        log_prefix = f"{step_number:02d}-{_safe_id(command['name'])}"
        stdout_path = step_dir / f"{log_prefix}.stdout.log"
        stderr_path = step_dir / f"{log_prefix}.stderr.log"
        step = {
            "name": command["name"],
            "argv": command["argv"],
            "cwd": str(step_dir),
            "status": "PASS",
        }
        step["stdout_log"] = str(stdout_path.relative_to(output))
        step["stderr_log"] = str(stderr_path.relative_to(output))
        try:
            metrics_file = step_dir / "metrics.json"
            if metrics_file.exists():
                metrics_file.unlink()
            with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
                proc = subprocess.Popen(
                    command["argv"],
                    cwd=step_dir,
                    env=base_env,
                    stdout=stdout,
                    stderr=stderr,
                    start_new_session=True,
                )
                try:
                    returncode = proc.wait(timeout=timeout)
                except KeyboardInterrupt:
                    _stop_process_group(proc)
                    step.update(
                        {"returncode": None, "status": "ERROR", "error": "interrupted"}
                    )
                    step["duration_seconds"] = round(time.monotonic() - started, 6)
                    result["steps"].append(step)
                    result["status"] = "ERROR"
                    result["interrupted"] = True
                    break
                except subprocess.TimeoutExpired:
                    _stop_process_group(proc)
                    raise
            step.update(
                {
                    "returncode": returncode,
                    "status": "PASS" if returncode == 0 else "ERROR",
                }
            )
        except subprocess.TimeoutExpired:
            step.update(
                {
                    "returncode": None,
                    "status": "ERROR",
                    "error": f"timed out after {timeout:g} seconds",
                }
            )
        except OSError as exc:
            step.update({"returncode": None, "status": "ERROR", "error": str(exc)})
        step["duration_seconds"] = round(time.monotonic() - started, 6)
        result["steps"].append(step)
        if step["status"] != "PASS":
            result["status"] = "ERROR"
            break
    if result["status"] == "PASS" and case.get("metrics"):
        metric_files = [
            case_dir / command.get("cwd", ".") / "metrics.json"
            for command in case["commands"]
        ]
        values: dict[str, Any] = {}
        for metric_file in metric_files:
            if metric_file.exists():
                try:
                    loaded = json.loads(metric_file.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError) as exc:
                    result["status"] = "ERROR"
                    result["metric_error"] = f"invalid {metric_file}: {exc}"
                    break
                if not isinstance(loaded, dict):
                    result["status"] = "ERROR"
                    result["metric_error"] = f"{metric_file} must contain an object"
                    break
                values.update(loaded)
        if result["status"] == "PASS":
            for check in case["metrics"]:
                actual = values.get(check["metric"])
                if (
                    isinstance(actual, bool)
                    or not isinstance(actual, (int, float))
                    or not math.isfinite(actual)
                ):
                    result["status"] = "ERROR"
                    result["metric_error"] = (
                        f"missing or non-numeric metric {check['metric']!r}"
                    )
                    break
                passed = abs(actual - check["expected"]) <= check["atol"]
                result["metrics"].append(
                    {**check, "actual": actual, "status": "PASS" if passed else "FAIL"}
                )
                if not passed:
                    result["status"] = "FAIL"
    return result


def _write_reports(output: Path, summary: dict[str, Any]) -> None:
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    cases = summary["cases"]
    suite = ET.Element(
        "testsuite",
        name=summary["suite"],
        tests=str(len(cases)),
        failures=str(sum(c["status"] == "FAIL" for c in cases)),
        errors=str(sum(c["status"] == "ERROR" for c in cases)),
        skipped=str(sum(c["status"] == "BLOCKED" for c in cases)),
    )
    for case in cases:
        node = ET.SubElement(
            suite, "testcase", name=case["id"], classname=summary["suite"]
        )
        if case["status"] == "FAIL":
            ET.SubElement(
                node, "failure", message="scientific metric outside tolerance"
            )
        elif case["status"] in ("ERROR", "BLOCKED"):
            ET.SubElement(
                node,
                "error" if case["status"] == "ERROR" else "skipped",
                message=case["status"],
            )
    ET.ElementTree(suite).write(
        output / "junit.xml", encoding="utf-8", xml_declaration=True
    )
    rows = []
    for case in cases:
        details = []
        for step in case.get("steps", []):
            stdout = html.escape(step["stdout_log"])
            stderr = html.escape(step["stderr_log"])
            details.append(
                f"<li>{html.escape(step['name'])}: {step['status']} (<a href='{stdout}'>stdout</a>, <a href='{stderr}'>stderr</a>)</li>"
            )
        for metric in case.get("metrics", []):
            details.append(
                f"<li>{html.escape(metric['metric'])}: {metric['actual']} ({metric['status']})</li>"
            )
        rows.append(
            f"<tr><td>{html.escape(case['id'])}</td><td>{html.escape(case['name'])}</td><td>{case['status']}</td><td><ul>{''.join(details)}</ul></td></tr>"
        )
    page = (
        "<!doctype html><meta charset='utf-8'><title>Validation report</title><h1>"
        + html.escape(summary["suite"])
        + "</h1><p>Status: "
        + summary["status"]
        + " · <a href='summary.json'>JSON summary</a></p><table><thead><tr><th>Case</th><th>Name</th><th>Status</th><th>Details</th></tr></thead><tbody>"
        + "\n".join(rows)
        + "</tbody></table>\n"
    )
    (output / "report.html").write_text(page, encoding="utf-8")


def run(
    suite_path: Path | None,
    case_path: Path | None,
    output: Path,
    environment_path: Path | None,
    resume: bool = False,
) -> dict[str, Any]:
    if suite_path:
        suite = load_suite(suite_path)
    else:
        case = validate_case(_read_json(case_path), str(case_path))
        suite = {"id": case["id"], "cases": [case]}
    environment = _environment(environment_path)
    output.mkdir(parents=True, exist_ok=True)
    provenance = _provenance(environment)
    previous: dict[str, Any] = {}
    summary_path = output / "summary.json"
    if resume and summary_path.exists():
        prior = _read_json(summary_path)
        previous = {
            case["id"]: case
            for case in prior.get("cases", [])
            if isinstance(case, dict)
        }
    results = []
    for case in suite["cases"]:
        case = _expand_case(case)
        tools, _ = _tool_fingerprints(case, environment)
        fingerprint = _fingerprint(case, provenance, tools)
        cached = previous.get(case["id"])
        if (
            cached
            and cached.get("status") == "PASS"
            and cached.get("fingerprint") == fingerprint
        ):
            case_result = cached
        else:
            try:
                case_result = _run_case(case, output, environment, provenance)
            except KeyboardInterrupt:
                case_result = {
                    "id": case["id"],
                    "name": case.get("name", case["id"]),
                    "status": "ERROR",
                    "error": "interrupted",
                    "interrupted": True,
                    "steps": [],
                    "metrics": [],
                }
        results.append(case_result)
        statuses = {item["status"] for item in results}
        current_status = next(
            (
                candidate
                for candidate in ("ERROR", "BLOCKED", "FAIL")
                if candidate in statuses
            ),
            "PASS",
        )
        _write_reports(
            output,
            {
                "suite": suite["id"],
                "status": current_status,
                "provenance": provenance,
                "cases": results,
            },
        )
        if case_result.get("interrupted"):
            for pending in suite["cases"][len(results) :]:
                results.append(
                    {
                        "id": pending["id"],
                        "name": pending.get("name", pending["id"]),
                        "status": "BLOCKED",
                        "reason": "run interrupted before case started",
                        "steps": [],
                        "metrics": [],
                    }
                )
            break
    statuses = {case["status"] for case in results}
    status = next(
        (
            candidate
            for candidate in ("ERROR", "BLOCKED", "FAIL")
            if candidate in statuses
        ),
        "PASS",
    )
    summary = {
        "suite": suite["id"],
        "status": status,
        "provenance": provenance,
        "cases": results,
        "interrupted": any(case.get("interrupted") for case in results),
    }
    _write_reports(output, summary)
    return summary


def report(run_path: Path) -> dict[str, Any]:
    path = run_path / "summary.json" if run_path.is_dir() else run_path
    summary = _read_json(path)
    if "status" not in summary or "cases" not in summary:
        raise ConfigError(f"{path}: not a validation run summary")
    return summary


def brief(summary: dict[str, Any], run_path: Path) -> dict[str, Any]:
    output = run_path if run_path.is_dir() else run_path.parent
    return {
        "status": summary["status"],
        "cases": [
            {"id": case["id"], "status": case["status"]}
            for case in summary.get("cases", [])
        ],
        "summary_path": str(output / "summary.json"),
        "report_path": str(output / "report.html"),
    }
