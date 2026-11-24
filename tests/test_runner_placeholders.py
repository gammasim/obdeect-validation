import json
import sys
from pathlib import Path

import pytest

from obdeect_validation import runner

PYTHON = sys.executable


@pytest.fixture
def clean_roots(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("VALIDATION_ROOT", raising=False)
    monkeypatch.delenv("SIMTEL_REFERENCE_SOURCE", raising=False)


def test_expand_placeholders_replaces_known_names_only() -> None:
    values = {"python": "/py", "validation": "/val", "reference": "/ref"}
    assert runner.expand_placeholders("{python} -m x", values) == "/py -m x"
    assert runner.expand_placeholders("{validation}/fixtures/a", values) == (
        "/val/fixtures/a"
    )
    assert runner.expand_placeholders("{reference}/cfg/{x}", values) == ("/ref/cfg/{x}")
    assert runner.expand_placeholders("{unknown}", values) == "{unknown}"
    assert runner.expand_placeholders("/validation/x", values) == "/validation/x"


def test_defaults_are_the_container_mounts(clean_roots: None) -> None:
    environment = {"env": {}}
    values = runner.placeholder_values(environment)
    assert values["validation"] == "/validation"
    assert values["reference"] == "/reference"
    assert values["python"] == sys.executable


def test_process_environment_sets_host_roots(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("VALIDATION_ROOT", "/host/obdeect-validation")
    monkeypatch.setenv("SIMTEL_REFERENCE_SOURCE", "/host/sim_telarray")
    values = runner.placeholder_values({"env": {}})
    assert values["validation"] == "/host/obdeect-validation"
    assert values["reference"] == "/host/sim_telarray"


def test_environment_file_overrides_process_environment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("VALIDATION_ROOT", "/process")
    environment = {"env": {"VALIDATION_ROOT": "/file"}}
    assert runner.placeholder_values(environment)["validation"] == "/file"


def test_empty_root_falls_back_to_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("VALIDATION_ROOT", "")
    assert runner.placeholder_values({"env": {}})["validation"] == "/validation"


def test_expand_case_rewrites_argv_and_inputs_without_mutating_input(
    clean_roots: None,
) -> None:
    case = {
        "id": "placeholder",
        "inputs": ["{validation}/inputs/a.iact", "{reference}/cfg/b.dat"],
        "commands": [
            {
                "name": "run",
                "argv": [
                    "{python}",
                    "-m",
                    "tool",
                    "--config",
                    "{validation}/fixtures/c.json",
                ],
            }
        ],
    }
    original = json.loads(json.dumps(case))
    expanded = runner._expand_case(case, {"env": {"VALIDATION_ROOT": "/v"}})
    assert case == original
    assert expanded["inputs"] == ["/v/inputs/a.iact", "/reference/cfg/b.dat"]
    assert expanded["commands"][0]["argv"] == [
        sys.executable,
        "-m",
        "tool",
        "--config",
        "/v/fixtures/c.json",
    ]


def test_expand_case_keeps_absent_inputs_absent(clean_roots: None) -> None:
    case = {"id": "plain", "commands": [{"name": "run", "argv": ["{python}", "-V"]}]}
    expanded = runner._expand_case(case, {"env": {}})
    assert "inputs" not in expanded


def test_run_resolves_placeholders_from_environment_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("VALIDATION_ROOT", raising=False)
    root = tmp_path / "validation-root"
    (root / "inputs").mkdir(parents=True)
    shower = root / "inputs" / "shower.iact"
    shower.write_bytes(b"x")
    marker = root / "marker.txt"
    case = {
        "id": "placeholder-run",
        "inputs": ["{validation}/inputs/shower.iact"],
        "commands": [
            {
                "name": "mark",
                "argv": [
                    "{python}",
                    "-c",
                    "import pathlib,sys; pathlib.Path(sys.argv[1]).write_text('ok')",
                    "{validation}/marker.txt",
                ],
            }
        ],
        "metrics": [],
    }
    case_file = tmp_path / "case.json"
    case_file.write_text(json.dumps(case), encoding="utf-8")
    environment = tmp_path / "environment.json"
    environment.write_text(
        json.dumps({"env": {"VALIDATION_ROOT": str(root)}}), encoding="utf-8"
    )
    result = runner.run(None, case_file, tmp_path / "run", environment, resume=False)
    case_result = result["cases"][0]
    assert case_result["status"] == "PASS", case_result
    assert case_result["missing_inputs"] == []
    assert marker.read_text(encoding="utf-8") == "ok"


def test_run_blocks_when_placeholder_input_is_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("VALIDATION_ROOT", raising=False)
    case = {
        "id": "missing-placeholder",
        "inputs": ["{validation}/inputs/absent.iact"],
        "commands": [{"name": "never", "argv": [PYTHON, "-c", "pass"]}],
        "metrics": [],
    }
    case_file = tmp_path / "case.json"
    case_file.write_text(json.dumps(case), encoding="utf-8")
    environment = tmp_path / "environment.json"
    environment.write_text(
        json.dumps({"env": {"VALIDATION_ROOT": str(tmp_path / "none")}}),
        encoding="utf-8",
    )
    result = runner.run(None, case_file, tmp_path / "run", environment, resume=False)
    case_result = result["cases"][0]
    assert case_result["status"] == "BLOCKED"
    assert case_result["steps"] == []
    assert case_result["missing_inputs"] == [
        str(tmp_path / "none" / "inputs" / "absent.iact")
    ]
