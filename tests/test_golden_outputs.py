"""Unit tests of the golden-output adapter: comparison rules, platform key, case file consistency."""

import json
import sys
from pathlib import Path

import pytest

from obdeect_validation.adapters import golden_outputs as golden

CASE = Path(__file__).resolve().parents[1] / "cases/determinism/golden-outputs.json"


def entry(sha, pipeline="p"):
    return {"sha256": sha, "rows": None, "sum_weights": None, "pipeline": pipeline}


def test_platform_key_format():
    key = golden.platform_key()
    assert key.startswith(f"{sys.platform}-")
    assert key.endswith(golden.platform.machine())


def test_compare_identical_outputs_gives_no_differences():
    base = {"a.pst": entry("1"), "b.json": entry("2")}
    details, mismatches, missing = golden.compare(base, dict(base))
    assert (mismatches, missing) == (0, 0)
    assert {d["verdict"] for d in details.values()} == {"match"}


def test_compare_counts_changed_missing_and_unexpected_outputs():
    base = {"a.pst": entry("1"), "b.pst": entry("2"), "c.pst": entry("3")}
    current = {"a.pst": entry("X"), "c.pst": entry("3"), "d.pst": entry("4")}
    details, mismatches, missing = golden.compare(base, current)
    assert mismatches == 1  # a.pst changed bit
    assert missing == 2  # b.pst not produced, d.pst not in baseline
    assert details["a.pst"]["verdict"] == "mismatch"
    assert details["b.pst"]["verdict"] == "missing"
    assert details["d.pst"]["verdict"] == "unexpected"


def test_stage_placeholders_expand_models_only():
    argv = ["obdeect", "ray-tracing", "compile", "--simulation-models", "{models}"]
    assert golden._expand(argv, "/m")[-1] == "/m"
    assert golden._expand(argv, "/m")[0] == "obdeect"


def test_pipeline_dependencies_exist():
    names = {p.name for p in golden.PIPELINES}
    for pipeline in golden.PIPELINES:
        assert set(pipeline.depends_on) <= names


def test_case_file_matches_adapter_metrics():
    case = json.loads(CASE.read_text(encoding="utf-8"))
    assert case["id"] == "determinism.golden-outputs"
    gated = {m["metric"]: (m["expected"], m["atol"]) for m in case["metrics"]}
    assert gated == {
        "outputs_mismatch": (0, 0),
        "outputs_missing": (0, 0),
        "pipeline_failures": (0, 0),
        "baseline_missing": (0, 0),
    }
    assert case["commands"][0]["argv"][-1] == "run"


def test_fixture_eventio_matches_recorded_sha256():
    path = golden.GOLDEN_DIR / "gamma_100gev_1216.eventio"
    if not path.is_file():
        pytest.skip("fixture not available")
    assert golden.common.sha256(path) == golden.EVENTIO_SHA256


def test_baseline_file_for_this_platform_has_all_outputs():
    path = golden.GOLDEN_DIR / f"{golden.platform_key()}.json"
    if not path.is_file():
        pytest.skip(f"no baseline for {golden.platform_key()}")
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["platform"] == golden.platform_key()
    assert len(data["outputs"]) == 23
