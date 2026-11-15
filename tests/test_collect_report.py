import importlib.util
import json
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "collect_report.py"
SPEC = importlib.util.spec_from_file_location("collect_report", SCRIPT)
collect_report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(collect_report)


def _report(tmp_path):
    report = tmp_path / "2026-01-01"
    case = report / "sources.flasher-dish"
    case.mkdir(parents=True)
    (case / "metrics.json").write_text(
        json.dumps({"photon_ratio": 1.0015, "photon_ratio_se": 0.0008})
    )
    (case / "stage1.png").write_bytes(b"png")
    (case / "extra.png").write_bytes(b"png")
    summary = {
        "suite": "full",
        "status": "FAIL",
        "cases": [
            {
                "id": "sources.flasher-dish",
                "name": "Dish flasher versus ff-1m",
                "status": "FAIL",
                "metrics": [
                    {"metric": "event_count_delta", "status": "PASS"},
                    {"metric": "t_mean_limit_exceeded", "status": "FAIL"},
                ],
            }
        ],
    }
    (report / "summary.json").write_text(json.dumps(summary))
    return report


def test_format_row_marks_missing_and_wrong_values():
    assert (
        collect_report.format_row("{a:.2f} ± {b:.2f}", {"a": 1, "b": 0.5})
        == "1.00 ± 0.50"
    )
    assert collect_report.format_row("{a:.2f}", {}) == collect_report.MISSING
    assert collect_report.format_row("{a:.2f}", {"a": "x"}) == collect_report.MISSING


def test_headline_drops_parts_with_missing_metrics():
    entry = {"headline": "photons {a:.1f}, camera {b:.1f}"}
    assert collect_report._headline(entry, {"a": 1.0}) == "photons 1.0"
    assert collect_report._headline(entry, {}) == collect_report.MISSING


def test_readme_has_step_sections_figures_and_failed_checks(tmp_path):
    report = _report(tmp_path)
    readme = collect_report.write_readme(report).read_text()
    assert "**Result: 1 FAIL**" in readme
    assert "## 1. Light sources" in readme
    assert "[Flasher on the dish vs ff-1m](#sources-flasher-dish)" in readme
    assert "photons 1.0015 ± 0.0008" in readme
    assert "failed: `t_mean_limit_exceeded`" in readme
    assert "(sources.flasher-dish/stage1.png)" in readme
    # Figures without a caption in the layout are still shown.
    assert "Other figures:" in readme and "(sources.flasher-dish/extra.png)" in readme


def test_rewrite_keeps_hand_written_notes(tmp_path):
    report = _report(tmp_path)
    readme = collect_report.write_readme(report)
    text = readme.read_text().replace(
        "<!-- notes:sources -->\n\n", "<!-- notes:sources -->\nff-1m finding.\n"
    )
    readme.write_text(text)
    rewritten = collect_report.write_readme(report).read_text()
    assert "<!-- notes:sources -->\nff-1m finding.\n<!-- /notes -->" in rewritten


def test_case_directory_missing_from_summary_is_report_only(tmp_path):
    report = _report(tmp_path)
    extra = report / "performance.stages"
    extra.mkdir()
    (extra / "metrics.json").write_text("{}")
    readme = collect_report.write_readme(report).read_text()
    assert "## 2. Performance" in readme
    assert "REPORT ONLY" in readme


def test_reading_guide_and_differences_only_for_current_reports(tmp_path):
    report = _report(tmp_path)
    layout = collect_report._read_json(collect_report.LAYOUT)
    readme = collect_report.render_readme(report, layout, {})
    assert (
        "## How to read this report" not in readme
    )  # 2026-01-01 predates layout["since"]
    assert "predates the reading guide" in readme
    current = report.rename(tmp_path / "2099-01-01")
    metrics = {"t_mean_delta": 0.001, "photon_ratio": 1.0, "photon_ratio_se": 0.001}
    (current / "sources.flasher-dish" / "metrics.json").write_text(json.dumps(metrics))
    readme = collect_report.render_readme(current, layout, {})
    assert "## How to read this report" in readme
    assert "> **Expected differences.**" in readme and "+0.001 ns" in readme
    # Paragraphs of the guide are separated, list items are not.
    assert (
        "**What PASS means.**"
        in readme.split("**What is compared.**")[1].split("\n\n")[1]
    )


def test_links_are_listed_not_embedded(tmp_path):
    report = _report(tmp_path)
    entry = {"links": {"label": "simtools plots", "files": ["*.png"]}}
    lines = collect_report._figures(report / "sources.flasher-dish", entry, report)
    assert lines[0].startswith(
        "simtools plots: [extra.png](sources.flasher-dish/extra.png)"
    )
    assert not any(line.startswith("![") for line in lines)
