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


def _run(tmp_path, name, case_ids):
    run = tmp_path / name
    for case_id in case_ids:
        (run / case_id).mkdir(parents=True)
        (run / case_id / "metrics.json").write_text("{}")
    summary = {"suite": "s", "status": "PASS", "cases": []}
    summary["cases"] = [{"id": c, "status": "PASS", "metrics": []} for c in case_ids]
    (run / "summary.json").write_text(json.dumps(summary))
    return run


def test_unlisted_case_directory_is_not_shown(tmp_path):
    report = _report(tmp_path)
    (report / "stale.case").mkdir()
    (report / "stale.case" / "metrics.json").write_text("{}")
    assert "stale.case" not in collect_report.write_readme(report).read_text()


def test_collect_replaces_stale_cases_and_keeps_notes(tmp_path):
    report = tmp_path / "report"
    collect_report.collect(_run(tmp_path, "old", ["sources.old", "stream.a"]), report)
    (report / "README.md").write_text("<!-- notes:intro -->\nkeep\n<!-- /notes -->\n")
    collect_report.collect(_run(tmp_path, "new", ["stream.a"]), report)
    assert not (report / "sources.old").exists()
    readme = collect_report.write_readme(report).read_text()
    assert "sources.old" not in readme and "keep" in readme


def test_add_records_report_only_cases(tmp_path):
    report = tmp_path / "report"
    collect_report.collect(_run(tmp_path, "full", ["stream.a"]), report)
    collect_report.collect(
        _run(tmp_path, "perf", ["performance.stages"]), report, add=True
    )
    assert (report / "stream.a").is_dir()  # adding keeps the main run
    listed = json.loads((report / collect_report.REPORT_ONLY).read_text())
    assert listed == {"cases": ["performance.stages"]}
    readme = collect_report.write_readme(report).read_text()
    assert "## 2. Performance" in readme and "REPORT ONLY" in readme
    assert json.loads((report / "summary.json").read_text())["suite"] == "s"


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


def test_curvature_cases_have_their_own_step_and_detail_values(tmp_path):
    report = tmp_path / "2099-01-01"
    shower = {
        "pe_ratio": {
            "measured_ratio": 1.25,
            "predicted_ratio": 1.24,
            "ratio_error": 0.01,
        }
    }
    cases = {
        "propagation.extinction-analytic": ({"sample_count": 16}, {}),
        "propagation.curvature-chapman": ({"sample_count": 30}, {}),
        "chain.curvature-high-zenith": (
            {"file_count": 2, "pe_ratio_vs_prediction_limit_exceeded": 0},
            {"showers": [shower, shower]},
        ),
    }
    for case_id, (metrics, details) in cases.items():
        (report / case_id).mkdir(parents=True)
        (report / case_id / "metrics.json").write_text(json.dumps(metrics))
        (report / case_id / "details.json").write_text(json.dumps(details))
    summary = {"suite": "s", "status": "PASS", "cases": []}
    summary["cases"] = [{"id": c, "status": "PASS", "metrics": []} for c in cases]
    (report / "summary.json").write_text(json.dumps(summary))
    readme = collect_report.write_readme(report).read_text()
    atmosphere = readme.split("## 1. Atmosphere (propagation)")[1].split("## 2.")[0]
    curvature = readme.split("## 2. Earth curvature (large zenith angles)")[1]
    assert "extinction-analytic" in atmosphere and "curvature-chapman" not in atmosphere
    assert "curvature-chapman" in curvature and "curvature-high-zenith" in curvature
    assert "65°: p.e. ratio 1.2500 vs predicted 1.2400" in readme
    # A details field that is absent gives the missing mark, not an error.
    assert collect_report.format_row("{showers[5][x]}", {"showers": []}) == "–"


def test_table_fills_row_keys_and_marks_missing_cells():
    table = {
        "row_label": "Stage",
        "rows": [
            {"key": "propagation", "label": "propagation"},
            {"key": "compile_model", "label": "compile model"},
            {"key": "to_simtel", "label": "to-simtel"},
        ],
        "columns": [
            {"label": "Wall s", "metric": "{row}_wall_s_median", "format": ".2f"},
            {"label": "Share", "metric": "{row}_cpu_share", "format": ".0%"},
            {"label": "Profiled", "metric": "{row}_profiled"},
        ],
    }
    metrics = {
        "propagation_wall_s_median": 2.2372,
        "propagation_cpu_share": 0.25,
        "propagation_profiled": True,
        "to_simtel_wall_s_median": "x",
    }
    assert collect_report._table(table, metrics) == [
        "| Stage | Wall s | Share | Profiled |",
        "|---|---|---|---|",
        "| propagation | 2.24 | 25% | True |",
        "",
    ]
    # Missing or wrong-typed values give the missing mark; a row with no value is left out.
    assert collect_report._table(table, {"to_simtel_cpu_share": 0.5}) == [
        "| Stage | Share |",
        "|---|---|",
        "| to-simtel | 50% |",
        "",
    ]  # columns without any value are left out
    assert collect_report._table(table, {}) == []
    assert collect_report._table(table, {"to_simtel_wall_s_median": "x"}) == []


def test_performance_stage_table_in_readme(tmp_path):
    report = tmp_path / "2099-01-01"
    case = report / "performance.stages"
    case.mkdir(parents=True)
    metrics = {
        "trace_stream_wall_s_median": 12.1,
        "trace_stream_cpu_s_median": 10.31,
        "trace_stream_cpu_share": 0.8,
        "simtel_direct_wall_s_median": 11.68,
        "trace_stream_ns_per_photon_cpu": 103.3,
        "simtel_direct_ns_per_photon_cpu": 73.4,
        "timing_unreliable": False,
    }
    (case / "metrics.json").write_text(json.dumps(metrics))
    summary = {"suite": "performance", "status": "PASS", "cases": []}
    summary["cases"] = [{"id": "performance.stages", "status": "REPORT ONLY"}]
    (report / "summary.json").write_text(json.dumps(summary))
    readme = collect_report.write_readme(report).read_text()
    # Columns without any value (photons/s, spread) are left out.
    assert "| Program | Wall s | CPU s | CPU ns per photon |" in readme
    assert (
        "| ray-tracing: trace-stream (telescope optics) | 12.10 | 10.31 | 103 |"
        in readme
    )
    assert "| propagation: run (atmosphere) |" not in readme  # no values, so no row
    assert "trace-stream 103 ns/photon, sim_telarray 73 ns/photon (CPU)" in readme


def test_add_renders_gated_budget_checks_and_blocked_cases(tmp_path):
    run = tmp_path / "perf"
    case = run / "performance.air-shower"
    case.mkdir(parents=True)
    (case / "metrics.json").write_text(json.dumps({"ratio_wall_threads1": 2.5}))
    checks = [
        {
            "metric": "reference_completed",
            "expected": 1,
            "atol": 0,
            "actual": 1,
            "status": "PASS",
        },
        {
            "metric": "candidate_envelope_wall_budget_exceeded",
            "expected": 0,
            "atol": 0,
            "actual": 1,
            "status": "FAIL",
        },
    ]
    summary = {
        "suite": "performance",
        "status": "FAIL",
        "cases": [
            {"id": "performance.air-shower", "status": "FAIL", "metrics": checks},
            {
                "id": "performance.air-shower-proton",
                "status": "BLOCKED",
                "missing_inputs": ["/validation/inputs/proton.iact"],
            },
        ],
    }
    (run / "summary.json").write_text(json.dumps(summary))
    report = tmp_path / "2099-01-01"
    report.mkdir()
    collect_report.collect(run, report, add=True)
    readme = collect_report.write_readme(report).read_text()
    # Checks in words; the bookkeeping checks (runs completed) in one line.
    label = "obdeect chain wall ≤ sim_telarray wall (trace-stream on all cores)"
    assert f"| {label} | 1 | FAIL |" in readme
    assert "All runs completed with their repetitions (1 checks)." in readme
    assert "air-shower-proton" in readme
    assert "BLOCKED" in readme and "/validation/inputs/proton.iact" in readme
    assert "REPORT ONLY" not in readme.split("## 1.")[0]


def test_performance_is_a_separate_section_at_the_bottom(tmp_path):
    run = _run(tmp_path, "full", ["stream.a", "sources.b"])
    report = tmp_path / "2099-01-01"
    collect_report.collect(run, report)
    perf = _run(
        tmp_path,
        "perf",
        ["performance.stages", "performance.profile", "performance.air-shower"],
    )
    (perf / "performance.air-shower" / "metrics.json").write_text(
        json.dumps(
            {"reference_wall_s": 20.0, "reference_cpu_s": 19.0, "timing_unreliable": 0}
        )
    )
    collect_report.collect(perf, report, add=True)
    readme = collect_report.write_readme(report).read_text()
    top, bottom = readme.split("## Scope and limits")
    # Validation overview and result count: no performance cases.
    assert (
        "performance." not in top.split("## How to read")[0].split("[3. Performance")[0]
    )
    assert "**Result: 2 PASS**" in top
    title = "3. Performance: obdeect chain against sim_telarray"
    assert f"[{title}](#step-performance) is at the bottom: 3 REPORT ONLY." in top
    assert "## 3. Performance" not in top and '<a id="step-performance">' in bottom
    # Own overview table, then the cases in the layout's case_order.
    assert "| Case | Headline | Verdict |" in bottom
    order = [
        bottom.index(f'<a id="performance-{name}"></a>')
        for name in ("air-shower", "stages", "profile")
    ]
    assert order == sorted(order)
    # Quantity rows without a value are left out in the performance section.
    assert "| sim_telarray: wall s, CPU s | 20.00, 19.00 |" in bottom
    assert "| Reference CPU s |" not in bottom


def test_performance_section_shows_the_host_of_its_first_case(tmp_path):
    run = _run(tmp_path, "perf", ["performance.stages", "performance.profile"])
    (run / "performance.profile" / "details.json").write_text(
        json.dumps(
            {"host": {"host": "Apple M3 Pro (6 performance + 6 efficiency cores)"}}
        )
    )
    (run / "performance.stages" / "details.json").write_text(
        json.dumps({"host_info": {"nproc": 12}})
    )
    report = tmp_path / "2099-01-01"
    collect_report.collect(run, report)
    readme = collect_report.write_readme(report).read_text()
    assert "Host: Apple M3 Pro (6 performance + 6 efficiency cores)" in readme
    assert readme.count("Host: ") == 1


def test_no_host_line_without_host_value(tmp_path):
    run = _run(tmp_path, "perf", ["performance.stages"])
    (run / "performance.stages" / "details.json").write_text(
        json.dumps({"host_info": {"nproc": 12, "host": None}})
    )
    report = tmp_path / "2099-01-01"
    collect_report.collect(run, report)
    assert "Host: " not in collect_report.write_readme(report).read_text()
