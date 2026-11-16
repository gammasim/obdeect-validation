"""Copy the small, reviewable part of a validation run into reports/ and write its README.md.

Keeps verdicts, metrics, details and figures; drops photon files, simulator
output, logs and generated model configurations. README.md is the entry point of
a report: overview table, then one section per validation step with each case's
verdict, headline numbers and figures (layout in scripts/report_layout.json).

Hand-written text goes into blocks of README.md and is kept when it is rewritten:

    <!-- notes:intro -->        run, versions, result in words
    <!-- notes:<step key> -->   findings for one step (stream, sources, ...)
    <!-- notes:end -->          limits, scope
    <!-- /notes -->             closes each block

A collection replaces the report: case directories and run files of an earlier collection
into the same directory are removed first (README.md and its notes are kept). Cases of a
separate run (e.g. performance) are added with --add; they are listed in report_only.json
and shown as REPORT ONLY. A report-only case keeps its verdict from the run (verdict.json in its
directory: status, gated checks with actual values and limits, what is missing); a case with gated
checks, or BLOCKED / ERROR / FAIL, shows that status instead of REPORT ONLY. The README shows only
the cases of summary.json and that list.

    python scripts/collect_report.py results/full-integrated reports/2026-11-13
    python scripts/collect_report.py --add results/perf reports/2026-11-13
    python scripts/collect_report.py --readme-only reports/2026-11-13
"""

import argparse
import fnmatch
import json
import re
import shutil
import string
from pathlib import Path

KEEP_NAMES = {"metrics.json", "details.json"}
KEEP_SUFFIXES = {".png", ".svg"}
TOP_LEVEL = ("summary.json", "junit.xml", "report.html")
REPORT_ONLY = "report_only.json"
VERDICT = "verdict.json"
VERDICT_KEYS = (
    "name",
    "status",
    "metrics",
    "missing",
    "missing_inputs",
    "missing_modules",
    "reason",
    "error",
    "metric_error",
)
MAX_BYTES = 2_000_000
LAYOUT = Path(__file__).with_name("report_layout.json")
NOTES = re.compile(r"<!-- notes:([\w-]+) -->\n(.*?)<!-- /notes -->", re.DOTALL)
MISSING = "–"


def collect(run: Path, report: Path, add: bool = False) -> list[Path]:
    """Copy a run into a report. Without add, the report is replaced (stale case directories
    and run files removed, README.md kept); with add, the run's cases are added as report-only
    cases and recorded in report_only.json."""
    copied = []
    report.mkdir(parents=True, exist_ok=True)
    if add:
        run_cases = {
            case["id"]: case
            for case in _read_json(run / "summary.json").get("cases", [])
            if isinstance(case, dict) and "id" in case
        }
        listed = set(_read_json(report / REPORT_ONLY).get("cases", []))
        listed.update(path.parent.name for path in run.glob("*/metrics.json"))
        listed.update(run_cases)
        (report / REPORT_ONLY).write_text(
            json.dumps({"cases": sorted(listed)}, indent=2) + "\n", encoding="utf-8"
        )
        for case_id, case in run_cases.items():
            target = report / case_id
            target.mkdir(parents=True, exist_ok=True)
            verdict = {key: case[key] for key in VERDICT_KEYS if key in case}
            (target / VERDICT).write_text(
                json.dumps(verdict, indent=2) + "\n", encoding="utf-8"
            )
            copied.append(target / VERDICT)
    else:
        for path in report.iterdir():
            if path.is_dir():
                shutil.rmtree(path)
            elif path.name in (*TOP_LEVEL, REPORT_ONLY):
                path.unlink()
        for name in TOP_LEVEL:
            if (run / name).is_file():
                shutil.copy2(run / name, report / name)
                copied.append(report / name)
    for path in sorted(run.rglob("*")):
        relative = path.relative_to(run)
        if (
            not path.is_file()
            or "cfg" in relative.parts
            or path.stat().st_size > MAX_BYTES
        ):
            continue
        if path.name in KEEP_NAMES or path.suffix in KEEP_SUFFIXES:
            target = report / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
            copied.append(target)
    return copied


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}


def _case_layout(layout: dict, case_id: str) -> dict:
    cases = layout["cases"]
    if case_id in cases:
        return cases[case_id]
    for pattern, entry in cases.items():
        if fnmatch.fnmatchcase(case_id, pattern):
            return entry
    return {}


def _step_of(layout: dict, case_id: str) -> dict:
    """Step of a case: a step listing the case in 'cases' (patterns) wins over id prefixes."""
    for step in layout["steps"]:
        if any(fnmatch.fnmatchcase(case_id, p) for p in step.get("cases", [])):
            return step
    for step in layout["steps"]:
        if "prefix" in step and case_id.startswith(step["prefix"]):
            return step
    return {"key": "other", "title": "Other cases", "intro": ""}


def format_row(template: str, values: dict) -> str:
    """Fill a format string with case values (metrics, and details fields such as
    {showers[0][zenith_deg]}); a missing or wrong-typed value gives '–'."""
    fields = [name for _, name, _, _ in string.Formatter().parse(template) if name]
    if any(values.get(re.split(r"[.\[]", name)[0]) is None for name in fields):
        return MISSING
    try:
        return template.format(**values)
    except ValueError, TypeError, KeyError, IndexError:
        return MISSING


def _values(case_dir: Path) -> dict:
    """Values for the layout templates: details.json fields, overridden by metrics.json."""
    return {
        **_read_json(case_dir / "details.json"),
        **_read_json(case_dir / "metrics.json"),
    }


def _anchor(case_id: str) -> str:
    return case_id.replace(".", "-")


def _figures(case_dir: Path, entry: dict, report: Path) -> list[str]:
    """Markdown of the case figures: listed ones with captions, files listed under
    'links' as links only, then any other figure."""
    lines, shown = [], set()

    def link(path: Path) -> str:
        return path.relative_to(report).as_posix()

    for figure in entry.get("figures", []):
        if (case_dir / figure["file"]).is_file():
            path = case_dir / figure["file"]
            shown.add(path)
            lines += [
                f"![{figure['caption']}]({link(path)})",
                "",
                f"*{figure['caption']}*",
                "",
            ]
    if "links" in entry:
        linked = sorted(
            {
                path
                for pattern in entry["links"]["files"]
                for path in case_dir.glob(pattern)
            }
        )
        shown.update(linked)
        if linked:
            names = " · ".join(
                f"[{p.relative_to(case_dir)}]({link(p)})" for p in linked
            )
            lines += [f"{entry['links']['label']}: {names}", ""]
    others = [
        path
        for path in sorted(case_dir.rglob("*"))
        if path.suffix in KEEP_SUFFIXES and path not in shown
    ]
    if others:
        lines += ["Other figures:", ""]
        lines += [f"![{path.stem}]({link(path)})" for path in others]
        lines.append("")
    return lines


def _checks(case: dict) -> str:
    metrics = case.get("metrics", [])
    passed = sum(1 for metric in metrics if metric.get("status") == "PASS")
    failed = [metric["metric"] for metric in metrics if metric.get("status") != "PASS"]
    text = f"{passed}/{len(metrics)} gated checks pass"
    if failed:
        text += "; failed: " + ", ".join(f"`{name}`" for name in failed)
    if case.get("missing") or case.get("missing_inputs") or case.get("missing_modules"):
        missing = (
            case.get("missing", [])
            + case.get("missing_inputs", [])
            + case.get("missing_modules", [])
        )
        text += "; missing: " + ", ".join(f"`{name}`" for name in missing)
    return text


def _notes(notes: dict, key: str) -> list[str]:
    return [
        f"<!-- notes:{key} -->",
        notes.get(key, "").rstrip("\n"),
        "<!-- /notes -->",
        "",
    ]


def _report_only_case(case_id: str, verdict: dict) -> dict:
    """A report-only case from its run verdict: its status, or REPORT ONLY when it has no gated
    checks and completed (blocked, failed and errored cases always show their status)."""
    gated = bool(verdict.get("metrics"))
    status = verdict.get("status", "UNKNOWN")
    if not gated and status not in ("BLOCKED", "ERROR", "FAIL"):
        status = "REPORT ONLY"
    return {**verdict, "id": case_id, "status": status, "report_only": True}


def _gated_table(case: dict) -> list[str]:
    """Markdown table of the gated checks of a report-only case: actual value and limit."""
    lines = ["| Gated check | Actual | Limit | Verdict |", "|---|---|---|---|"]
    for check in case["metrics"]:
        actual = check.get("actual")
        actual = (
            MISSING
            if actual is None
            else f"{actual:.4g}"
            if isinstance(actual, float)
            else str(actual)
        )
        limit = f"{check['expected']:g} ± {check['atol']:g}"
        lines.append(
            f"| `{check['metric']}` | {actual} | {limit} | {check.get('status', MISSING)} |"
        )
    return lines + [""]


def _headline(entry: dict, metrics: dict) -> str:
    """Comma-separated parts of the headline; parts with a missing metric are left out."""
    if "headline" in entry:
        parts = [format_row(part, metrics) for part in entry["headline"].split(", ")]
        return ", ".join(part for part in parts if part != MISSING) or MISSING
    rows = entry.get("rows", [])
    if not rows:
        return MISSING
    value = format_row(rows[0]["text"], metrics)
    return f"{rows[0]['label']}: {value}" if value != MISSING else MISSING


def _gated_headline(case: dict) -> str:
    """Headline of a case without a layout entry: its gated checks that pass, in words."""
    metrics = case.get("metrics", [])
    if not metrics:
        return MISSING
    passed = sum(1 for metric in metrics if metric.get("status") == "PASS")
    return f"{passed}/{len(metrics)} gated checks pass"


def _title(entry: dict, case: dict, step: dict) -> str:
    """Layout title; {suffix} is the case id without its step prefix (e.g. lst-7.0.0)."""
    if "title" in entry:
        return entry["title"].format(
            suffix=case["id"].removeprefix(step.get("prefix", ""))
        )
    return case.get("name") or case["id"]


def render_readme(report: Path, layout: dict, notes: dict) -> str:
    """Return README.md of a report directory from summary.json and the case directories."""
    summary = _read_json(report / "summary.json")
    cases = list(summary.get("cases", []))
    # Case directories from a separate run (e.g. performance) are shown as report only.
    # Cases of a separate run added with --add (report_only.json); other directories are ignored.
    known = {case["id"] for case in cases}
    cases += [
        _report_only_case(case_id, _read_json(report / case_id / VERDICT))
        for case_id in _read_json(report / REPORT_ONLY).get("cases", [])
        if case_id not in known and (report / case_id).is_dir()
    ]
    statuses = [case.get("status", "UNKNOWN") for case in cases]
    counts = ", ".join(
        f"{statuses.count(status)} {status}" for status in sorted(set(statuses))
    )
    lines = [f"# {layout['title']}, {report.name}", ""]
    lines += _notes(notes, "intro")
    result = (
        f"**Result: {counts or 'no cases'}** (suite `{summary.get('suite', MISSING)}`, "
        f"status {summary.get('status', MISSING)}). Each case below shows what is compared, "
        "its headline numbers and figures; all metrics are in the case's `metrics.json`, "
        "limits and settings in `details.json`."
    )
    lines += [result, ""]
    steps = []
    for case in cases:
        step = _step_of(layout, case["id"])
        if step not in steps:
            steps.append(step)
    order = [step["key"] for step in layout["steps"]]
    steps.sort(
        key=lambda step: (
            order.index(step["key"]) if step["key"] in order else len(order)
        )
    )

    lines += ["| Step | Case | Headline | Verdict |", "|---|---|---|---|"]
    for number, step in enumerate(steps, start=1):
        for case in cases:
            if _step_of(layout, case["id"]) != step:
                continue
            metrics = _values(report / case["id"])
            entry = _case_layout(layout, case["id"])
            lines.append(
                f"| {number}. {step['title']} | [{_title(entry, case, step)}](#{_anchor(case['id'])}) "
                f"| {_headline(entry, metrics) if entry else _gated_headline(case)} "
                f"| {case.get('status', MISSING)} |"
            )
    lines.append("")

    # The reading guide and the expected-differences texts describe the figures and the cases
    # from layout["since"] on; older reports get a short note instead.
    current = report.name >= layout.get("since", "")
    if current and layout.get("reading_guide"):
        lines += ["## How to read this report", ""]
        guide = layout["reading_guide"]
        for item, following in zip(guide, [*guide[1:], ""], strict=True):
            # Paragraphs are separated by a blank line; list items stay together.
            in_list = item.startswith("- ") and following.startswith("- ")
            lines += [item] if in_list else [item, ""]
    elif layout.get("before_since"):
        lines += [layout["before_since"], ""]

    for number, step in enumerate(steps, start=1):
        lines += [f"## {number}. {step['title']}", "", step.get("intro", ""), ""]
        lines += _notes(notes, step["key"])
        for case in cases:
            if _step_of(layout, case["id"]) != step:
                continue
            case_dir = report / case["id"]
            metrics = _values(case_dir)
            entry = _case_layout(layout, case["id"])
            lines += [
                f'<a id="{_anchor(case["id"])}"></a>',
                "",
                f"### {_title(entry, case, step)} — {case.get('status', MISSING)}",
                "",
            ]
            if case.get("name") and case["name"] != _title(entry, case, step):
                lines += [f"Case: {case['name']}", ""]
            if entry.get("description"):
                lines += [entry["description"], ""]
            if entry.get("rows"):
                lines += ["| Quantity | Value |", "|---|---|"]
                lines += [
                    f"| {row['label']} | {format_row(row['text'], metrics)} |"
                    for row in entry["rows"]
                ]
                lines.append("")
            if case.get("report_only") and case.get("metrics"):
                lines += _gated_table(case)
            files = [
                f"[{name}]({case['id']}/{name})"
                for name in ("metrics.json", "details.json")
                if (case_dir / name).is_file()
            ]
            lines += [f"`{case['id']}`: {_checks(case)}. " + " · ".join(files), ""]
            lines += _figures(case_dir, entry, report)
            if current and entry.get("differences"):
                text = format_row(entry["differences"], metrics)
                if text != MISSING:
                    lines += [f"> **Expected differences.** {text}", ""]
    lines += ["## Scope and limits", ""]
    lines += _notes(notes, "end")
    return "\n".join(lines).rstrip("\n") + "\n"


def read_notes(readme: Path) -> dict:
    """Hand-written blocks of an existing README.md, by key."""
    if not readme.is_file():
        return {}
    return {
        key: text for key, text in NOTES.findall(readme.read_text(encoding="utf-8"))
    }


def write_readme(report: Path, layout_path: Path = LAYOUT) -> Path:
    layout = _read_json(layout_path)
    readme = report / "README.md"
    readme.write_text(
        render_readme(report, layout, read_notes(readme)), encoding="utf-8"
    )
    return readme


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("paths", nargs="+", type=Path, help="RUN REPORT, or REPORT")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--add", action="store_true", help="add RUN as report-only cases")
    mode.add_argument(
        "--readme-only", action="store_true", help="rewrite REPORT/README.md"
    )
    args = parser.parse_args()
    if args.readme_only:
        print(write_readme(args.paths[0]))
        return
    if len(args.paths) != 2:
        parser.error("expected RUN REPORT")
    run, report = args.paths
    copied = collect(run, report, add=args.add)
    write_readme(report)
    total = sum(path.stat().st_size for path in copied)
    print(json.dumps({"files": len(copied), "bytes": total}))


if __name__ == "__main__":
    main()
