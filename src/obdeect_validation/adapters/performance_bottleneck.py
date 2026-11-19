"""Bottleneck of the production pipe: the member that takes the wall time of the envelope pipe.

The envelope pipe (performance_chain.envelope_pipeline) runs from-corsika | propagation |
trace-stream | to-simtel | sim_telarray (Bypass_Optics=3) for LSTN-01, with the inputs of the
air-shower envelope candidate. Within one pipe every member ends at about the pipe wall, so a
member's own wall says nothing. The first k members therefore run as prefix pipes (k = 1..5; the
last member of each prefix writes to /dev/null, sim_telarray writes its output file as usual). The
marginal wall of member k is prefix k wall minus prefix k-1 wall; the member with the largest
marginal wall is the bottleneck. Member CPU is taken from the full pipe.

The model is compiled once before the measurement (timed and recorded in details.json, not in the
metrics). Repetitions are the median. Metrics are report only (no budgets).

Limits: a marginal wall is a difference of medians of separate runs, so it can be negative within
noise (the sign is kept). Adding a downstream member can also block an upstream one through the
pipe, which the marginal wall then includes. trace-stream uses --threads (default 1; 0 = every core).

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.performance_bottleneck --config cfg/model/7.0.0/CTAO.cfg [--input SHOWERS.iact] [--events N] [--reps N] [--threads N]
"""

from __future__ import annotations

import argparse
import os
import shutil
import statistics
import tempfile
from pathlib import Path

from obdeect_validation.adapters import chain_optics as co
from obdeect_validation.adapters import common, performance, performance_chain

REPS = 3
THREADS = 1


def prefixes(stages):
    """The prefix pipes: the first k members of stages, for k = 1..len(stages)."""
    return [stages[:k] for k in range(1, len(stages) + 1)]


def run_prefixes(stages, reps, work, failures):
    """Run every prefix pipe reps times (repetitions outermost); pipe records per prefix length k.

    A prefix that raises is listed in failures (k -> messages) and its repetition is left out.
    """
    records = {k: [] for k in range(1, len(stages) + 1)}
    for _ in range(reps):
        for k, part in enumerate(prefixes(stages), start=1):
            logs = {name: work / f"prefix{k}-{name}.log" for name, _ in part}
            try:
                records[k].append(
                    performance_chain.run_pipeline(part, os.devnull, logs)
                )
            except (RuntimeError, OSError) as error:  # reported, not hidden
                failures.setdefault(k, []).append(str(error)[:400])
    return records


def median_or_none(values):
    values = list(values)
    return statistics.median(values) if values else None


def member_cpu(record, name):
    """CPU time of one member in a pipe record."""
    return next(child["cpu_s"] for child in record["children"] if child["name"] == name)


def bottleneck_metrics(names, records, reps):
    """Prefix walls and CPU, marginal walls, wall shares, bottleneck index and full-pipe spread.

    names: member names in pipe order. records: {k: pipe records of prefix k}. A prefix without
    all of its repetitions gives None for its wall and for the marginals that use it.
    """
    count = len(names)
    full = records[count]
    walls = {0: 0.0}
    metrics = {}
    for k in range(1, count + 1):
        done = records[k]
        metrics[f"prefix{k}_completed"] = int(len(done) == reps)
        walls[k] = median_or_none(r["wall_s"] for r in done)
        metrics[f"prefix{k}_wall_s"] = walls[k]
        metrics[f"prefix{k}_cpu_s"] = median_or_none(r["cpu_s"] for r in done)
    total = walls[count]
    marginal = {}
    for k, name in enumerate(names, start=1):
        marginal[k] = (
            None
            if walls[k] is None or walls[k - 1] is None
            else walls[k] - walls[k - 1]
        )
        metrics[f"bottleneck_{name}_marginal_wall_s"] = marginal[k]
        metrics[f"bottleneck_{name}_wall_share"] = (
            marginal[k] / total if marginal[k] is not None and total else None
        )
        metrics[f"bottleneck_{name}_cpu_s"] = median_or_none(
            member_cpu(r, name) for r in full
        )
    known = {k: value for k, value in marginal.items() if value is not None}
    metrics["bottleneck_stage_index"] = max(known, key=known.get) if known else None
    full_walls = [r["wall_s"] for r in full]
    metrics["bottleneck_wall_s_spread"] = (
        (max(full_walls) - min(full_walls)) / statistics.median(full_walls)
        if full_walls
        else None
    )
    return metrics


def _plot(path, names, metrics):
    """Marginal wall per member (bar) and member CPU in the full pipe (lighter bar), pipe order."""
    from obdeect_validation.adapters import plotting

    fig, axis = plotting.figure(1, 1, (9, 1.2 + 0.6 * len(names)))
    rows = range(len(names))
    colors = [plotting.CYCLE[i % len(plotting.CYCLE)] for i in rows]
    wall = [metrics[f"bottleneck_{name}_marginal_wall_s"] or 0.0 for name in names]
    cpu = [metrics[f"bottleneck_{name}_cpu_s"] or 0.0 for name in names]
    axis.barh(
        [r - 0.2 for r in rows],
        wall,
        height=0.4,
        color=colors,
        label="marginal wall",
    )
    axis.barh(
        [r + 0.2 for r in rows],
        cpu,
        height=0.4,
        color=colors,
        alpha=0.45,
        label="member CPU (full pipe)",
    )
    axis.set_yticks(list(rows), names)
    axis.invert_yaxis()
    axis.set_xlabel("seconds (median; marginal = prefix k wall minus prefix k-1 wall)")
    axis.grid(axis="y", visible=False)
    axis.legend(fontsize=8)
    plotting.save(fig, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", type=Path, required=True, help="sim_telarray config from simtools"
    )
    parser.add_argument("--input", type=Path, default=performance.INPUT)
    parser.add_argument(
        "--events", type=int, default=None, help="first N showers (default all)"
    )
    parser.add_argument("--reps", type=int, default=REPS)
    parser.add_argument(
        "--threads",
        type=int,
        default=THREADS,
        help="trace-stream threads (0 = every core)",
    )
    args = parser.parse_args()
    root = Path.cwd().resolve()
    work_dir = os.environ.get("PERF_WORK_DIR")
    work = Path(work_dir or tempfile.mkdtemp(prefix="perf-bottleneck-"))
    work.mkdir(parents=True, exist_ok=True)
    try:
        source = performance.select_input(
            args.input.resolve(), args.events, work / "showers-input.iact"
        )
        inputs, _ = performance_chain.prepare_inputs(work, args.config)
        headroom = performance_chain.headroom_for(inputs["qe_variation"])
        files = performance_chain.envelope_files(
            work / "envelope", inputs["qe"], headroom
        )
        sim_argv, _ = performance_chain.simtel_candidate(work, inputs, files)
        model = work / "lst.model"
        stages = performance_chain.envelope_pipeline(
            source,
            inputs["atmosphere"],
            files["envelope"],
            model,
            co.camera_map(work, inputs["rotation"]),
            sim_argv,
            args.threads,
        )
        names = [name for name, _ in stages]
        compile_record = performance_chain.run_pipeline(
            [("compile_model", performance_chain.compile_argv(inputs, model))],
            os.devnull,
            {"compile_model": work / "compile.log"},
        )
        failures: dict = {}
        records = run_prefixes(stages, args.reps, work, failures)
        metrics = bottleneck_metrics(names, records, args.reps)
        index = metrics["bottleneck_stage_index"]
        details = {
            "stages": names,
            "bottleneck_stage": names[index - 1] if index else None,
            "records": records,
            "compile": compile_record,
            "failures": failures,
            "host_info": performance.host_info(),
            "settings": {
                "reps": args.reps,
                "threads": args.threads,
                "input": str(source),
                "events": performance.shower_count(source),
                "config": str(inputs["config"]),
                "candidate_atmosphere": "obdeect propagation (reference atmosphere below 1000 nm); sim_telarray unit table",
                "headroom": headroom,
                "work_dir": str(work),
                "timing": "pipe wall from start to last exit; CPU per member (wait4 rusage)",
            },
        }
        _plot(root / "bottleneck.png", names, metrics)
        common.write_results(metrics, details, root)
    finally:
        if work_dir is None:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
