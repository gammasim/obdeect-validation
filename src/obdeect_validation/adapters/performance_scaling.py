"""Fixed and per-shower cost of the LSTN-01 air-shower chain stages, and trace-stream thread scaling (report only).

Measurement D: each stage of the chain of performance_chain.shower_candidate (threads=1, one file in
and one file out per stage, native programs) runs on the input of 0, N and 2N showers (the first
events of --input; 0 keeps the run header only). Wall and CPU times of each stage are fitted as
fixed + per_shower x showers by least squares. compile_model does not read the showers: its time is
reported as fixed (mean of the measurements, per_shower 0). fixed_share_at_N is the fixed part of the
fitted wall time at N showers.

Measurement G: trace-stream on the propagated stream of N showers (from D) with --threads 1, 2, 4, 8,
only the counts up to the host nproc. Median wall and CPU over --reps. speedup = wall(1) / wall(k),
efficiency = speedup / k.

Limits: sim_telarray stages are not measured. A stage that fails on 0 showers is recorded in details
and fitted from N and 2N only (a two-point fit is exact). A stage that fails on N or 2N is not
completed and has no fit metrics. The fit uses the actual shower counts of the subsets. There are no
budgets: all metrics are reported, *_completed = 1 when the measurement ran.

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.performance_scaling --config cfg/model/7.0.0/CTAO.cfg [--events N] [--reps N] [--input FILE]
"""

from __future__ import annotations

import argparse
import os
import shutil
import statistics
import tempfile
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import chain_optics as co
from obdeect_validation.adapters import common, performance, performance_chain
from obdeect_validation.adapters import obdeect_cli as cli

EVENTS = 10
REPS = 3
THREADS = (1, 2, 4, 8)
STAGES = ("from_corsika", "propagation", "compile_model", "trace_stream", "to_simtel")
#: Stages without a shower dependence: reported as fixed cost.
FIXED_STAGES = ("compile_model",)


def fit_line(xs, ys):
    """Least-squares (fixed, per_unit) of ys = fixed + per_unit * xs."""
    if len(set(xs)) < 2:
        raise ValueError("a fit needs at least two different shower counts")
    design = np.column_stack([np.ones(len(xs)), np.asarray(xs, dtype=float)])
    (fixed, per_unit), *_ = np.linalg.lstsq(
        design, np.asarray(ys, dtype=float), rcond=None
    )
    return float(fixed), float(per_unit)


def fixed_share(fixed, per_shower, showers):
    """Fixed part of the fitted time at `showers`; None when the fitted total is not positive."""
    total = fixed + per_shower * showers
    return fixed / total if total > 0 else None


def fit_stage(points, fixed_only=False):
    """Fit of one stage. points: (showers, wall_s, cpu_s) per input."""
    if fixed_only:
        return {
            "fixed_wall_s": statistics.fmean(p[1] for p in points),
            "per_shower_wall_s": 0.0,
            "fixed_cpu_s": statistics.fmean(p[2] for p in points),
            "per_shower_cpu_s": 0.0,
        }
    showers = [p[0] for p in points]
    fixed_wall, per_wall = fit_line(showers, [p[1] for p in points])
    fixed_cpu, per_cpu = fit_line(showers, [p[2] for p in points])
    return {
        "fixed_wall_s": fixed_wall,
        "per_shower_wall_s": per_wall,
        "fixed_cpu_s": fixed_cpu,
        "per_shower_cpu_s": per_cpu,
    }


def stage_metrics(stage, fit, showers_at_n):
    """<stage>_fixed_* / <stage>_per_shower_* metrics and the fixed share of the wall time at N."""
    metrics = {f"{stage}_{key}": value for key, value in fit.items()}
    share = fixed_share(fit["fixed_wall_s"], fit["per_shower_wall_s"], showers_at_n)
    if share is not None:
        metrics[f"{stage}_fixed_share_at_N"] = share
    return metrics


def median_of(records):
    """(median wall_s, median cpu_s) of run records."""
    return (
        statistics.median(r["wall_s"] for r in records),
        statistics.median(r["cpu_s"] for r in records),
    )


def thread_counts(nproc):
    """Thread counts of the scaling run: THREADS up to the host nproc (1 always included)."""
    return tuple(k for k in THREADS if k <= max(1, nproc))


def thread_metrics(rows, nproc):
    """Metrics of the thread scaling. rows: threads -> (median wall_s, median cpu_s)."""
    base = rows[1][0]
    metrics = {"host_nproc": nproc}
    for k, (wall, cpu) in sorted(rows.items()):
        speedup = base / wall
        metrics[f"trace_stream_threads{k}_wall_s"] = wall
        metrics[f"trace_stream_threads{k}_cpu_s"] = cpu
        metrics[f"trace_stream_threads{k}_speedup"] = speedup
        metrics[f"trace_stream_threads{k}_efficiency"] = speedup / k
    return metrics


def assemble(points, events, thread_rows, nproc):
    """Metrics and fits from the measured points.

    points: shower request -> (showers kept, point summary {stage: {wall_s, cpu_s, reps}}).
    The stage fits use every point where the stage ran; N and 2N are required.
    """
    metrics, fits = {}, {}
    needed = {events, 2 * events}
    for stage in STAGES:
        measured = [
            (kept, summary[stage]["wall_s"], summary[stage]["cpu_s"])
            for _, (kept, summary) in sorted(points.items())
            if stage in summary
        ]
        present = {
            requested for requested, (_, summary) in points.items() if stage in summary
        }
        if not needed <= present:
            metrics[f"{stage}_scaling_completed"] = 0
            continue
        fit = fit_stage(measured, fixed_only=stage in FIXED_STAGES)
        metrics.update(stage_metrics(stage, fit, points[events][0]))
        metrics[f"{stage}_scaling_completed"] = 1
        fits[stage] = {"points": measured, **fit}
    completed = set(thread_rows) == set(thread_counts(nproc))
    if completed:
        metrics.update(thread_metrics(thread_rows, nproc))
    metrics["host_nproc"] = nproc
    metrics["thread_scaling_completed"] = int(completed)
    return metrics, fits


def trace_argv(model, propagated, traced, threads):
    """trace-stream of performance_chain.shower_candidate with the given thread count."""
    return cli.trace_stream(
        model,
        propagated,
        traced,
        seed=co.SEED,
        threads=threads,
        telescopes=",".join(map(str, co.TELESCOPES)),
        pointing_zenith_deg=co.ZENITH_DEG,
        pointing_azimuth_deg=co.AZIMUTH_DEG,
        telescope_transmission="full",
    )


def chain_commands(work, iact, inputs):
    """(stage, argv, log, stdout_file) of the chain, in order. Each stage reads the previous output."""
    native, propagated = work / "showers.gst", work / "propagated.gst"
    model, traced = work / "lst.model", work / "traced.gst"
    camera_stream = work / "camera.iact"
    camera_map = co.camera_map(work, inputs["rotation"])
    return [
        (
            "from_corsika",
            cli.from_corsika(iact),
            work / "from-corsika.log",
            native,
        ),
        (
            "propagation",
            cli.propagation_run(
                native,
                propagated,
                transmission=inputs["atmosphere"],
                site_altitude=2156,
                pathlen=28,
                seed=123,
                mode="weighted",
                wavelength_range=(250, 700),
                telescopes=",".join(map(str, co.TELESCOPES)),
            ),
            work / "propagation.log",
            None,
        ),
        (
            "compile_model",
            performance_chain.compile_argv(inputs, model),
            work / "compile.log",
            None,
        ),
        (
            "trace_stream",
            trace_argv(model, propagated, traced, 1),
            work / "trace.log",
            None,
        ),
        (
            "to_simtel",
            cli.obdeect(
                "stream",
                "to-simtel",
                camera_map=camera_map,
                o=camera_stream,
                extra=(traced,),
            ),
            work / "to-simtel.log",
            None,
        ),
    ]


def run_chain_once(work, iact, inputs, records, failures):
    """One pass of the chain. Stops at the first failing stage; returns the stages that were skipped."""
    commands = chain_commands(work, iact, inputs)
    for position, (stage, argv, log, stdout) in enumerate(commands):
        try:
            records[stage].append(
                performance_chain.run_timed(argv, log, stdout_file=stdout)
            )
        except Exception as error:  # a failed stage is recorded; the chain stops here
            failures.setdefault(stage, str(error))
            return [name for name, *_ in commands[position + 1 :]]
    return []


def measure_points(runs, inputs, reps):
    """Run the chain reps times on each point; repetitions alternate between points.

    runs: name -> (work dir, iact). Rep 1 of every point, then rep 2 of every point, ...; a point
    whose stage fails stops repeating, the other points continue.
    Returns name -> (records, failures, skipped stages of its last pass).
    """
    state = {name: ({stage: [] for stage in STAGES}, {}, []) for name in runs}
    active = list(runs)
    for _ in range(reps):
        for name in active:
            work, iact = runs[name]
            records, failures, _ = state[name]
            skipped = run_chain_once(work, iact, inputs, records, failures)
            state[name] = (records, failures, skipped)
        active = [name for name in active if not state[name][1]]
    return state


def point_summary(records, failures):
    """Median, minimum wall and CPU per stage that ran and did not fail."""
    summary = {}
    for stage in STAGES:
        if records[stage] and stage not in failures:
            wall, cpu = median_of(records[stage])
            summary[stage] = {
                "wall_s": wall,
                "cpu_s": cpu,
                "wall_s_min": min(r["wall_s"] for r in records[stage]),
                "cpu_s_min": min(r["cpu_s"] for r in records[stage]),
                "reps": len(records[stage]),
            }
    return summary


def measure_threads(work, nproc, reps):
    """trace-stream on the propagated stream and model of work (N showers), per thread count.

    Repetitions alternate between thread counts: rep 1 of every count, then rep 2, ...; a count
    that fails stops repeating, the other counts continue.
    Returns (rows: threads -> (median wall_s, median cpu_s), failures: threads -> message).
    """
    rows, failures = {}, {}
    thread_dir = work / "threads"
    thread_dir.mkdir(parents=True, exist_ok=True)
    argvs = {
        k: trace_argv(
            work / "lst.model",
            work / "propagated.gst",
            thread_dir / f"traced-threads{k}.gst",
            k,
        )
        for k in thread_counts(nproc)
    }
    records = {k: [] for k in argvs}
    active = list(argvs)
    for _ in range(reps):
        for k in active:
            try:
                records[k].append(
                    performance_chain.run_timed(
                        argvs[k], thread_dir / f"threads{k}.log"
                    )
                )
            except Exception as error:
                failures[k] = str(error)
        active = [k for k in active if k not in failures]
    for k in argvs:
        if k not in failures:
            rows[k] = median_of(records[k])
    return rows, failures


def _plot(path, fits, thread_rows):
    """Left: wall time per stage against showers, with the fit lines. Right: speedup against threads."""
    from obdeect_validation.adapters import plotting

    fig, axes = plotting.figure(1, 2, (11, 3.8))
    for index, (stage, fit) in enumerate(fits.items()):
        color = plotting.CYCLE[index % len(plotting.CYCLE)]
        xs = [point[0] for point in fit["points"]]
        ys = [point[1] for point in fit["points"]]
        axes[0, 0].plot(xs, ys, "o", color=color, label=stage)
        grid = np.linspace(0, max(xs), 50)
        axes[0, 0].plot(
            grid,
            fit["fixed_wall_s"] + fit["per_shower_wall_s"] * grid,
            "-",
            color=color,
        )
    axes[0, 0].set_xlabel("showers")
    axes[0, 0].set_ylabel("wall time (s, median)")
    axes[0, 0].legend(fontsize=8)
    if 1 in thread_rows:
        ks = sorted(thread_rows)
        speedups = [thread_rows[1][0] / thread_rows[k][0] for k in ks]
        axes[0, 1].plot(ks, speedups, "o-", color=plotting.NATIVE, label="measured")
        axes[0, 1].plot(ks, ks, "--", color=plotting.LIMIT, label="ideal")
        axes[0, 1].set_xscale("log", base=2)
        axes[0, 1].set_xticks(ks, [str(k) for k in ks])
        axes[0, 1].set_xlabel("trace-stream threads")
        axes[0, 1].set_ylabel("speedup (wall at 1 thread / wall)")
        axes[0, 1].legend(fontsize=8)
    plotting.save(fig, path)


def run_scaling(source, work, config, events, reps):
    """Measurements D and G in work. Returns (metrics, details, fits, thread_rows)."""
    if performance.shower_count(source) < 2 * events:
        raise ValueError(
            f"{source} holds fewer than {2 * events} showers (2 x --events)"
        )
    inputs, _ = performance_chain.prepare_inputs(work, config)
    host = performance.host_info()
    nproc = host["nproc"]
    points, details_points, kept_by, runs = {}, {}, {}, {}
    for requested in (0, events, 2 * events):
        point = work / f"showers-{requested}"
        point.mkdir(parents=True, exist_ok=True)
        iact = point / "showers.iact"
        kept_by[requested] = performance.subset_showers(source, iact, requested)
        runs[requested] = (point, iact)
    results = measure_points(runs, inputs, reps)
    for requested, (records, failures, skipped) in results.items():
        summary = point_summary(records, failures)
        points[requested] = (kept_by[requested], summary)
        details_points[str(requested)] = {
            "showers": kept_by[requested],
            "stages": summary,
            "failures": failures,
            "skipped": skipped,
        }
    thread_rows, thread_failures = {}, {}
    if {"propagation", "compile_model"} <= set(points[events][1]):
        thread_rows, thread_failures = measure_threads(
            work / f"showers-{events}", nproc, reps
        )
    metrics, fits = assemble(points, events, thread_rows, nproc)
    details = {
        "host": host,
        "settings": {
            "events": events,
            "reps": reps,
            "input": str(source),
            "config": str(config),
            "work_dir": str(work),
            "threads_tested": list(thread_counts(nproc)),
            "stage_programs": {
                stage: performance_chain.native_argv(argv)[0]
                for stage, argv, _, _ in chain_commands(work, work / "x.iact", inputs)
            },
            "timing": "wait4 rusage of each child; wall time around the child",
            "fit": "least squares, fixed + per_shower x showers, wall and CPU separately",
        },
        "points": details_points,
        "fits": fits,
        "threads": {
            "rows": {
                str(k): {"wall_s": wall, "cpu_s": cpu}
                for k, (wall, cpu) in thread_rows.items()
            },
            "failures": {str(k): v for k, v in thread_failures.items()},
        },
    }
    return metrics, details, fits, thread_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="CTAO.cfg of the LSTN-01 configuration (cfg/model/7.0.0/CTAO.cfg)",
    )
    parser.add_argument("--input", type=Path, default=performance.INPUT)
    parser.add_argument(
        "--events",
        type=int,
        default=EVENTS,
        help="N: the stages run on 0, N and 2N showers",
    )
    parser.add_argument("--reps", type=int, default=REPS)
    args = parser.parse_args()
    if args.events < 1 or args.reps < 1:
        parser.error("--events and --reps must be at least 1")
    root = Path.cwd().resolve()
    # Stream files go to a container-local directory (PERF_WORK_DIR), not the bind-mounted results.
    env_dir = os.environ.get("PERF_WORK_DIR")
    base = Path(env_dir) if env_dir else Path(tempfile.mkdtemp(prefix="perf-scaling-"))
    work = base / "scaling-work"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True, exist_ok=True)
    try:
        metrics, details, fits, thread_rows = run_scaling(
            args.input.resolve(), work, args.config.resolve(), args.events, args.reps
        )
        _plot(root / "scaling.png", fits, thread_rows)
        common.write_results(metrics, details, root)
    finally:
        shutil.rmtree(work, ignore_errors=True)
        if not env_dir:
            shutil.rmtree(base, ignore_errors=True)


if __name__ == "__main__":
    main()
