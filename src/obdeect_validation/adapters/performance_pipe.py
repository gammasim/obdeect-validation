"""Stream pipe throughput: cat | tee and cat | obdeect stream to-simtel, budget >= 1 GB/s.

The native photon stream of the CORSIKA input (made once, not timed) is read from a file by cat and
piped into the member under test. Throughput = bytes of the native stream / pipe wall time (decimal
GB, 1e9 bytes/s); the pipe wall is the time from the start of the pipe to the last exit. The cat | tee
pipe is the floor of the pipe itself (tee writes a copy to disk, stdout goes to /dev/null).
Gates: cat_tee_pipe_throughput_budget_exceeded and to_simtel_pipe_throughput_budget_exceeded
(performance_budget.PIPE_LIMIT_GB_PER_S). Peak RSS and CPU of the pipe members are reported.

The pipe files go to a container-local directory (PERF_WORK_DIR, default a new temporary directory),
not to the case directory: on macOS a bind-mounted results directory is a slow file system and
would measure the mount, not the pipe.

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.performance_pipe [--input SHOWERS.iact] [--events N] [--reps N]
"""

from __future__ import annotations

import argparse
import os
import shutil
import statistics
import tempfile
import traceback
from pathlib import Path

from obdeect_validation.adapters import (
    common,
    performance,
    performance_budget,
    performance_chain,
)

REPS = 3


def pipes(native: Path, work: Path) -> dict[str, list[tuple[str, list]]]:
    """The measured pipes: (stage name, argv) in order; each stdout feeds the next stdin."""
    return {
        "cat_tee": [
            ("cat", ["cat", native]),
            ("tee", ["tee", work / "tee-copy.gst"]),
        ],
        "to_simtel": [
            ("cat", ["cat", native]),
            (
                "to_simtel",
                ["obdeect", "stream", "to-simtel", "-o", work / "roundtrip.iact"],
            ),
        ],
    }


def measure_pipe(name, stages, native_bytes, reps, work, metrics, details):
    """Run one pipe reps times; a failure is recorded and the other pipe still runs."""
    records = []
    try:
        for _ in range(reps):
            logs = {stage: work / f"{name}-{stage}.log" for stage, _ in stages}
            records.append(performance_chain.run_pipeline(stages, os.devnull, logs))
    except (RuntimeError, OSError) as error:  # reported, not hidden
        metrics[f"{name}_pipe_completed"] = 0
        details["failures"][name] = {
            "error": str(error)[:400],
            "traceback": traceback.format_exc()[-800:],
        }
        return None
    throughput = [
        performance_budget.throughput_gb_per_s(native_bytes, r["wall_s"])
        for r in records
    ]
    median = {
        key: statistics.median(record[key] for record in records)
        for key in ("wall_s", "cpu_s", "peak_rss_mb")
    }
    metrics.update(
        {
            f"{name}_pipe_completed": 1,
            f"{name}_pipe_reps": len(records),
            f"{name}_pipe_wall_s": median["wall_s"],
            f"{name}_pipe_cpu_s": median["cpu_s"],
            f"{name}_pipe_peak_rss_mb": median["peak_rss_mb"],
            f"{name}_pipe_bytes_in": native_bytes,
            f"{name}_pipe_gb_per_s": statistics.median(throughput),
        }
    )
    details["pipes"][name] = {
        "stages": [stage for stage, _ in stages],
        "records": records,
        "gb_per_s_per_rep": throughput,
    }
    return metrics[f"{name}_pipe_gb_per_s"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=performance.INPUT)
    parser.add_argument(
        "--events", type=int, default=None, help="first N showers (default all)"
    )
    parser.add_argument("--reps", type=int, default=REPS)
    args = parser.parse_args()
    root = Path.cwd().resolve()
    work = Path(
        os.environ.get("PERF_WORK_DIR") or tempfile.mkdtemp(prefix="perf-pipe-")
    )
    work.mkdir(parents=True, exist_ok=True)
    source = performance.select_input(
        args.input.resolve(), args.events, work / "showers-input.iact"
    )
    native = work / "showers.gst"
    common.command(
        ["obdeect", "stream", "from-corsika", source],
        work / "from-corsika.log",
        stdout_file=native,
    )
    native_bytes = native.stat().st_size
    metrics: dict = {}
    details: dict = {
        "host": performance.host_info(),
        "pipes": {},
        "failures": {},
        "settings": {
            "reps": args.reps,
            "input": str(source),
            "events": performance.shower_count(source),
            "work_dir": str(work),
            "native_stream_bytes": native_bytes,
            "limit_gb_per_s": performance_budget.PIPE_LIMIT_GB_PER_S,
            "gigabyte": "decimal, 1e9 bytes",
            "timing": "pipe wall from start to last exit; CPU and RSS per member (wait4, sampled VmHWM)",
        },
    }
    throughput = {}
    for name, stages in pipes(native, work).items():
        throughput[name] = measure_pipe(
            name, stages, native_bytes, args.reps, work, metrics, details
        )
    budgets, reasons = performance_budget.pipe_budgets(throughput)
    metrics.update(budgets)
    details["budget_reasons"] = reasons
    common.write_results(metrics, details, root)
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
