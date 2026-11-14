"""Per-stage throughput of the obdeect chain on the fixed CORSIKA IACT input (report only)."""

import argparse
import json
import os
import platform
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from eventio.iact import IACTFile

from obdeect_validation.adapters import chain_showers

REPS = 3
INPUT = Path("/validation/inputs/gamma-20deg-100-3000GeV.iact")
TRANSMISSION = Path("/reference/cfg/common/atm_trans_2156_1_3_2_0_0_0.1_0.1.dat")
MODELS = Path("/models/simulation-models")
RUSAGE = []


def run_timed(argv, log, cwd=None, stdout_file=None):
    """Run one child and keep its own wall time, CPU time and peak RSS (wait4 rusage)."""
    argv = list(map(str, argv))
    with (
        Path(log).open("wb") as log_stream,
        (
            Path(stdout_file).open("wb") if stdout_file else open(os.devnull, "wb")
        ) as output,
    ):
        start = time.perf_counter()
        process = subprocess.Popen(argv, cwd=cwd, stdout=output, stderr=log_stream)
        _, status, usage = os.wait4(process.pid, 0)
        wall = time.perf_counter() - start
    process.returncode = os.waitstatus_to_exitcode(status)
    rss_divisor = 1024 * 1024 if sys.platform == "darwin" else 1024
    record = {
        "wall_s": wall,
        "cpu_s": usage.ru_utime + usage.ru_stime,
        "peak_rss_mb": usage.ru_maxrss / rss_divisor,
        "returncode": process.returncode,
    }
    RUSAGE.append(record)
    if process.returncode:
        raise RuntimeError(f"{argv[0]} exited {process.returncode}; see {log}")
    return record


def stream_weight(path):
    """Photon count and bunch weight sum of a native photon stream, from obdeect's own summary."""
    text = subprocess.run(
        ["obdeect", "stream", "cat", "--summary", str(path)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    match = re.search(r"total: (\d+) photons, weight sum ([0-9.eE+-]+)", text)
    if not match:
        raise ValueError(f"no total in stream summary of {path}")
    return int(match[1]), float(match[2])


def iact_weight(path):
    """Sum of CORSIKA IACT bunch photon weights (the 'photons' field)."""
    total = 0.0
    with IACTFile(path) as stream:
        for event in stream:
            for photons in event.photon_bunches.values():
                total += float(np.sum(photons["photons"], dtype=np.float64))
    return total


def restrict_table(source, target, limit=1000.0):
    """Keep comment lines and rows whose wavelength [nm] is below limit."""
    kept = []
    for line in source.read_text().splitlines(keepends=True):
        fields = line.split()
        try:
            keep = (
                not fields or line.lstrip().startswith("#") or float(fields[0]) < limit
            )
        except ValueError:
            keep = True
        if keep:
            kept.append(line)
    target.write_text("".join(kept))
    return len(kept)


def cpu_model():
    try:
        text = Path("/proc/cpuinfo").read_text()
    except OSError:
        return platform.processor() or platform.machine()
    for key in ("model name", "Hardware", "cpu model", "CPU part"):
        match = re.search(rf"^{key}\s*:\s*(.+)$", text, re.MULTILINE)
        if match:
            return match[1].strip()
    return platform.machine()


def host_info():
    version = subprocess.run(
        ["obdeect", "--version"], capture_output=True, text=True, check=False
    ).stdout.strip()
    return {
        "nproc": len(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else os.cpu_count(),
        "cpu_model": cpu_model(),
        "machine": platform.machine(),
        "platform": platform.platform(),
        "obdeect_version": version,
    }


def measure(
    name, run_once, photons_in, input_file, output_file, metrics, details, reps
):
    start = len(RUSAGE)
    try:
        for _ in range(reps):
            run_once()
    except (
        Exception
    ) as error:  # a failed stage is recorded; the remaining stages still run
        details["failures"][name] = str(error)
        metrics[f"{name}_completed"] = 0
        return False
    records = RUSAGE[start:]
    median = {
        key: statistics.median(record[key] for record in records)
        for key in ("wall_s", "cpu_s", "peak_rss_mb")
    }
    metrics.update(
        {
            f"{name}_completed": 1,
            f"{name}_reps": len(records),
            f"{name}_wall_s_median": median["wall_s"],
            f"{name}_cpu_s_median": median["cpu_s"],
            f"{name}_peak_rss_mb": median["peak_rss_mb"],
            f"{name}_photons_in": photons_in,
            f"{name}_photons_per_s": photons_in / median["wall_s"],
            f"{name}_bytes_in": input_file.stat().st_size,
            f"{name}_bytes_out": output_file.stat().st_size,
        }
    )
    details["stages"][name] = {
        "input": str(input_file),
        "output": str(output_file),
        "records": records,
    }
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--reps", type=int, default=REPS)
    args = parser.parse_args()
    root = Path.cwd().resolve()
    work = root / "work"
    work.mkdir(exist_ok=True)
    source = args.input.resolve()
    chain_showers.command = (
        run_timed  # sim_telarray reference runs through the same timed runner
    )
    table = work / "transmission-below-1000nm.dat"
    restrict_table(TRANSMISSION, table)
    metrics = {}
    details = {
        "host": host_info(),
        "stages": {},
        "failures": {},
        "skips": {},
        "settings": {
            "reps": args.reps,
            "input": str(source),
            "transmission": str(TRANSMISSION),
            "wavelength_limit_nm": 1000,
            "propagation_seed": 123,
            "trace_seed": 1,
            "pointing": "zenith 20 deg, azimuth 180 deg (telescope_theta/phi of the chain config)",
            "threads": 1,
            "timing": "wait4 rusage of each child; wall time around the child",
            "limits": "conservation: relative delta <= 1e-6 (IACT stores weights as float32; stream sums in float64)",
        },
    }
    iact_in = iact_weight(source)
    done = {}

    done["from_corsika"] = measure(
        "from_corsika",
        lambda: run_timed(
            ["obdeect", "stream", "from-corsika", source],
            work / "from-corsika.log",
            stdout_file=work / "showers.gst",
        ),
        iact_in,
        source,
        work / "showers.gst",
        metrics,
        details,
        args.reps,
    )
    if done["from_corsika"]:
        stream_in = stream_weight(work / "showers.gst")[1]
        metrics["iact_input_vs_stream_weight_relative_delta"] = (
            abs(stream_in - iact_in) / iact_in
        )
    else:
        details["skips"]["to_simtel"] = details["skips"]["propagation"] = (
            "from_corsika failed; no input stream"
        )

    if done["from_corsika"]:
        done["to_simtel"] = measure(
            "to_simtel",
            lambda: run_timed(
                [
                    "obdeect",
                    "stream",
                    "to-simtel",
                    work / "showers.gst",
                    "-o",
                    work / "roundtrip.iact",
                ],
                work / "to-simtel.log",
            ),
            stream_in,
            work / "showers.gst",
            work / "roundtrip.iact",
            metrics,
            details,
            args.reps,
        )
    if done.get("to_simtel"):
        roundtrip_in = iact_weight(work / "roundtrip.iact")
        metrics["from_corsika_to_simtel_weight_relative_delta"] = (
            abs(roundtrip_in - stream_in) / stream_in
        )

    if done["from_corsika"]:
        done["propagation"] = measure(
            "propagation",
            lambda: run_timed(
                [
                    "obdeect",
                    "propagation",
                    "run",
                    "--transmission",
                    table,
                    "--site-altitude",
                    2156,
                    "--pathlen",
                    28,
                    "--seed",
                    123,
                    "--mode",
                    "weighted",
                    "--wavelength-range",
                    250,
                    700,
                    "-i",
                    work / "showers.gst",
                    "-o",
                    work / "propagated.gst",
                ],
                work / "propagation.log",
            ),
            stream_in,
            work / "showers.gst",
            work / "propagated.gst",
            metrics,
            details,
            args.reps,
        )

    done["trace_stream"] = False
    if done.get("propagation"):
        trace_in = stream_weight(work / "propagated.gst")[1]
        model = work / "lst.model"
        try:
            subprocess.run(
                [
                    "obdeect",
                    "ray-tracing",
                    "compile",
                    "--simulation-models",
                    str(MODELS),
                    "--model",
                    "LSTN-design",
                    "--version",
                    "7.0.0",
                    "--seed",
                    "1",
                    "--output",
                    str(model),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        except (subprocess.CalledProcessError, OSError) as error:
            details["skips"]["trace_stream"] = f"optical model compile failed: {error}"
        else:
            done["trace_stream"] = measure(
                "trace_stream",
                lambda: run_timed(
                    [
                        "obdeect",
                        "ray-tracing",
                        "trace-stream",
                        "--model",
                        model,
                        "--seed",
                        1,
                        "--threads",
                        1,
                        "--pointing-zenith-deg",
                        20,
                        "--pointing-azimuth-deg",
                        180,
                        "-i",
                        work / "propagated.gst",
                        "-o",
                        work / "traced.gst",
                    ],
                    work / "trace-stream.log",
                ),
                trace_in,
                work / "propagated.gst",
                work / "traced.gst",
                metrics,
                details,
                args.reps,
            )
    else:
        details["skips"].setdefault(
            "trace_stream", "propagation failed or skipped; no input stream"
        )

    done["simtel_direct"] = measure(
        "simtel_direct",
        lambda: chain_showers.simulate(work, source, "direct"),
        iact_in,
        source,
        work / "direct.simtel",
        metrics,
        details,
        args.reps,
    )

    completed = [name for name in done if done[name]]
    metrics["stages_failed"] = len(details["failures"])
    metrics["stages_skipped"] = len(details["skips"])
    metrics["stages_completed"] = len(completed)
    details["stage_order"] = list(done)

    if completed:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        values = [metrics[f"{name}_photons_per_s"] for name in completed]
        fig, axis = plt.subplots(figsize=(7, 3.6))
        bars = axis.barh(completed, values, color="#2f6f9f", height=0.6)
        axis.set_xscale("log")
        axis.set_xlabel("photons per second (log scale, median of repetitions)")
        axis.invert_yaxis()
        for bar, value in zip(bars, values):
            axis.text(
                value * 1.05,
                bar.get_y() + bar.get_height() / 2,
                f"{value:.3g}",
                va="center",
                fontsize=9,
            )
        axis.set_xlim(right=max(values) * 10)
        axis.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(root / "photons-per-second.png", dpi=150)
        plt.close(fig)

    (root / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    (root / "details.json").write_text(json.dumps(details, indent=2) + "\n")

    for path in (
        work / "showers.gst",
        work / "roundtrip.iact",
        work / "propagated.gst",
        work / "traced.gst",
        work / "direct.simtel",
        work / "direct.hdata",
        work / "roundtrip.hdata",
    ):
        path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
