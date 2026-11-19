"""Per-stage throughput of the obdeect chain on the fixed CORSIKA IACT input (report only).

Same setup as the air-shower case (performance_chain): LSTN-01 7.0.0 production atmosphere cut
below 1000 nm, QE envelope, one telescope, pointing of the chain config. Stages, in order:
from_corsika, to_simtel (round trip, conservation check), propagation (weighted), propagation_envelope
(thinned by the QE bound), compile_model (reported only), trace_stream (weighted input), simtel_direct
(sim_telarray Bypass_Optics=0 reference). The throughput stages are the timed ones in stages_completed.

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.performance --config cfg/model/7.0.0/CTAO.cfg [--events N] [--reps N]

Timing: each child is timed with performance_chain.run_timed (obdeect subcommands as native programs,
as in the air-shower case). Each stage records the load average before it starts; the host load
before and after the run and the spread of each stage's wall time (timing_unreliable) tell whether
the numbers are usable.
"""

import argparse
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import numpy as np
from eventio import EventIOFile
from eventio.iact import IACTFile

from obdeect_validation.adapters import common

REPS = 3
INPUT = common.VALIDATION / "inputs/gamma-20deg-100-3000GeV-100showers.iact"
RUSAGE = []
SPREAD_LIMIT = (
    0.05  # (max - min) / median of a stage's wall time above this: timing unreliable
)
LOAD_LIMIT = (
    1.0  # host load average (1 min) before the run above this: timing unreliable
)
CPU_SHARE_STAGES = ("from_corsika", "propagation_envelope", "trace_stream", "to_simtel")


def loadavg():
    """1-minute load average of the host, or None where it is not available."""
    try:
        return os.getloadavg()[0]
    except OSError, AttributeError:
        return None


PROC_VMHWM = re.compile(r"^VmHWM:\s+(\d+)\s+kB", re.MULTILINE)


def parse_vmhwm_kb(text):
    """Peak resident set [kB] from the text of /proc/<pid>/status, or None when absent."""
    match = PROC_VMHWM.search(text)
    return int(match[1]) if match else None


def read_vmhwm_kb(pid):
    try:
        return parse_vmhwm_kb(Path(f"/proc/{pid}/status").read_text())
    except OSError:  # the child has exited, or there is no /proc
        return None


class RssSampler:
    """Peak RSS of one running child, sampled from /proc until stop() (Linux only).

    wait4 ru_maxrss is not used on Linux: exec keeps the high-water mark of the forked parent,
    so every child of a large Python process reports the parent's peak (identical values for all
    stages). The sampled VmHWM is the child's own peak, a lower bound at the sampling interval.
    """

    def __init__(self, pid, read=None, interval=0.002):
        self.read = read or (lambda: read_vmhwm_kb(pid))
        self.interval = interval
        self.peak_kb = None
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        while True:
            kb = self.read()
            if kb is not None and (self.peak_kb is None or kb > self.peak_kb):
                self.peak_kb = kb
            if self._stop.wait(self.interval):
                return

    def start(self):
        self._thread.start()

    def stop(self):
        """Stop sampling and return the peak [kB], or None when no sample was taken."""
        self._stop.set()
        self._thread.join()
        return self.peak_kb


def run_timed(argv, log, cwd=None, stdout_file=None):
    """Run one child and keep its own wall time, CPU time and peak RSS.

    Peak RSS: sampled VmHWM on Linux, wait4 ru_maxrss elsewhere (see RssSampler). The record
    keeps the raw ru_maxrss as rusage_maxrss_mb for comparison.
    """
    argv = list(map(str, argv))
    with (
        Path(log).open("wb") as log_stream,
        Path(stdout_file).open("wb") if stdout_file else log_stream as output,
    ):
        load_before = loadavg()
        start = time.perf_counter()
        process = subprocess.Popen(argv, cwd=cwd, stdout=output, stderr=log_stream)
        sampler = RssSampler(process.pid) if sys.platform.startswith("linux") else None
        if sampler:
            sampler.start()
        _, status, usage = os.wait4(process.pid, 0)
        wall = time.perf_counter() - start
        sampled_kb = sampler.stop() if sampler else None
    process.returncode = os.waitstatus_to_exitcode(status)
    rss_divisor = 1024 * 1024 if sys.platform == "darwin" else 1024
    rusage_mb = usage.ru_maxrss / rss_divisor
    if sampled_kb is not None:
        peak_mb, source = sampled_kb / 1024, "proc-VmHWM-sampled"
    else:
        peak_mb, source = rusage_mb, "wait4-ru_maxrss"
    record = {
        "argv": argv,
        "wall_s": wall,
        "cpu_s": usage.ru_utime + usage.ru_stime,
        "peak_rss_mb": peak_mb,
        "peak_rss_source": source,
        "rusage_maxrss_mb": rusage_mb,
        "returncode": process.returncode,
        "loadavg_before": load_before,
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


EVENT_HEADER = (
    1202  # eventio type of a CORSIKA IACT event header; the event's blocks follow it
)
EVENT_TYPES = range(
    1202, 1210
)  # EventHeader .. EventEnd (1203-1209: the event's telescope blocks)
BLOCK_HEADER_BYTES = 16  # eventio block header in front of each top-level block payload


def _top_level_blocks(path):
    """(type, start, end) byte range of each top-level eventio block, header included."""
    with EventIOFile(path) as stream:
        for block in stream:
            start = block.address - BLOCK_HEADER_BYTES
            yield block.header.type, start, block.address + block.size


def shower_count(path):
    """Number of showers (event headers) of a CORSIKA IACT file."""
    return sum(
        1 for block_type, _, _ in _top_level_blocks(path) if block_type == EVENT_HEADER
    )


def subset_showers(source, target, events):
    """Copy the first `events` showers of an IACT file as raw eventio blocks; the run header and
    the blocks after the events (run end) are kept. Returns the number of showers written."""
    if events < 0:
        raise ValueError("number of showers must not be negative")
    kept = seen = 0
    with Path(source).open("rb") as raw, Path(target).open("wb") as out:
        for block_type, start, end in _top_level_blocks(source):
            if block_type == EVENT_HEADER:
                seen += 1
            if block_type not in EVENT_TYPES or seen <= events:
                kept += block_type == EVENT_HEADER
                raw.seek(start)
                out.write(raw.read(end - start))
    return kept


def select_input(source, events, target):
    """The input to measure: source itself (events None), or its first `events` showers in target."""
    if events is None:
        return source
    subset_showers(source, target, events)
    return target


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
        "host": os.environ.get("OBDEECT_VALIDATION_HOST"),
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
    wall = [record["wall_s"] for record in records]
    metrics.update(
        {
            f"{name}_completed": 1,
            f"{name}_reps": len(records),
            f"{name}_wall_s_median": median["wall_s"],
            f"{name}_wall_s_min": min(wall),
            f"{name}_wall_s_max": max(wall),
            f"{name}_wall_s_spread": (max(wall) - min(wall)) / median["wall_s"],
            f"{name}_cpu_s_median": median["cpu_s"],
            f"{name}_peak_rss_mb": median["peak_rss_mb"],
            f"{name}_bytes_in": input_file.stat().st_size,
            f"{name}_bytes_out": output_file.stat().st_size,
        }
    )
    if (
        photons_in is not None
    ):  # a stage without a photon count (compile_model) has no rate
        metrics.update(
            {
                f"{name}_photons_in": photons_in,
                f"{name}_photons_per_s": photons_in / median["wall_s"],
            }
        )
    details["stages"][name] = {
        "input": str(input_file),
        "output": str(output_file),
        "records": records,
    }
    return True


def cpu_metrics(metrics, stages):
    """Per stage: CPU share of the envelope chain (CPU_SHARE_STAGES) and CPU ns per photon.

    The share divides by the CPU medians of the CPU_SHARE_STAGES present in metrics.
    ns per photon uses the stage's photons_in (the rate convention of photons_per_s).
    """
    chain = [
        metrics[f"{name}_cpu_s_median"]
        for name in CPU_SHARE_STAGES
        if f"{name}_cpu_s_median" in metrics
    ]
    chain_cpu = sum(chain)
    out = {}
    for name in stages:
        cpu = metrics[f"{name}_cpu_s_median"]
        if chain_cpu > 0:
            out[f"{name}_cpu_share"] = cpu / chain_cpu
        photons = metrics.get(f"{name}_photons_in")
        if photons:
            out[f"{name}_ns_per_photon_cpu"] = cpu * 1e9 / photons
    return out


def timing_unreliable(metrics, stages):
    """1 when a stage's wall-time spread is above SPREAD_LIMIT or the host load before the run is
    above LOAD_LIMIT; a missing load average does not count."""
    spreads = [metrics.get(f"{name}_wall_s_spread") or 0.0 for name in stages]
    load = metrics.get("host_loadavg_before")
    busy = load is not None and load > LOAD_LIMIT
    return int(busy or any(spread > SPREAD_LIMIT for spread in spreads))


def _plot(path, stages, metrics):
    """Throughput per stage (one colour per stage), with the per-stage CPU and wall-time medians."""
    from obdeect_validation.adapters import plotting

    fig, axes = plotting.figure(1, 2, (11, 3.8))
    colors = [plotting.CYCLE[i % len(plotting.CYCLE)] for i in range(len(stages))]
    values = [metrics[f"{name}_photons_per_s"] for name in stages]
    bars = axes[0, 0].barh(stages, values, color=colors, height=0.6)
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_xlabel("photons per second (log scale, median of repetitions)")
    axes[0, 0].invert_yaxis()
    axes[0, 0].grid(axis="y", visible=False)
    for bar, value in zip(bars, values, strict=True):
        axes[0, 0].text(
            value * 1.05,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.3g}",
            va="center",
            fontsize=9,
        )
    axes[0, 0].set_xlim(right=max(values) * 10)
    wall = [metrics[f"{name}_wall_s_median"] for name in stages]
    cpu = [metrics[f"{name}_cpu_s_median"] for name in stages]
    rows = range(len(stages))
    axes[0, 1].barh(
        [r - 0.2 for r in rows], wall, height=0.4, color=colors, label="wall"
    )
    axes[0, 1].barh(
        [r + 0.2 for r in rows], cpu, height=0.4, color=colors, alpha=0.45, label="CPU"
    )
    axes[0, 1].set_yticks(list(rows), stages)
    axes[0, 1].invert_yaxis()
    axes[0, 1].set_xlabel("seconds (median)")
    axes[0, 1].grid(axis="y", visible=False)
    axes[0, 1].legend()
    plotting.save(fig, path)


def main():
    # Imported here: performance_chain imports this module at module level.
    from obdeect_validation.adapters import chain_optics as co
    from obdeect_validation.adapters import performance_chain as chain

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="sim_telarray CTAO.cfg of the LSTN-01 7.0.0 array (simtools generate-array-config)",
    )
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument(
        "--events",
        type=int,
        default=None,
        help="measure the first N showers (default all)",
    )
    parser.add_argument("--reps", type=int, default=REPS)
    args = parser.parse_args()
    root = Path.cwd().resolve()
    # Stream files go to a container-local directory (PERF_WORK_DIR), not the bind-mounted results.
    work_dir = os.environ.get("PERF_WORK_DIR")
    work = Path(work_dir or tempfile.mkdtemp(prefix="perf-stages-"))
    work.mkdir(parents=True, exist_ok=True)
    metrics = {}
    load_before = loadavg()
    if load_before is not None:
        metrics["host_loadavg_before"] = load_before
    source = select_input(
        args.input.resolve(), args.events, work / "showers-input.iact"
    )
    inputs, cfg_text = chain.prepare_inputs(work, args.config)
    headroom = chain.headroom_for(inputs["qe_variation"])
    (work / "envelope").mkdir(parents=True, exist_ok=True)
    envelope = chain.envelope_files(work / "envelope", inputs["qe"], headroom)
    details = {
        "host": host_info(),
        "stages": {},
        "failures": {},
        "skips": {},
        "settings": {
            "reps": args.reps,
            "input": str(source),
            "showers": shower_count(source),
            "work_dir": str(work),
            "config": str(inputs["config"]),
            "model": f"{co.MODEL} {co.VERSION}",
            "reference_atmosphere": co.cfg_value(cfg_text, "atmospheric_transmission"),
            "propagation_atmosphere": str(inputs["atmosphere"]),
            "wavelength_limit_nm": 1000,
            "propagation": "weighted 250-700 nm (as shower_candidate); envelope: thinned by the QE bound (as envelope_pipeline)",
            "propagation_seed": 123,
            "envelope_headroom": headroom,
            "envelope_bound_max": envelope["bound_max"],
            "telescopes": list(co.TELESCOPES),
            "pointing": f"zenith {co.ZENITH_DEG} deg, azimuth {co.AZIMUTH_DEG} deg",
            "trace_seed": co.SEED,
            "threads": 1,
            "simtel_direct": "performance_chain.shower_reference: sim_telarray Bypass_Optics=0, one telescope, reference atmosphere",
            "timing": "performance_chain.run_timed: wait4 rusage of each child (obdeect subcommands as native programs); wall time around the child; loadavg_before per record",
            "limits": "conservation: relative delta <= 1e-6 (IACT stores weights as float32; stream sums in float64)",
        },
    }
    iact_in = iact_weight(source)
    stream = work / "showers.gst"
    propagated = work / "propagated.gst"
    enveloped = work / "propagated-envelope.gst"
    model = work / "lst.model"
    traced = work / "traced.gst"
    simtel_dir = work / "simtel-direct"
    simtel_dir.mkdir(parents=True, exist_ok=True)
    run = chain.run_timed
    done = {}

    done["from_corsika"] = measure(
        "from_corsika",
        lambda: run(
            ["obdeect", "stream", "from-corsika", source],
            work / "from-corsika.log",
            stdout_file=stream,
        ),
        iact_in,
        source,
        stream,
        metrics,
        details,
        args.reps,
    )
    if done["from_corsika"]:
        stream_in = stream_weight(stream)[1]
        metrics["iact_input_vs_stream_weight_relative_delta"] = (
            abs(stream_in - iact_in) / iact_in
        )
        done["to_simtel"] = measure(
            "to_simtel",
            lambda: run(
                [
                    "obdeect",
                    "stream",
                    "to-simtel",
                    stream,
                    "-o",
                    work / "roundtrip.iact",
                ],
                work / "to-simtel.log",
            ),
            stream_in,
            stream,
            work / "roundtrip.iact",
            metrics,
            details,
            args.reps,
        )
    else:
        stream_in = None
        skip = "from_corsika failed; no input stream"
        details["skips"].update(
            to_simtel=skip, propagation=skip, propagation_envelope=skip
        )
    if done.get("to_simtel"):
        roundtrip_in = iact_weight(work / "roundtrip.iact")
        metrics["from_corsika_to_simtel_weight_relative_delta"] = (
            abs(roundtrip_in - stream_in) / stream_in
        )

    if done["from_corsika"]:
        done["propagation"] = measure(
            "propagation",
            lambda: run(
                [
                    "obdeect",
                    "propagation",
                    "run",
                    "--transmission",
                    inputs["atmosphere"],
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
                    "--telescopes",
                    ",".join(map(str, co.TELESCOPES)),
                    "-i",
                    stream,
                    "-o",
                    propagated,
                ],
                work / "propagation.log",
            ),
            stream_in,
            stream,
            propagated,
            metrics,
            details,
            args.reps,
        )
        done["propagation_envelope"] = measure(
            "propagation_envelope",
            lambda: run(
                [
                    "obdeect",
                    "propagation",
                    "run",
                    "--transmission",
                    inputs["atmosphere"],
                    "--site-altitude",
                    2156,
                    "--pathlen",
                    28,
                    "--seed",
                    123,
                    "--envelope",
                    f"{chain.ENVELOPE_FACTOR}={envelope['envelope']}",
                    "--telescopes",
                    ",".join(map(str, co.TELESCOPES)),
                    "-i",
                    stream,
                    "-o",
                    enveloped,
                ],
                work / "propagation-envelope.log",
            ),
            stream_in,
            stream,
            enveloped,
            metrics,
            details,
            args.reps,
        )
        if done["propagation_envelope"]:
            metrics["propagation_envelope_photons_out"] = stream_weight(enveloped)[0]

    done["compile_model"] = measure(
        "compile_model",
        lambda: run(chain.compile_argv(inputs, model), work / "compile.log"),
        None,
        inputs["config"],
        model,
        metrics,
        details,
        args.reps,
    )

    done["trace_stream"] = False
    if done.get("propagation") and done["compile_model"]:
        trace_in = stream_weight(propagated)[1]
        done["trace_stream"] = measure(
            "trace_stream",
            lambda: run(
                [
                    "obdeect",
                    "ray-tracing",
                    "trace-stream",
                    "--model",
                    model,
                    "--seed",
                    co.SEED,
                    "--threads",
                    1,
                    "--telescopes",
                    ",".join(map(str, co.TELESCOPES)),
                    "--pointing-zenith-deg",
                    co.ZENITH_DEG,
                    "--pointing-azimuth-deg",
                    co.AZIMUTH_DEG,
                    "--telescope-transmission",
                    "full",
                    "-i",
                    propagated,
                    "-o",
                    traced,
                ],
                work / "trace-stream.log",
            ),
            trace_in,
            propagated,
            traced,
            metrics,
            details,
            args.reps,
        )
    else:
        details["skips"].setdefault(
            "trace_stream",
            "propagation or model compile failed or skipped; no input stream",
        )

    _name, reference_run = chain.shower_reference(
        simtel_dir,
        source,
        inputs["config"],
        inputs["config_dir"],
        inputs["reference_atmosphere"],
    )[0]
    done["simtel_direct"] = measure(
        "simtel_direct",
        reference_run,
        iact_in,
        source,
        simtel_dir / "reference.simtel",
        metrics,
        details,
        args.reps,
    )

    # compile_model is reported, not a throughput stage (no photon rate, not in stages_completed).
    completed = [name for name in done if done[name] and name != "compile_model"]
    metrics["stages_failed"] = len(details["failures"])
    metrics["stages_skipped"] = len(details["skips"])
    metrics["stages_completed"] = len(completed)
    metrics.update(cpu_metrics(metrics, completed))
    load_after = loadavg()
    if load_after is not None:
        metrics["host_loadavg_after"] = load_after
    metrics["timing_unreliable"] = timing_unreliable(metrics, list(details["stages"]))
    details["stage_order"] = [name for name in done if name != "compile_model"]

    if completed:
        _plot(root / "photons-per-second.png", completed, metrics)

    common.write_results(metrics, details, root)

    if work_dir is None:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    # Run the canonical module: performance_chain imports it, and both must share one RUSAGE list.
    from obdeect_validation.adapters import performance

    performance.main()
