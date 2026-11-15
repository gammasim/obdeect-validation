"""Matched end-to-end timing: obdeect pipelines against their reference programs (report only).

Each case times a reference pipeline and an obdeect pipeline on the same physical input and
photon counts. The reference is CORSIKA IACT -> sim_telarray with full optics (air shower) or
LightEmission -> sim_telarray (sources). The candidate is the obdeect chain ending in the
sim_telarray camera stage (Bypass_Optics=1 for the air shower, as in cases/chain/optics-bypass1.json).
Every child process is timed with performance.run_timed. A pipeline total is the sum of its stage
wall and CPU times and the maximum of its stage peak RSS (stages run one after another).

The air-shower case also times an optimised candidate (candidate_envelope): from-corsika,
propagation (QE bound, thinned), trace-stream and to-simtel run as one pipe; sim_telarray then
reads a QE / bound table, so the product is the model QE. A physics check compares its camera
output with sim_telarray on the same input (unit atmosphere on both sides, report only).

Only ray-tracing trace-stream is multithreaded. The air-shower candidate therefore runs with
threads=1 and with the obdeect default (0 = every core); the source pipelines are single-threaded
and run once. Ratios (candidate / reference) are informational: no performance limit is defined.

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.performance_chain air-shower --config cfg/model/7.0.0/CTAO.cfg
    python -m obdeect_validation.adapters.performance_chain flasher-dish
    python -m obdeect_validation.adapters.performance_chain point-source
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import statistics
import subprocess
import sys
import threading
import time
import traceback
from pathlib import Path

import numpy as np
from eventio.iact import IACTFile

from obdeect_validation.adapters import chain_optics as co
from obdeect_validation.adapters import chain_showers as cs
from obdeect_validation.adapters import common, performance, source_camera
from obdeect_validation.adapters.performance import iact_weight, run_timed

REPS = 3
# Production atmosphere of the LSTN-01 model (CTAO.cfg, set only under #if TELESCOPE == 0). The
# reference and every candidate pass it explicitly and the case checks the sim_telarray log for it.
ATMOSPHERE_FILE = "atmospheric_transmission-2.0.0-OBS-North.dat"
ENVELOPE_FACTOR = "qe"
# An efficiency bound must be positive: bins with QE 0 keep this bound (their table value is 0).
ENVELOPE_FLOOR = 1e-6
# Table maximum of the sim_telarray QE / bound table. sim_telarray's first-cut probability is
# p2 * max_qe_rel (pixel qe_variation, 0.063 in LSTN-01), so a table at 1 saturates for upward
# pixel fluctuations and loses photons (clustered bias -6.6 % at headroom 1, 0 % at 0.8).
ENVELOPE_HEADROOM = 0.8
# Second node of each nm bin: the linear bound is constant on [n, n + 1) nm.
ENVELOPE_NODE_WIDTH = 0.999999
# CORSIKA Cherenkov band (EVTH words 96-97) on sim_telarray's integer nm bins.
WAVELENGTH_GRID = np.arange(100, 2001)
TRACE_SUMMARY = re.compile(r"trace-stream: \d+ events, (\d+) photons in, (\d+) out")
ALL_PHOTONS = re.compile(r"All photons:\s+(\d+)")
# sim_telarray for the candidate camera stage (Bypass_Optics=3 build); None: the production binary.
CANDIDATE_SIMTEL = os.environ.get("SIMTEL_CANDIDATE_BINARY") or None
THREAD_VARIANTS = (
    ("threads1", 1),
    ("threadsdefault", 0),
)  # 0 = obdeect default, every core
PHOTON_KEYS = ("wall_s", "cpu_s", "peak_rss_mb")
SKY_PROPAGATION = [
    "--no-extinction",
    "--site-altitude",
    "2150",
    "--seed",
    "1",
    "--wavelength-range",
    "300",
    "800",
]  # as in source_camera stage 2

SOURCES = {
    "flasher-dish": {
        "events": 2,
        "photons": 20_000_000,
        "candidate": lambda n, p: [
            "obdeect",
            "sources",
            "flasher",
            "--type",
            "dish",
            "--events",
            n,
            "--photons",
            p,
        ],
        "reference": lambda n, p: ["ff-1m", "--events", n, "--photons", p],
        "bypass": 1,  # as in source_camera stage 2
    },
    "point-source": {
        "events": 2,
        "photons": 10_000_000,
        "candidate": lambda n, p: [
            "obdeect",
            "sources",
            "point-source",
            "-N",
            n,
            "-n",
            p,
            "-z",
            "100000",
            "-a",
            "isotropic",
            "--seed",
            "17",
        ],
        "reference": lambda n, p: [
            "xyzls",
            "-N",
            n,
            "-n",
            p,
            "-z",
            "100000",
            "-a",
            "isotropic",
        ],
        "bypass": 0,  # photons reach the telescope sphere: full optics, as in source_camera
    },
}


def ratio(value, reference):
    """value / reference, or None when the reference is not positive."""
    return value / reference if reference and reference > 0 else None


def telescope_arguments(config_dir, atmosphere=None):
    """sim_telarray arguments of the matched one-telescope model (as in chain_optics)."""
    extra = [
        f"-I{config_dir}",
        "-DTELESCOPE=1",
        "-C",
        f"maximum_telescopes={len(co.TELESCOPES)}",
        "-C",
        "trigger_telescopes=1",  # one telescope: the default array trigger needs two
    ]
    if atmosphere:
        extra += ["-C", f"atmospheric_transmission={atmosphere}"]
    return extra


def summarise(children):
    """One stage (one or more timed children) as a single record."""
    if not children:
        raise RuntimeError("stage ran no timed child process")
    record = {
        "wall_s": sum(c["wall_s"] for c in children),
        "cpu_s": sum(c["cpu_s"] for c in children),
        "peak_rss_mb": max(c["peak_rss_mb"] for c in children),
        "argv": [c["argv"] for c in children],
    }
    members = [member for c in children for member in c.get("children", [])]
    if members:  # a pipeline record: its processes, each with its own wall, CPU and RSS
        record["children"] = members
    return record


def total(stages):
    """Pipeline total of one repetition: sum of wall and CPU, maximum of peak RSS."""
    return {
        "wall_s": sum(s["wall_s"] for s in stages),
        "cpu_s": sum(s["cpu_s"] for s in stages),
        "peak_rss_mb": max(s["peak_rss_mb"] for s in stages),
    }


def measure(steps, reps):
    """Run (name, step) pairs in order, reps times. Returns per-stage records and pipeline totals."""
    per_stage = {name: [] for name, _ in steps}
    totals = []
    for _ in range(reps):
        rep = []
        for name, step in steps:
            start = len(performance.RUSAGE)
            step()
            summary = summarise(performance.RUSAGE[start:])
            per_stage[name].append(summary)
            rep.append(summary)
        totals.append(total(rep))
    return per_stage, totals


def _median(records, key):
    return statistics.median(record[key] for record in records)


def _aggregate(prefix, records):
    return {
        f"{prefix}_wall_s": _median(records, "wall_s"),
        f"{prefix}_cpu_s": _median(records, "cpu_s"),
        f"{prefix}_peak_rss_mb": _median(records, "peak_rss_mb"),
    }


def pipeline_metrics(prefix, per_stage, totals):
    """Medians over repetitions: pipeline totals and each stage ('<prefix>_<stage>_wall_s', ...)."""
    metrics = {f"{prefix}_completed": 1, f"{prefix}_reps": len(totals)}
    metrics.update(_aggregate(prefix, totals))
    for name, records in per_stage.items():
        metrics.update(_aggregate(f"{prefix}_{name}", records))
    return metrics


def ratio_metrics(label, candidate, reference):
    """Candidate / reference for wall, CPU and peak RSS of the pipeline totals."""
    out = {}
    for key, name in (
        ("wall_s", "wall"),
        ("cpu_s", "cpu"),
        ("peak_rss_mb", "peak_rss"),
    ):
        value = ratio(candidate[key], reference[key])
        if value is not None:
            out[f"ratio_{name}_{label}"] = value
    return out


def _reference_sim(work, source, config, config_dir, atmosphere, name="reference"):
    """sim_telarray Bypass_Optics=0 on the CORSIKA input; its atmosphere is checked from the log."""
    output = cs.simulate(
        work,
        source,
        name,
        bypass=0,
        config=str(config),
        extra=telescope_arguments(config_dir, atmosphere),
        runner=run_timed,
    )
    common.check_atmosphere_log(work / f"{name}.log", ATMOSPHERE_FILE)
    return output


def shower_reference(work, source, config, config_dir, atmosphere=None):
    """CORSIKA IACT -> sim_telarray with full optics and the production atmosphere (one telescope)."""
    return [
        (
            "sim_telarray",
            lambda: _reference_sim(work, source, config, config_dir, atmosphere),
        )
    ]


def shower_candidate(work, source, inputs, threads):
    """CORSIKA IACT -> obdeect stream, atmosphere, optics, camera stream -> sim_telarray (Bypass_Optics=1)."""
    native, propagated = work / "showers.gst", work / "propagated.gst"
    model, traced = work / "lst.model", work / "traced.gst"
    camera_stream = work / "candidate.iact"
    camera_map = co.camera_map(work, inputs["rotation"])
    return [
        (
            "from_corsika",
            lambda: run_timed(
                ["obdeect", "stream", "from-corsika", source],
                work / "from-corsika.log",
                stdout_file=native,
            ),
        ),
        (
            "propagation",
            lambda: run_timed(
                [
                    "obdeect",
                    "propagation",
                    "run",
                    "--transmission",
                    inputs["table"],
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
                    native,
                    "-o",
                    propagated,
                ],
                work / "propagation.log",
            ),
        ),
        (
            "compile_model",
            lambda: run_timed(
                [
                    "obdeect",
                    "ray-tracing",
                    "compile",
                    "--simulation-models",
                    inputs["models"],
                    "--model",
                    co.MODEL,
                    "--version",
                    co.VERSION,
                    "--seed",
                    co.SEED,
                    "--output",
                    model,
                ],
                work / "compile.log",
            ),
        ),
        (
            "trace_stream",
            lambda: run_timed(
                [
                    "obdeect",
                    "ray-tracing",
                    "trace-stream",
                    "--model",
                    model,
                    "--seed",
                    co.SEED,
                    "--threads",
                    threads,
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
                work / "trace.log",
            ),
        ),
        (
            "to_simtel",
            lambda: run_timed(
                [
                    "obdeect",
                    "stream",
                    "to-simtel",
                    "--camera-map",
                    camera_map,
                    "-o",
                    camera_stream,
                    traced,
                ],
                work / "to-simtel.log",
            ),
        ),
        (
            "sim_telarray",
            lambda: cs.simulate(
                work,
                camera_stream,
                "candidate",
                bypass=1,
                config=str(inputs["config"]),
                extra=telescope_arguments(inputs["config_dir"], inputs["unit"]),
                runner=run_timed,
            ),
        ),
    ]


def source_candidate(work, spec):
    """obdeect sources -> propagation -> stream to-simtel -> sim_telarray (camera stage)."""
    native, propagated = work / "native.gst", work / "propagated.gst"
    candidate = work / "candidate.iact"
    argv = spec["candidate"](spec["events"], spec["photons"])
    return [
        ("sources", lambda: run_timed([*argv, "-o", native], work / "sources.log")),
        (
            "propagation",
            lambda: run_timed(
                [
                    "obdeect",
                    "propagation",
                    "run",
                    *SKY_PROPAGATION,
                    "-i",
                    native,
                    "-o",
                    propagated,
                ],
                work / "propagation.log",
            ),
        ),
        (
            "to_simtel",
            lambda: run_timed(
                ["obdeect", "stream", "to-simtel", "-o", candidate, propagated],
                work / "to-simtel.log",
            ),
        ),
        (
            "sim_telarray",
            lambda: source_camera._simtel(
                work, candidate, "candidate", spec["bypass"], runner=run_timed
            ),
        ),
    ]


def source_reference(work, spec):
    """LightEmission -> sim_telarray with the same photon counts and camera options."""
    reference = work / "reference.dat.gz"
    argv = source_camera._reference_argv(
        spec["reference"](spec["events"], spec["photons"]), work
    )
    return [
        (
            "light_emission",
            lambda: run_timed([*argv, "-o", reference], work / "reference.log"),
        ),
        (
            "sim_telarray",
            lambda: source_camera._simtel(
                work, reference, "reference", spec["bypass"], runner=run_timed
            ),
        ),
    ]


def _run(label, steps, reps, metrics, details):
    """Measure one pipeline; a failure is recorded and the other pipelines still run."""
    try:
        per_stage, totals = measure(steps, reps)
    except (RuntimeError, OSError, ValueError) as error:  # reported, not hidden
        metrics[f"{label}_completed"] = 0
        details["failures"][label] = {
            "error": str(error)[:400],
            "traceback": traceback.format_exc()[-800:],
        }
        return False
    metrics.update(pipeline_metrics(label, per_stage, totals))
    details["stages"][label] = per_stage
    details["pipelines"][label] = totals
    return True


def _ratios(metrics, reference_ok, candidates):
    """Candidate / reference ratios of the pipeline totals (only when both completed)."""
    if not reference_ok:
        return
    ref = {key: metrics[f"reference_{key}"] for key in PHOTON_KEYS}
    for label in candidates:
        if metrics.get(f"{label}_completed") != 1:
            continue
        cand = {key: metrics[f"{label}_{key}"] for key in PHOTON_KEYS}
        metrics.update(ratio_metrics(label.removeprefix("candidate_"), cand, ref))


def sim_telarray_qe(path, grid=WAVELENGTH_GRID):
    """QE as sim_telarray reads it: first two columns, linear in wavelength, at integer nm."""
    table = np.loadtxt(path, comments="#", usecols=(0, 1), ndmin=2)
    return np.interp(grid, table[:, 0], table[:, 1])


def bin_envelope(qe, headroom=ENVELOPE_HEADROOM):
    """Efficiency bound per nm bin (QE / headroom, at least ENVELOPE_FLOOR) and the table QE / bound.

    The bound is at least the QE in every bin, and the table is at most headroom (0.8), so
    bound times table is the QE again.
    """
    qe = np.asarray(qe, dtype=float)
    bound = np.where(qe > 0, qe / headroom, ENVELOPE_FLOOR)
    return bound, qe / bound


def write_envelope(path, grid, bound):
    """Two-column bound with nodes at n and n + 0.999999 nm, so it is constant on each nm bin.

    sim_telarray bins a photon by its wavelength truncated to integer nm. Constant bins make the
    thinning probability of a photon equal to the QE that sim_telarray applies to its bin.
    """
    lines = []
    for n, value in zip(grid, bound):
        lines.append(f"{n:.6f} {value:.9g}")
        lines.append(f"{n + ENVELOPE_NODE_WIDTH:.6f} {value:.9g}")
    Path(path).write_text("\n".join(lines) + "\n")


def write_scaled(path, grid, scaled):
    """sim_telarray quantum-efficiency table QE / bound at integer nm."""
    Path(path).write_text(
        "".join(f"{n} {value:.9g}\n" for n, value in zip(grid, scaled))
    )


def headroom_for(qe_variation, sigmas=6):
    """Table maximum that keeps the sim_telarray first-cut probability below 1.

    The first cut is p2 * max_qe_rel with max_qe_rel <= 1 + sigmas * qe_variation (per-pixel QE
    spread). 1 / (1 + sigmas * qe_variation) keeps p2 * max_qe_rel <= 1: 0.7257 for LSTN-01 (0.063).
    """
    return 1.0 / (1.0 + sigmas * float(qe_variation))


def envelope_files(work, qe_path, headroom=ENVELOPE_HEADROOM):
    """Efficiency bound for propagation, and the QE / bound table in a directory that holds
    a file of the config's name, so that sim_telarray reads it instead of the model's QE table."""
    bound, scaled = bin_envelope(sim_telarray_qe(qe_path), headroom=headroom)
    table_dir = work / "qe-table"
    table_dir.mkdir(parents=True, exist_ok=True)
    envelope = work / "qe-envelope.dat"
    write_envelope(envelope, WAVELENGTH_GRID, bound)
    write_scaled(table_dir / qe_path.name, WAVELENGTH_GRID, scaled)
    return {
        "envelope": envelope,
        "table_dir": table_dir,
        "bound_max": float(bound.max()),
        "headroom": headroom,
    }


def compile_argv(inputs, model):
    """obdeect ray-tracing compile of the LSTN-01 optical model (as in shower_candidate)."""
    return [
        "obdeect",
        "ray-tracing",
        "compile",
        "--simulation-models",
        inputs["models"],
        "--model",
        co.MODEL,
        "--version",
        co.VERSION,
        "--seed",
        co.SEED,
        "--output",
        model,
    ]


def envelope_pipeline(
    source,
    table,
    envelope,
    model,
    camera_map,
    sim_argv,
    threads,
    telescope=co.TELESCOPES[0],
):
    """One pipe: from-corsika | propagation | trace-stream | to-simtel | sim_telarray.

    envelope None gives the unthinned chain (weighted propagation, no scaled table). Only the
    telescope is propagated and traced. sim_argv reads standard input ("-") and is the last stage.
    """
    thinned = envelope is not None
    propagation = [
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
        "--telescopes",
        str(telescope),
    ]
    if thinned:
        propagation += ["--envelope", f"{ENVELOPE_FACTOR}={envelope}"]
    else:
        propagation += ["--mode", "weighted"]
    propagation += ["-i", "-", "-o", "-"]
    to_simtel = ["obdeect", "stream", "to-simtel", "--camera-map", camera_map]
    if thinned:
        to_simtel += ["--scaled-factor", ENVELOPE_FACTOR]
    to_simtel += ["-o", "-"]
    return [
        ("from_corsika", ["obdeect", "stream", "from-corsika", source]),
        ("propagation", propagation),
        (
            "trace_stream",
            [
                "obdeect",
                "ray-tracing",
                "trace-stream",
                "--model",
                model,
                "--seed",
                co.SEED,
                "--threads",
                threads,
                "--telescopes",
                str(telescope),
                "--pointing-zenith-deg",
                co.ZENITH_DEG,
                "--pointing-azimuth-deg",
                co.AZIMUTH_DEG,
                "--telescope-transmission",
                "full",
                "-i",
                "-",
                "-o",
                "-",
            ],
        ),
        ("to_simtel", to_simtel),
        ("sim_telarray", sim_argv),
    ]


def run_pipeline(stages, output, logs):
    """Run (name, argv) stages as one pipe: each stdout feeds the next stdin, the last writes output.

    All stages run at once. A stage's wall time runs from the start of the pipe to its exit; its
    CPU time and peak RSS are its own (wait4 rusage; sampled VmHWM on Linux). One pipeline record
    is appended to performance.RUSAGE: wall to the last exit, CPU summed, RSS maximal, and the
    stages in 'children'. A failed stage raises RuntimeError after the record is kept.
    """
    start = time.perf_counter()
    running = []
    previous = subprocess.DEVNULL
    for index, (name, argv) in enumerate(stages):
        last = index == len(stages) - 1
        argv = [str(arg) for arg in argv]
        with Path(logs[name]).open("wb") as log:
            if last:
                with Path(output).open("wb") as out:
                    proc = subprocess.Popen(
                        argv, stdin=previous, stdout=out, stderr=log
                    )
            else:
                proc = subprocess.Popen(
                    argv, stdin=previous, stdout=subprocess.PIPE, stderr=log
                )
        if previous is not subprocess.DEVNULL:
            previous.close()  # the parent keeps no copy of the pipe end
        sampler = (
            performance.RssSampler(proc.pid)
            if sys.platform.startswith("linux")
            else None
        )
        if sampler:
            sampler.start()
        running.append((name, argv, proc, sampler))
        previous = proc.stdout
    exits = {}

    def wait(name, proc, sampler):
        _, status, usage = os.wait4(proc.pid, 0)
        sampled_kb = sampler.stop() if sampler else None
        exits[name] = (time.perf_counter() - start, status, usage, sampled_kb)

    waiters = [
        threading.Thread(target=wait, args=(name, proc, sampler))
        for name, _, proc, sampler in running
    ]
    for waiter in waiters:
        waiter.start()
    for waiter in waiters:
        waiter.join()
    rss_divisor = 1024 * 1024 if sys.platform == "darwin" else 1024
    children = []
    for name, argv, proc, _ in running:
        wall, status, usage, sampled_kb = exits[name]
        proc.returncode = os.waitstatus_to_exitcode(status)
        children.append(
            {
                "name": name,
                "argv": argv,
                "wall_s": wall,
                "cpu_s": usage.ru_utime + usage.ru_stime,
                "peak_rss_mb": sampled_kb / 1024
                if sampled_kb is not None
                else usage.ru_maxrss / rss_divisor,
                "returncode": proc.returncode,
            }
        )
    record = {
        "argv": [child["argv"] for child in children],
        "wall_s": max(child["wall_s"] for child in children),
        "cpu_s": sum(child["cpu_s"] for child in children),
        "peak_rss_mb": max(child["peak_rss_mb"] for child in children),
        "children": children,
    }
    performance.RUSAGE.append(record)
    for child in children:
        if child["returncode"]:
            raise RuntimeError(
                f"{child['name']} exited {child['returncode']}; see {logs[child['name']]}"
            )
    return record


def simtel_candidate(work, inputs, files=None, name="candidate"):
    """sim_telarray argv of the candidate camera stage: stdin, Bypass_Optics=3 (camera only).

    The atmosphere is applied upstream by propagation, so sim_telarray gets the same config
    (CTAO.cfg) without an atmosphere override only through its mirror and QE settings.
    """
    extra = list(
        telescope_arguments(inputs["config_dir"], inputs["reference_atmosphere"])
    )
    if files is not None:
        extra = [f"-I{files['table_dir']}", *extra]
    argv, output = cs.simulate_argv(
        work,
        "-",
        name,
        bypass=3,
        config=str(inputs["config"]),
        extra=extra,
        binary=CANDIDATE_SIMTEL,
    )
    return argv, output


def envelope_candidate(work, source, inputs, files, threads=0):
    """Optimised candidate: photons thinned by the QE bound in propagation, one pipe to sim_telarray.

    Steps: compile of the optical model, then one pipe from-corsika | propagation | trace-stream |
    to-simtel | sim_telarray (Bypass_Optics=3). sim_telarray's statistics go to candidate.stdout.
    """
    model = work / "lst.model"
    sim_argv, _ = simtel_candidate(work, inputs, files)
    stages = envelope_pipeline(
        source,
        inputs["atmosphere"],
        files["envelope"],
        model,
        co.camera_map(work, inputs["rotation"]),
        sim_argv,
        threads,
    )
    logs = {name: work / f"{name}.log" for name, _ in stages}
    return [
        (
            "compile_model",
            lambda: run_timed(compile_argv(inputs, model), work / "compile.log"),
        ),
        (
            "pipeline",
            lambda: _checked_pipeline(stages, work / "candidate.stdout", logs),
        ),
    ]


def _checked_pipeline(stages, output, logs):
    """run_pipeline, then check the atmosphere that the piped sim_telarray read (its stdout)."""
    record = run_pipeline(stages, output, logs)
    # sim_telarray (the last stage) prints its setup on stderr, its statistics on stdout.
    last = stages[-1][0]
    text = Path(output).read_text(errors="replace") + Path(logs[last]).read_text(
        errors="replace"
    )
    common.check_atmosphere(text, ATMOSPHERE_FILE)
    return record


def pipeline_counts(work):
    """Photons into trace (telescope 1 rows, trace-stream summary) and into sim_telarray (All photons)."""
    counts = {}
    trace = work / "trace_stream.log"
    if trace.is_file():
        match = TRACE_SUMMARY.search(trace.read_text(errors="replace"))
        if match:
            counts.update(
                into_trace_photons=int(match[1]), out_of_trace_photons=int(match[2])
            )
    stdout = work / "candidate.stdout"
    if stdout.is_file():
        match = ALL_PHOTONS.search(stdout.read_text(errors="replace"))
        if match:
            counts["into_simtel_photons"] = int(match[1])
    return counts


def iact_telescope_weight(path, telescope=0):
    """Sum of the photon weights of one telescope (eventio index; CORSIKA telescope 1 is 0)."""
    total = 0.0
    with IACTFile(path) as stream:
        for event in stream:
            photons = event.photon_bunches.get(telescope)
            if photons is not None:
                total += float(np.sum(photons["photons"], dtype=np.float64))
    return total


def variant_counts(work, camera):
    """Photons into trace (telescope 1 rows in the trace-stream summary) and into sim_telarray (camera weight)."""
    counts = {"into_simtel_weight": iact_weight(camera)}
    for log in (work / "trace.log", work / "trace_stream.log"):
        if log.is_file():
            match = TRACE_SUMMARY.search(log.read_text(errors="replace"))
            if match:
                counts.update(
                    into_trace_photons=int(match[1]),
                    out_of_trace_photons=int(match[2]),
                )
    return counts


def child_metrics(prefix, summaries):
    """Median wall and CPU of each pipeline member over the repetitions: '<prefix>_<member>_wall_s'."""
    members = {}
    for summary in summaries:
        for child in summary.get("children", []):
            members.setdefault(child["name"], []).append(child)
    metrics = {}
    for name, records in members.items():
        metrics[f"{prefix}_{name}_wall_s"] = _median(records, "wall_s")
        metrics[f"{prefix}_{name}_cpu_s"] = _median(records, "cpu_s")
    return metrics


def physics_check(root, source, inputs, qe_path, headroom, metrics, details):
    """Camera comparison of the optimised chain with sim_telarray on the same input (report only).

    Normal atmosphere everywhere: reference = sim_telarray Bypass_Optics=0 with the config's
    atmosphere; candidates = obdeect propagation with the same transmission table. Controls:
    the unthinned pipe (weighted propagation, QE applied by sim_telarray).
    """
    work = root / "work" / "physics"
    work.mkdir(parents=True, exist_ok=True)
    reference = cs.camera_events(
        _reference_sim(
            work,
            source,
            inputs["config"],
            inputs["config_dir"],
            inputs["reference_atmosphere"],
            name="physics-reference",
        )
    )
    unthinned_dir = work / "unthinned"
    unthinned_dir.mkdir(exist_ok=True)
    unthinned = _pipe_camera(unthinned_dir, source, inputs, None)
    envelope_dir = work / "envelope"
    files = envelope_files(envelope_dir, qe_path, headroom)
    optimised = _pipe_camera(envelope_dir, source, inputs, files)
    details["counts"]["physics_envelope"] = pipeline_counts(envelope_dir)
    shower_ids = [str(event.event_number) for event in IACTFile(source)]
    pairs = {
        "unthinned": (reference, unthinned),
        "envelope": (reference, optimised),
        "envelope_vs_unthinned": (unthinned, optimised),
    }
    details["physics"] = {}
    for label, (left, right) in pairs.items():
        out = root / "physics" / label
        out.mkdir(parents=True, exist_ok=True)
        comparison = cs.compare_camera(left, right, out)
        triggers = co.trigger_difference(left, right, shower_ids)
        centroids = co.centroid_offsets(left, right)
        metrics.update({f"physics_{label}_{k}": v for k, v in comparison.items()})
        metrics.update({f"physics_{label}_{k}": v for k, v in centroids.items()})
        for tel, row in triggers.items():
            metrics[f"physics_{label}_trigger_difference_tel{tel}"] = row["difference"]
            metrics[f"physics_{label}_trigger_difference_se_tel{tel}"] = row["se"]
        details["physics"][label] = {"triggers": triggers}


def _pipe_camera(work, source, inputs, files):
    """Run one pipe (thinned when files is given) and return the camera events of its sim_telarray."""
    model = work / "lst.model"
    common.command(compile_argv(inputs, model), work / "compile.log")
    sim_argv, output = simtel_candidate(work, inputs, files)
    stages = envelope_pipeline(
        source,
        inputs["atmosphere"],
        files["envelope"] if files else None,
        model,
        co.camera_map(work, inputs["rotation"]),
        sim_argv,
        0,
    )
    logs = {name: work / f"{name}.log" for name, _ in stages}
    _checked_pipeline(stages, work / "candidate.stdout", logs)
    return cs.camera_events(output)


def run_air_shower(root, source, config, reps, metrics, details):
    config = config.resolve()
    config_dir = config.parent
    cfg_text = (config_dir / "CTAO-LSTN-01.cfg").read_text()
    inputs = {
        "config": config,
        "config_dir": config_dir,
        "rotation": co.camera_rotation_deg(config_dir, cfg_text),
        "unit": root / "unit-atmosphere.dat",
        "table": root / "transmission-below-1000nm.dat",
        "models": str(common.models_root()),
        "qe": config_dir / co.cfg_value(cfg_text, "quantum_efficiency"),
        "qe_variation": float(co.cfg_value(cfg_text, "qe_variation")),
        "atmosphere": root / "atmosphere-below-1000nm.dat",
        "reference_atmosphere": config_dir / ATMOSPHERE_FILE,
    }
    shutil.copyfile(co.UNIT_ATMOSPHERE, inputs["unit"])
    performance.restrict_table(performance.TRANSMISSION, inputs["table"])
    # Propagation needs integer wavelengths below 1000 nm, so the table is cut there; the reference
    # keeps the full table (its wavelengths are never above 1000 nm in the sim_telarray band).
    performance.restrict_table(inputs["reference_atmosphere"], inputs["atmosphere"])
    headroom = headroom_for(inputs["qe_variation"])
    details["settings"].update(
        {
            "input": str(source),
            "config": str(config),
            "reference_atmosphere": co.cfg_value(cfg_text, "atmospheric_transmission"),
            "candidate_atmosphere": "obdeect propagation (transmission table, 250-700 nm); sim_telarray unit table",
            "telescopes": list(co.TELESCOPES),
            "pointing": f"zenith {co.ZENITH_DEG} deg, azimuth {co.AZIMUTH_DEG} deg",
        }
    )
    reference_dir = root / "work" / "reference"
    reference_dir.mkdir(parents=True, exist_ok=True)
    reference_ok = _run(
        "reference",
        shower_reference(
            reference_dir, source, config, config_dir, inputs["reference_atmosphere"]
        ),
        reps,
        metrics,
        details,
    )
    candidates = []
    for name, threads in THREAD_VARIANTS:
        label = f"candidate_{name}"
        work = root / "work" / name
        work.mkdir(parents=True, exist_ok=True)
        if _run(
            label,
            shower_candidate(work, source, inputs, threads),
            reps,
            metrics,
            details,
        ):
            details["counts"][label] = variant_counts(work, work / "candidate.iact")
        details["threads"][label] = threads
        candidates.append(label)
    details["counts"]["reference"] = {
        "into_simtel_weight": iact_telescope_weight(source)
    }
    env_label = "candidate_envelope"
    env_work = root / "work" / "envelope"
    files = envelope_files(env_work, inputs["qe"], headroom)
    if _run(
        env_label,
        envelope_candidate(env_work, source, inputs, files, threads=0),
        reps,
        metrics,
        details,
    ):
        metrics.update(
            child_metrics(env_label, details["stages"][env_label]["pipeline"])
        )
        details["counts"][env_label] = pipeline_counts(env_work)
    metrics["candidate_envelope_headroom"] = headroom
    details["threads"][env_label] = 0
    candidates.append(env_label)
    try:
        physics_check(root, source, inputs, inputs["qe"], headroom, metrics, details)
        metrics["physics_completed"] = 1
    except (RuntimeError, OSError, ValueError) as error:  # reported, not hidden
        metrics["physics_completed"] = 0
        details["failures"]["physics"] = str(error)[:400]
    _ratios(metrics, reference_ok, candidates)


def run_source(root, name, reps, metrics, details):
    spec = SOURCES[name]
    reference_dir = root / "work" / "reference"
    candidate_dir = root / "work" / "candidate"
    reference_dir.mkdir(parents=True, exist_ok=True)
    candidate_dir.mkdir(parents=True, exist_ok=True)
    details["settings"].update(
        {
            "events": spec["events"],
            "photons_per_event": spec["photons"],
            "bypass_optics": spec["bypass"],
        }
    )
    reference_ok = _run(
        "reference", source_reference(reference_dir, spec), reps, metrics, details
    )
    _run(
        "candidate_threads1",
        source_candidate(candidate_dir, spec),
        reps,
        metrics,
        details,
    )
    details["threads"]["candidate_threads1"] = (
        "1 (all stages single-threaded; default is the same)"
    )
    _ratios(metrics, reference_ok, ["candidate_threads1"])


def _plot(path, stages_by_pipeline):
    """Stacked median wall times per pipeline, one colour per stage name, stages in sequence."""
    from obdeect_validation.adapters import plotting

    fig, axis = plotting.figure(1, 1, (9, 1.2 + 0.8 * len(stages_by_pipeline)))
    colours, seen = {}, set()
    for row, (label, stages) in enumerate(stages_by_pipeline.items()):
        left = 0.0
        for name, wall in stages.items():
            colour = colours.setdefault(
                name, plotting.CYCLE[len(colours) % len(plotting.CYCLE)]
            )
            axis.barh(
                row,
                wall,
                left=left,
                color=colour,
                height=0.6,
                label=None if name in seen else name,
            )
            seen.add(name)
            left += wall
    axis.set_yticks(range(len(stages_by_pipeline)), list(stages_by_pipeline))
    axis.invert_yaxis()
    axis.set_xlabel("wall seconds (median over repetitions; stages run in sequence)")
    axis.grid(axis="y", visible=False)
    axis.legend(fontsize=8)
    plotting.save(fig, path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", choices=["air-shower", *SOURCES])
    parser.add_argument(
        "--config", type=Path, help="sim_telarray config from simtools (air-shower)"
    )
    parser.add_argument("--input", type=Path, default=performance.INPUT)
    parser.add_argument("--reps", type=int, default=REPS)
    args = parser.parse_args()
    if args.case == "air-shower" and args.config is None:
        parser.error("air-shower needs --config")
    root = Path.cwd().resolve()
    metrics = {}
    details = {
        "case": args.case,
        "host": performance.host_info(),
        "settings": {
            "reps": args.reps,
            "timing": "wall and CPU: wait4 rusage per child; peak RSS: sampled /proc VmHWM on Linux",
        },
        "threads": {},
        "stages": {},
        "pipelines": {},
        "counts": {},
        "failures": {},
    }
    if args.case == "air-shower":
        run_air_shower(
            root, args.input.resolve(), args.config, args.reps, metrics, details
        )
    else:
        run_source(root, args.case, args.reps, metrics, details)
    walls = {
        label: {name: _median(records, "wall_s") for name, records in per_stage.items()}
        for label, per_stage in details["stages"].items()
    }
    if walls:
        _plot(root / "wall-by-stage.png", walls)
    common.write_results(metrics, details, root)
    shutil.rmtree(root / "work", ignore_errors=True)


if __name__ == "__main__":
    main()
