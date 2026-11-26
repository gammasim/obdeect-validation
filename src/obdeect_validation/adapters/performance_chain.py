"""Matched end-to-end timing: obdeect pipelines against their reference programs, with budgets.

Each case times a reference pipeline and an obdeect pipeline on the same physical input and
photon counts. The reference is CORSIKA IACT -> sim_telarray with full optics (air shower) or
LightEmission -> sim_telarray (sources). The candidate is the obdeect chain ending in the
sim_telarray camera stage (Bypass_Optics=1 for the air shower, as in cases/chain/optics-bypass1.json).
Every child process is timed with performance.run_timed. A pipeline total is the sum of its stage
wall and CPU times and the maximum of its stage peak RSS (stages run one after another).
The repetitions alternate between pipelines (rep 1 of each, then rep 2 of each), so a change in host
state affects all pipelines alike (run_interleaved).

The air-shower case also times an optimised candidate (candidate_envelope, obdeect default threads;
candidate_envelope_threads1 for the single-core comparison): from-corsika, propagation (QE bound,
thinned), trace-stream and to-simtel run as one pipe; sim_telarray then reads a QE / bound table, so
the product is the model QE. A physics check compares its camera output with sim_telarray on the
same input (unit atmosphere on both sides); its photoelectron bias is the equal-statistics condition
of the wall and CPU budgets (performance_budget). candidate_unthinned_threads1 is the same pipe
without thinning (report only: thinning and piping gains).

Only ray-tracing trace-stream is multithreaded. The air-shower candidate therefore runs with
threads=1 and with the obdeect default (0 = every core); the source pipelines are single-threaded
and run once. Budgets (performance_budget): *_budget_exceeded metrics, 0 = pass.

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.performance_chain air-shower --config cfg/model/7.0.0/CTAO.cfg [--events N] [--reps N]
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
import tempfile
import threading
import time
import traceback
from pathlib import Path

import numpy as np
from eventio.iact import IACTFile

from obdeect_validation.adapters import chain_optics as co
from obdeect_validation.adapters import chain_showers as cs
from obdeect_validation.adapters import (
    common,
    performance,
    performance_budget,
    source_camera,
)
from obdeect_validation.adapters import obdeect_cli as cli
from obdeect_validation.adapters.performance import iact_weight
from obdeect_validation.adapters.performance import run_timed as _run_timed

REPS = 3
# Production atmosphere of the LSTN-01 model (CTAO.cfg, set only under #if TELESCOPE == 0). The
# reference and every candidate pass it explicitly and the case checks the sim_telarray log for it.
ATMOSPHERE_FILE = "atmospheric_transmission-2.0.0-OBS-North.dat"
ENVELOPE_FACTOR = "qe"
# An efficiency bound must be positive: bins with QE 0 keep this bound (their table value is 0).
ENVELOPE_FLOOR = 1e-6
# The table maximum (headroom) is not a constant: it is derived from the model's qe_variation by
# headroom_for(). sim_telarray's first-cut probability is p2 * max_qe_rel with max_qe_rel up to
# 1 + 6 * qe_variation; a table maximum above 1 / (1 + 6 * qe_variation) saturates for upward
# pixel fluctuations and loses photons (clustered bias -6.6 % at headroom 1 in LSTN-01).
# Second node of each nm bin: the linear bound is constant on [n, n + 1) nm.
ENVELOPE_NODE_WIDTH = 0.999999
# CORSIKA Cherenkov band (EVTH words 96-97) on sim_telarray's integer nm bins.
WAVELENGTH_GRID = np.arange(100, 2001)
TRACE_SUMMARY = re.compile(r"trace-stream: \d+ events, (\d+) photons in, (\d+) out")
ALL_PHOTONS = re.compile(r"All photons:\s+(\d+)")
# sim_telarray for the candidate camera stage (Bypass_Optics=3 build); None: the production binary.
CANDIDATE_SIMTEL = os.environ.get("SIMTEL_CANDIDATE_BINARY") or None
# Extra sim_telarray arguments of the candidate camera stage only (chain_mode3 adds "-v": sim_telarray
# logs the QE table file only at verbose level; verbosity does not change the simulation).
CANDIDATE_OPTIONS = ()
THREAD_VARIANTS = (
    ("threads1", 1),
    ("threadsdefault", 0),
)  # 0 = obdeect default, every core
# Envelope candidate: the default-threads label is the legacy name (candidate_envelope_completed).
ENVELOPE_VARIANTS = (
    ("candidate_envelope", 0),
    ("candidate_envelope_threads1", 1),
)
PHOTON_KEYS = ("wall_s", "cpu_s", "peak_rss_mb")
UNTHINNED_LABEL = "candidate_unthinned_threads1"
# Timing reliability (report only, see timing_reliability): relative wall spread of a pipeline total
# over the repetitions, and the host 1-minute load average before the run.
SPREAD_LIMIT = 0.05
LOAD_LIMIT = 1.0
# The model compile runs once per repetition: left out of CPU shares and of the piping gain.
SHARE_EXCLUDED = ("compile_model",)
# obdeect wrapper subcommands that run as native programs in candidate pipelines. The wrapper adds
# ~0.1 s per call (light-start-up measurements); the native program takes the same arguments and
# gives byte-identical output (checked for propagation, to-simtel, from-corsika with the same seeds).
# trace-stream and the model compile have no native program here: they stay on the wrapper.
# obdeect sources run has no native program (the sources tools were replaced by run descriptions),
# so the sources stage uses the wrapper.
NATIVE_PROGRAMS = {
    ("propagation", "run"): ("obdeect.propagation", "obdeect-propagation-run"),
    ("stream", "to-simtel"): ("obdeect.stream", "obdeect-stream-to-simtel"),
    ("stream", "from-corsika"): ("obdeect.stream", "obdeect-stream-from-corsika"),
    ("ray-tracing", "trace-stream"): (
        "obdeect.ray_tracing",
        "obdeect-ray-tracing-trace-stream",
    ),
}


def native_argv(argv):
    """argv with an obdeect wrapper subcommand replaced by its native program (same arguments).

    Other commands, and subcommands without a native program, are returned unchanged.
    """
    argv = list(argv)
    if (
        len(argv) < 3
        or argv[0] != "obdeect"
        or (argv[1], argv[2]) not in NATIVE_PROGRAMS
    ):
        return argv
    from obdeect.stream.cli import program_path

    package, program = NATIVE_PROGRAMS[(argv[1], argv[2])]
    return [program_path(package, program), *argv[3:]]


def run_timed(argv, log, **kwargs):
    """performance.run_timed of the native form of argv (see native_argv)."""
    return _run_timed(native_argv(argv), log, **kwargs)


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
        "candidate": lambda n, p: cli.obdeect(
            "sources",
            "run",
            "flasher",
            type="dish",
            observation_level="2150",
            seed="1",
            atmosphere_id="1",
            position="0,0,16",
            camera_radius="1.2",
            spectrum="fixed:wavelength_nm=400",
            angular="isotropic",
            events=n,
            photons=p,
            extra=(
                "--set",
                "run.frame=sky",
                "--set",
                "events[0].emissions[0].to_plane_z0=false",
            ),
        ),
        "reference": lambda n, p: ["ff-1m", "--events", n, "--photons", p],
        "bypass": 1,  # as in source_camera stage 2
    },
    "point-source": {
        "events": 2,
        "photons": 10_000_000,
        "candidate": lambda n, p: cli.obdeect(
            "sources",
            "run",
            "point-source",
            observation_level="2150",
            seed="17",
            atmosphere_id="1",
            events=n,
            photons=p,
            position="0,0,1000",
            telescope_radius="15",
            spectrum="fixed:wavelength_nm=400",
            angular="isotropic",
        ),
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


#: sim_telarray skips a telescope event below min_photons (sim_telarray.c, "photons <
#: min_photons"), counting the photons of the bunches it reads. A thinned stream carries one photon
#: per surviving bunch, so the configured threshold (120 LST, 450 MST) is not comparable with the
#: reference's raw photon count. Both sides therefore run with 0 and simulate the same events.
MIN_PHOTONS = 0


def telescope_arguments(config_dir, atmosphere=None):
    """sim_telarray arguments of the matched one-telescope model (as in chain_optics)."""
    extra = [
        f"-I{config_dir}",
        "-DTELESCOPE=1",
        "-C",
        f"maximum_telescopes={len(co.TELESCOPES)}",
        "-C",
        "trigger_telescopes=1",  # one telescope: the default array trigger needs two
        "-C",
        f"min_photons={MIN_PHOTONS}",
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
    """Medians over repetitions, the range of wall and CPU, and the relative wall spread."""
    wall = [record["wall_s"] for record in records]
    cpu = [record["cpu_s"] for record in records]
    out = {
        f"{prefix}_wall_s": _median(records, "wall_s"),
        f"{prefix}_cpu_s": _median(records, "cpu_s"),
        f"{prefix}_peak_rss_mb": _median(records, "peak_rss_mb"),
        f"{prefix}_wall_s_min": min(wall),
        f"{prefix}_wall_s_max": max(wall),
        f"{prefix}_cpu_s_min": min(cpu),
        f"{prefix}_cpu_s_max": max(cpu),
    }
    spread = ratio(max(wall) - min(wall), out[f"{prefix}_wall_s"])
    if spread is not None:
        out[f"{prefix}_wall_s_spread"] = spread
    return out


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


def ratio_min_metrics(label, metrics):
    """Candidate / reference of the minimum wall and CPU over the repetitions (pipeline totals)."""
    out = {}
    for key, name in (("wall_s", "wall"), ("cpu_s", "cpu")):
        candidate = metrics.get(f"{label}_{key}_min")
        reference = metrics.get(f"reference_{key}_min")
        if candidate is None or reference is None:
            continue
        value = ratio(candidate, reference)
        if value is not None:
            out[f"ratio_{name}_{label.removeprefix('candidate_')}_min"] = value
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
                cli.from_corsika(source),
                work / "from-corsika.log",
                stdout_file=native,
            ),
        ),
        (
            "propagation",
            lambda: run_timed(
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
            ),
        ),
        (
            "compile_model",
            lambda: run_timed(compile_argv(inputs, model), work / "compile.log"),
        ),
        (
            "trace_stream",
            lambda: run_timed(
                cli.trace_stream(
                    model,
                    propagated,
                    traced,
                    seed=co.SEED,
                    threads=threads,
                    telescopes=",".join(map(str, co.TELESCOPES)),
                    pointing_zenith_deg=co.ZENITH_DEG,
                    pointing_azimuth_deg=co.AZIMUTH_DEG,
                    telescope_transmission="full",
                ),
                work / "trace.log",
            ),
        ),
        (
            "to_simtel",
            lambda: run_timed(
                cli.obdeect(
                    "stream",
                    "to-simtel",
                    camera_map=camera_map,
                    o=camera_stream,
                    extra=(traced,),
                ),
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
                cli.obdeect(
                    "propagation", "run", *SKY_PROPAGATION, i=native, o=propagated
                ),
                work / "propagation.log",
            ),
        ),
        (
            "to_simtel",
            lambda: run_timed(
                cli.obdeect("stream", "to-simtel", o=candidate, extra=(propagated,)),
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
    """LightEmission -> sim_telarray with the same photon counts and camera options.

    The photon file is uncompressed, as the obdeect side's stream: gzip would add its own time."""
    reference = work / "reference.dat"
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


def run_interleaved(runs, reps, metrics, details):
    """Measure several pipelines with their repetitions alternating: rep 1 of each, then rep 2 of each.

    runs: [(label, steps)]. A label that raises is recorded as failed and skipped in the later
    repetitions; the others still run. Returns {label: completed}.
    """
    per_stage = {label: {name: [] for name, _ in steps} for label, steps in runs}
    totals = {label: [] for label, _ in runs}
    failed = {}
    for _ in range(reps):
        for label, steps in runs:
            if label in failed:
                continue
            try:
                rep_stages, rep_totals = measure(steps, 1)
            except (RuntimeError, OSError, ValueError) as error:  # reported, not hidden
                failed[label] = (error, traceback.format_exc())
                continue
            for name, records in rep_stages.items():
                per_stage[label][name].extend(records)
            totals[label].extend(rep_totals)
    ok = {}
    for label, _ in runs:
        if label in failed:
            error, trace = failed[label]
            metrics[f"{label}_completed"] = 0
            details["failures"][label] = {
                "error": str(error)[:400],
                "traceback": trace[-800:],
            }
            ok[label] = False
            continue
        metrics.update(pipeline_metrics(label, per_stage[label], totals[label]))
        details["stages"][label] = per_stage[label]
        details["pipelines"][label] = totals[label]
        ok[label] = True
    return ok


def _run(label, steps, reps, metrics, details):
    """Measure one pipeline (run_interleaved with one label); a failure is recorded."""
    return run_interleaved([(label, steps)], reps, metrics, details)[label]


def _ratios(metrics, reference_ok, candidates):
    """Candidate / reference ratios of the pipeline totals, medians and minima (only when both completed)."""
    if not reference_ok:
        return
    ref = {key: metrics[f"reference_{key}"] for key in PHOTON_KEYS}
    for label in candidates:
        if metrics.get(f"{label}_completed") != 1:
            continue
        cand = {key: metrics[f"{label}_{key}"] for key in PHOTON_KEYS}
        metrics.update(ratio_metrics(label.removeprefix("candidate_"), cand, ref))
        metrics.update(ratio_min_metrics(label, metrics))


def sim_telarray_qe(path, grid=WAVELENGTH_GRID):
    """QE as sim_telarray reads it: first two columns, linear in wavelength, at integer nm."""
    table = np.loadtxt(path, comments="#", usecols=(0, 1), ndmin=2)
    return np.interp(grid, table[:, 0], table[:, 1])


def bin_envelope(qe, headroom):
    """Efficiency bound per nm bin (QE / headroom, at least ENVELOPE_FLOOR) and the table QE / bound.

    The bound is at least the QE in every bin, and the table is at most headroom (see headroom_for),
    so bound times table is the QE again.
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


def envelope_files(work, qe_path, headroom):
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
    return cli.compile_model(
        model,
        co.MODEL,
        co.VERSION,
        seed=co.SEED,
        simulation_models=inputs["models"],
    )


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

    envelope None gives the unthinned chain (weighted propagation, no scaled table). Propagation
    and trace-stream keep only the telescope; the photons of the others are not processed.
    sim_argv reads standard input ("-") and is the last stage.
    """
    thinned = envelope is not None
    propagation = cli.propagation_run(
        "-",
        "-",
        transmission=table,
        site_altitude=2156,
        pathlen=28,
        seed=123,
        telescopes=str(telescope),
        envelope=f"{ENVELOPE_FACTOR}={envelope}" if thinned else None,
        mode=None if thinned else "weighted",
    )
    to_simtel = cli.obdeect(
        "stream",
        "to-simtel",
        camera_map=camera_map,
        scaled_factor=ENVELOPE_FACTOR if thinned else None,
        o="-",
    )
    return [
        ("from_corsika", cli.from_corsika(source)),
        ("propagation", propagation),
        (
            "trace_stream",
            cli.trace_stream(
                model,
                "-",
                "-",
                seed=co.SEED,
                threads=threads,
                telescopes=str(telescope),
                pointing_zenith_deg=co.ZENITH_DEG,
                pointing_azimuth_deg=co.AZIMUTH_DEG,
                telescope_transmission="full",
            ),
        ),
        ("to_simtel", to_simtel),
        ("sim_telarray", sim_argv),
    ]


def _abort_pipeline(procs, samplers):
    """Kill and reap the children of a pipe whose setup failed; stop their samplers."""
    for proc in procs:
        proc.kill()
        proc.wait()
    for sampler in samplers:
        sampler.stop()


def wait_for_children(running, start: float) -> dict:
    """Reap the children of running ((name, argv, proc, sampler) tuples) concurrently.

    Returns {name: (wall seconds since start, wait status, rusage, sampled peak RSS kB or None)}.
    """
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
    return exits


def run_pipeline(stages, output, logs):
    """Run (name, argv) stages as one pipe: each stdout feeds the next stdin, the last writes output.

    All stages run at once. A stage's wall time runs from the start of the pipe to its exit; its
    CPU time and peak RSS are its own (wait4 rusage; sampled VmHWM on Linux). One pipeline record
    is appended to performance.RUSAGE: wall to the last exit, CPU summed, RSS maximal, and the
    stages in 'children'. A failed stage raises RuntimeError after the record is kept.
    """
    start = time.perf_counter()
    running = []
    procs, samplers = [], []  # everything started, for cleanup if the setup fails
    previous = subprocess.DEVNULL
    try:
        for index, (name, argv) in enumerate(stages):
            last = index == len(stages) - 1
            argv = native_argv([str(arg) for arg in argv])
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
            procs.append(proc)
            if previous is not subprocess.DEVNULL:
                previous.close()  # the parent keeps no copy of the pipe end
            sampler = (
                performance.RssSampler(proc.pid)
                if sys.platform.startswith("linux")
                else None
            )
            if sampler:
                sampler.start()
                samplers.append(sampler)
            running.append((name, argv, proc, sampler))
            previous = proc.stdout
    except BaseException:
        _abort_pipeline(procs, samplers)
        raise
    exits = wait_for_children(running, start)
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
    extra = [
        *CANDIDATE_OPTIONS,
        *telescope_arguments(inputs["config_dir"], inputs["reference_atmosphere"]),
    ]
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
    files None gives the unthinned pipe (weighted propagation, no scaled table).
    """
    model = work / "lst.model"
    sim_argv, _ = simtel_candidate(work, inputs, files)
    stages = envelope_pipeline(
        source,
        inputs["atmosphere"],
        files["envelope"] if files else None,
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
    """Median exit time and CPU of each pipeline member over the repetitions.

    The exit time runs from the start of the pipe to the member's exit ('<prefix>_<member>_exit_s'),
    not the member's own run time.
    """
    members = {}
    for summary in summaries:
        for child in summary.get("children", []):
            members.setdefault(child["name"], []).append(child)
    metrics = {}
    for name, records in members.items():
        metrics[f"{prefix}_{name}_exit_s"] = _median(records, "wall_s")
        metrics[f"{prefix}_{name}_cpu_s"] = _median(records, "cpu_s")
    return metrics


def timing_reliability(metrics, labels):
    """(1, reasons) when a pipeline total is too noisy to compare, else (0, reasons).

    Noisy: the relative wall spread over the repetitions above SPREAD_LIMIT, or a host load
    above LOAD_LIMIT before the run. Report only; no budget depends on it.
    """
    reasons = []
    for label in labels:
        spread = metrics.get(f"{label}_wall_s_spread")
        if spread is not None and spread > SPREAD_LIMIT:
            reasons.append(f"{label} wall spread {spread:.3f} > {SPREAD_LIMIT}")
    before = metrics.get("host_loadavg_before")
    if before is not None and before > LOAD_LIMIT:
        reasons.append(f"host_loadavg_before {before:.2f} > {LOAD_LIMIT}")
    details = {
        "spread_limit": SPREAD_LIMIT,
        "load_limit": LOAD_LIMIT,
        "reasons": reasons,
    }
    return int(bool(reasons)), details


def _shares(label, cpu, denominator):
    """'<label>_<stage>_cpu_share' = stage CPU / denominator; nothing when the denominator is not positive."""
    out = {}
    for name, value in cpu.items():
        share = ratio(value, denominator)
        if share is not None:
            out[f"{label}_{name}_cpu_share"] = share
    return out


def normalised_metrics(metrics, details):
    """Cost per shower and per telescope photon of each completed pipeline, and CPU shares of stages.

    Per shower: wall and CPU over settings.showers. Per photon: CPU in ns over telescope_photons.
    Shares: candidate_threads1 stages over its CPU without compile_model; envelope pipe members over
    the pipe stage (their CPU sums to it).
    """
    out = {}
    showers = details["settings"].get("showers")
    photons = metrics.get("telescope_photons")
    for label in details["pipelines"]:
        wall, cpu = metrics[f"{label}_wall_s"], metrics[f"{label}_cpu_s"]
        out[f"{label}_s_per_shower"] = ratio(wall, showers)
        out[f"{label}_cpu_s_per_shower"] = ratio(cpu, showers)
        out[f"{label}_cpu_ns_per_photon"] = ratio(cpu * 1e9, photons)
    out = {key: value for key, value in out.items() if value is not None}
    if "candidate_threads1" in details["stages"]:
        names = [
            name
            for name in details["stages"]["candidate_threads1"]
            if name not in SHARE_EXCLUDED
        ]
        cpu = {name: metrics[f"candidate_threads1_{name}_cpu_s"] for name in names}
        out.update(_shares("candidate_threads1", cpu, sum(cpu.values())))
    for label in ("candidate_envelope", "candidate_envelope_threads1"):
        if label not in details["stages"]:
            continue
        members = {
            child["name"]
            for summary in details["stages"][label]["pipeline"]
            for child in summary.get("children", [])
        }
        cpu = {name: metrics[f"{label}_{name}_cpu_s"] for name in sorted(members)}
        out.update(_shares(label, cpu, metrics[f"{label}_pipeline_cpu_s"]))
    return out


def reference_optics_metrics(metrics, details):
    """Optics cost of the reference, by difference of two sim_telarray runs (report only).

    The reference runs sim_telarray with full optics; candidate_threads1 runs the same events
    camera only (Bypass_Optics=1). Assumption: the camera work is the same in both runs, so the
    difference is the optics. Missing inputs give no metric.
    """
    out = {}
    for key in ("cpu_s", "wall_s"):
        reference = metrics.get(f"reference_sim_telarray_{key}")
        camera = metrics.get(f"candidate_threads1_sim_telarray_{key}")
        if reference is not None and camera is not None:
            out[f"reference_optics_{key}"] = reference - camera
    if "reference_optics_cpu_s" in out:
        out["reference_camera_cpu_s"] = metrics["candidate_threads1_sim_telarray_cpu_s"]
    trace = metrics.get("candidate_threads1_trace_stream_cpu_s")
    optics = out.get("reference_optics_cpu_s")
    if trace is None or optics is None:
        details["reference_optics"] = {
            "trace_stream_ratio": "omitted: inputs not measured"
        }
    elif ratio(trace, optics) is None:
        details["reference_optics"] = {
            "trace_stream_ratio": "omitted: reference_optics_cpu_s <= 0"
        }
    else:
        out["trace_stream_vs_reference_optics_cpu_ratio"] = ratio(trace, optics)
    return out


def _without_compile(metrics, label, key):
    """Pipeline wall or CPU of a label without its compile_model stage, or None if not measured."""
    total = metrics.get(f"{label}_{key}")
    compile_time = metrics.get(f"{label}_compile_model_{key}")
    if total is None or compile_time is None:
        return None
    return total - compile_time


def gain_metrics(metrics):
    """Thinning gain: unthinned / thinned envelope pipe (threads=1).

    Piping gain: separate stages (candidate_threads1) / unthinned pipe. Both run weighted
    propagation, so the separate stages are compared with the pipe. All without compile_model.
    """
    out = {}
    for name in ("cpu", "wall"):
        piped = _without_compile(metrics, UNTHINNED_LABEL, f"{name}_s")
        thinned = _without_compile(metrics, "candidate_envelope_threads1", f"{name}_s")
        if piped is not None and thinned is not None:
            out[f"thinning_gain_{name}"] = ratio(piped, thinned)
        separate = _without_compile(metrics, "candidate_threads1", f"{name}_s")
        if separate is not None and piped is not None:
            out[f"piping_gain_{name}"] = ratio(separate, piped)
    return {key: value for key, value in out.items() if value is not None}


def physics_check(
    root, source, inputs, qe_path, headroom, metrics, details, work_root=None
):
    """Camera comparison of the optimised chain with sim_telarray on the same input (report only).

    Normal atmosphere everywhere: reference = sim_telarray Bypass_Optics=0 with the config's
    atmosphere; candidates = obdeect propagation with the same atmosphere. Controls:
    the unthinned pipe (weighted propagation, QE applied by sim_telarray).
    """
    work = (work_root or root / "work") / "physics"
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
    shower_ids = co.camera_keys(source)
    pairs = {
        "unthinned": (reference, unthinned),
        "envelope": (reference, optimised),
        "envelope_vs_unthinned": (unthinned, optimised),
    }
    details["physics"] = {}
    for label, (left, right) in pairs.items():
        out = root / "physics" / label
        out.mkdir(parents=True, exist_ok=True)
        comparison = cs.compare_camera(left, right, out, co.CAMERA_NAMES)
        triggers = co.trigger_difference(left, right, shower_ids)
        centroids = co.centroid_offsets(left, right)
        metrics.update({f"physics_{label}_{k}": v for k, v in comparison.items()})
        metrics.update({f"physics_{label}_{k}": v for k, v in centroids.items()})
        for tel, row in triggers.items():
            metrics[f"physics_{label}_trigger_difference_tel{tel}"] = row["difference"]
            metrics[f"physics_{label}_trigger_difference_se_tel{tel}"] = row["se"]
        details["physics"][label] = {"triggers": triggers}
    return {
        "reference": reference,
        "unthinned": unthinned,
        "envelope": optimised,
        "envelope_dir": envelope_dir,
        "files": files,
        "work": work,
        "shower_ids": shower_ids,
    }


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


def _camera_only_lines(work):
    """Camera-only lines of the envelope pipe's sim_telarray (stdout and log), see performance_budget."""
    text = "".join(
        (work / name).read_text(errors="replace")
        for name in ("candidate.stdout", "sim_telarray.log")
        if (work / name).is_file()
    )
    return performance_budget.count_camera_only(text)


def prepare_inputs(root, config):
    """Inputs of the LSTN-01 air-shower chains from the simtools sim_telarray config (CTAO.cfg).

    Returns (inputs, cfg_text); the unit table and the cut atmosphere are written to root.
    The propagation atmosphere is the production table cut below 1000 nm (propagation needs
    integer wavelengths there); the reference keeps the full table.
    """
    config = config.resolve()
    config_dir = config.parent
    cfg_text = (config_dir / "CTAO-LSTN-01.cfg").read_text()
    inputs = {
        "config": config,
        "config_dir": config_dir,
        "rotation": co.camera_rotation_deg(config_dir, cfg_text),
        "unit": root / "unit-atmosphere.dat",
        "models": str(common.models_root()),
        "qe": config_dir / co.cfg_value(cfg_text, "quantum_efficiency"),
        "qe_variation": float(co.cfg_value(cfg_text, "qe_variation")),
        "atmosphere": root / "atmosphere-below-1000nm.dat",
        "reference_atmosphere": config_dir / ATMOSPHERE_FILE,
    }
    inputs["unit"] = common.write_zero_depth_table(root, co.UNIT_ATMOSPHERE)
    performance.restrict_table(inputs["reference_atmosphere"], inputs["atmosphere"])
    return inputs, cfg_text


def run_air_shower(root, work_root, source, config, reps, metrics, details):
    inputs, cfg_text = prepare_inputs(work_root, config)
    config_dir = inputs["config_dir"]
    config = inputs["config"]
    headroom = headroom_for(inputs["qe_variation"])
    details["settings"].update(
        {
            "input": str(source),
            "config": str(config),
            "reference_atmosphere": co.cfg_value(cfg_text, "atmospheric_transmission"),
            "candidate_atmosphere": "obdeect propagation (reference atmosphere below 1000 nm, 250-700 nm); sim_telarray unit table",
            "telescopes": list(co.TELESCOPES),
            "pointing": f"zenith {co.ZENITH_DEG} deg, azimuth {co.AZIMUTH_DEG} deg",
            "candidate_sim_telarray": CANDIDATE_SIMTEL
            or "production (SIMTEL_CANDIDATE_BINARY unset)",
        }
    )
    reference_dir = work_root / "reference"
    reference_dir.mkdir(parents=True, exist_ok=True)
    # Steps are built before any run (building runs nothing); repetitions alternate between labels.
    runs = [
        (
            "reference",
            shower_reference(
                reference_dir,
                source,
                config,
                config_dir,
                inputs["reference_atmosphere"],
            ),
        )
    ]
    for name, threads in THREAD_VARIANTS:
        work = work_root / name
        work.mkdir(parents=True, exist_ok=True)
        runs.append(
            (f"candidate_{name}", shower_candidate(work, source, inputs, threads))
        )
    files = envelope_files(work_root / "envelope", inputs["qe"], headroom)
    for env_label, threads in ENVELOPE_VARIANTS:
        env_work = work_root / env_label.removeprefix("candidate_")
        env_work.mkdir(parents=True, exist_ok=True)
        runs.append(
            (
                env_label,
                envelope_candidate(env_work, source, inputs, files, threads=threads),
            )
        )
    unthinned_work = work_root / "unthinned"
    unthinned_work.mkdir(parents=True, exist_ok=True)
    runs.append(
        (
            UNTHINNED_LABEL,
            envelope_candidate(unthinned_work, source, inputs, None, threads=1),
        )
    )
    ok = run_interleaved(runs, reps, metrics, details)
    reference_ok = ok["reference"]
    candidates = []
    for name, threads in THREAD_VARIANTS:
        label = f"candidate_{name}"
        work = work_root / name
        if ok[label]:
            details["counts"][label] = variant_counts(work, work / "candidate.iact")
        details["threads"][label] = threads
        candidates.append(label)
    details["counts"]["reference"] = {
        "into_simtel_weight": iact_telescope_weight(source)
    }
    metrics["telescope_photons"] = details["counts"]["reference"]["into_simtel_weight"]
    for env_label, threads in ENVELOPE_VARIANTS:
        env_work = work_root / env_label.removeprefix("candidate_")
        if ok[env_label]:
            metrics.update(
                child_metrics(env_label, details["stages"][env_label]["pipeline"])
            )
            details["counts"][env_label] = pipeline_counts(env_work)
            camera_only = _camera_only_lines(env_work)
            metrics[f"{env_label}_camera_only_lines"] = camera_only
            metrics[f"{env_label}_camera_only_missing"] = int(camera_only == 0)
        metrics["candidate_envelope_headroom"] = headroom
        details["threads"][env_label] = threads
        candidates.append(env_label)
    if ok[UNTHINNED_LABEL]:
        metrics.update(
            child_metrics(
                UNTHINNED_LABEL, details["stages"][UNTHINNED_LABEL]["pipeline"]
            )
        )
        details["counts"][UNTHINNED_LABEL] = pipeline_counts(unthinned_work)
    details["threads"][UNTHINNED_LABEL] = 1
    candidates.append(UNTHINNED_LABEL)
    try:
        physics_check(
            root,
            source,
            inputs,
            inputs["qe"],
            headroom,
            metrics,
            details,
            work_root=work_root,
        )
        metrics["physics_completed"] = 1
    except (RuntimeError, OSError, ValueError) as error:  # reported, not hidden
        metrics["physics_completed"] = 0
        details["failures"]["physics"] = str(error)[:400]
    _ratios(metrics, reference_ok, candidates)
    metrics.update(normalised_metrics(metrics, details))
    metrics.update(reference_optics_metrics(metrics, details))
    metrics.update(gain_metrics(metrics))
    budgets, reasons = performance_budget.shower_budgets(metrics)
    metrics.update(budgets)
    details["budget_reasons"] = reasons


def run_source(root, work_root, name, reps, metrics, details):
    spec = SOURCES[name]
    reference_dir = work_root / "reference"
    candidate_dir = work_root / "candidate"
    reference_dir.mkdir(parents=True, exist_ok=True)
    candidate_dir.mkdir(parents=True, exist_ok=True)
    details["settings"].update(
        {
            "events": spec["events"],
            "photons_per_event": spec["photons"],
            "bypass_optics": spec["bypass"],
        }
    )
    ok = run_interleaved(
        [
            ("reference", source_reference(reference_dir, spec)),
            ("candidate_threads1", source_candidate(candidate_dir, spec)),
        ],
        reps,
        metrics,
        details,
    )
    reference_ok = ok["reference"]
    details["threads"]["candidate_threads1"] = (
        "1 (all stages single-threaded; default is the same)"
    )
    _ratios(metrics, reference_ok, ["candidate_threads1"])
    budgets, reasons = performance_budget.source_budgets(metrics, "candidate_threads1")
    metrics.update(budgets)
    details["budget_reasons"] = reasons


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
    parser.add_argument(
        "--events",
        type=int,
        default=None,
        help="air-shower: first N showers (default all)",
    )
    parser.add_argument("--reps", type=int, default=REPS)
    args = parser.parse_args()
    if args.case == "air-shower" and args.config is None:
        parser.error("air-shower needs --config")
    root = Path.cwd().resolve()
    # Work files (pipes, photon streams, sim_telarray output) go to a container-local directory:
    # a bind-mounted results directory would measure the mount, not the chain.
    work_root = Path(
        os.environ.get("PERF_WORK_DIR") or tempfile.mkdtemp(prefix="perf-chain-")
    )
    work_root.mkdir(parents=True, exist_ok=True)
    showers = performance.select_input(
        args.input.resolve(), args.events, work_root / "showers-input.iact"
    )
    metrics = {}
    load_before = performance.loadavg()
    if load_before is not None:
        metrics["host_loadavg_before"] = load_before
    host = performance.host_info()
    if host["nproc"]:
        metrics["host_nproc"] = host["nproc"]
    details = {
        "case": args.case,
        "host": host,
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
        details["settings"]["showers"] = performance.shower_count(showers)
        run_air_shower(
            root, work_root, showers, args.config, args.reps, metrics, details
        )
    else:
        run_source(root, work_root, args.case, args.reps, metrics, details)
    load_after = performance.loadavg()
    if load_after is not None:
        metrics["host_loadavg_after"] = load_after
    metrics["timing_unreliable"], details["timing_reliability"] = timing_reliability(
        metrics, list(details["pipelines"])
    )
    walls = {
        label: {name: _median(records, "wall_s") for name, records in per_stage.items()}
        for label, per_stage in details["stages"].items()
    }
    if walls:
        _plot(root / "wall-by-stage.png", walls)
    common.write_results(metrics, details, root)
    shutil.rmtree(work_root, ignore_errors=True)


if __name__ == "__main__":
    main()
