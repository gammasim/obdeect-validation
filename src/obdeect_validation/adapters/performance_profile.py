"""In-process profile of the obdeect stages with Linux perf: where trace-stream and propagation spend time.

Each selected stage of the LSTN-01 air-shower chain (file mode, one thread, as in
performance_chain.shower_candidate) runs under `perf record -g -F 999`. Its report (perf report,
no children) goes to <stage>-perf-report.txt in the case directory; the top 15 symbols go to
details.json. The prerequisite stages (and the optical model compile) run once without profiling
to provide the stage inputs.

Fallbacks, recorded in details.json (profile.<stage>.mode and reason):
    perf missing or blocked (perf_event_paranoid, no PERFMON capability): perf stat with
    task-clock, cycles, instructions, cache-misses; if that fails too: plain timing.

Only native programs are profiled where they exist (performance_chain.native_argv). trace-stream
has no native program: its obdeect wrapper is profiled, so the report shows the wrapper's
interpreter and extension frames as well.

Metrics: <stage>_profile_completed (1 when the stage ran in any mode), <stage>_profiled (1 when a
perf record profile with samples was reported), <stage>_perf_samples, <stage>_top1_percent.

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.performance_profile --config cfg/model/7.0.0/CTAO.cfg [--events N] [--stage NAME]
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from obdeect_validation.adapters import chain_optics as co
from obdeect_validation.adapters import common, performance
from obdeect_validation.adapters.performance_chain import (
    compile_argv,
    envelope_files,
    headroom_for,
    native_argv,
    prepare_inputs,
)

STAGES = (
    "from_corsika",
    "propagation",
    "propagation_envelope",
    "trace_stream",
    "to_simtel",
)
# Stages that must run (unprofiled) before a stage can read its input file.
NEEDS = {
    "from_corsika": (),
    "propagation": ("from_corsika",),
    "propagation_envelope": ("from_corsika",),
    "trace_stream": ("propagation", "compile_model"),
    "to_simtel": ("trace_stream",),
    "compile_model": (),
}
PERF_RECORD = ["record", "-g", "-F", "999"]
PERF_STAT_EVENTS = "task-clock,cycles,instructions,cache-misses"
TOP_SYMBOLS = 15
# perf report line: overhead, command, shared object, [.]/[k] symbol type, symbol.
REPORT_LINE = re.compile(r"^\s*(\d+(?:\.\d+)?)%\s+(\S+)\s+(\S+)\s+\[(.)\]\s+(.*\S)\s*$")
# perf messages when the kernel refuses sampling (record) or counting (stat).
BLOCKED = re.compile(
    r"perf_event_paranoid|Permission denied|Operation not permitted|"
    r"Access to performance monitoring",
    re.IGNORECASE,
)
# Stage paths of the work directory, named as in performance_chain.shower_candidate.
FILES = {
    "native": "showers.gst",
    "propagated": "propagated.gst",
    "propagated_envelope": "propagated-envelope.gst",
    "model": "lst.model",
    "traced": "traced.gst",
    "camera_stream": "candidate.iact",
    "camera_map": "camera-map.txt",
}


def stage_argv(name, work, source, inputs, envelope):
    """(argv, stdout_file) of one stage; stdout_file is set only for the stage that writes stdout.

    Arguments are those of performance_chain.shower_candidate in file mode (threads 1); the
    envelope stage is the thinned propagation of performance_chain.envelope_pipeline.
    """
    paths = {key: work / value for key, value in FILES.items()}
    if name == "from_corsika":
        return ["obdeect", "stream", "from-corsika", source], paths["native"]
    if name == "propagation":
        return (
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
                "-i",
                paths["native"],
                "-o",
                paths["propagated"],
            ],
            None,
        )
    if name == "propagation_envelope":
        return (
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
                # ENVELOPE_FACTOR of performance_chain ("qe")
                "--envelope",
                f"qe={envelope['envelope']}",
                "-i",
                paths["native"],
                "-o",
                paths["propagated_envelope"],
            ],
            None,
        )
    if name == "compile_model":
        return compile_argv(inputs, paths["model"]), None
    if name == "trace_stream":
        return (
            [
                "obdeect",
                "ray-tracing",
                "trace-stream",
                "--model",
                paths["model"],
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
                paths["propagated"],
                "-o",
                paths["traced"],
            ],
            None,
        )
    if name == "to_simtel":
        return (
            [
                "obdeect",
                "stream",
                "to-simtel",
                "--camera-map",
                paths["camera_map"],
                "-o",
                paths["camera_stream"],
                paths["traced"],
            ],
            None,
        )
    raise ValueError(f"unknown stage {name}")


def write_camera_map(path, rotation):
    """Camera map of the candidate chain (format of chain_optics.camera_map, mirror class 1)."""
    path.write_text(
        "# telescope_id mirror_class camera_rotation_deg\n"
        + "".join(f"{tel} 1 {rotation}\n" for tel in co.TELESCOPES)
    )


def _execute(argv, log, stdout_file=None, prefix=()):
    """Run prefix + argv (native form of argv): stdout to stdout_file, stderr to log. Returns the exit code."""
    command = [*map(str, prefix), *map(str, native_argv(argv))]
    with Path(log).open("wb") as log_stream:
        if stdout_file:
            with Path(stdout_file).open("wb") as output:
                return subprocess.run(
                    command, stdout=output, stderr=log_stream, check=False
                ).returncode
        return subprocess.run(
            command, stdout=log_stream, stderr=log_stream, check=False
        ).returncode


def parse_report(text):
    """Symbols of a perf report (no children) in report order: overhead, command, shared object, symbol."""
    symbols = []
    for line in text.splitlines():
        match = REPORT_LINE.match(line)
        if match:
            symbols.append(
                {
                    "overhead_percent": float(match[1]),
                    "command": match[2],
                    "shared_object": match[3],
                    "symbol": match[5],
                }
            )
    return symbols


def count_samples(script_text):
    """Samples of perf script output: one unindented header line per sample, frames are indented."""
    return sum(
        1 for line in script_text.splitlines() if line.strip() and not line[0].isspace()
    )


def flamegraph(script_text, svg):
    """SVG flame graph of perf script output via stackcollapse-perf.pl and flamegraph.pl (FlameGraph)."""
    collapse = shutil.which("stackcollapse-perf.pl")
    render = shutil.which("flamegraph.pl")
    if not (collapse and render):
        return "stackcollapse-perf.pl or flamegraph.pl not on PATH; no flame graph"
    try:
        folded = subprocess.run(
            [collapse], input=script_text, capture_output=True, text=True, check=True
        ).stdout
        svg.write_text(
            subprocess.run(
                [render], input=folded, capture_output=True, text=True, check=True
            ).stdout
        )
    except subprocess.CalledProcessError as error:
        return f"{error.cmd[0]} exited {error.returncode}; no flame graph"
    return svg.name


def _profile_with_perf(perf, name, argv, stdout_file, work, case_dir, entry):
    """perf record, then perf report and perf script. Returns the symbols and sample count, or None."""
    data = work / f"{name}.perf.data"
    code = _execute(
        argv,
        work / f"{name}-record.log",
        stdout_file,
        prefix=[perf, *PERF_RECORD, "-o", data, "--"],
    )
    if code or not data.is_file():
        entry["reason"] = _reason("perf record", code, work / f"{name}-record.log")
        return None
    report = subprocess.run(
        [
            perf,
            "report",
            "--stdio",
            "--no-children",
            "--percent-limit",
            "1",
            "-i",
            data,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if report.returncode:
        entry["reason"] = f"perf report exited {report.returncode}"
        return None
    (case_dir / f"{name}-perf-report.txt").write_text(report.stdout)
    entry["report"] = f"{name}-perf-report.txt"
    script = subprocess.run(
        [perf, "script", "-i", data], capture_output=True, text=True, check=False
    )
    samples = count_samples(script.stdout) if script.returncode == 0 else 0
    if script.returncode == 0:
        entry["flamegraph"] = flamegraph(script.stdout, case_dir / f"{name}-flame.svg")
    return parse_report(report.stdout)[:TOP_SYMBOLS], samples


def _reason(tool, code, log):
    """Why a perf call failed: blocked by the kernel, or its exit code and the log tail."""
    text = Path(log).read_text(errors="replace")
    tail = " | ".join(text.strip().splitlines()[-2:])
    if BLOCKED.search(text):
        return f"{tool} blocked: {tail}"
    return f"{tool} exited {code}: {tail}"


def profile_stage(name, argv, stdout_file, perf, work, case_dir, metrics, details):
    """Profile one stage (perf record, else perf stat, else timing) and record its metrics and details."""
    entry = {
        "argv": [str(arg) for arg in native_argv(argv)],
        "mode": None,
        "reason": None,
    }
    details["profile"][name] = entry
    result = None
    if perf is None:
        entry["reason"] = "perf not found on PATH"
    else:
        result = _profile_with_perf(
            perf, name, argv, stdout_file, work, case_dir, entry
        )
    if result is not None:
        entry["mode"] = "perf-record"
        symbols, samples = result
        entry["top"] = symbols
        entry["samples"] = samples
        metrics[f"{name}_profiled"] = int(bool(symbols))
        metrics[f"{name}_perf_samples"] = samples
        if symbols:
            metrics[f"{name}_top1_percent"] = symbols[0]["overhead_percent"]
    else:
        metrics[f"{name}_profiled"] = 0
        metrics[f"{name}_perf_samples"] = 0
        if perf is not None:
            stat = case_dir / f"{name}-perf-stat.txt"
            code = _execute(
                argv,
                work / f"{name}-stat.log",
                stdout_file,
                prefix=[perf, "stat", "-e", PERF_STAT_EVENTS, "-o", stat, "--"],
            )
            if code == 0:
                entry["mode"] = "perf-stat"
                entry["report"] = stat.name
        if entry["mode"] is None:
            entry["mode"] = "timing"
            record = performance.run_timed(
                native_argv(argv), work / f"{name}-timing.log", stdout_file=stdout_file
            )
            entry["wall_s"] = record["wall_s"]
    metrics[f"{name}_profile_completed"] = 1
    return entry


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", type=Path, required=True, help="sim_telarray config from simtools"
    )
    parser.add_argument("--input", type=Path, default=performance.INPUT)
    parser.add_argument("--events", type=int, default=5, help="first N showers")
    parser.add_argument(
        "--stage", choices=[*STAGES, "all"], default="all", help="stage to profile"
    )
    args = parser.parse_args()
    case_dir = Path.cwd().resolve()
    # Work files (perf data, streams) go to a fresh directory under PERF_WORK_DIR or the temp
    # directory; only the reports stay in the case directory.
    base = Path(os.environ.get("PERF_WORK_DIR") or tempfile.gettempdir())
    base.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="perf-profile-", dir=base))
    perf = shutil.which("perf")
    metrics = {}
    details = {
        "case": "performance.profile",
        "host": performance.host_info(),
        "settings": {
            "events": args.events,
            "perf": perf or "not found",
            "perf_record": " ".join(PERF_RECORD),
            "perf_stat_events": PERF_STAT_EVENTS,
            "threads": 1,
        },
        "profile": {},
    }
    try:
        showers = performance.select_input(
            args.input.resolve(), args.events, work / "showers-input.iact"
        )
        details["settings"]["showers"] = performance.shower_count(showers)
        inputs, _ = prepare_inputs(work, args.config)
        details["settings"]["config"] = str(inputs["config"])
        envelope = envelope_files(
            work, inputs["qe"], headroom_for(inputs["qe_variation"])
        )
        write_camera_map(work / FILES["camera_map"], inputs["rotation"])
        selected = list(STAGES) if args.stage == "all" else [args.stage]
        done = set()

        def ensure(stage):
            """Run stage (unprofiled) and its prerequisites, unless they already produced their output."""
            if stage in done:
                return
            for dep in NEEDS[stage]:
                ensure(dep)
            argv, stdout_file = stage_argv(stage, work, showers, inputs, envelope)
            common.command(
                native_argv(argv),
                work / f"{stage}-prepare.log",
                stdout_file=stdout_file,
            )
            done.add(stage)

        for stage in selected:
            for dep in NEEDS[stage]:
                ensure(dep)
            argv, stdout_file = stage_argv(stage, work, showers, inputs, envelope)
            try:
                profile_stage(
                    stage, argv, stdout_file, perf, work, case_dir, metrics, details
                )
                done.add(stage)
            except RuntimeError as error:
                details["profile"][stage]["failure"] = str(error)
                metrics[f"{stage}_profile_completed"] = 0
    finally:
        common.write_results(metrics, details, case_dir)
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
