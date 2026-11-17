"""Multi-telescope air-shower timing and camera check on the North array (4 LSTs, 2 MSTs).

Input: the array shower file (CORSIKA telescopes 1-6 = LSTN-01..04, MSTN-03, MSTN-15, at their
simtools positions; 20 deg, 100-3000 GeV, CSCAT 25 km... see shower_matrix.ARRAY_NORTH_TELESCOPES).
The config is the simtools array config of the same six elements (CTAO.cfg with CTAO-<element>.cfg).

Reference (timed): one sim_telarray run with all six telescopes (no -DTELESCOPE, so sim_telarray
reads every telescope block of the config), Bypass_Optics=0 (full optics), production atmosphere.
Candidate (timed, gate) tee: the telescopes are grouped by QE envelope (LST 1-4, MST 5-6), because a
propagation takes one envelope. Per group one pipe: from-corsika | propagation --telescopes <group>
(QE envelope, thinned) | obdeect stream tee into one FIFO per telescope. Each telescope has a branch:
trace-stream (model of element k) | to-simtel | sim_telarray -DTELESCOPE=k (Bypass_Optics=3, camera
only). The groups run one after another, the branches of a group at once. The chain time is the sum
of the group pipes; the model compiles run before the timing and are left out (as in performance_budget).
Candidate (timed, report only) per_telescope: for each telescope k, one pipe as in
performance_chain.envelope_pipeline (from-corsika | propagation | trace-stream | to-simtel |
sim_telarray -DTELESCOPE=k). Each pipe reads CORSIKA and propagates again, so its CPU is overstated
against the tee. Its metrics carry the prefix candidate_per_telescope_.

Deviations from a single joint run, recorded in details.json:
  * No joint candidate run: the per-telescope camera streams are not merged into one sim_telarray
    input (no event merge tool), so the array trigger is not compared. Both sides use
    trigger_telescopes=1, so every telescope's photoelectrons are stored for the same events.
  * Two propagations, one per envelope group (details settings tee_groups): each reads the CORSIKA
    file once. The tee's standard output is discarded: the propagated stream is not stored.
  * Model compile uses --zenith-deg 20 for all six (the MSTs fail without it: axes_offsets).
  * min_photons = 0 on both sides (performance_chain.MIN_PHOTONS). sim_telarray skips a telescope
    event below min_photons (120 LST, 450 MST in the config) using the photons of its bunches; a
    thinned stream has one photon per surviving bunch, so the threshold would drop events on the
    candidate side only. min_photoelectrons (30 MST, 20 LST) is kept and recorded in details.

Limits (fixed here):
  * Per telescope, on the tee outputs: photoelectron bias within max(1 %, 5 sigma)
    (chain_showers.compare_camera, PE_BIAS_FLOOR and SIGMA of chain_showers). Counted in
    equal_statistics_budget_exceeded. The per-telescope candidate gets the same checks, reported only.
  * Tee chain wall and CPU (sum over the group pipes, compiles excluded) <= reference wall and CPU.
  * Every tee telescope's sim_telarray prints the camera-only line of Bypass_Optics=3.
  * Per telescope, sim_telarray's "Used photons" equals "All photons" (no photon is dropped before
    the camera), and the to-simtel bunch count (sim_telarray "All photons", the sum of the bunch
    photon counts it reads) equals the trace-stream output rows of that telescope. Metrics
    photon_used_ne_all_tel<k>, to_simtel_bunch_ne_trace_rows_tel<k>; 0 = pass.
  * A telescope with fewer than MIN_MATCHED_SAMPLES (20) matched photoelectron samples is not
    measured: its bias flag is 1 and the reason goes to details.

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.performance_array air-shower-array --config cfg/model/7.0.0/CTAO.cfg [--events N] [--reps N]
"""

from __future__ import annotations

import argparse
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import threading
import time
import traceback
from pathlib import Path

from obdeect_validation.adapters import (
    chain_optics as co,
)
from obdeect_validation.adapters import (
    chain_showers as cs,
)
from obdeect_validation.adapters import (
    common,
    performance,
)
from obdeect_validation.adapters import (
    performance_chain as pc,
)

INPUT = common.VALIDATION / "inputs/gamma-20deg-100-3000GeV-array-north.iact"
# (CORSIKA telescope number = sim_telarray TELESCOPE index, element name), as in the simtools config.
ARRAY = (
    (1, "LSTN-01"),
    (2, "LSTN-02"),
    (3, "LSTN-03"),
    (4, "LSTN-04"),
    (5, "MSTN-03"),
    (6, "MSTN-15"),
)
CONFIG_NAME = "CTAO.cfg"
ATMOSPHERE_FILE = pc.ATMOSPHERE_FILE
MIRROR_CLASS = 1  # as the single-LST camera map of performance_chain.camera_map
ZENITH_DEG, AZIMUTH_DEG = (
    co.ZENITH_DEG,
    co.AZIMUTH_DEG,
)  # pointing; CORSIKA PHIP 0 is 180 deg here
VERSION = co.VERSION
#: Matched photoelectron samples (event, telescope) below which a telescope's bias is not measured.
MIN_MATCHED_SAMPLES = 20
#: Metric prefix of the per-telescope candidate (report only).
REPORT_PREFIX = "candidate_per_telescope_"
TELESCOPE_OPTIONS = {
    "altitude": 2156,
    "telescope_theta": ZENITH_DEG,
    "telescope_phi": AZIMUTH_DEG,
}


def element_inputs(config_dir, element):
    """Camera rotation, QE table and QE variation of one element (from its CTAO-<element>.cfg)."""
    cfg_text = (config_dir / f"CTAO-{element}.cfg").read_text()
    return {
        "rotation": co.camera_rotation_deg(config_dir, cfg_text),
        "qe": config_dir / co.cfg_value(cfg_text, "quantum_efficiency"),
        "qe_variation": float(co.cfg_value(cfg_text, "qe_variation")),
        "min_photoelectrons": int(float(co.cfg_value(cfg_text, "min_photoelectrons"))),
    }


def camera_map(root, telescope, rotation):
    """to-simtel camera map with one telescope (id = CORSIKA telescope number)."""
    path = root / "camera-map.txt"
    path.write_text(
        "# telescope_id mirror_class camera_rotation_deg\n"
        f"{telescope} {MIRROR_CLASS} {rotation}\n"
    )
    return path


def compile_argv(models, element, model_path):
    """obdeect optical model of one element, with the pointing zenith (needed for the MSTs)."""
    return [
        "obdeect",
        "ray-tracing",
        "compile",
        "--simulation-models",
        models,
        "--model",
        element,
        "--version",
        VERSION,
        "--seed",
        co.SEED,
        "--zenith-deg",
        ZENITH_DEG,
        "--output",
        model_path,
    ]


def reference_options(atmosphere):
    """sim_telarray options of the reference: full optics, all six telescopes, min_photons as the candidate."""
    return {
        **TELESCOPE_OPTIONS,
        "maximum_telescopes": len(ARRAY),
        "Bypass_Optics": 0,
        "trigger_telescopes": 1,
        "min_photons": pc.MIN_PHOTONS,
        "atmospheric_transmission": atmosphere,
    }


def candidate_options(atmosphere_name):
    """sim_telarray options of the candidate camera stage (Bypass_Optics=3)."""
    return {
        **TELESCOPE_OPTIONS,
        "maximum_telescopes": len(ARRAY),
        "Bypass_Optics": 3,
        "trigger_telescopes": 1,
        "min_photons": pc.MIN_PHOTONS,
        "atmospheric_transmission": atmosphere_name,
    }


def reference_run(work, source, config, atmosphere):
    """One sim_telarray run with all six telescopes of the config, Bypass_Optics=0, timed."""
    output = common.sim_telarray(
        work,
        source,
        "reference",
        reference_options(atmosphere),
        config=str(config),
        defines=(),
        extra=(
            f"-I{config.parent}",
        ),  # CTAO.cfg includes InvalidTelescope.cfg and the telescopes
        runner=pc.run_timed,
    )
    common.check_atmosphere_log(work / "reference.log", ATMOSPHERE_FILE)
    return output


def candidate_stages(work, source, inputs, telescope, element, files):
    """Pipe of one telescope (performance_chain.envelope_pipeline with this telescope and model)."""
    model = work / "model"
    camera = camera_map(work, telescope, inputs["rotation"])
    table_dir = files["table_dir"]
    sim_options = candidate_options(inputs["atmosphere_name"])
    sim_argv, output = common.sim_telarray_argv(
        work,
        "-",
        "candidate",
        sim_options,
        config=str(inputs["config"]),
        defines=(f"-DTELESCOPE={telescope}",),
        extra=(f"-I{table_dir}", f"-I{inputs['config'].parent}"),
        binary=pc.CANDIDATE_SIMTEL,
    )
    stages = pc.envelope_pipeline(
        source,
        inputs["atmosphere"],
        files["envelope"],
        model,
        camera,
        sim_argv,
        0,
        telescope=telescope,
    )
    logs = {name: work / f"{name}.log" for name, _ in stages}
    return model, stages, logs, output


def telescope_steps(work_root, source, inputs, files, elements):
    """Timed steps of the per-telescope candidate: compile_tel<k> and pipe_tel<k>."""
    steps, outputs = [], {}
    for telescope, element in ARRAY:
        work = work_root / f"candidate-tel{telescope}"
        work.mkdir(parents=True, exist_ok=True)
        model, stages, logs, output = candidate_stages(
            work,
            source,
            {**inputs, **elements[element]},
            telescope,
            element,
            files[element],
        )
        outputs[telescope] = (work, output)
        steps.append(
            (
                f"compile_tel{telescope}",
                lambda m=model, e=element, w=work: pc.run_timed(
                    compile_argv(inputs["models"], e, m), w / "compile.log"
                ),
            )
        )
        steps.append(
            (
                f"pipe_tel{telescope}",
                lambda s=stages, o=work / "candidate.stdout", lg=logs: (
                    pc._checked_pipeline(s, o, lg)
                ),
            )
        )
    return steps, outputs


def tee_groups(elements):
    """Telescopes that share one QE envelope, in ARRAY order: one propagation and one tee each."""
    groups = {}
    for telescope, element in ARRAY:
        key = (str(elements[element]["qe"]), elements[element]["headroom"])
        groups.setdefault(key, []).append(telescope)
    return list(groups.values())


def group_name(telescopes):
    """Name of one tee group (step and settings), e.g. tee_tel1_4."""
    return f"tee_tel{telescopes[0]}_{telescopes[-1]}"


def with_telescopes(argv, telescopes):
    """Propagation argv with the telescope selection (--telescopes) before its input option."""
    index = argv.index("-i")
    return [
        *argv[:index],
        "--telescopes",
        ",".join(map(str, telescopes)),
        *argv[index:],
    ]


def with_input(argv, path):
    """argv with the input option ('-i') set to path."""
    argv = list(argv)
    argv[argv.index("-i") + 1] = str(path)
    return argv


def tee_branches(work_root, source, inputs, files, elements):
    """Per telescope of the tee candidate: work directory, model, stages and logs of its branch."""
    branches = {}
    for telescope, element in ARRAY:
        work = work_root / f"tee-tel{telescope}"
        work.mkdir(parents=True, exist_ok=True)
        model, stages, logs, output = candidate_stages(
            work,
            source,
            {**inputs, **elements[element]},
            telescope,
            element,
            files[element],
        )
        branches[telescope] = {
            "element": element,
            "work": work,
            "model": model,
            "stages": stages,
            "logs": logs,
            "output": output,
        }
    return branches


def compile_tee_models(branches, models):
    """Optical model of every telescope of the tee candidate (timed apart: not part of the chain)."""
    for telescope, _ in ARRAY:
        branch = branches[telescope]
        pc.run_timed(
            compile_argv(models, branch["element"], branch["model"]),
            branch["work"] / "compile.log",
        )


def _start_pipe(stages, output, logs, running):
    """Start (name, argv) stages as one pipe; the last one writes output (None: discarded)."""
    previous = subprocess.DEVNULL
    for index, (name, argv) in enumerate(stages):
        last = index == len(stages) - 1
        argv = pc.native_argv([str(arg) for arg in argv])
        with Path(logs[name]).open("wb") as log:
            if last and output is not None:
                with Path(output).open("wb") as out:
                    proc = subprocess.Popen(
                        argv, stdin=previous, stdout=out, stderr=log
                    )
            else:
                stdout = subprocess.DEVNULL if last else subprocess.PIPE
                proc = subprocess.Popen(argv, stdin=previous, stdout=stdout, stderr=log)
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


def run_fan_out(pipe, branches, logs):
    """Start pipe and every branch at once; time each process as performance_chain.run_pipeline.

    pipe: [(name, argv)] ending in the tee, whose standard output is discarded. branches:
    [(stages, output)], one pipe each; its first stage reads its FIFO through its argv. One pipeline
    record goes to performance.RUSAGE (wall to the last exit, CPU summed, RSS maximal, children).
    A failed stage raises RuntimeError after the record is kept.
    """
    start = time.perf_counter()
    running = []
    try:
        _start_pipe(pipe, None, logs, running)
        for stages, output in branches:
            _start_pipe(stages, output, logs, running)
    except BaseException:
        pc._abort_pipeline(
            [proc for _, _, proc, _ in running],
            [sampler for _, _, _, sampler in running if sampler],
        )
        raise
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


def tee_pipe(telescopes, branches):
    """One group: from-corsika | propagation --telescopes | tee, and a FIFO branch per telescope.

    The branch of telescope k is trace-stream (reads its FIFO) | to-simtel | sim_telarray -DTELESCOPE=k.
    Each branch's sim_telarray must read the atmosphere table (as performance_chain._checked_pipeline).
    """
    first = branches[telescopes[0]]
    source_stage, propagation_stage = first["stages"][:2]
    fifo_dir = Path(
        tempfile.mkdtemp(prefix="fifo-")
    )  # local: bind mounts may refuse FIFOs
    try:
        fifos = {t: fifo_dir / f"tel{t}.fifo" for t in telescopes}
        for fifo in fifos.values():
            os.mkfifo(fifo)
        pipe = [
            source_stage,
            (
                propagation_stage[0],
                with_telescopes(propagation_stage[1], telescopes),
            ),
            ("tee", ["obdeect", "stream", "tee", *fifos.values()]),
        ]
        logs = {
            "from_corsika": first["logs"]["from_corsika"],
            "propagation": first["logs"]["propagation"],
            "tee": first["work"] / "tee.log",
        }
        pipes = []
        for telescope in telescopes:
            branch = branches[telescope]
            stages = []
            for name, argv in branch["stages"][2:]:
                member = f"tel{telescope}_{name}"
                logs[member] = branch["logs"][name]
                if name == "trace_stream":
                    argv = with_input(argv, fifos[telescope])
                stages.append((member, argv))
            pipes.append((stages, branch["work"] / "candidate.stdout"))
        run_fan_out(pipe, pipes, logs)
        for telescope in telescopes:
            branch = branches[telescope]
            text = (branch["work"] / "candidate.stdout").read_text(
                errors="replace"
            ) + branch["logs"]["sim_telarray"].read_text(errors="replace")
            common.check_atmosphere(text, ATMOSPHERE_FILE)
    finally:
        shutil.rmtree(fifo_dir, ignore_errors=True)


def tee_candidate(work_root, source, inputs, files, elements, reps, metrics, details):
    """Tee candidate (gate): its branches, or None when it could not be measured (recorded)."""
    branches = tee_branches(work_root, source, inputs, files, elements)
    try:
        compile_tee_models(branches, inputs["models"])
    except (RuntimeError, OSError) as error:
        details["failures"]["candidate_tee_compile"] = str(error)[:400]
        metrics["candidate_tee_completed"] = 0
        return None
    steps = [
        (group_name(group), lambda g=group: tee_pipe(g, branches))
        for group in tee_groups(elements)
    ]
    ok = pc._run("candidate_tee", steps, reps, metrics, details)
    return branches if ok else None


def tee_member_cpu(stages):
    """Median over repetitions of each member's CPU, summed over the group pipes of a repetition."""
    groups = [name for name in stages if name.startswith("tee_tel")]
    per_rep = []
    for rep in range(len(stages[groups[0]])):
        cpu = {}
        for name in groups:
            for child in stages[name][rep].get("children", []):
                cpu[child["name"]] = cpu.get(child["name"], 0.0) + child["cpu_s"]
        per_rep.append(cpu)
    names = sorted({name for cpu in per_rep for name in cpu})
    return {
        name: statistics.median(cpu.get(name, 0.0) for cpu in per_rep) for name in names
    }


def tee_shares(total, members):
    """Metrics of each member's CPU and its share of the tee chain CPU."""
    out = {}
    for name, cpu in members.items():
        out[f"candidate_tee_{name}_cpu_s"] = cpu
        if total > 0:
            out[f"candidate_tee_{name}_cpu_share"] = cpu / total
    return out


def photon_statistics(stdout_text, telescope):
    """(All photons, Used photons) of one telescope from sim_telarray's statistics block.

    sim_telarray prints one column per telescope index (telescope - 1).
    """
    values = {}
    for label in ("All photons", "Used photons"):
        lines = [
            line
            for line in stdout_text.splitlines()
            if line.strip().startswith(f"{label}:")
        ]
        if len(lines) != 1:
            raise ValueError(
                f"expected one '{label}' line in sim_telarray output, found {len(lines)}"
            )
        columns = [float(x) for x in lines[0].split(":", 1)[1].split()]
        if telescope - 1 >= len(columns):
            raise ValueError(f"'{label}' line has no column for telescope {telescope}")
        values[label] = int(columns[telescope - 1])
    return values["All photons"], values["Used photons"]


def trace_out_rows(trace_text):
    """Output rows of the trace-stream summary (telescope of the run; see performance_chain)."""
    match = pc.TRACE_SUMMARY.search(trace_text)
    if match is None:
        raise ValueError("no trace-stream summary line")
    return int(match[2])


def photon_checks(outputs, metrics, details, prefix=""):
    """Per telescope: no photon dropped before the camera, and to-simtel bunches = trace rows.

    to-simtel bunch count is taken from sim_telarray's "All photons" (the sum of the bunch photon
    counts it reads; the to-simtel stream is piped and not stored). Each check is 0 = pass.
    prefix names the metrics and details of one candidate ('' for the gated tee candidate).
    """
    details[f"{prefix}photon_checks"] = {}
    failed = 0
    for telescope, _ in ARRAY:
        work, _ = outputs[telescope]
        key = f"tel{telescope}"
        all_photons, used = photon_statistics(
            (work / "candidate.stdout").read_text(errors="replace"), telescope
        )
        rows = trace_out_rows((work / "trace_stream.log").read_text(errors="replace"))
        flags = {
            f"{prefix}photon_used_ne_all_{key}": int(used != all_photons),
            f"{prefix}to_simtel_bunch_ne_trace_rows_{key}": int(all_photons != rows),
        }
        metrics.update(
            {
                f"{prefix}photon_all_{key}": all_photons,
                f"{prefix}photon_used_{key}": used,
                f"{prefix}to_simtel_bunches_{key}": all_photons,
                f"{prefix}trace_out_rows_{key}": rows,
                **flags,
            }
        )
        failed += sum(flags.values())
        details[f"{prefix}photon_checks"][key] = {
            "all": all_photons,
            "used": used,
            "trace_rows": rows,
        }
    metrics[f"{prefix}photon_checks_failed"] = failed


def measured(comparison, key, details, prefix=""):
    """A comparison with fewer than MIN_MATCHED_SAMPLES matched samples is not measured (flag 1)."""
    matched = comparison.get("matched_event_telescope_samples", 0)
    if matched >= MIN_MATCHED_SAMPLES:
        return comparison
    details["failures"][f"{prefix}physics_{key}"] = (
        f"not measured: {matched} matched samples < {MIN_MATCHED_SAMPLES}"
    )
    return {
        "matched_event_telescope_samples": matched,
        "photoelectron_bias_limit_exceeded": 1,
    }


def chain_medians(stages, key):
    """Median over repetitions of the sum of the telescope pipes (compile stages excluded)."""
    pipes = [name for name in stages if name.startswith("pipe_tel")]
    reps = len(stages[pipes[0]])
    return statistics.median(
        sum(stages[name][rep][key] for name in pipes) for rep in range(reps)
    )


def physics_check(root, reference_output, outputs, metrics, details, prefix=""):
    """Per telescope: reference camera events against the candidate events of that telescope.

    A telescope with no candidate events (or no matched samples) is not compared: its bias flag
    is set to 1 (exceeded) and the reason goes to details, so it cannot pass silently.
    prefix names the metrics and details of one candidate ('' for the gated tee candidate).
    """
    reference = cs.camera_events(reference_output)
    details[f"{prefix}physics"] = {}
    completed = 1
    for telescope, element in ARRAY:
        work, output = outputs[telescope]
        key = f"tel{telescope}"
        out = root / "physics" / f"{prefix}{key}-{element}"
        out.mkdir(parents=True, exist_ok=True)
        index = telescope - 1  # sim_telarray telescope index is 0-based
        try:
            candidate = cs.camera_events(output)
            comparison = measured(
                cs.compare_camera(reference, candidate, out, {index: element}),
                key,
                details,
                prefix,
            )
        except ValueError as error:  # no matched samples: exceeded, reported
            completed = 0
            metrics[f"{prefix}physics_{key}_photoelectron_bias_limit_exceeded"] = 1
            details["failures"][f"{prefix}physics_{key}"] = str(error)[:200]
            comparison = {}
        metrics.update({f"{prefix}physics_{key}_{k}": v for k, v in comparison.items()})
        metrics[f"{prefix}{key}_camera_only_lines"] = pc._camera_only_lines(work)
        metrics[f"{prefix}{key}_camera_only_missing"] = int(
            metrics[f"{prefix}{key}_camera_only_lines"] == 0
        )
        details[f"{prefix}physics"][key] = comparison
    metrics[f"{prefix}physics_completed"] = completed


def checks(root, reference_path, outputs, metrics, details, prefix):
    """Photon checks and camera comparison of one candidate; a failure is reported, not hidden."""
    try:
        photon_checks(outputs, metrics, details, prefix)
    except (RuntimeError, OSError, ValueError) as error:  # reported, not hidden
        details["failures"][f"{prefix}photon_checks"] = str(error)[:400]
    if reference_path is None:
        return
    try:
        physics_check(root, reference_path, outputs, metrics, details, prefix)
    except (RuntimeError, OSError, ValueError) as error:  # reported, not hidden
        metrics[f"{prefix}physics_completed"] = 0
        details["failures"][f"{prefix}physics"] = str(error)[:400]


def array_budgets(metrics, stages_wall_cpu):
    """Gates: chain wall and CPU against the reference, and the per-telescope photoelectron stats."""
    out, reasons = {}, {}
    for key, short in (("wall_s", "wall"), ("cpu_s", "cpu")):
        ref = f"reference_{key}"
        chain = stages_wall_cpu.get(key)
        if not isinstance(metrics.get(ref), (int, float)) or chain is None:
            reasons[short] = f"not measured: {ref}, chain {key}"
            out[f"array_chain_{short}_budget_exceeded"] = 1
            continue
        out[f"array_chain_{short}_ratio"] = chain / metrics[ref]
        out[f"array_chain_{short}_budget_exceeded"] = int(chain > metrics[ref])
    limits = [
        metrics.get(f"physics_tel{t}_photoelectron_bias_limit_exceeded")
        for t, _ in ARRAY
    ]
    if all(isinstance(v, int) for v in limits):
        out["equal_statistics_budget_exceeded"] = int(any(v != 0 for v in limits))
    else:
        reasons["equal_statistics"] = "not measured: physics per telescope"
        out["equal_statistics_budget_exceeded"] = 1
    failed = metrics.get("photon_checks_failed")
    if isinstance(failed, int):
        out["photon_gate_budget_exceeded"] = int(failed != 0)
    else:
        reasons["photon_gate"] = "not measured: photon checks"
        out["photon_gate_budget_exceeded"] = 1
    return out, reasons


def run_array(root, work_root, source, config, reps, metrics, details):
    config = config.resolve()
    config_dir = config.parent
    source = source.resolve()
    atmosphere_in = config_dir / ATMOSPHERE_FILE
    atmosphere = work_root / "atmosphere-below-1000nm.dat"
    performance.restrict_table(atmosphere_in, atmosphere)
    inputs = {
        "config": config,
        "models": str(common.models_root()),
        "atmosphere": atmosphere,
        "atmosphere_name": str(atmosphere_in),
    }
    elements = {element: element_inputs(config_dir, element) for _, element in ARRAY}
    files = {}
    for telescope, element in ARRAY:
        qe = elements[element]["qe"]
        headroom = pc.headroom_for(elements[element]["qe_variation"])
        files[element] = pc.envelope_files(
            work_root / f"envelope-{element}", qe, headroom
        )
        elements[element]["headroom"] = headroom
    details["settings"].update(
        {
            "input": str(source),
            "config": str(config),
            "telescopes": {str(t): e for t, e in ARRAY},
            "pointing": f"zenith {ZENITH_DEG} deg, azimuth {AZIMUTH_DEG} deg",
            "headroom": {e: elements[e]["headroom"] for _, e in ARRAY},
            "min_photons": {"reference": pc.MIN_PHOTONS, "candidate": pc.MIN_PHOTONS},
            "min_matched_samples": MIN_MATCHED_SAMPLES,
            "tee_groups": {group_name(g): g for g in tee_groups(elements)},
            "candidate_sim_telarray": pc.CANDIDATE_SIMTEL
            or "production (SIMTEL_CANDIDATE_BINARY unset)",
        }
    )
    reference_dir = work_root / "reference"
    reference_dir.mkdir(parents=True, exist_ok=True)
    reference_output = {}

    def reference_step():
        reference_output["path"] = reference_run(
            reference_dir, source, config, str(atmosphere_in)
        )

    reference_ok = pc._run(
        "reference", [("sim_telarray", reference_step)], reps, metrics, details
    )
    reference_path = reference_output["path"] if reference_ok else None

    # Report only: one pipe per telescope (each re-reads CORSIKA and propagates again).
    steps, outputs = telescope_steps(work_root, source, inputs, files, elements)
    per_ok = pc._run("candidate_per_telescope", steps, reps, metrics, details)
    if per_ok:
        checks(root, reference_path, outputs, metrics, details, REPORT_PREFIX)

    # Gate: tee candidate (one pipe per envelope group, FIFO branch per telescope).
    branches = tee_candidate(
        work_root, source, inputs, files, elements, reps, metrics, details
    )
    if branches is not None:
        metrics.update(
            tee_shares(
                metrics["candidate_tee_cpu_s"],
                tee_member_cpu(details["stages"]["candidate_tee"]),
            )
        )
        if per_ok:
            ratio = pc.ratio(
                metrics["candidate_tee_cpu_s"],
                chain_medians(details["stages"]["candidate_per_telescope"], "cpu_s"),
            )
            if ratio is not None:
                metrics["tee_vs_per_telescope_cpu_ratio"] = ratio
        tee_outputs = {t: (b["work"], b["output"]) for t, b in branches.items()}
        checks(root, reference_path, tee_outputs, metrics, details, "")
    chain = {
        "wall_s": metrics.get("candidate_tee_wall_s"),
        "cpu_s": metrics.get("candidate_tee_cpu_s"),
    }
    budgets, reasons = array_budgets(metrics, chain)
    metrics.update(budgets)
    details["budget_reasons"] = reasons
    details["telescope_elements"] = {
        e: {
            "rotation": v["rotation"],
            "qe_variation": v["qe_variation"],
            "headroom": v["headroom"],
            "min_photoelectrons": v["min_photoelectrons"],
        }
        for e, v in elements.items()
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", choices=["air-shower-array"])
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument(
        "--events", type=int, default=None, help="first N showers (default all)"
    )
    parser.add_argument("--reps", type=int, default=pc.REPS)
    args = parser.parse_args()
    root = Path.cwd().resolve()
    work_root = Path(
        os.environ.get("PERF_WORK_DIR") or tempfile.mkdtemp(prefix="perf-array-")
    )
    work_root.mkdir(parents=True, exist_ok=True)
    showers = performance.select_input(
        args.input.resolve(), args.events, work_root / "showers-input.iact"
    )
    metrics = {}
    load_before = pc.loadavg()
    if load_before is not None:
        metrics["host_loadavg_before"] = load_before
    details = {
        "case": args.case,
        "host": performance.host_info(),
        "settings": {
            "reps": args.reps,
            "timing": "wall and CPU: wait4 rusage per child; peak RSS: sampled /proc VmHWM on Linux",
            "showers": performance.shower_count(showers),
        },
        "threads": {},
        "stages": {},
        "pipelines": {},
        "counts": {},
        "failures": {},
    }
    try:
        run_array(root, work_root, showers, args.config, args.reps, metrics, details)
    except Exception:  # keep the failure in details, do not lose the metrics
        details["failures"]["run"] = traceback.format_exc()[-800:]
        raise
    finally:
        load_after = pc.loadavg()
        if load_after is not None:
            metrics["host_loadavg_after"] = load_after
        metrics["timing_unreliable"], details["timing_reliability"] = (
            pc.timing_reliability(metrics, list(details["pipelines"]))
        )
        walls = {
            label: {
                name: statistics.median(r["wall_s"] for r in records)
                for name, records in per_stage.items()
            }
            for label, per_stage in details["stages"].items()
        }
        if walls:
            pc._plot(root / "wall-by-stage.png", walls)
        common.write_results(metrics, details, root)
        if not details["failures"]:  # keep the logs of a failed run
            shutil.rmtree(work_root, ignore_errors=True)


if __name__ == "__main__":
    main()
