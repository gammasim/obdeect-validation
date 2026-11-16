"""Multi-telescope air-shower timing and camera check on the North array (4 LSTs, 2 MSTs).

Input: the array shower file (CORSIKA telescopes 1-6 = LSTN-01..04, MSTN-03, MSTN-15, at their
simtools positions; 20 deg, 100-3000 GeV, CSCAT 25 km... see shower_matrix.ARRAY_NORTH_TELESCOPES).
The config is the simtools array config of the same six elements (CTAO.cfg with CTAO-<element>.cfg).

Reference (timed): one sim_telarray run with all six telescopes (no -DTELESCOPE, so sim_telarray
reads every telescope block of the config), Bypass_Optics=0 (full optics), production atmosphere.
Candidate (timed): for each telescope k, one pipe as in performance_chain.envelope_pipeline:
from-corsika | propagation --telescopes k (QE envelope, thinned) | trace-stream (model of element k,
compiled once per telescope) | to-simtel | sim_telarray -DTELESCOPE=k (Bypass_Optics=3, camera only).
The chain time is the sum of the pipe stages over the telescopes; the compile stages are reported
and left out, as in performance_budget. The candidate runs the telescopes one after another.

Deviations from a single joint run, recorded in details.json:
  * No joint candidate run: the per-telescope camera streams are not merged into one sim_telarray
    input (no event merge tool), so the array trigger is not compared. Both sides use
    trigger_telescopes=1, so every telescope's photoelectrons are stored for the same events.
  * No tee: each telescope pipe reads the CORSIKA file and propagates again. The candidate CPU is
    overstated against a tee of the propagated stream (the design's fan-out).
  * Model compile uses --zenith-deg 20 for all six (the MSTs fail without it: axes_offsets).
  * min_photons = 0 on both sides (performance_chain.MIN_PHOTONS). sim_telarray skips a telescope
    event below min_photons (120 LST, 450 MST in the config) using the photons of its bunches; a
    thinned stream has one photon per surviving bunch, so the threshold would drop events on the
    candidate side only. min_photoelectrons (30 MST, 20 LST) is kept and recorded in details.

Limits (fixed here):
  * Per telescope: photoelectron bias within max(1 %, 5 sigma) (chain_showers.compare_camera,
    PE_BIAS_FLOOR and SIGMA of chain_showers). Counted in equal_statistics_budget_exceeded.
  * Chain wall and CPU (sum over telescope pipes, compiles excluded) <= reference wall and CPU.
  * Every telescope's sim_telarray prints the camera-only line of Bypass_Optics=3.
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
import tempfile
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
    """Timed steps of the candidate: per telescope, compile_tel<k> and pipe_tel<k>."""
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


def photon_checks(outputs, metrics, details):
    """Per telescope: no photon dropped before the camera, and to-simtel bunches = trace rows.

    to-simtel bunch count is taken from sim_telarray's "All photons" (the sum of the bunch photon
    counts it reads; the to-simtel stream is piped and not stored). Each check is 0 = pass.
    """
    details["photon_checks"] = {}
    failed = 0
    for telescope, _ in ARRAY:
        work, _ = outputs[telescope]
        key = f"tel{telescope}"
        all_photons, used = photon_statistics(
            (work / "candidate.stdout").read_text(errors="replace"), telescope
        )
        rows = trace_out_rows((work / "trace_stream.log").read_text(errors="replace"))
        flags = {
            f"photon_used_ne_all_{key}": int(used != all_photons),
            f"to_simtel_bunch_ne_trace_rows_{key}": int(all_photons != rows),
        }
        metrics.update(
            {
                f"photon_all_{key}": all_photons,
                f"photon_used_{key}": used,
                f"to_simtel_bunches_{key}": all_photons,
                f"trace_out_rows_{key}": rows,
                **flags,
            }
        )
        failed += sum(flags.values())
        details["photon_checks"][key] = {
            "all": all_photons,
            "used": used,
            "trace_rows": rows,
        }
    metrics["photon_checks_failed"] = failed


def measured(comparison, key, details):
    """A comparison with fewer than MIN_MATCHED_SAMPLES matched samples is not measured (flag 1)."""
    matched = comparison.get("matched_event_telescope_samples", 0)
    if matched >= MIN_MATCHED_SAMPLES:
        return comparison
    details["failures"][f"physics_{key}"] = (
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


def physics_check(root, reference_output, outputs, metrics, details):
    """Per telescope: reference camera events against the candidate events of that telescope.

    A telescope with no candidate events (or no matched samples) is not compared: its bias flag
    is set to 1 (exceeded) and the reason goes to details, so it cannot pass silently.
    """
    reference = cs.camera_events(reference_output)
    details["physics"] = {}
    completed = 1
    for telescope, element in ARRAY:
        work, output = outputs[telescope]
        key = f"tel{telescope}"
        out = root / "physics" / f"{key}-{element}"
        out.mkdir(parents=True, exist_ok=True)
        index = telescope - 1  # sim_telarray telescope index is 0-based
        try:
            candidate = cs.camera_events(output)
            comparison = measured(
                cs.compare_camera(reference, candidate, out, {index: element}),
                key,
                details,
            )
        except ValueError as error:  # no matched samples: exceeded, reported
            completed = 0
            metrics[f"physics_{key}_photoelectron_bias_limit_exceeded"] = 1
            details["failures"][f"physics_{key}"] = str(error)[:200]
            comparison = {}
        metrics.update({f"physics_{key}_{k}": v for k, v in comparison.items()})
        metrics[f"{key}_camera_only_lines"] = pc._camera_only_lines(work)
        metrics[f"{key}_camera_only_missing"] = int(
            metrics[f"{key}_camera_only_lines"] == 0
        )
        details["physics"][key] = comparison
    metrics["physics_completed"] = completed


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
    steps, outputs = telescope_steps(work_root, source, inputs, files, elements)
    candidate_ok = pc._run("candidate_envelope", steps, reps, metrics, details)
    if candidate_ok:
        stages = details["stages"]["candidate_envelope"]
        chain = {
            "wall_s": chain_medians(stages, "wall_s"),
            "cpu_s": chain_medians(stages, "cpu_s"),
        }
        metrics["array_chain_wall_s"] = chain["wall_s"]
        metrics["array_chain_cpu_s"] = chain["cpu_s"]
    else:
        chain = {}
    if candidate_ok:
        try:
            photon_checks(outputs, metrics, details)
        except (RuntimeError, OSError, ValueError) as error:  # reported, not hidden
            details["failures"]["photon_checks"] = str(error)[:400]
    if reference_ok and candidate_ok:
        try:
            physics_check(root, reference_output["path"], outputs, metrics, details)
        except (RuntimeError, OSError, ValueError) as error:  # reported, not hidden
            metrics["physics_completed"] = 0
            details["failures"]["physics"] = str(error)[:400]
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
