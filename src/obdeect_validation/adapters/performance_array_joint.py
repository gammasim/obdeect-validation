"""Joint array chain: all six telescopes of the North array in ONE obdeect pipe and ONE sim_telarray run.

Why: the tee candidate of performance_array runs one sim_telarray per telescope, so it writes six
files and has no array trigger. Here the streams of all telescopes stay in one stream:

    from-corsika | propagation --telescopes 1..6 (one QE envelope for all) |
    trace-stream --model k=<model of telescope k> ... | to-simtel (camera map of all six) |
    sim_telarray (Bypass_Optics=3, all six telescopes, the array trigger of sim_telarray)

One output file, one array trigger: the same sim_telarray code as the reference.

QE envelope: propagation takes one efficiency bound per process. The bound is the maximum over the
six telescopes of QE_k / headroom_k (per nm bin). The table that sim_telarray reads for telescope k
is QE_k / bound, so bound * table = QE_k (as in performance_chain.bin_envelope) and the table stays
below headroom_k. Telescopes with a lower QE than the maximum get more photons than they need
(thinned at the bound of the best telescope), which costs CPU time but not accuracy.

Variants (labels): mono = trigger_telescopes=1 and the array_triggers.dat of the simtools config
(any one telescope triggers), stereo = the LST line of that file changed to two LSTs
("Trigger 2 of 1, 2, 3, 4"). Each variant has a reference (one sim_telarray run, Bypass_Optics=0)
and a joint candidate with the same trigger settings; the stored events, the telescopes with data
and the trigger patterns are compared event by event, and the photoelectrons per telescope.

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.performance_array_joint air-shower-array-joint \
        --config cfg/model/7.0.0/CTAO.cfg [--events N] [--reps N] [--variants mono,stereo]
"""

from __future__ import annotations

import argparse
import os
import shutil
import statistics
import tempfile
import traceback
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import chain_optics as co
from obdeect_validation.adapters import chain_showers as cs
from obdeect_validation.adapters import common, performance
from obdeect_validation.adapters import performance_array as pa
from obdeect_validation.adapters import performance_chain as pc

ARRAY = pa.ARRAY
#: Photoelectron bias budget as in performance_array (checks of chain_showers).
VARIANTS = ("mono", "stereo")


def union_files(work, elements):
    """One efficiency bound for all telescopes and the QE / bound table of every QE file."""
    work.mkdir(parents=True, exist_ok=True)
    table_dir = work / "qe-table"
    table_dir.mkdir(parents=True, exist_ok=True)
    qes = {}
    bounds = []
    for _, element in ARRAY:
        info = elements[element]
        qe = pc.sim_telarray_qe(info["qe"])
        qes[info["qe"].name] = qe
        bounds.append(pc.bin_envelope(qe, info["headroom"])[0])
    bound = np.max(bounds, axis=0)
    envelope = work / "qe-envelope.dat"
    pc.write_envelope(envelope, pc.WAVELENGTH_GRID, bound)
    for name, qe in qes.items():
        pc.write_scaled(table_dir / name, pc.WAVELENGTH_GRID, qe / bound)
    return {
        "envelope": envelope,
        "table_dir": table_dir,
        "bound_max": float(bound.max()),
    }


def stereo_triggers(config_dir, path):
    """array_triggers.dat of the config with the LST line changed to a 2-telescope trigger."""
    lines = []
    changed = 0
    for line in (config_dir / "array_triggers.dat").read_text().splitlines():
        if line.strip().lower().startswith("trigger 1 of 1, 2, 3, 4 "):
            line = line.replace("Trigger 1 of", "Trigger 2 of", 1)
            changed += 1
        lines.append(line)
    if changed != 1:
        raise RuntimeError("array_triggers.dat: LST trigger line not found")
    path.write_text("\n".join(lines) + "\n")
    return path


def camera_map(path, elements):
    """to-simtel camera map of all six telescopes (id, mirror class, camera rotation)."""
    path.write_text(
        "# telescope_id mirror_class camera_rotation_deg\n"
        + "".join(
            f"{t} {pa.MIRROR_CLASS} {elements[e]['rotation']}\n" for t, e in ARRAY
        )
    )
    return path


def joint_stages(work, source, inputs, files, models, camera, sim_argv):
    """from-corsika | propagation | trace-stream | to-simtel | sim_telarray for all telescopes."""
    ids = ",".join(str(t) for t, _ in ARRAY)
    propagation = [
        "obdeect", "propagation", "run",
        "--transmission", inputs["atmosphere"],
        "--site-altitude", 2156,
        "--pathlen", 28,
        "--seed", 123,
        "--telescopes", ids,
        "--envelope", f"{pc.ENVELOPE_FACTOR}={files['envelope']}",
        "-i", "-", "-o", "-",
    ]  # fmt: skip
    trace = [
        "obdeect", "ray-tracing", "trace-stream",
        *[arg for t, _ in ARRAY for arg in ("--model", f"{t}={models[t]}")],
        "--seed", co.SEED,
        "--threads", 0,
        "--telescopes", ids,
        "--pointing-zenith-deg", co.ZENITH_DEG,
        "--pointing-azimuth-deg", co.AZIMUTH_DEG,
        "--telescope-transmission", "full",
        "-i", "-", "-o", "-",
    ]  # fmt: skip
    to_simtel = [
        "obdeect", "stream", "to-simtel",
        "--camera-map", camera,
        "--scaled-factor", pc.ENVELOPE_FACTOR,
        "--keep-telpos-height",  # sim_telarray's array trigger needs the real heights
        "-o", "-",
    ]  # fmt: skip
    return [
        ("from_corsika", ["obdeect", "stream", "from-corsika", source]),
        ("propagation", propagation),
        ("trace_stream", trace),
        ("to_simtel", to_simtel),
        ("sim_telarray", sim_argv),
    ]


def event_comparison(reference, candidate):
    """Event by event: stored events, telescopes with data and trigger patterns."""
    ref_keys, cand_keys = set(reference), set(candidate)
    common_keys = ref_keys & cand_keys
    tel_mismatch = pattern_mismatch = 0
    for key in common_keys:
        a, b = reference[key], candidate[key]
        if set(a.get("photoelectrons", {})) != set(b.get("photoelectrons", {})):
            tel_mismatch += 1
        if (
            a["trigger_information"]["trigger_pattern"]
            != b["trigger_information"]["trigger_pattern"]
        ):
            pattern_mismatch += 1
    return {
        "events_reference": len(ref_keys),
        "events_candidate": len(cand_keys),
        "events_common": len(common_keys),
        "events_only_reference": len(ref_keys - cand_keys),
        "events_only_candidate": len(cand_keys - ref_keys),
        "events_telescope_set_differs": tel_mismatch,
        "events_trigger_pattern_differs": pattern_mismatch,
        "telescope_images_reference": sum(
            len(e.get("photoelectrons", {})) for e in reference.values()
        ),
        "telescope_images_candidate": sum(
            len(e.get("photoelectrons", {})) for e in candidate.values()
        ),
    }


def restrict(events, index):
    """Events with only one telescope: photoelectrons are keyed by the 0-based index, camera events by the 1-based telescope id."""
    out = {}
    for key, event in events.items():
        copy = dict(event)
        copy["photoelectrons"] = {
            t: v for t, v in event.get("photoelectrons", {}).items() if int(t) == index
        }
        copy["telescope_events"] = {
            t: v
            for t, v in event.get("telescope_events", {}).items()
            if int(t) == index + 1
        }
        out[key] = copy
    return out


def physics_per_telescope(root, reference, candidate, metrics, details, prefix):
    """cs.compare_camera for each telescope of the joint file against the same telescope of the reference."""
    completed = 1
    for telescope, element in ARRAY:
        key = f"tel{telescope}"
        out = root / "physics" / f"{prefix}{key}-{element}"
        out.mkdir(parents=True, exist_ok=True)
        index = telescope - 1  # sim_telarray telescope index is 0-based
        try:
            comparison = pa.measured(
                cs.compare_camera(
                    restrict(reference, index),
                    restrict(candidate, index),
                    out,
                    {index: element},
                ),
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
    metrics[f"{prefix}physics_completed"] = completed


def run_joint(root, work_root, source, config, reps, variants, metrics, details):
    config = config.resolve()
    config_dir = config.parent
    atmosphere_in = config_dir / pa.ATMOSPHERE_FILE
    atmosphere = work_root / "atmosphere-below-1000nm.dat"
    performance.restrict_table(atmosphere_in, atmosphere)
    inputs = {
        "config": config,
        "models": str(common.models_root()),
        "atmosphere": atmosphere,
        "atmosphere_name": str(atmosphere_in),
    }
    elements = {e: pa.element_inputs(config_dir, e) for _, e in ARRAY}
    for _, element in ARRAY:
        elements[element]["headroom"] = pc.headroom_for(
            elements[element]["qe_variation"]
        )
    files = union_files(work_root / "envelope-union", elements)
    camera = camera_map(work_root / "camera-map.txt", elements)
    models = {}
    for telescope, element in ARRAY:
        models[telescope] = work_root / f"model-tel{telescope}"
        pc.run_timed(
            pa.compile_argv(inputs["models"], element, models[telescope]),
            work_root / f"compile-tel{telescope}.log",
        )
    trigger_files = {"mono": None}
    if "stereo" in variants:
        trigger_files["stereo"] = stereo_triggers(
            config_dir, work_root / "array_triggers_stereo.dat"
        )
    details["settings"].update(
        {
            "input": str(source),
            "config": str(config),
            "telescopes": {str(t): e for t, e in ARRAY},
            "headroom": {e: elements[e]["headroom"] for _, e in ARRAY},
            "envelope_bound_max": files["bound_max"],
            "variants": list(variants),
            "stereo_trigger": "Trigger 2 of 1, 2, 3, 4 (LSTs); the MST lines are unchanged",
            "candidate_sim_telarray": pc.CANDIDATE_SIMTEL
            or "production (SIMTEL_CANDIDATE_BINARY unset)",
        }
    )
    runs, outputs = [], {}
    for variant in variants:
        extra_options = {}
        if trigger_files.get(variant) is not None:
            extra_options["array_triggers"] = trigger_files[variant]
        ref_dir = work_root / f"reference-{variant}"
        ref_dir.mkdir(parents=True, exist_ok=True)
        cand_dir = work_root / f"joint-{variant}"
        cand_dir.mkdir(parents=True, exist_ok=True)

        def reference_step(v=variant, d=ref_dir, o=extra_options):
            options = {**pa.reference_options(str(atmosphere_in)), **o}
            outputs[f"reference_{v}"] = common.sim_telarray(
                d, source, "reference", options,
                config=str(config), defines=(), extra=(f"-I{config.parent}",),
                runner=pc.run_timed,
            )  # fmt: skip
            common.check_atmosphere_log(d / "reference.log", pa.ATMOSPHERE_FILE)

        sim_argv, cand_output = common.sim_telarray_argv(
            cand_dir, "-", "candidate",
            {**pa.candidate_options(inputs["atmosphere_name"]), **extra_options},
            config=str(config), defines=(),
            extra=(f"-I{files['table_dir']}", f"-I{config.parent}"),
            binary=pc.CANDIDATE_SIMTEL,
        )  # fmt: skip
        stages = joint_stages(cand_dir, source, inputs, files, models, camera, sim_argv)
        logs = {name: cand_dir / f"{name}.log" for name, _ in stages}
        outputs[f"joint_{variant}"] = (cand_dir, cand_output)
        if not os.environ.get("PERF_JOINT_NO_REFERENCE"):
            runs.append((f"reference_{variant}", [("sim_telarray", reference_step)]))
        runs.append(
            (
                f"candidate_joint_{variant}",
                [
                    (
                        "pipe",
                        lambda s=stages, o=cand_dir / "candidate.stdout", lg=logs: (
                            pc._checked_pipeline(s, o, lg)
                        ),
                    )
                ],
            )
        )
    ok = pc.run_interleaved(runs, reps, metrics, details)
    for variant in variants:
        ref_label, cand_label = f"reference_{variant}", f"candidate_joint_{variant}"
        ratio_wall = pc.ratio(
            metrics.get(f"{cand_label}_wall_s"), metrics.get(f"{ref_label}_wall_s")
        )
        ratio_cpu = pc.ratio(
            metrics.get(f"{cand_label}_cpu_s"), metrics.get(f"{ref_label}_cpu_s")
        )
        for name, value in (("wall", ratio_wall), ("cpu", ratio_cpu)):
            if value is not None:
                metrics[f"joint_{variant}_{name}_ratio"] = value
        if not (ok.get(ref_label) and ok[cand_label]):
            continue
        cand_dir, cand_output = outputs[f"joint_{variant}"]
        try:
            reference = cs.camera_events(outputs[f"reference_{variant}"])
            candidate = cs.camera_events(cand_output)
            for key, value in event_comparison(reference, candidate).items():
                metrics[f"joint_{variant}_{key}"] = value
            metrics[f"joint_{variant}_camera_only_lines"] = pc._camera_only_lines(
                cand_dir
            )
            physics_per_telescope(
                root, reference, candidate, metrics, details, f"joint_{variant}_"
            )
        except (RuntimeError, OSError, ValueError) as error:  # reported, not hidden
            details["failures"][f"joint_{variant}_checks"] = str(error)[:400]
            metrics[f"joint_{variant}_checks_completed"] = 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", choices=["air-shower-array-joint"])
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=pa.INPUT)
    parser.add_argument("--events", type=int, default=None)
    parser.add_argument("--reps", type=int, default=1)
    parser.add_argument("--variants", default=",".join(VARIANTS))
    args = parser.parse_args()
    variants = tuple(v for v in args.variants.split(",") if v)
    if not set(variants) <= set(VARIANTS):
        parser.error(f"variants must be among {VARIANTS}")
    root = Path.cwd().resolve()
    work_root = Path(
        os.environ.get("PERF_WORK_DIR") or tempfile.mkdtemp(prefix="perf-joint-")
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
            "showers": performance.shower_count(showers),
        },
        "threads": {},
        "stages": {},
        "pipelines": {},
        "counts": {},
        "failures": {},
    }
    try:
        run_joint(
            root, work_root, showers, args.config, args.reps, variants, metrics, details
        )
    except Exception:
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
        if not details["failures"] and not os.environ.get("PERF_JOINT_KEEP_WORK"):
            shutil.rmtree(work_root, ignore_errors=True)


if __name__ == "__main__":
    main()
