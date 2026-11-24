"""Fan-out: one stream copied to several consumers with obdeect stream tee (FIFO branches, concurrent).

(a) Same model on two branches. The atmosphere stream goes to four FIFOs: two byte copies (/bin/cat)
    and two LSTN-01 ray-tracing consumers. Checks: the byte copies equal the input (tee is byte
    transparent), the two traces are photon-identical, and their sim_telarray camera outputs (to-simtel
    with Bypass_Optics=1, the optics-bypass1 setup) are identical event by event.
(b) Two models from one stream. The atmosphere stream goes to a LSTN-01 and a MSTN-01 consumer. Each
    trace is compared with a separate full run (from-corsika, propagation, trace) of the same model.
    Same seeds, so identical by design (the RNG is keyed on seed, stage, event, telescope, photon).

Limits: exact equality for every stream and every camera event. Any difference is a failure. The camera
comparison uses the photo-electron arrays and the trigger flags of the sim_telarray output.
"""

import argparse
from pathlib import Path

from obdeect_validation.adapters import atmosphere_chain as atm
from obdeect_validation.adapters import chain_optics as optics
from obdeect_validation.adapters import chain_showers as cs
from obdeect_validation.adapters import common
from obdeect_validation.adapters import endtoend_common as e2e

LST, MST, VERSION = "LSTN-01", "MSTN-01", "7.0.0"
ZENITH_DEG, AZIMUTH_DEG, SEED = 20, 180, 5
PROPAGATION_SEED = 123
TELESCOPE = "1"  # the camera map of to-simtel covers telescope 1 only
MIN_PHOTONS = 1000


def trace_consumer(model, output, seed=SEED):
    """argv builder of a trace-stream consumer reading the FIFO and writing output."""

    def build(fifo):
        return [
            "obdeect",
            "ray-tracing",
            "trace-stream",
            "--model",
            model,
            "--seed",
            seed,
            "--threads",
            2,
            "--telescopes",
            TELESCOPE,
            "--pointing-zenith-deg",
            ZENITH_DEG,
            "--pointing-azimuth-deg",
            AZIMUTH_DEG,
            "-i",
            fifo,
            "-o",
            output,
        ]

    return build


def compile_model(root, name, zenith):
    model = root / f"{name}.json"
    common.command(
        [
            "obdeect",
            "ray-tracing",
            "compile",
            "--simulation-models",
            common.models_root(),
            "--model",
            name,
            "--version",
            VERSION,
            "--seed",
            1,
            "--zenith-deg",
            zenith,
            "--output",
            model,
        ],
        root / f"compile-{name}.log",
    )
    return model


def full_run(root, name, source, cut, band, model):
    """from-corsika, propagation and trace as files (the separate run for the two-model comparison)."""
    native, prop, out = (
        root / f"{name}-native.gst",
        root / f"{name}-prop.gst",
        root / f"{name}.gst",
    )
    stages = [
        e2e.run_stage(
            ["obdeect", "stream", "from-corsika", source],
            log=root / f"{name}-from-corsika.log",
            output=native,
        ),
        e2e.run_stage(
            [
                "obdeect",
                "propagation",
                "run",
                "--transmission",
                cut,
                "--pathlen",
                "0",
                "--site-altitude",
                atm.SITE_ALTITUDE_M,
                "--seed",
                PROPAGATION_SEED,
                "--mode",
                "weighted",
                "--wavelength-range",
                band[0],
                band[1],
                "-i",
                native,
                "-o",
                prop,
            ],
            log=root / f"{name}-propagation.log",
        ),
        e2e.run_stage(
            [
                "obdeect",
                "ray-tracing",
                "trace-stream",
                "--model",
                model,
                "--seed",
                SEED,
                "--threads",
                2,
                "--telescopes",
                TELESCOPE,
                "--pointing-zenith-deg",
                ZENITH_DEG,
                "--pointing-azimuth-deg",
                AZIMUTH_DEG,
                "-i",
                prop,
                "-o",
                out,
            ],
            log=root / f"{name}-trace.log",
        ),
    ]
    return out, stages


def camera_difference(left, right):
    """Event ids present on one side only, and events whose photo-electron arrays or triggers differ."""
    only = len(set(left) ^ set(right))
    differing = 0
    for key in sorted(set(left) & set(right)):
        a, b = left[key], right[key]
        pa = {
            tel: list(map(float, data["photoelectrons"]))
            for tel, data in a.get("photoelectrons", {}).items()
        }
        pb = {
            tel: list(map(float, data["photoelectrons"]))
            for tel, data in b.get("photoelectrons", {}).items()
        }
        if pa != pb or sorted(a.get("telescope_events", {})) != sorted(
            b.get("telescope_events", {})
        ):
            differing += 1
    return only, differing


def camera_output(root, traced, name, config, extra, camera_rotation):
    iact = root / f"{name}.iact"
    common.command(
        [
            "obdeect",
            "stream",
            "to-simtel",
            "--camera-map",
            optics.camera_map(root, camera_rotation),
            "-o",
            iact,
            traced,
        ],
        root / f"{name}-to-simtel.log",
    )
    return iact, cs.simulate(
        root, iact, name, bypass=1, config=str(config), extra=extra
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=common.VALIDATION / "inputs/gamma-20deg-20-100GeV.iact",
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="sim_telarray config generated by simtools for LSTN-01 7.0.0 (CTAO.cfg)",
    )
    args = parser.parse_args()
    root = Path.cwd().resolve()
    source = args.input.resolve()
    config = args.config.resolve()
    config_dir = config.parent
    cfg_text = (config_dir / "CTAO-LSTN-01.cfg").read_text()
    rotation = optics.camera_rotation_deg(config_dir, cfg_text)
    band = atm.cherenkov_band(
        source.with_suffix(".corsika.log").read_text(errors="replace")
    )
    cut = root / "cut-transmission.dat"
    cut.write_text(
        atm.cut_table(common.TRANSMISSION_TABLE.read_text(), atm.MAX_LAMBDA_NM)
    )
    unit = common.write_zero_depth_table(root, optics.UNIT_ATMOSPHERE)
    extra = [
        "-I" + str(config_dir),
        "-DTELESCOPE=1",
        "-C",
        "maximum_telescopes=1",
        "-C",
        f"atmospheric_transmission={unit}",
        "-C",
        "trigger_telescopes=1",
    ]

    metrics, details = (
        {},
        {
            "inputs": str(source),
            "seed": SEED,
            "propagation_seed": PROPAGATION_SEED,
            "models": [f"{LST} {VERSION}", f"{MST} {VERSION}"],
            "band_nm": band,
        },
    )
    stages = {}
    native = root / "native.gst"
    stages["from-corsika"] = e2e.run_stage(
        ["obdeect", "stream", "from-corsika", source],
        log=root / "from-corsika.log",
        output=native,
    )
    prop = root / "prop.gst"
    stages["propagation"] = e2e.run_stage(
        [
            "obdeect",
            "propagation",
            "run",
            "--transmission",
            cut,
            "--pathlen",
            "0",
            "--site-altitude",
            atm.SITE_ALTITUDE_M,
            "--seed",
            PROPAGATION_SEED,
            "--mode",
            "weighted",
            "--wavelength-range",
            band[0],
            band[1],
            "-i",
            native,
            "-o",
            prop,
        ],
        log=root / "propagation.log",
    )
    lst_model = compile_model(root, LST, ZENITH_DEG)
    mst_model = compile_model(root, MST, ZENITH_DEG)

    # (a) one stream, four FIFO branches: two byte copies, two LST traces.
    branches = [
        ("copy-a", lambda f: ["cat", f], root / "copy-a.gst"),
        ("copy-b", lambda f: ["cat", f], root / "copy-b.gst"),
        ("lst-a", trace_consumer(lst_model, root / "lst-a.out.gst"), None),
        ("lst-b", trace_consumer(lst_model, root / "lst-b.out.gst"), None),
    ]
    tee, consumers = e2e.fan_out(prop, branches, log_dir=root / "fanout-a", timeout=600)
    stages_a = [tee, *consumers]
    metrics["fanout_a_stage_failures"] = sum(e2e.classify(s) != "ok" for s in stages_a)
    details["fanout_a_verdicts"] = [e2e.describe(s) for s in stages_a]
    prop_sha = common.sha256(prop)
    metrics["tee_copies_not_byte_identical"] = sum(
        common.sha256(root / n) != prop_sha for n in ("copy-a.gst", "copy-b.gst")
    )
    metrics["tee_stdout_not_byte_identical"] = int(
        common.sha256(root / "fanout-a/tee-stdout.gst") != prop_sha
    )
    trace_a, trace_b = root / "lst-a.out.gst", root / "lst-b.out.gst"
    rows_a, rows_b = e2e.photon_rows(trace_a), e2e.photon_rows(trace_b)
    result = e2e.compare_rows(rows_a, rows_b)
    metrics["same_model_trace_photon_mismatch"] = result["photon_mismatches"]
    metrics["same_model_empty_comparison"] = int(
        result["photons_compared"] < MIN_PHOTONS
    )
    details["same_model_trace_bytes_equal"] = common.sha256(trace_a) == common.sha256(
        trace_b
    )

    iact_a, sim_a = camera_output(root, trace_a, "lst-a", config, extra, rotation)
    iact_b, sim_b = camera_output(root, trace_b, "lst-b", config, extra, rotation)
    metrics["same_model_simtel_input_bytes_differ"] = int(
        common.sha256(iact_a) != common.sha256(iact_b)
    )
    events_a, events_b = cs.camera_events(sim_a), cs.camera_events(sim_b)
    only, differing = camera_difference(events_a, events_b)
    metrics["same_model_camera_event_id_delta"] = only
    metrics["same_model_camera_event_mismatch"] = differing
    details["camera_events"] = {"lst-a": len(events_a), "lst-b": len(events_b)}

    # (b) one stream, LST and MST consumers; compare with separate full runs.
    (root / "fanout-b").mkdir(exist_ok=True)
    tee_b, consumers_b = e2e.fan_out(
        prop,
        [
            ("lst", trace_consumer(lst_model, root / "fanout-b/lst.out.gst"), None),
            ("mst", trace_consumer(mst_model, root / "fanout-b/mst.out.gst"), None),
        ],
        log_dir=root / "fanout-b",
        timeout=600,
    )
    stages_b = [tee_b, *consumers_b]
    metrics["fanout_b_stage_failures"] = sum(e2e.classify(s) != "ok" for s in stages_b)
    details["fanout_b_verdicts"] = [e2e.describe(s) for s in stages_b]
    for name, model in (("lst", lst_model), ("mst", mst_model)):
        tee_out = root / f"fanout-b/{name}.out.gst"
        separate, run_stages = full_run(root, f"sep-{name}", source, cut, band, model)
        stages_b.extend(run_stages)
        result = e2e.compare_rows(e2e.photon_rows(tee_out), e2e.photon_rows(separate))
        metrics[f"two_model_{name}_tee_vs_separate_photon_mismatch"] = result[
            "photon_mismatches"
        ]
        metrics[f"two_model_{name}_empty_comparison"] = int(
            result["photons_compared"] < MIN_PHOTONS
        )
        details[f"two_model_{name}_bytes_equal"] = common.sha256(
            tee_out
        ) == common.sha256(separate)
        details.setdefault("photons", {})[name] = result["photons_compared"]
    metrics["two_model_stage_failures"] = sum(
        e2e.classify(s) != "ok" for s in stages_b[len(consumers_b) + 1 :]
    )
    details["two_model_separate_verdicts"] = [
        e2e.describe(s) for s in stages_b[len(consumers_b) + 1 :]
    ]
    details["stage_verdicts"] = {name: e2e.describe(s) for name, s in stages.items()}
    common.write_results(metrics, details, root)
    e2e.summary_figure(
        root / "summary.png", metrics, "fan-out (tee branches, models from one stream)"
    )


if __name__ == "__main__":
    main()
