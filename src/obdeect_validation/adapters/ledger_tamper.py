"""Efficiency ledger end to end: tampered and untampered chains through obdeect stream to-simtel and
trace-stream, built with the obdeect.stream Python API.

Each must-reject case must exit with an error (exit status 1 or 2, a message, no traceback) that names the
factor. Each positive control must be accepted. Limits: zero accepted tampered chains, zero rejections
without the factor in the message, zero rejected positive controls, zero crashes or hangs.

Tampered chains (LSTN-01, telescope 1, 20 showers, zenith 20 deg):
  missing_atmosphere   from-corsika -> trace-stream -> to-simtel; no propagation stage, so no
                       atmosphere record (camera frame needs atmosphere applied upstream)
  camera_mirror_pending  camera-frame stream with mirror_reflectivity set back to pending
  ground_mirror_applied  ground stream with mirror_reflectivity already applied (sim_telarray applies it again)
  ground_tt_applied      ground stream with telescope_transmission already applied
  qe_envelope_undeclared ground stream with qe at envelope, no --scaled-factor
  qe_applied             ground stream with qe already applied
  trace_on_applied_mirror  ground stream with mirror_reflectivity applied, traced again
  angular_t0_undeclared  telescope_transmission_constant pending, no --allow-pending
  backward_state         qe envelope then pending in a second run record (state moves backwards)
  retrace_camera_frame   camera-frame stream traced again (frame rejected)
Positive controls: genuine camera chain; qe envelope with --scaled-factor qe; angular with --allow-pending.
"""

import argparse
from pathlib import Path

from obdeect_validation.adapters import atmosphere_chain as atm
from obdeect_validation.adapters import common
from obdeect_validation.adapters import endtoend_common as e2e

LST, VERSION, ZENITH_DEG, AZIMUTH_DEG, SEED = "LSTN-01", "7.0.0", 20, 180, 5
PROPAGATION_SEED = 123
SOURCE_NAME = "endtoend-ledger tamper"


def status_editor(name, state, *, add=False):
    """Edit function for rewrite(): set factor name to state in the efficiency_status records that name it.

    add=True also adds the factor to the first record that does not name it (a run record of a stream
    that never had the factor).
    """

    def edit(record):
        if record["type"] != "efficiency_status":
            return None
        metadata = dict(record["metadata"])
        named = name in str(metadata.get("factors", "")).split(",")
        if not named and not add:
            return None
        return [
            (
                "meta",
                "efficiency_status",
                e2e.with_status_factor(metadata, name, state, SOURCE_NAME),
            )
        ]

    return edit


def backward_editor():
    """Keep the run status and add a qe-only envelope record (valid; see write_backward)."""

    def edit(record):
        if record["type"] != "efficiency_status":
            return None
        metadata = dict(record["metadata"])
        envelope = {"factors": "qe", "qe.state": "envelope", "qe.source": SOURCE_NAME}
        return [
            ("meta", "efficiency_status", metadata),
            ("meta", "efficiency_status", envelope),
        ]

    return edit


def write_backward(source, destination):
    """The backward-state stream: qe envelope, then qe pending (a backward step).

    The obdeect Writer refuses to write an invalid status transition, so the stream is written up
    to the envelope record by the Writer and the pending record is spliced in as raw bytes after
    it (the tampering under test).
    """
    from obdeect_validation.adapters import robustness_rejections as raw

    valid = destination.with_name(destination.name + ".valid")
    e2e.rewrite(source, valid, backward_editor())
    records = raw.split_records(valid.read_bytes())
    status = raw.RECORD_TYPE["efficiency_status"]
    last = max(i for i, (kind, _) in enumerate(records) if kind == status)
    pending = {"factors": "qe", "qe.state": "pending", "qe.source": SOURCE_NAME}
    spliced = raw.raw_record(status, raw.encode_metadata(pending))
    data = [bytes(r) for _, r in records]
    data.insert(last + 1, spliced)
    destination.write_bytes(b"".join(data))
    valid.unlink()
    return destination


def tool_argv_to_simtel(camera_map, stream, output, flags=()):
    """to-simtel; camera_map None for ground-frame streams (telescope 1 is then a ground telescope)."""
    mapping = ["--camera-map", camera_map] if camera_map else []
    return ["obdeect", "stream", "to-simtel", *mapping, *flags, stream, "-o", output]


def tool_argv_trace(model, stream, output):
    return [
        "obdeect",
        "ray-tracing",
        "trace-stream",
        "--model",
        model,
        "--seed",
        SEED,
        "--telescopes",
        "1",
        "--pointing-zenith-deg",
        ZENITH_DEG,
        "--pointing-azimuth-deg",
        AZIMUTH_DEG,
        "-i",
        stream,
        "-o",
        output,
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=common.VALIDATION / "inputs/gamma-20deg-20-100GeV.iact",
    )
    args = parser.parse_args()
    root = Path.cwd().resolve()
    source = args.input.resolve()
    band = atm.cherenkov_band(
        source.with_suffix(".corsika.log").read_text(errors="replace")
    )
    cut = root / "cut-transmission.dat"
    cut.write_text(
        atm.cut_table(common.TRANSMISSION_TABLE.read_text(), atm.MAX_LAMBDA_NM)
    )
    camera_map = root / "camera-map.txt"
    camera_map.write_text("1 1 0.0\n")
    model = root / "model.json"
    common.command(
        [
            "obdeect",
            "ray-tracing",
            "compile",
            "--simulation-models",
            common.models_root(),
            "--model",
            LST,
            "--version",
            VERSION,
            "--seed",
            1,
            "--zenith-deg",
            ZENITH_DEG,
            "--output",
            model,
        ],
        root / "compile.log",
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
            123,
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
    # Streams with the chain's own records and with the tampered records.
    stages["trace-from-prop"] = e2e.run_stage(
        tool_argv_trace(model, prop, root / "lst.gst"), log=root / "trace-lst.log"
    )
    camera = root / "lst.gst"
    stages["trace-angular"] = e2e.run_stage(
        [
            *tool_argv_trace(model, prop, root / "lst-angular.gst"),
            "--telescope-transmission",
            "angular",
        ],
        log=root / "trace-angular.log",
    )
    angular = root / "lst-angular.gst"
    no_atmosphere = root / "nat-lst.gst"
    stages["trace-no-atmosphere"] = e2e.run_stage(
        tool_argv_trace(model, native, no_atmosphere), log=root / "trace-no-atm.log"
    )

    e2e.rewrite(
        camera,
        root / "t-camera-mirror-pending.gst",
        status_editor("mirror_reflectivity", "pending"),
    )
    e2e.rewrite(
        prop,
        root / "t-ground-mirror-applied.gst",
        status_editor("mirror_reflectivity", "applied", add=True),
    )
    e2e.rewrite(
        prop,
        root / "t-ground-tt-applied.gst",
        status_editor("telescope_transmission", "applied", add=True),
    )
    e2e.rewrite(
        prop, root / "t-qe-envelope.gst", status_editor("qe", "envelope", add=True)
    )
    e2e.rewrite(
        prop, root / "t-qe-applied.gst", status_editor("qe", "applied", add=True)
    )
    write_backward(prop, root / "t-backward.gst")
    e2e.rewrite(
        prop,
        root / "t-trace-applied-mirror.gst",
        status_editor("mirror_reflectivity", "applied", add=True),
    )

    # (name, argv, expected factor in the message, must_reject)
    checks = [
        (
            "missing_atmosphere",
            tool_argv_to_simtel(camera_map, no_atmosphere, root / "o.iact"),
            "atmosphere",
            True,
        ),
        (
            "camera_mirror_pending",
            tool_argv_to_simtel(
                camera_map, root / "t-camera-mirror-pending.gst", root / "o.iact"
            ),
            "mirror_reflectivity",
            True,
        ),
        (
            "ground_mirror_applied",
            tool_argv_to_simtel(
                None, root / "t-ground-mirror-applied.gst", root / "o.iact"
            ),
            "mirror_reflectivity",
            True,
        ),
        (
            "ground_tt_applied",
            tool_argv_to_simtel(
                None, root / "t-ground-tt-applied.gst", root / "o.iact"
            ),
            "telescope_transmission",
            True,
        ),
        (
            "qe_envelope_undeclared",
            tool_argv_to_simtel(None, root / "t-qe-envelope.gst", root / "o.iact"),
            "qe",
            True,
        ),
        (
            "qe_applied",
            tool_argv_to_simtel(None, root / "t-qe-applied.gst", root / "o.iact"),
            "qe",
            True,
        ),
        (
            "trace_on_applied_mirror",
            tool_argv_trace(model, root / "t-trace-applied-mirror.gst", root / "o.gst"),
            "mirror_reflectivity",
            True,
        ),
        (
            "angular_t0_undeclared",
            tool_argv_to_simtel(camera_map, angular, root / "o.iact"),
            "telescope_transmission_constant",
            True,
        ),
        (
            "backward_state",
            tool_argv_to_simtel(None, root / "t-backward.gst", root / "o.iact"),
            "qe",
            True,
        ),
        (
            "retrace_camera_frame",
            tool_argv_trace(model, camera, root / "o.gst"),
            # Either reason proves the stream was traced already: the run's efficiency status
            # (read first) has the optics factors applied, or the telescope frame is the camera.
            ("frame", "already applied"),
            True,
        ),
        (
            "control_camera_chain",
            tool_argv_to_simtel(camera_map, camera, root / "o.iact"),
            None,
            False,
        ),
        (
            "control_qe_scaled",
            tool_argv_to_simtel(
                None,
                root / "t-qe-envelope.gst",
                root / "o.iact",
                flags=("--scaled-factor", "qe"),
            ),
            None,
            False,
        ),
        (
            "control_angular_allow_pending",
            tool_argv_to_simtel(
                camera_map,
                angular,
                root / "o.iact",
                flags=("--allow-pending", "telescope_transmission_constant"),
            ),
            None,
            False,
        ),
    ]
    metrics, details = (
        {},
        {
            "stage_verdicts": {n: e2e.describe(s) for n, s in stages.items()},
            "checks": {},
        },
    )
    metrics["setup_stage_failures"] = sum(
        e2e.classify(s) != "ok" for s in stages.values()
    )
    crashes = 0
    for name, argv, factor, must_reject in checks:
        stage = e2e.run_stage(argv, log=root / f"{name}.log", timeout=120)
        verdict = e2e.classify(stage)
        details["checks"][name] = {"verdict": verdict, "message": e2e.describe(stage)}
        if verdict in ("crash", "hang"):
            crashes += 1
        if must_reject:
            metrics[f"{name}_accepted"] = int(verdict == "ok")
            reasons = factor if isinstance(factor, tuple) else (factor,)
            metrics[f"{name}_no_factor_in_message"] = int(
                verdict != "ok"
                and not any(reason in stage.stderr for reason in reasons)
            )
        else:
            metrics[f"{name}_rejected"] = int(verdict != "ok")
    metrics["crash_or_hang_total"] = crashes
    common.write_results(metrics, details, root)
    e2e.summary_figure(
        root / "summary.png", metrics, "efficiency ledger (tampered chains)"
    )


if __name__ == "__main__":
    main()
