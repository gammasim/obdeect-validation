"""Robustness of the stream pipe chain: obdeect stream cat | ray-tracing trace-stream | stream to-simtel,
fed with damaged or unusual streams (LSTN-01 on telescope 1; stream from 20 showers, zenith 20 deg).

Limits (exact, no tolerance):
  - no variant may crash (signal, exit status >= 128, Python traceback) or hang (timeout 300 s per chain);
  - a damaged input (truncated at half, cut inside the header, row count 70000 > 2^16, record magic
    overwritten, format_major 2) must not succeed: the chain must exit with an error;
  - an event without photons (event 3) must give exactly the photons of the reference run for the other
    events (photon by photon); a telescope without photons must deliver no photons;
  - a stream without events must give no photons; a consumer that closes its input early (head -c 4096)
    must end the pipeline without a hang.
"""

import argparse
import struct
from pathlib import Path

from obdeect_validation.adapters import atmosphere_chain as atm
from obdeect_validation.adapters import common
from obdeect_validation.adapters import endtoend_common as e2e

LST, VERSION, ZENITH_DEG, AZIMUTH_DEG, SEED = "LSTN-01", "7.0.0", 20, 180, 5
PROPAGATION_SEED = 123
HEADER = struct.Struct("<IHHQ")
MAGIC = 0x52545347
PHOTONS, EVENT = 5, 3
#: Variants whose input is damaged; success for them would be a silent wrong result.
DAMAGED = (
    "truncated_half",
    "truncated_header",
    "corrupt_rows",
    "corrupt_magic",
    "format_major",
)
EVENT_WITHOUT_PHOTONS = 3
MIN_PHOTONS = 1000


def records(data):
    """(record offset, type, payload offset, payload length) of each record; stops at a bad magic."""
    position = 0
    while position + HEADER.size <= len(data):
        magic, kind, _, length = HEADER.unpack_from(data, position)
        if magic != MAGIC:
            return
        yield position, kind, position + HEADER.size, length
        position += HEADER.size + length


def first_photons(data):
    """Record offset and payload offset of the first photons record."""
    for position, kind, payload, _ in records(data):
        if kind == PHOTONS:
            return position, payload
    raise ValueError("no photons record in stream")


def cut_half(data):
    del data[len(data) // 2 :]


def cut_header(data):
    del data[20:]


def rows_too_large(data):
    _, payload = first_photons(data)
    data[payload + 12 : payload + 16] = struct.pack(
        "<I", 70000
    )  # rows field, above 2^16


def record_magic(data):
    position, _ = first_photons(data)
    data[position : position + 4] = b"XXXX"


def format_major(data):
    key = (
        data.index(b"format_major") + len(b"format_major") + 1
    )  # skip the type byte (1 = int)
    data[key : key + 8] = struct.pack("<q", 2)


def variant(source, destination, change):
    data = bytearray(Path(source).read_bytes())
    change(data)
    Path(destination).write_bytes(bytes(data))
    return destination


def trace_argv(model, stream_in, stream_out):
    return [
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
        "1",
        "--pointing-zenith-deg",
        ZENITH_DEG,
        "--pointing-azimuth-deg",
        AZIMUTH_DEG,
        "-i",
        stream_in,
        "-o",
        stream_out,
    ]


def chain_argvs(model, camera_map):
    """cat | trace | to-simtel, every stage on stdin/stdout."""
    return [
        ["obdeect", "stream", "cat", "-"],
        trace_argv(model, "-", "-"),
        ["obdeect", "stream", "to-simtel", "--camera-map", camera_map, "-o", "-", "-"],
    ]


def verdict(stages):
    kinds = [e2e.classify(stage) for stage in stages]
    for word in ("hang", "crash", "error"):
        if word in kinds:
            return word
    return "ok"


def run_chain(root, name, stream, model, camera_map, timeout=300):
    stages = e2e.pipeline(
        chain_argvs(model, camera_map),
        stdin=stream,
        output=root / f"{name}.iact",
        log_dir=root / f"{name}.stages",
        timeout=timeout,
    )
    return stages


def run_trace(root, name, stream, model, timeout=300):
    """cat | trace-stream: the photons delivered for the photon comparisons."""
    out = root / f"{name}.trace.gst"
    stages = e2e.pipeline(
        [["obdeect", "stream", "cat", "-"], trace_argv(model, "-", "-")],
        stdin=stream,
        output=out,
        log_dir=root / f"{name}.trace-stages",
        timeout=timeout,
    )
    return stages, out


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
    native = root / "native.gst"
    prop = root / "prop.gst"
    setup = [
        e2e.run_stage(
            ["obdeect", "stream", "from-corsika", source],
            log=root / "from-corsika.log",
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
            log=root / "propagation.log",
        ),
    ]

    variants = {
        "truncated_half": variant(prop, root / "v-truncated-half.gst", cut_half),
        "truncated_header": variant(prop, root / "v-truncated-header.gst", cut_header),
        "corrupt_rows": variant(prop, root / "v-corrupt-rows.gst", rows_too_large),
        "corrupt_magic": variant(prop, root / "v-corrupt-magic.gst", record_magic),
        "format_major": variant(prop, root / "v-format-major.gst", format_major),
    }
    seen = {"event": 0}

    def drop_event_photons(record):
        if record["type"] == "event":
            seen["event"] += 1
        if record["type"] == "photons" and seen["event"] == EVENT_WITHOUT_PHOTONS:
            return []
        return None

    variants["empty_event"] = root / "v-empty-event.gst"
    e2e.rewrite(prop, variants["empty_event"], drop_event_photons)
    variants["telescope_no_photons"] = root / "v-telescope1-no-photons.gst"
    e2e.rewrite(
        prop,
        variants["telescope_no_photons"],
        lambda r: [] if r["type"] == "photons" and r["telescope_id"] == 1 else None,
    )
    variants["no_events"] = root / "v-no-events.gst"
    e2e.rewrite(
        prop,
        variants["no_events"],
        lambda r: (
            []
            if r["type"]
            in ("event", "telescope", "photons", "telescope_end", "event_end")
            else None
        ),
    )

    metrics = {"setup_stage_failures": sum(e2e.classify(s) != "ok" for s in setup)}
    details = {"variants": {}}
    reference_stages = run_chain(root, "reference", prop, model, camera_map)
    details["variants"]["reference"] = {
        "verdict": verdict(reference_stages),
        "stages": [e2e.describe(s) for s in reference_stages],
    }
    metrics["reference_chain_not_ok"] = int(verdict(reference_stages) != "ok")
    _, reference_trace = run_trace(root, "reference", prop, model)
    reference_rows = e2e.photon_rows(reference_trace)

    for name in (
        "truncated_half",
        "truncated_header",
        "corrupt_rows",
        "corrupt_magic",
        "format_major",
        "empty_event",
        "telescope_no_photons",
        "no_events",
    ):
        stages = run_chain(root, name, variants[name], model, camera_map)
        result = verdict(stages)
        details["variants"][name] = {
            "verdict": result,
            "stages": [e2e.describe(s) for s in stages],
        }
        metrics[f"{name}_crash_or_hang"] = int(result in ("crash", "hang"))
        if name in DAMAGED:
            metrics[f"{name}_silent_success"] = int(result == "ok")
        if name in ("empty_event", "telescope_no_photons", "no_events"):
            trace_stages, trace = run_trace(root, name, variants[name], model)
            trace_result = verdict(trace_stages)
            details["variants"][f"{name} (cat|trace)"] = {
                "verdict": trace_result,
                "stages": [e2e.describe(s) for s in trace_stages],
            }
            metrics[f"{name}_trace_crash_or_hang"] = int(
                trace_result in ("crash", "hang")
            )
            rows = e2e.photon_rows(trace)
            if name == "empty_event":
                comparison = e2e.compare_rows(
                    e2e.drop_events(reference_rows, [EVENT_WITHOUT_PHOTONS]), rows
                )
                metrics["empty_event_photon_mismatch"] = comparison["photon_mismatches"]
                metrics["empty_event_empty_comparison"] = int(
                    comparison["photons_compared"] < MIN_PHOTONS
                )
            else:
                metrics[f"{name}_photons_out"] = len(rows["x"]) if "x" in rows else 0

    early = e2e.pipeline(
        [["obdeect", "stream", "cat", "-"], ["head", "-c", "4096"]],
        stdin=prop,
        output=root / "early.bin",
        log_dir=root / "early-stages",
        timeout=120,
    )
    details["early_exit_verdicts"] = [e2e.describe(s) for s in early]
    metrics["early_exit_hang"] = int(any(s.timed_out for s in early))
    common.write_results(metrics, details, root)
    e2e.summary_figure(
        root / "summary.png", metrics, "robustness (damaged and unusual streams)"
    )


if __name__ == "__main__":
    main()
