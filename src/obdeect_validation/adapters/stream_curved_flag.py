"""CURVED flag of a flat and a curved 65 deg CORSIKA shower file, through the obdeect stream.

For each file: (1) the CORSIKA header flags read independently with eventio (EVTH word 79, bit
0x40 and bit 0x10 of word 77, see curved_showers.read_header_flags); (2) the stream converted with
`obdeect stream from-corsika` and its first run record (keys curvature_applied and
refraction_applied); (3) the same through the Python API (obdeect.stream.corsika_to_stream) when
it is importable; (4) the CORSIKA log next to the file (<stem>.corsika.log) for the CURVED pattern.

Metrics (all expected 0 unless stated): curved_file_curvature_applied (expected 1),
flat_file_curvature_applied (expected 0), word79_flag_mismatch (run key against word 79, per
conversion route), iact_curved_bit_mismatch (word 77 bit 0x40 against word 79),
refraction_bit_mismatch (refraction_applied against word 77 bit 0x10, per route),
run_key_missing_count (absent run keys or run records, per route), curved_log_pattern_mismatch
(curved log without the pattern, flat log with it, or a missing log).
"""

from __future__ import annotations

from pathlib import Path

from obdeect_validation.adapters import common, curved_showers
from obdeect_validation.adapters import obdeect_cli as cli
from obdeect_validation.adapters.common import command

CURVED_INPUT = "/validation/inputs/" + curved_showers.shower_name(
    65, "curved", 300, 5000
)
FLAT_INPUT = "/validation/inputs/" + curved_showers.shower_name(65, "flat", 300, 5000)
METRICS = (
    "curved_file_curvature_applied",
    "flat_file_curvature_applied",
    "word79_flag_mismatch",
    "iact_curved_bit_mismatch",
    "refraction_bit_mismatch",
    "run_key_missing_count",
    "curved_log_pattern_mismatch",
)


def run_metadata(stream_path):
    """Metadata of the first run record of a photon stream, or None."""
    from obdeect.stream import Reader

    for record in Reader(str(stream_path)):
        if record["type"] == "run":
            return dict(record["metadata"])
    return None


def routes_for(shower, root: Path, label: str) -> dict:
    """Run metadata per conversion route: the CLI, and the Python API when importable."""
    native = root / f"{label}.cli.gst"
    command(
        cli.from_corsika(shower),
        root / f"{label}.cli.from-corsika.log",
        stdout_file=native,
    )
    routes = {"cli": run_metadata(native)}
    try:
        from obdeect.stream import corsika_to_stream
    except ImportError:
        return routes
    converted = root / f"{label}.api.gst"
    corsika_to_stream(shower, converted)
    routes["api"] = run_metadata(converted)
    return routes


def evaluate(items: list[dict]) -> tuple[dict, dict]:
    """Metrics and details from per-file results.

    Each item has kind ("curved" or "flat"), flags (curved_showers.read_header_flags), routes (run
    metadata per route, None when no run record) and log_text (None when the log is missing).
    """
    metrics = {name: 0 for name in METRICS}
    details: dict = {"files": {}, "log_missing": []}
    for item in items:
        kind, flags = item["kind"], item["flags"]
        expected = bool(flags["corsika_curved_flag"])
        if bool(flags["iact_curved_bit"]) != expected:
            metrics["iact_curved_bit_mismatch"] += 1
        values = {}
        for route, metadata in item["routes"].items():
            if metadata is None or "curvature_applied" not in metadata:
                metrics["run_key_missing_count"] += 1
                values[route] = None
            else:
                values[route] = int(metadata["curvature_applied"])
                if bool(values[route]) != expected:
                    metrics["word79_flag_mismatch"] += 1
            if metadata is None or "refraction_applied" not in metadata:
                metrics["run_key_missing_count"] += 1
            elif bool(metadata["refraction_applied"]) != bool(flags["refraction_bit"]):
                metrics["refraction_bit_mismatch"] += 1
        cli_value = values.get("cli")
        if kind == "curved":
            metrics["curved_file_curvature_applied"] = (
                -1 if cli_value is None else cli_value
            )
        else:
            metrics["flat_file_curvature_applied"] = (
                -1 if cli_value is None else cli_value
            )
        log_text = item.get("log_text")
        if log_text is None:
            details["log_missing"].append(kind)
            metrics["curved_log_pattern_mismatch"] += 1
        elif curved_showers.log_is_curved(log_text) != (kind == "curved"):
            metrics["curved_log_pattern_mismatch"] += 1
        details["files"][kind] = {
            "flags": flags,
            "run_keys": {
                route: None
                if metadata is None
                else {
                    key: metadata.get(key)
                    for key in (
                        "curvature_applied",
                        "refraction_applied",
                        "atmosphere_id",
                    )
                }
                for route, metadata in item["routes"].items()
            },
            "log": item.get("log"),
        }
    return metrics, details


def analyse(root: Path, args) -> tuple[dict, dict]:
    root = root.resolve()
    items = []
    for kind, shower in (("curved", args.curved), ("flat", args.flat)):
        shower = shower.resolve()
        log = shower.with_suffix(".corsika.log")
        items.append(
            {
                "kind": kind,
                "flags": curved_showers.read_header_flags(shower),
                "routes": routes_for(shower, root, kind),
                "log_text": log.read_text(errors="replace") if log.is_file() else None,
                "log": str(log),
            }
        )
    return evaluate(items)


def add_arguments(parser):
    parser.add_argument("--curved", default=CURVED_INPUT, type=Path)
    parser.add_argument("--flat", default=FLAT_INPUT, type=Path)


def main(argv=None) -> None:
    common.results_main(
        analyse,
        __doc__.splitlines()[0],
        add_arguments=add_arguments,
        argv=argv,
    )


if __name__ == "__main__":
    main()
