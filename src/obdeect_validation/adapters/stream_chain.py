"""CORSIKA round trip through the obdeect stream, and stream command robustness variants.

Round trip: CORSIKA IACT -> from-corsika -> to-simtel, compared bunch by bunch and through
sim_telarray camera images (reuses chain_showers helpers).
Variants: the same stream through cat/tee/stdin pipes must give identical photons per
event and telescope (canonical digest of photons sorted by all columns, bit-exact).
"""

import hashlib
import subprocess
from pathlib import Path

import numpy as np
from obdeect.stream import Reader

from obdeect_validation.adapters import chain_showers as chain
from obdeect_validation.adapters import common
from obdeect_validation.adapters import obdeect_cli as cli


def run_stdin(argv, stdin_path, stdout_path, log):
    with (
        Path(stdin_path).open("rb") as stdin,
        Path(stdout_path).open("wb") as output,
        Path(log).open("wb") as err,
    ):
        result = subprocess.run(
            list(map(str, argv)), stdin=stdin, stdout=output, stderr=err, check=False
        )
    if result.returncode:
        raise RuntimeError(f"{argv[0]} exited {result.returncode}; see {log}")


def canonical(path):
    """Digest per 'event:telescope' of the photons, sorted by all columns."""
    digests, event, groups = {}, -1, {}

    def close():
        for telescope, parts in groups.items():
            columns = {
                name: np.concatenate([part[name] for part in parts])
                for name in parts[0]
            }
            order = np.lexsort(
                [columns[name] for name in sorted(columns, reverse=True)]
            )
            digest = hashlib.sha256()
            for name in sorted(columns):
                digest.update(np.ascontiguousarray(columns[name][order]).tobytes())
            digests[f"{event}:{telescope}"] = digest.hexdigest()

    for record in Reader(str(path)):
        if record["type"] == "event":
            close()
            groups, event = {}, event + 1
        elif record["type"] == "photons":
            groups.setdefault(record["telescope_id"], []).append(record["columns"])
    close()
    return digests


def differing(left, right):
    return len(set(left.items()) ^ set(right.items()))


def variants(root, source, native):
    """Stream commands with different input/output paths; each must reproduce the same photons."""
    reference = canonical(native)
    outputs = {
        "cat_file": root / "var-cat-file.gst",
        "cat_stdin": root / "var-cat-stdin.gst",
        "tee": root / "var-tee.gst",
        "from_corsika_stdin": root / "var-from-corsika-stdin.gst",
    }
    common.command(
        cli.obdeect("stream", "cat", native),
        root / "var-cat-file.log",
        stdout_file=outputs["cat_file"],
    )
    run_stdin(
        cli.obdeect("stream", "cat", "-"),
        native,
        outputs["cat_stdin"],
        root / "var-cat-stdin.log",
    )
    run_stdin(
        cli.obdeect("stream", "tee", outputs["tee"]),
        native,
        root / "var-tee-stdout.gst",
        root / "var-tee.log",
    )
    run_stdin(
        cli.from_corsika("-"),
        source,
        outputs["from_corsika_stdin"],
        root / "var-from-corsika.log",
    )
    metrics = {
        f"variant_{name}_photon_mismatch": differing(reference, canonical(path))
        for name, path in outputs.items()
    }
    metrics["variant_tee_stdout_photon_mismatch"] = differing(
        reference, canonical(root / "var-tee-stdout.gst")
    )
    to_simtel_stdin = root / "var-to-simtel-stdin.iact"
    run_stdin(
        cli.obdeect("stream", "to-simtel", o=to_simtel_stdin, extra=["-"]),
        native,
        "/dev/null",
        root / "var-to-simtel.log",
    )
    to_simtel_file = root / "var-to-simtel-file.iact"
    common.command(
        cli.to_simtel(native, to_simtel_file),
        root / "var-to-simtel-file.log",
    )
    photons, _ = chain.compare_photons(to_simtel_file, to_simtel_stdin)
    metrics["variant_to_simtel_stdin_key_delta"] = photons[
        "photon_event_telescope_key_delta"
    ]
    metrics["variant_to_simtel_stdin_bunch_delta"] = photons["bunch_count_delta"]
    metrics["variant_to_simtel_stdin_codec_bounds_exceeded"] = photons[
        "codec_bounds_exceeded"
    ]
    metrics["variant_to_simtel_stdin_bytes_differ"] = int(
        common.sha256(to_simtel_file) != common.sha256(to_simtel_stdin)
    )
    return metrics


def analyse(root: Path, args) -> tuple[dict, dict]:
    root = root.resolve()
    source = args.input.resolve()
    native, metrics, fields, a, b = chain.roundtrip(root, source)
    metrics.update(variants(root, source, native))
    return metrics, {
        "codec_per_field_normalized_error": fields,
        "reference_event_ids": list(a),
        "candidate_event_ids": list(b),
        "variant_commands": [
            "stream cat FILE",
            "stream cat -",
            "stream tee FILE",
            "stream from-corsika -",
            "stream to-simtel -o FILE -",
        ],
    }


def add_arguments(parser):
    parser.add_argument("--input", type=Path, required=True, help="CORSIKA IACT file")


def main():
    common.results_main(
        analyse,
        "CORSIKA IACT round trip and stream command variants",
        add_arguments=add_arguments,
    )


if __name__ == "__main__":
    main()
