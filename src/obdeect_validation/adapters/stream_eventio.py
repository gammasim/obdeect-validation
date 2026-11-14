"""Check stream-to-sim_telarray conversion with the independent eventio reader."""

import itertools
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
from eventio import EventIOFile
from eventio.exceptions import WrongType
from eventio.iact import IACTFile
from obdeect.stream import Reader, Writer

from obdeect_validation.adapters import common

FIXTURE = Path(__file__).with_name("compare-led-camera.txt")
HEADER_ORDER = [
    1200,
    1212,
    1216,
    1201,
]  # RUNH, INPUTCFG, ATMPROF, TELPOS as CORSIKA writes them


def _top_level_spans(path):
    """Byte span (start, end) and type of each top-level EventIO object."""
    spans = []
    with EventIOFile(str(path)) as file:
        for obj in file:
            header = obj.header
            start = header.content_address - header.header_size
            spans.append(
                (header.type, start, header.content_address + header.content_size)
            )
    return spans


def _object_types(path):
    return [span[0] for span in _top_level_spans(path)]


def _reorder_input_card(path, output_path):
    """Write a copy whose header blocks follow CORSIKA's order (HEADER_ORDER), rest unchanged.

    Diagnostic only: the converted file itself is not changed.
    """
    data = Path(path).read_bytes()
    spans = _top_level_spans(path)
    if (
        spans[0][1] != 0
        or spans[-1][2] != len(data)
        or any(a[2] != b[1] for a, b in itertools.pairwise(spans))
    ):
        raise ValueError("top-level EventIO objects do not cover the file contiguously")
    head = [
        next(i for i, span in enumerate(spans) if span[0] == block_type)
        for block_type in HEADER_ORDER
        if any(span[0] == block_type for span in spans)
    ]
    order = head + [i for i in range(len(spans)) if i not in head]
    Path(output_path).write_bytes(
        b"".join(data[spans[i][1] : spans[i][2]] for i in order)
    )


def _read_stream(path):
    counts = []
    for record in Reader(str(path)):
        if record["type"] == "event":
            counts.append(0.0)
        elif record["type"] == "photons":
            if not counts:
                raise ValueError("photon stream has photons before its first event")
            counts[-1] += float(np.sum(record["columns"]["w"], dtype=np.float64))
    return counts


def _read_eventio(path):
    counts = []
    with IACTFile(path) as file:
        for event in file:
            counts.append(
                sum(
                    float(np.sum(bunches["photons"], dtype=np.float64))
                    for bunches in event.photon_bunches.values()
                )
            )
    return counts


def _add_zem(input_path, output_path):
    """Add the required zero emission-height field to camera-frame LED photons."""
    with Writer(
        str(output_path), "obdeect validation stream interoperability"
    ) as writer:
        for record in Reader(str(input_path)):
            if record["type"] == "stream":
                continue
            if record["type"] == "photons":
                columns = dict(record["columns"])
                rows = len(columns["x"])
                columns["zem"] = np.zeros(rows, dtype=np.float32)
                writer.write_photons(
                    record["telescope_id"],
                    record["frame"],
                    record["first_index"],
                    columns,
                )
            elif record["type"] == "blob":
                writer.write_blob(
                    record["eventio_type"], record["eventio_id"], record["data"]
                )
            else:
                writer.write(record["type"], record["metadata"])


def main():
    obdeect = shutil.which("obdeect")
    config = os.environ.get("SIMTEL_CONFIG_PATH")
    profile = Path(config) / "atm_profile_model_1.dat" if config else None
    if obdeect is None:
        raise RuntimeError("obdeect must be installed")
    if profile is None or not profile.is_file():
        raise RuntimeError("SIMTEL_CONFIG_PATH must contain atm_profile_model_1.dat")

    root = Path.cwd()
    camera = root / "camera.txt"
    source = root / "sources.gst"
    prepared = root / "prepared.gst"
    converted = root / "converted.dat.gz"
    camera_map = root / "camera.map"
    shutil.copyfile(FIXTURE, camera)
    camera_map.write_text("1 1 0\n", encoding="utf-8")
    subprocess.run(
        [
            obdeect,
            "sources",
            "led",
            "--camera-file",
            str(camera),
            "--use",
            "0 1 2",
            "--require",
            "2",
            "--photons",
            "10",
            "--seed",
            "17",
            "-o",
            str(source),
        ],
        check=True,
    )
    _add_zem(source, prepared)
    subprocess.run(
        [
            obdeect,
            "stream",
            "to-simtel",
            "--camera-map",
            str(camera_map),
            "--atmprof",
            str(profile),
            str(prepared),
            "-o",
            str(converted),
        ],
        check=True,
    )

    source_counts = _read_stream(prepared)
    if not source_counts or sum(source_counts) <= 0:
        raise ValueError("obdeect stream contains no photon events")
    object_types = _object_types(converted)
    try:
        reference_counts = _read_eventio(converted)
        reader_error = None
    except (
        WrongType
    ) as exc:  # The reference reader rejects the emitted EventIO block order.
        reference_counts = []
        reader_error = f"{type(exc).__name__}: {exc}"
    # Diagnostic: the same photons with INPUTCFG moved to CORSIKA's position (after RUNH).
    reordered = root / "reordered.dat"
    _reorder_input_card(converted, reordered)
    reordered_counts = _read_eventio(reordered)
    deltas = common.event_deltas(source_counts, reference_counts)
    metrics = {
        "eventio_readable": int(reader_error is None),
        "source_event_count": len(source_counts),
        "eventio_event_count": len(reference_counts),
        "event_count_delta": abs(len(source_counts) - len(reference_counts)),
        "source_photon_count": sum(source_counts),
        "eventio_photon_count": sum(reference_counts),
        "photon_count_delta": abs(sum(source_counts) - sum(reference_counts)),
        "max_event_photon_delta": max(deltas),
        "reordered_event_count": len(reordered_counts),
        "reordered_photon_count": sum(reordered_counts),
        "reordered_photon_count_delta": abs(sum(source_counts) - sum(reordered_counts)),
    }
    provenance = {
        "eventio_reader_error": reader_error,
        "converted_object_types": object_types,
        "corsika_reference_order": "1200 RunHeader, 1212 InputCard, 1216 AtmProfile, 1201 TelDef",
        "reordered_object_types": _object_types(reordered),
        "seed": 17,
        "source_event_photons": source_counts,
        "eventio_event_photons": reference_counts,
        "camera_sha256": common.sha256(camera),
        "atmosphere_sha256": common.sha256(profile),
    }
    common.write_results(metrics, provenance)


if __name__ == "__main__":
    main()
