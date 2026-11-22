"""Photon directions over the full sphere through the photon stream (format 1.1).

Round trip: photons with upward, horizontal and downward cz are written to a .pst stream in each
frame (sky, ground, telescope, camera), read back, and every column is compared bit for bit with
the input (float32 bit patterns). Counts of the three cz signs must match the input.

Rejection: `obdeect stream to-simtel` cannot store the sign of cz, so it must refuse a frame with
an upward or horizontal photon, with the message of the stream. Per frame, a stream with one
upward and one horizontal photon must fail; a stream with two downward photons (the control, same
layout) must convert. The reference is the input itself; sim_telarray is not needed.
"""

from __future__ import annotations

import argparse
import importlib.metadata
from pathlib import Path

import numpy as np
from obdeect.stream import FRAMES, Reader, Writer, format_major, format_minor

from obdeect_validation.adapters import common, plotting

GRID_SIZE = 1000  # full-sphere directions per frame
UP, DOWN, FIRST_HORIZONTAL, HORIZONTAL = 0, 1, 2, 16  # exact directions in the grid
GOLDEN_ANGLE = np.pi * (3.0 - np.sqrt(5.0))
SEED = 20261122
ROUND_TRIP_FRAMES = ("sky", "ground", "telescope", "camera")
SIMTEL_FRAMES = ("ground", "telescope", "camera")
#: Frames whose cz sign sim_telarray cannot store: obdeect refuses upward and horizontal photons.
REJECTING_FRAMES = ("ground", "telescope")
#: The camera frame takes photons of any sign (they are at the camera entrance; only cx, cy are used).
ANY_SIGN_FRAME = "camera"
CAMERA_MAPS = {
    "telescope": "1 2 0.0\n",
    "camera": "1 1 0.0\n",
}  # dual mirror, single mirror
COLUMNS = ("x", "y", "z", "cx", "cy", "cz", "t", "zem", "w", "wavelength")
SIGNS = ("upward", "downward", "horizontal")
SIGN_COLORS = {
    "upward": plotting.CYCLE[3],
    "downward": plotting.CYCLE[2],
    "horizontal": plotting.CYCLE[5],
}
REJECTION = (
    "{frame}-frame photon must move downwards (cz < 0): "
    "sim_telarray does not store the sign of cz"
)


def directions() -> dict[str, np.ndarray]:
    """Unit vectors over the full sphere: a Fibonacci grid plus exact vertical and horizontal ones.

    Index UP is (0, 0, 1), DOWN is (0, 0, -1), and the HORIZONTAL indices from FIRST_HORIZONTAL on
    have cz = 0 exactly, in equal steps of azimuth. The rest of the grid has cz != 0.
    """
    index = np.arange(GRID_SIZE)
    cz = 1.0 - (2.0 * index + 1.0) / GRID_SIZE
    azimuth = index * GOLDEN_ANGLE
    radius = np.sqrt(1.0 - cz * cz)
    cx, cy = radius * np.cos(azimuth), radius * np.sin(azimuth)
    cx[UP], cy[UP], cz[UP] = 0.0, 0.0, 1.0
    cx[DOWN], cy[DOWN], cz[DOWN] = 0.0, 0.0, -1.0
    for k in range(HORIZONTAL):
        angle = 2.0 * np.pi * k / HORIZONTAL
        i = FIRST_HORIZONTAL + k
        cx[i], cy[i], cz[i] = np.cos(angle), np.sin(angle), 0.0
    return {
        "cx": cx.astype(np.float32),
        "cy": cy.astype(np.float32),
        "cz": cz.astype(np.float32),
    }


def sign_counts(cz) -> dict[str, int]:
    cz = np.asarray(cz)
    return {
        "upward": int(np.count_nonzero(cz > 0)),
        "downward": int(np.count_nonzero(cz < 0)),
        "horizontal": int(np.count_nonzero(cz == 0)),
    }


def photon_columns(directions_, seed: int) -> dict[str, np.ndarray]:
    """All stream columns for the directions; positions, times and energies are random float32."""
    rng = np.random.default_rng(seed)
    rows = len(directions_["cz"])

    def uniform(low: float, high: float) -> np.ndarray:
        return rng.uniform(low, high, rows).astype(np.float32)

    return {
        "x": uniform(-100.0, 100.0),
        "y": uniform(-100.0, 100.0),
        "z": uniform(-100.0, 100.0),
        "cx": directions_["cx"],
        "cy": directions_["cy"],
        "cz": directions_["cz"],
        "t": uniform(0.0, 50.0),
        "zem": uniform(2000.0, 20000.0),
        "w": np.ones(rows, dtype=np.float32),
        "wavelength": uniform(300.0, 600.0),
    }


def _run(telescopes: int) -> dict:
    return {
        "run_number": 1,
        "observation_level_m": 2200.0,
        "telescope_x_m": np.zeros(telescopes),
        "telescope_y_m": np.zeros(telescopes),
        "telescope_z_m": np.zeros(telescopes),
        "telescope_r_m": np.full(telescopes, 1.2),
    }


def write_round_trip(path: Path, directions_: dict) -> dict[int, tuple[str, dict]]:
    """One telescope per frame (ids 1..4); returns {telescope_id: (frame, columns)} as written."""
    written = {}
    with Writer(str(path), "obdeect validation: direction round trip") as writer:
        writer.write("run", _run(len(ROUND_TRIP_FRAMES)))
        writer.write("event", {"event_number": 1})
        for telescope_id, frame in enumerate(ROUND_TRIP_FRAMES, start=1):
            columns = photon_columns(directions_, SEED + telescope_id)
            writer.write(
                "telescope",
                {"telescope_id": telescope_id, "frame": frame, "array_index": 0},
            )
            writer.write_photons(telescope_id, frame, 0, columns)
            writer.write(
                "telescope_end",
                {"telescope_id": telescope_id, "photons": float(len(columns["cz"]))},
            )
            written[telescope_id] = (frame, columns)
        writer.write("event_end", {})
        writer.write("run_end", {})
    return written


def read_back(path: Path) -> tuple[dict, dict[int, list[dict]]]:
    """Stream metadata and the photon chunks of each telescope, as read."""
    header, chunks = {}, {}
    for record in Reader(str(path)):
        if record["type"] == "stream":
            header = dict(record["metadata"])
        elif record["type"] == "photons":
            chunks.setdefault(record["telescope_id"], []).append(record)
    return header, chunks


def compare(written: dict, chunks: dict) -> dict:
    """Bit mismatches, largest absolute difference and photon counts of the read-back chunks."""
    mismatches, largest, photons = 0, 0.0, 0
    per_frame = {}
    for telescope_id, (frame, expected) in written.items():
        parts = chunks.get(telescope_id, [])
        if any(part["frame"] != FRAMES[frame] for part in parts):
            raise ValueError(f"telescope {telescope_id}: frame differs from {frame}")
        actual = {
            name: np.concatenate([part["columns"][name] for part in parts])
            if all(name in part["columns"] for part in parts) and parts
            else np.empty(0, dtype=np.float32)
            for name in COLUMNS
        }
        frame_mismatch, frame_largest = 0, 0.0
        for name in COLUMNS:
            want = np.asarray(expected[name], dtype=np.float32)
            got = actual[name].astype(np.float32, copy=False)
            if got.shape != want.shape:  # a missing or short column: every row counts
                frame_mismatch += want.size
                continue
            differs = np.count_nonzero(got.view(np.uint32) != want.view(np.uint32))
            frame_mismatch += int(differs)
            if want.size:
                frame_largest = max(
                    frame_largest,
                    float(
                        np.max(np.abs(got.astype(np.float64) - want.astype(np.float64)))
                    ),
                )
        frame_photons = int(actual["x"].size)
        photons += frame_photons
        mismatches += frame_mismatch
        largest = max(largest, frame_largest)
        per_frame[frame] = {
            "photons": frame_photons,
            "bit_mismatches": frame_mismatch,
            "max_abs_difference": frame_largest,
        }
    return {
        "photon_count": photons,
        "bit_mismatch_count": mismatches,
        "max_abs_difference": largest,
        "per_frame": per_frame,
    }


def stream_with_photons(path: Path, frame: str, cx, cy, cz) -> Path:
    """Single-telescope stream (id 1) in `frame` with the given directions; positions are fixed."""
    rows = len(cz)

    def column(value):
        return np.full(rows, value, dtype=np.float32)

    columns = {
        "x": column(0.1),
        "y": column(0.2),
        "z": column(0.5),
        "cx": np.asarray(cx, dtype=np.float32),
        "cy": np.asarray(cy, dtype=np.float32),
        "cz": np.asarray(cz, dtype=np.float32),
        "t": column(3.0),
        "zem": column(2200.5),
        "w": column(1.0),
        "wavelength": column(400.0),
    }
    with Writer(str(path), "obdeect validation: to-simtel rejection") as writer:
        writer.write("run", _run(1))
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1, "frame": frame, "array_index": 0})
        writer.write_photons(1, frame, 0, columns)
        writer.write("telescope_end", {"telescope_id": 1, "photons": float(rows)})
        writer.write("event_end", {})
        writer.write("run_end", {})
    return path


def to_simtel(stream: Path, output: Path, camera_map: str | None, log: Path):
    """Run `obdeect stream to-simtel`; returns (converted, log text)."""
    argv = ["obdeect", "stream", "to-simtel"]
    if camera_map is not None:
        map_file = stream.with_suffix(".map")
        map_file.write_text(camera_map, encoding="utf-8")
        argv += ["--camera-map", map_file]
    argv += ["-o", output, stream]
    try:
        common.command(argv, log)
        converted = True
    except RuntimeError:
        converted = False
    return converted, Path(log).read_text(encoding="utf-8", errors="replace")


def rejection_checks(root: Path) -> tuple[dict, dict, dict]:
    """Ground and telescope: an upward plus horizontal stream is refused; camera: any sign converts.

    The downward control (same layout, cz < 0 only) must convert in every SIMTEL frame. The camera
    frame takes a mixed stream (upward, downward and horizontal photons) as it is.
    """
    rejected, accepted, camera_mixed = {}, {}, {}
    for frame in SIMTEL_FRAMES:
        camera_map = CAMERA_MAPS.get(frame)
        if frame == ANY_SIGN_FRAME:
            # Upward, downward and horizontal photons together: must convert.
            signs = stream_with_photons(
                root / f"{frame}-any-sign.pst",
                frame,
                cx=[0.6, 0.6, 0.0],
                cy=[0.0, 0.0, 1.0],
                cz=[0.8, -0.8, 0.0],
            )
            output = root / f"{frame}-any-sign.simtel"
            converted_all, text_all = to_simtel(
                signs, output, camera_map, root / f"{frame}-any-sign.log"
            )
            camera_mixed[frame] = {
                "converted": converted_all
                and output.is_file()
                and output.stat().st_size > 0,
                "log_tail": [line for line in text_all.strip().splitlines()[-2:]],
            }
        else:
            mixed = stream_with_photons(
                root / f"{frame}-mixed.pst",
                frame,
                cx=[0.6, 0.0],
                cy=[0.0, 1.0],
                cz=[0.8, 0.0],
            )
            converted, text = to_simtel(
                mixed,
                root / f"{frame}-mixed.simtel",
                camera_map,
                root / f"{frame}-mixed.log",
            )
            message = REJECTION.format(frame=frame)
            rejected[frame] = {
                "converted": converted,
                "message_found": message in text,
                "log_tail": [line for line in text.strip().splitlines()[-2:]],
            }
        control = stream_with_photons(
            root / f"{frame}-downward.pst",
            frame,
            cx=[0.6, -0.6],
            cy=[0.0, 0.0],
            cz=[-0.8, -0.8],
        )
        output = root / f"{frame}-downward.simtel"
        converted, text = to_simtel(
            control, output, camera_map, root / f"{frame}-downward.log"
        )
        accepted[frame] = {
            "converted": converted and output.is_file() and output.stat().st_size > 0,
            "output_bytes": output.stat().st_size if output.is_file() else 0,
        }
    return rejected, accepted, camera_mixed


def plot_sky(path: Path, chunks: dict) -> None:
    """cz against azimuth of the read-back sky-frame directions (obdeect markers), by sign."""
    readback = chunks[1][0]["columns"]
    cx, cy, cz = (
        np.asarray(readback[name], dtype=np.float64) for name in ("cx", "cy", "cz")
    )
    azimuth = np.degrees(np.arctan2(cy, cx)) % 360.0
    fig, axis = plotting.figure(1, 1, (7.0, 4.2))
    for sign in SIGNS:
        select = {
            "upward": cz > 0,
            "downward": cz < 0,
            "horizontal": cz == 0,
        }[sign]
        axis.plot(
            azimuth[select],
            cz[select],
            "o",
            markersize=3.5,
            markerfacecolor="none",
            color=SIGN_COLORS[sign],
            label=f"{sign} ({int(select.sum())})",
        )
    axis.set(
        xlabel="azimuth of direction [deg]",
        ylabel="cz (direction cosine, up positive)",
        ylim=(-1.05, 1.05),
    )
    axis.legend(
        title="read back from .pst (sky frame)",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
    )
    plotting.save(fig, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args()
    root = Path.cwd()
    grid = directions()
    expected_counts = sign_counts(grid["cz"])

    stream = root / "directions.pst"
    written = write_round_trip(stream, grid)
    header, chunks = read_back(stream)
    comparison = compare(written, chunks)
    readback_counts = {name: 0 for name in SIGNS}
    for parts in chunks.values():
        for part in parts:
            for name, count in sign_counts(part["columns"]["cz"]).items():
                readback_counts[name] += count
    expected_total = {
        name: count * len(written) for name, count in expected_counts.items()
    }

    rejected, accepted, camera_mixed = rejection_checks(root)
    plot_sky(root / "sky-map.png", chunks)

    metrics = {
        "photon_count": comparison["photon_count"],
        "bit_mismatch_count": comparison["bit_mismatch_count"],
        "max_abs_difference": comparison["max_abs_difference"],
        "format_major": int(header.get("format_major", -1)),
        "format_minor": int(header.get("format_minor", -1)),
        # Format 1.1 introduced cz of either sign; later minor versions keep it.
        "format_supports_directions": int(
            (int(header.get("format_major", -1)), int(header.get("format_minor", -1)))
            >= (1, 1)
        ),
        "cz_upward_count_delta": abs(
            readback_counts["upward"] - expected_total["upward"]
        ),
        "cz_downward_count_delta": abs(
            readback_counts["downward"] - expected_total["downward"]
        ),
        "cz_horizontal_count_delta": abs(
            readback_counts["horizontal"] - expected_total["horizontal"]
        ),
        "rejected_cases": sum(
            1
            for case in rejected.values()
            if not case["converted"] and case["message_found"]
        ),
        "camera_any_sign_accepted": sum(
            1 for case in camera_mixed.values() if case["converted"]
        ),
        "accepted_downward": sum(1 for case in accepted.values() if case["converted"]),
    }
    details = {
        "grid_size": GRID_SIZE,
        "exact_directions": {
            "up": [0, 0, 1],
            "down": [0, 0, -1],
            "horizontal_count": HORIZONTAL,
        },
        "seed": SEED,
        "stream_format_constants": {"major": format_major, "minor": format_minor},
        "obdeect_stream_version": importlib.metadata.version("obdeect-stream"),
        "input_sign_counts_per_frame": expected_counts,
        "readback_sign_counts_all_frames": readback_counts,
        "round_trip": comparison,
        "rejection_message": REJECTION,
        "rejection_by_frame": rejected,
        "camera_any_sign_by_frame": camera_mixed,
        "downward_control_by_frame": accepted,
        "command": "obdeect stream to-simtel [--camera-map FILE] -o OUT IN",
    }
    common.write_results(metrics, details, root)


if __name__ == "__main__":
    main()
