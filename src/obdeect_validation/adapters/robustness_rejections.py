"""Rejection of invalid or unsupported photon-stream input (robustness.stream-rejections).

Complements robustness.stream-damage (truncation, corrupt rows and magic, format major through the
pipe chain). This adapter feeds well-formed but invalid or unsupported input to the readers and
writers and requires an explicit error. AGENTS.md: "Reject invalid or unsupported input; never
silently discard fields."

Scenarios. Every gate is a count of violations with expected 0 and atol 0: a count is an exact
integer, so no tolerance applies. Each variant is a stream that must be rejected; each control is a
valid stream that must be accepted (a rejected control means the harness is wrong, not obdeect, and
is counted by control_rejected, also gate 0).

  1. Unknown format major (docs/stream/format.md: readers reject a different format_major). The
     Python reader must raise on the first call and on every later call. Gates
     reader_first_call_accepted and reader_continues_after_rejection (a later call returns a record).
  2. Record nesting (format.md, Nesting): photons before any telescope record, telescope_end without
     a telescope, run_end without a run, event inside an event. Consumers: the Python reader
     (full iteration) and obdeect stream cat. Gate nesting_violations_accepted.
  3. Values the format forbids: NaN in x, y, t or w; a direction that is not a unit vector (cz column
     -0.5 with small cx, cy); an upward cz (cz > 0) in a format 1.0 stream (format.md: 1.0 has cz
     negative; 1.1 allows upward, the control). Consumers: reader and cat. Gate
     invalid_values_accepted. format.md does not spell out finiteness of x, y, t, w; the NaN variants
     follow the brief (the photon chunk validation requires finite path, sphere rows finite).
  4. to-simtel with camera-frame (C) or telescope-frame (T) photons and no wavelength column.
     sim_telarray treats wavelength 0 as "sample a Cherenkov spectrum", so the converter must refuse.
     The control is the same T-frame stream with the wavelength column; it must be accepted, so the
     rejection is due to the missing column. Gate missing_wavelength_accepted.
  5. CORSIKA input with RUNH while a run is open (no RUNE between the runs). The control, the same two
     runs with RUNE closing the first, must be accepted. Gate corsika_run_reopen_accepted.
  6. Write failure: to-simtel, to-grisu and cat writing to /dev/full must exit non-zero. Each command
     first writes to a regular file (control, must exit 0). /dev/full exists only on Linux; elsewhere
     the gate is 0 with details "skipped: no /dev/full". That is a skip, not a pass. Gate
     write_failure_exit_zero.
  7. Efficiency status going backwards or repeated (docs/stream/efficiency-status.md: a factor only
     advances pending -> envelope -> applied; a repeated or backward state is invalid). Variants:
     run scope backward (qe pending, applied, pending), run scope repeated (qe pending, pending),
     telescope scope backward (qe applied, pending after the telescope record). Consumers: reader,
     cat and to-simtel. Control: atmosphere applied, qe pending, accepted by all three. Gate
     efficiency_backward_accepted.

Valid inputs are generated with the obdeect Writer and byte edits of the version fields. The invalid
nesting and efficiency-status streams are concatenations of raw records (the Writer refuses them):
header and stream_end, photon records and the metadata encoding are the same bytes the Writer writes. Photon columns use a fixed seed. The CORSIKA input is
EventIO bytes written here: sync 0xD41F8A37, type/version word, ident, length word, body (see
obdeect/stream/eventio.hpp and corsika.hpp for the block types 1200, 1202, 1209, 1210). Results do not
depend on the host except through /dev/full.
"""

from __future__ import annotations

import struct
import sys
import tempfile
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common
from obdeect_validation.adapters import endtoend_common as e2e

PRODUCER = "robustness-rejections"
TIMEOUT = 120.0
SEED = 20261123
EVENTIO_SYNC = 0xD41F8A37
FRAME_GROUND, FRAME_TELESCOPE, FRAME_CAMERA = 1, 2, 3
CAMERA_MAP_C = "1 1 0.0\n"  # camera frame, single mirror
CAMERA_MAP_T = "1 2 0.0\n"  # telescope frame, dual mirror
PHOTON_ROWS = 5
CAMERA_NAMES = ("x", "y", "cx", "cy", "t", "w", "zem", "wavelength")
TELESCOPE_NAMES = ("x", "y", "z", "cx", "cy", "t", "w", "zem", "wavelength")
GROUND_NAMES = CAMERA_NAMES

RUN = {
    "run_number": 1,
    "observation_level_m": 2150.0,
    "energy_min_gev": 0.01,
    "energy_max_gev": 100.0,
    "spectral_index": -2.7,
    "n_telescopes": 1,
    "telescope_x_m": np.array([0.0]),
    "telescope_y_m": np.array([0.0]),
    "telescope_z_m": np.array([0.0]),
    "telescope_r_m": np.array([5.0]),
    "wavelength_min_nm": 240.0,
    "wavelength_max_nm": 900.0,
}


def event_meta(number=1):
    return {
        "event_number": number,
        "primary_id": 1,
        "energy_gev": 0.1,
        "zenith_rad": 0.3,
        "azimuth_rad": 0.0,
        "first_interaction_m": 10000.0,
        "core_x_m": np.array([0.0]),
        "core_y_m": np.array([0.0]),
        "time_offset_ns": np.array([0.0]),
    }


def telescope_meta(tel=1, frame=FRAME_CAMERA):
    return {
        "telescope_id": tel,
        "frame": frame,
        "x_m": 0.0,
        "y_m": 0.0,
        "z_m": 0.0,
        "r_m": 5.0,
    }


def status(factor, state):
    """efficiency_status metadata updating one factor (each record names only its factors)."""
    return {"factors": factor, f"{factor}.state": state, f"{factor}.source": PRODUCER}


def chunk_columns(names, upward=False, cz_value=None):
    """Photon columns of PHOTON_ROWS rows with a unit direction, downward unless upward.

    cz_value replaces the cz column (for a non-unit direction). Fixed seed.
    """
    rng = np.random.default_rng(SEED)
    cx = rng.uniform(-0.05, 0.05, PHOTON_ROWS)
    cy = rng.uniform(-0.05, 0.05, PHOTON_ROWS)
    sign = 1.0 if upward else -1.0
    cz = sign * np.sqrt(1.0 - cx**2 - cy**2)
    if cz_value is not None:
        cz = np.full(PHOTON_ROWS, cz_value)
    values = {
        "x": rng.uniform(-2.0, 2.0, PHOTON_ROWS),
        "y": rng.uniform(-2.0, 2.0, PHOTON_ROWS),
        "z": np.full(PHOTON_ROWS, 0.5),
        "cx": cx,
        "cy": cy,
        "cz": cz,
        "t": rng.uniform(10.0, 20.0, PHOTON_ROWS),
        "zem": np.full(PHOTON_ROWS, 2160.0),
        "w": np.ones(PHOTON_ROWS),
        "wavelength": rng.uniform(300.0, 600.0, PHOTON_ROWS),
    }
    return {name: np.asarray(values[name], dtype=np.float32) for name in names}


def telescope_ops(columns, tel=1, frame=FRAME_CAMERA):
    return [
        ("telescope", telescope_meta(tel, frame)),
        ("photons", tel, frame, 0, columns),
        ("telescope_end", {"telescope_id": tel, "photons": float(len(columns["x"]))}),
    ]


def event_ops(columns, frame=FRAME_CAMERA, number=1):
    return (
        [("event", event_meta(number))]
        + telescope_ops(columns, frame=frame)
        + [("event_end", {})]
    )


def valid_ops(columns, frame=FRAME_CAMERA):
    return [("run", RUN)] + event_ops(columns, frame=frame) + [("run_end", {})]


def write_stream(path, ops):
    """Write ops with the obdeect Writer; it checks neither nesting nor photon values.

    An op is (record type, metadata) or ("photons", telescope, frame, first_index, columns).
    """
    from obdeect.stream import Writer

    with Writer(str(path), PRODUCER) as writer:
        for op in ops:
            if op[0] == "photons":
                _, tel, frame, first, columns = op
                writer.write_photons(tel, frame, first, columns)
            else:
                writer.write(op[0], op[1])
    return Path(path)


RECORD_TYPE = {
    "stream": 1,
    "run": 2,
    "event": 3,
    "telescope": 4,
    "photons": 5,
    "telescope_end": 6,
    "event_end": 7,
    "run_end": 8,
    "efficiency_status": 10,
}
RECORD_MAGIC = 0x52545347  # "GSTR"


def split_records(data):
    """(type, record bytes) of each record of stream bytes (header and payload)."""
    records, at = [], 0
    while at < len(data):
        size = int.from_bytes(data[at + 8 : at + 16], "little")
        records.append(
            (int.from_bytes(data[at + 4 : at + 6], "little"), data[at : at + 16 + size])
        )
        at += 16 + size
    return records


def encode_metadata(metadata):
    """Metadata payload: u32 count, then per entry u16 key length, key, u8 type and the value.

    Types: 1 int64, 2 float64, 3 string, 4 array of float64.
    """
    out = struct.pack("<I", len(metadata))
    for key, value in metadata.items():
        name = key.encode()
        out += struct.pack("<H", len(name)) + name
        if isinstance(value, str):
            text = value.encode()
            out += b"\x03" + struct.pack("<I", len(text)) + text
        elif isinstance(value, np.ndarray):
            values = np.asarray(value, dtype=float).ravel()
            out += (
                b"\x04"
                + struct.pack("<I", len(values))
                + struct.pack(f"<{len(values)}d", *values)
            )
        elif isinstance(value, float):
            out += b"\x02" + struct.pack("<d", value)
        else:
            out += b"\x01" + struct.pack("<q", int(value))
    return out


def raw_record(type_code, payload):
    """One record as raw bytes: magic, type, reserved flags, payload length, payload."""
    return struct.pack("<IHHQ", RECORD_MAGIC, type_code, 0, len(payload)) + payload


def photons_record(columns, tel, frame, first):
    """The bytes of one valid photons record, taken from a stream the Writer writes."""
    ops = [
        ("run", RUN),
        ("event", event_meta()),
        ("telescope", telescope_meta(tel, frame)),
        ("photons", tel, frame, first, columns),
        ("telescope_end", {"telescope_id": tel, "photons": float(len(columns["x"]))}),
        ("event_end", {}),
        ("run_end", {}),
    ]
    with tempfile.TemporaryDirectory() as directory:
        source = write_stream(Path(directory) / "photons.pst", ops)
        (record,) = [
            raw
            for kind, raw in split_records(source.read_bytes())
            if kind == RECORD_TYPE["photons"]
        ]
    return record


def write_raw_stream(path, ops):
    """Write ops as raw records: the stream record and stream_end of a valid Writer stream around
    them. An op is (record type, metadata) or ("photons", telescope, frame, first_index, columns).
    Nothing is checked: the streams are invalid on purpose."""
    with tempfile.TemporaryDirectory() as directory:
        source = write_stream(
            Path(directory) / "valid.pst", valid_ops(chunk_columns(CAMERA_NAMES))
        )
        records = split_records(source.read_bytes())
    header, end = records[0][1], records[-1][1]
    body = b""
    for op in ops:
        if op[0] == "photons":
            _, tel, frame, first, columns = op
            body += photons_record(columns, tel, frame, first)
        else:
            body += raw_record(RECORD_TYPE[op[0]], encode_metadata(op[1]))
    Path(path).write_bytes(header + body + end)
    return Path(path)


def set_int_field(path, key, value):
    """Overwrite an integer metadata value of the stream record (the width does not change).

    Entries are key length (u16), key, type byte (1 = int), i64 value; the type byte follows the key.
    """
    data = bytearray(Path(path).read_bytes())
    name = key.encode()
    at = data.index(name) + len(name)
    if data[at] != 1:
        raise ValueError(f"metadata '{key}' is not an integer entry")
    data[at + 1 : at + 9] = struct.pack("<q", value)
    Path(path).write_bytes(bytes(data))
    return Path(path)


def eventio_block(block_type, ident, body, *, top_level=True):
    """EventIO block: optional sync tag, type/version word (version 0), ident, length word, body."""
    body = bytes(body) + b"\0" * (-len(body) % 4)
    header = struct.pack("<II", block_type & 0xFFFF, ident) + struct.pack(
        "<I", len(body)
    )
    if top_level:
        header = struct.pack("<I", EVENTIO_SYNC) + header
    return header + body


def float_block(values):
    """Float array block body (RUNH, EVTH, EVTE, RUNE): int32 count, then float32 values."""
    return struct.pack("<i", len(values)) + struct.pack(f"<{len(values)}f", *values)


def corsika_run(run_number, *, events=1, run_end=True):
    """Blocks of one CORSIKA run: RUNH, per event EVTH and EVTE, RUNE when run_end."""
    runh = [0.0] * 273
    runh[1] = float(run_number)
    runh[5] = 215000.0  # observation level [cm]
    runh[15] = -2.7  # spectral slope
    runh[16], runh[17] = 0.01, 100.0  # energy range [GeV]
    blocks = [eventio_block(1200, 0, float_block(runh))]
    for number in range(1, events + 1):
        evth = [0.0] * 273
        evth[1] = float(number)
        evth[2] = 1.0  # gamma
        evth[3] = 0.1  # energy [GeV]
        evth[6] = 1.0e6  # first interaction [cm]
        evth[10] = 0.3  # zenith [rad]
        evth[95], evth[96] = 240.0, 900.0  # wavelength range [nm]
        blocks.append(eventio_block(1202, 0, float_block(evth)))
        blocks.append(eventio_block(1209, 0, float_block([0.0])))
    if run_end:
        blocks.append(eventio_block(1210, 0, float_block([0.0])))
    return blocks


# ---------------------------------------------------------------------------
# Consumers. Each returns (accepted, message): accepted means exit status 0 and no hang.


def accepted(stage):
    return stage.returncode == 0 and not stage.timed_out


def last_line(stage):
    lines = (stage.stderr or "").strip().splitlines()
    return lines[-1] if lines else f"exit status {stage.returncode}"


def run_command(argv, workdir, name, *, output=None):
    return e2e.run_stage(
        argv, log=workdir / f"{name}.log", output=output, timeout=TIMEOUT
    )


def reader_accepts(stream, workdir, name):
    """The Python reader iterates the whole stream without an exception (own process)."""
    code = "import sys\nfrom obdeect.stream import Reader\nprint(sum(1 for _ in Reader(sys.argv[1])))\n"
    stage = run_command(
        [sys.executable, "-c", code, str(stream)], workdir, f"{name}.reader"
    )
    return accepted(stage), last_line(stage)


def cat_accepts(stream, workdir, name):
    stage = run_command(
        ["obdeect", "stream", "cat", str(stream)],
        workdir,
        f"{name}.cat",
        output=workdir / f"{name}.cat.out",
    )
    return accepted(stage), last_line(stage)


def to_simtel_accepts(stream, camera_map, workdir, name):
    map_path = workdir / f"{name}.camera-map.txt"
    map_path.write_text(camera_map)
    stage = run_command(
        [
            "obdeect",
            "stream",
            "to-simtel",
            "--camera-map",
            str(map_path),
            "-o",
            str(workdir / f"{name}.iact"),
            str(stream),
        ],
        workdir,
        f"{name}.to-simtel",
    )
    return accepted(stage), last_line(stage)


# ---------------------------------------------------------------------------
# Scenarios. Each appends result rows: {scenario, name, consumer, role, accepted, message}.


def row(results, scenario, name, consumer, role, is_accepted, message=""):
    results.append(
        {
            "scenario": scenario,
            "name": name,
            "consumer": consumer,
            "role": role,
            "accepted": bool(is_accepted),
            "message": message,
        }
    )


def scenario_unknown_major(workdir, results):
    stream = write_stream(
        workdir / "major-2.pst", valid_ops(chunk_columns(CAMERA_NAMES))
    )
    set_int_field(stream, "format_major", 2)
    code = (
        "import sys\n"
        "from obdeect.stream import Reader\n"
        "r = Reader(sys.argv[1])\n"
        "out = []\n"
        "for _ in range(5):\n"
        "    try:\n"
        "        next(r)\n"
        "        out.append('record')\n"
        "    except StopIteration:\n"
        "        out.append('stop')\n"
        "    except Exception:\n"
        "        out.append('raise')\n"
        "print(' '.join(out))\n"
    )
    output = workdir / "major-2.calls.txt"
    stage = run_command(
        [sys.executable, "-c", code, str(stream)],
        workdir,
        "major-2.calls",
        output=output,
    )
    calls = output.read_text().split() if output.exists() else []
    first_rejected = bool(calls) and calls[0] == "raise"
    later_records = sum(1 for word in calls[1:] if word == "record")
    row(
        results,
        "unknown_major",
        "major-2",
        "reader",
        "variant",
        not first_rejected,
        f"calls: {' '.join(calls)}; {last_line(stage)}",
    )
    row(
        results,
        "unknown_major",
        "major-2-later-calls",
        "reader",
        "variant",
        later_records > 0,
        f"calls: {' '.join(calls)}",
    )
    return calls


def scenario_nesting(workdir, results):
    cols = chunk_columns(CAMERA_NAMES)
    variants = {
        "photons_before_telescope": [
            ("run", RUN),
            ("event", event_meta()),
            ("photons", 1, FRAME_CAMERA, 0, cols),
            ("telescope", telescope_meta()),
            ("telescope_end", {"telescope_id": 1, "photons": float(PHOTON_ROWS)}),
            ("event_end", {}),
            ("run_end", {}),
        ],
        "telescope_end_without_telescope": [
            ("run", RUN),
            ("event", event_meta()),
            ("telescope_end", {"telescope_id": 1, "photons": 0.0}),
            ("event_end", {}),
            ("run_end", {}),
        ],
        "run_end_without_run": valid_ops(cols) + [("run_end", {})],
        "event_inside_event": [
            ("run", RUN),
            ("event", event_meta(1)),
            ("event", event_meta(2)),
        ]
        + telescope_ops(cols)
        + [("event_end", {}), ("event_end", {}), ("run_end", {})],
    }
    for name, ops in variants.items():
        stream = write_raw_stream(workdir / f"nesting-{name}.pst", ops)
        ok, message = reader_accepts(stream, workdir, f"nesting-{name}")
        row(results, "nesting", name, "reader", "variant", ok, message)
        ok, message = cat_accepts(stream, workdir, f"nesting-{name}")
        row(results, "nesting", name, "cat", "variant", ok, message)


def scenario_invalid_values(workdir, results):
    cases = {}
    for column in ("x", "y", "t", "w"):
        cols = chunk_columns(CAMERA_NAMES)
        cols[column][2] = np.float32("nan")
        cases[f"nan-{column}"] = (valid_ops(cols), None)
    cz_names = CAMERA_NAMES + ("cz",)
    cases["nonunit-direction"] = (
        valid_ops(chunk_columns(cz_names, cz_value=-0.5)),
        None,
    )
    cases["upward-format-1-0"] = (valid_ops(chunk_columns(cz_names, upward=True)), 0)
    for name, (ops, minor) in cases.items():
        stream = write_stream(workdir / f"value-{name}.pst", ops)
        if minor is not None:
            set_int_field(stream, "format_minor", minor)
        ok, message = reader_accepts(stream, workdir, f"value-{name}")
        row(results, "invalid_values", name, "reader", "variant", ok, message)
        ok, message = cat_accepts(stream, workdir, f"value-{name}")
        row(results, "invalid_values", name, "cat", "variant", ok, message)
    # Control: the same upward photons are valid in format 1.1.
    control = write_stream(
        workdir / "control-upward-format-1-1.pst",
        valid_ops(chunk_columns(cz_names, upward=True)),
    )
    set_int_field(control, "format_minor", 1)
    for consumer, check in (("reader", reader_accepts), ("cat", cat_accepts)):
        ok, message = check(control, workdir, "control-upward-format-1-1")
        row(
            results,
            "invalid_values",
            "control-upward-format-1-1",
            consumer,
            "control",
            ok,
            message,
        )


def scenario_missing_wavelength(workdir, results):
    no_wl_c = valid_ops(chunk_columns(CAMERA_NAMES[:-1]))
    no_wl_t = valid_ops(chunk_columns(TELESCOPE_NAMES[:-1]), frame=FRAME_TELESCOPE)
    with_wl_t = valid_ops(chunk_columns(TELESCOPE_NAMES), frame=FRAME_TELESCOPE)
    cases = {
        "camera-no-wavelength": (no_wl_c, CAMERA_MAP_C, "variant"),
        "telescope-no-wavelength": (no_wl_t, CAMERA_MAP_T, "variant"),
        "control-telescope-with-wavelength": (with_wl_t, CAMERA_MAP_T, "control"),
    }
    for name, (ops, camera_map, role) in cases.items():
        stream = write_stream(workdir / f"wl-{name}.pst", ops)
        ok, message = to_simtel_accepts(stream, camera_map, workdir, f"wl-{name}")
        row(results, "missing_wavelength", name, "to-simtel", role, ok, message)


def scenario_corsika_reopen(workdir, results):
    reopened = corsika_run(1, run_end=False) + corsika_run(2)
    control = corsika_run(1) + corsika_run(2)
    for name, blocks, role in (
        ("runh-while-run-open", reopened, "variant"),
        ("control-run-closed-by-rune", control, "control"),
    ):
        path = workdir / f"{name}.iact"
        path.write_bytes(b"".join(blocks))
        stage = run_command(
            ["obdeect", "stream", "from-corsika", str(path)],
            workdir,
            f"{name}.from-corsika",
            output=workdir / f"{name}.pst",
        )
        row(
            results,
            "corsika_reopen",
            name,
            "from-corsika",
            role,
            accepted(stage),
            last_line(stage),
        )


def grisu_stream_ops(cols):
    """Ground-frame stream with atmosphere applied, as obdeect stream to-grisu requires."""
    return [
        ("run", RUN),
        ("efficiency_status", status("atmosphere", "applied")),
        ("event", event_meta()),
        ("telescope", telescope_meta(1, FRAME_GROUND)),
        ("photons", 1, FRAME_GROUND, 0, cols),
        ("telescope_end", {"telescope_id": 1, "photons": float(PHOTON_ROWS)}),
        ("event_end", {}),
        ("run_end", {}),
    ]


def scenario_write_failure(workdir, results):
    """Gate 6. Returns None when /dev/full does not exist (skip)."""
    dev_full = Path("/dev/full")
    if not dev_full.exists():
        return
    camera = write_stream(
        workdir / "wf-camera.pst", valid_ops(chunk_columns(CAMERA_NAMES))
    )
    ground = write_stream(
        workdir / "wf-ground.pst", grisu_stream_ops(chunk_columns(GROUND_NAMES))
    )
    camera_map = workdir / "wf-camera-map.txt"
    camera_map.write_text(CAMERA_MAP_C)
    commands = {
        "to-simtel": lambda out: [
            "obdeect",
            "stream",
            "to-simtel",
            "--camera-map",
            str(camera_map),
            "-o",
            out,
            str(camera),
        ],
        "to-grisu": lambda out: [
            "obdeect",
            "stream",
            "to-grisu",
            "-o",
            out,
            str(ground),
        ],
        "cat": lambda out: ["obdeect", "stream", "cat", str(camera)],
    }
    for name, argv_for in commands.items():
        control = run_command(
            argv_for(str(workdir / f"wf-{name}.out")), workdir, f"wf-{name}.control"
        )
        row(
            results,
            "write_failure",
            name,
            "control-file",
            "control",
            accepted(control),
            last_line(control),
        )
        stage = run_command(argv_for("/dev/full"), workdir, f"wf-{name}.dev-full")
        row(
            results,
            "write_failure",
            name,
            "/dev/full",
            "variant",
            accepted(stage),
            last_line(stage),
        )


def scenario_efficiency(workdir, results):
    cols = chunk_columns(CAMERA_NAMES)
    atmosphere = ("efficiency_status", status("atmosphere", "applied"))
    qe_pending = ("efficiency_status", status("qe", "pending"))
    run_head = [("run", RUN), atmosphere, qe_pending]
    tail = event_ops(cols) + [("run_end", {})]
    telescope_scope = [
        ("telescope", telescope_meta()),
        ("efficiency_status", status("qe", "applied")),
        ("efficiency_status", status("qe", "pending")),
        ("photons", 1, FRAME_CAMERA, 0, cols),
        ("telescope_end", {"telescope_id": 1, "photons": float(PHOTON_ROWS)}),
    ]
    variants = {
        "run-backward": [
            ("run", RUN),
            atmosphere,
            qe_pending,
            ("efficiency_status", status("qe", "applied")),
            ("efficiency_status", status("qe", "pending")),
        ]
        + tail,
        "run-repeated": [
            ("run", RUN),
            atmosphere,
            qe_pending,
            ("efficiency_status", status("qe", "pending")),
        ]
        + tail,
        "telescope-backward": run_head
        + [("event", event_meta())]
        + telescope_scope
        + [("event_end", {}), ("run_end", {})],
        "control-pending-only": run_head + tail,
    }
    for name, ops in variants.items():
        role = "control" if name.startswith("control") else "variant"
        stream = write_raw_stream(workdir / f"eff-{name}.pst", ops)
        ok, message = reader_accepts(stream, workdir, f"eff-{name}")
        row(results, "efficiency", name, "reader", role, ok, message)
        ok, message = cat_accepts(stream, workdir, f"eff-{name}")
        row(results, "efficiency", name, "cat", role, ok, message)
        ok, message = to_simtel_accepts(stream, CAMERA_MAP_C, workdir, f"eff-{name}")
        row(results, "efficiency", name, "to-simtel", role, ok, message)


#: Metric of each scenario that counts accepted variants.
SCENARIO_METRIC = {
    "nesting": "nesting_violations_accepted",
    "invalid_values": "invalid_values_accepted",
    "missing_wavelength": "missing_wavelength_accepted",
    "corsika_reopen": "corsika_run_reopen_accepted",
    "write_failure": "write_failure_exit_zero",
    "efficiency": "efficiency_backward_accepted",
}


def summarize(results):
    """Metrics from the result rows: accepted variants per scenario, rejected controls, unknown major."""
    metrics = {name: 0 for name in SCENARIO_METRIC.values()}
    metrics.update(
        reader_first_call_accepted=0,
        reader_continues_after_rejection=0,
        control_rejected=0,
    )
    for r in results:
        if r["role"] == "control":
            metrics["control_rejected"] += int(not r["accepted"])
        elif (
            r["role"] == "variant"
            and r["accepted"]
            and r["scenario"] in SCENARIO_METRIC
        ):
            metrics[SCENARIO_METRIC[r["scenario"]]] += 1
        if r["scenario"] == "unknown_major" and r["accepted"]:
            key = (
                "reader_first_call_accepted"
                if r["name"] == "major-2"
                else "reader_continues_after_rejection"
            )
            metrics[key] = 1
    return metrics


def main():
    workdir = Path.cwd()
    results: list[dict] = []
    details: dict = {"skipped": []}

    scenario_unknown_major(workdir, results)
    scenario_nesting(workdir, results)
    scenario_invalid_values(workdir, results)
    scenario_missing_wavelength(workdir, results)
    scenario_corsika_reopen(workdir, results)
    if scenario_write_failure(workdir, results) is None:
        details["skipped"].append("write_failure: skipped: no /dev/full")
    scenario_efficiency(workdir, results)

    metrics = summarize(results)
    if details["skipped"]:
        metrics["write_failure_exit_zero"] = 0
    details["accepted_variants"] = [
        {"scenario": r["scenario"], "name": r["name"], "consumer": r["consumer"]}
        for r in results
        if r["role"] == "variant" and r["accepted"]
    ]
    details["results"] = results
    common.write_results(metrics, details, workdir)


if __name__ == "__main__":
    main()
