"""obdeect photon streams in the adapters: read, write, and the float32 rounding unit.

Readers (all read the whole stream; validation inputs are small):
  read_stream(path)   run and run_end metadata, frame codes, telescope id per row, photon columns
  read_columns(path)  photon columns only
  endtoend_common.photon_rows(path) keeps event, telescope and row index per photon (exact comparison)
  readers.read_native(path) per-event weights and float64 columns (source comparisons)

Writer: write_stream(path, producer, run, blocks, ...) writes one run with one event; blocks are
(telescope_id, columns) pairs. Columns are written as given (the stream format takes float32).
The telescope record frame (frame) and the photon frame (photon_frame, e.g. "sky") may differ.

ulp32(values) is the float32 spacing at the magnitude of each value, floored at 1 (a metre, a ns):
the rounding unit of a float32 position or time.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Stream:
    """Content of a stream (read_stream)."""

    run: dict = field(default_factory=dict)
    run_end: dict = field(default_factory=dict)
    frames: set = field(default_factory=set)
    telescope: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    columns: dict = field(default_factory=dict)


def read_stream(path) -> Stream:
    """Run, run_end, frame codes, telescope id per row and concatenated photon columns of path."""
    from obdeect.stream import Reader

    content = Stream()
    parts: dict[str, list] = {}
    telescopes: list = []
    for record in Reader(str(path)):
        kind = record["type"]
        if kind == "run":
            content.run = dict(record["metadata"])
        elif kind == "run_end":
            content.run_end = dict(record["metadata"])
        elif kind == "photons":
            columns = record["columns"]
            if not columns:
                continue
            content.frames.add(int(record["frame"]))
            rows = len(next(iter(columns.values())))
            telescopes.append(np.full(rows, record["telescope_id"], dtype=np.int64))
            for name, values in columns.items():
                parts.setdefault(name, []).append(np.asarray(values))
    content.columns = {name: np.concatenate(v) for name, v in parts.items()}
    if telescopes:
        content.telescope = np.concatenate(telescopes)
    return content


def read_columns(path) -> dict:
    """Concatenated photon columns of path ({} for a stream without photons)."""
    return read_stream(path).columns


def write_stream(
    path,
    producer: str,
    run: dict,
    blocks,
    *,
    frame=0,
    photon_frame=None,
    telescope_extra: dict | None = None,
    event: dict | None = None,
    chunk_rows: int | None = None,
    telescope_end: bool = True,
):
    """Write one run with one event and return path.

    blocks: (telescope_id, columns) pairs. columns is a dict of arrays (written in chunks of
    chunk_rows, default one chunk) or an iterable of such dicts (one chunk each, e.g. a generator
    that draws the rows chunk by chunk). Per block: a telescope record {"telescope_id", "frame":
    frame, **telescope_extra}, the photon chunks (frame photon_frame, default frame; first_index =
    row offset) and (telescope_end) {"telescope_id", "photons": float(rows)}; then event_end and
    run_end, both empty.
    """
    from obdeect.stream import Writer

    photon_frame = frame if photon_frame is None else photon_frame
    with Writer(str(path), producer) as writer:
        writer.write("run", run)
        writer.write("event", event or {"event_number": 1})
        for telescope_id, columns in blocks:
            writer.write(
                "telescope",
                {
                    "telescope_id": telescope_id,
                    "frame": frame,
                    **(telescope_extra or {}),
                },
            )
            rows = 0
            for piece in _chunks(columns, chunk_rows):
                writer.write_photons(telescope_id, photon_frame, rows, piece)
                rows += len(next(iter(piece.values())))
            if telescope_end:
                writer.write(
                    "telescope_end",
                    {"telescope_id": telescope_id, "photons": float(rows)},
                )
        writer.write("event_end", {})
        writer.write("run_end", {})
    return path


def _chunks(columns, chunk_rows):
    if not isinstance(columns, dict):
        yield from columns
        return
    rows = len(next(iter(columns.values())))
    step = chunk_rows or rows
    for first in range(0, rows, max(step, 1)):
        yield {k: np.asarray(v)[first : first + step] for k, v in columns.items()}


def ulp32(values, floor: float = 1.0) -> np.ndarray:
    """float32 spacing at max(|values|, floor), as float64 (see the module docstring)."""
    magnitude = np.maximum(np.abs(np.asarray(values, dtype=np.float64)), floor)
    return np.spacing(magnitude.astype(np.float32)).astype(np.float64)
