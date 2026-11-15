"""Photon readers shared by the source adapters: obdeect streams and EventIO IACT files.

Both return (weight per event, columns). Columns are float64 arrays; positions in
metres, times in ns, wavelengths in nm, w the photon weight (bunch size).
"""

from __future__ import annotations

import numpy as np
from eventio.iact import IACTFile

NATIVE_KEYS = ("x", "y", "z", "cx", "cy", "cz", "t", "zem", "w", "wavelength")
# IACT field -> column; positions are stored in cm.
_IACT_FIELDS = (
    ("x", "x"),
    ("y", "y"),
    ("zem", "zem"),
    ("cx", "cx"),
    ("cy", "cy"),
    ("t", "time"),
    ("wavelength", "wavelength"),
    ("w", "photons"),
)
_CM_TO_M = {"x": 0.01, "y": 0.01, "zem": 0.01}


def _concatenate(chunks):
    return {k: np.concatenate(v) if v else np.empty(0) for k, v in chunks.items()}


def downward_cz(cx, cy):
    """cz of unit vectors pointing down (cz < 0) from cx and cy."""
    return -np.sqrt(np.clip(1.0 - cx**2 - cy**2, 0.0, None))


def read_native(path, keys=NATIVE_KEYS, require_events=True):
    """obdeect photon stream. Missing columns stay empty, except cz (from cx, cy)
    and z (zero for camera-frame streams).

    Photons before the first event record are an error unless require_events is False
    (then they count as one event; the per-event weights are not meaningful).
    """
    from obdeect.stream import Reader  # optional dependency, only for native streams

    counts, chunks = [], {k: [] for k in keys}
    for record in Reader(str(path)):
        if record["type"] == "event":
            counts.append(0.0)
        elif record["type"] == "photons":
            columns = dict(record["columns"])
            if "cz" in keys and "cz" not in columns and "cx" in columns:
                columns["cz"] = downward_cz(
                    np.asarray(columns["cx"], dtype=np.float64),
                    np.asarray(columns["cy"], dtype=np.float64),
                )
            weights = np.asarray(columns["w"], dtype=np.float64)
            if not counts:
                if require_events:
                    raise ValueError("photon stream has photons before its first event")
                counts.append(0.0)
            counts[-1] += float(np.sum(weights))
            for key in keys:
                if key in columns:
                    chunks[key].append(np.asarray(columns[key], dtype=np.float64))
                elif key == "z":
                    chunks[key].append(np.zeros(len(weights)))
    return counts, _concatenate(chunks)


def read_iact(path):
    """EventIO IACT photon bunches; t includes the event time offset, cz is derived."""
    counts, chunks = [], {k: [] for k, _ in _IACT_FIELDS}
    with IACTFile(str(path)) as stream:
        for event in stream:
            total = 0.0
            for bunches in event.photon_bunches.values():
                names = bunches.dtype.names or ()
                total += float(np.sum(bunches["photons"], dtype=np.float64))
                for key, field in _IACT_FIELDS:
                    if field not in names:
                        continue
                    values = np.asarray(bunches[field], dtype=np.float64)
                    values = values * _CM_TO_M.get(key, 1.0)
                    if key == "t":
                        values = values + float(event.time_offset)
                    chunks[key].append(values)
            counts.append(total)
    data = _concatenate(chunks)
    if data["cx"].size:
        data["cz"] = downward_cz(data["cx"], data["cy"])
    return counts, data
