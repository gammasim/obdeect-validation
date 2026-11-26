"""CORSIKA 3D photon bunches (1217) to sphere rows: from-corsika and to-simtel (plane conversion).

Case stream.corsika-photons3d. The input is the fixture fixtures/sphere/photons3d.iact (CORSIKA IACT
EventIO, one run, two telescopes, eleven bunches; provenance in photons3d_provenance.txt). The
reference is independent of obdeect: the EventIO structure is walked here from the bytes (documented
item header, containers without sync tags; pyeventio is used only to cross-check the top-level block
types and the telescope positions), the 1217 body is decoded here (i16 array, i16 telescope, f32
photons, i32 count, ten f32 per bunch, CORSIKA units cm, ns, nm), and the sphere rows are computed in
float64 from the documented rules (docs/stream/command-line.md, "CORSIKA 3D photon bunches (1217)"):

- E = A - dist d (emission point, A the arrival point relative to the sphere centre);
- E on or inside the sphere: row E, path 0, t = ctime - dist / c;
- otherwise the smallest root s >= 0 of |E + s d| = R gives the row E + s d, path s, t = ctime + (s -
  dist) / c; no such root: misses_sphere;
- zem = observation level + telescope z + E_z (a.s.l., m).

Rows are compared with the float32 values of `obdeect stream from-corsika` in the order of the bunches
(misses removed), column by column, in float32 units of the stored value (ulp, floored at 1). The
run key ground_representation is sphere, run_end misses_sphere equals the reference count.

to-simtel: the full fixture is refused (a photon emitted inside the sphere); a subset with the
downward on-sphere rows only converts, and each row moves to the plane z = 0 along its line (x, y at
the crossing, t + s/c with s = -z/cz, zem and cz unchanged), read back with pyeventio's IACT 1205
reader (long format: x, y, zem in cm, time in ns).
"""

from __future__ import annotations

import math
import struct
import subprocess
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common
from obdeect_validation.adapters import upward_common as uc

FIXTURE = uc.ROOT / "fixtures/sphere/photons3d.iact"
PROVENANCE = uc.ROOT / "fixtures/sphere/photons3d_provenance.txt"
SYNC = bytes.fromhex("378a1fd4")
C_M_NS = uc.C_M_NS
COLUMNS = ("x", "y", "z", "cx", "cy", "cz", "t", "zem", "w", "path", "wavelength")
CM = 0.01


def version() -> str:
    return common.package_version("obdeect-stream")


# ----------------------------------------------------------------- EventIO, independent reader


def _items(buf: bytes, start: int, end: int, top: bool) -> list[dict]:
    """Items of a byte range: top-level items carry the sync tag, container bodies do not.

    Item header (12 bytes; 16 with the length extension, flagged in bit 17 of the type word and bit 31
    of the length word): type/version word (type bits 0-15, version bits 20-31), ident, length word
    (length bits 0-29, bit 30 = only sub-items). Containers hold their children as the body.
    """
    out = []
    pos = start
    while pos < end:
        if top:
            if buf[pos : pos + 4] != SYNC:
                raise ValueError(f"no sync tag at byte {pos}")
            pos += 4
        type_version, ident, length_word = struct.unpack_from("<III", buf, pos)
        extended = bool(length_word & 0x80000000)
        if extended != bool(type_version & 0x20000):
            raise ValueError(f"inconsistent extension flags at byte {pos}")
        header = 16 if extended else 12
        length = length_word & 0x3FFFFFFF
        if extended:
            length |= struct.unpack_from("<I", buf, pos + 12)[0] << 30
        body = pos + header
        node = {
            "type": type_version & 0xFFFF,
            "version": (type_version >> 20) & 0xFFF,
            "ident": ident,
            "body": (body, body + length),
            "only_subobjects": bool(length_word & 0x40000000),
        }
        if node["only_subobjects"]:
            node["children"] = _items(buf, body, body + length, top=False)
        out.append(node)
        pos = body + length
    if pos != end:
        raise ValueError("item sizes do not fill the container")
    return out


def read_eventio(path=FIXTURE) -> list[dict]:
    """Top-level items of the file, with children for containers."""
    data = Path(path).read_bytes()
    return _items(data, 0, len(data), top=True)


def _block(items: list[dict], block_type: int) -> dict:
    found = [node for node in items if node["type"] == block_type]
    if len(found) != 1:
        raise ValueError(f"expected one block of type {block_type}, found {len(found)}")
    return found[0]


def read_run(path=FIXTURE) -> dict:
    """Observation level [m] (RUNH, runh[5] in cm), telescopes (TELPOS: int count, then x, y, z, r as
    float32 columns, cm) and the 1217 bunches per telescope (array, telescope, photons, bunches)."""
    data = Path(path).read_bytes()
    items = read_eventio(path)
    runh = _block(items, 1200)
    start, _ = runh["body"]
    observation_level_cm = struct.unpack_from("<f", data, start + 4 + 5 * 4)[0]
    telpos = _block(items, 1201)
    start, _ = telpos["body"]
    (count,) = struct.unpack_from("<i", data, start)
    columns = struct.unpack_from(f"<{4 * count}f", data, start + 4)
    telescopes = [
        {
            "x_cm": columns[i],
            "y_cm": columns[count + i],
            "z_cm": columns[2 * count + i],
            "r_cm": columns[3 * count + i],
        }
        for i in range(count)
    ]
    container = _block(items, 1204)
    bunches = []
    for child in container["children"]:
        if child["type"] != 1217:
            continue
        start, _ = child["body"]
        array, tel = struct.unpack_from("<hh", data, start)
        (photons,) = struct.unpack_from("<f", data, start + 4)
        (n,) = struct.unpack_from("<i", data, start + 8)
        values = np.frombuffer(
            data, dtype="<f4", count=10 * n, offset=start + 12
        ).reshape(n, 10)
        if child["body"][1] - child["body"][0] != 12 + 40 * n:
            raise ValueError("1217 block size does not match its count")
        bunches.append(
            {"array": array, "telescope": tel, "photons": photons, "values": values}
        )
    return {
        "observation_level_m": float(observation_level_cm) * CM,
        "telescopes": telescopes,
        "bunches": bunches,
        "blocks": [node["type"] for node in items],
    }


def pyeventio_cross_check(path=FIXTURE) -> dict:
    """Top-level block types and telescope positions as pyeventio reads them (independent check)."""
    from eventio import EventIOFile

    types, telescopes = [], None
    with EventIOFile(str(path)) as handle:
        for obj in handle:
            types.append(int(obj.eventio_type))
            if obj.eventio_type == 1201:
                telescopes = int(obj.n_telescopes)
    return {"types": types, "n_telescopes": telescopes}


# ----------------------------------------------------------------- reference rows


def reference_rows(run: dict) -> list[dict]:
    """Sphere rows of every bunch in float64 (documented rules); misses are kept with miss = True."""
    rows = []
    obs = run["observation_level_m"]
    for bunch in run["bunches"]:
        tel = run["telescopes"][bunch["telescope"]]
        radius = tel["r_cm"] * CM
        centre_z = tel["z_cm"] * CM
        for values in bunch["values"].astype(np.float64):
            x, y, z, cx, cy, cz, ctime, dist, photons, lam = values
            d = np.array([cx, cy, cz])
            dist_m = dist * CM
            arrival = np.array([x, y, z]) * CM
            emission = arrival - dist_m * d
            b = float(emission @ d)
            c2 = float(emission @ emission)
            inside = c2 <= radius * radius
            disc = b * b - (c2 - radius * radius)
            s = -b - math.sqrt(disc) if disc >= 0 else math.nan
            if not inside and not (s >= 0.0):
                rows.append({"miss": True, "telescope": bunch["telescope"] + 1})
                continue
            position, path = (emission, 0.0) if inside else (emission + s * d, s)
            rows.append(
                {
                    "miss": False,
                    "inside": bool(inside),
                    "telescope": bunch["telescope"] + 1,
                    "x": position[0],
                    "y": position[1],
                    "z": position[2],
                    "cx": cx,
                    "cy": cy,
                    "cz": cz,
                    "t": ctime + (path - dist_m) / C_M_NS,
                    "zem": obs + centre_z + emission[2],
                    "w": photons,
                    "path": path,
                    "wavelength": lam,
                }
            )
    return rows


def ulp_error(got, want) -> np.ndarray:
    """|got - want| in float32 spacings of want (spacing floored at 1 unit), per element."""
    want32 = np.asarray(want, dtype=np.float32)
    spacing = np.spacing(np.maximum(np.abs(want32), 1.0).astype(np.float32)).astype(
        np.float64
    )
    return (
        np.abs(np.asarray(got, dtype=np.float64) - want32.astype(np.float64)) / spacing
    )


def read_stream(path) -> dict:
    """from-corsika output: run metadata, run_end and every photon row in file order."""
    from obdeect.stream import Reader

    run, run_end, rows = {}, {}, []
    for record in Reader(str(path)):
        kind = record["type"]
        if kind == "run":
            run = dict(record["metadata"])
        elif kind == "run_end":
            run_end = dict(record["metadata"])
        elif kind == "photons":
            columns = record["columns"]
            for i in range(len(columns["x"])):
                row = {name: float(columns[name][i]) for name in columns}
                row["telescope"] = record["telescope_id"]
                rows.append(row)
    return {"run": run, "run_end": run_end, "rows": rows}


def compare_rows(reference: list[dict], rows: list[dict]) -> dict:
    """Column-wise ulp errors of obdeect's rows against the reference rows (misses removed)."""
    kept = [r for r in reference if not r["miss"]]
    worst = {name: 0.0 for name in COLUMNS}
    telescope_mismatch = 0
    if len(kept) != len(rows):
        return {
            "row_count_reference": len(kept),
            "row_count_obdeect": len(rows),
            "telescope_mismatch": len(kept) - len(rows),
            "max_ulp": {name: math.inf for name in COLUMNS},
        }
    for want, got in zip(kept, rows):
        telescope_mismatch += want["telescope"] != got["telescope"]
        for name in COLUMNS:
            error = float(ulp_error(got[name], want[name]))
            worst[name] = max(worst[name], error)
    return {
        "row_count_reference": len(kept),
        "row_count_obdeect": len(rows),
        "telescope_mismatch": telescope_mismatch,
        "max_ulp": worst,
    }


# ----------------------------------------------------------------- to-simtel


def downward_on_sphere(reference: list[dict]) -> np.ndarray:
    """Mask over the kept reference rows: emitted outside the sphere and moving downwards."""
    kept = [r for r in reference if not r["miss"]]
    return np.array([(not r["inside"]) and r["cz"] < 0 for r in kept])


def subset_stream(stream: dict, mask: np.ndarray, out: Path) -> Path:
    """obdeect stream of the masked rows (same run and event records), one chunk per telescope."""
    from obdeect.stream import Writer

    rows = [row for row, keep in zip(stream["rows"], mask) if keep]
    with Writer(str(out), "obdeect validation: sphere subset") as writer:
        writer.write("run", stream["run"])
        writer.write("event", {"event_number": 1})
        for tel in sorted({row["telescope"] for row in rows}):
            selected = [row for row in rows if row["telescope"] == tel]
            columns = {
                name: np.array([row[name] for row in selected], dtype=np.float32)
                for name in (
                    "x",
                    "y",
                    "z",
                    "cx",
                    "cy",
                    "cz",
                    "t",
                    "zem",
                    "w",
                    "path",
                    "wavelength",
                )
            }
            writer.write(
                "telescope", {"telescope_id": tel, "frame": "ground", "array_index": 0}
            )
            writer.write_photons(tel, "ground", 0, columns)
            writer.write(
                "telescope_end", {"telescope_id": tel, "photons": float(len(selected))}
            )
        writer.write("event_end", {})
        writer.write("run_end", {})
    return out


def simtel_rows(path: Path) -> list[dict]:
    """Bunches of a sim_telarray CORSIKA file (1205 long format, read with pyeventio), in SI units."""
    from eventio import EventIOFile

    rows = []
    with EventIOFile(str(path)) as handle:
        for top in handle:
            if top.eventio_type != 1204:
                continue
            for child in top:
                if child.eventio_type != 1205:
                    continue
                data = child.parse()[0]
                for row in data:
                    rows.append(
                        {
                            "telescope": int(child.telescope_id) + 1,
                            "x": float(row["x"]) * CM,
                            "y": float(row["y"]) * CM,
                            "cx": float(row["cx"]),
                            "cy": float(row["cy"]),
                            "t": float(row["time"]),
                            "zem": float(row["zem"]) * CM,
                            "photons": float(row["photons"]),
                            "wavelength": float(row["wavelength"]),
                        }
                    )
    return rows


def plane_reference(rows: list[dict]) -> list[dict]:
    """Plane conversion of downward sphere rows: x, y at z = 0 along the line, t + s / c with s = -z / cz."""
    out = []
    for row in rows:
        s = -row["z"] / row["cz"]
        out.append(
            {
                "telescope": row["telescope"],
                "x": row["x"] + s * row["cx"],
                "y": row["y"] + s * row["cy"],
                "t": row["t"] + s / C_M_NS,
                "zem": row["zem"],
                "photons": row["w"],
                "wavelength": row["wavelength"],
            }
        )
    return out


def status(argv, log: Path) -> int:
    """Exit status of argv; its stderr goes to log (no exception on failure)."""
    with Path(log).open("wb") as stream:
        return subprocess.run(
            [str(a) for a in argv],
            stdout=subprocess.DEVNULL,
            stderr=stream,
            check=False,
        ).returncode


def run_all(root: Path) -> tuple[dict, dict]:
    tool = common.tool("obdeect")
    pst = root / "photons3d.pst"
    common.command(
        [tool, "stream", "from-corsika", FIXTURE],
        root / "from-corsika.log",
        stdout_file=pst,
    )
    run = read_run()
    reference = reference_rows(run)
    stream = read_stream(pst)
    compared = compare_rows(reference, stream["rows"])
    run_end = stream["run_end"]
    misses_reference = sum(1 for r in reference if r["miss"])
    misses_obdeect = int(run_end.get("misses_sphere", -1))
    inside_reference = sum(1 for r in reference if not r["miss"] and r["inside"])

    # to-simtel: the full fixture is refused; the downward on-sphere subset converts to the plane.
    full_log = root / "to-simtel-full.log"
    full_status = status([tool, "stream", "to-simtel", pst], full_log)
    refused_full = int(
        full_status != 0 and "inside the telescope sphere" in full_log.read_text()
    )
    mask = downward_on_sphere(reference)
    subset = subset_stream(stream, mask, root / "subset.pst")
    converted = root / "subset.simtel"
    common.command(
        [tool, "stream", "to-simtel", subset],
        root / "to-simtel-subset.log",
        stdout_file=converted,
    )
    simtel = simtel_rows(converted)
    subset_rows = read_stream(subset)["rows"]
    expected = plane_reference(subset_rows)
    matched = len(simtel) == len(expected)
    pairs = list(zip(simtel, expected)) if matched else []
    plane = {
        key: max((abs(a[key] - b[key]) for a, b in pairs), default=math.inf)
        for key in ("x", "y", "t", "zem")
    }
    details = {
        "fixture": str(FIXTURE.relative_to(uc.ROOT)),
        "fixture_sha256": common.sha256(FIXTURE),
        "obdeect_stream_version": version(),
        "observation_level_m": run["observation_level_m"],
        "telescopes_m": [
            {k[:-3] + "_m": v * CM for k, v in tel.items()} for tel in run["telescopes"]
        ],
        "eventio_top_level_blocks": run["blocks"],
        "pyeventio_cross_check": pyeventio_cross_check(),
        "per_column_max_ulp": compared["max_ulp"],
        "to_simtel_full_status": full_status,
        "to_simtel_subset_rows": len(simtel),
        "to_simtel_plane_error": plane,
        "reference_rows": [
            {
                k: (float(v) if isinstance(v, (float, np.floating)) else v)
                for k, v in r.items()
            }
            for r in reference
        ],
    }
    metrics = {
        "row_count_obdeect": compared["row_count_obdeect"],
        "row_count_reference": compared["row_count_reference"],
        "row_count_mismatch": abs(
            compared["row_count_obdeect"] - compared["row_count_reference"]
        ),
        "telescope_mismatch": compared["telescope_mismatch"],
        "max_column_ulp": max(compared["max_ulp"].values()),
        "ground_representation_sphere": int(
            run_end is not None
            and stream["run"].get("ground_representation") == "sphere"
        ),
        "misses_sphere_obdeect": misses_obdeect,
        "misses_sphere_reference": misses_reference,
        "misses_sphere_mismatch": abs(misses_obdeect - misses_reference),
        "inside_rows_reference": inside_reference,
        "to_simtel_full_refused": refused_full,
        "to_simtel_subset_rows": len(simtel),
        "to_simtel_subset_row_mismatch": abs(len(simtel) - len(expected)),
        "to_simtel_max_position_error_m": max(plane["x"], plane["y"]),
        "to_simtel_max_time_error_ns": plane["t"],
        "to_simtel_max_zem_error_m": plane["zem"],
    }
    return metrics, details


def main() -> None:
    common.results_main(run_all, __doc__.splitlines()[0])


if __name__ == "__main__":
    main()
