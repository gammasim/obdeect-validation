"""Extinction of sphere rows: transmission against independent formulas (propagation.sphere-extinction).

Sky photons with one telescope (centre C = (0, 0, 30) m above the observation level 2000 m, sphere radius
12.5 m), transmission table fixtures/atm_trans_2000_1_10_0_0_2000.dat, no profile, no in-telescope path.
The output is weighted (w = p, store_p): p of each row is its column p. Rows are matched to the input photons by
their wavelength, a unique value per photon (the stages keep a wavelength they are given).

The reference is independent of obdeect code: upward_common.Transmission (log10-linear interpolation in
height, zero at H2, linear in wavelength between integer rows, dtau/dh of that interpolation), scipy quad.
For a sphere row: Q = C + r (the row point), E = Q - L d (emission point), h_E = zem, h_Q = observation level +
C_z + r_z, u = |cz| (flat). Spherical: Earth radius R_E = 6371315 m, centre (0, 0, -(R_E + H2)) in the sky frame,
h(X) = |X - centre| - R_E, n_Q = (Q - centre)/|Q - centre|, u = d . n_Q.
- flat, |u| >= 0.1: p = exp(-|tau(h_E) - tau(h_Q)| / |u|);
- flat or spherical, |u| < 0.1: p = exp(-int_0^L k(h(s)) ds), k = |dtau/dh| along the straight chord E -> Q;
- L = 0 (emitted inside the sphere): p = 1.

Sections (gates in cases/propagation/sphere-extinction.json):
(a) downward and upward steep rows, flat, plane formula;
(b) near-horizontal (|u| < 0.1) and exactly horizontal rows, flat and spherical, chord integral;
(c) continuity at |u| = 0.1: the two reference forms agree at u = 0.1 (depth difference, documented 1e-9), and
    obdeect's rows on both sides (|u| = 0.1 -/+ 2e-4, sightlines through the centre) match the reference;
(d) the same photons as plane rows (G plane through the telescope centre) and as sphere rows whose entry point
    lies at the telescope centre height (entry on the equator): transmissions agree to 1e-6 relative;
(e) photons emitted inside the sphere: p = 1.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from obdeect_validation.adapters import common
from obdeect_validation.adapters import upward_common as uc

SEED = 20261128
TELESCOPE_M = (0.0, 0.0, 30.0)  # sphere centre above the observation level [m]
RADIUS_M = 12.5
OBS_M = uc.SITE_M  # 2000 m: H2 of the transmission table
EARTH_RADIUS_M = 6371315.0
WAVELENGTH_BASE_NM = 400.0
WAVELENGTH_STEP_NM = 0.01
QUAD = {"limit": 2000, "epsabs": 1e-13, "epsrel": 1e-11}
NEAR_HORIZONTAL = 0.1


def version() -> str:
    try:
        return importlib.metadata.version("obdeect-propagation")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


_TABLE = None


def table() -> uc.Transmission:
    global _TABLE
    if _TABLE is None:
        _TABLE = uc.Transmission()
    return _TABLE


def unit(v) -> np.ndarray:
    v = np.asarray(v, dtype=np.float64)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def k_of(wl: float, h: float) -> float:
    return abs(float(table().dtau_dh(wl, h)))


def chord_depth(
    wl: float, h_start: float, h_end: float, length: float, geometry
) -> float:
    """int_0^L k ds along the straight chord; geometry(s) gives the altitude at distance s."""
    if length <= 0.0:
        return 0.0
    value, _ = quad(lambda s: k_of(wl, geometry(s)), 0.0, length, **QUAD)
    return value


def flat_chord_depth(wl, h_start, h_end, length) -> float:
    """Flat: h linear in s; the integral of k along the chord (equals the plane formula for any u)."""
    if length <= 0.0:
        return 0.0
    slope = (h_end - h_start) / length
    points = []
    for h_b in table().breakpoints(min(h_start, h_end), max(h_start, h_end)):
        points.append((h_b - h_start) / slope)
    value, _ = quad(
        lambda s: k_of(wl, h_start + slope * s),
        0.0,
        length,
        points=points or None,
        **QUAD,
    )
    return value


def plane_depth(wl, h_start, h_end, u) -> float:
    tau = table().tau(wl, np.array([h_start, h_end]))
    return abs(float(tau[0] - tau[1])) / abs(u)


def flat_row_depth(wl, zem, z_row_rel, cz, length) -> float:
    """Reference depth of a flat sphere row: plane formula for |u| >= 0.1, chord integral below."""
    h_e = float(zem)
    h_q = OBS_M + TELESCOPE_M[2] + float(z_row_rel)
    if length <= 0.0:
        return 0.0
    if abs(cz) >= NEAR_HORIZONTAL:
        return plane_depth(wl, h_e, h_q, cz)
    return flat_chord_depth(wl, h_e, h_q, length)


def spherical_chord_depth(wl, q_rel, direction, length) -> tuple[float, float]:
    """Spherical chord ending at Q (sky frame: Q = C + r, z above the observation level); returns (u, depth)."""
    centre = np.array([0.0, 0.0, -(EARTH_RADIUS_M + OBS_M)])
    q = np.array([TELESCOPE_M[0], TELESCOPE_M[1], TELESCOPE_M[2]]) + np.asarray(
        q_rel, dtype=np.float64
    )
    d = unit(direction)
    e = q - length * d

    def altitude(s):
        x = e + s * d - centre
        return float(np.linalg.norm(x)) - EARTH_RADIUS_M

    n_q = unit(q - centre)
    u = float(d @ n_q)
    depth = chord_depth(wl, altitude(0.0), altitude(length), length, altitude)
    return u, depth


def read_rows(path) -> dict:
    from obdeect.stream import Reader

    columns: dict[str, list] = {}
    for record in Reader(str(path)):
        if record["type"] == "photons":
            for name, values in record["columns"].items():
                columns.setdefault(name, []).append(np.asarray(values))
    return {k: np.concatenate(v) for k, v in columns.items()} if columns else {}


def ids_of(wavelength) -> np.ndarray:
    return np.rint(
        (np.asarray(wavelength, dtype=np.float64) - WAVELENGTH_BASE_NM)
        / WAVELENGTH_STEP_NM
    ).astype(np.int64)


def write_input(path: Path, columns: dict) -> Path:
    from obdeect.stream import Writer

    run = {
        "run_number": 1,
        "observation_level_m": OBS_M,
        "telescope_x_m": [TELESCOPE_M[0]],
        "telescope_y_m": [TELESCOPE_M[1]],
        "telescope_z_m": [TELESCOPE_M[2]],
        "telescope_r_m": [RADIUS_M],
        "wavelength_min_nm": 300.0,
        "wavelength_max_nm": 700.0,
        "curvature_applied": 0,
    }
    with Writer(str(path), "obdeect validation: sphere extinction") as writer:
        writer.write("run", run)
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1, "frame": 0, "array_index": 0})
        writer.write_photons(1, "sky", 0, columns)
        writer.write(
            "telescope_end", {"telescope_id": 1, "photons": float(len(columns["x"]))}
        )
        writer.write("event_end", {})
        writer.write("run_end", {})
    return path


def photon_sets(seed: int = SEED) -> dict:
    """Sky photons by group (float32 columns, unique wavelengths), and the group of each photon."""
    from obdeect_validation.adapters import sphere_geometry as sg

    rng = np.random.default_rng(seed)
    centre = np.asarray(TELESCOPE_M)
    rows, groups = [], []

    def add(group, origin, cx, cy, cz):
        rows.append((origin, cx, cy, cz))
        groups.append(group)

    for _ in range(40):  # (a) steep downward from high altitude, aimed at the sphere
        origin = np.array(
            [
                rng.uniform(-200, 200),
                rng.uniform(-200, 200),
                rng.uniform(500, 5000) + OBS_M - OBS_M,
            ]
        )
        cx, cy, cz = sg.direction_from_target(
            origin,
            centre + rng.uniform(-0.8 * RADIUS_M, 0.8 * RADIUS_M, 3) * [1, 1, 0.3],
        )
        add("a_steep", origin, cx, cy, cz)
    for _ in range(10):  # (a) steep upward from below the sphere
        origin = np.array(
            [rng.uniform(-40, 40), rng.uniform(-40, 40), rng.uniform(0.0, 15.0)]
        )
        cx, cy, cz = sg.direction_from_target(
            origin, centre + rng.uniform(-RADIUS_M, RADIUS_M, 3) * [1, 1, 0.3]
        )
        add("a_upward", origin, cx, cy, sg.stored_cz(cx, cy, False))
    for _ in range(50):  # (b) near-horizontal: 1 to 4 km away, at the telescope height
        azimuth = rng.uniform(0, 2 * math.pi)
        distance = rng.uniform(1000.0, 4000.0)
        origin = np.array(
            [
                distance * math.cos(azimuth),
                distance * math.sin(azimuth),
                centre[2] + rng.uniform(-3, 3),
            ]
        )
        cx, cy, cz = sg.direction_from_target(
            origin, centre + rng.uniform(-6, 6, 3) * [1, 1, 0.0]
        )
        add("b_near_horizontal", origin, cx, cy, cz)
    for _ in range(30):  # (b) exactly horizontal (cz = 0), at the centre height
        azimuth = rng.uniform(0, 2 * math.pi)
        distance = rng.uniform(1000.0, 4000.0)
        origin = np.array(
            [distance * math.cos(azimuth), distance * math.sin(azimuth), centre[2]]
        )
        add(
            "b_horizontal",
            origin,
            np.float32(-math.cos(azimuth)),
            np.float32(-math.sin(azimuth)),
            np.float32(0.0),
        )
    for cz_target in (
        -0.0998,
        -0.1002,
        -0.0999,
        -0.1001,
    ):  # (c) either side of |u| = 0.1, through the centre
        for _ in range(8):
            azimuth = rng.uniform(0, 2 * math.pi)
            sin_t = math.sqrt(1.0 - cz_target * cz_target)
            d = np.array(
                [sin_t * math.cos(azimuth), sin_t * math.sin(azimuth), cz_target]
            )
            origin = centre - 2500.0 * d
            add(
                "c_continuity",
                origin,
                np.float32(d[0]),
                np.float32(d[1]),
                np.float32(d[2]),
            )
    for _ in range(
        40
    ):  # (d) entry on the equator (Q at the centre height), steep direction
        phi = rng.uniform(0, 2 * math.pi)
        normal = np.array([math.cos(phi), math.sin(phi), 0.0])
        zen = rng.uniform(0.5, 0.95)  # cz magnitude
        horizontal = math.sqrt(1 - zen * zen)
        d = np.array(
            [-horizontal * normal[0], -horizontal * normal[1], -zen]
        )  # inward at the equator point
        target = centre + RADIUS_M * normal
        s = rng.uniform(1500.0, 4000.0)
        origin = target - s * d
        add("d_equator", origin, np.float32(d[0]), np.float32(d[1]), np.float32(d[2]))
    for _ in range(5):  # (e) emitted inside the sphere
        origin = centre + unit(rng.normal(size=3)) * RADIUS_M * 0.5
        azimuth, zenith = rng.uniform(0, 2 * math.pi), math.acos(rng.uniform(-1, 1))
        cx, cy = (
            math.sin(zenith) * math.cos(azimuth),
            math.sin(zenith) * math.sin(azimuth),
        )
        add(
            "e_inside",
            origin,
            np.float32(cx),
            np.float32(cy),
            np.float32(math.cos(zenith)),
        )
    n = len(rows)
    cols = {
        "x": np.array([r[0][0] for r in rows], dtype=np.float32),
        "y": np.array([r[0][1] for r in rows], dtype=np.float32),
        "z": np.array([r[0][2] for r in rows], dtype=np.float32),
        "cx": np.array([r[1] for r in rows], dtype=np.float32),
        "cy": np.array([r[2] for r in rows], dtype=np.float32),
        "cz": np.array([r[3] for r in rows], dtype=np.float32),
        "t": np.zeros(n, dtype=np.float32),
        "w": np.ones(n, dtype=np.float32),
        "wavelength": (WAVELENGTH_BASE_NM + WAVELENGTH_STEP_NM * np.arange(n)).astype(
            np.float32
        ),
    }
    return {"columns": cols, "groups": np.array(groups, dtype=object)}


def run_all(root: Path) -> tuple[dict, dict]:
    from obdeect.propagation import propagate

    sets = photon_sets()
    cols, groups = sets["columns"], sets["groups"]
    n = len(groups)
    inp = write_input(root / "sphere-extinction-in.pst", cols)
    transmission = str(uc.TRANSMISSION)
    flat = root / "sphere-flat.pst"
    propagate(
        inp,
        flat,
        transmission=transmission,
        site_altitude_m=OBS_M,
        seed=1,
        ground_representation="sphere",
        mode="weighted",
        store_p=True,
    )
    spherical_groups = np.isin(groups, ["b_near_horizontal", "b_horizontal"])
    spherical_in = root / "sphere-spherical-in.pst"
    sub = {k: v[spherical_groups] for k, v in cols.items()}
    write_input(spherical_in, sub)
    spherical = root / "sphere-spherical.pst"
    propagate(
        spherical_in,
        spherical,
        transmission=transmission,
        site_altitude_m=OBS_M,
        seed=1,
        ground_representation="sphere",
        mode="weighted",
        store_p=True,
        curvature="spherical",
    )
    plane = root / "plane-flat.pst"
    equator = np.isin(groups, ["d_equator"])
    plane_in = root / "plane-in.pst"
    write_input(plane_in, {k: v[equator] for k, v in cols.items()})
    propagate(
        plane_in,
        plane,
        transmission=transmission,
        site_altitude_m=OBS_M,
        seed=1,
        mode="weighted",
        store_p=True,
    )
    return {
        "flat": flat,
        "spherical": spherical,
        "plane": plane,
        "n": n,
        "groups": groups,
        "cols": cols,
        "spherical_mask": spherical_groups,
        "equator_mask": equator,
    }, {}


def rel(a, b) -> float:
    return (
        float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300))) if len(a) else 0.0
    )


def evaluate(paths: dict) -> tuple[dict, dict]:
    cols, groups = paths["cols"], paths["groups"]
    flat = read_rows(paths["flat"])
    spherical = read_rows(paths["spherical"])
    plane = read_rows(paths["plane"])
    wl_in = cols["wavelength"].astype(np.float64)
    index = {int(v): i for i, v in enumerate(ids_of(wl_in).tolist())}
    details = {
        "obdeect_propagation_version": version(),
        "groups": {g: int(np.count_nonzero(groups == g)) for g in sorted(set(groups))},
    }
    metrics = {}

    def compare(rows, mask_group, label, geometry):
        ids = ids_of(rows["wavelength"]).tolist()
        ref_values, obd_values, kept_idx = [], [], []
        for k, idx in enumerate(ids):
            i = index[int(idx)]
            if groups[i] not in mask_group:
                continue
            if geometry == "spherical":
                u, depth = spherical_chord_depth(
                    float(rows["wavelength"][k]),
                    np.array(
                        [rows["x"][k], rows["y"][k], rows["z"][k]], dtype=np.float64
                    ),
                    np.array(
                        [cols["cx"][i], cols["cy"][i], cols["cz"][i]], dtype=np.float64
                    ),
                    float(rows["path"][k]),
                )
                if abs(u) >= NEAR_HORIZONTAL:
                    continue
                ref = math.exp(-depth)
            else:
                ref = math.exp(
                    -flat_row_depth(
                        float(rows["wavelength"][k]),
                        float(rows["zem"][k]),
                        float(rows["z"][k]),
                        float(cols["cz"][i]),
                        float(rows["path"][k]),
                    )
                )
            ref_values.append(ref)
            obd_values.append(float(rows["p"][k]))
            kept_idx.append(i)
        ref_values = np.array(ref_values)
        obd_values = np.array(obd_values)
        out = {
            "count": len(ref_values),
            "max_rel_error": rel(obd_values, ref_values) if len(ref_values) else 0.0,
        }
        details[label] = out
        return out

    a = compare(flat, {"a_steep", "a_upward"}, "a_steep_flat", "flat")
    b = compare(flat, {"b_near_horizontal", "b_horizontal"}, "b_chord_flat", "flat")
    s = compare(
        spherical,
        {"b_near_horizontal", "b_horizontal"},
        "b_chord_spherical",
        "spherical",
    )
    c = compare(flat, {"c_continuity"}, "c_continuity_flat", "flat")
    e = compare(flat, {"e_inside"}, "e_inside", "flat")

    # geometry: every photon reaches the sphere (row present) exactly when its straight line meets it ahead
    centre = np.asarray(TELESCOPE_M, dtype=np.float64)
    n = len(cols["x"])
    hits = []
    for i in range(n):
        p0 = np.array([cols["x"][i], cols["y"][i], cols["z"][i]], dtype=np.float64)
        d0 = unit([cols["cx"][i], cols["cy"][i], cols["cz"][i]])
        offset = p0 - centre
        if offset @ offset <= RADIUS_M**2:
            hits.append(True)
            continue
        half = offset @ d0
        disc = half * half - (offset @ offset - RADIUS_M**2)
        hits.append(bool(disc >= 0.0 and -half - math.sqrt(disc) > 0.0))
    present = set(ids_of(flat["wavelength"]).tolist())
    obdeect_hit = np.array(
        [int(ids_of(cols["wavelength"][i : i + 1])[0]) in present for i in range(n)]
    )
    geometry_mismatch = int(np.count_nonzero(np.array(hits) != obdeect_hit))
    details["geometry_hits"] = {
        "reference": int(np.count_nonzero(hits)),
        "obdeect_rows": int(np.count_nonzero(obdeect_hit)),
    }

    # (c) the two reference forms agree at |u| = 0.1 (documented 1e-9 in depth)
    h_e, h_q = 2000.0 + 130.0, 2000.0 + 30.0 + 0.5
    length = (h_e - h_q) / 0.1
    form_plane = plane_depth(410.0, h_e, h_q, 0.1)
    form_chord = flat_chord_depth(410.0, h_e, h_q, length)
    # (d) plane rows and sphere rows of the same photons (wavelength identity)
    pl_ids = ids_of(plane["wavelength"]).tolist()
    sp_index = {int(v): k for k, v in enumerate(ids_of(flat["wavelength"]).tolist())}
    pairs = [(k, sp_index[int(v)]) for k, v in enumerate(pl_ids) if int(v) in sp_index]
    p_plane = np.array([plane["p"][k] for k, _ in pairs], dtype=np.float64)
    p_sphere = np.array([flat["p"][j] for _, j in pairs], dtype=np.float64)
    plane_ref = []
    for k, _ in pairs:
        i = index[int(ids_of(plane["wavelength"][k : k + 1])[0])]
        h_t = OBS_M + TELESCOPE_M[2]
        plane_ref.append(
            math.exp(
                -plane_depth(
                    float(plane["wavelength"][k]),
                    float(plane["zem"][k]),
                    h_t,
                    float(cols["cz"][i]),
                )
            )
        )
    details["d_equator"] = {
        "pairs": len(pairs),
        "plane_vs_sphere_max_rel": rel(p_sphere, p_plane) if pairs else math.inf,
        "plane_vs_reference_max_rel": rel(p_plane, np.array(plane_ref))
        if pairs
        else math.inf,
    }
    details["continuity_formula_depth_difference"] = abs(form_plane - form_chord)
    metrics = {
        "a_steep_count": a["count"],
        "a_steep_max_rel_error": a["max_rel_error"],
        "b_chord_flat_count": b["count"],
        "b_chord_flat_max_rel_error": b["max_rel_error"],
        "b_chord_spherical_count": s["count"],
        "b_chord_spherical_max_rel_error": s["max_rel_error"],
        "c_continuity_count": c["count"],
        "c_continuity_max_rel_error": c["max_rel_error"],
        "c_formula_depth_difference": abs(form_plane - form_chord),
        "d_pairs": len(pairs),
        "d_plane_vs_sphere_max_rel": details["d_equator"]["plane_vs_sphere_max_rel"],
        "d_plane_vs_reference_max_rel": details["d_equator"][
            "plane_vs_reference_max_rel"
        ],
        "geometry_mismatch_count": geometry_mismatch,
        "e_inside_count": e["count"],
        "e_inside_max_rel_error": e["max_rel_error"],
    }
    return metrics, details


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args()
    root = Path.cwd()
    paths, _ = run_all(root)
    metrics, details = evaluate(paths)
    common.write_results(metrics, details, root)


if __name__ == "__main__":
    main()
