"""Ground-frame sources with the telescope sphere: a muon-like line crossing the sphere, an isotropic point outside.

Case sources.ground-muon. One telescope at the origin (sphere radius R = 12.5 m), frame ground (sphere rows,
docs/sources/spec/run.md "Ground runs"), observation level 2000 m. Two emissions:

- line: a 120 m track at y = 2 m, z = 3 m (impact 3.6 m from the sphere centre) along +x, Cherenkov cone of
  1.2 degrees about the track, 1e5 photons, uniform along the track, beta = 1;
- point: isotropic light from (30, -40, 30) m (distance D = 57.4 m), 2e5 photons, no target.

Every reference is analytic and written here (numpy only). Each row is checked against the geometry it must
satisfy, with its emission point rebuilt from the row (P = entry - path d for rows with path > 0, P = position
for rows with path 0):
- inside rows (emitted inside the sphere): path 0, position = emission point on the track (|r| <= R);
- outside rows (entry points): |r| = R, path = |r - P|, time = (s_emit + path) / c with s_emit the distance
  of P along the track from its start (beta = 1), and zem = observation level + P_z (a.s.l.);
- the cone: the photon direction makes at most 1.2 degrees with the track;
- inside rows of the line: Binomial(N_line, chord / L) (uniform emission along the track; every photon
  emitted inside the sphere is a row); kept point rows: Binomial(N_point, (1 - sqrt(1 - (R/D)^2)) / 2);
- run_end misses_sphere equals the photons emitted minus the rows written.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common
from obdeect_validation.adapters import upward_common as uc

SEED = 20261126
OBSERVATION_M = 2000.0
RADIUS_M = 12.5
LINE_START = np.array([-60.0, 2.0, 3.0])
LINE_DIRECTION = np.array([1.0, 0.0, 0.0])
LINE_LENGTH_M = 120.0
CONE_DEG = 1.2
LINE_PHOTONS = 100_000
POINT_POSITION = np.array([30.0, -40.0, 30.0])
POINT_PHOTONS = 200_000
GEOMETRY_TOL_M = 1e-4
C_M_NS = uc.C_M_NS
EPS32 = 2.0**-24


def version() -> str:
    return common.package_version("obdeect-sources")


def description() -> dict:
    """Run file (dictionary form) of the two emissions, frame ground."""
    spectrum = {"kind": "fixed", "wavelength_nm": 400.0}
    return {
        "run": {"seed": SEED, "observation_level_m": OBSERVATION_M, "frame": "ground"},
        "telescopes": [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": RADIUS_M}],
        "events": [
            {
                "repeat": 1,
                "emissions": [
                    {
                        "kind": "line",
                        "photons": LINE_PHOTONS,
                        "start_m": LINE_START.tolist(),
                        "direction": LINE_DIRECTION.tolist(),
                        "length_m": LINE_LENGTH_M,
                        "angular": {"kind": "cone", "angle_deg": CONE_DEG},
                        "spectrum": spectrum,
                        "target": "none",
                    },
                    {
                        "kind": "point",
                        "photons": POINT_PHOTONS,
                        "position_m": POINT_POSITION.tolist(),
                        "spectrum": spectrum,
                        "target": "none",
                    },
                ],
            }
        ],
    }


def read_rows(path: Path) -> dict:
    from obdeect.stream import Reader

    columns = {}
    run_end = {}
    for record in Reader(str(path)):
        if record["type"] == "photons":
            for name, values in record["columns"].items():
                columns.setdefault(name, []).append(np.asarray(values))
        elif record["type"] == "run_end":
            run_end = dict(record["metadata"])
        elif record["type"] == "run":
            run_meta = dict(record["metadata"])
    data = {name: np.concatenate(parts) for name, parts in columns.items()}
    return {"columns": data, "run_end": run_end, "run": run_meta}


def chord_fraction() -> float:
    """Fraction of the line length inside the sphere (uniform emission along the track)."""
    u = LINE_DIRECTION / np.linalg.norm(LINE_DIRECTION)
    rel = LINE_START  # sphere centre at the origin
    b = rel @ u
    c = rel @ rel - RADIUS_M**2
    disc = b * b - c
    if disc <= 0.0:
        return 0.0
    s1, s2 = -b - math.sqrt(disc), -b + math.sqrt(disc)
    lo, hi = max(0.0, s1), min(LINE_LENGTH_M, s2)
    return max(hi - lo, 0.0) / LINE_LENGTH_M


def pull(count: int, n: int, p: float) -> float:
    return (count - n * p) / math.sqrt(n * p * (1.0 - p))


def float32_budget(value, extra=0.0):
    """One float32 spacing of the stored value (floored at 1 unit) plus the direction budget extra."""
    spacing = np.spacing(
        np.maximum(np.abs(np.asarray(value, dtype=np.float32)), 1.0).astype(np.float32)
    ).astype(np.float64)
    return spacing + extra


def check(data: dict, summary: dict) -> tuple[dict, dict]:
    c = data["columns"]
    pos = np.stack([c["x"], c["y"], c["z"]], axis=1).astype(np.float64)
    d = np.stack([c["cx"], c["cy"], c["cz"]], axis=1).astype(np.float64)
    d = d / np.linalg.norm(d, axis=1, keepdims=True)
    path = c["path"].astype(np.float64)
    t = c["t"].astype(np.float64)
    zem = c["zem"].astype(np.float64)
    emission = (
        pos - path[:, None] * d
    )  # rebuilt emission point (path 0: the position itself)
    u = LINE_DIRECTION / np.linalg.norm(LINE_DIRECTION)

    rel = emission - LINE_START
    along = rel @ u
    off = np.linalg.norm(rel - along[:, None] * u, axis=1)
    is_point = np.linalg.norm(emission - POINT_POSITION, axis=1) < 1e-3
    is_line = ~is_point
    inside = path == 0.0
    r = np.linalg.norm(pos, axis=1)

    # line rows
    line_inside = is_line & inside
    line_outside = is_line & ~inside
    line_geometry = {
        "on_track_max_m": float(np.max(off[is_line])) if is_line.any() else math.inf,
        "along_range_violations": int(
            np.count_nonzero(
                (along[is_line] < -GEOMETRY_TOL_M)
                | (along[is_line] > LINE_LENGTH_M + GEOMETRY_TOL_M)
            )
        ),
        "inside_radius_violations": int(
            np.count_nonzero(r[line_inside] > RADIUS_M + GEOMETRY_TOL_M)
        ),
        "outside_radius_error_m": float(np.max(np.abs(r[line_outside] - RADIUS_M)))
        if line_outside.any()
        else 0.0,
    }
    cos_cone = np.clip(d[is_line] @ u, -1.0, 1.0)
    cone_excess = (
        float(np.max(np.degrees(np.arccos(cos_cone)) - CONE_DEG))
        if is_line.any()
        else 0.0
    )
    # times and zem against the analytic emission (s_emit + path) / c and observation + P_z
    # point emission at t = 0 (default pulse); line emission at s_emit / (beta c) with beta = 1
    s_emit = np.where(is_point, 0.0, along)
    t_ref = (s_emit + path) / C_M_NS
    zem_ref = OBSERVATION_M + emission[:, 2]
    t_ratio = np.abs(t - t_ref) / float32_budget(
        t_ref, 2 * EPS32 * (s_emit + path) / C_M_NS
    )
    zem_ratio = np.abs(zem - zem_ref) / float32_budget(zem_ref)
    # counts
    n_line = LINE_PHOTONS
    n_point = POINT_PHOTONS
    chord = chord_fraction()
    inside_count = int(np.count_nonzero(line_inside))
    p_point = (
        1.0 - math.sqrt(1.0 - (RADIUS_M / np.linalg.norm(POINT_POSITION)) ** 2)
    ) / 2.0
    kept_point = int(np.count_nonzero(is_point))
    rows = int(c["x"].size)
    misses = int(summary.get("misses_sphere", -1))
    run_end_misses = int(data["run_end"].get("misses_sphere", -1))
    total = n_line + n_point
    metrics = {
        "rows": rows,
        "rows_plus_misses_minus_emitted": rows + run_end_misses - total,
        "misses_sphere_run_end": run_end_misses,
        "misses_sphere_summary_mismatch": abs(run_end_misses - misses),
        "ground_representation_sphere": int(
            data["run"].get("ground_representation") == "sphere"
        ),
        "line_rows_outside": int(np.count_nonzero(line_outside)),
        "line_rows_inside": inside_count,
        "line_inside_pull": pull(inside_count, n_line, chord) if 0 < chord < 1 else 0.0,
        "point_rows": kept_point,
        "point_kept_pull": pull(kept_point, n_point, p_point),
        "line_on_track_max_m": line_geometry["on_track_max_m"],
        "line_along_violations": line_geometry["along_range_violations"],
        "inside_radius_violations": line_geometry["inside_radius_violations"],
        "inside_path_nonzero": int(np.count_nonzero(inside & (path != 0.0))),
        "outside_path_zero": int(np.count_nonzero(~inside & (path <= 0.0))),
        "entry_radius_max_error_m": line_geometry["outside_radius_error_m"],
        "cone_excess_deg": max(cone_excess, 0.0),
        "time_budget_ratio_max": float(np.max(t_ratio)),
        "zem_budget_ratio_max": float(np.max(zem_ratio)),
        "point_position_max_error_m": float(
            np.max(np.linalg.norm(emission[is_point] - POINT_POSITION, axis=1))
        )
        if kept_point
        else 0.0,
    }
    details = {
        "obdeect_sources_version": version(),
        "chord_fraction": chord,
        "point_probability": p_point,
        "line_photons": n_line,
        "point_photons": n_point,
        "rows": rows,
        "misses_sphere": run_end_misses,
        "line_rows_inside": inside_count,
        "line_rows_outside": int(np.count_nonzero(line_outside)),
        "point_rows": kept_point,
        "emitted_total": total,
    }
    return metrics, details


def run_all(root: Path) -> tuple[dict, dict]:
    from obdeect.sources import run

    path = root / "ground-muon.pst"
    summary = run.write(description(), path)
    data = read_rows(path)
    metrics, details = check(data, summary)
    details["summary"] = {
        k: (v if isinstance(v, (int, float, str)) else str(v))
        for k, v in summary.items()
    }
    return metrics, details


def main() -> None:
    common.results_main(run_all, __doc__.splitlines()[0])


if __name__ == "__main__":
    main()
