"""Sphere representation of the atmosphere stage: sky photons to telescope-sphere rows.

Case propagation.sphere-geometry. One telescope (centre C = (0, 0, 20) m above the observation level,
sphere radius R = 12.5 m), photons in the sky frame in all directions: downward, horizontal, upward from
below the sphere, emitted inside it, emitted below the observation level, and aimed past it.

(a) No refraction (extinction off, no profile): the row is the straight-line entry point (smallest
    s > 0 with |P + s d - C| = R), its path s and time t0 + s / c; inside rows are the emission point with
    path 0 and time t0; below the observation level and misses are dropped with w = 0 and counted by
    reason. Reference: numpy, float64 from the stored float32 values, rounded to float32 at the end.
(b) Refraction (atmospheric profile of the propagation fixtures, plane-parallel layers): the entry point
    and the time are those of the 3D ray equation d(n d)/ds = grad n, dT/ds = n/c, with n = 1 + nu(h) from
    the documented interpolation (ln nu a natural cubic spline through the levels, tabulated on a 3 m grid,
    exponential inside a cell: upward_common.Profile). Integrated with scipy solve_ivp (DOP853, rtol 1e-12)
    from the emission to the sphere; the sphere is entered at the first event |X - C| = R. Accuracy gates:
    the hybrid rule and the integrator agree with this reference to 1e-4 m and 1e-4 ns (hybrid rows: steep
    rays, |cz| >= 0.1 at the telescope plane; integrator rows: near-horizontal rays); a turning ray (upward,
    apex about 285 m above its emission, back down on its descending branch) is run alone with the sphere
    centred on its descending crossing of the telescope plane: 5e-3 m and 5e-3 ns.
(c) Determinism: the refracted set gives identical rows for chunk sizes 1, 7 and all rows, threads 1 and 4,
    and a permuted input order (rows matched by photon index).

Direction convention (documented, refraction.md): |cz| = sqrt(1 - cx^2 - cy^2) from the stored cx, cy; the
stored cz is computed the same way and rounded to float32, so both conventions agree to float32 rounding.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from obdeect_validation.adapters import common, streams
from obdeect_validation.adapters import endtoend_common as e2e
from obdeect_validation.adapters import upward_common as uc

SEED = 20261125
TELESCOPE_M = (0.0, 0.0, 20.0)  # sphere centre above the observation level [m]
RADIUS_M = 12.5
C_M_NS = uc.C_M_NS
#: Two-phase ray search: coarse steps until SHELL_M outside the sphere, fine steps (FINE_STEP_M) inside.
SHELL_M = 100.0
FINE_STEP_M = 0.05
#: Path bound of the turning ray: apex about 870 m, back to the telescope plane near 560 km.
TURNING_MAX_PATH_M = 700000.0
EPS32 = 2.0**-24
ABS_TOL = 1e-10
REL_TOL = 1e-12
GATE_TIME_NS = 1e-4  # times [ns], plus one float32 spacing of the stored time
GATE_TURNING_NS = 5e-3  # turning ray time [ns], plus one float32 spacing
DROP_KEYS = ("below_ground", "misses_sphere")
#: Photon identity: a unique wavelength per input photon (the stage keeps a wavelength it is given);
#: the output drops photons that miss the sphere, so rows are matched by this value.
ID_BASE_NM = np.float32(300.0)
ID_STEP_NM = np.float32(0.01)


def version() -> str:
    return common.package_version("obdeect-propagation")


def unit(v) -> np.ndarray:
    v = np.asarray(v, dtype=np.float64)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def stored_cz(cx, cy, negative) -> np.float32:
    """cz from cx, cy (documented convention), rounded to float32; sign: negative = downward."""
    magnitude = math.sqrt(max(1.0 - float(cx) ** 2 - float(cy) ** 2, 0.0))
    return np.float32(-magnitude if negative else magnitude)


def direction_from_target(origin, target):
    """Unit vector from origin to target, as float32 (cx, cy, cz) with the cz convention above."""
    d = unit(np.asarray(target, float) - np.asarray(origin, float))
    cx, cy = np.float32(d[0]), np.float32(d[1])
    return cx, cy, stored_cz(cx, cy, d[2] < 0.0)


# ------------------------------------------------------------------ photon sets


def photon_set(seed: int = SEED) -> dict:
    """Sky-frame photons (float32 columns) of every class; each row carries its class name."""
    rng = np.random.default_rng(seed)
    centre = np.asarray(TELESCOPE_M)
    rows = []

    def add(name, origin, cx, cy, cz, t=None):
        rows.append(
            (name, origin, cx, cy, cz, rng.uniform(0.0, 50.0) if t is None else t)
        )

    # downward from above: aimed near the sphere (about half hit), some past it
    for k in range(60):
        origin = np.array(
            [rng.uniform(-300, 300), rng.uniform(-300, 300), rng.uniform(100, 3000)]
        )
        offset = rng.uniform(-2.2 * RADIUS_M, 2.2 * RADIUS_M, 3) * np.array([1, 1, 0.3])
        cx, cy, cz = direction_from_target(origin, centre + offset)
        add("downward", origin, cx, cy, cz)
    # horizontal (cz = 0 exactly): heights around the sphere centre, four exactly at the centre height
    for k in range(40):
        azimuth = rng.uniform(0, 2 * math.pi)
        origin = np.array(
            [
                rng.uniform(200, 800) * math.cos(azimuth),
                rng.uniform(200, 800) * math.sin(azimuth),
                0.0,
            ]
        )
        origin[2] = (
            centre[2] + rng.uniform(-1.6, 1.6) * RADIUS_M if k >= 4 else centre[2]
        )
        direction = np.array([-math.cos(azimuth), -math.sin(azimuth)])
        add(
            "horizontal",
            origin,
            np.float32(direction[0]),
            np.float32(direction[1]),
            np.float32(0.0),
        )
    # upward from below the sphere (observation level to 15 m), aimed up at it
    for k in range(60):
        origin = np.array(
            [rng.uniform(-40, 40), rng.uniform(-40, 40), rng.uniform(0.0, 15.0)]
        )
        offset = rng.uniform(-2.0 * RADIUS_M, 2.0 * RADIUS_M, 3)
        target = centre + offset
        target[2] = max(target[2], origin[2] + 1.0)
        cx, cy, cz = direction_from_target(origin, target)
        add("upward", origin, cx, cy, stored_cz(cx, cy, False))
    # emitted inside the sphere (any direction)
    for k in range(20):
        inside = (
            unit(rng.normal(size=3)) * RADIUS_M * 0.9 * rng.uniform(0, 1) ** (1 / 3)
        )
        origin = centre + inside
        azimuth, zenith = rng.uniform(0, 2 * math.pi), math.acos(rng.uniform(-1, 1))
        cx = np.float32(math.sin(zenith) * math.cos(azimuth))
        cy = np.float32(math.sin(zenith) * math.sin(azimuth))
        add("inside", origin, cx, cy, stored_cz(cx, cy, math.cos(zenith) < 0))
    # below the observation level, aimed at the sphere: below_ground
    for k in range(12):
        origin = np.array(
            [rng.uniform(-60, 60), rng.uniform(-60, 60), rng.uniform(-30.0, -1.0)]
        )
        cx, cy, cz = direction_from_target(origin, centre + rng.uniform(-5, 5, 3))
        add("below", origin, cx, cy, cz)
    # upward from high above, away from the sphere: behind the emission
    for k in range(12):
        origin = np.array(
            [rng.uniform(-200, 200), rng.uniform(-200, 200), rng.uniform(100, 500)]
        )
        cx, cy, cz = direction_from_target(
            origin, origin + np.array([rng.uniform(-1, 1), rng.uniform(-1, 1), 1.0])
        )
        add("behind", origin, cx, cy, cz)
    # steep upward rays from the ground: the sphere is above them but the line passes by
    for k in range(10):
        origin = np.array(
            [rng.uniform(-40, 40), rng.uniform(-40, 40), rng.uniform(0.0, 15.0)]
        )
        target = centre + np.array(
            [rng.uniform(3.0, 6.0) * RADIUS_M * rng.choice([-1, 1]), 0.0, 0.0]
        )
        cx, cy, cz = direction_from_target(origin, target)
        add("upward_aside", origin, cx, cy, stored_cz(cx, cy, False))
    # near-horizontal rays at the telescope altitude, a few km away: refraction bends them into the sphere
    for k in range(20):
        azimuth = rng.uniform(0, 2 * math.pi)
        distance = rng.uniform(2500.0, 4000.0)
        origin = np.array(
            [
                distance * math.cos(azimuth),
                distance * math.sin(azimuth),
                centre[2] + rng.uniform(-3.0, 3.0),
            ]
        )
        target = centre + rng.uniform(-6.0, 6.0, 3) * np.array([1, 1, 0.0])
        cx, cy, cz = direction_from_target(origin, target)
        add("near_horizontal", origin, cx, cy, cz)
    # steep upward from the ground with cz >= 0.1 (plane-mode hybrid rule): a set at 10 to 70 degrees
    for k in range(10):
        origin = np.array(
            [rng.uniform(-40, 40), rng.uniform(-40, 40), rng.uniform(0.0, 15.0)]
        )
        elevation = math.radians(rng.uniform(10.0, 70.0))
        azimuth = rng.uniform(0, 2 * math.pi)
        cx = np.float32(math.cos(elevation) * math.cos(azimuth))
        cy = np.float32(math.cos(elevation) * math.sin(azimuth))
        add("upward_steep", origin, cx, cy, stored_cz(cx, cy, False))
    # steep downward from high altitude (hybrid rule, plane-mode bending)
    for k in range(10):
        origin = np.array(
            [rng.uniform(-800, 800), rng.uniform(-800, 800), rng.uniform(300, 8000)]
        )
        offset = rng.uniform(-1.5 * RADIUS_M, 1.5 * RADIUS_M, 3) * np.array([1, 1, 0.0])
        cx, cy, cz = direction_from_target(origin, centre + offset)
        add("downward_steep", origin, cx, cy, cz)
    names = [r[0] for r in rows]
    cols = {
        "x": np.array([r[1][0] for r in rows], dtype=np.float32),
        "y": np.array([r[1][1] for r in rows], dtype=np.float32),
        "z": np.array([r[1][2] for r in rows], dtype=np.float32),
        "cx": np.array([r[2] for r in rows], dtype=np.float32),
        "cy": np.array([r[3] for r in rows], dtype=np.float32),
        "cz": np.array([r[4] for r in rows], dtype=np.float32),
        "t": np.array([r[5] for r in rows], dtype=np.float32),
        "w": np.ones(len(rows), dtype=np.float32),
        "wavelength": ID_BASE_NM + ID_STEP_NM * np.arange(len(rows), dtype=np.float32),
    }
    return {"columns": cols, "class": names}


def turning_photon() -> dict:
    """Upward ray from the ground at sin(theta) = 1 - 2e-5 (apex about 285 m), columns of one photon."""
    sin_e = np.float32(1.0 - 2e-5)
    cx = sin_e
    cy = np.float32(0.0)
    cz = stored_cz(cx, cy, False)
    return {
        "columns": {
            "x": np.array([0.0], dtype=np.float32),
            "y": np.array([0.0], dtype=np.float32),
            "z": np.array([0.0], dtype=np.float32),
            "cx": np.array([cx], dtype=np.float32),
            "cy": np.array([cy], dtype=np.float32),
            "cz": np.array([cz], dtype=np.float32),
            "t": np.array([0.0], dtype=np.float32),
            "w": np.ones(1, dtype=np.float32),
            "wavelength": np.full(1, ID_BASE_NM, dtype=np.float32),  # identity 0
        },
        "class": ["turning"],
    }


# ------------------------------------------------------------------ references


class RayReference:
    """Ray equation in the flat atmosphere (heights above the observation level, sky frame)."""

    def __init__(self, profile: uc.Profile, observation_m: float):
        self.profile = profile
        self.obs = observation_m
        self.grid = profile.grid
        self.lnnu = profile.lnnu
        self.slope = np.diff(self.lnnu) / np.diff(self.grid)

    def nu_and_slope(self, height_m: float):
        if height_m > self.profile.top:
            return 0.0, 0.0
        nu = float(self.profile.nu(height_m))
        i = int(
            np.clip(
                np.searchsorted(self.grid, height_m, side="right") - 1,
                0,
                len(self.slope) - 1,
            )
        )
        return nu, nu * float(self.slope[i])

    def rhs(self, s, y):
        x = y[:3]
        d = y[3:6]
        nu, dnu = self.nu_and_slope(self.obs + x[2])
        n = 1.0 + nu
        g = np.array([0.0, 0.0, dnu])
        dd = (g - (g @ d) * d) / n
        return np.concatenate([d, dd, [n / C_M_NS]])

    def trace(self, origin, direction, t0: float, centre, radius: float, s_max: float):
        """Integrate from origin until the sphere is entered, or until the ray is judged to miss.

        Two phases: coarse steps (50 m) to a shell SHELL_M outside the sphere, then fine steps (FINE_STEP_M)
        inside the shell, where the entry event is searched. A coarse step can cross a chord of the sphere
        (about 14 m) without a sign change, so the entry is never searched with coarse steps.
        Returns (entered, s, position relative to the centre, time) with time = t0 + integral n / c.
        """
        origin = np.asarray(origin, dtype=np.float64)
        d = unit(direction)
        centre = np.asarray(centre, dtype=np.float64)
        shell = radius + SHELL_M

        def outer(s, y):
            return float(np.linalg.norm(y[:3] - centre) - shell)

        def closest(s, y):
            return float((y[:3] - centre) @ y[3:6])

        def below(s, y):
            return float(y[2])

        def entry(s, y):
            return float((y[:3] - centre) @ (y[:3] - centre) - radius * radius)

        outer.terminal, outer.direction = True, -1
        closest.terminal, closest.direction = True, 1
        below.terminal, below.direction = True, -1
        entry.terminal, entry.direction = True, -1
        y0 = np.concatenate([origin, d, [t0]])
        if np.linalg.norm(origin - centre) > shell:
            coarse = solve_ivp(
                self.rhs,
                (0.0, s_max),
                y0,
                method="DOP853",
                rtol=REL_TOL,
                atol=ABS_TOL,
                events=[outer, closest, below],
                max_step=50.0,
            )
            if not coarse.t_events[0].size:
                return False, math.nan, np.full(3, math.nan), math.nan
            y0 = coarse.y_events[0][0]
            s_offset = float(coarse.t_events[0][0])
        else:
            s_offset = 0.0
        fine = solve_ivp(
            self.rhs,
            (0.0, s_max),
            y0,
            method="DOP853",
            rtol=REL_TOL,
            atol=ABS_TOL,
            events=[entry, closest, below],
            max_step=FINE_STEP_M,
        )
        if fine.t_events[0].size:
            y = fine.y_events[0][0]
            return (
                True,
                s_offset + float(fine.t_events[0][0]),
                y[:3] - centre,
                float(y[6]),
            )
        return False, math.nan, np.full(3, math.nan), math.nan


def reference_no_refraction(cols: dict) -> dict:
    """(a): straight lines, float64 from the stored values."""
    centre = np.asarray(TELESCOPE_M, dtype=np.float64)
    n = len(cols["x"])
    out = {
        k: np.zeros(n)
        for k in ("x", "y", "z", "cx", "cy", "cz", "t", "zem", "w", "path")
    }
    reason = np.array(
        [
            "",
        ]
        * n,
        dtype=object,
    )
    for i in range(n):
        P = np.array([cols["x"][i], cols["y"][i], cols["z"][i]], dtype=np.float64)
        stored = np.array(
            [cols["cx"][i], cols["cy"][i], cols["cz"][i]], dtype=np.float64
        )
        d = unit(
            stored
        )  # float32 components are unit to ~6e-8; the geometry uses the unit vector
        t0 = float(cols["t"][i])
        if P[2] < 0.0:
            reason[i] = "below_ground"
            continue
        rel = P - centre
        if rel @ rel <= RADIUS_M**2:
            pos, s, time = rel, 0.0, t0
        else:
            b = rel @ d
            disc = b * b - (rel @ rel - RADIUS_M**2)
            if disc < 0.0:
                reason[i] = "misses_sphere"
                continue
            s = -b - math.sqrt(disc)
            if s <= 0.0:
                reason[i] = "misses_sphere"
                continue
            if d[2] < 0.0 and P[2] / (-d[2]) < s:
                reason[i] = "below_ground"
                continue
            pos = rel + s * d
            time = t0 + s / C_M_NS
        out["x"][i], out["y"][i], out["z"][i] = pos
        out["cx"][i], out["cy"][i], out["cz"][i] = stored
        out["t"][i] = time
        out["zem"][i] = P[2] + uc.SITE_M
        out["w"][i] = 1.0
        out["path"][i] = s
        reason[i] = "ok"
    return {"values": out, "reason": reason}


def reference_refraction(cols: dict, profile: uc.Profile, centre, s_max: float):
    """(b): integrated rays; returns per-photon values and reasons (scipy solve_ivp, DOP853)."""
    ref = RayReference(profile, uc.SITE_M)
    c = np.asarray(centre, dtype=np.float64)
    n = len(cols["x"])
    out = {
        k: np.zeros(n)
        for k in ("x", "y", "z", "cx", "cy", "cz", "t", "zem", "w", "path")
    }
    reason = np.array(
        [
            "",
        ]
        * n,
        dtype=object,
    )
    for i in range(n):
        P = np.array([cols["x"][i], cols["y"][i], cols["z"][i]], dtype=np.float64)
        stored = np.array(
            [cols["cx"][i], cols["cy"][i], cols["cz"][i]], dtype=np.float64
        )
        d = unit(stored)
        t0 = float(cols["t"][i])
        if P[2] < 0.0:
            reason[i] = "below_ground"
            continue
        rel = P - c
        if rel @ rel <= RADIUS_M**2:
            pos, s, time = rel, 0.0, t0
            reason[i] = "ok"
        else:
            entered, s, pos, time = ref.trace(P, d, t0, c, RADIUS_M, s_max)
            if not entered:
                reason[i] = "misses_sphere"
                continue
            reason[i] = "ok"
        out["x"][i], out["y"][i], out["z"][i] = pos
        out["cx"][i], out["cy"][i], out["cz"][i] = stored
        out["t"][i] = time
        out["zem"][i] = P[2] + uc.SITE_M
        out["w"][i] = 1.0
        out["path"][i] = s
    return {"values": out, "reason": reason}


# ------------------------------------------------------------------ obdeect runs


def run_obdeect(
    root: Path,
    name: str,
    cols: dict,
    telescope=TELESCOPE_M,
    chunk_rows=None,
    indices=None,
    threads=1,
    profile=None,
):
    """Sky stream of cols (one telescope, id 1) through obdeect propagation (Python API, sphere rows)."""
    from obdeect.propagation import propagate

    stream = uc.write_stream(
        root / f"{name}-in.pst",
        [telescope],
        [(1, cols)],
        chunk_rows=chunk_rows,
        indices=None if indices is None else {1: indices},
        run_extra={"telescope_r_m": [RADIUS_M]},
    )
    out = root / f"{name}-out.pst"
    kwargs = {"profile": profile} if profile else {}
    summary = propagate(
        stream,
        out,
        extinction=False,
        site_altitude_m=uc.SITE_M,
        seed=1,
        ground_representation="sphere",
        threads=threads,
        **kwargs,
    )
    return summary, e2e.photon_rows(out), out


def match_by_index(rows: dict, n: int) -> dict:
    """Output rows placed at their input photon (identity = wavelength); `present` marks survivors.

    The output keeps only the photons that reach the sphere (weight > 0 in thinned and no-extinction
    runs), numbered consecutively, so the position in the output is not the input photon.
    """
    ids = np.rint(
        (np.asarray(rows["wavelength"], dtype=np.float64) - float(ID_BASE_NM))
        / float(ID_STEP_NM)
    ).astype(np.int64)
    if ids.size and (
        ids.min() < 0 or ids.max() >= n or len(np.unique(ids)) != ids.size
    ):
        raise ValueError("output photon identities are not a subset of the input")
    out = {"present": np.zeros(n, dtype=bool)}
    out["present"][ids] = True
    for key in OUTPUT_COLUMNS:
        full = np.zeros(n, dtype=np.float32)
        full[ids] = np.asarray(rows[key], dtype=np.float32)
        out[key] = full
    return out


# ------------------------------------------------------------------ comparison


OUTPUT_COLUMNS = ("x", "y", "z", "cx", "cy", "cz", "t", "zem", "w", "path")


def float32_budget(value, path, time: bool = False) -> np.ndarray:
    """Allowed error of a stored float32 value: one float32 spacing plus the direction budget.

    The stored direction is unit to 6e-8 (float32); a path s then moves the row by 2 eps s (eps = 2^-24).
    Times add the same budget divided by c. Spacing floored at 1 unit (a metre or a ns).
    """
    direction = 2.0 * EPS32 * np.asarray(path, dtype=np.float64)
    return streams.ulp32(value) + (direction / C_M_NS if time else direction)


def run_all(root: Path) -> tuple[dict, dict]:
    profile = uc.profile()
    photons = photon_set()
    cols, names = photons["columns"], np.array(photons["class"], dtype=object)
    n = len(names)
    details = {
        "obdeect_propagation_version": version(),
        "telescope_m": TELESCOPE_M,
        "radius_m": RADIUS_M,
    }

    # (a) no refraction: straight lines; ratios to the float32 budget (1 = at the budget)
    summary_a, rows_a_raw, _ = run_obdeect(root, "no-refraction", cols)
    rows_a = match_by_index(rows_a_raw, n)
    ref_a = reference_no_refraction(cols)
    kept_a = ref_a["reason"] == "ok"
    present_a = rows_a["present"]
    ratios = {}
    for name in ("x", "y", "z", "t", "zem", "path", "cx", "cy", "cz"):
        want = ref_a["values"][name][kept_a]
        got = rows_a[name][kept_a]
        path = ref_a["values"]["path"][kept_a]
        err = np.abs(got.astype(np.float64) - want)
        budget = (
            float32_budget(want, path, time=(name == "t"))
            if name in ("x", "y", "z", "t", "path")
            else np.maximum(float32_budget(want, path, False), 0.0)
        )
        ratios[name] = float(np.max(err / budget)) if want.size else 0.0
    inside_mask = names == "inside"
    drop_a = {
        k: int(np.count_nonzero(ref_a["reason"] == k))
        for k in ("below_ground", "misses_sphere")
    }
    details["no_refraction"] = {
        "classes": {c: int(np.count_nonzero(names == c)) for c in sorted(set(names))},
        "reference_drops": drop_a,
        "obdeect_drops": {
            k: int(summary_a["dropped_by_reason"].get(k, 0)) for k in DROP_KEYS
        },
        "budget_ratio_max": ratios,
        "inside_path_max_m": float(np.max(np.abs(rows_a["path"][inside_mask])))
        if inside_mask.any()
        else 0.0,
    }

    # (b) refraction (profile): 3D ray integrator reference; hybrid rule and integrator rows
    summary_b, rows_b_raw, _ = run_obdeect(root, "refraction", cols, profile=uc.PROFILE)
    rows_b = match_by_index(rows_b_raw, n)
    # Longest path needed: from the farthest origin to the sphere, with margin for bending.
    origins = np.stack([cols["x"], cols["y"], cols["z"]], axis=1).astype(np.float64)
    s_max = (
        1.2 * float(np.max(np.linalg.norm(origins - TELESCOPE_M, axis=1))) + RADIUS_M
    )
    ref_b = reference_refraction(cols, profile, TELESCOPE_M, s_max=s_max)
    kept_b = ref_b["reason"] == "ok"
    present_b = rows_b["present"]
    reason_mismatch = int(np.count_nonzero(kept_b != present_b))
    both = kept_b & present_b
    pos_err = np.full(n, np.nan)
    time_err = np.full(n, np.nan)
    time_ratio = np.full(n, np.nan)
    pos_err[both] = np.sqrt(
        sum(
            (rows_b[k][both].astype(np.float64) - ref_b["values"][k][both]) ** 2
            for k in ("x", "y", "z")
        )
    )
    time_err[both] = np.abs(
        rows_b["t"][both].astype(np.float64) - ref_b["values"]["t"][both]
    )
    # time gate: GATE_TIME_NS plus one float32 spacing of the stored time (its resolution at 1e4 ns)
    time_ratio[both] = time_err[both] / (
        GATE_TIME_NS + float32_budget(ref_b["values"]["t"][both], 0.0)
    )
    # Rule of docs/propagation/refraction.md ("Hybrid rule"): steep rays (|cz| >= 0.1 at the telescope plane)
    # whose emission lies farther than R from the plane use the plane-mode bending and the tangent line;
    # the others use the integrator. |cz| at the emission stands for |cz| at the plane (bending is small there).
    emission_far = np.abs(cols["z"].astype(np.float64) - TELESCOPE_M[2]) > RADIUS_M
    steep = np.abs(cols["cz"].astype(np.float64)) >= 0.1
    hybrid_mask = both & steep & emission_far
    integrator_mask = both & ~(steep & emission_far)

    def worst(values, mask):
        return float(np.max(values[mask])) if mask.any() else 0.0

    details["refraction"] = {
        "reference_drops": {
            k: int(np.count_nonzero(ref_b["reason"] == k))
            for k in ("below_ground", "misses_sphere")
        },
        "obdeect_drops": {
            k: int(summary_b["dropped_by_reason"].get(k, 0)) for k in DROP_KEYS
        },
        "hybrid_position_m_max": worst(pos_err, hybrid_mask),
        "hybrid_time_ns_max": worst(time_err, hybrid_mask),
        "integrator_position_m_max": worst(pos_err, integrator_mask),
        "integrator_time_ns_max": worst(time_err, integrator_mask),
        "integrator_time_gate_ratio_max": worst(time_ratio, integrator_mask),
        "hybrid_time_gate_ratio_max": worst(time_ratio, hybrid_mask),
        "hybrid_count": int(np.count_nonzero(hybrid_mask)),
        "integrator_count": int(np.count_nonzero(integrator_mask)),
        "compared_count": int(np.count_nonzero(both)),
        "worst_hybrid_photon": {
            "index": int(np.nanargmax(np.where(hybrid_mask, pos_err, np.nan)))
            if hybrid_mask.any()
            else -1,
        },
    }

    # turning ray: its own telescope on the descending branch
    turning = turning_photon()
    ref_t = turning_reference(turning["columns"])
    _summary_t, rows_t_raw, _ = run_obdeect(
        root,
        "turning",
        turning["columns"],
        telescope=ref_t["centre"],
        profile=uc.PROFILE,
    )
    rows_t = match_by_index(rows_t_raw, 1)
    turn_pos = float(
        np.sqrt(
            sum(
                (float(rows_t[k][0]) - ref_t["position"][i]) ** 2
                for i, k in enumerate(("x", "y", "z"))
            )
        )
    )
    turn_time = float(rows_t["t"][0]) - ref_t["time"]
    turn_budget = float(np.spacing(np.float32(ref_t["time"])))
    details["turning"] = {
        "centre_m": list(map(float, ref_t["centre"])),
        "apex_m_above_emission": ref_t["apex_m"],
        "path_m": ref_t["path_m"],
        "reference_entry_position": list(map(float, ref_t["position"])),
        "obdeect_entry_position": [float(rows_t[k][0]) for k in ("x", "y", "z")],
        "reference_time_ns": ref_t["time"],
        "obdeect_time_ns": float(rows_t["t"][0]),
        "time_float32_spacing_ns": turn_budget,
        "obdeect_present": bool(rows_t["present"][0]),
    }

    # (c) determinism of the refracted set: chunk sizes, threads, permuted order
    baseline = rows_b
    variants = {}
    for label, kwargs in {
        "chunk1": {"chunk_rows": 1},
        "chunk7": {"chunk_rows": 7},
        "threads4": {"threads": 4},
    }.items():
        _, raw, _ = run_obdeect(
            root, f"det-{label}", cols, profile=uc.PROFILE, **kwargs
        )
        variants[label] = match_by_index(raw, n)
    perm = np.random.default_rng(SEED + 7).permutation(n)
    permuted_cols = {k: v[perm] for k, v in cols.items()}
    _, raw_perm, _ = run_obdeect(
        root, "det-permuted", permuted_cols, profile=uc.PROFILE, indices=perm
    )
    variants["permuted"] = match_by_index(raw_perm, n)
    mismatch = {}
    for label, rows in variants.items():
        differ = int(np.count_nonzero(rows["present"] != baseline["present"]))
        for name in OUTPUT_COLUMNS:
            differ += int(
                np.count_nonzero(
                    rows[name].view(np.uint32) != baseline[name].view(np.uint32)
                )
            )
        mismatch[label] = differ
    details["determinism_bit_mismatches"] = mismatch

    metrics = {
        "no_refraction_rows": n,
        "no_refraction_kept_mismatch": int(np.count_nonzero(kept_a != present_a)),
        "no_refraction_below_ground_mismatch": abs(
            drop_a["below_ground"]
            - summary_a["dropped_by_reason"].get("below_ground", 0)
        ),
        "no_refraction_misses_mismatch": abs(
            drop_a["misses_sphere"]
            - summary_a["dropped_by_reason"].get("misses_sphere", 0)
        ),
        "no_refraction_budget_ratio_max": max(ratios.values()),
        "no_refraction_inside_path_m": details["no_refraction"]["inside_path_max_m"],
        "refraction_kept_mismatch": reason_mismatch,
        "refraction_hybrid_position_m": details["refraction"]["hybrid_position_m_max"],
        "refraction_hybrid_time_gate_ratio": details["refraction"][
            "hybrid_time_gate_ratio_max"
        ],
        "refraction_integrator_position_m": details["refraction"][
            "integrator_position_m_max"
        ],
        "refraction_integrator_time_gate_ratio": details["refraction"][
            "integrator_time_gate_ratio_max"
        ],
        "turning_position_m": turn_pos,
        "turning_time_gate_ratio": abs(turn_time)
        / (GATE_TURNING_NS + float(np.spacing(np.float32(ref_t["time"])))),
        "turning_reached_sphere": int(bool(rows_t["present"][0])),
        "determinism_mismatch_count": sum(mismatch.values()),
    }
    return metrics, details


def turning_reference(columns: dict) -> dict:
    """Upward ray from the ground; centre = its descending crossing of the telescope plane; entry and time."""
    profile = uc.profile()
    ref = RayReference(profile, uc.SITE_M)
    origin = np.array(
        [columns["x"][0], columns["y"][0], columns["z"][0]], dtype=np.float64
    )
    d = unit([columns["cx"][0], columns["cy"][0], columns["cz"][0]])

    def plane(s, y):
        return float(y[2] - TELESCOPE_M[2])

    plane.direction = -1
    plane.terminal = False
    y0 = np.concatenate([origin, d, [0.0]])
    sol = solve_ivp(
        ref.rhs,
        (0.0, TURNING_MAX_PATH_M),
        y0,
        method="DOP853",
        rtol=REL_TOL,
        atol=ABS_TOL,
        events=[plane],
        max_step=50.0,
    )
    crossing = sol.y_events[0][0]
    centre = np.array([crossing[0], 0.0, TELESCOPE_M[2]])
    apex = float(np.max(sol.y[2]) - 0.0)
    entered, s, position, time = ref.trace(
        origin, d, 0.0, centre, RADIUS_M, TURNING_MAX_PATH_M
    )
    return {
        "centre": centre,
        "entered": entered,
        "position": position,
        "time": time,
        "apex_m": apex,
        "path_m": s,
    }


def main() -> None:
    common.results_main(run_all, __doc__.splitlines()[0])


if __name__ == "__main__":
    main()
