"""Sources with upward directions: isotropic light moved to the plane z = 0, plane waves, laser beams.

Three parts. Every reference is analytic and written here (numpy, scipy), not taken from obdeect.

1. Isotropic point source at height H above the plane z = 0 (and at -H below it), N photons,
   ``to_plane_z0 = true``. Photons that move away from the plane are dropped and counted as
   ``not_on_plane``. Checks: kept fraction 1/2 (binomial pull); kept + not_on_plane = N exactly,
   also in the run summary and the ``run_end`` record; the radial density of the kept photons
   on the plane, dN/dA = (N / 4 pi) H / (r^2 + H^2)^(3/2), in equal-probability radial bins
   (chi2/ndf); arrival times t = sqrt(r^2 + H^2) / c to float32 rounding.
2. Upward plane wave (``towards_source.z < 0``, source below the disc): every photon has the
   requested direction, starts on the start plane z = start_z, and lies on the target disc
   (uniform in area). Rejected: a horizontal ``towards_source`` and a start plane on the wrong
   side of the disc, each with the message of obdeect.
3. Upward laser beam (``laser-beam`` preset): exact direction, start on the beam at the start
   height above the laser, arrival time offset path / c, cone bound for a divergent beam, tilt
   limit measured from the zenith. Rejected: a horizontal beam, a tilt above the limit, and a
   start height below the laser.

Radial geometry: a photon from (0, 0, H) reaches the plane at radius r with cos(theta) = H / d,
d = sqrt(r^2 + H^2). The fraction of the emitted photons in the annulus r1..r2 is
(c(r1) - c(r2)) / 2 with c(r) = H / sqrt(r^2 + H^2) (solid angle of the cone, hemisphere toward
the plane). The kept photons have P(r < R | kept) = 1 - c(R).
"""

from __future__ import annotations

import argparse
import importlib.metadata
import math
from pathlib import Path

import numpy as np
from obdeect.sources import run
from obdeect.stream import Reader
from scipy import stats

from obdeect_validation.adapters import common, plotting

C_M_NS = 0.299792458
SEED = 20261123
PHOTONS = 200_000
HEIGHT_M = 100.0
RADIAL_BINS = 20  # equal-probability bins of the kept photons
PLOT_EDGES_M = np.logspace(0.0, 7.0, 29)  # log-spaced for the density plot
DISC_RADIUS_M = 20.0
PLANE_PHOTONS = 10_000
PLANE_TOWARDS = (
    0.36,
    -0.48,
    -0.8,
)  # unit vector; photons travel along -towards: upward
PLANE_START_M = -50.0  # below the disc, on the source side
PLANE_REJECT_START_M = 50.0  # above the disc: wrong side for a source below
LASER_OBSERVATION_M = 2150.0
LASER_POSITION_M = (0.0, 0.0, 1000.0)  # sky frame: z above the observation level
LASER_ZENITH_DEG = 20.0
LASER_AZIMUTH_DEG = 30.0
LASER_START_M = 100000.0  # sky-frame height of the start on the upward beam
LASER_ROWS = 1000
LASER_CONE_FULL_DEG = 0.5
LASER_MAX_TILT_DEG = 20.0
LASER_TILT_OVER_DEG = 25.0
TIME_KEYS = ("x", "y", "z", "cx", "cy", "cz", "t")


# ------------------------------------------------------------------ reference formulas


def hemisphere_fraction(height: float, r1: float, r2: float) -> float:
    """Fraction of the emitted photons of a point at height H that reach the plane in r1..r2."""
    h = abs(height)

    def cosine(r: float) -> float:
        return h / math.sqrt(r * r + h * h) if math.isfinite(r) else 0.0

    return 0.5 * (cosine(r1) - cosine(r2))


def radial_density(height: float, n_emitted: int, radius: np.ndarray) -> np.ndarray:
    """dN/dA on the plane [1/m^2]: (N / 4 pi) |H| / (r^2 + H^2)^(3/2)."""
    h = abs(height)
    return n_emitted / (4.0 * np.pi) * h / (radius * radius + h * h) ** 1.5


def radial_edges(height: float, bins: int = RADIAL_BINS) -> np.ndarray:
    """Equal-probability radial edges of the kept photons (last edge infinite).

    With P(r < R | kept) = 1 - c(R), the edge k/bins solves c = 1 - k/bins.
    """
    h = abs(height)
    edges = np.empty(bins + 1)
    edges[0] = 0.0
    for k in range(1, bins + 1):
        c = 1.0 - k / bins
        edges[k] = h * math.sqrt(1.0 - c * c) / c if c > 0.0 else math.inf
    return edges


def arrival_time_reference(x, y, height: float) -> np.ndarray:
    """Arrival time at the plane of a photon emitted at (0, 0, height) at t = 0 [ns]."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    return np.sqrt(x * x + y * y + height * height) / C_M_NS


def angle_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Angle between rows of a and b [deg]; atan2 keeps tiny angles accurate."""
    cross = np.linalg.norm(np.cross(a, b), axis=-1)
    return np.degrees(np.arctan2(cross, np.sum(a * b, axis=-1)))


def beam_direction(zenith_deg: float, azimuth_deg: float) -> np.ndarray:
    """Direction of travel from zenith and astronomical azimuth (x north, y west, z up)."""
    zenith, azimuth = math.radians(zenith_deg), math.radians(azimuth_deg)
    return np.array(
        [
            math.sin(zenith) * math.cos(azimuth),
            -math.sin(zenith) * math.sin(azimuth),
            math.cos(zenith),
        ]
    )


# ------------------------------------------------------------------ stream reading


def read_photons(path) -> dict:
    """Photon columns (float32, as stored) of a stream, plus the run_end not_on_plane count."""
    columns = {name: [] for name in TIME_KEYS}
    run_end = None
    for record in Reader(str(path)):
        if record["type"] == "photons":
            for name in TIME_KEYS:
                columns[name].append(np.asarray(record["columns"][name]))
        elif record["type"] == "run_end":
            run_end = record["metadata"]
    data = {
        name: np.concatenate(parts) if parts else np.empty(0, dtype=np.float32)
        for name, parts in columns.items()
    }
    data["run_end_not_on_plane"] = (
        int(run_end.get("not_on_plane", -1)) if run_end else -1
    )
    return data


# ------------------------------------------------------------------ part 1: isotropic to plane


def isotropic_run(height: float, photons: int, seed: int) -> dict:
    """One isotropic point emission at (0, 0, height), moved to the plane z = 0."""
    return {
        "run": {"seed": seed, "observation_level_m": 0.0, "frame": "sky"},
        "telescopes": [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": 10.0}],
        "events": [
            {
                "repeat": 1,
                "emissions": [
                    {
                        "kind": "point",
                        "photons": float(photons),
                        "position_m": [0.0, 0.0, height],
                        "angular": {"kind": "isotropic"},
                        "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
                        "to_plane_z0": True,
                        "target": "none",
                    }
                ],
            }
        ],
    }


def radial_test(r: np.ndarray, height: float) -> dict:
    """Chi2 of the kept photons in equal-probability radial bins (expected: N_kept / bins each)."""
    edges = radial_edges(height)
    finite = np.where(np.isinf(edges), 1e30, edges)  # numpy needs finite edges
    observed, _ = np.histogram(r, bins=finite)
    expected = observed.sum() / RADIAL_BINS
    chi2 = float(np.sum((observed - expected) ** 2 / expected))
    ndf = RADIAL_BINS - 1
    return {
        "chi2": chi2,
        "ndf": ndf,
        "chi2_ndf": chi2 / ndf,
        "p_value": float(stats.chi2.sf(chi2, ndf)),
        "counts": observed.tolist(),
        "expected_per_bin": float(expected),
    }


def time_test(data: dict, height: float) -> tuple[float, float]:
    """Largest |t - reference| [ns] and the same in float32 units of the stored t (ulp)."""
    t = data["t"]
    error = np.abs(
        t.astype(np.float64) - arrival_time_reference(data["x"], data["y"], height)
    )
    ulp = np.spacing(t).astype(np.float64)
    return float(error.max()), float((error / ulp).max())


def isotropic_part(root: Path, photons: int = PHOTONS, seed: int = SEED) -> dict:
    """Both heights: kept fraction, photon balance, radial density and arrival times."""
    cases = {}
    for label, height in (("above", HEIGHT_M), ("below", -HEIGHT_M)):
        path = root / f"isotropic-{label}.pst"
        summary = run.write(isotropic_run(height, photons, seed), path)
        data = read_photons(path)
        kept = len(data["x"])
        r = np.hypot(data["x"].astype(np.float64), data["y"].astype(np.float64))
        time_ns, time_ulp = time_test(data, height)
        cases[label] = {
            "height_m": height,
            "emitted": photons,
            "kept_rows": kept,
            "summary": summary,
            "not_on_plane_run_end": data["run_end_not_on_plane"],
            "kept_fraction": kept / photons,
            "kept_fraction_pull": (kept - 0.5 * photons) / math.sqrt(0.25 * photons),
            "balance": kept + summary["not_on_plane"] - photons,
            "run_end_delta": data["run_end_not_on_plane"] - summary["not_on_plane"],
            "rows_delta": kept - summary["rows"],
            "weight_delta": summary["weight"] - kept,
            "on_plane_z_max_abs_m": float(np.abs(data["z"]).max()) if kept else 0.0,
            "radial": radial_test(r, height),
            "max_time_error_ns": time_ns,
            "max_time_error_ulp": time_ulp,
            "_r": r,
        }
    return cases


def density_comparison(r: np.ndarray, height: float, photons: int) -> dict:
    """Plot bins: obdeect counts per area against the analytic annulus fraction (pulls)."""
    edges = PLOT_EDGES_M
    observed, _ = np.histogram(r, bins=edges)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    expected = np.array(
        [
            photons * hemisphere_fraction(height, edges[i], edges[i + 1])
            for i in range(len(edges) - 1)
        ]
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        pulls = np.where(expected > 0, (observed - expected) / np.sqrt(expected), 0.0)
    return {
        "edges_m": edges,
        "area_m2": area,
        "counts": observed,
        "expected": expected,
        "density": observed / (photons * area),
        "pulls": pulls,
    }


# ------------------------------------------------------------------ part 2: upward plane wave


def plane_wave_run(
    start_z: float, towards, seed: int, photons: int = PLANE_PHOTONS
) -> dict:
    """Plane wave through a disc of radius DISC_RADIUS_M at z = 0, from the start plane z = start_z."""
    return {
        "run": {"seed": seed, "observation_level_m": 0.0, "frame": "sky"},
        "telescopes": [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": 10.0}],
        "events": [
            {
                "repeat": 1,
                "emissions": [
                    {
                        "kind": "plane_wave",
                        "photons": float(photons),
                        "towards_source": list(towards),
                        "disc_radius_m": DISC_RADIUS_M,
                        "disc_z_m": 0.0,
                        "start_z_m": start_z,
                        "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
                    }
                ],
            }
        ],
    }


def plane_wave_part(root: Path, seed: int = SEED, photons: int = PLANE_PHOTONS) -> dict:
    """Upward plane wave: direction, start plane, and hits of the target disc."""
    path = root / "plane-wave-upward.pst"
    summary = run.write(
        plane_wave_run(PLANE_START_M, PLANE_TOWARDS, seed, photons), path
    )
    data = read_photons(path)
    stored = np.stack([data[k].astype(np.float64) for k in ("cx", "cy", "cz")], axis=1)
    expected = -np.array(PLANE_TOWARDS, dtype=np.float64)  # direction of travel
    direction_error = angle_deg(stored, np.broadcast_to(expected, stored.shape))
    x = data["x"].astype(np.float64)
    y = data["y"].astype(np.float64)
    z = data["z"].astype(np.float64)
    cz = stored[:, 2]
    start_error = float(np.abs(z - PLANE_START_M).max())
    # Along the stored direction back to z = 0 (the disc); the photon starts at z = PLANE_START_M.
    steps = -z / cz
    tx, ty = x + steps * stored[:, 0], y + steps * stored[:, 1]
    tz = z + steps * cz
    radius = np.hypot(tx, ty)
    r2 = (radius / DISC_RADIUS_M) ** 2
    # Area-uniform disc: r^2 / R^2 is uniform on [0, 1], variance 1/12.
    pull = (float(r2.mean()) - 0.5) / math.sqrt(1.0 / 12.0 / len(r2))
    return {
        "photons": photons,
        "rows": len(x),
        "summary": summary,
        "direction_error_deg": float(direction_error.max()),
        "start_error_m": start_error,
        "target_z_max_abs_m": float(np.abs(tz).max()),
        "target_outside_count": int(np.count_nonzero(radius > DISC_RADIUS_M + 1e-4)),
        "target_r2_pull": float(pull),
        "targets_xy": np.stack([tx, ty], axis=1),
    }


def rejection(run_description, expected_text: str) -> dict:
    """Run the description; a ValueError must name expected_text. Other outcomes are not a pass."""
    try:
        run.write(run_description, "/dev/null")
    except ValueError as exc:
        message = str(exc)
        return {
            "rejected": True,
            "message": message,
            "message_found": expected_text in message,
        }
    return {"rejected": False, "message": "", "message_found": False}


def plane_wave_rejections() -> dict:
    horizontal = plane_wave_run(PLANE_START_M, (1.0, 0.0, 0.0), SEED, 10)
    wrong_side = plane_wave_run(PLANE_REJECT_START_M, PLANE_TOWARDS, SEED, 10)
    return {
        "horizontal_towards_source": rejection(horizontal, "towards_source (z = 0)"),
        "start_plane_wrong_side": rejection(wrong_side, "source side of the disc"),
    }


# ------------------------------------------------------------------ part 3: upward laser beam


def laser_argv(
    out: Path, zenith: float, azimuth: float, divergence: float, **extra
) -> list:
    """obdeect sources run laser-beam with the upward-beam options (rows, energy fixed)."""
    argv = [
        "obdeect",
        "sources",
        "run",
        "laser-beam",
        "--observation-level",
        str(LASER_OBSERVATION_M),
        "--laser-position",
        ",".join(str(v) for v in LASER_POSITION_M),
        "--beam-zenith",
        str(zenith),
        "--beam-azimuth",
        str(azimuth),
        "--divergence",
        str(divergence),
        "--pulse-energy",
        "1e-3",
        "--wavelength-nm",
        "355",
        "--telescope-radius",
        "12.5",
        "--rows",
        str(LASER_ROWS),
    ]
    for key, value in extra.items():
        argv += [f"--{key.replace('_', '-')}", str(value)]
    return argv + ["-o", str(out)]


def run_laser(root: Path, name: str, argv: list) -> dict:
    """Run one laser command; returns its outcome (ok, or the error text of the command)."""
    log = root / f"{name}.log"
    try:
        common.command(argv, log)
    except RuntimeError as exc:
        return {"ok": False, "message": str(exc)}
    return {"ok": True, "message": ""}


def laser_part(root: Path) -> dict:
    """Upward beam, zero and finite divergence; tilt limit and rejections."""
    out = root / "laser-upward.pst"
    outcome = run_laser(
        root,
        "laser-upward",
        laser_argv(
            out, LASER_ZENITH_DEG, LASER_AZIMUTH_DEG, 0.0, start_height=LASER_START_M
        ),
    )
    if not outcome["ok"]:
        raise RuntimeError(f"upward laser beam rejected: {outcome['message']}")
    data = read_photons(out)
    direction = beam_direction(LASER_ZENITH_DEG, LASER_AZIMUTH_DEG)
    laser = np.array(LASER_POSITION_M, dtype=np.float64)
    path = (LASER_START_M - laser[2]) / direction[2]
    start = laser + direction * path
    stored = np.stack([data[k].astype(np.float64) for k in ("cx", "cy", "cz")], axis=1)
    positions = np.stack([data[k].astype(np.float64) for k in ("x", "y", "z")], axis=1)
    time_ref = path / C_M_NS
    time_error = np.abs(data["t"].astype(np.float64) - time_ref)
    time_ulp = float((time_error / np.spacing(data["t"]).astype(np.float64)).max())

    cone_out = root / "laser-cone.pst"
    cone = run_laser(
        root,
        "laser-cone",
        laser_argv(cone_out, LASER_ZENITH_DEG, LASER_AZIMUTH_DEG, LASER_CONE_FULL_DEG),
    )
    cone_data = read_photons(cone_out) if cone["ok"] else None
    cone_angle = cone_start = float("nan")
    if cone_data is not None:
        cone_stored = np.stack(
            [cone_data[k].astype(np.float64) for k in ("cx", "cy", "cz")], axis=1
        )
        cone_angle = float(
            angle_deg(cone_stored, np.broadcast_to(direction, cone_stored.shape)).max()
        )
        cone_positions = np.stack(
            [cone_data[k].astype(np.float64) for k in ("x", "y", "z")], axis=1
        )
        cone_start = float(np.abs(cone_positions - laser).max())

    tilt_ok = run_laser(
        root,
        "laser-tilt-at-limit",
        laser_argv(
            root / "laser-tilt-at-limit.pst",
            LASER_ZENITH_DEG,
            LASER_AZIMUTH_DEG,
            0.0,
            max_tilt=LASER_MAX_TILT_DEG,
        ),
    )
    tilt_over = run_laser(
        root,
        "laser-tilt-over",
        laser_argv(
            root / "laser-tilt-over.pst",
            LASER_TILT_OVER_DEG,
            LASER_AZIMUTH_DEG,
            0.0,
            max_tilt=LASER_MAX_TILT_DEG,
        ),
    )
    horizontal = run_laser(
        root,
        "laser-horizontal",
        laser_argv(root / "laser-horizontal.pst", 90.0, LASER_AZIMUTH_DEG, 0.0),
    )
    below = run_laser(
        root,
        "laser-start-below",
        laser_argv(
            root / "laser-start-below.pst",
            LASER_ZENITH_DEG,
            LASER_AZIMUTH_DEG,
            0.0,
            start_height=500.0,
        ),
    )
    return {
        "reference_direction": direction.tolist(),
        "reference_start_m": start.tolist(),
        "reference_path_m": path,
        "rows": len(stored),
        "direction_error_deg": float(
            angle_deg(stored, np.broadcast_to(direction, stored.shape)).max()
        ),
        "start_error_m": float(np.abs(positions - start).max()),
        "start_height_error_m": float(np.abs(positions[:, 2] - LASER_START_M).max()),
        "time_error_ulp": time_ulp,
        "time_error_ns": float(time_error.max()),
        "cone_full_deg": LASER_CONE_FULL_DEG,
        "cone_max_angle_deg": cone_angle,
        "cone_excess_deg": max(0.0, cone_angle - 0.5 * LASER_CONE_FULL_DEG),
        "cone_start_error_m": cone_start,
        "cone_ok": cone["ok"],
        "tilt_at_limit_accepted": tilt_ok["ok"],
        "tilt_over_rejected": (not tilt_over["ok"])
        and "exceeds --max-tilt" in tilt_over["message"],
        "horizontal_rejected": (not horizontal["ok"])
        and "horizontal" in horizontal["message"],
        "start_below_rejected": (not below["ok"])
        and "must lie on the beam" in below["message"],
        "messages": {
            "tilt_over": tilt_over["message"],
            "horizontal": horizontal["message"],
            "start_below": below["message"],
        },
    }


# ------------------------------------------------------------------ plots


def plot_radial(path: Path, cases: dict, photons: int) -> None:
    """Radial density (top, log-log) and pulls per bin (bottom) of both heights."""
    fig, axes = plotting.figure(2, 1, (7.0, 6.2))
    top, bottom = axes[0, 0], axes[1, 0]
    colors = {"above": plotting.NATIVE, "below": plotting.CYCLE[2]}
    centres = np.sqrt(PLOT_EDGES_M[1:] * PLOT_EDGES_M[:-1])
    fine = np.logspace(0.0, 7.0, 400)
    for label, case in cases.items():
        comparison = density_comparison(case["_r"], case["height_m"], photons)
        top.plot(
            centres,
            comparison["density"],
            "o",
            markersize=4,
            markerfacecolor="none",
            color=colors[label],
            label=f"obdeect, source at z = {case['height_m']:+g} m",
        )
        bottom.plot(
            centres,
            comparison["pulls"],
            "o",
            markersize=3.5,
            color=colors[label],
            label=f"z = {case['height_m']:+g} m",
        )
    analytic = radial_density(
        HEIGHT_M, 1, fine
    )  # per emitted photon, as the density above
    top.plot(
        fine,
        analytic,
        "--",
        color=plotting.ANALYTIC,
        label="analytic (1/4pi) H/(r²+H²)^1.5",
    )
    top.set_xscale("log")
    top.set_yscale("log")
    top.set(ylabel="photons per m² per emitted photon")
    top.legend(loc="upper right", fontsize="small")
    bottom.axhline(0.0, color=plotting.ANALYTIC, linewidth=0.8)
    bottom.axhspan(-2.0, 2.0, color=plotting.LIMIT, alpha=0.12, linewidth=0)
    bottom.set_xscale("log")
    bottom.set(
        xlabel="radius on the plane r [m]", ylabel="pull (obdeect − analytic)/√analytic"
    )
    bottom.legend(loc="upper right", fontsize="small")
    plotting.save(fig, path)


def plot_plane_wave(path: Path, targets: np.ndarray) -> None:
    fig, axis = plotting.figure(1, 1, (4.8, 4.6))
    theta = np.linspace(0.0, 2.0 * np.pi, 200)
    axis.plot(
        DISC_RADIUS_M * np.cos(theta),
        DISC_RADIUS_M * np.sin(theta),
        "--",
        color=plotting.ANALYTIC,
        label="target disc (analytic)",
    )
    axis.plot(
        targets[:, 0],
        targets[:, 1],
        ".",
        markersize=2.0,
        color=plotting.NATIVE,
        label="obdeect targets (z = 0)",
    )
    axis.set_aspect("equal")
    axis.set(xlabel="x [m]", ylabel="y [m]")
    axis.legend(loc="upper right", fontsize="small")
    plotting.save(fig, path)


# ------------------------------------------------------------------ main


def version() -> str:
    try:
        return importlib.metadata.version("obdeect-sources")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def run_all(
    root: Path, photons: int = PHOTONS, plane_photons: int = PLANE_PHOTONS
) -> tuple[dict, dict, dict]:
    """All parts; returns (metrics, details, arrays for the plots). main() writes the outputs."""
    isotropic = isotropic_part(root, photons)
    plane = plane_wave_part(root, SEED, plane_photons)
    rejections = plane_wave_rejections()
    laser = laser_part(root)
    rejected_plane = sum(
        1 for case in rejections.values() if case["rejected"] and case["message_found"]
    )
    laser_rejected = sum(
        1
        for key in ("horizontal_rejected", "start_below_rejected", "tilt_over_rejected")
        if laser[key]
    )
    below, above = isotropic["below"], isotropic["above"]
    metrics = {
        "kept_fraction_pull_max_abs": max(
            abs(above["kept_fraction_pull"]), abs(below["kept_fraction_pull"])
        ),
        "photon_balance_error": (
            abs(above["balance"])
            + abs(below["balance"])
            + abs(above["run_end_delta"])
            + abs(below["run_end_delta"])
            + abs(above["rows_delta"])
            + abs(below["rows_delta"])
            + abs(above["weight_delta"])
            + abs(below["weight_delta"])
        ),
        "on_plane_z_max_abs_m": max(
            above["on_plane_z_max_abs_m"], below["on_plane_z_max_abs_m"]
        ),
        "radial_chi2_ndf": max(
            above["radial"]["chi2_ndf"], below["radial"]["chi2_ndf"]
        ),
        "max_time_error_ulp": max(
            above["max_time_error_ulp"], below["max_time_error_ulp"]
        ),
        "plane_wave_direction_error_deg": plane["direction_error_deg"],
        "plane_wave_start_error_m": plane["start_error_m"],
        "plane_wave_target_outside_count": plane["target_outside_count"],
        "plane_wave_target_r2_pull": abs(plane["target_r2_pull"]),
        "plane_wave_rows_delta": abs(plane["rows"] - PLANE_PHOTONS),
        "rejected_cases": rejected_plane,
        "laser_direction_error_deg": laser["direction_error_deg"],
        "laser_start_error_m": laser["start_error_m"],
        "laser_start_height_error_m": laser["start_height_error_m"],
        "laser_time_error_ulp": laser["time_error_ulp"],
        "laser_cone_excess_deg": laser["cone_excess_deg"],
        "laser_tilt_limit_ok": int(laser["tilt_at_limit_accepted"])
        + int(laser["tilt_over_rejected"]),
        "laser_rejected_cases": laser_rejected,
    }
    details = {
        "seed": SEED,
        "obdeect_sources_version": version(),
        "isotropic": {
            label: {
                key: value for key, value in case.items() if not key.startswith("_")
            }
            for label, case in isotropic.items()
        },
        "radial_bins": RADIAL_BINS,
        "plane_wave": {k: v for k, v in plane.items() if k != "targets_xy"},
        "plane_wave_rejections": rejections,
        "laser": laser,
        "references": {
            "kept_fraction": 0.5,
            "radial": "(N/4pi) H/(r^2+H^2)^1.5; equal-probability bins from 1 - H/sqrt(R^2+H^2)",
            "time": "sqrt(x^2+y^2+H^2)/c, c = 0.299792458 m/ns, float32 ulp of the stored t",
            "plane_wave_disc": "area-uniform, r^2/R^2 uniform on [0,1]",
            "laser_direction": "sin(zenith)cos(az), -sin(zenith)sin(az), cos(zenith) (x north, y west)",
        },
    }
    arrays = {
        "r_above": isotropic["above"]["_r"],
        "r_below": isotropic["below"]["_r"],
        "targets": plane["targets_xy"],
    }
    return metrics, details, arrays


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args()
    root = Path.cwd()
    metrics, details, arrays = run_all(root)
    cases = {
        "above": {"height_m": HEIGHT_M, "_r": arrays["r_above"]},
        "below": {"height_m": -HEIGHT_M, "_r": arrays["r_below"]},
    }
    plot_radial(root / "isotropic-radial.png", cases, PHOTONS)
    plot_plane_wave(root / "plane-wave-targets.png", arrays["targets"])
    common.write_results(metrics, details, root)


if __name__ == "__main__":
    main()
