"""Scattered laser light: obdeect single scattering (step 6b) versus the analytic single-scattering
integral, and ls-beam recorded for comparison.

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.scattered_light

Geometry: a laser above the telescope fires a beam downwards through the 15 deg field of view
(laser at (-1000, 0, 3000) m, zenith 108.5 deg, azimuth 0, telescope zenith 30 deg, radius 12.5 m).
Downward only: obdeect scattering assumes primaries travel towards the observation level
(scattering.hpp:555 sets cz = -sqrt(1 - cx^2 - cy^2) for every primary, and a cz column that disagrees
raises; scattering.hpp:365 returns for primaries at or below the observation level; the scattering
documentation, scattering.md:85, says points below the telescope and photons travel downwards only).
The upward beam of the laser-beam case (zenith 60 deg) therefore gives no scattered photons at all.

Output: the propagation CLI passes the primary photons through unchanged (scattering.md, "The primary
photon is unchanged") and writes the scattered photons as additional rows. Only the scattered rows are
compared here: the Python scatter() API returns them alone, and the CLI weighted output minus the
output without scattering gives the same weight sum (cross-checked in the run).

Medium: Rayleigh only, tau_aer = 0, single-scattering albedo 1, 355 nm, vertical optical depth
tau(h) = 0.00865 (lambda/um)^-4.05 exp(-h/8 km) (the Rayleigh form of the obdeect scattering
documentation, not a fitted model), tabulated at 300 log-spaced altitudes so that the table
interpolation is negligible. The transmission table is zero (T = 1, no extinction), which obdeect
requires with scattering.

Analytic expectation (per unit path length s along the beam, direction d):

    g(s) = W_up(s) beta(h) p_R(theta) Omega(D)     for h > 0 and D > R, else 0
    beta(h) = -dtau/dh = tau(h) / 8000 m           (sum over the Rayleigh vertical depth)
    p_R(theta) = 3/(16 pi) (1 + cos^2 theta),      theta between d and the direction to the telescope
    Omega(D) = 2 pi (1 - sqrt(1 - (R/D)^2))        solid angle of the telescope sphere from distance D
    W_up(s) = sum over primaries upstream of s of w_k exp(-(tau(h_s) - tau(h_k))/|cz|)
    (the primary survives to s with the slant Rayleigh depth between its height h_k and h_s, as in
    scattering.hpp: exp(-step j) with D measured from the emission point)

Expected scattered weight = integral of g(s) ds. Its first two angular moments and the mean
altitude of the scattering points are compared with the same moments of the obdeect output.
The obdeect run uses 32 scattering points per primary (its discretisation error is 3e-4, see the
scattering documentation). ls-beam is run on the same geometry and recorded (yield, altitude of
the emission points); it is not judged: ls-beam scales a brightness per metre with exp(-z/7.5 km)
and has no optical depth or phase function, so its yield is not the same physical quantity.

Limits (set before the first run): yield ratio max(1 %, 5 standard errors of the scattered
weights); mean cos(theta) and mean cos^2(theta) within max(0.01, 5 se) (the obdeect aim point is
on the sphere disc, the analytic direction on its centre: the angle differs by up to R/D, about
0.01 rad here); mean altitude within max(1 m, 5 se).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common, readers

LAMBDA_NM = 355
SITE_ALTITUDE_M = 2156.0
TAU_0 = 0.00865 * (LAMBDA_NM / 1000.0) ** -4.05
SCALE_HEIGHT_M = 8000.0
TELESCOPE_RADIUS_M = 12.5
TELESCOPE_CENTRE_M = np.zeros(
    3
)  # telescope at the origin of the sky frame (laser-beam convention)
SCATTERING_POINTS = 32
STEP_M = 0.5
YIELD_FLOOR = 0.01
ANGLE_FLOOR = 0.01
ALTITUDE_FLOOR_M = 1.0
SWEEP = {
    "laser_position": "-1000,0,3000",
    "laser_theta": "108.5",
    "laser_phi": "0",
    "telescope_theta": "30",
    "telescope_phi": "0",
    "telescope_radius": str(TELESCOPE_RADIUS_M),
    "fov": "15",
    "primary_per_m": "1000",
    "step": "0.1",
    "spectrum": str(LAMBDA_NM),
    "altitude": str(SITE_ALTITUDE_M),
    "seed": "17",
}
UNIT_TRANSMISSION = common.VALIDATION / "fixtures" / "unit-transmission.dat"


def write_scattering_table(path):
    """Rayleigh-only scattering table (see module docstring); altitudes log-spaced 0.1-60 km."""
    altitudes = np.geomspace(0.1, 60.0, 300)
    rows = ["# wavelength_nm altitude_km tau_mol tau_aer ssa g"]
    for h in altitudes:
        tau = TAU_0 * math.exp(-h / 8.0)
        rows.append(f"{LAMBDA_NM} {h:.8f} {tau:.10e} 0.0 1.0 0.65")
    Path(path).write_text("\n".join(rows) + "\n", encoding="utf-8")
    return Path(path)


def write_zero_transmission(path):
    """Copy of the zero optical-depth fixture without the 1000 nm row (obdeect wants < 1000 nm)."""
    lines = []
    for line in UNIT_TRANSMISSION.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or (line.split() and int(line.split()[0]) < 1000):
            lines.append(line)
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return Path(path)


def analytic_expectation(primary_s, primary_w, direction, origin_point):
    """Expected scattered weight and moments from the path integral (see module docstring).

    primary_s: position of each primary along the beam line, primary_w: their weights,
    direction: unit beam direction, origin_point: any point on the beam line.
    Returns a dict with the yield and the weighted moments of the scattering points.
    """
    s_min = float(np.min(primary_s))
    # Ground (z = 0, observation level) reached by the beam: s where z = 0.
    if direction[2] >= 0:
        raise ValueError("beam must point downwards")
    s_ground = -origin_point[2] / direction[2]
    grid = np.arange(s_min, s_ground, STEP_M) + 0.5 * STEP_M
    points = origin_point[None, :] + grid[:, None] * direction[None, :]
    z = points[:, 2]
    dist = np.linalg.norm(points - TELESCOPE_CENTRE_M, axis=1)
    allowed = (z > 0.0) & (dist > TELESCOPE_RADIUS_M)
    cos_theta = (direction[None, :] * (TELESCOPE_CENTRE_M - points)).sum(axis=1) / dist
    phase = 3.0 / (16.0 * math.pi) * (1.0 + cos_theta**2)
    omega = (
        2.0
        * math.pi
        * (1.0 - np.sqrt(np.clip(1.0 - (TELESCOPE_RADIUS_M / dist) ** 2, 0, None)))
    )
    altitude_km = (SITE_ALTITUDE_M + z) / 1000.0
    beta = TAU_0 * np.exp(-altitude_km / 8.0) / SCALE_HEIGHT_M
    # Primary survival to the scattering point: exp(-(tau(h_k) - tau(h_s)) / |cz|), the slant
    # Rayleigh depth between the primary's own height h_k and the point (obdeect applies the same
    # factor through exp(-step j), with the table's extinction). Only primaries upstream of s count.
    cz = abs(direction[2])
    tau_grid = TAU_0 * np.exp(-((SITE_ALTITUDE_M + z) / 1000.0) / 8.0)
    order = np.argsort(primary_s)
    s_sorted = np.asarray(primary_s)[order]
    w_sorted = np.asarray(primary_w, dtype=float)[order]
    h_sorted = origin_point[2] + s_sorted * direction[2]
    tau_prim = TAU_0 * np.exp(-((SITE_ALTITUDE_M + h_sorted) / 1000.0) / 8.0)
    # tau decreases with altitude: the depth from the primary (higher) down to the point (lower) is
    # tau(h_s) - tau(h_k) >= 0, so the survival factor is exp(-tau_s/cz) exp(+tau_k/cz).
    weighted_surv = np.cumsum(w_sorted * np.exp(tau_prim / cz))
    idx = np.searchsorted(s_sorted, grid, side="right")
    has = idx > 0
    weight_up = np.zeros(grid.size)
    weight_up[has] = weighted_surv[idx[has] - 1] * np.exp(-tau_grid[has] / cz)
    g = np.where(allowed, weight_up * beta * phase * omega, 0.0)
    total = float(np.sum(g) * STEP_M)
    if total <= 0:
        raise ValueError("analytic expectation is zero")
    mean_cos = float(np.sum(g * cos_theta) / np.sum(g))
    mean_cos2 = float(np.sum(g * cos_theta**2) / np.sum(g))
    mean_z = float(np.sum(g * z) / np.sum(g))
    return {
        "yield": total,
        "mean_cos_theta": mean_cos,
        "mean_cos2_theta": mean_cos2,
        "mean_altitude_above_obs_m": mean_z,
        "path_length_m": float(s_ground - s_min),
    }


def weighted_moments(weights, cos_theta, z):
    """Weighted mean cos(theta), cos^2(theta) and z with their standard errors (weights as samples)."""
    w = np.asarray(weights, dtype=float)
    total = w.sum()
    out = {}
    for name, values in (
        ("cos_theta", cos_theta),
        ("cos2_theta", cos_theta**2),
        ("z", z),
    ):
        mean = float(np.sum(w * values) / total)
        var = float(np.sum(w * (values - mean) ** 2) / total)
        out[name] = (mean, float(math.sqrt(var * np.sum(w**2)) / total))
    return out


def _propagate(work, primary, out, transmission, table, mode, scattering):
    """obdeect propagation run on the primaries, with or without single scattering (CLI)."""
    argv = [
        common.tool("obdeect"),
        "propagation",
        "run",
        "--transmission",
        str(transmission),
        "--site-altitude",
        SWEEP["altitude"],
        "--wavelength-range",
        "300",
        "800",
        "--seed",
        "1",
        "--mode",
        mode,
        "-i",
        str(primary),
        "-o",
        str(out),
    ]
    if scattering:
        argv += [
            "--scattering",
            str(table),
            "--scattering-points",
            str(SCATTERING_POINTS),
        ]
    common.run_tool(argv, out.with_suffix(".log"))


def run(work):
    cfg_dir = common.REFERENCE / "cfg/common"
    table = write_scattering_table(work / "scattering.dat")
    transmission = write_zero_transmission(work / "transmission.dat")
    primary = work / "primary.gst"
    sweep = SWEEP
    common.run_tool(
        [
            common.tool("obdeect"),
            "sources",
            "run",
            "laser-beam",
            "--observation-level",
            sweep["altitude"],
            "--atmosphere-id",
            "1",
            "--events",
            "1",
            "--seed",
            sweep["seed"],
            "--bunch-size",
            "3",
            "--photons-per-m",
            sweep["primary_per_m"],
            "--spectrum",
            f"fixed:wavelength_nm={sweep['spectrum']}",
            "--laser-position",
            sweep["laser_position"],
            "--laser-zenith",
            sweep["laser_theta"],
            "--laser-azimuth",
            sweep["laser_phi"],
            "--telescope-zenith",
            sweep["telescope_theta"],
            "--telescope-azimuth",
            sweep["telescope_phi"],
            "--telescope-radius",
            sweep["telescope_radius"],
            "--fov-deg",
            sweep["fov"],
            "-o",
            str(primary),
        ],
        work / "primary.log",
    )
    # Weighted output with and without scattering: the primaries pass through in both, so the
    # difference of their weight sums is the scattered weight (CLI cross-check of the API run).
    base, with_scatter = work / "cli-base.gst", work / "cli-scattered.gst"
    _propagate(work, primary, base, transmission, table, "weighted", scattering=False)
    _propagate(
        work, primary, with_scatter, transmission, table, "weighted", scattering=True
    )
    ls_out = work / "ls-beam.iact.gz"
    ls_text = common.run_tool(
        [
            common.tool("ls-beam"),
            f"-I{cfg_dir}",
            "--events",
            "1",
            "--brightness",
            "1000",
            "--bunchsize",
            "1",
            "--step",
            sweep["step"],
            "--spectrum",
            sweep["spectrum"],
            "--laser-position",
            sweep["laser_position"],
            "--laser-theta",
            sweep["laser_theta"],
            "--laser-phi",
            sweep["laser_phi"],
            "--telescope-theta",
            sweep["telescope_theta"],
            "--telescope-phi",
            sweep["telescope_phi"],
            "--telescope-radius",
            sweep["telescope_radius"],
            "--fov",
            sweep["fov"],
            "--altitude",
            sweep["altitude"],
            "-o",
            str(ls_out),
        ],
        work / "ls-beam.log",
    )
    return primary, base, with_scatter, table, transmission, ls_out, ls_text


def _sky_columns(path):
    """Primary photons of a sky-frame stream as float32 columns (python API input)."""
    from obdeect.stream import Reader

    cols = {k: [] for k in ("x", "y", "z", "cx", "cy", "t", "w")}
    for record in Reader(str(path)):
        if record["type"] == "photons":
            data = dict(record["columns"])
            for key, values in cols.items():
                values.append(np.asarray(data[key], dtype=np.float32))
    out = {key: np.concatenate(values) for key, values in cols.items()}
    out["wavelength"] = np.full(out["x"].size, LAMBDA_NM, dtype=np.float32)
    return out


def _weight_sum(path):
    from obdeect.stream import Reader

    total = 0.0
    for record in Reader(str(path)):
        if record["type"] == "photons":
            total += float(
                np.sum(np.asarray(dict(record["columns"])["w"], dtype=np.float64))
            )
    return total


API_SEEDS = (1, 2, 3, 4)


def main():
    from obdeect.propagation import ScatteringTable, TransmissionTable, scatter

    work = Path.cwd().resolve()
    primary, base, with_scatter, table, transmission, ls_out, _ = run(work)
    cols = _sky_columns(primary)
    cx_mean, cy_mean = float(np.mean(cols["cx"])), float(np.mean(cols["cy"]))
    direction = np.array(
        [cx_mean, cy_mean, -math.sqrt(max(0.0, 1.0 - cx_mean**2 - cy_mean**2))]
    )
    direction /= np.linalg.norm(direction)
    positions = np.stack([cols["x"], cols["y"], cols["z"]], axis=1).astype(float)
    origin_point = positions[0]
    primary_s = (positions - origin_point) @ direction
    primary_w = cols["w"].astype(float)
    expected = analytic_expectation(primary_s, primary_w, direction, origin_point)

    # Python API (the function of the unit tests): scattered photons only, one run per seed.
    scat_table = ScatteringTable(str(table))
    site = TransmissionTable(
        str(transmission), site_altitude_m=float(SWEEP["altitude"])
    )
    telescope = [(1, 0.0, 0.0, 0.0, TELESCOPE_RADIUS_M)]
    yields, first = [], None
    for seed in API_SEEDS:
        out = scatter(cols, scat_table, site, telescope, seed, 1, 0, SCATTERING_POINTS)[
            1
        ]
        yields.append(float(np.sum(out["w"], dtype=np.float64)))
        if first is None:
            first = out
    yields = np.asarray(yields)
    yield_mc = float(yields.mean())
    yield_se_rel = float(yields.std(ddof=1) / math.sqrt(yields.size) / yield_mc)
    w = first["w"].astype(float)
    cos_beam = (
        first["cx"].astype(float) * direction[0]
        + first["cy"].astype(float) * direction[1]
        + first["cz"].astype(float) * direction[2]
    )
    z = first["zem"].astype(float) - SITE_ALTITUDE_M
    moments = weighted_moments(w, cos_beam, z)
    cli_yield = _weight_sum(with_scatter) - _weight_sum(base)

    ratio = yield_mc / expected["yield"]
    ratio_limit = max(YIELD_FLOOR, 5.0 * yield_se_rel)
    cos_mean, cos_se = moments["cos_theta"]
    cos2_mean, cos2_se = moments["cos2_theta"]
    z_mean, z_se = moments["z"]
    metrics = {
        "primary_photon_weight": float(primary_w.sum()),
        "primary_path_length_m": expected["path_length_m"],
        "scattered_weight_mc": yield_mc,
        "scattered_weight_mc_se_rel": yield_se_rel,
        "scattered_weight_cli_difference": cli_yield,
        "scattered_weight_analytic": expected["yield"],
        "scattered_yield_ratio": ratio,
        "scattered_yield_limit_exceeded": int(abs(ratio - 1.0) > ratio_limit),
        "cos_theta_mean_mc": cos_mean,
        "cos_theta_mean_analytic": expected["mean_cos_theta"],
        "cos_theta_mean_delta": cos_mean - expected["mean_cos_theta"],
        "cos_theta_mean_delta_se": cos_se,
        "cos_theta_mean_limit_exceeded": int(
            abs(cos_mean - expected["mean_cos_theta"]) > max(ANGLE_FLOOR, 5 * cos_se)
        ),
        "cos2_theta_mean_mc": cos2_mean,
        "cos2_theta_mean_analytic": expected["mean_cos2_theta"],
        "cos2_theta_mean_delta": cos2_mean - expected["mean_cos2_theta"],
        "cos2_theta_mean_limit_exceeded": int(
            abs(cos2_mean - expected["mean_cos2_theta"]) > max(ANGLE_FLOOR, 5 * cos2_se)
        ),
        "altitude_mean_mc_m": z_mean,
        "altitude_mean_analytic_m": expected["mean_altitude_above_obs_m"],
        "altitude_mean_delta_m": z_mean - expected["mean_altitude_above_obs_m"],
        "altitude_mean_limit_exceeded": int(
            abs(z_mean - expected["mean_altitude_above_obs_m"])
            > max(ALTITUDE_FLOOR_M, 5 * z_se)
        ),
    }
    ls_counts, ls = readers.read_iact(ls_out)
    details = {
        "case": "laser beam scattered light: obdeect single scattering (python scatter, CLI cross-check) versus analytic integral; ls-beam recorded",
        "geometry": SWEEP,
        "medium": f"Rayleigh only, {LAMBDA_NM} nm, tau_0 = {TAU_0:.6f} at the ground, scale height 8 km, "
        "tau_aer = 0, ssa = 1, transmission table zero (T = 1)",
        "scattering_points_per_primary": SCATTERING_POINTS,
        "api_seeds": list(API_SEEDS),
        "api_yields": yields.tolist(),
        "analytic_expectation": expected,
        "output_definition": "scattered photons = rows of scatter() (Python API); CLI weighted run with minus "
        "without scattering gives the same weight sum; the primaries are passed through by propagation run "
        "(they are not scattered photons) and are excluded",
        "ls_beam_recorded": {
            "photon_weight": float(np.sum(ls["w"])) if ls["w"].size else 0.0,
            "events": len(ls_counts),
            "judged": False,
            "reason": "ls-beam scales a brightness per metre with exp(-z/7.5 km), no optical depth or phase "
            "function: a different physical model, not compared quantitatively",
        },
        "limits": {
            "yield": f"max({YIELD_FLOOR}, 5 se), se = std of 4 seeds / sqrt(4)",
            "angles": f"max({ANGLE_FLOOR}, 5 se) on the mean cos(theta) and cos^2(theta)",
            "altitude": f"max({ALTITUDE_FLOOR_M} m, 5 se) on the mean altitude",
        },
    }
    common.write_results(metrics, details)
    from obdeect_validation.adapters import plotting

    fig, axes = plotting.figure(1, 2, (10, 3.6))
    axes[0, 0].hist(z, weights=w, bins=40, color=plotting.NATIVE, label="obdeect")
    axes[0, 0].set(
        xlabel="scattering altitude above observation level [m]", ylabel="weight"
    )
    axes[0, 0].legend()
    axes[0, 1].hist(
        cos_beam,
        weights=w,
        bins=40,
        color=plotting.NATIVE,
        label="obdeect, angle to beam",
    )
    axes[0, 1].set(xlabel="cos(scattering angle)", ylabel="weight")
    axes[0, 1].legend()
    plotting.save(fig, Path("scattered-light.png"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
