"""Extinction of upward photons (plane-parallel and spherical), against independent references.

Case propagation.upward-extinction. Upward photons (cz > 0) below the telescope, weighted output
(w = T), telescope at 3 km above the observation level (plane z = 5 km a.s.l.), wavelengths random in
[400, 600] nm (given, not drawn). Checks:

1. Plane-parallel: T = exp(-|tau(lambda, h_t) - tau(lambda, h_e)| / |cz|), tau from the table
   (upward_common.Transmission: log10-linear in height, linear in wavelength, zero at H2).
2. Spherical: T = exp(-tau_slant), tau_slant = integral over [h_e, h_t] of |dtau/dh| / cos(theta(h)),
   with cos(theta(h)) = sqrt(1 - b^2 / (R + h)^2) for the straight line of impact parameter
   b = |(P - C) x d| about the Earth centre C (R = 6371315 m, the CORSIKA radius). The emission
   altitude is the spherical one, sqrt(rho^2 + (R + H + z)^2) - R. The integral is quad with the table
   heights as breakpoints. The flat limit (R = 1e12 m) must reproduce the plane-parallel value.
3. Downward rows are unchanged: telescope 2 (downward photons above the plane) is bit-identical in the
   mixed run and in a downward-only run, in plane and spherical modes.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import itertools
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from obdeect_validation.adapters import (
    common,
    plotting,
)
from obdeect_validation.adapters import (
    endtoend_common as e2e,
)
from obdeect_validation.adapters import (
    upward_common as uc,
)

SEED = 20261125
N_UP = 300
N_DOWN = 60
TELESCOPE_Z_M = 3000.0
EARTH_RADIUS_M = 6371315.0
FLAT_RADIUS_M = 1.0e12
WAVELENGTH_RANGE_NM = (400.0, 600.0)
MAX_ZENITH_DEG = 60.0  # plane-parallel and spherical samples
EMISSION_MAX_M = 2990.0
C_SITE_M = uc.SITE_M


def version() -> str:
    try:
        return importlib.metadata.version("obdeect-propagation")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def upward_sample(rng, n: int) -> dict:
    """Upward photons: random emission (x, y, z below the plane), zenith up to MAX_ZENITH_DEG."""
    zenith = rng.uniform(0.0, MAX_ZENITH_DEG, n)
    azimuth = rng.uniform(0.0, 360.0, n)
    th = np.radians(zenith)
    cx = np.sin(th) * np.cos(np.radians(azimuth))
    cy = np.sin(th) * np.sin(np.radians(azimuth))
    cz = np.cos(th)
    cols = {
        "x": rng.uniform(-400.0, 400.0, n),
        "y": rng.uniform(-400.0, 400.0, n),
        "z": rng.uniform(0.0, EMISSION_MAX_M, n),
        "cx": cx,
        "cy": cy,
        "cz": cz,
        "t": rng.uniform(0.0, 10.0, n),
        "wavelength": rng.uniform(*WAVELENGTH_RANGE_NM, n),
    }
    return {k: np.asarray(v, dtype=np.float32) for k, v in cols.items()}


def downward_sample(rng, n: int) -> dict:
    """Downward photons above the plane (z above TELESCOPE_Z_M) at moderate zenith."""
    zenith = rng.uniform(0.0, MAX_ZENITH_DEG, n)
    azimuth = rng.uniform(0.0, 360.0, n)
    th = np.radians(zenith)
    cols = {
        "x": rng.uniform(-400.0, 400.0, n),
        "y": rng.uniform(-400.0, 400.0, n),
        "z": rng.uniform(TELESCOPE_Z_M + 10.0, TELESCOPE_Z_M + 1500.0, n),
        "cx": np.sin(th) * np.cos(np.radians(azimuth)),
        "cy": np.sin(th) * np.sin(np.radians(azimuth)),
        "cz": -np.cos(th),
        "t": rng.uniform(0.0, 10.0, n),
        "wavelength": rng.uniform(*WAVELENGTH_RANGE_NM, n),
    }
    return {k: np.asarray(v, dtype=np.float32) for k, v in cols.items()}


def plane_reference(tab: uc.Transmission, cols: dict) -> np.ndarray:
    """exp(-|tau(h_t) - tau(h_e)| / |cz|) per photon (float64 from the stored values)."""
    h_t = uc.SITE_M + TELESCOPE_Z_M
    out = np.empty(len(cols["z"]))
    for i in range(len(out)):
        lam = float(cols["wavelength"][i])
        h_e = uc.SITE_M + float(cols["z"][i])
        cz = abs(float(cols["cz"][i]))
        depth = abs(float(tab.tau(lam, h_t)) - float(tab.tau(lam, h_e)))
        out[i] = math.exp(-depth / cz)
    return out


def _slant_depth(
    tab: uc.Transmission, lam: float, p: float, h_e: float, h_t: float, radius_m: float
) -> float:
    """Integral of |dtau/dh| / cos(theta(h)) over [h_e, h_t] with cos(theta) = sqrt(1 - (p / (R + h))^2)."""

    def integrand(h: float) -> float:
        cos_h = math.sqrt(max(1.0 - (p / (radius_m + h)) ** 2, 0.0))
        return abs(tab.dtau_dh(lam, h)) / cos_h

    edges = [h_e, *tab.breakpoints(h_e, h_t), h_t]
    return sum(
        quad(integrand, a, c, limit=400, epsabs=1e-13, epsrel=1e-11)[0]
        for a, c in itertools.pairwise(edges)
    )


def spherical_reference(
    tab: uc.Transmission, cols: dict, radius_m: float, convention: str = "telescope"
) -> np.ndarray:
    """exp(-tau_slant) for upward photons: slant depth over [h_e, h_t] along the straight line.

    cos(theta(h)) = sqrt(1 - p^2 / (R + h)^2) with p = (R + h_t) sin(theta_t) (convention "telescope",
    documented in docs/propagation/curvature.md: theta_t from the direction against the telescope's
    local vertical, sin(theta_t) = sqrt(1 - cz^2)). Convention "line": p = |(P - C) x d| about the
    Earth centre (the impact parameter of the actual line; differs by about offset / R). The emission
    altitude is the spherical one for both.
    """
    r_obs = radius_m + uc.SITE_M
    h_t = uc.SITE_M + TELESCOPE_Z_M
    out = np.empty(len(cols["z"]))
    for i in range(len(out)):
        lam = float(cols["wavelength"][i])
        x0, y0, z0 = (float(cols[k][i]) for k in ("x", "y", "z"))
        cx, cy, cz = (float(cols[k][i]) for k in ("cx", "cy", "cz"))
        p_vec = np.array([x0, y0, z0 + r_obs])  # P - C, C = (0, 0, -(R + H))
        d = np.array([cx, cy, cz])
        if convention == "line":
            p = float(np.linalg.norm(np.cross(p_vec, d)))
        else:
            p = (radius_m + h_t) * math.sqrt(max(1.0 - cz * cz, 0.0))
        h_e = float(np.linalg.norm(p_vec)) - radius_m
        if not p < radius_m + h_t:
            out[i] = 0.0
            continue

        depth = _slant_depth(tab, lam, p, h_e, h_t, radius_m)
        out[i] = math.exp(-depth)
    return out


def run(
    root: Path, name: str, blocks, curvature=None, radius=None
) -> tuple[Path, Path, dict]:
    from obdeect.propagation import propagate

    telescopes = [(0.0, 0.0, TELESCOPE_Z_M), (0.0, 0.0, TELESCOPE_Z_M)]
    stream = uc.write_stream(
        root / f"{name}-in.pst",
        telescopes,
        blocks,
        curvature_applied=0 if curvature else None,
        wavelength_nm=WAVELENGTH_RANGE_NM,
    )
    out = root / f"{name}-out.pst"
    options = {}
    if curvature:
        options["curvature"] = curvature
    if radius is not None:
        options["earth_radius_m"] = radius
    summary = propagate(
        stream,
        out,
        transmission=uc.TRANSMISSION,
        site_altitude_m=uc.SITE_M,
        seed=1,
        mode="weighted",
        **options,
    )
    return stream, out, summary


def analyse(root: Path) -> tuple[dict, dict]:
    tab = uc.Transmission()
    rng = np.random.default_rng(SEED)
    up = upward_sample(rng, N_UP)
    down = downward_sample(rng, N_DOWN)
    for cols in (up, down):
        cols["w"] = np.ones(len(cols["z"]), dtype=np.float32)
    blocks_mixed = [(1, up), (2, down)]
    blocks_down_only = [(2, down)]

    # plane-parallel
    _, out_mixed, summary_plane = run(root, "plane-mixed", blocks_mixed)
    _, out_down_only, _ = run(root, "plane-down-only", blocks_down_only)
    rows = uc.read_rows(out_mixed)
    up_rows = rows["telescope"] == 1
    ref_plane = plane_reference(tab, up)
    w_plane = rows["w"][up_rows].astype(np.float64)
    plane_err = float(np.max(np.abs(w_plane - ref_plane)))
    plane_mismatch = _mismatch_for(rows, uc.read_rows(out_down_only), telescope=2)

    # spherical, earth radius CORSIKA
    _, out_sph, summary_sph = run(
        root, "sph-mixed", blocks_mixed, curvature="spherical"
    )
    _, out_sph_down, _ = run(
        root, "sph-down-only", blocks_down_only, curvature="spherical"
    )
    rows_s = uc.read_rows(out_sph)
    ref_sph = spherical_reference(tab, up, EARTH_RADIUS_M)
    ref_line = spherical_reference(tab, up, EARTH_RADIUS_M, convention="line")
    w_sph = rows_s["w"][rows_s["telescope"] == 1].astype(np.float64)
    sph_err = float(np.max(np.abs(w_sph - ref_sph)))
    sph_mismatch = _mismatch_for(rows_s, uc.read_rows(out_sph_down), telescope=2)

    # flat limit: spherical with R -> infinity against plane-parallel (reference: plane formula)
    _, out_flat, _ = run(
        root, "flat-mixed", [(1, up)], curvature="spherical", radius=FLAT_RADIUS_M
    )
    rows_f = uc.read_rows(out_flat)
    w_flat = rows_f["w"][rows_f["telescope"] == 1].astype(np.float64)
    flat_err = float(np.max(np.abs(w_flat - ref_plane)))

    metrics = {
        "sample_count": int(N_UP),
        "plane_max_abs_error": plane_err,
        "spherical_max_abs_error": sph_err,
        "flat_limit_max_abs_error": flat_err,
        "spherical_line_convention_difference": float(
            np.max(np.abs(ref_line - ref_sph))
        ),
        "plane_downward_bitwise_mismatch": plane_mismatch,
        "spherical_downward_bitwise_mismatch": sph_mismatch,
        "dropped_upward": int(summary_plane["dropped"]),
    }
    details = {
        "seed": SEED,
        "obdeect_propagation_version": version(),
        "transmission": str(uc.TRANSMISSION),
        "telescope_z_m": TELESCOPE_Z_M,
        "earth_radius_m": EARTH_RADIUS_M,
        "summary_plane": summary_plane,
        "summary_spherical": summary_sph,
        "reference": "plane: exp(-|tau(h_t)-tau(h_e)|/|cz|); spherical: exp(-int |dtau/dh| / cos(theta) dh), "
        "cos(theta) = sqrt(1 - p^2/(R+h)^2), p from the telescope vertical (curvature.md); quad with "
        "table heights as breakpoints",
        "telescope_depth_convention": "tau(lambda, h_T) as in obdeect (pathlen 0: focus depth 0)",
        "_arrays": {
            "ref_plane": ref_plane,
            "obd_plane": w_plane,
            "ref_sph": ref_sph,
            "obd_sph": w_sph,
        },
    }
    return metrics, details


def _mismatch_for(rows_mixed: dict, rows_only: dict, telescope: int) -> int:
    left = {
        k: (v[rows_mixed["telescope"] == telescope] if not k.startswith("_") else v)
        for k, v in rows_mixed.items()
    }
    right = {
        k: (v[rows_only["telescope"] == telescope] if not k.startswith("_") else v)
        for k, v in rows_only.items()
    }
    return e2e.compare_rows(left, right)["photon_mismatches"]


def plot(path: Path, details: dict) -> None:
    arrays = details.pop("_arrays")
    fig, axes = plotting.figure(1, 2, (9.0, 3.8))
    axes[0, 0].plot(
        arrays["ref_plane"],
        arrays["obd_plane"],
        "o",
        markersize=3,
        color=plotting.NATIVE,
        label="obdeect (plane)",
    )
    axes[0, 0].plot([0, 1], [0, 1], "--", color=plotting.ANALYTIC, label="equal")
    axes[0, 0].set(xlabel="analytic transmission", ylabel="obdeect transmission")
    axes[0, 0].legend(fontsize="small")
    axes[0, 1].plot(
        arrays["ref_sph"],
        arrays["obd_sph"],
        "o",
        markersize=3,
        color=plotting.NATIVE,
        label="obdeect (spherical)",
    )
    axes[0, 1].plot([0, 1], [0, 1], "--", color=plotting.ANALYTIC, label="equal")
    axes[0, 1].set(xlabel="slant integral reference", ylabel="obdeect transmission")
    axes[0, 1].legend(fontsize="small")
    plotting.save(fig, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args()
    root = Path.cwd()
    metrics, details = analyse(root)
    plot(root / "upward-extinction.png", details)
    common.write_results(metrics, details, root)


if __name__ == "__main__":
    main()
