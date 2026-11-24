"""Flat-atmosphere refraction at high zenith angles: displacement and arrival time.

Case propagation.refraction-high-zenith. Photons are emitted at 5, 10 and 20 km a.s.l. at zenith
60, 70, 75, 80, 82, 85 and 87 deg (azimuth 0, observation level 2000 m, ECMWF profile fixture).

References
----------
1. Primary, independent of obdeect: numerical integration of the flat-atmosphere ray equation in
   depth form, u = h_e - h. Snell's law n sin(theta) = const, differentiated, gives
       d(theta)/du = (dn/dh) / n * tan(theta)
   and the displacement and time integrands are
       dx/du = tan(theta) - tan(theta_e),   dE/du = n/cos(theta) - 1/cos(theta_e)  [m]
   so that dx is the displacement against the straight line and the air excess is
   E / c; the vacuum path (h_e - obs)/(cos theta_e) is added to the time. The ODE is
   integrated by scipy solve_ivp (DOP853, rtol 1e-12, atol 1e-12 m) cell by cell. n - 1 is
   the model obdeect uses (docs/propagation/refraction.md, "Interpolation"): natural cubic
   spline of ln(n - 1) through the levels, linear in ln(n - 1) between nodes of at most 3 m
   (uc.Profile, the same as exponential inside each fine cell). Within a cell the integrand is
   smooth, so the step control is exact up to the stated tolerance. Cross-checked against
   scipy quad of the Snell integrals (upward_common.snell_integrals) in the details.
2. Secondary: CORSIKA IACT/atmo raybnd() (the compiled reference of the atmosphere-corsika
   case), run on the same profile and observation level. raybnd() is an approximation (mean
   angle 0.28 and fast tables), so its residuals are reported, not gated.

Obdeect candidate
-----------------
Stage candidate: a sky-frame stream with one photon per telescope, telescope at the straight
line's crossing of the observation plane (so the stored output x is the displacement alone), and
the arrival time referenced to the exact float32 rounding of the reference time, so the stored
time is the residual. obdeect propagation run, mode none, refraction_applied = 0 (the stage
uses the series `bend` and the quadrature `bend_quadrature` beyond its series limit).
Kernel candidate: obdeect.propagation.refract (double precision, series only); photons that
the series rejects are reported as not bent and are not gated there.

Gates (fixed before the run; no tolerance is chosen from a result)
-------------------------------------------------------------------
Displacement [m] and time [ns], per zenith angle:
- zenith <= 82 deg: 1e-6 m and 1e-6 ns. refraction.md:73 states the series truncation error is
  below 1e-6 m for zenith angles up to 82 deg; refraction.md:99-100 states agreement with a
  numerical integration to 1e-6 m and 1e-6 ns on the ECMWF table.
- zenith 85 and 87 deg: 4e-6 m and 1.4e-5 ns. refraction.md:247-248 states the only figures for
  the 87-89 deg band (quadrature sub-interval halving changes the displacement by at most 4e-6 m
  and the time by at most 1.4e-5 ns). The docs give no absolute bound at 83-86 deg; these figures
  are applied there as the nearest stated bound, and the gate is borrowed, not measured.
- Stage candidate only: the displacement gate is max(gate, ulp32(dx_ref)), the rounding of the
  float32 output column (one float32 spacing at |dx|, refraction.md and curvature.md note the
  float32 time/position column limits). The time gate needs no floor: the stored residual is
  near zero (sub-ps rounding).
Counts: displacement_violations and time_violations (stage candidate, expected 0),
kernel_displacement_violations and kernel_time_violations (kernel candidate, bent photons only,
expected 0). CORSIKA residuals and the quadrature cross-check are in details only.
"""

from __future__ import annotations

import itertools
import math
import os
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from obdeect_validation.adapters import common
from obdeect_validation.adapters import upward_common as uc

ZENITHS_DEG = (60.0, 70.0, 75.0, 80.0, 82.0, 85.0, 87.0)
HEIGHTS_M = (5000.0, 10000.0, 20000.0)
OBS_M = uc.SITE_M  # 2000 m a.s.l.
AZIMUTH_DEG = 0.0
C_M_NS = uc.C_M_NS
ODE_RTOL = 1e-12
ODE_ATOL = 1e-12  # m, for the displacement and the air excess integrals
GATE_LOW_ZENITH_M = 1e-6
GATE_LOW_ZENITH_NS = 1e-6
GATE_HIGH_ZENITH_M = 4e-6
GATE_HIGH_ZENITH_NS = 1.4e-5
HIGH_ZENITH_FROM_DEG = 85.0


def version() -> str:
    return common.package_version("obdeect-propagation")


def gates(zenith_deg: float) -> tuple[float, float]:
    """(displacement gate [m], time gate [ns]) for a zenith angle, as fixed in the module docstring."""
    if zenith_deg >= HIGH_ZENITH_FROM_DEG:
        return GATE_HIGH_ZENITH_M, GATE_HIGH_ZENITH_NS
    return GATE_LOW_ZENITH_M, GATE_LOW_ZENITH_NS


# ----------------------------------------------------------------- independent reference


def ray_reference(prof: uc.Profile, h_e: float, sin_e: float, obs: float = OBS_M):
    """Displacement [m] and air excess [ns] of a downward ray from h_e to obs (flat, Snell).

    Depth form of the ray equation, integrated by solve_ivp on each cell of the profile grid
    (the nodes of uc.Profile, where ln(n - 1) is linear). Returns (dx_m, excess_ns, theta_obs).
    dx is the displacement against the straight line along the emission direction, at obs,
    in the direction of (cx, cy)/sin_e (the same convention as the stage).
    """
    theta_e = math.asin(sin_e)
    tan_e = math.tan(theta_e)
    inv_cos_e = 1.0 / math.cos(theta_e)
    grid = prof.grid
    lnnu = prof.lnnu
    inner = grid[(grid > obs) & (grid < h_e)]
    points = sorted({h_e, *inner.tolist(), obs}, reverse=True)
    y = np.array([0.0, theta_e, 0.0])  # dx, theta, excess (m)
    for upper, lower in itertools.pairwise(points):
        mid = 0.5 * (upper + lower)
        k = min(
            max(int(np.searchsorted(grid, mid, side="right")) - 1, 0), len(grid) - 2
        )
        slope = (lnnu[k + 1] - lnnu[k]) / (grid[k + 1] - grid[k])

        def rhs(u, state, k=k, slope=slope):
            h = h_e - u
            nu = math.exp(lnnu[k] + slope * (h - grid[k]))
            theta = state[1]
            tan_t = math.tan(theta)
            n = 1.0 + nu
            return [
                tan_t - tan_e,
                slope * nu / n * tan_t,
                n / math.cos(theta) - inv_cos_e,
            ]

        sol = solve_ivp(
            rhs,
            (h_e - upper, h_e - lower),
            y,
            method="DOP853",
            rtol=ODE_RTOL,
            atol=[ODE_ATOL, 1e-14, ODE_ATOL],
        )
        y = sol.y[:, -1]
    return float(y[0]), float(y[2]) / C_M_NS, float(y[1])


def total_time_ns(
    h_e: float, sin_e: float, excess_ns: float, obs: float = OBS_M
) -> float:
    """Vacuum path along the emission direction plus the air excess [ns]."""
    cos_e = math.sqrt(1.0 - sin_e * sin_e)
    return (h_e - obs) / (C_M_NS * cos_e) + excess_ns


def quadrature_check(prof: uc.Profile, h_e: float, sin_e: float) -> tuple[float, float]:
    """Difference of the ODE reference from scipy quad of the Snell integrals (dx [m], time [ns])."""
    cos_e = math.sqrt(1.0 - sin_e * sin_e)
    dx, excess, _ = ray_reference(prof, h_e, sin_e)
    quad = uc.snell_integrals(prof, h_e, sin_e, OBS_M)
    dx_quad = quad["horizontal"] - (h_e - OBS_M) * sin_e / cos_e
    t_quad = quad["time_ns"] - (h_e - OBS_M) / (C_M_NS * cos_e)
    return dx - dx_quad, excess - t_quad


# ----------------------------------------------------------------- CORSIKA secondary reference


def corsika_reference(
    rays: list[tuple[float, float, float]], profile: Path
) -> list[tuple]:
    """CORSIKA raybnd for (h_e [m], u = sin, w = cos); returns (dx_m, dt_ns) per ray."""
    corsika_dir = Path(os.environ.get("CORSIKA_ATMO_DIR", "/reference/corsika"))
    support_dir = Path(os.environ.get("CORSIKA_SUPPORT_DIR", str(corsika_dir)))
    include_dirs = [corsika_dir] + [
        Path(p)
        for p in os.environ.get("CORSIKA_INCLUDE_DIRS", "").split(os.pathsep)
        if p
    ]
    wrapper = Path(__file__).with_name("refraction_high_zenith_corsika.c")
    executable = Path("corsika-refraction-high-zenith")
    include_args = [arg for directory in include_dirs for arg in ("-I", str(directory))]
    subprocess.run(
        [
            "cc",
            "-std=c11",
            "-O2",
            *include_args,
            str(wrapper),
            str(corsika_dir / "atmo.c"),
            str(support_dir / "fileopen.c"),
            str(support_dir / "straux.c"),
            "-lm",
            "-o",
            str(executable),
        ],
        check=True,
    )
    text = "".join(f"{h * 100.0:.17g} {u:.17g} {w:.17g}\n" for h, u, w in rays)
    run = subprocess.run(
        [str(executable.resolve()), str(profile), f"{OBS_M:.17g}"],
        input=text,
        check=True,
        capture_output=True,
        text=True,
    )
    out = []
    for line in run.stdout.splitlines():
        fields = line.split()
        if len(fields) != 4:  # CORSIKA prints its own banner on stdout
            continue
        try:
            values = [float(v) for v in fields]
        except ValueError:
            continue
        out.append((values[2] / 100.0, values[3]))
    if len(out) != len(rays):
        raise RuntimeError(f"CORSIKA returned {len(out)} rays, expected {len(rays)}")
    return out


# ----------------------------------------------------------------- obdeect candidates


def _columns(z, cx, cy, cz, t):
    f = np.float32
    return {
        "x": np.array([0.0], dtype=f),
        "y": np.array([0.0], dtype=f),
        "z": np.array([z], dtype=f),
        "cx": np.array([cx], dtype=f),
        "cy": np.array([cy], dtype=f),
        "cz": np.array([cz], dtype=f),
        "t": np.array([t], dtype=f),
        "w": np.ones(1, dtype=f),
        "wavelength": np.full(1, 400.0, dtype=f),
    }


def build_rays():
    """(h_e, zenith, cx, cy, cz_doc) per ray, cx and cy as stored float32 (the reference reads them back)."""
    rays = []
    az = math.radians(AZIMUTH_DEG)
    for h_e in HEIGHTS_M:
        for zenith in ZENITHS_DEG:
            th = math.radians(zenith)
            cx = np.float32(math.sin(th) * math.cos(az))
            cy = np.float32(math.sin(th) * math.sin(az))
            s = math.hypot(float(cx), float(cy))
            rays.append(
                {
                    "h_e": h_e,
                    "zenith": zenith,
                    "cx": float(cx),
                    "cy": float(cy),
                    "sin_e": s,
                    "cz_doc": math.sqrt(1.0 - s * s),
                }
            )
    return rays


def run_stage(rays, root: Path):
    from obdeect.propagation import propagate

    telescopes, blocks = [], []
    for i, ray in enumerate(rays, start=1):
        z = ray["h_e"] - OBS_M
        path = z / ray["cz_doc"]
        tx = ray["cx"] * path
        telescopes.append((tx, ray["cy"] * path, 0.0))
        # reference time, rounded to float32 for the input; the output is the residual
        t_ref = total_time_ns(ray["h_e"], ray["sin_e"], ray["excess_ns"])
        t_in = np.float32(-t_ref)
        ray["t_in"] = float(t_in)
        ray["t_ref"] = t_ref
        blocks.append((i, _columns(z, ray["cx"], ray["cy"], -ray["cz_doc"], t_in)))
    stream = uc.write_stream(
        root / "refraction-high-zenith-in.pst", telescopes, blocks, refraction_applied=0
    )
    summary = propagate(
        stream,
        root / "refraction-high-zenith-out.pst",
        extinction=False,
        site_altitude_m=OBS_M,
        seed=1,
        profile=uc.PROFILE,
    )
    rows = uc.read_rows(root / "refraction-high-zenith-out.pst")
    by = {int(t): k for k, t in enumerate(rows["telescope"])}
    for i, ray in enumerate(rays, start=1):
        if i not in by:
            ray["stage"] = None
            continue
        k = by[i]
        ray["stage"] = {
            "x": float(rows["x"][k]),
            "t": float(rows["t"][k]),
            "cx": float(rows["cx"][k]),
        }
    return summary


def kernel_run(rays):
    from obdeect.propagation import AtmosphereProfile, RefractionTables, refract

    tables = RefractionTables(AtmosphereProfile(str(uc.PROFILE)), OBS_M)
    h = np.array([r["h_e"] for r in rays])
    cx = np.array([r["cx"] for r in rays])
    cy = np.array([r["cy"] for r in rays])
    zeros = np.zeros(len(rays))
    _, _, x, _, t, bent = refract(tables, h, cx, cy, zeros, zeros, zeros)
    for r, xi, ti, bi in zip(rays, x, t, bent):
        r["kernel"] = {"x": float(xi), "t": float(ti), "bent": bool(bi)}


# ----------------------------------------------------------------- main


def main() -> None:
    root = Path.cwd()
    prof = uc.profile()
    rays = build_rays()
    for ray in rays:
        dx, excess, th_obs = ray_reference(prof, ray["h_e"], ray["sin_e"])
        ray["dx_ref"] = dx * (ray["cx"] / ray["sin_e"])  # direction of (cx, cy)/sin_e
        ray["excess_ns"] = excess
        ray["theta_obs"] = th_obs
        ray["quad_diff"] = quadrature_check(prof, ray["h_e"], ray["sin_e"])
    corsika = corsika_reference(
        [(ray["h_e"], ray["sin_e"], ray["cz_doc"]) for ray in rays], uc.PROFILE
    )
    for ray, (dx_corsika, dt_corsika) in zip(rays, corsika):
        ray["corsika"] = {
            "dx": dx_corsika * (ray["cx"] / ray["sin_e"]),
            "dt": dt_corsika,
        }
    summary = run_stage(rays, root)
    kernel_run(rays)

    rows = []
    for ray in rays:
        gate_dx, gate_t = gates(ray["zenith"])
        stage = ray["stage"]
        dx_ref = ray["dx_ref"]
        t_ref = ray["t_ref"]
        dx_err = stage["x"] - dx_ref
        t_err = (stage["t"] - ray["t_in"]) - t_ref
        dx_gate = max(gate_dx, float(uc.ulp32([dx_ref])[0]))
        rec = {
            "h_e_m": ray["h_e"],
            "zenith_deg": ray["zenith"],
            "dx_ref_m": dx_ref,
            "dx_stage_m": stage["x"],
            "dx_error_m": dx_err,
            "dx_gate_m": dx_gate,
            "dx_violation": bool(abs(dx_err) > dx_gate),
            "t_ref_ns": t_ref,
            "t_error_ns": t_err,
            "t_gate_ns": gate_t,
            "t_violation": bool(abs(t_err) > gate_t),
            "series_bent": ray["kernel"]["bent"],
            "theta_obs_ref_rad": ray["theta_obs"],
            "direction_cx_error": stage["cx"] - math.sin(ray["theta_obs"]),
            "quadrature_check_dx_m": ray["quad_diff"][0],
            "quadrature_check_ns": ray["quad_diff"][1],
            "corsika_dx_error_m": ray["corsika"]["dx"] - dx_ref,
            "corsika_t_error_ns": ray["corsika"]["dt"] - t_ref,
        }
        if ray["kernel"]["bent"]:
            k_dx_err = ray["kernel"]["x"] - dx_ref
            k_t_err = ray["kernel"]["t"] - t_ref
            rec["kernel_dx_error_m"] = k_dx_err
            rec["kernel_t_error_ns"] = k_t_err
            rec["kernel_dx_violation"] = bool(abs(k_dx_err) > gate_dx)
            rec["kernel_t_violation"] = bool(abs(k_t_err) > gate_t)
        rows.append(rec)

    metrics = {
        "sample_count": len(rows),
        "displacement_violations": sum(r["dx_violation"] for r in rows),
        "time_violations": sum(r["t_violation"] for r in rows),
        "kernel_displacement_violations": sum(
            r.get("kernel_dx_violation", False) for r in rows
        ),
        "kernel_time_violations": sum(r.get("kernel_t_violation", False) for r in rows),
    }
    per_zenith = {}
    for zen in ZENITHS_DEG:
        sel = [r for r in rows if r["zenith_deg"] == zen]
        per_zenith[str(zen)] = {
            "max_abs_dx_error_m": max(abs(r["dx_error_m"]) for r in sel),
            "max_abs_t_error_ns": max(abs(r["t_error_ns"]) for r in sel),
            "dx_gate_m": gates(zen)[0],
            "t_gate_ns": gates(zen)[1],
            "series_bent_count": sum(r["series_bent"] for r in sel),
            "max_corsika_abs_dx_error_m": max(
                abs(r["corsika_dx_error_m"]) for r in sel
            ),
            "max_corsika_abs_t_error_ns": max(
                abs(r["corsika_t_error_ns"]) for r in sel
            ),
        }
    details = {
        "obdeect_propagation_version": version(),
        "profile": str(uc.PROFILE),
        "observation_level_m": OBS_M,
        "heights_m": list(HEIGHTS_M),
        "zeniths_deg": list(ZENITHS_DEG),
        "azimuth_deg": AZIMUTH_DEG,
        "ode": {"method": "DOP853", "rtol": ODE_RTOL, "atol_m": ODE_ATOL},
        "stage_summary": summary,
        "per_zenith": per_zenith,
        "per_photon": rows,
        "max_quadrature_check_dx_m": max(abs(r["quadrature_check_dx_m"]) for r in rows),
        "max_quadrature_check_ns": max(abs(r["quadrature_check_ns"]) for r in rows),
        "gates": {
            "zenith_le_82": {
                "dx_m": GATE_LOW_ZENITH_M,
                "t_ns": GATE_LOW_ZENITH_NS,
                "source": "refraction.md:73, 99-100",
            },
            "zenith_ge_85": {
                "dx_m": GATE_HIGH_ZENITH_M,
                "t_ns": GATE_HIGH_ZENITH_NS,
                "source": "refraction.md:247-248 (87-89 deg, borrowed for 85 deg)",
            },
            "stage_dx_floor": "float32 ulp of |dx_ref| (uc.ulp32)",
        },
        "reference": "flat Snell ray equation, depth form, scipy solve_ivp DOP853 rtol 1e-12; "
        "CORSIKA IACT/atmo raybnd secondary (not gated)",
    }
    common.write_results(metrics, details)


if __name__ == "__main__":
    main()
