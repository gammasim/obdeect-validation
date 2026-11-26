"""Upward photons bent by the atmosphere (plane mode, refraction_applied = 0 with a profile).

Case propagation.upward-refraction. Reference: scipy quad of Snell's law n(h) sin(theta) = K between the
emission height and the telescope plane (see upward_common.snell_integrals), with n from the atmprof
table, ln(n - 1) a natural cubic spline. Checks:

1. Upward rays from several zeniths and emission heights: lateral displacement (the integral of
   tan(theta) minus the straight line), exit direction (cx, cy scaled by sin_t / sin_e, cz = +sqrt), and
   travel time (integral of n / cos(theta) / c).
2. Turning over: an upward ray with sin(theta_e) just above the critical value (sin = 1 at the plane)
   must be dropped and counted as turns_over; no output row.
3. Reversed downward ray: a photon from above the plane (height H_top) with the same Snell constant K
   and the reversed horizontal direction must arrive with the reversed exit direction of the upward ray
   (same K at the plane), and its own displacement and time must match the quadrature. obdeect's
   downward branch and upward branch are both exercised; the two references come from one integrator.

Positions are placed so that the output is small: the telescope of each photon is placed at the straight
line's crossing point (x = path cx), so the output x, y is the lateral shift alone (float32 rounding ~ 1e-8
m). The time input is -t_ref (rounded to float32), so the output time is the residual (rounding ~ 1e-11
ns). The compared quantities are therefore the bending itself, at the stated accuracy.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common, plotting
from obdeect_validation.adapters import upward_common as uc

ZENITHS_DEG = (30.0, 60.0, 75.0, 80.0, 85.0, 88.0)
EMISSION_Z_M = (0.0, 200.0, 400.0)  # sky z above the observation level, below the plane
PLANE_Z_M = 500.0
DOWN_Z_M = 1500.0  # reversed downward photons start here, above the plane
AZIMUTH_DEG = 20.0
MODERATE_MAX_DEG = 80.0  # zenith up to which the 1e-6 ns gate applies
GRAZING_TIME_NS = 1.4e-5  # zenith 85 and 88 deg: obdeect documents 1.4e-5 ns for sub-interval changes at 87 to 89 deg (refraction.md)
SHIFT_GATE_M = 1e-6
TIME_GATE_NS = 1e-6
DIRECTION_GATE = 1e-6


def version() -> str:
    return common.package_version("obdeect-propagation")


def upward_cases(prof: uc.Profile):
    """(zenith, z_e, cx, cy, cz) per upward photon, float32 stored."""
    out = []
    for zenith in ZENITHS_DEG:
        for z in EMISSION_Z_M:
            cx, cy, cz = uc.stored_direction(zenith, AZIMUTH_DEG)
            out.append({"zenith": zenith, "z": z, "cx": cx, "cy": cy, "cz": cz})
    return out


def upward_reference(prof, case) -> dict:
    """Expected shift, exit direction and time of an upward photon (float64 from the stored float32 values).

    Documented convention (refraction.md, "Upward photons"; sky_to_ground): |cz| = sqrt(1 - cx^2 - cy^2)
    from the stored cx, cy, also for upward rows. "stored_cz" gives the same quantities with |cz| from the
    stored cz column; it is informational only.
    """
    cx, cy, cz_col = (float(case[k]) for k in ("cx", "cy", "cz"))
    s_e = math.hypot(cx, cy)
    cz = math.sqrt(1.0 - s_e * s_e)
    h_e = uc.SITE_M + case["z"]
    h_t = uc.SITE_M + PLANE_Z_M
    dh = h_t - h_e
    snell = uc.snell_integrals(prof, h_e, s_e, h_t)

    def shifts(cz_used):
        scale = (snell["horizontal"] - dh * s_e / cz_used) / s_e
        return scale * cx, scale * cy

    shift_x, shift_y = shifts(cz)
    stored_x, stored_y = shifts(cz_col)
    return {
        "path": dh / cz,
        "shift_x": shift_x,
        "shift_y": shift_y,
        "shift_x_stored_cz": stored_x,
        "shift_y_stored_cz": stored_y,
        "time_ns": snell["time_ns"],
        "exit_cx": cx * snell["sin_b"] / s_e,
        "exit_cy": cy * snell["sin_b"] / s_e,
        "exit_cz": math.sqrt(1.0 - snell["sin_b"] ** 2),
        "sin_t": snell["sin_b"],
        "h_e": h_e,
    }


def gate_ratio(error_x: float, error_y: float, out_x: float, out_y: float) -> float:
    """Largest |error| / max(SHIFT_GATE_M, float32 spacing of the stored output position), per component."""
    ratios = []
    for error, out in ((error_x, out_x), (error_y, out_y)):
        ulp = float(np.spacing(np.float32(abs(out))))
        ratios.append(abs(error) / max(SHIFT_GATE_M, ulp))
    return max(ratios)


def build_upward(prof, root: Path):
    """Stream of upward photons (one telescope each, at the straight crossing point) and references."""
    cases = upward_cases(prof)
    telescopes, blocks, refs, t_inputs = [], [], [], []
    for i, case in enumerate(cases, start=1):
        ref = upward_reference(prof, case)
        t_in = np.float32(-ref["time_ns"])
        tx = float(case["cx"]) * ref["path"]
        ty = float(case["cy"]) * ref["path"]
        telescopes.append((tx, ty, PLANE_Z_M))
        blocks.append(
            (i, _columns(0.0, 0.0, case["z"], case["cx"], case["cy"], case["cz"], t_in))
        )
        refs.append(ref)
        t_inputs.append(float(t_in))
    # Turning over: sin(theta_e) above the critical value of this emission height.
    z_turn = 100.0
    h_e = uc.SITE_M + z_turn
    s_crit = uc.critical_sin(prof, h_e, uc.SITE_M + PLANE_Z_M)
    s_turn = s_crit + 0.5 * (
        1.0 - s_crit
    )  # below 1, but the plane is reached with sin > 1
    cx_t, cy_t, cz_t = _stored_from_sin(s_turn, AZIMUTH_DEG)
    turn_id = len(cases) + 1
    telescopes.append((0.0, 0.0, PLANE_Z_M))
    blocks.append(
        (turn_id, _columns(0.0, 0.0, z_turn, cx_t, cy_t, cz_t, np.float32(0.0)))
    )
    path = uc.write_stream(root / "upward-refraction-in.pst", telescopes, blocks)
    return (
        path,
        cases,
        refs,
        t_inputs,
        {"critical_sin": s_crit, "sin_turnover": s_turn, "turn_id": turn_id},
    )


def _columns(x, y, z, cx, cy, cz, t):
    f = np.float32
    return {
        "x": np.array([x], dtype=f),
        "y": np.array([y], dtype=f),
        "z": np.array([z], dtype=f),
        "cx": np.array([cx], dtype=f),
        "cy": np.array([cy], dtype=f),
        "cz": np.array([cz], dtype=f),
        "t": np.array([t], dtype=f),
        "w": np.ones(1, dtype=f),
        "wavelength": np.full(1, 400.0, dtype=f),
    }


def _stored_from_sin(sin_theta: float, azimuth_deg: float):
    az = math.radians(azimuth_deg)
    return (
        np.float32(sin_theta * math.cos(az)),
        np.float32(sin_theta * math.sin(az)),
        np.float32(math.sqrt(1.0 - sin_theta**2)),
    )


def build_reversed(prof, cases, root: Path):
    """Downward photons from DOWN_Z_M with the Snell constant of each upward ray, reversed direction.

    Convention of obdeect for downward rows (docs/propagation/propagation.md, sky.hpp): |cz| is
    sqrt(1 - cx^2 - cy^2) from the stored float32 cx, cy, not the stored cz column. The telescope of each
    photon is placed on the path of that convention. Two references are returned: "documented" (|cz| and the
    sine from the stored cx, cy) and "stored_cz" (|cz| from the stored cz column); their difference is the
    sensitivity of the lateral shift to the float32 rounding of the direction.
    """
    telescopes, blocks, refs, t_inputs = [], [], [], []
    h_top = uc.SITE_M + DOWN_Z_M
    h_t = uc.SITE_M + PLANE_Z_M
    dh = DOWN_Z_M - PLANE_Z_M
    for i, case in enumerate(cases, start=1):
        cx, cy, _cz = (float(case[k]) for k in ("cx", "cy", "cz"))
        s_e = math.hypot(cx, cy)
        k = (1.0 + float(prof.nu(uc.SITE_M + case["z"]))) * s_e
        s_top = k / (1.0 + float(prof.nu(h_top)))
        if not s_top < 1.0:
            raise RuntimeError("reversed ray has no real direction above the plane")
        cx_d = np.float32(-cx / s_e * s_top)
        cy_d = np.float32(-cy / s_e * s_top)
        cz_d = np.float32(-math.sqrt(1.0 - s_top * s_top))
        sin_doc = math.hypot(float(cx_d), float(cy_d))
        cz_doc = math.sqrt(max(1.0 - sin_doc * sin_doc, 0.0))
        cz_store = abs(float(cz_d))
        snell = uc.snell_integrals(prof, h_top, sin_doc, h_t)
        t_in = np.float32(-snell["time_ns"])
        path_doc = dh / cz_doc
        telescopes.append((float(cx_d) * path_doc, float(cy_d) * path_doc, PLANE_Z_M))
        blocks.append((i, _columns(0.0, 0.0, DOWN_Z_M, cx_d, cy_d, cz_d, t_in)))

        doc_scale = (snell["horizontal"] - dh * sin_doc / cz_doc) / sin_doc
        st_scale = (snell["horizontal"] - dh * sin_doc / cz_store) / sin_doc
        doc_x, doc_y = doc_scale * float(cx_d), doc_scale * float(cy_d)
        st_x, st_y = st_scale * float(cx_d), st_scale * float(cy_d)
        refs.append(
            {
                "shift_x": doc_x,
                "shift_y": doc_y,
                "shift_x_stored_cz": st_x,
                "shift_y_stored_cz": st_y,
                "time_ns": snell["time_ns"],
                "exit_cx": -cx * snell["sin_b"] / s_e,  # reversed upward exit
                "exit_cy": -cy * snell["sin_b"] / s_e,
                "exit_cz": -math.sqrt(1.0 - snell["sin_b"] ** 2),
            }
        )
        t_inputs.append(float(t_in))
    path = uc.write_stream(
        root / "upward-refraction-reversed-in.pst", telescopes, blocks
    )
    return path, refs, t_inputs


def analyse(root: Path) -> tuple[dict, dict]:
    from obdeect.propagation import propagate

    prof = uc.profile()
    stream, cases, refs, t_in, turn = build_upward(prof, root)
    out_up = root / "upward-refraction-out.pst"
    summary_up = propagate(
        stream,
        out_up,
        extinction=False,
        site_altitude_m=uc.SITE_M,
        seed=1,
        profile=uc.PROFILE,
    )
    rows_up = uc.read_rows(out_up)
    by_tel = {int(tel): idx for idx, tel in enumerate(rows_up["telescope"])}

    shift_err, shift_ratio, shift_stored_err, dir_err, time_err, zen_of = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for i, (case, ref) in enumerate(zip(cases, refs), start=1):
        idx = by_tel[i]
        x, y = float(rows_up["x"][idx]), float(rows_up["y"][idx])
        t = float(rows_up["t"][idx])
        exit_out = np.array(
            [rows_up["cx"][idx], rows_up["cy"][idx], rows_up["cz"][idx]],
            dtype=np.float64,
        )
        exit_ref = np.array([ref["exit_cx"], ref["exit_cy"], ref["exit_cz"]])
        ex, ey = x - ref["shift_x"], y - ref["shift_y"]
        shift_err.append(max(abs(ex), abs(ey)))
        shift_ratio.append(gate_ratio(ex, ey, x, y))
        shift_stored_err.append(
            max(abs(x - ref["shift_x_stored_cz"]), abs(y - ref["shift_y_stored_cz"]))
        )
        dir_err.append(float(np.max(np.abs(exit_out - exit_ref))))
        time_err.append(abs(t - (t_in[i - 1] + ref["time_ns"])))
        zen_of.append(case["zenith"])
    shift_err = np.array(shift_err)
    shift_ratio = np.array(shift_ratio)
    shift_stored_err = np.array(shift_stored_err)
    dir_err = np.array(dir_err)
    time_err = np.array(time_err)
    zen_of = np.array(zen_of)
    moderate = zen_of <= MODERATE_MAX_DEG

    turn_id = turn["turn_id"]
    turned = turn_id not in by_tel
    turn_counts = summary_up["dropped_by_reason"]

    down_path, refs_down, t_down = build_reversed(prof, cases, root)
    out_down = root / "upward-refraction-reversed-out.pst"
    summary_down = propagate(
        down_path,
        out_down,
        extinction=False,
        site_altitude_m=uc.SITE_M,
        seed=1,
        profile=uc.PROFILE,
    )
    rows_down = uc.read_rows(out_down)
    by_tel_d = {int(tel): idx for idx, tel in enumerate(rows_down["telescope"])}
    rev_shift, rev_ratio, rev_dir, rev_time = [], [], [], []
    rev_by_zenith = []
    zen_rev = [c["zenith"] for c in cases]
    for i, ref in enumerate(refs_down, start=1):
        idx = by_tel_d[i]
        x, y, t = (
            float(rows_down["x"][idx]),
            float(rows_down["y"][idx]),
            float(rows_down["t"][idx]),
        )
        out_d = np.array(
            [rows_down["cx"][idx], rows_down["cy"][idx], rows_down["cz"][idx]],
            dtype=np.float64,
        )
        ex, ey = x - ref["shift_x"], y - ref["shift_y"]
        rev_shift.append(max(abs(ex), abs(ey)))
        rev_ratio.append(gate_ratio(ex, ey, x, y))
        rev_by_zenith.append(
            {
                "zenith_deg": zen_rev[i - 1],
                "z_m": cases[i - 1]["z"],
                "shift_m": float(max(abs(ex), abs(ey))),
                "shift_gate_ratio": float(gate_ratio(ex, ey, x, y)),
                "shift_stored_cz_m": float(
                    max(
                        abs(x - ref["shift_x_stored_cz"]),
                        abs(y - ref["shift_y_stored_cz"]),
                    )
                ),
                "time_ns": float(abs(t - (t_down[i - 1] + ref["time_ns"]))),
            }
        )
        rev_dir.append(
            float(
                np.max(
                    np.abs(
                        out_d
                        - np.array([ref["exit_cx"], ref["exit_cy"], ref["exit_cz"]])
                    )
                )
            )
        )
        rev_time.append(abs(t - (t_down[i - 1] + ref["time_ns"])))

    metrics = {
        "sample_count": len(cases),
        "max_shift_error_m": float(shift_err.max()),
        "max_shift_gate_ratio": float(shift_ratio.max()),
        "max_direction_error": float(dir_err.max()),
        "max_time_error_ns_moderate": float(time_err[moderate].max()),
        "max_time_error_ns_grazing": float(time_err[~moderate].max()),
        "turns_over_reported": int(turned and turn_counts["turns_over"] == 1),
        "turns_over_other_drops": int(
            sum(v for k, v in turn_counts.items() if k != "turns_over")
        ),
        "reversed_max_shift_error_m": float(max(rev_shift)),
        "reversed_max_shift_gate_ratio": float(max(rev_ratio)),
        "reversed_max_direction_error": float(max(rev_dir)),
        "reversed_max_time_error_ns": float(max(rev_time)),
        "reversed_dropped": int(summary_down["dropped"]),
    }
    details = {
        "seed": 0,
        "obdeect_propagation_version": version(),
        "profile": str(uc.PROFILE),
        "zeniths_deg": list(ZENITHS_DEG),
        "emission_z_m": list(EMISSION_Z_M),
        "plane_z_m": PLANE_Z_M,
        "downward_start_z_m": DOWN_Z_M,
        "summary_upward": summary_up,
        "summary_reversed": summary_down,
        "turnover": turn,
        "upward_shift_error_by_zenith_m": {
            str(zen): float(np.max(shift_err[zen_of == zen])) for zen in ZENITHS_DEG
        },
        "upward_shift_gate_ratio_by_zenith": {
            str(zen): float(np.max(shift_ratio[zen_of == zen])) for zen in ZENITHS_DEG
        },
        "upward_shift_error_stored_cz_informational_by_zenith_m": {
            str(zen): float(np.max(shift_stored_err[zen_of == zen]))
            for zen in ZENITHS_DEG
        },
        "reversed_by_zenith": rev_by_zenith,
        "gates": {
            "shift_m": "max(1e-6 m, float32 spacing of the stored output position), per component",
            "time_ns_moderate": TIME_GATE_NS,
            "time_ns_grazing": GRAZING_TIME_NS,
            "direction": DIRECTION_GATE,
        },
        "time_error_ns_by_zenith": {
            str(z): float(time_err[zen_of == z].max()) for z in ZENITHS_DEG
        },
        "reference": "scipy quad of n(h) sin(theta) = K, n from the atmprof table, ln(n-1) natural cubic spline",
        "time_convention": "t_out - (t_in + integral of n / cos(theta) dh / c), t_in = -reference (float32)",
    }
    details["_arrays"] = {"zenith": zen_of, "time": time_err}
    plot(root / "upward-refraction.png", details)
    return metrics, details


def plot(path: Path, details: dict) -> None:
    arrays = details.pop("_arrays")
    fig, axis = plotting.figure(1, 1, (6.0, 4.0))
    axis.plot(
        arrays["zenith"],
        np.maximum(arrays["time"], 1e-12),
        "o",
        color=plotting.NATIVE,
        label="obdeect time error (ns)",
    )
    axis.axhline(
        TIME_GATE_NS,
        color=plotting.ANALYTIC,
        linestyle="--",
        linewidth=1,
        label="1e-6 ns gate",
    )
    axis.set(xlabel="zenith angle [deg]", ylabel="travel time error [ns]", yscale="log")
    axis.legend(loc="upper left", fontsize="small")
    plotting.save(fig, path)


if __name__ == "__main__":
    common.results_main(analyse, __doc__.splitlines()[0])
