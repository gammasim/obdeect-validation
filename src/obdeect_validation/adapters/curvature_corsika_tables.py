"""CORSIKA CURVED tables (DISTIP, TOFIP) against obdeect's spherical refraction, from a traced run.

Input. One trace file per primary zenith, written by a CURVED CORSIKA run that carries the output
patch fixtures/corsika-patches/curvature-tables-debug.patch and the environment variable
OBDEECT_CURVTAB (docs/curvature-corsika-tables.md, section 4). Names: curvature-tables-<zen>deg.dat
under VALIDATION_INPUTS (default /validation/inputs), zenith 65 and 75 deg. Missing inputs end the
case with exit code 3 and a message.

Trace format (first line 'CURVTAB 1'): 'CONST' lines (constants and flags of the run), 'ROW' lines
(table rows at the grid nodes, DISTEF and TOF) and 'BUNCH' lines (one per bunch that hits a
telescope, fixed width, see the patch). Units are CORSIKA's: cm, ns, cosines.

(a) Interpolation consistency: metrics with exact limits. Each bunch line carries the four table
    values of its cell for DISTEF and TOF. The bilinear interpolation of DISTIP and TOFIP is
    re-implemented in numpy with the same arithmetic as the Fortran (cell index INT(4 theta_deg),
    INT(4 h_km), the height clamp to the observation level, the two-step linear interpolation).
    The recomputed values are compared with the printed RDIST and TOFIP. The printed cell indices
    are compared with the ones recomputed from WEMIS and ZAPP. CARTIM = TEMIS*1e9 + TOFIP is
    checked the same way. The 'ROW' lines are compared with the cell corners where a corner is a
    printed grid node (details only: the node check depends on which cells the bunches fall in).

(b) Mapping comparison: details, plots, no pass or fail limits. The mapping between CORSIKA's
    apparent quantities (WEMIS, ZAPP) and obdeect's vacuum quantities (theta_s, h_e) is a
    HYPOTHESIS (design note, item 2), not verified. Under it:
      theta_s = arccos(WEMIS) (vacuum zenith at the telescope), h_e = ZEM (emission altitude).
    obdeect's shift and excess time at (theta_s, h_e) come from the numpy reference
    (refraction_curved_reference.ray_integrals) on the sphere R = C(1) with the telescope at the
    observation level. The public obdeect propagate() path is not used here: its stream helper
    places the telescope at 2000 m, the CORSIKA runs at 2156 m.
    Both sides are evaluated at obdeect's emission altitude h_e = ZEM, so the residuals test the
    whole hypothesis together with the altitude residual below.
    Residuals (sign and convention fixed here, not tuned to data):
      shift:        RDIST - [(h_e - h_obs) tan(theta_s) + shift_obdeect]          [m]
                    RDIST is the horizontal distance accumulated along the apparent ray
                    (DDIST = DHAPP tan(angle) in INRTAB; CERENK uses it as the full ground
                    offset), so it contains the straight part (h_e - h_obs) tan(theta_s), not
                    only the refraction shift. obdeect's shift is signed (the difference of the
                    horizontal travel of the refracted and the vacuum ray), so it is added as is.
      excess time:  [TOFIP - (h_e - h_obs) / (cos(theta_s) c)] - excess_obdeect   [ns]
                    the straight time to the tangent plane, the same reference as obdeect's
                    t = s/c + excess (s = (h_e - h_obs)/cos(theta_s)).
      altitude:     HE(ZAPP, theta_s) - ZEM, with HE the true altitude of the apparent geometry
                    (INRTAB's DHELP and HE formulas). Zero if H holds for the height part.
    The residuals are binned in zenith (theta_s) and height (h_e). Nothing is gated: the limits
    cannot be set honestly before the mapping is verified. The first real run decides it.
"""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common
from obdeect_validation.adapters import refraction_curved_reference as refraction

FORMAT_VERSION = "1"
ZENITHS_DEG = (65, 75)
TRACE_NAME = "curvature-tables-{zenith:g}deg.dat"

BUNCH_TAG = "BUNCH"
BUNCH_OFFSET = 34  # columns before the reals: 5 + I10 + I7 + I6 + I6
FIELD_WIDTH = 24  # 1X + E23.15
BUNCH_REALS = (
    "wemis",
    "zapp_cm",
    "zem_cm",
    "sthe",
    "rdist_cm",
    "tofip_ns",
    "temis_ns",
    "cartim_ns",
    "umean",
    "vmean",
    "wmean",
    "uemis2",
    "vemis2",
    "cos_cherenkov",
)
CELL_VALUES = ("distef_n2m2", "distef_n2m1", "distef_n1m2", "distef_n1m1")
CELL_TIMES = ("tof_n2m2", "tof_n2m1", "tof_n1m2", "tof_n1m1")
BUNCH_FIELDS = len(BUNCH_REALS) + len(CELL_VALUES) + len(CELL_TIMES)  # 22

ROW_THETA_DEG = (0, 30, 60, 70, 75, 80, 85, 89)
ROW_HEIGHT_KM = (3, 8, 20, 60)
REQUIRED_CONSTANTS = ("C1_CM", "C25_CMPS", "OBSLEV_CM", "ETADSN", "MODATM")

INTERPOLATION_LIMIT = (
    1e-6  # relative; the metric's limit lives in cases/propagation/...json
)
FORTRAN_FACT = 180.0 / math.pi  # FACT in DISTIP and TOFIP
C_M_PER_NS = 0.299792458
ZENITH_EDGES_DEG = (60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0)
HEIGHT_EDGES_KM = (0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 60.0, 113.0)

# Limits fixed before looking at results. Each entry: (rule, justification).
LIMITS = {
    "table_interpolation_max_relative_error": (
        "<= 1e-6 (relative to the largest corner of the cell)",
        (
            "The arithmetic is the same as DISTIP and TOFIP: the inputs are the printed cell "
            "values, printed with 15 significant digits (relative 1e-15), so the agreement "
            "should be near 1e-13. 1e-6 is a loose bound on that arithmetic; it does not test "
            "the physics of the table, which is not an input here."
        ),
    ),
    "cell_index_mismatch_count": (
        "== 0",
        (
            "The cell index INT(4 theta_deg), INT(4 h_km) recomputed from the printed WEMIS and "
            "ZAPP must equal the printed indices. A mismatch means a value was printed from a "
            "different state than the one used for the interpolation."
        ),
    ),
    "bunch_count_positive": (
        "== 1 (at least one bunch on a telescope in the trace)",
        "Without bunches none of the consistency checks is evaluated.",
    ),
    "file_count": (
        "== 2 (65 and 75 deg traces)",
        "Both traces are needed; the input check ends the case with exit code 3 otherwise.",
    ),
}

HYPOTHESIS = (
    "theta_s = arccos(WEMIS) (vacuum zenith at the telescope equals the CORSIKA emission-angle "
    "cosine); h_e = ZEM (true emission altitude). Design note item 2; not verified."
)


# ---- trace parsing ------------------------------------------------------------------------


def _constant(token: str):
    """Value of a CONST line: a float, an integer or a logical (T, F)."""
    if token in ("T", "F"):
        return token == "T"
    return float(token)


def parse_bunch(line: str, where: str) -> dict:
    """One BUNCH line by its fixed columns (see the patch format 9611)."""
    if len(line) < BUNCH_OFFSET + FIELD_WIDTH * BUNCH_FIELDS:
        raise ValueError(f"{where}: BUNCH line too short ({len(line)} columns)")
    reals = [
        float(
            line[BUNCH_OFFSET + FIELD_WIDTH * i : BUNCH_OFFSET + FIELD_WIDTH * (i + 1)]
        )
        for i in range(BUNCH_FIELDS)
    ]
    return {
        "event": int(line[5:15]),
        "index": int(line[15:22]),
        "n2": int(line[22:28]),
        "m2": int(line[28:34]),
        "values": reals,
    }


def parse_trace(path) -> dict:
    """Constants, table rows and bunches of one trace file (raises ValueError when malformed)."""
    where = str(path)
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].split()[:2] != ["CURVTAB", FORMAT_VERSION]:
        raise ValueError(f"{where}: not a CURVTAB {FORMAT_VERSION} trace")
    constants: dict = {}
    rows: dict = {}
    bunches: list[dict] = []
    for number, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue
        if line.startswith(BUNCH_TAG):
            bunches.append(parse_bunch(line, f"{where}:{number}"))
        elif line.startswith("CONST"):
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"{where}:{number}: CONST needs name and value")
            constants[fields[1]] = _constant(fields[2])
        elif line.startswith("ROW"):
            fields = line.split()
            if len(fields) != 5:
                raise ValueError(f"{where}:{number}: ROW needs N, M, DISTEF, TOF")
            rows[(int(fields[1]), int(fields[2]))] = (
                float(fields[3]),
                float(fields[4]),
            )
        else:
            raise ValueError(f"{where}:{number}: unknown line {line[:20]!r}")
    missing = [name for name in REQUIRED_CONSTANTS if name not in constants]
    if missing:
        raise ValueError(f"{where}: missing CONST {', '.join(missing)}")
    return {"constants": constants, "rows": rows, "bunches": bunches}


def bunch_arrays(bunches: list[dict]) -> dict:
    """Column arrays of the bunches: scalars by name, cell values as (N, 4) arrays."""
    values = np.array([bunch["values"] for bunch in bunches], dtype=np.float64).reshape(
        -1, BUNCH_FIELDS
    )
    columns = {name: values[:, i] for i, name in enumerate(BUNCH_REALS)}
    n_real = len(BUNCH_REALS)
    columns["cell_distef"] = values[:, n_real : n_real + 4]
    columns["cell_tof"] = values[:, n_real + 4 : n_real + 8]
    columns["n2"] = np.array([bunch["n2"] for bunch in bunches], dtype=np.int64)
    columns["m2"] = np.array([bunch["m2"] for bunch in bunches], dtype=np.int64)
    columns["event"] = np.array([bunch["event"] for bunch in bunches], dtype=np.int64)
    return columns


# ---- interpolation as DISTIP and TOFIP do it -------------------------------------------------


def cell_index(wemis, zapp_cm, obs_cm):
    """(theta*4, height*4, N2, M2, OBS*4) exactly as DISTIP and TOFIP form them.

    theta*4 = ACOS(WEMIS) * FACT * 4 (degrees times 4), height*4 = ZAPP * 1e-5 * 4 (km times 4).
    N2 and M2 are INT() of these (truncation, the values are positive).
    """
    wemis = np.asarray(wemis, dtype=np.float64)
    zapp_cm = np.asarray(zapp_cm, dtype=np.float64)
    theta4 = np.arccos(wemis) * FORTRAN_FACT * 4.0
    height4 = zapp_cm * 1.0e-5 * 4.0
    obs4 = float(obs_cm) * 1.0e-5 * 4.0
    n2 = np.trunc(theta4).astype(np.int64)
    m2 = np.trunc(height4).astype(np.int64)
    return theta4, height4, n2, m2, obs4


def bilinear(theta4, height4, n2, m2, obs4, d_n2m2, d_n2m1, d_n1m2, d_n1m1):
    """The two-step linear interpolation of DISTIP (and TOFIP) on the four cell corners.

    Corner names: first index is the theta cell (N2 or N1 = N2 + 1), second the height cell
    (M2 or M1 = M2 + 1). A height below the observation level uses HM2 = OBS in the weights,
    as the Fortran does (HM2 .LT. OBS).
    """
    hm1 = (m2 + 1).astype(np.float64)
    hm2 = m2.astype(np.float64)
    hm2 = np.where(hm2 < obs4, obs4, hm2)
    thn1 = (n2 + 1).astype(np.float64)
    thn2 = n2.astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        slope = (d_n2m1 - d_n2m2) / (hm1 - hm2)
        at_n2 = slope * (height4 - hm2) + d_n2m2
        slope = (d_n1m1 - d_n1m2) / (hm1 - hm2)
        at_n1 = slope * (height4 - hm2) + d_n1m2
        slope = (at_n1 - at_n2) / (thn1 - thn2)
        return slope * (theta4 - thn2) + at_n2


def bilinear_bruteforce(theta4, height4, n2, m2, d00, d01, d10, d11):
    """Reference: the same interpolation written with explicit cell weights (no clamp).

    Equal to bilinear() for heights at or above the observation level (where HM2 = M2 and the
    grid spacing is 1 in units of 0.25). Used only by the tests.
    """
    u = theta4 - n2
    v = height4 - m2
    return (1 - u) * (1 - v) * d00 + u * (1 - v) * d10 + (1 - u) * v * d01 + u * v * d11


def relative_error(recomputed, printed, corners) -> np.ndarray:
    """|recomputed - printed| relative to the largest corner of the cell (0 when all are 0)."""
    scale = np.max(np.abs(corners), axis=1)
    difference = np.abs(np.asarray(recomputed) - np.asarray(printed))
    with np.errstate(divide="ignore", invalid="ignore"):
        result = np.where(scale > 0.0, difference / scale, difference)
    return np.asarray(result, dtype=np.float64)


def check_interpolation(columns: dict, obs_cm: float) -> dict:
    """Consistency of RDIST, TOFIP, CARTIM and the cell indices with the printed corners."""
    if columns["n2"].size == 0:
        return {
            "distip_relative_error": np.zeros(0),
            "tofip_relative_error": np.zeros(0),
            "cartim_relative_error": np.zeros(0),
            "cell_index_mismatch": 0,
            "recomputed_n2": np.zeros(0, dtype=np.int64),
            "recomputed_m2": np.zeros(0, dtype=np.int64),
        }
    theta4, height4, n2, m2, obs4 = cell_index(
        columns["wemis"], columns["zapp_cm"], obs_cm
    )
    cells = columns["cell_distef"]
    times = columns["cell_tof"]
    distip = bilinear(theta4, height4, n2, m2, obs4, *cells.T)
    tofip = bilinear(theta4, height4, n2, m2, obs4, *times.T)
    mismatch = int(np.count_nonzero((n2 != columns["n2"]) | (m2 != columns["m2"])))
    cartim = columns["temis_ns"] + columns["tofip_ns"]
    cartim_error = np.abs(columns["cartim_ns"] - cartim) / np.maximum(
        np.abs(columns["cartim_ns"]), 1e-300
    )
    return {
        "distip_relative_error": relative_error(distip, columns["rdist_cm"], cells),
        "tofip_relative_error": relative_error(tofip, columns["tofip_ns"], times),
        "cartim_relative_error": cartim_error,
        "cell_index_mismatch": mismatch,
        "recomputed_n2": n2,
        "recomputed_m2": m2,
    }


def _relative_to(printed: float, reference: float) -> float:
    """|printed - reference| relative to |reference| (absolute when the reference is 0)."""
    difference = abs(printed - reference)
    return difference / abs(reference) if reference != 0.0 else difference


def node_check(columns: dict, rows: dict) -> tuple[int, float]:
    """Cell corners that are printed grid nodes, compared with the 'ROW' lines of the trace.

    Returns (number of corners compared, largest relative difference of DISTEF and TOF).
    Corner order as in the bunch line: (N2,M2), (N2,M1), (N1,M2), (N1,M1).
    """
    compared, worst = 0, 0.0
    corners = ((0, 0), (0, 1), (1, 0), (1, 1))  # (theta offset, height offset)
    for k in range(columns["n2"].size):
        n2, m2 = int(columns["n2"][k]), int(columns["m2"][k])
        for corner, (dn, dm) in enumerate(corners):
            node = (n2 + dn, m2 + dm)
            if node not in rows:
                continue
            distef, tof = rows[node]
            compared += 1
            worst = max(
                worst,
                _relative_to(float(columns["cell_distef"][k, corner]), distef),
                _relative_to(float(columns["cell_tof"][k, corner]), tof),
            )
    return compared, worst


# ---- mapping comparison (hypothesis, details only) ---------------------------------------


def apparent_altitude_m(zapp_m, theta_rad, obs_m, earth_radius_m) -> float:
    """HE of INRTAB: the true altitude of the straight line with apparent angle THEAP at HEAPP."""
    dhelp = (zapp_m - obs_m) * math.tan(theta_rad)
    return -earth_radius_m + math.sqrt(dhelp * dhelp + (earth_radius_m + zapp_m) ** 2)


def mapping_residuals(columns: dict, constants: dict, nu) -> list[dict]:
    """Per bunch: obdeect (theta_s, h_e) against the CORSIKA table values; see the module docstring."""
    earth_m = float(constants["C1_CM"]) * 1.0e-2
    obs_m = float(constants["OBSLEV_CM"]) * 1.0e-2
    out = []
    for k in range(columns["wemis"].size):
        theta_rad = math.acos(min(1.0, float(columns["wemis"][k])))
        theta_s = math.degrees(theta_rad)
        h_e = float(columns["zem_cm"][k]) * 1.0e-2
        zapp_m = float(columns["zapp_cm"][k]) * 1.0e-2
        record = {
            "event": int(columns["event"][k]),
            "bunch": int(columns["index"][k]) if "index" in columns else k + 1,
            "theta_s_deg": theta_s,
            "h_e_km": h_e * 1.0e-3,
            "zapp_km": zapp_m * 1.0e-3,
            "below_observation": not h_e > obs_m,
        }
        if record["below_observation"]:
            out.append(record)
            continue
        b = (earth_m + obs_m) * math.sin(theta_rad)
        obd = refraction.ray_integrals(nu, earth_m, obs_m, h_e, b)
        rdist_m = float(columns["rdist_cm"][k]) * 1.0e-2
        straight_horizontal_m = (h_e - obs_m) * math.tan(theta_rad)
        straight_s_m = (h_e - obs_m) / math.cos(theta_rad)
        excess_corsika = float(columns["tofip_ns"][k]) - straight_s_m / C_M_PER_NS
        record.update(
            {
                "rdist_m": rdist_m,
                "straight_horizontal_m": straight_horizontal_m,
                "shift_obdeect_m": obd["shift_m"],
                "shift_residual_m": rdist_m - (straight_horizontal_m + obd["shift_m"]),
                "excess_corsika_ns": excess_corsika,
                "excess_obdeect_ns": obd["excess_ns"],
                "excess_residual_ns": excess_corsika - obd["excess_ns"],
                "altitude_residual_m": apparent_altitude_m(
                    zapp_m, theta_rad, obs_m, earth_m
                )
                - h_e,
            }
        )
        out.append(record)
    return out


def bin_residuals(records: list[dict]) -> list[dict]:
    """Counts and mean, standard deviation and largest magnitude per (zenith, height) cell.

    Bins are half-open [low, high); a value outside the edges is not in any bin.
    """
    usable = [r for r in records if not r["below_observation"]]
    zenith = np.array([r["theta_s_deg"] for r in usable], dtype=np.float64)
    height = np.array([r["h_e_km"] for r in usable], dtype=np.float64)
    shift = np.array([r["shift_residual_m"] for r in usable], dtype=np.float64)
    excess = np.array([r["excess_residual_ns"] for r in usable], dtype=np.float64)
    zi = np.searchsorted(ZENITH_EDGES_DEG, zenith, side="right") - 1
    hi = np.searchsorted(HEIGHT_EDGES_KM, height, side="right") - 1
    cells = []
    for i in range(len(ZENITH_EDGES_DEG) - 1):
        for j in range(len(HEIGHT_EDGES_KM) - 1):
            mask = (zi == i) & (hi == j)
            count = int(np.count_nonzero(mask))
            cell = {
                "zenith_deg": [ZENITH_EDGES_DEG[i], ZENITH_EDGES_DEG[i + 1]],
                "height_km": [HEIGHT_EDGES_KM[j], HEIGHT_EDGES_KM[j + 1]],
                "count": count,
            }
            for label, values in (
                ("shift_residual_m", shift),
                ("excess_residual_ns", excess),
            ):
                selected = values[mask]
                cell[f"{label}_mean"] = float(selected.mean()) if count else None
                cell[f"{label}_std"] = (
                    float(selected.std(ddof=1)) if count > 1 else None
                )
                cell[f"{label}_max_abs"] = (
                    float(np.max(np.abs(selected))) if count else None
                )
            cells.append(cell)
    return cells


# ---- per file and per run ----------------------------------------------------------------


def trace_path(root: Path, zenith_deg) -> Path:
    return root / TRACE_NAME.format(zenith=zenith_deg)


def analyse_file(path, zenith_deg, nu) -> dict:
    trace = parse_trace(path)
    constants = trace["constants"]
    columns = bunch_arrays(trace["bunches"])
    obs_cm = float(constants["OBSLEV_CM"])
    check = check_interpolation(columns, obs_cm)
    compared, node_worst = node_check(columns, trace["rows"])
    records = mapping_residuals(
        {**columns, "index": np.array([b["index"] for b in trace["bunches"]])},
        constants,
        nu,
    )
    interpolation = {
        "distip_max_relative_error": float(np.max(check["distip_relative_error"]))
        if columns["n2"].size
        else 0.0,
        "tofip_max_relative_error": float(np.max(check["tofip_relative_error"]))
        if columns["n2"].size
        else 0.0,
    }
    return {
        "file": Path(path).name,
        "zenith_deg": float(zenith_deg),
        "constants": constants,
        "row_count": len(trace["rows"]),
        "bunch_count": int(columns["n2"].size),
        "cell_index_mismatch_count": check["cell_index_mismatch"],
        "interpolation": interpolation,
        "interpolation_max_relative_error": max(interpolation.values()),
        "cartim_max_relative_error": float(np.max(check["cartim_relative_error"]))
        if columns["n2"].size
        else 0.0,
        "node_checks": compared,
        "node_max_relative_error": node_worst,
        "below_observation_count": sum(1 for r in records if r["below_observation"]),
        "bins": bin_residuals(records),
        "bunches": records,
    }


def summarise(files: list[dict]) -> dict:
    """Metrics of the case (numerical consistency and counts only)."""
    bunch_count = sum(row["bunch_count"] for row in files)
    return {
        "file_count": len(files),
        "bunch_count_positive": 1 if bunch_count > 0 else 0,
        "cell_index_mismatch_count": sum(
            row["cell_index_mismatch_count"] for row in files
        ),
        "table_interpolation_max_relative_error": max(
            (row["interpolation_max_relative_error"] for row in files), default=0.0
        ),
    }


def _plot_interpolation(files, path) -> None:
    """log10 of the interpolation errors per trace (DISTIP and TOFIP), against the limit."""
    from obdeect_validation.adapters import plotting

    fig, axis = plotting.figure()
    samples = []
    for i, row in enumerate(files):
        errors = [
            row["interpolation"]["distip_max_relative_error"],
            row["interpolation"]["tofip_max_relative_error"],
        ]
        samples.append(
            (
                f"{row['zenith_deg']:g} deg",
                np.log10(np.maximum(errors, 1e-16)),
                plotting.CYCLE[i % len(plotting.CYCLE)],
            )
        )
    plotting.histograms(
        axis,
        samples,
        "log10 of the max relative error (DISTIP, TOFIP per trace)",
        bins=20,
        clip=None,
    )
    axis.axvline(math.log10(INTERPOLATION_LIMIT), color=plotting.LIMIT, ls="--")
    plotting.save(fig, path)


def _plot_mapping(files, path) -> None:
    """Residuals against emission altitude, one colour per trace (hypothesis, no limits)."""
    from obdeect_validation.adapters import plotting

    fig, axes = plotting.figure(rows=1, columns=3, size=(13.5, 3.8))
    axes = np.asarray(axes).ravel()
    keys = (
        ("shift_residual_m", "RDIST - (straight + shift) [m]"),
        ("excess_residual_ns", "excess time residual [ns]"),
        ("altitude_residual_m", "HE(ZAPP) - ZEM [m]"),
    )
    for i, row in enumerate(files):
        records = [r for r in row["bunches"] if not r["below_observation"]]
        if not records:
            continue
        color = plotting.CYCLE[i % len(plotting.CYCLE)]
        height = [r["h_e_km"] for r in records]
        for axis, (key, label) in zip(axes, keys):
            axis.scatter(
                height,
                [r[key] for r in records],
                s=8,
                color=color,
                label=f"{row['zenith_deg']:g} deg",
            )
            axis.set(xlabel="emission altitude h_e [km]", ylabel=label)
            axis.axhline(0.0, color=plotting.ANALYTIC, lw=0.8)
    axes[0].legend()
    plotting.save(fig, path)


def main() -> int:
    root = Path(os.environ.get("VALIDATION_INPUTS", common.VALIDATION / "inputs"))
    traces = [trace_path(root, zenith) for zenith in ZENITHS_DEG]
    missing = [str(path) for path in traces if not path.is_file()]
    if missing:
        print(
            "missing input(s): "
            + ", ".join(missing)
            + "\nrun the CURVED CORSIKA trace (docs/curvature-corsika-tables.md, section 4) "
            "and copy each curvature-tables.dat to its curvature-tables-<zen>deg.dat name",
            file=sys.stderr,
        )
        return 3
    profile_file = refraction.profile_path()
    if not profile_file.is_file():
        print(f"missing input: atmosphere profile {profile_file}", file=sys.stderr)
        return 3
    nu = refraction.NuProfile(profile_file)
    files = [
        analyse_file(path, zenith, nu) for path, zenith in zip(traces, ZENITHS_DEG)
    ]
    metrics = summarise(files)
    details = {
        "hypothesis": HYPOTHESIS,
        "conventions": {
            "interpolation": "DISTIP/TOFIP arithmetic on the printed cell corners",
            "relative_error": "|recomputed - printed| / max |corner| of the same cell",
            "excess_corsika": "TOFIP - (h_e - h_obs)/(cos(theta_s) c), h_e = ZEM (tangent plane, as obdeect)",
            "shift_residual": "RDIST - ((h_e - h_obs) tan(theta_s) + shift_obdeect), signed shift",
            "obdeect_side": "numpy reference (refraction_curved_reference.ray_integrals)",
        },
        "limits": LIMITS,
        "profile": str(profile_file),
        "row_nodes": {
            "theta_deg": list(ROW_THETA_DEG),
            "height_km": list(ROW_HEIGHT_KM),
        },
        "bin_edges": {
            "zenith_deg": list(ZENITH_EDGES_DEG),
            "height_km": list(HEIGHT_EDGES_KM),
        },
        "files": files,
    }
    common.write_results(metrics, details)
    _plot_interpolation(files, Path("interpolation_error.png"))
    _plot_mapping(files, Path("mapping_residuals.png"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
