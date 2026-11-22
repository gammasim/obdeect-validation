"""References and stream helpers for the upward-photon propagation cases (plane mode, phase 2).

Everything here is numpy and scipy, written independently of obdeect code:

- Profile: n - 1 from the atmprof table (fixtures/atmprof_ecmwf_north_winter_fixed.dat, a copy of
  obdeect's test data src/propagation/cpp/tests/data, "Do not edit" there; copied unchanged). ln(n - 1) is a natural cubic spline in
  height (scipy CubicSpline, bc_type "natural"), exactly as documented for obdeect's interpolation
  (docs/propagation/refraction.md, "Interpolation"); obdeect tabulates it on a 3 m grid, which differs
  by about 1e-8 relative.
- Snell integrals: along n(h) sin(theta) = K with K = n_e sin(theta_e), the horizontal displacement
  integral of tan(theta(h)) and the travel time integral of n(h) / cos(theta(h)) / c, by scipy quad.
- Transmission: tau(lambda, h) from the extinction table (log10-linear in height, zero at H2, linear in
  wavelength between integer rows, docs/propagation/extinction.md).

Direction convention: every photon is written with its float32 cx, cy, cz exactly as stored, and the
references read those values back (as in refraction_curved_reference). cz is always written, because
upward and horizontal photons need it.
"""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

from obdeect_validation.adapters import endtoend_common as e2e

ROOT = Path(__file__).resolve().parents[3]
PROFILE = ROOT / "fixtures/atmprof_ecmwf_north_winter_fixed.dat"
TRANSMISSION = ROOT / "fixtures/atm_trans_2000_1_10_0_0_2000.dat"
SITE_M = 2000.0
C_M_NS = 0.299792458
QUAD = {"limit": 800, "epsabs": 1e-12, "epsrel": 1e-11}


class Profile:
    """n - 1 against height [m a.s.l.] with the documented interpolation (docs/propagation/refraction.md).

    ln(n - 1) is a natural cubic spline through the levels (scipy CubicSpline). It is tabulated at nodes
    that include every level and are at most GRID_M apart; between nodes ln(n - 1) is linear (n - 1 is
    exponential inside a cell), zero above the table and the bottom value below it.
    """

    GRID_M = 3.0

    def __init__(self, path=PROFILE):
        rows = []
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            fields = line.split()
            if fields and not line.lstrip().startswith("#"):
                rows.append([float(v) for v in fields[:4]])
        table = np.asarray(rows, dtype=np.float64)
        levels = table[:, 0] * 1000.0
        spline = CubicSpline(levels, np.log(table[:, 3]), bc_type="natural")
        nodes = [levels[0]]
        for a_m, b_m in itertools.pairwise(levels):
            pieces = math.ceil((b_m - a_m) / self.GRID_M)
            nodes.extend(np.linspace(a_m, b_m, pieces + 1)[1:])
        self.grid = np.asarray(nodes)
        self.lnnu = spline(self.grid)
        self.altitude = levels
        self.top = float(levels[-1])

    def nu(self, height_m):
        h = np.asarray(height_m, dtype=np.float64)
        value = np.exp(np.interp(h, self.grid, self.lnnu))
        return np.where(h > self.top, 0.0, value)


PROFILE_REF = None


def profile() -> Profile:
    global PROFILE_REF
    if PROFILE_REF is None:
        PROFILE_REF = Profile()
    return PROFILE_REF


def snell_integrals(prof: Profile, h_a: float, sin_a: float, h_b: float) -> dict:
    """Ray with n(h) sin(theta) = K = n(h_a) sin_a, from h_a to h_b (either order).

    horizontal: integral of tan(theta(h)) over the height interval [m] (the lateral travel, a positive
    length along the horizontal direction of motion); time_ns: integral of n / cos(theta) / c;
    sin_b: sin(theta) at h_b (Snell). Turning (sin >= 1 inside the interval) is reported, not integrated.
    """
    k = (1.0 + float(prof.nu(h_a))) * sin_a
    lo, hi = min(h_a, h_b), max(h_a, h_b)

    def sin_h(h):
        return k / (1.0 + float(prof.nu(h)))

    def tan_h(h):
        s = sin_h(h)
        return s / math.sqrt(1.0 - s * s)

    def time_h(h):
        s = sin_h(h)
        return (1.0 + float(prof.nu(h))) / math.sqrt(1.0 - s * s)

    turns = sin_h(hi) >= 1.0 or sin_h(lo) >= 1.0
    if turns:
        return {
            "turns": True,
            "horizontal": math.nan,
            "time_ns": math.nan,
            "sin_b": sin_h(h_b),
        }
    horizontal = quad(tan_h, lo, hi, **QUAD)[0]
    time = quad(time_h, lo, hi, **QUAD)[0] / C_M_NS
    return {
        "turns": False,
        "horizontal": horizontal,
        "time_ns": time,
        "sin_b": sin_h(h_b),
    }


def critical_sin(prof: Profile, h_e: float, h_t: float) -> float:
    """sin(theta_e) of the upward ray at which sin(theta) reaches 1 at h_t: (1 + nu_t) / (1 + nu_e)."""
    return (1.0 + float(prof.nu(h_t))) / (1.0 + float(prof.nu(h_e)))


# ------------------------------------------------------------------------- transmission


class Transmission:
    """tau(lambda, h) from the extinction table (log10-linear in h, linear in wavelength)."""

    def __init__(self, path=TRANSMISSION):
        lines = Path(path).read_text(encoding="utf-8").splitlines()
        header = next(line for line in lines if "H2=" in line)
        values = header.split("H1=")[1].split()
        self.h1 = np.array([float(v) for v in values]) * 1000.0
        self.h2 = float(header.split("H2=")[1].split(",")[0]) * 1000.0
        self.rows = {}
        for line in lines:
            fields = line.split()
            if fields and not line.lstrip().startswith("#"):
                self.rows[round(float(fields[0]))] = np.array(
                    [float(v) for v in fields[1:]]
                )
        self.log_nodes = np.log10(np.concatenate([[self.h2], self.h1]))

    def tau_row(self, row: int, height_m) -> np.ndarray:
        """Optical depth of one integer-nm row at heights [m]; 0 at and below H2."""
        nodes = np.concatenate([[0.0], self.rows[row]])
        h = np.asarray(height_m, dtype=np.float64)
        log_h = np.log10(np.maximum(h, self.h2))
        value = np.interp(
            log_h, self.log_nodes, nodes
        )  # constant above the last height
        return np.where(h <= self.h2, 0.0, value)

    def tau(self, wavelength_nm: float, height_m) -> np.ndarray:
        n = math.floor(wavelength_nm)
        f = wavelength_nm - n
        return (1.0 - f) * self.tau_row(n, height_m) + f * self.tau_row(n + 1, height_m)

    def dtau_dh(self, wavelength_nm: float, height_m: float) -> float:
        """d tau / dh [1/m] at height_m (above H2): the log10-linear segment slope over h ln 10."""
        n = math.floor(wavelength_nm)
        f = wavelength_nm - n
        return (1.0 - f) * self._row_slope(n, height_m) + f * self._row_slope(
            n + 1, height_m
        )

    def _row_slope(self, row: int, height_m: float) -> float:
        h = float(height_m)
        if h <= self.h2 or h >= self.h1[-1]:
            return 0.0
        nodes = np.concatenate([[0.0], self.rows[row]])
        log_h = math.log10(h)
        i = int(np.searchsorted(self.log_nodes, log_h, side="right")) - 1
        i = min(max(i, 0), len(self.log_nodes) - 2)
        slope = (nodes[i + 1] - nodes[i]) / (self.log_nodes[i + 1] - self.log_nodes[i])
        return slope / (h * math.log(10.0))

    def breakpoints(self, low_m: float, high_m: float) -> list:
        """Table heights strictly inside (low, high), for the quadrature."""
        return [
            float(v) for v in np.concatenate([[self.h2], self.h1]) if low_m < v < high_m
        ]


# ------------------------------------------------------------------------- directions


def stored_direction(zenith_deg, azimuth_deg):
    """Unit direction of travel (upward, cz > 0) as float32 cx, cy, cz, as stored."""
    th = math.radians(zenith_deg)
    az = math.radians(azimuth_deg)
    return (
        np.float32(math.sin(th) * math.cos(az)),
        np.float32(math.sin(th) * math.sin(az)),
        np.float32(math.cos(th)),
    )


def path_to_plane(z_emission: float, plane_z: float, cz_abs: float) -> float:
    """Straight path [m] from the emission (sky z) to the plane, along a direction with |cz|."""
    return abs(plane_z - z_emission) / cz_abs


# ------------------------------------------------------------------------- streams


def write_stream(
    path,
    telescopes,
    blocks,
    *,
    plane_site_m=SITE_M,
    refraction_applied=0,
    curvature_applied=None,
    wavelength_nm=(300.0, 700.0),
    chunk_rows=None,
    indices=None,
    run_extra=None,
):
    """Write a sky-frame stream with one telescope block per entry of blocks.

    telescopes: list of (x, y, z) [m] of the telescope (z above the observation level), telescope ids
    1..n. blocks: list of (telescope_id, columns) with float32 columns x y z cx cy cz t w wavelength (any
    subset of cz/wavelength allowed). chunk_rows: rows per chunk (default: one chunk per block).
    indices: optional per-row photon index of each block (same length as the columns), written as one
    chunk per photon with first_index = that index (the photon keeps its identity when rows are
    shuffled). Without it the rows are numbered 0..n-1 in file order.
    """
    from obdeect.stream import Writer

    run = {
        "run_number": 1,
        "observation_level_m": plane_site_m,
        "telescope_x_m": [float(t[0]) for t in telescopes],
        "telescope_y_m": [float(t[1]) for t in telescopes],
        "telescope_z_m": [float(t[2]) for t in telescopes],
        "refraction_applied": refraction_applied,
        "wavelength_min_nm": wavelength_nm[0],
        "wavelength_max_nm": wavelength_nm[1],
    }
    if curvature_applied is not None:
        run["curvature_applied"] = curvature_applied
    run.update(run_extra or {})
    with Writer(str(path), "obdeect validation: upward propagation") as writer:
        writer.write("run", run)
        writer.write("event", {"event_number": 1})
        for telescope_id, columns in blocks:
            rows = len(columns["x"])
            writer.write("telescope", {"telescope_id": telescope_id, "frame": 0})
            if indices is not None:
                for row, index in enumerate(indices[telescope_id]):
                    piece = {
                        k: np.asarray(v)[row : row + 1] for k, v in columns.items()
                    }
                    writer.write_photons(telescope_id, 0, int(index), piece)
            else:
                step = chunk_rows or rows
                for first in range(0, rows, step):
                    piece = {
                        k: np.asarray(v)[first : first + step]
                        for k, v in columns.items()
                    }
                    writer.write_photons(telescope_id, 0, first, piece)
            writer.write(
                "telescope_end", {"telescope_id": telescope_id, "photons": float(rows)}
            )
        writer.write("event_end", {})
        writer.write("run_end", {})
    return path


def read_rows(path) -> dict:
    """Photon rows of a stream (photon_rows): event, telescope, index and float32 columns."""
    return e2e.photon_rows(path)


def ulp32(values) -> np.ndarray:
    """float32 spacing at the magnitude of each value, floored at 1 (a metre or a ns): the rounding unit."""
    magnitude = np.maximum(np.abs(np.asarray(values, dtype=np.float64)), 1.0)
    return np.spacing(magnitude.astype(np.float32)).astype(np.float64)
