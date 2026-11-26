"""Plane-parallel and curved optical depths of sim_telarray transmission tables (pure numpy).

Independent expectation for obdeect's curvature modes. Nothing here imports obdeect.

Model (the plane-parallel part is sim_telarray's atmospheric_transmission2):
  tau(lambda, h): vertical optical depth from the table. Zero at and below H2, linear in log10 h
  between the table heights, constant above the last height.
  Plane-parallel (what sim_telarray applies for every zenith angle, also for CURVED input):
      T_p = exp(-(|tau(zem) - tau(h_t)| / |cz| + tau_focus)),
      tau_focus = tau(h_t + pathlen) - tau(h_t).
  Curved (obdeect curvature='corsika'): with p = (R + h_t) sin(theta_t) and
      cos theta(h) = sqrt(1 - p^2 / (R + h)^2),   cos theta_t = |cz|,
      tau_slant = integral from h_t to zem of (dtau/dh) / cos theta(h) dh,
      T_c = exp(-(tau_slant + tau_focus)).
  Since cos theta(h) >= cos theta_t for h >= h_t, T_c >= T_p for every bunch and wavelength.

Quadrature: on each table segment dtau/dh = a_j / (h ln 10), so tau_slant is a sum over segments of
a_j(lambda) G_j, where G_j = integral of dh / (h cos theta(h)) does not depend on the wavelength.
G_j is evaluated by composite Simpson in u = sqrt(h - h_t) (h = h_t + u^2), which removes the
(h - h_t)^(-1/2) behaviour of grazing rays. For cos theta_t = 0 the integrand is finite.

Wavelength: sim_telarray and obdeect take the integer nm floor(lambda) of lambda drawn with density
proportional to 1/lambda^2 on the band. The band expectation sums the exact probability of each
integer bin instead of sampling it, so it carries no Monte Carlo error. Bins at or above
MAX_LAMBDA_NM have zero transmission (sim_telarray skips them; the obdeect table ends at 999 nm).

Units: heights and positions in m a.s.l., wavelengths in nm, optical depths dimensionless.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from obdeect_validation.adapters.atmosphere_chain import MAX_LAMBDA_NM, table_lines

EARTH_RADIUS_M = 6371315.0  # obdeect default (CORSIKA CURVED constant)
SIMPSON_INTERVALS = 256  # composite Simpson intervals per table segment (even)
CHUNK = 128  # bunches per block of the segment quadrature (memory)
LN10 = math.log(10.0)


@dataclass(frozen=True)
class Table:
    """Parsed transmission table: H2 and the table heights (m), wavelengths (nm), tau (nw, nh)."""

    h2_m: float
    heights_m: np.ndarray
    wavelengths: np.ndarray
    tau: np.ndarray

    @property
    def edges_m(self) -> np.ndarray:
        """Table nodes from H2 up to the last height (m); tau is zero at the first node."""
        return np.concatenate(([self.h2_m], self.heights_m))


def parse_table(text: str) -> Table:
    """Parse a sim_telarray transmission table ('# H2= ..., H1= ...' header and integer rows).

    Rows must be consecutive integer wavelengths with one value per height. Values are kept as
    written (the sentinel 99999 is a large optical depth, as sim_telarray reads it).
    """
    header, rows = table_lines(text)
    line = next((ln for ln in header if "H2=" in ln and "H1=" in ln), None)
    if line is None:
        raise ValueError("transmission table has no 'H2=' / 'H1=' header")
    h2_text, h1_text = line.split("H2=", 1)[1].split("H1=", 1)
    h2_km = float(h2_text.strip().split()[0].rstrip(","))
    heights_km = np.array([float(value) for value in h1_text.split()], dtype=float)
    if (
        heights_km.size == 0
        or np.any(np.diff(heights_km) <= 0.0)
        or heights_km[0] <= h2_km
    ):
        raise ValueError("table heights must be increasing and above H2")
    if not rows:
        raise ValueError("transmission table has no wavelength rows")
    wavelengths = np.array([wavelength for wavelength, _ in rows], dtype=int)
    if np.any(np.diff(wavelengths) != 1):
        raise ValueError("table wavelength rows must be consecutive integers (nm)")
    tau = np.array([[float(v) for v in fields[1:]] for _, fields in rows], dtype=float)
    if tau.shape[1] != heights_km.size:
        raise ValueError("table row length does not match the number of heights")
    return Table(
        h2_m=h2_km * 1000.0,
        heights_m=heights_km * 1000.0,
        wavelengths=wavelengths,
        tau=tau,
    )


def _log_nodes(table: Table) -> tuple[np.ndarray, np.ndarray]:
    """log10 of the table nodes (edges_m) and tau at the nodes, shape (n_wavelengths, n_nodes)."""
    x = np.log10(table.edges_m)
    nodes = np.concatenate([np.zeros((table.tau.shape[0], 1)), table.tau], axis=1)
    return x, nodes


def vertical_depth(table: Table, heights_m) -> np.ndarray:
    """Vertical optical depth tau(lambda, h) of the table, shape (n_wavelengths, n_heights)."""
    heights = np.atleast_1d(np.asarray(heights_m, dtype=float))
    if not np.all(np.isfinite(heights)) or np.any(heights < 0.0):
        raise ValueError("heights must be finite and not negative")
    x, nodes = _log_nodes(table)
    lh = np.log10(np.maximum(heights, 1.0))
    j = np.clip(np.searchsorted(x, lh, side="right") - 1, 0, x.size - 2)
    t = np.clip((lh - x[j]) / (x[j + 1] - x[j]), 0.0, 1.0)
    return nodes[:, j] + t * (nodes[:, j + 1] - nodes[:, j])


def segment_slopes(table: Table) -> np.ndarray:
    """a_j(lambda) = d tau / d log10(h) on each table segment, shape (n_wavelengths, n_heights)."""
    x, nodes = _log_nodes(table)
    return np.diff(nodes, axis=1) / np.diff(x)


def _as_bunches(zem_m, cz, telescope_altitude_m):
    zem = np.atleast_1d(np.asarray(zem_m, dtype=float)).ravel()
    c = np.atleast_1d(np.asarray(cz, dtype=float)).ravel()
    if zem.shape != c.shape:
        raise ValueError("zem and cz have different lengths")
    if not (np.all(np.isfinite(zem)) and np.all(np.isfinite(c))):
        raise ValueError("bunch emission heights and direction cosines must be finite")
    if np.any(c < 0.0) or np.any(c > 1.0):
        raise ValueError("|cz| must lie in [0, 1]")
    if np.any(zem <= telescope_altitude_m):
        raise ValueError("emission at or below the telescope altitude is not supported")
    return zem, c


def _geometry_block(edges, zem, c, h_t, radius, intervals):
    """G[b, j] = integral of dh / (h cos theta) over segment j clipped to [h_t, zem], one block."""
    top = zem[:, None]
    lo = np.clip(edges[:-1], h_t, top)
    hi = np.clip(edges[1:], h_t, top)
    su = np.sqrt(lo - h_t)
    sv = np.sqrt(hi - h_t)
    step = (sv - su) / intervals
    s = su[..., None] + step[..., None] * np.arange(intervals + 1)
    h = h_t + s * s
    cc = c[:, None, None]
    rel = 2.0 * radius + 2.0 * h_t + s * s  # 2R + h + h_t
    # cos^2 theta (R+h)^2 = (h - h_t)(2R + h + h_t) + cz^2 (R + h_t)^2 with h - h_t = s^2
    q = s * s * rel + cc * cc * (radius + h_t) ** 2
    grazing = cc > 0.0
    # integrand in s: 2 s (R + h) / (h sqrt(q)); s / sqrt(q) = 1 / sqrt(rel) when cz = 0
    ratio = np.where(
        grazing,
        s / np.sqrt(np.where(grazing, q, 1.0)),
        1.0 / np.sqrt(rel),
    )
    f = 2.0 * (radius + h) / h * ratio
    weights = np.full(intervals + 1, 2.0)
    weights[1::2] = 4.0
    weights[0] = weights[-1] = 1.0
    integral = step * (f @ weights) / 3.0
    return np.where(hi > lo, integral, 0.0)


def _geometry_integrals(edges, zem, c, h_t, radius, intervals):
    out = np.zeros((zem.size, edges.size - 1))
    for start in range(0, zem.size, CHUNK):
        part = slice(start, start + CHUNK)
        out[part] = _geometry_block(edges, zem[part], c[part], h_t, radius, intervals)
    return out


def curved_optical_depth(
    table: Table,
    zem_m,
    cz,
    telescope_altitude_m: float,
    radius_m: float = EARTH_RADIUS_M,
    intervals: int = SIMPSON_INTERVALS,
) -> np.ndarray:
    """Curved slant optical depth tau_slant(lambda) per bunch, shape (n_bunches, n_wavelengths).

    cz is |cz| of the bunch (cos of the zenith angle at the telescope, as obdeect takes it).
    Valid for zem above the telescope; cz = 0 is finite.
    """
    if intervals < 2 or intervals % 2:
        raise ValueError("intervals must be an even integer of at least 2")
    if not radius_m > 1.0e6:
        raise ValueError("Earth radius must be above 1e6 m")
    zem, c = _as_bunches(zem_m, cz, telescope_altitude_m)
    geometry = _geometry_integrals(
        table.edges_m, zem, c, telescope_altitude_m, radius_m, intervals
    )
    return geometry @ segment_slopes(table).T / LN10


def plane_optical_depth(
    table: Table, zem_m, cz, telescope_altitude_m: float
) -> np.ndarray:
    """Plane-parallel slant optical depth |tau(zem) - tau(h_t)| / |cz|, shape (bunches, wavelengths)."""
    zem, c = _as_bunches(zem_m, cz, telescope_altitude_m)
    if np.any(c <= 0.0):
        raise ValueError("plane-parallel airmass needs |cz| > 0")
    tau_t = vertical_depth(table, [telescope_altitude_m])[:, 0]
    tau_z = vertical_depth(table, zem)
    return (np.abs(tau_z - tau_t[:, None]) / c[None, :]).T


def focus_optical_depth(
    table: Table, telescope_altitude_m: float, pathlen_m: float
) -> np.ndarray:
    """Optical depth inside the telescope tau(h_t + pathlen) - tau(h_t), shape (wavelengths,)."""
    if not (np.isfinite(pathlen_m) and pathlen_m >= 0.0):
        raise ValueError("pathlen must be finite and not negative")
    top = vertical_depth(table, [telescope_altitude_m + pathlen_m])[:, 0]
    base = vertical_depth(table, [telescope_altitude_m])[:, 0]
    return top - base


def transmission_pair(
    table: Table,
    zem_m,
    cz,
    telescope_altitude_m: float,
    pathlen_m: float = 0.0,
    radius_m: float = EARTH_RADIUS_M,
) -> tuple[np.ndarray, np.ndarray]:
    """(curved, plane-parallel) transmission, each of shape (bunches, wavelengths)."""
    focus = focus_optical_depth(table, telescope_altitude_m, pathlen_m)
    curved = curved_optical_depth(table, zem_m, cz, telescope_altitude_m, radius_m)
    plane = plane_optical_depth(table, zem_m, cz, telescope_altitude_m)
    return np.exp(-(curved + focus)), np.exp(-(plane + focus))


def band_bins(band) -> tuple[np.ndarray, np.ndarray]:
    """Integer nm bins floor(lambda) of lambda ~ 1/lambda^2 on [lo, hi], and their probabilities."""
    lo, hi = float(band[0]), float(band[1])
    if not 0.0 < lo < hi:
        raise ValueError("wavelength band must satisfy 0 < lo < hi")
    n = np.arange(math.floor(lo), math.ceil(hi))
    a = np.maximum(n, lo)
    b = np.minimum(n + 1, hi)
    keep = b > a
    probability = (1.0 / a - 1.0 / b) / (1.0 / lo - 1.0 / hi)
    return n[keep], probability[keep]


def band_expectation(
    table: Table,
    band,
    zem_m,
    cz,
    telescope_altitude_m: float,
    pathlen_m: float = 0.0,
    radius_m: float = EARTH_RADIUS_M,
    chunk: int = 2048,
    weight=None,
) -> dict[str, np.ndarray]:
    """Expected transmission per bunch, averaged exactly over the band, and its variance.

    weight: optional detection weight w(nm) (callable on an array of integer nm, e.g. quantum
    efficiency times mirror reflectivity). The expectation is then the w-weighted mean of T:
    sum_n P_n w_n T_n / sum_n P_n w_n over all integer bins of the band. Bins at or above
    MAX_LAMBDA_NM count in the normalisation with their weight and contribute zero transmission;
    live_fraction reports the share of the normalisation they do not carry (1 when w vanishes there).

    Returns curved, curved_variance, plane, plane_variance (each shape (bunches,)) and live_fraction.
    The variance is the variance of w T over the wavelength draw of one bunch, normalised by the
    square of the normalisation (so that photon-weighted sums give the standard error).
    """
    zem, c = _as_bunches(zem_m, cz, telescope_altitude_m)
    bins, probability = band_bins(band)
    if weight is None:
        detection = np.ones(bins.size)
    else:
        detection = np.asarray(weight(bins), dtype=float)
        if detection.shape != bins.shape or not np.all(np.isfinite(detection)):
            raise ValueError("detection weight must be finite on every band bin")
        if np.any(detection < 0.0):
            raise ValueError("detection weight must not be negative")
    weighted = probability * detection
    norm = float(weighted.sum())
    if not norm > 0.0:
        raise ValueError("detection weight vanishes over the whole band")
    live = bins < MAX_LAMBDA_NM
    out = {
        key: np.zeros(zem.size)
        for key in ("curved", "curved_variance", "plane", "plane_variance")
    }
    out["live_fraction"] = float(weighted[live].sum() / norm)
    if not np.any(live):
        return out
    bins_live = bins[live]
    pw = weighted[live]
    pw2 = (probability * detection * detection)[live]
    rows = (
        np.clip(bins_live, table.wavelengths[0], table.wavelengths[-1])
        - table.wavelengths[0]
    )
    sub = Table(table.h2_m, table.heights_m, bins_live, table.tau[rows])
    for start in range(0, zem.size, chunk):
        part = slice(start, start + chunk)
        curved, plane = transmission_pair(
            sub, zem[part], c[part], telescope_altitude_m, pathlen_m, radius_m
        )
        for name, values in (("curved", curved), ("plane", plane)):
            total = values @ pw
            out[name][part] = total / norm
            out[f"{name}_variance"][part] = np.maximum(
                (values * values) @ pw2 - total * total, 0.0
            ) / (norm * norm)
    return out


def unit_variance(band, weight=None) -> float:
    """Variance (normalised as in band_expectation) of the wavelength draw of a unit transmission.

    This is the draw error of a count that carries the detection weight alone, e.g. a unit-
    transmission sim_telarray run. Zero when the weight is constant over the band.
    """
    bins, probability = band_bins(band)
    detection = (
        np.ones(bins.size) if weight is None else np.asarray(weight(bins), dtype=float)
    )
    norm = float((probability * detection).sum())
    if not norm > 0.0:
        raise ValueError("detection weight vanishes over the whole band")
    second = float((probability * detection * detection).sum())
    return max(second - norm * norm, 0.0) / (norm * norm)


def _photon_total(photons) -> tuple[np.ndarray, float]:
    """Photons as float array and their positive total."""
    photons = np.asarray(photons, dtype=float)
    total = float(photons.sum())
    if not total > 0.0:
        raise ValueError("photon total must be positive")
    return photons, total


def weighted_mean(values, photons) -> float:
    """Photon-weighted mean of per-bunch values."""
    values = np.asarray(values, dtype=float)
    photons, total = _photon_total(photons)
    return float(np.sum(photons * values) / total)


def weighted_standard_error(variances, photons) -> float:
    """Standard error of the photon-weighted mean from the per-bunch wavelength variances."""
    variances = np.asarray(variances, dtype=float)
    photons, total = _photon_total(photons)
    return float(math.sqrt(np.sum(photons * photons * variances)) / total)


def ratio_with_error(numerator, numerator_error, denominator, denominator_error):
    """Ratio of two quantities and its error from independent relative errors (first order)."""
    if not (numerator > 0.0 and denominator > 0.0):
        raise ValueError("ratio needs positive numerator and denominator")
    ratio = numerator / denominator
    relative = math.hypot(numerator_error / numerator, denominator_error / denominator)
    return ratio, ratio * relative
