"""Emission profiles of obdeect sources against their analytic distributions, and emission levels.

Part 1, profiles. Each profile is one point source at the origin without target (every photon
has weight 1, so the sample is the profile itself), N = 2e5 photons, axis +z. Every profile
is compared with its analytic distribution:

- angular, sample u = 1 - cos(theta) = 1 - cz (the profile is defined in u, see
  docs/sources/spec/distributions.md): isotropic (u uniform on [0, 2]), lambertian
  (pdf 1 - u on [0, 1]), filled_cone (u uniform on [0, u0]), gauss (intensity
  exp(-theta^2 / 2 sigma^2) per solid angle, truncated at min(5 sigma, pi); the reference is the
  exact Gaussian integral, the code uses a piecewise-linear table of it), rayleigh (intensity
  1 + cos^2 theta, normalised over u in [0, 2]), table (piecewise-linear intensity over the
  angle nodes, interpolated in u as the code does);
- spectral, sample the wavelength: uniform, inverse_square (pdf ~ 1 / lambda^2),
  table (piecewise-linear density);
- pulse, sample the emission time t: tophat (uniform), gauss, exponential (leading edge at
  mean - tau), complex (flat top + Gaussian + exponential tail, convolved), table
  (piecewise-linear rate).

Degenerate profiles (fixed angle "cone", fixed wavelength, fixed time) have no density; they are
checked for exact values instead. Every non-degenerate profile is also checked for azimuth
uniformity about the axis (angular profiles only).

Statistics. KS: scipy.stats.kstest against the analytic CDF, unweighted (all weights are 1 here,
asserted). Gate ``ks_fail_count`` = number of KS tests with p < 1e-4. Justification: 21 KS tests
(14 profile KS, 7 azimuth KS) with a fixed seed; at p < 1e-4 per test the family-wise false alarm
for a fresh seed is 21 x 1e-4 = 2.1e-3 < 0.3 %. Moments: sample mean and variance against the
analytic values; pull = difference / standard error, with SE(mean) = sigma / sqrt(N) and
SE(var) = sqrt((m4 - s^4) / N) from the sample fourth central moment. Gate ``moment_pull_max``
(the largest of the 2 x 14 pulls) <= 5 (k = 5 standard errors: P(|z| > 5) = 5.7e-7 per pull,
28 pulls: family-wise 1.6e-5). The metric is gated with expected 0 and atol 5, since pulls are
non-negative.

Control. The isotropic sample tested against the Lambertian CDF must be rejected (p < 1e-4):
``control_fail_detected`` = 1 shows that the test has the power to see a wrong profile.

Part 2, emission levels. The same point source (isotropic, at (0, 0, 40) m, uniform spectrum,
tophat pulse) written as a run with frame sky, telescope (T), camera (C) and ground (G), one
telescope at the origin of radius 12.5 m, same seed. The sources module has no frame transform:
the frame is a label, so C and T rows must equal the sky rows exactly (the requirement is
positions to 1e-5 m and directions to 1e-6). The G rows are the photons whose path enters the
telescope sphere: they must be the sky rows that hit the sphere, with unchanged directions,
entry points P + s d (smallest positive root), |r| = R, time t + s / c, and zem = observation
level + telescope z + P.z. Tolerances: positions 1e-5 m (about 10 float32 ulp at 12.5 m, and
direction rounding times the path of 28 m gives 2e-6 m), directions 1e-6 (float32 rounding of
unit components is 6e-8), zem 2 float32 ulp at 2190 m (4.9e-4 m), time 2 float32 ulp of the
value, row counts and frame codes exact.
"""

from __future__ import annotations

import math
import shutil
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.special import ndtr

from obdeect_validation.adapters import common

SEED = 20261123
PHOTONS = 200_000
P_THRESHOLD = 1e-4  # per KS test (see the module docstring)
MOMENT_PULL_LIMIT = 5.0
GRID_POINTS = 2_000_001
OBSERVATION_M = 2150.0
TELESCOPE_RADIUS_M = 12.5
SOURCE_POSITION_M = (0.0, 0.0, 40.0)
C_M_PER_NS = 0.299792458
COS_CONE_DEG = 30.0
DECIMAL_NM = 400.0


# ------------------------------------------------------------------ reference distributions


def grid_moments(lo: float, hi: float, pdf: Callable) -> tuple[float, float]:
    """Mean and variance of an unnormalised density on [lo, hi] (trapezoid on a fine grid)."""
    x = np.linspace(lo, hi, GRID_POINTS)
    f = pdf(x)
    z = np.trapezoid(f, x)
    mean = np.trapezoid(x * f, x) / z
    return float(mean), float(np.trapezoid((x - mean) ** 2 * f, x) / z)


def grid_cdf(lo: float, hi: float, pdf: Callable) -> Callable:
    """CDF of an unnormalised density on [lo, hi] by cumulative trapezoid, interpolated."""
    x = np.linspace(lo, hi, GRID_POINTS)
    f = pdf(x)
    cumulative = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))))
    cumulative /= cumulative[-1]

    def cdf(v):
        return np.interp(np.asarray(v, dtype=float), x, cumulative, left=0.0, right=1.0)

    return cdf


def piecewise_linear_cdf(nodes_x, nodes_y) -> Callable:
    """Exact CDF of the piecewise-linear density through (nodes_x, nodes_y) (zero outside)."""
    x = np.asarray(nodes_x, dtype=float)
    y = np.asarray(nodes_y, dtype=float)
    segment_area = 0.5 * (y[1:] + y[:-1]) * np.diff(x)
    cumulative = np.concatenate(([0.0], np.cumsum(segment_area)))
    total = cumulative[-1]

    def cdf(v):
        v = np.asarray(v, dtype=float)
        i = np.clip(np.searchsorted(x, v, side="right") - 1, 0, len(x) - 2)
        width = x[i + 1] - x[i]
        s = np.clip(v - x[i], 0.0, width)
        slope = (y[i + 1] - y[i]) / width
        value = cumulative[i] + y[i] * s + 0.5 * slope * s * s
        value = np.where(v <= x[0], 0.0, value)
        value = np.where(v >= x[-1], total, value)
        return value / total

    def pdf(v):
        return np.interp(np.asarray(v, dtype=float), x, y, left=0.0, right=0.0)

    return cdf, pdf


@dataclass(frozen=True)
class Reference:
    """Analytic distribution of one sampled quantity."""

    cdf: Callable
    mean: float
    var: float
    lo: float
    hi: float


def ref_isotropic() -> Reference:
    return Reference(
        lambda u: np.clip(np.asarray(u) / 2.0, 0.0, 1.0), 1.0, 1.0 / 3.0, 0.0, 2.0
    )


def ref_lambertian() -> Reference:
    def cdf(u):
        v = np.clip(np.asarray(u, dtype=float), 0.0, 1.0)
        return 2.0 * v - v * v

    # pdf (1 - u) / (1/2) on [0, 1]: mean 1/3, E[u^2] = 1/6, variance 1/18.
    return Reference(cdf, 1.0 / 3.0, 1.0 / 18.0, 0.0, 1.0)


def ref_filled_cone(angle_deg: float) -> Reference:
    u0 = 1.0 - math.cos(math.radians(angle_deg))
    return Reference(
        lambda u: np.clip(np.asarray(u, dtype=float) / u0, 0.0, 1.0),
        u0 / 2.0,
        u0 * u0 / 12.0,
        0.0,
        u0,
    )


def ref_gauss(sigma_deg: float) -> Reference:
    sigma = math.radians(sigma_deg)
    theta_max = min(5.0 * sigma, math.pi)
    u_max = 1.0 - math.cos(theta_max)

    def pdf(u):
        theta = 2.0 * np.arcsin(
            np.sqrt(np.clip(np.asarray(u, dtype=float) / 2.0, 0.0, 1.0))
        )
        return np.exp(-(theta * theta) / (2.0 * sigma * sigma))

    mean, var = grid_moments(0.0, u_max, pdf)
    return Reference(grid_cdf(0.0, u_max, pdf), mean, var, 0.0, u_max)


def ref_rayleigh() -> Reference:
    def cdf(u):
        v = np.clip(np.asarray(u, dtype=float), 0.0, 2.0)
        return (3.0 * v + 1.0 - (1.0 - v) ** 3) / 8.0

    def pdf(u):
        return 1.0 + (1.0 - np.asarray(u, dtype=float)) ** 2

    mean, var = grid_moments(0.0, 2.0, pdf)
    return Reference(cdf, mean, var, 0.0, 2.0)


def ref_angular_table(angles_deg, values) -> Reference:
    u_nodes = 1.0 - np.cos(np.radians(np.asarray(angles_deg, dtype=float)))
    cdf, pdf = piecewise_linear_cdf(u_nodes, values)
    mean, var = grid_moments(float(u_nodes[0]), float(u_nodes[-1]), pdf)
    return Reference(cdf, mean, var, float(u_nodes[0]), float(u_nodes[-1]))


def ref_uniform(low: float, high: float) -> Reference:
    """Uniform on [low, high] (the spectrum 'uniform' and the pulse 'tophat')."""
    return Reference(
        lambda v: np.clip((np.asarray(v, dtype=float) - low) / (high - low), 0.0, 1.0),
        0.5 * (low + high),
        (high - low) ** 2 / 12.0,
        low,
        high,
    )


def ref_inverse_square(low: float, high: float) -> Reference:
    """Density proportional to 1 / lambda^2 on [low, high]."""
    a, b = 1.0 / low, 1.0 / high

    def cdf(lam):
        v = np.clip(np.asarray(lam, dtype=float), low, high)
        return np.clip((a - 1.0 / v) / (a - b), 0.0, 1.0)

    mean = math.log(high / low) / (a - b)
    second = (high - low) / (a - b)
    return Reference(cdf, mean, second - mean * mean, low, high)


def ref_table_spectrum(wavelengths, density) -> Reference:
    cdf, pdf = piecewise_linear_cdf(wavelengths, density)
    mean, var = grid_moments(float(wavelengths[0]), float(wavelengths[-1]), pdf)
    return Reference(cdf, mean, var, float(wavelengths[0]), float(wavelengths[-1]))


def ref_gauss_pulse(mean: float, sigma: float) -> Reference:
    return Reference(
        lambda t: ndtr((np.asarray(t, dtype=float) - mean) / sigma),
        mean,
        sigma * sigma,
        mean - 10 * sigma,
        mean + 10 * sigma,
    )


def ref_exponential_pulse(mean: float, tau: float) -> Reference:
    """Leading edge at mean - tau; density exp(-(t - edge) / tau) / tau above it."""
    edge = mean - tau

    def cdf(t):
        t = np.asarray(t, dtype=float)
        return np.where(
            t >= edge, 1.0 - np.exp(-np.clip(t - edge, 0.0, None) / tau), 0.0
        )

    return Reference(cdf, mean, tau * tau, edge, edge + 40.0 * tau)


def complex_pulse_cdf(flat: float, rms: float, decay: float, delay: float) -> Callable:
    """CDF of delay + flat (U - 1/2) + rms N + decay (E - 1), U uniform, N normal, E ~ Exp(1).

    The flat-plus-Gaussian part has the closed CDF F(z) = (rms / flat) [G((z + a)/rms) -
    G((z - a)/rms)] with a = flat / 2 and G(w) = w Phi(w) + phi(w) (the integral of Phi). The
    exponential tail is integrated over its density (Gauss-Legendre on [-decay, 60 decay]),
    tabulated on a grid and interpolated.
    """
    half = 0.5 * flat

    def g(w):
        return w * ndtr(w) + np.exp(-0.5 * w * w) / math.sqrt(2.0 * math.pi)

    def flat_gauss_cdf(z):
        return (rms / flat) * (g((z + half) / rms) - g((z - half) / rms))

    nodes, weights = np.polynomial.legendre.leggauss(1500)
    lo, hi = -decay, 60.0 * decay
    y = 0.5 * (hi - lo) * nodes + 0.5 * (hi + lo)
    wy = 0.5 * (hi - lo) * weights
    density_y = np.exp(-(y + decay) / decay) / decay

    grid = np.linspace(delay - 30.0, delay + 60.0, 40_001)
    table = np.empty_like(grid)
    for i in range(0, grid.size, 2000):
        chunk = grid[i : i + 2000]
        table[i : i + 2000] = (
            flat_gauss_cdf(chunk[:, None] - delay - y[None, :])
            * density_y[None, :]
            * wy[None, :]
        ).sum(axis=1)

    def cdf(t):
        return np.interp(np.asarray(t, dtype=float), grid, table, left=0.0, right=1.0)

    return cdf


def ref_complex_pulse(flat: float, rms: float, decay: float, delay: float) -> Reference:
    var = flat * flat / 12.0 + rms * rms + decay * decay
    return Reference(
        complex_pulse_cdf(flat, rms, decay, delay),
        delay,
        var,
        delay - 30.0,
        delay + 60.0,
    )


# ------------------------------------------------------------------ profiles


@dataclass(frozen=True)
class Profile:
    name: str
    group: str  # "angular", "spectrum" or "pulse"
    emission: dict  # the emission (point, no target), without seed
    reference: (
        Reference | None
    )  # None for degenerate profiles (checked for exact values)
    fixed_value: float | None = (
        None  # degenerate: the exact value of the checked column
    )


ANGULAR_TABLE_DEG = [0.0, 20.0, 40.0, 70.0, 110.0, 150.0]
ANGULAR_TABLE_VALUE = [0.3, 1.0, 0.2, 0.8, 0.5, 0.1]
SPECTRUM_TABLE_NM = [300.0, 380.0, 420.0, 500.0, 600.0]
SPECTRUM_TABLE_VALUE = [0.2, 1.0, 0.5, 1.5, 0.0]
PULSE_TABLE_NS = [0.0, 1.0, 3.0, 6.0, 10.0]
PULSE_TABLE_RATE = [0.0, 2.0, 1.0, 0.5, 0.0]
FIXED_SPECTRUM = {"kind": "fixed", "wavelength_nm": DECIMAL_NM}
FIXED_PULSE = {"kind": "fixed", "time_ns": 0.0}


def point(angular=None, spectrum=None, pulse=None) -> dict:
    return {
        "kind": "point",
        "photons": float(PHOTONS),
        "position_m": list(SOURCE_POSITION_M),
        "angular": angular or {"kind": "isotropic"},
        "spectrum": spectrum or FIXED_SPECTRUM,
        "pulse": pulse or FIXED_PULSE,
        "target": "none",
    }


def profiles() -> list[Profile]:
    return [
        Profile(
            "ang_isotropic", "angular", point({"kind": "isotropic"}), ref_isotropic()
        ),
        Profile(
            "ang_lambertian", "angular", point({"kind": "lambertian"}), ref_lambertian()
        ),
        Profile(
            "ang_filled_cone",
            "angular",
            point({"kind": "filled_cone", "angle_deg": COS_CONE_DEG}),
            ref_filled_cone(COS_CONE_DEG),
        ),
        Profile(
            "ang_gauss",
            "angular",
            point({"kind": "gauss", "sigma_deg": 5.0}),
            ref_gauss(5.0),
        ),
        Profile("ang_rayleigh", "angular", point({"kind": "rayleigh"}), ref_rayleigh()),
        Profile(
            "ang_table",
            "angular",
            point(
                {
                    "kind": "table",
                    "angle_deg": ANGULAR_TABLE_DEG,
                    "value": ANGULAR_TABLE_VALUE,
                }
            ),
            ref_angular_table(ANGULAR_TABLE_DEG, ANGULAR_TABLE_VALUE),
        ),
        Profile(
            "ang_cone",
            "angular",
            point({"kind": "cone", "angle_deg": COS_CONE_DEG}),
            None,
            fixed_value=math.cos(math.radians(COS_CONE_DEG)),
        ),
        Profile(
            "spec_uniform",
            "spectrum",
            point(spectrum={"kind": "uniform", "low_nm": 300.0, "high_nm": 600.0}),
            ref_uniform(300.0, 600.0),
        ),
        Profile(
            "spec_inverse_square",
            "spectrum",
            point(
                spectrum={"kind": "inverse_square", "low_nm": 300.0, "high_nm": 600.0}
            ),
            ref_inverse_square(300.0, 600.0),
        ),
        Profile(
            "spec_table",
            "spectrum",
            point(
                spectrum={
                    "kind": "table",
                    "wavelength_nm": SPECTRUM_TABLE_NM,
                    "density": SPECTRUM_TABLE_VALUE,
                }
            ),
            ref_table_spectrum(SPECTRUM_TABLE_NM, SPECTRUM_TABLE_VALUE),
        ),
        Profile(
            "spec_fixed",
            "spectrum",
            point(spectrum={"kind": "fixed", "wavelength_nm": 450.0}),
            None,
            fixed_value=450.0,
        ),
        Profile(
            "pulse_tophat",
            "pulse",
            point(pulse={"kind": "tophat", "width_ns": 5.0, "delay_ns": 10.0}),
            ref_uniform(7.5, 12.5),
        ),
        Profile(
            "pulse_gauss",
            "pulse",
            point(pulse={"kind": "gauss", "sigma_ns": 1.5, "mean_ns": 3.0}),
            ref_gauss_pulse(3.0, 1.5),
        ),
        Profile(
            "pulse_exponential",
            "pulse",
            point(pulse={"kind": "exponential", "decay_ns": 1.5, "mean_ns": 2.0}),
            ref_exponential_pulse(2.0, 1.5),
        ),
        Profile(
            "pulse_complex",
            "pulse",
            point(
                pulse={
                    "kind": "complex",
                    "flat_ns": 4.0,
                    "gauss_rms_ns": 1.0,
                    "decay_ns": 2.0,
                    "delay_ns": 1.0,
                }
            ),
            ref_complex_pulse(4.0, 1.0, 2.0, 1.0),
        ),
        Profile(
            "pulse_table",
            "pulse",
            point(
                pulse={
                    "kind": "table",
                    "time_ns": PULSE_TABLE_NS,
                    "rate": PULSE_TABLE_RATE,
                }
            ),
            ref_table_spectrum(PULSE_TABLE_NS, PULSE_TABLE_RATE),
        ),
        Profile(
            "pulse_fixed",
            "pulse",
            point(pulse={"kind": "fixed", "time_ns": 2.5}),
            None,
            fixed_value=2.5,
        ),
    ]


# ------------------------------------------------------------------ statistics


def ks_p_value(sample: np.ndarray, cdf: Callable) -> tuple[float, float]:
    result = stats.kstest(sample, cdf)
    return float(result.statistic), float(result.pvalue)


def moment_pulls(sample: np.ndarray, reference: Reference) -> dict:
    n = sample.size
    mean = float(np.mean(sample))
    centred = sample - mean
    var = float(np.mean(centred**2))
    m4 = float(np.mean(centred**4))
    se_mean = math.sqrt(reference.var / n)
    se_var = math.sqrt(max(m4 - var * var, 0.0) / n)
    pull_mean = (mean - reference.mean) / se_mean
    pull_var = (var - reference.var) / se_var
    return {
        "sample_mean": mean,
        "reference_mean": reference.mean,
        "sample_var": var,
        "reference_var": reference.var,
        "mean_pull": pull_mean,
        "var_pull": pull_var,
        "moment_pull_max": max(abs(pull_mean), abs(pull_var)),
    }


def azimuth_sample(cx: np.ndarray, cy: np.ndarray) -> np.ndarray:
    """Azimuth about the +z axis in units of a full turn, in [0, 1)."""
    return np.mod(np.arctan2(cy.astype(float), cx.astype(float)) / (2.0 * math.pi), 1.0)


def sample_profile(profile: Profile, seed: int) -> dict:
    """Emit the profile at the sky level (directions unchanged)."""
    from obdeect.sources import source

    emission = {**profile.emission}
    photons = source.emit(emission, seed=seed)
    count = source.photon_count(emission, seed=seed)
    if count != PHOTONS or photons["w"].size != PHOTONS:
        raise RuntimeError(
            f"{profile.name}: expected {PHOTONS} photons, got {photons['w'].size}"
        )
    if not np.all(photons["w"] == 1.0):
        raise RuntimeError(
            f"{profile.name}: weights are not all 1 (the KS test is unweighted)"
        )
    return photons


def quantity_of(profile: Profile, photons: dict) -> np.ndarray:
    if profile.group == "angular":
        return 1.0 - photons["cz"].astype(float)
    if profile.group == "spectrum":
        return photons["wavelength"].astype(float)
    return photons["t"].astype(float)


def profile_part(index: int, profile: Profile) -> dict:
    seed = SEED + index
    photons = sample_profile(profile, seed)
    entry = {
        "name": profile.name,
        "group": profile.group,
        "emission": profile.emission,
        "seed": seed,
        "photons": PHOTONS,
    }
    if profile.reference is None:
        if profile.group == "angular":
            deviation = np.abs(photons["cz"].astype(float) - profile.fixed_value).max()
        elif profile.group == "spectrum":
            deviation = np.abs(
                photons["wavelength"].astype(float) - profile.fixed_value
            ).max()
        else:
            deviation = np.abs(photons["t"].astype(float) - profile.fixed_value).max()
        entry["degenerate"] = True
        entry["fixed_value"] = profile.fixed_value
        entry["max_dev"] = float(deviation)
    else:
        sample = quantity_of(profile, photons)
        statistic, p_value = ks_p_value(sample, profile.reference.cdf)
        entry.update(
            {
                "degenerate": False,
                "ks_statistic": statistic,
                "ks_p": p_value,
                **moment_pulls(sample, profile.reference),
            }
        )
    if profile.group == "angular":
        az_stat, az_p = ks_p_value(
            azimuth_sample(photons["cx"], photons["cy"]), _uniform_cdf
        )
        entry["azimuth_statistic"] = az_stat
        entry["azimuth_p"] = az_p
    # Other columns of the profile are fixed by the defaults (wavelength 400 nm, t = 0 ns).
    other = {}
    if profile.group != "spectrum":
        other["wavelength_dev_nm"] = float(
            np.abs(photons["wavelength"] - DECIMAL_NM).max()
        )
    if profile.group != "pulse":
        other["time_dev_ns"] = float(np.abs(photons["t"]).max())
    entry["fixed_other_columns"] = other
    return entry


def _uniform_cdf(v):
    return np.clip(np.asarray(v, dtype=float), 0.0, 1.0)


def control_part(seed: int = SEED) -> dict:
    """Isotropic sample against the Lambertian CDF: must be rejected."""
    photons = sample_profile(profiles()[0], seed)
    statistic, p_value = ks_p_value(
        1.0 - photons["cz"].astype(float), ref_lambertian().cdf
    )
    return {
        "statistic": statistic,
        "p_value": p_value,
        "rejected": p_value < P_THRESHOLD,
    }


# ------------------------------------------------------------------ part 2: emission levels


def level_emission() -> dict:
    return {
        "kind": "point",
        "photons": float(PHOTONS),
        "position_m": list(SOURCE_POSITION_M),
        "angular": {"kind": "isotropic"},
        "spectrum": {"kind": "uniform", "low_nm": 300.0, "high_nm": 600.0},
        "pulse": {"kind": "tophat", "width_ns": 5.0, "delay_ns": 10.0},
        "target": "none",
    }


def level_run(frame: str) -> dict:
    return {
        "run": {"seed": SEED, "observation_level_m": OBSERVATION_M, "frame": frame},
        "telescopes": [
            {"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": TELESCOPE_RADIUS_M}
        ],
        "events": [{"repeat": 1, "emissions": [level_emission()]}],
    }


def read_stream(path: Path) -> dict:
    """Photon columns and frame code of a stream (one telescope, one event)."""
    from obdeect.stream import Reader

    columns: dict[str, list] = {}
    frames = set()
    for record in Reader(str(path)):
        if record["type"] != "photons":
            continue
        frames.add(int(record["frame"]))
        for name, values in record["columns"].items():
            columns.setdefault(name, []).append(np.asarray(values))
    if len(frames) != 1:
        raise RuntimeError(f"{path}: expected one frame code, got {sorted(frames)}")
    data = {name: np.concatenate(parts) for name, parts in columns.items()}
    data["frame"] = frames.pop()
    return data


def ulp32(values: np.ndarray) -> np.ndarray:
    """Float32 spacing at each value (in double)."""
    return np.spacing(np.abs(values).astype(np.float32)).astype(float)


def ground_expected(sky: dict) -> dict:
    """Sky rows that hit the telescope sphere: entry point, path and time (smallest root)."""
    p = np.array(SOURCE_POSITION_M, dtype=float)
    d = np.stack([sky["cx"], sky["cy"], sky["cz"]]).astype(float)
    # The entry is solved with the normalised direction (|d| = 1 to float32 rounding, as the sources
    # step does); the stored cx, cy, cz are float32 and need not be exactly unit length.
    norm = np.linalg.norm(d, axis=0)
    d = d / np.where(norm > 0.0, norm, 1.0)
    b = d.T @ p
    c = float(p @ p) - TELESCOPE_RADIUS_M**2
    disc = b * b - c
    hit = (disc >= 0.0) & (-b > 0.0)
    s = (-b - np.sqrt(np.where(disc >= 0.0, disc, 0.0)))[hit]
    entry = p[None, :] + s[:, None] * d[:, hit].T
    return {
        "hit": hit,
        "s": s,
        "entry": entry,
        "t": sky["t"].astype(float)[hit] + s / C_M_PER_NS,
    }


def level_part(root: Path) -> tuple[dict, dict]:
    """Sky, T, C and G runs of the same source; returns (metrics, details)."""
    from obdeect.sources import config, run

    work = root / "levels"
    work.mkdir(exist_ok=True)
    data = {}
    summaries = {}
    try:
        for frame in ("sky", "telescope", "camera", "ground"):
            path = work / f"{frame}.pst"
            description = config.normalize(level_run(frame))
            summaries[frame] = run.write(description, path)
            data[frame] = read_stream(path)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    sky, telescope, camera, ground = (
        data[k] for k in ("sky", "telescope", "camera", "ground")
    )
    expected_frames = {"sky": 0, "telescope": 2, "camera": 3, "ground": 1}
    frame_codes_mismatch = sum(data[k]["frame"] != expected_frames[k] for k in data)
    rows = {k: data[k]["x"].size for k in ("sky", "telescope", "camera")}
    row_delta = max(abs(n - PHOTONS) for n in rows.values())

    def max_dev(a: dict, b: dict, names) -> float:
        return max(
            float(np.max(np.abs(a[n].astype(float) - b[n].astype(float))))
            for n in names
        )

    position = ("x", "y", "z")
    direction = ("cx", "cy", "cz")
    frame_position_dev = max(
        max_dev(telescope, camera, position), max_dev(sky, telescope, position)
    )
    frame_direction_dev = max(
        max_dev(telescope, camera, direction), max_dev(sky, telescope, direction)
    )
    frame_tw_dev = max(
        max_dev(telescope, camera, ("t", "w", "wavelength")),
        max_dev(sky, telescope, ("t", "w", "wavelength")),
    )
    sky_zem_dev = float(
        np.max(np.abs(sky["zem"].astype(float) - (OBSERVATION_M + sky["z"])))
    )

    expect = ground_expected(sky)
    hit = expect["hit"]
    ground_count_delta = int(ground["x"].size) - int(hit.sum())
    if ground_count_delta == 0:
        ground_direction_dev = max(
            float(np.max(np.abs(ground[n].astype(float) - sky[n][hit].astype(float))))
            for n in direction
        )
        g_pos = np.stack([ground["x"], ground["y"], ground["z"]]).astype(float).T
        ground_entry_dev = float(np.max(np.abs(g_pos - expect["entry"])))
        ground_radius_dev = float(
            np.max(np.abs(np.linalg.norm(g_pos, axis=1) - TELESCOPE_RADIUS_M))
        )
        ground_time_ulp = float(
            np.max(np.abs(ground["t"].astype(float) - expect["t"]) / ulp32(expect["t"]))
        )
        ground_zem_dev = float(
            np.max(
                np.abs(
                    ground["zem"].astype(float) - (OBSERVATION_M + SOURCE_POSITION_M[2])
                )
            )
        )
        ground_wavelength_dev = float(
            np.max(np.abs(ground["wavelength"].astype(float) - sky["wavelength"][hit]))
        )
    else:
        ground_direction_dev = ground_entry_dev = ground_radius_dev = float("inf")
        ground_time_ulp = ground_zem_dev = ground_wavelength_dev = float("inf")

    metrics = {
        "level_row_count_delta": row_delta,
        "frame_codes_mismatch": frame_codes_mismatch,
        "frame_position_max_dev_m": frame_position_dev,
        "frame_direction_max_dev": frame_direction_dev,
        "frame_time_wavelength_max_dev": frame_tw_dev,
        "sky_zem_max_dev_m": sky_zem_dev,
        "ground_count_delta": ground_count_delta,
        "ground_direction_max_dev": ground_direction_dev,
        "ground_entry_max_dev_m": ground_entry_dev,
        "ground_radius_max_dev_m": ground_radius_dev,
        "ground_time_max_ulp": ground_time_ulp,
        "ground_zem_max_dev_m": ground_zem_dev,
        "ground_wavelength_max_dev": ground_wavelength_dev,
    }
    details = {
        "rows": {k: int(v) for k, v in rows.items()}
        | {"ground": int(ground["x"].size)},
        "ground_expected_rows": int(hit.sum()),
        "ground_path_mean_m": float(np.mean(ground["path"]))
        if "path" in ground
        else None,
        "run_summaries": summaries,
        "note": (
            "sources have no T-to-C transform: the frame is a label, so T and C rows equal the "
            "sky rows; the telescope-to-camera conversion lives in obdeect stream to-simtel"
        ),
    }
    return metrics, details


# ------------------------------------------------------------------ run


def version() -> str:
    return common.package_version("obdeect")


def run_all(root: Path) -> tuple[dict, dict]:
    entries = [profile_part(i, p) for i, p in enumerate(profiles())]
    tested = [e for e in entries if not e["degenerate"]]
    azimuth = [e for e in entries if "azimuth_p" in e]
    ks_fail = [e["name"] for e in tested if e["ks_p"] < P_THRESHOLD] + [
        e["name"] for e in azimuth if e["azimuth_p"] < P_THRESHOLD
    ]
    moment_max = max(e["moment_pull_max"] for e in tested)
    degenerate = [e for e in entries if e["degenerate"]]
    fixed_dev = max(
        [e["max_dev"] for e in degenerate]
        + [max(e["fixed_other_columns"].values(), default=0.0) for e in entries]
    )
    control = control_part()
    level_metrics, level_details = level_part(root)
    metrics = {
        "ks_fail_count": len(ks_fail),
        "moment_pull_max": moment_max,
        "fixed_value_max_dev": fixed_dev,
        "control_fail_detected": int(control["rejected"]),
        **level_metrics,
    }
    details = {
        "version": version(),
        "photons": PHOTONS,
        "seed_base": SEED,
        "ks_p_threshold": P_THRESHOLD,
        "moment_pull_limit": MOMENT_PULL_LIMIT,
        "ks_tests": len(tested) + len(azimuth),
        "ks_failures": ks_fail,
        "profiles": entries,
        "control": control,
        "levels": level_details,
    }
    return metrics, details


def main() -> None:
    common.results_main(run_all, __doc__.splitlines()[0])


if __name__ == "__main__":
    main()
