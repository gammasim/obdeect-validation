"""Photon targets of obdeect point and line sources: light at the target, row waste, finiteness.

Targets checked: a disc (uniform-area sampling with weight I/I_max (h/d)^3), a sphere (isotropic
cone, weight 1), the automatic sphere of a run (threshold 1.02 R), and an isotropic line with a
sphere target that misses it, grazes it and crosses it. Every reference is written here from the
profile definitions (closed forms, scipy quadrature for the Gaussian), not taken from obdeect.

Physics of the reference. A point emits ``photons`` physical photons with directions distributed
as the profile ``I(theta)`` (normalised over the sphere). The light that reaches a target is
``N F``, with F the fraction of the emitted light inside the target's solid angle:

- disc on the axis of the point (half-angle alpha = atan(R/h) about the axis, h = distance of the
  point to the disc plane): isotropic F = (1 - cos alpha) / 2; lambertian (I = cos theta)
  F = sin^2 alpha = R^2 / (R^2 + h^2); filled cone (half-angle theta_c, I = 1 inside) F = (1 -
  cos min(alpha, theta_c)) / (1 - cos theta_c); Gaussian (I = exp(-theta^2 / 2 sigma^2) up to
  theta_max = min(5 sigma, pi), as the obdeect table) F from scipy quadrature.
- disc with the axis tilted by ``TILT_DEG`` from the direction to the disc centre, with
  ``TILT_DEG > theta_c + alpha`` for every disc: the cone does not reach the disc, F = 0 exactly.
- sphere (radius R, centre C, point at distance d): isotropic F = (1 - sqrt(1 - (R/d)^2)) / 2.

Checks and gates (each gate: tolerance = max(requirement, k x standard error), k stated):

1. ``weight_sum_pull_max``: |sum w - N F| / sigma over all disc and sphere configurations, with
   sigma^2 = (sum w^2 - (sum w)^2 / M) + (mean w)^2 / 4 (the sample variance of the M rows, plus the
   rounding of the row count, at most one row). Gate 5 (5 sigma): the weight sum is an unbiased
   estimator of the light at the target; the sample variance is the sampler's own noise, so the
   test is statistical (k = 5 as in the other sources validations).
2. ``zero_weight_rows``: rows with w == 0 written, summed over all configurations. Exact gate 0.
   Justification: a row with zero weight carries no light; it costs stream space, time and
   memory in every downstream stage (AGENTS: explicit loss accounting, no silent rows). A sampler
   that draws target points only where the light is writes none.
3. ``rows_per_hit_photon_max``: max over configurations with F > 0 of rows / (N F), the rows
   written per photon that reaches the target. Gate 1.5 (requirement: a sampler that emits only
   the photons that hit has ratio 1 up to the Poisson rounding of the count; 50 % margin for
   the rounding and for the solid-angle approximation of the count, documented as the
   threshold of acceptable waste). Configurations with F = 0 are not in the ratio (infinite); they
   are counted in ``zero_weight_rows``.
4. ``nonfinite_weights``: NaN or inf values in any column of any row, over all configurations
   including the line with a sphere target. Exact gate 0 (a NaN weight propagates into every sum).
5. ``sphere_threshold_jump``: number of sides (0 to 2) of the automatic-sphere threshold pair
   whose weight sum W disagrees with its own analytic sphere fraction, |W - N F(d)| <= 5 sigma,
   for d = 1.019 R and d = 1.021 R (isotropic point, ``target = "auto"``, one telescope of radius R,
   rows of weight 1; sigma = 0.5, the rounding bound of the row count, as in ``count_sigma``). The
   light must be continuous across the old 1.02 R threshold: the sphere target writes the light of
   the point at its own distance on either side of 1.02 R. Exact gate 0.

   Why the old definition (|W(1.019 R) - W(1.021 R)| <= 5 sigma with sigma^2 = 1/4 + 1/4) was wrong:
   the two distances have different true light (F = 0.4035 and 0.3988, a difference of about 940 of
   the 200000 photons, which is 1300 sigma), so a correct sampler fails it. It was a design error in
   the check (it compared two different physical quantities as if equal), not a tolerance to relax.
   The difference of the two sums is printed as a diagnostic only.
6. ``expected_count_pull_max``: |rows - expected_photons| / sigma over all emissions, with
   expected_photons the API count of candidate rows (before rounding) and sigma = 0.5, the fixed
   rounding bound of the count (``count_sigma``: the rounding moves the mean to a count by at most
   one row). The exact standard deviation sqrt(f (1 - f)) is not used (for f near 1 its pull is
   unbounded). Gate 5. When f = 0 the count is exact and the pull is 0 or infinite.

Extra (diagnostic, gate 0): ``sphere_nonisotropic_accepted``: point emissions with a sphere
target and a non-isotropic profile that obdeect accepts. The validation rejects them by design
(``a sphere target needs an isotropic point``), so the count is 0.

Configurations are written with ``PHOTONS`` photons each (the mean physical count), seeds fixed
per configuration. Angular profiles: ``cone`` is the fixed-angle profile (a delta on the cone);
the disc targets need a density, so the filled cone ``filled_cone`` is used for the cone cases.
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import numpy as np
from scipy import integrate

from obdeect_validation.adapters import common

SEED = 20261123
PHOTONS = 200_000  # mean physical photons per configuration (upper bound)
ROW_BUDGET = (
    3_000_000  # expected rows per configuration: photons are reduced above this
)
COLUMNS = ("x", "y", "z", "cx", "cy", "cz", "t", "w", "wavelength")
WAVELENGTH = {"kind": "fixed", "wavelength_nm": 400.0}

DISC_HEIGHT_M = 100.0  # distance of the point to the disc plane
DISC_RATIOS = (0.05, 0.2, 0.5, 1.0)  # R / h
TILT_DEG = 60.0  # axis tilt from the direction to the disc centre (narrow cones miss)
CONE_DEG = {"cone_11deg": 11.0, "cone_0.5deg": 0.5}
GAUSS_SIGMA_DEG = 5.0
DISC_PROFILES = ("isotropic", "lambertian", "cone_11deg", "cone_0.5deg", "gauss_5deg")
TILTED_PROFILES = ("cone_11deg", "cone_0.5deg")

SPHERE_RADIUS_M = 10.0
SPHERE_DISTANCES = (1.01, 1.03, 2.0, 10.0)  # d / R
THRESHOLD_PAIR = (1.019, 1.021)  # d / R around the 1.02 R automatic-target threshold
THRESHOLD_PHOTONS = 200_000

LINE_START_LENGTH = (
    ("line-miss", [-100.0, 30.0, 0.0], 200.0),  # y = 30 m > R: misses the sphere
    ("line-graze", [-100.0, 9.5, 0.0], 200.0),  # chord of 6.2 m: partly hits
    ("line-cross", [-100.0, 0.0, 0.0], 200.0),  # through the centre
)

GATE_PULL = 5.0  # k of the standard error, same as the other sources validations
GATE_ROWS_PER_HIT = 1.5
NONISOTROPIC_SPHERE = (
    {"kind": "lambertian"},
    {"kind": "filled_cone", "angle_deg": 11.0},
    {"kind": "filled_cone", "angle_deg": 0.5},
    {"kind": "gauss", "sigma_deg": GAUSS_SIGMA_DEG},
)


def angular_of(profile: str) -> dict:
    """Angular profile of a configuration name."""
    if profile == "isotropic":
        return {"kind": "isotropic"}
    if profile == "lambertian":
        return {"kind": "lambertian"}
    if profile in CONE_DEG:
        return {"kind": "filled_cone", "angle_deg": CONE_DEG[profile]}
    if profile == "gauss_5deg":
        return {"kind": "gauss", "sigma_deg": GAUSS_SIGMA_DEG}
    raise ValueError(f"unknown profile {profile}")


# ------------------------------------------------------------------ reference fractions


def gauss_cap_fraction(alpha_rad: float, sigma_deg: float = GAUSS_SIGMA_DEG) -> float:
    """Fraction of the Gaussian profile inside the cone of half-angle alpha (on axis).

    Intensity exp(-theta^2 / 2 sigma^2) for theta < theta_max = min(5 sigma, pi), zero beyond
    (the extent of the obdeect angular table); the normalisation is the integral over the table.
    """
    sigma = math.radians(sigma_deg)
    theta_max = min(5.0 * sigma, math.pi)

    def weight(theta: float) -> float:
        return math.exp(-0.5 * (theta / sigma) ** 2) * math.sin(theta)

    total = integrate.quad(weight, 0.0, theta_max, epsabs=0.0, epsrel=1e-12)[0]
    part = integrate.quad(
        weight, 0.0, min(alpha_rad, theta_max), epsabs=0.0, epsrel=1e-12
    )[0]
    return part / total


def cap_fraction(profile: str, alpha_rad: float) -> float:
    """Fraction of the emitted light of an on-axis profile inside the cone of half-angle alpha."""
    if profile == "isotropic":
        return 0.5 * (1.0 - math.cos(alpha_rad))
    if (
        profile == "lambertian"
    ):  # I = cos(theta) on the hemisphere: integral over the cap = sin^2
        return math.sin(min(alpha_rad, 0.5 * math.pi)) ** 2
    if profile in CONE_DEG:
        theta_c = math.radians(CONE_DEG[profile])
        return (1.0 - math.cos(min(alpha_rad, theta_c))) / (1.0 - math.cos(theta_c))
    if profile == "gauss_5deg":
        return gauss_cap_fraction(alpha_rad)
    raise ValueError(f"unknown profile {profile}")


def sphere_fraction(distance_over_radius: float) -> float:
    """Isotropic fraction of the emitted light that hits a sphere at distance d (d >= R)."""
    return 0.5 * (1.0 - math.sqrt(1.0 - distance_over_radius**-2))


# ------------------------------------------------------------------ configurations


def disc_configs() -> list[dict]:
    """Point sources aimed at a disc (centre origin, normal +z), at several R/h."""
    out = []
    for ratio in DISC_RATIOS:
        radius = ratio * DISC_HEIGHT_M
        alpha = math.atan(ratio)
        target = {
            "kind": "disc",
            "centre_m": [0.0, 0.0, 0.0],
            "normal": [0.0, 0.0, 1.0],
            "radius_m": radius,
        }
        for profile in DISC_PROFILES:
            out.append(
                _disc(
                    profile, ratio, alpha, target, axis=[0.0, 0.0, -1.0], tilted=False
                )
            )
        for profile in TILTED_PROFILES:
            theta_c = CONE_DEG[profile]
            if math.degrees(alpha) + theta_c >= TILT_DEG:
                raise ValueError("tilt too small: the cone may reach the disc")
            tilt = math.radians(TILT_DEG)
            axis = [math.sin(tilt), 0.0, -math.cos(tilt)]
            out.append(_disc(profile, ratio, alpha, target, axis=axis, tilted=True))
    return out


def _disc(profile, ratio, alpha, target, axis, tilted) -> dict:
    fraction = 0.0 if tilted else cap_fraction(profile, alpha)
    name = f"disc-{profile}-R_h{ratio:g}" + ("-tilted" if tilted else "")
    return {
        "id": name,
        "kind": "disc",
        "profile": profile,
        "fraction": fraction,
        "emission": {
            "kind": "point",
            "photons": PHOTONS,
            "position_m": [0.0, 0.0, DISC_HEIGHT_M],
            "axis": axis,
            "angular": angular_of(profile),
            "spectrum": WAVELENGTH,
            "target": target,
        },
    }


def sphere_configs() -> list[dict]:
    """Isotropic point sources aimed at a sphere (centre origin) at several distances."""
    out = []
    for factor in SPHERE_DISTANCES:
        out.append(
            {
                "id": f"sphere-isotropic-d{factor:g}R",
                "kind": "sphere",
                "profile": "isotropic",
                "fraction": sphere_fraction(factor),
                "emission": {
                    "kind": "point",
                    "photons": PHOTONS,
                    "position_m": [0.0, 0.0, factor * SPHERE_RADIUS_M],
                    "angular": {"kind": "isotropic"},
                    "spectrum": WAVELENGTH,
                    "target": {
                        "kind": "sphere",
                        "centre_m": [0.0, 0.0, 0.0],
                        "radius_m": SPHERE_RADIUS_M,
                    },
                },
            }
        )
    return out


def line_configs() -> list[dict]:
    """Isotropic line with a sphere target: misses, grazes, crosses (no light fraction is gated)."""
    out = []
    for name, start, length in LINE_START_LENGTH:
        out.append(
            {
                "id": name,
                "kind": "line",
                "profile": "isotropic",
                "fraction": None,
                "emission": {
                    "kind": "line",
                    "photons": PHOTONS,
                    "start_m": start,
                    "direction": [1.0, 0.0, 0.0],
                    "length_m": length,
                    "angular": {"kind": "isotropic"},
                    "spectrum": WAVELENGTH,
                    "target": {
                        "kind": "sphere",
                        "centre_m": [0.0, 0.0, 0.0],
                        "radius_m": SPHERE_RADIUS_M,
                    },
                },
            }
        )
    return out


# ------------------------------------------------------------------ statistics


def pull(value: float, expected: float, sigma: float) -> float:
    """|value - expected| / sigma; zero sigma: 0 if equal, infinity otherwise."""
    if sigma > 0.0:
        return abs(value - expected) / sigma
    return 0.0 if value == expected else math.inf


def count_sigma(mean: float) -> float:
    """Sigma of the count pull: the rounding of a mean to a count moves it by at most one row.

    The exact standard deviation sqrt(f (1 - f)) (f the fractional part) is not used: for f near 1
    the floor outcome has a pull sqrt(f / (1 - f)), which is unbounded. Sigma is the bound 0.5 of
    the rounding (|rows - mean| <= 1), so the pull is at most 2 for a correct count.
    """
    return 0.5


def tally(emission: dict, seed: int) -> dict:
    """Rows, weight sums, zero weights and non-finite values of one emission (streamed)."""
    from obdeect.sources import source

    rows = zeros = nonfinite = 0
    w_sum = w2_sum = 0.0
    for chunk in source.iter_chunks(emission, seed=seed):
        w = np.asarray(chunk["w"], dtype=np.float64)
        rows += int(w.size)
        zeros += int(np.count_nonzero(w == 0.0))
        w_sum += float(np.sum(w))
        w2_sum += float(np.sum(w * w))
        for name in COLUMNS:
            nonfinite += int(np.count_nonzero(~np.isfinite(chunk[name])))
    return {
        "rows": rows,
        "zeros": zeros,
        "nonfinite": nonfinite,
        "w_sum": w_sum,
        "w2_sum": w2_sum,
    }


def weight_sigma(t: dict) -> float:
    """sigma of the weight sum: sample variance of M rows plus the rounding of the count."""
    rows = t["rows"]
    if rows == 0:
        return 0.0
    mean_w = t["w_sum"] / rows
    spread = max(0.0, t["w2_sum"] - t["w_sum"] ** 2 / rows)
    return math.sqrt(spread + 0.25 * mean_w * mean_w)


def evaluate_config(config: dict, seed: int) -> dict:
    """Tallies of one configuration and its pulls, rows per hit and count pull."""
    from obdeect.sources import source

    emission = config["emission"]
    rows_per_photon = source.expected_photons({**emission, "photons": 1.0})
    if rows_per_photon > 0.0:  # cap the expected rows of the configuration (run time)
        emission["photons"] = min(
            float(emission["photons"]), ROW_BUDGET / rows_per_photon
        )
    t = tally(emission, seed)
    photons = float(emission["photons"])
    expected_rows = source.expected_photons(emission)
    expected_w = source.expected_weight(emission)
    fraction = config["fraction"]
    out = {
        "id": config["id"],
        "kind": config["kind"],
        "profile": config["profile"],
        "photons": photons,
        "fraction": fraction,
        "rows": t["rows"],
        "expected_rows_api": expected_rows,
        "expected_weight_api": expected_w,
        "weight_sum": t["w_sum"],
        "zero_weight_rows": t["zeros"],
        "nonfinite": t["nonfinite"],
        "count_pull": pull(t["rows"], expected_rows, count_sigma(expected_rows)),
    }
    if fraction is not None:
        light = photons * fraction
        sigma = weight_sigma(t)
        out["expected_light"] = light
        out["weight_sigma"] = sigma
        out["weight_pull"] = pull(t["w_sum"], light, sigma)
        out["rows_per_hit"] = t["rows"] / max(light, 1.0) if fraction > 0.0 else None
    return out


def threshold_jump(seed: int) -> dict:
    """Weight sum of a run with target auto at each side of 1.02 R, against its own light.

    Each side is checked against N F(d) of that side (docs: check 5). The rows have weight 1, so
    the sigma of the sum is the rounding of its count, 0.5 (``count_sigma``). The difference of the
    two sums (the old check) is returned as a diagnostic only: the two sides have different light.
    """
    from obdeect.sources import config, run

    sums = {}
    with tempfile.TemporaryDirectory() as tmp:
        for index, factor in enumerate(THRESHOLD_PAIR):
            description = {
                "run": {"seed": seed + index, "observation_level_m": 2150.0},
                "telescopes": [{"radius_m": SPHERE_RADIUS_M}],
                "events": [
                    {
                        "emissions": [
                            {
                                "kind": "point",
                                "photons": THRESHOLD_PHOTONS,
                                "position_m": [0.0, 0.0, factor * SPHERE_RADIUS_M],
                                "angular": {"kind": "isotropic"},
                                "spectrum": WAVELENGTH,
                                "target": "auto",
                            }
                        ]
                    }
                ],
            }
            summary = run.write(
                config.normalize(description), str(Path(tmp) / f"threshold{index}.pst")
            )
            light = THRESHOLD_PHOTONS * sphere_fraction(factor)
            sigma = count_sigma(light)
            sums[factor] = {
                "weight": float(summary["weight"]),
                "rows": int(summary["rows"]),
                "sphere_light": light,
                "pull": pull(float(summary["weight"]), light, sigma),
            }
    lo, hi = THRESHOLD_PAIR
    violations = sum(int(sums[f]["pull"] > GATE_PULL) for f in THRESHOLD_PAIR)
    return {
        "factors": list(THRESHOLD_PAIR),
        "weight_below": sums[lo]["weight"],
        "weight_above": sums[hi]["weight"],
        "sphere_light_below": sums[lo]["sphere_light"],
        "sphere_light_above": sums[hi]["sphere_light"],
        "pull_below": sums[lo]["pull"],
        "pull_above": sums[hi]["pull"],
        "sigma": count_sigma(sums[lo]["sphere_light"]),
        "violations": violations,
        "violation": violations,
        "difference_diagnostic": sums[lo]["weight"] - sums[hi]["weight"],
    }


def nonisotropic_sphere_accepted() -> int:
    """Non-isotropic point profiles that obdeect accepts with a sphere target (should be 0)."""
    from obdeect.sources import source

    accepted = 0
    for angular in NONISOTROPIC_SPHERE:
        emission = {
            "kind": "point",
            "photons": 1.0,
            "position_m": [0.0, 0.0, 20.0],
            "angular": angular,
            "spectrum": WAVELENGTH,
            "target": {"kind": "sphere", "centre_m": [0.0, 0.0, 0.0], "radius_m": 10.0},
        }
        try:
            source.expected_photons(emission)
        except ValueError:
            continue
        accepted += 1
    return accepted


# ------------------------------------------------------------------ run


def run_all() -> tuple[dict, dict]:
    """All configurations and the threshold pair; returns (metrics, details)."""
    configs = disc_configs() + sphere_configs() + line_configs()
    rows = [evaluate_config(c, SEED + i) for i, c in enumerate(configs)]
    threshold = threshold_jump(SEED + 1000)
    accepted = nonisotropic_sphere_accepted()

    targets = [r for r in rows if r["kind"] in ("disc", "sphere")]
    positive = [r for r in targets if r["fraction"] and r["fraction"] > 0.0]
    metrics = {
        "weight_sum_pull_max": max(r["weight_pull"] for r in targets),
        "zero_weight_rows": sum(r["zero_weight_rows"] for r in rows),
        "rows_per_hit_photon_max": max(r["rows_per_hit"] for r in positive),
        "nonfinite_weights": sum(r["nonfinite"] for r in rows),
        "sphere_threshold_jump": threshold["violation"],
        "expected_count_pull_max": max(r["count_pull"] for r in rows),
        "sphere_nonisotropic_accepted": accepted,
    }
    details = {
        "seed": SEED,
        "photons_per_configuration": PHOTONS,
        "gates": {
            "weight_sum_pull_max": {"expected": 0, "atol": GATE_PULL},
            "zero_weight_rows": {"expected": 0, "atol": 0},
            "rows_per_hit_photon_max": {"expected": 1, "atol": GATE_ROWS_PER_HIT - 1},
            "nonfinite_weights": {"expected": 0, "atol": 0},
            "sphere_threshold_jump": {"expected": 0, "atol": 0},
            "expected_count_pull_max": {"expected": 0, "atol": GATE_PULL},
            "sphere_nonisotropic_accepted": {"expected": 0, "atol": 0},
        },
        "configurations": rows,
        "threshold": threshold,
        "zero_weight_rows_by_configuration": {
            r["id"]: r["zero_weight_rows"] for r in rows if r["zero_weight_rows"]
        },
    }
    return metrics, details


def main() -> None:
    common.results_main(lambda _root: run_all(), __doc__.splitlines()[0])


if __name__ == "__main__":
    main()
