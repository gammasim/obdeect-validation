"""Substitute obdeect atmospheric propagation for sim_telarray's atmosphere and compare cameras.

Reference: sim_telarray with its own atmospheric transmission table on CORSIKA IACT photons.
Candidate: obdeect propagation of the same photons (transmission table cut below 1000 nm, weighted
mode), then sim_telarray with unit transmission (zero optical depth table).

The obdeect wavelength band must equal the band sim_telarray samples from (CORSIKA
CHERENKOV WAVELENGTH BAND, read from corsika.log). A narrower band changes the photon population
and biases the candidate high.

The transmission expectation (expected_transmission) is sim_telarray's atmospheric_transmission2:
it includes the telescope altitude term tau(h_t) and the in-telescope focus term
tau(h_t + pathlen) - tau(h_t), as obdeect does with --site-altitude and --pathlen.
"""

import argparse
import math
import re
from pathlib import Path

from obdeect_validation.adapters import common

MAX_LAMBDA_NM = 1000.0  # sim_telarray skips bunches at or above MAX_LAMBDA; obdeect table ends at 999 nm
SITE_ALTITUDE_M = 2156.0
PATHLEN_M = 28.0
# z of the telescopes above the observation level in the flat CORSIKA cards (0 cm, see
# chain_showers.FLAT_TELESCOPES). obdeect and sim_telarray put telescopes at level + z.
TELESCOPE_Z_M = 0.0
SEED = 123
SAMPLE_FRACTION = 0.02
BAND_PATTERN = re.compile(
    r"CHERENKOV WAVELENGTH BAND FROM\s+([\d.]+)\s+TO\s+([\d.]+)\s+NANOMETER"
)
STAGE_KEYS = {"Not absorbed/max. Q.E.": "not_absorbed", "Detected": "detected"}

# Limits fixed before looking at results. Each entry: (rule, justification).
LIMITS = {
    "photoelectron_bias": (
        "max(0.01, 5 se)",
        "1 % agreement of the photon chain; 5 se covers event sampling.",
    ),
    "arrival_time_bias": (
        "max(0.05 ns, 5 se)",
        "well below the 1 ns camera sample; 5 se covers event sampling.",
    ),
    "image_shape": (
        "mean reduced chi2 <= 2",
        "equal Poisson images give about 1; above 2 the shape differs.",
    ),
    "stage_counts": (
        "max(0.01, 5 sqrt(1/N_ref + 1/N_cand))",
        "1 % floor plus Poisson error of the stage counts.",
    ),
    "transmission_mean": (
        "max(0.01, 5 se)",
        (
            "obdeect weight loss must match the table expectation at 1 %. The expectation "
            "includes the telescope altitude and in-telescope focus terms (sim_telarray "
            "atmospheric_transmission2), as obdeect does."
        ),
    ),
    "trigger": (
        "0 mismatches",
        "trigger decisions must not change under unit-transmission substitution.",
    ),
}


def cherenkov_band(log_text):
    match = BAND_PATTERN.search(log_text)
    if not match:
        raise ValueError("no CHERENKOV WAVELENGTH BAND line in CORSIKA log")
    return float(match.group(1)), float(match.group(2))


def table_lines(text):
    """Split a sim_telarray transmission table into header lines and (wavelength, values) rows."""
    header, rows = [], []
    for line in text.splitlines():
        fields = line.split()
        if fields and fields[0].isdigit():
            rows.append((int(fields[0]), fields))
        else:
            header.append(line)
    return header, rows


def cut_table(text, max_nm):
    """Transmission table with rows at or above max_nm removed (obdeect accepts 0..999 nm)."""
    header, rows = table_lines(text)
    kept = [" ".join(fields) for wavelength, fields in rows if wavelength < max_nm]
    return "\n".join(header + kept) + "\n"


def zero_table(text):
    """Same table with all optical depths zero (unit transmission) and its bottom below the sources.

    The bottom is common.ZERO_DEPTH_BOTTOM_M, not the table's own H2 (2156 m): photons emitted below
    the bottom get zero transmission in sim_telarray.
    """
    return common.zero_depth_table(text)


def photon_statistics(log_text):
    """Per-telescope counts of sim_telarray's 'Photon statistics' block, keyed by stage name."""
    lines = log_text.splitlines()
    start = lines.index("Photon statistics:")
    stats = {}
    for line in lines[start + 1 :]:
        if not line.startswith("  "):
            break
        name, _, values = line.strip().partition(":")
        stats[name.strip()] = [int(value) for value in values.split()]
    return stats


def propagation_summary(log_text):
    match = re.search(
        r"photons in (\d+), out (\d+) \(weight ([\d.e+-]+)\), dropped (\d+)", log_text
    )
    if not match:
        raise ValueError("no propagation summary in propagate log")
    photons_in, photons_out, weight, dropped = match.groups()
    return {
        "photons_in": int(photons_in),
        "photons_out": int(photons_out),
        "weight_out": float(weight),
        "dropped": int(dropped),
    }


def limit_exceeded(bias, limit):
    return int(abs(bias) > limit)


def ratio_limit(reference, candidate):
    """Limit on the relative difference of two stage counts: 1 % floor or 5 Poisson sigma."""
    return max(0.01, 5 * math.sqrt(1 / reference + 1 / candidate))


def expected_transmission(
    samples,
    table_depth,
    band,
    rng_draws,
    telescope_altitude_m=SITE_ALTITUDE_M + TELESCOPE_Z_M,
    pathlen_m=PATHLEN_M,
):
    """Photon-weighted mean transmission for wavelengths drawn as 1/lambda^2 in band.

    Per photon, with airmass 1/|cz| (sim_telarray atmospheric_transmission2, obdeect --site-altitude
    and --pathlen):
        T = exp(-(|tau(zem) - tau_tel| / |cz| + tau_focus)),
        tau_tel = tau(h_t), tau_focus = tau(h_t + pathlen_m) - tau(h_t),
    with h_t = telescope_altitude_m (observation level + telescope z) and tau = table_depth(
    wavelength, height_m), the vertical optical depth above H2 (zero at and below H2).

    samples: list of (zem_m, cz_abs, photons). Wavelengths at or above MAX_LAMBDA give zero.
    Returns (mean, standard error of the weighted mean).
    """
    lo, hi = band
    telescope_depths = {}
    values = []
    total = 0.0
    for (zem_m, cz_abs, photons), u in zip(samples, rng_draws):
        wavelength = math.floor(1.0 / (1.0 / lo - u * (1.0 / lo - 1.0 / hi)))
        if wavelength >= MAX_LAMBDA_NM:
            probability = 0.0
        else:
            row = max(wavelength, 200)
            if row not in telescope_depths:
                tau_tel = table_depth(row, telescope_altitude_m)
                tau_focus = table_depth(row, telescope_altitude_m + pathlen_m) - tau_tel
                telescope_depths[row] = (tau_tel, tau_focus)
            tau_tel, tau_focus = telescope_depths[row]
            optical_depth = abs(table_depth(row, zem_m) - tau_tel) / cz_abs + tau_focus
            probability = math.exp(-optical_depth)
        values.append(photons * probability)
        total += photons
    mean = sum(values) / total
    deviation = [
        value - mean * photons for value, (_, _, photons) in zip(values, samples)
    ]
    standard_error = math.sqrt(sum(d * d for d in deviation)) / total
    return mean, standard_error


def read_bunches(path, fraction, seed):
    """Total photons of every bunch, and a random sample (zem_m, |cz|, photons) of the bunches."""
    import numpy as np
    from eventio.iact import IACTFile

    rng = np.random.default_rng(seed)
    total = 0.0
    sample = []
    with IACTFile(path) as stream:
        for event in stream:
            for bunches in event.photon_bunches.values():
                bunches = np.asarray(bunches)
                if len(bunches) == 0:
                    continue
                total += float(bunches["photons"].astype(float).sum())
                keep = rng.random(len(bunches)) < fraction
                cz = np.sqrt(
                    np.clip(
                        1.0
                        - bunches["cx"][keep].astype(float) ** 2
                        - bunches["cy"][keep].astype(float) ** 2,
                        1e-12,
                        None,
                    )
                )
                sample.extend(
                    zip(
                        bunches["zem"][keep].astype(float) / 100.0,
                        cz,
                        bunches["photons"][keep].astype(float),
                    )
                )
    return total, sample


def main():
    import numpy as np

    from obdeect_validation.adapters import chain_showers as cs
    from obdeect_validation.adapters.extinction_reference import table_reference

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--showers",
        type=Path,
        default=common.SHOWERS,
    )
    parser.add_argument(
        "--transmission",
        type=Path,
        default=common.TRANSMISSION_TABLE,
    )
    parser.add_argument("--site-altitude", type=float, default=SITE_ALTITUDE_M)
    parser.add_argument("--pathlen", type=float, default=PATHLEN_M)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()

    root = Path.cwd().resolve()
    source = args.showers.resolve()
    band = cherenkov_band(
        source.with_suffix(".corsika.log").read_text(errors="replace")
    )
    table_text = args.transmission.read_text(encoding="utf-8")
    cut_path = root / "transmission-below-1000nm.dat"
    cut_path.write_text(cut_table(table_text, MAX_LAMBDA_NM), encoding="utf-8")
    zero_path = root / "zero-transmission.dat"
    zero_path.write_text(zero_table(table_text), encoding="utf-8")

    native = root / "showers.gst"
    propagated = root / "propagated.gst"
    common.command(
        ["obdeect", "stream", "from-corsika", source],
        root / "from-corsika.log",
        stdout_file=native,
    )
    common.command(
        [
            "obdeect",
            "propagation",
            "run",
            "--transmission",
            cut_path,
            "--site-altitude",
            args.site_altitude,
            "--pathlen",
            args.pathlen,
            "--seed",
            args.seed,
            "--mode",
            "weighted",
            "--wavelength-interpolation",
            "floor",
            "--wavelength-range",
            band[0],
            band[1],
            "-i",
            native,
            "-o",
            propagated,
        ],
        root / "propagate.log",
    )
    candidate_input = root / "propagated.iact"
    common.command(
        ["obdeect", "stream", "to-simtel", propagated, "-o", candidate_input],
        root / "to-simtel.log",
    )

    summary = propagation_summary((root / "propagate.log").read_text(errors="replace"))
    linear = root / "propagated-linear.pst"
    common.command(
        [
            "obdeect",
            "propagation",
            "run",
            "--transmission",
            cut_path,
            "--site-altitude",
            args.site_altitude,
            "--pathlen",
            args.pathlen,
            "--seed",
            args.seed,
            "--mode",
            "weighted",
            "--wavelength-interpolation",
            "linear",
            "--wavelength-range",
            band[0],
            band[1],
            "-i",
            native,
            "-o",
            linear,
        ],
        root / "propagate-linear.log",
    )
    linear_summary = propagation_summary(
        (root / "propagate-linear.log").read_text(errors="replace")
    )
    input_weight, sample = read_bunches(source, SAMPLE_FRACTION, args.seed)
    depth = table_reference(cut_path)[1]
    draws = np.random.default_rng(args.seed).random(len(sample)).tolist()
    telescope_altitude = args.site_altitude + TELESCOPE_Z_M
    expected, expected_se = expected_transmission(
        sample,
        depth,
        band,
        draws,
        telescope_altitude_m=telescope_altitude,
        pathlen_m=args.pathlen,
    )
    observed = summary["weight_out"] / input_weight

    direct = cs.simulate(root, source, "direct", bypass=0)
    candidate = cs.simulate(
        root,
        candidate_input,
        "candidate",
        bypass=0,
        extra=["-C", f"atmospheric_transmission={zero_path}"],
    )
    stats_direct = photon_statistics((root / "direct.log").read_text(errors="replace"))
    stats_candidate = photon_statistics(
        (root / "candidate.log").read_text(errors="replace")
    )
    stage_metrics = {}
    for stage, key in STAGE_KEYS.items():
        reference = sum(stats_direct[stage])
        candidate_count = sum(stats_candidate[stage])
        bias = candidate_count / reference - 1.0
        stage_metrics[f"{key}_relative_bias"] = bias
        stage_metrics[f"{key}_limit_exceeded"] = limit_exceeded(
            bias, ratio_limit(reference, candidate_count)
        )

    a, b = cs.camera_events(direct), cs.camera_events(candidate)
    metrics = {
        "reference_camera_events": len(a),
        "candidate_camera_events": len(b),
        "camera_event_id_delta": len(set(a) ^ set(b)),
    }
    metrics.update(cs.compare_camera(a, b, root))
    metrics.update(stage_metrics)
    metrics.update(
        {
            "propagation_dropped_photons": summary["dropped"],
            "weight_ratio_observed": observed,
            "weight_ratio_expected": expected,
            "transmission_mean_difference": observed - expected,
            "linear_weight_ratio": linear_summary["weight_out"] / input_weight,
            "linear_minus_floor_weight_ratio": (
                linear_summary["weight_out"] - summary["weight_out"]
            )
            / input_weight,
            "transmission_mean_limit_exceeded": limit_exceeded(
                observed - expected, max(0.01, 5 * expected_se)
            ),
        }
    )
    details = {
        "band_nm": list(band),
        "site_altitude_m": args.site_altitude,
        "telescope_altitude_m": telescope_altitude,
        "pathlen_m": args.pathlen,
        "seed": args.seed,
        "input_photon_weight": input_weight,
        "propagation": summary,
        "expected_transmission_standard_error": expected_se,
        "camera_comparison_wavelength_interpolation": "floor (sim_telarray convention)",
        "linear_mode_weight_ratio": linear_summary["weight_out"] / input_weight,
        "linear_minus_floor_weight_ratio": (
            linear_summary["weight_out"] - summary["weight_out"]
        )
        / input_weight,
        "sample_bunches": len(sample),
        "limits": LIMITS,
        "direct_photon_statistics": stats_direct,
        "candidate_photon_statistics": stats_candidate,
    }
    common.write_results(
        {key: float(value) for key, value in metrics.items()}, details, root
    )


if __name__ == "__main__":
    main()
