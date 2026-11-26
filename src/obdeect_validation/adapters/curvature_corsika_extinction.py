"""Curved CORSIKA showers through obdeect's curvature='corsika' stage, against the curved expectation.

For each CURVED shower file (65 and 75 deg; inputs in /validation/inputs):
  * the EVTH CURVED flag (word 79 = 2) is checked and the CORSIKA Cherenkov band is read from the log;
  * obdeect converts the file (stream from-corsika) and propagates it with curvature='corsika'
    (weighted mode, acceptance stored), with the sim_telarray table cut below 1000 nm;
  * the photon-weighted mean transmission of obdeect (weight out / photons in) is compared with the
    expectation of the same bunches: the curved slant depth of curved_extinction, averaged exactly
    over the band. The plane-parallel expectation (what sim_telarray applies) and the predicted
    ratio curved / plane are reported too.

The only randomness in obdeect's mean is the wavelength draw, so its standard error follows from the
per-bunch wavelength variance of the expectation (curved_extinction.band_expectation).
Telescope altitude: 2156 m + z_t, z_t = the TELESCOPE height of the card layout (curved_showers.layout,
0 cm, the legacy flat layout); obdeect takes the same altitude from the stream's telescope_z_m.
Quantum efficiency and mirror reflectivity are not included: the expectation is of the transmission.
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import atmosphere_chain as chain
from obdeect_validation.adapters import common, curved_showers
from obdeect_validation.adapters import curved_extinction as ce
from obdeect_validation.adapters import obdeect_cli as cli

SITE_ALTITUDE_M = chain.SITE_ALTITUDE_M  # 2156 m, observation level of the CORSIKA runs
PATHLEN_M = chain.PATHLEN_M  # 28 m optical path inside the telescope
SEED = chain.SEED
SIGMA = 5.0
TRANSMISSION_FLOOR = 1.0e-3
RATIO_TOLERANCE = 1.0e-9

# Limits fixed before looking at results. Each entry: (rule, justification).
LIMITS = {
    "curved_transmission": (
        "|observed - curved expectation| <= max(1e-3, 5 se)",
        (
            "The only randomness of obdeect's weighted mean is the wavelength draw; se is the exact "
            "per-bunch wavelength variance of the expectation. 1e-3 is a numerical floor (about ten "
            "times the float32 accumulation and quadrature error), 5 se the sampling error."
        ),
    ),
    "predicted_ratio": (
        "curved / plane >= 1 at every zenith, non-decreasing with zenith",
        (
            "Exact per bunch and wavelength: cos theta(h) >= cos theta_t for h >= h_t, so the curved "
            "slant depth is at most the plane-parallel one and T_curved >= T_plane. The rise with "
            "zenith is the physics of the path; it is checked across the two files, not assumed."
        ),
    ),
    "dropped_photons": (
        "0",
        "Every bunch of a CURVED file has a finite curved transmission; a drop would change the mean.",
    ),
}


class MissingInputError(RuntimeError):
    """A shower input or its CORSIKA log is not in the input directory (exit code 3)."""


def curved_cases():
    """(zenith_deg, energy_min_GeV, energy_max_GeV) of the CURVED fixture showers."""
    return tuple(
        (zenith, low, high)
        for geometry, zenith, low, high in curved_showers.FIXTURE_CASES
        if geometry == "curved"
    )


def shower_file(input_dir, geometry, zenith, low, high) -> Path:
    """Path of an IACT shower file in input_dir, named as curved_showers names it."""
    return Path(input_dir) / curved_showers.shower_name(
        zenith, geometry, energy_min=low, energy_max=high
    )


def check_inputs(input_dir) -> None:
    """Raise MissingInputError unless every curved shower and its CORSIKA log are present."""
    missing = []
    for zenith, low, high in curved_cases():
        source = shower_file(input_dir, "curved", zenith, low, high)
        for path in (source, source.with_suffix(".corsika.log")):
            if not path.is_file():
                missing.append(str(path))
    if missing:
        raise MissingInputError(
            "missing shower input(s): "
            + ", ".join(missing)
            + "; generate with python -m obdeect_validation.adapters.curved_showers run "
            "--geometry curved --zenith Z --energy-min LO --energy-max HI --root WORKDIR, then copy "
            "WORKDIR/<name>.iact and WORKDIR/<name>.corsika.log to the input directory "
            "(fixtures/inputs.json)"
        )


def telescope_altitude_m() -> float:
    """Telescope altitude a.s.l.: site altitude plus the CORSIKA TELESCOPE height z_t (cm, card layout)."""
    telescopes = curved_showers.layout()["telescopes"]
    height_cm = telescopes[0][2]
    return SITE_ALTITUDE_M + height_cm / 100.0


def load_bunches(path):
    """Emission heights (m), |cz| and photon counts of all bunches of an IACT file."""
    zem, cz, photons = [], [], []
    for _, _, columns in curved_showers.read_bunches(path):
        if not {"zem", "cx", "cy", "photons"} <= columns.keys():
            raise ValueError(f"{path} lacks zem, cx, cy or photons columns")
        cx, cy = columns["cx"], columns["cy"]
        zem.append(columns["zem"])
        cz.append(np.sqrt(np.clip(1.0 - cx * cx - cy * cy, 1.0e-12, None)))
        photons.append(columns["photons"])
    if not zem:
        raise ValueError(f"no photon bunches in {path}")
    return np.concatenate(zem), np.concatenate(cz), np.concatenate(photons)


def expectation(
    zem, cz, photons, table, band, telescope_altitude, pathlen_m=PATHLEN_M, weight=None
) -> dict:
    """Photon-weighted expectations of the curved and plane-parallel transmission of the bunches.

    weight: optional detection weight of the wavelength (see curved_extinction.band_expectation);
    None gives the transmission alone, which is what obdeect's acceptance p is.
    """
    per_bunch = ce.band_expectation(
        table, band, zem, cz, telescope_altitude, pathlen_m, weight=weight
    )
    weights = np.asarray(photons, dtype=float)
    curved = ce.weighted_mean(per_bunch["curved"], weights)
    plane = ce.weighted_mean(per_bunch["plane"], weights)
    curved_se = ce.weighted_standard_error(per_bunch["curved_variance"], weights)
    plane_se = ce.weighted_standard_error(per_bunch["plane_variance"], weights)
    ratio, ratio_error = ce.ratio_with_error(curved, curved_se, plane, plane_se)
    return {
        "bunches": int(weights.size),
        "photon_total": float(weights.sum()),
        "live_fraction": per_bunch["live_fraction"],
        "curved": curved,
        "curved_standard_error": curved_se,
        "plane": plane,
        "plane_standard_error": plane_se,
        "predicted_ratio": ratio,
        "predicted_ratio_error": ratio_error,
    }


def ratio_violations(zeniths, ratios) -> int:
    """Predicted ratios below 1, plus decreases between neighbouring zenith angles (sorted)."""
    pairs = sorted(zip(zeniths, ratios))
    below = sum(1 for _, ratio in pairs if ratio < 1.0 - RATIO_TOLERANCE)
    falling = sum(1 for (_, a), (_, b) in itertools.pairwise(pairs) if b < a)
    return int(below + falling)


def tolerance(standard_error: float) -> float:
    """Allowed |observed - expected| of a transmission mean: max(floor, 5 se)."""
    return max(TRANSMISSION_FLOOR, SIGMA * standard_error)


def propagate(root, source, table_path, band, curvature, tag):
    """obdeect stream conversion and the atmosphere stage of one shower file.

    Returns (propagation summary, propagated stream path). Logs go to root/<tag>-*.log.
    """
    root = Path(root)
    native = root / f"{tag}-native.gst"
    propagated = root / f"{tag}-propagated.gst"
    common.command(
        cli.from_corsika(source), root / f"{tag}-from-corsika.log", stdout_file=native
    )
    log = root / f"{tag}-propagate.log"
    common.command(
        cli.propagation_run(
            native,
            propagated,
            transmission=table_path,
            site_altitude=SITE_ALTITUDE_M,
            pathlen=PATHLEN_M,
            seed=SEED,
            mode="weighted",
            store_p=True,
            curvature=curvature,
            wavelength_range=tuple(band),
        ),
        log,
    )
    return chain.propagation_summary(log.read_text(errors="replace")), propagated


def analyse_shower(root, source, zenith, table, table_path) -> dict:
    """Curved CORSIKA file through obdeect's corsika stage, with its expectations."""
    source = Path(source)
    flags = curved_showers.read_header_flags(source)
    if not flags["corsika_curved_flag"]:
        raise ValueError(
            f"{source} is not a CURVED CORSIKA file (EVTH word 79 = {flags['evth_word79']})"
        )
    band = chain.cherenkov_band(
        source.with_suffix(".corsika.log").read_text(errors="replace")
    )
    zem, cz, photons = load_bunches(source)
    h_t = telescope_altitude_m()
    expected = expectation(zem, cz, photons, table, band, h_t)
    work = Path(root) / f"gamma-{zenith:g}deg"
    work.mkdir(parents=True, exist_ok=True)
    summary, _ = propagate(work, source, table_path, band, "corsika", "curved")
    if summary["photons_in"] != zem.size:
        raise ValueError(
            f"obdeect read {summary['photons_in']} bunches, the file holds {zem.size}"
        )
    observed = summary["weight_out"] / expected["photon_total"]
    difference = observed - expected["curved"]
    allowed = tolerance(expected["curved_standard_error"])
    return {
        "zenith_deg": zenith,
        "source": source.name,
        "band_nm": list(band),
        "telescope_altitude_m": h_t,
        "bunches": expected["bunches"],
        "photon_total": expected["photon_total"],
        "propagation": summary,
        "observed_transmission": observed,
        "curved_expectation": expected["curved"],
        "curved_standard_error": expected["curved_standard_error"],
        "difference": difference,
        "tolerance": allowed,
        "limit_exceeded": int(abs(difference) > allowed),
        "plane_expectation": expected["plane"],
        "plane_standard_error": expected["plane_standard_error"],
        "predicted_ratio": expected["predicted_ratio"],
        "predicted_ratio_error": expected["predicted_ratio_error"],
        "observed_over_plane": observed / expected["plane"],
        "observed_over_plane_error": expected["curved_standard_error"]
        / expected["plane"],
        "curved_minus_plane_in_se": (expected["curved"] - expected["plane"])
        / expected["curved_standard_error"],
        "observed_minus_plane_in_se": (observed - expected["plane"])
        / expected["curved_standard_error"],
    }


def analyse_flat(root, source, table, table_path) -> dict:
    """Flat CORSIKA file at 65 deg through obdeect with plane-parallel geometry (information only)."""
    flags = curved_showers.read_header_flags(source)
    if flags["corsika_curved_flag"]:
        raise ValueError(
            f"{source} is flagged CURVED; the flat comparison needs a flat file"
        )
    band = chain.cherenkov_band(
        source.with_suffix(".corsika.log").read_text(errors="replace")
    )
    zem, cz, photons = load_bunches(source)
    expected = expectation(zem, cz, photons, table, band, telescope_altitude_m())
    work = Path(root) / "gamma-65deg-flat"
    work.mkdir(parents=True, exist_ok=True)
    summary, _ = propagate(work, source, table_path, band, "none", "flat")
    observed = summary["weight_out"] / expected["photon_total"]
    return {
        "source": source.name,
        "band_nm": list(band),
        "photon_total": expected["photon_total"],
        "transmitted_weight": summary["weight_out"],
        "transmitted_poisson_error": float(np.sqrt(summary["weight_out"])),
        "observed_transmission": observed,
        "plane_expectation": expected["plane"],
        "plane_standard_error": expected["plane_standard_error"],
        "difference_in_se": (observed - expected["plane"])
        / expected["plane_standard_error"],
        "dropped": summary["dropped"],
    }


def plot(path, rows) -> None:
    """Transmission ratio against zenith and the observed, curved and plane-parallel means."""
    from obdeect_validation.adapters import plotting

    if not rows:
        plotting.message(path, "no curved showers analysed")
        return
    zenith = [row["zenith_deg"] for row in rows]
    fig, axes = plotting.figure(1, 2, (10, 3.8))
    axes[0, 0].errorbar(
        zenith,
        [row["observed_over_plane"] for row in rows],
        yerr=[row["observed_over_plane_error"] for row in rows],
        fmt="o",
        color=plotting.NATIVE,
        label="obdeect / plane expectation",
    )
    axes[0, 0].plot(
        zenith,
        [row["predicted_ratio"] for row in rows],
        "--",
        color=plotting.ANALYTIC,
        label="predicted curved / plane",
    )
    axes[0, 0].set(xlabel="zenith angle [deg]", ylabel="transmission ratio")
    axes[0, 0].legend()
    axes[0, 1].errorbar(
        zenith,
        [row["observed_transmission"] for row in rows],
        yerr=[row["curved_standard_error"] for row in rows],
        fmt="o",
        color=plotting.NATIVE,
        label="obdeect (curved)",
    )
    axes[0, 1].plot(
        zenith,
        [row["curved_expectation"] for row in rows],
        "--",
        color=plotting.ANALYTIC,
        label="curved expectation",
    )
    axes[0, 1].plot(
        zenith,
        [row["plane_expectation"] for row in rows],
        "s",
        color=plotting.REFERENCE,
        label="plane-parallel expectation",
    )
    axes[0, 1].set(xlabel="zenith angle [deg]", ylabel="mean transmission")
    axes[0, 1].legend()
    plotting.save(fig, path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="curvature_corsika_extinction")
    parser.add_argument("--input-dir", type=Path, default=common.VALIDATION / "inputs")
    parser.add_argument("--transmission", type=Path, default=common.TRANSMISSION_TABLE)
    args = parser.parse_args(argv)
    root = Path.cwd().resolve()
    try:
        check_inputs(args.input_dir)
        full_text = args.transmission.read_text(encoding="utf-8")
        table = ce.parse_table(full_text)
        table_path = root / "transmission-below-1000nm.dat"
        table_path.write_text(
            chain.cut_table(full_text, chain.MAX_LAMBDA_NM), encoding="utf-8"
        )
        rows = [
            analyse_shower(
                root,
                shower_file(args.input_dir, "curved", zenith, low, high),
                zenith,
                table,
                table_path,
            )
            for zenith, low, high in curved_cases()
        ]
        flat_file = shower_file(args.input_dir, "flat", 65, 300, 5000)
        flat = None
        if flat_file.is_file():
            flat = analyse_flat(root, flat_file, table, table_path)
    except MissingInputError as exc:
        print(f"missing prerequisite: {exc}", file=sys.stderr)
        return 3
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"curvature_corsika_extinction: {exc}", file=sys.stderr)
        return 2

    metrics = {
        "file_count": len(rows),
        "propagation_dropped_photons": sum(
            row["propagation"]["dropped"] for row in rows
        ),
        "curved_transmission_limit_exceeded": sum(
            row["limit_exceeded"] for row in rows
        ),
        "predicted_ratio_monotonic_violations": ratio_violations(
            [row["zenith_deg"] for row in rows],
            [row["predicted_ratio"] for row in rows],
        ),
    }
    details = {
        "limits": LIMITS,
        "site_altitude_m": SITE_ALTITUDE_M,
        "pathlen_m": PATHLEN_M,
        "seed": SEED,
        "table": str(args.transmission),
        "table_cut_nm": chain.MAX_LAMBDA_NM,
        "showers": rows,
        "flat_65deg_information": flat
        or {
            "available": False,
            "reason": "flat 65 deg shower not in the input directory",
        },
        "notes": [
            "The expectation is the transmission only: no quantum efficiency or mirror weighting.",
            (
                "Telescope zenith for the curved slant is |cz| of each bunch (obdeect corsika mode uses "
                "the frame z axis as the local vertical)."
            ),
            (
                "Information only, not a metric: curved_minus_plane_in_se (per shower) gives the curved "
                "minus plane-parallel mean transmission in standard errors. A synthetic estimate at 75 deg "
                "(20 deg file emission heights, fixture table, 100 to 2000 nm) gives a difference of about "
                "0.2 % of the mean transmission, about 1 standard error at 1e6 bunches and 3 at 1e7, so no "
                "detection metric is set."
            ),
        ],
    }
    plot(root / "curvature-corsika-extinction.png", rows)
    common.write_results(
        {key: float(value) for key, value in metrics.items()}, details, root
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
