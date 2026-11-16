"""Curved CORSIKA showers through obdeect's corsika stage and sim_telarray, against the curved prediction.

Reference: sim_telarray with its own atmosphere on the CURVED CORSIKA photons (as atmosphere_chain
does), telescopes pointing at the shower zenith.
Candidate: obdeect propagation (curvature='corsika', weighted), to-simtel, then sim_telarray with the
zero-optical-depth table (unit transmission), as atmosphere_chain's candidate does.
Unit: sim_telarray with the zero-optical-depth table on the CORSIKA photons themselves. It gives the
detection count without atmosphere, so the plane-parallel check is a ratio and needs no absolute
normalisation of sim_telarray's quantum efficiency.

Checks per zenith angle (65 and 75 deg):
  1. Detected photo-electron ratio candidate / reference against the prediction from the same bunches.
     The prediction weights each wavelength by the detection weight of the LST-1 configuration:
     quantum efficiency (high-QE pixels) times mirror reflectivity. sim_telarray's 'Detected' count
     carries both, so an unweighted (transmission-only) prediction is biased high at large zenith
     angles (details: predicted_ratio_transmission_only).
  2. sim_telarray is plane-parallel: its 'Not absorbed/max. Q.E.' count of the reference over that of
     the unit run, against the plane-parallel prediction weighted by the detection weight. That stage
     draws against sim_telarray's element quantum_efficiency array, which is Q.E. x mirror
     reflectivity (sim_absorb.c; mc_aux.h), so the weight is the same as in check 1.
Camera event matching and trigger decisions are reported in the details only: the candidate is
curved by design, so trigger decisions of borderline events may change at high zenith.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import atmosphere_chain as chain
from obdeect_validation.adapters import chain_showers, common, curved_showers
from obdeect_validation.adapters import curvature_corsika_extinction as curved
from obdeect_validation.adapters import curved_extinction as ce
from obdeect_validation.adapters.curvature_corsika_extinction import MissingInputError

SITE_ALTITUDE_M = chain.SITE_ALTITUDE_M
SIGMA = 5.0
RATIO_FLOOR = 0.01
DETECTED = "Detected"
NOT_ABSORBED = "Not absorbed/max. Q.E."
SIMTEL_DEFINES = ("-DCTA_NORTH", "-DLST1", "-DNUM_TELESCOPES=3")
QE_TABLE = "cfg/CTA/qe_lst1_20200318_high+low.dat"  # LST-1 config: first column = high-QE pixels
MIRROR_TABLE = (
    "cfg/CTA/ref_LST_2020-04-23.dat"  # LST mirror reflectivity, column 2 (fraction)
)
# LST camera window (camera_filter in CTA-PROD6-LST.cfg). sim_config.c multiplies it into the
# element quantum_efficiency array, at the bin centre (rpolate_1d at lambda + 0.5).
WINDOW_TABLE = "cfg/CTA/transmission_lst_window_No7-10_ave.dat"

# Limits fixed before looking at results. Each entry: (rule, justification).
LIMITS = {
    "pe_ratio_prediction": (
        "|measured - predicted| <= max(0.01, 5 se_ratio), measured and predicted = detected p.e. ratio",
        (
            "The 1 % floor is that of the 20 deg chain case (atmosphere_chain). The prediction carries "
            "the detection weight (quantum efficiency x mirror reflectivity), the quantity sim_telarray "
            "counts as 'Detected'. se_ratio combines the Poisson errors of both detected counts with "
            "the wavelength-draw errors of the two programs, so 5 se covers the sampling."
        ),
    ),
    "simtel_plane_parallel": (
        "|N_ref / N_unit - E_qe(plane-parallel)| <= max(0.01, 5 se), N = 'Not absorbed/max. Q.E.'",
        (
            "sim_telarray's own atmosphere against the plane-parallel expectation of the same bunches, "
            "weighted by the quantum efficiency that this stage counts (the ratio to the unit run removes "
            "the normalisation). The curved-to-plane difference of the mean transmission at 75 deg is "
            "about 0.2 % (synthetic estimate, curvature_corsika_extinction notes), below the 1 % floor: "
            "the check bounds a plane-parallel sim_telarray at the 1 % level and does not resolve the "
            "curved correction."
        ),
    ),
}


def read_two_column(path) -> np.ndarray:
    """(wavelength nm, value) rows of a sim_telarray table; comments and extra columns are ignored."""
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        fields = line.split("#", 1)[0].split()
        if len(fields) < 2:
            continue
        try:
            rows.append((float(fields[0]), float(fields[1])))
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"no wavelength rows in {path}")
    table = np.array(sorted(rows), dtype=float)
    if np.any(np.diff(table[:, 0]) <= 0.0):
        raise ValueError(f"wavelengths of {path} must be strictly increasing")
    return table


class DetectionWeights:
    """Wavelength weights of the sim_telarray detection stages, from the LST-1 configuration.

    quantum: quantum efficiency (zero outside the table, which ends at 800 nm with zero QE).
    detection: element efficiency of sim_telarray (sim_config.c, sim_absorb.c): quantum efficiency
    x mirror reflectivity x camera window transmission (evaluated at lambda + 0.5, as sim_config.c
    does; zero below 259 nm). Every factor that sim_telarray multiplies into el->quantum_efficiency
    before the photon draw is included. The funnel efficiency is applied after the draw (per pixel
    type, in the photoelectron stage), so it is not in this weight; the 'Not absorbed' stage does
    not reach it. The mirror table must cover every wavelength where the quantum efficiency is
    non-zero; outside the quantum table the weight is zero by definition, so a wider Cherenkov band
    (CORSIKA: 100 to 2000 nm) is accepted there.
    """

    def __init__(self, reference=common.REFERENCE):
        reference = Path(reference)
        self._qe = read_two_column(reference / QE_TABLE)
        self._mirror = read_two_column(reference / MIRROR_TABLE)
        self._window = read_two_column(reference / WINDOW_TABLE)
        if np.any(self._mirror[:, 1] < 0.0) or np.any(self._mirror[:, 1] > 1.0):
            raise ValueError("mirror reflectivity must be a fraction between 0 and 1")

    def quantum(self, nm):
        nm = np.asarray(nm, dtype=float)
        return np.interp(nm, self._qe[:, 0], self._qe[:, 1], left=0.0, right=0.0)

    def detection(self, nm):
        nm = np.asarray(nm, dtype=float)
        quantum = self.quantum(nm)
        active = nm[quantum > 0.0]
        if active.size and (
            active.min() < self._mirror[0, 0] or active.max() > self._mirror[-1, 0]
        ):
            raise ValueError(
                "mirror reflectivity table does not cover the wavelength band where the "
                "quantum efficiency is non-zero"
            )
        reflectivity = np.interp(nm, self._mirror[:, 0], self._mirror[:, 1])
        window = np.interp(
            nm + 0.5, self._window[:, 0], self._window[:, 1], left=0.0, right=0.0
        )
        return quantum * reflectivity * window


def run_simtel(work, input_file, name, zenith_deg, extra=()):
    """sim_telarray on one input with the telescopes pointing at the shower zenith (as chain_showers)."""
    return common.sim_telarray(
        work,
        input_file,
        name,
        {
            "maximum_telescopes": 3,
            "altitude": SITE_ALTITUDE_M,
            "telescope_theta": zenith_deg,
            "telescope_phi": 180,
            "Bypass_Optics": 0,
        },
        defines=SIMTEL_DEFINES,
        extra=tuple(extra),
    )


def stage_total(stats, stage) -> float:
    """Photons (or p.e.) of one 'Photon statistics' stage summed over telescopes."""
    if stage not in stats:
        raise ValueError(f"no '{stage}' line in the sim_telarray photon statistics")
    return float(sum(stats[stage]))


def event_id_delta(left, right) -> int:
    """Camera events present in only one of the two runs."""
    return len(set(left) ^ set(right))


def trigger_mismatches(left, right) -> int:
    """Events present in both runs whose trigger pattern differs (informational)."""
    common_ids = set(left) & set(right)
    return sum(
        1
        for key in common_ids
        if left[key].get("trigger_information", {}).get("trigger_pattern")
        != right[key].get("trigger_information", {}).get("trigger_pattern")
    )


def ratio_check(measured, predicted, predicted_error, reference, candidate) -> dict:
    """Ratio of two sim_telarray counts against a prediction; limit from all error sources.

    predicted_error: standard error of the prediction (draws of both programs); the Poisson errors of
    the two counts are added in quadrature.
    """
    relative = math.hypot(
        predicted_error / predicted,
        math.sqrt(1.0 / reference + 1.0 / candidate),
    )
    error = predicted * relative
    limit = max(RATIO_FLOOR, SIGMA * error)
    difference = measured - predicted
    return {
        "measured_ratio": measured,
        "predicted_ratio": predicted,
        "ratio_error": error,
        "difference": difference,
        "limit": limit,
        "limit_exceeded": int(abs(difference) > limit),
        "relative_difference": difference / predicted,
    }


def analyse(root, source, zenith, table, cut_path, zero_path, weights) -> dict:
    """Reference, unit and candidate sim_telarray runs of one curved shower file, with the checks."""
    source = Path(source)
    flags = curved_showers.read_header_flags(source)
    if not flags["corsika_curved_flag"]:
        raise ValueError(f"{source} is not a CURVED CORSIKA file")
    band = chain.cherenkov_band(
        source.with_suffix(".corsika.log").read_text(errors="replace")
    )
    zem, cz, photons = curved.load_bunches(source)
    h_t = curved.telescope_altitude_m()
    transmission = curved.expectation(zem, cz, photons, table, band, h_t)
    detected = curved.expectation(
        zem, cz, photons, table, band, h_t, weight=weights.detection
    )
    # The 'Not absorbed/max. Q.E.' stage draws against sim_telarray's element quantum_efficiency
    # array, which is Q.E. x mirror reflectivity (sim_absorb.c), so its prediction uses the
    # detection weight as well.
    quantum = curved.expectation(
        zem, cz, photons, table, band, h_t, weight=weights.detection
    )
    for name, result in (("detection", detected), ("quantum", quantum)):
        if abs(result["live_fraction"] - 1.0) > 1e-9:
            raise ValueError(f"{name} weight does not vanish at and above 1000 nm")
    unit_standard_error = ce.weighted_standard_error(
        np.full(photons.size, ce.unit_variance(band, weights.detection)), photons
    )

    work = Path(root) / f"gamma-{zenith:g}deg"
    work.mkdir(parents=True, exist_ok=True)
    reference = run_simtel(work, source, "reference", zenith)
    run_simtel(
        work,
        source,
        "unit",
        zenith,
        extra=["-C", f"atmospheric_transmission={zero_path}"],
    )
    summary, propagated = curved.propagate(
        work, source, cut_path, band, "corsika", "candidate"
    )
    candidate_input = work / "candidate.iact"
    common.command(
        ["obdeect", "stream", "to-simtel", propagated, "-o", candidate_input],
        work / "to-simtel.log",
    )
    candidate = run_simtel(
        work,
        candidate_input,
        "candidate",
        zenith,
        extra=["-C", f"atmospheric_transmission={zero_path}"],
    )

    def stats(name):
        return chain.photon_statistics(
            (work / f"{name}.log").read_text(errors="replace")
        )

    stats_reference = stats("reference")
    stats_unit = stats("unit")
    stats_candidate = stats("candidate")
    detected_reference = stage_total(stats_reference, DETECTED)
    detected_candidate = stage_total(stats_candidate, DETECTED)
    not_absorbed_reference = stage_total(stats_reference, NOT_ABSORBED)
    not_absorbed_unit = stage_total(stats_unit, NOT_ABSORBED)
    not_absorbed_candidate = stage_total(stats_candidate, NOT_ABSORBED)
    if not (detected_reference > 0.0 and not_absorbed_unit > 0.0):
        raise ValueError(f"sim_telarray detected no photons for {source.name}")

    pe = ratio_check(
        detected_candidate / detected_reference,
        detected["predicted_ratio"],
        detected["predicted_ratio_error"],
        detected_reference,
        detected_candidate,
    )
    plane_error = math.hypot(
        quantum["plane_standard_error"] / quantum["plane"],
        unit_standard_error,
    )
    simtel_plane = ratio_check(
        not_absorbed_reference / not_absorbed_unit,
        quantum["plane"],
        quantum["plane"] * plane_error,
        not_absorbed_reference,
        not_absorbed_unit,
    )
    cameras_reference = chain_showers.camera_events(reference)
    cameras_candidate = chain_showers.camera_events(candidate)
    return {
        "zenith_deg": zenith,
        "source": source.name,
        "band_nm": list(band),
        "telescope_altitude_m": h_t,
        "bunches": transmission["bunches"],
        "photon_total": transmission["photon_total"],
        "propagation": summary,
        "transmission": {
            "curved": transmission["curved"],
            "curved_standard_error": transmission["curved_standard_error"],
            "plane": transmission["plane"],
            "plane_standard_error": transmission["plane_standard_error"],
        },
        "predicted_ratio_transmission_only": transmission["predicted_ratio"],
        "predicted_ratio_detection_weighted": detected["predicted_ratio"],
        "predicted_ratio_detection_error": detected["predicted_ratio_error"],
        "detected_reference": detected_reference,
        "detected_candidate": detected_candidate,
        "detected_unit": stage_total(stats_unit, DETECTED),
        "not_absorbed_reference": not_absorbed_reference,
        "not_absorbed_unit": not_absorbed_unit,
        "not_absorbed_candidate": not_absorbed_candidate,
        "candidate_over_unit_not_absorbed_predicted": quantum["curved"],
        "candidate_over_unit_not_absorbed_measured": not_absorbed_candidate
        / not_absorbed_unit,
        "pe_ratio": pe,
        "simtel_plane": simtel_plane,
        "camera_events_reference": len(cameras_reference),
        "camera_events_candidate": len(cameras_candidate),
        "camera_event_id_delta": event_id_delta(cameras_reference, cameras_candidate),
        "trigger_pattern_mismatches": trigger_mismatches(
            cameras_reference, cameras_candidate
        ),
    }


def plot(path, rows) -> None:
    from obdeect_validation.adapters import plotting

    if not rows:
        plotting.message(path, "no curved showers analysed")
        return
    zenith = [row["zenith_deg"] for row in rows]
    fig, axes = plotting.figure(1, 2, (10, 3.8))
    pe = [row["pe_ratio"] for row in rows]
    axes[0, 0].errorbar(
        zenith,
        [item["measured_ratio"] for item in pe],
        yerr=[item["ratio_error"] for item in pe],
        fmt="o",
        color=plotting.NATIVE,
        label="measured (candidate / reference)",
    )
    axes[0, 0].plot(
        zenith,
        [item["predicted_ratio"] for item in pe],
        "--",
        color=plotting.ANALYTIC,
        label="predicted curved / plane (detection-weighted)",
    )
    axes[0, 0].plot(
        zenith,
        [row["predicted_ratio_transmission_only"] for row in rows],
        ":",
        color=plotting.LIMIT,
        label="transmission only (not used)",
    )
    axes[0, 0].set(xlabel="zenith angle [deg]", ylabel="detected p.e. ratio")
    axes[0, 0].legend()
    plane = [row["simtel_plane"] for row in rows]
    axes[0, 1].errorbar(
        zenith,
        [item["measured_ratio"] for item in plane],
        yerr=[item["ratio_error"] for item in plane],
        fmt="o",
        color=plotting.REFERENCE,
        label="sim_telarray not absorbed (reference / unit)",
    )
    axes[0, 1].plot(
        zenith,
        [item["predicted_ratio"] for item in plane],
        "--",
        color=plotting.ANALYTIC,
        label="plane-parallel prediction (Q.E. weighted)",
    )
    axes[0, 1].set(xlabel="zenith angle [deg]", ylabel="ratio")
    axes[0, 1].legend()
    plotting.save(fig, path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="chain_curvature")
    parser.add_argument("--input-dir", type=Path, default=common.VALIDATION / "inputs")
    parser.add_argument("--transmission", type=Path, default=common.TRANSMISSION_TABLE)
    args = parser.parse_args(argv)
    root = Path.cwd().resolve()
    try:
        curved.check_inputs(args.input_dir)
        full_text = args.transmission.read_text(encoding="utf-8")
        table = ce.parse_table(full_text)
        weights = DetectionWeights()
        cut_path = root / "transmission-below-1000nm.dat"
        cut_path.write_text(
            chain.cut_table(full_text, chain.MAX_LAMBDA_NM), encoding="utf-8"
        )
        zero_path = root / "zero-transmission.dat"
        zero_path.write_text(chain.zero_table(full_text), encoding="utf-8")
        rows = [
            analyse(
                root,
                curved.shower_file(args.input_dir, "curved", zenith, low, high),
                zenith,
                table,
                cut_path,
                zero_path,
                weights,
            )
            for zenith, low, high in curved.curved_cases()
        ]
    except MissingInputError as exc:
        print(f"missing prerequisite: {exc}", file=sys.stderr)
        return 3
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"chain_curvature: {exc}", file=sys.stderr)
        return 2

    metrics = {
        "file_count": len(rows),
        "pe_ratio_vs_prediction_limit_exceeded": sum(
            row["pe_ratio"]["limit_exceeded"] for row in rows
        ),
        "simtel_plane_parallel_limit_exceeded": sum(
            row["simtel_plane"]["limit_exceeded"] for row in rows
        ),
        "propagation_dropped_photons": sum(
            row["propagation"]["dropped"] for row in rows
        ),
    }
    details = {
        "limits": LIMITS,
        "site_altitude_m": SITE_ALTITUDE_M,
        "pathlen_m": chain.PATHLEN_M,
        "seed": chain.SEED,
        "simtel_binary": str(common.SIMTEL_BINARY),
        "transmission_table": str(args.transmission),
        "quantum_efficiency_table": QE_TABLE,
        "mirror_reflectivity_table": MIRROR_TABLE,
        "showers": rows,
        "notes": [
            "Detected = sum of the 'Detected' stage over telescopes (photo-electron stage).",
            (
                "Not absorbed/max. Q.E. counts the quantum efficiency once and no mirror reflectivity; "
                "the unit run (zero optical depth) gives its normalisation."
            ),
            "Camera events and triggers are informational (candidate is curved by design).",
            "The light-guide efficiency table is averaged over wavelength and is not a wavelength weight.",
        ],
    }
    plot(root / "chain-curvature-high-zenith.png", rows)
    common.write_results(
        {key: float(value) for key, value in metrics.items()}, details, root
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
