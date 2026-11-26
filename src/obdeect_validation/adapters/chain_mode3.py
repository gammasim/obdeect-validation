"""Functional pass/fail of the production chain with QE envelope thinning and Bypass_Optics=3.

Chain: CORSIKA IACT -> obdeect stream from-corsika -> propagation (QE envelope thinning, one telescope)
-> trace-stream -> stream to-simtel -> sim_telarray Bypass_Optics=3 (camera only, QE table QE / envelope).
Reference: direct sim_telarray Bypass_Optics=0 (full optics) on the same CORSIKA input. Both sides use
LSTN-01 7.0.0 (CTAO.cfg from simtools) and the production atmosphere atmospheric_transmission-2.0.0-OBS-North.dat.

The chain runs through performance_chain.physics_check (the envelope pipe, and its unthinned control).
This adapter adds the pass/fail metrics. Every limit is stated here; none was set to make a run pass.

Metrics (each judged metric is a count of failed checks, expected 0):
- camera_event_id_delta: showers present on one side only.
- photoelectron_bias_limit_exceeded: relative p.e. bias sum(cand - ref) / sum(ref) over telescope events,
  with cluster-robust standard error (clusters = showers). Limit: |bias| <= max(1 %, 3 se). The 1 % floor
  is the camera limit of the validation plan (step 5); the 3 se term is the statistical limit.
- image_shape_limit_exceeded: mean per-pixel chi2 / ndf of the paired images (compare_camera) <= 2
  (chain_showers.SHAPE_CHI2_LIMIT).
- arrival_time_bias_limit_exceeded: mean pixel-time difference (cand - ref) over showers, cluster standard
  error over showers. Limit: |mean| <= max(0.05 ns, 3 se). The floor is chain_showers.TIME_BIAS_FLOOR_NS.
- trigger_pattern_mismatches: showers whose sim_telarray trigger pattern differs (expected 0).
- thinned_photons_limit_exceeded: rows that propagation keeps after envelope thinning (trace-stream
  "photons in" of telescope 1) against the expectation sum_i w'_i E(floor(lambda_i)), where w' are the
  weights of the unthinned weighted propagation (the same input and seed) and E the envelope. Limit:
  |observed - expected| <= 3 sigma, sigma^2 = sum_i w'_i E_i. Thinning keeps a bunch of weight w with
  probability p = p_atm E, so Var = sum w p (1 - p) <= sum w' E: this sigma is an upper bound.
- atmosphere_log_failures: sim_telarray logs (reference and candidate) without the production atmosphere
  line "Read atmospheric transmission data from file ...". CTAO.cfg sets the atmosphere only under
  #if TELESCOPE == 0, so the line is the evidence.
- qe_table_failures: (a) the QE table name of the config is read by both runs; (b) the scaled table is
  found first on the sim_telarray search path (-I table_dir before the config directory); (c) the scaled
  table holds the headroom as its maximum. sim_telarray logs the configured file name only, not its
  directory, so (b) and (c) are checked from the command line and the file, not from the log.
- bypass_mode_failures: the candidate log has the Bypass_Optics=3 line of sim_telarray (camera only).

Event selection: min_photons = 0 on both sides (performance_chain.MIN_PHOTONS, via telescope_arguments).
sim_telarray otherwise skips an event below the config's min_photons (120 for LSTN-01), counted on the
photons of the bunches it reads; a thinned stream has one photon per bunch, so the candidate would drop
events the reference keeps. min_photoelectrons (20 for LSTN-01) is kept and recorded in details.

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.chain_mode3 --config cfg/model/7.0.0/CTAO.cfg \
        --candidate-binary /opt/mode3/sim_telarray
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from obdeect_validation.adapters import chain_optics as co
from obdeect_validation.adapters import chain_showers as cs
from obdeect_validation.adapters import common
from obdeect_validation.adapters import obdeect_cli as cli
from obdeect_validation.adapters import performance_chain as pc

PE_BIAS_FLOOR = 0.01  # relative; camera limit of plan step 5
SIGMA = 3.0  # statistical limit in standard errors
TIME_BIAS_FLOOR_NS = cs.TIME_BIAS_FLOOR_NS
BYPASS3_LINE = "Bypass_Optics=3 (camera only)"
#: Comparison outputs replaced by the clustered metrics of this module.
REPLACED = (
    "photoelectron_relative_bias",
    "photoelectron_relative_bias_se",
    "photoelectron_bias_limit_exceeded",
    "mean_arrival_time_difference_ns",
    "arrival_time_bias_limit_exceeded",
)


def shower_rows(reference, candidate, keys):
    """Per shower: (sum of cand - ref p.e., sum of ref p.e., mean pixel-time difference or None).

    Telescope events are matched by telescope index within each shower; the time difference is the
    mean over the matched telescope events that have pixel times on both sides.
    """
    rows = []
    for key in keys:
        ref_tel = reference[key].get("photoelectrons", {})
        cand_tel = candidate[key].get("photoelectrons", {})
        diff = ref_sum = 0.0
        times = []
        for tel in set(ref_tel) & set(cand_tel):
            a, b = ref_tel[tel], cand_tel[tel]
            diff += float(b["n_pe"]) - float(a["n_pe"])
            ref_sum += float(a["n_pe"])
            if len(a["time"]) and len(b["time"]):
                times.append(float(np.mean(b["time"]) - np.mean(a["time"])))
        rows.append((diff, ref_sum, float(np.mean(times)) if times else None))
    return rows


def clustered_pe_bias(rows):
    """Ratio bias sum(D) / sum(R) and its cluster-robust standard error (delta method, clusters = showers).

    D_k and R_k are the summed differences and reference p.e. of shower k. The residual
    D_k - bias R_k gives se = sqrt(K / (K - 1) * sum r_k^2) / sum R_k.
    """
    diff = np.array([row[0] for row in rows], dtype=float)
    ref = np.array([row[1] for row in rows], dtype=float)
    clusters = len(rows)
    if clusters < 2 or ref.sum() <= 0:
        raise ValueError("p.e. bias needs at least two showers with reference p.e.")
    bias = diff.sum() / ref.sum()
    residual = diff - bias * ref
    se = np.sqrt(clusters / (clusters - 1) * np.sum(residual**2)) / ref.sum()
    return float(bias), float(se)


def clustered_time_bias(rows):
    """Mean over showers of the pixel-time difference, and its standard error over showers."""
    means = [row[2] for row in rows if row[2] is not None]
    if len(means) < 2:
        raise ValueError(
            "arrival-time bias needs at least two showers with pixel times"
        )
    return float(np.mean(means)), float(np.std(means, ddof=1) / np.sqrt(len(means)))


def bias_limit_exceeded(value, se, floor):
    """1 if |value| is above max(floor, SIGMA * se)."""
    return int(abs(value) > max(floor, SIGMA * se))


def thinned_expectation(weights, wavelengths, envelope):
    """Expected surviving rows sum(w' E(floor(lambda))) and the sigma bound sqrt(sum(w' E)).

    weights are the bunch weights of the unthinned weighted propagation (w' = w p_atm), envelope the
    two-column table written by performance_chain.write_envelope (constant on each nm bin).
    """
    nodes = np.asarray(envelope, dtype=float)
    efficiency = np.interp(np.floor(wavelengths), nodes[:, 0], nodes[:, 1])
    expected = float(np.sum(weights * efficiency))
    return expected, float(np.sqrt(np.sum(weights * efficiency)))


def telescope_photons(path, telescope):
    """Bunch weights and wavelengths of one telescope in an obdeect photon stream."""
    from obdeect.stream import Reader  # optional dependency, only for native streams

    weights, wavelengths = [], []
    current = None
    for record in Reader(str(path)):
        if record["type"] == "telescope":
            current = int(record["metadata"]["telescope_id"])
        elif record["type"] == "photons":
            tel = int(record.get("telescope_id", current))
            if tel == telescope:
                columns = dict(record["columns"])
                weights.append(np.asarray(columns["w"], dtype=float))
                wavelengths.append(np.asarray(columns["wavelength"], dtype=float))
    if not weights:
        raise RuntimeError(f"no photons of telescope {telescope} in {path}")
    return np.concatenate(weights), np.concatenate(wavelengths)


def weighted_reference(root, source, inputs):
    """Weights and wavelengths of telescope 1 without thinned output (same input, atmosphere, seed).

    The propagation covers all telescopes (the image's obdeect has no --telescopes option); the
    thinning and wavelength draws are keyed by telescope, so telescope 1 is the same as in the pipe.
    """
    root.mkdir(parents=True, exist_ok=True)
    native, weighted = root / "showers.gst", root / "weighted.gst"
    common.command(
        cli.from_corsika(source), root / "from-corsika.log", stdout_file=native
    )
    common.command(
        cli.propagation_run(
            native,
            weighted,
            transmission=inputs["atmosphere"],
            site_altitude=2156,
            pathlen=28,
            seed=123,
            mode="weighted",
        ),
        root / "propagation.log",
    )
    return telescope_photons(weighted, co.TELESCOPES[0])


def _failed(check):
    """1 if check() raises RuntimeError (a log without the expected line), else 0."""
    try:
        check()
    except RuntimeError:
        return 1
    return 0


def configuration_failures(inputs, files, simtel_argv, texts):
    """Failure counts of the log, table and mode checks (0 = all as intended).

    texts: {'reference': text, 'candidate': text} of the sim_telarray outputs.
    """
    qe_line = f"Reading quantum efficiency from file {inputs['qe'].name}."
    atmosphere = sum(
        _failed(lambda text=text: common.check_atmosphere(text, pc.ATMOSPHERE_FILE))
        for text in texts.values()
    )
    qe = int(qe_line not in texts["candidate"])
    table = Path(files["table_dir"]) / inputs["qe"].name
    order = [str(arg) for arg in simtel_argv]
    table_first = (
        f"-I{files['table_dir']}" in order
        and f"-I{inputs['config_dir']}" in order
        and order.index(f"-I{files['table_dir']}")
        < order.index(f"-I{inputs['config_dir']}")
    )
    scaled = (
        table.is_file()
        and abs(float(np.loadtxt(table, usecols=1).max()) - files["headroom"]) < 1e-6
    )
    qe += int(not table_first) + int(not scaled)
    bypass = int(BYPASS3_LINE not in texts["candidate"])
    return {"atmosphere": atmosphere, "qe_table": qe, "bypass_mode": bypass}


def _add_arguments(parser):
    cli.input_argument(parser, default=common.SHOWERS)
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="sim_telarray config generated by simtools for LSTN-01 7.0.0 (CTAO.cfg)",
    )
    parser.add_argument(
        "--candidate-binary",
        required=True,
        help="sim_telarray with Bypass_Optics=3 (camera only) for the candidate camera stage",
    )


def analyse(root, args):
    """Functional checks of the candidate chain against the reference: (metrics, details)."""
    root = Path(root).resolve()
    source = args.input.resolve()
    pc.CANDIDATE_SIMTEL = str(args.candidate_binary)
    pc.CANDIDATE_OPTIONS = (
        "-v",
    )  # sim_telarray prints the QE table file only when verbose
    inputs, cfg_text = pc.prepare_inputs(root, args.config)
    headroom = pc.headroom_for(inputs["qe_variation"])
    metrics = {"candidate_pipeline_complete": 0}
    details = {
        "settings": {
            "input": str(source),
            "config": str(inputs["config"]),
            "candidate_binary": str(args.candidate_binary),
            "reference_atmosphere": pc.ATMOSPHERE_FILE,
            "candidate_atmosphere": "obdeect propagation, production table below 1000 nm (sim_telarray camera only)",
            "qe_table": inputs["qe"].name,
            "qe_variation": inputs["qe_variation"],
            "headroom": headroom,
            "telescopes": list(co.TELESCOPES),
            "pointing": f"zenith {co.ZENITH_DEG} deg, azimuth {co.AZIMUTH_DEG} deg",
            "min_photons": {"reference": pc.MIN_PHOTONS, "candidate": pc.MIN_PHOTONS},
            "min_photoelectrons": int(
                float(co.cfg_value(cfg_text, "min_photoelectrons"))
            ),
        },
        "factors": co.factor_record(cfg_text, inputs["config_dir"]),
        "counts": {},
        "physics": {},
        "physics_limits": {
            "photoelectron_relative_bias": f"max({PE_BIAS_FLOOR}, {SIGMA:g} se), clustered by shower",
            "arrival_time_bias_ns": f"max({TIME_BIAS_FLOOR_NS}, {SIGMA:g} se over showers)",
            "image_reduced_chi2": f"<= {cs.SHAPE_CHI2_LIMIT:g}",
            "thinned_photons": f"|observed - expected| <= {SIGMA:g} sigma",
        },
    }
    try:
        run = pc.physics_check(
            root, source, inputs, inputs["qe"], headroom, {}, details
        )
        metrics["candidate_pipeline_complete"] = 1
        reference, candidate = run["reference"], run["envelope"]
        keys = sorted(set(reference) & set(candidate), key=int)
        metrics["reference_camera_events"] = len(reference)
        metrics["candidate_camera_events"] = len(candidate)
        metrics["camera_event_id_delta"] = len(set(reference) ^ set(candidate))
        comparison = cs.compare_camera(reference, candidate, root, co.CAMERA_NAMES)
        for key in REPLACED:
            comparison.pop(key, None)
        metrics.update(comparison)
        rows = shower_rows(reference, candidate, keys)
        bias, bias_se = clustered_pe_bias(rows)
        time_bias, time_se = clustered_time_bias(rows)
        metrics.update(
            {
                "photoelectron_relative_bias": bias,
                "photoelectron_relative_bias_se": bias_se,
                "photoelectron_bias_limit_exceeded": bias_limit_exceeded(
                    bias, bias_se, PE_BIAS_FLOOR
                ),
                "arrival_time_bias_ns": time_bias,
                "arrival_time_bias_se_ns": time_se,
                "arrival_time_bias_limit_exceeded": bias_limit_exceeded(
                    time_bias, time_se, TIME_BIAS_FLOOR_NS
                ),
                "showers_clustered": len(rows),
            }
        )
        weights, wavelengths = weighted_reference(root / "weighted", source, inputs)
        expected, sigma = thinned_expectation(
            weights, wavelengths, np.loadtxt(run["files"]["envelope"])
        )
        observed = pc.pipeline_counts(run["envelope_dir"]).get("into_trace_photons")
        if observed is None:
            raise RuntimeError(
                "no trace-stream photon count in the envelope pipeline log"
            )
        metrics.update(
            {
                "thinned_photons_observed": observed,
                "thinned_photons_expected": expected,
                "thinned_photons_sigma": sigma,
                "thinned_photons_limit_exceeded": int(
                    abs(observed - expected) > SIGMA * sigma
                ),
            }
        )
        details["counts"]["thinned_photons"] = {
            "observed": observed,
            "expected": expected,
            "sigma_bound": sigma,
        }
        texts = {
            "reference": (run["work"] / "physics-reference.log").read_text(
                errors="replace"
            ),
            "candidate": (run["envelope_dir"] / "candidate.stdout").read_text(
                errors="replace"
            )
            + (run["envelope_dir"] / "sim_telarray.log").read_text(errors="replace"),
        }
        argv, _ = pc.simtel_candidate(run["envelope_dir"], inputs, run["files"])
        failures = configuration_failures(inputs, run["files"], argv, texts)
        metrics.update(
            {
                "atmosphere_log_failures": failures["atmosphere"],
                "qe_table_failures": failures["qe_table"],
                "bypass_mode_failures": failures["bypass_mode"],
            }
        )
        details["settings"]["candidate_command"] = argv
    except (RuntimeError, OSError, ValueError) as error:  # reported, not hidden
        details["blocked"] = str(error)[:800]
    return metrics, details


def main():
    common.results_main(analyse, "", add_arguments=_add_arguments)


if __name__ == "__main__":
    main()
