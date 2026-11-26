"""simtools simulation chain against the simtools direct chain, per event.

Reference: ``simtools-simulate-prod --simulation_software corsika sim_telarray`` (sim_telarray applies
all efficiency factors to the CORSIKA photons). Candidate: the same run with the obdeect chain
(``corsika corsika_to_stream obdeect_atmosphere obdeect_optics stream_to_simtel sim_telarray``).
Both use the same CORSIKA input and seeds (50 gamma showers, 100-3000 GeV, LSTN-01 7.0.0, North,
zenith 20 deg, azimuth 180 deg), so events are matched by event id.

Limits (shared with ``chain_showers``): photoelectron relative bias within max(1 %, 5 standard
errors); mean image reduced chi2 at most 2; mean arrival-time bias within max(0.05 ns, 5 standard
errors). Triggers: paired per-shower trigger difference is flagged when |mean| > 5 standard errors
(``chain_optics.trigger_difference``). A trigger mismatch for a shower near threshold is expected
between two chains with different random sequences; it counts only through that 5-sigma limit.

The candidate must read the atmospheric table that the chain prescribes: for bypass 1 the
zero-depth table (obdeect applies the extinction, sim_telarray none); for bypass 3 the site table
is not read by sim_telarray for the atmosphere (no extinction stage in sim_telarray), so the
expected name is given by the case.
"""

import gzip
from pathlib import Path

from obdeect_validation.adapters import chain_optics as co
from obdeect_validation.adapters import chain_showers as cs
from obdeect_validation.adapters import common

#: Legend label of the one telescope of the simtools runs (sim_telarray index 1).
CAMERA_NAMES = {1: "LSTN-01 (model 7.0.0)"}


def _log_text(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as handle:
        return handle.read()


def atmosphere_check(log_path, expected):
    """Return 0 if the log reads exactly the expected atmospheric table, else 1."""
    try:
        common.check_atmosphere(_log_text(log_path), expected)
    except RuntimeError:
        return 1
    return 0


def generated_shower_ids(reduced_path):
    """Event ids of all generated showers (SHOWERS table of the reference reduced event file)."""
    import h5py

    with h5py.File(reduced_path, "r") as handle:
        return [str(int(event_id)) for event_id in handle["SHOWERS"]["event_id"]]


def run_comparison(
    reference_path,
    candidate_path,
    candidate_log,
    expected_atmosphere,
    root,
    reduced_path,
):
    """Compare two sim_telarray outputs of matched runs; return (metrics, details)."""
    metrics = {"candidate_pipeline_complete": 0}
    details = {"reference": str(reference_path), "candidate": str(candidate_path)}
    if not Path(candidate_path).exists():
        details["candidate_blocked"] = "candidate sim_telarray output missing"
        return metrics, details
    metrics["candidate_pipeline_complete"] = 1
    a = cs.camera_events(reference_path)
    b = cs.camera_events(candidate_path)
    metrics.update(
        {
            "reference_camera_events": len(a),
            "candidate_camera_events": len(b),
            # showers triggered in one chain only: recorded, not a limit (threshold statistics)
            "trigger_set_mismatches": len(set(a) ^ set(b)),
        }
    )
    metrics.update(cs.compare_camera(a, b, Path(root), CAMERA_NAMES))
    shower_ids = generated_shower_ids(reduced_path)
    triggers = co.trigger_difference(a, b, shower_ids)
    details["trigger_probability_per_telescope"] = triggers
    metrics["trigger_probability_limit_exceeded"] = int(
        any(r["flag"] for r in triggers.values())
    )
    metrics["atmosphere_table_limit_exceeded"] = atmosphere_check(
        candidate_log, expected_atmosphere
    )
    details["expected_atmosphere"] = expected_atmosphere
    return metrics, details


def _add_arguments(parser):
    parser.add_argument(
        "--reference", type=Path, required=True, help="direct chain sim_telarray output"
    )
    parser.add_argument(
        "--candidate",
        type=Path,
        required=True,
        help="candidate chain sim_telarray output",
    )
    parser.add_argument(
        "--candidate-log",
        type=Path,
        required=True,
        help="candidate sim_telarray stderr/log",
    )
    parser.add_argument(
        "--expected-atmosphere", required=True, help="table name sim_telarray must read"
    )
    parser.add_argument(
        "--reduced",
        type=Path,
        required=True,
        help="reduced event file of the reference run",
    )


def analyse(root, args):
    return run_comparison(
        args.reference.resolve(),
        args.candidate.resolve(),
        args.candidate_log.resolve(),
        args.expected_atmosphere,
        root,
        args.reduced.resolve(),
    )


def main():
    common.results_main(analyse, __doc__.splitlines()[0], add_arguments=_add_arguments)


if __name__ == "__main__":
    main()
