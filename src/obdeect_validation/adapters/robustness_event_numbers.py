"""Event numbers of the photon-stream chain: the same photons in ten events with event numbers from 1
to the largest allowed value, through propagation, ray tracing (native CLI and Python API) and
stream to-simtel. LSTN-01 on telescope 1, pointing zenith 20 deg, azimuth 180 deg.

Input: one point source 1000 m above the observation level, placed on the optical axis of the
pointing (x = -1000 sin 20 deg, z = 1000 cos 20 deg; the focal spot is at the camera centre), 1e10
photons into 4 pi in bunches of 10, one telescope of radius 12 m. It gives 36002 sky-frame rows
(about 3.6e5 photons). The rows are written once and copied into each event by the Writer with only
the event number changed; so every event carries the same photons.

Event numbers: 1, 2, 255, 256, 300, 65535, 65536, 100000, 2^31-1 and 2^32-1. The stream format
stores event_number as an integer, the propagation stage accepts up to 2^32-1 (u32), and the
ray tracer needs event numbers below 2^32 (bunch ids). 2^32-1 is therefore the largest allowed
value. The random numbers are keyed on (seed, event, telescope, photon), so each event must give
different camera photons.

Stages per event: propagation run, ray-tracing trace-stream (native CLI), ray_tracing.stream.trace_stream
(Python API, in a subprocess), stream to-simtel (on the native trace output). A stage whose input is
missing (its upstream stage failed) is not run and is listed in details.json as skipped.

Metrics, with gates (tolerance exact (0) for counts; see each justification):
  - stage_failures: executed stage x event runs whose classification is not ok (crash, hang, error).
    Gate exact 0: an event number must never make a stage fail.
  - events_identical_output: pairs of events (i < j) whose traced camera photons agree bit for bit in
    every column, row by row, with at least one photon. Gate exact 0: the random numbers are keyed on
    the event, so different events must give different photons. Pairs need native output of both events.
  - detected_fraction_pull_max: max over events e with native output of
    |f_e - f_1| / sqrt(var_e + var_1), f_e = (camera photons) / (photons into the trace) and
    var_e = f_e (1 - f_e) / n_e, n_e = photons into the trace of event e (the sum of the propagated
    weights, i.e. photons; each camera row has weight 1, so the camera photons are counted as weights).
    Binomial variances: each photon of the input is detected or not, independently, given the optics;
    the draws are independent across events. Baseline event 1. Gate 5 sigma (tolerance 5, expected 0):
    the detected fraction is the same for every event. There are nine comparisons against event 1; a
    two-sided 5 sigma gate has a two-sided normal tail of 5.7e-7 each, so the family-wise false-alarm
    rate is about 9 x 5.7e-7 = 5e-6 (normal approximation of the pull). The gate tests sameness, not a
    physics value, so no reference is involved. If event 1 or any event has no output, the metric is
    PULL_MISSING, which fails the gate.
  - native_python_mismatch: events where both native and Python traces produced output whose camera
    photon counts differ. Gate exact 0: the two interfaces are meant to be equivalent. Bit identity of
    the two outputs is recorded in details.json, not gated.
  - events_without_camera_photons: events without native camera photons (a guard: with empty outputs the
    identity and pull metrics would pass vacuously). Gate exact 0.

Usage: python -m obdeect_validation.adapters.robustness_event_numbers   (writes metrics.json and
details.json in the current directory).
"""

from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common
from obdeect_validation.adapters import endtoend_common as e2e
from obdeect_validation.adapters import obdeect_cli as cli

LST = "LSTN-01"
VERSION = "7.0.0"
COMPILE_SEED = 11
TRACE_SEED = 7
PROPAGATION_SEED = 5
TRACE_THREADS = 2
TELESCOPE = 1
POINTING_ZENITH_DEG = 20.0
POINTING_AZIMUTH_DEG = 180.0
SITE_ALTITUDE_M = 2000.0
SOURCE_DISTANCE_M = 1000.0
SOURCE_PHOTONS = 1e10
SOURCE_BUNCH = 10
SOURCE_RADIUS_M = 12.0
SOURCE_SEED = 3
TRANSMISSION = "fixtures/atm_trans_2000_1_10_0_0_2000.dat"
EVENTS = (1, 2, 255, 256, 300, 65535, 65536, 100000, 2**31 - 1, 2**32 - 1)
BASELINE = EVENTS[0]
#: Gate of detected_fraction_pull_max (sigma).
PULL_GATE = 5.0
#: Value of the pull metric when it cannot be computed (missing output); above PULL_GATE, so it fails.
PULL_MISSING = 1.0e6
#: Python API call of trace_stream, run in a subprocess (argv: model, input, output).
TRACE_SCRIPT = f"""
import sys
from obdeect.ray_tracing import OpticalModel
from obdeect.ray_tracing.stream import trace_stream
model, stream_in, stream_out = sys.argv[1:4]
counts = trace_stream(OpticalModel.load(model), stream_in, stream_out, seed={TRACE_SEED},
                      threads={TRACE_THREADS}, telescopes=[{TELESCOPE}],
                      pointing=({POINTING_ZENITH_DEG}, {POINTING_AZIMUTH_DEG}))
print(counts)
"""


def source_position_m() -> tuple[float, float, float]:
    """Position of the point source on the optical axis (sky frame, x north, z above the level)."""
    zenith = math.radians(POINTING_ZENITH_DEG)
    azimuth = math.radians(POINTING_AZIMUTH_DEG)
    # Along the pointing direction (zenith, azimuth in the array frame, x north, y west). With
    # azimuth 180 deg this is x = -1000 sin 20 deg: the focal spot falls on the camera centre
    # (checked by the camera centroid in details.json).
    return (
        SOURCE_DISTANCE_M * math.sin(zenith) * math.cos(azimuth),
        SOURCE_DISTANCE_M * math.sin(zenith) * math.sin(azimuth),
        SOURCE_DISTANCE_M * math.cos(zenith),
    )


def binomial_fraction(detected: float, trials: float) -> tuple[float, float]:
    """(f, var) of a detected fraction f = detected / trials with binomial variance f (1 - f) / trials."""
    if trials <= 0:
        raise ValueError("binomial_fraction: no input photons")
    fraction = detected / trials
    return fraction, fraction * (1.0 - fraction) / trials


def pull(
    fraction: float, variance: float, reference: float, reference_variance: float
) -> float:
    """|f - f_ref| / sqrt(var + var_ref). Identical fractions with zero variance give 0; a difference
    with zero variance gives PULL_MISSING (infinitely significant)."""
    spread = math.sqrt(variance + reference_variance)
    difference = abs(fraction - reference)
    if spread == 0.0:
        return 0.0 if difference == 0.0 else PULL_MISSING
    return difference / spread


def detected_fraction_pull_max(
    fractions: dict[int, tuple[float, float] | None], baseline: int = BASELINE
) -> tuple[float, dict[int, float | None]]:
    """Largest pull of the events against the baseline; (value, pull per event).

    fractions maps event -> (f, var) or None (no output). A missing baseline or a missing event
    gives PULL_MISSING as the maximum.
    """
    reference = fractions.get(baseline)
    pulls: dict[int, float | None] = {}
    for event, value in fractions.items():
        if event == baseline:
            continue
        if value is None or reference is None:
            pulls[event] = None
        else:
            pulls[event] = pull(value[0], value[1], reference[0], reference[1])
    if reference is None or any(v is None for v in pulls.values()):
        return PULL_MISSING, pulls
    return max(pulls.values(), default=0.0), pulls


def identical_pairs(rows: dict[int, dict | None]) -> list[tuple[int, int]]:
    """Pairs of events whose camera photons are bit-identical (with at least one photon)."""
    found = []
    for left, right in itertools.combinations(sorted(rows), 2):
        a, b = rows[left], rows[right]
        if a is None or b is None or "x" not in a or "x" not in b:
            continue
        if len(a["x"]) == 0:
            continue
        result = e2e.compare_rows(a, b)
        if result["photon_mismatches"] == 0 and result["column_sets_equal"] == 1:
            found.append((left, right))
    return found


def centroid(rows: dict | None) -> list[float] | None:
    """Mean camera position (x, y) [m] of the camera photons; None without photons."""
    if rows is None or "x" not in rows or len(rows["x"]) == 0:
        return None
    return [
        float(np.mean(rows["x"], dtype=np.float64)),
        float(np.mean(rows["y"], dtype=np.float64)),
    ]


def count_photons(rows: dict | None) -> int:
    return 0 if rows is None or "x" not in rows else len(rows["x"])


def weight_sum(rows: dict) -> float:
    return float(np.sum(rows["w"], dtype=np.float64)) if "w" in rows else 0.0


def run_event(
    root: Path, event: int, base: Path, model: Path, camera_map: Path
) -> dict:
    """Propagation, native trace, Python trace and to-simtel of one event; returns its record."""
    name = f"event-{event}"
    record = {"event_number": event, "stages": {}, "skipped": []}
    stream = root / f"{name}.pst"
    e2e.rewrite(base, stream, lambda r: _event_edit(r, event))
    prop = root / f"{name}.prop.pst"
    stage = e2e.run_stage(
        cli.propagation_run(
            stream,
            prop,
            transmission=common.VALIDATION / TRANSMISSION,
            pathlen="0",
            site_altitude=SITE_ALTITUDE_M,
            wavelength_range=(300, 600),
            seed=PROPAGATION_SEED,
            mode="weighted",
        ),
        log=root / f"{name}.propagation.log",
    )
    record["stages"]["propagation"] = e2e.describe(stage)
    if e2e.classify(stage) != "ok":
        record["skipped"] = ["trace-native", "trace-python", "to-simtel"]
        return record

    native = root / f"{name}.trace-native.pst"
    stage = e2e.run_stage(
        cli.trace_stream(
            model,
            prop,
            native,
            seed=TRACE_SEED,
            threads=TRACE_THREADS,
            telescopes=str(TELESCOPE),
            pointing_zenith_deg=POINTING_ZENITH_DEG,
            pointing_azimuth_deg=POINTING_AZIMUTH_DEG,
        ),
        log=root / f"{name}.trace-native.log",
    )
    record["stages"]["trace-native"] = e2e.describe(stage)
    native_ok = e2e.classify(stage) == "ok" and native.is_file()

    python_out = root / f"{name}.trace-python.pst"
    stage = e2e.run_stage(
        [sys.executable, "-c", TRACE_SCRIPT, model, prop, python_out],
        log=root / f"{name}.trace-python.log",
    )
    record["stages"]["trace-python"] = e2e.describe(stage)
    python_ok = e2e.classify(stage) == "ok" and python_out.is_file()

    prop_rows = e2e.photon_rows(prop)
    record["input_photons"] = weight_sum(prop_rows)
    record["_native"] = e2e.photon_rows(native) if native_ok else None
    record["_python"] = e2e.photon_rows(python_out) if python_ok else None

    if native_ok:
        simtel = root / f"{name}.simtel.iact"
        stage = e2e.run_stage(
            cli.obdeect(
                "stream", "to-simtel", camera_map=camera_map, o=simtel, extra=[native]
            ),
            log=root / f"{name}.to-simtel.log",
        )
        record["stages"]["to-simtel"] = e2e.describe(stage)
    else:
        record["skipped"].append("to-simtel")
    return record


def _event_edit(record: dict, event: int):
    if record["type"] == "event":
        return [("meta", "event", {"event_number": event})]
    return None


def analyse(root: Path) -> tuple[dict, dict]:
    camera_map = root / "camera-map.txt"
    camera_map.write_text("1 1 0.0\n")
    model = root / "model.json"
    common.command(
        cli.compile_model(
            model, LST, VERSION, seed=COMPILE_SEED, zenith_deg=POINTING_ZENITH_DEG
        ),
        root / "compile.log",
    )

    base = root / "base.pst"
    x, y, z = source_position_m()
    common.command(
        cli.obdeect(
            "sources",
            "run",
            "point-source",
            position=f"{x:.6f},{y:.6f},{z:.6f}",
            photons=f"{SOURCE_PHOTONS:g}",
            bunch_size=SOURCE_BUNCH,
            telescope_radius=SOURCE_RADIUS_M,
            spectrum="fixed:wavelength_nm=400",
            observation_level=SITE_ALTITUDE_M,
            atmosphere_id=1,
            seed=SOURCE_SEED,
            o=base,
        ),
        root / "point-source.log",
    )

    records = [run_event(root, event, base, model, camera_map) for event in EVENTS]

    # e2e.describe starts with the classification word: "ok rc=0 ..." for a successful stage.
    stage_failures = sum(
        1
        for record in records
        for stage in record["stages"].values()
        if not stage.startswith("ok ")
    )

    native = {r["event_number"]: r["_native"] for r in records}
    python = {r["event_number"]: r["_python"] for r in records}
    fractions: dict[int, tuple[float, float] | None] = {}
    for record in records:
        event = record["event_number"]
        rows = record["_native"]
        if rows is None or record["input_photons"] <= 0:
            fractions[event] = None
            continue
        fractions[event] = binomial_fraction(weight_sum(rows), record["input_photons"])
    pull_max, pulls = detected_fraction_pull_max(fractions)
    identical = identical_pairs(native)

    both = [
        event
        for event in EVENTS
        if native[event] is not None and python[event] is not None
    ]
    python_mismatch = [
        event
        for event in both
        if count_photons(native[event]) != count_photons(python[event])
    ]
    bit_identical_native_python = [
        event
        for event in both
        if e2e.compare_rows(native[event], python[event])["photon_mismatches"] == 0
    ]

    metrics = {
        "stage_failures": stage_failures,
        "events_identical_output": len(identical),
        "detected_fraction_pull_max": float(pull_max),
        "native_python_mismatch": len(python_mismatch),
        "events_without_camera_photons": sum(
            count_photons(native[event]) == 0 for event in EVENTS
        ),
    }
    details = {
        "events": [
            {
                "event_number": record["event_number"],
                "stages": record["stages"],
                "skipped": record["skipped"],
                "input_photons": record["input_photons"],
                "camera_photons_native": count_photons(record["_native"]),
                "camera_photons_python": count_photons(record["_python"]),
                "detected_fraction": fractions[record["event_number"]][0]
                if fractions[record["event_number"]] is not None
                else None,
                "camera_centroid_m_native": centroid(record["_native"]),
                "pull_vs_baseline": pulls.get(record["event_number"]),
                "native_python_bit_identical": record["event_number"]
                in bit_identical_native_python,
            }
            for record in records
        ],
        "identical_pairs": [list(pair) for pair in identical],
        "native_python_count_mismatch_events": python_mismatch,
        "baseline_event": BASELINE,
        "pull_gate_sigma": PULL_GATE,
        "pull_missing_value": PULL_MISSING,
        "source_position_m": [x, y, z],
    }
    return metrics, details


def main():
    common.results_main(analyse, __doc__.split("\n")[0])


if __name__ == "__main__":
    main()
