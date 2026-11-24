"""Determinism of mixed-direction photons through the atmosphere stage, photon by photon.

Case determinism.mixed-directions. One telescope, N sky-frame photons in every direction class (upward and
downward reaching the plane, upward and downward dropped, horizontal), refraction with a profile,
extinction (weighted), wavelengths not assigned (drawn per photon from the Philox stream keyed on its
index, so the photon index is part of the result). Compared photon by photon (event, telescope, index and
every column bit for bit, endtoend_common.compare_rows):

- CLI `obdeect propagation run` with --threads 1 and --threads 4 (same input);
- Python API propagate with threads 4 and input chunks of 500 rows;
- Python API with one input chunk per photon (chunk size 1);
- shuffled input: the rows are permuted and each keeps its original photon index (one chunk per
  photon with first_index = original index). The output rows of the survivors are then in the shuffled
  order; they are mapped back to the original index (survivors are known from the input classification
  below) and compared with the natural-order output.

Survivors are the photons that reach the plane: upward with 0 < z < plane, downward with z > plane.
Negative control: another seed must change photons (otherwise the comparisons prove nothing).
Limits: zero mismatched photons, zero survivor-count errors, identical summaries.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common
from obdeect_validation.adapters import endtoend_common as e2e
from obdeect_validation.adapters import upward_common as uc

SEED = 20261127
PLANE_Z_M = 500.0
N = 4000
CHUNK_ROWS = 500
MIN_COMPARED = 1000


def version() -> str:
    return common.package_version("obdeect-propagation")


def population():
    rng = np.random.default_rng(SEED)
    classes = rng.integers(
        0, 6, N
    )  # 0,1 upward/downward reaching, 2,3 dropped (above/below), 4 horizontal,
    cols = {
        k: np.empty(N, dtype=np.float32) for k in ("x", "y", "z", "cx", "cy", "cz", "t")
    }
    for i, kind in enumerate(classes):
        zenith = rng.uniform(1.0, 60.0)
        az = rng.uniform(0.0, 360.0)
        th = math.radians(zenith)
        cx, cy = (
            math.sin(th) * math.cos(math.radians(az)),
            math.sin(th) * math.sin(math.radians(az)),
        )
        cz = math.cos(th)
        if kind == 0:  # upward, below the plane
            z, cz_v = rng.uniform(0.0, PLANE_Z_M - 1.0), cz
        elif kind == 1:  # downward, above the plane
            z, cz_v = rng.uniform(PLANE_Z_M + 1.0, PLANE_Z_M + 800.0), -cz
        elif kind == 2:  # upward, above the plane: moving away
            z, cz_v = rng.uniform(PLANE_Z_M + 1.0, PLANE_Z_M + 800.0), cz
        elif kind == 3:  # downward, below the plane
            z, cz_v = rng.uniform(0.0, PLANE_Z_M - 1.0), -cz
        else:  # horizontal
            z, cz_v, cx, cy = (
                rng.uniform(0.0, PLANE_Z_M),
                0.0,
                math.cos(math.radians(az)),
                math.sin(math.radians(az)),
            )
        vals = {
            "x": rng.uniform(-300, 300),
            "y": rng.uniform(-300, 300),
            "z": z,
            "cx": cx,
            "cy": cy,
            "cz": cz_v,
            "t": rng.uniform(0, 5),
        }
        for k, v in vals.items():
            cols[k][i] = v
    cols["w"] = np.ones(N, dtype=np.float32)
    cols["wavelength"] = np.zeros(
        N, dtype=np.float32
    )  # not assigned: drawn per photon index
    reach = np.zeros(N, dtype=bool)
    z = cols["z"].astype(np.float64)
    cz = cols["cz"].astype(np.float64)
    reach = ((cz > 0) & (z > 0) & (z < PLANE_Z_M)) | ((cz < 0) & (z > PLANE_Z_M))
    return cols, reach


def cli(root: Path, stream: Path, out: Path, threads: int, seed: int) -> None:
    argv = [
        "obdeect",
        "propagation",
        "run",
        "--transmission",
        uc.TRANSMISSION,
        "--pathlen",
        "0",
        "--site-altitude",
        uc.SITE_M,
        "--profile",
        uc.PROFILE,
        "--seed",
        seed,
        "--mode",
        "weighted",
        "--threads",
        threads,
        "-i",
        stream,
        "-o",
        out,
    ]
    common.command(argv, root / f"cli-{out.stem}.log")


def api(stream: Path, out: Path, threads: int, seed: int) -> dict:
    from obdeect.propagation import propagate

    return propagate(
        stream,
        out,
        transmission=uc.TRANSMISSION,
        site_altitude_m=uc.SITE_M,
        seed=seed,
        profile=uc.PROFILE,
        mode="weighted",
        threads=threads,
        pathlen_m=0.0,  # as the CLI run of the same case (--pathlen 0)
    )


def _subset(rows: dict, mask) -> dict:
    return {k: (v[mask] if not k.startswith("_") else v) for k, v in rows.items()}


def compare(left: dict, right: dict) -> int:
    result = e2e.compare_rows(left, right)
    return result["photon_mismatches"] + (1 - result["column_sets_equal"])


def analyse(root: Path) -> tuple[dict, dict]:
    cols, reach = population()
    n_reach = int(reach.sum())
    tel = [(0.0, 0.0, PLANE_Z_M)]
    natural = uc.write_stream(
        root / "det-natural.pst",
        tel,
        [(1, cols)],
        chunk_rows=N,
        wavelength_nm=(300.0, 700.0),
    )
    chunked = uc.write_stream(
        root / "det-chunked.pst",
        tel,
        [(1, cols)],
        chunk_rows=CHUNK_ROWS,
        wavelength_nm=(300.0, 700.0),
    )
    per_photon = uc.write_stream(
        root / "det-per-photon.pst",
        tel,
        [(1, cols)],
        indices={1: np.arange(N)},
        wavelength_nm=(300.0, 700.0),
    )
    perm = np.random.default_rng(SEED + 1).permutation(N)
    shuffled_cols = {k: v[perm] for k, v in cols.items()}
    shuffled = uc.write_stream(
        root / "det-shuffled.pst",
        tel,
        [(1, shuffled_cols)],
        indices={1: perm},
        wavelength_nm=(300.0, 700.0),
    )

    outs = {}
    outs["cli1"] = root / "out-cli1.pst"
    cli(root, natural, outs["cli1"], 1, 7)
    outs["cli4"] = root / "out-cli4.pst"
    cli(root, natural, outs["cli4"], 4, 7)
    summaries = {}
    outs["api4_chunked"] = root / "out-api4-chunked.pst"
    summaries["api4_chunked"] = api(chunked, outs["api4_chunked"], 4, 7)
    outs["api1_per_photon"] = root / "out-api1-per-photon.pst"
    summaries["api1_per_photon"] = api(per_photon, outs["api1_per_photon"], 1, 7)
    outs["api_shuffled"] = root / "out-api-shuffled.pst"
    summaries["api_shuffled"] = api(shuffled, outs["api_shuffled"], 1, 7)
    outs["control_seed"] = root / "out-control-seed.pst"
    summaries["control_seed"] = api(natural, outs["control_seed"], 1, 8)
    summaries["api1_natural"] = api(natural, root / "out-api1-natural.pst", 1, 7)

    rows = {k: e2e.photon_rows(v) for k, v in outs.items()}
    base = rows["cli1"]
    survivor_errors = {
        k: int(len(rows[k]["x"]) != n_reach) for k in rows if k != "control_seed"
    }

    mismatches = {}
    mismatches["cli4_vs_cli1"] = compare(base, rows["cli4"])
    mismatches["api4_chunked_vs_cli1"] = compare(base, rows["api4_chunked"])
    mismatches["api1_per_photon_vs_cli1"] = compare(base, rows["api1_per_photon"])
    # Shuffled: output row j is survivor number j in the shuffled input order; map it to its photon.
    survivors_shuffled = [int(i) for i in perm if reach[i]]
    survivors_natural = [int(i) for i in range(N) if reach[i]]
    shuf = rows["api_shuffled"]
    if len(shuf["x"]) == n_reach:
        order = np.argsort(
            np.asarray(survivors_shuffled)
        )  # rows ordered by original index
        # rows of shuffled output sorted by their original index, then compared with natural order
        shuf_sorted = {
            k: (v[order] if not k.startswith("_") else v) for k, v in shuf.items()
        }
        shuf_sorted["index"] = np.asarray(survivors_shuffled, dtype=np.int64)[order]
        natural_rows = dict(base)
        natural_rows["index"] = np.asarray(survivors_natural, dtype=np.int64)
        # Same photons in the same order; the index column holds the photon index in both.
        shuf_sorted["event"] = natural_rows["event"]
        shuf_sorted["telescope"] = natural_rows["telescope"]
        mismatches["shuffled_vs_natural_by_index"] = compare(natural_rows, shuf_sorted)
    else:
        mismatches["shuffled_vs_natural_by_index"] = n_reach
    seed_changed = compare(base, rows["control_seed"])

    summary_equal = int(
        all(
            summaries[k] == summaries["api1_natural"]
            for k in ("api4_chunked", "api1_per_photon")
        )
    )
    metrics = {
        "sample_count": int(N),
        "survivor_count": n_reach,
        "mismatch_count": int(sum(mismatches.values())),
        "empty_comparison": int(len(base["x"]) < MIN_COMPARED),
        "survivor_count_error": int(sum(survivor_errors.values())),
        "summary_mismatch": int(1 - summary_equal),
        "control_seed_change_detected": int(seed_changed > 0),
        "control_seed_mismatch_count": int(seed_changed),
    }
    details = {
        "seed": SEED,
        "obdeect_propagation_version": version(),
        "photons": N,
        "reaching_photons": n_reach,
        "mismatches": mismatches,
        "summaries": summaries,
        "classes": "0 upward reach, 1 downward reach, 2 upward above, 3 downward below, 4 horizontal, 5 unused (downward reach)",
        "reference": "survivors from the input classification (upward 0 < z < plane, downward z > plane)",
    }
    return metrics, details


def main() -> None:
    common.results_main(analyse, __doc__.splitlines()[0])


if __name__ == "__main__":
    main()
