"""Curvature mode rules of the atmosphere stage, through the CLI and the Python API.

Each guard case is a tiny synthetic photon stream (ground frame G or sky frame) with a run key
curvature_applied (absent, 0 or 1) and one curvature setting. The documented rule decides whether
the run must succeed (exit status 0, no exception) or fail (nonzero exit status, ValueError or
RuntimeError). A failure must name the offending key in its message (substring check on stderr
of the CLI and on the text of the Python exception).

Rules (docs/propagation/curvature.md, src/propagation/cpp/include/obdeect/propagation/propagator.hpp):
  mode none        rejects curvature_applied = 1
  mode corsika     needs curvature_applied = 1 and G-frame photons; rejects scattering
  mode spherical   needs curvature_applied = 0 (key present); rejects 1; rejects scattering
  earth_radius_m   must be finite and above 1e6 m
  curvature value  none, corsika or spherical
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common
from obdeect_validation.adapters.refraction_curved_reference import (
    SITE_ALTITUDE_M,
    profile_path,
)

FIXTURES = Path(__file__).resolve().parents[3] / "fixtures"
TRANSMISSION = FIXTURES / "atm_trans_2000_1_10_0_0_2000.dat"
SCATTERING = FIXTURES / "scattering_exponential.dat"
SEED = 1

# (id, frame of the photons, run key curvature_applied or None, curvature setting or None,
#  extra, expected outcome, substring that the error message must contain).
# extra: profile, scattering, earth_radius (float) or None.
CASES = (
    {
        "id": "none-rejects-curvature-applied",
        "frame": "ground",
        "key": 1,
        "curvature": None,
        "extra": None,
        "expect": "error",
        "substring": "curvature_applied",
    },
    {
        "id": "corsika-requires-key",
        "frame": "ground",
        "key": None,
        "curvature": "corsika",
        "extra": None,
        "expect": "error",
        "substring": "curvature_applied",
    },
    {
        "id": "corsika-rejects-zero",
        "frame": "ground",
        "key": 0,
        "curvature": "corsika",
        "extra": None,
        "expect": "error",
        "substring": "curvature_applied",
    },
    {
        "id": "corsika-accepts-one",
        "frame": "ground",
        "key": 1,
        "curvature": "corsika",
        "extra": None,
        "expect": "ok",
        "substring": None,
    },
    {
        "id": "spherical-rejects-one",
        "frame": "sky",
        "key": 1,
        "curvature": "spherical",
        "extra": None,
        "expect": "error",
        "substring": "curvature_applied",
    },
    {
        "id": "spherical-accepts-zero",
        "frame": "sky",
        "key": 0,
        "curvature": "spherical",
        "extra": None,
        "expect": "ok",
        "substring": None,
    },
    {
        "id": "spherical-accepts-profile",
        "frame": "sky",
        "key": 0,
        "curvature": "spherical",
        "extra": {"profile": True},
        "expect": "ok",
        "substring": None,
    },
    {
        "id": "spherical-rejects-scattering",
        "frame": "sky",
        "key": 0,
        "curvature": "spherical",
        "extra": {"scattering": True},
        "expect": "error",
        "substring": "scattering",
    },
    {
        "id": "corsika-rejects-scattering",
        "frame": "ground",
        "key": 1,
        "curvature": "corsika",
        "extra": {"scattering": True},
        "expect": "error",
        "substring": "scattering",
    },
    {
        "id": "corsika-rejects-sky-frame",
        "frame": "sky",
        "key": 1,
        "curvature": "corsika",
        "extra": None,
        "expect": "error",
        "substring": "curvature",
    },
    {
        "id": "unknown-curvature-value",
        "frame": "ground",
        "key": 1,
        "curvature": "flat",
        "extra": None,
        "expect": "error",
        "substring": "curvature",
    },
    {
        "id": "earth-radius-below-1e6",
        "frame": "sky",
        "key": 0,
        "curvature": "spherical",
        "extra": {"earth_radius": 1.0e5},
        "expect": "error",
        "substring": "earth",
    },
    {
        "id": "earth-radius-nan",
        "frame": "sky",
        "key": 0,
        "curvature": "spherical",
        "extra": {"earth_radius": float("nan")},
        "expect": "error",
        "substring": "earth",
    },
)


def _run_metadata(key) -> dict:
    run = {
        "run_number": 1,
        "observation_level_m": SITE_ALTITUDE_M,
        "telescope_x_m": [0.0],
        "telescope_y_m": [0.0],
        "telescope_z_m": [0.0],
        "telescope_r_m": [12.0],
        "refraction_applied": 0,
        "wavelength_min_nm": 250.0,
        "wavelength_max_nm": 700.0,
    }
    if key is not None:
        run["curvature_applied"] = key
    return run


def write_stream(path, frame: str, key, n: int = 40, seed: int = 7) -> None:
    """Tiny stream: n photons per telescope (ground frame G or sky frame), one telescope."""
    from obdeect.stream import Writer

    rng = np.random.default_rng(seed)
    writer = Writer(str(path), "obdeect-validation")
    writer.write("run", _run_metadata(key))
    writer.write("event", {"event_number": 1})
    writer.write(
        "telescope", {"telescope_id": 1, "frame": 1 if frame == "ground" else 0}
    )
    cx = rng.uniform(-0.3, 0.3, n)
    cy = rng.uniform(-0.3, 0.3, n)
    if frame == "ground":
        columns = {
            "x": rng.uniform(-300, 300, n),
            "y": rng.uniform(-300, 300, n),
            "cx": cx,
            "cy": cy,
            "t": rng.uniform(0, 10, n),
            "zem": rng.uniform(2000, 8000, n),
            "w": np.ones(n),
        }
    else:
        columns = {
            "x": rng.uniform(-3000, 3000, n),
            "y": rng.uniform(-3000, 3000, n),
            "z": rng.uniform(3000, 15000, n),
            "cx": cx,
            "cy": cy,
            "cz": -np.sqrt(1.0 - cx**2 - cy**2),
            "t": rng.uniform(0, 10, n),
            "w": np.ones(n),
        }
    writer.write_photons(
        1,
        1 if frame == "ground" else 0,
        0,
        {
            name: np.asarray(values, dtype=np.float32)
            for name, values in columns.items()
        },
    )
    writer.write("telescope_end", {"telescope_id": 1, "photons": float(n)})
    writer.write("event_end", {})
    writer.write("run_end", {})
    writer.close()


def _cli_command(case: dict, stream: Path, output: Path) -> list[str]:
    obdeect = shutil.which("obdeect")
    head = (
        [obdeect, "propagation", "run"]
        if obdeect
        else [sys.executable, "-m", "obdeect.stream.cli", "propagation", "run"]
    )
    extra = case["extra"] or {}
    argv = [
        *head,
        "--transmission",
        str(TRANSMISSION),
        "--site-altitude",
        str(SITE_ALTITUDE_M),
        "--seed",
        str(SEED),
    ]
    if case["curvature"] is not None:
        argv += ["--curvature", case["curvature"]]
    if extra.get("profile"):
        argv += ["--profile", str(profile_path())]
    if extra.get("scattering"):
        argv += ["--scattering", str(SCATTERING), "--scattering-points", "8"]
    if "earth_radius" in extra:
        argv += ["--earth-radius-m", repr(extra["earth_radius"])]
    argv += ["-i", str(stream), "-o", str(output)]
    return [str(arg) for arg in argv]


def run_cli(case: dict, stream: Path, output: Path) -> dict:
    argv = _cli_command(case, stream, output)
    result = subprocess.run(argv, capture_output=True, text=True, check=False)
    return {
        "exit_status": result.returncode,
        "stderr": result.stderr.strip(),
        "outcome": "ok" if result.returncode == 0 else "error",
        "argv": argv,
    }


def run_api(case: dict, stream: Path, output: Path) -> dict:
    from obdeect.propagation import propagate

    extra = case["extra"] or {}
    kwargs = {}
    if case["curvature"] is not None:
        kwargs["curvature"] = case["curvature"]
    if extra.get("profile"):
        kwargs["profile"] = profile_path()
    if extra.get("scattering"):
        kwargs["scattering"] = SCATTERING
        kwargs["scattering_points"] = 8
    if "earth_radius" in extra:
        kwargs["earth_radius_m"] = extra["earth_radius"]
    try:
        propagate(
            stream,
            output,
            transmission=TRANSMISSION,
            site_altitude_m=SITE_ALTITUDE_M,
            seed=SEED,
            **kwargs,
        )
    except (ValueError, RuntimeError) as error:
        return {"outcome": "error", "message": str(error)}
    return {"outcome": "ok", "message": ""}


def main() -> None:
    root = Path.cwd()
    streams_dir = root / "streams"
    streams_dir.mkdir(exist_ok=True)
    rows = []
    unexpected = 0
    mismatch = 0
    for case in CASES:
        stream = streams_dir / f"{case['id']}.pst"
        write_stream(stream, case["frame"], case["key"])
        output_cli = streams_dir / f"{case['id']}-cli.pst"
        output_api = streams_dir / f"{case['id']}-api.pst"
        cli = run_cli(case, stream, output_cli)
        api = run_api(case, stream, output_api)
        for channel, result, text in (
            ("cli", cli, cli["stderr"]),
            ("api", api, api["message"]),
        ):
            if result["outcome"] != case["expect"]:
                unexpected += 1
            elif (
                case["expect"] == "error"
                and case["substring"].lower() not in text.lower()
            ):
                mismatch += 1
        rows.append(
            {
                "id": case["id"],
                "frame": case["frame"],
                "curvature_applied": case["key"],
                "curvature": case["curvature"],
                # Floats such as NaN are recorded as text: strict JSON has no NaN.
                "extra": {
                    key: repr(value) if isinstance(value, float) else value
                    for key, value in (case["extra"] or {}).items()
                },
                "expect": case["expect"],
                "substring": case["substring"],
                "cli_exit_status": cli["exit_status"],
                "cli_outcome": cli["outcome"],
                "cli_stderr": cli["stderr"],
                "api_outcome": api["outcome"],
                "api_message": api["message"],
            }
        )
    metrics = {
        "guard_case_count": len(CASES),
        "guard_unexpected_outcome_count": unexpected,
        "guard_message_mismatch_count": mismatch,
    }
    common.write_results(
        metrics, {"cases": rows, "limits": "documented curvature rules"}
    )


if __name__ == "__main__":
    main()
