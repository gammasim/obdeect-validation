"""Per-stage CPU cost per photon against a stored per-host baseline (case performance.regression).

Five stages. Every child is timed with performance.run_timed (os.wait4: ru_utime + ru_stime of the
child, threads 1). A stage's CPU is the sum over its children; the cost is CPU / photons.

  propagation           propagation run, extinction + refraction (fixture atmosphere and profile),
                        5e6 sky-frame rows
  propagation_scattering  propagation run --scattering (fixture table), 2e5 sky-frame rows at the
                        table's wavelengths (345-395 nm)
  trace_stream          ray-tracing trace-stream of compiled LSTN-01, MSTN-01, SSTS-01 models (7.0.0,
                        zenith 20 deg), 2e6 sphere rows: per telescope, half on-axis and half 2 deg
                        off-axis, spread over the dish of the model (its primary radius)
  stream                stream cat (to stdout) and stream to-simtel of 1e7 ground-plane rows; the
                        stage CPU is the sum of both
  flasher_dish          sources run flasher --type dish (disc target), 2e9 photons into 4 pi; the
                        photons are the rows written (only photons that reach the disc are written)

Inputs are generated with the obdeect Python Writer (untimed) in a work directory and read back with
`obdeect stream cat --summary` for the photon counts. Models are compiled untimed from simulation-models.

Modes:
  baseline --output FILE   measure and write the host baseline (fixtures/performance/<host key>.json)
  compare [--baseline F]   measure, compare with the baseline of this host, write metrics.json and
                           details.json of the case (default baseline: fixtures/performance/<host key>)

Host key: platform (darwin, linux), machine and CPU brand, lower case, runs of non-alphanumerics
replaced by '-' (macOS: sysctl machdep.cpu.brand_string; Linux: /proc/cpuinfo 'model name').

Repetitions: REPS (5). The stage's best time per photon is the minimum over the repetitions. The
relative spread of a stage is (max - min) / min over its repetitions (SPREAD_K below).

Gates (metrics of the case):
  regressions = number of stages with ratio > 1 + max(FLOOR, K x spread), with ratio = new best /
    baseline best and spread = larger of the new and the baseline relative spread. FLOOR = 0.05: a
    5 % slowdown of the best time is the smallest change the gate is meant to report (a stated
    requirement). K = 3: a 3-sigma-like allowance for the repetition scatter of the stage, so that the
    tolerance is max(0.05, 3 x spread) as the case states. Expected 0: the case compares the code it
    runs with a baseline made from the same code on this host (the ratios are ~1 by construction), so
    any ratio above the tolerance would be a real change, and a change is what the gate must catch.
  baseline_missing = number of measured stages with no entry in the baseline of this host. Expected
    0: the baseline must hold every stage (a missing baseline cannot be compared, so it is not a pass).
  stages_failed = number of stages that did not run to completion. Expected 0.
  stages_completed = number of stages with all repetitions measured. Expected 5.

Usage:
    python -m obdeect_validation.adapters.performance_regression baseline --output FILE
    python -m obdeect_validation.adapters.performance_regression compare [--baseline FILE]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common, performance
from obdeect_validation.adapters import obdeect_cli as cli

REPS = 5
FLOOR = 0.05
SPREAD_K = 3.0
CHUNK_ROWS = 65536  # rows per photon chunk of the Python Writer (its limit is 2**16)
MODEL_VERSION = "7.0.0"
MODEL_ZENITH_DEG = 20.0
MODEL_SEED = 11
MODEL_NAMES = ("LSTN-01", "MSTN-01", "SSTS-01")
POINTING_ZENITH_DEG = 20.0
POINTING_AZIMUTH_DEG = 0.0
OFF_AXIS_DEG = 2.0
SPHERE_RADIUS_M = 35.0  # encloses every compiled model (the check of trace-stream)
SKY_HEIGHT_M = 1000.0  # sky-frame sources above the observation level
SKY_DISC_M = 15.0  # radius of the sky-frame beam on the telescope
OBSERVATION_LEVEL_M = 2000.0  # site of the fixture atmosphere (H2 = 2000 m)
SPHERE_PER_TELESCOPE = 666_666
GROUND_ROWS = 10_000_000
PROPAGATION_ROWS = 5_000_000
SCATTERING_ROWS = 200_000
FLASHER_PHOTONS = 2_000_000_000
FLASHER_CAMERA_RADIUS_M = 1.2
SEED_PROPAGATION = 1
SEED_TRACE = 3
SCATTERING_WAVELENGTHS_NM = (
    345.0,
    395.0,
)  # inside scattering_exponential.dat (340-400 nm)
WAVELENGTHS_NM = (300.0, 600.0)
TRANSMISSION = "fixtures/atm_trans_2000_1_10_0_0_2000.dat"
PROFILE = "fixtures/atmprof_ecmwf_north_winter_fixed.dat"
SCATTERING = "fixtures/scattering_exponential.dat"
BASELINE_DIR = "fixtures/performance"
STAGES = (
    "propagation",
    "propagation_scattering",
    "trace_stream",
    "stream",
    "flasher_dish",
)


@dataclass(frozen=True)
class Child:
    """One timed command of a stage: argv, its log file and, optionally, the stdout file."""

    argv: tuple[str, ...]
    log: str
    stdout: str | None = None


@dataclass(frozen=True)
class Stage:
    """A stage: its children, and the photon count they process (fixed or read back)."""

    name: str
    children: tuple[Child, ...]
    photons: int


# --- host -------------------------------------------------------------------------------------


def cpu_brand() -> str:
    """CPU brand string: sysctl on macOS, /proc/cpuinfo 'model name' on Linux, else the processor."""
    if sys.platform == "darwin":
        try:
            out = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                capture_output=True,
                text=True,
                check=True,
            ).stdout.strip()
            if out:
                return out
        except OSError, subprocess.CalledProcessError:
            pass
    try:
        match = re.search(
            r"^model name\s*:\s*(.+)$",
            Path("/proc/cpuinfo").read_text(),
            re.MULTILINE,
        )
        if match:
            return match[1].strip()
    except OSError:
        pass
    return platform.processor() or platform.machine()


def host_key(system: str, machine: str, brand: str) -> str:
    """Stable file-name key of a host: platform, machine and CPU brand, sanitised."""
    text = "-".join(part for part in (system, machine, brand) if part)
    key = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")
    if not key:
        raise ValueError("host key is empty: platform, machine and CPU brand unknown")
    return key


def current_host_key() -> str:
    return host_key(platform.system(), platform.machine(), cpu_brand())


def default_baseline(key: str) -> Path:
    return common.VALIDATION / BASELINE_DIR / f"{key}.json"


# --- statistics -------------------------------------------------------------------------------


def relative_spread(values: list[float]) -> float:
    """(max - min) / min over the repetitions; 0 for a single value."""
    low = min(values)
    return (max(values) - low) / low if low > 0 else math.inf


def per_photon(cpu_s: list[float], photons: int) -> list[float]:
    """CPU nanoseconds per photon of each repetition."""
    if photons <= 0:
        raise ValueError("photon count must be positive")
    return [value * 1e9 / photons for value in cpu_s]


def summarise(cpu_s: list[float], photons: int) -> dict:
    """Best (minimum) time per photon, the relative spread and every repetition."""
    values = per_photon(cpu_s, photons)
    return {
        "photons": photons,
        "cpu_s": cpu_s,
        "ns_per_photon": values,
        "best_ns_per_photon": min(values),
        "spread": relative_spread(values),
    }


def classify(new: dict, baseline: dict) -> dict:
    """Ratio of the best times, the tolerance and whether the stage regressed (see module docstring)."""
    ratio = new["best_ns_per_photon"] / baseline["best_ns_per_photon"]
    spread = max(new["spread"], baseline["spread"])
    threshold = 1.0 + max(FLOOR, SPREAD_K * spread)
    return {
        "ratio": ratio,
        "threshold": threshold,
        "spread_used": spread,
        "regressed": int(ratio > threshold),
    }


def compare(measured: dict, baseline_stages: dict) -> tuple[dict, dict]:
    """Metrics and per-stage details of the measured stages against the baseline stages."""
    metrics = {"regressions": 0, "baseline_missing": 0}
    details = {}
    for name, new in measured.items():
        entry = {"new": new}
        base = baseline_stages.get(name)
        if base is None:
            metrics["baseline_missing"] += 1
            entry["baseline"] = None
            entry["regressed"] = None
        else:
            entry["baseline"] = base
            result = classify(new, base)
            entry.update(result)
            metrics["regressions"] += result["regressed"]
            metrics[f"{name}_ratio"] = result["ratio"]
        metrics[f"{name}_ns_per_photon"] = new["best_ns_per_photon"]
        metrics[f"{name}_spread"] = new["spread"]
        details[name] = entry
    return metrics, details


# --- inputs (untimed) -------------------------------------------------------------------------


def _chunks(total: int):
    first = 0
    while first < total:
        rows = min(CHUNK_ROWS, total - first)
        yield first, rows
        first += rows


def _axes(zenith_deg: float, azimuth_deg: float = 0.0):
    """Unit vector towards the source (s) and two unit vectors spanning the plane perpendicular to it."""
    th, ph = math.radians(zenith_deg), math.radians(azimuth_deg)
    s = np.array(
        [math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)]
    )
    e1 = np.array(
        [math.cos(th) * math.cos(ph), math.cos(th) * math.sin(ph), -math.sin(th)]
    )
    e2 = np.cross(s, e1)
    return s, e1, e2


def _disc(rng, rows: int, radius_m: float):
    """Radii and azimuths of points uniform in a disc of radius radius_m."""
    r = radius_m * np.sqrt(rng.uniform(0.0, 1.0, rows))
    phi = rng.uniform(0.0, 2.0 * math.pi, rows)
    return r, phi


def _float32(columns: dict) -> dict:
    return {k: np.asarray(v, dtype=np.float32) for k, v in columns.items()}


def _run_record(extra: dict, radius_m: float) -> dict:
    return {
        "run_number": 1,
        "observation_level_m": OBSERVATION_LEVEL_M,
        "telescope_x_m": [0.0],
        "telescope_y_m": [0.0],
        "telescope_z_m": [0.0],
        "telescope_r_m": [radius_m],
        "wavelength_min_nm": WAVELENGTHS_NM[0],
        "wavelength_max_nm": WAVELENGTHS_NM[1],
        **extra,
    }


def write_sky(path, rows: int, seed: int, wavelengths_nm) -> Path:
    """Sky-frame beam: photons above the telescope at SKY_HEIGHT_M, moving straight down, uniform
    over a disc of SKY_DISC_M. Refraction not yet applied (refraction_applied = 0)."""
    from obdeect.stream import Writer

    rng = np.random.default_rng(seed)
    s, e1, e2 = _axes(0.0)
    d = -s
    run = _run_record({"refraction_applied": 0}, SKY_DISC_M)
    with Writer(str(path), "obdeect validation: performance regression") as writer:
        writer.write("run", run)
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1, "frame": 0, "array_index": 0})
        for first, count in _chunks(rows):
            r, phi = _disc(rng, count, SKY_DISC_M)
            q = r[:, None] * (np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2)
            pos = q + SKY_HEIGHT_M * s
            columns = _float32(
                {
                    "x": pos[:, 0],
                    "y": pos[:, 1],
                    "z": pos[:, 2],
                    "cx": np.full(count, d[0]),
                    "cy": np.full(count, d[1]),
                    "cz": np.full(count, d[2]),
                    "t": np.zeros(count),
                    "w": np.ones(count),
                    "wavelength": rng.uniform(*wavelengths_nm, count),
                }
            )
            writer.write_photons(1, "sky", first, columns)
        writer.write("telescope_end", {"telescope_id": 1, "photons": float(rows)})
        writer.write("event_end", {})
        writer.write("run_end", {})
    return Path(path)


def _sphere_entry(rng, rows: int, zenith_deg: float, radius_m: float):
    """Entry points of photons moving along -s through a disc of radius radius_m perpendicular to s.

    The line through a point q of that disc (q . s = 0) enters the sphere of radius SPHERE_RADIUS_M at
    p = q + sqrt(R^2 - |q|^2) s, and moves along d = -s (downward for zenith < 90 deg)."""
    s, e1, e2 = _axes(zenith_deg)
    r, phi = _disc(rng, rows, radius_m)
    q = r[:, None] * (np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2)
    p = q + np.sqrt(SPHERE_RADIUS_M**2 - np.sum(q * q, axis=1))[:, None] * s
    return p, np.broadcast_to(-s, p.shape)


def write_sphere(path, per_telescope: int, dish_radii_m: dict, seed: int) -> Path:
    """Ground-frame sphere rows (ground_representation = sphere), telescopes 1..n, each with half of its
    photons on the axis and half at OFF_AXIS_DEG, over the dish radius of that telescope."""
    from obdeect.stream import Writer

    ids = sorted(dish_radii_m)
    rng = np.random.default_rng(seed)
    run = {
        "run_number": 1,
        "observation_level_m": OBSERVATION_LEVEL_M,
        "ground_representation": "sphere",
        "telescope_x_m": [0.0] * len(ids),
        "telescope_y_m": [0.0] * len(ids),
        "telescope_z_m": [0.0] * len(ids),
        "telescope_r_m": [SPHERE_RADIUS_M] * len(ids),
        "wavelength_min_nm": WAVELENGTHS_NM[0],
        "wavelength_max_nm": WAVELENGTHS_NM[1],
    }
    with Writer(str(path), "obdeect validation: performance regression") as writer:
        writer.write("run", run)
        writer.write("event", {"event_number": 1})
        for index, telescope in enumerate(ids):
            writer.write(
                "telescope",
                {"telescope_id": telescope, "frame": 1, "array_index": index},
            )
            first = 0
            half = per_telescope // 2
            for zenith, rows in (
                (POINTING_ZENITH_DEG, half),
                (POINTING_ZENITH_DEG + OFF_AXIS_DEG, per_telescope - half),
            ):
                for _, count in _chunks(rows):
                    p, d = _sphere_entry(rng, count, zenith, dish_radii_m[telescope])
                    columns = _float32(
                        {
                            "x": p[:, 0],
                            "y": p[:, 1],
                            "z": p[:, 2],
                            "cx": d[:, 0],
                            "cy": d[:, 1],
                            "cz": d[:, 2],
                            "t": np.zeros(count),
                            "zem": np.full(count, OBSERVATION_LEVEL_M + SKY_HEIGHT_M),
                            "w": np.ones(count),
                            "wavelength": rng.uniform(*WAVELENGTHS_NM, count),
                        }
                    )
                    writer.write_photons(telescope, "ground", first, columns)
                    first += count
            writer.write(
                "telescope_end", {"telescope_id": telescope, "photons": float(first)}
            )
        writer.write("event_end", {})
        writer.write("run_end", {})
    return Path(path)


def write_ground(path, rows: int, seed: int) -> Path:
    """Ground-frame rows on the plane z = 0, moving straight down, uniform over SKY_DISC_M (telescope 1)."""
    from obdeect.stream import Writer

    rng = np.random.default_rng(seed)
    run = _run_record({}, SKY_DISC_M)
    with Writer(str(path), "obdeect validation: performance regression") as writer:
        writer.write("run", run)
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1, "frame": 1, "array_index": 0})
        for first, count in _chunks(rows):
            r, phi = _disc(rng, count, SKY_DISC_M)
            columns = _float32(
                {
                    "x": r * np.cos(phi),
                    "y": r * np.sin(phi),
                    "cx": np.zeros(count),
                    "cy": np.zeros(count),
                    "cz": np.full(count, -1.0),
                    "t": np.zeros(count),
                    "zem": np.full(count, OBSERVATION_LEVEL_M + SKY_HEIGHT_M),
                    "w": np.ones(count),
                    "wavelength": rng.uniform(*WAVELENGTHS_NM, count),
                }
            )
            writer.write_photons(1, "ground", first, columns)
        writer.write("telescope_end", {"telescope_id": 1, "photons": float(rows)})
        writer.write("event_end", {})
        writer.write("run_end", {})
    return Path(path)


def dish_radius_m(document: dict) -> float:
    """Radius of the primary of a compiled optical model: facet centres plus circumradius, or the
    outer radius of an asphere."""
    primary = document["optics"]["primary"]
    if "facets" in primary:
        radius = 0.0
        for facet in primary["facets"]:
            cx, cy = float(facet["centre_m"][0]), float(facet["centre_m"][1])
            diameter = float(facet["diameter_m"])
            circum = (
                diameter / math.sqrt(3.0)
                if str(facet["shape"]).startswith("hexagon")
                else 0.5 * diameter
            )
            radius = max(radius, math.hypot(cx, cy) + circum)
        return radius
    if "asphere" in primary:
        return float(primary["asphere"]["outer_radius_m"])
    raise ValueError("optical model primary has neither facets nor an asphere")


# --- stages -----------------------------------------------------------------------------------


def _fixture(name: str) -> str:
    path = common.VALIDATION / name
    if not path.is_file():
        raise RuntimeError(f"fixture missing: {path}")
    return str(path)


def _compile_models(work: Path) -> dict[str, Path]:
    """Compile the optical models untimed; returns the model file of each name and its dish radius."""
    models = {}
    for name in MODEL_NAMES:
        output = work / f"model-{name}.json"
        common.command(
            cli.compile_model(
                str(output),
                name,
                MODEL_VERSION,
                seed=str(MODEL_SEED),
                zenith_deg=str(MODEL_ZENITH_DEG),
                simulation_models=str(common.models_root()),
            ),
            work / f"compile-{name}.log",
        )
        models[name] = output
    return models


def flasher_argv(output: Path) -> tuple[str, ...]:
    """sources run flasher: a dish flasher in front of the camera, a disc target of FLASHER_CAMERA_RADIUS_M."""
    return tuple(
        cli.obdeect(
            "sources",
            "run",
            "flasher",
            type="dish",
            observation_level=str(int(OBSERVATION_LEVEL_M)),
            seed="1",
            atmosphere_id="1",
            position="0,0,16",
            camera_radius=str(FLASHER_CAMERA_RADIUS_M),
            spectrum="fixed:wavelength_nm=400",
            angular="isotropic",
            events="1",
            photons=str(FLASHER_PHOTONS),
            bunch_size="1",
            extra=(
                "--set",
                "run.frame=sky",
                "--set",
                "events[0].emissions[0].to_plane_z0=false",
                "-o",
                str(output),
            ),
        )
    )


def prepare(work: Path) -> dict[str, Stage]:
    """Write the inputs and compile the models (untimed), and return the stages."""
    work.mkdir(parents=True, exist_ok=True)
    transmission = _fixture(TRANSMISSION)
    profile = _fixture(PROFILE)
    scattering = _fixture(SCATTERING)
    models = _compile_models(work)
    dish = {}
    for index, name in enumerate(MODEL_NAMES, start=1):
        dish[index] = dish_radius_m(
            json.loads((work / f"model-{name}.json").read_text())
        )
    sky_wide = write_sky(
        work / "sky-propagation.pst", PROPAGATION_ROWS, 101, WAVELENGTHS_NM
    )
    sky_scatter = write_sky(
        work / "sky-scattering.pst", SCATTERING_ROWS, 102, SCATTERING_WAVELENGTHS_NM
    )
    sphere = write_sphere(work / "sphere.pst", SPHERE_PER_TELESCOPE, dish, 103)
    ground = write_ground(work / "ground.pst", GROUND_ROWS, 104)
    sphere_rows = performance.stream_weight(sphere)[0]
    ground_rows = performance.stream_weight(ground)[0]
    sky_rows = performance.stream_weight(sky_wide)[0]
    scatter_rows = performance.stream_weight(sky_scatter)[0]
    if sphere_rows != SPHERE_PER_TELESCOPE * len(dish):
        raise RuntimeError(f"sphere input has {sphere_rows} rows")
    if ground_rows != GROUND_ROWS:
        raise RuntimeError(f"ground input has {ground_rows} rows")
    if sky_rows != PROPAGATION_ROWS or scatter_rows != SCATTERING_ROWS:
        raise RuntimeError("sky input has an unexpected number of rows")
    flasher_out = work / "flasher-untimed.pst"
    common.command(flasher_argv(flasher_out), work / "flasher-untimed.log")
    flasher_rows = performance.stream_weight(flasher_out)[0]

    model_args = []
    for index, name in enumerate(MODEL_NAMES, start=1):
        model_args += ["--model", f"{index}={models[name]}"]
    telescopes = ",".join(str(i) for i in dish)
    propagation_out = work / "propagated.pst"
    scattered_out = work / "scattered.pst"
    camera_out = work / "camera.pst"
    simtel_out = work / "simtel.iact"
    cat_out = work / "cat.pst"
    return {
        "propagation": Stage(
            "propagation",
            (
                Child(
                    tuple(
                        cli.propagation_run(
                            str(sky_wide),
                            str(propagation_out),
                            transmission=transmission,
                            pathlen="0",
                            profile=profile,
                            site_altitude=str(int(OBSERVATION_LEVEL_M)),
                            seed=str(SEED_PROPAGATION),
                            wavelength_range=(
                                str(WAVELENGTHS_NM[0]),
                                str(WAVELENGTHS_NM[1]),
                            ),
                            threads="1",
                        )
                    ),
                    "propagation.log",
                ),
            ),
            sky_rows,
        ),
        "propagation_scattering": Stage(
            "propagation_scattering",
            (
                Child(
                    tuple(
                        cli.propagation_run(
                            str(sky_scatter),
                            str(scattered_out),
                            scattering=scattering,
                            transmission=transmission,
                            pathlen="0",
                            site_altitude=str(int(OBSERVATION_LEVEL_M)),
                            seed=str(SEED_PROPAGATION),
                            threads="1",
                        )
                    ),
                    "scattering.log",
                ),
            ),
            scatter_rows,
        ),
        "trace_stream": Stage(
            "trace_stream",
            (
                Child(
                    tuple(
                        cli.obdeect(
                            "ray-tracing",
                            "trace-stream",
                            *model_args,
                            seed=str(SEED_TRACE),
                            threads="1",
                            telescopes=telescopes,
                            pointing_zenith_deg=str(POINTING_ZENITH_DEG),
                            pointing_azimuth_deg=str(POINTING_AZIMUTH_DEG),
                            telescope_transmission="full",
                            i=str(sphere),
                            o=str(camera_out),
                        )
                    ),
                    "trace-stream.log",
                ),
            ),
            sphere_rows,
        ),
        "stream": Stage(
            "stream",
            (
                Child(
                    tuple(cli.obdeect("stream", "cat", str(ground))),
                    "cat.log",
                    stdout=str(cat_out),
                ),
                Child(
                    tuple(
                        cli.obdeect(
                            "stream",
                            "to-simtel",
                            o=str(simtel_out),
                            extra=(str(ground),),
                        )
                    ),
                    "to-simtel.log",
                ),
            ),
            ground_rows,
        ),
        "flasher_dish": Stage(
            "flasher_dish",
            (
                Child(
                    flasher_argv(work / "flasher-timed.pst"),
                    "flasher.log",
                ),
            ),
            flasher_rows,
        ),
    }


def run_stage(stage: Stage, work: Path) -> float:
    """CPU seconds of one pass of the stage (sum of its children; wall is not used)."""
    cpu = 0.0
    for child in stage.children:
        record = performance.run_timed(
            list(child.argv),
            work / child.log,
            stdout_file=(work / child.stdout) if child.stdout else None,
        )
        cpu += record["cpu_s"]
    return cpu


def measure(work: Path, reps: int) -> tuple[dict, dict]:
    """Time every stage reps times, interleaved (repetition 1 of every stage, then 2, ...).

    Returns (measured stages, details). A stage that fails is recorded in details['failures'] and
    left out of the measured stages."""
    stages = prepare(work)
    samples: dict[str, list[float]] = {name: [] for name in stages}
    failures: dict[str, str] = {}
    loads = {"before": performance.loadavg()}
    for _ in range(reps):
        for name, stage in stages.items():
            if name in failures:
                continue
            try:
                samples[name].append(run_stage(stage, work))
            except (RuntimeError, OSError) as error:
                failures[name] = str(error)[:400]
    loads["after"] = performance.loadavg()
    measured = {
        name: summarise(samples[name], stages[name].photons)
        for name in stages
        if name not in failures and len(samples[name]) == reps
    }
    details = {
        "stages": {
            name: {
                "photons": stage.photons,
                "commands": [list(c.argv) for c in stage.children],
            }
            for name, stage in stages.items()
        },
        "failures": failures,
        "loadavg": loads,
    }
    return measured, details


def run_host() -> dict:
    """Host record of the baseline and the case."""
    info = performance.host_info()
    info["cpu_brand"] = cpu_brand()
    info["host_key"] = current_host_key()
    return info


def _work_dir() -> Path:
    configured = os.environ.get("PERF_WORK_DIR")
    if configured:
        return Path(configured)
    return Path(tempfile.mkdtemp(prefix="perf-regression-"))


def _run(reps: int) -> tuple[dict, dict, Path]:
    work = _work_dir()
    measured, details = measure(work, reps)
    return measured, details, work


def cmd_baseline(args) -> int:
    work = _work_dir()
    measured, details = measure(work, args.reps)
    if details["failures"] or len(measured) != len(STAGES):
        raise SystemExit(
            f"baseline not written: {details['failures'] or 'stages missing'}"
        )
    host = run_host()
    record = {
        "schema": 1,
        "host_key": host["host_key"],
        "host": host,
        "reps": args.reps,
        "stages": measured,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    common.write_json(output, record)
    shutil.rmtree(work, ignore_errors=True)
    print(f"baseline {host['host_key']}: {output}")
    return 0


def cmd_compare(args) -> int:
    key = current_host_key()
    path = Path(args.baseline) if args.baseline else default_baseline(key)
    measured, details, work = _run(args.reps)
    baseline_stages: dict = {}
    if path.is_file():
        record = json.loads(path.read_text())
        if record.get("host_key") != key:
            raise SystemExit(
                f"baseline {path} is for host {record.get('host_key')}, not {key}"
            )
        baseline_stages = record.get("stages", {})
    metrics, per_stage = compare(measured, baseline_stages)
    metrics["stages_completed"] = len(measured)
    metrics["stages_failed"] = len(details["failures"])
    details.update(
        {
            "case": "performance.regression",
            "host": run_host(),
            "baseline": str(path),
            "baseline_found": path.is_file(),
            "settings": {"reps": args.reps, "floor": FLOOR, "spread_k": SPREAD_K},
            "stages_compared": per_stage,
        }
    )
    common.write_results(metrics, details, Path.cwd())
    if not details["failures"]:
        shutil.rmtree(work, ignore_errors=True)
    for name, entry in per_stage.items():
        ratio = entry.get("ratio")
        text = (
            "no baseline"
            if ratio is None
            else f"ratio {ratio:.3f} limit {entry['threshold']:.3f}"
        )
        print(f"{name}: {entry['new']['best_ns_per_photon']:.1f} ns/photon, {text}")
    print(
        f"regressions {metrics['regressions']}, baseline_missing {metrics['baseline_missing']}"
    )
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="mode", required=True)
    base = sub.add_parser("baseline", help="write the baseline of this host")
    base.add_argument("--output", required=True, type=Path)
    base.add_argument("--reps", type=int, default=REPS)
    comp = sub.add_parser("compare", help="compare with the baseline of this host")
    comp.add_argument("--baseline", type=Path, default=None)
    comp.add_argument("--reps", type=int, default=REPS)
    args = parser.parse_args(argv)
    try:
        if args.mode == "baseline":
            return cmd_baseline(args)
        return cmd_compare(args)
    except (RuntimeError, OSError, ValueError, KeyError) as error:
        print(f"performance.regression: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
