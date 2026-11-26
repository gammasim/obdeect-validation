"""Golden outputs: every output byte of a fixed set of obdeect pipelines must equal a frozen baseline.

Guards refactors and performance changes in obdeect that must not change any output bit. The
baseline is produced by the frozen baseline obdeect (main at c125628) with the `baseline` command and
stored per platform in fixtures/golden/<platform-key>.json, platform key f"{sys.platform}-{machine}"
(libm differs between macOS and Linux, so outputs are not comparable across platforms).

Pipelines (threads 1, fixed seeds, all inputs staged from fixtures/ or the sim_telarray tables
$SIMTEL_REFERENCE_SOURCE/cfg/common; the model compile reads $VALIDATION_MODELS):
  compile      optical models LSTN-01, MSTN-01, SSTS-01 (compiled JSON)
  lst-chain    point source (sky, 2 events) -> propagation (2156 m transmission table, refraction
               with a profile, dispersion, single scattering) -> trace-stream LSTN-01 -> to-simtel
  mirror-models the same propagated light through MSTN-01 (single mirror) and SSTS-01 (dual mirror)
               -> trace-stream -> to-simtel
  ground-sphere sphere representation (frame ground): emission inside the telescope sphere -> trace-stream
               LSTN-01; emission outside, sphere-entry rows written (not traced, see below)
  curvature-dispersion propagation with spherical curvature and dispersion; flat propagation with
               dispersion -> propagation to-sky
  corsika      stream from-corsika on a small CORSIKA IACT file (plane, 1205 blocks, 100 GeV, 2200 m)
               -> propagation at 2200 m -> to-grisu (to-grisu needs a plane run with energy_gev)
  events-1-2   point source, two events (event numbers 1 and 2), no extinction -> trace-stream LSTN-01
Every output file of every stage is recorded: sha256, size, and for photon streams the row count and
the sum of weights (obdeect.stream Reader; float64 sums, math.fsum over the photon records).

Limits (exact, no tolerance; a refactor must not change any output bit):
  - outputs_mismatch = 0: the sha256 of every output equals the baseline. Justification: the
    AGENTS.md rule "results must not depend on ..." extends to the bit level for refactors and
    performance changes; any change is either intended (then regenerate the baseline with a reason)
    or a regression.
  - outputs_missing = 0: every baseline output is produced and no output is unexpected (the
    outputs of a failed pipeline count as missing). Exact: a count of files.
  - pipeline_failures = 0: every stage exits 0. Exact: a count of stages that fail.
  - baseline_missing = 0: a platform without a baseline fails the gate (no silent pass). Exact.
The metrics of a run without baseline: outputs_mismatch and outputs_missing are 0 (nothing compared),
baseline_missing is 1, so the case fails.

Known limitation: the sphere-entry rows of the outside emission are not traced, because
`ray-tracing trace-stream` rejects rows that lie about 4 mm outside the telescope sphere ("sphere row
lies outside the telescope sphere", observed with radius 30 m). The file itself is kept and compared.

Inputs: fixtures/golden/gamma_100gev_1216.eventio (copied from pyeventio tests/resources, 35 KB,
sha256 287be666...0efd0a0e, checked here). Derived inputs are staged into the work directory and
recorded with their sha256 in the baseline: the atmospheric transmission tables (cut below 1000 nm, as
the other atmosphere cases do; 2156 m and 2200 m), the atmosphere profile and the scattering table
from fixtures/.

Usage (in the case output directory; writes metrics.json and details.json there):
  python -m obdeect_validation.adapters.golden_outputs run [--baseline FILE] [--workdir DIR]
  python -m obdeect_validation.adapters.golden_outputs baseline --output fixtures/golden/<key>.json
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from obdeect_validation.adapters import atmosphere_chain as atm
from obdeect_validation.adapters import common
from obdeect_validation.adapters import obdeect_cli as cli

EVENTIO_SHA256 = "287be6668e14c9da931fabf89b8573e2902cba1561b1e34ee1a6acab1efd0a0e"
GOLDEN_DIR = common.VALIDATION / "fixtures/golden"
FIXTURES = common.VALIDATION / "fixtures"
TRANSMISSION_2200 = "cfg/common/atm_trans_2200_1_10_0_0_2200.dat"
TRANSMISSION_2156 = "cfg/common/atm_trans_2156_1_3_2_0_0_0.1_0.1.dat"
MODEL_VERSION = "7.0.0"


@dataclass(frozen=True)
class Stage:
    """One obdeect command: argv with placeholders, the output file it writes and, for
    commands that write to standard output, the file that receives it."""

    argv: tuple
    output: str | None = None
    stdout: str | None = None


@dataclass(frozen=True)
class Pipeline:
    name: str
    stages: tuple
    depends_on: tuple = ()


def _trace(model, source, output):
    return Stage(
        tuple(
            cli.trace_stream(
                f"{model}.json",
                source,
                output,
                seed="7",
                threads="1",
                pointing_zenith_deg="20",
                pointing_azimuth_deg="180",
            )
        ),
        output,
    )


def _to_simtel(source, camera_map, output):
    return Stage(
        tuple(
            cli.obdeect(
                "stream",
                "to-simtel",
                camera_map=camera_map,
                atmprof="atmprof.dat",
                o=output,
                extra=[source],
            )
        ),
        output,
    )


def _point_source(
    position, photons, seed, output, telescope_radius, events=2, frame=None
):
    argv = cli.obdeect(
        "sources",
        "run",
        "point-source",
        observation_level="2156",
        spectrum="fixed:wavelength_nm=400",
        position=position,
        photons=photons,
        events=str(events),
        telescope_radius=telescope_radius,
        seed=str(seed),
        set=f"run.frame={frame}" if frame else None,
        o=output,
    )
    return Stage(tuple(argv), output)


def _propagation(source, output, *, site, transmission, seed, extra=()):
    return Stage(
        tuple(
            cli.obdeect(
                "propagation",
                "run",
                transmission=transmission,
                pathlen="0",
                site_altitude=site,
                seed=str(seed),
                threads="1",
                extra=[*extra, "-i", source, "-o", output],
            )
        ),
        output,
    )


PROFILE = ("--profile", "atmprof.dat", "--dispersion-reference-nm", "500")

PIPELINES = (
    Pipeline(
        "compile",
        tuple(
            Stage(
                tuple(
                    cli.compile_model(
                        f"{model}.json",
                        model,
                        MODEL_VERSION,
                        seed="11",
                        zenith_deg="20",
                        simulation_models="{models}",
                    )
                ),
                f"{model}.json",
            )
            for model in ("LSTN-01", "MSTN-01", "SSTS-01")
        ),
    ),
    Pipeline(
        "lst-chain",
        (
            _point_source("0,0,2000", "1e9", 1, "point.pst", "12.5"),
            _propagation(
                "point.pst",
                "prop.pst",
                site="2156",
                transmission="transmission-2156.dat",
                seed=3,
                extra=(*PROFILE, "--scattering", "scattering.dat"),
            ),
            _trace("LSTN-01", "prop.pst", "lst-trace.pst"),
            _to_simtel("lst-trace.pst", "cam-single.txt", "lst.simtel"),
        ),
        depends_on=("compile",),
    ),
    Pipeline(
        "mirror-models",
        (
            _trace("MSTN-01", "prop.pst", "mstn-trace.pst"),
            _to_simtel("mstn-trace.pst", "cam-single.txt", "mstn.simtel"),
            _trace("SSTS-01", "prop.pst", "ssts-trace.pst"),
            _to_simtel("ssts-trace.pst", "cam-dual.txt", "ssts.simtel"),
        ),
        depends_on=("compile", "lst-chain"),
    ),
    Pipeline(
        "ground-sphere",
        (
            _point_source(
                "0,0,10", "1e5", 2, "ground-inside.pst", "30", frame="ground"
            ),
            _trace("LSTN-01", "ground-inside.pst", "ground-trace.pst"),
            _point_source(
                "0,0,2000", "1e8", 4, "ground-sphere.pst", "30", frame="ground"
            ),
        ),
        depends_on=("compile",),
    ),
    Pipeline(
        "curvature-dispersion",
        (
            _propagation(
                "point.pst",
                "curved.pst",
                site="2156",
                transmission="transmission-2156.dat",
                seed=3,
                extra=(*PROFILE, "--curvature", "spherical"),
            ),
            _propagation(
                "point.pst",
                "flat-dispersion.pst",
                site="2156",
                transmission="transmission-2156.dat",
                seed=3,
                extra=PROFILE,
            ),
            Stage(
                tuple(
                    cli.obdeect(
                        "propagation",
                        "to-sky",
                        profile="atmprof.dat",
                        i="flat-dispersion.pst",
                        o="sky-from-ground.pst",
                    )
                ),
                "sky-from-ground.pst",
            ),
        ),
        depends_on=("lst-chain",),
    ),
    Pipeline(
        "corsika",
        (
            Stage(
                tuple(cli.from_corsika("corsika.eventio")),
                "corsika.pst",
                stdout="corsika.pst",
            ),
            _propagation(
                "corsika.pst",
                "corsika-prop.pst",
                site="2200",
                transmission="transmission-2200.dat",
                seed=3,
            ),
            Stage(
                tuple(
                    cli.obdeect(
                        "stream",
                        "to-grisu",
                        o="corsika.grisu",
                        extra=["corsika-prop.pst"],
                    )
                ),
                "corsika.grisu",
            ),
        ),
    ),
    Pipeline(
        "events-1-2",
        (
            _point_source("0,0,2000", "1e8", 9, "events12.pst", "12.5"),
            Stage(
                tuple(
                    cli.propagation_run(
                        "events12.pst",
                        "events12-prop.pst",
                        no_extinction=True,
                        site_altitude="2156",
                        seed="9",
                        threads="1",
                    )
                ),
                "events12-prop.pst",
            ),
            _trace("LSTN-01", "events12-prop.pst", "events12-trace.pst"),
        ),
    ),
)


def platform_key() -> str:
    return f"{sys.platform}-{platform.machine()}"


def _expand(argv, models):
    return [str(arg).replace("{models}", str(models)) for arg in argv]


def _cut(source: Path, destination: Path) -> None:
    destination.write_text(
        atm.cut_table(Path(source).read_text(encoding="utf-8"), atm.MAX_LAMBDA_NM),
        encoding="utf-8",
    )


def stage_inputs(workdir: Path) -> dict:
    """Copy the fixtures and derive the transmission tables into workdir; return their provenance."""
    reference = common.REFERENCE
    staged = {}

    def copy(source: Path, name: str, label: str):
        shutil.copyfile(source, workdir / name)
        staged[name] = {"source": label, "sha256": common.sha256(source)}

    eventio = GOLDEN_DIR / "gamma_100gev_1216.eventio"
    if common.sha256(eventio) != EVENTIO_SHA256:
        raise RuntimeError(f"{eventio}: sha256 differs from the fixture record")
    copy(eventio, "corsika.eventio", "fixtures/golden/gamma_100gev_1216.eventio")
    copy(
        FIXTURES / "atmprof_ecmwf_north_winter_fixed.dat",
        "atmprof.dat",
        "fixtures/atmprof_ecmwf_north_winter_fixed.dat",
    )
    copy(
        FIXTURES / "scattering_exponential.dat",
        "scattering.dat",
        "fixtures/scattering_exponential.dat",
    )
    for name, table in (
        ("transmission-2156.dat", TRANSMISSION_2156),
        ("transmission-2200.dat", TRANSMISSION_2200),
    ):
        source = reference / table
        _cut(source, workdir / name)
        staged[name] = {
            "source": f"$SIMTEL_REFERENCE_SOURCE/{table}, cut below {atm.MAX_LAMBDA_NM:g} nm",
            "sha256": common.sha256(source),
        }
    for name, mirror_class in (("cam-single.txt", 1), ("cam-dual.txt", 2)):
        (workdir / name).write_text(f"1 {mirror_class} 0.0\n", encoding="utf-8")
        staged[name] = {
            "source": "written by golden_outputs",
            "sha256": common.sha256(workdir / name),
        }
    return staged


def stream_stats(path: Path) -> tuple[int, float]:
    """Photon rows and sum of weights of an obdeect photon stream (obdeect.stream Reader)."""
    import numpy as np
    from obdeect.stream import Reader

    rows, sums = 0, []
    for record in Reader(str(path)):
        if record["type"] == "photons":
            weights = np.asarray(dict(record["columns"])["w"], dtype=np.float64)
            rows += int(weights.size)
            sums.append(float(np.sum(weights, dtype=np.float64)))
    return rows, math.fsum(sums)


def describe_output(path: Path, stage_argv) -> dict:
    entry = {
        "sha256": common.sha256(path),
        "bytes": path.stat().st_size,
        "rows": None,
        "sum_weights": None,
        "command": stage_argv[0:4],
    }
    if path.suffix == ".pst":
        entry["rows"], entry["sum_weights"] = stream_stats(path)
    return entry


def obdeect_provenance() -> dict:
    version = (
        subprocess.run(
            [common.tool("obdeect"), "--version"],
            capture_output=True,
            text=True,
            check=True,
        )
        .stdout.strip()
        .splitlines()
    )
    checkout = Path(
        os.environ.get("OBDEECT_CHECKOUT", common.VALIDATION.parent / "obdeect")
    )
    try:
        commit = subprocess.run(
            ["git", "-C", str(checkout), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except OSError, subprocess.CalledProcessError:
        commit = "unknown"
    return {"obdeect_version": version, "obdeect_commit": commit}


def run_pipelines(workdir: Path, models) -> dict:
    """Run every pipeline in workdir. Returns outputs, pipeline status and the stage argvs."""
    logs = workdir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    outputs, status, failed = {}, {}, set()
    for pipeline in PIPELINES:
        started = time.monotonic()
        error = None
        if set(pipeline.depends_on) & failed:
            error = f"dependency failed: {', '.join(sorted(set(pipeline.depends_on) & failed))}"
        for index, stage in enumerate(pipeline.stages if error is None else ()):
            argv = _expand(stage.argv, models)
            log = logs / f"{pipeline.name}-{index}.log"
            stdout_file = workdir / stage.stdout if stage.stdout else None
            argv[0] = common.tool("obdeect")
            try:
                common.command(argv, log, cwd=workdir, stdout_file=stdout_file)
                if stage.output:
                    path = workdir / stage.output
                    if not path.is_file():
                        raise RuntimeError(f"{stage.output} was not written")
                    outputs[stage.output] = {
                        **describe_output(path, stage.argv),
                        "pipeline": pipeline.name,
                    }
            except (RuntimeError, OSError, ValueError) as exc:
                error = f"stage {index}: {exc}"
                break
        if error:
            failed.add(pipeline.name)
        status[pipeline.name] = {
            "status": "FAILED" if error else "OK",
            "error": error,
            "seconds": round(time.monotonic() - started, 3),
            "stages": [list(stage.argv) for stage in pipeline.stages],
            "depends_on": list(pipeline.depends_on),
        }
    return {"outputs": outputs, "pipelines": status, "failures": len(failed)}


def compare(baseline_outputs: dict, current_outputs: dict) -> tuple[dict, int, int]:
    """Per-output verdicts; returns (details, mismatches, missing)."""
    details, mismatches, missing = {}, 0, 0
    for name in sorted(set(baseline_outputs) | set(current_outputs)):
        base = baseline_outputs.get(name)
        new = current_outputs.get(name)
        if base is None:
            verdict = "unexpected"
        elif new is None:
            verdict = "missing"
        elif base["sha256"] != new["sha256"]:
            verdict = "mismatch"
        else:
            verdict = "match"
        if verdict in ("unexpected", "missing"):
            missing += 1
        elif verdict == "mismatch":
            mismatches += 1
        details[name] = {
            "verdict": verdict,
            "pipeline": (new or base or {}).get("pipeline"),
            "baseline_sha256": base and base["sha256"],
            "sha256": new and new["sha256"],
            "baseline_rows": base and base["rows"],
            "rows": new and new["rows"],
            "baseline_sum_weights": base and base["sum_weights"],
            "sum_weights": new and new["sum_weights"],
        }
    return details, mismatches, missing


def _workdir(path: str | None, name: str) -> Path:
    workdir = Path(path).resolve() if path else (Path.cwd() / name).resolve()
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    return workdir


def run_command(args) -> int:
    started = time.monotonic()
    key = platform_key()
    baseline_path = Path(args.baseline) if args.baseline else GOLDEN_DIR / f"{key}.json"
    workdir = _workdir(args.workdir, "golden-work")
    inputs = stage_inputs(workdir)
    result = run_pipelines(workdir, common.models_root())
    details = {
        "platform": key,
        "baseline_file": str(baseline_path),
        "current": obdeect_provenance(),
        "inputs": inputs,
        "pipelines": result["pipelines"],
    }
    metrics = {"pipeline_failures": result["failures"]}
    if baseline_path.is_file():
        baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
        comparison, mismatches, missing = compare(
            baseline["outputs"], result["outputs"]
        )
        details["baseline"] = baseline["provenance"]
        details["outputs"] = comparison
        metrics.update(
            outputs_mismatch=mismatches, outputs_missing=missing, baseline_missing=0
        )
    else:
        details["outputs"] = {
            name: {"verdict": "not_compared", **entry}
            for name, entry in result["outputs"].items()
        }
        metrics.update(outputs_mismatch=0, outputs_missing=0, baseline_missing=1)
        details["note"] = f"no baseline for {key}: {baseline_path} does not exist"
    details["runtime_seconds"] = round(time.monotonic() - started, 3)
    Path("details.json").write_text(
        json.dumps(details, indent=2, sort_keys=True) + "\n"
    )
    Path("metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n"
    )
    return 0


def baseline_command(args) -> int:
    started = time.monotonic()
    key = platform_key()
    workdir = _workdir(args.workdir, "golden-baseline-work")
    inputs = stage_inputs(workdir)
    result = run_pipelines(workdir, common.models_root())
    if result["failures"]:
        failed = {
            k: v["error"] for k, v in result["pipelines"].items() if v["status"] != "OK"
        }
        print(f"baseline not written: pipelines failed: {failed}", file=sys.stderr)
        return 2
    data = {
        "platform": key,
        "platform_detail": {
            "sys.platform": sys.platform,
            "machine": platform.machine(),
            "platform": platform.platform(),
            "python": sys.version.split()[0],
        },
        "provenance": {
            **obdeect_provenance(),
            "date": datetime.date.today().isoformat(),
            "command": [
                "python",
                "-m",
                "obdeect_validation.adapters.golden_outputs",
                "baseline",
            ],
            "runtime_seconds": round(time.monotonic() - started, 3),
            "inputs": inputs,
            "pipelines": {
                name: {"depends_on": value["depends_on"], "stages": value["stages"]}
                for name, value in result["pipelines"].items()
            },
            "placeholders": {
                "{models}": "$VALIDATION_MODELS (simulation-models checkout)",
                "work directory": "stage commands run inside a fresh work directory",
            },
        },
        "outputs": result["outputs"],
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"baseline for {key}: {len(result['outputs'])} outputs -> {output}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser(
        "run", help="compare the outputs with the baseline of this platform"
    )
    run.add_argument(
        "--baseline", help="baseline JSON (default fixtures/golden/<platform>.json)"
    )
    run.add_argument("--workdir", help="work directory (default ./golden-work)")
    run.set_defaults(func=run_command)
    base = commands.add_parser(
        "baseline", help="write the baseline JSON of this platform"
    )
    base.add_argument(
        "--output", required=True, help="e.g. fixtures/golden/<platform>.json"
    )
    base.add_argument(
        "--workdir", help="work directory (default ./golden-baseline-work)"
    )
    base.set_defaults(func=baseline_command)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
