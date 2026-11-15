"""Helpers shared by the adapters: tool locations, command execution, sim_telarray, outputs."""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

REFERENCE = Path(os.environ.get("SIMTEL_REFERENCE_SOURCE", "/reference"))
SIMTEL_BINARY = os.environ.get(
    "SIMTEL_BINARY", "/workdir/simulation_software/sim_telarray/bin/sim_telarray"
)
SIMTEL_SEED = 20261113
VALIDATION = Path(os.environ.get("VALIDATION_ROOT", "/validation"))
SHOWERS = VALIDATION / "inputs/gamma-20deg-100-3000GeV.iact"
TRANSMISSION_TABLE = REFERENCE / "cfg/common/atm_trans_2156_1_3_2_0_0_0.1_0.1.dat"


def reference_config_dir() -> Path:
    """Directory with the sim_telarray atmosphere and configuration files."""
    return Path(os.environ.get("SIMTEL_CONFIG_PATH", REFERENCE / "cfg/common"))


def models_root() -> Path:
    """simulation-models checkout (the directory that holds model_parameters/).

    VALIDATION_MODELS overrides; the container mounts the checkout at /models, either
    as the repository (with a simulation-models/ subdirectory) or as the checkout itself.
    """
    configured = os.environ.get("VALIDATION_MODELS")
    if configured:
        return Path(configured)
    mounted = Path("/models")
    if (mounted / "model_parameters").is_dir():
        return mounted
    return mounted / "simulation-models"


def sha256(path) -> str:
    """SHA-256 of a file, read in blocks (CORSIKA data files are large)."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"{name} not found on PATH")
    return path


def command(argv, log, cwd=None, stdout_file=None) -> None:
    """Run argv; stdout and stderr go to log (stdout to stdout_file when given).

    Raises RuntimeError with the log tail when the command fails.
    """
    argv = [str(arg) for arg in argv]
    with Path(log).open("wb") as log_stream:
        if stdout_file:
            with Path(stdout_file).open("wb") as output:
                process = subprocess.run(
                    argv, cwd=cwd, stdout=output, stderr=log_stream, check=False
                )
        else:
            process = subprocess.run(
                argv,
                cwd=cwd,
                stdout=log_stream,
                stderr=subprocess.STDOUT,
                check=False,
            )
    if process.returncode:
        tail = " | ".join(
            Path(log).read_text(errors="replace").strip().splitlines()[-3:]
        )
        raise RuntimeError(
            f"{argv[0]} exited {process.returncode}; see {log}"
            + (f": {tail}" if tail else "")
        )


def run_tool(argv, log) -> str:
    """Run an executable in the current directory and return its log text."""
    command(argv, log)
    return Path(log).read_text(encoding="utf-8", errors="replace")


def sim_telarray_argv(
    root: Path,
    input_file,
    name: str,
    options: dict[str, Any],
    config: str | Path | None = None,
    defines=("-DCTA_NORTH", "-DLST1"),
    extra=(),
    binary: str | None = None,
):
    """argv and output path of sim_telarray on input_file (input_file "-" reads standard input).

    options become -C key=value; extra holds further sim_telarray arguments.
    """
    output = root / f"{name}.simtel"
    options = {
        "histogram_file": root / f"{name}.hdata",
        "output_file": output,
        "random_seed": SIMTEL_SEED,
        "save_photons": 2,
        "store_photoelectrons": 0,
        **options,
    }
    argv = [
        binary or SIMTEL_BINARY,
        "-c",
        config or REFERENCE / "cfg/CTA/CTA-PROD6-LST.cfg",
        f"-I{REFERENCE / 'cfg/common'}",
        f"-I{REFERENCE / 'cfg/CTA'}",
        *defines,
        "-W",
        "random_state=none",
    ]
    for key, value in options.items():
        argv += ["-C", f"{key}={value}"]
    argv += [*extra, input_file]
    return [str(arg) for arg in argv], output


def sim_telarray(
    root: Path,
    input_file,
    name: str,
    options: dict[str, Any],
    config: str | Path | None = None,
    defines=("-DCTA_NORTH", "-DLST1"),
    extra=(),
    runner=command,
    binary: str | None = None,
) -> Path:
    """Run the released sim_telarray on input_file and return <root>/<name>.simtel.

    The command line is kept in <name>.command.json.
    """
    output = root / f"{name}.simtel"
    output.unlink(missing_ok=True)
    (root / f"{name}.hdata").unlink(missing_ok=True)
    argv, output = sim_telarray_argv(
        root, input_file, name, options, config, defines, extra, binary
    )
    (root / f"{name}.command.json").write_text(json.dumps(argv, indent=2) + "\n")
    runner(argv, root / f"{name}.log", cwd=REFERENCE)
    if not output.is_file():
        raise RuntimeError(f"no sim_telarray output: {root / name}.log")
    return output


def write_json(path, value) -> None:
    """Strict JSON (NaN and infinity are errors, not invalid output)."""
    Path(path).write_text(
        json.dumps(value, indent=2, allow_nan=False, default=_plain) + "\n",
        encoding="utf-8",
    )


def write_results(metrics: dict, details: dict | None = None, root: Path | None = None):
    """Write metrics.json and details.json into root (default: the current directory)."""
    root = root or Path.cwd()
    write_json(root / "metrics.json", metrics)
    if details is not None:
        write_json(root / "details.json", details)


def _plain(value):
    """JSON form of numpy scalars and arrays and of paths."""
    if hasattr(value, "tolist"):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"not JSON serialisable: {type(value).__name__}")


def event_deltas(left, right) -> list[float]:
    """Per-event |left - right| photon totals; an event present on one side only counts in full."""
    return [
        abs(left[i] - right[i])
        if i < len(left) and i < len(right)
        else (left[i] if i < len(left) else right[i])
        for i in range(max(len(left), len(right)))
    ]


ATMOSPHERE_READ = re.compile(r"Read atmospheric transmission data from file (\S+)")


def check_atmosphere(text, expected):
    """Raise unless the sim_telarray output read exactly the atmospheric table `expected` (by file name).

    sim_telarray loads its atmosphere from the configuration only under some #if branches, and falls
    back to a built-in table silently otherwise; the read line is the evidence.
    """
    files = [Path(match[1]).name for match in ATMOSPHERE_READ.finditer(text)]
    if not files:
        raise RuntimeError(
            "sim_telarray output has no 'Read atmospheric transmission' line"
        )
    wrong = sorted({name for name in files if name != Path(expected).name})
    if wrong:
        raise RuntimeError(
            f"sim_telarray read atmospheric table(s) {wrong}, expected {Path(expected).name}"
        )
    return files


def check_atmosphere_log(path, expected):
    """check_atmosphere on the text of a log file."""
    return check_atmosphere(Path(path).read_text(errors="replace"), expected)
