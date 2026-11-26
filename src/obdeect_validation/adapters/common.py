"""Helpers shared by the adapters: tool locations, command execution, sim_telarray, outputs."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import math
import os
import re
import shutil
import subprocess
from contextlib import nullcontext
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
UNIT_TRANSMISSION = VALIDATION / "fixtures/unit-transmission.dat"
#: Bottom (H2, m a.s.l.) of zero-optical-depth tables. Below the lowest emission height of every
#: validation case: 2147.01 m (South CORSIKA showers, site 2147 m) and 2147.5 m (ground sources).
ZERO_DEPTH_BOTTOM_M = 2140.0
_H2_HEADER = re.compile(r"H2=\s*[0-9.]+")


def package_version(distribution: str) -> str:
    """Installed version of a distribution (for example obdeect-propagation), or 'unknown'."""
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def mean_std(values: list[float]) -> tuple[float, float]:
    """Mean and sample standard deviation (0 for fewer than two values)."""
    mean = sum(values) / len(values)
    if len(values) < 2:
        return mean, 0.0
    return mean, math.sqrt(sum((v - mean) ** 2 for v in values) / (len(values) - 1))


def obdeect_supports_upward() -> bool:
    """True if the installed obdeect takes photons in all directions (stream format 1.1 or later).

    Unit tests of the upward-photon cases skip on an older obdeect (obdeect PR 80).
    """
    try:
        from obdeect.stream import format_major, format_minor
    except ImportError:
        return False
    return (format_major, format_minor) >= (1, 1)


def obdeect_has(module: str) -> bool:
    """True if the installed obdeect provides module (for example obdeect.ray_tracing)."""
    try:
        importlib.import_module(module)
    except ImportError:
        return False
    return True


def obdeect_supports_sphere() -> bool:
    """True if the installed obdeect has the sphere representation (stream format 1.2 or later).

    Unit tests of the sphere cases skip on an older obdeect (obdeect plan, phases 5 and 6).
    """
    try:
        from obdeect.stream import format_major, format_minor
    except ImportError:
        return False
    return (format_major, format_minor) >= (1, 2)


def zero_depth_table(text: str, bottom_m: float = ZERO_DEPTH_BOTTOM_M) -> str:
    """Sim_telarray transmission table with zero optical depth and its bottom (H2) at bottom_m.

    sim_telarray returns zero transmission for emission below H2 (atmospheric_transmission2,
    zem < h2), so a table whose bottom lies above a source's emission heights removes those photons.
    The header H2 is replaced; wavelength rows keep their wavelength and get zero in every column.
    """
    if not _H2_HEADER.search(text):
        raise ValueError("transmission table has no H2 header")
    lines = []
    for line in text.splitlines():
        fields = line.split()
        if fields and fields[0].isdigit():
            lines.append(" ".join([fields[0]] + ["0"] * (len(fields) - 1)))
        else:
            lines.append(_H2_HEADER.sub(f"H2= {bottom_m / 1000.0:.3f}", line, count=1))
    return "\n".join(lines) + "\n"


def write_zero_depth_table(
    root: Path,
    source: Path = UNIT_TRANSMISSION,
    bottom_m: float = ZERO_DEPTH_BOTTOM_M,
    name: str = "unit-atmosphere.dat",
) -> Path:
    """Write zero_depth_table(source) to root/name and return the path (for -C atmospheric_transmission)."""
    path = Path(root) / name
    path.write_text(
        zero_depth_table(Path(source).read_text(encoding="utf-8"), bottom_m),
        encoding="utf-8",
    )
    return path


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


def command(
    argv, log, cwd=None, stdout_file=None, *, check=True, stdout=None, stdin_file=None
) -> int:
    """Run argv; stdout and stderr go to log (stdout to stdout_file, or to stdout, e.g.
    subprocess.DEVNULL, when given); stdin_file is read as standard input. Returns the exit status.

    check: raise RuntimeError with the log tail when the command fails. Every adapter runs its
    commands through this function, endtoend_common.pipeline (pipes, timeouts) or
    performance.run_timed (CPU time and peak RSS).
    """
    argv = [str(arg) for arg in argv]
    with (
        Path(log).open("wb") as log_stream,
        Path(stdin_file).open("rb") if stdin_file else nullcontext() as stdin,
    ):
        if stdout_file:
            with Path(stdout_file).open("wb") as output:
                process = subprocess.run(
                    argv,
                    cwd=cwd,
                    stdin=stdin,
                    stdout=output,
                    stderr=log_stream,
                    check=False,
                )
        elif stdout is not None:
            process = subprocess.run(
                argv,
                cwd=cwd,
                stdin=stdin,
                stdout=stdout,
                stderr=log_stream,
                check=False,
            )
        else:
            process = subprocess.run(
                argv,
                cwd=cwd,
                stdin=stdin,
                stdout=log_stream,
                stderr=subprocess.STDOUT,
                check=False,
            )
    if process.returncode and check:
        tail = " | ".join(
            Path(log).read_text(errors="replace").strip().splitlines()[-3:]
        )
        raise RuntimeError(
            f"{argv[0]} exited {process.returncode}; see {log}"
            + (f": {tail}" if tail else "")
        )
    return process.returncode


def run_tool(argv, log) -> str:
    """Run an executable in the current directory and return its log text."""
    command(argv, log)
    return Path(log).read_text(encoding="utf-8", errors="replace")


# sim_telarray log lines that mean input was dropped although the run exits 0.
SIMTEL_DROPPED_INPUT = (
    "I/O block is skipped",
    "Cannot extend I/O buffer",
    "Error reading",
)


def simtel_dropped_input(log_text: str) -> list[str]:
    """Lines of a sim_telarray log that report skipped or unreadable input blocks."""
    return [
        line.strip()
        for line in log_text.splitlines()
        if any(marker in line for marker in SIMTEL_DROPPED_INPUT)
    ]


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
    min_photons = 0 unless an adapter sets it: sim_telarray otherwise skips a telescope event
    with fewer photons than the configuration's min_photons (120 for LSTN-01), so the two sides
    of a comparison must use the same setting to simulate the same event set.
    """
    output = root / f"{name}.simtel"
    options = {
        "histogram_file": root / f"{name}.hdata",
        "output_file": output,
        "random_seed": SIMTEL_SEED,
        "save_photons": 2,
        "store_photoelectrons": 0,
        "min_photons": 0,
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
    log = root / f"{name}.log"
    runner(argv, log, cwd=REFERENCE)
    dropped = (
        simtel_dropped_input(log.read_text(errors="replace")) if log.is_file() else []
    )
    if dropped:
        raise RuntimeError(f"sim_telarray dropped input ({dropped[0]}); see {log}")
    if not output.is_file():
        raise RuntimeError(f"no sim_telarray output: {log}")
    return output


def write_json(path, value) -> None:
    """Strict JSON (NaN and infinity are errors, not invalid output)."""
    Path(path).write_text(
        json.dumps(value, indent=2, allow_nan=False, default=_plain) + "\n",
        encoding="utf-8",
    )


def results_main(analyse, description: str, add_arguments=None, argv=None) -> None:
    """The command line of an adapter: analyse the current directory and write its results there.

    analyse(root) returns (metrics, details); both are written to root (metrics.json, details.json).
    With add_arguments(parser), the adapter takes options and analyse is called as analyse(root, args).
    """
    parser = argparse.ArgumentParser(description=description)
    if add_arguments is not None:
        add_arguments(parser)
    args = parser.parse_args(argv)
    root = Path.cwd()
    metrics, details = analyse(root) if add_arguments is None else analyse(root, args)
    write_results(metrics, details, root)


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
