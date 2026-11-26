"""Shared primitives of the performance cases: timing of child processes, host description, input
selection, work directories, command-line options."""

import contextlib
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import numpy as np
from eventio import EventIOFile
from eventio.iact import IACTFile

from obdeect_validation.adapters import common
from obdeect_validation.adapters import obdeect_cli as cli

INPUT = common.VALIDATION / "inputs/gamma-20deg-100-3000GeV-100showers.iact"

RUSAGE = []


def loadavg():
    """1-minute load average of the host, or None where it is not available."""
    try:
        return os.getloadavg()[0]
    except OSError, AttributeError:
        return None


PROC_VMHWM = re.compile(r"^VmHWM:\s+(\d+)\s+kB", re.MULTILINE)


def parse_vmhwm_kb(text):
    """Peak resident set [kB] from the text of /proc/<pid>/status, or None when absent."""
    match = PROC_VMHWM.search(text)
    return int(match[1]) if match else None


def read_vmhwm_kb(pid):
    try:
        return parse_vmhwm_kb(Path(f"/proc/{pid}/status").read_text())
    except OSError:  # the child has exited, or there is no /proc
        return None


class RssSampler:
    """Peak RSS of one running child, sampled from /proc until stop() (Linux only).

    wait4 ru_maxrss is not used on Linux: exec keeps the high-water mark of the forked parent,
    so every child of a large Python process reports the parent's peak (identical values for all
    stages). The sampled VmHWM is the child's own peak, a lower bound at the sampling interval.
    """

    def __init__(self, pid, read=None, interval=0.002):
        self.read = read or (lambda: read_vmhwm_kb(pid))
        self.interval = interval
        self.peak_kb = None
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        while True:
            kb = self.read()
            if kb is not None and (self.peak_kb is None or kb > self.peak_kb):
                self.peak_kb = kb
            if self._stop.wait(self.interval):
                return

    def start(self):
        self._thread.start()

    def stop(self):
        """Stop sampling and return the peak [kB], or None when no sample was taken."""
        self._stop.set()
        self._thread.join()
        return self.peak_kb


def run_timed(argv, log, cwd=None, stdout_file=None):
    """Run one child and keep its own wall time, CPU time and peak RSS.

    Peak RSS: sampled VmHWM on Linux, wait4 ru_maxrss elsewhere (see RssSampler). The record
    keeps the raw ru_maxrss as rusage_maxrss_mb for comparison.
    """
    argv = list(map(str, argv))
    with (
        Path(log).open("wb") as log_stream,
        Path(stdout_file).open("wb") if stdout_file else log_stream as output,
    ):
        load_before = loadavg()
        start = time.perf_counter()
        process = subprocess.Popen(argv, cwd=cwd, stdout=output, stderr=log_stream)
        sampler = RssSampler(process.pid) if sys.platform.startswith("linux") else None
        if sampler:
            sampler.start()
        _, status, usage = os.wait4(process.pid, 0)
        wall = time.perf_counter() - start
        sampled_kb = sampler.stop() if sampler else None
    process.returncode = os.waitstatus_to_exitcode(status)
    rss_divisor = 1024 * 1024 if sys.platform == "darwin" else 1024
    rusage_mb = usage.ru_maxrss / rss_divisor
    if sampled_kb is not None:
        peak_mb, source = sampled_kb / 1024, "proc-VmHWM-sampled"
    else:
        peak_mb, source = rusage_mb, "wait4-ru_maxrss"
    record = {
        "argv": argv,
        "wall_s": wall,
        "cpu_s": usage.ru_utime + usage.ru_stime,
        "peak_rss_mb": peak_mb,
        "peak_rss_source": source,
        "rusage_maxrss_mb": rusage_mb,
        "returncode": process.returncode,
        "loadavg_before": load_before,
    }
    RUSAGE.append(record)
    if process.returncode:
        raise RuntimeError(f"{argv[0]} exited {process.returncode}; see {log}")
    return record


def stream_weight(path):
    """Photon count and bunch weight sum of a native photon stream, from obdeect's own summary."""
    text = subprocess.run(
        cli.obdeect("stream", "cat", summary=True, extra=(str(path),)),
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    match = re.search(r"total: (\d+) photons, weight sum ([0-9.eE+-]+)", text)
    if not match:
        raise ValueError(f"no total in stream summary of {path}")
    return int(match[1]), float(match[2])


def iact_weight(path):
    """Sum of CORSIKA IACT bunch photon weights (the 'photons' field)."""
    total = 0.0
    with IACTFile(path) as stream:
        for event in stream:
            for photons in event.photon_bunches.values():
                total += float(np.sum(photons["photons"], dtype=np.float64))
    return total


EVENT_HEADER = (
    1202  # eventio type of a CORSIKA IACT event header; the event's blocks follow it
)
EVENT_TYPES = range(
    1202, 1210
)  # EventHeader .. EventEnd (1203-1209: the event's telescope blocks)
BLOCK_HEADER_BYTES = 16  # eventio block header in front of each top-level block payload


def _top_level_blocks(path):
    """(type, start, end) byte range of each top-level eventio block, header included."""
    with EventIOFile(path) as stream:
        for block in stream:
            start = block.address - BLOCK_HEADER_BYTES
            yield block.header.type, start, block.address + block.size


def shower_count(path):
    """Number of showers (event headers) of a CORSIKA IACT file."""
    return sum(
        1 for block_type, _, _ in _top_level_blocks(path) if block_type == EVENT_HEADER
    )


def subset_showers(source, target, events):
    """Copy the first `events` showers of an IACT file as raw eventio blocks; the run header and
    the blocks after the events (run end) are kept. Returns the number of showers written."""
    if events < 0:
        raise ValueError("number of showers must not be negative")
    kept = seen = 0
    with Path(source).open("rb") as raw, Path(target).open("wb") as out:
        for block_type, start, end in _top_level_blocks(source):
            if block_type == EVENT_HEADER:
                seen += 1
            if block_type not in EVENT_TYPES or seen <= events:
                kept += block_type == EVENT_HEADER
                raw.seek(start)
                out.write(raw.read(end - start))
    return kept


def select_input(source, events, target):
    """The input to measure: source itself (events None), or its first `events` showers in target."""
    if events is None:
        return source
    subset_showers(source, target, events)
    return target


def restrict_table(source, target, limit=1000.0):
    """Keep comment lines and rows whose wavelength [nm] is below limit."""
    kept = []
    for line in source.read_text().splitlines(keepends=True):
        fields = line.split()
        try:
            keep = (
                not fields or line.lstrip().startswith("#") or float(fields[0]) < limit
            )
        except ValueError:
            keep = True
        if keep:
            kept.append(line)
    target.write_text("".join(kept))
    return len(kept)


def cpu_model():
    try:
        text = Path("/proc/cpuinfo").read_text()
    except OSError:
        return platform.processor() or platform.machine()
    for key in ("model name", "Hardware", "cpu model", "CPU part"):
        match = re.search(rf"^{key}\s*:\s*(.+)$", text, re.MULTILINE)
        if match:
            return match[1].strip()
    return platform.machine()


def host_info():
    version = subprocess.run(
        cli.obdeect("--version"), capture_output=True, text=True, check=False
    ).stdout.strip()
    return {
        "nproc": len(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else os.cpu_count(),
        "cpu_model": cpu_model(),
        "machine": platform.machine(),
        "platform": platform.platform(),
        "obdeect_version": version,
        "host": os.environ.get("OBDEECT_VALIDATION_HOST"),
    }


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


@contextlib.contextmanager
def work_dir(prefix, name):
    """Work directory <base>/<name>, base = PERF_WORK_DIR when set, else a new temporary directory.

    The directory is emptied on entry and removed on exit; a temporary base is removed as well.
    """
    configured = os.environ.get("PERF_WORK_DIR")
    base = Path(configured) if configured else Path(tempfile.mkdtemp(prefix=prefix))
    work = base / name
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True, exist_ok=True)
    try:
        yield work
    finally:
        shutil.rmtree(work, ignore_errors=True)
        if not configured:
            shutil.rmtree(base, ignore_errors=True)


def add_arguments(
    parser,
    *,
    config_help=None,
    input_default=INPUT,
    events_default=None,
    events_help=None,
    reps_default,
):
    """Options shared by the performance cases: --config (required), --input, --events, --reps."""
    parser.add_argument("--config", type=Path, required=True, help=config_help)
    parser.add_argument("--input", type=Path, default=input_default)
    parser.add_argument("--events", type=int, default=events_default, help=events_help)
    parser.add_argument("--reps", type=int, default=reps_default)
