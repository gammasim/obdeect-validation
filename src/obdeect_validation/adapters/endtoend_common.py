"""Helpers for the end-to-end stream cases: pipelines with timeouts, photon-by-photon comparison,
stream rewriting for tampered and damaged inputs.

Comparison is exact. A photon is equal to another only if its event, telescope, row index and every
column agree bit for bit (NaN payloads included). Chunk boundaries are not compared: a stream rewritten
with other chunk sizes must give the same rows in the same order.

Classification of a command result (classify):
  hang      the stage did not exit within the timeout (all stages are killed)
  crash     killed by a signal, exit status >= 128, or a Python traceback on stderr
  error     exit status 1 or 2 with a message (a clean, reported error)
  ok        exit status 0
"""

from __future__ import annotations

import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

#: Columns compared bit for bit (names of the photon stream format).
ROW_KEYS = ("event", "telescope", "index")


@dataclass
class Stage:
    argv: list
    returncode: int | None = None
    stderr: str = ""
    timed_out: bool = False
    seconds: float = 0.0
    log: str = ""
    extra: dict = field(default_factory=dict)


def classify(stage: Stage) -> str:
    """hang, crash, error or ok (see the module docstring)."""
    if stage.timed_out:
        return "hang"
    code = stage.returncode
    if (
        code is None
        or code < 0
        or code >= 128
        or "Traceback (most recent call last)" in stage.stderr
    ):
        return "crash"
    if code == 0:
        return "ok"
    if code in (1, 2):
        return "error"
    return "crash"


def pipeline(
    argvs,
    *,
    output=None,
    stdin=None,
    log_dir,
    timeout=120.0,
) -> list[Stage]:
    """Run argvs as a shell pipeline (stage i writes to stage i+1); the last stage writes to output.

    stdin: file for the first stage (default: empty input). output: file for the last stage (default:
    discarded). Each stage's stderr goes to log_dir/stage<i>.log. A stage that does not finish within
    timeout seconds, and the whole pipeline with it, is killed; the stage is then reported as hang.
    The pipe read end is closed in the parent, so a consumer that exits early delivers SIGPIPE upstream.
    """
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    stages, processes, handles = [], [], []
    try:
        first_in = open(stdin, "rb") if stdin else subprocess.DEVNULL  # noqa: SIM115
        if stdin:
            handles.append(first_in)
        last_out = open(output, "wb") if output else subprocess.DEVNULL  # noqa: SIM115
        if output:
            handles.append(last_out)
        previous = None
        for index, argv in enumerate(argvs):
            argv = [str(arg) for arg in argv]
            log = log_dir / f"stage{index}.log"
            log_handle = log.open("wb")
            handles.append(log_handle)
            stage_in = first_in if index == 0 else previous.stdout
            stage_out = last_out if index == len(argvs) - 1 else subprocess.PIPE
            started = time.monotonic()
            process = subprocess.Popen(
                argv, stdin=stage_in, stdout=stage_out, stderr=log_handle
            )
            if index > 0:
                previous.stdout.close()
            processes.append(process)
            stages.append(Stage(argv=argv, log=str(log)))
            stages[-1].extra["started"] = started
            previous = process
        deadline = time.monotonic() + timeout
        hung = False
        for process in processes:
            remaining = max(0.0, deadline - time.monotonic())
            try:
                process.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                hung = True
                break
        running = [process.poll() is None for process in processes]
        if hung:
            for process in processes:
                if process.poll() is None:
                    process.kill()
                process.wait()
        for process, stage, still_running in zip(processes, stages, running):
            stage.returncode = process.returncode
            stage.timed_out = still_running
            stage.seconds = time.monotonic() - stage.extra.pop("started")
            stage.stderr = Path(stage.log).read_text(errors="replace")
        return stages
    finally:
        for handle in handles:
            handle.close()


def run_stage(argv, *, log, output=None, stdin=None, timeout=120.0) -> Stage:
    """pipeline with one stage."""
    log_path = Path(log)
    stages = pipeline(
        [argv],
        output=output,
        stdin=stdin,
        log_dir=log_path.parent / f"{log_path.stem}.stages",
        timeout=timeout,
    )
    stage = stages[0]
    stage.log = str(log_path)
    log_path.write_text(stage.stderr)
    return stage


def describe(stage: Stage) -> str:
    """One line for details.json: verdict, exit status and the last stderr line."""
    lines = [line for line in stage.stderr.strip().splitlines() if line.strip()]
    tail = lines[-1] if lines else ""
    return f"{classify(stage)} rc={stage.returncode} {tail}".strip()


# ---------------------------------------------------------------------------------------------
# Photon-by-photon comparison


def photon_rows(path) -> dict:
    """Every photon of a stream in file order: event (index among event records, from 1), telescope,
    row index (first_index + offset within the telescope and event), and the columns present.

    Chunk boundaries are dropped; the row index is the one of the format (per telescope and event).
    """
    from obdeect.stream import Reader

    event = 0
    parts: dict[str, list] = {key: [] for key in ROW_KEYS}
    columns: dict[str, list] = {}
    names = None
    chunks = 0
    for record in Reader(str(path)):
        kind = record["type"]
        if kind == "event":
            event += 1
        elif kind == "photons" and record["columns"]:
            chunks += 1
            rows = len(record["columns"]["x"])
            parts["event"].append(np.full(rows, event, dtype=np.int64))
            parts["telescope"].append(
                np.full(rows, record["telescope_id"], dtype=np.int64)
            )
            parts["index"].append(
                np.arange(
                    record["first_index"], record["first_index"] + rows, dtype=np.int64
                )
            )
            present = tuple(sorted(record["columns"]))
            if names is None:
                names = present
                columns = {name: [] for name in names}
            if present != names:
                raise ValueError(
                    f"{path}: photon chunks with different columns {present} and {names}"
                )
            for name in names:
                columns[name].append(
                    np.asarray(record["columns"][name], dtype=np.float32)
                )
    result = {key: _concat(values, np.int64) for key, values in parts.items()}
    for name in names or ():
        result[name] = _concat(columns[name], np.float32)
    result["_chunks"] = chunks
    return result


def _concat(values, dtype):
    if not values:
        return np.zeros(0, dtype=dtype)
    return np.concatenate(values)


def compare_rows(left: dict, right: dict) -> dict:
    """Photons that differ between two photon_rows results, by key or bit pattern of any column.

    Returns photons_compared (the larger count), photon_mismatches (rows that differ, plus the rows
    present on one side only) and column_sets_equal.
    """
    names = sorted(name for name in left if not name.startswith("_"))
    other = sorted(name for name in right if not name.startswith("_"))
    n_left, n_right = (
        len(left["x"]) if "x" in left else 0,
        len(right["x"]) if "x" in right else 0,
    )
    common = min(n_left, n_right)
    mismatch = abs(n_left - n_right)
    if common:
        differs = np.zeros(common, dtype=bool)
        for name in ROW_KEYS:
            differs |= left[name][:common] != right[name][:common]
        for name in names:
            if name in ROW_KEYS or name not in right:
                continue
            a = np.asarray(left[name][:common], dtype=np.float32).view(np.uint32)
            b = np.asarray(right[name][:common], dtype=np.float32).view(np.uint32)
            differs |= a != b
        mismatch += int(differs.sum())
    return {
        "photons_compared": int(max(n_left, n_right)),
        "photon_mismatches": int(mismatch),
        "column_sets_equal": int(names == other),
    }


def compare_streams(left, right) -> dict:
    """compare_rows of two stream files."""
    return compare_rows(photon_rows(left), photon_rows(right))


def drop_events(rows: dict, events) -> dict:
    """photon_rows result without the photons of the given event numbers."""
    keep = ~np.isin(rows["event"], list(events))
    return {
        name: values[keep] if not name.startswith("_") else values
        for name, values in rows.items()
    }


# ---------------------------------------------------------------------------------------------
# Stream rewriting


def rewrite(source, destination, edit=None, *, split_rows=None):
    """Copy a stream record by record through the Python Writer.

    edit(record) returns None to keep the record unchanged, or a list of replacement items:
      ("meta", record_type, metadata)            metadata record (run, event, telescope, ...)
      ("photons", telescope_id, frame, first, columns)
      ("blob", eventio_type, eventio_id, data)
      an empty list drops the record.
    split_rows: write photon chunks of at most this many rows (rechunking; rows keep their indices).
    The stream and stream_end records are written by the Writer. A failure while writing aborts the
    output without stream_end (the reader then reports it as truncated).
    """
    from obdeect.stream import Reader, Writer

    with Writer(str(destination), "endtoend-rewrite") as writer:
        for record in Reader(str(source)):
            kind = record["type"]
            if kind in ("stream", "stream_end"):
                continue
            items = edit(record) if edit else None
            if items is None:
                items = _unchanged(record, split_rows)
            for item in items:
                _write_item(writer, item, split_rows)


def _unchanged(record, split_rows):
    kind = record["type"]
    if kind == "photons":
        columns = record["columns"]
        if not columns:
            return [
                (
                    "photons",
                    record["telescope_id"],
                    record["frame"],
                    record["first_index"],
                    {},
                )
            ]
        return [
            (
                "photons",
                record["telescope_id"],
                record["frame"],
                record["first_index"],
                columns,
            )
        ]
    if kind == "blob":
        return [("blob", record["eventio_type"], record["eventio_id"], record["data"])]
    return [("meta", kind, dict(record.get("metadata", {})))]


def _write_item(writer, item, split_rows):
    kind = item[0]
    if kind == "meta":
        writer.write(item[1], item[2])
    elif kind == "blob":
        writer.write_blob(item[1], item[2], item[3])
    elif kind == "photons":
        _, telescope, frame, first, columns = item
        rows = len(next(iter(columns.values()))) if columns else 0
        if rows == 0:
            return
        step = split_rows or rows
        for start in range(0, rows, step):
            writer.write_photons(
                telescope,
                frame,
                first + start,
                {
                    name: np.ascontiguousarray(
                        values[start : start + step], dtype=np.float32
                    )
                    for name, values in columns.items()
                },
            )
    else:
        raise ValueError(f"unknown item {kind}")


def efficiency_records(path) -> list[dict]:
    """Metadata of every efficiency_status record, in file order, with its scope (run or telescope)."""
    from obdeect.stream import Reader

    scope, found = "run", []
    for record in Reader(str(path)):
        if record["type"] == "telescope":
            scope = f"telescope {record['metadata'].get('telescope_id')}"
        elif record["type"] == "event":
            scope = "run"
        elif record["type"] == "efficiency_status":
            found.append({"scope": scope, "metadata": dict(record["metadata"])})
    return found


def with_status_factor(metadata: dict, name: str, state: str, source: str) -> dict:
    """Copy of an efficiency status with factor name set to state (factors list extended)."""
    updated = dict(metadata)
    factors = [item for item in str(updated.get("factors", "")).split(",") if item]
    if name not in factors:
        factors.append(name)
    updated["factors"] = ",".join(factors)
    updated[f"{name}.state"] = state
    updated[f"{name}.source"] = source
    return updated


def without_factor(metadata: dict, name: str) -> dict:
    """Copy of an efficiency status without the factor name."""
    updated = {k: v for k, v in metadata.items() if not k.startswith(f"{name}.")}
    factors = [
        item
        for item in str(updated.get("factors", "")).split(",")
        if item and item != name
    ]
    updated["factors"] = ",".join(factors)
    return updated


def write_text_list(path, lines) -> Path:
    path = Path(path)
    path.write_text("\n".join(lines) + "\n")
    return path


def fan_out(source, branches, *, log_dir, timeout=300.0):
    """Copy source with obdeect stream tee to one FIFO per branch while every branch consumes its FIFO.

    branches: list of (name, argv_of, stdout): argv_of(fifo) returns the consumer command line (its input
    is the FIFO path), and
    stdout is the file for its standard output (None discards it). The consumers run concurrently; the
    tee's own standard output goes to log_dir/tee-stdout.gst. Returns (tee Stage, [branch Stage]). A
    branch still running at the timeout is killed and reported as hang; the tee is then also reported.
    """
    import os
    import tempfile

    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    # FIFOs live in a local temporary directory: bind-mounted host filesystems may refuse them.
    fifo_dir = Path(tempfile.mkdtemp(prefix="fifo-"))
    fifos = [fifo_dir / f"{name}.fifo" for name, _, _ in branches]
    running, handles = [], []
    try:
        for fifo in fifos:
            os.mkfifo(fifo)
        # Consumers first: a FIFO opened for reading waits for the writer, which tee opens below.
        for (name, argv_of, stdout), fifo in zip(branches, fifos):
            argv = [str(arg) for arg in argv_of(fifo)]
            log = (log_dir / f"{name}.log").open("wb")
            out = open(stdout, "wb") if stdout else subprocess.DEVNULL  # noqa: SIM115
            handles += [log] + ([out] if stdout else [])
            process = subprocess.Popen(argv, stdout=out, stderr=log)
            running.append((name, argv, process, time.monotonic()))
        tee = pipeline(
            [["obdeect", "stream", "tee", *fifos]],
            stdin=source,
            output=log_dir / "tee-stdout.gst",
            log_dir=log_dir / "tee",
            timeout=timeout,
        )[0]
        deadline = time.monotonic() + timeout
        stages = []
        for name, argv, process, started in running:
            timed_out = False
            try:
                process.wait(timeout=max(0.0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                timed_out = True
                process.kill()
                process.wait()
            stage = Stage(
                argv=argv,
                returncode=process.returncode,
                timed_out=timed_out,
                seconds=time.monotonic() - started,
                log=str(log_dir / f"{name}.log"),
            )
            stage.stderr = Path(stage.log).read_text(errors="replace")
            stages.append(stage)
        return tee, stages
    finally:
        for handle in handles:
            handle.close()
        for fifo in fifos:
            fifo.unlink(missing_ok=True)
        fifo_dir.rmdir()


def summary_figure(path, metrics: dict, title: str) -> None:
    """Horizontal bars of the metric values (zero: no difference, no failure), one figure per case."""
    from obdeect_validation.adapters import plotting

    names = list(metrics)
    fig, axis = plotting.figure(size=(9, 0.32 * len(names) + 1.2))
    axis.barh(
        range(len(names)),
        [float(metrics[name]) for name in names],
        color=plotting.case_color(title),
    )
    axis.set_yticks(range(len(names)))
    axis.set_yticklabels(names, fontsize=7)
    axis.invert_yaxis()
    axis.set_xlabel("value (expected: see the case file)")
    axis.set_title(title, fontsize=9)
    plotting.save(fig, path)
