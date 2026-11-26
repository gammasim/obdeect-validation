"""Determinism of the photon stages (propagation and ray tracing), photon by photon.

Limits: none beyond exact equality. Every comparison must give zero mismatched photons (event,
telescope, row index and every column bit for bit). Compared: thread counts 1, 4 and 0 (all cores);
input chunk sizes 1024 and 65536 rows (the output chunk size is fixed by the format, 2^16 rows, and
the CLI has no option for it); the file chain, the pipe chain (from-corsika | propagation | trace-stream
with '-' stages) and the fused in-process Python calls (propagate then trace_stream). Propagation
threads are only in the Python API (propagate(threads=)); the propagation CLI has none.
"""

import os
from pathlib import Path

from obdeect_validation.adapters import atmosphere_chain as atm
from obdeect_validation.adapters import common
from obdeect_validation.adapters import endtoend_common as e2e
from obdeect_validation.adapters import obdeect_cli as cli

MODEL, VERSION, ZENITH_DEG, AZIMUTH_DEG = "LSTN-01", "7.0.0", 20, 180
SEED = 5
PROPAGATION_SEED = 123
MIN_PHOTONS = 1000


def trace_argv(model, threads, stream_in=None, stream_out=None):
    """ray-tracing trace-stream with the chain's options; -i and -o only for file streams."""
    streams = {} if stream_in is None else {"i": stream_in, "o": stream_out}
    return cli.obdeect(
        "ray-tracing",
        "trace-stream",
        model=model,
        seed=SEED,
        threads=threads,
        pointing_zenith_deg=ZENITH_DEG,
        pointing_azimuth_deg=AZIMUTH_DEG,
        **streams,
    )


def compare(metrics, details, name, left, right):
    result = e2e.compare_rows(left, right)
    metrics[f"{name}_photon_mismatch"] = result["photon_mismatches"]
    metrics[f"{name}_column_set_mismatch"] = 1 - result["column_sets_equal"]
    metrics[f"{name}_empty_comparison"] = int(result["photons_compared"] < MIN_PHOTONS)
    details.setdefault("comparisons", {})[name] = result
    return result


def add_input(parser):
    parser.add_argument(
        "--input",
        type=Path,
        default=common.VALIDATION / "inputs/gamma-20deg-20-100GeV.iact",
    )


def analyse(root, args):
    source = args.input.resolve()
    band = atm.cherenkov_band(
        source.with_suffix(".corsika.log").read_text(errors="replace")
    )
    cut = root / "cut-transmission.dat"
    cut.write_text(
        atm.cut_table(common.TRANSMISSION_TABLE.read_text(), atm.MAX_LAMBDA_NM)
    )
    metrics, details = (
        {},
        {
            "inputs": str(source),
            "seed": SEED,
            "propagation_seed": PROPAGATION_SEED,
            "model": f"{MODEL} {VERSION}",
            "band_nm": band,
        },
    )
    stages = {}

    native = root / "native.gst"
    stages["from-corsika"] = e2e.run_stage(
        cli.from_corsika(source),
        log=root / "from-corsika.log",
        output=native,
    )
    model = root / "model.json"
    common.command(
        cli.compile_model(model, MODEL, VERSION, zenith_deg=ZENITH_DEG),
        root / "compile.log",
    )

    # File chain.
    prop = root / "prop-file.gst"
    stages["propagation-file"] = e2e.run_stage(
        cli.weighted_propagation(cut, band, native, prop, seed=PROPAGATION_SEED),
        log=root / "prop-file.log",
    )
    trace = {}
    for threads in (1, 4, 0):
        out = root / f"trace-file-t{threads}.gst"
        stages[f"trace-file-t{threads}"] = e2e.run_stage(
            trace_argv(model, threads, prop, out),
            log=root / f"trace-file-t{threads}.log",
        )
        trace[threads] = out
    compare(
        metrics,
        details,
        "threads_4_vs_1",
        e2e.photon_rows(trace[1]),
        e2e.photon_rows(trace[4]),
    )
    compare(
        metrics,
        details,
        "threads_0_vs_1",
        e2e.photon_rows(trace[1]),
        e2e.photon_rows(trace[0]),
    )

    # Negative control: another run seed must show up as differing photons, else the checks above prove nothing.
    control = root / "trace-control-seed.gst"
    control_argv = trace_argv(model, 4, prop, control)
    control_argv[control_argv.index("--seed") + 1] = str(SEED + 1)
    stages["trace-control"] = e2e.run_stage(
        control_argv, log=root / "trace-control.log"
    )
    details["control_photon_mismatch"] = e2e.compare_rows(
        e2e.photon_rows(trace[1]), e2e.photon_rows(control)
    )["photon_mismatches"]
    metrics["control_seed_change_detected"] = int(
        details["control_photon_mismatch"] > 0
    )

    # Input chunk sizes: the same photons rewritten in chunks of 1024 and 65536 rows.
    for rows in (1024, 65536):
        rechunked = root / f"prop-chunk{rows}.gst"
        e2e.rewrite(prop, rechunked, split_rows=rows)
        out = root / f"trace-chunk{rows}.gst"
        stages[f"trace-chunk{rows}"] = e2e.run_stage(
            trace_argv(model, 4, rechunked, out), log=root / f"trace-chunk{rows}.log"
        )
        compare(
            metrics,
            details,
            f"input_chunk{rows}_vs_original",
            e2e.photon_rows(trace[1]),
            e2e.photon_rows(out),
        )
        details.setdefault("input_chunk_records", {})[rows] = e2e.photon_rows(
            rechunked
        )["_chunks"]
    details["input_chunk_records"]["original"] = e2e.photon_rows(prop)["_chunks"]

    # Pipes: from-corsika | propagation | trace-stream, every stage on stdin/stdout.
    pipe_out = root / "trace-pipe.gst"
    pipe = e2e.pipeline(
        [
            cli.from_corsika(source),
            cli.weighted_propagation(cut, band, "-", "-", seed=PROPAGATION_SEED),
            trace_argv(model, 4),
        ],
        output=pipe_out,
        log_dir=root / "pipe-stages",
        timeout=600,
    )
    metrics["pipe_stage_failures"] = sum(e2e.classify(stage) != "ok" for stage in pipe)
    details["pipe_stage_verdicts"] = [e2e.describe(stage) for stage in pipe]
    compare(
        metrics,
        details,
        "pipe_vs_file_threads_4",
        e2e.photon_rows(trace[4]),
        e2e.photon_rows(pipe_out),
    )

    # Fused in process: propagate and trace_stream from the Python API.
    from obdeect.propagation.propagation import propagate
    from obdeect.ray_tracing.model import OpticalModel
    from obdeect.ray_tracing.stream import trace_stream

    prop_threads = {}
    all_cores = min(
        256, os.cpu_count() or 1
    )  # propagate accepts 1 to 256 threads (0 is an error)
    for threads in (1, 4, all_cores):
        out = root / f"prop-api-t{threads}.gst"
        propagate(
            native,
            out,
            transmission=cut,
            site_altitude_m=atm.SITE_ALTITUDE_M,
            seed=PROPAGATION_SEED,
            mode="weighted",
            wavelength_range=band,
            threads=threads,
            pathlen_m=0.0,  # as the CLI run of the same case (--pathlen 0)
        )
        prop_threads[threads] = out
    compare(
        metrics,
        details,
        "propagation_api_threads_4_vs_file",
        e2e.photon_rows(prop),
        e2e.photon_rows(prop_threads[4]),
    )
    compare(
        metrics,
        details,
        "propagation_api_threads_all_vs_file",
        e2e.photon_rows(prop),
        e2e.photon_rows(prop_threads[all_cores]),
    )
    details["all_cores"] = all_cores
    optical = OpticalModel.load(model)
    for threads in (1, 4, 0):
        out = root / f"trace-api-t{threads}.gst"
        trace_stream(
            optical,
            prop_threads[1],
            out,
            seed=SEED,
            threads=threads,
            telescopes=None,
            pointing=(ZENITH_DEG, AZIMUTH_DEG),
        )
        compare(
            metrics,
            details,
            f"fused_api_threads_{threads}_vs_file",
            e2e.photon_rows(trace[1]),
            e2e.photon_rows(out),
        )
    details["fused_vs_pipe_bytes_equal"] = common.sha256(trace[4]) == common.sha256(
        pipe_out
    )
    details["stage_verdicts"] = {
        name: e2e.describe(stage) for name, stage in stages.items()
    }
    metrics["file_stage_failures"] = sum(
        e2e.classify(stage) != "ok" for stage in stages.values()
    )
    e2e.summary_figure(
        root / "summary.png", metrics, "determinism (threads, chunks, pipes, fused API)"
    )
    return metrics, details


def main():
    common.results_main(analyse, __doc__.split("\n")[0], add_arguments=add_input)


if __name__ == "__main__":
    main()
