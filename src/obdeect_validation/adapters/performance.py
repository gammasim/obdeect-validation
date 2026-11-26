"""Per-stage throughput of the obdeect chain on the fixed CORSIKA IACT input (report only).

Same setup as the air-shower case (performance_chain): LSTN-01 7.0.0 production atmosphere cut
below 1000 nm, QE envelope, one telescope, pointing of the chain config. Stages, in order:
from_corsika, to_simtel (round trip, conservation check), propagation (weighted), propagation_envelope
(thinned by the QE bound), compile_model (reported only), trace_stream (weighted input), simtel_direct
(sim_telarray Bypass_Optics=0 reference). The throughput stages are the timed ones in stages_completed.

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.performance --config cfg/model/7.0.0/CTAO.cfg [--events N] [--reps N]

Timing: each child is timed with performance_chain.run_timed (obdeect subcommands as native programs,
as in the air-shower case). Each stage records the load average before it starts; the host load
before and after the run and the spread of each stage's wall time (timing_unreliable) tell whether
the numbers are usable.
"""

import argparse
import os
import shutil
import statistics
import tempfile
from pathlib import Path

from obdeect_validation.adapters import chain_optics as co
from obdeect_validation.adapters import common
from obdeect_validation.adapters import obdeect_cli as cli
from obdeect_validation.adapters import performance_chain as chain
from obdeect_validation.adapters import performance_common as pcm

REPS = 3
SPREAD_LIMIT = (
    0.05  # (max - min) / median of a stage's wall time above this: timing unreliable
)
LOAD_LIMIT = (
    1.0  # host load average (1 min) before the run above this: timing unreliable
)
CPU_SHARE_STAGES = ("from_corsika", "propagation_envelope", "trace_stream", "to_simtel")


def measure(
    name, run_once, photons_in, input_file, output_file, metrics, details, reps
):
    start = len(pcm.RUSAGE)
    try:
        for _ in range(reps):
            run_once()
    except (
        Exception  # noqa: BLE001
    ) as error:  # a failed stage is recorded; the remaining stages still run
        details["failures"][name] = str(error)
        metrics[f"{name}_completed"] = 0
        return False
    records = pcm.RUSAGE[start:]
    median = {
        key: statistics.median(record[key] for record in records)
        for key in ("wall_s", "cpu_s", "peak_rss_mb")
    }
    wall = [record["wall_s"] for record in records]
    metrics.update(
        {
            f"{name}_completed": 1,
            f"{name}_reps": len(records),
            f"{name}_wall_s_median": median["wall_s"],
            f"{name}_wall_s_min": min(wall),
            f"{name}_wall_s_max": max(wall),
            f"{name}_wall_s_spread": (max(wall) - min(wall)) / median["wall_s"],
            f"{name}_cpu_s_median": median["cpu_s"],
            f"{name}_peak_rss_mb": median["peak_rss_mb"],
            f"{name}_bytes_in": input_file.stat().st_size,
            f"{name}_bytes_out": output_file.stat().st_size,
        }
    )
    if (
        photons_in is not None
    ):  # a stage without a photon count (compile_model) has no rate
        metrics.update(
            {
                f"{name}_photons_in": photons_in,
                f"{name}_photons_per_s": photons_in / median["wall_s"],
            }
        )
    details["stages"][name] = {
        "input": str(input_file),
        "output": str(output_file),
        "records": records,
    }
    return True


def cpu_metrics(metrics, stages):
    """Per stage: CPU share of the envelope chain (CPU_SHARE_STAGES) and CPU ns per photon.

    The share divides by the CPU medians of the CPU_SHARE_STAGES present in metrics.
    ns per photon uses the stage's photons_in (the rate convention of photons_per_s).
    """
    chain = [
        metrics[f"{name}_cpu_s_median"]
        for name in CPU_SHARE_STAGES
        if f"{name}_cpu_s_median" in metrics
    ]
    chain_cpu = sum(chain)
    out = {}
    for name in stages:
        cpu = metrics[f"{name}_cpu_s_median"]
        if chain_cpu > 0:
            out[f"{name}_cpu_share"] = cpu / chain_cpu
        photons = metrics.get(f"{name}_photons_in")
        if photons:
            out[f"{name}_ns_per_photon_cpu"] = cpu * 1e9 / photons
    return out


def timing_unreliable(metrics, stages):
    """1 when a stage's wall-time spread is above SPREAD_LIMIT or the host load before the run is
    above LOAD_LIMIT; a missing load average does not count."""
    spreads = [metrics.get(f"{name}_wall_s_spread") or 0.0 for name in stages]
    load = metrics.get("host_loadavg_before")
    busy = load is not None and load > LOAD_LIMIT
    return int(busy or any(spread > SPREAD_LIMIT for spread in spreads))


def _plot(path, stages, metrics):
    """Throughput per stage (one colour per stage), with the per-stage CPU and wall-time medians."""
    from obdeect_validation.adapters import plotting

    fig, axes = plotting.figure(1, 2, (11, 3.8))
    colors = [plotting.CYCLE[i % len(plotting.CYCLE)] for i in range(len(stages))]
    values = [metrics[f"{name}_photons_per_s"] for name in stages]
    bars = axes[0, 0].barh(stages, values, color=colors, height=0.6)
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_xlabel("photons per second (log scale, median of repetitions)")
    axes[0, 0].invert_yaxis()
    axes[0, 0].grid(axis="y", visible=False)
    for bar, value in zip(bars, values, strict=True):
        axes[0, 0].text(
            value * 1.05,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.3g}",
            va="center",
            fontsize=9,
        )
    axes[0, 0].set_xlim(right=max(values) * 10)
    wall = [metrics[f"{name}_wall_s_median"] for name in stages]
    cpu = [metrics[f"{name}_cpu_s_median"] for name in stages]
    rows = range(len(stages))
    axes[0, 1].barh(
        [r - 0.2 for r in rows], wall, height=0.4, color=colors, label="wall"
    )
    axes[0, 1].barh(
        [r + 0.2 for r in rows], cpu, height=0.4, color=colors, alpha=0.45, label="CPU"
    )
    axes[0, 1].set_yticks(list(rows), stages)
    axes[0, 1].invert_yaxis()
    axes[0, 1].set_xlabel("seconds (median)")
    axes[0, 1].grid(axis="y", visible=False)
    axes[0, 1].legend()
    plotting.save(fig, path)


def main():
    parser = argparse.ArgumentParser()
    pcm.add_arguments(
        parser,
        config_help="sim_telarray CTAO.cfg of the LSTN-01 7.0.0 array (simtools generate-array-config)",
        input_default=pcm.INPUT,
        events_help="measure the first N showers (default all)",
        reps_default=REPS,
    )
    args = parser.parse_args()
    root = Path.cwd().resolve()
    # Stream files go to a container-local directory (PERF_WORK_DIR), not the bind-mounted results.
    work_dir = os.environ.get("PERF_WORK_DIR")
    work = Path(work_dir or tempfile.mkdtemp(prefix="perf-stages-"))
    work.mkdir(parents=True, exist_ok=True)
    metrics = {}
    load_before = pcm.loadavg()
    if load_before is not None:
        metrics["host_loadavg_before"] = load_before
    source = pcm.select_input(
        args.input.resolve(), args.events, work / "showers-input.iact"
    )
    inputs, cfg_text = chain.prepare_inputs(work, args.config)
    headroom = chain.headroom_for(inputs["qe_variation"])
    (work / "envelope").mkdir(parents=True, exist_ok=True)
    envelope = chain.envelope_files(work / "envelope", inputs["qe"], headroom)
    details = {
        "host": pcm.host_info(),
        "stages": {},
        "failures": {},
        "skips": {},
        "settings": {
            "reps": args.reps,
            "input": str(source),
            "showers": pcm.shower_count(source),
            "work_dir": str(work),
            "config": str(inputs["config"]),
            "model": f"{co.MODEL} {co.VERSION}",
            "reference_atmosphere": co.cfg_value(cfg_text, "atmospheric_transmission"),
            "propagation_atmosphere": str(inputs["atmosphere"]),
            "wavelength_limit_nm": 1000,
            "propagation": "weighted 250-700 nm (as shower_candidate); envelope: thinned by the QE bound (as envelope_pipeline)",
            "propagation_seed": 123,
            "envelope_headroom": headroom,
            "envelope_bound_max": envelope["bound_max"],
            "telescopes": list(co.TELESCOPES),
            "pointing": f"zenith {co.ZENITH_DEG} deg, azimuth {co.AZIMUTH_DEG} deg",
            "trace_seed": co.SEED,
            "threads": 1,
            "simtel_direct": "performance_chain.shower_reference: sim_telarray Bypass_Optics=0, one telescope, reference atmosphere",
            "timing": "performance_chain.run_timed: wait4 rusage of each child (obdeect subcommands as native programs); wall time around the child; loadavg_before per record",
            "limits": "conservation: relative delta <= 1e-6 (IACT stores weights as float32; stream sums in float64)",
        },
    }
    iact_in = pcm.iact_weight(source)
    stream = work / "showers.gst"
    propagated = work / "propagated.gst"
    enveloped = work / "propagated-envelope.gst"
    model = work / "lst.model"
    traced = work / "traced.gst"
    simtel_dir = work / "simtel-direct"
    simtel_dir.mkdir(parents=True, exist_ok=True)
    run = chain.run_timed
    done = {}

    done["from_corsika"] = measure(
        "from_corsika",
        lambda: run(
            ["obdeect", "stream", "from-corsika", source],
            work / "from-corsika.log",
            stdout_file=stream,
        ),
        iact_in,
        source,
        stream,
        metrics,
        details,
        args.reps,
    )
    if done["from_corsika"]:
        stream_in = pcm.stream_weight(stream)[1]
        metrics["iact_input_vs_stream_weight_relative_delta"] = (
            abs(stream_in - iact_in) / iact_in
        )
        done["to_simtel"] = measure(
            "to_simtel",
            lambda: run(
                cli.to_simtel(stream, work / "roundtrip.iact"),
                work / "to-simtel.log",
            ),
            stream_in,
            stream,
            work / "roundtrip.iact",
            metrics,
            details,
            args.reps,
        )
    else:
        stream_in = None
        skip = "from_corsika failed; no input stream"
        details["skips"].update(
            to_simtel=skip, propagation=skip, propagation_envelope=skip
        )
    if done.get("to_simtel"):
        roundtrip_in = pcm.iact_weight(work / "roundtrip.iact")
        metrics["from_corsika_to_simtel_weight_relative_delta"] = (
            abs(roundtrip_in - stream_in) / stream_in
        )

    if done["from_corsika"]:
        done["propagation"] = measure(
            "propagation",
            lambda: run(
                cli.propagation_run(
                    stream,
                    propagated,
                    transmission=inputs["atmosphere"],
                    site_altitude=2156,
                    pathlen=28,
                    seed=123,
                    mode="weighted",
                    wavelength_range=(250, 700),
                    telescopes=",".join(map(str, co.TELESCOPES)),
                ),
                work / "propagation.log",
            ),
            stream_in,
            stream,
            propagated,
            metrics,
            details,
            args.reps,
        )
        done["propagation_envelope"] = measure(
            "propagation_envelope",
            lambda: run(
                cli.propagation_run(
                    stream,
                    enveloped,
                    transmission=inputs["atmosphere"],
                    site_altitude=2156,
                    pathlen=28,
                    seed=123,
                    envelope=f"{chain.ENVELOPE_FACTOR}={envelope['envelope']}",
                    telescopes=",".join(map(str, co.TELESCOPES)),
                ),
                work / "propagation-envelope.log",
            ),
            stream_in,
            stream,
            enveloped,
            metrics,
            details,
            args.reps,
        )
        if done["propagation_envelope"]:
            metrics["propagation_envelope_photons_out"] = pcm.stream_weight(enveloped)[
                0
            ]

    done["compile_model"] = measure(
        "compile_model",
        lambda: run(chain.compile_argv(inputs, model), work / "compile.log"),
        None,
        inputs["config"],
        model,
        metrics,
        details,
        args.reps,
    )

    done["trace_stream"] = False
    if done.get("propagation") and done["compile_model"]:
        trace_in = pcm.stream_weight(propagated)[1]
        done["trace_stream"] = measure(
            "trace_stream",
            lambda: run(
                cli.trace_stream(
                    model,
                    propagated,
                    traced,
                    seed=co.SEED,
                    threads=1,
                    telescopes=",".join(map(str, co.TELESCOPES)),
                    pointing_zenith_deg=co.ZENITH_DEG,
                    pointing_azimuth_deg=co.AZIMUTH_DEG,
                    telescope_transmission="full",
                ),
                work / "trace-stream.log",
            ),
            trace_in,
            propagated,
            traced,
            metrics,
            details,
            args.reps,
        )
    else:
        details["skips"].setdefault(
            "trace_stream",
            "propagation or model compile failed or skipped; no input stream",
        )

    _name, reference_run = chain.shower_reference(
        simtel_dir,
        source,
        inputs["config"],
        inputs["config_dir"],
        inputs["reference_atmosphere"],
    )[0]
    done["simtel_direct"] = measure(
        "simtel_direct",
        reference_run,
        iact_in,
        source,
        simtel_dir / "reference.simtel",
        metrics,
        details,
        args.reps,
    )

    # compile_model is reported, not a throughput stage (no photon rate, not in stages_completed).
    completed = [name for name in done if done[name] and name != "compile_model"]
    metrics["stages_failed"] = len(details["failures"])
    metrics["stages_skipped"] = len(details["skips"])
    metrics["stages_completed"] = len(completed)
    metrics.update(cpu_metrics(metrics, completed))
    load_after = pcm.loadavg()
    if load_after is not None:
        metrics["host_loadavg_after"] = load_after
    metrics["timing_unreliable"] = timing_unreliable(metrics, list(details["stages"]))
    details["stage_order"] = [name for name in done if name != "compile_model"]

    if completed:
        _plot(root / "photons-per-second.png", completed, metrics)

    common.write_results(metrics, details, root)

    if work_dir is None:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
