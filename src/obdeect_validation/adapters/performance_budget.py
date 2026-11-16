"""Performance budgets: gates on the timing metrics of the performance cases.

Every gate is an integer metric '<name>_budget_exceeded': 0 = within budget (PASS), 1 = exceeded
(FAIL). A budget whose inputs were not measured is exceeded (1) and its reason is recorded in
details.json: an unmeasured quantity cannot pass. Peak RSS is reported, not gated.

Limits (fixed here, from the plan's "Performance budget per shower"):
  * from-corsika + propagation CPU <= 10 % of the reference sim_telarray CPU time (per shower).
  * stream pipe throughput >= 1 GB/s (decimal, 1e9 bytes/s) for to-simtel and for a cat | tee pipe.
  * whole chain (envelope candidate) wall and CPU <= reference, at equal photoelectron statistics
    (the physics check of the air-shower case: photoelectron bias within max(1 %, 5 sigma)).
    The model compile (compile_model stage, repeated every repetition) is excluded; it is reported.
  * light-source chains (flasher-dish, point-source): candidate wall <= reference wall.
"""

from __future__ import annotations

import re

CPU_FRACTION_LIMIT = 0.10
PIPE_LIMIT_GB_PER_S = 1.0
GIGABYTE = 1e9

# Pipe members of the envelope candidate (performance_chain.envelope_pipeline names).
# sim_telarray (Bypass_Optics=3) line that says the candidate is camera-only: the optics are applied
# by obdeect, so the timed chain is the production chain (obdeect optics + sim_telarray camera).
CAMERA_ONLY_LINE = re.compile(r"Bypass_Optics=3 \(camera only\) for telescope")
PIPE_MEMBERS = ("from_corsika", "propagation")
THREAD_LABELS = ("candidate_envelope_threads1", "candidate_envelope")


def _measured(metrics: dict, keys: list[str]) -> bool:
    return all(isinstance(metrics.get(key), (int, float)) for key in keys)


def _gate(out: dict, reasons: dict, name: str, exceeded: bool) -> None:
    """Integer flag '<name>_budget_exceeded'; a budget with a recorded reason is always exceeded."""
    out[f"{name}_budget_exceeded"] = int(bool(exceeded) or name in reasons)


def _cpu_gate(
    out: dict,
    reasons: dict,
    name: str,
    cpu_keys: list[str],
    ref_key: str,
    metrics: dict,
) -> None:
    """CPU of cpu_keys summed, as a fraction of ref_key; exceeded above CPU_FRACTION_LIMIT."""
    if not _measured(metrics, [*cpu_keys, ref_key]) or metrics[ref_key] <= 0:
        reasons[name] = f"not measured: {', '.join([*cpu_keys, ref_key])}"
        _gate(out, reasons, name, True)
        return
    fraction = sum(metrics[key] for key in cpu_keys) / metrics[ref_key]
    out[f"{name}_fraction"] = fraction
    _gate(out, reasons, name, fraction > CPU_FRACTION_LIMIT)


def shower_budgets(metrics: dict) -> tuple[dict, dict]:
    """Budget metrics of the air-shower case, and the reasons of any budget not measured.

    Gates (names as in metrics.json):
      from_corsika_propagation_cpu_budget_exceeded           separate stages, threads=1
      envelope_from_corsika_propagation_cpu_budget_exceeded  pipe members of the envelope candidate
      candidate_envelope[_threads1]_{wall,cpu}_budget_exceeded  whole chain vs reference
      equal_statistics_budget_exceeded                       photoelectron bias of the envelope chain
    """
    out: dict = {}
    reasons: dict = {}
    ref_cpu = "reference_sim_telarray_cpu_s"
    _cpu_gate(
        out,
        reasons,
        "from_corsika_propagation_cpu",
        [
            "candidate_threads1_from_corsika_cpu_s",
            "candidate_threads1_propagation_cpu_s",
        ],
        ref_cpu,
        metrics,
    )
    _cpu_gate(
        out,
        reasons,
        "envelope_from_corsika_propagation_cpu",
        [f"candidate_envelope_{name}_cpu_s" for name in PIPE_MEMBERS],
        ref_cpu,
        metrics,
    )
    for label in THREAD_LABELS:
        for key, short in (("wall_s", "wall"), ("cpu_s", "cpu")):
            ref, cand = f"reference_{key}", f"{label}_{key}"
            compile_key = f"{label}_compile_model_{key}"
            budget = f"{label}_{short}"
            if _measured(metrics, [ref, cand, compile_key]) and metrics[ref] > 0:
                # The model compile runs once per repetition, not per shower: it is reported as
                # its own stage (<label>_compile_model_*) and left out of the chain comparison.
                chain = metrics[cand] - metrics[compile_key]
                out[f"{budget}_ratio"] = chain / metrics[ref]
                _gate(out, reasons, budget, chain > metrics[ref])
            else:
                reasons[budget] = f"not measured: {ref}, {cand}, {compile_key}"
                _gate(out, reasons, budget, True)
    stats = "physics_envelope_photoelectron_bias_limit_exceeded"
    if _measured(metrics, [stats]):
        _gate(out, reasons, "equal_statistics", metrics[stats] != 0)
    else:
        reasons["equal_statistics"] = f"not measured: {stats}"
        _gate(out, reasons, "equal_statistics", True)
    return out, reasons


def source_budgets(
    metrics: dict, label: str = "candidate_threads1"
) -> tuple[dict, dict]:
    """Light-source chains: candidate wall time <= reference wall time (flasher-dish, point-source)."""
    out: dict = {}
    reasons: dict = {}
    name = f"{label}_wall"
    if _measured(metrics, [f"{label}_wall_s", "reference_wall_s"]):
        out[f"{name}_ratio"] = metrics[f"{label}_wall_s"] / metrics["reference_wall_s"]
        _gate(
            out, reasons, name, metrics[f"{label}_wall_s"] > metrics["reference_wall_s"]
        )
    else:
        reasons[name] = f"not measured: {label}_wall_s, reference_wall_s"
        _gate(out, reasons, name, True)
    return out, reasons


def pipe_budgets(throughputs: dict[str, float | None]) -> tuple[dict, dict]:
    """Pipe throughput [GB/s] per member: exceeded when below PIPE_LIMIT_GB_PER_S (or not measured)."""
    out: dict = {}
    reasons: dict = {}
    for name, value in throughputs.items():
        budget = f"{name}_pipe_throughput"
        if value is None:
            reasons[budget] = "not measured"
            _gate(out, reasons, budget, True)
            continue
        out[f"{name}_pipe_gb_per_s_limit"] = PIPE_LIMIT_GB_PER_S
        _gate(out, reasons, budget, value < PIPE_LIMIT_GB_PER_S)
    return out, reasons


def throughput_gb_per_s(size_bytes: float, wall_s: float) -> float | None:
    """Decimal GB per second of size_bytes moved in wall_s, or None for a non-positive time."""
    if wall_s is None or wall_s <= 0:
        return None
    return size_bytes / wall_s / GIGABYTE


def count_camera_only(text: str) -> int:
    """Number of sim_telarray 'Bypass_Optics=3 (camera only)' lines in its output text."""
    return len(CAMERA_ONLY_LINE.findall(text))
