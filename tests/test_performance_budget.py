"""Budget gates of the performance cases (performance_budget) and the shower subset helper."""

import json
from pathlib import Path

import pytest

from obdeect_validation.adapters import performance, performance_budget


def shower_metrics(**overrides):
    """A passing set of air-shower metrics (medians, seconds); overrides change single keys."""
    metrics = {
        "reference_wall_s": 100.0,
        "reference_cpu_s": 100.0,
        "reference_sim_telarray_cpu_s": 100.0,
        "candidate_threads1_from_corsika_cpu_s": 2.0,
        "candidate_threads1_propagation_cpu_s": 3.0,
        "candidate_envelope_wall_s": 60.0,
        "candidate_envelope_cpu_s": 70.0,
        "candidate_envelope_from_corsika_cpu_s": 1.0,
        "candidate_envelope_propagation_cpu_s": 2.0,
        "candidate_envelope_threads1_wall_s": 80.0,
        "candidate_envelope_threads1_cpu_s": 90.0,
        "candidate_envelope_compile_model_wall_s": 0.5,
        "candidate_envelope_compile_model_cpu_s": 0.4,
        "candidate_envelope_threads1_compile_model_wall_s": 0.5,
        "candidate_envelope_threads1_compile_model_cpu_s": 0.4,
        "physics_envelope_photoelectron_bias_limit_exceeded": 0,
    }
    metrics.update(overrides)
    return metrics


def test_passing_shower_budgets_are_zero():
    budgets, reasons = performance_budget.shower_budgets(shower_metrics())
    assert reasons == {}
    flags = {k: v for k, v in budgets.items() if k.endswith("_budget_exceeded")}
    assert flags and set(flags.values()) == {0}
    assert budgets["from_corsika_propagation_cpu_fraction"] == pytest.approx(0.05)


def test_cpu_budget_is_exceeded_above_ten_percent():
    budgets, _ = performance_budget.shower_budgets(
        shower_metrics(candidate_threads1_propagation_cpu_s=9.0)
    )
    assert budgets["from_corsika_propagation_cpu_budget_exceeded"] == 1
    assert budgets["envelope_from_corsika_propagation_cpu_budget_exceeded"] == 0


def test_chain_wall_and_cpu_gates_compare_with_reference():
    budgets, _ = performance_budget.shower_budgets(
        shower_metrics(
            candidate_envelope_threads1_cpu_s=101.0, candidate_envelope_cpu_s=99.0
        )
    )
    assert budgets["candidate_envelope_threads1_cpu_budget_exceeded"] == 1
    assert budgets["candidate_envelope_cpu_budget_exceeded"] == 0
    assert budgets["candidate_envelope_wall_budget_exceeded"] == 0
    assert budgets["candidate_envelope_threads1_wall_budget_exceeded"] == 0


def test_model_compile_is_excluded_from_the_chain_budget():
    # 100.4 s with the compile step, 100.0 s without: only the compile-free total is compared (not above 100).
    budgets, _ = performance_budget.shower_budgets(
        shower_metrics(candidate_envelope_threads1_cpu_s=100.4)
    )
    assert budgets["candidate_envelope_threads1_cpu_budget_exceeded"] == 0
    assert budgets["candidate_envelope_threads1_cpu_ratio"] == pytest.approx(1.0)


def test_equal_statistics_flag_gates_the_chain_budgets():
    budgets, _ = performance_budget.shower_budgets(
        shower_metrics(physics_envelope_photoelectron_bias_limit_exceeded=1)
    )
    assert budgets["equal_statistics_budget_exceeded"] == 1


def test_unmeasured_budget_is_exceeded_with_reason():
    metrics = shower_metrics()
    del metrics["candidate_threads1_propagation_cpu_s"]
    budgets, reasons = performance_budget.shower_budgets(metrics)
    assert budgets["from_corsika_propagation_cpu_budget_exceeded"] == 1
    assert (
        "candidate_threads1_propagation_cpu_s"
        in reasons["from_corsika_propagation_cpu"]
    )


def test_zero_reference_cpu_is_not_a_pass():
    budgets, reasons = performance_budget.shower_budgets(
        shower_metrics(reference_sim_telarray_cpu_s=0.0)
    )
    assert budgets["from_corsika_propagation_cpu_budget_exceeded"] == 1
    assert "from_corsika_propagation_cpu" in reasons


def test_source_wall_budget():
    metrics = {"reference_wall_s": 10.0, "candidate_threads1_wall_s": 9.0}
    budgets, reasons = performance_budget.source_budgets(metrics)
    assert budgets["candidate_threads1_wall_budget_exceeded"] == 0
    assert budgets["candidate_threads1_wall_ratio"] == pytest.approx(0.9)
    metrics["candidate_threads1_wall_s"] = 11.0
    budgets, _ = performance_budget.source_budgets(metrics)
    assert budgets["candidate_threads1_wall_budget_exceeded"] == 1
    budgets, reasons = performance_budget.source_budgets({})
    assert budgets["candidate_threads1_wall_budget_exceeded"] == 1 and reasons


def test_pipe_budgets_are_at_one_gigabyte_per_second():
    budgets, reasons = performance_budget.pipe_budgets(
        {"cat_tee": 1.2, "to_simtel": 0.4}
    )
    assert budgets["cat_tee_pipe_throughput_budget_exceeded"] == 0
    assert budgets["to_simtel_pipe_throughput_budget_exceeded"] == 1
    assert reasons == {}
    budgets, reasons = performance_budget.pipe_budgets({"cat_tee": None})
    assert budgets["cat_tee_pipe_throughput_budget_exceeded"] == 1
    assert "cat_tee_pipe_throughput" in reasons


def test_throughput_is_decimal_gigabytes_per_second():
    assert performance_budget.throughput_gb_per_s(2e9, 1.0) == pytest.approx(2.0)
    assert performance_budget.throughput_gb_per_s(1e9, 0.0) is None


def fake_blocks(monkeypatch, layout):
    """Replace the eventio reader by (type, start, end) blocks laid out in the given order."""
    monkeypatch.setattr(performance, "_top_level_blocks", lambda path: iter(layout))


def test_subset_keeps_header_first_events_and_run_end(
    tmp_path: Path, monkeypatch
) -> None:
    source = tmp_path / "in.iact"
    data = bytes(range(64))
    source.write_bytes(data)
    layout = [
        (1200, 0, 8),  # run header
        (1202, 8, 16),  # event 1
        (1204, 16, 24),
        (1202, 24, 32),  # event 2
        (1204, 32, 40),
        (1202, 40, 48),  # event 3
        (1209, 48, 56),
        (1210, 56, 64),  # run end
    ]
    fake_blocks(monkeypatch, layout)
    target = tmp_path / "out.iact"
    assert performance.subset_showers(source, target, 2) == 2
    assert target.read_bytes() == data[0:40] + data[56:64]


def test_subset_zero_showers_and_negative_count(tmp_path: Path, monkeypatch) -> None:
    source = tmp_path / "in.iact"
    source.write_bytes(bytes(32))
    fake_blocks(monkeypatch, [(1200, 0, 8), (1202, 8, 16), (1210, 16, 32)])
    target = tmp_path / "out.iact"
    assert performance.subset_showers(source, target, 0) == 0
    assert target.read_bytes() == bytes(8) + bytes(16)
    with pytest.raises(ValueError):
        performance.subset_showers(source, target, -1)


def test_select_input_returns_source_without_events(tmp_path: Path) -> None:
    source = tmp_path / "in.iact"
    assert performance.select_input(source, None, tmp_path / "x") == source


def test_metrics_file_is_strict_json(tmp_path: Path) -> None:
    budgets, _ = performance_budget.shower_budgets(shower_metrics())
    performance.common.write_results(budgets, {}, tmp_path)
    assert json.loads((tmp_path / "metrics.json").read_text())


def test_camera_only_line_is_counted():
    text = (
        "Bypass_Optics=3 (camera only) for telescope 1: photon bunches are on the camera plane;\n"
        "Bypass_Optics=1 for telescope 1\n"
    )
    assert performance_budget.count_camera_only(text) == 1
    assert performance_budget.count_camera_only("nothing") == 0
