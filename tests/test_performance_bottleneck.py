import sys

import pytest

from obdeect_validation.adapters import performance_bottleneck as pb

NAMES = ["a", "b", "c"]


def pipe_record(wall, cpus):
    return {
        "wall_s": wall,
        "cpu_s": sum(cpus.values()),
        "children": [{"name": n, "cpu_s": c} for n, c in cpus.items()],
    }


def records_for(walls, cpus, reps=1):
    """Records of prefixes 1..len(walls): walls[k-1] is the wall of prefix k, repeated reps times."""
    return {
        k: [pipe_record(wall, cpus) for _ in range(reps)]
        for k, wall in enumerate(walls, start=1)
    }


def test_prefixes_are_the_first_k_members() -> None:
    stages = [("a", ["x"]), ("b", ["y"]), ("c", ["z"])]
    assert [len(part) for part in pb.prefixes(stages)] == [1, 2, 3]
    assert pb.prefixes(stages)[1] == stages[:2]


def test_marginal_walls_shares_and_bottleneck() -> None:
    cpus = {"a": 0.5, "b": 2.0, "c": 0.1}
    metrics = pb.bottleneck_metrics(NAMES, records_for([1.0, 3.0, 3.5], cpus), reps=1)
    assert metrics["prefix1_completed"] == metrics["prefix3_completed"] == 1
    assert metrics["prefix2_wall_s"] == 3.0
    assert metrics["prefix2_cpu_s"] == pytest.approx(2.6)
    assert metrics["bottleneck_a_marginal_wall_s"] == pytest.approx(1.0)
    assert metrics["bottleneck_b_marginal_wall_s"] == pytest.approx(2.0)
    assert metrics["bottleneck_c_marginal_wall_s"] == pytest.approx(0.5)
    assert metrics["bottleneck_b_wall_share"] == pytest.approx(2.0 / 3.5)
    assert metrics["bottleneck_b_cpu_s"] == 2.0
    assert metrics["bottleneck_stage_index"] == 2


def test_negative_marginal_keeps_its_sign() -> None:
    metrics = pb.bottleneck_metrics(
        NAMES, records_for([1.0, 3.0, 2.5], {"a": 1, "b": 1, "c": 1}), reps=1
    )
    assert metrics["bottleneck_c_marginal_wall_s"] == pytest.approx(-0.5)
    assert metrics["bottleneck_c_wall_share"] == pytest.approx(-0.5 / 2.5)
    assert metrics["bottleneck_stage_index"] == 2


def test_wall_is_median_and_spread_uses_full_pipe() -> None:
    cpus = {"a": 1.0, "b": 1.0, "c": 1.0}
    records = {
        1: [pipe_record(1.0, cpus)] * 3,
        2: [pipe_record(2.0, cpus)] * 3,
        3: [pipe_record(w, cpus) for w in (2.0, 3.0, 4.0)],
    }
    metrics = pb.bottleneck_metrics(NAMES, records, reps=3)
    assert metrics["prefix3_wall_s"] == 3.0
    assert metrics["bottleneck_wall_s_spread"] == pytest.approx((4.0 - 2.0) / 3.0)


def test_incomplete_prefix_gives_none_and_is_not_scored() -> None:
    cpus = {"a": 1.0, "b": 1.0, "c": 1.0}
    records = records_for([1.0, 3.0], cpus)
    records[3] = []
    metrics = pb.bottleneck_metrics(NAMES, records, reps=1)
    assert metrics["prefix3_completed"] == 0
    assert metrics["prefix3_wall_s"] is None
    assert metrics["bottleneck_c_marginal_wall_s"] is None
    assert metrics["bottleneck_c_wall_share"] is None
    assert metrics["bottleneck_wall_s_spread"] is None
    assert metrics["bottleneck_stage_index"] == 2


def test_run_prefixes_runs_shell_pipes_and_records_failures(tmp_path) -> None:
    python = sys.executable
    relay = [python, "-c", "import sys; sys.stdout.write(sys.stdin.read())"]
    stages = [
        ("produce", [python, "-c", "import sys; sys.stdout.write('x' * 5000)"]),
        ("relay", relay),
        ("last", relay),
    ]
    failures = {}
    records = pb.run_prefixes(stages, 2, tmp_path, failures)
    assert failures == {}
    assert {k: len(v) for k, v in records.items()} == {1: 2, 2: 2, 3: 2}
    assert [c["name"] for c in records[3][0]["children"]] == [
        "produce",
        "relay",
        "last",
    ]
    assert all(r["wall_s"] > 0 for r in records[3])


def test_run_prefixes_records_a_failing_member(tmp_path) -> None:
    python = sys.executable
    stages = [
        ("produce", [python, "-c", "pass"]),
        ("broken", [python, "-c", "import sys; sys.exit(3)"]),
    ]
    failures = {}
    records = pb.run_prefixes(stages, 1, tmp_path, failures)
    assert len(records[1]) == 1
    assert records[2] == []
    assert list(failures) == [2]
    assert "broken exited 3" in failures[2][0]
