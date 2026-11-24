import math

import numpy as np

from obdeect_validation.adapters import loss_accounting as la


class Rng:
    """Deterministic stand-in for numpy's Generator with the two methods the populations use."""

    def __init__(self, seed=1):
        self._rng = np.random.default_rng(seed)

    def uniform(self, low=0.0, high=1.0):
        return float(self._rng.uniform(low, high))

    def choice(self, values):
        return values[int(self._rng.integers(len(values)))]


def _summary(by_reason, photons_in, photons_out):
    dropped_by_reason = {r: 0 for r in la.REASONS}
    dropped_by_reason.update(by_reason)
    return {
        "photons_in": photons_in,
        "photons_out": photons_out,
        "dropped": sum(dropped_by_reason.values()),
        "dropped_by_reason": dropped_by_reason,
        "thinned_out": 0.0,
        "weight_out": float(photons_out),
    }


def test_summary_line_parses_counts():
    text = "obdeect propagation run: photons in 160, out 120 (weight 120), dropped 40\n"
    match = la.SUMMARY_LINE.search(text)
    assert [int(v) for v in match.groups()] == [160, 120, 40]


def test_zero_table_bottom_at_site_and_no_row_from_1000_nm(tmp_path):
    table = la.zero_table(tmp_path / "zero.dat")
    lines = table.read_text().splitlines()
    header = next(line for line in lines if "H2=" in line)
    assert "H2= 2.000" in header
    rows = [line.split() for line in lines if line and not line.startswith("#")]
    assert rows and all(int(r[0]) < 1000 for r in rows)
    assert all(float(v) == 0.0 for r in rows for v in r[1:])


def test_directions_are_unit_vectors():
    for zenith in (0.0, 30.0, 90.0):
        cx, cy, cz = la._direction(zenith, 77.0)
        assert math.isclose(cx * cx + cy * cy + cz * cz, 1.0, rel_tol=1e-12)


def test_plane_population_expected_reasons_are_the_constructed_ones():
    categories = la.plane_population(Rng(3))
    counts, reaching, total = la.expected_reasons(categories)
    assert total == 7 * la.N
    assert reaching == 2 * la.N
    assert counts["moving_away"] == 2 * la.N
    assert counts["horizontal"] == la.N
    assert counts["below_ground"] == 2 * la.N
    assert sum(counts.values()) + reaching == total


def test_horizontal_category_has_no_vertical_component():
    categories = {name: rows for name, rows, _ in la.plane_population(Rng(4))}
    horizontal = categories["horizontal"]
    assert all(r["cz"] == 0.0 and r["cx"] == 1.0 for r in horizontal)


def test_sphere_hits_pass_within_the_sphere_radius():
    centre = (0.0, 0.0, 5.0)
    categories = {name: rows for name, rows, _ in la.sphere_population(Rng(5), centre)}
    for row in categories["hit_down"]:
        u = np.array([row["cx"], row["cy"], row["cz"]])
        p = np.array([row["x"], row["y"], row["z"]])
        c = np.array(centre)
        # distance of C from the line through p along u
        offset = np.linalg.norm((c - p) - np.dot(c - p, u) * u)
        assert offset <= la.TELESCOPE_R_M
    for row in categories["side_down"]:
        assert abs(row["x"]) > la.TELESCOPE_R_M


def test_population_metrics_zero_for_a_correct_summary():
    categories = la.plane_population(Rng(6))
    counts, reaching, total = la.expected_reasons(categories)
    summary = _summary(counts, total, reaching)
    metrics, _ = la.population_metrics("p", summary, categories, la.PLANE_EXPECTED)
    assert metrics["p_balance_error"] == 0
    assert metrics["p_dropped_sum_error"] == 0
    assert metrics["p_reason_error"] == 0
    assert metrics["p_reaching_error"] == 0
    assert metrics["p_input_count_error"] == 0
    assert metrics["p_triggered_missing"] == 0


def test_population_metrics_flag_a_missing_drop():
    categories = la.plane_population(Rng(6))
    counts, reaching, total = la.expected_reasons(categories)
    counts["moving_away"] -= 1  # one photon the summary does not report
    summary = _summary(counts, total, reaching)
    summary["photons_out"] += 1  # ... so it is counted as written
    metrics, _ = la.population_metrics("p", summary, categories, la.PLANE_EXPECTED)
    assert metrics["p_balance_error"] == 0
    assert metrics["p_reaching_error"] == 1
    assert metrics["p_reason_error"] == 1


def test_population_metrics_flag_unreported_loss():
    categories = la.plane_population(Rng(6))
    counts, reaching, total = la.expected_reasons(categories)
    summary = _summary(
        counts, total, reaching - 2
    )  # two photons neither written nor dropped
    metrics, _ = la.population_metrics("p", summary, categories, la.PLANE_EXPECTED)
    assert metrics["p_balance_error"] == 2


def test_triggered_missing_counts_zero_reasons():
    categories = la.plane_population(Rng(7))
    counts, reaching, total = la.expected_reasons(categories)
    counts["horizontal"] = 0
    summary = _summary(counts, total, reaching + la.N)
    metrics, _ = la.population_metrics("p", summary, categories, la.PLANE_EXPECTED)
    assert metrics["p_triggered_missing"] == 1


def test_scattered_unaccounted_is_exact_difference():
    assert la.scattered_unaccounted(960, 0, 957, 0) == 3
    assert la.scattered_unaccounted(960, 2, 955, 3) == 0
    assert la.scattered_unaccounted(10, 0, 10, 0) == 0


def test_balance_error_on_summary_dict():
    base = {"photons_in": 250, "photons_out": 100, "dropped": 150, "thinned_out": 0.0}
    assert la.balance_error(base) == 0
    assert la.balance_error({**base, "photons_out": 98}) == 2


def test_balance_error_subtracts_thinned_out():
    # 40 primaries removed by the acceptance draw (thinned_out) are accounted for, not lost.
    summary = {
        "photons_in": 250,
        "photons_out": 60,
        "dropped": 150,
        "thinned_out": 40.0,
    }
    assert la.balance_error(summary) == 0
    summary["thinned_out"] = 39.6  # rounded to the nearest photon
    assert la.balance_error(summary) == 0


def test_balance_error_fails_without_thinned_out():
    import pytest

    with pytest.raises(RuntimeError, match="thinned_out"):
        la.balance_error({"photons_in": 250, "photons_out": 100, "dropped": 150})


def test_scatter_metrics_use_reported_scattered_dropped(tmp_path):
    import pytest

    base_summary = _summary({}, 10, 10)
    scatter_summary = _summary({}, 10, 10)
    scatter_summary["scattered_dropped"] = 3
    scatter_summary["scattered_dropped_by_reason"] = {r: 0 for r in la.REASONS}
    scatter_summary["scattered_thinned_out"] = 0.0
    empty = tmp_path / "empty.pst"
    empty.write_bytes(b"")
    # 10 scattered rows generated, 7 written (photons_out 17 against 10 without scattering), 3
    # reported as dropped. The stream helpers are stubbed: only the arithmetic of the gate is tested.
    original_rows = la.e2e.photon_rows
    original_weight = la.stream_weight
    try:
        la.e2e.photon_rows = lambda path: {"x": np.zeros(0, dtype=np.float32)}
        la.stream_weight = lambda path: (0.0, 0)
        scatter_summary["photons_out"] = 17
        metrics, details = la.scatter_metrics(
            "s", base_summary, scatter_summary, 10, empty, empty, 0
        )
    finally:
        la.e2e.photon_rows = original_rows
        la.stream_weight = original_weight
    assert details["reported_scattered_dropped"] == 3
    assert metrics["s_unaccounted"] == 10 - 0 - 7 - 3

    del scatter_summary["scattered_dropped"]
    la.e2e.photon_rows = lambda path: {"x": np.zeros(0, dtype=np.float32)}
    la.stream_weight = lambda path: (0.0, 0)
    try:
        with pytest.raises(RuntimeError, match="scattered_dropped"):
            la.scatter_metrics("s", base_summary, scatter_summary, 10, empty, empty, 0)
    finally:
        la.e2e.photon_rows = original_rows
        la.stream_weight = original_weight
