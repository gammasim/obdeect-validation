"""Tests of the simtools flasher chain pixel comparison (no simulation needed)."""

import numpy as np
import pytest

from obdeect_validation.adapters import simtools_flasher_chain as sfc


def _samples(rng, n_events, mean, scatter=0.0, pixels=200):
    """Events with pixel means ``mean`` and a relative event-to-event scatter."""
    base = np.full(pixels, float(mean))
    return base * (1.0 + scatter * rng.standard_normal((n_events, pixels)))


def test_identical_samples_pass():
    rng = np.random.default_rng(1)
    events = _samples(rng, 20, 1000.0, scatter=0.02)
    metrics = sfc.pixel_comparison(events, events)

    assert metrics["pixels_tested"] == 200
    assert metrics["pixel_ratio_limit_exceeded"] == 0
    assert metrics["image_reduced_chi2"] == pytest.approx(0.0, abs=1e-12)


def test_consistent_samples_are_within_limits():
    rng = np.random.default_rng(2)
    reference = _samples(rng, 40, 1000.0, scatter=0.03)
    candidate = _samples(rng, 40, 1000.0, scatter=0.03)
    metrics = sfc.pixel_comparison(reference, candidate)

    assert metrics["pixel_ratio_limit_exceeded"] == 0
    assert 0.3 < metrics["image_reduced_chi2"] < 2.0


def test_one_deviating_pixel_is_counted():
    rng = np.random.default_rng(3)
    reference = _samples(rng, 40, 1000.0, scatter=0.03)
    candidate = _samples(rng, 40, 1000.0, scatter=0.03)
    candidate[:, 0] *= 1.5  # shape difference far beyond the event scatter
    metrics = sfc.pixel_comparison(reference, candidate)

    assert metrics["pixel_ratio_limit_exceeded"] == 1
    assert metrics["image_reduced_chi2"] > sfc.REDUCED_CHI2_LIMIT


def test_pixels_below_minimum_are_not_tested():
    rng = np.random.default_rng(4)
    reference = _samples(rng, 10, 1000.0, scatter=0.02, pixels=3)
    reference[:, 0] = 1.0  # 1 p.e. per event
    candidate = reference.copy()
    metrics = sfc.pixel_comparison(reference, candidate)

    assert metrics["pixels_tested"] == 2


def test_shares_remove_the_total_scale():
    rng = np.random.default_rng(5)
    reference = _samples(rng, 20, 1000.0, scatter=0.02)
    metrics = sfc.pixel_comparison(reference, 3.0 * reference)

    assert metrics["pixel_ratio_limit_exceeded"] == 0
    assert metrics["image_reduced_chi2"] == pytest.approx(0.0, abs=1e-9)


def test_one_event_is_rejected():
    with pytest.raises(ValueError, match="at least two events"):
        sfc.pixel_comparison(np.ones((1, 5)), np.ones((1, 5)))


def _moment_events(rng, n_events, centroid_x=0.0, width=140.0):
    """Synthetic pixel grid with a Gaussian spot: rows of p.e. per event, pixel x and y [mm]."""
    grid = np.arange(-10, 11) * 50.0
    x, y = (g.ravel() for g in np.meshgrid(grid, grid))
    events = []
    for _ in range(n_events):
        pe = 1000.0 * np.exp(
            -0.5 * ((x - centroid_x) ** 2 + y**2) / (width / 2.355) ** 2
        )
        events.append(rng.poisson(pe).astype(float))
    return np.vstack(events), x, y


def test_event_moments_of_a_centred_gaussian_spot():
    x = np.arange(-10, 11) * 50.0
    xx, yy = (g.ravel() for g in np.meshgrid(x, x))
    pe = np.exp(-0.5 * (xx**2 + yy**2) / 100.0**2)
    row = sfc.event_moments(np.vstack([pe, pe]), xx, yy)[0]
    assert row[1] == pytest.approx(0.0, abs=1e-9)
    assert row[2] == pytest.approx(0.0, abs=1e-9)
    assert row[3] == pytest.approx(row[4], rel=1e-6)
    assert row[3] == pytest.approx(100.0, rel=0.02)


def test_moment_size_gate_has_no_standard_error_condition():
    rng = np.random.default_rng(2)
    ref, x, y = _moment_events(rng, 40, width=140.0)
    cand, _, _ = _moment_events(rng, 40, width=140.0 * 1.005)
    ref_m = sfc.event_moments(ref, x, y)
    cand_m = sfc.event_moments(cand, x, y)
    rel = abs(cand_m[:, 3].mean() - ref_m[:, 3].mean()) / ref_m[:, 3].mean()
    # A 0.5 % width difference is inside the 1 % gate, however many standard errors it is.
    assert rel < sfc.MOMENT_RELATIVE
    assert sfc.MOMENT_RELATIVE == pytest.approx(0.01)


def test_nearest_pitch_of_a_regular_grid():
    x = np.arange(-3, 4) * 50.0
    xx, yy = (g.ravel() for g in np.meshgrid(x, x))
    assert sfc._nearest_pitch_mm(xx, yy) == pytest.approx(50.0)
