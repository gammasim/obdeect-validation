"""simtools flasher chain (obdeect photon stream) against the default light-emission chain.

Reference: ``simtools-simulate-flasher`` with the default chain (LightEmission program of
sim_telarray, then sim_telarray). Candidate: the same flasher with the obdeect chain
(``obdeect_flasher stream_to_simtel sim_telarray:bypass_optics=1``). Both use the light source
LSFS-design, telescope LSTS-01, South, model 7.0.0, and the same flasher photon yield per event.
The two samples are independent (different random seeds), so events are not paired.

Both runs must use the same instrument seed (``--instrument_seed 1``): the instrument seed sets
the per-pixel realisation (gains, mirror deviations). With different instrument seeds the pixel
shapes differ by more than the statistics (two default runs with different seeds: reduced chi2
about 11, 230 pixels beyond 5 sigma; with the same seed: chi2 about 1, none). The simulation seeds
differ, so the samples are independent.

Limits (standard errors from the event-to-event scatter of each sample):
- Total photoelectrons per event: candidate/reference ratio within max(1 %, 5 standard errors),
  as in ``source_camera.camera_metrics`` (flasher-dish case).
- Pixel shape: per pixel, the share of the event total (pixel p.e. over the event p.e.). The
  difference of the mean shares is compared with its standard error. Reduced chi2 over the tested
  pixels at most 2 (the shares sum to one, which removes one degree of freedom).
- Per pixel: |difference| above 5 standard errors (``SIGMA``, as in ``chain_showers``) is
  counted; the expected count is 0. Only pixels with at least ``MIN_PIXEL_PE`` mean p.e. per
  event in both samples are tested.

Moment gate (``--gate moments``, used by the illuminator case): the illuminator spot is
reflected by the mirrors, and the mirror misalignment realisation differs between sim_telarray
(its instrument seed) and obdeect (its model seed). The per-pixel test is therefore not a valid
gate for this case. Reference-vs-reference check (South ILLS-01/LSTS-01, 40 events, instrument
seeds 1 and 2, same simulation settings): per-pixel reduced chi2 70.4 with 74 of 138 pixels
beyond 5 sigma; centroid x differs by 0.73 +- 0.21 mm, width by 1.03 +- 0.09 mm (0.7 %). The
reference-vs-candidate per-pixel chi2 (28.5) is of the same order, so it is recorded only.
The moment gate, per event:
- centroid x and y (pixel units of ``MOMENT_CENTROID_PITCH``): |mean difference| at most
  0.1 pixel pitch and within 5 standard errors (``SIGMA``);
- length and width (sqrt of the eigenvalues of the second moments): |relative mean difference|
  at most 1 %. No standard-error condition: the reference-vs-reference width differs by
  1.03 +- 0.09 mm (0.7 %), i.e. 12 standard errors, so the event scatter does not bound the
  run-to-run (mirror realisation) spread of the size;
- total p.e. ratio: as above.
"""

from pathlib import Path

import numpy as np
from eventio import SimTelFile

from obdeect_validation.adapters import common, source_camera

SIGMA = 5.0  # pixel limit in standard errors of the mean share difference
MIN_PIXEL_PE = 10.0  # mean p.e. per event in both samples for a pixel to be tested
REDUCED_CHI2_LIMIT = 2.0
TELESCOPE = 0  # sim_telarray index of the one telescope of the flasher runs
MOMENT_CENTROID_PITCH = 0.1  # centroid tolerance in pixel pitches
MOMENT_RELATIVE = 0.01  # length and width tolerance (relative)


def event_images(path, telescope=TELESCOPE):
    """Return the photoelectrons per pixel of each triggered event: array (events, pixels)."""
    images = []
    with SimTelFile(str(path)) as stream:
        for event in stream:
            pe = event.get("photoelectrons", {}).get(telescope)
            if pe is not None:
                images.append(np.asarray(pe["photoelectrons"], dtype=np.float64))
    if not images:
        raise ValueError(f"no triggered events of telescope {telescope} in {path}")
    return np.vstack(images)


def _mean_and_error(values):
    """Mean over events and its standard error (event-to-event scatter)."""
    if values.shape[0] < 2:
        raise ValueError("pixel comparison needs at least two events per sample")
    return values.mean(axis=0), values.std(axis=0, ddof=1) / np.sqrt(values.shape[0])


def pixel_comparison(reference_events, candidate_events, min_pe=MIN_PIXEL_PE):
    """
    Compare the pixel shapes of two samples of events and return the metrics.

    Parameters
    ----------
    reference_events, candidate_events : array
        Photoelectrons per pixel, one row per event (independent samples).
    min_pe : float
        Minimum mean photoelectrons per event in both samples for a pixel to be tested.

    Returns
    -------
    dict
        Tested pixels, pixels beyond the limit, the largest significance and the reduced chi2.
    """
    reference = np.asarray(reference_events, dtype=float)
    candidate = np.asarray(candidate_events, dtype=float)
    ref_share = reference / reference.sum(axis=1, keepdims=True)
    cand_share = candidate / candidate.sum(axis=1, keepdims=True)
    ref_mean, ref_error = _mean_and_error(ref_share)
    cand_mean, cand_error = _mean_and_error(cand_share)
    tested = (
        (reference.mean(axis=0) >= min_pe)
        & (candidate.mean(axis=0) >= min_pe)
        & (ref_error > 0)
        & (cand_error > 0)
    )
    significance = np.abs(cand_mean - ref_mean)[tested] / np.hypot(
        ref_error[tested], cand_error[tested]
    )
    ndof = max(int(tested.sum()) - 1, 1)
    return {
        "pixels_tested": int(tested.sum()),
        "pixel_ratio_limit_exceeded": int(np.count_nonzero(significance > SIGMA)),
        "pixel_max_significance": float(significance.max())
        if significance.size
        else 0.0,
        "image_reduced_chi2": float(np.sum(significance**2) / ndof),
    }


def camera_positions_mm(path, telescope=TELESCOPE):
    """Pixel x and y [mm] of one telescope (sim_telarray camera settings, metres)."""
    with SimTelFile(str(path)) as stream:
        for _ in stream:
            break
        description = stream.telescope_descriptions[
            next(iter(stream.telescope_descriptions))
        ]
        settings = description["camera_settings"]
        return (
            np.asarray(settings["pixel_x"], dtype=np.float64) * 1000.0,
            np.asarray(settings["pixel_y"], dtype=np.float64) * 1000.0,
        )


def event_moments(images, x, y):
    """Per event: total p.e., centroid x, y, length, width [mm] (rows of the array)."""
    rows = []
    for pe in images:
        total = pe.sum()
        if total <= 0:
            continue
        cx = (pe * x).sum() / total
        cy = (pe * y).sum() / total
        dx, dy = x - cx, y - cy
        sxx = (pe * dx * dx).sum() / total
        syy = (pe * dy * dy).sum() / total
        sxy = (pe * dx * dy).sum() / total
        eigenvalues = np.linalg.eigvalsh(np.array([[sxx, sxy], [sxy, syy]]))
        rows.append((total, cx, cy, np.sqrt(eigenvalues[1]), np.sqrt(eigenvalues[0])))
    if len(rows) < 2:
        raise ValueError(
            "moment comparison needs at least two triggered events per sample"
        )
    return np.array(rows)


def moment_comparison(reference_path, candidate_path):
    """Moment gate (see module docstring): counts of exceeded limits and the mean differences."""
    x, y = camera_positions_mm(reference_path)
    x_cand, y_cand = camera_positions_mm(candidate_path)
    if not (np.allclose(x, x_cand) and np.allclose(y, y_cand)):
        raise ValueError(
            "reference and candidate cameras differ: moment gate not applicable"
        )
    nearest = _nearest_pitch_mm(x, y)
    ref = event_moments(event_images(reference_path), x, y)
    cand = event_moments(event_images(candidate_path), x, y)
    ref_mean, ref_error = _mean_and_error(ref)
    cand_mean, cand_error = _mean_and_error(cand)
    difference = cand_mean - ref_mean
    error = np.hypot(ref_error, cand_error)
    centroid_ok = (np.abs(difference[1:3]) <= MOMENT_CENTROID_PITCH * nearest) & (
        np.abs(difference[1:3]) <= SIGMA * error[1:3]
    )
    size_ok = np.abs(difference[3:5]) <= MOMENT_RELATIVE * ref_mean[3:5]
    pe_limit = max(0.01 * ref_mean[0], SIGMA * error[0])
    return {
        "centroid_limit_exceeded": int(np.count_nonzero(~centroid_ok)),
        "width_length_limit_exceeded": int(np.count_nonzero(~size_ok)),
        "pe_total_difference_within_limit": bool(abs(difference[0]) <= pe_limit),
        "centroid_pitch_mm": nearest,
        "mean_reference": dict(zip(MOMENT_NAMES, ref_mean.tolist(), strict=True)),
        "mean_candidate": dict(zip(MOMENT_NAMES, cand_mean.tolist(), strict=True)),
        "difference": dict(zip(MOMENT_NAMES, difference.tolist(), strict=True)),
        "error": dict(zip(MOMENT_NAMES, error.tolist(), strict=True)),
    }


MOMENT_NAMES = ("pe_total", "centroid_x_mm", "centroid_y_mm", "length_mm", "width_mm")


def _nearest_pitch_mm(x, y):
    """Median nearest-neighbour distance of the pixels [mm] (the pixel pitch)."""
    points = np.column_stack([x, y])
    distances = np.sqrt(((points[:, None, :] - points[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(distances, np.inf)
    return float(np.median(distances.min(axis=1)))


def camera_comparison(reference_path, candidate_path):
    """Return (metrics, details) of the comparison of two sim_telarray outputs."""
    reference = source_camera.camera_totals(reference_path)
    candidate = source_camera.camera_totals(candidate_path)
    camera = source_camera.camera_metrics(reference, candidate)
    pixels = pixel_comparison(
        event_images(reference_path), event_images(candidate_path)
    )
    reference_events = reference[TELESCOPE][0].size
    candidate_events = candidate[TELESCOPE][0].size
    metrics = {
        "candidate_pipeline_complete": 1,
        "reference_event_count": int(reference_events),
        "candidate_event_count": int(candidate_events),
        "event_count_delta": int(candidate_events - reference_events),
        "pe_ratio_limit_exceeded": camera["pe_ratio_limit_exceeded"],
        "image_chi2_limit_exceeded": int(
            pixels["image_reduced_chi2"] > REDUCED_CHI2_LIMIT
        ),
        "pixel_ratio_limit_exceeded": pixels["pixel_ratio_limit_exceeded"],
    }
    details = {
        "camera": camera,
        "pixels": pixels,
        "limits": {
            "pe_ratio": "max(1 %, 5 standard errors)",
            "image_reduced_chi2": REDUCED_CHI2_LIMIT,
            "pixel_significance": SIGMA,
            "min_pixel_pe_per_event": MIN_PIXEL_PE,
        },
        "reference": str(reference_path),
        "candidate": str(candidate_path),
    }
    return metrics, details


def moment_camera_comparison(reference_path, candidate_path):
    """Return (metrics, details) with the moment gate; per-pixel chi2 is recorded only."""
    reference = source_camera.camera_totals(reference_path)
    candidate = source_camera.camera_totals(candidate_path)
    camera = source_camera.camera_metrics(reference, candidate)
    moments = moment_comparison(reference_path, candidate_path)
    pixels = pixel_comparison(
        event_images(reference_path), event_images(candidate_path)
    )
    reference_events = reference[TELESCOPE][0].size
    candidate_events = candidate[TELESCOPE][0].size
    metrics = {
        "candidate_pipeline_complete": 1,
        "reference_event_count": int(reference_events),
        "candidate_event_count": int(candidate_events),
        "event_count_delta": int(candidate_events - reference_events),
        "pe_ratio_limit_exceeded": camera["pe_ratio_limit_exceeded"],
        "centroid_limit_exceeded": moments["centroid_limit_exceeded"],
        "width_length_limit_exceeded": moments["width_length_limit_exceeded"],
    }
    details = {
        "camera": camera,
        "moments": moments,
        "pixels_recorded_only": pixels,
        "limits": {
            "pe_ratio": "max(1 %, 5 standard errors)",
            "centroid": f"{MOMENT_CENTROID_PITCH} pixel pitch and 5 standard errors",
            "width_length": f"{MOMENT_RELATIVE:.0%} relative (no standard-error condition)",
        },
        "reference": str(reference_path),
        "candidate": str(candidate_path),
    }
    return metrics, details


def _add_arguments(parser):
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--gate", choices=("pixels", "moments"), default="pixels")


def analyse(root, args):
    """Metrics and details of the comparison; a candidate that did not produce output is blocked."""
    if not args.candidate.exists():
        return (
            {"candidate_pipeline_complete": 0},
            {"candidate_blocked": "candidate sim_telarray output missing"},
        )
    if args.gate == "moments":
        return moment_camera_comparison(
            args.reference.resolve(), args.candidate.resolve()
        )
    return camera_comparison(args.reference.resolve(), args.candidate.resolve())


def main():
    """Compare the two flasher outputs and write metrics.json and details.json."""
    common.results_main(analyse, __doc__, add_arguments=_add_arguments)


if __name__ == "__main__":
    main()
