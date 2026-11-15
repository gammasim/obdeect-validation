"""Calibration sources versus LightEmission: emission (stage 1) and camera (stage 2).

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.source_camera CASE

Stage 1 compares photon streams at the same plane: photons per event, cx/cy means and KS
distance, arrival time and wavelength moments. Native times are emission times; a photon at
height z (m) with direction cz < 0 reaches the plane z = 0 after -z/cz/c, which is added to the
native times before comparison (time-origin correction, zero for camera-frame output).

Stage 2 propagates the native sky stream, converts it with obdeect stream to-simtel and runs
sim_telarray on both inputs. The two samples are independent (no pairing): photoelectrons per
event as a mean ratio with standard error, and the summed image shape as a per-pixel chi2 with
Poisson variances. If sim_telarray gives no triggered events for either input, stage 2 is BLOCKED
and no camera metric is written.
"""

from __future__ import annotations

import os
import shutil
import sys
import traceback
from pathlib import Path

import numpy as np
from eventio import SimTelFile

from obdeect_validation.adapters import common, readers

C_M_NS = 0.299792458
TIME_RADIUS_CM = 1500.0
MIN_CAMERA_EVENTS = (
    1  # with one event the standard error is the Poisson error of the total
)
KEYS = ("cx", "cy", "cz", "t", "z", "w", "wavelength", "x", "y")

# Matched parameter sets. Flasher cases follow obdeect/validation/sources/compare_*.sh; point
# source follows the xyzls defaults (-r 1500, -s 400, -h 2150 m, -d 0,0,-1); beam scan follows the
# xyls defaults with -r 123 (obdeect beam-scan defaults are the same fractions of 123 cm).
CASES = {
    "flasher-dish": {
        "native": ["obdeect", "sources", "flasher", "--type", "dish", "--events", "20"],
        "reference": ["ff-1m", "--events", "20"],
        "stage2": {
            "path": "propagate",
            "bypass": 1,
            "native": [
                "obdeect",
                "sources",
                "flasher",
                "--type",
                "dish",
                "--events",
                "5",
                "--photons",
                "150000000",
            ],
            "reference": ["ff-1m", "--events", "5", "--photons", "150000000"],
            "photon_factor": 60,
        },
        "frame": {
            "native": "sky (emission at flasher, z = 16 m)",
            "reference": "observation plane z = 0",
        },
        "time_correction": "geometric",
        "native_frame": "sky",
    },
    "flasher-secondary": {
        "native": [
            "obdeect",
            "sources",
            "flasher",
            "--type",
            "secondary",
            "--photons",
            "1e5",
            "--events",
            "1",
        ],
        "reference": ["ff-gct", "--photons", "1e5", "--events", "1"],
        "frame": {
            "native": "telescope (T)",
            "reference": "ff-gct photon bunches (no cz)",
        },
        "time_correction": "none",
        # Replaces the time keys of emission_metrics (see ff_gct_time_metrics).
        "time_model": "ff_gct",
        "time_rule": "ff-gct model: t + |z/cz|/c for all photons (impact at the signed path); "
        "|cz| >= 0.3 on both sides; ref times + time_offset",
        "native_frame": "telescope",
        "stage2": {
            "path": "direct",
            # ff-gct documentation: Bypass_Optics=1 skips the primary only; the secondary
            # reflection stays (Bypass_Optics=2 would skip both mirrors).
            "bypass": 1,
            "cfg": "CTA/CTA-ULTRA6-SST-GCT-S.cfg",
            "altitude": 2100,
            "atmospheric_transmission": "atm_trans_2100_1_10_0_0_2100.dat",
            "camera_map": True,
            "defines": [],
            "native": [
                "obdeect",
                "sources",
                "flasher",
                "--type",
                "secondary",
                "--photons",
                "1000000",
                "--events",
                "1",
            ],
            # 1e6 per flasher: ~57000 p.e. per event. At 1e7 the 1.28 GB event block exceeds the
            # default sim_telarray input buffer (CORSIKA_IO_BUFFER, MAX_BUNCHES), which skips the block.
            "reference": ["ff-gct", "--photons", "1000000", "--events", "1"],
            "photon_factor": 10,
        },
    },
    "point-source": {
        "native": [
            "obdeect",
            "sources",
            "point-source",
            "-N",
            "30",
            "-n",
            "10000000",
            "-z",
            "100000",
            "-a",
            "isotropic",
            "--seed",
            "17",
        ],
        "reference": [
            "xyzls",
            "-N",
            "30",
            "-n",
            "10000000",
            "-z",
            "100000",
            "-a",
            "isotropic",
        ],
        "frame": {
            "native": "sky (emission at source, z = 1000 m)",
            "reference": "plane through telescope sphere centre",
        },
        "time_correction": "geometric_air",
        "native_frame": "sky",
        # Native emits into a cone of 1.02 R (deliberate margin), xyzls into R: compare the photons
        # crossing the plane z = 0 within the same radius, inside both cones.
        "acceptance_radius_cm": 1400.0,
        "stage2": {
            # Photons reach the telescope sphere, not the camera: full optics (simtools: bypass only for flat fielding).
            "path": "propagate",
            "bypass": 0,
            "native": [
                "obdeect",
                "sources",
                "point-source",
                "-N",
                "2",
                "-n",
                "3000000000",
                "-z",
                "100000",
                "-a",
                "isotropic",
                "--seed",
                "17",
            ],
            "reference": [
                "xyzls",
                "-N",
                "2",
                "-n",
                "3000000000",
                "-z",
                "100000",
                "-a",
                "isotropic",
            ],
            "photon_factor": 300,
        },
    },
    "beam-scan": {
        "native": ["obdeect", "sources", "beam-scan", "-r", "123", "--seed", "1"],
        "reference": ["xyls", "-r", "123"],
        "stage2": {
            "blocked": "no intensity option: obdeect beam-scan and xyls emit a fixed 100 photons per "
            "position, so the flat-field intensity cannot be matched between the two"
        },
        "frame": {"native": "sky (beam at z = 6.15 cm)", "reference": "plane z = 0"},
        "time_correction": "geometric",
        "native_frame": "sky",
    },
}


def _read_native(path):
    counts, columns = readers.read_native(path, KEYS)
    return counts, columns


def _read_reference(path):
    counts, columns = readers.read_iact(path)
    # Bunch distance [cm] from the telescope axis, for the camera acceptance cut.
    columns["r"] = 100.0 * np.hypot(columns["x"], columns["y"])
    return counts, columns


def _weighted_mean(x, w):
    total = float(np.sum(w))
    mean = float(np.sum(w * x) / total)
    se = float(np.sqrt(np.sum(w**2 * (x - mean) ** 2)) / total)
    return mean, se


def _weighted_rms(x, w):
    mean, _ = _weighted_mean(x, w)
    return float(np.sqrt(np.sum(w * (x - mean) ** 2) / np.sum(w)))


def _ks(a, b):
    from scipy.stats import ks_2samp

    # Degenerate samples (beam at cx = 0) differ only by round-off, e.g. 0 and 6e-17; snap them.
    a, b = (np.where(np.abs(x) < 1e-12, 0.0, x) for x in (a, b))
    result = ks_2samp(a, b)
    return float(result.statistic), float(result.pvalue)


def _flag(value, limit):
    return int(bool(value > limit))


def accepted_ratio(native, n_events, ref, r_events, radius_cm):
    """Photons per event crossing z = 0 within radius_cm of the telescope axis, native / reference."""
    cz = native["cz"]
    path = np.where(cz < 0, -native["z"] / np.where(cz < 0, cz, -1.0), np.inf)  # m
    radius = 100.0 * np.hypot(
        native["x"] + path * native["cx"], native["y"] + path * native["cy"]
    )
    wn = native["w"][np.isfinite(radius) & (radius < radius_cm)]
    wr = ref["w"][ref["r"] < radius_cm]
    if not wn.size or not wr.size:
        raise ValueError("no photons within the acceptance radius")
    ratio = float((wn.sum() / n_events) / (wr.sum() / r_events))
    # Effective photon numbers (bunches carry weights).
    se = ratio * np.sqrt(np.sum(wn**2) / wn.sum() ** 2 + np.sum(wr**2) / wr.sum() ** 2)
    return ratio, float(se)


ATMPROF_FILE = "atmprof1.dat"  # xyzls default -A 1
# xyzls default -h: observation level of the point-source plane.
ATMPROF_OBSLEV_M = 2150.0


def atmprof_path():
    config = os.environ.get("SIMTEL_CONFIG_PATH", common.REFERENCE / "cfg/common")
    return Path(config) / ATMPROF_FILE


def load_atmprof(path):
    """Heights [m] and n-1 of a sim_telarray atmprof table (columns: alt km, rho, thick, n-1)."""
    rows = [
        line.split()
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    table = np.array([[float(r[0]), float(r[3])] for r in rows])
    return table[:, 0] * 1000.0, table[:, 1]


def index_integral_m(heights_m, index_minus_one, low_m, high_m):
    """Integral of (n-1) dh [m] between two heights; n-1 is log-linear between table levels."""
    h = np.asarray(heights_m, dtype=np.float64)
    y = np.asarray(index_minus_one, dtype=np.float64)
    if not low_m <= high_m or low_m < h[0] or high_m > h[-1]:
        raise ValueError(
            "index integral: heights must lie inside the table, low <= high"
        )
    total = 0.0
    for i in range(h.size - 1):
        a, b = max(low_m, h[i]), min(high_m, h[i + 1])
        if b <= a:
            continue
        k = np.log(y[i + 1] / y[i]) / (h[i + 1] - h[i])
        ya, yb = y[i] * np.exp(k * (a - h[i])), y[i] * np.exp(k * (b - h[i]))
        total += (b - a) * ya if abs(k) * (b - a) < 1e-12 else (yb - ya) / k
    return float(total)


def air_delay_ns(z_m, cz, heights_m, index_minus_one, obslev_m):
    """Extra light time of straight paths in air compared with vacuum, for photons heading to the
    plane (cz < 0) from height z_m above it: (1/c) integral of (n-1) dl, dl = dh / |cz|."""
    z, cz = np.asarray(z_m, dtype=np.float64), np.asarray(cz, dtype=np.float64)
    delay = np.zeros(cz.shape)
    toward = cz < 0
    for z_value in np.unique(z[toward]):
        sel = toward & (z == z_value)
        integral = index_integral_m(
            heights_m, index_minus_one, obslev_m, obslev_m + float(z_value)
        )
        delay[sel] = integral / (C_M_NS * -cz[sel])
    return delay


# flasher-secondary (ff-gct) time model. ff-gct writes every photon at its impact on the plane
# z = 0, including photons heading away from it (impact on the backward extension), with time
# t_em + |z/cz|/c (atmo.c raybnd_: dz / |w| / c). Its raybnd_ also adds a refraction term
# -(t s)^2 / c^3 that does not scale with the 0.5 m path: about -4 ns at |cz| = 0.04, about
# -0.01 ns at |cz| = 0.3. Both samples are therefore compared for |cz| >= FF_GCT_MIN_ABS_CZ.
FF_GCT_MIN_ABS_CZ = 0.3


def _time_block(native_t, native_w, ref_t, ref_w, tol=0.05):
    """Weighted mean and rms of two time samples with their limit flags (as in emission_metrics)."""
    mn, sn = _weighted_mean(native_t, native_w)
    mr, sr = _weighted_mean(ref_t, ref_w)
    se = float(np.hypot(sn, sr))
    rn, rr = _weighted_rms(native_t, native_w), _weighted_rms(ref_t, ref_w)
    se_rms = float(max(rn, rr) / np.sqrt(2 * max(len(native_t), len(ref_t))))
    return {
        "t_mean_native": mn,
        "t_mean_reference": mr,
        "t_mean_delta": mn - mr,
        "t_mean_delta_se": se,
        "t_mean_limit_exceeded": _flag(abs(mn - mr), max(tol, 5 * se)),
        "t_rms_native": rn,
        "t_rms_reference": rr,
        "t_rms_delta": rn - rr,
        "t_rms_limit_exceeded": _flag(abs(rn - rr), max(tol, 5 * se_rms)),
    }


def ff_gct_times(native, ref, min_abs_cz=FF_GCT_MIN_ABS_CZ):
    """Arrival times of the ff-gct time model: obdeect times on the plane z = 0 and the photons
    compared on both sides (impact within TIME_RADIUS_CM, |cz| >= min_abs_cz).

    Returns t_native, native_ok, ref_ok and the native impact-within-radius mask.
    """
    cz = native["cz"]
    if np.any(cz == 0):
        raise ValueError("native photons parallel to the plane z = 0")
    # Impact on the plane along the signed path (as ff-gct), arrival time |path| / c (as ff-gct).
    path = -native["z"] / cz  # m, negative for photons heading away from the plane
    t_native = native["t"] + np.abs(path) / C_M_NS
    xp, yp = native["x"] + path * native["cx"], native["y"] + path * native["cy"]
    inside = 100.0 * np.hypot(xp, yp) < TIME_RADIUS_CM
    native_ok = inside & (np.abs(cz) >= min_abs_cz)
    ref_ok = (ref["r"] < TIME_RADIUS_CM) & (np.abs(ref["cz"]) >= min_abs_cz)
    return t_native, native_ok, ref_ok, inside


def ff_gct_time_metrics(native, ref, min_abs_cz=FF_GCT_MIN_ABS_CZ):
    """Time keys of stage 1 for the ff-gct time model; both samples use the same selection:
    impact within TIME_RADIUS_CM of the axis and |cz| >= min_abs_cz (cz of the reference from cx, cy).
    """
    t_native, native_ok, ref_ok, inside = ff_gct_times(native, ref, min_abs_cz)
    if not native_ok.any() or not ref_ok.any():
        raise ValueError("no photons left for the time comparison")
    metrics = {
        "time_correction_ns_mean": float(np.mean((t_native - native["t"])[native_ok])),
        "native_photons_moving_away_fraction": float(np.mean(native["cz"] >= 0)),
        "t_native_outside_radius_fraction": float(1.0 - np.mean(inside)),
        "t_native_excluded_fraction": float(1.0 - np.mean(native_ok)),
        "t_reference_outside_radius_fraction": float(
            1.0 - np.mean(ref["r"] < TIME_RADIUS_CM)
        ),
        "t_reference_excluded_fraction": float(1.0 - np.mean(ref_ok)),
        "t_min_abs_cz": float(min_abs_cz),
    }
    metrics.update(
        _time_block(
            t_native[native_ok],
            native["w"][native_ok],
            ref["t"][ref_ok],
            ref["w"][ref_ok],
        )
    )
    return metrics


def emission_times(native, ref, time_correction):
    """Arrival times compared in stage 1: obdeect times moved to the reference plane
    (time_correction 'geometric', 'geometric_air' or 'none') and the photons compared.

    Returns t_native, native_ok, ref_ok and the native 'heading for the plane' mask (None
    without a time correction).
    """
    t_native = native["t"].copy()
    native_ok = np.ones(t_native.size, dtype=bool)
    toward = None
    if time_correction in ("geometric", "geometric_air"):
        cz = native["cz"]
        if np.any(cz == 0):
            raise ValueError("native photons parallel to the plane z = 0")
        # Only photons heading for the plane (cz < 0) have a flight path to it (about half of the
        # flasher-secondary photons do not); their time moments are not compared.
        toward = cz < 0
        path = np.where(toward, -native["z"] / np.where(toward, cz, -1.0), 0.0)  # m
        t_native = t_native + path / C_M_NS
        if time_correction == "geometric_air":
            # xyzls adds the light time in air (raybnd_); obdeect emits at the source, no air delay.
            heights, index_minus_one = load_atmprof(atmprof_path())
            t_native = t_native + air_delay_ns(
                native["z"], cz, heights, index_minus_one, ATMPROF_OBSLEV_M
            )
        # Same camera acceptance as the reference: arrival position within TIME_RADIUS_CM.
        xp, yp = native["x"] + path * native["cx"], native["y"] + path * native["cy"]
        native_ok = toward & (100.0 * np.hypot(xp, yp) < TIME_RADIUS_CM)
    # Reference bunches within the fiducial radius of 1500 cm (xyzls -r default). Bunches beyond
    # it (flasher-secondary: 3.6 %, incl. times near -2e8 ns) are not camera photons; their
    # fraction is reported, not compared.
    ref_ok = ref["r"] < TIME_RADIUS_CM
    return t_native, native_ok, ref_ok, toward


def emission_metrics(
    native_counts, native, ref_counts, ref, time_correction, acceptance_radius_cm=None
):
    """Stage 1 numbers and integer limit flags. time_correction is 'geometric' (vacuum path),
    'geometric_air' (vacuum path plus air delay of the atmprof table) or 'none'."""
    n_ev, r_ev = len(native_counts), len(ref_counts)
    metrics = {
        "native_event_count": n_ev,
        "reference_event_count": r_ev,
        "event_count_delta": abs(n_ev - r_ev),
    }
    nc, rc = np.asarray(native_counts), np.asarray(ref_counts)
    if nc.size == 0 or rc.size == 0 or nc.mean() == 0 or rc.mean() == 0:
        raise ValueError("stage 1 needs photons in both streams")
    ratio = float(nc.mean() / rc.mean())
    # Event-to-event scatter; with fewer than two events, the Poisson error of the photon totals.
    if min(nc.size, rc.size) > 1:
        se_native = nc.std(ddof=1) / np.sqrt(nc.size) / nc.mean()
        se_reference = rc.std(ddof=1) / np.sqrt(rc.size) / rc.mean()
    else:
        se_native, se_reference = 1 / np.sqrt(nc.sum()), 1 / np.sqrt(rc.sum())
    se_ratio = float(ratio * np.hypot(se_native, se_reference))
    metrics.update(
        {
            "photons_per_event_native": float(nc.mean()),
            "photons_per_event_reference": float(rc.mean()),
            "photon_ratio": ratio,
            "photon_ratio_se": se_ratio,
            "photon_ratio_limit_exceeded": _flag(
                abs(ratio - 1), max(0.005, 5 * se_ratio)
            ),
        }
    )
    if acceptance_radius_cm:
        accepted, se_accepted = accepted_ratio(
            native, n_ev, ref, r_ev, acceptance_radius_cm
        )
        metrics.update(
            {
                "photon_ratio_accepted": accepted,
                "photon_ratio_accepted_se": se_accepted,
                "photon_ratio_accepted_limit_exceeded": _flag(
                    abs(accepted - 1), max(0.005, 5 * se_accepted)
                ),
            }
        )
    wn, wr = native["w"], ref["w"]
    t_native, native_ok, good, toward = emission_times(native, ref, time_correction)
    if toward is not None:
        metrics["native_photons_moving_away_fraction"] = float(np.mean(~toward))
        metrics["t_native_outside_radius_fraction"] = float(1.0 - np.mean(native_ok))
    metrics["time_correction_ns_mean"] = float(np.mean(t_native - native["t"]))
    for key in ("cx", "cy"):
        mn, sn = _weighted_mean(native[key], wn)
        mr, sr = _weighted_mean(ref[key], wr)
        se = float(np.hypot(sn, sr))
        d_ks, p_ks = _ks(native[key], ref[key])
        metrics.update(
            {
                f"{key}_native": mn,
                f"{key}_reference": mr,
                f"{key}_delta": mn - mr,
                f"{key}_delta_sigma": float((mn - mr) / se) if se else 0.0,
                f"{key}_mean_limit_exceeded": _flag(
                    abs(mn - mr) / se if se else 0.0, 5
                ),
                f"{key}_ks_distance": d_ks,
                f"{key}_ks_pvalue": p_ks,
                f"{key}_ks_limit_exceeded": int(d_ks >= 0.01 and p_ks <= 1e-3),
            }
        )
    metrics["t_reference_outside_radius_fraction"] = (
        float(1.0 - np.mean(good)) if good.size else 0.0
    )
    for key, tol in (("t", 0.05), ("wavelength", 0.01)):
        if key == "t":
            xn, wn_t = t_native[native_ok], wn[native_ok]
            xr, wr_ = ref["t"][good], wr[good]
        else:
            xn, xr, wr_ = native["wavelength"], ref["wavelength"], wr
        wn_k = wn_t if key == "t" else wn
        mn, sn = _weighted_mean(xn, wn_k)
        mr, sr = _weighted_mean(xr, wr_)
        se = float(np.hypot(sn, sr))
        rn, rr = _weighted_rms(xn, wn_k), _weighted_rms(xr, wr_)
        se_rms = float(max(rn, rr) / np.sqrt(2 * max(len(xn), len(xr))))
        metrics.update(
            {
                f"{key}_mean_native": mn,
                f"{key}_mean_reference": mr,
                f"{key}_mean_delta": mn - mr,
                f"{key}_mean_delta_se": se,
                f"{key}_mean_limit_exceeded": _flag(abs(mn - mr), max(tol, 5 * se)),
                f"{key}_rms_native": rn,
                f"{key}_rms_reference": rr,
                f"{key}_rms_delta": rn - rr,
                f"{key}_rms_limit_exceeded": _flag(abs(rn - rr), max(tol, 5 * se_rms)),
            }
        )
    return metrics


def camera_totals(path):
    """Per telescope: (p.e. per triggered event, summed image) of a sim_telarray file."""
    totals, images = {}, {}
    with SimTelFile(str(path)) as stream:
        for event in stream:
            for tel, pe in event.get("photoelectrons", {}).items():
                image = np.asarray(pe["photoelectrons"], dtype=np.float64)
                totals.setdefault(tel, []).append(float(pe["n_pe"]))
                images[tel] = images.get(tel, 0) + image
    return {tel: (np.asarray(totals[tel]), images[tel]) for tel in totals}


def _relative_mean_error(per_event, mean):
    """Relative standard error of a mean over events; with one event, the Poisson error of the
    total photoelectron count (the event-to-event scatter is then unknown)."""
    if per_event.size > 1:
        return per_event.std(ddof=1) / np.sqrt(per_event.size) / mean
    return float(1.0 / np.sqrt(per_event.sum()))


def camera_metrics(reference, candidate):
    """Independent samples: mean p.e. ratio with standard error and summed image chi2 per pixel."""
    shared = sorted(set(reference) & set(candidate))
    if not shared:
        raise ValueError("no telescope with photoelectrons in both samples")
    tel = shared[0]
    na, ima = reference[tel]
    nb, imb = candidate[tel]
    if min(na.size, nb.size) < MIN_CAMERA_EVENTS:
        raise ValueError(
            f"camera comparison needs at least {MIN_CAMERA_EVENTS} triggered events per sample"
        )
    mean_a, mean_b = na.mean(), nb.mean()
    ratio = float(mean_b / mean_a)
    se = float(
        ratio
        * np.hypot(
            _relative_mean_error(na, mean_a),
            _relative_mean_error(nb, mean_b),
        )
    )
    # Pixel shapes normalised to unit sum; variance of a normalised Poisson count.
    pa, pb = ima / ima.sum(), imb / imb.sum()
    var = ima / ima.sum() ** 2 + imb / imb.sum() ** 2
    used = var > 0
    chi2 = float(np.sum((pa[used] - pb[used]) ** 2 / var[used]))
    ndof = max(int(np.count_nonzero(used)) - 1, 1)
    reduced = chi2 / ndof
    return {
        "telescope": int(tel),
        "reference_events": int(na.size),
        "candidate_events": int(nb.size),
        "pe_per_event_reference": float(mean_a),
        "pe_per_event_candidate": float(mean_b),
        "pe_ratio": ratio,
        "pe_ratio_se": se,
        "pe_ratio_limit_exceeded": _flag(abs(ratio - 1), max(0.01, 5 * se)),
        "image_reduced_chi2": reduced,
        "image_chi2_limit_exceeded": _flag(reduced, 2.0),
    }


def _plot(native, ref, output, times):
    """Stage 1 histograms; times = (native, reference) arrival times as compared by the metrics
    (obdeect times moved to the reference plane, same photon selection on both sides)."""
    from obdeect_validation.adapters import plotting

    columns = {
        "cx": (native["cx"], ref["cx"], "cx"),
        "cy": (native["cy"], ref["cy"], "cy"),
        "t": (*times, "arrival time at the reference plane [ns]"),
        "wavelength": (native["wavelength"], ref["wavelength"], "wavelength [nm]"),
    }
    fig, axes = plotting.figure(2, 2, (9, 7))
    for axis, (values_native, values_ref, label) in zip(axes.flat, columns.values()):
        plotting.histograms(
            axis,
            [
                ("obdeect", values_native, plotting.NATIVE),
                ("LightEmission", values_ref, plotting.REFERENCE),
            ],
            label,
            clip=(0.001, 0.999),
        )
    plotting.save(fig, output)


def compared_times(native, ref, spec):
    """Native and reference arrival times as compared by the stage 1 metrics of a case."""
    if "time_model" in spec:
        t_native, native_ok, ref_ok, _ = ff_gct_times(native, ref)
    else:
        t_native, native_ok, ref_ok, _ = emission_times(
            native, ref, spec["time_correction"]
        )
    return t_native[native_ok], ref["t"][ref_ok]


def _simtel(
    root,
    input_file,
    name,
    bypass,
    cfg="CTA/CTA-PROD6-LST.cfg",
    defines=("-DCTA_NORTH", "-DLST1"),
    altitude=2150,
    atmospheric_transmission=None,
    runner=common.command,
):
    """sim_telarray with the options of simtools' light-emission step (released binary)."""
    return common.sim_telarray(
        root,
        input_file,
        name,
        {
            "altitude": altitude,
            **(
                {"atmospheric_transmission": atmospheric_transmission}
                if atmospheric_transmission
                else {}
            ),
            "TRIGGER_TELESCOPES": 1,
            "TELTRIG_MIN_SIGSUM": 2,
            "PULSE_ANALYSIS": -30,
            "MAXIMUM_TELESCOPES": 1,
            "telescope_theta": 0,
            "telescope_phi": 0,
            "Bypass_Optics": bypass,
        },
        config=common.REFERENCE / "cfg" / cfg,
        defines=(*defines, "-DNUM_TELESCOPES=1"),
        runner=runner,
    )


def _camera_stage(root, native_gst, reference_file, stage2):
    """Stage 2: candidate and reference through the same sim_telarray options."""
    candidate = root / "candidate.iact"
    if stage2["path"] == "propagate":
        propagated = root / "propagated.gst"
        common.command(
            [
                "obdeect",
                "propagation",
                "run",
                "--no-extinction",
                "--site-altitude",
                "2150",
                "--seed",
                "1",
                "--wavelength-range",
                "300",
                "800",
                "-i",
                native_gst,
                "-o",
                propagated,
            ],
            root / "propagation.log",
        )
        common.command(
            ["obdeect", "stream", "to-simtel", "-o", candidate, propagated],
            root / "to-simtel.log",
        )
    else:
        argv = ["obdeect", "stream", "to-simtel", "-o", candidate]
        if stage2.get("camera_map"):
            # Telescope 1 of the reference program as a dual-mirror (class 2) telescope in the telescope frame.
            (root / "camera-map.txt").write_text("1 2 0.0\n", encoding="utf-8")
            argv += ["--camera-map", root / "camera-map.txt"]
        common.command([*argv, native_gst], root / "to-simtel.log")
    cfg = stage2.get("cfg", "CTA/CTA-PROD6-LST.cfg")
    defines = stage2.get("defines", ["-DCTA_NORTH", "-DLST1"])
    site = {
        "altitude": stage2.get("altitude", 2150),
        "atmospheric_transmission": stage2.get("atmospheric_transmission"),
    }
    ref_out = _simtel(
        root, reference_file, "reference", stage2["bypass"], cfg, defines, **site
    )
    cand_out = _simtel(
        root, candidate, "candidate", stage2["bypass"], cfg, defines, **site
    )
    a, b = camera_totals(ref_out), camera_totals(cand_out)
    if (
        not a
        or not b
        or any(
            len(v[0]) < MIN_CAMERA_EVENTS for v in list(a.values()) + list(b.values())
        )
    ):
        return None, {
            "reference_triggered_telescopes": {
                str(k): int(len(v[0])) for k, v in a.items()
            },
            "candidate_triggered_telescopes": {
                str(k): int(len(v[0])) for k, v in b.items()
            },
        }
    return camera_metrics(a, b), {}


def _reference_argv(argv, root):
    """LightEmission command with the configuration directory the program needs."""
    config = common.reference_config_dir()
    argv = list(argv)
    if argv[0] in ("ff-1m", "ff-gct", "xyzls"):
        argv.insert(1, f"-I{config}")
    elif argv[0] == "xyls":
        # xyls has no configuration option and reads the atmosphere table from the working directory.
        shutil.copyfile(
            Path(config) / "atm_profile_model_1.dat", root / "atm_profile_model_1.dat"
        )
    return argv


def main(case):
    spec = CASES[case]
    root = Path.cwd()
    native_gst, reference_file = root / "native.gst", root / "reference.dat.gz"
    common.command([*spec["native"], "-o", native_gst], root / "native.log")
    common.command(
        [*_reference_argv(spec["reference"], root), "-o", reference_file],
        root / "reference.log",
    )
    native_counts, native = _read_native(native_gst)
    ref_counts, ref = _read_reference(reference_file)
    metrics = emission_metrics(
        native_counts,
        native,
        ref_counts,
        ref,
        spec["time_correction"],
        spec.get("acceptance_radius_cm"),
    )
    if "time_model" in spec:
        metrics.update(ff_gct_time_metrics(native, ref))
    _plot(native, ref, root / "stage1.png", compared_times(native, ref, spec))
    details = {
        "case": case,
        "native_command": spec["native"],
        "reference_command": spec["reference"],
        "frames": spec["frame"],
        "time_correction": spec["time_correction"],
        "time_correction_rule": (
            spec["time_rule"]
            if "time_rule" in spec
            else "t_arrival = t + (-z/cz)/c + (1/c) int(n-1) dh/|cz| over atmprof1 at 2150 m (z in m, cz<0); ref times + time_offset"
            if spec["time_correction"] == "geometric_air"
            else "t_arrival = t + (-z/cz)/c (z in m, cz<0); ref times + time_offset"
        ),
        "stage2": {"status": "BLOCKED", "reason": "not attempted"},
    }
    stage2 = spec["stage2"]
    try:
        if "blocked" in stage2:
            raise RuntimeError(stage2["blocked"])
        native2, reference2 = (
            root / "stage2-native.gst",
            root / "stage2-reference.dat.gz",
        )
        common.command([*stage2["native"], "-o", native2], root / "stage2-native.log")
        common.command(
            [*_reference_argv(stage2["reference"], root), "-o", reference2],
            root / "stage2-reference.log",
        )
        camera, reason = _camera_stage(root, native2, reference2, stage2)
        if camera is None:
            details["stage2"] = {
                "status": "BLOCKED",
                "reason": f"fewer than {MIN_CAMERA_EVENTS} triggered events per telescope in sim_telarray for one input",
                **reason,
            }
        else:
            details["stage2"] = {"status": "OK"}
            metrics.update(
                {
                    f"camera_{k}": v
                    for k, v in camera.items()
                    if isinstance(v, (int, float))
                }
            )
    except Exception as exc:  # stage 2 is reported, not hidden
        details["stage2"] = {
            "status": "BLOCKED",
            "reason": str(exc)[:400],
            "traceback": traceback.format_exc()[-800:],
        }
    if "native" in stage2:
        details["stage2"].update(
            {
                "native_command": stage2["native"],
                "reference_command": stage2["reference"],
                "photon_factor": stage2["photon_factor"],
                "photon_factor_note": "same factor applied to native and reference; "
                "stage 1 uses the matched stage 1 counts",
            }
        )
    common.write_results(metrics, details, root)


if __name__ == "__main__":
    main(sys.argv[1])
