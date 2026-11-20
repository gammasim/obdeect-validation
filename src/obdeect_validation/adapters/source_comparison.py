"""Compare obdeect light sources with sim_telarray LightEmission and analytic models.

Each invocation writes metrics.json, details.json and one PNG plot into the current
case directory. Cases: muon-track (fake-muon), night-sky (nsbls), laser-beam (ls-beam).
The flasher, point-source and beam-scan cases are in source_camera.
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common, readers
from obdeect_validation.adapters.chain_showers import camera_events, simulate
from obdeect_validation.adapters.common import REFERENCE, run_tool, tool

C_M_NS = 0.299792458
ALPHA = 0.0072973526
KS_ALPHA = 1e-3
RATIO_FLOOR = 0.01
FWHM_TO_SIGMA = 2.35482004503095

# Physical configuration shared by obdeect and sim_telarray (see details.json for the mapping).
MUON = {
    "radius_m": 15.0,
    "altitude_m": 2150.0,
    "height_m": 3000.0,
    "theta_deg": 0.0,
    "phi_deg": 0.0,
    "wavelength_nm": (300.0, 600.0),
    "events": 10,
    "seed": 17,
    "index_step_m": 100.0,
}
# Old night-sky jitter: uniform in +-pulse_width/2.35482 about 0 with pulse_width = FWHM = 2.35482 ns,
# i.e. a top hat of 1 ns centred at 0.
NIGHT_PULSE = "tophat:width_ns=1"
NIGHT = {
    "radius_m": 15.0,
    "altitude_m": 2150.0,
    "viewcone_deg": 10.0,
    "spectrum_nm": 400.0,
    "band_nm": (300.0, 650.0),
    # 10 x the nsbls help-text rate (0.24): at 0.24 the camera never triggers (0 events in 4).
    # The same configured rate is used for both programs; yield and rate are checked against it.
    "rate_per_sr_ns_cm2": 2.4,
    "pulse_fwhm_ns": FWHM_TO_SIGMA,
    "events": 4,
    "bunch_size": 2,
    "seed": 17,
}
LASER = {
    "laser_position_m": (1000.0, 0.0, 0.0),
    "laser_theta_deg": 60.0,
    "laser_phi_deg": 180.0,
    "telescope_theta_deg": 30.0,
    "telescope_phi_deg": 0.0,
    "telescope_radius_m": 12.5,
    "fov_deg": 15.0,
    "spectrum_nm": 355.0,
    "primary_per_m": 1000.0,
    "step_m": 0.1,
    "brightness_per_m": 1000.0,
    "seed": 17,
    "altitude_m": 2150.0,
}


# ---------------------------------------------------------------- statistics helpers


def ks_statistic(a, b):
    """Two-sample Kolmogorov-Smirnov distance."""
    a, b = np.sort(np.asarray(a, dtype=float)), np.sort(np.asarray(b, dtype=float))
    grid = np.concatenate((a, b))
    cdf_a = np.searchsorted(a, grid, side="right") / a.size
    cdf_b = np.searchsorted(b, grid, side="right") / b.size
    return float(np.max(np.abs(cdf_a - cdf_b)))


def ks_one_sample(a, cdf):
    """One-sample Kolmogorov-Smirnov distance against an analytic CDF."""
    a = np.sort(np.asarray(a, dtype=float))
    n = a.size
    f = cdf(a)
    return float(
        max(np.max(np.arange(1, n + 1) / n - f), np.max(f - np.arange(0, n) / n))
    )


def ks_critical(n, m=None, alpha=KS_ALPHA):
    """Critical KS distance at significance alpha (asymptotic); one or two samples."""
    c = math.sqrt(-0.5 * math.log(alpha / 2.0))
    return c / math.sqrt(n) if m is None else c * math.sqrt((n + m) / (n * m))


def ratio_limit(ratio_se):
    return max(RATIO_FLOOR, 5.0 * ratio_se)


def mean_se(values):
    values = np.asarray(values, dtype=float)
    return float(np.mean(values)), float(
        np.std(values, ddof=1) / math.sqrt(values.size)
    )


# ---------------------------------------------------------------- physics helpers


def frank_tamm_per_metre(n_minus_one, wl_min_nm, wl_max_nm, beta=1.0, charge=1.0):
    """Cherenkov photons per metre of track in [wl_min, wl_max] at constant index, into 4 pi."""
    cos_c = 1.0 / ((1.0 + n_minus_one) * beta)
    return (
        2.0
        * math.pi
        * ALPHA
        * charge**2
        * (1.0 - cos_c**2)
        * (1.0 / wl_min_nm - 1.0 / wl_max_nm)
        * 1e9
    )


def cherenkov_cos(n_minus_one, beta=1.0):
    return 1.0 / ((1.0 + n_minus_one) * beta)


def night_sky_photons(rate_per_sr_ns_cm2, viewcone_deg, radius_m, window_ns):
    """Photons per event that give the configured rate on the telescope sphere."""
    omega = 2.0 * math.pi * (1.0 - math.cos(math.radians(viewcone_deg)))
    area_cm2 = math.pi * (100.0 * radius_m) ** 2
    return rate_per_sr_ns_cm2 * omega * area_cm2 * window_ns


def observation_crossing(p, d):
    """Point where each ray (from p along unit d) crosses the observation plane z = 0."""
    return p + (-p[:, 2] / d[:, 2])[:, None] * d


def impact_parameter(q, d):
    """Distance from the origin to each ray; any point on the ray gives the same value."""
    return np.linalg.norm(q - np.sum(q * d, axis=1)[:, None] * d, axis=1)


def closest_approach(points, d):
    """Closest-approach distance of the ray line through each point along unit d: |P x d|.

    Equal for every point on the same ray, so sphere entry points and observation-plane
    crossings give the same value for the same ray.
    """
    return np.linalg.norm(np.cross(points, d), axis=1)


def centre_plane_time(t, points, d, speed_m_per_ns=C_M_NS):
    """Arrival time at the plane through the origin perpendicular to d: t - (X . d) / c.

    t is the time at the point X on the ray; X may be any point on the ray.
    """
    return t - np.sum(points * d, axis=1) / speed_m_per_ns


# ---------------------------------------------------------------- data readers


def read_native(path):
    """obdeect photon stream: positions m, times ns, wavelengths nm, weights (missing columns empty)."""
    # Only the columns are used here, not the per-event weights.
    return readers.read_native(path, require_events=False)[1]


def read_reference(path):
    """EventIO IACT photon bunches: positions cm -> m, zem absolute m, time includes the event offset."""
    return readers.read_iact(path)[1]


def positions_n(data):
    return np.stack([data["x"], data["y"], data["z"]], axis=1)


def unit_vectors(data):
    return np.stack([data["cx"], data["cy"], data["cz"]], axis=1)


def fake_muon_argv(output, events):
    return [
        tool("fake-muon"),
        f"-I{REFERENCE / 'cfg/common'}",
        "-N",
        str(events),
        "-h",
        str(MUON["height_m"]),
        "--theta",
        str(MUON["theta_deg"]),
        "--phi",
        str(MUON["phi_deg"]),
        "-x",
        "0",
        "-y",
        "0",
        "-r",
        str(MUON["radius_m"]),
        "-a",
        str(MUON["altitude_m"]),
        "-l",
        str(MUON["wavelength_nm"][0]),
        "-L",
        str(MUON["wavelength_nm"][1]),
        "-A",
        "0",
        "-b",
        "1",
        "-o",
        output,
    ]


def index_table(work):
    """Refractive index above the observation level from the sim_telarray atmosphere (fake-muon printout)."""
    heights = np.arange(
        0.0, MUON["height_m"] + 0.5 * MUON["index_step_m"], MUON["index_step_m"]
    )
    values = []
    for height in heights:
        argv = [
            tool("fake-muon"),
            f"-I{REFERENCE / 'cfg/common'}",
            "-N",
            "1",
            "-h",
            str(int(height)),
            "-A",
            "0",
            "-o",
            str(work / "index-probe.iact"),
        ]
        text = run_tool(argv, work / "index-probe.log")
        values.append(
            float(re.search(r"at altitude level [0-9.]+ m is ([0-9.]+)", text).group(1))
            - 1.0
        )
    table = np.column_stack((heights, np.array(values)))
    np.savetxt(work / "index.txt", table, fmt=["%.1f", "%.9e"])
    return table


def muon_expected_hits(table, step_m=0.5, azimuths=3600):
    """Expected Cherenkov photons per event on the telescope sphere (numerical, same geometry)."""
    radius = MUON["radius_m"]
    wl1, wl2 = MUON["wavelength_nm"]
    theta, phi = math.radians(MUON["theta_deg"]), math.radians(MUON["phi_deg"])
    direction = np.array(
        [
            -math.sin(theta) * math.cos(phi),
            -math.sin(theta) * math.sin(phi),
            -math.cos(theta),
        ]
    )
    e1 = np.array(
        [
            math.cos(theta) * math.cos(phi),
            math.cos(theta) * math.sin(phi),
            -math.sin(theta),
        ]
    )
    e2 = np.array([-math.sin(phi), math.cos(phi), 0.0])
    length = MUON["height_m"] / math.cos(theta)
    s = (np.arange(int(round(length / step_m))) + 0.5) * step_m
    height = s * math.cos(theta)
    n_minus_one = np.interp(height, table[:, 0], table[:, 1])
    cos_c = cherenkov_cos(n_minus_one)
    sin_c = np.sqrt(1.0 - cos_c**2)
    dn_ds = 2.0 * math.pi * ALPHA * sin_c**2 * (1.0 / wl1 - 1.0 / wl2) * 1e9
    psi = (np.arange(azimuths) + 0.5) * 2.0 * math.pi / azimuths
    total = 0.0
    for chunk in np.array_split(np.arange(s.size), max(1, s.size // 400)):
        p = -s[chunk][:, None] * direction[None, :]
        d = cos_c[chunk][:, None, None] * direction[None, None, :] + sin_c[chunk][
            :, None, None
        ] * (np.cos(psi)[None, :, None] * e1 + np.sin(psi)[None, :, None] * e2)
        pd = np.sum(p[:, None, :] * d, axis=2)
        pp = np.sum(p * p, axis=1)[:, None]
        disc = pd * pd - (pp - radius * radius)
        hit = (disc >= 0.0) & (-pd + np.sqrt(np.clip(disc, 0.0, None)) > 0.0)
        total += float(np.sum(dn_ds[chunk] * step_m * hit.mean(axis=1)))
    return total


def muon_air_integral_m(table, heights_m, step_m=0.5):
    """Integral of (n-1) dh [m] from the observation level to each height [m]; table rows (h [m], n-1)."""
    heights_m = np.asarray(heights_m, dtype=float)
    grid = np.arange(0.0, float(heights_m.max()) + step_m, step_m)
    n_minus_one = np.interp(grid, table[:, 0], table[:, 1])
    cumulative = np.concatenate(
        ([0.0], np.cumsum(0.5 * (n_minus_one[1:] + n_minus_one[:-1]) * np.diff(grid)))
    )
    return np.interp(heights_m, grid, cumulative)


def muon_arrival_expected_ns(heights_m, abs_cz, table):
    """Arrival time [ns] at the telescope plane z = 0 of Cherenkov photons from heights h [m].

    Emission time is -h/c (track time zero at impact, beta = 1). The straight path to z = 0 has length
    h/|cz| in vacuum; the air excess is (1/|cz|) integral (n-1) dh, so t = (-h + (h + air) / |cz|) / c.
    """
    h = np.asarray(heights_m, dtype=float)
    cz = np.asarray(abs_cz, dtype=float)
    air = muon_air_integral_m(table, h)
    return (-h + (h + air) / cz) / C_M_NS


def camera_pe(path):
    """Photo-electrons per triggered event, summed over telescopes."""
    values = []
    for event in camera_events(path).values():
        pe = event.get("photoelectrons", {})
        values.append(float(sum(pe[t]["n_pe"] for t in pe)))
    return np.asarray(values, dtype=float)


def camera_comparison(native_pe, reference_pe):
    """Independent samples: relative difference of mean p.e. per event with its standard error."""
    if native_pe.size < 2 or reference_pe.size < 2:
        raise ValueError(
            "camera comparison needs at least two triggered events per sample"
        )
    mean_a, mean_b = native_pe.mean(), reference_pe.mean()
    se_a = native_pe.std(ddof=1) / math.sqrt(native_pe.size)
    se_b = reference_pe.std(ddof=1) / math.sqrt(reference_pe.size)
    ratio = mean_a / mean_b
    ratio_se = math.sqrt((se_a / mean_b) ** 2 + (mean_a * se_b / mean_b**2) ** 2)
    return {
        "camera_events_native": int(native_pe.size),
        "camera_events_reference": int(reference_pe.size),
        "camera_pe_mean_native": float(mean_a),
        "camera_pe_mean_reference": float(mean_b),
        "camera_pe_ratio": float(ratio),
        "camera_pe_ratio_se": float(ratio_se),
        "camera_pe_ratio_limit_exceeded": int(abs(ratio - 1.0) > ratio_limit(ratio_se)),
    }


def flag(metrics, name, value, limit):
    """Store a raw value and an integer flag set when the value exceeds the limit."""
    metrics[name] = float(value)
    metrics[f"{name}_limit_exceeded"] = int(abs(value) > limit)


def ks_flags(metrics, name, statistic, critical):
    metrics[f"{name}_ks"] = float(statistic)
    metrics[f"{name}_ks_critical"] = float(critical)
    metrics[f"{name}_ks_limit_exceeded"] = int(statistic > critical)


def write_outputs(metrics, details, plot, plot_path, rows=2, columns=2):
    """metrics.json, details.json and the figure; plot(axes) draws on rows x columns axes."""
    from obdeect_validation.adapters import plotting

    metrics = {
        k: (int(v) if isinstance(v, (bool, np.bool_)) else v)
        for k, v in metrics.items()
    }
    common.write_results(metrics, details)
    fig, axes = plotting.figure(rows, columns, (5 * columns, 3.6 * rows))
    plot(axes)
    plotting.save(fig, plot_path)


def time_difference(native_times, reference_times):
    """Mean native - reference time and its standard error [ns]."""
    mean = float(np.mean(native_times) - np.mean(reference_times))
    se = math.sqrt(
        np.var(native_times) / len(native_times)
        + np.var(reference_times) / len(reference_times)
    )
    return mean, se


def propagate_to_simtel(
    work, native_file, altitude_m, wavelength_nm, seed, profile=None
):
    """obdeect propagation without extinction, then stream to-simtel; returns the sim_telarray input."""
    cfg_dir = REFERENCE / "cfg/common"
    propagated, camera_input = work / "native_prop.gst", work / "native.iact"
    run_tool(
        [
            tool("obdeect"),
            "propagation",
            "run",
            "--no-extinction",
            "--site-altitude",
            str(altitude_m),
            *(["--profile", str(profile)] if profile else []),
            "--wavelength-range",
            *map(str, wavelength_nm),
            "--seed",
            str(seed),
            "-i",
            str(native_file),
            "-o",
            str(propagated),
        ],
        work / "propagate.log",
    )
    run_tool(
        [
            tool("obdeect"),
            "stream",
            "to-simtel",
            "--atmprof",
            str(cfg_dir / "atmprof1.dat"),
            "-o",
            str(camera_input),
            str(propagated),
        ],
        work / "to-simtel.log",
    )
    return camera_input


# ---------------------------------------------------------------- muon-track


def muon_track(work):
    table = index_table(work)
    native_file = work / "native_sources.gst"
    run_tool(
        [
            tool("obdeect"),
            "sources",
            "run",
            "muon",
            "--observation-level",
            str(MUON["altitude_m"]),
            "--seed",
            str(MUON["seed"]),
            "--events",
            str(MUON["events"]),
            "--primary",
            "muon",
            "--zenith",
            str(MUON["theta_deg"]),
            "--azimuth",
            str(MUON["phi_deg"]),
            "--start-height",
            str(MUON["height_m"]),
            "--impact-x",
            "0",
            "--impact-y",
            "0",
            "--telescope-radius",
            str(MUON["radius_m"]),
            "--wavelength-min",
            str(MUON["wavelength_nm"][0]),
            "--wavelength-max",
            str(MUON["wavelength_nm"][1]),
            "--index-profile",
            str(work / "index.txt"),
            "-o",
            str(native_file),
        ],
        work / "native-sources.log",
    )
    camera_native = propagate_to_simtel(
        work,
        native_file,
        MUON["altitude_m"],
        MUON["wavelength_nm"],
        MUON["seed"],
        # Same profile as fake-muon: adds the air delay and refraction that the IACT time includes.
        profile=REFERENCE / "cfg/common/atmprof1.dat",
    )
    prop_file = work / "native_prop.gst"
    reference_file = work / "reference.iact.gz"
    run_tool(
        fake_muon_argv(str(reference_file), MUON["events"]), work / "fake-muon.log"
    )

    native, reference = read_native(native_file), read_reference(reference_file)
    n_events = MUON["events"]
    dirs_n, dirs_r = unit_vectors(native), unit_vectors(reference)
    wl1, wl2 = MUON["wavelength_nm"]
    metrics, details = {}, {}

    # Yield: native and reference counts from the same geometry; analytic from the index table.
    native_total, reference_total = (
        float(native["w"].sum()),
        float(reference["w"].sum()),
    )
    analytic = muon_expected_hits(table)
    yield_ratio = native_total / reference_total
    yield_se = math.sqrt(1.0 / native_total + 1.0 / reference_total)
    metrics.update(
        {
            "native_photons_per_event": native_total / n_events,
            "reference_photons_per_event": reference_total / n_events,
            "analytic_photons_per_event": analytic,
            "yield_ratio_native_reference": yield_ratio,
            "yield_ratio_se": yield_se,
            "yield_ratio_limit_exceeded": int(
                abs(yield_ratio - 1) > ratio_limit(yield_se)
            ),
            "native_analytic_ratio": (native_total / n_events) / analytic,
            "native_analytic_limit_exceeded": int(
                abs((native_total / n_events) / analytic - 1)
                > ratio_limit(1 / math.sqrt(native_total))
            ),
            "reference_analytic_ratio": (reference_total / n_events) / analytic,
            "reference_analytic_limit_exceeded": int(
                abs((reference_total / n_events) / analytic - 1)
                > ratio_limit(1 / math.sqrt(reference_total))
            ),
        }
    )

    # Cherenkov angle: photon direction against the track, cos(theta_c) = 1/(n beta) at emission height.
    cos_native = -dirs_n[:, 2]
    cos_reference = -dirs_r[:, 2]
    h_native = native["z"]
    h_reference = reference["zem"] - MUON["altitude_m"]
    res_native = cos_native - cherenkov_cos(
        np.interp(h_native, table[:, 0], table[:, 1])
    )
    res_reference = cos_reference - cherenkov_cos(
        np.interp(h_reference, table[:, 0], table[:, 1])
    )
    for label, res in (("native", res_native), ("reference", res_reference)):
        mean, se = mean_se(res)
        flag(metrics, f"cherenkov_cos_residual_{label}_mean", mean, max(1e-6, 5 * se))
        metrics[f"cherenkov_cos_residual_{label}_max_abs"] = float(np.max(np.abs(res)))

    # Wavelength: 1/lambda^2 spectrum in the band (analytic CDF) and native against reference.
    def wl_cdf(x):
        return (1.0 / wl1 - 1.0 / x) / (1.0 / wl1 - 1.0 / wl2)

    ks_flags(
        metrics,
        "wavelength_native_vs_analytic",
        ks_one_sample(native["wavelength"], wl_cdf),
        ks_critical(native["wavelength"].size),
    )
    ks_flags(
        metrics,
        "wavelength_reference_vs_analytic",
        ks_one_sample(reference["wavelength"], wl_cdf),
        ks_critical(reference["wavelength"].size),
    )
    ks_flags(
        metrics,
        "wavelength_native_vs_reference",
        ks_statistic(native["wavelength"], reference["wavelength"]),
        ks_critical(native["wavelength"].size, reference["wavelength"].size),
    )

    # Azimuth around the track (uniform for Cherenkov light) and direction against reference.
    azimuth_n = (np.arctan2(dirs_n[:, 1], dirs_n[:, 0]) + math.pi) / (2 * math.pi)
    azimuth_r = (np.arctan2(dirs_r[:, 1], dirs_r[:, 0]) + math.pi) / (2 * math.pi)
    ks_flags(
        metrics,
        "azimuth_native_vs_analytic",
        ks_one_sample(azimuth_n, lambda x: x),
        ks_critical(azimuth_n.size),
    )
    ks_flags(
        metrics,
        "azimuth_reference_vs_analytic",
        ks_one_sample(azimuth_r, lambda x: x),
        ks_critical(azimuth_r.size),
    )
    # cz spans only ~2e-5 around cos(theta_c) and is stored as float32 (6e-8 resolution): the KS distance is
    # reported but not limited; the mean is limited at the float32 quantum and the 5 se.
    metrics["cz_native_vs_reference_ks"] = ks_statistic(dirs_n[:, 2], dirs_r[:, 2])
    mean_cz_diff = float(np.mean(dirs_n[:, 2]) - np.mean(dirs_r[:, 2]))
    se_cz = math.sqrt(
        np.var(dirs_n[:, 2]) / dirs_n.shape[0] + np.var(dirs_r[:, 2]) / dirs_r.shape[0]
    )
    flag(metrics, "cz_mean_difference", mean_cz_diff, max(6e-8, 5 * se_cz))

    # Positions: IACT x, y are the crossing of each photon with the observation plane z = 0 (checked here:
    # the crossing of the native rays reproduces the reference x, y radii). Native rays start at the emission point.
    crossing = observation_crossing(positions_n(native), dirs_n)
    radial_n = np.hypot(crossing[:, 0], crossing[:, 1])
    radial_r = np.hypot(reference["x"], reference["y"])
    ks_flags(
        metrics,
        "obs_plane_xy_radius_native_vs_reference",
        ks_statistic(radial_n, radial_r),
        ks_critical(radial_n.size, radial_r.size),
    )

    # Time: arrival at the telescope plane z = 0 (IACT time = emission time + optical path to z = 0, air
    # delay included). Native arrival is the propagated stream (--profile); the expectation is analytic.
    t_native = read_native(prop_file)["t"]
    t_expected_native = muon_arrival_expected_ns(native["z"], -dirs_n[:, 2], table)
    t_expected_reference = muon_arrival_expected_ns(
        reference["zem"] - MUON["altitude_m"], -dirs_r[:, 2], table
    )
    metrics["time_native_mean_ns"] = float(np.mean(t_native))
    metrics["time_native_sd_ns"] = float(np.std(t_native))
    metrics["time_reference_mean_ns"] = float(np.mean(reference["t"]))
    metrics["time_reference_sd_ns"] = float(np.std(reference["t"]))
    se_diff = math.sqrt(
        np.var(t_native) / t_native.size + np.var(reference["t"]) / reference["t"].size
    )
    flag(
        metrics,
        "time_mean_difference_ns",
        float(np.mean(t_native) - np.mean(reference["t"])),
        max(0.05, 5 * se_diff),
    )
    # Judged: obdeect against the analytic arrival time (same geometry and n(h)).
    se_expected = math.sqrt(
        np.var(t_native) / t_native.size
        + np.var(t_expected_native) / t_expected_native.size
    )
    flag(
        metrics,
        "time_native_vs_analytic_mean_difference_ns",
        float(np.mean(t_native) - np.mean(t_expected_native)),
        max(0.05, 5 * se_expected),
    )
    ks_flags(
        metrics,
        "time_native_vs_analytic",
        ks_statistic(t_native, t_expected_native),
        ks_critical(t_native.size, t_expected_native.size),
    )
    # Recorded, not judged: sim_telarray against the analytic arrival and native against reference KS.
    metrics["time_reference_vs_analytic_mean_difference_ns"] = float(
        np.mean(reference["t"]) - np.mean(t_expected_reference)
    )
    metrics["time_reference_vs_analytic_ks"] = float(
        ks_statistic(reference["t"], t_expected_reference)
    )
    metrics["time_reference_vs_analytic_ks_critical"] = float(
        ks_critical(reference["t"].size, t_expected_reference.size)
    )
    metrics["time_native_vs_reference_ks"] = float(
        ks_statistic(t_native, reference["t"])
    )
    metrics["time_native_vs_reference_ks_critical"] = float(
        ks_critical(t_native.size, reference["t"].size)
    )
    # Sources stage only (before propagation, no air delay): emission-time diagnostics.
    t_centre = centre_plane_time(native["t"], positions_n(native), dirs_n)
    metrics["time_centre_native_mean_ns"] = float(np.mean(t_centre))
    metrics["time_centre_native_sd_ns"] = float(np.std(t_centre))
    t_obs = native["t"] + (-positions_n(native)[:, 2] / dirs_n[:, 2]) / C_M_NS
    metrics["time_obs_plane_mean_difference_ns"] = float(
        np.mean(t_obs) - np.mean(reference["t"])
    )
    metrics["time_obs_plane_sd_ns"] = float(np.std(t_obs))

    # Camera level: independent sim_telarray runs of both streams.
    camera_native_pe = camera_pe(
        simulate(work, camera_native, "camera-native", extra=CAMERA_MUON)
    )
    camera_reference_pe = camera_pe(
        simulate(work, reference_file, "camera-reference", extra=CAMERA_MUON)
    )
    metrics.update(camera_comparison(camera_native_pe, camera_reference_pe))

    details.update(
        {
            "case": "muon-track versus fake-muon (LightEmission)",
            "mapping": {
                "geometry": "telescope sphere radius 15 m at origin; track from 3000 m above observation level to impact (0,0), vertical",
                "index": "constant index replaced by sim_telarray profile: fake-muon -A 0 only uses the profile; "
                "obdeect gets the same n(h) through --index-profile (index.txt from fake-muon printout)",
                "spectrum": "Cherenkov 1/lambda^2, 300-600 nm in both",
                "photons": "fake-muon -b 1 (single photons); obdeect muon-track emits single photons",
                "positions": "IACT x, y are the observation-plane crossing (z=0) of each ray; native rays crossed at z=0 (compared in obs_plane_xy_radius_native_vs_reference)",
                "time_origin": "emission time -h/c (track time zero at impact, beta=1); arrival at z=0 for both codes",
                "time_convention": "native: propagate --profile atmprof1.dat (air delay and refraction) and read t; "
                "reference: IACT time (raybnd: air delay included) + event time_offset; expected: "
                "t_em + (h + int_0^h (n-1) dh)/(|cz| c), straight path, n(h) from index.txt",
                "time_reference_finding": "sim_telarray IACT time is 0.005 ns above the straight-path expectation "
                "on average and differs at KS 0.02 (see time_reference_vs_analytic_*); obdeect matches the "
                "expectation (time_native_vs_analytic_*). Not judged: the reference residual is recorded only.",
                "time_earlier_failure": "the harness propagated without a profile (no air delay, 0.26 ns mean) and "
                "compared the sphere-centre plane; both were harness errors, not obdeect errors",
                "sim_telarray": "propagate --no-extinction --profile atmprof1.dat then obdeect stream to-simtel; both with the same chain settings, telescope_theta=0, trigger_telescopes=1",
            },
            "physics_limits": {
                "yield": "max(0.01, 5 se) (task rule)",
                "cherenkov_cos": "max(1e-6, 5 se): float32 directions and 5e-8 index table precision",
                "time_mean": "max(0.05 ns, 5 se): float32 timestamps",
                "ks": "two-sample / one-sample KS at alpha=1e-3 (statistical)",
                "camera": "max(0.01, 5 se) on the relative p.e. per event ratio (task rule)",
            },
            "limit_setting": "Set before the first run: yield, Cherenkov angle, wavelength, azimuth, camera ratio. "
            "Set after exploratory runs: the cz check (KS replaced by a mean check, cz is float32 at "
            "6e-8 resolution against a 2e-5 spread) and the position convention (observation plane). "
            "Time: judged against the analytic arrival; the sim_telarray-vs-analytic residual is recorded, not limited.",
            "analytic": "numerical integral of dN/ds x (fraction of Cherenkov azimuths hitting the sphere) over the track",
            "reproducibility": "fake-muon has no seed option: the reference sample changes between runs at the 1e-4 level in yield; the native run is seeded (17)",
            "index_table": table.tolist(),
            "camera_pe_native": camera_native_pe.tolist(),
            "camera_pe_reference": camera_reference_pe.tolist(),
        }
    )

    def plot(axes):
        from obdeect_validation.adapters import plotting as pl

        span = np.linspace(wl1, wl2, 200)
        axes[0, 0].plot(
            span,
            1.0 / span**2 / (1.0 / wl1 - 1.0 / wl2),
            color=pl.ANALYTIC,
            linestyle="--",
            label="1/lambda^2",
        )
        pl.histograms(
            axes[0, 0],
            [
                ("obdeect", native["wavelength"], pl.NATIVE),
                ("fake-muon", reference["wavelength"], pl.REFERENCE),
            ],
            "wavelength [nm]",
            clip=None,
        )
        pl.histograms(
            axes[0, 1],
            [
                ("obdeect", res_native, pl.NATIVE),
                ("fake-muon", res_reference, pl.REFERENCE),
            ],
            "cos(angle to track) - 1/(n beta)",
        )
        pl.histograms(
            axes[1, 0],
            [
                ("obdeect, z=0", t_native, pl.NATIVE),
                ("analytic, z=0", t_expected_native, pl.ANALYTIC),
                ("fake-muon", reference["t"], pl.REFERENCE),
            ],
            "time [ns]",
            bins=80,
        )
        pl.camera_bars(axes[1, 1], camera_native_pe, camera_reference_pe, "fake-muon")

    write_outputs(metrics, details, plot, Path("muon-track.png"))
    return metrics


CAMERA_MUON = ["-C", "telescope_theta=0", "-C", "trigger_telescopes=1"]
# The LST discriminator (99999) never fires for the 1 ns night-sky pulse; a zero threshold stores every
# event so that p.e. per event can be compared. p.e. do not depend on the threshold.
CAMERA_NIGHT = [
    "-C",
    "telescope_theta=0",
    "-C",
    "trigger_telescopes=1",
    "-C",
    "discriminator_threshold=0",
]


# ---------------------------------------------------------------- night-sky


def night_sky(work):
    cfg_dir = REFERENCE / "cfg/common"
    window_ns = NIGHT["pulse_fwhm_ns"] / FWHM_TO_SIGMA
    photons = round(
        night_sky_photons(
            NIGHT["rate_per_sr_ns_cm2"],
            NIGHT["viewcone_deg"],
            NIGHT["radius_m"],
            window_ns,
        )
    )
    native_file = work / "native_sources.gst"
    run_tool(
        [
            tool("obdeect"),
            "sources",
            "run",
            "night-sky",
            "--observation-level",
            str(NIGHT["altitude_m"]),
            "--atmosphere-id",
            "1",
            "--seed",
            str(NIGHT["seed"]),
            "--events",
            str(NIGHT["events"]),
            "--photons",
            str(photons),
            "--spectrum",
            f"fixed:wavelength_nm={NIGHT['spectrum_nm']}",
            "--pulse",
            NIGHT_PULSE,
            "--bunch-size",
            str(NIGHT["bunch_size"]),
            "--viewcone-deg",
            str(NIGHT["viewcone_deg"]),
            "--telescope-position",
            "0,0,0",
            "--telescope-radius",
            str(NIGHT["radius_m"]),
            "-o",
            str(native_file),
        ],
        work / "native-sources.log",
    )
    reference_file = work / "reference.iact.gz"
    run_tool(
        [
            tool("nsbls"),
            f"-I{cfg_dir}",
            "--events",
            str(NIGHT["events"]),
            "--photons",
            str(photons),
            "--spectrum",
            str(NIGHT["spectrum_nm"]),
            "--pulse-width",
            str(NIGHT["pulse_fwhm_ns"]),
            "-b",
            str(NIGHT["bunch_size"]),
            "-r",
            str(NIGHT["radius_m"]),
            "-v",
            str(NIGHT["viewcone_deg"]),
            "--altitude",
            str(NIGHT["altitude_m"]),
            "-A",
            "1",
            "--run",
            "1",
            "-o",
            str(reference_file),
        ],
        work / "nsbls.log",
    )
    camera_native = propagate_to_simtel(
        work, native_file, NIGHT["altitude_m"], NIGHT["band_nm"], NIGHT["seed"]
    )

    native, reference = read_native(native_file), read_reference(reference_file)
    n_events = NIGHT["events"]
    dirs_n, dirs_r = unit_vectors(native), unit_vectors(reference)
    radius, viewcone = NIGHT["radius_m"], math.radians(NIGHT["viewcone_deg"])
    omega = 2.0 * math.pi * (1.0 - math.cos(viewcone))
    area_cm2 = math.pi * (100.0 * radius) ** 2
    metrics, details = {}, {}

    # Yield and configured rate.
    native_total, reference_total = (
        float(native["w"].sum()),
        float(reference["w"].sum()),
    )
    expected = float(photons * n_events)
    rate_native = native_total / n_events / (omega * area_cm2 * window_ns)
    rate_reference = reference_total / n_events / (omega * area_cm2 * window_ns)
    metrics.update(
        {
            "photons_per_event_configured": photons,
            "native_photons_per_event": native_total / n_events,
            "reference_photons_per_event": reference_total / n_events,
            "yield_ratio_native_reference": native_total / reference_total,
            "yield_ratio_native_reference_limit_exceeded": int(
                abs(native_total / reference_total - 1) > ratio_limit(0.0)
            ),
            "rate_native_per_sr_ns_cm2": rate_native,
            "rate_reference_per_sr_ns_cm2": rate_reference,
            "rate_configured_per_sr_ns_cm2": NIGHT["rate_per_sr_ns_cm2"],
        }
    )
    flag(
        metrics,
        "rate_native_relative_error",
        rate_native / NIGHT["rate_per_sr_ns_cm2"] - 1,
        ratio_limit(0.0),
    )
    flag(
        metrics,
        "rate_reference_relative_error",
        rate_reference / NIGHT["rate_per_sr_ns_cm2"] - 1,
        ratio_limit(0.0),
    )
    flag(
        metrics,
        "yield_native_expected_relative_error",
        native_total / expected - 1,
        ratio_limit(0.0),
    )
    flag(
        metrics,
        "yield_reference_expected_relative_error",
        reference_total / expected - 1,
        ratio_limit(0.0),
    )

    # Direction: uniform in solid angle inside the cone (cos of angle to the axis uniform in [cos v, 1]).
    cos_n = -dirs_n[:, 2]
    cos_r = -dirs_r[:, 2]
    uniform_n = (1.0 - cos_n) / (1.0 - math.cos(viewcone))
    uniform_r = (1.0 - cos_r) / (1.0 - math.cos(viewcone))
    ks_flags(
        metrics,
        "cone_uniform_native_vs_analytic",
        ks_one_sample(uniform_n, lambda x: x),
        ks_critical(uniform_n.size),
    )
    ks_flags(
        metrics,
        "cone_uniform_reference_vs_analytic",
        ks_one_sample(uniform_r, lambda x: x),
        ks_critical(uniform_r.size),
    )
    ks_flags(
        metrics,
        "cz_native_vs_reference",
        ks_statistic(dirs_n[:, 2], dirs_r[:, 2]),
        ks_critical(dirs_n.shape[0], dirs_r.shape[0]),
    )

    # Impact density on the projected sphere: an isotropic flux is uniform in projected area,
    # i.e. rho^2 = sin^2(angle between entry point and -d) is uniform in [0, 1].
    # One definition for both: closest-approach distance of each ray line to the origin, |P x d|.
    # Native points are the sphere entry points; reference points are the observation-plane
    # crossings (x, y, 0). Both are points on the same ray, so the value does not depend on the point.
    entry_n = np.stack([native["x"], native["y"], native["z"]], 1)
    entry_r = np.stack(
        [reference["x"], reference["y"], np.zeros(reference["x"].size)], 1
    )
    b_n = closest_approach(entry_n, dirs_n)
    b_r = closest_approach(entry_r, dirs_r)
    rho2_n, rho2_r = (b_n / radius) ** 2, (b_r / radius) ** 2
    ks_flags(
        metrics,
        "projected_area_uniform_native_vs_analytic",
        ks_one_sample(np.clip(rho2_n, 0, 1), lambda x: x),
        ks_critical(rho2_n.size),
    )
    ks_flags(
        metrics,
        "projected_area_uniform_reference_vs_analytic",
        ks_one_sample(np.clip(rho2_r, 0, 1), lambda x: x),
        ks_critical(rho2_r.size),
    )
    ks_flags(
        metrics,
        "impact_parameter_native_vs_reference",
        ks_statistic(b_n, b_r),
        ks_critical(b_n.size, b_r.size),
    )
    # Analytic CDF of b/R = ri for an area-uniform disc (ri^2 uniform): P(ri <= r) = r^2.
    ks_flags(
        metrics,
        "impact_parameter_native_vs_sampler",
        ks_one_sample(np.clip(b_n / radius, 0, 1), lambda r: r * r),
        ks_critical(b_n.size),
    )

    # Time: arrival at the sphere-centre plane (perpendicular to the ray through the origin) for both,
    # t - (X . d)/c with X the sample point. The pulse is uniform over window_ns there.
    # Reference time is the raw IACT time (event offset included), i.e. at the observation plane.
    t_centre = centre_plane_time(native["t"], entry_n, dirs_n)
    t_ref_obs = reference["t"]
    t_ref_centre = centre_plane_time(t_ref_obs, entry_r, dirs_r)

    def pulse_cdf(x):
        return np.clip(x / window_ns + 0.5, 0.0, 1.0)

    ks_flags(
        metrics,
        "time_native_vs_pulse",
        ks_one_sample(t_centre, pulse_cdf),
        ks_critical(t_centre.size),
    )
    ks_flags(
        metrics,
        "time_reference_vs_pulse",
        ks_one_sample(t_ref_centre, pulse_cdf),
        ks_critical(t_ref_centre.size),
    )
    ks_flags(
        metrics,
        "time_native_vs_reference",
        ks_statistic(t_centre, t_ref_centre),
        ks_critical(t_centre.size, t_ref_centre.size),
    )
    mean_diff, se_diff = time_difference(t_centre, t_ref_centre)
    flag(metrics, "time_mean_difference_ns", mean_diff, max(0.05, 5 * se_diff))
    metrics["time_reference_sd_ns"] = float(np.std(t_ref_centre))
    metrics["time_reference_observation_plane_sd_ns"] = float(np.std(t_ref_obs))
    metrics["time_native_sd_ns"] = float(np.std(t_centre))

    # Wavelength is fixed at the configured value in both.
    metrics["wavelength_mismatch_count"] = int(
        np.sum(np.abs(native["wavelength"] - NIGHT["spectrum_nm"]) > 1e-3)
        + np.sum(np.abs(reference["wavelength"] - NIGHT["spectrum_nm"]) > 1e-3)
    )

    # Camera level: at the configured rate (10x help-text rate, same on both sides). A sample with fewer
    # than two triggered events is recorded as BLOCKED, not as a pass.
    camera_native_pe = camera_pe(
        simulate(work, camera_native, "camera-native", extra=CAMERA_NIGHT)
    )
    camera_reference_pe = camera_pe(
        simulate(work, reference_file, "camera-reference", extra=CAMERA_NIGHT)
    )
    camera_status = "run"
    try:
        metrics.update(camera_comparison(camera_native_pe, camera_reference_pe))
    except ValueError as exc:
        camera_status = f"BLOCKED: {exc}; sim_telarray triggered events: native {camera_native_pe.size}, reference {camera_reference_pe.size}"

    details.update(
        {
            "case": "night-sky versus nsbls (LightEmission)",
            "impact_parameter": "closest-approach distance of the ray line to the origin, |P x d| for unit d. "
            "Native: point on the sphere (entry point). Reference: observation-plane crossing (x, y, 0). Same ray, same value.",
            "mapping": {
                "geometry": "telescope sphere radius 15 m at origin, viewing cone 10 deg about (0,0,-1), altitude 2150 m (both)",
                "spectrum": "fixed 400 nm (both); inside the 300-650 nm band of the configured rate",
                "rate": "configured 2.4 ph/(sr ns cm^2) (10 x the nsbls help-text rate 0.24, raised so that the camera triggers; "
                "same on both sides) -> photons per event = rate x Omega x pi R^2 x window",
                "window": "uniform pulse of FWHM 2.35482 ns, i.e. width 1 ns (both)",
                "bunch": "bunch size 2 photons (both)",
                "time_convention": "both at the sphere-centre plane (through the origin, perpendicular to the ray): "
                "t - (X.d)/c with X the sample point (native: sphere point; reference: observation-plane crossing "
                "(x, y, 0) with the IACT time including the event time_offset); c = vacuum speed (native uses c too)",
                "atmosphere": "nsbls -A 1 and obdeect -A 1 (recorded; profile for air speed only)",
                "sim_telarray": "propagate --no-extinction then obdeect stream to-simtel; telescope_theta=0, trigger_telescopes=1",
            },
            "physics_limits": {
                "yield_and_rate": "max(0.01, 5 se) with se=0: fixed photon count per event in both",
                "cone_and_direction": "KS alpha=1e-3",
                "projected_area": "for a fixed direction an isotropic flux is uniform over the disc perpendicular to d, "
                "so b^2 is uniform on [0, R^2] (checked with a numpy toy MC: KS 0.002 at 2e5 rays); KS alpha=1e-3. "
                "Impact parameter native vs sampler: analytic CDF of b/R is r^2",
                "time_mean": "max(0.05 ns, 5 se)",
                "camera": "max(0.01, 5 se) on the relative p.e. per event ratio (task rule)",
            },
            "reproducibility": "nsbls has no seed option: the reference sample is not reproducible; the native run is seeded (17)",
            "limit_setting": "Limits set before the first run and not changed. Same definitions for native and reference "
            "(closest-approach distance, centre-plane time). The projected-area expectation is the physical one (b^2 uniform), "
            "which is the obdeect sampler after the 2026-11-14 fix (ri = sqrt(U), docs/sources/spec/background.md).",
            "judged_and_recorded": "Judged: yield, rate, cone, cz, camera p.e. ratio, projected area and time of obdeect "
            "against physics, impact parameter of obdeect against the sampler. Recorded only (values in metrics.json, "
            "no limit): impact_parameter_native_vs_reference, projected_area_uniform_reference_vs_analytic, "
            "time_reference_vs_pulse, time_native_vs_reference, time_mean_difference_ns. The nsbls values are not "
            "judged because nsbls has the two issues in reference_findings.",
            "reference_findings": [
                (
                    "(a) nsbls samples ri = 1 - sqrt(U) (copied from sim_telarray nsbls; the same sampler as obdeect "
                    "before 2026-11-14). An isotropic flux needs ri = sqrt(U): with 1 - sqrt(U), b^2 is not uniform "
                    "(toy MC KS 0.50), so projected_area_uniform_reference_vs_analytic fails and the impact-parameter "
                    "distribution is wrong. Recorded only."
                ),
                (
                    "(b) nsbls (IactLightEmission.cc ~l.1395) rotates the sphere point with az = d.Azimuth() instead of "
                    "the azimuth of -d. The z axis then maps to (cx, cy, -cz), not to -d, so about 0.03 % of the entry "
                    "points are on the wrong side (away from the source). Impact parameters and centre-plane arrival "
                    "times are then wrong: the dt time correction is wrong. A numpy replica reproduces the nsbls "
                    "impact-parameter quantiles and the centre-plane time SD (3.6 ns); with az = atan2(-d.y, -d.x) "
                    "the SD is 0. Recorded only."
                ),
            ],
            "camera_status_note": "2026-11-13 run at 1x rate: 0 triggered events in both samples (BLOCKED)",
            "camera_status": camera_status,
            "camera_pe_native": camera_native_pe.tolist(),
            "camera_pe_reference": camera_reference_pe.tolist(),
        }
    )

    def plot(axes):
        from obdeect_validation.adapters import plotting as pl

        sides = (("obdeect", pl.NATIVE), ("nsbls", pl.REFERENCE))
        pl.histograms(
            axes[0, 0],
            [
                (sides[0][0], uniform_n, pl.NATIVE),
                (sides[1][0], uniform_r, pl.REFERENCE),
            ],
            "(1-cos theta)/(1-cos v)",
            bins=50,
            clip=None,
        )
        pl.histograms(
            axes[0, 1],
            [(sides[0][0], rho2_n, pl.NATIVE), (sides[1][0], rho2_r, pl.REFERENCE)],
            "rho^2 (projected area fraction)",
            bins=50,
            clip=None,
        )
        axes[0, 1].axhline(1.0, color=pl.ANALYTIC, linestyle="--", label="uniform")
        axes[0, 1].legend()
        pl.histograms(
            axes[1, 0],
            [
                ("obdeect, centre plane", t_centre, pl.NATIVE),
                ("nsbls, centre plane", t_ref_centre, pl.REFERENCE),
            ],
            "time [ns]",
        )
        pl.camera_bars(axes[1, 1], camera_native_pe, camera_reference_pe, "nsbls")

    write_outputs(metrics, details, plot, Path("night-sky.png"))
    return metrics


# ---------------------------------------------------------------- laser-beam (geometry only)


def _ls_beam_point(label, point, text):
    """Position [m] ls-beam prints for the beam crossing the field of view."""
    match = re.search(label + point, text)
    if match is None:
        raise RuntimeError(
            f"ls-beam printed no field-of-view point for {label.strip()!r}"
        )
    return np.array([float(v) for v in match.groups()]) / 100.0


def laser_beam(work):
    cfg_dir = REFERENCE / "cfg/common"
    native_file, reference_file = (
        work / "native_sources.gst",
        work / "reference.iact.gz",
    )
    run_tool(
        [
            tool("obdeect"),
            "sources",
            "run",
            "laser-beam",
            "--observation-level",
            str(LASER["altitude_m"]),
            "--atmosphere-id",
            "1",
            "--events",
            "1",
            "--seed",
            str(LASER["seed"]),
            "--bunch-size",
            "3",
            "--photons-per-m",
            str(LASER["primary_per_m"]),
            "--spectrum",
            f"fixed:wavelength_nm={LASER['spectrum_nm']}",
            "--laser-position",
            ",".join(str(v) for v in LASER["laser_position_m"]),
            "--laser-zenith",
            str(LASER["laser_theta_deg"]),
            "--laser-azimuth",
            str(LASER["laser_phi_deg"]),
            "--telescope-zenith",
            str(LASER["telescope_theta_deg"]),
            "--telescope-azimuth",
            str(LASER["telescope_phi_deg"]),
            "--telescope-radius",
            str(LASER["telescope_radius_m"]),
            "--fov-deg",
            str(LASER["fov_deg"]),
            "-o",
            str(native_file),
        ],
        work / "native-sources.log",
    )
    ls_text = run_tool(
        [
            tool("ls-beam"),
            f"-I{cfg_dir}",
            "--events",
            "1",
            "--brightness",
            str(LASER["brightness_per_m"]),
            "--bunchsize",
            "1",
            "--step",
            str(LASER["step_m"]),
            "--spectrum",
            str(LASER["spectrum_nm"]),
            "--laser-position",
            ",".join(str(v) for v in LASER["laser_position_m"]),
            "--laser-theta",
            str(LASER["laser_theta_deg"]),
            "--laser-phi",
            str(LASER["laser_phi_deg"]),
            "--telescope-theta",
            str(LASER["telescope_theta_deg"]),
            "--telescope-phi",
            str(LASER["telescope_phi_deg"]),
            "--telescope-radius",
            str(LASER["telescope_radius_m"]),
            "--fov",
            str(LASER["fov_deg"]),
            "-o",
            str(reference_file),
        ],
        work / "ls-beam.log",
    )
    point = r"SpaceVect: x=([-0-9.e+]+), y=([-0-9.e+]+), z=([-0-9.e+]+) cm"
    # "Laser exits enters FoV" is ls-beam's own wording.
    ls_entry = _ls_beam_point(r"beam enters FoV at \{ ", point, ls_text)
    ls_exit = _ls_beam_point(r"Laser exits enters FoV at \{ ", point, ls_text)

    native = read_native(native_file)
    positions = positions_n(native)
    direction = unit_vectors(native)
    metrics, details = {}, {}
    order = np.argsort(positions[:, 2])
    native_lo, native_hi = positions[order[0]], positions[order[-1]]
    # ls-beam entry is the lower point (beam goes upward), exit the upper one.
    lo_ref, hi_ref = (
        (ls_entry, ls_exit) if ls_entry[2] < ls_exit[2] else (ls_exit, ls_entry)
    )
    endpoint_lo = float(np.linalg.norm(native_lo - lo_ref))
    endpoint_hi = float(np.linalg.norm(native_hi - hi_ref))
    flag(metrics, "fov_entry_offset_m", endpoint_lo, LASER["step_m"])
    flag(metrics, "fov_exit_offset_m", endpoint_hi, LASER["step_m"])
    metrics["fov_path_native_m"] = float(np.linalg.norm(native_hi - native_lo))
    metrics["fov_path_reference_m"] = float(np.linalg.norm(hi_ref - lo_ref))
    flag(
        metrics,
        "fov_path_offset_m",
        metrics["fov_path_native_m"] - metrics["fov_path_reference_m"],
        LASER["step_m"],
    )

    # The primaries form one straight beam: all share one direction and lie on the same line.
    beam_dir = direction[0]
    metrics["primary_direction_spread"] = float(
        np.max(np.linalg.norm(direction - beam_dir, axis=1))
    )
    off_line = np.linalg.norm(np.cross(positions - native_lo, beam_dir), axis=1)
    metrics["primary_off_line_max_m"] = float(np.max(off_line))
    # ls-beam entry and exit lie on the primary line (distance does not depend on the 0.1 m sampling).
    ref_off_line = np.linalg.norm(
        np.cross(np.stack([lo_ref, hi_ref]) - native_lo, beam_dir), axis=1
    )
    flag(
        metrics,
        "ls_beam_fov_points_off_primary_line_max_m",
        float(np.max(ref_off_line)),
        1e-2,
    )
    metrics["primary_count"] = int(native["w"].size)

    # ls-beam emission altitudes lie on the same beam segment.
    zem = read_reference(reference_file)["zem"] - LASER["altitude_m"]
    lo_z, hi_z = min(lo_ref[2], hi_ref[2]), max(lo_ref[2], hi_ref[2])
    metrics["ls_beam_emission_count"] = int(zem.size)
    metrics["ls_beam_outside_segment_count"] = int(
        np.sum((zem < lo_z - 1e-3) | (zem > hi_z + 1e-3))
    )
    details.update(
        {
            "case": "laser-beam versus ls-beam: geometry only (BLOCKED for yield, direction, time, camera)",
            "blocked": {
                "reason": "obdeect laser-beam emits primary beam photons; ls-beam emits light scattered from the beam "
                "(Rayleigh/aerosol). obdeect propagation has no scattering yet (plan step 6b: "
                "beta_s(355 nm) aerosol, single-scattering albedo and phase asymmetry missing), so scattered "
                "yield, direction, wavelength, time and camera p.e. cannot be compared.",
            },
            "mapping": {
                "geometry": "laser (1000,0,0) m, theta 60, phi 180 (astronomical), telescope theta 30 phi 0, "
                "radius 12.5 m, FoV 15 deg, step 0.1 m (both)",
                "spectrum": "355 nm (both; not used in the geometry check)",
                "compared": "FoV entry and exit printed by ls-beam: endpoint offsets to the primary extremes (one step) and "
                "perpendicular distance to the primary line (1 cm); emission altitudes inside the primary segment",
            },
            "physics_limits": {
                "endpoints": "one step (0.1 m): primary sampling at the ends",
                "off_line": "1 cm: ls-beam writes cm with float precision",
            },
            "limit_setting": "Limits set before the first run.",
            "ls_beam_entry_m": ls_entry.tolist(),
            "ls_beam_exit_m": ls_exit.tolist(),
            "native_endpoints_m": [native_lo.tolist(), native_hi.tolist()],
        }
    )

    def plot(axes):
        from obdeect_validation.adapters import plotting as pl

        axes[0, 0].scatter(
            positions[:, 0],
            positions[:, 2],
            s=3,
            color=pl.NATIVE,
            label="obdeect primaries",
        )
        axes[0, 0].scatter(
            [lo_ref[0], hi_ref[0]],
            [lo_ref[2], hi_ref[2]],
            marker="x",
            s=60,
            color=pl.REFERENCE,
            label="ls-beam FoV",
        )
        axes[0, 0].set(xlabel="x [m]", ylabel="z [m]")
        axes[0, 0].legend()
        pl.histograms(
            axes[0, 1],
            [
                ("ls-beam emission altitude", zem, pl.REFERENCE),
                ("obdeect primaries", positions[:, 2], pl.NATIVE),
            ],
            "z [m]",
            density=False,
            clip=None,
        )

    write_outputs(metrics, details, plot, Path("laser-beam.png"), rows=1)
    return metrics


# ---------------------------------------------------------------- entry point


def main(case):
    work = Path.cwd().resolve()
    if case == "muon-track":
        muon_track(work)
    elif case == "night-sky":
        night_sky(work)
    elif case == "laser-beam":
        laser_beam(work)
    else:
        raise ValueError(f"unknown source case: {case}")


if __name__ == "__main__":
    main(sys.argv[1])
