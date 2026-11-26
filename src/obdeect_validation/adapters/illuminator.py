"""Illuminators ILLN-01 and ILLS-01 versus LightEmission xyzls: emission and camera stage.

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.illuminator ILLN-01

Geometry as simtools builds it (simtel/simulator_light_emission.py): the illuminator sits at its
ground position plus illuminator_tower_height; the telescope is the LST of the same site at its
ground position; the beam points from the illuminator to the telescope; the telescope zenith and
azimuth follow _calibration_pointing_direction. Each illuminator runs at three distances along
the line of sight: nominal, 2x and 0.5x (illuminator moved on the line to the telescope).

Reference: xyzls with the arguments simtools passes (-x -y -z in cm, -d pointing, -n photons per
event, -N events, -s wavelength, -p pulse, -a angular distribution, -h observation level, -A
atmosphere). The telescope is at the origin (position relative to the telescope, radius -r
TELESCOPE_RADIUS_CM). Native: obdeect sources point-source with the same arguments. The illuminator design values (ILLN-design / ILLS-design: 300 nm, Gauss pulse
of width 0, angular distribution gauss of width 0, which is isotropic) are fixed here and checked
against the model files in the metrics.

Atmosphere: both sides use atmprof1 (the xyzls -A 1 profile). The native stream is propagated by
obdeect with that profile (--profile: bent rays, air delay), and this propagated stream is used for
stage 1 (arrival time at the plane z = 0 as obdeect gives it, photons landing within the radius)
and for stage 2. A straight-line air delay is not used: at the South geometry (zenith about 87 deg)
it differs from the bent-ray travel time by about 0.8 ns, and the bent-ray time agrees with an
independent numerical integral of the Snell's-law ray (plane-parallel, atmprof1).

Stage 1 compares the two photon streams as source_camera does (photons per event, direction,
arrival time, wavelength). Stage 2 converts the propagated stream with obdeect stream to-simtel
and runs both streams through sim_telarray with the telescope pointed at the illuminator
(telescope_theta and telescope_phi from the geometry, Bypass_Optics=0). The
sim_telarray configuration is the LST-1 North one for both sites: only position and pointing
differ between the sites. The two camera samples are independent (no pairing).

Limits (set before the first run, not changed to pass):
  stage 1 and camera p.e. ratio: as in source_camera (photon ratio max(0.5 %, 5 se); mean
  direction 5 sigma; KS shape: fail when the two-sample p-value is below 1e-4 (source_camera _ks_flag); mean time max(0.05 ns, 5 se); wavelength
  max(0.01 nm, 5 se); camera p.e. ratio max(1 %, 5 se); image shape reduced chi2 <= 2).
  Per-pixel p.e. map: pixels with at least PIXEL_MIN_MEAN_PE p.e. per event in the reference;
  R_i = mean native / mean reference with the event-to-event standard errors of both means;
  reduced chi2 of (R_i - 1)/sigma_i must be <= 2.

sim_telarray stores every event (discriminator threshold 0, see TRIGGER). If sim_telarray gives fewer than two triggered events for either input, stage 2 is BLOCKED and
no camera metric is written.
"""

from __future__ import annotations

import json
import math
import traceback
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common, source_camera
from obdeect_validation.adapters import obdeect_cli as cli

ILLUMINATORS = {
    "ILLN-01": {"site": "North", "telescope": "LSTN-01"},
    "ILLS-01": {"site": "South", "telescope": "LSTS-01"},
}
DISTANCE_SCALES = {"nominal": 1.0, "x2": 2.0, "half": 0.5}
# 200 events: at 2x distance the illuminator gives about 60 p.e. per event, and the trigger
# (TRIGGER) fires in about one event in ten; 20 events left the reference with no triggered event.
EVENTS = 200
PHOTONS_PER_EVENT = 100_000_000
SEED = 17
WAVELENGTH_NM = 300
PULSE = "Gauss:0.0"
ANGULAR = "isotropic"  # simtools: gauss with width 0 is isotropic
OBS_LEVEL_M = 2150.0
TELESCOPE_RADIUS_CM = (
    1250.0  # telescope_sphere_radius 12.5 m (LSTN-design, LSTS-design)
)
PIXEL_MIN_MEAN_PE = 1.0  # per event, reference
PIXEL_CHI2_LIMIT = 2.0
# The sim_telarray analog-sum trigger writes only the triggered events (the .simtel file holds no
# event without a trigger), so at 2x distance (about 60 p.e. per event) the pixel map was built from
# 9 and 3 events. The forced trigger (FORCE_FAKE_TRIGGER=1, FAKE_TRIGGER=0) makes every event
# count as triggered; the p.e. list is written for every telescope event either way. The trigger
# is not part of the photon-to-p.e. model under test, and intensity is unchanged.
TRIGGER = {
    "discriminator_threshold": 0,
    "TRIGGER_TELESCOPES": 1,
    "FAKE_TRIGGER": 0,  # forced trigger (sim_telarray common/sim_telarray.c): every telescope
    "FORCE_FAKE_TRIGGER": 1,  # event counts as triggered, so every event and its p.e. list is written
    "PULSE_ANALYSIS": -30,
    "MAXIMUM_TELESCOPES": 1,
    "Bypass_Optics": 0,
}


# ---------------------------------------------------------------- geometry


def _model_file(instrument, parameter):
    """Model parameter JSON of an instrument (latest version of the file)."""
    files = sorted(
        (common.models_root() / "model_parameters" / instrument / parameter).glob(
            "*.json"
        )
    )
    if not files:
        raise FileNotFoundError(f"no model file for {instrument} {parameter}")
    return json.loads(files[-1].read_text(encoding="utf-8"))


def _ground_position(instrument, parameter="array_element_position_ground"):
    """Ground position [m] of an array element or calibration device."""
    return np.asarray(_model_file(instrument, parameter)["value"], dtype=float)


def illuminator_position(instrument):
    """Illuminator emission position [m]: ground position plus tower height (simtools)."""
    tower = float(_model_file(instrument, "illuminator_tower_height")["value"])
    return _ground_position(instrument) + np.array([0.0, 0.0, tower])


def pointing(illuminator_m, telescope_m):
    """Pointing as simtools computes it: unit vector from illuminator to telescope (light
    travels that way) and the telescope zenith and azimuth in degrees (simtools rounds to 6)."""
    direction = np.asarray(telescope_m, dtype=float) - np.asarray(
        illuminator_m, dtype=float
    )
    norm = float(np.linalg.norm(direction))
    unit = direction / norm
    theta = 180.0 - math.degrees(math.acos(direction[2] / norm))
    phi = 180.0 - math.degrees(math.atan2(direction[1], direction[0]))
    return unit, theta, phi


def scaled_illuminator(illuminator_m, telescope_m, scale):
    """Illuminator moved along the line of sight to the telescope: T + scale (P - T) [m]."""
    telescope_m = np.asarray(telescope_m, dtype=float)
    return telescope_m + scale * (np.asarray(illuminator_m, dtype=float) - telescope_m)


def geometry(instrument, scale):
    """Dict with the distance, pointing and xyzls/point-source arguments of one configuration."""
    spec = ILLUMINATORS[instrument]
    telescope = _ground_position(spec["telescope"])
    nominal = illuminator_position(instrument)
    source = scaled_illuminator(nominal, telescope, scale)
    unit, theta, phi = pointing(source, telescope)
    return {
        "illuminator_m": source,
        "telescope_m": telescope,
        "distance_m": float(np.linalg.norm(telescope - source)),
        "pointing": unit,
        "telescope_theta_deg": theta,
        "telescope_phi_deg": phi,
        # Telescope at the origin (as the point-source case): source relative to the telescope, cm.
        "source_cm": (source - telescope) * 100.0,
    }


def source_argv(geo, events, photons):
    """Arguments shared by xyzls and obdeect point-source (cm, simtools values)."""
    x, y, z = geo["source_cm"]
    return [
        "-R",
        "1",
        "-r",
        f"{TELESCOPE_RADIUS_CM:.1f}",
        "-N",
        str(events),
        "-x",
        f"{x:.6f}",
        "-y",
        f"{y:.6f}",
        "-z",
        f"{z:.6f}",
        "-d",
        ",".join(f"{v:.8f}" for v in geo["pointing"]),
        "-n",
        str(photons),
        "-b",
        "1",
        "-s",
        str(WAVELENGTH_NM),
        "-p",
        PULSE,
        "-a",
        ANGULAR,
        "-h",
        str(OBS_LEVEL_M),
        "-A",
        "1",
    ]


# ---------------------------------------------------------------- pixel map


def camera_event_images(path):
    """Per telescope: array (events, pixels) of photoelectrons of a sim_telarray file."""
    from eventio import SimTelFile

    images = {}
    with SimTelFile(str(path)) as stream:
        for event in stream:
            for tel, pe in event.get("photoelectrons", {}).items():
                images.setdefault(tel, []).append(
                    np.asarray(pe["photoelectrons"], dtype=float)
                )
    return {tel: np.vstack(values) for tel, values in images.items()}


def pixel_map_metrics(reference, candidate):
    """Per-pixel p.e. ratio of the mean image per event with the event-to-event errors.

    reference, candidate: arrays (events, pixels). Pixels with at least PIXEL_MIN_MEAN_PE mean
    p.e. per event in the reference. R_i = mean_cand_i / mean_ref_i, sigma_i from the standard
    error of both means (sample variance over events; this includes the noise of the sim_telarray
    runs). Returns the reduced chi2 of (R_i - 1) / sigma_i, the pixel count and the number of pixels
    beyond 5 sigma.
    """
    n_ref, n_cand = reference.shape[0], candidate.shape[0]
    if min(n_ref, n_cand) < 2:
        raise ValueError("pixel map needs at least two events per sample")
    mean_r, mean_c = reference.mean(axis=0), candidate.mean(axis=0)
    var_r, var_c = reference.var(axis=0, ddof=1), candidate.var(axis=0, ddof=1)
    ok = (mean_r >= PIXEL_MIN_MEAN_PE) & (mean_c > 0)
    if not ok.any():
        raise ValueError("no pixel with enough p.e. in the reference sample")
    ratio = mean_c[ok] / mean_r[ok]
    sigma = ratio * np.sqrt(
        var_c[ok] / (n_cand * mean_c[ok] ** 2) + var_r[ok] / (n_ref * mean_r[ok] ** 2)
    )
    good = sigma > 0
    if not good.any():
        raise ValueError("pixel map: no pixel with a positive error")
    z = (ratio[good] - 1.0) / sigma[good]
    chi2 = float(np.mean(z**2))
    return {
        "pixel_map_pixels_used": int(z.size),
        "pixel_map_reduced_chi2": chi2,
        "pixel_map_reduced_chi2_limit_exceeded": int(chi2 > PIXEL_CHI2_LIMIT),
        "pixel_map_pixels_beyond_5sigma": int(np.sum(np.abs(z) > 5.0)),
        "pixel_map_ratio_median": float(np.median(ratio[good])),
    }


# ---------------------------------------------------------------- timing


# Pass/fail timing: obdeect against the exact bent-ray integral, per photon (plane-parallel Snell's
# law with the atmprof1 profile). Recorded only: xyzls (LightEmission, raybnd_ in sim_telarray
# corsika/atmo.c) against the same integral, and its model offset -K/cz^3. K = 6.5e-4 ns from a
# scan at cz = 0.10, 0.056, 0.035 (K = 6.5, 6.2, 6.8 e-4 ns; the difference is about -K/cz^3,
# consistent with the raybnd_ term -(t s)^2/c^3). Grazing emission (cz < 0.1) is where it matters.
RAYBND_K_NS = 6.5e-4


def bent_excess_ns(
    z_m, cz, heights_m, index_minus_one, obs_m=OBS_LEVEL_M, points=20001
):
    """Exact travel-time excess [ns] of a downward ray from height obs_m + z_m to obs_m, relative
    to the vacuum path z_m / cz: (1/c) [integral of n(h) / cos(theta(h)) dh - z_m / cz] with
    n(h) sin(theta(h)) constant (Snell's law, plane-parallel layers). cz > 0 is |cos| of the
    emission zenith angle."""
    h = np.linspace(obs_m, obs_m + z_m, points)
    n = 1.0 + np.exp(np.interp(h, heights_m, np.log(index_minus_one)))
    n_emit = 1.0 + float(
        np.exp(np.interp(obs_m + z_m, heights_m, np.log(index_minus_one)))
    )
    sin_emit = math.sqrt(1.0 - cz * cz)
    cos_h = np.sqrt(1.0 - (sin_emit * n_emit / n) ** 2)
    total = float(np.trapezoid(n / cos_h, h) - z_m / cz)
    return total / source_camera.C_M_NS


def bent_excess_per_photon(z_m, cz_abs, heights_m, index_minus_one):
    """bent_excess_ns for an array of |cz|, interpolated on a grid over the photon range."""
    cz_abs = np.asarray(cz_abs, dtype=np.float64)
    grid = np.linspace(cz_abs.min(), cz_abs.max(), 201)
    values = np.array(
        [bent_excess_ns(z_m, c, heights_m, index_minus_one) for c in grid]
    )
    return np.interp(cz_abs, grid, values)


def cone_ratio_expected(distance_m, radius_m=TELESCOPE_RADIUS_CM / 100.0):
    """Expected native / reference photon ratio: obdeect emits into the cone that contains the
    1.02 R sphere, xyzls into the cone of the sphere R (solid angles at the source distance)."""
    rho = radius_m / distance_m

    def solid(k):
        return 1.0 - math.sqrt(1.0 - (k * rho) ** 2)

    return solid(1.02) / solid(1.0)


def timing_metrics(geo, emission, native, ref, metrics):
    """Timing keys of one distance: pass/fail obdeect vs exact (bent-ray integral) per photon;
    recorded xyzls vs exact on the reference's landing cut and on the straight landing cut,
    and the -K/cz^3 model offset. native: bent (propagated) columns; emission: emission columns."""
    heights, index_minus_one = source_camera.load_atmprof(source_camera.atmprof_path())
    cz_e, z_e = emission["cz"], emission["z"]
    if np.any(cz_e >= 0) or not np.allclose(z_e, z_e[0]):
        raise ValueError("timing: photons must go down from one source height")
    z_m = float(z_e[0])
    path = -z_e / cz_e  # m, straight to z = 0
    geo_t = path / source_camera.C_M_NS
    excess = bent_excess_per_photon(z_m, -cz_e, heights, index_minus_one)
    t_exact = geo_t + excess
    w = native["w"]
    sel_b = 100.0 * np.hypot(native["x"], native["y"]) < source_camera.TIME_RADIUS_CM
    landing_s = 100.0 * np.hypot(
        emission["x"] + path * emission["cx"], emission["y"] + path * emission["cy"]
    )
    sel_s = landing_s < source_camera.TIME_RADIUS_CM
    block = source_camera._time_block(
        native["t"][sel_b], w[sel_b], t_exact[sel_b], w[sel_b]
    )
    out = {
        "t_mean_native": block["t_mean_native"],
        "t_mean_exact": block["t_mean_reference"],
        "t_mean_delta": block["t_mean_delta"],
        "t_mean_delta_se": block["t_mean_delta_se"],
        "t_mean_limit_exceeded": block["t_mean_limit_exceeded"],
        "t_rms_native": block["t_rms_native"],
        "t_rms_exact": block["t_rms_reference"],
        "t_rms_delta": block["t_rms_delta"],
        "t_rms_limit_exceeded": block["t_rms_limit_exceeded"],
        "t_xyzls_mean_delta_recorded": metrics["t_mean_delta"],
    }
    rw = ref["w"][ref["r"] < source_camera.TIME_RADIUS_CM]
    rt = ref["t"][ref["r"] < source_camera.TIME_RADIUS_CM]
    xyzls_mean = float(np.sum(rt * rw) / np.sum(rw))
    pred = -RAYBND_K_NS / np.abs(cz_e) ** 3
    out.update(
        {
            "t_xyzls_minus_exact_recorded": xyzls_mean
            - float(np.average(t_exact[sel_b], weights=w[sel_b])),
            "t_xyzls_minus_exact_straight_cut_recorded": xyzls_mean
            - float(np.average(t_exact[sel_s], weights=w[sel_s])),
            "t_xyzls_predicted_offset_recorded": float(
                np.average(pred[sel_b], weights=w[sel_b])
            ),
            "t_xyzls_predicted_offset_straight_cut_recorded": float(
                np.average(pred[sel_s], weights=w[sel_s])
            ),
        }
    )
    distance = float(np.linalg.norm(geo["source_cm"]) / 100.0)
    expected = cone_ratio_expected(distance)
    ratio = metrics["photon_ratio"]
    se = metrics["photon_ratio_se"]
    out.update(
        {
            "photon_ratio_cone_expected": expected,
            "photon_ratio_over_cone_expected": ratio / expected,
            "photon_ratio_limit_exceeded": source_camera._flag(
                abs(ratio / expected - 1), max(0.005, 5 * se / expected)
            ),
        }
    )
    return out


# ---------------------------------------------------------------- one configuration


def propagate_bent(root, native_gst):
    """obdeect propagation of the native stream with the atmosphere (bent rays, air delay), no
    extinction; photons end on the plane z = 0 with their arrival time (xyzls -A 1 counterpart)."""
    propagated = root / "native-bent.gst"
    common.command(
        cli.propagation_run(
            native_gst,
            propagated,
            no_extinction=True,
            site_altitude=str(OBS_LEVEL_M),
            wavelength_range=("300", "800"),
            seed="1",
            profile=str(common.reference_config_dir() / "atmprof1.dat"),
        ),
        root / "propagate.log",
    )
    return propagated


def _camera_stage(root, propagated, reference_iact, geo, label):
    """Stage 2: propagated native stream through to-simtel, reference IACT as written by xyzls;
    both through sim_telarray with the telescope pointed at the illuminator."""
    candidate = root / "stage2-candidate.iact"
    common.command(
        cli.obdeect(
            "stream",
            "to-simtel",
            atmprof=str(common.reference_config_dir() / "atmprof1.dat"),
            o=candidate,
            extra=[propagated],
        ),
        root / "stage2-to-simtel.log",
    )
    options = {
        "altitude": OBS_LEVEL_M,
        "telescope_theta": f"{geo['telescope_theta_deg']:.6f}",
        "telescope_phi": f"{geo['telescope_phi_deg']:.6f}",
        **TRIGGER,
    }
    ref_out = common.sim_telarray(
        root,
        reference_iact,
        f"{label}-reference",
        options,
        config=common.REFERENCE / "cfg/CTA/CTA-PROD6-LST.cfg",
        defines=("-DCTA_NORTH", "-DLST1", "-DNUM_TELESCOPES=1"),
    )
    cand_out = common.sim_telarray(
        root,
        candidate,
        f"{label}-candidate",
        options,
        config=common.REFERENCE / "cfg/CTA/CTA-PROD6-LST.cfg",
        defines=("-DCTA_NORTH", "-DLST1", "-DNUM_TELESCOPES=1"),
    )
    a, b = source_camera.camera_totals(ref_out), source_camera.camera_totals(cand_out)
    images_ref, images_cand = (
        camera_event_images(ref_out),
        camera_event_images(cand_out),
    )
    counts = {
        "reference_triggered_telescopes": {str(k): len(v[0]) for k, v in a.items()},
        "candidate_triggered_telescopes": {str(k): len(v[0]) for k, v in b.items()},
    }
    if (
        not a
        or not b
        or any(
            len(v[0]) < source_camera.MIN_CAMERA_EVENTS
            for v in (*a.values(), *b.values())
        )
    ):
        return None, counts
    camera = source_camera.camera_metrics(a, b)
    tel = camera["telescope"]
    pixel = pixel_map_metrics(images_ref[tel], images_cand[tel])
    return {**camera, **pixel}, counts


def native_argv(geo, events, photons, output):
    """obdeect sources run point-source with the arguments of source_argv (SI units, m), -o output."""
    x, y, z = (v / 100.0 for v in geo["source_cm"])
    return cli.obdeect(
        "sources",
        "run",
        "point-source",
        observation_level=f"{OBS_LEVEL_M:.1f}",
        seed=str(SEED),
        atmosphere_id="1",
        events=str(events),
        photons=str(photons),
        position=f"{x:.8f},{y:.8f},{z:.8f}",
        direction=",".join(f"{v:.8f}" for v in geo["pointing"]),
        telescope_radius=f"{TELESCOPE_RADIUS_CM / 100.0:.4f}",
        bunch_size="1",
        spectrum=f"fixed:wavelength_nm={WAVELENGTH_NM}",
        angular=ANGULAR,
        o=output,
    )


def run_configuration(root, instrument, label, scale):
    """Emission and camera comparison of one illuminator distance; returns (metrics, details)."""
    root.mkdir(parents=True, exist_ok=True)
    geo = geometry(instrument, scale)
    native_gst, reference_iact = root / "native.gst", root / "reference.iact.gz"
    common.command(
        native_argv(geo, EVENTS, PHOTONS_PER_EVENT, native_gst),
        root / "native.log",
    )
    common.command(
        [
            "xyzls",
            f"-I{common.reference_config_dir()}",
            *source_argv(geo, EVENTS, PHOTONS_PER_EVENT),
            "-o",
            reference_iact,
        ],
        root / "reference.log",
    )
    bent_gst = propagate_bent(root, native_gst)
    native_counts, native = source_camera._read_native(bent_gst)
    ref_counts, ref = source_camera._read_reference(reference_iact)
    # Both cones contain the telescope sphere of radius TELESCOPE_RADIUS_CM (obdeect uses 1.02 R,
    # as in the point-source case): photons are compared inside R.
    metrics = source_camera.emission_metrics(
        native_counts,
        native,
        ref_counts,
        ref,
        "bent",
        acceptance_radius_cm=TELESCOPE_RADIUS_CM,
    )
    emission = source_camera._read_native(native_gst)[
        1
    ]  # emission z, cz (same photon order)
    metrics.update(timing_metrics(geo, emission, native, ref, metrics))
    details = {
        "distance_m": geo["distance_m"],
        "telescope_theta_deg": geo["telescope_theta_deg"],
        "telescope_phi_deg": geo["telescope_phi_deg"],
        "illuminator_m": geo["illuminator_m"].tolist(),
        "telescope_m": geo["telescope_m"].tolist(),
        "stage2": {"status": "BLOCKED", "reason": "not attempted"},
        "timing": "pass/fail: obdeect vs exact bent-ray integral (plane-parallel Snell, atmprof1); "
        "recorded: xyzls vs exact and the -K/cz^3 model offset (K = RAYBND_K_NS)",
    }
    source_camera._plot(
        native,
        ref,
        root / f"stage1-{label}.png",
        source_camera.compared_times(native, ref, {"time_correction": "bent"}),
    )
    try:
        camera, counts = _camera_stage(root, bent_gst, reference_iact, geo, label)
        if camera is None:
            details["stage2"] = {
                "status": "BLOCKED",
                "reason": f"fewer than {source_camera.MIN_CAMERA_EVENTS} triggered events per telescope",
                **counts,
            }
        else:
            details["stage2"] = {"status": "OK", **counts}
            metrics.update(
                {
                    f"camera_{k}": v
                    for k, v in camera.items()
                    if isinstance(v, (int, float))
                }
            )
    except Exception as exc:  # noqa: BLE001 stage 2 is reported, not hidden
        details["stage2"] = {
            "status": "BLOCKED",
            "reason": str(exc)[:400],
            "traceback": traceback.format_exc()[-800:],
        }
    return metrics, details


def analyse(root: Path, args) -> tuple[dict, dict]:
    """Both distances of the illuminator args.instrument in root (see common.results_main)."""
    instrument = args.instrument
    if instrument not in ILLUMINATORS:
        raise ValueError(f"unknown illuminator {instrument}")
    metrics, details = (
        {},
        {"case": f"{instrument} versus xyzls (LightEmission)", "configurations": {}},
    )
    for label, scale in DISTANCE_SCALES.items():
        sub = root / label
        m, d = run_configuration(sub, instrument, label, scale)
        for key, value in m.items():
            metrics[f"{label}_{key}"] = value
        details["configurations"][label] = d
    details["mapping"] = {
        "geometry": "simtools: illuminator = ground position + illuminator_tower_height; telescope = "
        f"{ILLUMINATORS[instrument]['telescope']} ground position; pointing from illuminator to telescope",
        "source": f"xyzls / point-source: {WAVELENGTH_NM} nm, pulse {PULSE}, angular {ANGULAR}, "
        f"{PHOTONS_PER_EVENT} photons per event, {EVENTS} events, observation level {OBS_LEVEL_M} m",
        "design_values": "ILLN-design / ILLS-design: flasher wavelength 300 nm, pulse Gauss 0, "
        "angular distribution gauss width 0 (isotropic)",
        "camera": "sim_telarray LST-1 North configuration, Bypass_Optics=0, telescope pointed at the "
        "illuminator; the native stream, propagated with atmprof1 (bent rays, no extinction), goes "
        "through to-simtel",
        "time": "pass/fail: obdeect arrival time at z = 0 (propagation with atmprof1) vs the exact "
        "bent-ray integral per photon; xyzls times recorded only (grazing raybnd_ term)",
        "limits": "see the module docstring",
    }
    return metrics, details


def _add_arguments(parser) -> None:
    parser.add_argument("instrument", help=f"illuminator: {', '.join(ILLUMINATORS)}")


def main() -> None:
    common.results_main(analyse, __doc__.splitlines()[0], _add_arguments)


if __name__ == "__main__":
    main()
