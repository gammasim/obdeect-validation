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

import json
import os
import shutil
import subprocess
import sys
import traceback
from pathlib import Path

import numpy as np
from eventio import SimTelFile
from eventio.iact import IACTFile
from obdeect.stream import Reader

from obdeect_validation.adapters import chain_showers

C_M_NS = 0.299792458
TIME_RADIUS_CM = 1500.0
KEYS = ("cx", "cy", "cz", "t", "z", "w", "wavelength", "x", "y")

# Matched parameter sets. Flasher cases follow obdeect/validation/sources/compare_*.sh; point
# source follows the xyzls defaults (-r 1500, -s 400, -h 2150 m, -d 0,0,-1); beam scan follows the
# xyls defaults with -r 123 (obdeect beam-scan defaults are the same fractions of 123 cm).
CASES = {
    "flasher-dish": {
        "native": ["obdeect", "sources", "flasher", "--type", "dish", "--events", "20"],
        "reference": ["ff-1m", "--events", "20"],
        "stage2": {
            "path": "propagate", "bypass": 1,
            "native": ["obdeect", "sources", "flasher", "--type", "dish", "--events", "5", "--photons", "150000000"],
            "reference": ["ff-1m", "--events", "5", "--photons", "150000000"],
            "photon_factor": 60,
        },
        "frame": {"native": "sky (emission at flasher, z = 16 m)", "reference": "observation plane z = 0"},
        "time_correction": "geometric",
        "native_frame": "sky",
    },
    "flasher-secondary": {
        "native": ["obdeect", "sources", "flasher", "--type", "secondary", "--photons", "1e5", "--events", "1"],
        "reference": ["ff-gct", "--photons", "1e5", "--events", "1"],
        "frame": {"native": "telescope (T)", "reference": "ff-gct photon bunches (no cz)"},
        "time_correction": "geometric",
        "native_frame": "telescope",
        "stage2": {
            "path": "direct", "bypass": 2, "cfg": "CTA/CTA-PROD6-SST.cfg", "camera_map": True,
            "defines": ["-DCTA_NORTH", "-DSST1"],
            "native": ["obdeect", "sources", "flasher", "--type", "secondary", "--photons", "10000000", "--events", "1"],
            "reference": ["ff-gct", "--photons", "10000000", "--events", "1"],
            "photon_factor": 100,
        },
    },
    "point-source": {
        "native": ["obdeect", "sources", "point-source", "-N", "30", "-n", "10000000", "-z", "100000",
                   "-a", "isotropic", "--seed", "17"],
        "reference": ["xyzls", "-N", "30", "-n", "10000000", "-z", "100000", "-a", "isotropic"],
        "frame": {"native": "sky (emission at source, z = 1000 m)", "reference": "plane through telescope sphere centre"},
        "time_correction": "geometric",
        "native_frame": "sky",
        # Native emits into a cone of 1.02 R (deliberate margin), xyzls into R: compare the photons
        # crossing the plane z = 0 within the same radius, inside both cones.
        "acceptance_radius_cm": 1400.0,
        "stage2": {
            # Photons reach the telescope sphere, not the camera: full optics (simtools: bypass only for flat fielding).
            "path": "propagate", "bypass": 0,
            "native": ["obdeect", "sources", "point-source", "-N", "2", "-n", "3000000000", "-z", "100000",
                       "-a", "isotropic", "--seed", "17"],
            "reference": ["xyzls", "-N", "2", "-n", "3000000000", "-z", "100000", "-a", "isotropic"],
            "photon_factor": 300,
        },
    },
    "beam-scan": {
        "native": ["obdeect", "sources", "beam-scan", "-r", "123", "--seed", "1"],
        "reference": ["xyls", "-r", "123"],
        "stage2": {"blocked": "no intensity option: obdeect beam-scan and xyls emit a fixed 100 photons per "
                              "position, so the flat-field intensity cannot be matched between the two"},
        "frame": {"native": "sky (beam at z = 6.15 cm)", "reference": "plane z = 0"},
        "time_correction": "geometric",
        "native_frame": "sky",
    },
}


def _read_native(path):
    counts, cols = [], {k: [] for k in KEYS}
    for record in Reader(str(path)):
        if record["type"] == "event":
            counts.append(0.0)
        elif record["type"] == "photons":
            columns = record["columns"]
            if not counts:
                raise ValueError("photons before the first event record")
            if "cz" not in columns:
                columns = {**columns, "cz": -np.sqrt(np.clip(1 - columns["cx"] ** 2 - columns["cy"] ** 2, 0, None))}
            w = np.asarray(columns["w"], dtype=np.float64)
            counts[-1] += float(np.sum(w))
            for key in KEYS:
                if key in columns:
                    cols[key].append(np.asarray(columns[key], dtype=np.float64))
                elif key == "z":
                    cols[key].append(np.zeros(len(w)))
    return counts, {k: np.concatenate(v) if v else np.empty(0) for k, v in cols.items()}


def _read_reference(path):
    counts, cols = [], {k: [] for k in ("cx", "cy", "t", "w", "wavelength", "r")}
    with IACTFile(str(path)) as stream:
        for event in stream:
            total = 0.0
            for bunches in event.photon_bunches.values():
                photons = np.asarray(bunches["photons"], dtype=np.float64)
                total += float(np.sum(photons))
                cols["cx"].append(np.asarray(bunches["cx"], dtype=np.float64))
                cols["cy"].append(np.asarray(bunches["cy"], dtype=np.float64))
                # Bunch position [cm] from the telescope axis, for the camera acceptance cut.
                cols["r"].append(np.hypot(np.asarray(bunches["x"], dtype=np.float64),
                                          np.asarray(bunches["y"], dtype=np.float64)))
                # Reference times are relative to the event time offset; restore the origin.
                cols["t"].append(np.asarray(bunches["time"], dtype=np.float64) + float(event.time_offset))
                cols["w"].append(photons)
                cols["wavelength"].append(np.asarray(bunches["wavelength"], dtype=np.float64))
            counts.append(total)
    out = {k: np.concatenate(v) if v else np.empty(0) for k, v in cols.items()}
    out["cz"] = -np.sqrt(np.clip(1.0 - out["cx"] ** 2 - out["cy"] ** 2, 0.0, None))
    return counts, out


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
    radius = 100.0 * np.hypot(native["x"] + path * native["cx"], native["y"] + path * native["cy"])
    wn = native["w"][np.isfinite(radius) & (radius < radius_cm)]
    wr = ref["w"][ref["r"] < radius_cm]
    if not wn.size or not wr.size:
        raise ValueError("no photons within the acceptance radius")
    ratio = float((wn.sum() / n_events) / (wr.sum() / r_events))
    # Effective photon numbers (bunches carry weights).
    se = ratio * np.sqrt(np.sum(wn**2) / wn.sum() ** 2 + np.sum(wr**2) / wr.sum() ** 2)
    return ratio, float(se)


def emission_metrics(native_counts, native, ref_counts, ref, time_correction, acceptance_radius_cm=None):
    """Stage 1 numbers and integer limit flags. time_correction is 'geometric' or 'none'."""
    n_ev, r_ev = len(native_counts), len(ref_counts)
    metrics = {"native_event_count": n_ev, "reference_event_count": r_ev,
               "event_count_delta": abs(n_ev - r_ev)}
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
    metrics.update({"photons_per_event_native": float(nc.mean()),
                    "photons_per_event_reference": float(rc.mean()),
                    "photon_ratio": ratio, "photon_ratio_se": se_ratio,
                    "photon_ratio_limit_exceeded": _flag(abs(ratio - 1), max(0.005, 5 * se_ratio))})
    if acceptance_radius_cm:
        accepted, se_accepted = accepted_ratio(native, n_ev, ref, r_ev, acceptance_radius_cm)
        metrics.update({"photon_ratio_accepted": accepted, "photon_ratio_accepted_se": se_accepted,
                        "photon_ratio_accepted_limit_exceeded":
                            _flag(abs(accepted - 1), max(0.005, 5 * se_accepted))})
    wn, wr = native["w"], ref["w"]
    t_native = native["t"].copy()
    native_ok = np.ones(t_native.size, dtype=bool)
    if time_correction == "geometric":
        cz = native["cz"]
        if np.any(cz == 0):
            raise ValueError("native photons parallel to the plane z = 0")
        # Only photons heading for the plane (cz < 0) have a flight path to it (about half of the
        # flasher-secondary photons do not); their time moments are not compared.
        toward = cz < 0
        path = np.where(toward, -native["z"] / np.where(toward, cz, -1.0), 0.0)  # m
        t_native = t_native + path / C_M_NS
        # Same camera acceptance as the reference: arrival position within TIME_RADIUS_CM.
        xp, yp = native["x"] + path * native["cx"], native["y"] + path * native["cy"]
        native_ok = toward & (100.0 * np.hypot(xp, yp) < TIME_RADIUS_CM)
        metrics["native_photons_moving_away_fraction"] = float(np.mean(~toward))
        metrics["t_native_outside_radius_fraction"] = float(1.0 - np.mean(native_ok))
    metrics["time_correction_ns_mean"] = float(np.mean(t_native - native["t"]))
    for key in ("cx", "cy"):
        mn, sn = _weighted_mean(native[key], wn)
        mr, sr = _weighted_mean(ref[key], wr)
        se = float(np.hypot(sn, sr))
        d_ks, p_ks = _ks(native[key], ref[key])
        metrics.update({f"{key}_native": mn, f"{key}_reference": mr, f"{key}_delta": mn - mr,
                        f"{key}_delta_sigma": float((mn - mr) / se) if se else 0.0,
                        f"{key}_mean_limit_exceeded": _flag(abs(mn - mr) / se if se else 0.0, 5),
                        f"{key}_ks_distance": d_ks, f"{key}_ks_pvalue": p_ks,
                        f"{key}_ks_limit_exceeded": int(d_ks >= 0.01 and p_ks <= 1e-3)})
    # Time moments use reference bunches within the fiducial radius of 1500 cm (xyzls -r default).
    # Bunches beyond it (flasher-secondary: 3.6 %, incl. times near -2e8 ns) are not camera photons;
    # their fraction is reported, not compared.
    good = ref["r"] < TIME_RADIUS_CM
    metrics["t_reference_outside_radius_fraction"] = float(1.0 - np.mean(good)) if good.size else 0.0
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
        metrics.update({f"{key}_mean_native": mn, f"{key}_mean_reference": mr,
                        f"{key}_mean_delta": mn - mr, f"{key}_mean_delta_se": se,
                        f"{key}_mean_limit_exceeded": _flag(abs(mn - mr), max(tol, 5 * se)),
                        f"{key}_rms_native": rn, f"{key}_rms_reference": rr,
                        f"{key}_rms_delta": rn - rr,
                        f"{key}_rms_limit_exceeded": _flag(abs(rn - rr), max(tol, 5 * se_rms))})
    return metrics


def camera_events(path):
    """Per-event photoelectron totals and images: {telescope: (n_pe list, image sum, n_pixels)}."""
    totals, images = {}, {}
    with SimTelFile(str(path)) as stream:
        for event in stream:
            for tel, pe in event.get("photoelectrons", {}).items():
                image = np.asarray(pe["photoelectrons"], dtype=np.float64)
                totals.setdefault(tel, []).append(float(pe["n_pe"]))
                images[tel] = images.get(tel, 0) + image
    return {tel: (np.asarray(totals[tel]), images[tel]) for tel in totals}


def camera_metrics(reference, candidate):
    """Independent samples: mean p.e. ratio with standard error and summed image chi2 per pixel."""
    common = sorted(set(reference) & set(candidate))
    if not common:
        raise ValueError("no telescope with photoelectrons in both samples")
    tel = common[0]
    na, ima = reference[tel]
    nb, imb = candidate[tel]
    mean_a, mean_b = na.mean(), nb.mean()
    ratio = float(mean_b / mean_a)
    se = float(ratio * np.hypot(na.std(ddof=1) / np.sqrt(na.size) / mean_a,
                                nb.std(ddof=1) / np.sqrt(nb.size) / mean_b))
    # Pixel shapes normalised to unit sum; variance of a normalised Poisson count.
    pa, pb = ima / ima.sum(), imb / imb.sum()
    var = ima / ima.sum() ** 2 + imb / imb.sum() ** 2
    used = var > 0
    chi2 = float(np.sum((pa[used] - pb[used]) ** 2 / var[used]))
    ndof = max(int(np.count_nonzero(used)) - 1, 1)
    reduced = chi2 / ndof
    return {"telescope": int(tel), "reference_events": int(na.size), "candidate_events": int(nb.size),
            "pe_per_event_reference": float(mean_a), "pe_per_event_candidate": float(mean_b),
            "pe_ratio": ratio, "pe_ratio_se": se,
            "pe_ratio_limit_exceeded": _flag(abs(ratio - 1), max(0.01, 5 * se)),
            "image_reduced_chi2": reduced, "image_chi2_limit_exceeded": _flag(reduced, 2.0)}


def _plot(native, ref, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    for ax, key in zip(axes.flat, ("cx", "cy", "t", "wavelength")):
        a, b = native[key], ref[key]
        if a.size and b.size:
            lo, hi = np.quantile(np.concatenate((a, b)), [0.001, 0.999])
            hi = hi if hi > lo else lo + 1e-6
            ax.hist(a, bins=60, range=(lo, hi), histtype="step", density=True, label="obdeect")
            ax.hist(b, bins=60, range=(lo, hi), histtype="step", density=True, label="LightEmission")
        ax.set(xlabel=key)
    axes[0, 0].legend()
    fig.tight_layout()
    fig.savefig(output, dpi=120)
    plt.close(fig)


def _run(argv, log):
    with Path(log).open("wb") as handle:
        result = subprocess.run([str(a) for a in argv], stdout=handle, stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        tail = Path(log).read_text(errors="replace").strip().splitlines()[-3:]
        raise RuntimeError(f"{argv[0]} exited {result.returncode}: {' | '.join(tail)}")


def _simtel(root, input_file, name, bypass, cfg="CTA/CTA-PROD6-LST.cfg", defines=("-DCTA_NORTH", "-DLST1")):
    """sim_telarray with the options of simtools' light-emission step (released binary)."""
    binary = os.environ.get("SIMTEL_BINARY", "/workdir/simulation_software/sim_telarray/bin/sim_telarray")
    reference = Path(os.environ.get("SIMTEL_REFERENCE_SOURCE", "/reference"))
    output = root / f"{name}.simtel"
    output.unlink(missing_ok=True)
    options = {"altitude": 2150, "TRIGGER_TELESCOPES": 1, "TELTRIG_MIN_SIGSUM": 2, "PULSE_ANALYSIS": -30,
               "MAXIMUM_TELESCOPES": 1, "telescope_theta": 0, "telescope_phi": 0, "Bypass_Optics": bypass,
               "output_file": output, "histogram_file": root / f"{name}.hdata",
               "random_seed": 20261113, "save_photons": 2, "store_photoelectrons": 0}
    argv = [binary, f"-I{reference / 'cfg/common'}", f"-I{reference / 'cfg/CTA'}",
            "-c", reference / "cfg" / cfg, *defines, "-DNUM_TELESCOPES=1", "-W", "random_state=none"]
    for key, value in options.items():
        argv += ["-C", f"{key}={value}"]
    argv.append(input_file)
    chain_showers.command(argv, root / f"{name}.log", cwd=reference)
    if not output.is_file():
        raise RuntimeError(f"no sim_telarray output: {name}.log")
    return output


def _camera_stage(root, native_gst, reference_file, stage2):
    """Stage 2: candidate and reference through the same sim_telarray options."""
    candidate = root / "candidate.iact"
    if stage2["path"] == "propagate":
        propagated = root / "propagated.gst"
        _run(["obdeect", "propagation", "run", "--no-extinction", "--site-altitude", "2150", "--seed", "1",
              "--wavelength-range", "300", "800", "-i", native_gst, "-o", propagated], root / "propagation.log")
        _run(["obdeect", "stream", "to-simtel", "-o", candidate, propagated], root / "to-simtel.log")
    else:
        argv = ["obdeect", "stream", "to-simtel", "-o", candidate]
        if stage2.get("camera_map"):
            # Telescope 1 of the reference program as a dual-mirror (class 2) telescope in the telescope frame.
            (root / "camera-map.txt").write_text("1 2 0.0\n", encoding="utf-8")
            argv += ["--camera-map", root / "camera-map.txt"]
        _run([*argv, native_gst], root / "to-simtel.log")
    cfg = stage2.get("cfg", "CTA/CTA-PROD6-LST.cfg")
    defines = stage2.get("defines", ["-DCTA_NORTH", "-DLST1"])
    ref_out = _simtel(root, reference_file, "reference", stage2["bypass"], cfg, defines)
    cand_out = _simtel(root, candidate, "candidate", stage2["bypass"], cfg, defines)
    a, b = camera_events(ref_out), camera_events(cand_out)
    if not a or not b or any(len(v[0]) == 0 for v in list(a.values()) + list(b.values())):
        return None, {"reference_triggered_telescopes": {str(k): int(len(v[0])) for k, v in a.items()},
                      "candidate_triggered_telescopes": {str(k): int(len(v[0])) for k, v in b.items()}}
    return camera_metrics(a, b), {}


def _reference_argv(argv, root):
    """LightEmission command with the configuration directory the program needs."""
    config = os.environ.get("SIMTEL_CONFIG_PATH", "/reference/cfg/common")
    argv = list(argv)
    if argv[0] in ("ff-1m", "ff-gct", "xyzls"):
        argv.insert(1, f"-I{config}")
    elif argv[0] == "xyls":
        # xyls has no configuration option and reads the atmosphere table from the working directory.
        shutil.copyfile(Path(config) / "atm_profile_model_1.dat", root / "atm_profile_model_1.dat")
    return argv


def main(case):
    spec = CASES[case]
    root = Path.cwd()
    native_gst, reference_file = root / "native.gst", root / "reference.dat.gz"
    _run([*spec["native"], "-o", native_gst], root / "native.log")
    _run([*_reference_argv(spec["reference"], root), "-o", reference_file], root / "reference.log")
    native_counts, native = _read_native(native_gst)
    ref_counts, ref = _read_reference(reference_file)
    metrics = emission_metrics(native_counts, native, ref_counts, ref, spec["time_correction"],
                               spec.get("acceptance_radius_cm"))
    _plot(native, ref, root / "stage1.png")
    details = {"case": case, "native_command": spec["native"], "reference_command": spec["reference"],
               "frames": spec["frame"], "time_correction": spec["time_correction"],
               "time_correction_rule": "t_arrival = t + (-z/cz)/c (z in m, cz<0); ref times + time_offset",
               "stage2": {"status": "BLOCKED", "reason": "not attempted"}}
    stage2 = spec["stage2"]
    try:
        if "blocked" in stage2:
            raise RuntimeError(stage2["blocked"])
        native2, reference2 = root / "stage2-native.gst", root / "stage2-reference.dat.gz"
        _run([*stage2["native"], "-o", native2], root / "stage2-native.log")
        _run([*_reference_argv(stage2["reference"], root), "-o", reference2], root / "stage2-reference.log")
        camera, reason = _camera_stage(root, native2, reference2, stage2)
        if camera is None:
            details["stage2"] = {"status": "BLOCKED",
                                 "reason": "no triggered telescope events in sim_telarray for one input",
                                 **reason}
        else:
            details["stage2"] = {"status": "OK"}
            metrics.update({f"camera_{k}": v for k, v in camera.items() if isinstance(v, (int, float))})
    except Exception as exc:  # stage 2 is reported, not hidden
        details["stage2"] = {"status": "BLOCKED", "reason": str(exc)[:400],
                             "traceback": traceback.format_exc()[-800:]}
    if "native" in stage2:
        details["stage2"].update({"native_command": stage2["native"], "reference_command": stage2["reference"],
                                  "photon_factor": stage2["photon_factor"],
                                  "photon_factor_note": "same factor applied to native and reference; "
                                                        "stage 1 uses the matched stage 1 counts"})
    (root / "metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (root / "details.json").write_text(json.dumps(details, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1])
