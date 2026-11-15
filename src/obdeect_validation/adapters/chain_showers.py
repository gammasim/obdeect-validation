"""Real CORSIKA showers through the stream codec and sim_telarray camera."""

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
from eventio import SimTelFile
from eventio.iact import IACTFile

from obdeect_validation.adapters import common
from obdeect_validation.adapters.common import command  # noqa: F401  (re-exported for the adapters)

PE_BIAS_FLOOR = 0.01  # mean p.e. per telescope-event within 1 % (plan step 5)
TIME_BIAS_FLOOR_NS = 0.05
SHAPE_CHI2_LIMIT = 2.0
SIGMA = 5.0  # limit = max(floor, SIGMA * standard error)


FLAT_TELESCOPES = (
    (0.0, 0.0, 0.0, 1200.0),
    (5000.0, 0.0, 0.0, 1200.0),
    (0.0, 5000.0, 0.0, 1200.0),
)  # (x, y, z, radius) in cm, as the flat showers of this chain have always used
FLAT_SCATTER_CM = 10000.0  # CSCAT radius in cm
SEED_OFFSETS = (0, 12, 26, 44)  # SEED lines are seed_base + offset: 17, 29, 43, 61


def card_text(
    *,
    showers,
    zenith,
    energy_min,
    energy_max,
    data_dir,
    direct,
    telfil,
    telescopes=FLAT_TELESCOPES,
    scatter_cm=FLAT_SCATTER_CM,
    seed_base=17,
    tstart=False,
):
    """CORSIKA input card of the IACT showers; the values are written as given.

    data_dir is the DATDIR, direct the DIRECT prefix (with trailing slash) and telfil the IACT
    output name below it. Telescope positions, radii and the CSCAT radius are in cm (CORSIKA
    units). With the defaults the card is the one generate() has always written.
    """
    seeds = "".join(f"SEED {seed_base + offset} 0 0\n" for offset in SEED_OFFSETS)
    scatter = "".join(f"TELESCOPE {x} {y} {z} {r}\n" for x, y, z, r in telescopes)
    start = "TSTART T\n" if tstart else ""
    return f"""RUNNR 1
EVTNR 1
NSHOW {showers}
PRMPAR 1
ERANGE {energy_min} {energy_max}
ESLOPE -2.0
THETAP {zenith} {zenith}
PHIP 0.0 0.0
{seeds}OBSLEV 215600.0
ATMOSPHERE 1 Y
{start}MAGNET 20.0 40.0
ARRANG 0.0
HADFLG 0 0 0 0 0 2
ECUTS 0.3 0.3 0.02 0.02
ELMFLG F F
CSCAT 5 {scatter_cm} 0.0
DATDIR {data_dir}/
DIRECT {direct}
{scatter}TELFIL {direct}{telfil}
IACT SPLIT_AUTO 0
IACT IO_BUFFER 800MB
IACT MAX_BUNCHES 1000000
EXIT
"""


def prepare_data(root):
    """Link the CORSIKA data files into root/data and copy the atmosphere; return SHA-256 sums."""
    root.mkdir(parents=True, exist_ok=True)
    data = root / "data"
    data.mkdir(exist_ok=True)
    source = Path(os.environ.get("CORSIKA_DATA_DIR", "/corsika-data"))
    names = [
        "NUCNUCCS",
        "UrQMD-1.3.1-xs.dat",
        "QGSDAT01",
        "qgsdat-II-04",
        "sectnu-II-04",
    ]
    names += [p.name for p in source.glob("EGSDAT6_*")]
    hashes = {}
    for name in names:
        original = source / name
        if not original.is_file():
            raise FileNotFoundError(original)
        destination = data / {
            "qgsdat-II-04": "QGSDAT-II-04",
            "sectnu-II-04": "SECTNU-II-04",
        }.get(name, name)
        if not destination.exists():
            destination.symlink_to(original)
        hashes[name] = common.sha256(original)
    shutil.copyfile(common.REFERENCE / "cfg/common/atmprof1.dat", root / "atmprof1.dat")
    return hashes


def generate(
    root, showers=20, zenith=20, energy_min=100, energy_max=3000, geometry="flat"
):
    """Run CORSIKA in root for IACT showers and return root/showers.iact.

    geometry "flat" keeps the card, executable and environment default this chain has always
    used. geometry "curved" uses the compact-format layout of curved_showers and the curved
    executable (CORSIKA_BINARY still overrides the executable).
    """
    if geometry == "flat":
        layout = {"telescopes": FLAT_TELESCOPES, "scatter_cm": FLAT_SCATTER_CM}
    elif geometry == "curved":
        from obdeect_validation.adapters import curved_showers

        layout = curved_showers.layout(zenith)
    else:
        raise ValueError(f"unknown CORSIKA geometry: {geometry!r}")
    hashes = prepare_data(root)
    card = card_text(
        showers=showers,
        zenith=zenith,
        energy_min=energy_min,
        energy_max=energy_max,
        data_dir=root / "data",
        direct=f"{root}/",
        telfil="showers.iact",
        **layout,
    )
    (root / "corsika.card").write_text(card)
    corsika = os.environ.get("CORSIKA_BINARY") or (
        "/workdir/simulation_software/corsika7/corsika_qgs2_urqmd_flat"
        if geometry == "flat"
        else str(curved_showers.corsika_executable("curved"))
    )
    with (
        (root / "corsika.card").open("rb") as stdin,
        (root / "corsika.log").open("wb") as log,
    ):
        result = subprocess.run(
            [corsika], cwd=root, stdin=stdin, stdout=log, stderr=subprocess.STDOUT
        )
    if result.returncode or not (root / "showers.iact").is_file():
        raise RuntimeError(f"CORSIKA failed: see {root}/corsika.log")
    (root / "input-hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")
    return root / "showers.iact"


def bunches(path):
    result = {}
    with IACTFile(path) as stream:
        for number, event in enumerate(stream):
            for telescope, photons in event.photon_bunches.items():
                result[number, telescope] = np.array(photons, copy=True)
    return result


def compare_photons(reference, candidate):
    left, right = bunches(reference), bunches(candidate)
    metrics = {
        "photon_event_telescope_key_delta": len(set(left) ^ set(right)),
        "bunch_count_delta": 0,
        "max_codec_normalized_error": 0.0,
    }
    per_field = {}
    for key in left.keys() & right.keys():
        a, b = left[key], right[key]
        metrics["bunch_count_delta"] += abs(len(a) - len(b))
        if len(a) != len(b):
            continue
        for field in a.dtype.names:
            av, bv = a[field].astype(float), b[field].astype(float)
            if not np.isfinite(av).all() or not np.isfinite(bv).all():
                raise ValueError(f"nonfinite {field}")
            # EventIO cm positions undergo a float32 cm -> m -> cm round trip.
            error = np.max(
                np.abs(av - bv) / (1e-6 * np.maximum(np.abs(av), 1.0)), initial=0
            )
            per_field[field] = max(per_field.get(field, 0), float(error))
    metrics["max_codec_normalized_error"] = max(per_field.values(), default=0)
    metrics["codec_bounds_exceeded"] = sum(value > 1 for value in per_field.values())
    return metrics, per_field


def simulate_argv(
    root, input_file, name, bypass=0, config=None, extra=None, binary=None
):
    """argv and output of the three-telescope sim_telarray run (see simulate)."""
    return common.sim_telarray_argv(
        root,
        input_file,
        name,
        {
            "maximum_telescopes": 3,
            "altitude": 2156,
            "telescope_theta": 20,
            "telescope_phi": 180,
            "Bypass_Optics": bypass,
        },
        config=config,
        defines=("-DCTA_NORTH", "-DLST1", "-DNUM_TELESCOPES=3"),
        extra=extra or (),
        binary=binary,
    )


def simulate(
    root,
    input_file,
    name,
    bypass=0,
    config=None,
    extra=None,
    runner=command,
    binary=None,
):
    """Three-telescope sim_telarray run of the CORSIKA showers (zenith 20 deg, azimuth 180 deg)."""
    return common.sim_telarray(
        root,
        input_file,
        name,
        {
            "maximum_telescopes": 3,
            "altitude": 2156,
            "telescope_theta": 20,
            "telescope_phi": 180,
            "Bypass_Optics": bypass,
        },
        config=config,
        defines=("-DCTA_NORTH", "-DLST1", "-DNUM_TELESCOPES=3"),
        extra=extra or (),
        runner=runner,
        binary=binary,
    )


def camera_events(path):
    events = {}
    with SimTelFile(path) as stream:
        for event in stream:
            events[str(event["event_id"])] = event
    return events


def bias_exceeded(value, standard_error, floor):
    """Integer flag: |value| above max(floor, SIGMA * standard_error)."""
    return int(abs(value) > max(floor, SIGMA * standard_error))


def _standard_error(values):
    values = np.asarray(values, dtype=float)
    if values.size < 2:
        raise ValueError("camera comparison needs at least two matched samples")
    return float(np.std(values, ddof=1) / np.sqrt(values.size))


def compare_camera(left, right, output):
    """Paired per-telescope comparison of two sim_telarray runs; writes the CSV and the figure."""
    if not left or not right:
        raise ValueError("camera comparison requires nonempty triggered event samples")
    keys = sorted(set(left) & set(right), key=int)
    pairs, images, times, adc, shape_chi2 = [], [], [], [], []
    trigger_changes = 0
    for key in keys:
        a, b = left[key], right[key]
        trigger_changes += (
            a["trigger_information"]["trigger_pattern"]
            != b["trigger_information"]["trigger_pattern"]
        )
        for tel in set(a.get("photoelectrons", {})) & set(b.get("photoelectrons", {})):
            ap, bp = a["photoelectrons"][tel], b["photoelectrons"][tel]
            pairs.append((int(key), int(tel), float(ap["n_pe"]), float(bp["n_pe"])))
            av, bv = (
                np.asarray(ap["photoelectrons"], dtype=float),
                np.asarray(bp["photoelectrons"], dtype=float),
            )
            images.append(float(np.abs(av - bv).sum() / max(av.sum(), 1)))
            mask = (av + bv) > 0
            if mask.any():
                shape_chi2.append(
                    float(np.mean((av[mask] - bv[mask]) ** 2 / (av[mask] + bv[mask])))
                )
            if len(ap["time"]) and len(bp["time"]):
                times.append(float(np.mean(bp["time"]) - np.mean(ap["time"])))
        for tel in set(a["telescope_events"]) & set(b["telescope_events"]):
            av = a["telescope_events"][tel]["adc_samples"].astype(float)
            bv = b["telescope_events"][tel]["adc_samples"].astype(float)
            adc.append(float(np.mean(bv - av)))
    if not pairs:
        raise ValueError("no matched photoelectron samples")
    if not shape_chi2 or not times:
        raise ValueError("no matched images or pixel times to compare")
    array = np.asarray(pairs)
    difference = array[:, 3] - array[:, 2]
    bias = float(difference.sum() / array[:, 2].sum())
    se = _standard_error(difference) / float(np.mean(array[:, 2]))
    np.savetxt(
        output / "paired-photoelectrons.csv",
        array,
        delimiter=",",
        header="event,telescope,reference,candidate",
        comments="",
    )
    _plot_camera(array, difference, output / "camera-comparison.png")
    mean_time = float(np.mean(times))
    return {
        "matched_event_telescope_samples": len(pairs),
        "photoelectron_relative_bias": bias,
        "photoelectron_relative_bias_se": se,
        "photoelectron_bias_limit_exceeded": bias_exceeded(bias, se, PE_BIAS_FLOOR),
        "mean_image_reduced_chi2": float(np.mean(shape_chi2)),
        "image_shape_limit_exceeded": int(np.mean(shape_chi2) > SHAPE_CHI2_LIMIT),
        "mean_arrival_time_difference_ns": mean_time,
        "arrival_time_bias_limit_exceeded": bias_exceeded(
            mean_time, _standard_error(times), TIME_BIAS_FLOOR_NS
        ),
        "max_image_relative_l1_difference": max(images),
        "mean_image_relative_l1_difference": float(np.mean(images)),
        "max_mean_arrival_time_difference_ns": max(map(abs, times)),
        "mean_adc_sample_difference": float(np.mean(adc)),
        "trigger_pattern_mismatches": int(trigger_changes),
    }


def _plot_camera(array, difference, path):
    """Paired p.e. (colour per telescope) against the identity line, and the relative differences."""
    from obdeect_validation.adapters import plotting

    fig, axes = plotting.figure(1, 2, (10, 4))
    for i, tel in enumerate(sorted({int(t) for t in array[:, 1]})):
        rows = array[:, 1] == tel
        axes[0, 0].scatter(
            array[rows, 2],
            array[rows, 3],
            s=14,
            color=plotting.CYCLE[i % len(plotting.CYCLE)],
            alpha=0.8,
            label=f"telescope {tel}",
        )
    top = float(max(array[:, 2].max(), array[:, 3].max()))
    axes[0, 0].plot(
        [0, top], [0, top], color=plotting.ANALYTIC, linewidth=1, label="identity"
    )
    axes[0, 0].set(xlabel="sim_telarray reference p.e.", ylabel="candidate p.e.")
    axes[0, 0].legend()
    relative = difference / np.maximum(array[:, 2], 1.0)
    axes[0, 1].hist(relative, bins=30, color=plotting.NATIVE, alpha=0.85)
    axes[0, 1].axvline(0, color=plotting.ANALYTIC, linewidth=1)
    axes[0, 1].set(
        xlabel="(candidate - reference) / reference p.e.",
        ylabel="event/telescope samples",
    )
    plotting.save(fig, path)


def roundtrip(root, source):
    """CORSIKA IACT -> stream -> sim_telarray input, compared bunch by bunch and in the camera.

    Returns (native stream, metrics, per-field codec errors, reference events, candidate events).
    """
    native, converted = root / "showers.gst", root / "roundtrip.iact"
    command(
        ["obdeect", "stream", "from-corsika", source],
        root / "from-corsika.log",
        stdout_file=native,
    )
    command(
        ["obdeect", "stream", "to-simtel", native, "-o", converted],
        root / "to-simtel.log",
    )
    metrics, fields = compare_photons(source, converted)
    a = camera_events(simulate(root, source, "direct"))
    b = camera_events(simulate(root, converted, "roundtrip"))
    metrics.update(
        {
            "reference_camera_events": len(a),
            "candidate_camera_events": len(b),
            "camera_event_id_delta": len(set(a) ^ set(b)),
        }
    )
    metrics.update(compare_camera(a, b, root))
    return native, metrics, fields, a, b


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--showers", type=int, default=20)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--energy-min", type=float, default=100)
    parser.add_argument("--energy-max", type=float, default=3000)
    args = parser.parse_args()
    root = Path.cwd().resolve()
    source = (
        args.input.resolve()
        if args.input
        else generate(
            root, args.showers, energy_min=args.energy_min, energy_max=args.energy_max
        )
    )
    _, metrics, fields, a, b = roundtrip(root, source)
    common.write_results(
        metrics,
        {
            "codec_per_field_normalized_error": fields,
            "reference_event_ids": list(a),
            "candidate_event_ids": list(b),
            "camera_event_keys": list(next(iter(a.values())).keys()) if a else [],
        },
    )


if __name__ == "__main__":
    main()
