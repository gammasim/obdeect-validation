"""Compare the native LED stream writer with sim_telarray pixled."""

import os
import shutil
import subprocess
from importlib.metadata import version
from pathlib import Path

import numpy as np
from eventio.iact import IACTFile
from obdeect.stream import Reader

from obdeect_validation.adapters import common

FIXTURE = Path(__file__).with_name("compare-led-camera.txt")


def read_stream(path):
    event_counts = []
    x_weight = 0.0
    y_weight = 0.0
    total_weight = 0.0
    for record in Reader(str(path)):
        if record["type"] == "event":
            event_counts.append(0.0)
        elif record["type"] == "photons":
            if not event_counts:
                raise ValueError("photon stream has photons before its first event")
            columns = record["columns"]
            weights = columns["w"].astype(np.float64)
            event_counts[-1] += float(np.sum(weights))
            x_weight += float(np.dot(columns["x"].astype(np.float64), weights))
            y_weight += float(np.dot(columns["y"].astype(np.float64), weights))
            total_weight += float(np.sum(weights))
    if total_weight == 0:
        raise ValueError("obdeect produced no photons")
    return (
        event_counts,
        total_weight,
        100.0 * x_weight / total_weight,
        100.0 * y_weight / total_weight,
    )


def read_reference(path):
    event_counts = []
    x_weight = 0.0
    y_weight = 0.0
    total_weight = 0.0
    with IACTFile(path) as file:
        for event in file:
            event_weight = 0.0
            for bunches in event.photon_bunches.values():
                weights = np.asarray(bunches["photons"], dtype=np.float64)
                event_weight += float(np.sum(weights))
                x_weight += float(
                    np.dot(np.asarray(bunches["x"], dtype=np.float64), weights)
                )
                y_weight += float(
                    np.dot(np.asarray(bunches["y"], dtype=np.float64), weights)
                )
                total_weight += float(np.sum(weights))
            event_counts.append(event_weight)
    if total_weight == 0:
        raise ValueError("pixled produced no photons")
    return event_counts, total_weight, x_weight / total_weight, y_weight / total_weight


def main():
    obdeect = shutil.which("obdeect")
    sim_telarray = os.environ.get("SIMTOOLS_SIM_TELARRAY_PATH")
    pixled = shutil.which("pixled")
    if pixled is None and sim_telarray:
        candidate = Path(sim_telarray) / "LightEmission/pixled"
        if candidate.is_file():
            pixled = str(candidate)
    config = os.environ.get("SIMTEL_CONFIG_PATH")
    if obdeect is None or pixled is None:
        raise RuntimeError("obdeect and pixled must be installed")
    if not config or not Path(config).is_dir():
        raise RuntimeError(
            "SIMTEL_CONFIG_PATH must name the sim_telarray configuration directory"
        )
    if not (Path(config) / "atm_profile_model_1.dat").is_file():
        raise RuntimeError("SIMTEL_CONFIG_PATH must contain atm_profile_model_1.dat")

    camera = Path.cwd() / "camera.txt"
    native_output = Path.cwd() / "sources.gst"
    reference_output = Path.cwd() / "reference.dat.gz"
    shutil.copyfile(FIXTURE, camera)
    seed = 17
    shared = ["--use", "0 1 2", "--require", "2", "--photons", "10"]
    subprocess.run(
        [
            obdeect,
            "sources",
            "run",
            "led",
            "--camera-file",
            str(camera),
            *shared,
            "--seed",
            str(seed),
            "--observation-level",
            "2150",
            "--atmosphere-id",
            "1",
            "--spectrum",
            "fixed:wavelength_nm=400",
            "--camera-radius",
            "1.23",
            "-o",
            str(native_output),
        ],
        check=True,
    )
    process = subprocess.run(
        [
            pixled,
            "-I" + config,
            "-c",
            str(camera),
            *shared,
            "-o",
            str(reference_output),
        ],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    (Path.cwd() / "pixled.log").write_text(process.stdout, encoding="utf-8")
    process.check_returncode()
    native_metrics = read_stream(native_output)
    reference_metrics = read_reference(reference_output)
    event_deltas = common.event_deltas(native_metrics[0], reference_metrics[0])
    metrics = {
        "native_event_count": len(native_metrics[0]),
        "reference_event_count": len(reference_metrics[0]),
        "event_count_delta": abs(len(native_metrics[0]) - len(reference_metrics[0])),
        "native_photon_count": native_metrics[1],
        "reference_photon_count": reference_metrics[1],
        "photon_count_delta": abs(native_metrics[1] - reference_metrics[1]),
        "max_event_photon_delta": max(event_deltas),
        "position_centroid_delta_cm": float(
            np.hypot(
                native_metrics[2] - reference_metrics[2],
                native_metrics[3] - reference_metrics[3],
            )
        ),
    }
    profile = Path(config) / "atm_profile_model_1.dat"
    provenance = {
        "seed": seed,
        "events_native_photons": native_metrics[0],
        "events_reference_photons": reference_metrics[0],
        "camera_sha256": common.sha256(camera),
        "atmosphere_sha256": common.sha256(profile),
        "pixled_sha256": common.sha256(Path(pixled)),
        "obdeect_sources_version": version("obdeect-sources"),
        "obdeect_stream_version": version("obdeect-stream"),
    }
    common.write_results(metrics, provenance)


if __name__ == "__main__":
    main()
