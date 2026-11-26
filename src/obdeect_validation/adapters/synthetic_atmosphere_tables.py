"""Cross-check synthetic simtools atmosphere tables against obdeect and analytic depth."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common


def analyse(root: Path) -> tuple[dict, None]:
    import astropy.units as u
    from obdeect.propagation import TransmissionTable, in_telescope_path_m
    from obdeect.ray_tracing import OpticalModel
    from simtools.atmosphere.optical_tables import (
        AtmosphericOpticalProperties,
        get_backend,
        write_obdeect_transmission_table,
        write_optical_properties,
    )

    # A small deterministic ECSV grid with exponential vertical optical depth.
    wavelength = np.array([400.0, 600.0])
    altitude = np.array([1000.0, 2000.0, 4000.0, 8000.0])
    amplitude = np.array([0.36, 0.24])
    scale_height = 5000.0
    tau = amplitude[:, None] * np.exp(-altitude[None, :] / scale_height)
    properties = AtmosphericOpticalProperties(
        wavelength * u.nm,
        altitude * u.m,
        0.7 * tau,
        0.3 * tau,
        np.zeros_like(tau),
        np.full_like(tau, 0.9),
        np.full_like(tau, 0.2),
        metadata={"backend": "synthetic-exponential", "version": "1"},
    )
    source = root / "synthetic-properties.ecsv"
    write_optical_properties(properties, source)
    generated = get_backend("ecsv", file_name=source).compute(
        wavelength * u.nm, altitude * u.m
    )
    if not np.allclose(generated.total_optical_depth(), tau, rtol=1e-12, atol=1e-14):
        raise AssertionError(
            "simtools ECSV backend changed the synthetic optical depths"
        )

    observation_level = 1000.0
    transmission_file = root / "synthetic-transmission.dat"
    write_obdeect_transmission_table(
        generated, transmission_file, observation_level * u.m
    )
    table = TransmissionTable(str(transmission_file), observation_level)

    # obdeect's table interpolates depth linearly in log altitude. Compute that
    # independently from the source grid, including telescope path inside optics.
    def depth(wl: float, height: float) -> float:
        """Independent log-height interpolation of generated table depths."""
        row = int(np.where(wavelength == wl)[0][0])
        heights = np.concatenate(([observation_level], altitude[1:]))
        values = np.concatenate(([0.0], tau[row, 0] - tau[row, 1:]))
        return float(np.interp(np.log(height), np.log(heights), values))

    samples = []
    for wl in wavelength:
        for height in (2000.0, 4000.0, 8000.0):
            vertical = math.exp(-depth(wl, height))
            for zenith in (0.0, 40.0, 60.0):
                airmass = 1.0 / math.cos(math.radians(zenith))
                expected = vertical**airmass
                actual = float(table.transmission(wl, height, airmass))
                samples.append((actual, expected))

    pathlen_m = 28.0
    depths = table.telescope_depths(observation_level, pathlen_m)
    expected_focus_depth = depth(400.0, observation_level + pathlen_m)
    expected_telescope_depth = 0.0
    error = max(abs(actual - expected) for actual, expected in samples)
    telescope_error = max(
        abs(float(depths[0][400]) - expected_telescope_depth),
        abs(float(depths[1][400]) - expected_focus_depth),
    )
    if error > 2e-6 or telescope_error > 2e-6:
        raise AssertionError(
            f"synthetic atmosphere mismatch: transmission={error:g}, telescope={telescope_error:g}"
        )

    # These constructors and expected lengths mirror the focused path-length
    # unit tests in obdeect, while keeping the validation case self-contained.
    def record(name, value, unit=None, *, file=False):
        return {
            "parameter": name,
            "parameter_version": "1.0.0",
            "value": value,
            "unit": unit,
            "file": file,
        }

    single_records = {
        "mirror_class": record("mirror_class", 0),
        "focal_length": record("focal_length", 800.0, "cm"),
        "mirror_offset": record("mirror_offset", 0.0, "cm"),
        "dish_shape_length": record("dish_shape_length", 800.0, "cm"),
        "parabolic_dish": record("parabolic_dish", False),
        "mirror_focal_length": record("mirror_focal_length", 0.0, "cm"),
        "mirror_list": record("mirror_list", "mirrors.dat", file=True),
        "camera_body_diameter": record("camera_body_diameter", 60.0, "cm"),
        "camera_body_shape": record("camera_body_shape", 2),
        "camera_depth": record("camera_depth", 0.0, "cm"),
    }
    mirror_list = """# x y diameter focal shape z (z = 0: automatic height)
   0.0     0.0  100.0  0.0  0  0.0
 125.0     0.0  100.0  0.0  0  0.0
-125.0     0.0  100.0  0.0  0  0.0
"""
    dual_records = {
        "mirror_class": record("mirror_class", 2),
        "focal_length": record("focal_length", 250.0, "cm"),
        "mirror_offset": record("mirror_offset", 0.0, "cm"),
        "primary_mirror_parameters": record("primary_mirror_parameters", [0.0], ["cm"]),
        "primary_mirror_diameter": record("primary_mirror_diameter", 400.0, "cm"),
        "primary_mirror_hole_diameter": record(
            "primary_mirror_hole_diameter", 0.0, "cm"
        ),
        "secondary_mirror_parameters": record(
            "secondary_mirror_parameters", [300.0], ["cm"]
        ),
        "secondary_mirror_diameter": record("secondary_mirror_diameter", 100.0, "cm"),
        "secondary_mirror_hole_diameter": record(
            "secondary_mirror_hole_diameter", 0.0, "cm"
        ),
        "focal_surface_parameters": record("focal_surface_parameters", [250.0], ["cm"]),
        "camera_body_diameter": record("camera_body_diameter", 50.0, "cm"),
        "camera_body_shape": record("camera_body_shape", 2),
        "camera_depth": record("camera_depth", 0.0, "cm"),
    }
    models = (
        (
            "single",
            OpticalModel.from_records(
                single_records,
                {"mirror_list": mirror_list}.__getitem__,
                configuration_seed=1,
                propagation_group_index=1.0,
            ),
            8.0,
        ),
        (
            "dual",
            OpticalModel.from_records(
                dual_records,
                {}.__getitem__,
                configuration_seed=1,
                propagation_group_index=1.0,
            ),
            3.5,
        ),
    )
    path_errors = []
    path_transmission_errors = []
    model_samples = []
    for model_name, model, analytic_path_m in models:
        model_path_m = float(in_telescope_path_m(model))
        path_errors.append(abs(model_path_m - analytic_path_m))
        od_tel_values, od_focus_values = table.telescope_depths(
            observation_level, model_path_m
        )
        for wl in wavelength:
            airmass = 1.0 / math.cos(math.radians(40.0))
            actual = float(
                table.transmission2(
                    wl,
                    8000.0,
                    airmass,
                    od_tel_values[int(wl)],
                    od_focus_values[int(wl)],
                )
            )
            tel_depth = depth(wl, observation_level)
            focus_depth = depth(wl, observation_level + model_path_m)
            shower_depth = depth(wl, 8000.0)
            expected = math.exp(-abs(shower_depth - tel_depth) * airmass - focus_depth)
            path_transmission_errors.append(abs(actual - expected))
            model_samples.append(
                {
                    "model": model_name,
                    "path_length_m": model_path_m,
                    "path_transmission": actual,
                    "path_transmission_reference": expected,
                }
            )
    path_error = max(path_errors)
    path_transmission_error = max(path_transmission_errors)
    if path_error > 1e-9 or path_transmission_error > 2e-6:
        raise AssertionError(
            f"optical path mismatch: length={path_error:g}, transmission={path_transmission_error:g}"
        )

    common.write_json(
        root / "samples.json",
        [
            *(
                {"candidate": actual, "reference": expected}
                for actual, expected in samples
            ),
            *model_samples,
        ],
    )
    return {
        "wavelength_count": len(wavelength),
        "height_count": 3,
        "zenith_count": 3,
        "transmission_max_absolute_error": error,
        "telescope_path_depth_max_absolute_error": telescope_error,
        "optical_model_count": len(models),
        "optical_model_path_max_absolute_error_m": path_error,
        "path_extinction_max_absolute_error": path_transmission_error,
    }, None


def main() -> None:
    common.results_main(analyse, __doc__)


if __name__ == "__main__":
    main()
