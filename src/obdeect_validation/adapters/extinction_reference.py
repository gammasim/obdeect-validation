"""Check real sim_telarray extinction tables and weighted/thinned efficiencies."""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common


def table_reference(path: Path):
    """Read the documented sim_telarray table and interpolate optical depth independently."""
    lines = path.read_text(encoding="utf-8").splitlines()
    header = next(line for line in lines if "H2=" in line and "H1=" in line)
    h2_text, h1_text = header.split("H2=", 1)[1].split("H1=", 1)
    h2_km = float(h2_text.strip().split()[0].rstrip(","))
    heights_km = np.asarray(
        [float(value) for value in h1_text.split()], dtype=np.float64
    )
    rows = []
    for line in lines:
        fields = line.split()
        if fields and fields[0].isdigit() and len(fields) == len(heights_km) + 1:
            rows.append([int(fields[0]), *map(float, fields[1:])])
    data = np.asarray(rows, dtype=np.float64)
    wavelengths = data[:, 0].astype(np.int64)
    tau = data[:, 1:]

    def depth(wavelength, height_m):
        row = int(
            np.clip(int(wavelength), wavelengths[0], wavelengths[-1]) - wavelengths[0]
        )
        if height_m <= h2_km * 1000.0:
            return 0.0
        log_grid = np.log10(heights_km * 1000.0)
        reference_grid = np.concatenate(([np.log10(h2_km * 1000.0)], log_grid))
        reference_tau = np.concatenate(([0.0], tau[row]))
        return float(np.interp(np.log10(height_m), reference_grid, reference_tau))

    return h2_km * 1000.0, depth


def main() -> None:
    from obdeect.propagation import TransmissionTable, apply_extinction

    table_path = Path(
        os.environ.get(
            "EXTINCTION_TABLE",
            "/validation-config/atm_trans_2158_1_10_0_0_0.dat",
        )
    )
    raw_lines = table_path.read_text(encoding="utf-8").splitlines()
    ignored_rows = 0
    normalized_lines = []
    for line in raw_lines:
        fields = line.split()
        if fields and fields[0].isdigit() and int(fields[0]) >= 1000:
            ignored_rows += 1
        else:
            normalized_lines.append(line)
    # The sim_telarray table also has a 1000 nm row; the photon API supports integer 0..999 nm.
    normalized_path = Path("simtel_atm_transmission_0_999nm.dat")
    normalized_path.write_text("\n".join(normalized_lines) + "\n", encoding="utf-8")
    table = TransmissionTable(str(normalized_path), 2158.0)
    site_m, depth = table_reference(normalized_path)
    wavelengths = [300, 350, 400, 450, 500, 550, 600, 700, 850]
    heights = [3000.0, 4000.0, 6000.0, 12000.0, 20000.0]
    zeniths = [0.0, 20.0, 40.0, 60.0]
    elevations = [site_m, site_m + 22.0]
    paths_m = [0.0, 5.0]
    cases = [
        (wavelength, height, zenith, elevation, path)
        for wavelength in wavelengths
        for height in heights
        for zenith in zeniths
        for elevation in elevations
        for path in paths_m
    ]
    repetitions = 512
    group_ids = np.repeat(np.arange(len(cases), dtype=np.float32), repetitions)
    wavelength_col = np.repeat([case[0] for case in cases], repetitions).astype(
        np.float32
    )
    height_col = np.repeat([case[1] for case in cases], repetitions).astype(np.float32)
    zenith_col = np.repeat([case[2] for case in cases], repetitions)
    cx = np.sin(np.deg2rad(zenith_col)).astype(np.float32)
    cy = np.zeros_like(cx)
    columns = {
        "x": np.zeros_like(cx),
        "y": np.zeros_like(cx),
        "cx": cx,
        "cy": cy,
        "t": group_ids,
        "zem": height_col,
        "w": np.ones_like(cx),
        "wavelength": wavelength_col,
    }

    expected_by_group = []
    for wavelength, height, zenith, elevation, path in cases:
        tau_z = depth(wavelength, height)
        tau_tel = depth(wavelength, elevation)
        tau_focus = depth(wavelength, elevation + path) - tau_tel
        expected_by_group.append(
            math.exp(-abs(tau_z - tau_tel) / math.cos(math.radians(zenith)) - tau_focus)
        )
    expected = np.repeat(np.asarray(expected_by_group), repetitions)

    options = dict(
        table=table,
        envelope_wavelength=[200.0, 999.0],
        envelope_value=[1.0, 1.0],
        seed=20261113,
        event=7,
        telescope=3,
        first_index=0,
        wavelength_range=(300.0, 850.0),
    )
    actual_p = []
    actual_w = []
    expected_actual = []
    thinned_count = 0
    for elevation in elevations:
        for path in paths_m:
            selected_groups = [
                i
                for i, case in enumerate(cases)
                if case[3] == elevation and case[4] == path
            ]
            rows = np.isin(group_ids, selected_groups)
            subset = {name: value[rows] for name, value in columns.items()}
            weighted = apply_extinction(
                subset,
                weighted=True,
                telescope_altitude_m=elevation,
                pathlen_m=path,
                **options,
            )
            thinned = apply_extinction(
                subset,
                weighted=False,
                telescope_altitude_m=elevation,
                pathlen_m=path,
                **options,
            )
            actual_p.append(weighted["p"].astype(np.float64))
            actual_w.append(weighted["w"].astype(np.float64))
            expected_actual.append(expected[rows])
            thinned_count += len(thinned["x"])
            if not np.all(thinned["w"] == 1.0):
                raise AssertionError("thinned photons must carry unit weight")
    candidate_p = np.concatenate(actual_p)
    candidate_w = np.concatenate(actual_w)
    expected = np.concatenate(expected_actual)
    weighted_error = float(np.max(np.abs(candidate_p - expected)))
    weight_error = float(np.max(np.abs(candidate_w - expected)))
    expected_count = float(expected.sum())
    variance = float(np.sum(expected * (1.0 - expected)))
    thinning_sigma = abs(thinned_count - expected_count) / math.sqrt(variance)

    # Check telescope elevation/path corrections in the deterministic weighted calculation.
    telescope_error = 0.0
    for path in paths_m:
        actual = table.telescope_depths(site_m + 22.0, path)
        for wavelength in wavelengths:
            reference_tel = depth(wavelength, site_m + 22.0)
            reference_focus = depth(wavelength, site_m + 22.0 + path) - reference_tel
            telescope_error = max(
                telescope_error,
                abs(actual[0][wavelength] - reference_tel),
                abs(actual[1][wavelength] - reference_focus),
            )

    metrics = {
        "wavelength_count": len(wavelengths),
        "height_count": len(heights),
        "zenith_count": len(zeniths),
        "telescope_elevation_count": len(elevations),
        "path_length_count": len(paths_m),
        "source_rows_ignored_above_api_wavelength_range": ignored_rows,
        "weighted_photon_count": len(candidate_w),
        "thinned_photon_count": thinned_count,
        "weighted_table_max_absolute_error": weighted_error,
        "weighted_column_max_absolute_error": weight_error,
        "telescope_depth_max_absolute_error": telescope_error,
        "thinning_count_residual_sigma": thinning_sigma,
    }
    common.write_json(
        "extinction-samples.json",
        {
            "source_table": str(table_path),
            "table_used": str(normalized_path),
            "ignored_rows_above_999_nm": ignored_rows,
            "site_altitude_m": site_m,
            "cases": [
                {
                    "wavelength_nm": w,
                    "emission_height_m": h,
                    "zenith_deg": z,
                    "telescope_altitude_m": e,
                    "telescope_path_m": p,
                    "transmission": expected_by_group[i],
                }
                for i, (w, h, z, e, p) in enumerate(cases)
            ],
        },
    )
    _plot(
        expected,
        candidate_p,
        [c[0] for c in cases],
        repetitions,
        Path("extinction-weighted.png"),
    )
    common.write_results(metrics)


def _plot(expected, weighted, wavelengths, repetitions, path):
    """obdeect weighted transmission against the table expectation, coloured by wavelength."""
    from obdeect_validation.adapters import plotting

    # Rows follow the (elevation, path) grouping, so wavelengths repeat per group.
    colors = np.resize(np.repeat(wavelengths, repetitions), expected.size).astype(float)
    fig, axes = plotting.figure(1, 2, (10, 3.8))
    points = axes[0, 0].scatter(expected, weighted, c=colors, s=6, cmap="viridis")
    axes[0, 0].plot([0, 1], [0, 1], color=plotting.ANALYTIC, linewidth=1)
    axes[0, 0].set(xlabel="table expectation", ylabel="obdeect weighted p")
    fig.colorbar(points, ax=axes[0, 0], label="wavelength [nm]")
    axes[0, 1].hist(weighted - expected, bins=50, color=plotting.NATIVE, alpha=0.85)
    axes[0, 1].set(xlabel="obdeect - expectation", ylabel="photons")
    plotting.save(fig, path)


if __name__ == "__main__":
    main()
