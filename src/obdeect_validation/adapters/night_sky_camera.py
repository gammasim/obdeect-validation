"""Night-sky camera p.e. of obdeect against an analytic expectation (no nsbls).

Usage (inside the container, in an empty case directory):
    python -m obdeect_validation.adapters.night_sky_camera

Native: obdeect sources night-sky (the configuration of source_comparison.NIGHT: rate 2.4 ph /
(sr ns cm^2), 10 deg view cone, 1 ns window, 400 nm), propagated without extinction, converted
with obdeect stream to-simtel and run through sim_telarray with the night-sky camera options
(source_comparison.CAMERA_NIGHT). The atmosphere is zero on both sides (see zero_extinction_table).
The analytic expectation uses the same photon rate and the sim_telarray configuration of the LST-1
camera; no reference program is run.

Analytic expectation of the mean p.e. per event:

    N_pe = rate * dt * A_mirror * Omega_eff * R * D * LG * W * QE_high * T_tel

    rate        2.4 ph / (sr ns cm^2)                       (source_comparison.NIGHT)
    dt          1 ns, the pulse window of the configuration
    A_mirror    mirror area traced by sim_telarray with the 348 cm camera body and frame explicit,
                377.71 m^2 (CTA-PROD6-LST.cfg, comment at telescope_transmission). T_tel = 366.0 / 377.71
                is defined relative to this area. The projected area without the body, 388.38 m^2, would
                count the body shadow twice (earlier version: 2.8 % too high).
    Omega_eff   sum_i a_i q_i / f^2: the sky photons from the 10 deg cone that reach the focal plane
                land uniformly per unit focal-plane area; a_i is the funnel area of pixel i (hexagon,
                flat-to-flat 4.9 cm from the camera file), q_i its relative QE (camera file, mean 0.983),
                f = effective focal length 2930.565 cm (CTA-PROD6-LST.cfg, optical ray-tracing)
    R           mirror reflectivity at 400 nm (ref_LST1_2022_04_01.dat)
    D           mirror_degraded_reflection = 0.80 (CTA-PROD6-LST.cfg). sim_telarray multiplies the
                reflectivity by D in the quantum-efficiency table (sim_absorb.c, rpt_qe_ref); it is not
                only an NSB factor. Earlier versions left D out (0.80 too high in the expectation).
    LG          light-guide (funnel) efficiency averaged over the angles at which the camera receives
                the light: a facet at radius rho with focal length f sends its light to the focal plane
                at slope rho/f, so the funnel sees atan(rho/f) (sim_imaging.c, cathode_hit with the
                funnel angle table). Area-weighted over the 198 facets of the LST-1 mirror list. At normal
                incidence LG = 0.7535; the average is 0.8156 (1.082 x normal incidence). The earlier
                normal-incidence value was too low.
    W           camera window transmission at 400 nm (transmission_lst_window_No7-10_ave.dat)
    QE_high     PMT quantum efficiency at 400 nm, high-QE column (qe_lst1_20200318_high+low.dat)
    T_tel       telescope transmission (shadowing by the camera body, CTA-PROD6-LST.cfg)

Not in the expectation, on purpose: the atmosphere (zero, see below), the per-pixel gain and QE
random variations (zero mean), the photon sphere radius (cancels: the photon count per event is
rate x Omega x pi R^2 x dt and only the mirror-hitting fraction enters).

Atmosphere. The unit fixture (fixtures/unit-transmission.dat) has its bottom at 2.156 km. sim_telarray
returns zero transmission for emission below that bottom (atmospheric_transmission2, zem < h2), and the
source emits between 2147.5 m and 2165 m. Before this fix 16.6 % of the photon weight was removed by
the sim (photon statistics "incident photons on sphere after accounting for atm. extinction": 0.837 of
the photons on the sphere). zero_extinction_table moves the bottom to OBS_ALTITUDE_M = 2140 m and the
run uses the same observation level (-C altitude), so all emission lies above the bottom.

Check against the sim_telarray photon statistics (run of 2026-11-15, 4 events): the first-selection
efficiency per photon, sum(p1 p2) / sum(p1) = 0.26355, equals R * D * W * QE_high * T_tel = 0.26355
(LG enters later, in releff). The mirror-trace fraction of the photons is 0.529, i.e. 374 m^2 on the
15 m disk, within 1 % of A_mirror. Run with telescope_transmission = 1 (diagnostic, not in the case):
the sim gives 8931 +- 20 p.e./event, 3.5 % above the 8622 of the case, as 1/T_tel predicts. The
photon statistics are stored in details.json as photon_statistics.

Limit: |sim / expected - 1| <= max(5 %, 5 standard errors), set before the first run. The
expectation has no fitted parameter. The 5 % covers the geometry approximations that are not
modelled here: the hexagon footprint instead of the exact facet shadowing, the facet-to-focus slope
rho/f, and the focal-plane mapping with f = 2930.565 cm. The limit is the same as in the first version.

If sim_telarray gives fewer than two events, the case is BLOCKED and no camera metric is judged.
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common
from obdeect_validation.adapters import obdeect_cli as cli
from obdeect_validation.adapters import source_comparison as sc

CFG = common.REFERENCE / "cfg"
MIRROR_LIST = CFG / "CTA" / "mirror_CTA-N-LST1_v2019-03-31_rotated.dat"
CAMERA_FILE = CFG / "CTA" / "camera_CTA-LST-1_analogsum21_v2020-04-14.dat"
LST_CFG = CFG / "CTA" / "CTA-PROD6-LST.cfg"
LG_TABLE = CFG / "CTA" / "CTA-LST_lightguide_eff_2020-04-12_average.dat"
WAVELENGTH_NM = sc.NIGHT["spectrum_nm"]
ANALYTIC_LIMIT = 0.05
UNIT_TRANSMISSION = common.VALIDATION / "fixtures" / "unit-transmission.dat"
# Observation level of the camera run: below the lowest emission height of the night-sky photons (2147.5 m).
OBS_ALTITUDE_M = 2140.0


def _table_rows(path, column):
    """(x, column) rows of a table, sorted by x; '#' comments are skipped."""
    rows = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split("#")[0].split()
        if len(fields) > column and _is_number(fields[0]):
            rows.append((float(fields[0]), float(fields[column])))
    return np.array(sorted(rows))


def _table_value(path, column, wavelength):
    """Linear interpolation of one column of a two-column-or-more table (skips '#' lines)."""
    table = _table_rows(path, column)
    return float(np.interp(wavelength, table[:, 0], table[:, 1]))


def _is_number(text):
    try:
        float(text)
    except ValueError:
        return False
    return True


def _cfg_number(name):
    """Value of 'name = number' in the LST configuration (first match, comment stripped)."""
    text = LST_CFG.read_text(encoding="utf-8", errors="replace")
    match = re.search(rf"^\s*{re.escape(name)}\s*=\s*([-+0-9.eE]+)", text, re.MULTILINE)
    if match is None:
        raise ValueError(f"{name} not found in {LST_CFG.name}")
    return float(match.group(1))


def traced_area_cm2():
    """Mirror area traced by sim_telarray with the 348 cm camera body and frame (377.71 m^2)."""
    text = LST_CFG.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"([0-9.]+)\+-[0-9.]+ m\^2 with the 348 cm", text)
    if match is None:
        raise ValueError(
            "traced mirror area (348 cm camera body) not found in the LST configuration"
        )
    return float(match.group(1)) * 1e4


def projected_area_cm2():
    """Projected mirror area quoted by the configuration (388.38 m^2 in the mirror_list comment)."""
    text = LST_CFG.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"([0-9.]+)\s*m\^2 projected", text)
    if match is None:
        raise ValueError("projected mirror area not found in the LST configuration")
    return float(match.group(1)) * 1e4


def mirror_facets():
    """(radius cm, hexagon area cm^2, focal length cm) of the facets of the LST-1 mirror list."""
    facets = []
    for line in MIRROR_LIST.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split("#")[0].split()
        if len(fields) >= 4 and all(_is_number(f) for f in fields[:4]):
            x, y, flat, focal = (float(f) for f in fields[:4])
            facets.append((math.hypot(x, y), math.sqrt(3.0) / 2.0 * flat**2, focal))
    return np.asarray(facets)


def funnel_efficiency_average():
    """Light-guide efficiency averaged over the incidence angles of the light on the camera.

    A facet at radius rho with focal length f sends its light to the focal plane at slope rho/f,
    so the funnel sees the angle atan(rho/f) (sim_imaging.c, funnel angle table at the slope of the
    photon). Area-weighted over the facets. Returns (average, normal incidence).
    """
    facets = mirror_facets()
    table = _table_rows(LG_TABLE, 1)
    angle = np.degrees(np.arctan(facets[:, 0] / facets[:, 2]))
    eff = np.interp(angle, table[:, 0], table[:, 1])
    average = float(np.sum(facets[:, 1] * eff) / np.sum(facets[:, 1]))
    return average, float(np.interp(0.0, table[:, 0], table[:, 1]))


def hexagon_sum_area_cm2():
    """Sum of hexagon areas (flat-to-flat diameter, column 3 of the mirror list), cm^2."""
    total = 0.0
    for line in MIRROR_LIST.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split("#")[0].split()
        if len(fields) >= 3 and _is_number(fields[2]):
            total += math.sqrt(3.0) / 2.0 * float(fields[2]) ** 2
    return total


def funnel_area_and_qe():
    """Funnel area (hexagon, flat-to-flat from the camera file, cm^2) and the relative QE of each
    enabled pixel of the LST-1 camera layout (camera_pixel_layout model, columns enabled, relative_qe)."""
    funnel = None
    for line in CAMERA_FILE.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("PixType"):
            fields = line.split()
            # PixType 1 0 0 2.500 3 4.900 5.100 "file": funnel shape 3 (hexagon), diameter 4.9 cm
            funnel = (
                math.sqrt(3.0) / 2.0 * float(fields[6]) ** 2
                if fields[5] == "3"
                else None
            )
            break
    if funnel is None:
        raise ValueError("camera file: hexagon funnel not found")
    layout = sorted(
        (
            common.models_root()
            / "model_parameters"
            / "LSTN-01"
            / "camera_pixel_layout"
        ).glob("*.ecsv")
    )
    if not layout:
        raise FileNotFoundError("LSTN-01 camera_pixel_layout model not found")
    enabled, relative_qe = [], []
    for line in layout[-1].read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if line.startswith("#") or len(fields) < 10 or not _is_number(fields[0]):
            continue
        enabled.append(float(fields[8]))
        relative_qe.append(float(fields[9]))
    enabled = np.asarray(enabled)
    relative_qe = np.asarray(relative_qe)
    if enabled.size != 1855:
        raise ValueError(f"camera layout has {enabled.size} pixels, expected 1855")
    return funnel, enabled, relative_qe


def analytic_expectation():
    """Analytic mean p.e. per event and its factors (see the module docstring)."""
    funnel_cm2, enabled, relative_qe = funnel_area_and_qe()
    focal_cm = _cfg_number("effective_focal_length")
    telescope_transmission = _cfg_number("telescope_transmission")
    mirror_cm2 = traced_area_cm2()
    reflectivity = _table_value(
        CFG / "CTA" / "ref_LST1_2022_04_01.dat", 1, WAVELENGTH_NM
    )
    mirror_degradation = _cfg_number("mirror_degraded_reflection")
    lightguide, lightguide_normal = funnel_efficiency_average()
    window = _table_value(
        CFG / "CTA" / "transmission_lst_window_No7-10_ave.dat", 1, WAVELENGTH_NM
    )
    qe_high = _table_value(
        CFG / "CTA" / "qe_lst1_20200318_high+low.dat", 1, WAVELENGTH_NM
    )
    omega_eff = float(np.sum(funnel_cm2 * enabled * relative_qe) / focal_cm**2)
    window_ns = sc.NIGHT["pulse_fwhm_ns"] / sc.FWHM_TO_SIGMA
    hits = sc.NIGHT["rate_per_sr_ns_cm2"] * window_ns * mirror_cm2 * omega_eff
    efficiency = (
        reflectivity
        * mirror_degradation
        * lightguide
        * window
        * qe_high
        * telescope_transmission
    )
    return {
        "expected_pe_per_event": float(hits * efficiency),
        "expected_photons_on_camera_per_event": float(hits),
        "efficiency_reflectivity": reflectivity,
        "efficiency_mirror_degradation": mirror_degradation,
        "efficiency_lightguide": lightguide,
        "efficiency_lightguide_normal_incidence": lightguide_normal,
        "efficiency_window": window,
        "efficiency_qe_high": qe_high,
        "efficiency_telescope_transmission": telescope_transmission,
        "efficiency_total": float(efficiency),
        "mirror_area_traced_m2": mirror_cm2 * 1e-4,
        "mirror_area_projected_m2": projected_area_cm2() * 1e-4,
        "mirror_area_hexagon_sum_m2": hexagon_sum_area_cm2() * 1e-4,
        "omega_eff_sr": omega_eff,
        "focal_length_cm": focal_cm,
        "funnel_area_cm2": funnel_cm2,
        "pixels_enabled": int(enabled.sum()),
        "relative_qe_mean": float(np.mean(relative_qe)),
    }


def zero_extinction_table(work):
    """Unit transmission table with its bottom moved to OBS_ALTITUDE_M (shared zero-depth helper).

    The fixture starts at H2 = 2.156 km. sim_telarray returns zero transmission for emission below
    the bottom of the table, so the photons emitted between 2147.5 m and 2156 m were removed.
    """
    return common.write_zero_depth_table(
        work,
        UNIT_TRANSMISSION,
        bottom_m=OBS_ALTITUDE_M,
        name=f"unit-transmission-{int(OBS_ALTITUDE_M)}m.dat",
    )


def photon_statistics(log):
    """Photon statistics of sim_telarray for telescope 1 (sums over the events of the run)."""
    text = Path(log).read_text(encoding="utf-8", errors="replace")
    if "Photon statistics:" not in text:
        return {}
    block = text.split("Photon statistics:", 1)[1].split("Trigger statistics:", 1)[0]
    stats = {}
    for line in block.splitlines():
        if ":" in line:
            label, values = line.split(":", 1)
            numbers = [float(v) for v in values.split() if _is_number(v)]
            if numbers:
                stats[label.strip()] = numbers[0]
    return stats


def run(work):
    """Simulate the night-sky camera input; returns the p.e. per event of the sim_telarray run."""
    native_file = work / "native_sources.gst"
    radius, viewcone = sc.NIGHT["radius_m"], sc.NIGHT["viewcone_deg"]
    window_ns = sc.NIGHT["pulse_fwhm_ns"] / sc.FWHM_TO_SIGMA
    photons = round(
        sc.night_sky_photons(
            sc.NIGHT["rate_per_sr_ns_cm2"], viewcone, radius, window_ns
        )
    )
    argv = cli.obdeect(
        "sources",
        "run",
        "night-sky",
        observation_level=str(sc.NIGHT["altitude_m"]),
        atmosphere_id="1",
        seed=str(sc.NIGHT["seed"]),
        events=str(sc.NIGHT["events"]),
        photons=str(photons),
        spectrum=f"fixed:wavelength_nm={sc.NIGHT['spectrum_nm']}",
        pulse=sc.NIGHT_PULSE,
        bunch_size=str(sc.NIGHT["bunch_size"]),
        viewcone_deg=str(viewcone),
        telescope_position="0,0,0",
        telescope_radius=str(radius),
        executable=common.tool("obdeect"),
        o=str(native_file),
    )
    common.run_tool(argv, work / "native-sources.log")
    camera_input = sc.propagate_to_simtel(
        work, native_file, sc.NIGHT["altitude_m"], sc.NIGHT["band_nm"], sc.NIGHT["seed"]
    )
    # Zero optical depth on both sides: the expectation has no atmosphere factor. The observation
    # level is set below all emission heights (later -C options override the defaults of the run).
    table = zero_extinction_table(work)
    extra = [
        *sc.CAMERA_NIGHT,
        "-C",
        f"altitude={int(OBS_ALTITUDE_M)}",
        "-C",
        f"atmospheric_transmission={table}",
    ]
    output = sc.simulate(work, camera_input, "camera-native", extra=extra)
    return sc.camera_pe(output), photons


def main() -> None:
    common.results_main(_analyse, __doc__.splitlines()[0])


def _analyse(root: Path):
    work = root.resolve()
    values, photons = run(work)
    expected = analytic_expectation()
    metrics = {
        "camera_events": int(values.size),
        "photons_per_event_configured": photons,
        "camera_pe_per_event_expected": expected["expected_pe_per_event"],
        "expected_photons_on_camera_per_event": expected[
            "expected_photons_on_camera_per_event"
        ],
        "efficiency_total": expected["efficiency_total"],
        "mirror_area_hexagon_relative_difference": expected[
            "mirror_area_hexagon_sum_m2"
        ]
        / expected["mirror_area_projected_m2"]
        - 1.0,
    }
    details = {
        "case": "night-sky camera (analytic expectation, no nsbls)",
        "expectation": expected,
        "photon_statistics": photon_statistics(work / "camera-native.log"),
        "limit": f"|sim/expected - 1| <= max({ANALYTIC_LIMIT}, 5 se)",
        "camera_pe_native": values.tolist(),
    }
    if values.size < 2:
        details["status"] = "BLOCKED: fewer than two sim_telarray events"
    else:
        mean = float(values.mean())
        se = float(values.std(ddof=1) / math.sqrt(values.size))
        ratio = mean / expected["expected_pe_per_event"]
        ratio_se = se / expected["expected_pe_per_event"]
        limit = max(ANALYTIC_LIMIT, 5.0 * ratio_se)
        metrics.update(
            {
                "camera_pe_per_event_native": mean,
                "camera_pe_per_event_native_se": se,
                "camera_pe_ratio_to_expected": ratio,
                "camera_pe_ratio_to_expected_se": ratio_se,
                "camera_pe_ratio_limit_exceeded": int(abs(ratio - 1.0) > limit),
            }
        )
        details["status"] = "run"
    _plot(values, expected, root / "night-sky-camera.png")
    return metrics, details


def _plot(values, expected, path) -> None:
    from obdeect_validation.adapters import plotting

    fig, axis = plotting.figure(1, 1, (6, 3.6))
    if values.size:
        axis.hist(
            values,
            bins=12,
            color=plotting.NATIVE,
            alpha=0.8,
            label="obdeect, sim_telarray",
        )
    axis.axvline(
        expected["expected_pe_per_event"],
        color=plotting.ANALYTIC,
        linestyle="--",
        label="analytic expectation",
    )
    axis.set(xlabel="p.e. per event", ylabel="events")
    axis.legend()
    plotting.save(fig, path)


if __name__ == "__main__":
    sys.exit(main())
