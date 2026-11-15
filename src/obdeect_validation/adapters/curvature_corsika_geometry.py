"""Check CORSIKA's CURVED-option assumptions on real IACT files (independent of obdeect).

obdeect's corsika mode takes CURVED-CORSIKA photons as already curved: emission altitudes zem
above sea level, directions and positions in the detector frame, times and refraction applied by
CORSIKA. This adapter checks those assumptions on the two curved shower files (65 and 75 deg) with numpy and the
EventIO reader only. No obdeect code is used.

Per file:
  (a) Emission altitude. The photon-weighted median zem of the file is placed on the primary
      axis (zenith from EVTH word 11, the curved_showers reader) and its slant depth is computed
      in the spherical atmosphere of the atmosphere profile (Earth radius R, straight axis through
      the core, array reference at the observation level). The curved slant depth must lie in
      XMAX_BAND. The plane-parallel value thick(zem)/cos(zenith) is reported for contrast.
  (b) Directions. The photon-weighted mean zenith of the bunches within CORE_RADIUS_M of the core
      must equal the primary zenith within ZENITH_TOLERANCE_DEG.
  (c) Bunch format. Every photon block (EventIO type 1205) must be in the long format, i.e. its
      header version is not 1000 + n (the compact format, io_simtel.c). The long format is the
      default of corsika7/bernlohr/iact.c and is also requested by TELFIL '+' (curved_showers).
  (d) Header and log. The CURVED log banner and the EVTH flags (word 79 and word 77 bit 0x40)
      must agree and mark the file as curved.

Inputs are the shower files named by curved_showers.shower_name under VALIDATION/inputs, each with
its <name>.corsika.log. A missing input ends the case with exit code 3 and a message.
"""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path

import numpy as np
from eventio import EventIOFile
from eventio.iact import IACTFile
from numpy.polynomial.legendre import leggauss

from obdeect_validation.adapters import common, curved_showers

EARTH_RADIUS_M = 6371315.0
# Observation level of the cards (OBSLEV 215600 cm); the run header of each file is checked too.
SITE_ALTITUDE_M = 2156.0
GAUSS_ORDER = (
    32  # nodes per profile segment; the density is exponential within a segment
)

# Limits fixed before looking at results. Each entry: (rule, justification).
XMAX_BAND_G_CM2 = (200.0, 500.0)
ZENITH_TOLERANCE_DEG = 0.5
CORE_RADIUS_M = 150.0
PHOTON_BLOCK_TYPE = (
    1205  # EventIO type of an IACT photon bunch block (eventio iact.objects)
)
COMPACT_VERSION_BASE = (
    1000  # header version of a compact block is 1000 + n (io_simtel.c)
)

LIMITS = {
    "xmax_band_g_cm2": (
        "200 <= X_slant(median zem) <= 500 g/cm2",
        (
            "Shower maximum of a gamma shower: X_max ~ X0 (ln(E/E_c) - 0.5) with X0 = 36.7 g/cm2 "
            "and E_c = 86 MeV (air). For 300 GeV, 1 TeV, 5 TeV and 10 TeV this gives about 280, "
            "325, 385 and 410 g/cm2 (without the -0.5 term: 300 to 430). The band [200, 500] leaves "
            "a margin of about 30 % below the lowest and about 15 to 20 % above the highest "
            "prediction, covering the spread of the photon emission profile around X_max. It is a "
            "plausibility window, not a precision test. The plane-parallel value is not used for "
            "the metric."
        ),
    ),
    "zenith_tolerance_deg": (
        "|mean zenith(|x_core| <= 150 m) - theta0| <= 0.5 deg",
        (
            "Refraction bends the arrival direction by about (n-1) tan(theta) = 2.8e-4 * 3.73 = "
            "0.06 deg at 75 deg. The Cherenkov cone (about 1 deg half-angle) is symmetric about the "
            "axis, so its effect on the mean zenith is second order, alpha^2 cot(theta) / 2 (0.004 "
            "deg at 65 deg, 0.002 deg at 75 deg). 0.5 deg is about 8 times the refraction effect at "
            "75 deg."
        ),
    ),
    "core_radius_m": (
        "150 m",
        (
            "Fixed a priori (not tuned to results): the bunches of the core region, where the "
            "light pool is dense and the arrival directions are most sensitive to the shower axis. "
            "Impact points are corrected to the ground (sphere-centre plane offset), see "
            "ground_impact_m."
        ),
    ),
    "long_format_violation_count": (
        "0 compact photon blocks (EventIO version 1000 + n) in each file, at least one block",
        (
            "The long bunch format has float32 fields without the int16 range limits of the compact "
            "format (bernlohr README section d), so the limits no longer apply. The format is "
            "checked on the file itself, not inferred from the TELFIL name: the block version is "
            "0 for the long and 1000 + n for the compact format (io_simtel.c)."
        ),
    ),
    "header_flags": (
        "word 79 == 2, word 77 bit 0x40 set, CURVED log banner present",
        "A CURVED run must carry all three; a standard binary must not pass as curved.",
    ),
}


# ---- atmosphere profile -------------------------------------------------------------------


def read_profile(path) -> dict:
    """Atmosphere profile (sim_telarray atmprof format): altitude [km], density [g/cm3], thickness.

    Returns arrays h_m (altitudes in m above sea level), rho (g/cm3) and thick_table (g/cm2, as
    tabulated). Rows need at least three columns; '#' lines are comments.
    """
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        fields = text.split()
        if len(fields) < 3:
            raise ValueError(f"atmosphere profile row needs three columns: {line!r}")
        rows.append([float(value) for value in fields[:3]])
    table = np.asarray(rows, dtype=np.float64)
    if table.shape[0] < 2:
        raise ValueError(f"atmosphere profile {path} has fewer than two rows")
    h_m = table[:, 0] * 1000.0
    rho = table[:, 1]
    if np.any(np.diff(h_m) <= 0.0) or np.any(rho <= 0.0):
        raise ValueError(
            f"atmosphere profile {path}: altitudes must rise, densities be positive"
        )
    return {"h_m": h_m, "rho": rho, "thick_table": table[:, 2]}


def _segments(profile):
    """Log-linear density segments (h0, h1, rate per m, rho at h0): exact for exponential layers."""
    h, rho = profile["h_m"], profile["rho"]
    rate = np.log(rho[1:] / rho[:-1]) / np.diff(h)
    return list(zip(h[:-1], h[1:], rate, rho[:-1]))


def _column(profile, h_low, geometry=None) -> float:
    """Integral of rho(h) * geometry(h) dh from h_low to the profile top, in g/cm2.

    geometry is a function of altitude (metres) or None for the vertical column. Each segment
    uses Gauss-Legendre nodes; the density is exponential within a segment.
    """
    bottom, top = profile["h_m"][0], profile["h_m"][-1]
    if h_low < bottom:
        raise ValueError(
            f"altitude {h_low} m lies below the atmosphere profile ({bottom} m)"
        )
    if h_low >= top:
        return 0.0
    nodes, weights = leggauss(GAUSS_ORDER)
    total = 0.0
    for h0, h1, rate, rho0 in _segments(profile):
        lo = max(h_low, h0)
        if lo >= h1:
            continue
        half = 0.5 * (h1 - lo)
        altitude = lo + half * (nodes + 1.0)
        density = rho0 * np.exp(rate * (altitude - h0))
        factor = 1.0 if geometry is None else geometry(altitude)
        total += float(np.sum(weights * half * density * factor))
    return total * 100.0  # g/cm3 times m -> g/cm2


def thick_from_density(profile, h_m: float) -> float:
    """Vertical column above altitude h_m in g/cm2, integrated from the density profile."""
    return _column(profile, float(h_m))


def slant_depth_curved(
    profile,
    h_m: float,
    theta0_deg: float,
    earth_radius_m: float = EARTH_RADIUS_M,
    site_altitude_m: float = SITE_ALTITUDE_M,
) -> float:
    """Slant depth in g/cm2 from altitude h_m up along the primary axis, spherical atmosphere.

    The axis is the straight line through the core, which lies at site_altitude_m above the
    Earth's surface (sphere of radius earth_radius_m, centre below the core). Points on it are at
    radius r with (s + a cos z)^2 = r^2 - a^2 sin^2 z, where a = R + h_site, z = theta0 and s the
    path length from the core. Hence ds/dh = r / sqrt(r^2 - a^2 sin^2 z) and

        X(h) = integral from h to the top of rho(h') ds/dh' dh'.

    The vertical case (z = 0) gives thick(h); z -> 90 deg diverges only through the factor above.
    """
    h_m = float(h_m)
    zenith = math.radians(float(theta0_deg))
    if not 0.0 <= zenith < math.pi / 2.0:
        raise ValueError(f"primary zenith must be in [0, 90) deg, got {theta0_deg}")
    radius = float(earth_radius_m)
    site = radius + float(site_altitude_m)
    sin_z, cos_z = math.sin(zenith), math.cos(zenith)
    if sin_z == 0.0:
        return _column(profile, h_m)

    def geometry(altitude):
        # r^2 - a^2 sin^2 = (r - a)(r + a) + a^2 cos^2 avoids cancellation for large R.
        r = radius + altitude
        denominator = (altitude - site_altitude_m) * (r + site) + (site * cos_z) ** 2
        if np.any(denominator <= 0.0):
            raise ValueError(
                f"altitude {h_m} m is not on the primary axis for zenith {theta0_deg}"
            )
        return r / np.sqrt(denominator)

    return _column(profile, h_m, geometry)


def slant_depth_plane(profile, h_m: float, theta0_deg: float) -> float:
    """Plane-parallel slant depth thick(h)/cos(theta0) in g/cm2 (straight lines, flat Earth)."""
    zenith = math.radians(float(theta0_deg))
    if not 0.0 <= zenith < math.pi / 2.0:
        raise ValueError(f"primary zenith must be in [0, 90) deg, got {theta0_deg}")
    return _column(profile, float(h_m)) / math.cos(zenith)


# ---- pure helpers -------------------------------------------------------------------------


def photon_weighted_median(values, weights) -> float:
    """Smallest value whose cumulative weight reaches half of the total weight (no interpolation)."""
    values = np.asarray(values, dtype=np.float64).ravel()
    weights = np.asarray(weights, dtype=np.float64).ravel()
    if values.size == 0 or values.size != weights.size:
        raise ValueError("median needs equal, non-empty value and weight arrays")
    if np.any(weights < 0.0) or not np.sum(weights) > 0.0:
        raise ValueError("weights must be non-negative with a positive sum")
    order = np.argsort(values, kind="stable")
    cumulative = np.cumsum(weights[order])
    index = int(np.searchsorted(cumulative, 0.5 * cumulative[-1], side="left"))
    return float(values[order][index])


def mean_zenith_deg(cx, cy, weights) -> float:
    """Photon-weighted mean zenith angle in degrees of directions from direction cosines cx, cy.

    The downward cosine magnitude is sqrt(1 - cx^2 - cy^2); the zenith is arctan2(|c_t|, |c_z|),
    which is stable at large zenith angles.
    """
    cx = np.asarray(cx, dtype=np.float64).ravel()
    cy = np.asarray(cy, dtype=np.float64).ravel()
    weights = np.asarray(weights, dtype=np.float64).ravel()
    if cx.size == 0:
        raise ValueError("mean zenith needs at least one direction")
    transverse = np.hypot(cx, cy)
    longitudinal = np.sqrt(np.clip(1.0 - cx**2 - cy**2, 0.0, None))
    zenith = np.degrees(np.arctan2(transverse, longitudinal))
    total = float(np.sum(weights))
    if not total > 0.0:
        raise ValueError("mean zenith needs a positive total weight")
    return float(np.sum(weights * zenith) / total)


def ground_impact_m(x_rel_m, y_rel_m, cx, cy, z0_m, x_tel_m, y_tel_m):
    """Impact points at the observation level in the array frame, from the bunch positions.

    CORSIKA stores positions in the horizontal plane through the detector sphere centre (height
    z0 above the observation level): x_plane = x_ground - (cx/|cz|) z0 relative to the detector.
    Adding the detector position and the slope correction gives the ground impact point.
    """
    x_rel_m = np.asarray(x_rel_m, dtype=np.float64)
    y_rel_m = np.asarray(y_rel_m, dtype=np.float64)
    cx = np.asarray(cx, dtype=np.float64)
    cy = np.asarray(cy, dtype=np.float64)
    longitudinal = np.sqrt(np.clip(1.0 - cx**2 - cy**2, 1e-12, None))
    x_ground = x_rel_m + x_tel_m + cx / longitudinal * z0_m
    y_ground = y_rel_m + y_tel_m + cy / longitudinal * z0_m
    return x_ground, y_ground


def within_radius(x_m, y_m, centre_x_m, centre_y_m, radius_m) -> np.ndarray:
    """Boolean mask of points within radius_m (inclusive) of the centre."""
    x_m = np.asarray(x_m, dtype=np.float64)
    y_m = np.asarray(y_m, dtype=np.float64)
    return (x_m - centre_x_m) ** 2 + (y_m - centre_y_m) ** 2 <= radius_m**2


def is_compact_version(version) -> bool:
    """True if an IACT photon block header version marks the compact bunch format (1000 + n)."""
    return int(version) // COMPACT_VERSION_BASE == 1


def compact_block_count(versions) -> int:
    """Number of photon blocks, given their header versions, that use the compact format."""
    return sum(1 for version in versions if is_compact_version(version))


def _objects(parent):
    """All EventIO objects below parent (an iterable of objects), depth first."""
    for obj in parent:
        yield obj
        if obj.header.only_subobjects:
            yield from _objects(obj)


def photon_block_versions(path) -> list[int]:
    """Header versions of all photon bunch blocks (type 1205) of an IACT file, in file order.

    Raises RuntimeError when the file holds no photon block, so that an empty file cannot pass.
    """
    with EventIOFile(str(path)) as stream:
        versions = [
            int(obj.header.version)
            for obj in _objects(stream)
            if obj.header.type == PHOTON_BLOCK_TYPE
        ]
    if not versions:
        raise RuntimeError(f"{path}: no photon bunch block (type {PHOTON_BLOCK_TYPE})")
    return versions


# ---- reading IACT files --------------------------------------------------------------------


def read_shower_file(path) -> dict:
    """Bunches and header of one IACT file, in metres, ns and degrees.

    Returns a dict with per-bunch arrays (zem_m, photons, cx, cy, x_rel_m, y_rel_m, time_ns,
    x_ground_m, y_ground_m, core_distance_m, event), per-event arrays (zenith_deg, azimuth_deg,
    core_x_m, core_y_m) and the observation level. Telescope ids are 1-based indices into the
    telescope positions of the run (TELESCOPE rows, see eventio TelescopeDefinition).
    """
    columns = {
        key: []
        for key in (
            "zem_m",
            "photons",
            "cx",
            "cy",
            "x_rel_m",
            "y_rel_m",
            "time_ns",
            "x_ground_m",
            "y_ground_m",
            "core_x_m",
            "core_y_m",
            "event",
        )
    }
    events = {"zenith_deg": [], "azimuth_deg": [], "core_x_m": [], "core_y_m": []}
    with IACTFile(str(path)) as stream:
        run = stream.header
        site_m = (
            float(np.asarray(run["observation_height"], dtype=np.float64)[0]) * 0.01
        )
        positions = stream.telescope_positions
        for number, event in enumerate(stream):
            events["zenith_deg"].append(math.degrees(float(event.header["zenith"])))
            events["azimuth_deg"].append(math.degrees(float(event.header["azimuth"])))
            core_x, core_y = float(event.impact_x) * 0.01, float(event.impact_y) * 0.01
            events["core_x_m"].append(core_x)
            events["core_y_m"].append(core_y)
            for telescope, bunches in event.photon_bunches.items():
                if len(bunches) == 0:
                    continue
                index = int(telescope) - 1
                if not 0 <= index < len(positions):
                    raise ValueError(
                        f"telescope id {telescope} has no TELESCOPE position"
                    )
                x_tel = float(positions["x"][index]) * 0.01
                y_tel = float(positions["y"][index]) * 0.01
                z0 = float(positions["z"][index]) * 0.01
                x_rel = np.asarray(bunches["x"], dtype=np.float64) * 0.01
                y_rel = np.asarray(bunches["y"], dtype=np.float64) * 0.01
                cx = np.asarray(bunches["cx"], dtype=np.float64)
                cy = np.asarray(bunches["cy"], dtype=np.float64)
                x_g, y_g = ground_impact_m(x_rel, y_rel, cx, cy, z0, x_tel, y_tel)
                columns["zem_m"].append(
                    np.asarray(bunches["zem"], dtype=np.float64) * 0.01
                )
                columns["photons"].append(
                    np.asarray(bunches["photons"], dtype=np.float64)
                )
                columns["cx"].append(cx)
                columns["cy"].append(cy)
                columns["x_rel_m"].append(x_rel)
                columns["y_rel_m"].append(y_rel)
                columns["time_ns"].append(np.asarray(bunches["time"], dtype=np.float64))
                columns["x_ground_m"].append(x_g)
                columns["y_ground_m"].append(y_g)
                columns["core_x_m"].append(np.full(len(x_rel), core_x))
                columns["core_y_m"].append(np.full(len(x_rel), core_y))
                columns["event"].append(np.full(len(x_rel), number))
    result = {key: _concatenate(values) for key, values in columns.items()}
    result.update(
        {key: np.asarray(values, dtype=np.float64) for key, values in events.items()}
    )
    result["observation_level_m"] = site_m
    return result


def _concatenate(chunks):
    return np.concatenate(chunks) if chunks else np.empty(0, dtype=np.float64)


# ---- per-file analysis -------------------------------------------------------------------


def analyse_file(
    name: str,
    shower: dict,
    profile,
    log_text: str,
    flags: dict,
    long_format_violations: int,
) -> dict:
    """Checks (a) to (d) for one file, from the output of read_shower_file."""
    zeniths = shower["zenith_deg"]
    if zeniths.size == 0 or shower["zem_m"].size == 0:
        raise RuntimeError(f"{name}: no events or no photon bunches")
    theta0 = float(np.median(zeniths))
    zenith_spread = float(np.max(zeniths) - np.min(zeniths))
    zem_median = photon_weighted_median(shower["zem_m"], shower["photons"])
    curved = slant_depth_curved(
        profile, zem_median, theta0, site_altitude_m=shower["observation_level_m"]
    )
    plane = slant_depth_plane(profile, zem_median, theta0)

    near = within_radius(
        shower["x_ground_m"],
        shower["y_ground_m"],
        shower["core_x_m"],
        shower["core_y_m"],
        CORE_RADIUS_M,
    )
    if not np.any(near):
        raise RuntimeError(f"{name}: no bunches within {CORE_RADIUS_M} m of the core")
    mean_zenith = mean_zenith_deg(
        shower["cx"][near], shower["cy"][near], shower["photons"][near]
    )
    per_event = []
    for number in np.unique(shower["event"]):
        mask = shower["event"] == number
        per_event.append(
            {
                "event": int(number),
                "zem_median_m": photon_weighted_median(
                    shower["zem_m"][mask], shower["photons"][mask]
                ),
            }
        )
    return {
        "file": name,
        "events": int(zeniths.size),
        "bunches": int(shower["zem_m"].size),
        "photons": float(np.sum(shower["photons"])),
        "primary_zenith_deg": theta0,
        "zenith_spread_deg": zenith_spread,
        "observation_level_m": shower["observation_level_m"],
        "zem_median_m": zem_median,
        "slant_depth_curved_g_cm2": curved,
        "slant_depth_plane_g_cm2": plane,
        "curved_within_band": bool(XMAX_BAND_G_CM2[0] <= curved <= XMAX_BAND_G_CM2[1]),
        "near_core_bunches": int(np.count_nonzero(near)),
        "mean_zenith_near_core_deg": mean_zenith,
        "mean_zenith_error_deg": mean_zenith - theta0,
        "long_format_violations": long_format_violations,
        "curved_log": curved_showers.log_is_curved(log_text),
        "header": flags,
        "per_event_zem_median_m": [entry["zem_median_m"] for entry in per_event],
    }


def _sources() -> list[tuple[str, int, int, int]]:
    """(name, zenith, energy min, energy max) of the curved files, from curved_showers."""
    out = []
    for geometry, zenith, low, high in curved_showers.FIXTURE_CASES:
        if geometry == "curved":
            out.append(
                (
                    curved_showers.shower_name(zenith, geometry, low, high),
                    zenith,
                    low,
                    high,
                )
            )
    return out


def check_inputs(root: Path) -> list[str]:
    """Names of the missing shower files or logs in root; empty when all inputs exist."""
    missing = []
    for name, *_ in _sources():
        iact = root / name
        log = iact.with_suffix(".corsika.log")
        for path in (iact, log):
            if not path.is_file():
                missing.append(str(path))
    return missing


def profile_path() -> Path:
    configured = os.environ.get("ATMOSPHERE_PROFILE")
    if configured:
        return Path(configured)
    return common.REFERENCE / "cfg/common/atmprof1.dat"


def _plot_slant(rows, path) -> None:
    from obdeect_validation.adapters import plotting

    fig, axis = plotting.figure()
    low, high = XMAX_BAND_G_CM2
    axis.axhspan(low, high, color=plotting.LIMIT, alpha=0.2, label="X_max band")
    rows = sorted(rows, key=lambda row: row["primary_zenith_deg"])
    zeniths = [row["primary_zenith_deg"] for row in rows]
    axis.plot(
        zeniths,
        [row["slant_depth_curved_g_cm2"] for row in rows],
        "o-",
        color=plotting.NATIVE,
        label="curved (spherical)",
    )
    axis.plot(
        zeniths,
        [row["slant_depth_plane_g_cm2"] for row in rows],
        "s--",
        color=plotting.ANALYTIC,
        label="plane-parallel",
    )
    axis.set(
        xlabel="zenith angle theta0 [deg]",
        ylabel="slant depth of median zem [g/cm2]",
        yscale="log",
    )
    axis.legend()
    plotting.save(fig, path)


def _plot_zem(rows, path) -> None:
    from obdeect_validation.adapters import plotting

    fig, axis = plotting.figure()
    samples = [
        (
            f"{row['primary_zenith_deg']:g} deg",
            np.asarray(row["per_event_zem_median_m"]) / 1000.0,
            plotting.CYCLE[i % len(plotting.CYCLE)],
        )
        for i, row in enumerate(sorted(rows, key=lambda r: r["primary_zenith_deg"]))
    ]
    plotting.histograms(
        axis,
        samples,
        "per-event photon-weighted median zem [km a.s.l.]",
        bins=30,
        clip=None,
    )
    plotting.save(fig, path)


def main() -> int:
    root = Path(os.environ.get("VALIDATION_INPUTS", common.VALIDATION / "inputs"))
    missing = check_inputs(root)
    if missing:
        print(
            "missing input(s): " + ", ".join(missing) + "\n"
            "generate them with python -m obdeect_validation.adapters.curved_showers run "
            "--geometry curved --zenith ZENITH --energy-min LO --energy-max HI --root WORKDIR "
            "and copy WORKDIR/<name>.iact and WORKDIR/<name>.corsika.log to inputs/",
            file=sys.stderr,
        )
        return 3
    profile_file = profile_path()
    if not profile_file.is_file():
        print(f"missing input: atmosphere profile {profile_file}", file=sys.stderr)
        return 3
    profile = read_profile(profile_file)

    rows = []
    for name, _zenith, _low, _high in _sources():
        iact = root / name
        log_text = iact.with_suffix(".corsika.log").read_text(errors="replace")
        flags = curved_showers.read_header_flags(iact)
        shower = read_shower_file(iact)
        long_format_violations = compact_block_count(photon_block_versions(iact))
        rows.append(
            analyse_file(name, shower, profile, log_text, flags, long_format_violations)
        )

    outside = [row["file"] for row in rows if not row["curved_within_band"]]
    metrics = {
        "file_count": len(rows),
        "median_slant_depth_outside_range_count": len(outside),
        "mean_zenith_error_max_deg": max(
            abs(row["mean_zenith_error_deg"]) for row in rows
        ),
        "long_format_violation_count": sum(
            row["long_format_violations"] for row in rows
        ),
        "header_flag_mismatch_count": sum(
            1
            for row in rows
            if not (
                row["header"]["corsika_curved_flag"]
                and row["header"]["iact_curved_bit"]
            )
        ),
        "curved_log_missing_count": sum(1 for row in rows if not row["curved_log"]),
    }
    metrics["mean_zenith_error_max_deg"] = float(metrics["mean_zenith_error_max_deg"])
    details = {
        "limits": LIMITS,
        "atmosphere_profile": str(profile_file),
        "earth_radius_m": EARTH_RADIUS_M,
        "files": rows,
        "outside_band": outside,
    }
    _plot_slant(rows, Path("slant-depth-vs-zenith.png"))
    _plot_zem(rows, Path("zem-histograms.png"))
    common.write_results(metrics, details)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
