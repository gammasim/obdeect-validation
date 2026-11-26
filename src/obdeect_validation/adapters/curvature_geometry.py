"""Spherical geometry of sky-frame photons (geometry only, no extinction).

Input: sky-frame photons (x, y, z above the observation level, direction cosines, t = 0, w = 1),
written once per telescope block. Run keys: curvature_applied = 0 (required by curvature
"spherical"), observation level H = 2000 m a.s.l.

Runs (obdeect.propagation.propagate with extinction=False):
- spherical, earth_radius_m = R = 6371315 m (the case under test),
- none (plane-parallel, straight lines, the flat reference),
- spherical, earth_radius_m = 1e12 m (flat limit, must reproduce none).

Telescopes: lateral distance d in (0, 2, 5, 10) km along x, altitude z in (0, 3.7, 9.9) m: 12
telescopes. The telescope plane of the stream is z = telescope z (sky), G positions are relative
to the telescope centre.

Analytic expectations (independent of obdeect):

(a) Emission altitude a.s.l. of a sky point (x, y, z):
        zem = sqrt(x^2 + y^2 + (z + R + H)^2) - R
    (the sphere of radius R with its centre at z = -(R + H), the observation level at z = 0).
    Flat: zem = z + H. Shift of zem for a lateral distance d: sqrt(d^2 + L^2) - L with
    L = R + z + H, which is d^2 / (2R) to a relative (z + H) / R, about 1e-3.

(b) Vacuum arrival: the straight line from the emission point to the plane z = z_t. With
    cz_mag = sqrt(1 - cx^2 - cy^2) (the stream definition of cz), the path is
    s = (z - z_t) / cz_mag, the arrival is (x + cx s, y + cy s), the G position subtracts the
    telescope (x_t, y_t), and the time is t + s / c.

(c) Photons with cz > 0 (upwards) emitted above a telescope plane (all sky heights here) never
    reach it (moving_away): each one is dropped once per telescope, so the dropped count is
    n_up x n_telescopes.

(d) The flat limit R -> 1e12 m reproduces the plane-parallel run (curvature none) to float32
    precision, because the residual of the spherical geometry is d^2/(2R) and the telescope
    tilt (at most lateral / R) is not applied to the arrival plane.

(e) The lateral shift at d = 5 km is d^2 / (2R) within 1 %.

Output precision: the stream columns are float32. Positions up to 2^15 m, times up to 2^15 ns
and zem up to 2^15 m are stored with at most half an ulp of rounding (about 1e-3 m, ns or m).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common, endtoend_common, streams

EARTH_RADIUS_M = 6371315.0
FLAT_RADIUS_M = 1e12
SITE_M = 2000.0
LIGHT_SPEED_M_PER_S = 299792458.0
LATERAL_KM = (0.0, 2.0, 5.0, 10.0)
TELESCOPE_Z_M = (0.0, 3.7, 9.9)
EMISSION_LATERAL_M = (
    (0.0, 0.0),
    (2000.0, 0.0),
    (5000.0, 0.0),
    (-5000.0, 0.0),
    (0.0, 5000.0),
    (7000.0, -7000.0),
)
SKY_HEIGHT_M = (
    1000.0,
    3000.0,
)  # above the observation level: emission at 3 and 5 km a.s.l.
# (zenith deg, azimuth deg) of downward photons; one upward direction (30 deg, 0) is added.
DOWNWARD = (
    (0.0, 0.0),
    (30.0, 0.0),
    (30.0, 90.0),
    (30.0, 180.0),
    (30.0, 270.0),
    (60.0, 0.0),
    (60.0, 90.0),
    (60.0, 180.0),
    (60.0, 270.0),
)
UPWARD = ((30.0, 0.0),)
LATERAL_SHIFT_KM = 5.0

LIMITS = {
    "max_zem_error_m": (
        5e-3,
        (
            "zem is stored as float32: half an ulp is 2.4e-4 m at 5 km a.s.l. (the highest emission "
            "altitude used). The limit is 20 times that and about 1/400 of the curvature shift it "
            "tests (1.96 m at 5 km lateral), so it resolves the sphere and not the storage."
        ),
    ),
    "max_position_error_m": (
        1e-3,
        (
            "G positions are float32 with magnitudes below 2.5e4 m: half an ulp is at most 9.8e-4 m "
            "there, so the limit is the float32 storage, not a physical tolerance."
        ),
    ),
    "max_time_error_ns": (
        2e-3,
        (
            "times are float32 with magnitudes below 2^15 ns (about 2.7e4 ns here): half an ulp is "
            "9.8e-4 ns; the limit is twice that."
        ),
    ),
    "flat_limit_max_difference": (
        1e-3,
        (
            "difference of the R = 1e12 m spherical run and the plane-parallel run in x, y and zem "
            "[m]: one float32 ulp at 8 to 16 km (9.8e-4 m). The physical 1/R residual is below "
            "1e-6 m for these geometries."
        ),
    ),
    "dropped_count_mismatch": (
        0,
        "exact count: n_upward x n_telescopes photons must be dropped, none others.",
    ),
    "lateral_shift_relative_error": (
        1e-2,
        (
            "relative error of the zem shift at 5 km lateral distance against d^2 / (2R). The "
            "closed form differs from d^2 / (2R) by (z + H) / R, about 1e-3 at these heights, so "
            "1 % is ten times that."
        ),
    ),
    "sample_count": (
        1296,
        "108 downward photons x 12 telescopes; every one reaches its plane (exact count).",
    ),
}


def telescope_table():
    """(telescope_id, x, y, z) of the 12 telescopes, lateral distance outer, altitude inner."""
    rows = []
    telescope_id = 1
    for lateral_km in LATERAL_KM:
        for z in TELESCOPE_Z_M:
            rows.append((telescope_id, lateral_km * 1000.0, 0.0, z))
            telescope_id += 1
    return rows


def photon_table():
    """Sky photons: dicts with x, y, z (above obs level), cx, cy, upward flag, emission index."""
    photons = []
    for x, y in EMISSION_LATERAL_M:
        for z in SKY_HEIGHT_M:
            for zenith, azimuth in DOWNWARD:
                photons.append(_direction(x, y, z, zenith, azimuth, upward=False))
            for zenith, azimuth in UPWARD:
                photons.append(_direction(x, y, z, zenith, azimuth, upward=True))
    return photons


def _direction(x, y, z, zenith, azimuth, upward):
    sin_t = math.sin(math.radians(zenith))
    cos_t = math.cos(math.radians(zenith))
    cx = np.float32(sin_t * math.cos(math.radians(azimuth)))
    cy = np.float32(sin_t * math.sin(math.radians(azimuth)))
    cz = np.float32(cos_t if upward else -cos_t)
    return {
        "x": np.float32(x),
        "y": np.float32(y),
        "z": np.float32(z),
        "cx": cx,
        "cy": cy,
        "cz": cz,
        "upward": upward,
        "zenith_deg": zenith,
    }


def cz_magnitude(cx, cy):
    """sqrt(1 - cx^2 - cy^2) in float64 from the stored float32 direction (stream definition)."""
    cx = np.asarray(cx, dtype=np.float64)
    cy = np.asarray(cy, dtype=np.float64)
    return np.sqrt(np.clip(1.0 - cx * cx - cy * cy, 0.0, None))


def emission_altitude_sphere(x, y, z, observation_m=SITE_M, radius_m=EARTH_RADIUS_M):
    """(a) zem a.s.l. = sqrt(x^2 + y^2 + (z + R + H)^2) - R, float64."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    return np.sqrt(x * x + y * y + (z + radius_m + observation_m) ** 2) - radius_m


def straight_arrival(x, y, z, cx, cy, z_plane):
    """(b) arrival (x, y) on the plane z = z_plane, the path s and the time in ns (no profile).

    Inputs are the stored float32 values; the path is (z - z_plane) / cz_mag and the distance
    along the unit direction.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    cx = np.asarray(cx, dtype=np.float64)
    cy = np.asarray(cy, dtype=np.float64)
    s = (z - z_plane) / cz_magnitude(cx, cy)
    return x + cx * s, y + cy * s, s / LIGHT_SPEED_M_PER_S * 1e9


def lateral_shift_closed_form(
    d_m, z_sky_m, observation_m=SITE_M, radius_m=EARTH_RADIUS_M
):
    """sqrt(d^2 + L^2) - L with L = R + z + H: the closed-form shift of zem at lateral d."""
    level = radius_m + z_sky_m + observation_m
    return math.sqrt(d_m * d_m + level * level) - level


def _write_input(path: Path, curvature_applied: int, photons, telescopes) -> None:
    n = len(photons)
    arrays = {
        "x": np.array([p["x"] for p in photons], dtype=np.float32),
        "y": np.array([p["y"] for p in photons], dtype=np.float32),
        "z": np.array([p["z"] for p in photons], dtype=np.float32),
        "cx": np.array([p["cx"] for p in photons], dtype=np.float32),
        "cy": np.array([p["cy"] for p in photons], dtype=np.float32),
        "cz": np.array([p["cz"] for p in photons], dtype=np.float32),
        "t": np.zeros(n, dtype=np.float32),
        "w": np.ones(n, dtype=np.float32),
    }
    run = {
        "run_number": 1,
        "observation_level_m": SITE_M,
        "telescope_x_m": [t[1] for t in telescopes],
        "telescope_y_m": [t[2] for t in telescopes],
        "telescope_z_m": [t[3] for t in telescopes],
        "telescope_r_m": [12.0] * len(telescopes),
        "refraction_applied": 0,
        "wavelength_min_nm": 300.0,
        "wavelength_max_nm": 600.0,
        "curvature_applied": curvature_applied,
    }
    blocks = [(telescope_id, arrays) for telescope_id, _, _, _ in telescopes]
    streams.write_stream(path, "obdeect-validation", run, blocks)


def _read_output(path: Path):
    """Per telescope id: arrays x, y, t, zem (float64) of the kept photons, in stream order."""
    rows = endtoend_common.photon_rows(path)
    blocks = {}
    for telescope in dict.fromkeys(rows["telescope"].tolist()):
        on_telescope = rows["telescope"] == telescope
        blocks[telescope] = {
            key: np.asarray(rows[key][on_telescope], dtype=np.float64)
            for key in ("x", "y", "t", "zem")
        }
    return blocks


def _run(propagation, source, output, radius, curvature, **kwargs):
    # The Earth radius belongs to a curved run only (the API rejects it with curvature "none").
    if curvature != "none":
        kwargs["earth_radius_m"] = radius
    return propagation.propagate(
        str(source),
        str(output),
        site_altitude_m=SITE_M,
        seed=20261115,
        extinction=False,
        curvature=curvature,
        **kwargs,
    )


def analyse(root: Path) -> tuple[dict, dict]:
    from obdeect import propagation

    telescopes = telescope_table()
    photons = photon_table()
    reaching = [p for p in photons if not p["upward"]]
    upward = [p for p in photons if p["upward"]]
    _write_input(root / "in-sky.pst", 0, photons, telescopes)

    summary = _run(
        propagation,
        root / "in-sky.pst",
        root / "spherical.pst",
        EARTH_RADIUS_M,
        "spherical",
    )
    flat_summary = _run(
        propagation, root / "in-sky.pst", root / "flat.pst", FLAT_RADIUS_M, "spherical"
    )
    plane_summary = _run(
        propagation, root / "in-sky.pst", root / "plane.pst", EARTH_RADIUS_M, "none"
    )
    curved = _read_output(root / "spherical.pst")
    flat = _read_output(root / "flat.pst")
    plane = _read_output(root / "plane.pst")

    # Column arrays of the reaching photons (input order), float64 from the stored float32.
    x_in = np.array([p["x"] for p in reaching], dtype=np.float64)
    y_in = np.array([p["y"] for p in reaching], dtype=np.float64)
    z_in = np.array([p["z"] for p in reaching], dtype=np.float64)
    cx_in = np.array([p["cx"] for p in reaching], dtype=np.float64)
    cy_in = np.array([p["cy"] for p in reaching], dtype=np.float64)
    zem_expected = emission_altitude_sphere(x_in, y_in, z_in)
    zem_flat = z_in + SITE_M

    position_errors, time_errors, zem_errors, flat_differences = [], [], [], []
    for telescope_id, x_t, y_t, z_t in telescopes:
        block = curved[telescope_id]
        if block["x"].size != len(reaching):
            raise RuntimeError(
                f"telescope {telescope_id}: {block['x'].size} photons kept, "
                f"expected {len(reaching)}"
            )
        x_arrival, y_arrival, t_arrival = straight_arrival(
            x_in, y_in, z_in, cx_in, cy_in, z_t
        )
        position_errors.append(block["x"] - (x_arrival - x_t))
        position_errors.append(block["y"] - (y_arrival - y_t))
        time_errors.append(block["t"] - t_arrival)
        zem_errors.append(block["zem"] - zem_expected)
        # Flat limit: R = 1e12 m against plane-parallel (curvature none) on the same photons.
        for key in ("x", "y", "zem"):
            flat_differences.append(
                np.abs(flat[telescope_id][key] - plane[telescope_id][key])
            )
    dropped_expected = len(upward) * len(telescopes)

    # Lateral shift of zem: emission points at radial distance 5 km (zem does not depend on the
    # telescope), the closed form and d^2 / (2R).
    radial = np.hypot(x_in, y_in)
    shift_index = np.flatnonzero(np.isclose(radial, LATERAL_SHIFT_KM * 1000.0))
    first = curved[telescopes[0][0]]
    shift_obdeect = first["zem"][shift_index] - (z_in[shift_index] + SITE_M)
    shift_expected = (LATERAL_SHIFT_KM * 1000.0) ** 2 / (2.0 * EARTH_RADIUS_M)
    shift_relative = np.abs(shift_obdeect / shift_expected - 1.0)
    shift_closed = np.array(
        [
            lateral_shift_closed_form(LATERAL_SHIFT_KM * 1000.0, z)
            for z in z_in[shift_index]
        ]
    )

    metrics = {
        "sample_count": len(reaching) * len(telescopes),
        "max_zem_error_m": float(max(np.max(np.abs(e)) for e in zem_errors)),
        "max_position_error_m": float(max(np.max(np.abs(e)) for e in position_errors)),
        "max_time_error_ns": float(max(np.max(np.abs(e)) for e in time_errors)),
        "flat_limit_max_difference": float(max(np.max(d) for d in flat_differences)),
        "dropped_count_mismatch": abs(int(summary["dropped"]) - dropped_expected),
        "lateral_shift_relative_error": float(np.max(shift_relative)),
    }
    details = {
        "limits": LIMITS,
        "telescopes": [
            {"telescope_id": t[0], "lateral_m": t[1], "z_m": t[3]} for t in telescopes
        ],
        "photons": {
            "sky_total": len(photons),
            "reaching": len(reaching),
            "upward": len(upward),
            "dropped_expected": dropped_expected,
        },
        "run_summaries": {
            "spherical": summary,
            "flat_1e12": flat_summary,
            "none": plane_summary,
        },
        "flat_limit_difference_m": "max over x, y, zem of |R=1e12 spherical - none| (metric)",
        "lateral_shift": [
            {
                "x_m": float(x_in[i]),
                "y_m": float(y_in[i]),
                "z_sky_m": float(z_in[i]),
                "shift_obdeect_m": float(shift_obdeect[k]),
                "d2_over_2R_m": shift_expected,
                "closed_form_m": float(shift_closed[k]),
            }
            for k, i in enumerate(shift_index)
        ],
        "time_reference": "t_out = s / c, s = (z - z_t) / sqrt(1 - cx^2 - cy^2)",
    }
    _plot(
        x_in,
        y_in,
        z_in,
        zem_expected,
        zem_flat,
        first["zem"],
        root / "curvature-geometry.png",
    )
    return metrics, details


def _plot(x_in, y_in, z_in, zem_expected, zem_flat, zem_obdeect, path: Path) -> None:
    """Emission altitude shift against lateral distance: obdeect, the sphere and the flat value."""
    from obdeect_validation.adapters import plotting

    fig, axis = plotting.figure()
    radial = np.hypot(x_in, y_in) / 1000.0
    axis.plot(
        radial,
        zem_expected - zem_flat,
        "o",
        color=plotting.NATIVE,
        markerfacecolor="none",
        label="sphere (closed form)",
    )
    axis.plot(
        radial,
        zem_obdeect - zem_flat,
        "x",
        color=plotting.REFERENCE,
        label="obdeect",
    )
    d = np.linspace(0.0, max(radial.max(), 1e-9), 101)
    axis.plot(
        d,
        (d * 1000.0) ** 2 / (2.0 * EARTH_RADIUS_M),
        "--",
        color=plotting.ANALYTIC,
        label="d^2 / (2R)",
    )
    axis.set(xlabel="lateral distance of emission [km]", ylabel="zem - flat zem [m]")
    axis.legend()
    plotting.save(fig, path)


if __name__ == "__main__":
    common.results_main(analyse, __doc__.splitlines()[0])
