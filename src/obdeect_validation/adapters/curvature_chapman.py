"""Spherical slant depth of an exponential atmosphere, against the Chapman integral.

Atmosphere: density rho(h) = (A / H) exp(-(h - H2) / H) with the vertical optical depth above
the site A = 1, scale height H = 8 km and site (observation level) H2 = 2 km. The synthetic
transmission table samples the exact cumulative depth tau(h) = A (1 - exp(-(h - H2) / H)) at
99 heights (log-spaced, up to 100 km) for wavelengths 300 to 600 nm.

Photons: G-frame photons (one per zenith angle and emission height), propagated with
curvature="corsika" (curvature_applied = 1) and earth radius 6371315 m. The acceptance p must
equal exp(-tau_slant), with tau_slant the slant depth of the straight line from the telescope
(h_t = H2) to the emission altitude z.

Direction: the stream stores cx, cy (float32) and no cz column, so the zenith cosine is the
stream definition cos(theta_t) = sqrt(1 - cx^2 - cy^2), evaluated in float64 from the stored
values. Near the horizon this is sensitive to the float32 rounding of cx (relative error of
about 4e-4 in cos at 89.5 deg); the references use the same stored values, so the comparison
tests the kernel and not the input precision. A separately rounded cz column would give a second,
inconsistent zenith angle; this is why none is written.

References (no obdeect, no table interpolation in the analytic one):

- Analytic: the column of an exponential atmosphere along a straight line is the Chapman
  function. With x = (R + h) / H and c the cosine of the zenith angle of the ray (measured
  from the local vertical, away from the Earth),

      Ch(x, c) = x * integral_0^inf exp(-x (sqrt(1 + 2 y c + y^2) - 1)) dy,

  and the slant depth between the telescope and the emission point on the same ray is

      tau = A [Ch(x_t, c_t) - exp(-(z - h_t) / H) Ch(x_e, c_e)],

  with c_t = cos(theta_t) the zenith cosine at the telescope, p = (R + h_t) sin(theta_t) the
  ray invariant and c_e = sqrt(1 - (p / (R + z))^2). There is no elementary
  closed form for general zenith angle; the integral is evaluated with composite Simpson.
  Closed forms checked in the unit tests: Ch(x, 1) = 1 and Ch(x, 0) = x e^x K_1(x).

- Table form: the same integral as obdeect defines it for the table, with tau piecewise
  linear in log10 h between the table heights (composite Simpson, 10^6 intervals, in the
  substituted variable h = h_t + u^2 that removes the grazing singularity). It isolates the
  curvature kernel from the table interpolation.

- Plane-parallel: (tau(z) - tau(h_t)) / cos(theta_t) on the table, the flat limit.

The table-only error (table form against the exact Chapman depth) is measured and written
to details.json. The Chapman limit follows from it (see LIMITS).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common

EARTH_RADIUS_M = 6371315.0
SCALE_HEIGHT_M = 8000.0
SITE_M = 2000.0
VERTICAL_DEPTH = 1.0
TOP_M = 100000.0
N_HEIGHTS = 99
WAVELENGTHS_NM = tuple(range(300, 601))
ZENITH_DEG = (0, 30, 45, 60, 70, 75, 80, 85, 88, 89.5)
EMISSION_KM = (5, 15, 40)
FLAT_RADIUS_M = 1e12
PLANE_ZENITH_DEG = 60.0  # the monotonic check holds from this zenith angle
FLAT_MAX_ZENITH_DEG = 60.0  # see LIMITS["flat_limit_max_relative_error"]
REFERENCE_INTERVALS = 10**6
CHAPMAN_INTERVALS = 2 * 10**6
CHAPMAN_UPPER_Y = 2.0

LIMITS = {
    "max_transmission_relative_error_vs_reference": (
        1e-5,
        "the table-form reference is converged to about 5e-6 in tau (Simpson with 10^6 "
        "intervals agrees with 8-point Gauss-Legendre per table segment to 4e-6), which is the "
        "tolerance in transmission; the float32 acceptance column adds 6e-8 relative.",
    ),
    "max_transmission_relative_error_vs_chapman": (
        2.8e-3,
        "triangle inequality: |T_obdeect/T_chapman - 1| <= 1e-5 (vs reference) + the table-only "
        "error measured in this adapter (2.711e-3 at 89.5 deg, 5 km; the tabulated optical depth "
        "is log-linear in h), rounded up to 2.8e-3.",
    ),
    "flat_limit_max_relative_error": (
        1e-6,
        "R = 1e12 m against the plane-parallel result on the table, for zenith <= 60 deg. "
        "The physical 1/R residual of the spherical geometry is (H_emission/R)/cos^2(theta) "
        "relative in cos^2: about 1.5e-7 at 60 deg, 1e-3 at 89.5 deg, so the flat limit is "
        "tested only where it is reached (the residual at larger zenith is in details.json).",
    ),
    "curved_below_plane_violations": (
        0,
        "curvature shortens the slant path, so the curved transmission is not below the "
        "plane-parallel one from 60 deg on. A violation counts only below plane * (1 - 1e-6), "
        "the float32 noise of the acceptance column.",
    ),
    "sample_count": (30, "10 zenith angles x 3 emission heights; exact count."),
}


def table_heights_km() -> np.ndarray:
    """H1 heights [km], log-spaced from just above H2 to the top, rounded as written."""
    k = np.arange(1, N_HEIGHTS + 1)
    return np.round((SITE_M / 1000.0) * (TOP_M / SITE_M) ** (k / N_HEIGHTS), 6)


def exact_depth(height_m):
    """Exact cumulative optical depth from the site up to height_m (exponential atmosphere)."""
    return VERTICAL_DEPTH * (
        1.0 - np.exp(-(np.asarray(height_m, dtype=float) - SITE_M) / SCALE_HEIGHT_M)
    )


def write_table(path: Path) -> None:
    """Extinction table in the documented format: '# H2= ..., H1= ...' and one row per nm."""
    heights_km = table_heights_km()
    heights_m = heights_km * 1000.0
    tau = exact_depth(heights_m)
    lines = [
        "# Synthetic exponential atmosphere (validation case propagation.curvature-chapman)",
        f"# tau(h) = {VERTICAL_DEPTH:g} (1 - exp(-(h - H2) / {SCALE_HEIGHT_M / 1000:g} km)); "
        "the same depth at every wavelength",
        f"# H2= {SITE_M / 1000:.3f}, H1= " + " ".join(f"{h:.6f}" for h in heights_km),
    ]
    row = " ".join(f"{t:.12f}" for t in tau)
    for wavelength in WAVELENGTHS_NM:
        lines.append(f"{wavelength} {row}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_table(path: Path):
    """(heights [m], tau row at 400 nm [1]) read back from the file the adapter wrote."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    header = next(line for line in lines if line.startswith("# H2="))
    h2_km = float(header.split("H2=", 1)[1].split("H1=", 1)[0].strip().rstrip(","))
    h1_km = [float(v) for v in header.split("H1=", 1)[1].split()]
    row = next(line.split() for line in lines if line.split()[:1] == ["400"])
    tau = [float(v) for v in row[1:]]
    if len(tau) != len(h1_km):
        raise ValueError("table row and height list differ")
    heights = np.array([h2_km * 1000.0, *(h * 1000.0 for h in h1_km)])
    return heights, np.array([0.0, *tau])


class PiecewiseTable:
    """tau linear in log10 h between the table heights (the table's own form), zero slope above."""

    def __init__(self, heights_m, tau):
        self.log_h = np.log10(np.asarray(heights_m, dtype=float))
        self.tau = np.asarray(tau, dtype=float)
        self.slope = np.diff(self.tau) / np.diff(self.log_h)
        self.top_m = float(heights_m[-1])

    def depth(self, height_m):
        return np.interp(np.log10(height_m), self.log_h, self.tau)

    def derivative(self, height_m):
        h = np.asarray(height_m, dtype=float)
        k = np.clip(
            np.searchsorted(self.log_h, np.log10(h), side="right") - 1, 0, len(self.slope) - 1
        )
        value = self.slope[k] / (h * math.log(10.0))
        return np.where(h >= self.top_m, 0.0, value)


def cos_theta_squared(h, h_t, cos_t):
    """cos^2 of the angle to the local vertical at h on the ray with p = (R + h_t) sin(theta_t).

    cos^2 = 1 - sin^2(theta_t) ((R + h_t) / (R + h))^2 = cos^2(theta_t) + sin^2(theta_t)
    (1 - rho^2), with 1 - rho = (h - h_t) / (R + h) evaluated without cancellation.
    """
    sin2_t = 1.0 - cos_t * cos_t
    one_minus_rho = (h - h_t) / (EARTH_RADIUS_M + h)
    one_plus_rho = 1.0 + (EARTH_RADIUS_M + h_t) / (EARTH_RADIUS_M + h)
    return cos_t * cos_t + sin2_t * one_minus_rho * one_plus_rho


def simpson(f, a, b, intervals):
    """Composite Simpson rule with an even number of intervals."""
    if intervals % 2:
        raise ValueError("Simpson needs an even number of intervals")
    x = np.linspace(a, b, intervals + 1)
    weights = np.ones(intervals + 1)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    return (b - a) / (3.0 * intervals) * float(np.sum(weights * f(x)))


def slant_depth_table(table: PiecewiseTable, h_t, z, cos_t, intervals=REFERENCE_INTERVALS):
    """int_{h_t}^{z} d tau/dh / cos(theta(h)) dh with h = h_t + u^2 (the table form)."""
    upper = math.sqrt(z - h_t)

    def integrand(u):
        h = h_t + u * u
        return table.derivative(h) * 2.0 * u / np.sqrt(cos_theta_squared(h, h_t, cos_t))

    return simpson(integrand, 0.0, upper, intervals)


def chapman(x, c, intervals=CHAPMAN_INTERVALS, upper_y=CHAPMAN_UPPER_Y):
    """Chapman function Ch(x, c) = x int_0^inf exp(-x (sqrt(1 + 2 y c + y^2) - 1)) dy."""

    def integrand(y):
        q = 2.0 * y * c + y * y
        return np.exp(-x * q / (np.sqrt(1.0 + q) + 1.0))

    return x * simpson(integrand, 0.0, upper_y, intervals)


def slant_depth_chapman(h_t, z, cos_t):
    """Exact slant depth of the exponential atmosphere between h_t and z on one straight line."""
    sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
    p = (EARTH_RADIUS_M + h_t) * sin_t
    r_e = EARTH_RADIUS_M + z
    cos_e = math.sqrt(1.0 - (p / r_e) ** 2)
    x_t = (EARTH_RADIUS_M + h_t) / SCALE_HEIGHT_M
    x_e = r_e / SCALE_HEIGHT_M
    return VERTICAL_DEPTH * (
        chapman(x_t, cos_t) - math.exp(-(z - h_t) / SCALE_HEIGHT_M) * chapman(x_e, cos_e)
    )


def plane_parallel_depth(table: PiecewiseTable, h_t, z, cos_t):
    """Flat-limit slant depth on the table: (tau(z) - tau(h_t)) / cos(theta_t)."""
    return (table.depth(z) - table.depth(h_t)) / cos_t


def sample_geometry():
    """(zenith deg, emission height m a.s.l.) of the 30 samples, in order."""
    return [(zenith, km * 1000.0) for zenith in ZENITH_DEG for km in EMISSION_KM]


def stored_direction(zenith_deg: float):
    """(cx, cy) float32 of a photon arriving at zenith_deg, from the exact sine."""
    return np.float32(math.sin(math.radians(zenith_deg))), np.float32(0.0)


def zenith_cosine(cx, cy):
    """cos(theta_t) = sqrt(1 - cx^2 - cy^2) in float64 from the stored float32 direction."""
    cx = np.asarray(cx, dtype=np.float64)
    cy = np.asarray(cy, dtype=np.float64)
    return np.sqrt(np.clip(1.0 - cx * cx - cy * cy, 0.0, None))


def _write_input(path: Path, curvature_applied: int, cx, cy, z_emission) -> None:
    import obdeect.stream as stream

    writer = stream.Writer(str(path), "obdeect-validation")
    writer.write(
        "run",
        {
            "run_number": 1,
            "observation_level_m": SITE_M,
            "telescope_x_m": [0.0],
            "telescope_y_m": [0.0],
            "telescope_z_m": [0.0],
            "telescope_r_m": [12.0],
            "refraction_applied": 0,
            "wavelength_min_nm": float(WAVELENGTHS_NM[0]),
            "wavelength_max_nm": float(WAVELENGTHS_NM[-1]),
            "curvature_applied": curvature_applied,
        },
    )
    writer.write("event", {"event_number": 1})
    writer.write("telescope", {"telescope_id": 1, "frame": 1})
    n = len(cx)
    arrays = {
        "x": np.zeros(n, dtype=np.float32),
        "y": np.zeros(n, dtype=np.float32),
        "cx": np.asarray(cx, dtype=np.float32),
        "cy": np.asarray(cy, dtype=np.float32),
        "t": np.zeros(n, dtype=np.float32),
        "zem": np.asarray(z_emission, dtype=np.float32),
        "w": np.ones(n, dtype=np.float32),
    }
    writer.write_photons(1, 1, 0, arrays)
    writer.write("telescope_end", {"telescope_id": 1, "photons": float(n)})
    writer.write("event_end", {})
    writer.write("run_end", {})
    writer.close()


def _read_output(path: Path):
    """Acceptance p and weight w of the single telescope block, as float64."""
    from obdeect.stream import Reader

    p, w, zem = [], [], []
    for record in Reader(str(path)):
        if record["type"] == "photons":
            columns = record["columns"]
            p.append(np.asarray(columns["p"], dtype=np.float64))
            w.append(np.asarray(columns["w"], dtype=np.float64))
            zem.append(np.asarray(columns["zem"], dtype=np.float64))
    return np.concatenate(p), np.concatenate(w), np.concatenate(zem)


def main() -> None:
    import obdeect.propagation as propagation

    table_path = Path("curvature-table.dat")
    write_table(table_path)
    heights, tau_row = read_table(table_path)
    table = PiecewiseTable(heights, tau_row)

    samples = sample_geometry()
    directions = [stored_direction(zenith) for zenith, _ in samples]
    cx = np.array([d[0] for d in directions], dtype=np.float32)
    cy = np.array([d[1] for d in directions], dtype=np.float32)
    cos_stored = zenith_cosine(cx, cy)
    z_emission = np.array([z for _, z in samples], dtype=np.float32)

    _write_input(Path("in-curved.pst"), 1, cx, cy, z_emission)
    kwargs = dict(
        transmission=table_path,
        site_altitude_m=SITE_M,
        seed=20261115,
        mode="weighted",
        store_p=True,
        curvature="corsika",
    )
    curved = Path("curved.pst")
    propagation.propagate("in-curved.pst", curved, earth_radius_m=EARTH_RADIUS_M, **kwargs)
    p_curved, w_curved, zem_out = _read_output(curved)
    if p_curved.size != len(samples):
        raise RuntimeError(f"expected {len(samples)} output photons, got {p_curved.size}")
    if not np.array_equal(zem_out.astype(np.float32), z_emission):
        raise RuntimeError("output photons are not in input order")

    # Flat limit: the zenith <= 60 deg photons only. Weighted output drops w = 0 rows, and at
    # 89.5 deg the plane-parallel depth of 40 km is above 100 (T = 0 by the kernel's cutoff),
    # so the flat run must not contain it; a missing row would shift the comparison.
    flat_index = [i for i, (zenith, _) in enumerate(samples) if zenith <= FLAT_MAX_ZENITH_DEG]
    _write_input(Path("in-flat.pst"), 1, cx[flat_index], cy[flat_index], z_emission[flat_index])
    flat = Path("flat.pst")
    propagation.propagate("in-flat.pst", flat, earth_radius_m=FLAT_RADIUS_M, **kwargs)
    p_flat_subset, _, zem_flat = _read_output(flat)
    if p_flat_subset.size != len(flat_index) or not np.array_equal(
        zem_flat.astype(np.float32), z_emission[flat_index]
    ):
        raise RuntimeError("flat-limit run did not keep every photon in order")
    p_flat = np.full(len(samples), np.nan)
    p_flat[flat_index] = p_flat_subset

    rows = []
    for i, (zenith, z) in enumerate(samples):
        cos_t = float(cos_stored[i])  # the stream definition from the stored cx, cy
        tau_ref = slant_depth_table(table, SITE_M, z, cos_t)
        tau_chap = slant_depth_chapman(SITE_M, z, cos_t)
        tau_plane = plane_parallel_depth(table, SITE_M, z, cos_t)
        rows.append({
            "zenith_deg": zenith,
            "emission_m": z,
            "cos_theta_t": cos_t,
            "tau_reference_table": tau_ref,
            "tau_chapman": tau_chap,
            "tau_plane_table": tau_plane,
            "transmission_curved": float(p_curved[i]),
            "transmission_flat": None if np.isnan(p_flat[i]) else float(p_flat[i]),
            "transmission_reference_table": math.exp(-tau_ref),
            "transmission_chapman": math.exp(-tau_chap),
            "transmission_plane": math.exp(-tau_plane),
            "weight": float(w_curved[i]),
        })

    rel_ref = [
        abs(r["transmission_curved"] / r["transmission_reference_table"] - 1.0) for r in rows
    ]
    rel_chap = [abs(r["transmission_curved"] / r["transmission_chapman"] - 1.0) for r in rows]
    flat_rows = [r for r in rows if r["transmission_flat"] is not None]
    rel_flat = [abs(r["transmission_flat"] / r["transmission_plane"] - 1.0) for r in flat_rows]
    violations = sum(
        1
        for r in rows
        if r["zenith_deg"] >= PLANE_ZENITH_DEG
        and r["transmission_curved"] < r["transmission_plane"] * (1.0 - 1e-6)
    )

    hs = np.linspace(SITE_M, 40000.0, 400001)
    table_error = np.abs(table.depth(hs) - exact_depth(hs))
    table_transmission = max(
        abs(math.exp(-r["tau_reference_table"]) / math.exp(-r["tau_chapman"]) - 1.0) for r in rows
    )

    metrics = {
        "sample_count": len(rows),
        "max_transmission_relative_error_vs_reference": max(rel_ref),
        "max_transmission_relative_error_vs_chapman": max(rel_chap),
        "flat_limit_max_relative_error": max(rel_flat),
        "curved_below_plane_violations": violations,
    }
    details = {
        "limits": LIMITS,
        "table": {
            "file": str(table_path),
            "heights_km": heights[1:] / 1000.0,
            "wavelengths_nm": [WAVELENGTHS_NM[0], WAVELENGTHS_NM[-1]],
            "tau_exact_model": "A (1 - exp(-(h - H2)/H)), A = 1, H = 8 km, H2 = 2 km",
            "interpolation_max_abs_error_tau_2_40km": float(np.max(table_error)),
            "interpolation_height_of_max_error_m": float(hs[int(np.argmax(table_error))]),
            "max_transmission_error_table_vs_chapman": table_transmission,
        },
        "flat_residual_estimate": {
            "formula": "(z - h_t) / R / cos^2(theta_t): physical 1/R residual of the spherical "
            "geometry relative to the plane-parallel depth, order of magnitude only",
            "zenith_deg": list(ZENITH_DEG),
            "value_at_40km": [
                (40000.0 - SITE_M) / EARTH_RADIUS_M / math.cos(math.radians(z)) ** 2
                for z in ZENITH_DEG
            ],
        },
        "samples": rows,
        "reference": {
            "earth_radius_m": EARTH_RADIUS_M,
            "table_intervals": REFERENCE_INTERVALS,
            "chapman_intervals": CHAPMAN_INTERVALS,
            "chapman_upper_y": CHAPMAN_UPPER_Y,
        },
    }
    common.write_results(metrics, details)
    _plot(rows, Path("curvature-chapman.png"))


def _plot(rows, path: Path) -> None:
    """Transmission against zenith (curved, plane-parallel, Chapman) and the curved error."""
    from obdeect_validation.adapters import plotting

    fig, axes = plotting.figure(1, 2, (10, 3.8))
    heights = sorted({r["emission_m"] for r in rows})
    for i, height in enumerate(heights):
        color = plotting.CYCLE[i % len(plotting.CYCLE)]
        subset = sorted(
            (r for r in rows if r["emission_m"] == height), key=lambda r: r["zenith_deg"]
        )
        zenith = [r["zenith_deg"] for r in subset]
        label = f"{height / 1000:g} km"
        axes[0, 0].plot(
            zenith,
            [r["transmission_chapman"] for r in subset],
            "--",
            color=plotting.ANALYTIC,
            linewidth=1.0,
        )
        axes[0, 0].plot(
            zenith, [r["transmission_plane"] for r in subset], ":", color=plotting.LIMIT
        )
        axes[0, 0].plot(
            zenith,
            [r["transmission_curved"] for r in subset],
            "o",
            color=color,
            label=label,
            markerfacecolor="none",
        )
        axes[0, 1].plot(
            zenith,
            [r["transmission_curved"] / r["transmission_chapman"] - 1.0 for r in subset],
            "o-",
            color=color,
            label=label,
        )
    axes[0, 0].set(xlabel="zenith angle [deg]", ylabel="transmission", yscale="log")
    axes[0, 0].legend(title="emission (dashed: Chapman, dotted: plane-parallel)")
    axes[0, 1].set(xlabel="zenith angle [deg]", ylabel="obdeect / Chapman - 1")
    axes[0, 1].legend(title="emission")
    plotting.save(fig, path)


if __name__ == "__main__":
    main()
