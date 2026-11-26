"""Curved-ray refraction in spherical shells: independent numpy ray tracer against obdeect.

Reference (this module, numpy only, no obdeect code): Bouguer's law n r sin(theta) = K in
spherical shells about the Earth centre, n(r) = 1 + nu(h), nu from the atmprof table with a
natural cubic spline of ln(nu) in altitude (obdeect interpolates ln(n-1) the same way, on a
3 m grid with log-linear cells, which differs from the exact spline by about 1e-8 relative).

Candidate: obdeect `propagate` with curvature="spherical" (and "none" / R -> infinity for the
plane-parallel limit) on sky-frame photons written with obdeect.stream.Writer.

Geometry (sky frame, telescope at the origin, Earth centre C = (0, 0, -(R + h_site))): a photon
emitted at altitude h_e from P with direction d = (cx, 0, cz), cz < 0. Its vacuum straight line
passes the telescope level r_t = R + h_site at zenith angle theta_s = asin(b / r_t), with
b = |(P - C) x d| and K = n_e b. Photons are stored in float32, so the reference reads back the
stored values (P, cx) and uses them everywhere; cz is not written and is derived by obdeect, so
the stored direction is exactly unit.

Integrals over r from r_t to r_e (r = r_t + u^2, composite Simpson in u):
    shift   = r_t int (tan theta_ref - tan theta_vac) / r dr             [m]
    excess  = (1/c) int (n / cos theta_ref - 1 / cos theta_vac) dr      [ns]
    straight= (1/c) int nu / cos theta_vac dr                            [ns]
Time = s/c + excess, with s the straight path from P to the telescope plane.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common

C_M_PER_NS = 0.299792458
EARTH_RADIUS_M = 6371315.0
PLANE_RADIUS_M = 1.0e12
SITE_ALTITUDE_M = 2000.0
SEED = 1
HEIGHTS_KM = (3.0, 8.0, 20.0, 60.0)
ZENITHS_DEG = (0.0, 20.0, 40.0, 60.0, 75.0, 80.0, 85.0)
DISPERSION_REFERENCE_NM = 500.0
PRODUCTION_INTERVALS = 2**16
COARSE_INTERVALS = 2**15
WORKSPACE = Path(__file__).resolve().parents[4]


def ciddor_refractivity_e8(wavelength_nm: float) -> float:
    """Ciddor (1996) standard dry-air refractivity, in units of 1e-8."""
    s2 = (1000.0 / wavelength_nm) ** 2
    return 5792105.0 / (238.0185 - s2) + 167917.0 / (57.362 - s2)


def ciddor_group_factor(wavelength_nm: float) -> float:
    """Group refractivity relative to the profile's reference-wavelength phase refractivity."""
    s2 = (1000.0 / wavelength_nm) ** 2
    phase = ciddor_refractivity_e8(wavelength_nm)
    group = phase + 2.0 * s2 * (
        5792105.0 / (238.0185 - s2) ** 2 + 167917.0 / (57.362 - s2) ** 2
    )
    return group / phase


# Limits fixed before looking at results. Each entry: (rule, justification).
LIMITS = {
    "position": (
        "max(5 mm, 1e-3 |shift_ref|) for theta_s <= 80 deg; 0.05 m at 85 deg",
        (
            "obdeect docs/propagation/curvature.md 'Accuracy and limits': table and altitude "
            "interpolation about 1e-3 m at 80 deg and 0.02 to 0.03 m at 85 deg (emission 20 km); "
            "the limits keep a margin of about 2 at 85 deg and 5 at 80 deg."
        ),
    ),
    "time": (
        "0.01 ns + half float32 ulp of the output t + 5e-6 |excess_ref|",
        (
            "docs: 0.06 to 0.12 ns on a 3.6e4 ns delay at 85 deg (2 to 3e-6 relative); 5e-6 of "
            "the excess keeps that margin. The float32 time column rounds to half an ulp (up to "
            "0.13 ns near 2e6 ns), which is added explicitly. The 0.01 ns floor is the absolute "
            "floor of the docs' test tolerance."
        ),
    ),
    "direction": (
        "|cx_out / cx_in - (1+nu_e)/(1+nu_t)| <= 1e-6",
        (
            "Snell's law is exact at the telescope level; float32 storage has a relative "
            "precision of 6e-8, so 1e-6 leaves about 16 ulp."
        ),
    ),
    "straight_delay": (
        "time limit (above) on the air delay with refraction_applied = 1",
        "same accuracy argument as the time limit; the air delay is a plain integral of nu.",
    ),
    "plane_limit": (
        (
            "the position and time limits (above), normalised: R -> infinity spherical against "
            "curvature none and against this tracer with R -> infinity"
        ),
        "the spherical mode must reduce to the plane-parallel result in the flat limit.",
    ),
    "convergence": (
        "relative change <= 1e-8 when the Simpson step is halved",
        (
            "composite Simpson has error O(h^4); a halved step must change the reference by far "
            "less than the candidate tolerances, otherwise the reference itself is not converged."
        ),
    ),
}


def profile_path() -> Path:
    """Atmprof table: ATMOSPHERE_PROFILE, then the sim_telarray cfg file, then obdeect's test data."""
    configured = os.environ.get("ATMOSPHERE_PROFILE")
    if configured:
        return Path(configured)
    reference = Path(os.environ.get("SIMTEL_REFERENCE_SOURCE", common.REFERENCE))
    candidate = reference / "cfg/common/atmprof1.dat"
    if candidate.is_file():
        return candidate
    relative = Path(
        "obdeect/src/propagation/cpp/tests/data/atmprof_ecmwf_north_winter_fixed.dat"
    )
    for checkout in (WORKSPACE, Path(__file__).resolve().parents[3]):
        candidate = checkout / relative
        if candidate.is_file():
            return candidate
    return WORKSPACE / relative


def load_table(path) -> tuple[np.ndarray, np.ndarray]:
    """Altitudes [m a.s.l.] and n-1 from an atmprof table (column 1 km, column 4 n-1)."""
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        rows.append([float(value) for value in fields[:4]])
    table = np.asarray(rows, dtype=np.float64)
    if table.ndim != 2 or table.shape[0] < 3:
        raise ValueError(f"no atmprof levels in {path}")
    altitude = table[:, 0] * 1000.0
    nu = table[:, 3]
    if np.any(np.diff(altitude) <= 0) or np.any(nu <= 0):
        raise ValueError("atmprof altitudes must ascend and n-1 must be positive")
    return altitude, nu


class LogSpline:
    """Natural cubic spline of ln(nu) in altitude; ln(nu) is exact at the table levels."""

    def __init__(self, x, nu):
        self.x = np.asarray(x, dtype=np.float64)
        self.y = np.log(np.asarray(nu, dtype=np.float64))
        n = self.x.size
        h = np.diff(self.x)
        system = np.zeros((n, n))
        rhs = np.zeros(n)
        system[0, 0] = system[-1, -1] = 1.0
        for i in range(1, n - 1):
            system[i, i - 1] = h[i - 1] / 6.0
            system[i, i] = (h[i - 1] + h[i]) / 3.0
            system[i, i + 1] = h[i] / 6.0
            rhs[i] = (self.y[i + 1] - self.y[i]) / h[i] - (
                self.y[i] - self.y[i - 1]
            ) / h[i - 1]
        self.m = np.linalg.solve(system, rhs)
        self.h = h

    def __call__(self, altitude):
        """ln(nu) at the altitudes (linear extension of the bottom value below the table)."""
        q = np.clip(np.asarray(altitude, dtype=np.float64), self.x[0], self.x[-1])
        i = np.clip(np.searchsorted(self.x, q, side="right") - 1, 0, self.x.size - 2)
        hi = self.h[i]
        a = (self.x[i + 1] - q) / hi
        b = (q - self.x[i]) / hi
        return (
            a * self.y[i]
            + b * self.y[i + 1]
            + ((a**3 - a) * self.m[i] + (b**3 - b) * self.m[i + 1]) * hi * hi / 6.0
        )


class NuProfile:
    """nu(h) = n - 1 from the spline of ln(nu); zero above the table top, constant below it."""

    def __init__(self, path):
        self.altitude, nu = load_table(path)
        self.spline = LogSpline(self.altitude, nu)
        self.top_m = float(self.altitude[-1])

    def __call__(self, altitude_m):
        h = np.asarray(altitude_m, dtype=np.float64)
        value = np.exp(self.spline(h))
        return np.where(h > self.top_m, 0.0, value)


def simpson_weights(intervals: int) -> np.ndarray:
    if intervals < 2 or intervals % 2:
        raise ValueError("Simpson needs an even number of intervals >= 2")
    weights = np.ones(intervals + 1)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    return weights * (1.0 / (3.0 * intervals))


def ray_integrals(nu, R, h_t, h_e, b, intervals=PRODUCTION_INTERVALS) -> dict:
    """shift [m], excess and straight air delay [ns], and the straight path integral [m].

    nu: callable n-1 at altitudes [m]; R Earth radius [m]; h_t telescope level; h_e emission
    altitude; b impact parameter of the vacuum line about the Earth centre [m].
    """
    delta_e = h_e - h_t
    if delta_e <= 0.0:
        raise ValueError("emission must be above the telescope level")
    r_t = R + h_t
    a0 = (r_t - b) * (r_t + b)  # r_t^2 - b^2, without cancellation
    n_e = 1.0 + float(nu(np.array([h_e]))[0])
    kk = n_e * b
    upper = math.sqrt(delta_e)
    u = np.linspace(0.0, upper, intervals + 1)
    delta = u * u
    r = r_t + delta
    nu_r = nu(h_t + delta)
    n = 1.0 + nu_r
    p = n * r
    q = np.sqrt(a0 + delta * (2.0 * r_t + delta))  # r cos(theta_vac)
    sq = np.sqrt((p - kk) * (p + kk))  # p cos(theta_ref)
    jacobian = 2.0 * u  # dr = 2 u du
    weights = simpson_weights(intervals) * upper
    shift = r_t * (kk / sq - b / q) / r
    excess = n * p / sq - r / q
    straight = nu_r * r / q
    path = r / q
    return {
        "shift_m": float(np.sum(weights * shift * jacobian)),
        "excess_ns": float(np.sum(weights * excess * jacobian)) / C_M_PER_NS,
        "straight_ns": float(np.sum(weights * straight * jacobian)) / C_M_PER_NS,
        "path_m": float(np.sum(weights * path * jacobian)),
        "nu_e": float(nu_e_value(nu, h_e)),
    }


def nu_e_value(nu, h_e) -> float:
    return float(nu(np.array([h_e]))[0])


def emission_point(
    h_e_m: float, theta_deg: float, R: float, h_t: float = SITE_ALTITUDE_M
):
    """Point P at altitude h_e whose vacuum line passes the telescope; direction (cx, 0, cz).

    Exact closed form: |P - C|^2 = s^2 + 2 r_t s cos(theta) + r_t^2 = r_e^2 with P = O - s d.
    The root is written in the cancellation-free form s = X / (r_t c + sqrt(r_t^2 c^2 + X)).
    """
    theta = math.radians(theta_deg)
    c = math.cos(theta)
    r_t = R + h_t
    x_excess = (h_e_m - h_t) * (2.0 * R + h_e_m + h_t)  # r_e^2 - r_t^2
    s = x_excess / (r_t * c + math.sqrt(r_t * r_t * c * c + x_excess))
    point = np.array([-s * math.sin(theta), 0.0, s * c])
    direction = np.array([math.sin(theta), 0.0, -c])
    return point, direction


def stored_sample(point, direction) -> dict:
    """Float32 values as written to the stream (cz omitted: obdeect derives it from cx, cy)."""
    stored = {
        "x": np.float32(point[0]),
        "y": np.float32(point[1]),
        "z": np.float32(point[2]),
        "cx": np.float32(direction[0]),
        "cy": np.float32(direction[1]),
    }
    return {name: float(value) for name, value in stored.items()}


def reference_sample(nu, R, h_t, stored, intervals=PRODUCTION_INTERVALS) -> dict:
    """Reference for one stored sample: positions, times, direction and convergence inputs."""
    point = np.array([stored["x"], stored["y"], stored["z"]])
    cx, cy = stored["cx"], stored["cy"]
    cz = -math.sqrt(1.0 - cx * cx - cy * cy)
    d = np.array([cx, cy, cz])
    centre = np.array([0.0, 0.0, -(R + h_t)])
    pc = point - centre
    r_e = float(np.linalg.norm(pc))
    h_e = r_e - R
    b = float(np.linalg.norm(np.cross(pc, d)))
    sin_s = b / (R + h_t)
    theta_s_deg = math.degrees(math.asin(min(sin_s, 1.0)))
    s_plane = point[2] / -cz
    x_s = point[0] + s_plane * cx
    y_s = point[1] + s_plane * cy
    main = ray_integrals(nu, R, h_t, h_e, b, intervals)
    nu_t = float(nu(np.array([h_t]))[0])
    scale = (1.0 + main["nu_e"]) / (1.0 + nu_t)
    horizontal = math.hypot(cx, cy)
    hx = cx / horizontal if horizontal > 0.0 else 0.0
    hy = cy / horizontal if horizontal > 0.0 else 0.0
    return {
        "h_e_m": h_e,
        "theta_s_deg": theta_s_deg,
        "s_plane_m": s_plane,
        "path_integral_m": main["path_m"],
        "x_straight_m": x_s,
        "y_straight_m": y_s,
        "shift_m": main["shift_m"],
        "x_ref_m": x_s + main["shift_m"] * hx,
        "y_ref_m": y_s + main["shift_m"] * hy,
        "excess_ns": main["excess_ns"],
        "straight_ns": main["straight_ns"],
        "t_ref_ns": s_plane / C_M_PER_NS + main["excess_ns"],
        "t_straight_ref_ns": s_plane / C_M_PER_NS + main["straight_ns"],
        "direction_scale_ref": scale,
        "cx_in": cx,
    }


def position_limit(theta_deg: float, shift_m: float) -> float:
    if theta_deg >= 84.0:
        return 0.05
    return max(0.005, 1e-3 * abs(shift_m))


def time_limit(t_out: float, excess_ns: float) -> float:
    half_ulp = float(np.spacing(np.float32(abs(t_out)))) / 2.0
    return 0.01 + half_ulp + 5e-6 * abs(excess_ns)


def convergence(nu, R, h_t, stored_samples) -> float:
    """Max relative change of the reference integrals when the Simpson step is halved."""
    worst = 0.0
    for stored in stored_samples:
        fine = reference_sample(nu, R, h_t, stored, PRODUCTION_INTERVALS)
        coarse = reference_sample(nu, R, h_t, stored, COARSE_INTERVALS)
        for key, floor in (
            ("shift_m", 1e-9),
            ("excess_ns", 1e-9),
            ("straight_ns", 1e-9),
        ):
            change = abs(fine[key] - coarse[key]) / max(abs(fine[key]), floor)
            worst = max(worst, change)
    return worst


def _plot(samples, path) -> None:
    """obdeect minus reference against emission height, one colour per zenith angle."""
    from obdeect_validation.adapters import plotting

    panels = (
        ("x_out_m", "x_ref_m", "position x [m]"),
        ("t_out_ns", "t_ref_ns", "travel time [ns]"),
        ("scale_out", "scale_ref", "direction scale cx_out/cx_in"),
    )
    fig, axes = plotting.figure(1, 3, (13, 3.6))
    zeniths = sorted({round(s["theta_s_deg"], 1) for s in samples})
    for i, zenith in enumerate(zeniths):
        rows = sorted(
            (s for s in samples if round(s["theta_s_deg"], 1) == zenith),
            key=lambda s: s["h_e_m"],
        )
        color = plotting.CYCLE[i % len(plotting.CYCLE)]
        for column, (candidate, reference, label) in enumerate(panels):
            points = [
                (s["h_e_m"] / 1000.0, s[candidate] - s[reference])
                for s in rows
                if s[candidate] is not None and s[reference] is not None
            ]
            if not points:
                continue
            axes[0, column].plot(
                [x for x, _ in points],
                [y for _, y in points],
                "o-",
                color=color,
                label=f"{zenith:g} deg",
            )
            axes[0, column].set(xlabel="emission height [km]", ylabel=label)
    for axis in axes[0]:
        axis.axhline(0.0, color=plotting.ANALYTIC, linestyle="--", linewidth=0.8)
    axes[0, 0].legend(title="zenith at telescope")
    fig.suptitle(
        "obdeect - independent ray tracer (spherical, R = 6371315 m)", fontsize="small"
    )
    plotting.save(fig, path)


def write_stream(
    path, geometry, run_extra: dict, wavelength_nm=DISPERSION_REFERENCE_NM
) -> None:
    """One event per sample: one photon of telescope 1 in the sky frame, no cz column."""
    run = {
        "run_number": 1,
        "observation_level_m": SITE_ALTITUDE_M,
        "telescope_x_m": [0.0],
        "telescope_y_m": [0.0],
        "telescope_z_m": [0.0],
        "refraction_applied": 0,
        "curvature_applied": 0,
        "wavelength_min_nm": wavelength_nm,
        "wavelength_max_nm": wavelength_nm,
        **run_extra,
    }
    from obdeect.stream import Writer

    writer = Writer(str(path), "obdeect-validation")
    writer.write("run", run)
    for number, sample in enumerate(geometry, start=1):
        writer.write("event", {"event_number": number})
        writer.write("telescope", {"telescope_id": 1, "frame": 0})
        writer.write_photons(
            1,
            0,
            0,
            {
                "x": np.array([sample["x"]], dtype=np.float32),
                "y": np.array([sample["y"]], dtype=np.float32),
                "z": np.array([sample["z"]], dtype=np.float32),
                "cx": np.array([sample["cx"]], dtype=np.float32),
                "cy": np.array([sample["cy"]], dtype=np.float32),
                "t": np.zeros(1, dtype=np.float32),
                "w": np.ones(1, dtype=np.float32),
                "wavelength": np.array([wavelength_nm], dtype=np.float32),
            },
        )
        writer.write("telescope_end", {"telescope_id": 1, "photons": 1.0})
        writer.write("event_end", {})
    writer.write("run_end", {})
    writer.close()


def read_outputs(path) -> list[dict]:
    """Per event: the single photon of telescope 1 (G frame) and the emission altitude."""
    from obdeect.stream import Reader

    outputs: list[dict] = []
    for record in Reader(str(path)):
        if record["type"] == "event":
            outputs.append({})
        elif record["type"] == "photons" and record["telescope_id"] == 1:
            columns = record["columns"]
            if len(columns["x"]) != 1:
                raise RuntimeError("expected one output photon per event")
            outputs[-1] = {name: float(columns[name][0]) for name in columns}
    if len(outputs) != 0 and any(not entry for entry in outputs):
        raise RuntimeError("an event has no output photon")
    return outputs


def run_obdeect(
    input_path, output_path, curvature: str, earth_radius_m, profile, dispersion=False
) -> list[dict]:
    from obdeect.propagation import propagate

    propagate(
        input_path,
        output_path,
        site_altitude_m=SITE_ALTITUDE_M,
        seed=SEED,
        profile=profile,
        extinction=False,
        mode="weighted",
        curvature=curvature,
        earth_radius_m=earth_radius_m,
        dispersion_reference_nm=DISPERSION_REFERENCE_NM if dispersion else None,
    )
    return read_outputs(output_path)


def build_geometry(R: float) -> tuple[list[dict], list[dict]]:
    """Samples (labels and stored values) for one Earth radius."""
    samples, stored = [], []
    for height_km in HEIGHTS_KM:
        for zenith in ZENITHS_DEG:
            point, direction = emission_point(height_km * 1000.0, zenith, R)
            values = stored_sample(point, direction)
            label = {"height_km": height_km, "zenith_deg": zenith}
            samples.append(label)
            stored.append(values)
    return samples, stored


def analyse(root: Path) -> tuple[dict, dict]:
    profile = profile_path()
    nu_earth = NuProfile(profile)
    nu_flat = nu_earth

    samples, stored_earth = build_geometry(EARTH_RADIUS_M)
    _, stored_flat = build_geometry(PLANE_RADIUS_M)

    write_stream(root / "earth-spherical.pst", stored_earth, {})
    write_stream(
        root / "earth-spherical-refracted.pst", stored_earth, {"refraction_applied": 1}
    )
    write_stream(root / "flat-spherical.pst", stored_flat, {})
    write_stream(root / "flat-none.pst", stored_flat, {})

    out_earth = run_obdeect(
        root / "earth-spherical.pst",
        root / "earth-spherical-out.pst",
        "spherical",
        EARTH_RADIUS_M,
        profile,
    )
    out_straight = run_obdeect(
        root / "earth-spherical-refracted.pst",
        root / "earth-straight-out.pst",
        "spherical",
        EARTH_RADIUS_M,
        profile,
    )
    out_dispersion = run_obdeect(
        root / "earth-spherical.pst",
        root / "earth-dispersion-out.pst",
        "spherical",
        EARTH_RADIUS_M,
        profile,
        dispersion=True,
    )
    dispersed_outputs = {DISPERSION_REFERENCE_NM: out_dispersion}
    for wavelength_nm in (300.0, 700.0):
        input_path = root / f"earth-dispersion-{int(wavelength_nm)}.pst"
        write_stream(input_path, stored_earth, {}, wavelength_nm=wavelength_nm)
        dispersed_outputs[wavelength_nm] = run_obdeect(
            input_path,
            root / f"earth-dispersion-{int(wavelength_nm)}-out.pst",
            "spherical",
            EARTH_RADIUS_M,
            profile,
            dispersion=True,
        )
    out_flat = run_obdeect(
        root / "flat-spherical.pst",
        root / "flat-spherical-out.pst",
        "spherical",
        PLANE_RADIUS_M,
        profile,
    )
    out_none = run_obdeect(
        root / "flat-none.pst", root / "flat-none-out.pst", "none", None, profile
    )

    rows, plane_rows = [], []
    scale_errors, time_norm, straight_norm, plane_norm = [], [], [], []
    group_delay_norm = []
    dispersion_time_norm = []
    dispersion_position_norm = []
    x_norm, shift_abs = [], []
    for index, (label, stored) in enumerate(zip(samples, stored_earth)):
        ref = reference_sample(nu_earth, EARTH_RADIUS_M, SITE_ALTITUDE_M, stored)
        out = out_earth[index]
        out_s = out_straight[index]
        out_d = out_dispersion[index]
        dx = out["x"] - ref["x_ref_m"]
        tol_x = position_limit(ref["theta_s_deg"], ref["shift_m"])
        dt = out["t"] - ref["t_ref_ns"]
        tol_t = time_limit(out["t"], ref["excess_ns"])
        dt_s = out_s["t"] - ref["t_straight_ref_ns"]
        tol_s = time_limit(out_s["t"], ref["straight_ns"])
        group_delay_expected = (
            ciddor_group_factor(DISPERSION_REFERENCE_NM) - 1.0
        ) * ref["straight_ns"]
        group_delay_error = (out_d["t"] - out["t"]) - group_delay_expected
        group_delay_norm.append(
            abs(group_delay_error)
            / (time_limit(out_d["t"], ref["straight_ns"]) + tol_t + tol_s)
        )
        for wavelength_nm, outputs in dispersed_outputs.items():
            if wavelength_nm == DISPERSION_REFERENCE_NM:
                continue
            factor = ciddor_refractivity_e8(wavelength_nm) / ciddor_refractivity_e8(
                DISPERSION_REFERENCE_NM
            )
            group_factor = ciddor_group_factor(wavelength_nm) * factor

            def scaled_profile(height, k=factor):
                return k * nu_earth(height)

            scaled_ref = reference_sample(
                scaled_profile, EARTH_RADIUS_M, SITE_ALTITUDE_M, stored
            )
            time_expected = (
                scaled_ref["t_ref_ns"] + (group_factor - factor) * ref["straight_ns"]
            )
            time_error = outputs[index]["t"] - time_expected
            time_tolerance = max(
                time_limit(outputs[index]["t"], scaled_ref["excess_ns"]),
                1e-3 * abs(scaled_ref["excess_ns"]),
            )
            dispersion_time_norm.append(abs(time_error) / time_tolerance)
            position_error = outputs[index]["x"] - scaled_ref["x_ref_m"]
            position_tolerance = max(
                position_limit(scaled_ref["theta_s_deg"], scaled_ref["shift_m"]),
                1e-3 * abs(scaled_ref["shift_m"]),
            )
            dispersion_position_norm.append(abs(position_error) / position_tolerance)
        x_norm.append(abs(dx) / tol_x)
        shift_abs.append(abs(dx))
        time_norm.append(abs(dt) / tol_t)
        straight_norm.append(abs(dt_s) / tol_s)
        # The direction scale is undefined for a vertical photon (cx = 0): recorded as None.
        scale_out = None
        if stored["cx"] != 0.0:
            scale_out = out["cx"] / stored["cx"]
            scale_errors.append(abs(scale_out - ref["direction_scale_ref"]))
        rows.append(
            {
                **label,
                "theta_s_deg": ref["theta_s_deg"],
                "h_e_m": ref["h_e_m"],
                "s_plane_m": ref["s_plane_m"],
                "path_integral_m": ref["path_integral_m"],
                "shift_ref_m": ref["shift_m"],
                "x_out_m": out["x"],
                "x_ref_m": ref["x_ref_m"],
                "dx_m": dx,
                "tol_x_m": tol_x,
                "excess_ref_ns": ref["excess_ns"],
                "t_out_ns": out["t"],
                "t_ref_ns": ref["t_ref_ns"],
                "dt_ns": dt,
                "tol_t_ns": tol_t,
                "straight_ref_ns": ref["straight_ns"],
                "t_straight_out_ns": out_s["t"],
                "t_straight_ref_ns": ref["t_straight_ref_ns"],
                "straight_dt_ns": dt_s,
                "dispersion_reference_nm": DISPERSION_REFERENCE_NM,
                "group_delay_expected_ns": group_delay_expected,
                "group_delay_observed_ns": out_d["t"] - out["t"],
                "group_delay_error_ns": group_delay_error,
                "dispersion_wavelengths_checked_nm": [300.0, 500.0, 700.0],
                "x_straight_out_m": out_s["x"],
                "x_straight_ref_m": ref["x_straight_m"],
                "cx_in": stored["cx"],
                "cx_out": out["cx"],
                "cx_out_straight": out_s["cx"],
                "scale_out": scale_out,
                "scale_ref": ref["direction_scale_ref"],
                "zem_out_m": out["zem"],
                "zem_error_m": out["zem"] - (ref["h_e_m"]),
            }
        )
        plane = out_flat[index]
        plane_none = out_none[index]
        ref_flat = reference_sample(
            nu_flat, PLANE_RADIUS_M, SITE_ALTITUDE_M, stored_flat[index]
        )
        tol_flat_x = position_limit(ref_flat["theta_s_deg"], ref_flat["shift_m"])
        tol_flat_t = time_limit(plane["t"], ref_flat["excess_ns"])
        plane_rows.append(
            {
                **label,
                "theta_s_deg": ref_flat["theta_s_deg"],
                "x_flat_out_m": plane["x"],
                "x_none_out_m": plane_none["x"],
                "x_tracer_m": ref_flat["x_ref_m"],
                "t_flat_out_ns": plane["t"],
                "t_none_out_ns": plane_none["t"],
                "t_tracer_ns": ref_flat["t_ref_ns"],
                "diff_none_x_m": plane["x"] - plane_none["x"],
                "diff_none_t_ns": plane["t"] - plane_none["t"],
                "diff_tracer_x_m": plane["x"] - ref_flat["x_ref_m"],
                "diff_tracer_t_ns": plane["t"] - ref_flat["t_ref_ns"],
            }
        )
        plane_norm.extend(
            [
                abs(plane["x"] - plane_none["x"]) / tol_flat_x,
                abs(plane["t"] - plane_none["t"]) / tol_flat_t,
                abs(plane["x"] - ref_flat["x_ref_m"]) / tol_flat_x,
                abs(plane["t"] - ref_flat["t_ref_ns"]) / tol_flat_t,
            ]
        )

    convergence_change = max(
        convergence(nu_earth, EARTH_RADIUS_M, SITE_ALTITUDE_M, stored_earth),
        convergence(nu_flat, PLANE_RADIUS_M, SITE_ALTITUDE_M, stored_flat),
    )
    _plot(rows, root / "refraction-curved-ray.png")
    metrics = {
        "sample_count": len(rows),
        "max_shift_error_m": max(shift_abs),
        "max_shift_error_normalised": max(x_norm),
        "max_time_error_ns_normalised": max(time_norm),
        "max_direction_scale_error": max(scale_errors),
        "straight_delay_max_error_ns_normalised": max(straight_norm),
        "reference_wavelength_group_delay_max_error_normalised": max(group_delay_norm),
        "dispersion_time_max_error_normalised": max(dispersion_time_norm),
        "dispersion_position_max_error_normalised": max(dispersion_position_norm),
        "plane_limit_max_difference": max(plane_norm),
        "reference_convergence_relative_change": convergence_change,
    }
    details = {
        "profile": str(profile),
        "site_altitude_m": SITE_ALTITUDE_M,
        "earth_radius_m": EARTH_RADIUS_M,
        "plane_radius_m": PLANE_RADIUS_M,
        "intervals": {"production": PRODUCTION_INTERVALS, "coarse": COARSE_INTERVALS},
        "limits": LIMITS,
        "samples": rows,
        "plane_limit": plane_rows,
    }
    return metrics, details


if __name__ == "__main__":
    common.results_main(analyse, __doc__.splitlines()[0])
