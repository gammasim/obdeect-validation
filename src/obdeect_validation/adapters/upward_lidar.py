"""Upward laser light scattered once into an off-axis telescope, against the numerical lidar equation.

Case propagation.upward-lidar. Sky-frame primaries (w = 1, wavelength 355 nm) from a laser emitting
upwards at 20 deg from the zenith (upward configurations) or vertically downwards (downward control),
scattered by obdeect's single scattering (Python API, one primary per call, first_index = primary
number, so the draws of each primary are fixed) into one telescope sphere (radius 12.5 m, centre 30 m
above the observation level, off-axis at 100 m and 1 km horizontal distance).

Reference (numpy and scipy only, no obdeect code). Per unit path length s of the primary, with
z(s) the height above the observation level and h = site + z the altitude:

    Y = integral over s of  T0(s) * [beta_mol(h) p_R(mu) + omega beta_aer(h) p_HG(mu)] * Omega(D) * T(P)  ds

- T0(s) = exp(-|tau(h0) - tau(h)| / |cz|): the primary's survival, tau the scattering table's vertical
  optical depth (molecular + aerosol), log10-altitude linear, zero slope above the top node;
- beta = -d tau / dh (the extinction per metre of path equals the vertical one for a slant path);
- mu is the cosine between the primary direction and the direction from P to the aim point, averaged
  over the aim point on the disc of radius R perpendicular to the line of sight through the telescope
  centre C (obdeect draws its aim point uniformly on that disc; the expectation of its weight is the
  disc average, not the centre value; the difference is reported in details);
- Omega(D) = 2 pi (1 - sqrt(1 - (R/D)^2)), D = |C - P|;
- T(P) = exp(-|tau_T(h_P) - tau_T(h_C)| D / |h_P - h_C|), tau_T the transmission table (atm_trans).

Emission: obdeect does not scatter primaries at or below the observation level (scattering.hpp, the
primary check), so the upward laser starts 1 m above it. The reference starts there too; the
0 to 1 m segment of a laser at the ground would add 0.46 % of the yield (not simulated by obdeect).

Path: upward primaries from the emission to the top node of the scattering table (60 km a.s.l.),
downward primaries to the observation level. Arrival time at the sphere centre (no air profile, so
vacuum): t = (s + D) / c. Each time bin of arrival time gives the reference by the same integral over
the s-interval where s + D lies in the bin (s + D is monotonic).

Checks: total yield pull per configuration; per-bin chi2 of the binned time distribution (the bins
with a positive reference; chi2/ndf is recorded, the gate is the chi2 p-value flag, see
chi2_p_below_limit) and largest per-bin pull; propagation of the scattered photons (sky to G, extinction off, since T is already
in their weight): photons reaching the telescope plane with cz > 0 from below it and with cz < 0 from
above it, and the moving_away count, exact against the scattering points.
"""

from __future__ import annotations

import argparse
import dataclasses
import math
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common, plotting
from obdeect_validation.adapters import upward_common as uc

SEED = 20261122
WAVELENGTH_NM = 355
SCATTERING = uc.ROOT / "fixtures/scattering_exponential.dat"
TRANSMISSION = uc.TRANSMISSION
SITE_M = uc.SITE_M
RADIUS_M = 12.5
TELESCOPE_Z_M = 30.0
EMIT_Z_M = 1.0  # obdeect scatters primaries only above the observation level
DOWN_START_Z_M = 48000.0
TOP_Z_M = 60000.0 - SITE_M  # top node of the scattering table (60 km a.s.l.)
SCATTER_POINTS = 128
CHUNK_ROWS = 20000  # photon chunks are limited to 2^16 rows
GRID_POINTS = 200000
DISC_RADIAL = 12
DISC_AZIMUTHAL = 24
DEPOLARIZATION = 0.0  # the fixture has no depolarization line
LN10 = math.log(10.0)
C_M_NS = uc.C_M_NS

UP_EDGES_NS = [
    0.0,
    500.0,
    1000.0,
    2000.0,
    4000.0,
    8000.0,
    16000.0,
    32000.0,
    64000.0,
    128000.0,
    256000.0,
]
DOWN_EDGES_NS = [159950.0 + 50.0 * k for k in range(12)]


@dataclasses.dataclass(frozen=True)
class Config:
    key: str
    upward: bool
    tilt_deg: float
    distance_m: float
    start_z_m: float
    n_primaries: int
    edges_ns: tuple

    def direction(self) -> tuple[np.float32, np.float32, np.float32]:
        """Stored float32 cx, cy, cz (cz signed) of the primary direction."""
        sin_t = math.sin(math.radians(self.tilt_deg))
        cos_t = math.cos(math.radians(self.tilt_deg))
        cz = cos_t if self.upward else -cos_t
        return np.float32(sin_t), np.float32(0.0), np.float32(cz)

    def telescope(self) -> tuple[float, float, float]:
        return (self.distance_m, 0.0, TELESCOPE_Z_M)


CONFIGS = (
    Config("up-100m", True, 20.0, 100.0, EMIT_Z_M, 4000, tuple(UP_EDGES_NS)),
    Config("up-1km", True, 20.0, 1000.0, EMIT_Z_M, 4000, tuple(UP_EDGES_NS)),
    Config("down-100m", False, 0.0, 100.0, DOWN_START_Z_M, 4000, tuple(DOWN_EDGES_NS)),
)


def version() -> str:
    return common.package_version("obdeect-propagation")


# ------------------------------------------------------------------------- references


class Scattering:
    """Scattering table at one wavelength: tau_mol, tau_aer, ssa, g against altitude a.s.l.

    Each quantity is linear in log10(altitude) between the nodes and constant beyond the ends
    (np.interp), as documented in docs/propagation/scattering.md; the extinction per metre is the
    exact slope of that interpolant.
    """

    def __init__(self, path=SCATTERING, wavelength_nm=WAVELENGTH_NM):
        rows = []
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if fields[0] in ("angles_deg", "phase", "depolarization"):
                raise ValueError("phase rows are not used by this reference")
            if int(fields[0]) == wavelength_nm:
                rows.append([float(v) for v in fields[1:]])
        table = np.asarray(rows, dtype=np.float64)
        if table.size == 0:
            raise ValueError(f"no rows at {wavelength_nm} nm")
        self.h_m = table[:, 0] * 1000.0
        self.logh = np.log10(self.h_m)
        if np.any(np.diff(self.logh) <= 0.0):
            raise ValueError("altitudes must be strictly ascending")
        self.tau_mol = table[:, 1]
        self.tau_aer = table[:, 2]
        self.ssa = table[:, 3]
        self.g = table[:, 4]
        self.bottom_m = self.h_m[0]
        self.top_m = self.h_m[-1]

    def _value(self, values, h):
        return np.interp(np.log10(h), self.logh, values)

    def _slope(self, values, h):
        """d value / d h [1/m] of the log10-linear interpolant; zero outside the table."""
        x = np.log10(h)
        i = np.clip(
            np.searchsorted(self.logh, x, side="right") - 1, 0, self.logh.size - 2
        )
        slope = (values[i + 1] - values[i]) / (self.logh[i + 1] - self.logh[i])
        inside = (h > self.bottom_m) & (h < self.top_m)
        return np.where(inside, slope / (h * LN10), 0.0)

    def tau(self, h):
        return self._value(self.tau_mol, h) + self._value(self.tau_aer, h)

    def state(self, h):
        return {
            "beta_mol": -self._slope(self.tau_mol, h),
            "beta_aer": -self._slope(self.tau_aer, h),
            "ssa": self._value(self.ssa, h),
            "g": self._value(self.g, h),
        }


def rayleigh_phase(mu):
    gamma = DEPOLARIZATION / (2.0 - DEPOLARIZATION)
    return (
        3.0 / (8.0 * math.pi * (2.0 + gamma)) * ((1.0 + gamma) + (1.0 - gamma) * mu**2)
    )


def henyey_greenstein_phase(mu, g):
    return (1.0 - g**2) / (4.0 * math.pi * (1.0 + g**2 - 2.0 * g * mu) ** 1.5)


def disc_nodes():
    """Gauss-Legendre in radius (weight r) times the uniform azimuth: nodes and weights of the unit disc."""
    x, w = np.polynomial.legendre.leggauss(DISC_RADIAL)
    r = 0.5 * (x + 1.0)
    wr = 0.5 * w * r
    phi = 2.0 * math.pi * (np.arange(DISC_AZIMUTHAL) + 0.5) / DISC_AZIMUTHAL
    rr = np.repeat(r, DISC_AZIMUTHAL)
    pp = np.tile(phi, DISC_RADIAL)
    weights = np.repeat(wr, DISC_AZIMUTHAL) / (wr.sum() * DISC_AZIMUTHAL)
    return rr * np.cos(pp), rr * np.sin(pp), weights


class Transmission:
    """Telescope path transmission T(P -> C) from the transmission table (numpy reference)."""

    def __init__(self, wavelength_nm=WAVELENGTH_NM):
        self.table = uc.Transmission(TRANSMISSION)
        self.wavelength_nm = wavelength_nm

    def path(self, z_point, z_centre, distance):
        """exp(-|tau(h_P) - tau(h_C)| * D / |h_P - h_C|), heights above the observation level."""
        tau_p = self.table.tau(self.wavelength_nm, SITE_M + z_point)
        tau_c = self.table.tau(self.wavelength_nm, SITE_M + z_centre)
        dz = np.abs(z_point - z_centre)
        od = np.abs(tau_p - tau_c)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(dz > 0.0, od / np.where(dz > 0.0, dz, 1.0), 0.0)
        return np.exp(-ratio * distance)


def _geometry(cfg: Config, s, scat: Scattering, tel: Transmission, disc=True):
    cx, cy, cz = (float(v) for v in cfg.direction())
    cz_abs = math.sqrt(max(0.0, 1.0 - cx * cx - cy * cy))
    d = np.array([cx, cy, math.copysign(cz_abs, cz)])
    P = np.stack([d[0] * s, d[1] * s, cfg.start_z_m + d[2] * s], axis=1)
    C = np.asarray(cfg.telescope())
    rel = C[None, :] - P
    dist = np.linalg.norm(rel, axis=1)
    u = rel / dist[:, None]
    tau0 = scat.tau(np.array([SITE_M + cfg.start_z_m]))[0]
    tau_s = scat.tau(SITE_M + P[:, 2])
    attenuation = np.exp(-np.abs(tau_s - tau0) / cz_abs)
    st = scat.state(SITE_M + P[:, 2])
    if disc:
        e1 = np.cross(u, np.array([0.0, 1.0, 0.0]))
        e1 /= np.linalg.norm(e1, axis=1)[:, None]
        e2 = np.cross(u, e1)
        dxs, dys, weights = disc_nodes()
        phase = np.zeros(s.size)
        for a, b, w in zip(dxs, dys, weights):
            v = rel + RADIUS_M * (a * e1 + b * e2)
            mu = (v @ d) / np.linalg.norm(v, axis=1)
            phase += w * (
                st["beta_mol"] * rayleigh_phase(mu)
                + st["ssa"] * st["beta_aer"] * henyey_greenstein_phase(mu, st["g"])
            )
    else:
        mu = (rel @ d) / dist
        phase = st["beta_mol"] * rayleigh_phase(mu) + st["ssa"] * st[
            "beta_aer"
        ] * henyey_greenstein_phase(mu, st["g"])
    ratio = RADIUS_M / dist
    solid = 2.0 * math.pi * (1.0 - np.sqrt(np.clip(1.0 - ratio**2, 0.0, None)))
    transmitted = tel.path(P[:, 2], C[2], dist)
    integrand = attenuation * phase * solid * transmitted
    arrival_m = s + dist
    return integrand, arrival_m, P


def path_end(cfg: Config) -> float:
    vz = float(cfg.direction()[2])
    z_end = TOP_Z_M if cfg.upward else 0.0
    return (z_end - cfg.start_z_m) / vz


def reference(
    cfg: Config, scat: Scattering, tel: Transmission, grid=GRID_POINTS, disc=True
):
    """Reference yield and per-time-bin yields (per primary of weight 1) on a trapezoid grid."""
    s_end = path_end(cfg)
    # Uniform grid plus the nodes of both tables (the extinction per metre and the telescope path
    # jump there): trapezoid on smooth pieces.
    vz = float(cfg.direction()[2])
    z_nodes = np.concatenate([scat.h_m, tel.table.h1]) - SITE_M
    s_nodes = (z_nodes - cfg.start_z_m) / vz
    s_nodes = s_nodes[(s_nodes > 0.0) & (s_nodes < s_end)]
    s = np.union1d(np.linspace(0.0, s_end, grid + 1), s_nodes)
    f, arrival_m, _ = _geometry(cfg, s, scat, tel, disc=disc)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(s))])
    arrival_m = np.maximum.accumulate(arrival_m)
    edges_m = np.asarray(cfg.edges_ns) * C_M_NS
    s_edges = np.interp(edges_m, arrival_m, s)
    bins = np.interp(s_edges, s, cum)
    return {
        "total": float(cum[-1]),
        "bins": np.diff(bins),
        "s": s,
        "density": f,
        "arrival_m": arrival_m,
        "cum": cum,
    }


# ------------------------------------------------------------------------- obdeect


def obdeect_run(cfg: Config, scat_path: Path, tr_path: Path, n: int | None = None):
    """Scatter one primary per call (first_index = primary number); per-primary totals and time bins."""
    from obdeect.propagation import ScatteringTable, TransmissionTable, scatter

    table = ScatteringTable(str(scat_path))
    transmission = TransmissionTable(str(tr_path), site_altitude_m=SITE_M)
    n = cfg.n_primaries if n is None else n
    cx, cy, cz = cfg.direction()
    xt, yt, zt = cfg.telescope()
    edges = np.asarray(cfg.edges_ns)
    totals = np.zeros(n)
    binned = np.zeros((n, edges.size - 1))
    columns = {k: [] for k in ("x", "y", "z", "cx", "cy", "cz", "t", "w", "wavelength")}
    for i in range(n):
        one = {
            "x": np.zeros(1, np.float32),
            "y": np.zeros(1, np.float32),
            "z": np.full(1, cfg.start_z_m, np.float32),
            "cx": np.full(1, cx, np.float32),
            "cy": np.full(1, cy, np.float32),
            "cz": np.full(1, cz, np.float32),
            "t": np.zeros(1, np.float32),
            "w": np.ones(1, np.float32),
            "wavelength": np.full(1, WAVELENGTH_NM, np.float32),
        }
        out = scatter(
            one,
            table,
            transmission,
            [(1, xt, yt, zt, RADIUS_M)],
            SEED,
            1,
            i,
            SCATTER_POINTS,
        )[1]
        w = np.asarray(out["w"], dtype=np.float64)
        dist = np.sqrt(
            (np.asarray(out["x"], np.float64) - xt) ** 2
            + (np.asarray(out["y"], np.float64) - yt) ** 2
            + (np.asarray(out["z"], np.float64) - zt) ** 2
        )
        arrival_ns = np.asarray(out["t"], np.float64) + dist / C_M_NS
        totals[i] = w.sum()
        binned[i] = np.histogram(arrival_ns, bins=edges, weights=w)[0]
        for name in columns:
            columns[name].append(np.asarray(out[name]))
    stacked = {
        k: np.concatenate(v) if v else np.zeros(0, np.float32)
        for k, v in columns.items()
    }
    return {"totals": totals, "binned": binned, "rows": stacked}


def count_check(cfg: Config, rows: dict, root: Path) -> dict:
    """Scattered photons through the sky to G stage (no extinction): reaching counts by cz sign."""
    from obdeect.propagation import propagate

    zt = TELESCOPE_Z_M
    z = rows["z"].astype(np.float64)
    cz = rows["cz"].astype(np.float64)
    up_reach = int(np.sum((cz > 0) & (z < zt)))
    down_reach = int(np.sum((cz < 0) & (z > zt)))
    away = int(z.size - up_reach - down_reach)
    stream = uc.write_stream(
        root / f"{cfg.key}-scattered.pst",
        [cfg.telescope()],
        [
            (
                1,
                {
                    k: rows[k]
                    for k in ("x", "y", "z", "cx", "cy", "cz", "t", "w", "wavelength")
                },
            )
        ],
        wavelength_nm=(WAVELENGTH_NM, WAVELENGTH_NM),
        chunk_rows=CHUNK_ROWS,
    )
    out_path = root / f"{cfg.key}-plane.pst"
    summary = propagate(
        stream,
        out_path,
        site_altitude_m=SITE_M,
        seed=1,
        mode="weighted",
        extinction=False,
    )
    out = uc.read_rows(out_path)
    on_plane = out["telescope"] == 1
    out_cz = out["cz"][on_plane].astype(np.float64)
    return {
        "expected_up": up_reach,
        "expected_down": down_reach,
        "expected_away": away,
        "obdeect_up": int(np.sum(out_cz > 0)),
        "obdeect_down": int(np.sum(out_cz < 0)),
        "obdeect_horizontal": int(np.sum(out_cz == 0)),
        "dropped_by_reason": {
            k: int(v) for k, v in summary["dropped_by_reason"].items()
        },
        "photons_in": int(summary["photons_in"]),
        "photons_out": int(summary["photons_out"]),
    }


# ------------------------------------------------------------------------- analysis


CHI2_P_LIMIT = (
    1e-4  # significance of the time-shape flag, the same level as the KS gate
)


def chi2_bins(obs, obs_err, ref):
    """chi2, the used bins (positive reference and error) and their pulls."""
    used = np.flatnonzero((ref > 0.0) & (obs_err > 0.0))
    pulls = (obs[used] - ref[used]) / obs_err[used]
    return float(np.sum(pulls**2)), used, pulls


def chi2_p_value(chi2, ndf):
    """Upper-tail p-value of chi2 with ndf degrees of freedom (the reference is fixed: no fit)."""
    from scipy.stats import chi2 as chi2_law

    return float(chi2_law.sf(chi2, ndf))


def chi2_p_below_limit(chi2, ndf):
    """1 if the time-shape chi2 is significant at CHI2_P_LIMIT, else 0.

    Justification (same rule as the KS shape gate of source_camera): if obdeect and the reference
    agree, the p-value is uniform on (0, 1), so a flag has a false-alarm rate of CHI2_P_LIMIT per
    configuration; the case has three configurations, about 3e-4 in all. The reference has no free
    parameters, so ndf is the number of used bins. Only a large chi2 (obdeect and the reference
    disagree) is flagged; a small chi2 would mean overestimated errors and is not a model error.
    """
    return int(chi2_p_value(chi2, ndf) < CHI2_P_LIMIT)


def analyse(root: Path, n_override: int | None = None, grid: int = GRID_POINTS):
    scat = Scattering()
    tel = Transmission()
    scat_path = SCATTERING
    tr_path = TRANSMISSION
    metrics = {}
    details = {
        "seed": SEED,
        "obdeect_propagation_version": version(),
        "scattering_table": str(SCATTERING),
        "transmission": str(TRANSMISSION),
        "wavelength_nm": WAVELENGTH_NM,
        "site_altitude_m": SITE_M,
        "telescope_radius_m": RADIUS_M,
        "telescope_z_m": TELESCOPE_Z_M,
        "emit_z_m": EMIT_Z_M,
        "scatter_points": SCATTER_POINTS,
        "reference": "trapezoid on %d points of the path (see module docstring)" % grid,
        "configurations": {},
    }
    arrays = {}
    for cfg in CONFIGS:
        k = cfg.key.replace("-", "_")
        n = cfg.n_primaries if n_override is None else n_override
        ref = reference(cfg, scat, tel, grid=grid)
        ref_fine = reference(cfg, scat, tel, grid=2 * grid)
        ref_centre = reference(cfg, scat, tel, grid=grid, disc=False)
        quad_rel = abs(ref["total"] - ref_fine["total"]) / ref["total"]
        obd = obdeect_run(cfg, scat_path, tr_path, n=n)
        totals = obd["totals"]
        mean = float(totals.mean())
        err = float(totals.std(ddof=1) / math.sqrt(n))
        pull = (mean - ref["total"]) / err
        obs_bin = obd["binned"].mean(axis=0)
        err_bin = obd["binned"].std(axis=0, ddof=1) / math.sqrt(n)
        chi2, used, pulls = chi2_bins(obs_bin, err_bin, ref["bins"])
        ndf = int(used.size)
        empty = ref["bins"] == 0.0
        empty_content = float(np.sum(np.abs(obs_bin[empty])))
        counts = count_check(cfg, obd["rows"], root)
        metrics[f"{k}_total_yield_pull"] = float(pull)
        metrics[f"{k}_time_chi2_ndf"] = chi2 / ndf if ndf else math.nan
        metrics[f"{k}_time_chi2_p_below_limit"] = (
            chi2_p_below_limit(chi2, ndf) if ndf else math.nan
        )
        metrics[f"{k}_time_max_abs_pull"] = (
            float(np.max(np.abs(pulls))) if ndf else math.nan
        )
        metrics[f"{k}_plane_up_count_error"] = (
            counts["obdeect_up"] - counts["expected_up"]
        )
        metrics[f"{k}_plane_down_count_error"] = (
            counts["obdeect_down"] - counts["expected_down"]
        )
        metrics[f"{k}_moving_away_error"] = (
            counts["dropped_by_reason"].get("moving_away", 0) - counts["expected_away"]
        )
        metrics[f"{k}_reference_quadrature_relative"] = float(quad_rel)
        details["configurations"][cfg.key] = {
            "primaries": n,
            "upward": cfg.upward,
            "tilt_deg": cfg.tilt_deg,
            "distance_m": cfg.distance_m,
            "start_z_m": cfg.start_z_m,
            "direction_stored": [float(v) for v in cfg.direction()],
            "edges_ns": list(cfg.edges_ns),
            "yield_obdeect": mean,
            "yield_obdeect_error": err,
            "yield_reference": ref["total"],
            "yield_reference_fine_grid": ref_fine["total"],
            "yield_reference_centre_aim": ref_centre["total"],
            "aim_disc_effect_relative": (ref["total"] - ref_centre["total"])
            / ref["total"],
            "bins_obdeect": obs_bin.tolist(),
            "bins_obdeect_error": err_bin.tolist(),
            "bins_reference": ref["bins"].tolist(),
            "empty_bin_content": empty_content,
            "bins_used": used.tolist(),
            "bins_ndf": ndf,
            "bins_chi2": chi2,
            "bins_chi2_ndf": chi2 / ndf if ndf else math.nan,
            "bins_chi2_p_value": chi2_p_value(chi2, ndf) if ndf else math.nan,
            "bins_pulls": pulls.tolist(),
            "count_check": counts,
        }
        arrays[cfg.key] = {
            "edges": np.asarray(cfg.edges_ns),
            "obs": obs_bin,
            "obs_err": err_bin,
            "ref_bins": ref["bins"],
            "arrival_us": ref["arrival_m"] / C_M_NS / 1000.0,
            "density": ref["density"],
            "s": ref["s"],
        }
    metrics["empty_bin_content"] = float(
        sum(details["configurations"][c.key]["empty_bin_content"] for c in CONFIGS)
    )
    metrics["sample_count"] = int(
        sum(c.n_primaries if n_override is None else n_override for c in CONFIGS)
    )
    metrics["reference_quadrature_relative"] = max(
        metrics[f"{c.key.replace('-', '_')}_reference_quadrature_relative"]
        for c in CONFIGS
    )
    details["_arrays"] = arrays
    return metrics, details


def plot(path: Path, details: dict) -> None:
    arrays = details.pop("_arrays")
    fig, axes = plotting.figure(2, len(CONFIGS), (4.0 * len(CONFIGS) + 1.0, 6.4))
    for col, cfg in enumerate(CONFIGS):
        a = arrays[cfg.key]
        conf = details["configurations"][cfg.key]
        edges_us = a["edges"] / 1000.0
        width = np.diff(edges_us)
        centres = 0.5 * (edges_us[1:] + edges_us[:-1])
        top = axes[0, col]
        top.errorbar(
            centres,
            a["obs"] / width,
            yerr=a["obs_err"] / width,
            fmt="o",
            ms=3.5,
            color=plotting.NATIVE,
            label="obdeect",
        )
        top.stairs(
            a["ref_bins"] / width,
            edges_us,
            color="black",
            linestyle="--",
            label="reference",
        )
        if cfg.upward:
            top.set_xscale("symlog", linthresh=1.0)
        top.set(
            xlabel="arrival time at sphere centre [us]",
            ylabel="scattered photons per primary per us",
            title=cfg.key,
        )
        top.legend(fontsize="small")
        bottom = axes[1, col]
        used = np.asarray(conf["bins_used"], dtype=int)
        bottom.bar(
            centres[used],
            np.asarray(conf["bins_pulls"]),
            width=0.8 * width[used],
            color=plotting.NATIVE,
        )
        if cfg.upward:
            bottom.set_xscale("symlog", linthresh=1.0)
        bottom.axhline(0.0, color="black", linestyle="--", linewidth=0.8)
        for limit in (-4.0, 4.0):
            bottom.axhline(limit, color="black", linestyle=":", linewidth=0.8)
        bottom.set(xlabel="arrival time bin centre [us]", ylabel="pull")
    plotting.save(fig, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args()
    root = Path.cwd()
    metrics, details = analyse(root)
    plot(root / "upward-lidar.png", details)
    common.write_results(metrics, details, root)


if __name__ == "__main__":
    main()
