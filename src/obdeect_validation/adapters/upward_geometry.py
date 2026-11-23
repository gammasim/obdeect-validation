"""Upward photons in plane mode: straight-line crossing point and time from emission below the plane.

Case propagation.upward-geometry. No refraction, no extinction (extinction=False, no profile): the
reference is the straight line from the emission point (x0, y0, z0) with the stored direction to the
telescope plane z = h_t, t = t0 + path / c, path = (h_t - z0) / |cz|. Directions run from vertical
(0.5 deg) to near-horizontal (89.5 deg) and all azimuths; emission heights from the observation level
to just below the plane. Errors are in units of the float32 spacing of the output (0.5 ulp is the
rounding of a stored value). Also: cx, cy, cz are kept bit for bit, and the G -> sky conversion
(to_sky) returns the emission positions and times to float32 rounding.
"""

from __future__ import annotations

import argparse
import importlib.metadata
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common, plotting
from obdeect_validation.adapters import upward_common as uc

SEED = 20261124
ZENITHS_DEG = (0.5, 10.0, 30.0, 60.0, 80.0, 89.0, 89.5)
AZIMUTHS_DEG = (0.0, 60.0, 120.0, 240.0, 300.0)
EMISSION_M = (
    0.0,
    150.0,
    420.0,
    499.0,
)  # sky z above the observation level, below the plane
PLANE_Z_M = 500.0  # telescope plane above the observation level
TELESCOPE = (0.0, 0.0, PLANE_Z_M)


def version() -> str:
    try:
        return importlib.metadata.version("obdeect-propagation")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def samples(seed: int = SEED) -> dict:
    """Upward photons (one per zenith, azimuth and emission height), float32 as stored."""
    rng = np.random.default_rng(seed)
    rows = []
    for zenith in ZENITHS_DEG:
        for azimuth in AZIMUTHS_DEG:
            cx, cy, cz = uc.stored_direction(zenith, azimuth)
            for z in EMISSION_M:
                rows.append(
                    (
                        np.float32(rng.uniform(-500.0, 500.0)),
                        np.float32(rng.uniform(-500.0, 500.0)),
                        np.float32(z),
                        cx,
                        cy,
                        cz,
                        np.float32(rng.uniform(0.0, 50.0)),
                        zenith,
                    )
                )
    names = ("x", "y", "z", "cx", "cy", "cz", "t", "zenith")
    return {name: np.array([row[i] for row in rows]) for i, name in enumerate(names)}


def columns_of(data: dict) -> dict:
    keys = ("x", "y", "z", "cx", "cy", "cz", "t")
    out = {k: np.asarray(data[k], dtype=np.float32) for k in keys}
    out["w"] = np.ones(len(data["x"]), dtype=np.float32)
    out["wavelength"] = np.full(len(data["x"]), 400.0, dtype=np.float32)
    return out


def reference(data: dict) -> dict:
    """Straight-line crossing and time from the stored values (float64)."""
    x0 = data["x"].astype(np.float64)
    y0 = data["y"].astype(np.float64)
    z0 = data["z"].astype(np.float64)
    cx = data["cx"].astype(np.float64)
    cy = data["cy"].astype(np.float64)
    t0 = data["t"].astype(np.float64)
    # Documented convention (refraction.md, "Upward photons"): |cz| = sqrt(1 - cx^2 - cy^2) from the stored
    # cx, cy, for upward rows too; the cz column gives the sign only.
    cz = np.sqrt(1.0 - cx * cx - cy * cy)
    path = uc.path_to_plane(z0, PLANE_Z_M, cz)
    return {
        "x": x0 + path * cx - TELESCOPE[0],
        "y": y0 + path * cy - TELESCOPE[1],
        "t": t0 + path / uc.C_M_NS,
        "path": path,
    }


def run_propagation(root: Path, data: dict):
    from obdeect.propagation import propagate

    telescopes = [TELESCOPE]
    stream = uc.write_stream(
        root / "upward-geometry-in.pst",
        telescopes,
        [(1, columns_of(data))],
        chunk_rows=64,
    )
    out = root / "upward-geometry-out.pst"
    summary = propagate(
        stream, out, extinction=False, site_altitude_m=uc.SITE_M, seed=1
    )
    rows = uc.read_rows(out)
    return summary, rows, out


def roundtrip(root: Path, out: Path, data: dict, g_rows: dict) -> np.ndarray:
    """G -> sky with to_sky; each error divided by its float32 error budget (1 = at the rounding bound).

    to_sky back-projects with s = |zem - h_t| / |cz| from the stored emission altitude zem. The budget
    of the path s: half the spacing of zem over |cz|, plus the change from |cz| = sqrt(1 - cx^2 - cy^2)
    against the stored cz; x, y, t are then the path times |cx|, |cy| or 1 / c, plus one spacing of each
    stored value.
    """
    from obdeect.propagation import to_sky

    back = root / "upward-geometry-back.pst"
    to_sky(out, back)
    rows = uc.read_rows(back)
    back_order = np.argsort(rows["index"])
    g_order = np.argsort(g_rows["index"])
    x_back = rows["x"][back_order].astype(np.float64)
    y_back = rows["y"][back_order].astype(np.float64)
    t_back = rows["t"][back_order].astype(np.float64)
    x_g = g_rows["x"][g_order].astype(np.float64)
    y_g = g_rows["y"][g_order].astype(np.float64)
    t_g = g_rows["t"][g_order].astype(np.float64)
    zem = g_rows["zem"][g_order].astype(np.float64)
    cx = g_rows["cx"][g_order].astype(np.float64)
    cy = g_rows["cy"][g_order].astype(np.float64)
    cz = np.sqrt(np.maximum(1.0 - cx * cx - cy * cy, 0.0))
    zem_half = 0.5 * uc.ulp32(zem)
    x_in = data["x"].astype(np.float64)
    y_in = data["y"].astype(np.float64)
    t_in = data["t"].astype(np.float64)
    # Both stages take |cz| from the stored cx, cy (the same convention), so no magnitude mismatch enters.
    path_budget = zem_half / cz
    x_budget = path_budget * np.abs(cx) + uc.ulp32(x_g) + uc.ulp32(x_in)
    y_budget = path_budget * np.abs(cy) + uc.ulp32(y_g) + uc.ulp32(y_in)
    t_budget = path_budget / uc.C_M_NS + uc.ulp32(t_g) + uc.ulp32(t_in)
    return np.maximum.reduce(
        [
            np.abs(x_back + TELESCOPE[0] - x_in) / x_budget,
            np.abs(y_back + TELESCOPE[1] - y_in) / y_budget,
            np.abs(t_back - t_in) / t_budget,
        ]
    )


def analyse(root: Path) -> tuple[dict, dict]:
    data = samples()
    ref = reference(data)
    summary, rows, out = run_propagation(root, data)
    assert summary["photons_in"] == len(data["x"]), summary
    order = np.argsort(rows["index"])
    x = rows["x"][order].astype(np.float64)
    y = rows["y"][order].astype(np.float64)
    t = rows["t"][order].astype(np.float64)
    # Dropped photons would shift the rows: none is expected, so the order is the input order.
    pos_err_ulp = np.maximum(
        np.abs(x - ref["x"]) / uc.ulp32(ref["x"]),
        np.abs(y - ref["y"]) / uc.ulp32(ref["y"]),
    )
    time_err_ulp = np.abs(t - ref["t"]) / uc.ulp32(ref["t"])
    dir_err = np.max(
        np.abs(
            np.stack([rows["cx"], rows["cy"], rows["cz"]], axis=1).astype(np.float64)
            - np.stack([data["cx"], data["cy"], data["cz"]], axis=1).astype(np.float64)
        )
    )
    rt_ulp = roundtrip(root, out, data, rows)
    per_zenith = {}
    for zenith in ZENITHS_DEG:
        mask = data["zenith"] == zenith
        per_zenith[str(zenith)] = {
            "position_error_ulp": float(pos_err_ulp[mask].max()),
            "time_error_ulp": float(time_err_ulp[mask].max()),
            "path_m_max": float(ref["path"][mask].max()),
        }
    metrics = {
        "sample_count": len(data["x"]),
        "dropped": int(summary["dropped"]),
        "max_position_error_ulp": float(pos_err_ulp.max()),
        "max_time_error_ulp": float(time_err_ulp.max()),
        "max_direction_error": float(dir_err),
        "round_trip_max_error_over_budget": float(rt_ulp.max()),
    }
    details = {
        "seed": SEED,
        "obdeect_propagation_version": version(),
        "zeniths_deg": list(ZENITHS_DEG),
        "azimuths_deg": list(AZIMUTHS_DEG),
        "emission_z_m": list(EMISSION_M),
        "plane_z_m": PLANE_Z_M,
        "summary": summary,
        "per_zenith": per_zenith,
        "reference": "straight line, path = (plane - z0) / |cz| with the stored float32 values; float64",
        "error_unit": "float32 spacing at max(|value|, 1) (m or ns) of the stored value; round trip: ratio to the float32 error budget (<= 1 passes)",
        "longest_path_m": float(ref["path"].max()),
    }
    details["_arrays"] = {
        "zenith": data["zenith"],
        "pos": pos_err_ulp,
        "x": ref["x"],
        "y": ref["y"],
    }
    return metrics, details


def plot(path: Path, details: dict) -> None:
    arrays = details.pop("_arrays")
    fig, axis = plotting.figure(1, 1, (6.0, 4.0))
    axis.plot(
        np.abs(arrays["zenith"]) * 0 + arrays["zenith"],
        arrays["pos"],
        "o",
        markersize=4,
        color=plotting.NATIVE,
        label="obdeect, max over x and y",
    )
    axis.axhline(
        0.5,
        color=plotting.ANALYTIC,
        linestyle="--",
        linewidth=1,
        label="0.5 ulp (rounding)",
    )
    axis.set(
        xlabel="zenith angle [deg]",
        ylabel="position error [float32 ulp]",
        yscale="symlog",
    )
    axis.legend(loc="upper left", fontsize="small")
    plotting.save(fig, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args()
    root = Path.cwd()
    metrics, details = analyse(root)
    plot(root / "upward-geometry.png", details)
    common.write_results(metrics, details, root)


if __name__ == "__main__":
    main()
