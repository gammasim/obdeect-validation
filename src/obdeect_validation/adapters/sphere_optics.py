"""Ray tracing from the telescope sphere: plane rows and sphere entry rows of the same downward photons.

Case optics.sphere-vs-plane. A parallel beam from a star at zenith 20 deg, azimuth 0 (downward photons, the
telescope pointing the same way: the beam is aimed along the optical axis reversed), is propagated twice from
the sky frame with the same random stream: plane mode (G rows on the plane z = 0 through the telescope centre)
and sphere mode (rows at the entry point of each photon's line into the telescope sphere). Both go through
`obdeect ray-tracing trace-stream` (the stream path) with the same optical model and pointing.

Compared, per photon (identity: a unique wavelength per photon, which the stages keep):
- the detected sets are identical;
- camera-frame positions and times agree within float32 rounding of the start point (budget: four float32
  spacings of the stored value, two roundings per path; the two paths start from different rows of one line);
- camera-frame directions agree to float32 rounding.

Models: the segmented fixture of the obdeect tests (fixtures/optical-models/segmented.json, always) and LSTN-01
(simulation models 7.0.0, compiled by obdeect) when a simulation-models checkout is present (VALIDATION_MODELS or
/models, as the other optics cases); without it the LST part is not compared and says so in details.json.

Also: horizontal sphere rows (cz = 0) of the tilted fixture telescope are traced without rejection
(photons_in of trace-stream = number of rows); and the enclosure check rejects a sphere of 0.5 m (the message
names the element) and accepts the sphere of the fixture (9 m, the camera front is 8.52 m from the centre).
"""

from __future__ import annotations

import argparse
import importlib.metadata
import math
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common
from obdeect_validation.adapters import optics_directions as od
from obdeect_validation.adapters import upward_common as uc

SEED = 20261127
ZENITH_DEG = 20.0
AZIMUTH_DEG = 0.0
START_M = 2000.0
FIXTURE = uc.ROOT / "fixtures/optical-models/segmented.json"
FIXTURE_RADIUS_M = 9.0  # camera front at 8.521 m from the centre
FIXTURE_DISK_M = 2.1
LST_RADIUS_M = 35.0  # LST focal length 28 m, camera front beyond it
LST_DISK_M = 12.0
PHOTONS_FIXTURE = 20_000
PHOTONS_LST = 40_000
WAVELENGTH_BASE_NM = 400.0
WAVELENGTH_STEP_NM = 0.01
EPS32 = 2.0**-24


def version() -> str:
    try:
        return importlib.metadata.version("obdeect-ray-tracing")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def star_axis() -> np.ndarray:
    th, ph = math.radians(ZENITH_DEG), math.radians(AZIMUTH_DEG)
    return np.array(
        [math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)]
    )


def parallel_beam(n: int, disk_m: float, seed: int) -> dict:
    """Downward parallel photons from the star direction, start on a disc of radius disk_m at START_M."""
    s = star_axis()
    e1 = (
        np.array([-s[1], s[0], 0.0]) if abs(s[2]) < 0.999 else np.array([1.0, 0.0, 0.0])
    )
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(s, e1)
    rng = np.random.default_rng(seed)
    uv = np.empty((0, 2))
    while len(uv) < n:
        cand = rng.uniform(-1, 1, (2 * n, 2))
        uv = np.concatenate([uv, cand[np.sum(cand * cand, axis=1) <= 1.0]])
    uv = uv[:n]
    pos = START_M * s[None, :] + disk_m * (
        uv[:, :1] * e1[None, :] + uv[:, 1:] * e2[None, :]
    )
    d = -s
    columns = {
        "x": pos[:, 0],
        "y": pos[:, 1],
        "z": pos[:, 2],
        "cx": np.full(n, d[0]),
        "cy": np.full(n, d[1]),
        "cz": np.full(n, d[2]),
        "t": np.zeros(n),
        "w": np.ones(n),
        "wavelength": WAVELENGTH_BASE_NM + WAVELENGTH_STEP_NM * np.arange(n),
    }
    return {k: np.asarray(v, dtype=np.float32) for k, v in columns.items()}


def horizontal_rows(n: int, seed: int) -> dict:
    """Horizontal photons (cz = 0) from the side, heights within the fixture sphere, identities as above."""
    rng = np.random.default_rng(seed)
    z = rng.uniform(-6.0, 6.0, n)
    azimuth = rng.uniform(0.0, 2 * math.pi, n)
    x = 200.0 * np.cos(azimuth)
    y = 200.0 * np.sin(azimuth)
    cx, cy = -np.cos(azimuth), -np.sin(azimuth)
    columns = {
        "x": x,
        "y": y,
        "z": z,
        "cx": cx,
        "cy": cy,
        "cz": np.zeros(n),
        "t": np.zeros(n),
        "w": np.ones(n),
        "wavelength": WAVELENGTH_BASE_NM + WAVELENGTH_STEP_NM * np.arange(n),
    }
    return {k: np.asarray(v, dtype=np.float32) for k, v in columns.items()}


def _write(path: Path, columns: dict, radius_m: float) -> Path:
    from obdeect.stream import Writer

    run = {
        "run_number": 1,
        "observation_level_m": 0.0,
        "telescope_x_m": [0.0],
        "telescope_y_m": [0.0],
        "telescope_z_m": [0.0],
        "telescope_r_m": [radius_m],
        "wavelength_min_nm": 250.0,
        "wavelength_max_nm": 800.0,
    }
    with Writer(str(path), "obdeect validation: sphere optics") as writer:
        writer.write("run", run)
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1, "frame": 0, "array_index": 0})
        writer.write_photons(1, "sky", 0, columns)
        writer.write(
            "telescope_end", {"telescope_id": 1, "photons": float(len(columns["x"]))}
        )
        writer.write("event_end", {})
        writer.write("run_end", {})
    return path


def propagate_pair(
    root: Path, name: str, columns: dict, radius_m: float
) -> tuple[Path, Path]:
    from obdeect.propagation import propagate

    inp = _write(root / f"{name}-sky.pst", columns, radius_m)
    plane = root / f"{name}-plane.pst"
    sphere = root / f"{name}-sphere.pst"
    propagate(inp, plane, extinction=False, site_altitude_m=0.0, seed=1)
    propagate(
        inp,
        sphere,
        extinction=False,
        site_altitude_m=0.0,
        seed=1,
        ground_representation="sphere",
    )
    return plane, sphere


def read_camera(path: Path) -> dict:
    from obdeect.stream import Reader

    columns: dict[str, list] = {}
    for record in Reader(str(path)):
        if record["type"] == "photons":
            for name, values in record["columns"].items():
                columns.setdefault(name, []).append(np.asarray(values))
    return {
        k: np.concatenate(v) if v else np.empty(0, dtype=np.float32)
        for k, v in columns.items()
    }


def trace(model, grid: Path, out: Path, seed: int, pointing) -> dict:
    from obdeect.ray_tracing.stream import trace_stream

    return trace_stream(model, str(grid), str(out), seed=seed, pointing=pointing)


def ids_of(wavelength) -> np.ndarray:
    return np.rint(
        (np.asarray(wavelength, dtype=np.float64) - WAVELENGTH_BASE_NM)
        / WAVELENGTH_STEP_NM
    ).astype(np.int64)


#: Each path rounds the value twice (the stream row, then the camera row): four float32 spacings.
BUDGET_SPACINGS = 4.0


def budget(values, factor: float = BUDGET_SPACINGS) -> np.ndarray:
    v = np.asarray(values, dtype=np.float32)
    return factor * np.spacing(np.maximum(np.abs(v), 1.0).astype(np.float32)).astype(
        np.float64
    )


def compare_camera(plane: dict, sphere: dict) -> dict:
    ida, idb = ids_of(plane["wavelength"]), ids_of(sphere["wavelength"])
    sa, sb = set(ida.tolist()), set(idb.tolist())
    common_ids = sorted(sa & sb)
    mismatch = len(sa ^ sb)
    if not common_ids:
        return {
            "detected_plane": len(sa),
            "detected_sphere": len(sb),
            "set_mismatch": mismatch,
        }
    ia = {v: i for i, v in enumerate(ida.tolist())}
    ib = {v: i for i, v in enumerate(idb.tolist())}
    pa = np.array([ia[v] for v in common_ids])
    pb = np.array([ib[v] for v in common_ids])
    out = {
        "detected_plane": len(sa),
        "detected_sphere": len(sb),
        "set_mismatch": mismatch,
    }
    worst_pos = 0.0
    pos_ratio = 0.0
    for key in ("x", "y"):
        a = plane[key][pa].astype(np.float64)
        b = sphere[key][pb].astype(np.float64)
        worst_pos = max(worst_pos, float(np.max(np.abs(a - b))))
        pos_ratio = max(pos_ratio, float(np.max(np.abs(a - b) / budget(a))))
    dt = np.abs(plane["t"][pa].astype(np.float64) - sphere["t"][pb].astype(np.float64))
    t_ratio = float(np.max(dt / budget(plane["t"][pa].astype(np.float64))))
    dirs = 0.0
    for key in ("cx", "cy", "cz"):
        dirs = max(
            dirs,
            float(
                np.max(
                    np.abs(
                        plane[key][pa].astype(np.float64)
                        - sphere[key][pb].astype(np.float64)
                    )
                )
            ),
        )
    out.update(
        {
            "max_position_difference_m": worst_pos,
            "position_budget_ratio": pos_ratio,
            "max_time_difference_ns": float(np.max(dt)),
            "time_budget_ratio": t_ratio,
            "max_direction_difference": dirs,
            "zem_difference_m": float(
                np.max(
                    np.abs(
                        plane["zem"][pa].astype(np.float64)
                        - sphere["zem"][pb].astype(np.float64)
                    )
                )
            ),
        }
    )
    return out


def fixture_model():
    from obdeect.ray_tracing import OpticalModel

    return OpticalModel.load(str(FIXTURE))


def lst_model():
    from obdeect.ray_tracing import OpticalModel

    models = common.models_root()
    if not (models / "productions" / od.MODEL_VERSION).is_dir():
        return None
    return OpticalModel.from_simulation_models(
        models,
        "LSTN-01",
        od.MODEL_VERSION,
        configuration_seed=od.CONFIGURATION_SEED,
        zenith_deg=ZENITH_DEG,
    )


def time_convention(
    model, plane_grid: Path, sphere_grid: Path, plane_cam: dict, sphere_cam: dict
) -> dict:
    """Expected camera-time offset of plane rows against sphere rows, per photon, and the residual.

    Convention (obdeect ray tracer, stream_trace.hpp, plane rows): a plane row is moved along its line to
    the ground-frame input plane z = plane_z (vacuum, t + s |d| / c), and the traced path from that start
    point is multiplied by the group index n_g (trace.hpp). A sphere row starts at its entry point E with its
    own time, and the same path is multiplied by n_g. With P = the input-plane start of the plane row and
    D = (P - E) . d (signed along the ray), the camera times differ by (1 - n_g) D / c, from the straight-line
    geometry of the common ray. The residual is the measured difference minus this offset.
    """
    from obdeect.ray_tracing.stream import minimum_input_plane_height_m

    n_g = float(model.propagation_group_index)
    plane_z = float(minimum_input_plane_height_m(model, ZENITH_DEG, AZIMUTH_DEG))
    grid_plane = read_camera(plane_grid)
    grid_sphere = read_camera(sphere_grid)

    def index(cols):
        return {v: i for i, v in enumerate(ids_of(cols["wavelength"]).tolist())}

    ip, isph = index(grid_plane), index(grid_sphere)
    cp, cs = index(plane_cam), index(sphere_cam)
    ids = sorted(set(ip) & set(isph) & set(cp) & set(cs))
    if not ids:
        return {"count": 0}
    gp = np.array([ip[v] for v in ids])
    gs = np.array([isph[v] for v in ids])
    camp = np.array([cp[v] for v in ids])
    cams = np.array([cs[v] for v in ids])
    cz = grid_plane["cz"][gp].astype(np.float64)
    s = (
        plane_z / cz
    )  # plane row z = 0 in the ground frame; move to the input plane z = plane_z
    start = np.stack(
        [
            grid_plane["x"][gp].astype(np.float64)
            + s * grid_plane["cx"][gp].astype(np.float64),
            grid_plane["y"][gp].astype(np.float64)
            + s * grid_plane["cy"][gp].astype(np.float64),
            np.full(len(ids), plane_z),
        ],
        axis=1,
    )
    entry = np.stack(
        [grid_sphere[k][gs].astype(np.float64) for k in ("x", "y", "z")], axis=1
    )
    direction = np.stack(
        [grid_sphere[k][gs].astype(np.float64) for k in ("cx", "cy", "cz")], axis=1
    )
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)
    along = np.sum((start - entry) * direction, axis=1)  # D, signed along the ray [m]
    expected = (1.0 - n_g) * along / uc.C_M_NS
    t_plane = plane_cam["t"][camp].astype(np.float64)
    t_sphere = sphere_cam["t"][cams].astype(np.float64)
    measured = t_plane - t_sphere
    residual = measured - expected
    ratio = np.abs(residual) / budget(t_plane)
    return {
        "count": len(ids),
        "n_g": n_g,
        "input_plane_z_m": plane_z,
        "path_difference_max_m": float(np.max(np.abs(along))),
        "offset_max_ns": float(np.max(np.abs(expected))),
        "measured_max_ns": float(np.max(np.abs(measured))),
        "residual_max_ns": float(np.max(np.abs(residual))),
        "residual_budget_ratio": float(np.max(ratio)),
    }


def pair_part(root: Path, name: str, model, columns: dict, radius_m: float) -> dict:
    plane_grid, sphere_grid = propagate_pair(root, name, columns, radius_m)
    trace(
        model,
        plane_grid,
        root / f"{name}-plane-cam.pst",
        seed=5,
        pointing=(ZENITH_DEG, AZIMUTH_DEG),
    )
    trace(
        model,
        sphere_grid,
        root / f"{name}-sphere-cam.pst",
        seed=5,
        pointing=(ZENITH_DEG, AZIMUTH_DEG),
    )
    plane_cam = read_camera(root / f"{name}-plane-cam.pst")
    sphere_cam = read_camera(root / f"{name}-sphere-cam.pst")
    result = compare_camera(plane_cam, sphere_cam)
    result["time_convention"] = time_convention(
        model, plane_grid, sphere_grid, plane_cam, sphere_cam
    )
    result["photons"] = len(columns["x"])
    return result


def run_all(root: Path) -> tuple[dict, dict]:
    fixture = fixture_model()
    details: dict = {
        "obdeect_ray_tracing_version": version(),
        "zenith_deg": ZENITH_DEG,
        "azimuth_deg": AZIMUTH_DEG,
    }

    fx = pair_part(
        root,
        "fixture",
        fixture,
        parallel_beam(PHOTONS_FIXTURE, FIXTURE_DISK_M, SEED),
        FIXTURE_RADIUS_M,
    )
    details["fixture"] = fx

    # horizontal sphere rows on the tilted fixture telescope: traced, not rejected
    horiz = horizontal_rows(2000, SEED + 1)
    _hplane, hsphere = propagate_pair(root, "horizontal", horiz, FIXTURE_RADIUS_M)
    hrows = read_camera(hsphere)
    hsummary = trace(
        fixture,
        hsphere,
        root / "horizontal-cam.pst",
        seed=5,
        pointing=(ZENITH_DEG, AZIMUTH_DEG),
    )
    details["horizontal"] = {
        "sphere_rows": len(hrows["x"]),
        "trace_summary": hsummary,
    }

    # the fixture sphere is accepted: its trace above returned (the enclosure check passed)
    fixture_radius_accepted = True
    # enclosure: a 0.5 m sphere is rejected (message), the fixture sphere is accepted
    small_grid = _write(
        root / "small-sphere.pst", parallel_beam(200, FIXTURE_DISK_M, SEED + 2), 0.5
    )
    from obdeect.propagation import propagate
    from obdeect.ray_tracing.errors import ModelError

    small = root / "small-sphere-sphere.pst"
    propagate(
        small_grid,
        small,
        extinction=False,
        site_altitude_m=0.0,
        seed=1,
        ground_representation="sphere",
    )
    rejected, message = 0, ""
    try:
        trace(
            fixture,
            small,
            root / "small-cam.pst",
            seed=5,
            pointing=(ZENITH_DEG, AZIMUTH_DEG),
        )
    except ModelError as error:  # the enclosure message of the tracer
        rejected, message = 1, str(error)
    details["enclosure"] = {
        "small_radius_m": 0.5,
        "rejected": rejected,
        "message": message,
    }

    lst = lst_model()
    lst_part = None
    if lst is not None:
        lst_part = pair_part(
            root,
            "lst",
            lst,
            parallel_beam(PHOTONS_LST, LST_DISK_M, SEED + 3),
            LST_RADIUS_M,
        )
    details["lst"] = (
        lst_part
        if lst_part is not None
        else {"compared": False, "reason": "no simulation models checkout"}
    )

    metrics = {
        "fixture_detected_set_mismatch": fx["set_mismatch"],
        "fixture_detected_plane": fx["detected_plane"],
        "fixture_detected_sphere": fx["detected_sphere"],
        "fixture_position_budget_ratio": fx.get("position_budget_ratio", math.inf),
        "fixture_time_residual_budget_ratio": fx["time_convention"].get(
            "residual_budget_ratio", math.inf
        ),
        "fixture_time_offset_max_ns": fx["time_convention"].get(
            "offset_max_ns", math.inf
        ),
        "fixture_time_residual_max_ns": fx["time_convention"].get(
            "residual_max_ns", math.inf
        ),
        "fixture_time_raw_budget_ratio": fx.get("time_budget_ratio", math.inf),
        "fixture_max_position_difference_m": fx.get(
            "max_position_difference_m", math.inf
        ),
        "fixture_max_time_difference_ns": fx.get("max_time_difference_ns", math.inf),
        "fixture_max_direction_difference": fx.get(
            "max_direction_difference", math.inf
        ),
        "lst_compared": int(lst_part is not None),
        "lst_detected_set_mismatch": (lst_part or {}).get("set_mismatch", 0),
        "lst_position_budget_ratio": (lst_part or {}).get("position_budget_ratio", 0.0),
        "lst_time_residual_budget_ratio": (
            (lst_part or {}).get("time_convention") or {}
        ).get("residual_budget_ratio", 0.0),
        "lst_time_offset_max_ns": ((lst_part or {}).get("time_convention") or {}).get(
            "offset_max_ns", 0.0
        ),
        "lst_time_residual_max_ns": ((lst_part or {}).get("time_convention") or {}).get(
            "residual_max_ns", 0.0
        ),
        "lst_time_raw_budget_ratio": (lst_part or {}).get("time_budget_ratio", 0.0),
        "lst_max_position_difference_m": (lst_part or {}).get(
            "max_position_difference_m", 0.0
        ),
        "lst_max_time_difference_ns": (lst_part or {}).get(
            "max_time_difference_ns", 0.0
        ),
        "horizontal_rows_traced": int(hsummary["photons_in"] == len(hrows["x"])),
        "horizontal_sphere_rows": len(hrows["x"]),
        "enclosure_small_rejected": rejected,
        "enclosure_message_names_element": int(
            "must enclose the optical model" in message
        ),
        "enclosure_fixture_accepted": int(
            fx["detected_sphere"] > 0 and fixture_radius_accepted
        ),
    }
    return metrics, details


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args()
    metrics, details = run_all(Path.cwd())
    common.write_results(metrics, details, Path.cwd())


if __name__ == "__main__":
    main()
