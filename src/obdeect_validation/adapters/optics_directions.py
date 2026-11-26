"""Telescope optics with light from below: uniform beams from below, a downward control, determinism.

Three parts on the optical models LSTN-01, MSTN-01 and SSTS-01 (simulation models 7.0.0, compiled by
obdeect) and on the segmented test model of the obdeect fixtures (``fixtures/optical-models``; it
needs no simulation models).

1. Blocked fraction. A uniform parallel beam from below, along the telescope axis and at 20 and
   45 degrees, on a disc larger than the primary. A photon that meets the back of a primary facet or
   segment is blocked there (status ``blocked``, terminal surface a primary id). The reference is the
   area of the primary projected along the beam onto the start plane, divided by the disc area,
   computed here from the mirror geometry of the optical model (numpy only):

   - facets: spherical caps (sphere of radius R about C + R n) over their tangent-plane apertures
     (hexagons, flat-to-flat = diameter, in the orientation of the model, or circles). A beam meets
     a cap when the first hit of its sphere lies in the aperture. Shadows of facets can overlap
     (the dishes at large inclinations), so the union of the shadows is taken: scrambled Sobol
     points on the start disc, each tested against the facets whose shadow contains it;
   - segments of the aspheric primary: z = vertex + sum_i c_i (r / scale)^(2 i) over the segment
     apertures; a segment is the hexagon (flat-to-flat = diameter) in the tangent plane at its
     centre, a point of the asphere mapped to that plane by its offset from the centre (radial part
     scaled by sqrt(1 + slope^2)), clipped to the asphere's inner and outer radii; area element
     |(-z_x, -z_y, 1) . d| dx dy.

   The pulls are binomial. Obscurers must lie above the primary (checked), so they never shade it.
2. Downward control. The same beam from above, traced alone and mixed with the upward beam: the
   downward rows are bit-identical in every output.
3. Determinism. The mixed beam traced with one thread, four threads, a thread pool, in chunks of
   several sizes and in a permuted row order: identical per photon (matched by photon id).

Only the geometry is read from the models; the photons and their first hits are traced by obdeect.
Fixture provenance: ``fixtures/optical-models/segmented.json`` is a copy of
``src/ray_tracing/python/tests/data/optical-models/segmented.json`` of the obdeect checkout at commit
f4ea1bb (7 circular facets, one plate obscurer; the model's own sha256 is inside the file).
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

try:  # the traced runs need obdeect; the reference helpers do not
    from obdeect.ray_tracing import STATUS_NAMES, OpticalModel, ThreadPool
except ImportError:  # pragma: no cover - exercised by the reference-only environment
    STATUS_NAMES = OpticalModel = ThreadPool = None

from obdeect_validation.adapters import common, plotting

MODEL_VERSION = "7.0.0"
MODELS = (("LSTN-01", "LST"), ("MSTN-01", "MST"), ("SSTS-01", "SST"))
CONFIGURATION_SEED = 11  # mirror errors, as the simtools optics cases
ZENITH_DEG = 20.0
FIXTURE = "fixtures/optical-models/segmented.json"
INCLINATIONS_DEG = (0.0, 20.0, 45.0)
SEED = 20261124
WAVELENGTH_NM = 400.0
PHOTONS = 500_000
DETERMINISM_PHOTONS = 100_000
START_MARGIN_M = 1.0  # start plane this far below the lowest point of the primary
DISK_SAFETY = 1.02  # the start disc is this much larger than the shadow of the primary
FACET_POINTS_LOG2 = (
    22  # scrambled Sobol points on the start disc (union of facet shadows)
)
SEGMENT_GRID = 400  # cells across the box of a segment
PULL_LIMIT = 4.0
SQRT3_2 = math.sqrt(3.0) / 2.0
STATUS_COLOURS = {
    "mirror back (primary)": plotting.NATIVE,
    "obscurer": plotting.CYCLE[5],
    "non_sequential": plotting.CYCLE[2],
    "direct_focal": plotting.CYCLE[3],
    "missed_primary": plotting.LIMIT,
}


FRACTION_COLOURS = {
    "LST": plotting.CASE_COLORS["lst"],
    "MST": plotting.CASE_COLORS["mst"],
    "SST": plotting.CASE_COLORS["sst"],
    "fixture": plotting.CYCLE[6],
}


class AdapterError(ValueError):
    """Geometry or configuration the reference cannot handle."""


# ------------------------------------------------------------------ mirror geometry (reference)


def unit(vector) -> np.ndarray:
    value = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(value))
    if not math.isfinite(norm) or norm == 0.0:
        raise AdapterError(f"zero or undefined vector {vector!r}")
    return value / norm


def hexagon(u, v, apothem: float):
    """Regular hexagon with flats at distance `apothem` from its centre (flat-to-flat 2 apothem)."""
    return (
        (np.abs(u) <= apothem)
        & (np.abs(0.5 * u + SQRT3_2 * v) <= apothem)
        & (np.abs(0.5 * u - SQRT3_2 * v) <= apothem)
    )


def facet_frame(facet: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unit normal, tangent and binormal (n x t) of a facet; the aperture lies in t, binormal."""
    normal = unit(facet["normal"])
    tangent = unit(
        np.asarray(facet["tangent"], float) - np.dot(facet["tangent"], normal) * normal
    )
    return normal, tangent, np.cross(normal, tangent)


def facet_aperture(facet: dict[str, Any], u, v):
    """Points (u, v) of the tangent plane inside the aperture of a facet."""
    apothem = 0.5 * float(facet["diameter_m"])
    shape = facet["shape"]
    if shape == "hexagon_flat_x":
        return hexagon(u, v, apothem)
    if (
        shape == "hexagon_flat_y"
    ):  # flats normal to the binormal: the same hexagon, axes swapped
        return hexagon(v, u, apothem)
    if shape == "circle":
        return u * u + v * v <= apothem * apothem
    raise AdapterError(f"unsupported facet shape {shape!r}")


def facet_hits(
    facet: dict[str, Any], x0, y0, z0: float, direction: np.ndarray
) -> np.ndarray:
    """Which start points (x0, y0, z0) meet this facet from below.

    The facet is the cap of its sphere (centre C + R n) inside its aperture; a ray from a start
    point below the facet meets the sphere first on the cap side. The hit must lie in the aperture.
    """
    normal, tangent, binormal = facet_frame(facet)
    centre = np.asarray(facet["centre_m"], float)
    radius = float(facet.get("radius_of_curvature_m") or 0.0)
    origin = np.stack([x0, y0, np.full_like(x0, z0, dtype=float)], axis=-1)
    if radius > 0.0:
        w = origin - (centre + radius * normal)
        along = w @ direction
        discriminant = along * along - ((w * w).sum(axis=-1) - radius * radius)
        reached = discriminant >= 0.0
        step = -along - np.sqrt(np.where(reached, discriminant, 0.0))
        reached &= step > 0.0
    else:
        denominator = float(np.dot(direction, normal))
        step = ((centre - origin) @ normal) / denominator
        reached = step > 0.0
    point = origin + step[..., None] * direction
    relative = point - centre
    return reached & facet_aperture(facet, relative @ tangent, relative @ binormal)


def facet_shadow_disc(
    facet: dict[str, Any], z0: float, direction: np.ndarray
) -> tuple[float, float, float]:
    """Disc (centre x, centre y, radius) that contains the shadow of a facet on the start plane.

    The shadow of a point P of the facet is P minus (P_z - z0) d / d_z. Relative to the projected
    centre C', its horizontal offset is at most |P - C|_xy + |P_z - C_z| tan(theta), and for a cap of
    circumradius rho: |P - C| <= rho (plus the sag), |P_z - C_z| <= rho sqrt(1 - n_z^2) (plus sag).
    """
    centre = np.asarray(facet["centre_m"], float)
    normal = unit(facet["normal"])
    reach = math.hypot(direction[0], direction[1]) / abs(direction[2])
    cx = centre[0] - (centre[2] - z0) * direction[0] / direction[2]
    cy = centre[1] - (centre[2] - z0) * direction[1] / direction[2]
    circum = (
        float(facet["diameter_m"]) / math.sqrt(3.0)
        if facet["shape"].startswith("hexagon")
        else 0.5 * float(facet["diameter_m"])
    )
    radius_of_curvature = float(facet.get("radius_of_curvature_m") or 0.0)
    sag = (
        circum * circum / (2.0 * radius_of_curvature)
        if radius_of_curvature > 0
        else 0.0
    )
    tilt = math.sqrt(max(0.0, 1.0 - normal[2] ** 2))
    delta = circum + (circum * tilt + sag) * reach + 0.02
    return cx, cy, delta


def facet_union_fraction(
    facets: list[dict[str, Any]],
    z0: float,
    inclination_deg: float,
    radius: float,
    log2_points: int = FACET_POINTS_LOG2,
) -> tuple[float, float]:
    """Fraction of the start disc whose beam meets at least one facet from below.

    Scrambled Sobol points (2^log2_points) fill the start disc; each facet is tested on the points
    of its shadow only (a disc around the projected facet centre). Overlapping shadows count once.
    Returns (union fraction, sum of the single shadows) with the same points.
    """
    from scipy.stats import qmc

    t = math.radians(inclination_deg)
    direction = np.array([math.sin(t), 0.0, math.cos(t)])
    unit_square = qmc.Sobol(d=2, scramble=True, seed=SEED).random_base2(log2_points)
    x = (2.0 * unit_square[:, 0] - 1.0) * radius
    y = (2.0 * unit_square[:, 1] - 1.0) * radius
    inside = x * x + y * y <= radius * radius
    x, y = x[inside], y[inside]
    order = np.argsort(x)
    xs = x[order]
    covered = np.zeros(len(x), dtype=bool)
    single = 0
    for facet in facets:
        cx, cy, delta = facet_shadow_disc(facet, z0, direction)
        if math.hypot(cx, cy) + delta > radius:
            raise AdapterError(
                f"the shadow of facet {facet['id']} reaches the edge of the disc"
            )
        low, high = np.searchsorted(xs, [cx - delta, cx + delta])
        candidates = order[low:high]
        candidates = candidates[np.abs(y[candidates] - cy) <= delta]
        hit = facet_hits(facet, x[candidates], y[candidates], z0, direction)
        covered[candidates[hit]] = True
        single += int(np.count_nonzero(hit))
    return float(np.count_nonzero(covered)) / len(x), single / len(x)


def asphere_height(asphere: dict[str, Any], x, y):
    """Sag of an even polynomial asphere: vertex + sum_i c_i (r / scale)^(2 i)."""
    scale = float(asphere["radial_scale_m"])
    q = (np.asarray(x, float) ** 2 + np.asarray(y, float) ** 2) / scale**2
    value = np.zeros_like(q)
    for coefficient in reversed(asphere["coefficients_m"]):
        value = value * q + coefficient
    return float(asphere["vertex_z_m"]) + value


def asphere_gradient(asphere: dict[str, Any], x, y):
    """(dz/dx, dz/dy) of the asphere: 2 x / scale^2 * sum_{i>=1} i c_i q^(i-1), q = r^2 / scale^2."""
    scale = float(asphere["radial_scale_m"])
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    q = (x * x + y * y) / scale**2
    slope = np.zeros_like(q)
    coefficients = asphere["coefficients_m"]
    for i in range(len(coefficients) - 1, 0, -1):
        slope = slope * q + i * coefficients[i]
    return 2.0 * x / scale**2 * slope, 2.0 * y / scale**2 * slope


def segment_aperture(asphere: dict[str, Any], segment: dict[str, Any], x, y):
    """Mask of the points (x, y) of the asphere that belong to a segment (the model's rule).

    The segment is the hexagon (flat-to-flat = diameter) in the tangent plane at the centre of the
    segment, with the rotation of the segment. A point is mapped to that plane from its offset to
    the centre: the radial offset, including the sag, is measured in the tangent plane (scaled by
    sqrt(1 + slope^2) for the radial part); the azimuthal offset is not scaled.
    """
    if segment["shape"] != "hexagon":
        raise AdapterError(f"unsupported segment shape {segment['shape']!r}")
    cx, cy = (float(c) for c in segment["centre_xy_m"])
    radius = math.hypot(cx, cy)
    er = np.array([cx, cy]) / radius
    ephi = np.array([-er[1], er[0]])
    scale = float(asphere["radial_scale_m"])
    q_centre = (radius / scale) ** 2
    slope_sum = sum(
        2 * i * c * q_centre ** (i - 1)
        for i, c in enumerate(asphere["coefficients_m"])
        if i
    )
    slope = radius / scale**2 * slope_sum
    cosine = 1.0 / math.sqrt(1.0 + slope * slope)
    sine = slope * cosine
    dx, dy = x - cx, y - cy
    dz = asphere_height(asphere, x, y) - float(asphere_height(asphere, cx, cy))
    radial = (dx * er[0] + dy * er[1]) * cosine + dz * sine
    azimuthal = -dx * er[1] + dy * er[0]
    px = radial * er[0] + azimuthal * ephi[0]
    py = radial * er[1] + azimuthal * ephi[1]
    angle = math.radians(float(segment.get("rotation_deg", 0.0)))
    u = px * math.cos(angle) + py * math.sin(angle)
    v = -px * math.sin(angle) + py * math.cos(angle)
    inside = hexagon(u, v, 0.5 * float(segment["diameter_m"]))
    r = np.hypot(x, y)
    return inside & (r >= asphere["inner_radius_m"]) & (r <= asphere["outer_radius_m"])


def segment_projected_area(
    asphere: dict[str, Any], segment: dict[str, Any], direction
) -> float:
    """Area of one segment of the asphere projected along `direction` onto z = const [m^2]."""
    cx, cy = (float(c) for c in segment["centre_xy_m"])
    half = 1.25 * 0.5 * float(segment["diameter_m"])
    step = 2.0 * half / SEGMENT_GRID
    axis = -half + step * (np.arange(SEGMENT_GRID) + 0.5)
    xx, yy = (g.ravel() for g in np.meshgrid(cx + axis, cy + axis, indexing="ij"))
    inside = segment_aperture(asphere, segment, xx, yy).reshape(
        SEGMENT_GRID, SEGMENT_GRID
    )
    if inside[0].any() or inside[-1].any() or inside[:, 0].any() or inside[:, -1].any():
        raise AdapterError("a segment reaches the edge of its grid box")
    inside = inside.ravel()
    fx, fy = asphere_gradient(asphere, xx[inside], yy[inside])
    along = -fx * direction[0] - fy * direction[1] + direction[2]
    return float(np.sum(np.abs(along))) * step * step / abs(direction[2])


def primary_summary(document: dict[str, Any]) -> dict[str, Any]:
    """Kind, ids, items and height bounds of the primary."""
    primary = document["optics"]["primary"]
    if "facets" in primary:
        facets = primary["facets"]
        circumradius = [
            float(f["diameter_m"]) / math.sqrt(3.0)
            if f["shape"].startswith("hexagon")
            else 0.5 * float(f["diameter_m"])
            for f in facets
        ]
        centres_z = [float(f["centre_m"][2]) for f in facets]
        return {
            "kind": "facets",
            "ids": {int(f["id"]) for f in facets},
            "items": facets,
            "circumradius_m": circumradius,
            "lowest_z_m": min(centres_z),
            "highest_z_m": max(centres_z) + 0.05,  # margin for the sag of the caps
        }
    if "asphere" in primary:
        asphere = primary["asphere"]
        radii = np.linspace(asphere["inner_radius_m"], asphere["outer_radius_m"], 2001)
        heights = asphere_height(asphere, radii, np.zeros_like(radii))
        return {
            "kind": "asphere",
            "ids": {int(s["id"]) for s in primary["segments"]},
            "items": primary["segments"],
            "asphere": asphere,
            "outer_radius_m": float(asphere["outer_radius_m"]),
            "lowest_z_m": float(heights.min()),
            "highest_z_m": float(heights.max()),
        }
    raise AdapterError("primary has neither facets nor an asphere")


def obscurer_z_range(obscurer: dict[str, Any]) -> tuple[float, float]:
    if obscurer["shape"] == "plate":
        z = float(obscurer["centre_m"][2])
        return z, z
    if obscurer["shape"] in ("frustum", "cylinder"):
        ends = (float(obscurer["start_m"][2]), float(obscurer["end_m"][2]))
        return min(ends), max(ends)
    raise AdapterError(
        f"cannot bound obscurer {obscurer['name']} ({obscurer['shape']})"
    )


@dataclass(frozen=True)
class Geometry:
    """Reference geometry of one optical model, from its compiled document."""

    primary: dict[str, Any]
    obscurer_ids: frozenset[int]
    start_z_m: float  # start plane of upward photons, below every primary point
    above_z_m: float  # start plane of downward photons, above every optic

    @classmethod
    def of(cls, document: dict[str, Any]) -> Geometry:
        optics = document["optics"]
        primary = primary_summary(document)
        for obscurer in optics["obscurers"]:
            lowest, _ = obscurer_z_range(obscurer)
            if lowest <= primary["highest_z_m"]:
                raise AdapterError(
                    f"obscurer {obscurer['name']} reaches the height of the primary: upward "
                    "photons would not meet the primary first"
                )
        focal = optics["focal_surface"]
        focal_z = (
            float(focal["plane_z_m"])
            if "plane_z_m" in focal
            else float(focal["asphere"]["vertex_z_m"])
        )
        tops = [focal_z, *(obscurer_z_range(o)[1] for o in optics["obscurers"])]
        return cls(
            primary=primary,
            obscurer_ids=frozenset(int(o["id"]) for o in optics["obscurers"]),
            start_z_m=primary["lowest_z_m"] - START_MARGIN_M,
            above_z_m=max(tops) + START_MARGIN_M,
        )

    def shadow_reach(self, inclination_deg: float) -> float:
        """Largest distance from the axis of the shadow of the primary on the start plane [m]."""
        if self.primary["kind"] == "facets":
            t = math.radians(inclination_deg)
            direction = np.array([math.sin(t), 0.0, math.cos(t)])
            reach = 0.0
            for facet in self.primary["items"]:
                cx, cy, delta = facet_shadow_disc(facet, self.start_z_m, direction)
                reach = max(reach, math.hypot(cx, cy) + delta)
            return reach
        tan = math.tan(math.radians(inclination_deg))
        height = self.primary["highest_z_m"] - self.start_z_m
        return self.primary["outer_radius_m"] + height * abs(tan)

    def disk_radius_m(self, inclination_deg: float) -> float:
        """Radius of the start disc, centred on the axis, that contains the whole shadow."""
        return DISK_SAFETY * self.shadow_reach(inclination_deg) + 0.05

    def reference_fraction(self, inclination_deg: float) -> tuple[float, float, float]:
        """(fraction of the start disc blocked by the primary, disc radius [m], single-shadow sum).

        Facets: the union of their shadows (grid). Segments of the asphere do not overlap in
        projection (the asphere is convex), so their shadows add up (exact quadrature).
        """
        radius = self.disk_radius_m(inclination_deg)
        if self.primary["kind"] == "facets":
            fraction, sum_fraction = facet_union_fraction(
                self.primary["items"], self.start_z_m, inclination_deg, radius
            )
        else:
            t = math.radians(inclination_deg)
            direction = np.array([math.sin(t), 0.0, math.cos(t)])
            area = sum(
                segment_projected_area(self.primary["asphere"], s, direction)
                for s in self.primary["items"]
            )
            fraction = sum_fraction = area / (math.pi * radius * radius)
        if not 0.0 < fraction < 1.0:
            raise AdapterError(f"reference fraction {fraction:.4f} outside (0, 1)")
        return fraction, radius, sum_fraction


# ------------------------------------------------------------------ beams and tracing


def direction_of(inclination_deg: float, downward: bool = False) -> np.ndarray:
    """Unit direction of travel: upward (+z) or downward, tilted by the inclination towards +x."""
    t = math.radians(inclination_deg)
    return np.array([math.sin(t), 0.0, -math.cos(t) if downward else math.cos(t)])


def disk_photons(
    rng, count: int, radius: float, z: float, direction, first_id: int
) -> dict:
    """Photons uniform in area on the disc of `radius` at height z, all travelling along direction."""
    r = radius * np.sqrt(rng.uniform(size=count))
    phi = rng.uniform(0.0, 2.0 * np.pi, size=count)
    position = np.stack([r * np.cos(phi), r * np.sin(phi), np.full(count, z)], axis=1)
    return {
        "position_m": position,
        "direction": np.tile(np.asarray(direction, float), (count, 1)),
        "wavelength_nm": np.full(count, WAVELENGTH_NM),
        "time_ns": np.zeros(count),
        "weight": np.ones(count),
        "photon_id": np.arange(first_id, first_id + count, dtype=np.uint64),
    }


def concatenate(parts: list[dict]) -> dict:
    return {key: np.concatenate([p[key] for p in parts]) for key in parts[0]}


def take(photons: dict, index) -> dict:
    return {key: value[index] for key, value in photons.items()}


def differing_rows(a: dict, b: dict) -> np.ndarray:
    """Rows where any output differs; the comparison is on the bits, so NaN equals NaN."""
    rows = len(next(iter(a.values())))
    differs = np.zeros(rows, dtype=bool)
    for key, left in a.items():
        left = np.ascontiguousarray(left)
        right = np.ascontiguousarray(b[key])
        if len(right) != rows:
            raise AdapterError(f"output {key} has {len(right)} rows, expected {rows}")
        differs |= np.any(
            left.view(np.uint8).reshape(rows, -1)
            != right.view(np.uint8).reshape(rows, -1),
            axis=1,
        )
    return differs


def status_code(name: str) -> int:
    return STATUS_NAMES.index(name)


def status_counts(status: np.ndarray) -> dict[str, int]:
    return {
        name: int(np.count_nonzero(status == i)) for i, name in enumerate(STATUS_NAMES)
    }


def binomial_pull(count: int, photons: int, probability: float) -> float:
    return (count - photons * probability) / math.sqrt(
        photons * probability * (1.0 - probability)
    )


def trace(model: OpticalModel, photons: dict, **options) -> dict:
    return model.trace(photons, seed=SEED, **options)


def beam_part(
    model, geometry: Geometry, inclination: float, count: int, seed: int
) -> tuple:
    """Upward beam at one inclination: blocked fraction against the reference.

    Returns the record (JSON) and the start positions and statuses for the figure.
    """
    fraction, radius, shadow_sum = geometry.reference_fraction(inclination)
    photons = disk_photons(
        np.random.default_rng(seed),
        count,
        radius,
        geometry.start_z_m,
        direction_of(inclination),
        first_id=0,
    )
    out = trace(model, photons)
    counts = status_counts(out["status"])
    primary = (out["status"] == status_code("blocked")) & np.isin(
        out["terminal_surface"], sorted(geometry.primary["ids"])
    )
    obscurer = (out["status"] == status_code("blocked")) & np.isin(
        out["terminal_surface"], sorted(geometry.obscurer_ids)
    )
    unknown = (out["status"] == status_code("blocked")) & ~primary & ~obscurer
    record = {
        "inclination_deg": inclination,
        "disk_radius_m": radius,
        "start_z_m": geometry.start_z_m,
        "photons": count,
        "status_counts": counts,
        "status_sum_error": abs(sum(counts.values()) - count),
        "detected": counts["detected"],
        "blocked_primary": int(np.count_nonzero(primary)),
        "blocked_obscurer": int(np.count_nonzero(obscurer)),
        "blocked_unknown_surface": int(np.count_nonzero(unknown)),
        "reference_fraction": fraction,
        "reference_area_on_plane_m2": fraction * math.pi * radius * radius,
        "shadow_sum_fraction": shadow_sum,
        "overlap_fraction": shadow_sum - fraction,
        "obdeect_fraction": int(np.count_nonzero(primary)) / count,
        "obdeect_fraction_error": math.sqrt(fraction * (1.0 - fraction) / count),
        "mirror_back_pull": binomial_pull(
            int(np.count_nonzero(primary)), count, fraction
        ),
    }
    figure = {
        "positions": photons["position_m"],
        "status": out["status"],
        "primary": primary,
        "inclination_deg": inclination,
        "start_z_m": geometry.start_z_m,
    }
    return record, figure


def downward_part(
    model, geometry: Geometry, inclination: float, count: int, seed: int
) -> dict:
    """Downward control: the downward rows of a mixed trace equal a downward-only trace."""
    rng = np.random.default_rng(seed)
    radius = geometry.disk_radius_m(inclination)
    up = disk_photons(
        rng, count, radius, geometry.start_z_m, direction_of(inclination), 0
    )
    down = disk_photons(
        rng, count, radius, geometry.above_z_m, direction_of(inclination, True), count
    )
    alone = trace(model, down)
    mixed = trace(model, concatenate([up, down]))
    mixed_down = {key: value[count:] for key, value in mixed.items()}
    counts = status_counts(alone["status"])
    return {
        "inclination_deg": inclination,
        "downward_photons": count,
        "downward_mismatch_rows": int(
            np.count_nonzero(differing_rows(mixed_down, alone))
        ),
        "status_counts": counts,
        "status_sum_error": abs(sum(counts.values()) - count),
    }


def determinism_part(model, geometry: Geometry, count: int, seed: int) -> dict:
    """Mixed beam at 20 degrees: threads, a pool, chunks of several sizes, a permuted order."""
    rng = np.random.default_rng(seed)
    radius = geometry.disk_radius_m(20.0)
    half = count // 2
    up = disk_photons(rng, half, radius, geometry.start_z_m, direction_of(20.0), 0)
    down = disk_photons(
        rng, count - half, radius, geometry.above_z_m, direction_of(20.0, True), half
    )
    photons = concatenate([up, down])
    total = len(photons["photon_id"])
    base = trace(model, photons, threads=1)
    variants: dict[str, dict] = {
        "threads_4": trace(model, photons, threads=4),
        "pool_3": trace(model, photons, pool=ThreadPool(3)),
    }
    bounds = [0, 997, 5000, 5001, 17000, total]
    chunks = [
        trace(model, take(photons, slice(lo, hi)), threads=2)
        for lo, hi in itertools.pairwise(bounds)
        if hi > lo
    ]
    variants["chunks"] = {key: np.concatenate([c[key] for c in chunks]) for key in base}
    permutation = rng.permutation(total)
    permuted = trace(model, take(photons, permutation))
    inverse = np.argsort(permutation)
    variants["row_order"] = {key: value[inverse] for key, value in permuted.items()}
    mismatches = {
        name: int(np.count_nonzero(differing_rows(base, out)))
        for name, out in variants.items()
    }
    return {
        "photons": total,
        "mismatch_rows": mismatches,
        "determinism_mismatch_count": sum(mismatches.values()),
        "status_counts": status_counts(base["status"]),
    }


def model_part(
    name: str, label: str, model, photons: int, determinism_photons: int
) -> dict:
    geometry = Geometry.of(model.document)
    beams, figures = [], {}
    for index, inclination in enumerate(INCLINATIONS_DEG):
        record, figure = beam_part(model, geometry, inclination, photons, SEED + index)
        beams.append(record)
        figures[inclination] = figure
    downward = [
        downward_part(model, geometry, inclination, photons, SEED + 100 + index)
        for index, inclination in enumerate(INCLINATIONS_DEG)
    ]
    section = {
        "name": name,
        "label": label,
        "sha256": model.sha256,
        "primary_kind": geometry.primary["kind"],
        "primary_surfaces": len(geometry.primary["ids"]),
        "obscurer_ids": sorted(geometry.obscurer_ids),
        "beams": beams,
        "downward": downward,
        "determinism": determinism_part(
            model, geometry, determinism_photons, SEED + 200
        ),
    }
    return {"section": section, "figures": figures, "geometry": geometry}


def run_all(
    models_root: Path,
    fixture: Path,
    photons: int = PHOTONS,
    determinism_photons: int = DETERMINISM_PHOTONS,
) -> tuple[dict, dict, dict]:
    """All parts for the three optical models and the fixture: metrics, details, figure arrays."""
    parts = []
    for name, label in MODELS:
        model = OpticalModel.from_simulation_models(
            models_root,
            name,
            MODEL_VERSION,
            configuration_seed=CONFIGURATION_SEED,
            zenith_deg=ZENITH_DEG,
        )
        parts.append(model_part(name, label, model, photons, determinism_photons))
    parts.append(
        model_part(
            "segmented-fixture",
            "fixture",
            OpticalModel.load(fixture),
            photons,
            determinism_photons,
        )
    )
    sections = [p["section"] for p in parts]
    metrics = summarise(sections)
    details = {
        "photons": photons,
        "determinism_photons": determinism_photons,
        "seed": SEED,
    }
    details["models"] = sections
    lst = next(p for p in parts if p["section"]["label"] == "LST")
    arrays = {
        "geometry": lst["geometry"],
        "label": "LST",
        **lst["figures"][INCLINATIONS_DEG[-1]],
    }
    return metrics, details, arrays


def summarise(sections: list[dict]) -> dict:
    real = [s for s in sections if s["label"] != "fixture"]
    fixture = [s for s in sections if s["label"] == "fixture"]
    beams = [b for s in sections for b in s["beams"]]
    real_pulls = [b["mirror_back_pull"] for s in real for b in s["beams"]]
    fixture_pulls = [b["mirror_back_pull"] for s in fixture for b in s["beams"]]
    return {
        "detected_count": sum(b["detected"] for b in beams),
        "status_sum_error": sum(b["status_sum_error"] for b in beams)
        + sum(d["status_sum_error"] for s in sections for d in s["downward"])
        + sum(
            abs(
                sum(s["determinism"]["status_counts"].values())
                - s["determinism"]["photons"]
            )
            for s in sections
        ),
        "mirror_back_cases": len(real_pulls),
        "mirror_back_pull_max_abs": max((abs(p) for p in real_pulls), default=0.0),
        "fixture_mirror_back_cases": len(fixture_pulls),
        "fixture_mirror_back_pull_max_abs": max(
            (abs(p) for p in fixture_pulls), default=0.0
        ),
        "blocked_unknown_surface_count": sum(
            b["blocked_unknown_surface"] for b in beams
        ),
        "downward_mismatch_count": sum(
            d["downward_mismatch_rows"] for s in sections for d in s["downward"]
        ),
        "determinism_mismatch_count": sum(
            s["determinism"]["determinism_mismatch_count"] for s in sections
        ),
        "obscurer_blocked_fraction_max": max(
            b["blocked_obscurer"] / b["photons"] for b in beams
        ),
    }


# ------------------------------------------------------------------ plots


def plot_start_map(path: Path, arrays: dict) -> None:
    """Start positions of the upward beam on the start plane, coloured by the first hit."""
    geometry: Geometry = arrays["geometry"]
    positions, status, primary = (
        arrays["positions"],
        arrays["status"],
        arrays["primary"],
    )
    obscurer = (status == status_code("blocked")) & ~primary
    masks = {
        "mirror back (primary)": primary,
        "obscurer": obscurer,
        "non_sequential": status == STATUS_NAMES.index("non_sequential"),
        "direct_focal": status == STATUS_NAMES.index("direct_focal"),
        "missed_primary": status == STATUS_NAMES.index("missed_primary"),
    }
    fig, axis = plotting.figure(1, 1, (5.8, 5.4))
    keep = np.arange(0, len(positions), max(1, len(positions) // 60_000))
    for label, mask in masks.items():
        selected = keep[mask[keep]]
        axis.plot(
            positions[selected, 0],
            positions[selected, 1],
            ".",
            markersize=1.2,
            color=STATUS_COLOURS[label],
            label=f"{label} ({int(mask.sum()) / len(mask):.1%})",
        )
    tan = math.tan(math.radians(arrays["inclination_deg"]))
    rim = geometry.primary["highest_z_m"] - arrays["start_z_m"]
    if geometry.primary["kind"] == "asphere":
        outer = geometry.primary["outer_radius_m"]
    else:
        outer = max(
            math.hypot(f["centre_m"][0], f["centre_m"][1])
            + 0.5 * float(f["diameter_m"])
            for f in geometry.primary["items"]
        )
    theta = np.linspace(0.0, 2.0 * np.pi, 400)
    axis.plot(
        outer * np.cos(theta) - rim * tan,
        outer * np.sin(theta),
        "--",
        color=plotting.ANALYTIC,
        linewidth=0.9,
        label="rim of the primary, projected",
    )
    axis.set_aspect("equal")
    axis.set(
        xlabel="x on the start plane [m]",
        ylabel="y on the start plane [m]",
        title=f"{arrays['label']} at {arrays['inclination_deg']:.0f}° from the axis",
    )
    axis.legend(loc="upper right", fontsize="x-small", markerscale=6)
    plotting.save(fig, path)


def plot_fractions(path: Path, sections: list[dict]) -> None:
    """Blocked fraction per inclination (obdeect with binomial errors, reference dashed); pulls."""
    fig, axes = plotting.figure(2, 1, (6.2, 6.0))
    top, bottom = axes[0, 0], axes[1, 0]
    for index, section in enumerate(sections):
        colour = FRACTION_COLOURS.get(
            section["label"], plotting.CYCLE[index % len(plotting.CYCLE)]
        )
        beams = section["beams"]
        x = [b["inclination_deg"] for b in beams]
        top.errorbar(
            x,
            [b["obdeect_fraction"] for b in beams],
            yerr=[b["obdeect_fraction_error"] for b in beams],
            fmt="o",
            color=colour,
            markersize=4,
            capsize=2,
            label=f"obdeect {section['label']}",
        )
        top.plot(
            x,
            [b["reference_fraction"] for b in beams],
            "--",
            color=plotting.ANALYTIC,
            linewidth=0.9,
            label="reference (projected primary)" if index == 0 else None,
        )
        bottom.plot(
            x,
            [b["mirror_back_pull"] for b in beams],
            "o",
            color=colour,
            markersize=4,
            label=section["label"],
        )
    top.set(ylabel="blocked by the primary back (fraction)")
    top.legend(loc="best", fontsize="small")
    bottom.axhspan(
        -PULL_LIMIT, PULL_LIMIT, color=plotting.LIMIT, alpha=0.12, linewidth=0
    )
    bottom.axhline(0.0, color=plotting.ANALYTIC, linewidth=0.8)
    bottom.set(
        xlabel="inclination from the telescope axis [deg]", ylabel="binomial pull"
    )
    bottom.legend(loc="best", fontsize="small")
    plotting.save(fig, path)


# ------------------------------------------------------------------ main


def main() -> None:
    common.results_main(_analyse, __doc__.splitlines()[0])


def _analyse(root: Path):
    models = common.models_root()
    if not (models / "productions" / MODEL_VERSION).is_dir():
        raise AdapterError(f"no simulation models {MODEL_VERSION} under {models}")
    metrics, details, arrays = run_all(models, common.VALIDATION / FIXTURE)
    plot_start_map(root / "upward-start-map.png", arrays)
    plot_fractions(root / "upward-blocked-fraction.png", details["models"])
    return metrics, details


if __name__ == "__main__":
    main()
