"""Unit tests of the upward-light optics validation (blocked fraction, downward control, determinism).

The analytic and geometric references are checked first. The obdeect runs need obdeect.ray_tracing
and are skipped otherwise; the real-model run is skipped unless the simulation models are present
(VALIDATION_MODELS or /models).
"""

import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest

OpticalModel = pytest.importorskip("obdeect.ray_tracing").OpticalModel
oa = pytest.importorskip("obdeect_validation.adapters.optics_directions")
common = pytest.importorskip("obdeect_validation.adapters.common")
FIXTURE = Path(__file__).resolve().parents[1] / oa.FIXTURE


def flat_facet(
    centre, normal, diameter, shape="hexagon_flat_x", radius=0.0, tangent=None
):
    if tangent is None:
        tangent = [1.0, 0.0, 0.0] if abs(normal[0]) < 0.9 else [0.0, 1.0, 0.0]
    return {
        "id": 0,
        "centre_m": list(centre),
        "normal": list(normal),
        "tangent": list(tangent),
        "shape": shape,
        "diameter_m": diameter,
        "radius_of_curvature_m": radius,
    }


def test_hexagon_flat_to_flat_is_the_diameter():
    apothem = 0.5
    grid = np.linspace(-1.0, 1.0, 1201)
    u, v = np.meshgrid(grid, grid, indexing="ij")
    area = np.count_nonzero(oa.hexagon(u, v, apothem)) * (grid[1] - grid[0]) ** 2
    assert area == pytest.approx(2.0 * math.sqrt(3.0) * apothem**2, rel=2e-3)


def test_flat_y_hexagon_is_the_flat_x_hexagon_turned_by_90_degrees():
    apothem = 0.6
    facet_x = {"shape": "hexagon_flat_x", "diameter_m": 2 * apothem}
    facet_y = {"shape": "hexagon_flat_y", "diameter_m": 2 * apothem}

    def inside(facet, u, v):
        return bool(oa.facet_aperture(facet, np.array([u]), np.array([v]))[0])

    # Flats of flat-x are normal to u; flats of flat-y are normal to v.
    assert not inside(facet_x, 1.05 * apothem, 0.0)
    assert inside(facet_y, 1.05 * apothem, 0.0)
    assert inside(facet_x, 0.0, 1.05 * apothem)
    assert not inside(facet_y, 0.0, 1.05 * apothem)


def test_unknown_facet_shape_is_rejected():
    with pytest.raises(oa.AdapterError):
        oa.facet_aperture(
            {"shape": "square", "diameter_m": 1.0}, np.zeros(1), np.zeros(1)
        )


def test_flat_facet_shadow_on_axis_and_tilted_matches_the_projection():
    # A flat hexagon of flat-to-flat 1 m: area sqrt(3)/2; along d its shadow on z = const is
    # A |n . d| / |d_z|; for a horizontal facet at any inclination the shadow area is A.
    facet = flat_facet([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0)
    radius = 3.0
    fraction, single = oa.facet_union_fraction(
        [facet], -1.0, 30.0, radius, log2_points=20
    )
    area = math.sqrt(3.0) / 2.0
    assert fraction * math.pi * radius**2 == pytest.approx(area, rel=5e-3)
    assert single == pytest.approx(fraction, rel=1e-12)


def test_tilted_flat_facet_projects_with_cosine_over_cosine():
    normal = [math.sin(math.radians(20.0)), 0.0, math.cos(math.radians(20.0))]
    facet = flat_facet([0.0, 0.0, 0.0], normal, 1.0, tangent=[0.0, 1.0, 0.0])
    radius = 3.0
    inclination = 35.0
    fraction, _ = oa.facet_union_fraction(
        [facet], -1.0, inclination, radius, log2_points=20
    )
    t = math.radians(inclination)
    direction = np.array([math.sin(t), 0.0, math.cos(t)])
    expected = math.sqrt(3.0) / 2.0 * abs(np.dot(normal, direction)) / direction[2]
    assert fraction * math.pi * radius**2 == pytest.approx(expected, rel=5e-3)


def test_disjoint_shadows_add_up_and_identical_shadows_count_once():
    a = flat_facet([-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.6)
    b = flat_facet([1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.6)
    b["id"] = 1
    fraction, single = oa.facet_union_fraction([a, b], -1.0, 0.0, 3.0, log2_points=18)
    assert fraction == pytest.approx(single, rel=1e-12)
    assert single == pytest.approx(
        2.0 * (2.0 * math.sqrt(3.0) * 0.09) / (math.pi * 9.0), rel=1e-2
    )
    twin = flat_facet([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.6)
    twin2 = dict(twin, id=1)
    fraction, single = oa.facet_union_fraction(
        [twin, twin2], -1.0, 0.0, 3.0, log2_points=18
    )
    assert single == pytest.approx(2.0 * fraction, rel=1e-12)


def test_facet_shadow_disc_contains_every_shadow_point():
    rng = np.random.default_rng(3)
    normal = np.array([0.2, 0.1, math.sqrt(1 - 0.05)])
    facet = flat_facet([3.0, -2.0, 1.2], normal.tolist(), 1.5, radius=40.0)
    t = math.radians(45.0)
    direction = np.array([math.sin(t), 0.0, math.cos(t)])
    cx, cy, delta = oa.facet_shadow_disc(facet, -1.0, direction)
    # Sample the cap (points of the sphere inside the aperture) and project them to z = -1.
    normal_u, tangent, binormal = oa.facet_frame(facet)
    u = rng.uniform(-1.5, 1.5, 200000)
    v = rng.uniform(-1.5, 1.5, 200000)
    keep = oa.facet_aperture(facet, u, v)
    u, v = u[keep], v[keep]
    h = 40.0 - np.sqrt(40.0**2 - u * u - v * v)
    point = np.asarray(facet["centre_m"]) + np.outer(u, tangent) + np.outer(v, binormal)
    point += np.outer(h, normal_u)
    shadow = point[:, :2] - np.outer((point[:, 2] + 1.0), direction[:2]) / direction[2]
    distance = np.hypot(shadow[:, 0] - cx, shadow[:, 1] - cy)
    assert distance.max() <= delta


def test_asphere_segment_on_a_flat_asphere_projects_to_its_hexagon():
    asphere = {
        "vertex_z_m": 0.4,
        "inner_radius_m": 0.0,
        "outer_radius_m": 5.0,
        "radial_scale_m": 1.0,
        "coefficients_m": [0.0],
    }
    segment = {
        "id": 1,
        "shape": "hexagon",
        "centre_xy_m": [1.0, 0.5],
        "diameter_m": 0.8,
        "rotation_deg": 0.0,
    }
    apothem = 0.4
    hex_area = 2.0 * math.sqrt(3.0) * apothem**2
    assert oa.segment_projected_area(
        asphere, segment, np.array([0.0, 0.0, 1.0])
    ) == pytest.approx(hex_area, rel=1e-3)
    t = math.radians(40.0)
    direction = np.array([math.sin(t), 0.0, math.cos(t)])
    assert oa.segment_projected_area(asphere, segment, direction) == pytest.approx(
        hex_area, rel=1e-3
    )


def test_sloped_asphere_segment_mapping_is_the_tangent_plane_projection():
    # The radial tangent coordinate of a point is D . e_r' with the unit tangent
    # e_r' = (e_r, slope) / sqrt(1 + slope^2), D the offset from the segment centre on the asphere;
    # the azimuthal coordinate is D . e_phi. The aperture (rotation 0) is the hexagon in these.
    asphere = {
        "vertex_z_m": 0.0,
        "inner_radius_m": 0.0,
        "outer_radius_m": 4.0,
        "radial_scale_m": 0.5,
        "coefficients_m": [0.0, 0.3, -0.05],
    }
    segment = {
        "id": 1,
        "shape": "hexagon",
        "centre_xy_m": [1.6, -0.9],
        "diameter_m": 0.9,
        "rotation_deg": 0.0,
    }
    cx, cy = segment["centre_xy_m"]
    radius = math.hypot(cx, cy)
    er = np.array([cx, cy]) / radius
    ephi = np.array([-er[1], er[0]])
    q = (radius / asphere["radial_scale_m"]) ** 2
    coefficients = asphere["coefficients_m"]
    slope = (
        radius
        / asphere["radial_scale_m"] ** 2
        * sum(2 * i * c * q ** (i - 1) for i, c in enumerate(coefficients) if i)
    )
    axis = np.linspace(-0.6, 0.6, 301)
    x, y = (g.ravel() for g in np.meshgrid(cx + axis, cy + axis, indexing="ij"))
    inside = oa.segment_aperture(asphere, segment, x, y)
    z = oa.asphere_height(asphere, x, y)
    offset = np.stack([x - cx, y - cy, z - oa.asphere_height(asphere, cx, cy)], axis=1)
    radial = (offset[:, :2] @ er + slope * offset[:, 2]) / math.sqrt(
        1.0 + slope * slope
    )
    azimuthal = offset[:, :2] @ ephi
    in_plane = np.outer(radial, er) + np.outer(azimuthal, ephi)
    expected = oa.hexagon(in_plane[:, 0], in_plane[:, 1], 0.45) & (
        np.hypot(x, y) <= 4.0
    )
    assert np.count_nonzero(inside != expected) <= 5
    assert 0.0 < inside.mean() < 1.0


def test_bitwise_row_comparison_counts_changes_and_equal_nans():
    a = {"x": np.array([1.0, np.nan, 3.0]), "k": np.array([1, 2, 3], dtype=np.uint64)}
    b = {"x": np.array([1.0, np.nan, 3.0]), "k": np.array([1, 2, 3], dtype=np.uint64)}
    assert not oa.differing_rows(a, b).any()
    b["x"][2] = np.nextafter(3.0, 4.0)
    b["k"][0] = 9
    assert oa.differing_rows(a, b).tolist() == [True, False, True]


def test_binomial_pull_definition():
    assert oa.binomial_pull(500, 1000, 0.5) == pytest.approx(0.0)
    assert oa.binomial_pull(540, 1000, 0.5) == pytest.approx(40 / math.sqrt(250.0))


def test_obscurer_below_the_primary_is_rejected():
    document = json.loads(FIXTURE.read_text(encoding="utf-8"))
    document["optics"]["obscurers"][0]["centre_m"][2] = -0.5
    with pytest.raises(oa.AdapterError, match="reaches the height of the primary"):
        oa.Geometry.of(document)


def test_fixture_blocked_fraction_matches_the_reference():
    model = OpticalModel.load(FIXTURE)
    geometry = oa.Geometry.of(model.document)
    for inclination in oa.INCLINATIONS_DEG:
        record, _ = oa.beam_part(model, geometry, inclination, 200_000, seed=5)
        assert record["detected"] == 0
        assert record["status_sum_error"] == 0
        assert record["blocked_unknown_surface"] == 0
        assert abs(record["mirror_back_pull"]) < oa.PULL_LIMIT


def test_fixture_downward_rows_are_identical_to_downward_only_trace():
    model = OpticalModel.load(FIXTURE)
    geometry = oa.Geometry.of(model.document)
    record = oa.downward_part(model, geometry, 20.0, 20_000, seed=9)
    assert record["downward_mismatch_rows"] == 0
    assert record["status_sum_error"] == 0


def test_fixture_determinism_across_threads_pools_chunks_and_order():
    model = OpticalModel.load(FIXTURE)
    geometry = oa.Geometry.of(model.document)
    record = oa.determinism_part(model, geometry, 30_000, seed=11)
    assert record["determinism_mismatch_count"] == 0
    assert sum(record["status_counts"].values()) == record["photons"]


def test_summary_flags_detected_and_unknown_surfaces():
    fake_beam = {
        "detected": 3,
        "status_sum_error": 0,
        "mirror_back_pull": 0.5,
        "blocked_unknown_surface": 0,
        "blocked_obscurer": 1,
        "photons": 10,
    }
    section = {
        "label": "LST",
        "beams": [fake_beam],
        "downward": [{"downward_mismatch_rows": 0, "status_sum_error": 0}],
        "determinism": {
            "determinism_mismatch_count": 0,
            "photons": 4,
            "status_counts": {"detected": 4},
        },
    }
    metrics = oa.summarise([section])
    assert metrics["detected_count"] == 3
    assert metrics["mirror_back_cases"] == 1


@pytest.mark.skipif(
    not (common.models_root() / "productions" / oa.MODEL_VERSION).is_dir(),
    reason="simulation models not available (set VALIDATION_MODELS)",
)
def test_lst_blocked_fraction_matches_the_reference():
    model = OpticalModel.from_simulation_models(
        common.models_root(),
        "LSTN-01",
        oa.MODEL_VERSION,
        configuration_seed=oa.CONFIGURATION_SEED,
        zenith_deg=oa.ZENITH_DEG,
    )
    geometry = oa.Geometry.of(model.document)
    record, _ = oa.beam_part(model, geometry, 20.0, 200_000, seed=13)
    assert record["detected"] == 0
    assert record["blocked_unknown_surface"] == 0
    assert abs(record["mirror_back_pull"]) < oa.PULL_LIMIT


def test_copy_of_document_is_not_changed_by_the_geometry_reference():
    model = OpticalModel.load(FIXTURE)
    before = copy.deepcopy(model.document)
    oa.Geometry.of(model.document).reference_fraction(20.0)
    assert model.document == before
