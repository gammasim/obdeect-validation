"""Pure helpers of the CORSIKA CURVED geometry case, on synthetic data (no IACT file needed).

Slant depths are checked against the flat limit (R -> 1e12 m), the vertical column, the ordering
curved < plane-parallel, a direct path integral along the axis, and the weighted median, mean
zenith, core-selection and bunch-format helpers. The input-missing path is checked with a
temporary inputs folder.
"""

import math
import os
from types import SimpleNamespace

import numpy as np
import pytest

from obdeect_validation.adapters import curvature_corsika_geometry as geometry_case

EARTH = geometry_case.EARTH_RADIUS_M
SITE = geometry_case.SITE_ALTITUDE_M


@pytest.fixture(scope="module")
def profile(tmp_path_factory):
    """Synthetic tropical-like profile: exponential layers with thickness from the density."""
    h_km = np.array([0.0, 5.0, 10.0, 20.0, 30.0, 50.0, 80.0, 120.0])
    rho = np.array([1.2e-3, 7.2e-4, 4.2e-4, 9.5e-5, 1.8e-5, 1.0e-6, 1.8e-8, 2.0e-10])
    path = tmp_path_factory.mktemp("atmosphere") / "profile.dat"
    lines = ["# synthetic profile", "# alt rho thick n-1"]
    lines += [f"{h:10.3f} {r:14.6e} {0.0:14.6e} {0.0:14.6e}" for h, r in zip(h_km, rho)]
    path.write_text("\n".join(lines) + "\n")
    return geometry_case.read_profile(path)


def _analytic_thick(profile, h_m):
    """Vertical column from the closed form of each exponential segment (independent check)."""
    total = 0.0
    for h0, h1, rate, rho0 in geometry_case._segments(profile):
        if h1 <= h_m:
            continue
        lo = max(h_m, h0)
        rho_lo = rho0 * math.exp(rate * (lo - h0))
        rho_hi = rho0 * math.exp(rate * (h1 - h0))
        total += (rho_lo - rho_hi) / (-rate)
    return total * 100.0


def _direct_path_slant(profile, h_m, theta_deg, earth_radius=EARTH, points=400001):
    """Slant depth by integrating rho along the path s (independent of the change of variables)."""
    zenith = math.radians(theta_deg)
    radius_site = earth_radius + SITE
    cos_z, sin_z = math.cos(zenith), math.sin(zenith)

    def path_length(altitude):
        r = earth_radius + altitude
        return math.sqrt(r * r - (radius_site * sin_z) ** 2) - radius_site * cos_z

    s_start = path_length(h_m)
    s_end = path_length(profile["h_m"][-1])
    s = np.linspace(s_start, s_end, points)
    altitude = (
        np.sqrt(radius_site**2 + 2 * radius_site * s * cos_z + s * s) - earth_radius
    )
    density = np.exp(np.interp(altitude, profile["h_m"], np.log(profile["rho"])))
    return float(np.trapezoid(density, s)) * 100.0


# ---- slant depth ---------------------------------------------------------------------------


@pytest.mark.parametrize("zenith", [0.0, 30.0, 60.0, 75.0])
@pytest.mark.parametrize("altitude", [SITE, 10000.0, 20000.0])
def test_flat_limit_reproduces_plane_parallel(profile, zenith, altitude):
    # The residual grows with h / (R cos^2 z) and is only checked where it is small (up to 75 deg).
    curved = geometry_case.slant_depth_curved(
        profile, altitude, zenith, earth_radius_m=1e12
    )
    plane = geometry_case.slant_depth_plane(profile, altitude, zenith)
    assert curved == pytest.approx(plane, rel=1e-6)


@pytest.mark.parametrize("zenith", [30.0, 65.0, 75.0])
@pytest.mark.parametrize("altitude", [SITE + 1.0, 8000.0, 15000.0, 25000.0])
def test_curved_slant_depth_below_plane_parallel(profile, zenith, altitude):
    # The axis rises away from the Earth's surface: ds/dh <= 1/cos z for h >= site altitude.
    curved = geometry_case.slant_depth_curved(profile, altitude, zenith)
    plane = geometry_case.slant_depth_plane(profile, altitude, zenith)
    assert 0.0 < curved < plane


def test_vertical_case_equals_thickness_exactly(profile):
    for altitude in (SITE, 10000.0, 25000.0, 60000.0):
        curved = geometry_case.slant_depth_curved(profile, altitude, 0.0)
        assert curved == geometry_case.thick_from_density(profile, altitude)
        assert curved == pytest.approx(_analytic_thick(profile, altitude), rel=1e-10)


def test_slant_depth_decreases_with_altitude(profile):
    altitudes = np.linspace(SITE, 100000.0, 60)
    values = [geometry_case.slant_depth_curved(profile, h, 75.0) for h in altitudes]
    assert np.all(np.diff(values) < 0.0)


@pytest.mark.parametrize("zenith", [45.0, 75.0])
def test_curved_slant_depth_matches_direct_path_integral(profile, zenith):
    altitude = 12000.0
    direct = _direct_path_slant(profile, altitude, zenith)
    assert geometry_case.slant_depth_curved(profile, altitude, zenith) == pytest.approx(
        direct, rel=1e-6
    )


def test_zenith_must_be_below_horizon(profile):
    with pytest.raises(ValueError, match="zenith"):
        geometry_case.slant_depth_curved(profile, 10000.0, 90.0)
    with pytest.raises(ValueError, match="zenith"):
        geometry_case.slant_depth_plane(profile, 10000.0, -1.0)


def test_altitude_below_profile_is_rejected(profile):
    with pytest.raises(ValueError, match="below the atmosphere"):
        geometry_case.slant_depth_curved(profile, -10.0, 0.0)


def test_profile_reader_rejects_non_monotonic_altitudes(tmp_path):
    path = tmp_path / "bad.dat"
    path.write_text("0.0 1e-3 1000.0 0.0\n0.0 9e-4 900.0 0.0\n")
    with pytest.raises(ValueError, match="altitudes must rise"):
        geometry_case.read_profile(path)


# ---- helpers on arrays ---------------------------------------------------------------------


def test_photon_weighted_median_uses_cumulative_weight():
    assert geometry_case.photon_weighted_median([1.0, 2.0, 3.0], [1.0, 1.0, 8.0]) == 3.0
    assert geometry_case.photon_weighted_median([4.0, 1.0, 3.0, 2.0], [1.0] * 4) == 2.0
    assert geometry_case.photon_weighted_median([5.0, 7.0], [0.0, 1.0]) == 7.0


def test_photon_weighted_median_rejects_bad_weights():
    with pytest.raises(ValueError):
        geometry_case.photon_weighted_median([1.0, 2.0], [0.0, 0.0])
    with pytest.raises(ValueError):
        geometry_case.photon_weighted_median([1.0], [1.0, 2.0])


@pytest.mark.parametrize("zenith", [65.0, 75.0])
def test_mean_zenith_of_cherenkov_cone_is_axis_zenith(zenith):
    """Photon directions on a 1 deg cone about the axis (Cherenkov emission) average to the axis.

    The residual is second order, about alpha^2 cot(theta) / 2: 0.004 deg at 65 deg and 0.002 deg
    at 75 deg for alpha = 1 deg (a cone about a near-vertical axis would instead average to alpha).
    """
    theta = math.radians(zenith)
    alpha = math.radians(1.0)
    axis = np.array([math.sin(theta), 0.0, math.cos(theta)])
    e1 = np.array([math.cos(theta), 0.0, -math.sin(theta)])
    e2 = np.array([0.0, 1.0, 0.0])
    phi = np.linspace(0.0, 2.0 * math.pi, 3600, endpoint=False)
    direction = math.cos(alpha) * axis[:, None] + math.sin(alpha) * (
        e1[:, None] * np.cos(phi) + e2[:, None] * np.sin(phi)
    )
    cx, cy = direction[0], direction[1]
    mean = geometry_case.mean_zenith_deg(cx, cy, np.ones_like(cx))
    assert mean == pytest.approx(zenith, abs=0.01)


def test_mean_zenith_known_values_and_weights():
    cx = np.array([math.sin(math.radians(58.0)), math.sin(math.radians(62.0))])
    cy = np.zeros(2)
    assert geometry_case.mean_zenith_deg(cx, cy, [1.0, 1.0]) == pytest.approx(
        60.0, abs=1e-12
    )
    assert geometry_case.mean_zenith_deg(cx, cy, [3.0, 1.0]) == pytest.approx(
        59.0, abs=1e-12
    )
    straight = geometry_case.mean_zenith_deg([0.0], [0.0], [5.0])
    assert straight == pytest.approx(0.0, abs=1e-12)


def test_within_radius_is_inclusive_and_centred():
    x = np.array([149.0, 150.0, 151.0])
    y = np.zeros(3)
    mask = geometry_case.within_radius(x, y, 0.0, 0.0, 150.0)
    assert mask.tolist() == [True, True, False]
    shifted = geometry_case.within_radius(x + 10.0, y, 10.0, 0.0, 150.0)
    assert shifted.tolist() == [True, True, False]


def test_ground_impact_adds_detector_position_and_slope_shift():
    # straight down: ground impact = detector position + bunch offset
    x, y = geometry_case.ground_impact_m([1.0], [2.0], [0.0], [0.0], 10.0, 300.0, -50.0)
    assert x[0] == pytest.approx(301.0) and y[0] == pytest.approx(-48.0)
    # inclined: shift = (c / |cz|) * z0, with cz = sqrt(1 - c^2)
    c = 0.5
    x, _ = geometry_case.ground_impact_m([0.0], [0.0], [c], [0.0], 10.0, 0.0, 0.0)
    assert x[0] == pytest.approx(c / math.sqrt(1.0 - c * c) * 10.0)


class _Stream(list):
    """Stand-in for EventIOFile: a context manager over its top-level objects."""

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class _Container:
    """Stand-in for an EventIO container: iterating yields its children."""

    def __init__(self, type_, version, children=()):
        self.header = SimpleNamespace(
            type=type_, version=version, only_subobjects=bool(children)
        )
        self._children = list(children)

    def __iter__(self):
        return iter(self._children)


def test_photon_blocks_are_found_below_containers(monkeypatch):
    tree = [
        _Container(1202, 0),
        _Container(
            1203,
            0,
            [
                _Container(
                    1204,
                    0,
                    [
                        _Container(1205, 0),
                        _Container(1205, 1000),
                    ],
                ),
                _Container(1204, 0, [_Container(1205, 1001)]),
            ],
        ),
    ]
    monkeypatch.setattr(geometry_case, "EventIOFile", lambda path: _Stream(tree))
    assert geometry_case.photon_block_versions("unused.iact") == [0, 1000, 1001]
    monkeypatch.setattr(
        geometry_case, "EventIOFile", lambda path: _Stream([_Container(1202, 0)])
    )
    with pytest.raises(RuntimeError, match="no photon bunch block"):
        geometry_case.photon_block_versions("unused.iact")


@pytest.mark.skipif(
    not os.environ.get("OBDEECT_TEST_IACT"),
    reason="set OBDEECT_TEST_IACT to a flat IACT file (long format from the legacy card)",
)
def test_real_file_photon_blocks_are_long_format():
    versions = geometry_case.photon_block_versions(os.environ["OBDEECT_TEST_IACT"])
    assert versions
    assert geometry_case.compact_block_count(versions) == 0


def test_compact_block_count_uses_the_eventio_version():
    """Long blocks have version 0, compact blocks 1000 + n (io_simtel.c)."""
    assert geometry_case.is_compact_version(1000)
    assert geometry_case.is_compact_version(1001)
    assert not geometry_case.is_compact_version(0)
    assert not geometry_case.is_compact_version(999)
    assert geometry_case.compact_block_count([0, 1000, 0, 1005]) == 2
    assert geometry_case.compact_block_count([0, 0]) == 0


# ---- inputs --------------------------------------------------------------------------------


def test_missing_inputs_are_listed_and_exit_code_is_3(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("VALIDATION_INPUTS", str(tmp_path))
    missing = geometry_case.check_inputs(tmp_path)
    assert len(missing) == 4  # two shower files and two logs
    assert all("curved" in name for name in missing)
    assert geometry_case.main() == 3
    assert "missing input" in capsys.readouterr().err


def test_curved_sources_are_the_two_zenith_angles():
    names = [name for name, *_ in geometry_case._sources()]
    assert names == [
        "gamma-65deg-300-5000GeV-curved.iact",
        "gamma-75deg-1000-10000GeV-curved.iact",
    ]
