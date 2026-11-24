"""Azimuth convention of the obdeect sources presets against the CORSIKA convention of the stream.

Two subcommands (two cases share this module).

``presets``: every preset that takes a zenith and azimuth (muon, shower, laser-beam, lidar,
beam-scan, night-sky) is run at zenith 30 deg and azimuth 0, 90, 180 and 270 deg. The direction
of motion of the light is recovered as the weighted mean of the photon direction cosines
(cx, cy, cz) of the output stream and compared with the expected direction below.

``corsika``: a CORSIKA IACT file with known THETAP and PHIP (from its .corsika.log) is converted
with ``obdeect stream from-corsika``. The event azimuth_rad and zenith_rad of the stream, and the
mean photon direction, must match PHIP and THETAP. The same mean is also taken from the EventIO
IACT bunches directly (independent of obdeect) and compared with the stream.

Convention of the presets (the reference of the ``presets`` part)
-----------------------------------------------------------------
Sky frame (docs/sources/spec/conventions.md, docs/stream/format.md): origin at the telescope,
x north, y west, z up. Zenith and azimuth are astronomical: zenith from up, azimuth clockwise
from north towards east, so the unit vector of zenith z and azimuth a is

    astro(z, a) = (sin z cos a, -sin z sin a, cos z)

Rule (one for every preset): sky sources (muon, shower, night-sky) give the direction the light
COMES FROM, as seen by the observer; the light travels along -astro(z, a). Beam devices (laser-
beam, lidar, beam-scan) give the direction the beam POINTS TO; the light travels along astro(z, a).
The disagreement counts compare the physical quantity each role names (source direction for sky
sources, pointing for beams), see ``pair_disagrees``.

Documented convention of the stream (the CORSIKA reference of the ``corsika`` part)
-----------------------------------------------------------------------------------
A direction of motion is the unit vector

    u(theta, phi) = (sin theta cos phi, sin theta sin phi, -cos theta),

with theta the angle of the motion from straight down (CORSIKA THETAP; 0 is vertical) and phi
the azimuth of the motion measured from north (x) towards west (y), i.e. counter-clockwise seen
from above. This is the CORSIKA PHI (EVTH word 12, stream event key ``azimuth_rad``) and the
convention of the photon direction cosines of the stream (grisu.hpp: ``az_c = atan2(cy, cx)``
for photons).

Evidence for this convention (see the corsika case):

- The stream direction cosines are those of the motion: PHIP 180 (gamma-20deg-az180) gives a
  mean cx of -0.340 = sin(20 deg) cos(180 deg) (motion to the south). The CORSIKA log of that
  run prints "Azimuth=180.00 deg (S->E) which is 360.00 deg (N->E)".
- PHIP 0 (gamma-20deg-20-100GeV, and the card named "South") gives cx = +0.343: motion to the
  north, so the shower arrives from the south, as simtools maps it (simtools
  ``utils/geometry.py:geographic_to_corsika_azimuth``: PHI = 180 - az, az clockwise from north
  for the arrival direction). A PHI = 90 event (motion west) is not available here; the west sign
  of y rests on this simtools mapping and on the frame of obdeect grisu.md (x north, y west).

Presets and their roles (zenith = the zenith given to the preset, theta = angle from down):

- muon, shower (role ``source``): ``--zenith``, ``--azimuth`` give the direction the light comes
  from. Expected travel: -astro(zenith, azimuth).
- night-sky (role ``source``): ``--telescope-zenith``, ``--telescope-azimuth`` give the pointing,
  the sky direction the photons come from. Expected travel: -astro(zenith, azimuth).
- laser-beam, lidar, beam-scan (role ``beam``): ``--beam-zenith``, ``--beam-azimuth`` give the
  direction the beam points to (zenith 180 is down). The preset is run with zenith 150 (30 deg
  from down). Expected travel: astro(150, azimuth).

point-source takes a Cartesian ``--direction`` (no azimuth) and flasher takes a position angle
of its flashers; neither takes a compass azimuth, so neither is in this case.

Before the presets followed this rule the case measured 16 violations: the muon and shower travel
was rotated (sign of the azimuth and a 180 deg offset), and night-sky was compared with the
reverse of its pointing. The beams were already correct for the rule.

Gates and their justification
-----------------------------
- ``convention_violations`` (expected 0, atol 0). Count of (preset, azimuth) pairs whose measured
  travel direction differs from the expected one by more than
  ``max(ANGLE_TOL_RAD, K * se)``. ``ANGLE_TOL_RAD = 1e-3`` rad (0.06 deg) is the requirement: a
  convention error is an angle of at least 2 sin(30 deg) = 1.0 rad (rotation by 180 deg) or
  2 sin(30 deg) sin(|a|) for a sign error, far above it, while the direction of a direct
  emission (laser, beam) is exact to 1e-5 rad. ``se`` is the standard error of the weighted mean
  of the photon directions, ``sqrt(sum w_i^2 |v_i - u|^2) / sum w_i``; ``K = 5``.
- ``presets_disagree_90`` (expected 0, atol 0). Number of preset pairs at azimuth 90 deg whose
  named directions differ by more than ``max(ANGLE_TOL_RAD, K * sqrt(se_i^2 + se_j^2))``. The
  named direction is the physical quantity the role names: the source direction (-travel) of the
  sky sources muon, shower and night-sky, the pointing (travel) of the beams laser-beam, lidar and
  beam-scan. Pairs are taken within one role only: a muon and a laser at the same azimuth point in
  opposite directions by design, so comparing them would test nothing.
- ``presets_disagree_all`` (expected 0, atol 0): the same count summed over the four azimuths.
- The muon and shower emit Cherenkov light on a cone around the track (half angle about 0.024
  rad). The mean of the cone photons recovers the track axis, and the standard error above
  includes the cone spread. The shower uses a fluorescence yield of 1e-12 (the input must be
  finite and > 0): its isotropic fluorescence photons are then negligible, and the photons are
  the Cherenkov cones of the shower steps.
- Seeds are fixed (``--seed 20261123``); the inputs are the run descriptions of the presets.

Corsika subcommand gates (the stream and the EventIO reader must agree with the CORSIKA log):

- ``azimuth_error_deg_max``, ``zenith_error_deg_max`` (expected 0, atol 0.001 deg): EVTH words are
  float32; an ulp at pi is 1.4e-5 deg, so 0.001 deg is about 70 ulp.
- ``stream_vs_eventio_delta_rad_max`` (expected 0, atol 1e-6 rad): the stream and the EventIO
  reader read the same bunches (float32 in the stream, float64 in the reader).
- ``direction_error_rad_max`` (expected 0, atol 0.05 rad). Requirement: the mean photon direction
  must follow the CORSIKA travel direction to 0.05 rad (2.9 deg), and it must separate the travel
  hypothesis from the arrival hypothesis, which differ by 2 sin(theta) = 0.68 rad at theta = 20
  deg. A shower mean is not exactly the axis (lateral and geometric effects of a finite shower),
  which is why the requirement is not 1e-3 rad. ``K * se`` is added (negligible at 1e6 photons).
  Honest note: the mean photon directions were inspected before this threshold was written
  (cx = -0.3399 against -0.3420 for PHIP 180); the threshold is fixed from the separation of the
  two hypotheses, not from that residual.
- ``arrival_closer_files`` (expected 0, atol 0): number of files whose mean direction is closer to
  the arrival hypothesis (-u) than to the travel hypothesis u.

Inputs (CORSIKA IACT files and logs, not in Git) are listed in the case JSON; the case is
blocked when they are missing. Photons are read with ``readers.read_native`` (the obdeect
stream) and ``readers.read_iact`` (EventIO).
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common, readers

#: Zenith of the directions in degrees (the zenith of the motion from straight down).
ZENITH_DEG = 30.0
AZIMUTHS_DEG = (0.0, 90.0, 180.0, 270.0)
#: Requirement on the direction of a preset, rad (see the docstring).
ANGLE_TOL_RAD = 1e-3
#: Requirement on the mean direction of a CORSIKA shower, rad (see the docstring).
CORSIKA_TOL_RAD = 0.05
#: Multiple of the standard error added to a tolerance.
K_SE = 5.0
SEED = 20261123
#: Fluorescence yield that isolates the Cherenkov cones of the shower (finite and > 0 is required).
SHOWER_FLUORESCENCE_YIELD = 1e-12
#: Angle between two directions below which they count as equal (rad, numerical only).
_EPS = 1e-12

#: Role of a preset: sky source (zenith/azimuth = direction the light comes from) or beam device
#: (zenith/azimuth = direction the beam points to). See the docstring.
ROLE_SOURCE = "source"
ROLE_BEAM = "beam"


def corsika_direction(theta_deg: float, phi_deg: float) -> np.ndarray:
    """Unit vector of the motion with polar angle theta from straight down and CORSIKA PHI.

    Sky frame: x north, y west, z up. phi is measured from north towards west.
    """
    theta = math.radians(theta_deg)
    phi = math.radians(phi_deg)
    return np.array(
        [
            math.sin(theta) * math.cos(phi),
            math.sin(theta) * math.sin(phi),
            -math.cos(theta),
        ]
    )


def angle_rad(a: np.ndarray, b: np.ndarray) -> float:
    """Angle between two vectors [rad]; atan2 keeps tiny angles accurate."""
    cross = float(np.linalg.norm(np.cross(a, b)))
    return math.atan2(cross, float(np.dot(a, b)))


def mean_direction(cx, cy, cz, weights) -> tuple[np.ndarray, float, int]:
    """Weighted mean direction of photons (unit vector), its standard error [rad], and the rows.

    The standard error is sqrt(sum w_i^2 |v_i - u|^2) / sum w_i for unit vectors v_i. It is zero
    for a single direction and grows with the spread of the photons around the mean.
    """
    v = np.column_stack(
        [np.asarray(cx, float), np.asarray(cy, float), np.asarray(cz, float)]
    )
    w = np.asarray(weights, float)
    if v.shape[0] == 0 or not np.isfinite(w).all() or w.sum() <= 0.0:
        raise ValueError("no photons with positive weight to take a direction from")
    mean = (w[:, None] * v).sum(axis=0) / w.sum()
    u = mean / np.linalg.norm(mean)
    spread = np.linalg.norm(v - u, axis=1)
    se = math.sqrt(float(np.sum((w * spread) ** 2))) / float(w.sum())
    return u, se, int(v.shape[0])


def tolerance(se_rad: float, base: float = ANGLE_TOL_RAD) -> float:
    """Gate on an angle: the requirement or K standard errors, whichever is larger."""
    return max(base, K_SE * se_rad)


# ------------------------------------------------------------------ presets


def astro_unit(zenith_deg: float, azimuth_deg: float) -> np.ndarray:
    """Unit vector of a zenith and a clockwise-from-north azimuth (x north, y west, z up)."""
    zenith = math.radians(zenith_deg)
    azimuth = math.radians(azimuth_deg)
    return np.array(
        [
            math.sin(zenith) * math.cos(azimuth),
            -math.sin(zenith) * math.sin(azimuth),
            math.cos(zenith),
        ]
    )


class Preset:
    """A preset run: how to call it and what its zenith/azimuth mean (see the docstring)."""

    def __init__(self, name, zenith_input, role, argv_builder):
        self.name = name
        self.zenith_input = zenith_input  # zenith given to the preset [deg]
        self.role = role
        self.argv_builder = argv_builder

    def expected(self, azimuth_deg: float) -> np.ndarray:
        """Expected direction of motion of the light for azimuth_deg (the rule of the docstring)."""
        direction = astro_unit(self.zenith_input, azimuth_deg)
        return -direction if self.role == ROLE_SOURCE else direction

    def named(self, travel: np.ndarray) -> np.ndarray:
        """The quantity the role names: source direction (-travel) or pointing (travel)."""
        return -travel if self.role == ROLE_SOURCE else travel

    def argv(self, azimuth_deg: float, out: Path) -> list:
        return self.argv_builder(self.zenith_input, azimuth_deg) + ["-o", str(out)]


def _muon_argv(zenith, azimuth):
    return [
        "sources", "run", "muon", "--observation-level", "2150",
        "--wavelength-min", "300", "--wavelength-max", "600", "--n-minus-one", "2.9e-4",
        "--telescope-radius", "12.5", "--start-height", "10000",
        "--zenith", str(zenith), "--azimuth", str(azimuth), "--events", "1",
        "--seed", str(SEED),
    ]  # fmt: skip


def _shower_argv(zenith, azimuth):
    return [
        "sources", "run", "shower", "--observation-level", "2150",
        "--wavelength-min", "300", "--wavelength-max", "600", "--n-minus-one", "2.9e-4",
        "--telescope-radius", "12.5", "--start-height", "10000",
        "--zenith", str(zenith), "--azimuth", str(azimuth),
        "--step", "100", "--nmax", "100", "--xmax", "600", "--x0", "-50", "--gh-lambda", "70",
        "--atmosphere-depth0", "800", "--scale-height", "8000",
        "--fluorescence-yield", str(SHOWER_FLUORESCENCE_YIELD),
        "--fluorescence-spectrum", "fixed:wavelength_nm=337",
        "--seed", str(SEED),
    ]  # fmt: skip


def _laser_argv(zenith, azimuth):
    return [
        "sources", "run", "laser-beam", "--observation-level", "2150",
        "--laser-position", "200,0,3000", "--beam-zenith", str(zenith),
        "--beam-azimuth", str(azimuth), "--divergence", "0.001", "--pulse-energy", "1e-3",
        "--wavelength-nm", "355", "--telescope-radius", "12.5", "--rows", "200",
        "--allow-direct-beam", "--seed", str(SEED),
    ]  # fmt: skip


def _lidar_argv(zenith, azimuth):
    return [
        "sources", "run", "lidar", "--observation-level", "2000",
        "--satellite-altitude", "400000", "--satellite-offset", "1000,0",
        "--beam-zenith", str(zenith), "--beam-azimuth", str(azimuth), "--divergence", "0.01",
        "--pulse-energy", "1e-3", "--wavelength-nm", "355", "--telescope-radius", "12.5",
        "--start-height", "100000", "--rows", "128", "--allow-direct-beam", "--seed", str(SEED),
    ]  # fmt: skip


def _beam_scan_argv(zenith, azimuth):
    return [
        "sources", "run", "beam-scan", "--observation-level", "2150",
        "--spectrum", "fixed:wavelength_nm=400", "--photons-per-position", "1e4",
        "--camera-radius", "1.5", "--scan-start", "-1", "-1", "--scan-extent-x", "1",
        "--scan-extent-y", "1", "--scan-step", "1", "--height", "100",
        "--beam-zenith", str(zenith), "--beam-azimuth", str(azimuth), "--seed", str(SEED),
    ]  # fmt: skip


def _night_sky_argv(zenith, azimuth):
    return [
        "sources", "run", "night-sky", "--observation-level", "2150",
        "--spectrum", "uniform:low_nm=300,high_nm=600", "--photons", "1e5",
        "--viewcone-deg", "2", "--telescope-position", "0,0,0", "--telescope-radius", "12.5",
        "--telescope-zenith", str(zenith), "--telescope-azimuth", str(azimuth),
        "--seed", str(SEED),
    ]  # fmt: skip


PRESETS = (
    Preset("muon", ZENITH_DEG, ROLE_SOURCE, _muon_argv),
    Preset("shower", ZENITH_DEG, ROLE_SOURCE, _shower_argv),
    # Beam zenith counts from up (180 is down): 150 is 30 deg from down.
    Preset("laser-beam", 180.0 - ZENITH_DEG, ROLE_BEAM, _laser_argv),
    Preset("lidar", 180.0 - ZENITH_DEG, ROLE_BEAM, _lidar_argv),
    Preset("beam-scan", 180.0 - ZENITH_DEG, ROLE_BEAM, _beam_scan_argv),
    Preset("night-sky", ZENITH_DEG, ROLE_SOURCE, _night_sky_argv),
)


def measure_preset(preset: Preset, azimuth_deg: float, root: Path) -> dict:
    """Run one preset at one azimuth, read its photons and return the measured direction."""
    tag = f"{preset.name}-az{int(azimuth_deg)}"
    out = root / f"{tag}.pst"
    common.command(
        [common.tool("obdeect"), *preset.argv(azimuth_deg, out)], root / f"{tag}.log"
    )
    _, data = readers.read_native(
        out, keys=("cx", "cy", "cz", "w"), require_events=False
    )
    u, se, rows = mean_direction(data["cx"], data["cy"], data["cz"], data["w"])
    out.unlink()  # the stream is large; the measurement and the log are kept
    expected = preset.expected(azimuth_deg)
    angle = angle_rad(u, expected)
    return {
        "preset": preset.name,
        "azimuth_deg": azimuth_deg,
        "zenith_input_deg": preset.zenith_input,
        "role": preset.role,
        "measured": u.tolist(),
        "named": preset.named(u).tolist(),
        "expected": expected.tolist(),
        "angle_deg": math.degrees(angle),
        "se_rad": se,
        "tolerance_rad": tolerance(se),
        "rows": rows,
        "violation": angle > tolerance(se),
    }


def pair_disagrees(first: dict, second: dict) -> bool:
    """True if two named directions (same azimuth and role) differ by more than their tolerance."""
    angle = angle_rad(np.array(first["named"]), np.array(second["named"]))
    limit = max(ANGLE_TOL_RAD, K_SE * math.hypot(first["se_rad"], second["se_rad"]))
    return angle > limit


def presets_part(root: Path) -> dict:
    """Measure every preset at every azimuth; count the convention violations and disagreements."""
    rows = []
    for azimuth in AZIMUTHS_DEG:
        for preset in PRESETS:
            rows.append(measure_preset(preset, azimuth, root))
    violations = [r for r in rows if r["violation"]]
    disagree = {}
    for azimuth in AZIMUTHS_DEG:
        at = [r for r in rows if r["azimuth_deg"] == azimuth]
        count = 0
        pairs = []
        for i, first in enumerate(at):
            for second in at[i + 1 :]:
                if first["role"] != second["role"]:
                    continue  # different physical quantities (source direction vs pointing)
                if pair_disagrees(first, second):
                    count += 1
                    pairs.append([first["preset"], second["preset"]])
        disagree[azimuth] = {"count": count, "pairs": pairs}
    return {
        "rows": rows,
        "convention_violations": len(violations),
        "violations": [[r["preset"], r["azimuth_deg"]] for r in violations],
        "presets_disagree_90": disagree[90.0]["count"],
        "presets_disagree_all": sum(d["count"] for d in disagree.values()),
        "disagreements": {str(int(a)): d for a, d in disagree.items()},
    }


# ------------------------------------------------------------------ corsika


def parse_corsika_angles(log_path: Path) -> tuple[float, float]:
    """THETAP and PHIP (first value of each line) from a CORSIKA log that echoes its input."""
    text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    values = {}
    for key in ("THETAP", "PHIP"):
        match = re.search(rf"^{key}\s+([-+0-9.eE]+)", text, re.MULTILINE)
        if match is None:
            raise ValueError(f"{log_path}: no {key} line")
        values[key] = float(match.group(1))
    return values["THETAP"], values["PHIP"]


def _wrap_deg(value: float) -> float:
    """Signed difference of two azimuths in degrees, in [-180, 180)."""
    return (value + 180.0) % 360.0 - 180.0


def corsika_file(path: Path, root: Path) -> dict:
    """Stream and EventIO directions of one CORSIKA IACT file against its log."""
    path = Path(path)
    theta, phi = parse_corsika_angles(
        path.with_name(path.name.replace(".iact", ".corsika.log"))
    )
    stem = path.name.removesuffix(".iact")
    stream = root / f"{stem}.pst"
    common.command(
        [common.tool("obdeect"), "stream", "from-corsika", str(path)],
        root / f"{stem}-from-corsika.log",
        stdout_file=stream,
    )
    events = [r["metadata"] for r in readers_records(stream) if r["type"] == "event"]
    if not events:
        raise ValueError(f"{stem}: the stream has no event record")
    event = events[0]
    stream_az_deg = math.degrees(float(event["azimuth_rad"]))
    stream_zen_deg = math.degrees(float(event["zenith_rad"]))
    _, data = readers.read_native(
        stream, keys=("cx", "cy", "cz", "w"), require_events=False
    )
    u_stream, _se_stream, rows_stream = mean_direction(
        data["cx"], data["cy"], data["cz"], data["w"]
    )
    _, iact = readers.read_iact(path)
    u_iact, se_iact, rows_iact = mean_direction(
        iact["cx"], iact["cy"], iact["cz"], iact["w"]
    )
    stream.unlink()

    travel = corsika_direction(theta, phi)
    arrival = -travel
    tol = tolerance(se_iact, CORSIKA_TOL_RAD)
    err_travel = angle_rad(u_iact, travel)
    err_arrival = angle_rad(u_iact, arrival)
    return {
        "file": path.name,
        "thetap_deg": theta,
        "phip_deg": phi,
        "stream_azimuth_deg": stream_az_deg,
        "stream_zenith_deg": stream_zen_deg,
        "azimuth_error_deg": abs(_wrap_deg(stream_az_deg - phi)),
        "zenith_error_deg": abs(stream_zen_deg - theta),
        "mean_direction_stream": u_stream.tolist(),
        "mean_direction_eventio": u_iact.tolist(),
        "se_eventio_rad": se_iact,
        "stream_vs_eventio_delta_rad": angle_rad(u_stream, u_iact),
        "rows_stream": rows_stream,
        "rows_eventio": rows_iact,
        "tolerance_rad": tol,
        "angle_to_travel_rad": err_travel,
        "angle_to_arrival_rad": err_arrival,
        "direction_error_rad": err_travel,
        "arrival_closer": err_arrival < err_travel,
        "violation": err_travel > tol,
    }


def readers_records(path):
    """Records of an obdeect stream (the Reader of the obdeect package)."""
    from obdeect.stream import Reader

    return list(Reader(str(path)))


def corsika_part(inputs: list[Path], root: Path) -> dict:
    files = [corsika_file(path, root) for path in inputs]
    return {
        "files": files,
        "azimuth_error_deg_max": max(f["azimuth_error_deg"] for f in files),
        "zenith_error_deg_max": max(f["zenith_error_deg"] for f in files),
        "stream_vs_eventio_delta_rad_max": max(
            f["stream_vs_eventio_delta_rad"] for f in files
        ),
        "direction_error_rad_max": max(f["direction_error_rad"] for f in files),
        "arrival_closer_files": sum(1 for f in files if f["arrival_closer"]),
        "direction_violations": sum(1 for f in files if f["violation"]),
    }


# ------------------------------------------------------------------ main


def version() -> str:
    return common.package_version("obdeect-sources")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="part", required=True)
    sub.add_parser("presets", help="presets at zenith 30 deg and four azimuths")
    corsika = sub.add_parser("corsika", help="CORSIKA IACT files against their logs")
    corsika.add_argument("inputs", nargs="+", type=Path, metavar="FILE.iact")
    args = parser.parse_args()
    root = Path.cwd()
    if args.part == "presets":
        part = presets_part(root)
        metrics = {
            "convention_violations": part["convention_violations"],
            "presets_disagree_90": part["presets_disagree_90"],
            "presets_disagree_all": part["presets_disagree_all"],
        }
        details = {
            "seed": SEED,
            "obdeect_sources_version": version(),
            "zenith_deg": ZENITH_DEG,
            "azimuths_deg": list(AZIMUTHS_DEG),
            "angle_tolerance_rad": ANGLE_TOL_RAD,
            "k_standard_errors": K_SE,
            "convention": "sources: travel = -astro(z, a); beams: travel = astro(z, a); "
            "astro(z, a) = (sin z cos a, -sin z sin a, cos z), a clockwise from north",
            "violations": part["violations"],
            "rows": part["rows"],
            "disagreements": part["disagreements"],
        }
    else:
        part = corsika_part(args.inputs, root)
        metrics = {
            "azimuth_error_deg_max": part["azimuth_error_deg_max"],
            "zenith_error_deg_max": part["zenith_error_deg_max"],
            "stream_vs_eventio_delta_rad_max": part["stream_vs_eventio_delta_rad_max"],
            "direction_error_rad_max": part["direction_error_rad_max"],
            "arrival_closer_files": part["arrival_closer_files"],
        }
        details = {
            "obdeect_sources_version": version(),
            "corsika_tolerance_rad": CORSIKA_TOL_RAD,
            "files": part["files"],
        }
    common.write_results(metrics, details, root)


if __name__ == "__main__":
    main()
