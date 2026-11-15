"""CORSIKA IACT inputs of the physics parity matrix: other zeniths, azimuths, protons and a South site.

The shower files of this module extend the checked-in 20 deg North gamma file (see
chain_showers.card_text and fixtures/corsika/gamma-20deg-100-3000GeV.card). Each file has the
same telescope layout: three full-size 1200 cm spheres at (0, 0), (5000, 0) and (0, 5000) cm
(the LST, MST and SST of the array origin in the existing cards), CSCAT 5 with 10000 cm, the long
bunch format is not used (plain TELFIL name, as the 20 deg card), seeds seed_base + 0, 12, 26, 44.

Per shower file: NSHOW primary showers, each thrown CSCAT 5 times (5 IACT events per shower, as
in the existing file), ESLOPE -2, no thinning. Each file has its own seed_base (no seed is shared
with another file of this module or with the curved 65 and 75 deg files, which use 17).

Sites. North: observation level 2156 m, CORSIKA ATMOSPHERE 1 (the atmprof1.dat of
sim_telarray cfg/common, as in the existing cards). South: observation level 2147 m (OBS-South
corsika_observation_level), ATMOSPHERE 99 with the South profile from
sim_telarray cfg/CTA/atmprof_ecmwf_south_winter_fixed.dat (the file the OBS-South
atmospheric_profile model parameter was generated from; CORSIKA reads it as atmprof99.dat in the
run directory). MAGNET and ARRANG are the North values of the existing cards; they are not
changed for South.

Run. The shower files are produced in the simtools-prod container, which has the flat CORSIKA
binary and the high-energy model qgs2:

  python -m obdeect_validation.adapters.shower_matrix run --name gamma-40deg-100-3000GeV \
      --root WORKDIR --data-source CORSIKA_DATA --reference SIM_TELARRAY_CHECKOUT [--cpus 3]

The run directory WORKDIR receives the card, the atmosphere file, the CORSIKA data (copied) and the
output; the shower file and its log are copied to the name the cases read, <name>.iact and
<name>.corsika.log. The log is checked before the file is accepted: the flat binary must not write
the CURVED banner, the atmosphere line must name the site profile, and the output must exist.

Checks. check_file reads a shower file with eventio and compares each event header with the
card: primary particle, energy window, zenith and azimuth windows (the THETAP and PHIP ranges),
observation level, and the number of distinct showers. Limits are the card values with a
float32 tolerance of 1e-3 (header words are single precision).
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from obdeect_validation.adapters import chain_showers, common

IMAGE = "ghcr.io/gammasim/simtools-prod:v0.39.0-v78010-v2025-11-30-rc-generic"
FLAT_EXECUTABLE = "/workdir/simulation_software/corsika7/corsika_qgs2_urqmd_flat"
CURVED_PATTERN = "CURVED VERSION WITH SLIDING PLANAR ATMOSPHERE"
PROFILE_DIRECTORY = "cfg/common"
SOUTH_PROFILE = "cfg/CTA/atmprof_ecmwf_south_winter_fixed.dat"
HEADER_TOLERANCE = 1e-3  # float32 header words (degrees, GeV, cm)
EVENTS_PER_SHOWER = (
    5  # CSCAT 5: each shower is reused five times, as in the 20 deg file
)
SEED_OFFSETS = chain_showers.SEED_OFFSETS
CORSIKA_DATA_FILES = (
    "NUCNUCCS",
    "UrQMD-1.3.1-xs.dat",
    "QGSDAT01",
    "qgsdat-II-04",
    "sectnu-II-04",
)
CORSIKA_DATA_RENAME = {
    "qgsdat-II-04": "QGSDAT-II-04",
    "sectnu-II-04": "SECTNU-II-04",
}


@dataclass(frozen=True)
class Site:
    """Observation level (cm), CORSIKA atmosphere model, profile source and run-directory name."""

    name: str
    obslev_cm: float
    atmosphere: int
    profile: str
    profile_name: str


SITES = {
    "North": Site(
        "North", 215600.0, 1, f"{PROFILE_DIRECTORY}/atmprof1.dat", "atmprof1.dat"
    ),
    "South": Site("South", 214700.0, 99, SOUTH_PROFILE, "atmprof99.dat"),
}


@dataclass(frozen=True)
class Shower:
    """One shower file: primary, zenith and azimuth windows (deg), energy window (GeV)."""

    name: str
    primary: int  # CORSIKA particle id: 1 gamma, 14 proton
    zenith: tuple[float, float]  # THETAP min max
    azimuth: tuple[float, float]  # PHIP min max
    energy: tuple[float, float]  # ERANGE min max
    site: str
    seed_base: int
    showers: int = 20
    telescopes: tuple = chain_showers.FLAT_TELESCOPES  # (x, y, z, radius) cm
    scatter_cm: float = chain_showers.FLAT_SCATTER_CM  # CSCAT radius (cm)

    @property
    def particle(self) -> str:
        return {1: "gamma", 14: "proton"}[self.primary]


# Seeds: each file has its own seed_base (the 20 deg North gamma file uses 17).
PROTON_ZENITH = (19.0, 21.0)

# North array of the simulation models (OBS-North array_layouts "CTAO-North-4-LSTs-1-MSTs" elements
# LSTN-01..04 and MSTN-03, plus MSTN-15 of "1mst"). Positions are array_element_position_ground
# (m, relative to the array reference point) written in cm as simtools writes the TELESCOPE lines
# (get_corsika_telescope_list). Sphere radii: LST 1300 cm (mirror reach 1229 cm, LSTN-01 mirror_list
# 2.0.0; the LST design value is 1250 cm); MST 915 cm (telescope_sphere_radius of MSTx-NectarCam;
# mirror reach 681 cm). The MSTs are 210 m (MSTN-03) and 0 m (MSTN-15) from the reference point.
ARRAY_NORTH_TELESCOPES = (
    (-7091.0, -5235.0, 4500.0, 1300.0),  # LSTN-01
    (-3528.0, 6592.0, 3250.0, 1300.0),  # LSTN-02
    (7531.0, 5021.0, 2820.0, 1300.0),  # LSTN-03
    (3088.0, -6378.0, 3280.0, 1300.0),  # LSTN-04
    (2686.0, 20854.0, 1140.0, 915.0),  # MSTN-03
    (464.0, -71.0, 2950.0, 915.0),  # MSTN-15
)
# CSCAT radius (cm) covering the array: cores within 250 m, the farthest MST at 210 m.
ARRAY_NORTH_SCATTER_CM = 25000.0

SPECS = {
    shower.name: shower
    for shower in (
        # The checked-in 20 deg North file (same card, seeds 17); listed for the input check only.
        Shower(
            "gamma-20deg-100-3000GeV", 1, (20, 20), (0, 0), (100, 3000), "North", 17
        ),
        # Single LST at the origin with a 1300 cm sphere (the 1200 cm sphere is too small for the
        # 1229 cm mirror reach); CSCAT 10000 cm as the 20 deg card.
        Shower(
            "gamma-20deg-100-3000GeV-r1300",
            1,
            (20, 20),
            (0, 0),
            (100, 3000),
            "North",
            717,
            telescopes=((0.0, 0.0, 0.0, 1300.0),),
        ),
        # Four LSTs and two MSTs of the North array at their real positions, CSCAT over the array
        # area, so several telescopes see each event (20 showers x 5 scatter reuses).
        Shower(
            "gamma-20deg-100-3000GeV-array-north",
            1,
            (20, 20),
            (0, 0),
            (100, 3000),
            "North",
            617,
            telescopes=ARRAY_NORTH_TELESCOPES,
            scatter_cm=ARRAY_NORTH_SCATTER_CM,
        ),
        Shower(
            "gamma-40deg-100-3000GeV", 1, (40, 40), (0, 0), (100, 3000), "North", 117
        ),
        Shower(
            "gamma-60deg-300-5000GeV", 1, (60, 60), (0, 0), (300, 5000), "North", 217
        ),
        Shower(
            "gamma-20deg-az180-100-3000GeV",
            1,
            (20, 20),
            (180, 180),
            (100, 3000),
            "North",
            317,
        ),
        # 45 showers: a proton file is about 40 MB per shower (5 showers: 199 MB), so 200 showers
        # would be about 8 GB; 45 keeps the file below 2 GB.
        Shower(
            "proton-20deg-300-10000GeV",
            14,
            PROTON_ZENITH,
            (0, 0),
            (300, 10000),
            "North",
            417,
            showers=45,
        ),
        Shower(
            "gamma-20deg-South-100-3000GeV",
            1,
            (20, 20),
            (0, 0),
            (100, 3000),
            "South",
            517,
        ),
    )
}


def _number(value: float) -> str:
    return f"{value:g}"


def card_text(shower: Shower, data_dir: str = "./data/", direct: str = "./") -> str:
    """CORSIKA input card of one shower file, in the style of the 20 deg North card.

    The telescope lines, CSCAT and the fixed flags are chain_showers.card_text's; only the primary,
    energy, angle, site and seed lines differ. TELFIL is the plain name (no '+' prefix).
    """
    site = SITES[shower.site]
    seeds = "".join(
        f"SEED {shower.seed_base + offset} 0 0\n" for offset in SEED_OFFSETS
    )
    scatter = "".join(
        f"TELESCOPE {x:g} {y:g} {z:g} {r:g}\n" for x, y, z, r in shower.telescopes
    )
    return f"""RUNNR 1
EVTNR 1
NSHOW {shower.showers}
PRMPAR {shower.primary}
ERANGE {_number(shower.energy[0])} {_number(shower.energy[1])}
ESLOPE -2.0
THETAP {_number(shower.zenith[0])} {_number(shower.zenith[1])}
PHIP {_number(shower.azimuth[0])} {_number(shower.azimuth[1])}
{seeds}OBSLEV {site.obslev_cm:g}
ATMOSPHERE {site.atmosphere} Y
MAGNET 20.0 40.0
ARRANG 0.0
HADFLG 0 0 0 0 0 2
ECUTS 0.3 0.3 0.02 0.02
ELMFLG F F
CSCAT 5 {shower.scatter_cm:g} 0.0
DATDIR {data_dir}
DIRECT {direct}
{scatter}TELFIL {shower.name}.iact
IACT SPLIT_AUTO 0
IACT IO_BUFFER 800MB
IACT MAX_BUNCHES 1000000
EXIT
"""


def stage(root: Path, shower: Shower, data_source: Path, reference: Path) -> Path:
    """Write the card, copy the CORSIKA data and the site atmosphere into root; return the card.

    data_source holds the CORSIKA data files (NUCNUCCS, QGSDAT01, ...); reference is a
    sim_telarray checkout (cfg/common/atmprof1.dat, cfg/CTA/atmprof_ecmwf_south_winter_fixed.dat).
    """
    site = SITES[shower.site]
    data = root / "data"
    data.mkdir(parents=True, exist_ok=True)
    originals = [data_source / name for name in CORSIKA_DATA_FILES]
    originals += sorted(data_source.glob("EGSDAT6_*"))
    for original in originals:
        if not original.is_file():
            raise FileNotFoundError(original)
        target = data / CORSIKA_DATA_RENAME.get(original.name, original.name)
        if not target.exists():
            shutil.copyfile(original, target)
    shutil.copyfile(reference / site.profile, root / site.profile_name)
    card = root / f"{shower.name}.card"
    card.write_text(card_text(shower))
    return card


def check_log(text: str, shower: Shower) -> None:
    """Raise ValueError unless the CORSIKA log is a flat run with the site atmosphere."""
    site = SITES[shower.site]
    if CURVED_PATTERN in text:
        raise ValueError("log has the CURVED banner; the flat binary was expected")
    if f"Atmospheric profile {site.atmosphere} with" not in text or (
        f"read from file {site.profile_name}" not in text
    ):
        raise ValueError(
            f"log lacks the atmosphere line of {site.name} ({site.profile_name})"
        )


def run(
    root: Path,
    shower: Shower,
    data_source: Path,
    reference: Path,
    image: str = IMAGE,
    cpus: int = 3,
    memory: str = "4g",
) -> Path:
    """Run CORSIKA (flat, qgs2) in the container for one shower file; return root/<name>.iact."""
    root = Path(root).resolve()
    stage(root, shower, data_source.resolve(), reference.resolve())
    card = f"{shower.name}.card"
    script = f"{FLAT_EXECUTABLE} < {card} > corsika.log 2>&1"
    result = subprocess.run(
        [
            "podman",
            "run",
            "--rm",
            "--cpus",
            str(cpus),
            "--memory",
            memory,
            "--volume",
            f"{root}:/work:Z",
            "--workdir",
            "/work",
            "--entrypoint",
            "/bin/bash",
            image,
            "-c",
            script,
        ],
        check=False,
    )
    log = root / "corsika.log"
    if result.returncode:
        raise RuntimeError(f"CORSIKA container exited {result.returncode}; see {log}")
    text = log.read_text(errors="replace")
    check_log(text, shower)
    iact = root / f"{shower.name}.iact"
    if not iact.is_file():
        raise RuntimeError(f"CORSIKA wrote no {iact}; see {log}")
    shutil.copyfile(log, root / f"{shower.name}.corsika.log")
    return iact


def header_mismatches(
    headers, shower: Shower, tolerance: float = HEADER_TOLERANCE
) -> dict:
    """Count event headers that disagree with the card of shower.

    headers: EVTH tuples of the IACT events (element i is word i + 1 of the EVTH block): word 2
    shower number, word 3 particle id, word 4 energy (GeV), word 11 zenith and word 12 azimuth
    (rad), word 19 the observation levels (cm).
    """
    site = SITES[shower.site]
    zen_lo, zen_hi = shower.zenith
    az_centre = sum(shower.azimuth) / 2
    az_half = (shower.azimuth[1] - shower.azimuth[0]) / 2
    e_lo, e_hi = shower.energy
    counts = {
        "particle": 0,
        "energy": 0,
        "zenith": 0,
        "azimuth": 0,
        "observation_level": 0,
    }
    for header in headers:
        if round(float(header[2])) != shower.primary:
            counts["particle"] += 1
        energy = float(header[3])
        if not (e_lo - tolerance * e_lo <= energy <= e_hi + tolerance * e_hi):
            counts["energy"] += 1
        theta = math.degrees(float(header[10]))
        if not (zen_lo - tolerance <= theta <= zen_hi + tolerance):
            counts["zenith"] += 1
        # Azimuth difference on the circle, in degrees, against the centre of the PHIP window.
        phi = math.degrees(float(header[11]))
        delta = (phi - az_centre + 180.0) % 360.0 - 180.0
        if abs(delta) > az_half + tolerance:
            counts["azimuth"] += 1
        levels = [float(v) for v in header[18] if float(v) > 0]
        if not levels or abs(levels[0] - site.obslev_cm) > 1.0:
            counts["observation_level"] += 1
    return counts


def check_file(path: Path, shower: Shower) -> dict:
    """Read a shower file with eventio and compare its headers and counts with shower's card."""
    from eventio.iact import IACTFile

    headers, telescopes, shower_numbers = [], [], set()
    events_with_two = 0
    with IACTFile(str(path)) as stream:
        for event in stream:
            headers.append(event.header)
            shower_numbers.add(round(float(event.header[1])))
            hit = len(event.photon_bunches)
            telescopes.append(hit)
            events_with_two += int(hit >= 2)
    mismatches = header_mismatches(headers, shower)
    return {
        "readable": 1,
        "events": len(headers),
        "distinct_showers": len(shower_numbers),
        "events_per_shower_expected": EVENTS_PER_SHOWER,
        "shower_count_delta": len(shower_numbers) - shower.showers,
        "event_count_delta": len(headers) - EVENTS_PER_SHOWER * shower.showers,
        "events_with_two_or_more_telescopes": events_with_two,
        "events_with_telescope_hits": sum(1 for n in telescopes if n),
        **{f"header_{key}_mismatches": value for key, value in mismatches.items()},
    }


def _log_mismatch(path: Path, shower: Shower) -> int:
    """1 when the <name>.corsika.log next to a shower file is missing or fails check_log."""
    try:
        check_log(path.with_suffix(".corsika.log").read_text(errors="replace"), shower)
    except (OSError, ValueError):
        return 1
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="shower_matrix")
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="print the shower files of the matrix")
    listing.add_argument("--json", action="store_true")
    cards = commands.add_parser("card", help="print the card of one shower file")
    cards.add_argument("name", choices=sorted(SPECS))
    gen = commands.add_parser(
        "run", help="run CORSIKA in the container for one shower file"
    )
    gen.add_argument("--name", choices=sorted(SPECS), required=True)
    gen.add_argument("--root", type=Path, required=True)
    gen.add_argument("--data-source", type=Path, required=True)
    gen.add_argument("--reference", type=Path, required=True)
    gen.add_argument("--image", default=IMAGE)
    gen.add_argument("--cpus", type=int, default=3)
    check = commands.add_parser("check", help="check shower files against their cards")
    check.add_argument("paths", nargs="+", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "list":
            rows = {name: vars(shower) for name, shower in SPECS.items()}
            print(
                json.dumps(rows, indent=2, default=list)
                if args.json
                else "\n".join(rows)
            )
            return 0
        if args.command == "card":
            print(card_text(SPECS[args.name]), end="")
            return 0
        if args.command == "run":
            iact = run(
                args.root,
                SPECS[args.name],
                args.data_source,
                args.reference,
                image=args.image,
                cpus=args.cpus,
            )
            print(iact)
            return 0
        metrics, details = {}, {}
        for path in args.paths:
            name = path.name.removesuffix(".iact")
            result = check_file(path, SPECS[name])
            result["log_mismatch"] = _log_mismatch(path, SPECS[name])
            details[name] = result
            key = name.replace("-", "_")
            metrics[f"{key}_readable"] = result["readable"]
            metrics[f"{key}_shower_count_delta"] = result["shower_count_delta"]
            metrics[f"{key}_event_count_delta"] = result["event_count_delta"]
            metrics[f"{key}_header_mismatches"] = sum(
                value for field, value in result.items() if field.startswith("header_")
            )
            metrics[f"{key}_log_mismatch"] = result["log_mismatch"]
        common.write_results(metrics, details)
        return 0
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        print(f"shower_matrix: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
