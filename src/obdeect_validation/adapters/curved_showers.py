"""CORSIKA showers for the curvature cases: flat and CURVED atmosphere, cards, flags and logs.

Geometry. A CORSIKA binary is either the standard (flat) program or the CURVED version (sliding
planar atmosphere, meant for large zenith angles). simtools picks the executable from
$SIMTOOLS_CORSIKA_PATH/build_opts.yml by atmosphere geometry and uses the curved one from
65 deg zenith (simtools defaults.py CURVED_ATMOSPHERE_MIN_ZENITH_ANGLE_DEG). simtools does not
write TSTART for the curved binary (corsika_config.py, _corsika_configuration_interaction_flags).
A curved run writes CURVED_LOG_PATTERN into its log.

Header flags (CORSIKA EVTH, 1-based Fortran words; the file stores them 0-based):
  word 79  CURVED flag: 0 standard, 2 curved (corsika.tex, EVTH table)
  word 77  integer flags; bit 0x40 = CURVED option compiled in, bit 0x10 = refraction
           (ATMEXT refraction, corsika.tex footnote). sim_telarray reads the same word.

Compact bunch format limits (CORSIKA bernlohr/README, section d). The IACT bunches are written
in the compact format (16 bytes per bunch) unless TELFIL is prefixed with '+'. Two limits matter:
  (2) for a detector sphere of radius R, sec(z) * R < 32.7 m (impact points are taken in the
      horizontal plane through the sphere centre);
  (3) times within 3.27 us of the reference time, and for large zenith angles horizontal core
      distances of about 1000 m.
The cards therefore follow these rules:
  * The sphere radius, not the core scatter, sets limit (2): R = min(1200 cm, 0.8 * 32.7 m *
    cos z) with a 20 % margin. At 65 deg R = 11.05 m, at 75 deg 6.76 m, at 80 deg 4.54 m and at
    85 deg 2.27 m. Smaller spheres collect fewer photons; this is the price of the compact format
    without '+'.
  * The telescopes sit at (0, 0), (300 m, 0) and (0, 300 m) with Z = R, the CORSIKA
    recommendation that a sphere is not cut by the observation level (corsika.tex, TELESCOPE).
  * CSCAT scatters the core within 300 m, so every core-telescope distance is below 600 m and
    limit (3) holds. CSCAT and TELESCOPE values are in cm in CORSIKA; the constants below are cm.
    (The 800 quoted for CSCAT in the brief would be 8 m in CORSIKA units.)
The flat cards of this module use the same telescope layout as the curved ones (only TSTART and
the executable differ). generate() in chain_showers keeps the legacy flat layout.

Cards. NSHOW 20 showers per file, PRMPAR 1 (gamma), OBSLEV 215600 cm, ATMOSPHERE 1 Y (refraction on),
fixed seeds 17, 29, 43, 61 (SEED lines, seed_base offsets). Energy windows: the threshold rises
with zenith (longer slant path, more absorption, wider light pool), so 300 to 5000 GeV is used at
65 deg and 1 to 10 TeV at 75, 80 and 85 deg. These windows are the brief's choice; whether the
telescopes trigger in each window is to be checked after generation (no CORSIKA run was possible
here). Placeholder DATDIR ("<corsika data>") in the checked-in cards; run_corsika substitutes the
real data directory.

Units and frames: positions in cm (CORSIKA), returned bunch positions in m, times in ns, wavelengths
in nm, zem in m above sea level (with CURVED the emission altitude is the altitude a.s.l.).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from eventio import EventIOFile
from eventio.iact import IACTFile

from obdeect_validation.adapters import chain_showers

CURVED_LOG_PATTERN = "CURVED VERSION WITH SLIDING PLANAR ATMOSPHERE"
GEOMETRIES = ("flat", "curved")
DEFAULT_CORSIKA_PATH = "/workdir/simulation_software/corsika7"
PREFERRED_HE_MODELS = ("qgs2", "qgs3", "epos")
DATA_PLACEHOLDER = "<corsika data>"
FIXTURE_CASES = (
    # (geometry, zenith_deg, energy_min_GeV, energy_max_GeV): the checked-in cards
    ("curved", 65, 300, 5000),
    ("flat", 65, 300, 5000),
    ("curved", 75, 1000, 10000),
    ("curved", 80, 1000, 10000),
    ("curved", 85, 1000, 10000),
)

EVTH_TYPE = 1202
EVTH_WORDS = 273
EVTH_FLAGS_INDEX = 76  # word 77
EVTH_CURVED_INDEX = 78  # word 79
IACT_CURVED_BIT = 0x40
REFRACTION_BIT = 0x10

SPHERE_LIMIT_M = 32.7  # compact format limit (2)
SPHERE_MARGIN = 0.8
MAX_RADIUS_CM = 1200.0  # radius of the legacy flat layout
CURVED_TELESCOPES_XY_CM = ((0.0, 0.0), (30000.0, 0.0), (0.0, 30000.0))
CURVED_SCATTER_CM = 30000.0


class MissingPrerequisite(RuntimeError):
    """A CORSIKA build, executable or input this case needs is not installed (exit code 3)."""


def _geometry(geometry: str) -> str:
    value = str(geometry).lower()
    if value not in GEOMETRIES:
        raise ValueError(
            f"unknown CORSIKA geometry {geometry!r}; expected flat or curved"
        )
    return value


def log_is_curved(log_text: str) -> bool:
    """True if a CORSIKA log was written by the CURVED version."""
    return CURVED_LOG_PATTERN in log_text


def shower_name(zenith_deg, geometry, energy_min=300, energy_max=5000) -> str:
    """File name of the IACT output, e.g. gamma-65deg-300-5000GeV-curved.iact."""
    return f"gamma-{zenith_deg:g}deg-{energy_min:g}-{energy_max:g}GeV-{_geometry(geometry)}.iact"


def telescope_radius_cm(zenith_deg) -> float:
    """Sphere radius in cm that keeps sec(z) * R below the compact-format limit (2)."""
    limit_cm = (
        SPHERE_MARGIN * SPHERE_LIMIT_M * 100.0 * math.cos(math.radians(zenith_deg))
    )
    return min(MAX_RADIUS_CM, float(math.floor(limit_cm)))


def layout(zenith_deg) -> dict:
    """Telescope and scatter keywords of the compact-format cards (cm), for chain_showers.card_text."""
    radius = telescope_radius_cm(zenith_deg)
    if radius / 100.0 / math.cos(math.radians(zenith_deg)) >= SPHERE_LIMIT_M:
        raise ValueError(
            f"no compact-format telescope radius at zenith {zenith_deg} deg"
        )
    return {
        "telescopes": tuple((x, y, radius, radius) for x, y in CURVED_TELESCOPES_XY_CM),
        "scatter_cm": CURVED_SCATTER_CM,
    }


def decode_flags(word77: float, word79: float) -> dict:
    """CURVED and refraction flags from EVTH word 77 and word 79 (floats as stored)."""
    flags = round(float(word77))
    curved = round(float(word79))
    return {
        "evth_word77_int": flags,
        "evth_word79": curved,
        "iact_curved_bit": bool(flags & IACT_CURVED_BIT),
        "corsika_curved_flag": curved == 2,
        "refraction_bit": bool(flags & REFRACTION_BIT),
    }


def _evth_words(iact_path) -> np.ndarray:
    """The 273 float words of the first EVTH block (word n is element n - 1)."""
    with EventIOFile(str(iact_path)) as file:
        for obj in file:
            if obj.header.type != EVTH_TYPE:
                continue
            obj.seek(0)
            data = obj.read()
            count = int(np.frombuffer(data[:4], dtype="<i4")[0])
            if count != EVTH_WORDS:
                raise ValueError(
                    f"EVTH block of {iact_path} has {count} words, expected 273"
                )
            return np.frombuffer(data[4 : 4 + 4 * count], dtype="<f4").astype(
                np.float64
            )
    raise ValueError(f"no EVTH block in {iact_path}")


def read_header_flags(iact_path) -> dict:
    """Curvature and refraction flags of the first shower of an IACT file (see decode_flags)."""
    words = _evth_words(iact_path)
    return decode_flags(words[EVTH_FLAGS_INDEX], words[EVTH_CURVED_INDEX])


def read_bunches(iact_path):
    """Yield (event_number, telescope_id, columns) for each event and telescope.

    Columns (numpy arrays): x, y, zem in m; cx, cy; time in ns (including the event time offset);
    photons (bunch size); wavelength in nm. Positions are converted from CORSIKA cm.
    """
    with IACTFile(str(iact_path)) as stream:
        for number, event in enumerate(stream):
            for telescope, bunches in event.photon_bunches.items():
                names = bunches.dtype.names or ()
                columns = {}
                for key in ("x", "y", "zem"):
                    if key in names:
                        columns[key] = np.asarray(bunches[key], dtype=np.float64) * 0.01
                for key in ("cx", "cy", "photons", "wavelength"):
                    if key in names:
                        columns[key] = np.asarray(bunches[key], dtype=np.float64)
                if "time" in names:
                    columns["time"] = np.asarray(
                        bunches["time"], dtype=np.float64
                    ) + float(event.time_offset)
                yield number, int(telescope), columns


def _manifest_variants(manifest: Path) -> list[dict]:
    """Variant entries of build_opts.yml: PyYAML when installed, else the simtools layout."""
    text = manifest.read_text(encoding="utf-8")
    try:
        import yaml
    except ImportError:
        return _variants_without_yaml(text)
    data = yaml.safe_load(text)
    variants = data.get("variant") if isinstance(data, dict) else None
    return [entry for entry in variants or [] if isinstance(entry, dict)]


def _variants_without_yaml(text: str) -> list[dict]:
    """Parse the 'variant:' list of build_opts.yml with the standard library.

    Handles what simtools writes: a top-level 'variant:' key, then '  - key: value' items with
    further 'key: value' lines indented by four spaces. Other top-level keys end the list.
    """
    variants, current, in_list = [], None, False
    for raw in text.splitlines():
        line = raw.split(" #", 1)[0].rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip())
        stripped = line.strip()
        if indent == 0:
            in_list, current = stripped == "variant:", None
            continue
        if not in_list:
            continue
        if stripped.startswith("- "):
            current = {}
            variants.append(current)
            stripped = stripped[2:].strip()
        if current is None or ":" not in stripped:
            continue
        key, _, value = stripped.partition(":")
        current[key.strip()] = value.strip().strip("\"'")
    return variants


def corsika_executable(geometry: str, he_model: str | None = None) -> Path:
    """Installed CORSIKA executable for an atmosphere geometry (flat or curved).

    The root is $SIMTOOLS_CORSIKA_PATH (default /workdir/simulation_software/corsika7) and the
    variants come from its build_opts.yml. With he_model the variant must use that high-energy
    model; otherwise CORSIKA_HE_MODEL (if set) is tried, then qgs2, qgs3 and epos.

    Raises MissingPrerequisite when the manifest, a matching variant or its executable is absent.
    """
    geometry = _geometry(geometry)
    root = Path(os.environ.get("SIMTOOLS_CORSIKA_PATH") or DEFAULT_CORSIKA_PATH)
    manifest = root / "build_opts.yml"
    if not manifest.is_file():
        raise MissingPrerequisite(f"CORSIKA build manifest not found: {manifest}")
    variants = [
        entry
        for entry in _manifest_variants(manifest)
        if str(entry.get("atmosphere_geometry", "")).lower() == geometry
    ]
    if he_model:
        wanted = [he_model]
    else:
        wanted = [os.environ.get("CORSIKA_HE_MODEL", ""), *PREFERRED_HE_MODELS]
        wanted = [model for model in wanted if model]
    for model in wanted:
        for entry in variants:
            if str(entry.get("he_hadronic_model", "")).lower() != model.lower():
                continue
            executable = root / str(entry.get("executable", ""))
            if not (executable.is_file() and os.access(executable, os.X_OK)):
                raise MissingPrerequisite(
                    f"CORSIKA variant is not an executable: {executable}"
                )
            return executable
    available = sorted(
        f"{entry.get('he_hadronic_model')}/{entry.get('atmosphere_geometry')}"
        for entry in _manifest_variants(manifest)
    )
    raise MissingPrerequisite(
        f"no installed {geometry} CORSIKA variant for HE model {', '.join(wanted)} in "
        f"{manifest}; available: {', '.join(available) or 'none'}"
    )


def write_card(
    root,
    geometry,
    zenith_deg,
    *,
    showers=20,
    energy_min=300,
    energy_max=5000,
    seed_base=17,
    data_dir=DATA_PLACEHOLDER,
    direct="./",
) -> Path:
    """Write the CORSIKA card of one shower file as root/<shower stem>.card and return it.

    TSTART T is written for flat geometry only (simtools omits it for the curved binary).
    """
    geometry = _geometry(geometry)
    name = shower_name(
        zenith_deg, geometry, energy_min=energy_min, energy_max=energy_max
    )
    text = chain_showers.card_text(
        showers=showers,
        zenith=zenith_deg,
        energy_min=energy_min,
        energy_max=energy_max,
        data_dir=data_dir,
        direct=direct,
        telfil=name,
        seed_base=seed_base,
        tstart=geometry == "flat",
        **layout(zenith_deg),
    )
    path = Path(root) / f"{name.removesuffix('.iact')}.card"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def _execute(executable, card, log, cwd) -> None:
    """Run CORSIKA with the card on stdin; stdout and stderr go to log."""
    with Path(card).open("rb") as stdin, Path(log).open("wb") as output:
        result = subprocess.run(
            [str(executable)],
            cwd=cwd,
            stdin=stdin,
            stdout=output,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode:
        raise RuntimeError(f"{executable} exited {result.returncode}; see {log}")


def run_corsika(
    root,
    geometry,
    zenith_deg,
    *,
    showers=20,
    energy_min=300,
    energy_max=5000,
    seed_base=17,
    runner=_execute,
) -> Path:
    """Run CORSIKA for one shower file in root and return the .iact path.

    The log is root/corsika.log and is copied next to the shower file as <stem>.corsika.log
    (the name the stream case reads). A curved run must contain CURVED_LOG_PATTERN and a flat
    run must not, so a flat binary cannot pass as curved.
    """
    geometry = _geometry(geometry)
    root = Path(root).resolve()
    executable = corsika_executable(geometry)
    hashes = chain_showers.prepare_data(root)
    write_card(
        root,
        geometry,
        zenith_deg,
        showers=showers,
        energy_min=energy_min,
        energy_max=energy_max,
        seed_base=seed_base,
        data_dir=root / "data",
        direct=f"{root}/",
    )
    name = shower_name(
        zenith_deg, geometry, energy_min=energy_min, energy_max=energy_max
    )
    card = root / f"{name.removesuffix('.iact')}.card"
    iact = root / name
    log = root / "corsika.log"
    iact.unlink(missing_ok=True)
    runner(executable, card, log, root)
    text = log.read_text(errors="replace") if log.is_file() else ""
    curved = log_is_curved(text)
    if curved != (geometry == "curved"):
        raise RuntimeError(
            f"{executable} {'lacks' if geometry == 'curved' else 'has'} the CURVED log pattern "
            f"for geometry {geometry}; see {log}"
        )
    if not iact.is_file():
        raise RuntimeError(f"CORSIKA wrote no {iact}; see {log}")
    shutil.copyfile(log, iact.with_suffix(".corsika.log"))
    (root / "input-hashes.json").write_text(
        json.dumps({"executable": str(executable), "data": hashes}, indent=2) + "\n"
    )
    return iact


def main(argv=None) -> int:
    """Write the checked-in cards (cards) or run one shower file (run).

    Exit codes follow the validation harness: 0 ok, 2 error, 3 missing prerequisite.
    """
    parser = argparse.ArgumentParser(prog="curved_showers")
    commands = parser.add_subparsers(dest="command", required=True)
    cards = commands.add_parser(
        "cards", help="write the placeholder cards of the fixtures"
    )
    cards.add_argument("--output", type=Path, required=True)
    shower = commands.add_parser("run", help="run CORSIKA for one shower file")
    shower.add_argument("--geometry", choices=GEOMETRIES, required=True)
    shower.add_argument("--zenith", type=float, required=True)
    shower.add_argument("--energy-min", type=float, required=True)
    shower.add_argument("--energy-max", type=float, required=True)
    shower.add_argument("--showers", type=int, default=20)
    shower.add_argument("--seed-base", type=int, default=17)
    shower.add_argument("--root", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "cards":
            for geometry, zenith, low, high in FIXTURE_CASES:
                write_card(
                    args.output, geometry, zenith, energy_min=low, energy_max=high
                )
            return 0
        iact = run_corsika(
            args.root,
            args.geometry,
            args.zenith,
            showers=args.showers,
            energy_min=args.energy_min,
            energy_max=args.energy_max,
            seed_base=args.seed_base,
        )
        print(iact)
        return 0
    except MissingPrerequisite as exc:
        print(f"missing prerequisite: {exc}", file=sys.stderr)
        return 3
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"curved_showers: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
