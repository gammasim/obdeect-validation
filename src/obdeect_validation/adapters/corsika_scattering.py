"""Exercise to-sky and single scattering on real CORSIKA shower photons.

The scattering yield is reported as a diagnostic. This case has no independent shower-scattering
reference, so it checks the frame conversion and records output without treating the yield as a
physics agreement test.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common, curved_showers
from obdeect_validation.adapters import obdeect_cli as cli

SUMMARY = re.compile(
    r"photons in (\d+), out (\d+) \(weight ([\d.e+-]+)\), dropped (\d+)"
)
TO_SKY_SUMMARY = re.compile(r"photons in (\d+), out (\d+), dropped (\d+)")


def write_rayleigh_table(path: Path) -> Path:
    """Write a Rayleigh-only table covering CORSIKA's usual optical band."""
    wavelengths = np.arange(200, 1001)
    heights = np.geomspace(0.1, 60.0, 120)
    rows = ["# wavelength_nm altitude_km tau_mol tau_aer ssa g"]
    for wavelength in wavelengths:
        tau0 = 0.00865 * (wavelength / 1000.0) ** -4.05
        for height in heights:
            rows.append(
                f"{wavelength} {height:.8g} {tau0 * np.exp(-height / 8):.10e} 0 1 0.65"
            )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _summary(pattern, path):
    match = pattern.search(Path(path).read_text(errors="replace"))
    if not match:
        raise ValueError(f"no photon summary in {path}")
    return tuple(float(value) for value in match.groups())


def telescope_totals(path: Path) -> tuple[dict[int, int], dict[int, float]]:
    """Photon-row counts and weights keyed by telescope id."""
    from obdeect.stream import Reader

    rows, weights = {}, {}
    for record in Reader(str(path)):
        if record["type"] != "photons":
            continue
        telescope = int(record["telescope_id"])
        columns = record["columns"]
        rows[telescope] = rows.get(telescope, 0) + len(columns["w"])
        weights[telescope] = weights.get(telescope, 0.0) + float(
            np.asarray(columns["w"], dtype=np.float64).sum()
        )
    return rows, weights


def run(source: Path, root: Path) -> dict:
    from obdeect_validation.adapters.endtoend_common import rewrite

    source = source.resolve()
    ground, sky = root / "shower.pst", root / "shower-sky.pst"
    common.command(
        cli.from_corsika(source), root / "from-corsika.log", stdout_file=ground
    )
    profile = common.REFERENCE / "cfg/common/atmprof1.dat"
    flags = curved_showers.read_header_flags(source)
    use_profile = profile.is_file() and flags["refraction_bit"]
    common.command(
        cli.obdeect(
            "propagation",
            "to-sky",
            profile=profile if use_profile else None,
            i=ground,
            o=sky,
        ),
        root / "to-sky.log",
    )
    to_sky_in, to_sky_out, to_sky_dropped = _summary(
        TO_SKY_SUMMARY, root / "to-sky.log"
    )

    table = write_rayleigh_table(root / "rayleigh-scattering.dat")
    base, scattered = root / "base.pst", root / "scattered.pst"
    chunked_sky = root / "shower-sky-rechunked.pst"
    rewrite(sky, chunked_sky, split_rows=257)

    def propagation(i, o, *, scattering=False, threads=None):
        """Weighted propagation of telescopes 2 and 3, optionally with single scattering."""
        scatter = {}
        if scattering:
            scatter = {
                "scattering": table,
                "scattering_points": "32",
                "scattering_primaries": "1",
            }
        return cli.obdeect(
            "propagation",
            "run",
            transmission=common.UNIT_TRANSMISSION,
            pathlen="0",
            site_altitude="2156",
            seed="17",
            mode="weighted",
            wavelength_range=("200", "1000"),
            profile=profile if use_profile else None,
            i=i,
            **scatter,
            telescopes="2,3",
            threads=threads,
            o=o,
        )

    common.command(propagation(sky, base), root / "base.log")
    common.command(propagation(sky, scattered, scattering=True), root / "scattered.log")
    _, base_out, _, base_dropped = _summary(SUMMARY, root / "base.log")
    in_count, out_count, _, dropped = _summary(SUMMARY, root / "scattered.log")
    _, base_weights = telescope_totals(base)
    selected_rows, selected_weights = telescope_totals(scattered)
    extra_weight = sum(selected_weights.values()) - sum(base_weights.values())

    invariant = root / "scattered-rechunked.pst"
    invariant_log = root / "scattered-rechunked.log"
    common.command(
        propagation(chunked_sky, invariant, scattering=True, threads="2"),
        invariant_log,
    )
    invariant_rows, invariant_weights = telescope_totals(invariant)
    row_delta = sum(
        abs(selected_rows.get(t, 0) - invariant_rows.get(t, 0)) for t in (2, 3)
    )
    weight_delta = max(
        (
            abs(selected_weights.get(t, 0.0) - invariant_weights.get(t, 0.0))
            for t in (2, 3)
        ),
        default=0.0,
    )
    return {
        "to_sky_photons_in": int(to_sky_in),
        "to_sky_photons_out": int(to_sky_out),
        "to_sky_dropped": int(to_sky_dropped),
        "primary_photons_in": int(in_count),
        "primary_photons_out_without_scattering": int(base_out),
        "photons_out_with_scattering": int(out_count),
        "scattered_weight_diagnostic": extra_weight,
        "scattered_weight_positive_diagnostic": int(extra_weight > 0),
        "selected_output_telescope_count": len(selected_rows),
        "selected_output_has_primary_telescope": int(1 in selected_rows),
        "selection_chunk_thread_row_delta": row_delta,
        "selection_chunk_thread_max_weight_difference": weight_delta,
        "scattering_dropped": int(dropped),
        "base_dropped": int(base_dropped),
        "profile_used": int(use_profile),
    }


def _add_arguments(parser):
    cli.input_argument(parser, default=common.SHOWERS)


def analyse(root, args):
    """Scattering run of the input shower file in root: (result, details)."""
    result = run(args.input, Path(root).resolve())
    return result, {"interpretation": "scattered yield is diagnostic only"}


def main():
    common.results_main(analyse, "", add_arguments=_add_arguments)


if __name__ == "__main__":
    main()
