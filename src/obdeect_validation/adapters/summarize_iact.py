"""Summary of the photon bunches of an EventIO IACT file from a reference light-source program.

Prints events, photons (sum of bunch sizes), bunches, and the mean wavelength and the mean and
r.m.s. of the bunch impact positions at the observation level [m], over all telescopes.
Moved from obdeect/validation/sources/summarize_iact.py (which summarised per telescope).

Usage: python -m obdeect_validation.adapters.summarize_iact showers.iact.gz
"""

import sys

import numpy as np

from obdeect_validation.adapters import readers


def summarise(path):
    counts, data = readers.read_iact(path)
    summary = {
        "events": len(counts),
        "photons": float(np.sum(counts)) if counts else 0.0,
        "bunches": int(data["x"].size),
    }
    if summary["bunches"]:
        summary.update(
            wavelength_mean_nm=float(np.mean(data["wavelength"])),
            x_mean_m=float(np.mean(data["x"])),
            y_mean_m=float(np.mean(data["y"])),
            x_rms_m=float(np.std(data["x"])),
            y_rms_m=float(np.std(data["y"])),
        )
    return summary


def format_summary(summary):
    if "wavelength_mean_nm" not in summary:
        return (
            f"events {summary['events']}, photons {summary['photons']:.1f}, bunches 0"
        )
    return (
        f"events {summary['events']}, photons {summary['photons']:.1f}, "
        f"bunches {summary['bunches']}, wavelength mean {summary['wavelength_mean_nm']:.3f} nm, "
        f"impact mean ({summary['x_mean_m']:.3f}, {summary['y_mean_m']:.3f}) m, "
        f"rms ({summary['x_rms_m']:.3f}, {summary['y_rms_m']:.3f}) m"
    )


def main(path):
    print(format_summary(summarise(path)))


if __name__ == "__main__":
    main(sys.argv[1])
