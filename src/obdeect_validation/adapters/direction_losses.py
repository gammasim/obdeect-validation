"""Loss accounting of a mixed-direction population in plane mode (phase 2, upward photons).

Case propagation.direction-losses. Eight categories of N photons each, one telescope, the observation
level 2 km, telescope plane 500 m above it, profile given (refraction_applied = 0), extinction off:

    category                    z [m]   direction                          expected
    downward_above_plane        800     downward                             reaches
    downward_below_plane        300     downward                             moving_away
    upward_below_plane          100     upward                               reaches
    upward_above_plane          700     upward                               moving_away
    horizontal                  200     cz = 0                               horizontal
    below_ground_downward       -10     downward                             below_ground
    below_ground_upward         -10     upward                               below_ground
    upward_turnover             100     upward, sin(theta) past the turn     turns_over

The expected reason of each category is written here, independently of the code. Checks: the exact count
per reason and in total, photons_in = photons_out + dropped, dropped = sum of dropped_by_reason, and the
output rows of the reaching categories (zem identifies them) are exactly the expected ones.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from obdeect_validation.adapters import common
from obdeect_validation.adapters import upward_common as uc

N = 25
PLANE_Z_M = 500.0
REASONS = (
    "moving_away",
    "horizontal",
    "below_ground",
    "turns_over",
    "not_bendable",
    "misses_level",
)
# category: (z, cx, cy, cz sign, expected reason or None for reaching)
CATEGORIES = {
    "downward_above_plane": (800.0, "down", None),
    "downward_below_plane": (300.0, "down", "moving_away"),
    "upward_below_plane": (100.0, "up", None),
    "upward_above_plane": (700.0, "up", "moving_away"),
    "horizontal": (200.0, "horizontal", "horizontal"),
    "below_ground_downward": (-10.0, "down", "below_ground"),
    "below_ground_upward": (-10.0, "up", "below_ground"),
    "upward_turnover": (100.0, "turnover", "turns_over"),
}
ZENITH_DEG = 40.0
AZIMUTH_DEG = 70.0


def version() -> str:
    return common.package_version("obdeect-propagation")


def population():
    """Columns of all categories (float32), the category of each row, and the expected reason."""
    rng = np.random.default_rng(20261126)
    cols = {k: [] for k in ("x", "y", "z", "cx", "cy", "cz", "t")}
    names, expected = [], []
    prof = uc.profile()
    turn_sin = uc.critical_sin(prof, uc.SITE_M + 100.0, uc.SITE_M + PLANE_Z_M)
    turn_sin = turn_sin + 0.5 * (
        1.0 - turn_sin
    )  # below 1: the plane is reached with sin > 1
    for name, (z, kind, reason) in CATEGORIES.items():
        for _ in range(N):
            if kind == "horizontal":
                cx, cy, cz = 1.0, 0.0, 0.0
            elif kind == "turnover":
                az = math.radians(rng.uniform(0.0, 360.0))
                cx, cy, cz = (
                    turn_sin * math.cos(az),
                    turn_sin * math.sin(az),
                    math.sqrt(1 - turn_sin**2),
                )
            else:
                th = math.radians(ZENITH_DEG)
                az = math.radians(AZIMUTH_DEG + rng.uniform(-30, 30))
                cx, cy = math.sin(th) * math.cos(az), math.sin(th) * math.sin(az)
                cz = math.cos(th) if kind == "up" else -math.cos(th)
            values = {
                "x": rng.uniform(-200, 200),
                "y": rng.uniform(-200, 200),
                "z": z,
                "cx": cx,
                "cy": cy,
                "cz": cz,
                "t": rng.uniform(0, 5),
            }
            for k, v in values.items():
                cols[k].append(v)
            names.append(name)
            expected.append(reason)
    arrays = {k: np.asarray(v, dtype=np.float32) for k, v in cols.items()}
    arrays["w"] = np.ones(len(names), dtype=np.float32)
    arrays["wavelength"] = np.full(len(names), 400.0, dtype=np.float32)
    return arrays, np.asarray(names), expected, turn_sin


def analyse(root: Path) -> tuple[dict, dict]:
    from obdeect.propagation import propagate

    cols, names, expected, turn_sin = population()
    stream = uc.write_stream(
        root / "direction-losses-in.pst",
        [(0.0, 0.0, PLANE_Z_M)],
        [(1, cols)],
        refraction_applied=0,
        chunk_rows=64,
    )
    summary = propagate(
        stream,
        root / "direction-losses-out.pst",
        extinction=False,
        site_altitude_m=uc.SITE_M,
        seed=1,
        profile=uc.PROFILE,
    )
    # Expected counts from the categories.
    expected_reason = {r: 0 for r in REASONS}
    for reason in expected:
        if reason is not None:
            expected_reason[reason] += 1
    expected_dropped = sum(expected_reason.values())
    expected_out = len(names) - expected_dropped
    reason_errors = sum(
        abs(summary["dropped_by_reason"][r] - expected_reason[r]) for r in REASONS
    )
    # Reaching rows: zem identifies them (above the observation level by z; the reaching ones: 2100 m,
    # 2800 m). Count the categories in the output.
    from obdeect_validation.adapters import endtoend_common as e2e

    rows = e2e.photon_rows(root / "direction-losses-out.pst")
    zem = rows["zem"]
    up_count = int(np.count_nonzero(np.isclose(zem, uc.SITE_M + 100.0, atol=1e-3)))
    down_count = int(np.count_nonzero(np.isclose(zem, uc.SITE_M + 800.0, atol=1e-3)))
    metrics = {
        "sample_count": len(names),
        "reason_count_error": int(reason_errors),
        "photons_out_error": int(summary["photons_out"] - expected_out),
        "balance_error": int(
            summary["photons_in"] - summary["photons_out"] - summary["dropped"]
        ),
        "dropped_sum_error": int(
            summary["dropped"] - sum(summary["dropped_by_reason"].values())
        ),
        "reaching_upward_count_error": int(up_count - N),
        "reaching_downward_count_error": int(down_count - N),
        "not_bendable_or_misses_level": int(
            summary["dropped_by_reason"]["not_bendable"]
            + summary["dropped_by_reason"]["misses_level"]
        ),
    }
    details = {
        "seed": 20261126,
        "obdeect_propagation_version": version(),
        "categories": {
            k: {"z_m": v[0], "expected": v[2]} for k, v in CATEGORIES.items()
        },
        "per_category_count": N,
        "turnover_sin": turn_sin,
        "summary": summary,
        "expected_reason": expected_reason,
        "expected_photons_out": expected_out,
        "reference": "category classification written in the adapter (see its docstring)",
    }
    return metrics, details


def main() -> None:
    common.results_main(analyse, __doc__.splitlines()[0])


if __name__ == "__main__":
    main()
