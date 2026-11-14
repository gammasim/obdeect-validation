"""Compare dlufd propagation with CORSIKA's compiled IACT/atmo implementation."""

from __future__ import annotations

import json
import math
import os
import subprocess
from pathlib import Path

import numpy as np


def main() -> None:
    from obdeect.propagation import AtmosphereProfile, RefractionTables, refract

    profile_path = Path(
        os.environ.get("ATMOSPHERE_PROFILE", "/validation-config/atmprof1.dat")
    )
    corsika_dir = Path(os.environ.get("CORSIKA_ATMO_DIR", "/reference/corsika"))
    # fileopen.c and straux.c live in the hessioxxx source tree, not in CORSIKA's atmo directory.
    support_dir = Path(os.environ.get("CORSIKA_SUPPORT_DIR", str(corsika_dir)))
    include_dirs = [corsika_dir] + [
        Path(p)
        for p in os.environ.get("CORSIKA_INCLUDE_DIRS", "").split(os.pathsep)
        if p
    ]
    wrapper = Path(__file__).with_name("corsika_atmosphere_reference.c")
    executable = Path("corsika-atmosphere-reference")
    include_args = [arg for directory in include_dirs for arg in ("-I", str(directory))]
    subprocess.run(
        [
            "cc",
            "-std=c11",
            "-O2",
            *include_args,
            str(wrapper),
            str(corsika_dir / "atmo.c"),
            str(support_dir / "fileopen.c"),
            str(support_dir / "straux.c"),
            "-lm",
            "-o",
            str(executable),
        ],
        check=True,
    )
    reference = subprocess.run(
        [str(executable.resolve()), str(profile_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    corsika_rows = []
    for line in reference.stdout.splitlines():
        fields = line.split()
        if len(fields) == 6:
            corsika_rows.append(list(map(float, fields)))
    if len(corsika_rows) != 24:
        raise RuntimeError(
            f"CORSIKA reference returned {len(corsika_rows)} samples, expected 24"
        )

    profile = AtmosphereProfile(profile_path)
    tables = RefractionTables(profile, 2158.0)
    heights = np.array([row[0] for row in corsika_rows], dtype=np.float64)
    zeniths = np.array([row[1] for row in corsika_rows], dtype=np.float64)
    theta = np.deg2rad(zeniths)
    cx = np.sin(theta)
    cy = np.zeros_like(cx)
    zeros = np.zeros_like(cx)
    bent_cx, _, x_out, _, times, bent = refract(
        tables, heights, cx, cy, zeros, zeros, zeros
    )

    samples = []
    for i, ref in enumerate(corsika_rows):
        height_m, zenith_deg, ref_cx, ref_cz, ref_dx, ref_dt = ref
        samples.append(
            {
                "height_m": height_m,
                "zenith_deg": zenith_deg,
                "corsika_cx": ref_cx,
                "candidate_cx": float(bent_cx[i]),
                "corsika_cz": ref_cz,
                "candidate_cz": math.sqrt(1.0 - float(bent_cx[i]) ** 2),
                "corsika_dx_m": ref_dx / 100.0,
                "candidate_dx_m": float(x_out[i]),
                "corsika_dt_ns": ref_dt,
                "candidate_dt_ns": float(times[i]),
                "bent": bool(bent[i]),
            }
        )
    Path("samples.json").write_text(
        json.dumps(samples, indent=2) + "\n", encoding="utf-8"
    )
    dx_error = max(abs(s["candidate_dx_m"] - s["corsika_dx_m"]) for s in samples)
    cx_error = max(abs(s["candidate_cx"] - s["corsika_cx"]) for s in samples)
    dt_error = max(abs(s["candidate_dt_ns"] - s["corsika_dt_ns"]) for s in samples)
    metrics = {
        "sample_count": len(samples),
        "bent_count": sum(s["bent"] for s in samples),
        "max_direction_cosine_absolute_error": cx_error,
        "max_displacement_absolute_error_m": dx_error,
        "max_travel_time_absolute_error_ns": dt_error,
    }
    Path("metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
