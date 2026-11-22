"""Independent stage checks for the plane horizon and scattered-light timing."""

from pathlib import Path

import numpy as np
import pytest

from obdeect_validation.adapters.refraction_curved_reference import (
    NuProfile,
    profile_path,
)

ROOT = Path(__file__).resolve().parents[1]
SITE = 2000.0
C = 0.299792458


def _stream(path, *, angles, refracted=False, telescopes=1):
    from obdeect.stream import Writer

    with Writer(str(path), "validation") as writer:
        writer.write(
            "run",
            {
                "run_number": 1,
                "observation_level_m": SITE,
                "telescope_x_m": [0.0, 200.0][:telescopes],
                "telescope_y_m": [0.0] * telescopes,
                "telescope_z_m": [0.0] * telescopes,
                "telescope_r_m": [12.0] * telescopes,
                "refraction_applied": int(refracted),
                "wavelength_min_nm": 355.0,
                "wavelength_max_nm": 355.0,
            },
        )
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1, "frame": 0})
        n = len(angles)
        columns = {
            "x": np.zeros(n),
            "y": np.zeros(n),
            "z": np.full(n, 10000.0),
            "cx": np.sin(np.radians(angles)),
            "cy": np.zeros(n),
            "t": np.zeros(n),
            "w": np.ones(n),
            "wavelength": np.full(n, 355.0),
        }
        writer.write_photons(
            1, 0, 0, {k: np.asarray(v, dtype=np.float32) for k, v in columns.items()}
        )
        writer.write("telescope_end", {"telescope_id": 1, "photons": float(n)})
        writer.write("event_end", {})
        writer.write("run_end", {})


def _photons(path):
    from obdeect.stream import Reader

    chunks = [r["columns"] for r in Reader(str(path)) if r["type"] == "photons"]
    assert chunks
    return {key: np.concatenate([c[key] for c in chunks]) for key in chunks[0]}


def test_plane_horizon_against_independent_snell_integral(tmp_path):
    propagate = pytest.importorskip("obdeect.propagation").propagate

    angles = np.array([84.0, 86.0, 87.0, 88.0, 89.0])
    source, output = tmp_path / "sky.pst", tmp_path / "plane.pst"
    _stream(source, angles=angles)
    propagate(
        source,
        output,
        site_altitude_m=SITE,
        seed=1,
        profile=profile_path(),
        extinction=False,
        mode="weighted",
    )
    actual = _photons(output)
    assert len(actual["x"]) == len(angles)
    profile = NuProfile(profile_path())
    heights = np.linspace(SITE, SITE + 10000.0, 200001)
    nu = profile(heights)
    for index, angle in enumerate(angles):
        sin_e = float(np.float32(np.sin(np.radians(angle))))
        sin_h = sin_e * (1.0 + nu[-1]) / (1.0 + nu)
        cos_h = np.sqrt(1.0 - sin_h**2)
        expected_x = np.trapezoid(sin_h / cos_h, heights)
        expected_t = np.trapezoid((1.0 + nu) / cos_h, heights) / C
        assert actual["x"][index] == pytest.approx(expected_x, abs=0.2)
        assert actual["t"][index] == pytest.approx(expected_t, abs=0.8)


def test_scattered_photon_timing_includes_both_air_paths(tmp_path):
    propagate = pytest.importorskip("obdeect.propagation").propagate

    source = tmp_path / "sky.pst"
    _stream(source, angles=[0.0], refracted=True, telescopes=2)
    outputs = []
    for use_profile in (False, True):
        output = tmp_path / f"scatter-{use_profile}.pst"
        propagate(
            source,
            output,
            site_altitude_m=SITE,
            seed=17,
            mode="weighted",
            transmission=ROOT / "fixtures/atm_trans_2000_1_10_0_0_2000.dat",
            scattering=ROOT / "fixtures/scattering_exponential.dat",
            scattering_points=16,
            scattering_primaries=[1],
            telescopes=[2],
            profile=profile_path() if use_profile else None,
        )
        outputs.append(_photons(output))
    vacuum, air = outputs
    assert len(vacuum["t"]) > 0
    for key in ("x", "y", "cx", "cy", "zem", "w"):
        np.testing.assert_array_equal(vacuum[key], air[key])
    profile = NuProfile(profile_path())
    for i, height in enumerate(vacuum["zem"]):
        up = np.linspace(float(height), SITE + 10000.0, 10001)
        down = np.linspace(SITE, float(height), 10001)
        cos_out = np.sqrt(
            1.0 - float(vacuum["cx"][i]) ** 2 - float(vacuum["cy"][i]) ** 2
        )
        expected = (
            np.trapezoid(profile(up), up) + np.trapezoid(profile(down), down) / cos_out
        ) / C
        assert float(air["t"][i] - vacuum["t"][i]) == pytest.approx(expected, abs=0.1)
