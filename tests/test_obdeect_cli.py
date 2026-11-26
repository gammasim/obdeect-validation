from pathlib import Path

import numpy as np
import pytest

from obdeect_validation.adapters import common, streams
from obdeect_validation.adapters import obdeect_cli as cli


def test_obdeect_maps_options_to_flags_in_keyword_order():
    assert cli.obdeect(
        "stream",
        "to-simtel",
        "-",
        o="out.simtel",
        camera_map=Path("map.txt"),
        keep_telpos_height=True,
        atmprof=None,
        scaled_factor=False,
        extra=["--set", "a=1"],
    ) == [
        "obdeect",
        "stream",
        "to-simtel",
        "-",
        "-o",
        "out.simtel",
        "--camera-map",
        Path("map.txt"),
        "--keep-telpos-height",
        "--set",
        "a=1",
    ]


def test_weighted_propagation_is_the_chain_command_line():
    assert cli.weighted_propagation(
        "cut.dat", (300.0, 600.0), "in.gst", "out.gst", seed=123
    ) == [
        "obdeect",
        "propagation",
        "run",
        "--transmission",
        "cut.dat",
        "--pathlen",
        "0",
        "--site-altitude",
        2156.0,
        "--seed",
        123,
        "--mode",
        "weighted",
        "--wavelength-range",
        300.0,
        600.0,
        "-i",
        "in.gst",
        "-o",
        "out.gst",
    ]


def test_compile_model_and_trace_stream():
    assert cli.compile_model("m.json", zenith_deg=20, simulation_models="/models") == [
        "obdeect",
        "ray-tracing",
        "compile",
        "--simulation-models",
        "/models",
        "--model",
        "LSTN-01",
        "--version",
        "7.0.0",
        "--seed",
        1,
        "--zenith-deg",
        20,
        "--output",
        "m.json",
    ]
    assert cli.trace_stream("m.json", "-", "-", seed=5, threads=1) == [
        "obdeect",
        "ray-tracing",
        "trace-stream",
        "--model",
        "m.json",
        "--seed",
        5,
        "--threads",
        1,
        "-i",
        "-",
        "-o",
        "-",
    ]


def test_prepare_chain_writes_band_cut_table_and_camera_map(tmp_path, monkeypatch):
    source = tmp_path / "showers.iact"
    source.with_suffix(".corsika.log").write_text(
        " CHERENKOV WAVELENGTH BAND FROM  240.0 TO  700.0 NANOMETER\n"
    )
    table = tmp_path / "table.dat"
    table.write_text("# H2= 2.156\n 999 0.1 0.2\n 1000 0.1 0.2\n")
    monkeypatch.setattr(common, "TRANSMISSION_TABLE", table)
    chain = cli.prepare_chain(tmp_path, source, model=None)
    assert chain.band == (240.0, 700.0)
    assert chain.model is None
    assert chain.camera_map.read_text() == "1 1 0.0\n"
    assert "1000" not in chain.cut.read_text()


def test_ulp32_floors_at_one():
    np.testing.assert_array_equal(
        streams.ulp32([0.0, 1.0, 2.0]),
        np.spacing(np.float32([1.0, 1.0, 2.0])).astype(np.float64),
    )


@pytest.mark.skipif(
    not common.obdeect_has("obdeect.stream"), reason="obdeect.stream is not installed"
)
def test_write_and_read_stream_round_trip(tmp_path):
    columns = {
        name: np.arange(5, dtype=np.float32)
        for name in ("x", "y", "z", "cx", "cy", "cz", "t", "w")
    }
    columns["cz"] = -np.ones(5, dtype=np.float32)
    columns["cx"] = columns["cy"] = np.zeros(5, dtype=np.float32)
    run = {
        "run_number": 1,
        "observation_level_m": 2000.0,
        "telescope_x_m": [0.0],
        "telescope_y_m": [0.0],
        "telescope_z_m": [0.0],
        "refraction_applied": 0,
        "wavelength_min_nm": 300.0,
        "wavelength_max_nm": 700.0,
    }
    path = streams.write_stream(
        tmp_path / "s.pst", "test", run, [(1, columns)], chunk_rows=2
    )
    content = streams.read_stream(path)
    np.testing.assert_array_equal(content.columns["t"], columns["t"])
    assert content.telescope.tolist() == [1] * 5
    assert content.frames == {0}
