"""End-to-end chain: light source -> atmosphere -> optics -> sim_telarray input.

    obdeect sources run point-source | obdeect propagation run | obdeect ray-tracing trace-stream
        | obdeect stream to-simtel

Moved from obdeect/validation/test_chain.py. Skipped unless the four obdeect components are
installed (pip install obdeect). Checks that the components work together through the photon
stream, not physics accuracy (each component has its own tests in obdeect).
The transmission table is a copy of obdeect src/propagation/cpp/tests/data/
atm_trans_2000_1_10_0_0_2000.dat (fixtures/).
"""

import subprocess
import sys
from pathlib import Path

import pytest

obdeect_stream = pytest.importorskip("obdeect.stream")
OpticalModel = pytest.importorskip("obdeect.ray_tracing").OpticalModel
for _module in ("obdeect.propagation", "obdeect.sources"):
    pytest.importorskip(_module)

REPOSITORY = Path(__file__).resolve().parents[1]
TRANSMISSION = REPOSITORY / "fixtures" / "atm_trans_2000_1_10_0_0_2000.dat"
OBDEECT = [sys.executable, "-m", "obdeect.stream.cli"]
CAMERA_FRAME = 3

SOURCE = [
    "sources",
    "run",
    "point-source",
    "--position",
    "0,0,1000",
    "--photons",
    "5e10",
    "--bunch-size",
    "10",
    "--telescope-radius",
    "3",
    "--spectrum",
    "fixed:wavelength_nm=400",
    "--observation-level",
    "2000",
    "--atmosphere-id",
    "1",
    "--seed",
    "3",
]
ATMOSPHERE = [
    "propagation",
    "run",
    "--transmission",
    str(TRANSMISSION),
    "--site-altitude",
]
ATMOSPHERE += ["2000", "--wavelength-range", "300", "600", "--seed", "5"]
OPTICS = ["ray-tracing", "trace-stream", "--seed", "7", "--pointing-zenith-deg", "0"]
OPTICS += ["--pointing-azimuth-deg", "0"]

# Synthetic telescope (not a CTAO model): seven circular 1 m mirrors on a Davies-Cotton dish,
# f = 8 m. Same inputs as the obdeect ray_tracing test telescope.
MIRROR_LIST = """# x y diameter focal shape
   0.0     0.0  100.0  0.0  0
 125.0     0.0  100.0  0.0  0
-125.0     0.0  100.0  0.0  0
  62.5   108.3  100.0  0.0  0
 -62.5   108.3  100.0  0.0  0
  62.5  -108.3  100.0  0.0  0
 -62.5  -108.3  100.0  0.0  0
"""


def _record(name, value, unit=None, *, file=False):
    return {
        "parameter": name,
        "parameter_version": "1.0.0",
        "value": value,
        "unit": unit,
        "file": file,
    }


def single_mirror_model_inputs():
    """Records and file contents of the synthetic telescope."""
    records = {
        "mirror_class": _record("mirror_class", 0),
        "focal_length": _record("focal_length", 800.0, "cm"),
        "mirror_offset": _record("mirror_offset", -50.0, "cm"),
        "dish_shape_length": _record("dish_shape_length", 800.0, "cm"),
        "parabolic_dish": _record("parabolic_dish", False),
        "mirror_focal_length": _record("mirror_focal_length", 0.0, "cm"),
        "mirror_list": _record("mirror_list", "mirrors.dat", file=True),
        "focus_offset": _record(
            "focus_offset", [0.0, 20.0, 0.0, 0.0], ["cm", "deg", "null", "null"]
        ),
        "camera_body_diameter": _record("camera_body_diameter", 60.0, "cm"),
        "camera_body_shape": _record("camera_body_shape", 2),
        "camera_depth": _record("camera_depth", 0.0, "cm"),
        "telescope_transmission": _record(
            "telescope_transmission", [0.9, 1, 0.1, 5.0, 2.0, 0.0]
        ),
    }
    return records, {"mirror_list": MIRROR_LIST}


def obdeect(*args, stdin=None):
    return subprocess.run(
        [*OBDEECT, *map(str, args)], input=stdin, capture_output=True, check=True
    ).stdout


@pytest.fixture(scope="module")
def model(tmp_path_factory):
    records, files = single_mirror_model_inputs()
    path = tmp_path_factory.mktemp("model") / "model.json"
    OpticalModel.from_records(records, files.__getitem__, configuration_seed=1).save(
        path
    )
    return path


@pytest.fixture(scope="module")
def chain(model):
    """The stream after each step."""
    source = obdeect(*SOURCE)
    atmosphere = obdeect(*ATMOSPHERE, "-i", "-", "-o", "-", stdin=source)
    camera = obdeect(*OPTICS, "--model", model, stdin=atmosphere)
    return {"source": source, "atmosphere": atmosphere, "camera": camera}


def records(data, tmp_path, name):
    path = tmp_path / f"{name}.pst"
    path.write_bytes(data)
    return list(obdeect_stream.Reader(str(path)))


def photon_weight(recs):
    return sum(float(r["columns"]["w"].sum()) for r in recs if r["type"] == "photons")


def test_each_step_writes_a_complete_stream(chain, tmp_path):
    for name, data in chain.items():
        types = [r["type"] for r in records(data, tmp_path, name)]
        assert types[0] == "stream" and types[-1] == "run_end", name
        assert "photons" in types, name


def test_photons_reach_the_camera_frame(chain, tmp_path):
    camera = records(chain["camera"], tmp_path, "camera")
    telescopes = [r for r in camera if r["type"] == "telescope"]
    assert telescopes and all(
        t["metadata"]["frame"] == CAMERA_FRAME for t in telescopes
    )
    source = photon_weight(records(chain["source"], tmp_path, "source"))
    atmosphere = photon_weight(records(chain["atmosphere"], tmp_path, "atmosphere"))
    detected = photon_weight(camera)
    # Extinction and the optics each remove light, but never all of it.
    assert 0 < detected < atmosphere < source


def test_efficiency_status_records_the_atmosphere(chain, tmp_path):
    status = [
        r["metadata"]
        for r in records(chain["camera"], tmp_path, "c")
        if r["type"] == "efficiency_status"
    ]
    assert any(s.get("atmosphere.state") == "applied" for s in status)


def test_chain_is_deterministic(chain, model):
    source = obdeect(*SOURCE)
    assert source == chain["source"]
    atmosphere = obdeect(*ATMOSPHERE, "-i", "-", "-o", "-", stdin=source)
    assert atmosphere == chain["atmosphere"]
    assert obdeect(*OPTICS, "--model", model, stdin=atmosphere) == chain["camera"]


def test_camera_photons_convert_to_sim_telarray_input(chain, tmp_path):
    camera_map = tmp_path / "cam.txt"
    camera_map.write_text("1 1 0.0\n")
    simtel = tmp_path / "chain.simtel"
    obdeect(
        "stream",
        "to-simtel",
        "--camera-map",
        camera_map,
        "-o",
        simtel,
        "-",
        stdin=chain["camera"],
    )
    assert simtel.stat().st_size > len(chain["camera"]) // 2
