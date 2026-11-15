"""Curved-atmosphere helpers: log pattern, cards, executable selection, flags, run enforcement.

Real IACT files are only read when OBDEECT_TEST_IACT points to a flat shower file (the header
flags and the bunch reader are then checked against it). Everything else uses synthetic input.
"""

import json
import math
import os
import struct
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

curved_showers = pytest.importorskip("obdeect_validation.adapters.curved_showers")
chain_showers = pytest.importorskip("obdeect_validation.adapters.chain_showers")
stream_curved_flag = pytest.importorskip(
    "obdeect_validation.adapters.stream_curved_flag"
)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
CURVED_LOG = "CORSIKA 7.8010\n" + curved_showers.CURVED_LOG_PATTERN + "\nRUN ENDED\n"
FLAT_LOG = "CORSIKA 7.8010\nstandard atmosphere\nRUN ENDED\n"
FLAT_WORD77 = 1211.0  # word 77 of a real flat run: refraction bit set, CURVED bit clear
CURVED_WORD77 = 1104.0  # 0x450: CURVED bit 0x40, refraction 0x10, atmosphere 1


# ---- log pattern and names -------------------------------------------------------------


def test_log_is_curved_matches_only_the_curved_banner():
    assert curved_showers.log_is_curved(CURVED_LOG)
    assert not curved_showers.log_is_curved(FLAT_LOG)
    assert not curved_showers.log_is_curved(curved_showers.CURVED_LOG_PATTERN.lower())


def test_shower_name_format():
    assert (
        curved_showers.shower_name(65, "curved", energy_min=300, energy_max=5000)
        == "gamma-65deg-300-5000GeV-curved.iact"
    )
    assert (
        curved_showers.shower_name(85, "flat", energy_min=1000, energy_max=10000)
        == "gamma-85deg-1000-10000GeV-flat.iact"
    )


def test_unknown_geometry_is_rejected():
    with pytest.raises(ValueError, match="geometry"):
        curved_showers.shower_name(65, "spherical")


# ---- compact-format layout --------------------------------------------------------------


@pytest.mark.parametrize("zenith", [65, 75, 80, 85])
def test_layout_respects_compact_format_limits(zenith):
    layout = curved_showers.layout(zenith)
    for x, y, z, radius in layout["telescopes"]:
        # limit (2): sec(z) * R < 32.7 m for the sphere radius R (cm -> m)
        assert radius / 100.0 / math.cos(math.radians(zenith)) < 32.7
        assert z == radius  # CORSIKA recommendation Z >= R
    # limit (3): cores lie within the scatter radius of the origin and telescopes within the
    # largest distance from the origin, so core-telescope distances are bounded by their sum
    farthest = max(math.hypot(x, y) for x, y, *_ in layout["telescopes"])
    assert (layout["scatter_cm"] + farthest) / 100.0 < 1000.0


def test_telescope_radius_shrinks_with_zenith():
    radii = [curved_showers.telescope_radius_cm(z) for z in (65, 75, 80, 85)]
    assert radii == sorted(radii, reverse=True)
    assert curved_showers.telescope_radius_cm(65) == pytest.approx(1105.0)
    assert curved_showers.telescope_radius_cm(85) == pytest.approx(227.0)


# ---- cards ------------------------------------------------------------------------------


def test_curved_card_has_no_tstart_and_flat_card_has(tmp_path):
    curved = curved_showers.write_card(
        tmp_path / "c", "curved", 65, energy_min=300, energy_max=5000
    )
    flat = curved_showers.write_card(
        tmp_path / "f", "flat", 65, energy_min=300, energy_max=5000
    )
    assert "TSTART" not in curved.read_text()
    assert "TSTART T" in flat.read_text()
    assert curved.name == "gamma-65deg-300-5000GeV-curved.card"


def test_flat_and_curved_cards_differ_only_in_tstart_and_telfil(tmp_path):
    kwargs = {"energy_min": 300, "energy_max": 5000}
    curved = curved_showers.write_card(
        tmp_path / "c", "curved", 65, **kwargs
    ).read_text()
    flat = curved_showers.write_card(tmp_path / "f", "flat", 65, **kwargs).read_text()

    def body(text):
        return [
            line
            for line in text.splitlines()
            if not line.startswith(("TSTART", "TELFIL"))
        ]

    assert body(flat) == body(curved)


def test_card_records_seeds_and_placeholder(tmp_path):
    text = curved_showers.write_card(
        tmp_path, "curved", 80, energy_min=1000, energy_max=10000
    ).read_text()
    for seed in (17, 29, 43, 61):
        assert f"SEED {seed} 0 0\n" in text
    assert "NSHOW 20\n" in text
    assert "THETAP 80 80\n" in text
    assert "DATDIR <corsika data>/\n" in text
    assert "ATMOSPHERE 1 Y\n" in text
    assert "TELFIL ./gamma-80deg-1000-10000GeV-curved.iact\n" in text


def test_legacy_flat_card_is_unchanged():
    """The card generate() has always written, byte for byte (default arguments)."""
    root = Path("/x/root")
    text = chain_showers.card_text(
        showers=20,
        zenith=20,
        energy_min=100,
        energy_max=3000,
        data_dir=root / "data",
        direct=f"{root}/",
        telfil="showers.iact",
    )
    assert text == (
        "RUNNR 1\nEVTNR 1\nNSHOW 20\nPRMPAR 1\nERANGE 100 3000\nESLOPE -2.0\n"
        "THETAP 20 20\nPHIP 0.0 0.0\nSEED 17 0 0\nSEED 29 0 0\nSEED 43 0 0\nSEED 61 0 0\n"
        "OBSLEV 215600.0\nATMOSPHERE 1 Y\nMAGNET 20.0 40.0\nARRANG 0.0\n"
        "HADFLG 0 0 0 0 0 2\nECUTS 0.3 0.3 0.02 0.02\nELMFLG F F\nCSCAT 5 10000.0 0.0\n"
        "DATDIR /x/root/data/\nDIRECT /x/root/\n"
        "TELESCOPE 0.0 0.0 0.0 1200.0\nTELESCOPE 5000.0 0.0 0.0 1200.0\n"
        "TELESCOPE 0.0 5000.0 0.0 1200.0\nTELFIL /x/root/showers.iact\n"
        "IACT SPLIT_AUTO 0\nIACT IO_BUFFER 800MB\nIACT MAX_BUNCHES 1000000\nEXIT\n"
    )


@pytest.mark.parametrize(
    ("geometry", "zenith", "energy_min", "energy_max"),
    curved_showers.FIXTURE_CASES,
)
def test_checked_in_cards_equal_write_card(
    tmp_path, geometry, zenith, energy_min, energy_max
):
    generated = curved_showers.write_card(
        tmp_path, geometry, zenith, energy_min=energy_min, energy_max=energy_max
    )
    checked_in = FIXTURES / "corsika" / generated.name
    assert checked_in.read_text() == generated.read_text()


def test_inputs_json_lists_the_new_showers_with_null_hashes():
    inputs = json.loads((FIXTURES / "inputs.json").read_text())["files"]
    for geometry, zenith, energy_min, energy_max in curved_showers.FIXTURE_CASES:
        name = curved_showers.shower_name(zenith, geometry, energy_min, energy_max)
        entry = inputs[name]
        assert (FIXTURES / entry["card"]).is_file()
        assert entry["sha256"] is None or len(entry["sha256"]) == 64


# ---- CORSIKA executable selection ------------------------------------------------------

BUILD_OPTS = """corsika_version: 78010
variant:
  - executable: corsika_qgs2_urqmd_flat
    config: config_qgs2_urqmd_flat
    atmosphere_geometry: flat
    he_hadronic_model: qgs2
    le_hadronic_model: urqmd
  - executable: corsika_qgs2_urqmd_curved
    config: config_qgs2_urqmd_curved
    atmosphere_geometry: curved
    he_hadronic_model: qgs2
    le_hadronic_model: urqmd
  - executable: corsika_epos_urqmd_flat
    config: config_epos_urqmd_flat
    atmosphere_geometry: flat
    he_hadronic_model: epos  # comment
    le_hadronic_model: urqmd
  - executable: corsika_qgs3_urqmd_curved
    config: config_qgs3_urqmd_curved
    atmosphere_geometry: curved
    he_hadronic_model: qgs3
    le_hadronic_model: urqmd
simtel_tag: "x"
"""


@pytest.fixture
def corsika_root(tmp_path, monkeypatch):
    """Fake CORSIKA installation with an executable for every variant of BUILD_OPTS."""
    root = tmp_path / "corsika7"
    root.mkdir()
    (root / "build_opts.yml").write_text(BUILD_OPTS)
    for name in (
        "corsika_qgs2_urqmd_flat",
        "corsika_qgs2_urqmd_curved",
        "corsika_epos_urqmd_flat",
        "corsika_qgs3_urqmd_curved",
    ):
        executable = root / name
        executable.write_text("#!/bin/sh\n")
        executable.chmod(0o755)
    monkeypatch.setenv("SIMTOOLS_CORSIKA_PATH", str(root))
    monkeypatch.delenv("CORSIKA_HE_MODEL", raising=False)
    return root


def test_curved_executable_prefers_qgs2(corsika_root):
    assert (
        curved_showers.corsika_executable("curved")
        == corsika_root / "corsika_qgs2_urqmd_curved"
    )


def test_flat_executable_prefers_qgs2_then_falls_back_to_epos(
    corsika_root, monkeypatch
):
    assert (
        curved_showers.corsika_executable("flat")
        == corsika_root / "corsika_qgs2_urqmd_flat"
    )
    monkeypatch.setenv("CORSIKA_HE_MODEL", "epos")
    assert (
        curved_showers.corsika_executable("flat")
        == corsika_root / "corsika_epos_urqmd_flat"
    )


def test_environment_model_is_preferred(corsika_root, monkeypatch):
    monkeypatch.setenv("CORSIKA_HE_MODEL", "qgs3")
    assert (
        curved_showers.corsika_executable("curved")
        == corsika_root / "corsika_qgs3_urqmd_curved"
    )


def test_explicit_model_must_exist(corsika_root):
    assert curved_showers.corsika_executable("curved", "qgs3") == corsika_root / (
        "corsika_qgs3_urqmd_curved"
    )
    with pytest.raises(curved_showers.MissingPrerequisite, match="epos"):
        curved_showers.corsika_executable("curved", "epos")


def test_missing_manifest_is_a_missing_prerequisite(tmp_path, monkeypatch):
    monkeypatch.setenv("SIMTOOLS_CORSIKA_PATH", str(tmp_path))
    with pytest.raises(curved_showers.MissingPrerequisite, match="build_opts"):
        curved_showers.corsika_executable("curved")


def test_missing_executable_is_a_missing_prerequisite(corsika_root):
    (corsika_root / "corsika_qgs2_urqmd_curved").unlink()
    with pytest.raises(curved_showers.MissingPrerequisite, match="not an executable"):
        curved_showers.corsika_executable("curved")


def test_no_matching_geometry_lists_available_variants(corsika_root, monkeypatch):
    (corsika_root / "build_opts.yml").write_text(
        "variant:\n  - executable: corsika_qgs2_urqmd_flat\n    config: c\n"
        "    atmosphere_geometry: flat\n    he_hadronic_model: qgs2\n    le_hadronic_model: urqmd\n"
    )
    with pytest.raises(
        curved_showers.MissingPrerequisite, match="available: qgs2/flat"
    ):
        curved_showers.corsika_executable("curved")


def test_standard_library_parser_gives_the_same_variants(corsika_root, monkeypatch):
    monkeypatch.setitem(sys.modules, "yaml", None)  # import yaml fails
    assert (
        curved_showers.corsika_executable("curved")
        == corsika_root / "corsika_qgs2_urqmd_curved"
    )
    assert (
        curved_showers.corsika_executable("flat")
        == corsika_root / "corsika_qgs2_urqmd_flat"
    )
    variants = curved_showers._variants_without_yaml(BUILD_OPTS)
    assert [v["executable"] for v in variants] == [
        "corsika_qgs2_urqmd_flat",
        "corsika_qgs2_urqmd_curved",
        "corsika_epos_urqmd_flat",
        "corsika_qgs3_urqmd_curved",
    ]
    assert variants[2]["he_hadronic_model"] == "epos"


def test_main_maps_missing_prerequisite_to_exit_code_3(tmp_path, monkeypatch):
    monkeypatch.setenv("SIMTOOLS_CORSIKA_PATH", str(tmp_path / "absent"))
    code = curved_showers.main(
        [
            "run",
            "--geometry",
            "curved",
            "--zenith",
            "65",
            "--energy-min",
            "300",
            "--energy-max",
            "5000",
            "--root",
            str(tmp_path / "work"),
        ]
    )
    assert code == 3


# ---- header flags -----------------------------------------------------------------------


def test_decode_flags_flat_run_word_values():
    flags = curved_showers.decode_flags(FLAT_WORD77, 0.0)
    assert flags == {
        "evth_word77_int": 1211,
        "evth_word79": 0,
        "iact_curved_bit": False,
        "corsika_curved_flag": False,
        "refraction_bit": True,
    }


def test_decode_flags_curved_run_word_values():
    flags = curved_showers.decode_flags(CURVED_WORD77, 2.0)
    assert flags["iact_curved_bit"] is True
    assert flags["corsika_curved_flag"] is True
    assert flags["refraction_bit"] is True


def test_decode_flags_reports_disagreement_between_indicators():
    flags = curved_showers.decode_flags(
        16.0, 2.0
    )  # refraction only, but word 79 says curved
    assert flags["iact_curved_bit"] is False
    assert flags["corsika_curved_flag"] is True


class _FakeObject:
    def __init__(self, type_, data):
        self.header = SimpleNamespace(type=type_)
        self._data = data
        self._pos = 0

    def seek(self, position):
        self._pos = position

    def read(self):
        return self._data


class _FakeFile:
    def __init__(self, objects):
        self._objects = objects

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def __iter__(self):
        return iter(self._objects)


def _evth_block(words):
    return struct.pack("<i", 273) + np.asarray(words, dtype="<f4").tobytes()


def test_read_header_flags_uses_the_first_evth(monkeypatch):
    words = np.zeros(273, dtype=np.float32)
    words[76], words[78] = CURVED_WORD77, 2.0
    later = np.zeros(273, dtype=np.float32)
    later[76] = 0.0
    objects = [
        _FakeObject(1200, b""),
        _FakeObject(curved_showers.EVTH_TYPE, _evth_block(words)),
        _FakeObject(curved_showers.EVTH_TYPE, _evth_block(later)),
    ]
    monkeypatch.setattr(curved_showers, "EventIOFile", lambda path: _FakeFile(objects))
    flags = curved_showers.read_header_flags("unused.iact")
    assert flags["corsika_curved_flag"] is True
    assert flags["iact_curved_bit"] is True


def test_read_header_flags_rejects_missing_or_short_evth(monkeypatch):
    monkeypatch.setattr(curved_showers, "EventIOFile", lambda path: _FakeFile([]))
    with pytest.raises(ValueError, match="no EVTH"):
        curved_showers.read_header_flags("unused.iact")
    short = struct.pack("<i", 10) + b"\0" * 40
    monkeypatch.setattr(
        curved_showers,
        "EventIOFile",
        lambda path: _FakeFile([_FakeObject(curved_showers.EVTH_TYPE, short)]),
    )
    with pytest.raises(ValueError, match="10 words"):
        curved_showers.read_header_flags("unused.iact")


SAMPLE = os.environ.get("OBDEECT_TEST_IACT")


@pytest.mark.skipif(not SAMPLE, reason="set OBDEECT_TEST_IACT to a flat IACT file")
def test_real_flat_file_header_and_bunches():
    flags = curved_showers.read_header_flags(SAMPLE)
    assert flags["corsika_curved_flag"] is False
    assert flags["iact_curved_bit"] is False
    events = list(curved_showers.read_bunches(SAMPLE))
    assert events
    _, _, columns = events[0]
    assert {"x", "y", "zem", "cx", "cy", "time", "photons", "wavelength"} <= set(
        columns
    )
    assert np.all(np.abs(columns["x"]) < 1.0e4)  # metres, not cm


# ---- CORSIKA runs and the curved log check ----------------------------------------------


@pytest.fixture
def run_setup(corsika_root, monkeypatch, tmp_path):
    monkeypatch.setattr(chain_showers, "prepare_data", lambda root: {})
    return tmp_path / "run"


def _fake_runner(log_text):
    calls = []

    def runner(executable, card, log, cwd):
        calls.append((Path(executable).name, Path(card).name))
        Path(log).write_text(log_text)
        shower = Path(cwd) / Path(card).name.replace(".card", ".iact")
        shower.write_bytes(b"iact")

    return runner, calls


def test_run_corsika_curved_checks_pattern_and_copies_log(corsika_root, run_setup):
    runner, calls = _fake_runner(CURVED_LOG)
    iact = curved_showers.run_corsika(
        run_setup, "curved", 65, energy_min=300, energy_max=5000, runner=runner
    )
    assert iact == run_setup.resolve() / "gamma-65deg-300-5000GeV-curved.iact"
    assert calls == [
        ("corsika_qgs2_urqmd_curved", "gamma-65deg-300-5000GeV-curved.card")
    ]
    assert iact.with_suffix(".corsika.log").read_text() == CURVED_LOG
    card = (run_setup / "gamma-65deg-300-5000GeV-curved.card").read_text()
    assert f"DATDIR {run_setup.resolve() / 'data'}/" in card
    assert "TSTART" not in card


def test_run_corsika_curved_rejects_flat_log(corsika_root, run_setup):
    runner, _ = _fake_runner(FLAT_LOG)
    with pytest.raises(RuntimeError, match="lacks the CURVED log pattern"):
        curved_showers.run_corsika(
            run_setup, "curved", 65, energy_min=300, energy_max=5000, runner=runner
        )


def test_run_corsika_flat_rejects_curved_log(corsika_root, run_setup):
    runner, _ = _fake_runner(CURVED_LOG)
    with pytest.raises(RuntimeError, match="has the CURVED log pattern"):
        curved_showers.run_corsika(
            run_setup, "flat", 65, energy_min=300, energy_max=5000, runner=runner
        )


def test_run_corsika_without_output_fails(corsika_root, run_setup):
    def runner(executable, card, log, cwd):
        Path(log).write_text(CURVED_LOG)

    with pytest.raises(RuntimeError, match="wrote no"):
        curved_showers.run_corsika(
            run_setup, "curved", 65, energy_min=300, energy_max=5000, runner=runner
        )


# ---- stream case helpers ----------------------------------------------------------------


def _flags(kind):
    if kind == "curved":
        return curved_showers.decode_flags(CURVED_WORD77, 2.0)
    return curved_showers.decode_flags(FLAT_WORD77, 0.0)


def _metadata(kind, **overrides):
    value = 1 if kind == "curved" else 0
    metadata = {"curvature_applied": value, "refraction_applied": 1}
    metadata.update(overrides)
    return metadata


def _items(**overrides):
    curved = {
        "kind": "curved",
        "flags": _flags("curved"),
        "routes": {"cli": _metadata("curved"), "api": _metadata("curved")},
        "log_text": CURVED_LOG,
        "log": "curved.corsika.log",
    }
    flat = {
        "kind": "flat",
        "flags": _flags("flat"),
        "routes": {"cli": _metadata("flat"), "api": _metadata("flat")},
        "log_text": FLAT_LOG,
        "log": "flat.corsika.log",
    }
    curved.update(overrides.get("curved", {}))
    flat.update(overrides.get("flat", {}))
    return [curved, flat]


def test_evaluate_all_consistent_gives_expected_values():
    metrics, details = stream_curved_flag.evaluate(_items())
    assert metrics == {
        "curved_file_curvature_applied": 1,
        "flat_file_curvature_applied": 0,
        "word79_flag_mismatch": 0,
        "iact_curved_bit_mismatch": 0,
        "refraction_bit_mismatch": 0,
        "run_key_missing_count": 0,
        "curved_log_pattern_mismatch": 0,
    }
    assert details["log_missing"] == []
    assert details["files"]["curved"]["run_keys"]["api"]["curvature_applied"] == 1


def test_evaluate_counts_wrong_run_key_against_word79():
    items = _items(curved={"routes": {"cli": _metadata("curved", curvature_applied=0)}})
    metrics, _ = stream_curved_flag.evaluate(items)
    assert metrics["curved_file_curvature_applied"] == 0
    assert metrics["word79_flag_mismatch"] == 1


def test_evaluate_counts_missing_run_key_once_per_route():
    items = _items(flat={"routes": {"cli": {"refraction_applied": 1}, "api": None}})
    metrics, _ = stream_curved_flag.evaluate(items)
    assert (
        metrics["run_key_missing_count"] == 3
    )  # cli curvature; api record (both keys)
    assert metrics["word79_flag_mismatch"] == 0


def test_evaluate_counts_refraction_and_iact_bit_disagreement():
    flags = dict(_flags("curved"), refraction_bit=False)
    items = _items(curved={"flags": flags, "routes": {"cli": _metadata("curved")}})
    metrics, _ = stream_curved_flag.evaluate(items)
    assert metrics["refraction_bit_mismatch"] == 1
    flags = dict(_flags("curved"), iact_curved_bit=False)
    metrics, _ = stream_curved_flag.evaluate(
        _items(curved={"flags": flags, "routes": {"cli": _metadata("curved")}})
    )
    assert metrics["iact_curved_bit_mismatch"] == 1


def test_evaluate_log_pattern_enforced_and_missing_log_recorded():
    metrics, _ = stream_curved_flag.evaluate(_items(flat={"log_text": CURVED_LOG}))
    assert metrics["curved_log_pattern_mismatch"] == 1
    metrics, details = stream_curved_flag.evaluate(
        _items(curved={"log_text": None, "log": "missing"})
    )
    assert metrics["curved_log_pattern_mismatch"] == 0
    assert details["log_missing"] == ["curved"]


def test_routes_for_skips_api_when_not_importable(monkeypatch, tmp_path):
    calls = []

    def fake_command(argv, log, stdout_file=None, **kwargs):
        calls.append(list(argv))
        Path(stdout_file).write_bytes(b"")

    monkeypatch.setattr(stream_curved_flag, "command", fake_command)
    monkeypatch.setattr(
        stream_curved_flag, "run_metadata", lambda path: {"curvature_applied": 1}
    )
    monkeypatch.setitem(sys.modules, "obdeect.stream", None)
    routes = stream_curved_flag.routes_for(Path("a.iact"), tmp_path, "curved")
    assert set(routes) == {"cli"}
    assert calls == [["obdeect", "stream", "from-corsika", Path("a.iact")]]


def test_routes_for_uses_api_when_importable(monkeypatch, tmp_path):
    stream = pytest.importorskip("obdeect.stream")
    monkeypatch.setattr(
        stream_curved_flag, "command", lambda argv, log, stdout_file=None, **kw: None
    )
    monkeypatch.setattr(
        stream_curved_flag, "run_metadata", lambda path: {"path": str(path)}
    )
    converted = []
    monkeypatch.setattr(
        stream,
        "corsika_to_stream",
        lambda src, dst: converted.append((src, dst)),
        raising=False,
    )
    routes = stream_curved_flag.routes_for(Path("a.iact"), tmp_path, "flat")
    assert set(routes) == {"cli", "api"}
    assert converted == [(Path("a.iact"), tmp_path / "flat.api.gst")]
