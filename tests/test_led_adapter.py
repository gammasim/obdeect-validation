import importlib.util
import json
import math
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[1]


class Values(list):
    def astype(self, _dtype):
        return self


def load_adapter(monkeypatch):
    numpy = types.ModuleType("numpy")
    numpy.float64 = float
    numpy.sum = sum
    numpy.dot = lambda left, right: sum(x * y for x, y in zip(left, right, strict=True))
    numpy.asarray = lambda values, dtype=None: Values(values)
    numpy.hypot = math.hypot

    eventio = types.ModuleType("eventio")
    eventio.__path__ = []
    iact = types.ModuleType("eventio.iact")
    iact.IACTFile = lambda _path: None
    obdeect = types.ModuleType("obdeect")
    obdeect.__path__ = []
    stream = types.ModuleType("obdeect.stream")
    stream.Reader = lambda _path: iter(())
    for name, module in (
        ("numpy", numpy),
        ("eventio", eventio),
        ("eventio.iact", iact),
        ("obdeect", obdeect),
        ("obdeect.stream", stream),
    ):
        monkeypatch.setitem(sys.modules, name, module)

    path = ROOT / "src" / "obdeect_validation" / "adapters" / "led_reference.py"
    spec = importlib.util.spec_from_file_location("led_reference_test_module", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module, stream, iact


def test_empty_native_output_is_an_error(monkeypatch, tmp_path):
    adapter, stream, _iact = load_adapter(monkeypatch)
    stream.Reader = lambda _path: iter([{"type": "event"}])
    with pytest.raises(ValueError, match="obdeect produced no photons"):
        adapter.read_stream(tmp_path / "empty.gst")


def test_empty_reference_output_is_an_error(monkeypatch, tmp_path):
    adapter, _stream, iact = load_adapter(monkeypatch)

    class EmptyFile:
        def __enter__(self):
            return iter(())

        def __exit__(self, *_args):
            return False

    adapter.IACTFile = lambda _path: EmptyFile()
    with pytest.raises(ValueError, match="pixled produced no photons"):
        adapter.read_reference(tmp_path / "empty.dat.gz")


def test_mismatched_event_counts_do_not_write_metrics(monkeypatch, tmp_path):
    adapter, _stream, _iact = load_adapter(monkeypatch)
    monkeypatch.chdir(tmp_path)
    config = tmp_path / "simtel-config"
    config.mkdir()
    profile = config / "atm_profile_model_1.dat"
    profile.write_text("fixture\n", encoding="utf-8")
    monkeypatch.setenv("SIMTEL_CONFIG_PATH", str(config))
    pixled = tmp_path / "pixled"
    pixled.write_text("fixture binary placeholder\n", encoding="utf-8")
    monkeypatch.setattr(adapter.shutil, "which", lambda name: str(pixled) if name == "pixled" else name)
    monkeypatch.setattr(
        adapter.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(stdout="pixled log", returncode=0, check_returncode=lambda: None),
    )
    monkeypatch.setattr(adapter, "version", lambda _name: "0.0.test")
    monkeypatch.setattr(adapter, "read_stream", lambda _path: ([10.0], 10.0, 0.0, 0.0))
    monkeypatch.setattr(adapter, "read_reference", lambda _path: ([5.0, 5.0], 10.0, 0.0, 0.0))

    adapter.main()
    metrics = json.loads((tmp_path / "metrics.json").read_text(encoding="utf-8"))
    assert metrics["event_count_delta"] == 1
    assert metrics["max_event_photon_delta"] == 5
