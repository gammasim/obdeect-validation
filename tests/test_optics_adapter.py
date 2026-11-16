import json
import math
from pathlib import Path
from types import SimpleNamespace

import pytest

astropy = pytest.importorskip("astropy")
import astropy.units as u  # noqa: E402
from astropy.table import QTable  # noqa: E402

from obdeect_validation.adapters import simtools_validate_optics as adapter  # noqa: E402

LIMITS = {"psf_cm": 0.02, "eff_area": 0.01, "eff_flen": 0.002}


def table_file(path: Path, *, include_units: bool = True, rows: int = 2) -> Path:
    results = path / "results"
    results.mkdir(parents=True, exist_ok=True)
    table = QTable()
    if include_units:
        table["off_x"] = [0, 60][:rows] * u.arcmin
        table["off_y"] = [0, 0][:rows] * u.deg
        table["psf_cm"] = [10, 20][:rows] * u.mm
        table["eff_area"] = [10000, 20000][:rows] * u.cm**2
        table["eff_flen"] = [float("nan"), 28][:rows] * u.m
    else:
        for name in ("off_x", "off_y", "psf_cm", "eff_area", "eff_flen"):
            table[name] = [1.0] * rows
    target = results / "result.ecsv"
    table.write(target, format="ascii.ecsv")
    return target


def side(psf, area, flen, flen_err=None, off_x=None):
    n = len(psf)
    return {
        "off_x": off_x or [float(i) for i in range(n)],
        "off_y": [0.0] * n,
        "psf_cm": psf,
        "eff_area": area,
        "eff_flen": flen,
        "eff_flen_err": flen_err
        or [None if flen[i] is None else 0.0 for i in range(n)],
    }


def test_table_converts_units_and_preserves_on_axis_undefined_flen(
    tmp_path: Path,
) -> None:
    table_file(tmp_path)
    values = adapter._table(tmp_path)
    assert values["off_x"] == pytest.approx([0, 1])
    assert values["off_y"] == pytest.approx([0, 0])
    assert values["psf_cm"] == pytest.approx([1, 2])
    assert values["eff_area"] == pytest.approx([1, 2])
    assert values["eff_flen"] == [None, 2800]
    assert values["eff_flen_err"] == [None, None]


def test_table_rejects_missing_units(tmp_path: Path) -> None:
    table_file(tmp_path, include_units=False)
    with pytest.raises(adapter.AdapterError, match="has no unit"):
        adapter._table(tmp_path)


def test_table_rejects_nonfinite_offsets_and_infinite_focal_length(
    tmp_path: Path,
) -> None:
    target = table_file(tmp_path)
    table = QTable.read(target, format="ascii.ecsv")
    table["off_x"][0] = float("nan") * u.arcmin
    table.write(target, format="ascii.ecsv", overwrite=True)
    with pytest.raises(adapter.AdapterError, match="undefined off_x"):
        adapter._table(tmp_path)

    table["off_x"][0] = 0 * u.arcmin
    table["eff_flen"][1] = float("inf") * u.m
    table.write(target, format="ascii.ecsv", overwrite=True)
    with pytest.raises(adapter.AdapterError, match="undefined eff_flen"):
        adapter._table(tmp_path)


def test_table_rejects_missing_and_empty_outputs(tmp_path: Path) -> None:
    with pytest.raises(adapter.AdapterError, match="found 0"):
        adapter._table(tmp_path)
    table_file(tmp_path, rows=0)
    with pytest.raises(adapter.AdapterError, match="has no rows"):
        adapter._table(tmp_path)


def test_rows_signed_delta_and_stat_threshold() -> None:
    sim = side([1.0, 2.0], [1.0, 2.0], [None, 28.0], flen_err=[None, 0.0])
    obd = side([1.0, 2.1], [1.0, 2.0], [None, 28.0 * 1.001], flen_err=[None, 0.0])
    rows = adapter._rows(sim, obd, LIMITS)
    assert rows[1]["psf_cm"]["delta"] == pytest.approx(0.05)
    assert rows[1]["psf_cm"]["exceeded"] is True  # 5% > 2%
    assert rows[1]["eff_flen"]["delta"] == pytest.approx(0.001)
    assert rows[1]["eff_flen"]["exceeded"] is False  # 0.1% < 0.2%
    assert rows[0]["eff_flen"] is None


def test_rows_statistical_error_widens_flen_threshold() -> None:
    # 0.5% deviation exceeds the 0.2% limit, but 3 * rel_se (se = 0.3 / 28 per backend) covers it.
    sim = side([1.0, 1.0], [1.0, 1.0], [None, 28.0], flen_err=[None, 0.3])
    obd = side([1.0, 1.0], [1.0, 1.0], [None, 28.14], flen_err=[None, 0.3])
    row = adapter._rows(sim, obd, LIMITS)[1]["eff_flen"]
    assert row["rel_se"] == pytest.approx(math.hypot(0.3, 0.3 * 28.14 / 28.0) / 28.0)
    assert row["threshold"] == pytest.approx(3 * row["rel_se"])
    assert row["exceeded"] is False


def test_rows_reject_one_sided_undefined_and_zero_reference() -> None:
    sim = side([1.0, 1.0], [1.0, 1.0], [None, 28.0])
    obd = side([1.0, 1.0], [1.0, 1.0], [1.0, 28.0])
    with pytest.raises(adapter.AdapterError, match="only one backend"):
        adapter._rows(sim, obd, LIMITS)
    zero = side([0.0, 1.0], [1.0, 1.0], [None, 28.0])
    with pytest.raises(adapter.AdapterError, match="reference psf_cm is zero"):
        adapter._rows(zero, side([1.0, 1.0], [1.0, 1.0], [None, 28.0]), LIMITS)


def test_summary_reports_worst_offset_and_sign() -> None:
    sim = side(
        [1.0, 2.0, 4.0], [1.0, 1.0, 1.0], [None, 28.0, 28.0], flen_err=[None, 0.0, 0.0]
    )
    obd = side(
        [1.0, 2.2, 3.6], [1.0, 1.0, 1.0], [None, 28.0, 28.0], flen_err=[None, 0.0, 0.0]
    )
    summary = adapter._summary(adapter._rows(sim, obd, LIMITS))
    assert summary["psf_cm"]["max_abs_delta"] == pytest.approx(0.1)
    assert summary["psf_cm"]["signed_delta_at_max"] == pytest.approx(
        0.1
    )  # first of the two equal maxima
    assert [entry["angle_deg"] for entry in summary["psf_cm"]["exceeded_offsets"]] == [
        1.0,
        2.0,
    ]
    assert summary["eff_area"]["exceeded_offsets"] == []


def test_execute_writes_application_config_outputs_and_compares_paired_offsets(
    monkeypatch, tmp_path: Path
) -> None:
    pytest.importorskip("matplotlib")
    config = tmp_path / "input.json"
    config.write_text(
        json.dumps(
            {
                "simtools_config": {
                    "application": "simtools-validate-optics",
                    "configuration": {"model_version": "7.0.0", "telescope": "LSTN-01"},
                },
                "limits": LIMITS,
            }
        ),
        encoding="utf-8",
    )
    configs = []

    def fake_run(command, **_kwargs):
        written = json.loads(Path(command[2]).read_text(encoding="utf-8"))
        configs.append(written)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(adapter.subprocess, "run", fake_run)
    monkeypatch.setattr(
        adapter,
        "_table",
        lambda path: {
            "off_x": [0.0, 1.0],
            "off_y": [0.0, 0.0],
            "psf_cm": [1.0, 2.0],
            "eff_area": [3.0, 4.0],
            "eff_flen": [None, 28.0],
            "eff_flen_err": [None, 0.0],
        },
    )
    metrics = adapter.execute(config, tmp_path / "output")
    assert [
        entry["applications"][0]["configuration"]["ray_tracing_backend"]
        for entry in configs
    ] == [
        "sim_telarray",
        "obdeect",
    ]
    assert all(
        entry["applications"][0]["configuration"]["model_version"] == "7.0.0"
        for entry in configs
    )
    assert metrics["n_offsets"] == 2
    assert metrics["max_rel_eff_flen"] == 0
    assert metrics["exceeded_limits"] == 0
    assert (
        json.loads((tmp_path / "output" / "metrics.json").read_text(encoding="utf-8"))
        == metrics
    )
    details = json.loads(
        (tmp_path / "output" / "details.json").read_text(encoding="utf-8")
    )
    assert details["status"].startswith(
        "approved 2026-11-15, from obdeect-validation-org docs/phase1.md"
    )
    assert (tmp_path / "output" / "optics_relative_difference.png").stat().st_size > 0


def test_execute_rejects_mismatched_offset_pairing(monkeypatch, tmp_path: Path) -> None:
    config = tmp_path / "input.json"
    config.write_text(
        json.dumps(
            {
                "simtools_config": {
                    "application": "simtools-validate-optics",
                    "configuration": {"model_version": "7.0.0"},
                },
                "limits": LIMITS,
            }
        ),
        encoding="utf-8",
    )
    result_tables = iter(
        [
            {
                "off_x": [0.0],
                "off_y": [0.0],
                "psf_cm": [1.0],
                "eff_area": [1.0],
                "eff_flen": [2.0],
                "eff_flen_err": [0.0],
            },
            {
                "off_x": [1.0],
                "off_y": [0.0],
                "psf_cm": [1.0],
                "eff_area": [1.0],
                "eff_flen": [2.0],
                "eff_flen_err": [0.0],
            },
        ]
    )
    monkeypatch.setattr(
        adapter.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=0),
    )
    monkeypatch.setattr(adapter, "_table", lambda _path: next(result_tables))
    with pytest.raises(adapter.AdapterError, match="offset rows do not align"):
        adapter.execute(config, tmp_path / "output")
