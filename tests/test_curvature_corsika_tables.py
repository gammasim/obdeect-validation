"""CORSIKA table trace: parser, DISTIP/TOFIP arithmetic, residual binning and the output patch.

Synthetic traces use the Fortran field widths of the patch formats (BUNCH 9611, CONST 9602/9603/9604,
ROW 9605). The patch is dry-run against a copy of corsika.F when the source is available (CORSIKA_SOURCE,
or the container path); otherwise that test is skipped.
"""

import json
import math
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from obdeect_validation import runner
from obdeect_validation.adapters import curvature_corsika_tables as tables

REPO = Path(__file__).resolve().parents[1]
PATCH = REPO / "fixtures" / "corsika-patches" / "curvature-tables-debug.patch"
CASE = REPO / "cases" / "propagation" / "curvature-corsika-tables.json"
OBS_CM = 215600.0
EARTH_CM = 637131500.0
C_CM_PER_S = 2.99792458e10


# ---- Fortran-format writers for synthetic traces -------------------------------------------


def fortran_e(value: float) -> str:
    """An E23.15 field: 0.ddd...(15 digits)E+ee, right-justified in 23 columns."""
    if value == 0.0:
        mantissa, exponent = 0.0, 1
    else:
        exponent = math.floor(math.log10(abs(value))) + 1
        mantissa = value / 10.0**exponent
        if abs(mantissa) >= 1.0:  # rounding at a decade boundary
            mantissa /= 10.0
            exponent += 1
    return f"{mantissa:.15f}E{exponent:+03d}".rjust(23)


def bunch_line(event, index, n2, m2, reals) -> str:
    """BUNCH line: 'BUNCH', I10, I7, I6, I6 and 22 fields of 1X, E23.15."""
    assert len(reals) == tables.BUNCH_FIELDS
    head = f"BUNCH{event:10d}{index:7d}{n2:6d}{m2:6d}"
    return head + "".join(" " + fortran_e(value) for value in reals)


def const_line(name, value) -> str:
    return f"CONST {name:>12s} {fortran_e(float(value))}"


def row_line(n, m, distef, tof) -> str:
    return f"ROW{n:6d}{m:6d} {fortran_e(distef)} {fortran_e(tof)}"


def table_value(n, m):
    """Synthetic table at grid node (N, M): DISTEF in cm and TOF in ns."""
    return 1000.0 + 3.0 * n + 0.5 * m * m, 1.0e4 + 7.0 * n + 11.0 * m


def bilinear_reference(theta4, height4, n2, m2):
    """Brute-force bilinear interpolation of table_value (explicit weights, above OBS)."""
    u, v = theta4 - n2, height4 - m2
    corners = [
        (n2, m2, (1 - u) * (1 - v)),
        (n2 + 1, m2, u * (1 - v)),
        (n2, m2 + 1, (1 - u) * v),
        (n2 + 1, m2 + 1, u * v),
    ]
    distance = sum(w * table_value(n, m)[0] for n, m, w in corners)
    time = sum(w * table_value(n, m)[1] for n, m, w in corners)
    return distance, time


def synthetic_trace(path, bunches, rdist_scale=1.0, tofip_offset=0.0, nodes=True):
    """Write a trace; bunches = [(wemis, zapp_cm, zem_cm, temis_s)]; values from table_value."""
    lines = ["CURVTAB 1 OBDEECT VALIDATION TRACE (NOT CORSIKA OUTPUT)"]
    for name, value in (
        ("C1_CM", EARTH_CM),
        ("C25_CMPS", C_CM_PER_S),
        ("OBSLEV_CM", OBS_CM),
        ("ETADSN", 2.8e-4),
        ("MODATM", 1),
        ("IATMOX", 1),
    ):
        lines.append(const_line(name, value))
    lines.append(f"CONST {'FREFRX':>12s} {'T':>2s}")
    if nodes:
        for theta in tables.ROW_THETA_DEG:
            for height in tables.ROW_HEIGHT_KM:
                n, m = 4 * theta, 4 * height
                lines.append(row_line(n, m, *table_value(n, m)))
    for number, (wemis, zapp_cm, zem_cm, temis_s) in enumerate(bunches, start=1):
        theta4 = math.acos(wemis) * (180.0 / math.pi) * 4.0
        height4 = zapp_cm * 1.0e-5 * 4.0
        n2, m2 = int(theta4), int(height4)
        distance, time = bilinear_reference(theta4, height4, n2, m2)
        corners_d = [
            table_value(n2, m2)[0],
            table_value(n2, m2 + 1)[0],
            table_value(n2 + 1, m2)[0],
            table_value(n2 + 1, m2 + 1)[0],
        ]
        corners_t = [
            table_value(n2, m2)[1],
            table_value(n2, m2 + 1)[1],
            table_value(n2 + 1, m2)[1],
            table_value(n2 + 1, m2 + 1)[1],
        ]
        tofip = time + tofip_offset
        temis_ns = temis_s * 1.0e9
        reals = [
            wemis,
            zapp_cm,
            zem_cm,
            math.sqrt(1 - wemis * wemis),
            distance * rdist_scale,
            tofip,
            temis_ns,
            temis_ns + tofip,
            0.0,
            0.0,
            wemis,
            math.sqrt(1 - wemis * wemis),
            0.0,
            1.0 / 1.0003,
            *corners_d,
            *corners_t,
        ]
        lines.append(bunch_line(1, number, n2, m2, reals))
    path.write_text("\n".join(lines) + "\n")
    return path


def zero_nu(altitude):
    """No refraction: n - 1 = 0 everywhere (the vacuum limit of the mapping)."""
    return np.zeros_like(np.asarray(altitude, dtype=np.float64))


# ---- parser -------------------------------------------------------------------------------


def test_fortran_e_matches_e23_15_layout():
    field = tables.BUNCH_OFFSET  # the first real starts after 34 columns
    assert field == 34
    assert fortran_e(1.5e6) == "0.150000000000000E+07".rjust(23)
    assert float(fortran_e(-2.25e-3)) == pytest.approx(-2.25e-3, rel=1e-14)
    assert len(fortran_e(-2.25e-3)) == 23


def test_parse_trace_reads_constants_rows_and_bunches(tmp_path):
    path = synthetic_trace(
        tmp_path / "t.dat", [(math.cos(math.radians(70.125)), 2.0125e6, 1.5e6, 1e-5)]
    )
    trace = tables.parse_trace(path)
    constants = trace["constants"]
    assert constants["C1_CM"] == pytest.approx(EARTH_CM, rel=1e-14)
    assert constants["C25_CMPS"] == pytest.approx(C_CM_PER_S, rel=1e-14)
    assert constants["OBSLEV_CM"] == OBS_CM
    assert constants["MODATM"] == 1.0
    assert constants["FREFRX"] is True
    assert len(trace["rows"]) == len(tables.ROW_THETA_DEG) * len(tables.ROW_HEIGHT_KM)
    assert trace["rows"][(280, 80)] == pytest.approx(table_value(280, 80), rel=1e-14)
    bunch = trace["bunches"][0]
    assert (bunch["event"], bunch["index"], bunch["n2"], bunch["m2"]) == (1, 1, 280, 80)
    assert bunch["values"][0] == pytest.approx(
        math.cos(math.radians(70.125)), rel=1e-14
    )
    _, printed_time = bilinear_reference(280.5, 80.5, 280, 80)
    assert bunch["values"][5] == pytest.approx(printed_time, rel=1e-14)


def test_parse_rejects_wrong_header(tmp_path):
    path = tmp_path / "bad.dat"
    path.write_text("CURVTAB 2 something\n")
    with pytest.raises(ValueError, match="not a CURVTAB 1 trace"):
        tables.parse_trace(path)


def test_parse_rejects_missing_constant(tmp_path):
    path = tmp_path / "missing.dat"
    lines = [
        "CURVTAB 1 x",
        const_line("C1_CM", EARTH_CM),
        const_line("OBSLEV_CM", OBS_CM),
    ]
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="missing CONST"):
        tables.parse_trace(path)


def test_parse_rejects_short_bunch_and_unknown_line(tmp_path):
    short = tmp_path / "short.dat"
    short.write_text("CURVTAB 1 x\nBUNCH         1     1   280    80 0.1\n")
    with pytest.raises(ValueError, match="too short"):
        tables.parse_trace(short)
    unknown = tmp_path / "unknown.dat"
    unknown.write_text("CURVTAB 1 x\nNOPE 1 2 3\n")
    with pytest.raises(ValueError, match="unknown line"):
        tables.parse_trace(unknown)


# ---- interpolation arithmetic --------------------------------------------------------------


def test_cell_index_follows_distip_truncation():
    wemis = math.cos(math.radians(70.125))
    theta4, height4, n2, m2, obs4 = tables.cell_index([wemis], [2.0125e6], OBS_CM)
    assert theta4[0] == pytest.approx(280.5, rel=1e-13)
    assert height4[0] == pytest.approx(80.5, rel=1e-13)
    assert (int(n2[0]), int(m2[0])) == (280, 80)
    assert obs4 == pytest.approx(8.624, rel=1e-14)


def test_bilinear_hand_computed_example():
    """Corners 100, 110 (N2,M1), 130 (N1,M2), 150 (N1,M1) at theta 280.5 and height 80.5.

    Height weights 0.5: at N2 100 + 0.5*10 = 105; at N1 130 + 0.5*20 = 140; theta weight 0.5:
    105 + 0.5*35 = 122.5. Time corners 1000, 1200, 1100, 1500 give 1200.
    """
    theta4, height4, n2, m2 = (
        np.array([280.5]),
        np.array([80.5]),
        np.array([280]),
        np.array([80]),
    )
    corners = [np.array([v]) for v in (100.0, 110.0, 130.0, 150.0)]
    distip = tables.bilinear(theta4, height4, n2, m2, 8.624, *corners)
    corners = [np.array([v]) for v in (1000.0, 1200.0, 1100.0, 1500.0)]
    tofip = tables.bilinear(theta4, height4, n2, m2, 8.624, *corners)
    assert distip[0] == pytest.approx(122.5, abs=1e-12)
    assert tofip[0] == pytest.approx(1200.0, abs=1e-12)


def test_bilinear_height_clamp_below_observation_level():
    """Heights below the observation level: HM2 = OBS in the weights, as in DISTIP (hand-evaluated).

    m2 = 0 gives hm2 = 8.624 and hm1 = 1, so the height slope is 10 / (1 - 8.624) = -1.311647
    at both theta cells. at N2: -1.311647 * (0.2 - 8.624) + 10 = 10 + 11.049319 = 21.049319;
    at N1: 41.049319; theta weight 0.5: 21.049319 + 0.5 * 20 = 31.049319.
    """
    slope = 10.0 / (1.0 - 8.624)
    at_n2 = slope * (0.2 - 8.624) + 10.0
    at_n1 = slope * (0.2 - 8.624) + 30.0
    expected = 0.5 * (at_n1 - at_n2) + at_n2
    value = tables.bilinear(
        np.array([280.5]),
        np.array([0.2]),
        np.array([280]),
        np.array([0]),
        8.624,
        np.array([10.0]),
        np.array([20.0]),
        np.array([30.0]),
        np.array([40.0]),
    )
    assert value[0] == pytest.approx(expected, abs=1e-12)
    assert value[0] == pytest.approx(31.0493179, abs=1e-7)


def test_bilinear_matches_bruteforce_reference():
    """Above the observation level the two-step interpolation equals explicit bilinear weights."""
    rng = np.random.default_rng(20261115)
    count = 400
    theta4 = rng.uniform(0.0, 356.0, count)
    height4 = rng.uniform(9.0, 440.0, count)
    n2 = np.trunc(theta4).astype(np.int64)
    m2 = np.trunc(height4).astype(np.int64)
    corners = rng.uniform(-1.0e6, 1.0e6, (4, count))
    mine = tables.bilinear(theta4, height4, n2, m2, 8.624, *corners)
    reference = tables.bilinear_bruteforce(
        theta4, height4, n2, m2, corners[0], corners[1], corners[2], corners[3]
    )
    scale = np.max(np.abs(corners), axis=0)
    assert np.max(np.abs(mine - reference) / scale) < 1e-12


def test_consistent_trace_has_no_interpolation_error(tmp_path):
    wemis = math.cos(math.radians(70.125))
    path = synthetic_trace(tmp_path / "ok.dat", [(wemis, 2.0125e6, 1.5e6, 1.0e-5)])
    trace = tables.parse_trace(path)
    columns = tables.bunch_arrays(trace["bunches"])
    check = tables.check_interpolation(columns, OBS_CM)
    assert check["cell_index_mismatch"] == 0
    assert float(np.max(check["distip_relative_error"])) < 1e-12
    assert float(np.max(check["tofip_relative_error"])) < 1e-12
    assert float(np.max(check["cartim_relative_error"])) < 1e-12


def test_altered_printed_distip_is_detected(tmp_path):
    """A printed RDIST 1e-4 off the table interpolation is far above the limit of 1e-6."""
    wemis = math.cos(math.radians(70.125))
    path = synthetic_trace(
        tmp_path / "bad.dat", [(wemis, 2.0125e6, 1.5e6, 1.0e-5)], rdist_scale=1.0 + 1e-4
    )
    columns = tables.bunch_arrays(tables.parse_trace(path)["bunches"])
    check = tables.check_interpolation(columns, OBS_CM)
    assert float(np.max(check["distip_relative_error"])) > 1e-5


def test_cell_index_mismatch_is_counted(tmp_path):
    wemis = math.cos(math.radians(70.125))
    path = synthetic_trace(tmp_path / "idx.dat", [(wemis, 2.0125e6, 1.5e6, 1.0e-5)])
    columns = tables.bunch_arrays(tables.parse_trace(path)["bunches"])
    columns["n2"] = columns["n2"] + 1  # printed index disagrees with the recomputed one
    assert tables.check_interpolation(columns, OBS_CM)["cell_index_mismatch"] == 1


def test_node_check_compares_corners_with_rows(tmp_path):
    # Theta 70 deg (N2 = 280) and height 20 km (M2 = 80) are printed nodes; the offsets keep the
    # truncation away from integer boundaries.
    wemis = math.cos(math.radians(70.0001))
    path = synthetic_trace(tmp_path / "node.dat", [(wemis, 2.00001e6, 1.5e6, 1.0e-5)])
    trace = tables.parse_trace(path)
    columns = tables.bunch_arrays(trace["bunches"])
    compared, worst = tables.node_check(columns, trace["rows"])
    assert compared >= 1
    assert worst < 1e-12


# ---- mapping residuals and binning --------------------------------------------------------


def test_mapping_residuals_vanish_in_the_vacuum_limit(tmp_path):
    """Without refraction, RDIST is the straight horizontal part and TOFIP the straight time.

    obdeect's shift and excess are zero for n - 1 = 0, so both residuals must vanish.
    """
    theta = math.radians(70.125)
    wemis = math.cos(theta)
    zapp_cm, zem_cm = 2.0125e6, 1.5e6
    h_e_m, obs_m = zem_cm * 1e-2, OBS_CM * 1e-2
    path = synthetic_trace(tmp_path / "vac.dat", [(wemis, zapp_cm, zem_cm, 1e-5)])
    trace = tables.parse_trace(path)
    columns = tables.bunch_arrays(trace["bunches"])
    columns["rdist_cm"] = np.array([(h_e_m - obs_m) * math.tan(theta) * 1e2])
    columns["tofip_ns"] = np.array([(h_e_m - obs_m) / wemis / tables.C_M_PER_NS])
    records = tables.mapping_residuals(
        {**columns, "index": np.array([1])},
        trace["constants"],
        zero_nu,
    )
    assert records[0]["shift_residual_m"] == pytest.approx(0.0, abs=1e-9)
    assert records[0]["excess_residual_ns"] == pytest.approx(0.0, abs=1e-9)


def test_bin_residuals_half_open_edges_and_counts():
    records = [
        {
            "theta_s_deg": 65.0,
            "h_e_km": 12.0,
            "below_observation": False,
            "shift_residual_m": 0.1,
            "excess_residual_ns": 1.0,
        },
        {
            "theta_s_deg": 69.9,
            "h_e_km": 14.9,
            "below_observation": False,
            "shift_residual_m": 0.3,
            "excess_residual_ns": 3.0,
        },
        {
            "theta_s_deg": 70.0,
            "h_e_km": 12.0,
            "below_observation": False,
            "shift_residual_m": -0.5,
            "excess_residual_ns": -2.0,
        },
        {
            "theta_s_deg": 90.0,
            "h_e_km": 12.0,
            "below_observation": False,
            "shift_residual_m": 9.0,
            "excess_residual_ns": 9.0,
        },
        {
            "theta_s_deg": 66.0,
            "h_e_km": 1.0,
            "below_observation": True,
            "shift_residual_m": 9.0,
            "excess_residual_ns": 9.0,
        },
    ]
    cells = tables.bin_residuals(records)
    by_key = {(tuple(c["zenith_deg"]), tuple(c["height_km"])): c for c in cells}
    first = by_key[((65.0, 70.0), (10.0, 15.0))]
    assert first["count"] == 2
    assert first["shift_residual_m_mean"] == pytest.approx(0.2)
    assert first["excess_residual_ns_max_abs"] == pytest.approx(3.0)
    second = by_key[((70.0, 75.0), (10.0, 15.0))]
    assert second["count"] == 1
    assert second["shift_residual_m_std"] is None
    assert (
        sum(c["count"] for c in cells) == 3
    )  # 90 deg and the below-observation bunch are in no bin


def test_metrics_and_limits_are_consistent_with_the_case():
    case = runner.validate_case(json.loads(CASE.read_text()))
    produced = tables.summarise(
        [
            {
                "bunch_count": 3,
                "cell_index_mismatch_count": 0,
                "interpolation_max_relative_error": 0.0,
            }
        ]
    )
    assert {check["metric"] for check in case["metrics"]} <= set(produced)
    for check in case["metrics"]:
        assert check["metric"] in tables.LIMITS
    assert case["metrics"][0] == {"metric": "file_count", "expected": 2, "atol": 0}
    assert set(tables.summarise([])) == set(produced)
    assert "requires" not in case


def test_bunch_format_matches_the_patch():
    text = PATCH.read_text()
    assert "'BUNCH',I10,I7,2I6,22(1X,E23.15)" in text
    assert tables.BUNCH_FIELDS == 22
    assert "OBDEECT_CURVTAB" in text and "curvature-tables.dat" in text


def test_analyse_file_reports_bunches_and_metrics(tmp_path):
    wemis = math.cos(math.radians(70.125))
    path = synthetic_trace(
        tmp_path / "curvature-tables-65deg.dat",
        [
            (wemis, 2.0125e6, 1.5e6, 1e-5),
            (math.cos(math.radians(72.6)), 1.8125e6, 1.2e6, 2e-5),
        ],
    )
    row = tables.analyse_file(path, 65, zero_nu)
    assert row["bunch_count"] == 2
    assert row["row_count"] == 32
    assert row["cell_index_mismatch_count"] == 0
    assert row["interpolation_max_relative_error"] < 1e-12
    assert len(row["bunches"]) == 2
    metrics = tables.summarise([row])
    assert metrics["bunch_count_positive"] == 1
    assert metrics["file_count"] == 1


# ---- inputs and the patch -----------------------------------------------------------------


def test_main_exits_3_when_a_trace_is_missing(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("VALIDATION_INPUTS", str(tmp_path))
    assert tables.main() == 3
    assert "missing input" in capsys.readouterr().err


def _corsika_source():
    candidates = []
    if os.environ.get("CORSIKA_SOURCE"):
        candidates.append(Path(os.environ["CORSIKA_SOURCE"]))
    candidates.append(Path("/workdir/simulation_software/corsika7/src/corsika.F"))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def test_patch_applies_to_corsika_source(tmp_path):
    source = _corsika_source()
    if source is None:
        pytest.skip("corsika.F not available (set CORSIKA_SOURCE to check the patch)")
    patch = shutil.which("patch")
    if patch is None:
        pytest.skip("patch not installed")
    shutil.copy(source, tmp_path / "corsika.F")
    result = subprocess.run(
        [patch, "--dry-run", "-p0", "-i", str(PATCH)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
