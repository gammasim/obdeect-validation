"""Parameters of the optics chain adapter: pointing, camera class, Bypass_Optics, default command text."""

import importlib
import math

import pytest

pytest.importorskip("eventio")
chain_optics = importlib.import_module("obdeect_validation.adapters.chain_optics")


def test_deg_text_keeps_integer_pointing_of_the_default_chain():
    assert chain_optics.deg_text(20.0) == "20"
    assert chain_optics.deg_text(180) == "180"
    assert chain_optics.deg_text(0.0) == "0"
    assert chain_optics.deg_text(22.5) == "22.5"
    assert chain_optics.deg_text(math.degrees(0.34906584)) == "20"


def test_camera_class_maps_sim_telarray_mirror_class():
    single = "mirror_class = 0\n"
    dual = "mirror_class = 2\n"
    assert chain_optics.camera_class(single) == 1
    assert chain_optics.camera_class(dual) == 2
    assert chain_optics.camera_class("camera_config_file = x\n") == 1


def test_camera_map_text_of_the_default_chain_is_unchanged(tmp_path):
    path = chain_optics.camera_map(tmp_path, 10.893)
    assert path.read_text() == (
        "# telescope_id mirror_class camera_rotation_deg\n1 1 10.893\n"
    )
    dual = chain_optics.camera_map(tmp_path, 0.0, mirror_class=2)
    assert dual.read_text().splitlines()[1] == "1 2 0.0"


def test_bypass_optics_is_read_from_the_conversion_log(tmp_path):
    log = tmp_path / "to-simtel.log"
    log.write_text(
        "obdeect stream to-simtel: telescope 1: camera frame, mirror class 2, Bypass_Optics=2, "
        "camera rotation 0 deg\n"
    )
    assert chain_optics.bypass_from_log(log) == 2
    log.write_text("no value here\n")
    with pytest.raises(RuntimeError):
        chain_optics.bypass_from_log(log)


def _header(**overrides):
    summary = {
        "zenith_min_deg": 20.0,
        "zenith_max_deg": 20.0,
        "zenith_mean_deg": 20.0,
        "n_showers": 3,
        "telescope_azimuth_deg": 180.0,
        "header_azimuth_deg": 0.0,
    }
    summary.update(overrides)
    return summary


def test_pointing_derives_zenith_mean_and_n_to_e_azimuth(monkeypatch):
    monkeypatch.setattr(chain_optics, "header_pointing", lambda source: _header())
    zenith, azimuth, header = chain_optics.pointing("x.iact")
    assert (zenith, azimuth) == ("20", "180")
    assert header["n_showers"] == 3


def test_pointing_azimuth_follows_header_plus_180(monkeypatch):
    monkeypatch.setattr(
        chain_optics,
        "header_pointing",
        lambda source: _header(telescope_azimuth_deg=0.0, header_azimuth_deg=180.0),
    )
    assert chain_optics.pointing("x.iact")[1] == "0"
    with pytest.raises(ValueError, match="N->E"):
        chain_optics.pointing("x.iact", azimuth=180.0)


def test_pointing_zenith_must_lie_in_the_shower_range(monkeypatch):
    monkeypatch.setattr(
        chain_optics,
        "header_pointing",
        lambda source: _header(
            zenith_min_deg=19.2, zenith_max_deg=21.0, zenith_mean_deg=20.11
        ),
    )
    assert chain_optics.pointing("x.iact", zenith=20.0)[0] == "20"
    with pytest.raises(ValueError, match="outside the shower zenith range"):
        chain_optics.pointing("x.iact", zenith=40.0)


def test_factor_record_keeps_the_bypass_one_text_by_default(tmp_path):
    (tmp_path / "camera.dat").write_text("Rotate 10.893\n")
    cfg = (
        "camera_config_file = camera.dat\n"
        "mirror_degraded_reflection = 1.0\n"
        "mirror_class = 0\n"
    )
    record = chain_optics.factor_record(cfg, tmp_path)
    assert record["mirror_degraded_reflection"]["candidate"].endswith(
        "sim_telarray Bypass_Optics=1 only warns"
    )
    assert record["camera_rotation_deg"] == pytest.approx(10.893)
    dual = chain_optics.factor_record(cfg, tmp_path, bypass=2)
    assert dual["mirror_degraded_reflection"]["candidate"].endswith(
        "sim_telarray Bypass_Optics=2 only warns"
    )


class _FakeStream:
    def __init__(self, numbers):
        self.numbers = numbers

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def __iter__(self):
        return iter([type("Event", (), {"event_number": n})() for n in self.numbers])


def test_camera_keys_are_100_times_shower_plus_reuse(monkeypatch):
    # two showers with three reuses, one with two (file order, reuses consecutive)
    numbers = [1, 1, 1, 2, 2, 2, 7, 7]
    monkeypatch.setattr(chain_optics, "IACTFile", lambda source: _FakeStream(numbers))
    assert chain_optics.camera_keys("x.iact") == [
        "100",
        "101",
        "102",
        "200",
        "201",
        "202",
        "700",
        "701",
    ]


def test_mcnemar_exact_p_values():
    assert chain_optics.mcnemar_exact_p(0, 0) == 1.0
    assert chain_optics.mcnemar_exact_p(1, 0) == 1.0
    assert chain_optics.mcnemar_exact_p(5, 5) == 1.0
    assert chain_optics.mcnemar_exact_p(0, 10) == pytest.approx(2 / 1024)
    assert chain_optics.mcnemar_exact_p(17, 0) < chain_optics.TRIGGER_P_LIMIT


def test_trigger_test_on_synthetic_mapping(monkeypatch):
    monkeypatch.setattr(chain_optics, "TELESCOPES", (1,))
    monkeypatch.setattr(
        chain_optics, "IACTFile", lambda source: _FakeStream([1, 1, 2, 2])
    )
    keys = chain_optics.camera_keys("x.iact")
    reference = {k: {"telescope_events": {"1": {}}} for k in ("100", "101", "200")}
    candidate = {
        k: {"telescope_events": {"1": {}}} for k in ("100", "101", "200", "201")
    }
    row = chain_optics.trigger_difference(reference, candidate, keys)[1]
    # 201 triggers only in the candidate; 200 and 100, 101 agree; keys cover 4 showers
    assert (row["both"], row["candidate_only"], row["reference_only"]) == (3, 1, 0)
    assert row["difference"] == pytest.approx(0.25)
    assert row["flag"] == 0
    # a systematic excess of candidate-only triggers over ten showers is flagged
    many = [1] * 5 + [2] * 5
    monkeypatch.setattr(chain_optics, "IACTFile", lambda source: _FakeStream(many))
    keys = chain_optics.camera_keys("x.iact")
    reference = {k: {"telescope_events": {}} for k in keys}
    candidate = {k: {"telescope_events": {"1": {}}} for k in keys}
    row = chain_optics.trigger_difference(reference, candidate, keys)[1]
    assert row["candidate_only"] == 10 and row["flag"] == 1


def test_triggered_accepts_integer_and_string_telescope_keys():
    # sim_telarray writes telescope_events keyed by the 1-based telescope number (integer)
    events = {"100": {"telescope_events": {1: {}}}, "101": {"telescope_events": {}}}
    assert chain_optics.triggered(events, "100", 1)
    assert chain_optics.triggered({"100": {"telescope_events": {"1": {}}}}, "100", 1)
    assert not chain_optics.triggered(events, "101", 1)
    assert not chain_optics.triggered(events, "999", 1)
