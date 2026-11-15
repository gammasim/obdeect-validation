"""Cards and header checks of the matrix shower files (no CORSIKA run needed)."""

import math
from pathlib import Path

import pytest

from obdeect_validation.adapters import chain_showers
from obdeect_validation.adapters import shower_matrix as sm

FIXTURE_CARD = (
    Path(__file__).resolve().parents[1]
    / "fixtures/corsika/gamma-20deg-100-3000GeV.card"
)


def _keys(text):
    return [line.split(" ", 1)[0] for line in text.splitlines() if line.strip()]


def _values(text, key):
    return [line.split() for line in text.splitlines() if line.startswith(key + " ")]


def _header(
    shower, *, energy=None, zenith=None, azimuth=None, primary=None, obslev=None
):
    """EVTH tuple as eventio returns it (element i is word i + 1)."""
    theta = math.radians(sum(shower.zenith) / 2 if zenith is None else zenith)
    phi = math.radians(sum(shower.azimuth) / 2 if azimuth is None else azimuth)
    head = [0.0] * 30
    head[0] = b"EVTH"
    head[2] = float(shower.primary if primary is None else primary)
    head[3] = float(sum(shower.energy) / 2 if energy is None else energy)
    head[10] = theta
    head[11] = phi
    level = sm.SITES[shower.site].obslev_cm if obslev is None else obslev
    head[18] = [float(level)]
    return tuple(head)


def test_specs_cover_the_matrix_points():
    assert set(sm.SPECS) == {
        "gamma-20deg-100-3000GeV",
        "gamma-40deg-100-3000GeV",
        "gamma-60deg-300-5000GeV",
        "gamma-20deg-az180-100-3000GeV",
        "proton-20deg-300-10000GeV",
        "gamma-20deg-South-100-3000GeV",
        "gamma-20deg-100-3000GeV-r1300",
        "gamma-20deg-100-3000GeV-array-north",
    }


def test_seeds_are_explicit_and_unique_per_file():
    bases = [shower.seed_base for shower in sm.SPECS.values()]
    assert len(set(bases)) == len(bases)
    for shower in sm.SPECS.values():
        seeds = _values(sm.card_text(shower), "SEED")
        assert len(seeds) == 4
        assert [int(line[1]) for line in seeds] == [
            shower.seed_base + offset for offset in sm.SEED_OFFSETS
        ]


def test_telescopes_match_the_checked_in_card():
    existing = _values(FIXTURE_CARD.read_text(), "TELESCOPE")
    for shower in sm.SPECS.values():
        if shower.telescopes != chain_showers.FLAT_TELESCOPES:
            continue  # the 1300 cm and array files are checked below
        generated = _values(sm.card_text(shower), "TELESCOPE")
        assert [[float(v) for v in line[1:]] for line in generated] == [
            [float(v) for v in line[1:]] for line in existing
        ]


@pytest.mark.parametrize(
    ("name", "primary", "zenith", "azimuth", "obslev", "atmosphere"),
    [
        ("gamma-40deg-100-3000GeV", "1", "40 40", "0 0", "215600", "1 Y"),
        ("gamma-60deg-300-5000GeV", "1", "60 60", "0 0", "215600", "1 Y"),
        ("gamma-20deg-az180-100-3000GeV", "1", "20 20", "180 180", "215600", "1 Y"),
        ("proton-20deg-300-10000GeV", "14", "19 21", "0 0", "215600", "1 Y"),
        ("gamma-20deg-South-100-3000GeV", "1", "20 20", "0 0", "214700", "99 Y"),
    ],
)
def test_card_lines(name, primary, zenith, azimuth, obslev, atmosphere):
    text = sm.card_text(sm.SPECS[name])
    assert _values(text, "PRMPAR") == [["PRMPAR", primary]]
    assert _values(text, "THETAP") == [["THETAP", *zenith.split()]]
    assert _values(text, "PHIP") == [["PHIP", *azimuth.split()]]
    assert _values(text, "OBSLEV") == [["OBSLEV", obslev]]
    assert _values(text, "ATMOSPHERE") == [["ATMOSPHERE", *atmosphere.split()]]
    assert f"TELFIL {name}.iact" in text
    assert "NSHOW" in text and "CSCAT 5 10000" in text


def test_card_keys_are_the_checked_in_card_keys():
    generated = set(_keys(sm.card_text(sm.SPECS["gamma-40deg-100-3000GeV"])))
    existing = set(_keys(FIXTURE_CARD.read_text()))
    assert existing <= generated


def test_site_profiles():
    assert sm.SITES["North"].atmosphere == 1
    assert sm.SITES["South"].atmosphere == 99
    assert sm.SITES["South"].profile_name == "atmprof99.dat"
    assert sm.SITES["South"].obslev_cm == 214700.0


def test_check_log_accepts_flat_site_profile():
    shower = sm.SPECS["gamma-20deg-South-100-3000GeV"]
    text = "Atmospheric profile 99 with 55 levels read from file atmprof99.dat\n"
    sm.check_log(text, shower)


def test_check_log_rejects_wrong_profile_and_curved_banner():
    shower = sm.SPECS["gamma-20deg-South-100-3000GeV"]
    with pytest.raises(ValueError, match="atmosphere line"):
        sm.check_log(
            "Atmospheric profile 1 with 50 levels read from file atmprof1.dat", shower
        )
    with pytest.raises(ValueError, match="CURVED"):
        sm.check_log(sm.CURVED_PATTERN + "\n" + "Atmospheric profile 99 with x", shower)


def test_header_check_passes_matching_headers():
    shower = sm.SPECS["gamma-60deg-300-5000GeV"]
    headers = [_header(shower) for _ in range(3)]
    assert set(sm.header_mismatches(headers, shower).values()) == {0}


def test_header_check_counts_each_mismatch():
    shower = sm.SPECS["gamma-20deg-az180-100-3000GeV"]
    bad = [
        _header(shower, primary=14),
        _header(shower, energy=50.0),
        _header(shower, zenith=30.0),
        _header(shower, azimuth=0.0),
        _header(shower, obslev=214700.0),
    ]
    counts = sm.header_mismatches(bad, shower)
    assert counts == {
        "particle": 1,
        "energy": 1,
        "zenith": 1,
        "azimuth": 1,
        "observation_level": 1,
    }


def test_azimuth_window_is_circular():
    shower = sm.SPECS["gamma-20deg-az180-100-3000GeV"]
    assert (
        sm.header_mismatches([_header(shower, azimuth=179.9999)], shower)["azimuth"]
        == 0
    )
    shower = sm.SPECS["gamma-40deg-100-3000GeV"]
    assert (
        sm.header_mismatches([_header(shower, azimuth=359.9999)], shower)["azimuth"]
        == 0
    )


def test_proton_zenith_window_accepts_cone():
    shower = sm.SPECS["proton-20deg-300-10000GeV"]
    headers = [_header(shower, zenith=z) for z in (19.0, 20.0, 21.0)]
    assert sm.header_mismatches(headers, shower)["zenith"] == 0
    assert sm.header_mismatches([_header(shower, zenith=22.0)], shower)["zenith"] == 1


def test_r1300_card_is_one_lst_sphere_of_1300_cm():
    text = sm.card_text(sm.SPECS["gamma-20deg-100-3000GeV-r1300"])
    assert _values(text, "TELESCOPE") == [["TELESCOPE", "0", "0", "0", "1300"]]
    assert _values(text, "CSCAT") == [["CSCAT", "5", "10000", "0.0"]]


def test_array_card_lists_lsts_and_msts_in_simtools_units():
    shower = sm.SPECS["gamma-20deg-100-3000GeV-array-north"]
    lines = _values(sm.card_text(shower), "TELESCOPE")
    assert len(lines) == 6
    assert [[float(v) for v in line[1:]] for line in lines] == [
        list(row) for row in sm.ARRAY_NORTH_TELESCOPES
    ]
    # LSTN-01 position (-70.91, -52.35, 45.0) m as written by simtools (cm)
    assert lines[0][1:4] == ["-7091", "-5235", "4500"]
    assert _values(sm.card_text(shower), "CSCAT") == [["CSCAT", "5", "25000", "0.0"]]


def test_array_spheres_contain_the_mirror_reach():
    # LST mirror reach 1228.8 cm (LSTN-01 mirror_list 2.0.0); MST reach 680.8 cm (MSTx-NectarCam).
    radii = {row[3] for row in sm.ARRAY_NORTH_TELESCOPES}
    assert radii == {1300.0, 915.0}
    lsts, msts = sm.ARRAY_NORTH_TELESCOPES[:4], sm.ARRAY_NORTH_TELESCOPES[4:]
    assert all(row[3] >= 1228.8 for row in lsts)
    assert all(row[3] >= 680.8 for row in msts)
    # every telescope lies inside the CSCAT disc with its sphere
    for x, y, _z, r in sm.ARRAY_NORTH_TELESCOPES:
        assert (x * x + y * y) ** 0.5 + r <= sm.ARRAY_NORTH_SCATTER_CM
