"""Zero-optical-depth transmission table: bottom (H2) below every emission height of the cases."""

import pytest

from obdeect_validation.adapters import common

TABLE = """#============================
# End altitude: 2.156 km
#============================

# H2= 2.156, H1=    2.206     2.256     2.356
200 0.1 0.2 0.3
999 0.4 0.5 0.6
1000 0.7 0.8 0.9
"""

LOWEST_CORSIKA_EMISSION_M = (
    2147.012  # South site, gamma-20deg-South-100-3000GeV.iact (zem scan)
)
LOWEST_GROUND_SOURCE_M = 2147.5  # night-sky source at 2150 m (lowest emission 2147.5 m, nsb-camera docstring)


def test_bottom_is_below_all_validation_emission_heights():
    assert common.ZERO_DEPTH_BOTTOM_M < LOWEST_CORSIKA_EMISSION_M
    assert common.ZERO_DEPTH_BOTTOM_M < LOWEST_GROUND_SOURCE_M


def test_zero_depth_table_sets_bottom_and_zeroes_every_column():
    lines = common.zero_depth_table(TABLE).splitlines()
    assert lines[0] == "#============================"
    assert lines[4] == "# H2= 2.140, H1=    2.206     2.256     2.356"
    assert lines[-3:] == ["200 0 0 0", "999 0 0 0", "1000 0 0 0"]


def test_zero_depth_table_keeps_comment_lines_and_row_count():
    source = TABLE.splitlines()
    result = common.zero_depth_table(TABLE).splitlines()
    assert len(result) == len(source)
    assert result[1] == source[1]


def test_zero_depth_table_requires_h2_header():
    with pytest.raises(ValueError, match="H2"):
        common.zero_depth_table("200 0.1 0.2\n")


def test_write_zero_depth_table_uses_source_and_bottom(tmp_path):
    source = tmp_path / "table.dat"
    source.write_text(TABLE, encoding="utf-8")
    path = common.write_zero_depth_table(
        tmp_path, source, bottom_m=2100.0, name="zero.dat"
    )
    assert path == tmp_path / "zero.dat"
    text = path.read_text(encoding="utf-8")
    assert "H2= 2.100" in text
    assert "H2= 2.156" not in text
    assert text.splitlines()[-1] == "1000 0 0 0"


def test_write_zero_depth_table_default_bottom_is_zero_depth_bottom(tmp_path):
    source = tmp_path / "table.dat"
    source.write_text(TABLE, encoding="utf-8")
    path = common.write_zero_depth_table(tmp_path, source)
    assert path.name == "unit-atmosphere.dat"
    assert f"H2= {common.ZERO_DEPTH_BOTTOM_M / 1000.0:.3f}" in path.read_text(
        encoding="utf-8"
    )
