# ABOUTME: `lunaNMR series` processes spectra in series order - by titration point or delay.
# ABOUTME: Real DNAJA1 titration: a cascade seeded in name order tracks every peak from the wrong point.

"""
The CLI sorted the files by name, splitting digit runs, so 0o5 < 0o25 < 0o125 and the
cascade ran 0 -> 0.5 -> 0.25 -> 0.125 eq: each point seeded from the wrong previous one,
territories and linewidth caps sized from it. The GUI sorts by the parsed value.
"""

import json
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


def test_titration_points_are_ordered_by_value():
    from lunaNMR.cli import _discover_spectra, _series_order
    files = _series_order(_discover_spectra(str(DNAJA1)), 'titration')
    assert [Path(f).name for f in files] == [
        "2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o125.ft", "2D_DNAJA1_HSPA8_0o25.ft",
        "2D_DNAJA1_HSPA8_0o5.ft", "2D_DNAJA1_HSPA8_1o0.ft", "2D_DNAJA1_HSPA8_2o0.ft",
        "2D_DNAJA1_HSPA8_3o0.ft"]


def test_delays_are_ordered_by_time_whatever_their_unit():
    from lunaNMR.cli import _series_order
    names = ["T1_1s.ft", "T1_500ms.ft", "T1_20ms.ft"]
    assert _series_order(names, 'time') == ["T1_20ms.ft", "T1_500ms.ft", "T1_1s.ft"]


def test_a_spectrum_without_a_value_is_kept_after_the_others():
    from lunaNMR.cli import _series_order
    assert _series_order(["T1_500ms.ft", "reference.ft", "T1_20ms.ft"], 'time') == [
        "T1_20ms.ft", "T1_500ms.ft", "reference.ft"]


def test_the_series_run_cascades_in_titration_order(tmp_path):
    from lunaNMR.cli import main
    (tmp_path / "spectra").mkdir()
    for name in ("2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o5.ft", "2D_DNAJA1_HSPA8_0o25.ft"):
        (tmp_path / "spectra" / name).symlink_to(DNAJA1 / name)
    assert main(["series", "--spectra", str(tmp_path / "spectra"),
                 "--peaks", str(DNAJA1 / "DNAJA1.txt"), "--mode", "titration",
                 "--peak-source", "cascade", "--out", str(tmp_path / "out")]) == 0
    columns = json.loads((tmp_path / "out" / "series_metadata.json").read_text())["columns"]
    assert [c["value"] for c in columns] == [0.0, 0.25, 0.5]
