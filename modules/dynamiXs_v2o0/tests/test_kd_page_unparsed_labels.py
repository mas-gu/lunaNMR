# ABOUTME: The Kd page leaves [L] empty for a spectrum whose label names no titration point.
# ABOUTME: Real DNAJA1 spectra under index names, run through `series`, opened in the real page.

"""
`series` names a spectrum it cannot read a titration point from after its file stem, and
records `value: null` for it. The page prefilled [L] with the stem's trailing digits, so
`sample_2` became 2 - for a spectrum acquired at 3 eq. `kd` refuses such labels as
concentrations (kd_input.labels_without_points); the page judges them the same way.
"""

import contextlib
import io
import json
import os
import sys
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_MOD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MOD))
sys.path.insert(0, str(_MOD / "dynamiXs_Kd"))

DNAJA1 = Path(__file__).resolve().parents[4] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
# Spectrum acquired at each point -> the index name `series` cannot read a point from.
INDEX_NAMES = {"0o0": "sample_0.ft", "0o5": "sample_1.ft", "3o0": "sample_2.ft"}
LABELS = ["sample_0", "sample_1", "sample_2"]
EQUIVALENTS = [0.0, 0.5, 3.0]

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def index_named_tidy(tmp_path_factory):
    from lunaNMR.cli import main
    root = tmp_path_factory.mktemp("index_named")
    (root / "spectra").mkdir()
    for point, name in INDEX_NAMES.items():
        (root / "spectra" / name).symlink_to(next(DNAJA1.glob(f"2D_*_{point}.ft")))
    with contextlib.redirect_stdout(io.StringIO()):
        code = main(["series", "--spectra", str(root / "spectra"),
                     "--peaks", str(DNAJA1 / "DNAJA1.txt"), "--mode", "titration",
                     "--out", str(root / "out"), "--format", "json"])
    assert code == 0
    meta = json.loads((root / "out" / "series_metadata.json").read_text())
    assert meta["n_value_unparsed"] == 3, "the premise: series could not parse these names"
    return root / "out" / "series_analysis_tidy.csv"


def _page():
    from kd_titration_page import KdTitrationPage
    return KdTitrationPage(types.SimpleNamespace(show_main_menu=lambda: None,
                                                 current_dir=""))


def _empty(page):
    return [page._conc_is_empty(spin) for spin in page.conc_spins]


class TestLabelsThatNameNoPointLeaveTheirRowEmpty:

    def test_no_row_is_prefilled(self, app, index_named_tidy):
        page = _page()
        page._set_input_file(str(index_named_tidy))
        assert len(page.conc_spins) == 3
        assert _empty(page) == [True, True, True]

    def test_the_page_names_them_and_says_what_to_do(self, app, index_named_tidy):
        page = _page()
        page._set_input_file(str(index_named_tidy))
        assert page.result_label.isVisibleTo(page)
        text = page.result_label.text()
        for label in LABELS:
            assert label in text
        assert "typed or loaded" in text
        assert all(label in page.log_text.toPlainText() for label in LABELS)

    def test_run_is_refused_on_the_page(self, app, index_named_tidy):
        page = _page()
        page._set_input_file(str(index_named_tidy))
        page.obs_combo.setCurrentIndex(2)
        page.conc_units_combo.setCurrentIndex(1)          # equivalents
        page._start_analysis()
        assert page.worker is None
        assert not (index_named_tidy.parent / "kd_analysis").exists()
        assert "sample_2" in page.result_label.text()
        assert "sample_2" in page.log_text.toPlainText().splitlines()[-1]

    def test_one_row_left_empty_still_refuses(self, app, index_named_tidy):
        page = _page()
        page._set_input_file(str(index_named_tidy))
        page.conc_spins[0].setValue(0.0)
        page.conc_spins[1].setValue(0.5)
        page.conc_units_combo.setCurrentIndex(1)
        page._start_analysis()
        assert page.worker is None
        assert "sample_2" in page.result_label.text()
        assert "sample_0" not in page.result_label.text()

    def test_typed_concentrations_fit(self, app, index_named_tidy, tmp_path):
        page = _page()
        page._set_input_file(str(index_named_tidy))
        page._set_output_dir(str(tmp_path / "kd"))
        for spin, value in zip(page.conc_spins, EQUIVALENTS):
            spin.setValue(value)
        assert _empty(page) == [False, False, False]
        page.obs_combo.setCurrentIndex(2)
        page.conc_units_combo.setCurrentIndex(1)
        page._start_analysis()
        assert page.worker is not None
        assert page.worker.wait(300_000), "the fit did not finish"
        app.processEvents()
        meta = json.loads(Path(page.last_json_file).read_text())["metadata"]
        assert meta["concentrations"] == [0.0, 25.0, 150.0]
        assert meta["conc_units"] == "absolute"

    def test_a_params_file_for_this_input_fills_them(self, app, index_named_tidy, tmp_path):
        csv = tmp_path / index_named_tidy.name
        csv.write_bytes(index_named_tidy.read_bytes())
        (tmp_path / "series_metadata.json").write_bytes(
            (index_named_tidy.parent / "series_metadata.json").read_bytes())
        (tmp_path / "x_kd_params.json").write_text(json.dumps(
            {"input": csv.name, "points": [0.0, 1.0, 2.0],
             "concentrations": EQUIVALENTS, "conc_units": "equivalents"}))
        page = _page()
        page._set_input_file(str(csv))
        assert [s.value() for s in page.conc_spins] == EQUIVALENTS
        assert not page.result_label.isVisibleTo(page)
