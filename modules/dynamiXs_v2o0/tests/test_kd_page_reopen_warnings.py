# ABOUTME: A reopened Kd fit states the verdict its run stated: the quality warnings stored
# ABOUTME: in its fit JSON appear on the page, and a stale outcome from another fit does not.

import json
import os
import shutil
import sys
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_MOD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MOD))
sys.path.insert(0, str(_MOD / "dynamiXs_Kd"))

DNAJA1 = (Path(__file__).resolve().parents[4] / "data_test" / "data_binding" / "tidy"
          / "DNAJA1_HSPA8_tidy.csv")

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _page():
    from kd_titration_page import KdTitrationPage
    return KdTitrationPage(types.SimpleNamespace(show_main_menu=lambda: None,
                                                 current_dir=""))


def _saved_fit(app, folder, equivalents):
    """Fit DNAJA1 on a page with the real worker; return the page and the fit as a saved
    analysis."""
    folder.mkdir()
    csv = folder / DNAJA1.name
    shutil.copy(DNAJA1, csv)
    page = _page()
    page._set_input_file(str(csv))
    if equivalents:
        page.conc_units_combo.setCurrentIndex(1)
    page._start_analysis()
    assert page.worker.wait(300_000), "the fit did not finish"
    app.processEvents()
    return page, {'fit_data_path': page.last_json_file, 'meta': page.get_session_state()}


def test_the_run_and_the_reopen_state_the_same_verdict(app, tmp_path):
    """Labels 0-3 equivalents read as micromolar: a Kd outside what the titration
    resolves. The run says so; reopening the saved fit must too."""
    _, entry = _saved_fit(app, tmp_path / "labels", equivalents=False)
    warnings_ = json.loads(Path(entry['fit_data_path']).read_text())[
        'metadata']['quality_warnings']
    assert warnings_

    page = _page()
    page.open_saved_analysis(entry, name="labels")

    assert page.result_label.isVisibleTo(page)
    assert "not usable" in page.result_label.text()
    assert warnings_[0] in page.result_label.text()


def test_a_clean_fit_does_not_inherit_another_fits_verdict(app, tmp_path):
    """The page that just ran an unusable fit opens a clean one: the warning shown is
    the run's, not the reopened fit's."""
    page, _ = _saved_fit(app, tmp_path / "labels", equivalents=False)
    assert page.result_label.isVisibleTo(page)
    _, good = _saved_fit(app, tmp_path / "equivalents", equivalents=True)

    page.open_saved_analysis(good, name="equivalents")

    assert not page.result_label.isVisibleTo(page)
