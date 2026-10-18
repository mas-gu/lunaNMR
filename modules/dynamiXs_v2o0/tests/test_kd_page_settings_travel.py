# ABOUTME: The Kd page carries a dataset's settings the way the CLI does: units with their
# ABOUTME: concentrations, and stored thresholds into the fit. Real DNAJA1 titration, real worker.

import contextlib
import io
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

TIDY = Path(__file__).resolve().parents[4] / "data_test" / "data_binding" / "tidy"
DNAJA1 = TIDY / "DNAJA1_HSPA8_tidy.csv"
EQUIVALENTS = [0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 3.0]

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture
def csv(tmp_path):
    path = tmp_path / "DNAJA1_HSPA8_tidy.csv"
    shutil.copy(DNAJA1, path)
    return path


def _page():
    from kd_titration_page import KdTitrationPage
    return KdTitrationPage(types.SimpleNamespace(show_main_menu=lambda: None,
                                                 current_dir=""))


def _run(app, page):
    page._start_analysis()
    assert page.worker.wait(300_000), "the fit did not finish"
    app.processEvents()
    return json.loads(Path(page.last_json_file).read_text())


def _cli_kd(csv, *args):
    """`kd` on `csv`, writing where it writes by default for a dataset folder."""
    from lunaNMR.cli import main
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = main(["kd", "--format", "json", "--p0", "50", "--input", str(csv),
                     "--out", str(csv.parent / "kd_analysis"), *args])
    assert code == 0, out.getvalue()
    summary = json.loads(out.getvalue())
    return json.loads(Path(summary.get("result", summary)["json_file"]).read_text())


class TestTheUnitsTravelWithTheConcentrations:
    """The page's state holds the [L] values; without their units a reopened page reads
    equivalents as micromolar and refits a Kd wrong by the factor [P]0."""

    def test_a_restored_page_keeps_equivalents(self, app, csv):
        page = _page()
        page._set_input_file(str(csv))
        page.conc_units_combo.setCurrentIndex(1)

        reopened = _page()
        reopened.restore_session_state(page.get_session_state())
        assert reopened._conc_units() == "equivalents"

    def test_a_restored_page_keeps_absolute(self, app, csv):
        page = _page()
        page._set_input_file(str(csv))

        reopened = _page()
        reopened.conc_units_combo.setCurrentIndex(1)
        reopened.restore_session_state(page.get_session_state())
        assert reopened._conc_units() == "absolute"

    def test_a_refit_after_reopening_gives_the_same_kd(self, app, csv):
        page = _page()
        page._set_input_file(str(csv))
        page.obs_combo.setCurrentIndex(2)                 # intensity only
        page.conc_units_combo.setCurrentIndex(1)
        first = _run(app, page)["global"]["intensity"]["Kd"]

        reopened = _page()
        reopened.restore_session_state(page.get_session_state())
        again = _run(app, reopened)["global"]["intensity"]["Kd"]
        assert again == pytest.approx(first)
