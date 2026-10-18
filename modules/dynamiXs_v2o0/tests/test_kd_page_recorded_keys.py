# ABOUTME: The Kd page applies only the settings a params source records, leaving every other
# ABOUTME: control as the user set it. Real DNAJA1 titration, files written by the real CLI.

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
    path = tmp_path / DNAJA1.name
    shutil.copy(DNAJA1, path)
    return path


def _page():
    from kd_titration_page import KdTitrationPage
    return KdTitrationPage(types.SimpleNamespace(show_main_menu=lambda: None,
                                                 current_dir=""))


def _cli_survey(csv):
    """`kd --survey` with non-default P0/alpha/observable; leaves its settings file beside
    the input. Returns that file."""
    from lunaNMR.cli import main
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        code = main(["kd", "--survey", "--format", "json", "--input", str(csv),
                     "--out", str(csv.parent / "survey"), "--prefix", "DNAJA1",
                     "--p0", "30", "--alpha", "0.2", "--observable", "intensity",
                     "--conc-units", "equivalents", "--kd-outlier-z", "2"])
    assert code == 0, out.getvalue()
    return csv.parent / "DNAJA1_kd_params.json"


def _set_controls(page):
    """Every control a settings file does not record, away from its default."""
    page.p0_spin.setValue(20.0)
    page.alpha_spin.setValue(0.3)
    page.obs_combo.setCurrentIndex(1)                 # CSP only
    page.intvalue_combo.setCurrentText("volume")
    page.boot_spin.setValue(100)


def _assert_controls_kept(page):
    assert page.p0_spin.value() == 20.0
    assert page.alpha_spin.value() == 0.3
    assert page.obs_combo.currentIndex() == 1
    assert page.intvalue_combo.currentText() == "volume"
    assert page.boot_spin.value() == 100


class TestASettingsFileLeavesWhatItDoesNotRecord:

    def test_the_source_found_on_load(self, app, csv):
        _cli_survey(csv)
        page = _page()
        _set_controls(page)
        page._set_input_file(str(csv))
        _assert_controls_kept(page)
        assert page._conc_units() == "equivalents"
        assert page.get_session_state()["thresholds"]["kd_outlier_z"] == 2.0

    def test_a_file_loaded_by_hand(self, app, csv):
        page = _page()
        page._set_input_file(str(csv))                    # no settings file yet
        path = _cli_survey(csv)
        _set_controls(page)
        page._load_params_file(str(path))
        _assert_controls_kept(page)
        assert page._conc_units() == "equivalents"

    def test_thresholds_it_does_not_record_are_kept(self, app, csv, tmp_path):
        page = _page()
        page._set_input_file(str(csv))
        page.thresholds["ref_max_ratio"] = 3.0
        path = tmp_path / "only_z.json"
        path.write_text(json.dumps({"kd_outlier_z": 2.0}))
        page._load_params_file(str(path))
        assert page.thresholds["ref_max_ratio"] == 3.0
        assert page.thresholds["kd_outlier_z"] == 2.0


class TestAFitJsonStillCarriesItsRecord:
    """A *_kd_fit_data.json records the fit's P0, alpha and observables in its metadata:
    loading one applies them, because it records them."""

    def test_its_recorded_controls_are_applied(self, app, csv):
        from lunaNMR.cli import main
        out = io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
            code = main(["kd", "--format", "json", "--input", str(csv),
                         "--out", str(csv.parent / "fit"), "--p0", "30", "--alpha", "0.2",
                         "--observable", "intensity", "--conc-units", "equivalents"])
        assert code == 0, out.getvalue()
        summary = json.loads(out.getvalue())
        fit_json = summary.get("result", summary)["json_file"]

        page = _page()
        page._set_input_file(str(csv))
        _set_controls(page)
        page._load_params_file(fit_json)
        assert page.p0_spin.value() == 30.0
        assert page.alpha_spin.value() == 0.2
        assert page.obs_combo.currentIndex() == 2          # intensity only
        assert page.intvalue_combo.currentText() == "height"
        assert page.boot_spin.value() == 0
