# ABOUTME: The Kd page fits only concentrations whose units someone stated, as `kd` does:
# ABOUTME: unstated by default, refused on the page, reset per input. Real titrations, real worker.

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
ERDJ6 = TIDY / "ERDJ6_HSPA8_tidy.csv"
EQUIVALENTS = [0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 3.0]

pytestmark = pytest.mark.skipif(not (DNAJA1.exists() and ERDJ6.exists()),
                                reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _copy(src, folder):
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / src.name
    shutil.copy(src, path)
    return path


def _page():
    from kd_titration_page import KdTitrationPage
    return KdTitrationPage(types.SimpleNamespace(show_main_menu=lambda: None,
                                                 current_dir=""))


class TestTheUnitsStartUnstated:

    def test_the_selector_says_to_state_them(self, app):
        page = _page()
        assert "state the units" in page.conc_units_combo.currentText()
        assert page._conc_units() is None

    def test_the_point_labels_still_fill_the_concentrations(self, app, tmp_path):
        page = _page()
        page._set_input_file(str(_copy(DNAJA1, tmp_path)))
        assert [s.value() for s in page.conc_spins] == EQUIVALENTS
        assert page._conc_units() is None


class TestARunWithUnstatedUnitsIsRefusedOnThePage:
    """Labels 0-3 eq fitted as micromolar gave a 'reliable' intensity Kd of 0.233
    instead of 11.67 - wrong by the factor [P]0, reported in green."""

    def test_no_fit_starts(self, app, tmp_path):
        page = _page()
        page._set_input_file(str(_copy(DNAJA1, tmp_path)))
        page.obs_combo.setCurrentIndex(2)                 # intensity only
        page._start_analysis()
        assert page.worker is None
        assert not (tmp_path / "kd_analysis").exists()

    def test_the_page_says_why(self, app, tmp_path):
        page = _page()
        page._set_input_file(str(_copy(DNAJA1, tmp_path)))
        page._start_analysis()
        assert page.result_label.isVisibleTo(page)
        assert "units" in page.result_label.text()
        assert "units" in page.log_text.toPlainText().splitlines()[-1]

    def test_stating_them_lets_it_run(self, app, tmp_path):
        page = _page()
        page._set_input_file(str(_copy(DNAJA1, tmp_path)))
        page.obs_combo.setCurrentIndex(2)
        page.conc_units_combo.setCurrentIndex(1)          # equivalents
        page._start_analysis()
        assert page.worker.wait(300_000), "the fit did not finish"
        app.processEvents()
        fit = json.loads(Path(page.last_json_file).read_text())
        assert fit["global"]["intensity"]["Kd"] == pytest.approx(11.67, rel=0.01)


class TestUnitsComeOnlyFromThisInput:

    def test_a_new_input_resets_them_with_its_concentrations(self, app, tmp_path):
        page = _page()
        page._set_input_file(str(_copy(DNAJA1, tmp_path / "a")))
        page.conc_units_combo.setCurrentIndex(1)
        page._set_input_file(str(_copy(ERDJ6, tmp_path / "b")))
        assert page._conc_units() is None

    def test_a_params_file_for_this_input_states_them(self, app, tmp_path):
        csv = _copy(DNAJA1, tmp_path)
        data = tmp_path / "kd_analysis" / "data"
        data.mkdir(parents=True)
        (data / "x_kd_params.json").write_text(json.dumps(
            {"input": csv.name, "points": EQUIVALENTS, "conc_units": "equivalents"}))
        page = _page()
        page._set_input_file(str(csv))
        assert page._conc_units() == "equivalents"

    def test_a_params_file_silent_on_units_leaves_them_unstated(self, app, tmp_path):
        csv = _copy(DNAJA1, tmp_path)
        (tmp_path / "kd_params.json").write_text(json.dumps({"protein_conc": 30.0}))
        page = _page()
        page._set_input_file(str(csv))
        assert page._conc_units() is None

    def test_the_log_names_the_protein_concentration_it_kept(self, app, tmp_path):
        """[P]0 is form state like alpha and the observable; the log says it carried."""
        page = _page()
        page._set_input_file(str(_copy(DNAJA1, tmp_path / "a")))
        page.p0_spin.setValue(20.0)
        page._set_input_file(str(_copy(ERDJ6, tmp_path / "b")))
        assert page.p0_spin.value() == 20.0
        log = page.log_text.toPlainText().splitlines()
        assert any("[P]0" in line and "20" in line and "kept" in line for line in log)

    def test_a_first_input_names_no_kept_value(self, app, tmp_path):
        page = _page()
        page._set_input_file(str(_copy(DNAJA1, tmp_path)))
        assert "kept" not in page.log_text.toPlainText()


class TestAParamsFileRecordsNoUnitsNobodyChose:
    """A file recording the point labels as 'absolute' became a statement the CLI then
    trusted, defeating its refusal to fit labels with no stated unit."""

    def test_unstated_units_write_neither_units_nor_concentrations(self, app, tmp_path):
        from kd_params import dump_params_json
        page = _page()
        page._set_input_file(str(_copy(DNAJA1, tmp_path)))
        path = tmp_path / "kd_params.json"
        dump_params_json(path, page._gather_params_dict())
        written = json.loads(path.read_text())
        assert "conc_units" not in written
        assert written["concentrations"] is None

    def test_the_cli_still_refuses_after_such_a_file(self, app, tmp_path):
        from kd_params import dump_params_json
        from lunaNMR.cli import main
        csv = _copy(DNAJA1, tmp_path)
        page = _page()
        page._set_input_file(str(csv))
        dump_params_json(tmp_path / "kd_params.json", page._gather_params_dict())

        out = io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
            code = main(["kd", "--format", "json", "--p0", "50", "--input", str(csv),
                         "--out", str(tmp_path / "kd_analysis")])
        assert code != 0
        assert "No concentrations" in out.getvalue()

    def test_stated_units_are_written(self, app, tmp_path):
        from kd_params import dump_params_json
        page = _page()
        page._set_input_file(str(_copy(DNAJA1, tmp_path)))
        page.conc_units_combo.setCurrentIndex(1)
        path = tmp_path / "kd_params.json"
        dump_params_json(path, page._gather_params_dict())
        written = json.loads(path.read_text())
        assert written["conc_units"] == "equivalents"
        assert written["concentrations"] == EQUIVALENTS


class TestTheSessionKeepsThemUnstated:

    def test_unstated_round_trips(self, app, tmp_path):
        page = _page()
        page._set_input_file(str(_copy(DNAJA1, tmp_path)))
        state = json.loads(json.dumps(page.get_session_state()))

        reopened = _page()
        reopened.conc_units_combo.setCurrentIndex(1)
        reopened.restore_session_state(state)
        assert reopened._conc_units() is None

    def test_a_cli_fit_bundled_headless_reopens_in_its_units(self, app, tmp_path):
        """The fit JSON records the L it used, in absolute units; a page state built
        from it must say so rather than leave them unstated."""
        from lunaNMR.utils.project_builder import kd_page_state
        page = _page()
        page._set_input_file(str(_copy(DNAJA1, tmp_path)))
        page.obs_combo.setCurrentIndex(2)
        page.conc_units_combo.setCurrentIndex(1)
        page._start_analysis()
        assert page.worker.wait(300_000), "the fit did not finish"
        app.processEvents()
        fit = json.loads(Path(page.last_json_file).read_text())

        reopened = _page()
        reopened.restore_session_state(kd_page_state(fit, page.last_json_file))
        assert reopened._conc_units() == "absolute"
        assert [s.value() for s in reopened.conc_spins] == [50 * p for p in EQUIVALENTS]
