# ABOUTME: The GUI T1rho -> T2 conversion converts the fields in use, and asks a carrier of those only.
# ABOUTME: Real JA1 data; one decay matrix serves as both series, as only the bookkeeping is tested.

import os
import sys
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
MATRIX = RUN / "series" / "600_T1" / "peak_intensity_matrix.csv"
PEAKS = RUN.parent / "600_WT" / "T1" / "600_T1_assignment.txt"

pytestmark = pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def fits(tmp_path_factory):
    from lunaNMR.cli import main
    out = tmp_path_factory.mktemp("fits")
    assert main(["dynamixs", "t1rho", "--input", str(MATRIX), "--t1", str(MATRIX),
                 "--peaks", str(PEAKS), "--omega1", "2000", "--carrier", "118",
                 "--theta", "45", "--field-freq", "600", "--out", str(out),
                 "--prefix", "f"]) == 0
    return out


@pytest.fixture
def page(fits, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication(sys.argv)
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    import dynamiXs_gui
    dialog = DynamiXsDialog(parent=None, main_window=types.SimpleNamespace())
    p = dialog.t1t2_page
    p.errors = []
    monkeypatch.setattr(dynamiXs_gui, "show_error", lambda parent, t, m: p.errors.append(m))
    p.output_dir = str(tmp_path)
    p.peak_list_file = str(PEAKS)
    p._toggle_field2()
    for field in ("field1", "field2"):
        p.session_results[field] = {
            "frequency_mhz": 600, "t2": None,
            "t1": {"results_file": str(fits / "f_T1_fit_results.txt")},
            "t1rho": {"results_file": str(fits / "f_T1rho_fit_results.txt")}}
    yield p
    dialog.deleteLater()


def test_a_disabled_field_is_neither_converted_nor_asked_for(page, tmp_path):
    """Field 2 disabled keeps its session results - and must not block field 1."""
    page._toggle_field2()                       # disable field 2 again
    page.field1_carrier_spin.setValue(118)
    page._calculate_t2_from_t1rho()
    assert not page.errors, page.errors
    assert (tmp_path / f"{page._run_prefix()}_field1_T2_from_T1rho.csv").exists()
    assert not (tmp_path / f"{page._run_prefix()}_field2_T2_from_T1rho.csv").exists()


def test_an_enabled_field_without_its_carrier_is_refused(page, tmp_path):
    page.field1_carrier_spin.setValue(118)
    page._calculate_t2_from_t1rho()
    assert page.errors and "Field 2" in page.errors[0]
    assert not list(tmp_path.glob("*T2_from_T1rho*"))


def test_a_measured_t2_is_kept_and_the_derived_one_written_beside_it(page, tmp_path):
    """A fitted T2 is the more direct number; the conversion does not replace it, and
    `project build` does the same with CLI runs."""
    page._toggle_field2()
    measured = {"results_file": str(tmp_path / "measured_T2_fit_results.txt")}
    page.session_results["field1"]["t2"] = dict(measured)
    page.field1_carrier_spin.setValue(118)
    page._calculate_t2_from_t1rho()
    assert not page.errors, page.errors
    assert page.session_results["field1"]["t2"] == measured
    assert (tmp_path / f"{page._run_prefix()}_field1_T2_from_T1rho.csv").exists()
    assert "kept the measured T2" in page.results_text.toPlainText()


def test_a_derived_t2_is_replaced_by_a_new_conversion(page, tmp_path):
    page._toggle_field2()
    page.session_results["field1"]["t2"] = {"results_file": "old", "derived_from": "t1rho"}
    page.field1_carrier_spin.setValue(118)
    page._calculate_t2_from_t1rho()
    assert page.session_results["field1"]["t2"]["results_file"].endswith(
        f"{page._run_prefix()}_field1_T2_from_T1rho_fit_results.txt")
