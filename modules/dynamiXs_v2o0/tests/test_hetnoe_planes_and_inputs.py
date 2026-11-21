# ABOUTME: One rule names a hetNOE plane from a column or file name, for the hetNOE page and
# ABOUTME: the model-free Auto-Load alike; matrix dummy rows are no residues; no unused field.

import os
import sys
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def window(app):
    from dynamiXs_gui import DynamiXsMainWindow
    window = DynamiXsMainWindow()
    yield window
    window.deleteLater()


@pytest.fixture
def hetnoe_page(window):
    from hetnoe_page import HetNOEPage
    page = HetNOEPage(window)
    yield page
    page.deleteLater()


@pytest.mark.parametrize("name, plane", [
    ("hetnoe_600_saturated", "sat"), ("f600_noe_sat.csv", "sat"), ("NOE_SAT", "sat"),
    ("hetnoe_600_unsaturated", "unsat"), ("f600_noe_unsat.csv", "unsat"),
    ("NOE_UNSAT.csv", "unsat"), ("f600_NOSAT.csv", "unsat"),
    ("600_T1", None), ("reference", None),
])
def test_a_name_states_its_plane(name, plane):
    from dynamiXs_integrated.data_converters import hetnoe_plane
    assert hetnoe_plane(name) == plane


def _matrix(path, columns):
    rows = [{"Peak_Number": n, "Assignment": name, "Reference_X": 8.0, "Reference_Y": 120.0,
             **{column: 1.0e6 * (n + i) for i, column in enumerate(columns)}}
            for n, name in enumerate(("3LysH", "8TyrH", "dummy_1", "Dummy_2"), start=1)]
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


@pytest.mark.parametrize("columns, sat, unsat", [
    (["hetnoe_600_saturated", "hetnoe_600_unsaturated"],
     "hetnoe_600_saturated", "hetnoe_600_unsaturated"),
    (["hetnoe_600_nosat", "hetnoe_600_sat"], "hetnoe_600_sat", "hetnoe_600_nosat"),
])
def test_the_hetnoe_page_preselects_its_planes_by_that_rule(hetnoe_page, tmp_path,
                                                             columns, sat, unsat):
    hetnoe_page._set_matrix(_matrix(tmp_path / "matrix.csv", columns))
    assert hetnoe_page.sat_combo.currentText() == sat
    assert hetnoe_page.unsat_combo.currentText() == unsat


def test_auto_load_puts_each_plane_file_in_its_own_slot(window, tmp_path, monkeypatch):
    import dynamiXs_gui
    for name in ("f600_T1.csv", "f600_T2.csv", "f600_noe_sat.csv", "f600_noe_unsat.csv"):
        (tmp_path / name).write_text("3LysH,1.0\n")
    monkeypatch.setattr(dynamiXs_gui, "open_directory_dialog", lambda *a, **k: str(tmp_path))
    page = window.integrated_page

    page._auto_load_folder("field1")

    assert os.path.basename(page.field1_noe_sat_file) == "f600_noe_sat.csv"
    assert os.path.basename(page.field1_noe_unsat_file) == "f600_noe_unsat.csv"
    assert page.field1_noe_unsat_display.text() == "f600_noe_unsat.csv"
    assert os.path.basename(page.field1_t1_file) == "f600_T1.csv"
    assert os.path.basename(page.field1_t2_file) == "f600_T2.csv"


def test_a_matrix_dummy_row_is_not_a_residue(tmp_path):
    from dynamiXs_integrated.data_converters import intensities_from_matrix
    matrix = _matrix(tmp_path / "matrix.csv", ["hetnoe_600_saturated", "hetnoe_600_unsaturated"])
    assert intensities_from_matrix(matrix, "hetnoe_600_saturated") == {"3LysH": 1.0e6,
                                                                       "8TyrH": 2.0e6}


def test_the_hetnoe_page_computes_no_hetnoe_for_dummy_rows(hetnoe_page, tmp_path):
    hetnoe_page._set_matrix(_matrix(tmp_path / "matrix.csv",
                                    ["hetnoe_600_saturated", "hetnoe_600_unsaturated"]))
    hetnoe_page.output_dir = str(tmp_path / "out")
    hetnoe_page._run()
    results = pd.read_csv(tmp_path / "out" / "field1_hetnoe.csv")
    assert sorted(results["Residue"]) == ["3LysH", "8TyrH"]


def test_the_hetnoe_page_asks_no_field_frequency(hetnoe_page):
    """hetNOE is a ratio of two intensities: no field enters it."""
    from PySide6.QtWidgets import QLabel
    labels = [label.text() for label in hetnoe_page.findChildren(QLabel)]
    assert not any("frequency" in text.lower() for text in labels)
    assert "field_freq" not in hetnoe_page.get_session_state()
    # a state that carries a field frequency restores, the value unread
    hetnoe_page.restore_session_state({"prefix": "f600", "field_freq": 700.0})
    assert hetnoe_page.prefix_edit.text() == "f600"
