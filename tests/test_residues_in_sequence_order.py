# ABOUTME: Residues are ordered by one rule everywhere: by residue number, a name without a
# ABOUTME: number last. lunaNMR and the DynamiXs/Kd modules each carry the rule; they must agree.

import json
import os
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "modules" / "dynamiXs_v2o0"))

RUN = ROOT.parent / "data_test" / "data_JA1_relaxation" / "run_20260825_1405"

# What a peak list, a matrix read without a dtype and a detection can each name a row.
NAMES = ["10ValH", "TrpSide", "A17", 45.0, "3LysH", "", "argSide", "17.0", "dummy1", "K14",
         "ArgSide", 8, "dummy_3", "3.LYS"]
IN_ORDER = ["3.LYS", "3LysH", 8, "10ValH", "K14", "17.0", "A17", 45.0,
            "", "ArgSide", "TrpSide", "argSide", "dummy1", "dummy_3"]


def _number(name):
    found = re.search(r"\d+", str(name))
    return int(found.group()) if found else float("inf")


class TestTheTwoCopiesAgree:
    def test_lunanmr_orders_by_number_and_puts_unnumbered_names_last(self):
        from lunaNMR.utils.residue_order import residue_sort_key
        assert sorted(NAMES, key=residue_sort_key) == IN_ORDER

    def test_the_dynamixs_copy_gives_the_same_order(self):
        from lunaNMR.utils.residue_order import residue_sort_key
        from residue_order import residue_sort_key as dynamixs_sort_key
        assert sorted(NAMES, key=dynamixs_sort_key) == sorted(NAMES, key=residue_sort_key)
        assert [dynamixs_sort_key(name) for name in NAMES] == [
            residue_sort_key(name) for name in NAMES]


def test_the_r2_table_from_t1rho_is_in_sequence_order():
    from lunaNMR.utils.t1rho_calculator import r2_table_from_fits
    names = ["TrpSide", "10ValH", "3LysH", "ArgSide", "21GlyH"]
    t1 = {name: {"value": 650.0 + i, "error": 1.0} for i, name in enumerate(names)}
    t1rho = {name: {"value": 80.0 + i, "error": 0.5} for i, name in enumerate(names)}
    peaks = pd.DataFrame([{"residue": name, "N_ppm": 118.0 + i} for i, name in enumerate(names)])

    table = r2_table_from_fits(t1, t1rho, peaks, omega1_hz=2000.0, carrier_ppm=118.0,
                               theta_deg=90.0, spec_freq_mhz=600.0)

    assert list(table["residue"]) == ["3LysH", "10ValH", "21GlyH", "ArgSide", "TrpSide"]
    assert list(table["R1rho"]) == pytest.approx([1000.0 / 82.0, 1000.0 / 81.0, 1000.0 / 84.0,
                                                  1000.0 / 83.0, 1000.0 / 80.0])


def test_diagnose_names_duplicate_assignments_in_sequence_order(tmp_path):
    from lunaNMR.validation.spectra_check import load_peak_list
    rows = [("10ValH", 8.20, 120.1), ("3LysH", 8.18, 126.5), ("10ValH", 8.25, 121.0),
            ("21GlyH", 8.40, 109.3), ("3LysH", 8.10, 125.9)]
    peak_list = tmp_path / "peaks.txt"
    peak_list.write_text("Assignment, Position_X, Position_Y\n"
                         + "".join(f"{name}, {x:.4f}, {y:.4f}\n" for name, x, y in rows))

    findings = load_peak_list(peak_list)["findings"]

    assert [f.message for f in findings] == [
        "peaks.txt: duplicate assignments ['3LysH', '10ValH']"]


class TestResultsBrowser:
    @pytest.fixture(scope="class")
    def app(self):
        from PySide6.QtWidgets import QApplication
        return QApplication.instance() or QApplication(sys.argv)

    @pytest.fixture
    def browser(self, app, dnaja1_titration_series):
        from lunaNMR.gui.dialogs.results_browser_dialog import ResultsBrowserDialog
        from lunaNMR.utils.project_manager import batch_results_from_dict
        batch = batch_results_from_dict(
            json.loads((dnaja1_titration_series / "batch_results.json").read_text()))
        dialog = ResultsBrowserDialog(batch_results=batch)
        yield dialog
        dialog.deleteLater()

    @pytest.mark.parametrize("combo", ["peak_combo", "evolution_peak_combo"])
    def test_the_peak_selectors_list_residues_in_sequence_order(self, browser, combo):
        box = getattr(browser, combo)
        listed = [box.itemText(i) for i in range(box.count())]
        dummies = [name for name in listed if name.lower().startswith('dummy')]
        residues = listed[:len(listed) - len(dummies)]
        numbers = [_number(name) for name in residues]

        assert len(residues) > 20
        assert numbers == sorted(numbers)
        assert listed.index("3LysH") < listed.index("17AlaH")
        assert listed[len(residues):] == dummies


@pytest.mark.skipif(not RUN.exists(), reason="data_test not present")
def test_model_free_names_the_residues_dropped_at_the_merge_in_sequence_order(tmp_path, capsys):
    """10ValH is in no hetNOE table of the real 600 MHz inputs; 3LysH is taken out of one."""
    from lunaNMR.cli import main
    noe = RUN / "results" / "hetnoe"
    saturated = tmp_path / "noe_sat.csv"
    saturated.write_text("".join(
        line for line in (noe / "f600_noe_sat.csv").read_text().splitlines(keepends=True)
        if not line.startswith("3LysH,")))

    assert main(["dynamixs", "modelfree",
                 "--f1-t1", str(RUN / "series" / "600_T1" / "peak_intensity_matrix.csv"),
                 "--f1-t2", str(RUN / "series" / "600_T2" / "peak_intensity_matrix.csv"),
                 "--f1-noe-sat", str(saturated),
                 "--f1-noe-unsat", str(noe / "f600_noe_unsat.csv"),
                 "--field1-freq", "600", "--out", str(tmp_path / "mf"), "--prefix", "mf",
                 "--n-monte-carlo", "2", "--format", "json"]) == 0

    excluded = json.loads(capsys.readouterr().out)["excluded_residues"]["field1"]
    assert excluded[:2] == ["3LysH", "10ValH"]
    assert [_number(name) for name in excluded] == sorted(_number(name) for name in excluded)
