# ABOUTME: The spectral density page's session records the results table its run wrote,
# ABOUTME: single- and dual-field. Real page and analysis on a small synthetic rate table.

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))


def _rates(path, r2):
    """A spectral-density input table of six residues."""
    rng = np.random.default_rng(1)
    n = 6
    pd.DataFrame({
        "Residue": [f"{i}AlaH" for i in range(3, 3 + n)],
        "R1": 1.4 + 0.02 * rng.standard_normal(n), "R1err": 0.03,
        "R2": r2 + 0.1 * rng.standard_normal(n), "R2err": 0.2,
        "hetNOE": 0.80 + 0.01 * rng.standard_normal(n), "hetNOEerr": 0.02,
    }).to_csv(path, index=False)
    return str(path)


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def page(app, tmp_path):
    from dynamiXs_gui import DynamiXsMainWindow
    window = DynamiXsMainWindow()
    window.current_dir = str(tmp_path / "out")
    os.makedirs(window.current_dir)
    yield window.spectral_page
    if window.spectral_page.worker is not None:
        window.spectral_page.worker.wait()
    window.deleteLater()


def _run(app, page, timeout_s=300):
    page._start_analysis()
    deadline = time.monotonic() + timeout_s
    while not page.session_results["analysis_complete"]:
        assert time.monotonic() < deadline, page.results_text.toPlainText()
        assert "Error during analysis" not in page.results_text.toPlainText(), \
            page.results_text.toPlainText()
        app.processEvents()
        time.sleep(0.02)
    return page.session_results


def test_a_single_field_run_records_its_results_table(app, page, tmp_path):
    page.input_file1 = _rates(tmp_path / "f1.csv", 11.0)
    page.prefix_edit.setText("sd")

    session = _run(app, page)

    assert session["results_file"] == os.path.join(page.output_dir, "sd_results.csv")
    assert len(pd.read_csv(session["results_file"])) == 6


def test_a_dual_field_run_records_its_detailed_table(app, page, tmp_path):
    page.dual_field_radio.setChecked(True)
    page.input_file1 = _rates(tmp_path / "f1.csv", 11.0)
    page.input_file2 = _rates(tmp_path / "f2.csv", 12.5)
    page.prefix_edit.setText("sd")

    session = _run(app, page)

    assert session["results_file"] == os.path.join(page.output_dir, "sd_detailed.csv")
    assert os.path.exists(session["results_file"])
