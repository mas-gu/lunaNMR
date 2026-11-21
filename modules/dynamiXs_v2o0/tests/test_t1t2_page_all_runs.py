# ABOUTME: The T1/T2 page's "All" selections fit T1 then T2 for each field they name.
# ABOUTME: Real page and real fitter on small synthetic decay matrices written per test.

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

DELAYS_MS = [0, 50, 100, 200, 400, 800, 1200]


def _matrix(path, time_constant_ms):
    """A series intensity matrix of four residues decaying with `time_constant_ms`."""
    rows = []
    for number, name in enumerate(("3LysH", "8TyrH", "10ValH", "11LeuH"), start=1):
        amplitude = 1.0e6 * number
        row = {"Peak_Number": number, "Assignment": name,
               "Reference_X": 8.0 + 0.1 * number, "Reference_Y": 118.0 + number}
        for delay in DELAYS_MS:
            row[str(delay)] = amplitude * np.exp(-delay / (time_constant_ms + 10 * number))
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def page(app, tmp_path):
    from dynamiXs_gui import DynamiXsMainWindow
    window = DynamiXsMainWindow()
    page = window.t1t2_page
    page.output_dir = str(tmp_path / "out")
    os.makedirs(page.output_dir)
    page._toggle_field2()
    for field in ("field1", "field2"):
        setattr(page, f"{field}_t1_file", _matrix(tmp_path / f"{field}_T1.csv", 700.0))
        setattr(page, f"{field}_t2_file", _matrix(tmp_path / f"{field}_T2.csv", 90.0))
    yield page
    _run_to_the_end(app, page)
    window.deleteLater()


def _run_to_the_end(app, page, timeout_s=300):
    """Process events until the page has no fit running and none waiting."""
    deadline = time.monotonic() + timeout_s
    while page.worker is not None or getattr(page, "_fits_waiting", []):
        assert time.monotonic() < deadline, "the page's fits did not finish"
        app.processEvents()
        time.sleep(0.02)
    app.processEvents()


def _fitted(page):
    """The (field, experiment) pairs the page's session holds a fit for."""
    return {(field, exp) for field, data in page.session_results.items() if data
            for exp in ("t1", "t2", "t1rho") if data.get(exp)}


@pytest.mark.parametrize("selection, expected", [
    ("Field 1 - All", {("field1", "t1"), ("field1", "t2")}),
    ("Field 2 - All", {("field2", "t1"), ("field2", "t2")}),
    ("Both Fields - All", {("field1", "t1"), ("field1", "t2"),
                           ("field2", "t1"), ("field2", "t2")}),
])
def test_all_fits_t1_then_t2_for_its_fields(app, page, selection, expected):
    page.analysis_combo.setCurrentText(selection)
    assert page.analysis_combo.currentText() == selection

    page._start_analysis()
    _run_to_the_end(app, page)

    assert _fitted(page) == expected
    for field, exp in expected:
        entry = page.session_results[field][exp]
        assert entry["n_fitted"] == 4
        assert os.path.exists(entry["results_file"])
    log = page.results_text.toPlainText()
    for field in {field for field, _ in expected}:
        assert log.index(f"{field}_t1 stored") < log.index(f"{field}_t2 stored")


def test_all_fits_only_the_experiments_that_have_data(app, page):
    page.field1_t2_file = None
    page.analysis_combo.setCurrentText("Field 1 - All")

    page._start_analysis()
    _run_to_the_end(app, page)

    assert _fitted(page) == {("field1", "t1")}
