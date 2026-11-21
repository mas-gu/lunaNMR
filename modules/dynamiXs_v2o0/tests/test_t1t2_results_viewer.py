# ABOUTME: The T1/T2 results viewer opened on a session exports its plot, shows a results
# ABOUTME: file loaded by hand, and labels times in the units each fit was made in.

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))

RESIDUES = ("3LysH", "8TyrH", "10ValH", "11LeuH", "12GlyH")


def _results_file(path, experiment, base):
    lines = [f"Residue\tA\t{experiment}\tC\tA_err\t{experiment}_err\tC_err\tSuccess\tWindowRatio"]
    for n, name in enumerate(RESIDUES):
        value = base * (1 + 0.02 * n)
        lines.append(f"{name}\t1.000000e+06\t{value:.6e}\t0.000000e+00\t1.000000e+03\t"
                     f"{value * 0.02:.6e}\t0.000000e+00\tYes\t3.000")
    path.write_text("\n".join(lines) + "\n")
    return str(path)


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def boxes(monkeypatch):
    """The viewer's message boxes, recorded as (kind, title) instead of shown."""
    import visualization.t1t2_results_viewer as module
    shown = []
    for kind in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(module, kind,
                            lambda parent, title, text, kind=kind: shown.append((kind, title)))
    return shown


@pytest.fixture
def session_viewer(app, tmp_path, boxes):
    """A viewer opened as the T1/T2 page opens it: a T1 in seconds and a T2 in ms."""
    from visualization.t1t2_results_viewer import T1T2ResultsViewer
    viewer = T1T2ResultsViewer(session_results={
        "field1_t1": {"results_file": _results_file(tmp_path / "T1_fit_results.txt", "T1", 0.7),
                      "frequency_mhz": 600.0, "experiment_type": "T1", "time_units": "s"},
        "field1_t2": {"results_file": _results_file(tmp_path / "T2_fit_results.txt", "T2", 90.0),
                      "frequency_mhz": 600.0, "experiment_type": "T2", "time_units": "ms"},
    })
    yield viewer
    viewer.close()


def test_export_plot_writes_the_session_plot(session_viewer, tmp_path, boxes, monkeypatch):
    import visualization.t1t2_results_viewer as module
    target = tmp_path / "exported.pdf"
    monkeypatch.setattr(module.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(target), "PDF files (*.pdf)")))

    session_viewer._export_plot()

    assert boxes == [("show_info", "Export Complete")]
    assert target.exists() and target.stat().st_size > 0


def test_a_results_file_loaded_by_hand_replaces_the_session_plot(session_viewer, tmp_path, boxes):
    assert len(session_viewer.figure.axes) == 2
    other = _results_file(tmp_path / "other_fit_results.txt", "T2", 40.0)

    session_viewer._load_results_file(other)

    (ax,) = session_viewer.figure.axes
    assert ax.get_title() == "T2 Relaxation Times vs Residue Sequence"
    assert session_viewer.file_label.text() == "other_fit_results.txt"
    assert "Total residues: 5" in session_viewer.stats_label.text()
    assert "Mean: 41.60" in session_viewer.stats_label.text()


def test_each_experiment_is_labelled_in_its_own_units(session_viewer):
    labels = [ax.get_ylabel() for ax in session_viewer.figure.axes]
    assert labels == ["T1 (s)", "T2 (ms)"]
    stats = session_viewer.stats_label.text().splitlines()
    assert stats[0].startswith("field1_t1:") and stats[0].endswith(" s")
    assert stats[1].startswith("field1_t2:") and stats[1].endswith(" ms")
    assert session_viewer.t1_ymax_spin.suffix() == " s"
    assert session_viewer.t2_ymax_spin.suffix() == " ms"


def test_a_session_that_records_no_units_is_labelled_in_ms(app, tmp_path, boxes):
    from visualization.t1t2_results_viewer import T1T2ResultsViewer
    viewer = T1T2ResultsViewer(session_results={
        "field1_t1": {"results_file": _results_file(tmp_path / "T1_fit_results.txt", "T1", 700.0),
                      "frequency_mhz": 600.0, "experiment_type": "T1"}})
    try:
        assert [ax.get_ylabel() for ax in viewer.figure.axes] == ["T1 (ms)"]
        assert viewer.t1_ymax_spin.suffix() == " ms"
    finally:
        viewer.close()


def test_the_page_opens_the_viewer_in_the_units_of_its_fit(app, tmp_path, boxes):
    """A series in seconds, fitted on the T1/T2 page with the selector on s."""
    import numpy as np
    import pandas as pd
    from dynamiXs_gui import DynamiXsMainWindow
    delays_s = [0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2]
    rows = [{"Assignment": name,
             **{f"{d:g}": 1.0e6 * np.exp(-d / (0.7 + 0.01 * n)) for d in delays_s}}
            for n, name in enumerate(RESIDUES)]
    matrix = tmp_path / "T1_seconds.csv"
    pd.DataFrame(rows).to_csv(matrix, index=False)
    window = DynamiXsMainWindow()
    try:
        page = window.t1t2_page
        page.output_dir = str(tmp_path / "out")
        os.makedirs(page.output_dir)
        page.field1_t1_file = str(matrix)
        page.units_combo.setCurrentText("s")
        page.analysis_combo.setCurrentText("Field 1 - T1")
        page._start_analysis()
        page.worker.wait()
        app.processEvents()
        assert page.session_results["field1"]["t1"], page.results_text.toPlainText()

        page._display_results()

        (ax,) = page.results_viewer.figure.axes
        assert ax.get_ylabel() == "T1 (s)"
        page.results_viewer.close()
    finally:
        window.deleteLater()
