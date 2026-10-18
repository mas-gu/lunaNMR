# ABOUTME: The 'Auto-add dummy peaks' switch reaches a series, from the dialog and from `series --params`.
# ABOUTME: Real DNAJA1 series: with it off the outputs carry no dummy_ rows, with it on (the default) they do.

"""
The series builds its own integrator, whose auto_add_dummy_peaks was always on: a user
who unticked the Expert box still got dummy_ peaks from the first spectrum's detection,
fitted at every point and written to the matrices and the tidy CSV. The CLI had no
switch at all.
"""

import contextlib
import io
import json
import sys
import types
import warnings
from pathlib import Path

import pandas as pd
import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")

SPECTRA = ("2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o125.ft")


def _spectra_folder(tmp_path_factory, name):
    folder = tmp_path_factory.mktemp(name)
    for spectrum in SPECTRA:
        (folder / spectrum).symlink_to(DNAJA1 / spectrum)
    return folder


def _dummy_rows(tidy_csv):
    tidy = pd.read_csv(tidy_csv)
    return int(tidy["assignment"].astype(str).str.startswith("dummy_").sum())


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def _dialog_run(app, folder, auto_add_dummy_peaks):
    """A cascade series with no Fit All Peaks first, so spectrum 1 is detected."""
    from numba.core.errors import NumbaWarning
    from PySide6.QtWidgets import QDoubleSpinBox
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    from lunaNMR.utils.file_manager import NMRFileManager
    contour = QDoubleSpinBox()
    contour.setDecimals(5)
    contour.setValue(0.05)
    with contextlib.redirect_stdout(io.StringIO()):
        peak_list = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
    main_window = types.SimpleNamespace(
        current_nmr_folder=str(folder), search_window_x=0.07, search_window_y=0.07,
        use_parallel_processing=False, saved_series={}, contour_min_spin=contour,
        auto_add_dummy_peaks=auto_add_dummy_peaks,
        integrator=types.SimpleNamespace(fitted_peaks=None, peak_list=peak_list))
    dialog = SeriesIntegrationDialog(main_window=main_window)
    dialog._prompt_series_name = lambda: "x"   # the user's answer to the name prompt
    dialog.cascade_radio.setChecked(True)
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", NumbaWarning)
        dialog._start_processing()
        while dialog.thread.isRunning():
            app.processEvents()
            dialog.thread.wait(50)
        app.processEvents()
    dialog.outcome_box.close()
    dialog.deleteLater()
    tidy, = folder.glob("series_results_*/series_analysis_tidy.csv")
    return tidy


@pytest.fixture(scope="module")
def dialog_tidy(app, tmp_path_factory):
    return {on: _dialog_run(app, _spectra_folder(tmp_path_factory, f"dialog_{on}"), on)
            for on in (False, True)}


def _cli_run(tmp_path_factory, auto_add_dummy_peaks):
    from numba.core.errors import NumbaWarning
    from lunaNMR.cli import main
    folder = _spectra_folder(tmp_path_factory, f"cli_{auto_add_dummy_peaks}")
    args = ["series", "--spectra", str(folder), "--peaks", str(DNAJA1 / "DNAJA1.txt"),
            "--out", str(folder / "out"), "--format", "json"]
    if auto_add_dummy_peaks is not None:
        params = folder / "params.json"
        params.write_text(json.dumps(
            {"gui_params": {"auto_add_dummy_peaks": auto_add_dummy_peaks}}))
        args += ["--params", str(params)]
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", NumbaWarning)
        assert main(args) == 0
    return folder / "out" / "series_analysis_tidy.csv"


@pytest.fixture(scope="module")
def cli_tidy(tmp_path_factory):
    return {on: _cli_run(tmp_path_factory, on) for on in (False, None)}


def test_the_dialog_with_it_off_writes_no_dummies(dialog_tidy):
    assert _dummy_rows(dialog_tidy[False]) == 0


def test_the_dialog_with_it_on_writes_them(dialog_tidy):
    assert _dummy_rows(dialog_tidy[True]) > 0


def test_the_cli_with_it_off_writes_no_dummies(cli_tidy):
    assert _dummy_rows(cli_tidy[False]) == 0


def test_the_cli_adds_them_by_default(cli_tidy):
    assert _dummy_rows(cli_tidy[None]) > 0


def test_the_cli_default_is_on():
    from lunaNMR.cli import _default_series_params
    assert _default_series_params()['gui_params']['auto_add_dummy_peaks'] is True
