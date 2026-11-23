# ABOUTME: The 2D Voigt view states the R-squared and the quality label the fit it draws carries.
# ABOUTME: Real DNAJA1 fits of both kinds: as the fitter returns them, and as a series' per-peak rows.

"""
Two records describe one fit. The fitter's result, which the main window holds after
Fit Spectrum, names its label `fitting_quality`; a series' per-peak row names it
`quality`. Both carry `r_squared`. The view must state the same R-squared and label for
either, where it stated "Quality: Unknown, R2 = 0.0000" for every series row.
"""

import copy
import json
import math
import os
import sys
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRUM = "2D_DNAJA1_HSPA8_0o5.ft"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture(scope="module")
def stored(dnaja1_titration_series):
    return json.loads((dnaja1_titration_series / "batch_results.json").read_text())


@pytest.fixture(scope="module")
def written(dnaja1_titration_series):
    """What the series wrote to its tidy CSV for a peak of the second spectrum whose fit
    is neither the best nor a failure."""
    folder = dnaja1_titration_series
    columns = json.loads((folder / "series_metadata.json").read_text())['columns']
    label = next(str(c['column']) for c in columns if c['spectrum'] == SPECTRUM)
    tidy = pd.read_csv(folder / "series_analysis_tidy.csv", dtype={'spectrum_name': str})
    rows = tidy[(tidy['spectrum_name'] == label) & tidy['quality'].isin(['Good', 'Fair'])]
    return rows.iloc[0]


def _rows(stored, kind):
    return stored['results'][SPECTRUM][kind]


@pytest.fixture(scope="module")
def integrator(stored):
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    it = EnhancedVoigtIntegrator()
    it.load_nmr_file(os.path.join(stored['metadata']['data_folder'], SPECTRUM))
    return it


def _fit(stored, kind, assignment, integrator):
    """The stored record of `assignment` of one kind, with its region and surfaces."""
    from lunaNMR.utils.project_manager import rebuild_fit_arrays
    fit = copy.deepcopy(next(row for row in _rows(stored, kind)
                             if row.get('assignment') == assignment))
    assert rebuild_fit_arrays(fit, integrator)
    return fit


def _shown(plotter):
    """Every text the plotter's figure draws."""
    texts = [text.get_text() for ax in plotter.figure.axes for text in ax.texts]
    texts += [ax.get_title() for ax in plotter.figure.axes]
    texts += [text.get_text() for text in plotter.figure.texts]
    return "\n".join(texts)


def test_the_two_records_name_the_label_differently_and_the_r_squared_alike(stored, written):
    fitted = next(r for r in _rows(stored, 'fitted_results') if r['assignment'] == written['assignment'])
    row = next(r for r in _rows(stored, 'integration_results') if r['assignment'] == written['assignment'])
    assert fitted['fitting_quality'] == row['quality'] == written['quality']
    assert 'fitting_quality' not in row
    assert fitted['r_squared'] == row['r_squared']
    assert math.isclose(row['r_squared'], written['r_squared'], abs_tol=1e-4)


@pytest.mark.parametrize("kind", ["fitted_results", "integration_results"])
def test_the_view_states_the_r_squared_and_label_of_the_fit(app, stored, written, integrator, kind):
    from lunaNMR.gui.components.voigt_analysis_plotter import VoigtAnalysisPlotter
    fit = _fit(stored, kind, written['assignment'], integrator)
    plotter = VoigtAnalysisPlotter()
    try:
        plotter.plot_voigt_analysis(fit)
        shown = _shown(plotter)
        assert f"Quality: {written['quality']}" in shown
        assert f"R² = {fit['r_squared']:.4f}" in shown
        assert f"(R²={fit['r_squared']:.3f})" in shown
    finally:
        plotter.deleteLater()


def test_the_overlay_states_them_for_a_series_fit(app, dnaja1_titration_series, written,
                                                  open_overlay_viewer):
    viewer = open_overlay_viewer(dnaja1_titration_series)
    try:
        spectrum = next(i for i, spec in enumerate(viewer.spectra) if spec['name'] == SPECTRUM)
        peaks = viewer._get_fitted_peaks(spectrum)
        viewer.selected_spectrum_index = spectrum
        viewer.selected_peak_index = next(i for i, peak in enumerate(peaks)
                                          if peak['assignment'] == written['assignment'])
        viewer._plot_voigt_analysis()
        fit = peaks[viewer.selected_peak_index]
        shown = _shown(viewer.voigt_2d_plotter)
        assert f"Quality: {written['quality']}" in shown
        assert f"R² = {fit['r_squared']:.4f}" in shown
        parameters = viewer.peak_params_text.toPlainText()
        assert f"Overall Quality: {written['quality']}" in parameters
        assert f"R² = {fit['r_squared']:.4f}" in parameters
    finally:
        viewer.deleteLater()


def test_the_main_window_states_them_for_its_own_fit(app, tmp_path):
    """Fit Spectrum on ten peaks, then the 2D Voigt tab of the first successful fit."""
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow, PeakDetectionWorker, VoigtFittingWorker
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    (tmp_path / "ten.txt").write_text("\n".join(lines[:11]) + "\n")
    (tmp_path / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    set_ps2d_config('15N')
    window = LunaNMRMainWindow()
    try:
        window.load_data(str(tmp_path / SPECTRUM), str(tmp_path / "ten.txt"))
        window.sync_parameters_to_integrator()
        PeakDetectionWorker(window.integrator, window.workflow_mode).run()
        window.sync_parameters_to_integrator()
        worker = VoigtFittingWorker(window.integrator, window.param_manager,
                                    {'use_parallel': False, 'use_voigt_fitting': True,
                                     'skip_series_params': True})
        worker.fitting_complete.connect(lambda summary: window.on_fitting_complete(summary, None))
        worker.run()
        fit = window.current_voigt_result
        assert fit['fitting_quality'] in ('Excellent', 'Good', 'Fair', 'Poor')
        shown = _shown(window.voigt_2d_plotter)
        assert f"Quality: {fit['fitting_quality']}" in shown
        assert f"R² = {fit['r_squared']:.4f}" in shown
    finally:
        window.close()
        window.deleteLater()
        set_ps2d_config('15N')
