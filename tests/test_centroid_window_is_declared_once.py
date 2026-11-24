# ABOUTME: The centroid window (how far a detected maximum's sub-point position may lie from its point) has one value.
# ABOUTME: The parameter manager, the CLI series defaults, the main window and an integrator given no parameters agree.

import contextlib
import io
import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

SPECTRUM = (Path(__file__).resolve().parents[2] / "data_test" / "data_binding"
            / "DNAJA1_HSPA8" / "2D_DNAJA1_0o0.ft")


def _window(gui_params):
    return gui_params['centroid_window_x_ppm'], gui_params['centroid_window_y_ppm']


def test_the_parameter_manager_starts_with_it():
    from lunaNMR.core.ps2d_config import CENTROID_WINDOW_PPM
    from lunaNMR.utils.parameter_manager import NMRParameterManager
    assert _window(NMRParameterManager().get_integrator_parameters()['gui_params']) == CENTROID_WINDOW_PPM


def test_a_cli_series_starts_with_it():
    from lunaNMR.cli import _default_series_params
    from lunaNMR.core.ps2d_config import CENTROID_WINDOW_PPM
    assert _window(_default_series_params()['gui_params']) == CENTROID_WINDOW_PPM


def test_the_main_window_starts_with_it():
    from PySide6.QtWidgets import QApplication
    from lunaNMR.core.ps2d_config import CENTROID_WINDOW_PPM, set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    app = QApplication.instance() or QApplication(sys.argv)
    set_ps2d_config('15N')
    window = LunaNMRMainWindow()
    try:
        assert (window.centroid_window_x_ppm, window.centroid_window_y_ppm) == CENTROID_WINDOW_PPM
    finally:
        window.close()
        window.deleteLater()
        app.processEvents()
        set_ps2d_config('15N')


@pytest.mark.skipif(not SPECTRUM.exists(), reason="data_test not present")
def test_an_integrator_given_no_parameters_detects_as_one_given_it():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.ps2d_config import CENTROID_WINDOW_PPM
    integrator = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        integrator.load_nmr_file(str(SPECTRUM))
    integrator.noise_level = integrator._estimate_noise_level_advanced()

    def positions(gui_params):
        integrator.gui_params = gui_params
        return [(p['x_ppm'], p['y_ppm'])
                for p in integrator._detect_peaks_by_threshold(10 * integrator.noise_level)]

    assert positions({}) == positions({'centroid_window_x_ppm': CENTROID_WINDOW_PPM[0],
                                       'centroid_window_y_ppm': CENTROID_WINDOW_PPM[1]})
