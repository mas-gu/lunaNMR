# ABOUTME: Startup, Reset, a nucleus switch and the CLI all start detection at the same search window.
# ABOUTME: Real LunaNMRMainWindow offscreen with its Expert panel open, and the CLI's series defaults.

"""
The search window starts at 0.070 ppm in both dimensions: a wider one in the indirect
dimension lets a low-S/N peak fit to noise. Reset put the Expert panel back to
0.01/0.04, and a nucleus switch took 0.05/0.02 from the PS2D config, so the window was
0.070 only until the user touched either.
"""

import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def window(app, monkeypatch):
    """The main window with its Expert panel built (its modal exec answered at once)."""
    from PySide6.QtWidgets import QDialog
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    monkeypatch.setattr(QDialog, "exec", lambda self: 0)
    w = LunaNMRMainWindow()
    w.on_expert_mode()
    yield w
    app.processEvents()
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _windows(window):
    """The search window as the state, the Expert panel and the integrator each hold it."""
    return {
        'state': (window.search_window_x, window.search_window_y),
        'panel': (window.expert_search_x_spin.value(), window.expert_search_y_spin.value()),
        'integrator': (window.integrator.search_window_x, window.integrator.search_window_y),
    }


def _all_at(value):
    return {'state': (value, value), 'panel': (value, value), 'integrator': (value, value)}


def test_startup_uses_the_default(window):
    from lunaNMR.core.ps2d_config import SEARCH_WINDOW_PPM
    window.sync_parameters_to_integrator()
    assert SEARCH_WINDOW_PPM == 0.070
    assert _windows(window) == _all_at(SEARCH_WINDOW_PPM)


def test_reset_puts_the_default_back(window, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    from lunaNMR.core.ps2d_config import SEARCH_WINDOW_PPM
    window.expert_search_x_spin.setValue(0.03)
    window.expert_search_y_spin.setValue(0.30)
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.Yes))

    window.on_reset_results()

    assert _windows(window) == _all_at(SEARCH_WINDOW_PPM)


@pytest.mark.parametrize("nucleus", ['13C', '15N'])
def test_a_nucleus_switch_puts_the_default_back(window, nucleus):
    from lunaNMR.core.ps2d_config import SEARCH_WINDOW_PPM
    window.expert_search_x_spin.setValue(0.03)
    window.expert_search_y_spin.setValue(0.30)

    window._on_nucleus_type_change(nucleus)

    assert _windows(window) == _all_at(SEARCH_WINDOW_PPM)


def test_the_cli_series_uses_the_default():
    from lunaNMR.cli import _default_series_params
    from lunaNMR.core.ps2d_config import SEARCH_WINDOW_PPM
    detection = _default_series_params()['detection_params']
    assert (detection['search_window_x'], detection['search_window_y']) == (SEARCH_WINDOW_PPM,) * 2


def test_the_integrator_starts_at_the_default():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.ps2d_config import SEARCH_WINDOW_PPM
    integrator = EnhancedVoigtIntegrator()
    assert (integrator.search_window_x, integrator.search_window_y) == (SEARCH_WINDOW_PPM,) * 2


def test_the_parameter_manager_starts_at_the_default():
    from lunaNMR.core.ps2d_config import SEARCH_WINDOW_PPM
    from lunaNMR.utils.parameter_manager import NMRParameterManager
    detection = NMRParameterManager().get_integrator_parameters()['detection_params']
    assert (detection['search_window_x'], detection['search_window_y']) == (SEARCH_WINDOW_PPM,) * 2


def test_the_preflight_check_matches_the_default():
    """spectra_check keeps its own literal so that `diagnose` does not import lunaNMR.core."""
    from lunaNMR.core.ps2d_config import SEARCH_WINDOW_PPM
    from lunaNMR.validation.spectra_check import SEARCH_WINDOW
    assert SEARCH_WINDOW == SEARCH_WINDOW_PPM
