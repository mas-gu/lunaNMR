# ABOUTME: Startup, Reset, a nucleus switch and the CLI all start detection at the same search window.
# ABOUTME: Real LunaNMRMainWindow offscreen with its Expert panel open, and the CLI's series defaults.

"""
The search window starts at 0.070 ppm in 1H and at the indirect window of the
spectrum's nucleus: 0.3 ppm in 15N, which reaches a peak one 15N point (0.074 ppm at
600 MHz) off its reference, and the same in 13C linewidths (radF1 0.15 against 0.4).
Reset put the Expert panel back to 0.01/0.04, and a nucleus switch took 0.05/0.02 from
the PS2D config, so the window was the default only until the user touched either.
"""

import contextlib
import io
import os
import sys
import warnings
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

METHYL = Path(__file__).resolve().parents[2] / "data_test" / "methyl" / "sampleNBD"


def _default(nucleus):
    """(1H, indirect) default search window of a nucleus."""
    from lunaNMR.core.ps2d_config import SEARCH_WINDOW_PPM, PS2DConfig
    return SEARCH_WINDOW_PPM, PS2DConfig(nucleus).search_window_f1


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
    return {'state': value, 'panel': value, 'integrator': value}


def test_the_defaults():
    from lunaNMR.core.ps2d_config import PS2DConfig
    assert _default('15N') == (0.070, 0.3)
    # the 15N window in 13C linewidths, to the Expert panel's 0.01 ppm step
    proportional = 0.3 * PS2DConfig('13C').radF1 / PS2DConfig('15N').radF1
    assert _default('13C')[1] == pytest.approx(proportional, abs=0.005)


def test_startup_uses_the_default(window):
    window.sync_parameters_to_integrator()
    assert _windows(window) == _all_at(_default('15N'))


def test_reset_puts_the_default_back(window, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    window._on_nucleus_type_change('13C')
    window.expert_search_x_spin.setValue(0.03)
    window.expert_search_y_spin.setValue(0.50)
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.Yes))

    window.on_reset_results()

    assert _windows(window) == _all_at(_default('15N'))


@pytest.mark.parametrize("nucleus", ['13C', '15N'])
def test_a_nucleus_switch_puts_the_default_back(window, nucleus):
    window.expert_search_x_spin.setValue(0.03)
    window.expert_search_y_spin.setValue(0.50)

    window._on_nucleus_type_change(nucleus)

    assert _windows(window) == _all_at(_default(nucleus))


def test_the_cli_series_leaves_the_indirect_window_to_the_nucleus():
    from lunaNMR.cli import _default_series_params
    from lunaNMR.core.ps2d_config import SEARCH_WINDOW_PPM
    detection = _default_series_params()['detection_params']
    assert (detection['search_window_x'], detection['search_window_y']) == (SEARCH_WINDOW_PPM, None)


@pytest.mark.skipif(not METHYL.exists(), reason="data_test not present")
def test_a_cli_series_detects_with_the_window_of_its_nucleus(tmp_path):
    """A 13C series is detected in 13C linewidths although the CLI does not know the
    nucleus when it builds its parameters."""
    from numba.core.errors import NumbaWarning
    from lunaNMR.cli import _default_series_params
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    peaks = pd.read_csv(METHYL / "001_NBDADP_NBD.txt", skipinitialspace=True)
    peaks.columns = peaks.columns.str.strip()
    set_ps2d_config('15N')
    processor = MultiSpectrumProcessor(_default_series_params())
    try:
        with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore", NumbaWarning)
            processor.process_nmr_series([str(METHYL / "001_T2_NBD_3ms.ft")],
                                         peaks[['Assignment', 'Position_X', 'Position_Y']].head(3),
                                         str(tmp_path / "out"), peak_source_mode="reference")
        window = (processor.integrator.search_window_x, processor.integrator.search_window_y)
    finally:
        set_ps2d_config('15N')
    assert window == _default('13C')


def test_the_integrator_starts_at_the_default():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    integrator = EnhancedVoigtIntegrator()
    assert (integrator.search_window_x, integrator.search_window_y) == _default('15N')


def test_the_parameter_manager_starts_at_the_default():
    from lunaNMR.utils.parameter_manager import NMRParameterManager
    detection = NMRParameterManager().get_integrator_parameters()['detection_params']
    assert (detection['search_window_x'], detection['search_window_y']) == _default('15N')


def test_the_preflight_check_matches_the_default():
    """spectra_check keeps its own literal so that `diagnose` does not import lunaNMR.core.
    Its registration is a 15N check."""
    from lunaNMR.validation.spectra_check import SEARCH_WINDOW
    assert SEARCH_WINDOW == _default('15N')
