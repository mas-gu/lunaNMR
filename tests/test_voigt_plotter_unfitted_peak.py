# ABOUTME: Tests what the Voigt analysis views show for a peak whose 2D fit failed, on a real failed fit.
# ABOUTME: Every fit is 2D: the views must say there is no fit to display and why, never name a 1D fit.

import contextlib
import io
import sys
from pathlib import Path

import pytest

from PySide6.QtWidgets import QApplication

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture(scope="module")
def failed_fit():
    """The failure the integrator returns for a peak it is given no peak list to place."""
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.ps2d_config import set_ps2d_config
    set_ps2d_config('15N')
    integrator = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        integrator.load_nmr_file(str(DNAJA1 / "2D_DNAJA1_0o0.ft"))
        fit = integrator.fit_peak_voigt_2d(8.225, 126.271, '3LysH')
    assert fit['success'] is False and fit['error']
    return fit


def _shown(plotter):
    """Every text the plotter's figure draws."""
    texts = [text.get_text() for ax in plotter.figure.axes for text in ax.texts]
    texts += [ax.get_title() for ax in plotter.figure.axes]
    texts += [text.get_text() for text in plotter.figure.texts]
    return "\n".join(texts)


@pytest.mark.parametrize("view", ["plot_voigt_analysis", "plot_voigt_analysis_3d"])
def test_a_failed_fit_is_shown_as_not_fitted_with_its_reason(app, failed_fit, view):
    from lunaNMR.gui.components.voigt_analysis_plotter import VoigtAnalysisPlotter
    plotter = VoigtAnalysisPlotter()
    try:
        getattr(plotter, view)(failed_fit)
        shown = _shown(plotter)
        assert "3LysH" in shown
        assert "1D" not in shown and "Dimension Fit" not in shown
        assert "no 2D fit context" in " ".join(shown.split())
    finally:
        plotter.deleteLater()
