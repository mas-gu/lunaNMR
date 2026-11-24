# ABOUTME: A Voigt view asked for its placeholder or its no-data message after a fit shows the message, not the fit.
# ABOUTME: Real VoigtAnalysisPlotter offscreen, on a real DNAJA1 fit drawn in the 2D view and in the 3D view.

"""
Drawing a fit clears the figure and builds axes of its own. The placeholder and the 2D
no-data message were written on the axes the widget was built with, which the figure then
no longer held: the message never appeared and the fit stayed drawn.
"""

import contextlib
import io
import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")

VIEWS = pytest.mark.parametrize("view", ["plot_voigt_analysis", "plot_voigt_analysis_3d"])


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture(scope="module")
def fit():
    """The 2D fit of one isolated peak of the first DNAJA1 spectrum."""
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.parallel_voigt_processor import ParallelVoigtProcessor
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.utils.file_manager import NMRFileManager
    set_ps2d_config('15N')
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(DNAJA1 / "2D_DNAJA1_0o0.ft"))
        it.noise_level = it._estimate_noise_level_advanced()
        it.peak_list = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
        it.gui_params = {'use_ps2d_multi_peak': True, 'use_adaptive_optimization': False}
        context = ParallelVoigtProcessor(it.enhanced_fitter, in_process=True)._peaks_context(it.peak_list)
        listed = it.peak_list.set_index('Assignment').loc['3LysH']
        result = it.fit_peak_voigt_2d(float(listed['Position_X']), float(listed['Position_Y']),
                                      '3LysH', all_peaks_context=context)
    assert '2d_simultaneous' in result['method'] and result['fitted_2d_surface'] is not None
    return result


@pytest.fixture
def plotter(app):
    from lunaNMR.gui.components.voigt_analysis_plotter import VoigtAnalysisPlotter
    widget = VoigtAnalysisPlotter()
    yield widget
    widget.deleteLater()


def _messages(plotter):
    """The texts written in the axes the figure holds."""
    return [text.get_text() for ax in plotter.figure.axes for text in ax.texts]


def _layout(plotter):
    """Where each axes of the figure sits."""
    return [tuple(round(edge, 6) for edge in ax.get_position().bounds) for ax in plotter.figure.axes]


def _draws_a_fit(plotter):
    """Whether the figure holds a title, data, or an axes that is not switched off."""
    title = plotter.figure._suptitle
    return bool((title is not None and title.get_text())
                or any(ax.has_data() or ax.axison for ax in plotter.figure.axes))


@VIEWS
def test_the_fit_is_drawn_before_anything_replaces_it(plotter, fit, view):
    getattr(plotter, view)(fit)

    assert _draws_a_fit(plotter)
    assert "3LysH" in plotter.figure._suptitle.get_text()


@VIEWS
def test_the_placeholder_replaces_a_drawn_fit(plotter, fit, view):
    plotter.show_placeholder("Select a peak to view Voigt analysis")
    before_any_fit = _messages(plotter), _layout(plotter)
    getattr(plotter, view)(fit)

    plotter.show_placeholder("Select a peak to view Voigt analysis")

    assert (_messages(plotter), _layout(plotter)) == before_any_fit
    assert not _draws_a_fit(plotter)
    assert plotter.current_result is None


@VIEWS
def test_the_no_data_message_replaces_a_drawn_fit(plotter, fit, view):
    getattr(plotter, view)(None)
    before_any_fit = _messages(plotter)
    getattr(plotter, view)(fit)

    getattr(plotter, view)(None)

    assert len(_messages(plotter)) == 1 and "No Voigt analysis data" in _messages(plotter)[0]
    assert _messages(plotter) == before_any_fit
    assert not _draws_a_fit(plotter)


@VIEWS
def test_a_fit_is_drawn_again_after_the_placeholder(plotter, fit, view):
    getattr(plotter, view)(fit)
    plotter.show_placeholder("Select a peak to view Voigt analysis")

    getattr(plotter, view)(fit)

    assert _draws_a_fit(plotter)
    assert "Select a peak to view Voigt analysis" not in _messages(plotter)


def test_a_new_view_shows_its_placeholder_as_before_any_fit(plotter):
    plotter.show_placeholder("Select a peak to view Voigt analysis")

    assert _messages(plotter) == ["Select a peak to view Voigt analysis"] * 4
    assert not _draws_a_fit(plotter)
