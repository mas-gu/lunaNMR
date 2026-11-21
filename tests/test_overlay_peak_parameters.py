# ABOUTME: Tests the overlay viewer's Peak Parameters text on a real series fit.
# ABOUTME: The text must not name a colour per peak that the fit plots do not draw it in.

import sys

import pytest

from PySide6.QtWidgets import QApplication


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


def test_peak_parameters_name_no_colour_the_plots_do_not_draw(app, dnaja1_titration_series,
                                                             open_overlay_viewer):
    """The fit plots colour a peak by assignment from the plotter's own palette; the
    parameter text is plain and carries no colour of its own."""
    from matplotlib.colors import CSS4_COLORS
    viewer = open_overlay_viewer(dnaja1_titration_series)
    try:
        fit = next(p for p in viewer._get_fitted_peaks(0) if len(p.get('all_peaks') or []) > 1)
        viewer.update_peak_parameters(fit)
        text = viewer.peak_params_text.toPlainText()
        headers = [line for line in text.splitlines()
                   if any(line.startswith(str(p['assignment'])) for p in fit['all_peaks'])]
        assert len(headers) == len(fit['all_peaks'])
        assert not [line for line in headers if any(f"({name})" in line for name in CSS4_COLORS)]
    finally:
        viewer.deleteLater()
