# ABOUTME: The Expert Mode dialog offers only controls the engine reads, and its tooltips say what the engine does.
# ABOUTME: Real LunaNMRMainWindow offscreen; the dialog is built and its widgets inspected.

"""No engine code reads custom initial linewidths or a detection rectangle in pixels, and
a detection started from the window always takes the contour threshold over a peak-list
S/N threshold: the dialog holds no control for any of the three. Tooltips are checked
against core_integrator._add_dummy_peaks (10th percentile, one PS2D radius) and against
the series, which never runs the adaptive radius search."""

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
    """A window whose Expert Mode dialog has been opened."""
    from PySide6.QtWidgets import QDialog
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    w = LunaNMRMainWindow()
    monkeypatch.setattr(QDialog, "exec", lambda self: 0)
    w.on_expert_mode()
    yield w
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _texts(window, widget_type):
    return [widget.text() for widget in window._expert_dialog.findChildren(widget_type)]


def test_the_dialog_offers_no_custom_initial_linewidths(window):
    from PySide6.QtWidgets import QCheckBox, QLabel
    assert "Custom Initial Linewidths" not in _texts(window, QCheckBox)
    assert not {"Lor:", "Gauss:"} & set(_texts(window, QLabel))
    for name in ('use_custom_linewidths', 'lw_lorentz_1h', 'lw_gauss_1h', 'lw_lorentz_15n',
                 'lw_gauss_15n', 'custom_lw_visible', 'expert_custom_lw_check'):
        assert not hasattr(window, name), name


def test_the_dialog_offers_no_detection_rectangle_in_pixels(window):
    from PySide6.QtWidgets import QLabel
    assert not [text for text in _texts(window, QLabel) if "pixels" in text]
    for name in ('detection_square_size', 'detection_rectangle_y', 'detection_square_ppm_x',
                 'detection_square_ppm_y', 'expert_det_x_spin', 'expert_det_y_spin'):
        assert not hasattr(window, name), name


def test_the_dialog_offers_no_peak_list_sn_threshold(window):
    from PySide6.QtWidgets import QLabel
    assert not [text for text in _texts(window, QLabel) if "S/N" in text]
    for name in ('peak_list_sn_threshold', 'expert_peak_list_sn_spin'):
        assert not hasattr(window, name), name


def test_the_integrator_keeps_its_own_peak_list_sn_threshold(window):
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    window.sync_parameters_to_integrator()
    assert (window.integrator.peak_list_sn_threshold
            == EnhancedVoigtIntegrator().peak_list_sn_threshold)


def test_the_dialog_has_no_group_without_a_control(window):
    from PySide6.QtWidgets import (QAbstractButton, QAbstractSpinBox, QComboBox, QGroupBox)
    empty = [group.title() for group in window._expert_dialog.findChildren(QGroupBox)
             if not group.findChildren(QAbstractButton) and not group.findChildren(QAbstractSpinBox)
             and not group.findChildren(QComboBox)]
    assert empty == []


def test_the_dummy_peak_tooltip_states_the_filters_detection_applies(window):
    tooltip = window.expert_auto_add_dummy_cb.toolTip()
    assert "10th percentile" in tooltip and "one PS2D radius" in tooltip
    assert "3rd percentile" not in tooltip and "3×" not in tooltip


def test_no_series_tooltip_promises_an_adaptive_search_per_spectrum(window):
    for box in (window.expert_rerun_adaptive_check, window.expert_lock_clusters_check):
        tooltip = box.toolTip()
        assert "grid search" not in tooltip and "adaptive optimization" not in tooltip, box.text()
        assert "linewidth statistics" in tooltip, box.text()
