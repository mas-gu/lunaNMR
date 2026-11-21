# ABOUTME: Every value the Expert panel lets the user set for a validated parameter lies inside
# ABOUTME: the parameter manager's range for it, so a permitted value never logs a range warning.

import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

SPIN_OF_PARAMETER = {
    'search_window_x': 'expert_search_x_spin',
    'search_window_y': 'expert_search_y_spin',
    'centroid_window_x_ppm': 'expert_centroid_x_spin',
    'centroid_window_y_ppm': 'expert_centroid_y_spin',
}


@pytest.fixture
def window(monkeypatch):
    from PySide6.QtWidgets import QApplication, QDialog
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    app = QApplication.instance() or QApplication(sys.argv)
    set_ps2d_config('15N')
    w = LunaNMRMainWindow()
    monkeypatch.setattr(QDialog, "exec", lambda self: 0)
    w.on_expert_mode()
    yield w
    app.processEvents()
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


@pytest.mark.parametrize("parameter", sorted(SPIN_OF_PARAMETER))
def test_the_spin_box_range_is_inside_the_validated_range(window, parameter):
    spin = getattr(window, SPIN_OF_PARAMETER[parameter])
    definition = window.param_manager.parameter_definitions[parameter]
    assert definition['min'] <= spin.minimum(), parameter
    assert spin.maximum() <= definition['max'], parameter
