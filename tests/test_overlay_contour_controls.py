# ABOUTME: Tests the overlay viewer's contour controls on a real series: each control offered
# ABOUTME: in the Contour group must change the contour levels the overlay draws.

import sys

import pytest

from PySide6.QtWidgets import QAbstractSpinBox, QApplication, QGroupBox, QLabel


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def viewer(app, dnaja1_titration_series, open_overlay_viewer):
    viewer = open_overlay_viewer(dnaja1_titration_series)
    yield viewer
    viewer.deleteLater()


def _drawn_levels(viewer):
    viewer._update_plot()
    return [tuple(c.levels) for c in viewer.plot_widget.axes.collections if hasattr(c, 'levels')]


def test_every_contour_control_changes_the_contours_drawn(viewer):
    group = next(g for g in viewer.findChildren(QGroupBox) if g.title() == "Contour")
    names = [label.text() for label in group.findChildren(QLabel)]
    controls = group.findChildren(QAbstractSpinBox)
    assert controls and len(names) == len(controls)
    without_effect = []
    for name, control in zip(names, controls):
        before = _drawn_levels(viewer)
        assert before
        if control.value() == control.maximum():
            control.stepDown()
        else:
            control.stepUp()
        if _drawn_levels(viewer) == before:
            without_effect.append(name)
    assert not without_effect
