# ABOUTME: Tests grab-to-move hooks on the shared NMR navigation handler.
# ABOUTME: Left-press on a peak grabs+drags+places; empty space still pans.

import os
import sys
from pathlib import Path
from unittest.mock import Mock

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from PySide6.QtWidgets import QApplication
from lunaNMR.gui.components.nmr_navigation_handler import NMRNavigationHandler


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


class _Event:
    def __init__(self, button, xdata, ydata, ax, x=0, y=0):
        self.button = button
        self.xdata, self.ydata = xdata, ydata
        self.inaxes = ax
        self.x, self.y = x, y


def _handler():
    h = NMRNavigationHandler()
    h._ax = Mock()            # events carry inaxes=this; get_xlim/ylim auto-stub
    h._canvas = Mock()        # setCursor is a no-op
    return h


def test_left_click_without_drag_fires_on_left_click(app):
    h = _handler()
    clicks = []
    h.on_left_click = lambda x, y: clicks.append((x, y))
    # press and release at the same pixel = a click, not a pan
    h._on_mouse_press(_Event(1, 8.0, 120.0, h._ax, x=100, y=100))
    assert h._pan_active
    h._on_mouse_release(_Event(1, 8.0, 120.0, h._ax, x=101, y=101))
    assert clicks == [(8.0, 120.0)] and not h._pan_active


def test_left_drag_pans_and_does_not_fire_click(app):
    h = _handler()
    clicks = []
    h.on_left_click = lambda x, y: clicks.append((x, y))
    h._on_mouse_press(_Event(1, 8.0, 120.0, h._ax, x=100, y=100))
    # release far from press (a drag/pan) → no click
    h._on_mouse_release(_Event(1, 8.5, 121.0, h._ax, x=160, y=140))
    assert clicks == []


def test_no_click_hook_left_press_pans_unchanged(app):
    h = _handler()                                # on_left_click is None (other dialogs)
    h._on_mouse_press(_Event(1, 5.0, 100.0, h._ax, x=10, y=10))
    assert h._pan_active
    h._on_mouse_release(_Event(1, 5.0, 100.0, h._ax, x=10, y=10))
    assert not h._pan_active
