# ABOUTME: SpectrumPlotter emits peak_left_click_requested on a no-drag left click.
# ABOUTME: This drives the main-window two-click select->move editing.

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from PySide6.QtWidgets import QApplication
from lunaNMR.gui.components.spectrum_plotter import SpectrumPlotter


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


def test_left_click_emits_signal(app):
    p = SpectrumPlotter()
    try:
        # the nav handler's left-click hook is wired to the plotter's emitter
        assert p._nav_handler.on_left_click == p._on_left_click_request
        got = []
        p.peak_left_click_requested.connect(lambda x, y: got.append((x, y)))
        p._nav_handler.on_left_click(8.0, 120.0)   # a no-drag left click
        assert got == [(8.0, 120.0)]
    finally:
        p.deleteLater()


def test_middle_click_emits_select_signal(app):
    p = SpectrumPlotter()
    try:
        assert p._nav_handler.on_peak_select == p._on_peak_select_request
        got = []
        p.peak_select_requested.connect(lambda x, y: got.append((x, y)))
        p._nav_handler.on_peak_select(7.5, 119.0)   # a no-drag middle click
        assert got == [(7.5, 119.0)]
    finally:
        p.deleteLater()
