# ABOUTME: What the methyl T2 page says about itself matches what it fits and opens: the
# ABOUTME: shared-amplitude model, and a summary viewer showing one quantity at a time.

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def page(app):
    from dynamiXs_gui import DynamiXsMainWindow
    window = DynamiXsMainWindow()
    yield window.methyl_t2_page
    window.deleteLater()


def test_the_summary_button_says_what_it_opens(page):
    """The summary viewer draws one quantity across residues, chosen among five."""
    from visualization.methyl_t2_fit_viewer import MethylT2FitViewer
    assert len(MethylT2FitViewer.SUMMARY_QUANTITIES) == 5
    tooltip = page.summary_btn.toolTip()
    for quantity in ("T2_avg", "ΔT2", "T2a", "T2b", "η"):
        assert quantity in tooltip
    assert "one quantity at a time" in tooltip
    assert "pre-enabled" not in tooltip


def test_the_page_describes_the_shared_amplitude_model():
    import methyl_t2_page
    source = Path(methyl_t2_page.__file__).read_text()
    about = source.splitlines()[0]
    assert "0.5*A*[exp(-t/T2a) + exp(-t/T2b)]" in about
    for absent in ("A_a", "A_b", "+ C"):
        assert absent not in about
        assert absent not in methyl_t2_page.MethylT2FittingPage.__doc__
