# ABOUTME: The DynamiXs Fit Viewer opens on its fits without a modal box announcing what it loaded.
# ABOUTME: Real JA1 600 MHz T1 matrix -> CLI t1t2 fit JSON -> FitViewer, offscreen.

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

MATRIX = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
          / "run_20260825_1405" / "series" / "600_T1" / "peak_intensity_matrix.csv")

pytestmark = pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture(scope="module")
def json_file(tmp_path_factory):
    from lunaNMR.cli import main
    out = tmp_path_factory.mktemp("t1")
    assert main(["dynamixs", "t1t2", "--input", str(MATRIX), "--out", str(out),
                 "--exp", "T1", "--field-freq", "600", "--prefix", "600_T1"]) == 0
    found = sorted(out.rglob("*_fit_data.json"))
    assert len(found) == 1
    return found[0]


def test_loading_fits_opens_no_message_box(app, json_file, monkeypatch):
    import lunaNMR.gui.dialogs.dynamixs_dialog  # noqa: F401 - puts the module on sys.path
    import visualization.fit_viewer as fit_viewer
    shown = []
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(fit_viewer, box, lambda *a, _box=box, **k: shown.append((_box, a[1:])))
    viewer = fit_viewer.FitViewer(json_files=[str(json_file)])
    try:
        assert viewer.data
        assert shown == []
    finally:
        viewer.close()
        viewer.deleteLater()
