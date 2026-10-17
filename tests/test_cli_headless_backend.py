# ABOUTME: The CLI must never open a plot window, whatever matplotlib backend its process has.
# ABOUTME: In a process on the Qt backend, an unguarded fit call blocks in plt.show() forever.

"""
The fitters show their figures when the backend is interactive. From a shell the backend
is non-interactive and nothing shows; in a process that already uses Qt - a test session
after a GUI test, or an in-process caller such as an MCP server - the backend is QtAgg,
and `plt.show()` starts an event loop that never returns: the call hangs.
"""

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

MATRIX = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
          / "run_20260825_1405" / "series" / "600_T1" / "peak_intensity_matrix.csv")


@pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")
def test_a_fit_run_in_a_qt_process_returns(tmp_path):
    import matplotlib
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication(sys.argv)
    matplotlib.use("QtAgg")
    from lunaNMR.cli import main
    assert main(["dynamixs", "t1t2", "--input", str(MATRIX), "--out", str(tmp_path),
                 "--exp", "T1", "--field-freq", "600", "--prefix", "t"]) == 0
    assert matplotlib.get_backend().lower() == "agg"
