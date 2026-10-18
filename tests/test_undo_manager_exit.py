# ABOUTME: A process that moved a peak exits without Qt errors, and the undo signals still reach listeners.
# ABOUTME: Real PeakUndoManager on a stand-alone object; the exit check runs in a fresh interpreter.

import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

PUSH_THEN_EXIT = """
import sys
from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QUndoCommand
from lunaNMR.utils.undo_manager import PeakUndoManager
app = QApplication(sys.argv)
manager = PeakUndoManager(main_window=None)
manager.push(QUndoCommand("move"))
print("pushed")
"""


def test_exit_after_an_undoable_edit_reports_no_error():
    run = subprocess.run([sys.executable, "-c", PUSH_THEN_EXIT], cwd=REPO_ROOT,
                         capture_output=True, text=True, timeout=120,
                         env={**os.environ, "QT_QPA_PLATFORM": "offscreen"})

    assert run.returncode == 0, run.stderr
    assert run.stdout.strip() == "pushed"
    assert "Signal source has been deleted" not in run.stderr


def test_the_undo_and_redo_availability_reaches_listeners():
    from PySide6.QtWidgets import QApplication
    from PySide6.QtGui import QUndoCommand
    from lunaNMR.utils.undo_manager import PeakUndoManager
    QApplication.instance() or QApplication(sys.argv)
    manager = PeakUndoManager(main_window=None)
    can_undo, can_redo = [], []
    manager.can_undo_changed.connect(can_undo.append)
    manager.can_redo_changed.connect(can_redo.append)

    manager.push(QUndoCommand("move"))
    manager.undo_stack.undo()

    assert can_undo == [True, False]
    assert can_redo == [False, True]
