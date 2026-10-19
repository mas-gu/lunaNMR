# ABOUTME: Tests that launch_lunaNMR.py builds and shows the main window through its own
# ABOUTME: start-up path, and that the launcher and modules package import only what exists.
import ast
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

# Run in a child process: the launcher forces the QtAgg backend and ends in os._exit.
# QApplication.exec is replaced so the event loop returns once the main window is up.
_LAUNCH = """
import sys
from PySide6.QtWidgets import QApplication

def _exec_until_main_window(self):
    for _ in range(200):
        self.processEvents()
        shown = [w for w in self.topLevelWidgets()
                 if type(w).__name__ == "LunaNMRMainWindow" and w.isVisible()]
        if shown:
            import logging
            root = logging.getLogger()
            formats = [h.formatter._fmt for h in root.handlers if h.formatter]
            print("MAIN_WINDOW_SHOWN", flush=True)
            print(f"ROOT_LOGGING {logging.getLevelName(root.level)} {formats}", flush=True)
            for w in self.topLevelWidgets():
                w.close()
            return 0
    return 1

QApplication.exec = _exec_until_main_window
sys.argv = ["launch_lunaNMR.py"]
import runpy
runpy.run_path("launch_lunaNMR.py", run_name="__main__")
"""


def test_launcher_shows_the_main_window():
    env = {**os.environ, "QT_QPA_PLATFORM": "offscreen"}
    result = subprocess.run([sys.executable, "-c", _LAUNCH], cwd=str(REPO_ROOT),
                            capture_output=True, text=True, env=env, timeout=300)
    assert "MAIN_WINDOW_SHOWN" in result.stdout, result.stdout[-2000:] + result.stderr[-2000:]
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]


def test_the_launcher_sends_info_logging_to_the_terminal():
    """The terminal the app is started from shows the INFO records of lunaNMR, in the
    logging module's basic format."""
    import logging
    env = {**os.environ, "QT_QPA_PLATFORM": "offscreen"}
    result = subprocess.run([sys.executable, "-c", _LAUNCH], cwd=str(REPO_ROOT),
                            capture_output=True, text=True, env=env, timeout=300)
    assert f"ROOT_LOGGING INFO {[logging.BASIC_FORMAT]}" in result.stdout, \
        result.stdout[-2000:] + result.stderr[-2000:]


def test_importing_the_utilities_configures_no_logging():
    """Configuring the root logger is the application's choice, not a side effect of an
    import: a host importing lunaNMR keeps the logging it had."""
    probe = ("import logging, lunaNMR.utils, lunaNMR.utils.global_optimization_manager; "
             "print(len(logging.getLogger().handlers))")
    result = subprocess.run([sys.executable, "-c", probe], cwd=str(REPO_ROOT),
                            capture_output=True, text=True, timeout=300)
    assert result.stdout.strip() == "0", result.stdout[-2000:] + result.stderr[-2000:]


def _absolute_imports(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            yield from (alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            prefix = "modules." if node.level and path.parent.name == "modules" else ""
            yield prefix + node.module


@pytest.mark.parametrize("source", ["launch_lunaNMR.py", "modules/__init__.py"])
def test_every_module_named_exists(source, monkeypatch):
    monkeypatch.syspath_prepend(str(REPO_ROOT))
    missing = [name for name in _absolute_imports(REPO_ROOT / source)
               if importlib.util.find_spec(name) is None]
    assert not missing, f"{source} imports modules that do not exist: {missing}"
