# ABOUTME: Importing lunaNMR.utils loads no tkinter: the GUI is PySide6, and the CLI runs with no display.
# ABOUTME: Checked in a fresh interpreter, since any earlier test may already have imported tkinter.

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_importing_utils_loads_no_tkinter():
    code = ("import sys\n"
            "import lunaNMR.utils\n"
            "print(sorted(m for m in sys.modules if 'tkinter' in m))\n")
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=ROOT)
    assert p.returncode == 0, p.stderr[-600:]
    assert p.stdout.strip().splitlines()[-1] == "[]"
