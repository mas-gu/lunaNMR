# ABOUTME: The launcher starts from a checkout whatever its folder is called, including
# ABOUTME: `lunaNMR`, the name `git clone` gives it.

import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TOP_LEVEL = [name for name in ("launch_lunaNMR.py", "__init__.py", "lunaNMR", "modules")
             if (REPO / name).exists()]


def _checkout_named(tmp_path, name):
    checkout = tmp_path / name
    checkout.mkdir()
    for entry in TOP_LEVEL:
        (checkout / entry).symlink_to(REPO / entry)
    return checkout


def _import_the_gui_as_the_launcher_does(checkout):
    script = (
        "import importlib.util, sys\n"
        f"spec = importlib.util.spec_from_file_location('launcher', r'{checkout / 'launch_lunaNMR.py'}')\n"
        "launcher = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(launcher)\n"
        "launcher.setup_paths()\n"
        "from lunaNMR.gui.main_window import main\n"
        "import lunaNMR\n"
        "print(lunaNMR.__file__)\n"
    )
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen")
    env.pop("PYTHONPATH", None)
    return subprocess.run([sys.executable, "-c", script], cwd=str(checkout.parent), env=env,
                          capture_output=True, text=True, timeout=300)


def test_a_checkout_named_lunaNMR_imports_its_own_gui(tmp_path):
    checkout = _checkout_named(tmp_path, "lunaNMR")
    done = _import_the_gui_as_the_launcher_does(checkout)
    assert done.returncode == 0, done.stderr[-2000:]
    assert Path(done.stdout.strip().splitlines()[-1]) == checkout / "lunaNMR" / "__init__.py"


def test_a_checkout_under_any_other_name_imports_its_own_gui(tmp_path):
    checkout = _checkout_named(tmp_path, "somewhere_else")
    done = _import_the_gui_as_the_launcher_does(checkout)
    assert done.returncode == 0, done.stderr[-2000:]
    assert Path(done.stdout.strip().splitlines()[-1]) == checkout / "lunaNMR" / "__init__.py"
