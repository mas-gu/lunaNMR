# ABOUTME: The module README's statement on matplotlib backends matches the code: no fitter
# ABOUTME: or plot script selects a backend, and the README does not say one does.

import re
from pathlib import Path

_DIR = Path(__file__).resolve().parent.parent


def test_the_readme_says_what_the_code_does_about_backends():
    selecting = [str(path.relative_to(_DIR)) for path in _DIR.rglob("*.py")
                 if "tests" not in path.parts
                 and re.search(r"^\s*(matplotlib|mpl)\.use\(", path.read_text(), re.M)]
    assert selecting == []
    readme = (_DIR / "README.md").read_text()
    assert "select matplotlib's non-interactive Agg backend" not in readme
    assert "select no matplotlib backend" in readme
