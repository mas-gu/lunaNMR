# ABOUTME: Every place lunaNMR states its own version must carry the same number.
# ABOUTME: The release version is at least 1.1.0, the first with the MCP server.

import re
from pathlib import Path

import lunaNMR
from lunaNMR.batch_processing import __version__ as batch_version

ROOT = Path(__file__).resolve().parents[1]


def _first_match(relative_path, pattern):
    return re.search(pattern, (ROOT / relative_path).read_text(encoding="utf-8")).group(1)


def _version_tuple(version):
    return tuple(int(part) for part in version.split("."))


def test_every_stated_version_is_the_same_number():
    release = _first_match("pyproject.toml", r'(?m)^version = "([\d.]+)"')
    stated = {
        "lunaNMR.__version__": lunaNMR.__version__,
        "batch_processing.__version__": batch_version,
        "top-level __init__.py": _first_match("__init__.py", r'__version__ = "([\d.]+)"'),
        "batch --version text": _first_match(
            "lunaNMR/batch_processing/cli_interface.py", r"LunaNMR Batch Processing v([\d.]+)"),
        "About dialog": _first_match("lunaNMR/gui/main_window.py", r"Version ([\d.]+) \(Qt Port\)"),
        "config default": _first_match("lunaNMR/utils/config_manager.py", r'"version": "([\d.]+)"'),
    }
    assert {name: version for name, version in stated.items() if version != release} == {}


def test_the_release_includes_the_mcp_server():
    assert _version_tuple(lunaNMR.__version__) >= (1, 1, 0)
