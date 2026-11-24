# ABOUTME: The documentation states the one supported way to install lunaNMR: from its
# ABOUTME: clone, run in place or installed in editable mode, never a regular install.

import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
INSTALL_DOCUMENTS = ["README.md", "AGENTS.md", "docs/INSTALLATION.md", "docs/README.md",
                     "docs/QUICK_REFERENCE.md", "docs/CLI.md", "docs/CLI_AGENT.md",
                     "docs/MCP.md"]
# `pip install .`, `pip install ".[mcp]"`, `pip install lunaNMR`: the project itself,
# installed as a copy.
_REGULAR_INSTALL = re.compile(
    r"pip3? install\s+(?!-e\b|--editable\b)(?:\.(?:\s|$)|[\"']\.\[|lunaNMR\b)")


def _tracked_documents():
    listed = subprocess.run(["git", "ls-files", "*.md", "*.html"], cwd=REPO,
                            capture_output=True, text=True)
    if listed.returncode != 0:
        pytest.skip("not a git checkout")
    return [REPO / name for name in listed.stdout.split()]


def _command_lines(text):
    """The lines of a document's fenced code blocks: what a reader is told to run."""
    fenced = False
    for line in text.splitlines():
        if line.lstrip().startswith("```"):
            fenced = not fenced
        elif fenced:
            yield line


@pytest.mark.parametrize("name", INSTALL_DOCUMENTS)
def test_it_names_the_editable_install_and_why(name):
    text = (REPO / name).read_text(encoding="utf-8")
    assert "pip install -e" in text
    assert "editable" in text
    assert re.search(r"`pip install \.`[^.\n]*\b(not|leaves|without|does not)\b", text), (
        f"{name} does not say that a regular `pip install .` is not enough")
    assert "`modules/`" in text or "`modules`" in text


def test_the_installation_guide_gives_both_supported_ways():
    text = (REPO / "docs" / "INSTALLATION.md").read_text(encoding="utf-8")
    commands = list(_command_lines(text))
    assert any(line.startswith("pip install -e .") for line in commands)
    assert "python3 launch_lunaNMR.py" in text
    assert "not supported" in text


def test_no_document_tells_the_reader_to_run_a_regular_install():
    offered = [f"{doc.relative_to(REPO)}: {line.strip()}"
               for doc in _tracked_documents()
               for line in _command_lines(doc.read_text(encoding="utf-8"))
               if _REGULAR_INSTALL.search(line)]
    assert offered == []


@pytest.mark.parametrize("line,regular", [
    ("pip install .", True),
    ('pip install ".[mcp]"', True),
    ("pip install lunaNMR", True),
    ("pip install -e .", False),
    ('pip install -e ".[mcp]" && lunanmr-mcp --check', False),
    ("pip install --editable .", False),
    ("pip install -r requirements.txt", False),
    ("pip install .venv/wheels/numba.whl", False),
])
def test_a_regular_install_line_is_told_from_an_editable_one(line, regular):
    assert bool(_REGULAR_INSTALL.search(line)) is regular
