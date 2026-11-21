# ABOUTME: The shipped documentation names only folders and files a `git clone` of the
# ABOUTME: published repository contains: its root is `lunaNMR/`, with no parent above it.

import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
ABSENT_FROM_A_CLONE = ("lunaNMR_v1o0", "CLAUDE.md", "JOURNAL.md")


def _tracked_documents():
    listed = subprocess.run(["git", "ls-files", "*.md", "*.html"], cwd=REPO,
                            capture_output=True, text=True)
    if listed.returncode != 0:
        pytest.skip("not a git checkout")
    return [REPO / name for name in listed.stdout.split()]


@pytest.mark.parametrize("name", ABSENT_FROM_A_CLONE)
def test_no_document_names_what_a_clone_does_not_hold(name):
    naming = [f"{doc.relative_to(REPO)}:{number}"
              for doc in _tracked_documents()
              for number, line in enumerate(doc.read_text(encoding="utf-8").splitlines(), 1)
              if name in line]
    assert naming == []
