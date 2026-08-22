# ABOUTME: Guards every test root against a test class or function name defined
# ABOUTME: twice in one module, where the second binding silently discards the first.
"""A duplicate top-level name in a test module costs coverage with no error.

Python binds the later definition over the earlier one, so pytest only ever
collects the survivor. `test_kd_models.py` carried two `TestObservables`
classes this way, and the two binding-model tests in the first one had never
run. Parsed with `ast` rather than imported, so this stays independent of
whether a module's own dependencies are installed.
"""

import ast
from collections import Counter
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

TEST_ROOTS = (
    REPO_ROOT / "tests",
    REPO_ROOT / "modules" / "dynamiXs_v2o0" / "tests",
    REPO_ROOT / "modules" / "integration_1d_v1o0" / "tests",
)


def _test_modules():
    for root in TEST_ROOTS:
        if root.is_dir():
            yield from sorted(root.glob("test_*.py"))


def _duplicate_top_level_names(path):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    collected = [
        node.name
        for node in tree.body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
        and (node.name.startswith("Test") or node.name.startswith("test_"))
    ]
    return sorted(name for name, count in Counter(collected).items() if count > 1)


@pytest.mark.parametrize("module", _test_modules(), ids=lambda p: p.name)
def test_module_defines_each_test_name_once(module):
    duplicates = _duplicate_top_level_names(module)
    assert not duplicates, (
        f"{module.relative_to(REPO_ROOT)} defines {duplicates} more than once; "
        "the later definition shadows the earlier, so those tests never run"
    )
