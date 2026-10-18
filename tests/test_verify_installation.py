# ABOUTME: Tests that the installation verifier reports success on a working checkout
# ABOUTME: and that every module and launch command it names actually exists.
"""The verifier is the first thing a new user runs, per README.md and docs/INSTALLATION.md.

It checked the flat v0.9 module layout (`core_integrator`, `spectrum_browser`,
`main_gui`, …) long after v1.0 moved everything into packages, so it reported
`INSTALLATION INCOMPLETE` on a healthy tree. These tests pin the module list to
the real layout and pin the exit status to the health of the checkout.
"""

import ast
import importlib
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "lunaNMR" / "validation" / "verify_installation.py"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _run(cwd):
    env = {**os.environ, "QT_QPA_PLATFORM": "offscreen"}
    return subprocess.run([sys.executable, str(SCRIPT)],
                          cwd=str(cwd), capture_output=True, text=True, env=env)


@pytest.fixture(scope="module")
def verifier():
    return importlib.import_module("lunaNMR.validation.verify_installation")


def test_every_internal_module_it_checks_is_importable(verifier):
    """The module list must track the real package layout, not a past one."""
    unimportable = []
    for module_name, _description in verifier.INTERNAL_MODULES:
        try:
            importlib.import_module(module_name)
        except ImportError as exc:
            unimportable.append(f"{module_name} ({exc})")
    assert not unimportable, (
        "verify_installation checks modules that cannot be imported: "
        + ", ".join(unimportable)
    )


def test_advertised_launch_command_exists(verifier):
    assert (REPO_ROOT / verifier.LAUNCH_SCRIPT).is_file(), (
        f"verifier tells users to run {verifier.LAUNCH_SCRIPT}, which does not exist"
    )


def test_reports_success_on_this_checkout():
    result = _run(REPO_ROOT)
    assert result.returncode == 0, (
        "verifier failed on a working checkout:\n" + result.stdout[-2000:]
    )
    assert "INSTALLATION COMPLETE" in result.stdout


def test_runs_from_any_working_directory(tmp_path):
    """README invocation puts only the script's own directory on sys.path."""
    result = _run(tmp_path)
    assert result.returncode == 0, (
        "verifier depends on the working directory:\n" + result.stdout[-2000:]
    )


def test_reports_missing_dependency(tmp_path):
    """A genuinely broken environment must still be reported as incomplete."""
    stub = tmp_path / "nmrglue.py"
    stub.write_text("raise ImportError('simulated missing install')\n", encoding="utf-8")
    env = {**os.environ, "QT_QPA_PLATFORM": "offscreen", "PYTHONPATH": str(tmp_path)}
    result = subprocess.run([sys.executable, str(SCRIPT)], cwd=str(REPO_ROOT),
                            capture_output=True, text=True, env=env)
    assert result.returncode == 1
    assert "INSTALLATION INCOMPLETE" in result.stdout
    assert "nmrglue" in result.stdout


def _declared_in_pyproject():
    tomllib = pytest.importorskip("tomllib")
    with open(REPO_ROOT / "pyproject.toml", "rb") as fh:
        deps = tomllib.load(fh)["project"]["dependencies"]
    return {re.split(r"[<>=!~ ;\[]", dep, maxsplit=1)[0] for dep in deps}


def _declared_in_requirements():
    names = set()
    for line in (REPO_ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            names.add(re.split(r"[<>=!~ ;\[]", line, maxsplit=1)[0])
    return names


def _third_party_imports():
    """Top-level names imported anywhere in the shipped code, minus stdlib and local files."""
    sources = [*(REPO_ROOT / "lunaNMR").rglob("*.py"), *(REPO_ROOT / "modules").rglob("*.py"),
               REPO_ROOT / "launch_lunaNMR.py"]
    local = {p.stem for p in sources} | {p.parent.name for p in sources}
    found = set()
    for path in sources:
        tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                found.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                found.add(node.module.split(".")[0])
    return found - set(sys.stdlib_module_names) - local - {"pytest"}


def test_dependency_declarations_agree(verifier):
    required = set(verifier.REQUIRED_DEPENDENCIES)
    assert _declared_in_pyproject() == required
    assert _declared_in_requirements() == required | set(verifier.OPTIONAL_DEPENDENCIES)


def test_every_imported_library_is_declared(verifier):
    declared = set(verifier.REQUIRED_DEPENDENCIES) | set(verifier.OPTIONAL_DEPENDENCIES)
    undeclared = _third_party_imports() - declared
    assert not undeclared, f"imported but not declared: {sorted(undeclared)}"


def test_every_declared_library_is_imported(verifier):
    declared = set(verifier.REQUIRED_DEPENDENCIES) | set(verifier.OPTIONAL_DEPENDENCIES)
    unused = declared - _third_party_imports()
    assert not unused, f"declared but imported nowhere: {sorted(unused)}"
