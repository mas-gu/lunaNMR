# ABOUTME: The `lunanmr-mcp --check` self-test that tells a user whether their install can run the server.
# ABOUTME: Real imports and a real CLI call; a missing module is named, never a bare traceback.

import subprocess
import sys

import pytest

pytest.importorskip("mcp")

from lunaNMR.mcp.check import REQUIRED_MODULES, run_checks


def test_every_check_passes_in_a_working_environment():
    results = run_checks()
    assert all(ok for _, ok, _ in results), results
    names = [name for name, _, _ in results]
    assert "cli" in names
    assert set(REQUIRED_MODULES) <= set(names)


def test_a_missing_module_is_reported_by_name_and_fails():
    results = run_checks(modules=["numpy", "no_such_module_xyz"])
    by_name = {name: (ok, detail) for name, ok, detail in results}
    assert by_name["numpy"][0] is True
    assert by_name["no_such_module_xyz"][0] is False
    assert "no_such_module_xyz" in by_name["no_such_module_xyz"][1]


def test_the_command_line_check_needs_no_root_and_exits_zero():
    completed = subprocess.run([sys.executable, "-m", "lunaNMR.mcp.server", "--check"],
                               capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "ok   cli" in completed.stdout


def test_the_server_without_a_root_or_check_is_a_usage_error():
    completed = subprocess.run([sys.executable, "-m", "lunaNMR.mcp.server"],
                               capture_output=True, text=True)
    assert completed.returncode == 2
    assert "--root" in completed.stderr


def test_the_command_line_check_names_a_missing_mcp_sdk_without_a_traceback():
    """An install without the `mcp` extra: `None` in sys.modules makes `import mcp` fail in
    the child as it does when the package is absent."""
    completed = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.modules['mcp'] = None; "
         "from lunaNMR.mcp.server import main; main(['--check'])"],
        capture_output=True, text=True)
    assert completed.returncode == 1, completed.stdout + completed.stderr
    assert "FAIL mcp: cannot import mcp" in completed.stdout
    assert "ok   cli" in completed.stdout
    assert "Traceback" not in completed.stderr
