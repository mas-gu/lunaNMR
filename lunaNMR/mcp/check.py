# ABOUTME: Install self-test for the MCP server: required modules import and the CLI answers with JSON.
# ABOUTME: Failed installs are the likeliest support request, so each problem is named in one line.

import importlib
import json
import subprocess
import sys
import tempfile

REQUIRED_MODULES = ['numpy', 'scipy', 'pandas', 'matplotlib', 'nmrglue', 'numba', 'lmfit', 'mcp']


def _module_check(name):
    try:
        module = importlib.import_module(name)
    except Exception as error:
        return name, False, f'cannot import {name}: {error}'
    return name, True, getattr(module, '__version__', 'imported')


def _cli_check():
    """Run a real command on an empty folder: it must fail cleanly, with a JSON envelope."""
    with tempfile.TemporaryDirectory() as empty:
        completed = subprocess.run(
            [sys.executable, '-m', 'lunaNMR', 'diagnose', empty, '--format', 'json'],
            capture_output=True, text=True)
    try:
        json.loads(completed.stdout)
    except ValueError:
        tail = ' '.join(completed.stderr.strip().splitlines()[-3:])
        return 'cli', False, f'`python -m lunaNMR` printed no JSON: {tail}'
    return 'cli', True, sys.executable


def run_checks(modules=REQUIRED_MODULES):
    return [_module_check(name) for name in modules] + [_cli_check()]


def report(results):
    """One line per check; the exit status the command should return."""
    for name, ok, detail in results:
        print(f"{'ok  ' if ok else 'FAIL'} {name}: {detail}")
    return 0 if all(ok for _, ok, _ in results) else 1
