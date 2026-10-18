# ABOUTME: The relaxation fitters import as dynamiXs_T1_T2.<module> with only the module root on
# ABOUTME: sys.path, as the DynamiXs dialog and its workers import them.

import subprocess
import sys
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parent.parent


@pytest.mark.parametrize("module", ["fit_Tx_NMRRE", "fitmulti__Tx_NMRRE", "fit_methyl_T2"])
def test_a_fitter_imports_with_only_the_module_root_on_the_path(module):
    code = ("import sys; "
            f"sys.path[:] = [p for p in sys.path if 'dynamiXs' not in p] + [{str(_DIR)!r}]; "
            f"from dynamiXs_T1_T2 import {module}")
    run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         cwd=str(_DIR.parent))
    assert run.returncode == 0, run.stderr
