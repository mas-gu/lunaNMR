# ABOUTME: Importing the relaxation modules leaves the process environment as it found it.
# ABOUTME: Thread counts are the host's to set; a module setting them is inherited by every worker.

"""
Four density modules and one T1/T2 fitter set OMP_NUM_THREADS and MKL_NUM_THREADS to 80% of
the CPUs when imported, so any host that imported them - an MCP server, a GUI - had its
thread counts changed, and passed the change to every process it spawned.
"""

import subprocess
import sys
from pathlib import Path

import pytest

DYNAMIXS = Path(__file__).resolve().parents[1] / "modules" / "dynamiXs_v2o0"
MODULES = ["dynamiXs_density_functions.ZZ_multi_density",
           "dynamiXs_density_functions.ZZ_multi_2fields_density",
           "fitmulti__Tx_NMRRE"]


@pytest.mark.parametrize("module", MODULES)
def test_importing_it_sets_no_thread_count(module):
    code = (f"import os, sys, importlib\n"
            f"sys.path[:0] = [{str(DYNAMIXS)!r}, {str(DYNAMIXS / 'dynamiXs_T1_T2')!r}]\n"
            "import matplotlib; matplotlib.use('Agg')\n"
            "for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.pop(k, None)\n"
            f"importlib.import_module({module!r})\n"
            "print([k for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS') if k in os.environ])\n")
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-600:]
    assert p.stdout.strip().splitlines()[-1] == "[]"
