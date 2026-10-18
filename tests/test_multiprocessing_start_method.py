# ABOUTME: Importing lunaNMR's parallel modules leaves the process's multiprocessing start method unset.
# ABOUTME: Each pool asks for a spawn context itself; a parallel series run still completes.

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
MODULE_DIRS = [REPO_ROOT / "modules" / "dynamiXs_v2o0" / "dynamiXs_T1_T2",
               REPO_ROOT / "modules" / "dynamiXs_v2o0"]

_MODULES = [
    "lunaNMR.core.parallel_voigt_processor",
    "lunaNMR.core.adaptive_optimizer",
    "lunaNMR.processors.parallel_fitting",
    "fitmulti__Tx_NMRRE",
    # by package name, as every caller imports them
    "dynamiXs_density_functions.ZZ_multi_density",
    "dynamiXs_density_functions.ZZ_multi_density_087",
    "dynamiXs_density_functions.ZZ_multi_2fields_density",
    "dynamiXs_density_functions.ZZ_multi_2fields_density_087",
]


@pytest.mark.parametrize("module", _MODULES)
def test_importing_leaves_the_start_method_unset(module):
    """A host that imports the series engine kept whatever start method the import
    chose for it; on Linux that turned the host's own fork pools into spawn ones."""
    code = ("import multiprocessing as mp, importlib, sys\n"
            f"sys.path[:0] = {[str(d) for d in MODULE_DIRS]!r}\n"
            "import matplotlib; matplotlib.use('Agg')\n"
            f"importlib.import_module({module!r})\n"
            "print(mp.get_start_method(allow_none=True))\n")
    p = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT),
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-800:]
    assert p.stdout.strip().splitlines()[-1] == "None"


@pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")
def test_a_parallel_series_run_completes(tmp_path):
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    for name in ("2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o125.ft", "2D_DNAJA1_HSPA8_0o25.ft"):
        os.symlink(DNAJA1 / name, spectra / name)
    p = subprocess.run([sys.executable, "-m", "lunaNMR", "series", "--spectra", str(spectra),
                        "--peaks", str(DNAJA1 / "DNAJA1.txt"), "--out", str(tmp_path / "out"),
                        "--mode", "titration", "--parallel", "--format", "json"],
                       cwd=str(REPO_ROOT), capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-1500:]
    payload = json.loads(p.stdout)
    assert payload["ok"] is True
    assert payload["spectra_fitted"] == 3
    assert "Parallel cluster execution summary" in p.stderr, p.stderr[-1500:]
    assert "Parallel cluster execution failed" not in p.stderr
    assert "Falling back" not in p.stderr, p.stderr[-1500:]
