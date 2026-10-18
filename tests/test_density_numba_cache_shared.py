# ABOUTME: `dynamixs density` must still run after `dynamixs modelfree` has written the numba cache.
# ABOUTME: The cache records the module name it was compiled under, so both must import one name.

"""
The spectral-density kernels are `@jit(cache=True, forceobj=True)`. An object-mode cache
entry is rebuilt by importing the module under the name it was compiled in. The
model-free pipeline imported the density module by its bare name and `density` by its
package name, so once one had run, the other died with
`ModuleNotFoundError: No module named 'ZZ_multi_density_087'` - in every later process,
since the cache lives on disk.
"""

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
NOE = RUN / "results" / "hetnoe"
BASIC = RUN / "modelfree" / "f600" / "f600_spectral_density_basic.csv"

pytestmark = pytest.mark.skipif(not RUN.exists(), reason="data_test not present")


def _cli(args, cache_dir):
    env = dict(os.environ, NUMBA_CACHE_DIR=str(cache_dir), QT_QPA_PLATFORM="offscreen")
    return subprocess.run([sys.executable, '-m', 'lunaNMR'] + args, cwd=str(REPO_ROOT),
                          env=env, capture_output=True, text=True)


def test_density_runs_on_the_cache_modelfree_wrote(tmp_path):
    cache = tmp_path / "numba_cache"
    series = RUN / "series"
    mf = _cli(["dynamixs", "modelfree",
               "--f1-t1", str(series / "600_T1" / "peak_intensity_matrix.csv"),
               "--f1-t2", str(series / "600_T2" / "peak_intensity_matrix.csv"),
               "--f1-noe-sat", str(NOE / "f600_noe_sat.csv"),
               "--f1-noe-unsat", str(NOE / "f600_noe_unsat.csv"),
               "--field1-freq", "600", "--out", str(tmp_path / "mf"), "--prefix", "mf",
               "--n-monte-carlo", "2"], cache)
    assert mf.returncode == 0, mf.stderr[-2000:]
    assert any(cache.rglob("*.nbi")), "model-free wrote no numba cache"

    frame = pd.read_csv(BASIC)[['Residue', 'R1', 'R1_err', 'R2', 'R2_err',
                                'hetNOE', 'hetNOE_err']]
    frame.columns = ['Residue', 'R1', 'R1err', 'R2', 'R2err', 'hetNOE', 'hetNOEerr']
    frame.to_csv(tmp_path / "table.csv", index=False)
    density = _cli(["dynamixs", "density", "--input", str(tmp_path / "table.csv"),
                    "--field1-freq", "600", "--no-parallel", "--no-plot",
                    "--out", str(tmp_path / "dens")], cache)
    assert density.returncode == 0, density.stderr[-2000:]
