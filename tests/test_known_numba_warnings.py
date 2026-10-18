# ABOUTME: The numba warnings lunaNMR's own JIT kernels are known to raise are silenced where they arise.
# ABOUTME: Each check runs in a fresh process, since numba warns only on a function's first compile.

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
# The density modules are imported by package name, as every caller imports them.
DYNAMIXS_DIR = REPO_ROOT / "modules" / "dynamiXs_v2o0"


def _recorded(body):
    """The numba warnings `body` raised in a fresh process, and whether it left the
    process's warning filters as it found them."""
    code = ("import sys, warnings\n"
            f"sys.path.insert(0, {str(DYNAMIXS_DIR)!r})\n"
            "import matplotlib; matplotlib.use('Agg')\n"
            "import numpy as np\n"
            # What these libraries add to the filters on import is theirs, not lunaNMR's.
            "import numba, pandas, scipy.optimize, scipy.special, scipy.stats\n"
            "import matplotlib.pyplot\n"
            "from numba.core.errors import NumbaWarning\n"
            "with warnings.catch_warnings(record=True) as caught:\n"
            "    warnings.simplefilter('always')\n"
            "    before = list(warnings.filters)\n"
            + "".join(f"    {line}\n" for line in body.splitlines()) +
            "    unchanged = list(warnings.filters) == before\n"
            "print(repr([str(w.message).strip()[:80] for w in caught\n"
            "            if issubclass(w.category, NumbaWarning)]))\n"
            "print(unchanged)\n")
    p = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT),
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-800:]
    numba_warnings, unchanged = p.stdout.strip().splitlines()[-2:]
    return numba_warnings, unchanged == "True"


def test_the_2d_voigt_kernel_compiles_quietly():
    """multi_voigt_profile_2d compiles in object mode and cannot be cached; numba says
    both on its first call in every process, including every spawn worker."""
    numba_warnings, unchanged = _recorded(
        "from lunaNMR.core.ps2d_style_fitter import multi_voigt_profile_2d\n"
        "f1, f2 = np.meshgrid(np.linspace(118, 122, 20), np.linspace(7.8, 8.2, 30))\n"
        "params = np.array([120.0, 0.1, 0.1, 8.0, 0.01, 0.01, 1.0, 0.0])\n"
        "surface = multi_voigt_profile_2d(f1, f2, params, 1)\n"
        "assert surface.shape == f1.shape and surface.max() > 0\n")
    assert numba_warnings == "[]"
    assert unchanged


def test_the_density_kernels_raise_no_deprecation():
    numba_warnings, unchanged = _recorded(
        "from dynamiXs_density_functions import ZZ_multi_density, ZZ_multi_density_087\n"
        "from dynamiXs_density_functions import (ZZ_multi_2fields_density,\n"
        "                                        ZZ_multi_2fields_density_087)\n")
    assert numba_warnings == "[]"
    assert unchanged
