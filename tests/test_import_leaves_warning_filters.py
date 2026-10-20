# ABOUTME: Importing lunaNMR's engine and the density modules leaves warnings.filters unchanged.
# ABOUTME: A module-level filterwarnings('ignore') silenced every warning in the importing process.

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DENSITY_DIR = REPO_ROOT / "modules" / "dynamiXs_v2o0" / "dynamiXs_density_functions"

_MODULES = [
    "lunaNMR.core.core_integrator",
    "lunaNMR.core.enhanced_peak_picker",
    "lunaNMR.core.enhanced_voigt_fitter",
    "lunaNMR.utils.file_manager",
    "lunaNMR.integrators.inplace_advanced_nmr_integrator",
    # by package name, as every caller imports them
    "dynamiXs_density_functions.ZZ_multi_density",
    "dynamiXs_density_functions.ZZ_multi_2fields_density",
]

# The third-party libraries those modules import, loaded before the snapshot: what they
# do to the filters on import is theirs, not lunaNMR's.
_THIRD_PARTY = ("numpy", "scipy.optimize", "scipy.signal", "scipy.special", "pandas",
                "matplotlib.pyplot", "lmfit", "nmrglue", "numba")


@pytest.mark.parametrize("module", _MODULES)
def test_importing_leaves_the_filters_unchanged(module):
    code = ("import importlib, sys, warnings\n"
            f"sys.path[:0] = [{str(DENSITY_DIR)!r}, {str(DENSITY_DIR.parent)!r}]\n"
            "import matplotlib; matplotlib.use('Agg')\n"
            f"for name in {_THIRD_PARTY!r}:\n"
            "    try:\n"
            "        importlib.import_module(name)\n"
            "    except ImportError:\n"
            "        pass\n"
            "before = list(warnings.filters)\n"
            f"importlib.import_module({module!r})\n"
            "added = [f for f in warnings.filters if f not in before]\n"
            "print(repr(added))\n")
    p = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT),
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-800:]
    assert p.stdout.strip().splitlines()[-1] == "[]"
