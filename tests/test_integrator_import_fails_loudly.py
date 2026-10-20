# ABOUTME: A broken install stops the integrator at import, rather than running it with a piece missing.
# ABOUTME: Each module it is built on is made unimportable in a fresh interpreter and the import must raise.

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

BLOCK_AND_IMPORT = """
import importlib.abc, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        if name == {blocked!r}:
            raise ImportError('blocked: ' + name)
        return None
sys.meta_path.insert(0, Block())
try:
    import lunaNMR.core.core_integrator
except ImportError as e:
    print('IMPORT_ERROR', e)
else:
    print('IMPORTED')
"""


@pytest.mark.parametrize('blocked', [
    'lunaNMR.integrators.inplace_advanced_nmr_integrator',
    'lunaNMR.core.enhanced_voigt_fitter',
    'lunaNMR.core.ps2d_2d_fitter',
])
def test_the_integrator_does_not_import_without(blocked):
    out = subprocess.run([sys.executable, '-c', BLOCK_AND_IMPORT.format(blocked=blocked)],
                         cwd=ROOT, capture_output=True, text=True, timeout=120)
    assert out.stdout.strip().splitlines()[-1] == f'IMPORT_ERROR blocked: {blocked}'
