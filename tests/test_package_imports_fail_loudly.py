# ABOUTME: A missing piece stops the batch package at import, and lunaNMR.utils loads a submodule only
# ABOUTME: when one of its names is used, so importing output_manager never pulls in Qt or the project code.

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

BLOCK = """
import importlib.abc, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        if name == {blocked!r}:
            raise ImportError('blocked: ' + name)
        return None
sys.meta_path.insert(0, Block())
try:
    import {module}
except ImportError as e:
    print('IMPORT_ERROR', e)
else:
    print('IMPORTED')
"""


def _run(code):
    out = subprocess.run([sys.executable, '-c', code], cwd=ROOT, capture_output=True, text=True,
                         timeout=120)
    assert out.returncode == 0, out.stderr
    return out.stdout.strip().splitlines()[-1]


@pytest.mark.parametrize('module, blocked', [
    ('lunaNMR.batch_processing.batch_processor', 'lunaNMR.core.core_integrator'),
    ('lunaNMR.batch_processing.cli_interface', 'lunaNMR.batch_processing.batch_processor'),
    ('lunaNMR.batch_processing', 'lunaNMR.batch_processing.config_manager'),
])
def test_the_batch_package_does_not_import_without(module, blocked):
    assert _run(BLOCK.format(module=module, blocked=blocked)) == f'IMPORT_ERROR blocked: {blocked}'


def test_output_manager_loads_no_other_utils_module():
    out = _run("import sys, lunaNMR.utils.output_manager; "
               "print(sorted(m for m in sys.modules if m.startswith('lunaNMR.utils.')))")
    assert out == "['lunaNMR.utils.output_manager']"


def test_a_utils_name_resolves_from_the_package():
    out = _run("import lunaNMR.utils as u; print(u.ProjectManager.__name__, u.log_info.__name__)")
    assert out == "ProjectManager log_info"


def test_an_unknown_utils_name_is_an_attribute_error():
    assert _run("import lunaNMR.utils as u; print(hasattr(u, 'ConfigManager'))") == "False"


def test_the_peak_picker_holds_no_fitter():
    from lunaNMR.core.enhanced_peak_picker import EnhancedPeakPicker
    assert not hasattr(EnhancedPeakPicker(), 'fitter')
