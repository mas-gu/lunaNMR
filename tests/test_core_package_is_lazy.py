# ABOUTME: lunaNMR.core loads its engines only when one is named, so importing a light submodule
# ABOUTME: (ps2d_config) cannot drag in core_integrator and close an import cycle through lunaNMR.utils.

import subprocess
import sys


def _run(code):
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert p.returncode == 0, p.stderr
    return p.stdout.strip()


def test_a_light_submodule_does_not_load_the_integrator():
    out = _run("import sys, lunaNMR.core.ps2d_config; "
               "print('lunaNMR.core.core_integrator' in sys.modules)")
    assert out == "False"


def test_the_engines_resolve_from_the_package():
    out = _run("import lunaNMR.core as c; "
               "print(c.EnhancedVoigtIntegrator.__name__, c.EnhancedVoigtFitter.__name__, "
               "c.EnhancedPeakPicker.__name__)")
    assert out == "EnhancedVoigtIntegrator EnhancedVoigtFitter EnhancedPeakPicker"


def test_every_exported_name_exists():
    out = _run("import lunaNMR.core as c; "
               "print([n for n in c.__all__ if not hasattr(c, n)])")
    assert out == "[]"


def test_an_unknown_name_is_an_attribute_error():
    out = _run("import lunaNMR.core as c; print(hasattr(c, 'CoreIntegrator'))")
    assert out == "False"
