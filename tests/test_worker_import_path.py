# ABOUTME: The import path a parallel worker is handed is the folder that HOLDS the lunaNMR
# ABOUTME: package, so no subpackage of it (mcp, utils, core, gui) becomes a top-level module.

import contextlib
import io
import os
from pathlib import Path

import pytest

SPECTRA = Path(__file__).resolve().parents[1] / "data_example" / "2DHSQC"
REFERENCE_POINT = "600_T1_0o0.ft"

pytestmark = pytest.mark.skipif(not (SPECTRA / REFERENCE_POINT).exists(),
                                reason="data_example/2DHSQC not present")


def test_the_worker_path_holds_the_package_and_is_not_the_package():
    import lunaNMR
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.parallel_voigt_processor import ParallelVoigtProcessor
    from lunaNMR.core.ps2d_config import set_ps2d_config
    set_ps2d_config('15N')
    integrator = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        integrator.load_nmr_file(str(SPECTRA / REFERENCE_POINT))
    fitter = integrator.enhanced_fitter
    fitter.nmr_data, fitter.ppm_x_axis, fitter.ppm_y_axis = (
        integrator.nmr_data, integrator.ppm_x_axis, integrator.ppm_y_axis)
    integrator.gui_params = {'use_ps2d_multi_peak': True, 'use_adaptive_optimization': False}
    processor = ParallelVoigtProcessor(fitter, enable_adaptive=False)
    try:
        shared = processor._prepare_shared_context([])
    finally:
        processor._cleanup_shared_memory()
    package = os.path.dirname(os.path.abspath(lunaNMR.__file__))
    assert os.path.samefile(shared['lunaNMR_path'], os.path.dirname(package))
    assert not os.path.isdir(os.path.join(shared['lunaNMR_path'], 'mcp'))
