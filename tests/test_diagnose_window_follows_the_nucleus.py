# ABOUTME: `diagnose` judges capture and registration with the search window the pipeline uses
# ABOUTME: for the spectrum's nucleus (0.3 ppm in 15N, 0.11 ppm in 13C), read from the PS2D config.

import contextlib
import io
from pathlib import Path

import pytest

from lunaNMR.core.ps2d_config import PS2DConfig
from lunaNMR.validation import spectra_check

HSQC_15N = Path(__file__).resolve().parents[1] / "data_example" / "2DHSQC" / "600_T1_0o0.ft"
METHYL_13C = (Path(__file__).resolve().parents[2] / "data_test" / "methyl" / "sampleNBD"
              / "001_T2_NBD_18ms.ft")


def _integrator_nucleus(path):
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    integrator = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        assert integrator.load_nmr_file(str(path))
    return integrator.auto_detected_nucleus


@pytest.mark.skipif(not HSQC_15N.exists(), reason="data_example/2DHSQC not present")
def test_a_15n_spectrum_is_judged_with_the_15n_window():
    config = PS2DConfig('15N')
    assert _integrator_nucleus(HSQC_15N) == '15N'
    assert spectra_check.search_window(HSQC_15N) == (config.search_window_f2, config.search_window_f1)


@pytest.mark.skipif(not METHYL_13C.exists(), reason="data_test not present")
def test_a_13c_spectrum_is_judged_with_the_13c_window():
    config = PS2DConfig('13C')
    assert _integrator_nucleus(METHYL_13C) == '13C'
    assert spectra_check.search_window(METHYL_13C) == (config.search_window_f2, config.search_window_f1)
    assert spectra_check.search_window(METHYL_13C)[1] < PS2DConfig('15N').search_window_f1


def test_the_module_holds_no_window_of_its_own():
    assert not hasattr(spectra_check, 'SEARCH_WINDOW')
