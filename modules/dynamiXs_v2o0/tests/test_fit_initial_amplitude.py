# ABOUTME: A relaxation fit must not depend on the absolute intensity scale of the spectra.
# ABOUTME: A fixed initial amplitude silently fails datasets recorded at a different gain.

import sys
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_DIR))

from fit_Tx_NMRRE import fit_single_residue

TRUE_T2 = 120.0
DELAYS = np.array([8, 17, 34, 51, 68, 85, 102, 136, 170, 204], float)


def _decay(amplitude):
    return amplitude * np.exp(-DELAYS / TRUE_T2)


class TestAmplitudeScaleInvariance:
    """The same decay recorded at any receiver gain must give the same T2."""

    @pytest.mark.parametrize('amplitude', [1e2, 1e6, 1e10, 3e10])
    def test_derived_initials_recover_t2_at_any_scale(self, amplitude):
        fit = fit_single_residue(DELAYS, _decay(amplitude), 'A1', n_bootstrap=0)
        assert fit is not None, f'no fit at amplitude {amplitude:g}'
        assert fit['t2'] == pytest.approx(TRUE_T2, rel=0.02)

    def test_two_scales_of_the_same_decay_agree(self):
        """Two fields of one real dataset differ ~30,000x in intensity scale; the fitted T2
        must not move because of that."""
        small = fit_single_residue(DELAYS, _decay(1e6), 'A1', n_bootstrap=0)
        large = fit_single_residue(DELAYS, _decay(3e10), 'A1', n_bootstrap=0)
        assert abs(small['t2'] - large['t2']) < 0.01


class TestPipelineDefaults:
    """The integrated pipeline must not override the data-driven initials with a constant."""

    def test_initial_amplitudes_default_to_data_driven(self):
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        from dynamiXs_integrated.integrated_analysis import IntegratedAnalysisParameters
        p = IntegratedAnalysisParameters()
        assert p.t1_initial_amplitude is None, 'a fixed amplitude breaks other intensity scales'
        assert p.t2_initial_amplitude is None
