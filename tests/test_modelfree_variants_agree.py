# ABOUTME: The single- and dual-field model-free variants predict the same relaxation at one field.
# ABOUTME: Same NOE (from the calculated T1), same J frequencies and same J errors, plain and 0.87 pairs.

"""
Two inconsistencies made the same residue fit differently depending on the variant: the
single-field NOE was predicted from the MEASURED T1, the dual-field one from the T1 the
model predicts; and the dual-field 0.87 variant evaluated J at 0.87wH +/- wN where the
single-field one, and the definition of the variant, use the true wH +/- wN and replace
only J(wH). A prediction is now a function of the parameters alone.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "modules" / "dynamiXs_v2o0"))

# A JA1 600 MHz residue: measured R1, R2, NOE and the run's fixed tc, plus a trial point.
R1, R2, NOE, TC = 1.62, 6.31, 0.71, 6.61e-9
TRIAL = (0.72, 3e-11, 0.4)


def _pair(factor):
    from dynamiXs_density_functions.ZZ_multi_density import ReducedSpectralDensityAnalysis
    from dynamiXs_density_functions.ZZ_multi_2fields_density import DualFieldSpectralDensityAnalysis
    return (ReducedSpectralDensityAnalysis(spectrometer_frequency=600.0, rNH=1.015e-10,
                                           csaN=-172.0e-6, omega_h_factor=factor),
            DualFieldSpectralDensityAnalysis(field1_freq=600.0, field2_freq=700.0, rNH=1.015e-10,
                                             csaN=-172.0e-6, omega_h_factor=factor))


@pytest.mark.parametrize("factor", [1.0, 0.87])
def test_one_field_of_a_dual_fit_predicts_what_the_single_fit_does(factor):
    single, dual = _pair(factor)
    s2, te, rex = TRIAL
    a = single._single_field_residuals(s2, te, rex, 1 / R1, 1 / R2, NOE, TC)
    b = dual._field_residuals(1, s2, te, rex, 1 / R1, 1 / R2, NOE, TC)
    assert a == pytest.approx(b, rel=1e-12)


@pytest.mark.parametrize("factor", [1.0, 0.87])
def test_the_predicted_noe_does_not_depend_on_the_measured_t1(factor):
    single, _ = _pair(factor)
    s2, te, rex = TRIAL
    noe_dev = [single._single_field_residuals(s2, te, rex, 1 / r1, 1 / R2, NOE, TC)[2]
               for r1 in (R1, 1.5 * R1)]
    assert noe_dev[0] == pytest.approx(noe_dev[1], rel=1e-12)


ERRORS = dict(r1_err=0.05, r2_err=0.3, noe_err=0.03)


@pytest.mark.parametrize("factor", [1.0, 0.87])
def test_both_variants_propagate_the_same_spectral_density_errors(factor):
    """First-order propagation, sqrt(sum((dJ/dx * sigma_x)^2)), in both: the single-field
    J(wN) error was once halved, so one residue reported two different errors."""
    single, dual = _pair(factor)
    e = ERRORS
    a = single.propagate_errors(R1, e['r1_err'], R2, e['r2_err'], NOE, e['noe_err'])
    b = dual.propagate_dual_field_errors(R1, e['r1_err'], R2, e['r2_err'], NOE, e['noe_err'],
                                         R1, e['r1_err'], R2, e['r2_err'], NOE, e['noe_err'])['Field1']
    for key in ('J0_err', 'JwN_err'):
        assert a[key] == pytest.approx(b[key], rel=1e-12, abs=0), key
    jwh_a, = [k for k in a if k.startswith('JwH') and k.endswith('_err')]
    jwh_b, = [k for k in b if k.startswith('JwH') and k.endswith('_err')]
    assert a[jwh_a] == pytest.approx(b[jwh_b], rel=1e-12, abs=0)
