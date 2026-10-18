# ABOUTME: The single- and dual-field model-free variants predict the same relaxation at one field.
# ABOUTME: Same NOE (from the calculated T1) and same J frequencies, for the plain and the 0.87 pair.

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


def _pair(single_name, dual_name):
    import importlib
    single = importlib.import_module(f"dynamiXs_density_functions.{single_name}")
    dual = importlib.import_module(f"dynamiXs_density_functions.{dual_name}")
    return (single.ReducedSpectralDensityAnalysis(spectrometer_frequency=600.0,
                                                  rNH=1.015e-10, csaN=-172.0e-6),
            dual.DualFieldSpectralDensityAnalysis(field1_freq=600.0, field2_freq=700.0,
                                                  rNH=1.015e-10, csaN=-172.0e-6))


@pytest.mark.parametrize("single_name, dual_name", [
    ("ZZ_multi_density", "ZZ_multi_2fields_density"),
    ("ZZ_multi_density_087", "ZZ_multi_2fields_density_087")])
def test_one_field_of_a_dual_fit_predicts_what_the_single_fit_does(single_name, dual_name):
    single, dual = _pair(single_name, dual_name)
    s2, te, rex = TRIAL
    a = single._single_field_residuals(s2, te, rex, 1 / R1, 1 / R2, NOE, TC)
    b = dual._field_residuals(1, s2, te, rex, 1 / R1, 1 / R2, NOE, TC)
    assert a == pytest.approx(b, rel=1e-12)


@pytest.mark.parametrize("single_name", ["ZZ_multi_density", "ZZ_multi_density_087"])
def test_the_predicted_noe_does_not_depend_on_the_measured_t1(single_name):
    single, _ = _pair(single_name, "ZZ_multi_2fields_density")
    s2, te, rex = TRIAL
    noe_dev = [single._single_field_residuals(s2, te, rex, 1 / r1, 1 / R2, NOE, TC)[2]
               for r1 in (R1, 1.5 * R1)]
    assert noe_dev[0] == pytest.approx(noe_dev[1], rel=1e-12)
