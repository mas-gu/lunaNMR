# ABOUTME: A dual-field model-free fit asked for fewer than ten Monte-Carlo refits returns no
# ABOUTME: errors (NaN), and the line it prints says that rather than naming errors it does not return.

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# R1, R2 (s-1) and hetNOE of an ordered backbone amide at 600 and 700 MHz, 2% / 2% / 0.02 errors
FIELD1 = dict(r1=1.60, r2=10.5, noe=0.78)
FIELD2 = dict(r1=1.35, r2=11.8, noe=0.81)
TC = 7.0e-9


def test_too_few_refits_says_the_errors_are_not_estimated(capsys):
    from dynamiXs_density_functions.ZZ_multi_2fields_density import DualFieldSpectralDensityAnalysis
    analysis = DualFieldSpectralDensityAnalysis(field1_freq=600.0, field2_freq=700.0)

    fit = analysis.fit_dual_field_model_free(
        FIELD1['r1'], FIELD1['r2'], FIELD1['noe'], FIELD2['r1'], FIELD2['r2'], FIELD2['noe'], TC,
        r1_f1_err=0.02 * FIELD1['r1'], r2_f1_err=0.02 * FIELD1['r2'], noe_f1_err=0.02,
        r1_f2_err=0.02 * FIELD2['r1'], r2_f2_err=0.02 * FIELD2['r2'], noe_f2_err=0.02,
        n_monte_carlo=5)

    assert all(np.isnan(fit[key]) for key in ('S2_err', 'te_err', 'Rex_f1_err', 'Rex_f2_err'))
    printed = capsys.readouterr().out
    assert "5/5 dual-field MC fits succeeded" in printed
    assert "ensemble" not in printed
    assert "errors are not estimated (NaN)" in printed
