# ABOUTME: Model-free fits reach the minimum of their objective, deterministically, single and dual field,
# ABOUTME: J(wH) at 1.0 or 0.87 wH. Real JA1 R1/R2/NOE at 600/700 MHz from a real modelfree run.

"""
The fit was a Monte-Carlo search whose steps MULTIPLY each parameter, so a member with
Rex = 0 keeps Rex = 0: on 13ValH an exact solution has Rex 0.31 s-1, and four seeds in five
reported Rex 0 with S2 0.687 instead of 0.671. Its only early exit was an exact fit, so
every fit the data could not match exactly ran to the iteration cap, and a dual-field fit
always did. The search is now a bounded least-squares fit from a grid of starts.
"""

import random
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
TRAPPED = "13ValH"

pytestmark = pytest.mark.skipif(not RUN.exists(), reason="data_test not present")


def _basic(tmp_path_factory, field):
    """The per-residue R1/R2/NOE and fixed tc a real single-field modelfree run used."""
    from lunaNMR.cli import main
    out = tmp_path_factory.mktemp(f"mf{field}")
    noe = RUN / "results" / "hetnoe"
    assert main(["dynamixs", "modelfree",
                 "--f1-t1", str(RUN / "series" / f"{field}_T1" / "peak_intensity_matrix.csv"),
                 "--f1-t2", str(RUN / "series" / f"{field}_T2" / "peak_intensity_matrix.csv"),
                 "--f1-noe-sat", str(noe / f"f{field}_noe_sat.csv"),
                 "--f1-noe-unsat", str(noe / f"f{field}_noe_unsat.csv"),
                 "--field1-freq", str(field), "--out", str(out), "--prefix", "mf",
                 "--n-monte-carlo", "2"]) == 0
    return pd.read_csv(out / "mf_spectral_density_basic.csv").set_index("Residue")


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    return {600: _basic(tmp_path_factory, 600), 700: _basic(tmp_path_factory, 700)}


def _single(factor, field):
    from dynamiXs_density_functions.ZZ_multi_density import ReducedSpectralDensityAnalysis
    return ReducedSpectralDensityAnalysis(spectrometer_frequency=float(field), rNH=1.015e-10,
                                          csaN=-172.0e-6, omega_h_factor=factor)


def _dual(factor):
    from dynamiXs_density_functions.ZZ_multi_2fields_density import DualFieldSpectralDensityAnalysis
    return DualFieldSpectralDensityAnalysis(field1_freq=600.0, field2_freq=700.0, rNH=1.015e-10,
                                            csaN=-172.0e-6, omega_h_factor=factor)


def _fit_single(factor, inputs, residue):
    row = inputs[600].loc[residue]
    return _single(factor, 600)._fit_single_dataset(
        row.R1, row.R2, row.hetNOE, row.tc * 1e-9, 10_000)


def _fit_dual(factor, inputs, residue):
    a, b = inputs[600].loc[residue], inputs[700].loc[residue]
    return _dual(factor)._fit_dual_field_dataset(
        a.R1, a.R2, a.hetNOE, b.R1, b.R2, b.hetNOE, a.tc * 1e-9, 10_000)


@pytest.fixture(autouse=True)
def modules_path():
    import sys
    path = str(Path(__file__).resolve().parents[1] / "modules" / "dynamiXs_v2o0")
    if path not in sys.path:
        sys.path.insert(0, path)


def test_a_residue_that_needs_a_small_rex_gets_it(inputs):
    fit = _fit_single(0.87, inputs, TRAPPED)
    assert fit["Rex"] == pytest.approx(0.31, abs=0.02)
    assert fit["S2"] == pytest.approx(0.6716, abs=0.002)
    assert fit["chi2"] < 1e-8
    assert fit["at_iteration_cap"] is False


@pytest.mark.parametrize("factor", [1.0, 0.87])
def test_a_single_field_fit_does_not_depend_on_the_seed(inputs, factor):
    fits = []
    for seed in (0, 1):
        random.seed(seed)
        np.random.seed(seed)
        fits.append(_fit_single(factor, inputs, TRAPPED))
    assert fits[0]["S2"] == fits[1]["S2"]
    assert fits[0]["Rex"] == fits[1]["Rex"]


@pytest.mark.parametrize("factor", [1.0, 0.87])
def test_a_dual_field_fit_converges_and_does_not_depend_on_the_seed(inputs, factor):
    fits = []
    for seed in (0, 1):
        random.seed(seed)
        np.random.seed(seed)
        fits.append(_fit_dual(factor, inputs, TRAPPED))
    assert fits[0]["at_iteration_cap"] is False
    assert np.isfinite(fits[0]["chi2"])
    assert fits[0]["S2"] == fits[1]["S2"]
    assert fits[0]["Rex_f1"] == fits[1]["Rex_f1"]


def test_without_monte_carlo_the_errors_are_not_invented(inputs):
    """The ensemble's spread was ~1e-9: copies of one member, not an uncertainty."""
    fit = _fit_single(0.87, inputs, TRAPPED)
    assert np.isnan(fit["S2_err"])
