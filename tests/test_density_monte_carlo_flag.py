# ABOUTME: Dual-field spectral density runs Monte-Carlo error refits only when asked, as single-field does.
# ABOUTME: Real JA1 600/700 MHz R1/R2/NOE for four residues, sequential and multiprocessing paths.

"""
analyze_dual_field_csv took use_monte_carlo_errors and never read it: every dual-field
`dynamixs density` run, and the GUI's dual spectral density page with Monte-Carlo off,
refitted each residue n_monte_carlo times (1000 from the CLI) with unseeded sampling.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "modules" / "dynamiXs_v2o0"))

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405" / "modelfree")
SOURCES = {'f600': RUN / "f600" / "f600_spectral_density_basic.csv",
           'f700': RUN / "f700_raw" / "f700_raw_spectral_density_basic.csv"}

pytestmark = pytest.mark.skipif(not RUN.exists(), reason="data_test not present")


def _analyzer(factor):
    from dynamiXs_density_functions.ZZ_multi_2fields_density import DualFieldSpectralDensityAnalysis
    return DualFieldSpectralDensityAnalysis(field1_freq=600.0, field2_freq=700.0, rNH=1.015e-10,
                                            csaN=-172.0e-6, omega_h_factor=factor)


@pytest.fixture
def tables(tmp_path):
    paths = []
    for name, source in SOURCES.items():
        frame = pd.read_csv(source)[['Residue', 'R1', 'R1_err', 'R2', 'R2_err',
                                     'hetNOE', 'hetNOE_err']].head(4)
        frame.columns = ['Residue', 'R1', 'R1err', 'R2', 'R2err', 'hetNOE', 'hetNOEerr']
        frame.to_csv(tmp_path / f"{name}.csv", index=False)
        paths.append(str(tmp_path / f"{name}.csv"))
    return paths


@pytest.mark.parametrize("factor", [1.0, 0.87])
@pytest.mark.parametrize("parallel", [False, True])
def test_no_monte_carlo_refits_unless_asked(factor, parallel, tables, capsys):
    analyzer = _analyzer(factor)
    df = analyzer.analyze_dual_field_csv(tables[0], tables[1], use_monte_carlo_errors=False,
                                         n_monte_carlo=12, use_multiprocessing=parallel,
                                         n_processes=2)
    assert len(df) == 4
    assert df['S2_err'].isna().all()
    assert 'Monte Carlo' not in capsys.readouterr().out


@pytest.mark.parametrize("factor", [1.0, 0.87])
def test_monte_carlo_refits_give_errors_when_asked(factor, tables):
    analyzer = _analyzer(factor)
    df = analyzer.analyze_dual_field_csv(tables[0], tables[1], use_monte_carlo_errors=True,
                                         n_monte_carlo=12, use_multiprocessing=False)
    assert np.isfinite(df['S2_err']).all()
