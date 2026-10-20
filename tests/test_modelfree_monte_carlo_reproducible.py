# ABOUTME: Model-free Monte-Carlo errors are the same on every run of the same input, sequential or parallel.
# ABOUTME: Real JA1 600/700 MHz R1/R2/NOE for four residues, single- and dual-field analyses.

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
N_SAMPLES = 12

pytestmark = pytest.mark.skipif(not RUN.exists(), reason="data_test not present")


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


def _single(tables, parallel):
    from dynamiXs_density_functions.ZZ_multi_density import ReducedSpectralDensityAnalysis
    analyzer = ReducedSpectralDensityAnalysis(spectrometer_frequency=600.0)
    return analyzer.analyze_csv(tables[0], use_monte_carlo_errors=True, n_monte_carlo=N_SAMPLES,
                                use_multiprocessing=parallel, n_processes=2)


def _dual(tables, parallel):
    from dynamiXs_density_functions.ZZ_multi_2fields_density import DualFieldSpectralDensityAnalysis
    analyzer = DualFieldSpectralDensityAnalysis(field1_freq=600.0, field2_freq=700.0)
    return analyzer.analyze_dual_field_csv(tables[0], tables[1], use_monte_carlo_errors=True,
                                           n_monte_carlo=N_SAMPLES, use_multiprocessing=parallel,
                                           n_processes=2)


def _errors(df):
    columns = [c for c in df.columns if c.endswith('_err') and c.split('_')[0] in ('S2', 'te', 'Rex')]
    assert 'S2_err' in columns
    frame = df.set_index('Residue')[columns]
    assert np.isfinite(frame['S2_err']).all()
    return frame


@pytest.mark.parametrize("analysis", [_single, _dual], ids=["single", "dual"])
def test_two_runs_give_the_same_errors(analysis, tables):
    first = _errors(analysis(tables, parallel=False))
    second = _errors(analysis(tables, parallel=False))
    pd.testing.assert_frame_equal(first, second, check_exact=True)


@pytest.mark.parametrize("analysis", [_single, _dual], ids=["single", "dual"])
def test_parallel_gives_the_sequential_errors(analysis, tables):
    sequential = _errors(analysis(tables, parallel=False))
    parallel = _errors(analysis(tables, parallel=True))
    pd.testing.assert_frame_equal(sequential.sort_index(), parallel.sort_index(), check_exact=True)


@pytest.mark.parametrize("analysis", [_single, _dual], ids=["single", "dual"])
def test_each_residue_draws_its_own_samples(analysis, tables):
    """One seed shared by every residue would give every residue the same standard
    normal draws, so the sampling error of the error bar would be the same factor on all."""
    errors = _errors(analysis(tables, parallel=False))
    from dynamiXs_density_functions.spectral_density_physics import monte_carlo_rng
    draws = [monte_carlo_rng(1.0, 2.0, 0.5).standard_normal(4),
             monte_carlo_rng(1.1, 2.0, 0.5).standard_normal(4)]
    assert not np.array_equal(*draws)
    assert errors['S2_err'].nunique() == len(errors)
