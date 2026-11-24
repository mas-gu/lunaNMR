# ABOUTME: The DynamiXs pages, their worker parameters and the CLI start a density or a model-free
# ABOUTME: analysis with the same 15N CSA, J(wH) variant and Monte-Carlo sample count, declared once.

import json
import os
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
BASIC = RUN / "modelfree" / "f600" / "f600_spectral_density_basic.csv"

pytestmark = pytest.mark.skipif(not BASIC.exists(), reason="data_test not present")

CSA_PPM = -172.0
DENSITY_SAMPLES = 1000
MODELFREE_SAMPLES = 50


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def pages(app):
    """The real Spectral Density and Model-free pages, as a fresh dialog shows them."""
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    dialog = DynamiXsDialog(parent=None, main_window=types.SimpleNamespace())
    yield dialog.spectral_page, dialog.integrated_page
    dialog.deleteLater()


@pytest.fixture(scope="module")
def workers():
    import lunaNMR.gui.dialogs.dynamixs_dialog  # noqa: F401 - puts the module on sys.path
    from workers import IntegratedAnalysisParams, SpectralDensityParams
    return SpectralDensityParams(input_file="table.csv"), IntegratedAnalysisParams()


def _record(folder, prefix):
    return json.loads((folder / f"{prefix}_dynamixs.json").read_text())


@pytest.fixture(scope="module")
def density_run(tmp_path_factory):
    """`dynamixs density` with none of the three settings given: its record and results."""
    from lunaNMR.cli import main
    folder = tmp_path_factory.mktemp("density")
    table = pd.read_csv(BASIC)[['Residue', 'R1', 'R1_err', 'R2', 'R2_err',
                                'hetNOE', 'hetNOE_err']].head(6)
    table.columns = ['Residue', 'R1', 'R1err', 'R2', 'R2err', 'hetNOE', 'hetNOEerr']
    table.to_csv(folder / "table.csv", index=False)
    assert main(["dynamixs", "density", "--input", str(folder / "table.csv"),
                 "--out", str(folder / "out"), "--no-plot", "--no-parallel"]) == 0
    record = _record(folder / "out", "spectral_density")
    return record['args'], pd.read_csv(folder / "out" / "spectral_density_results.csv")


@pytest.fixture(scope="module")
def modelfree_run(tmp_path_factory):
    """`dynamixs modelfree` with none of the three settings given: its record."""
    from lunaNMR.cli import main
    folder = tmp_path_factory.mktemp("modelfree")
    noe = RUN / "results" / "hetnoe"
    assert main(["dynamixs", "modelfree",
                 "--f1-t1", str(RUN / "series" / "600_T1" / "peak_intensity_matrix.csv"),
                 "--f1-t2", str(RUN / "series" / "600_T2" / "peak_intensity_matrix.csv"),
                 "--f1-noe-sat", str(noe / "f600_noe_sat.csv"),
                 "--f1-noe-unsat", str(noe / "f600_noe_unsat.csv"),
                 "--field1-freq", "600", "--out", str(folder), "--prefix", "mf"]) == 0
    record = _record(folder, "mf")
    return record['args'], pd.read_csv(folder / "mf_spectral_density_basic.csv")


def _checked_method(integrated_page):
    (method,) = [value for radio, value in integrated_page.method_radios if radio.isChecked()]
    return method


def test_the_csa_is_minus_172_ppm_everywhere(pages, workers, density_run, modelfree_run):
    spectral_page, integrated_page = pages
    assert spectral_page.csa_spin.value() == CSA_PPM
    assert integrated_page.csa_spin.value() == CSA_PPM
    assert workers[0].csa == CSA_PPM
    assert workers[1].csa == CSA_PPM
    assert density_run[0]['csa'] == CSA_PPM
    assert modelfree_run[0]['csa'] == CSA_PPM


def test_the_pipeline_and_the_analyses_start_from_the_same_csa(workers):
    from dynamiXs_integrated.integrated_analysis import IntegratedAnalysisParameters
    from dynamiXs_density_functions.ZZ_multi_density import ReducedSpectralDensityAnalysis
    from dynamiXs_density_functions.ZZ_multi_2fields_density import DualFieldSpectralDensityAnalysis
    assert IntegratedAnalysisParameters().csaN_ppm == CSA_PPM
    assert ReducedSpectralDensityAnalysis().csaN == CSA_PPM * 1e-6
    assert DualFieldSpectralDensityAnalysis().csaN == CSA_PPM * 1e-6


def test_j_is_evaluated_at_087_omega_h_everywhere(pages, workers, density_run, modelfree_run):
    spectral_page, integrated_page = pages
    assert spectral_page.j087_radio.isChecked() and not spectral_page.jwh_radio.isChecked()
    assert _checked_method(integrated_page).endswith("_087")
    assert workers[0].use_087 is True
    assert workers[1].analysis_method.endswith("_087")
    assert density_run[0]['use_087'] is True
    assert 'JwH_087' in density_run[1].columns and 'JwH' not in density_run[1].columns
    assert modelfree_run[0]['method'] == "single_087"
    assert 'JwH_087' in modelfree_run[1].columns


def test_the_analyses_evaluate_j_at_087_omega_h_unless_told(workers):
    from dynamiXs_density_functions.ZZ_multi_density import ReducedSpectralDensityAnalysis
    from dynamiXs_density_functions.ZZ_multi_2fields_density import DualFieldSpectralDensityAnalysis
    assert ReducedSpectralDensityAnalysis().omega_h_factor == 0.87
    assert DualFieldSpectralDensityAnalysis().omega_h_factor == 0.87


def test_a_density_analysis_draws_1000_samples_everywhere(pages, workers, density_run):
    assert pages[0].mc_spin.value() == DENSITY_SAMPLES
    assert workers[0].n_samples == DENSITY_SAMPLES
    assert density_run[0]['n_samples'] == DENSITY_SAMPLES


def test_a_model_free_analysis_draws_50_samples_everywhere(pages, workers, modelfree_run):
    from dynamiXs_integrated.integrated_analysis import IntegratedAnalysisParameters
    assert pages[1].mc_spin.value() == MODELFREE_SAMPLES
    assert workers[1].n_monte_carlo == MODELFREE_SAMPLES
    assert IntegratedAnalysisParameters().monte_carlo_iterations == MODELFREE_SAMPLES
    assert modelfree_run[0]['n_monte_carlo'] == MODELFREE_SAMPLES


def test_each_value_is_declared_once(workers):
    from dynamiXs_density_functions import analysis_defaults
    assert analysis_defaults.CSA_N_PPM == CSA_PPM
    assert analysis_defaults.OMEGA_H_FACTOR == 0.87
    assert analysis_defaults.DENSITY_MONTE_CARLO_SAMPLES == DENSITY_SAMPLES
    assert analysis_defaults.MODELFREE_MONTE_CARLO_SAMPLES == MODELFREE_SAMPLES


def _help(*command):
    from lunaNMR.cli import build_parser
    parser = build_parser()
    for name in command:
        parser = parser._subparsers._group_actions[0].choices[name]
    return " ".join(parser.format_help().split())


def test_the_help_states_the_defaults_the_commands_use(workers):
    density, modelfree = _help("dynamixs", "density"), _help("dynamixs", "modelfree")
    for text in (density, modelfree):
        assert "15N CSA in ppm (default: -172)" in text
    assert "(default: 0.87)" in density
    assert "Monte-Carlo samples (default: 1000)" in density
    assert "(default: 087)" in modelfree
    assert "Monte-Carlo samples (default: 50)" in modelfree


def test_a_page_restores_what_a_run_recorded(pages):
    spectral_page, integrated_page = pages
    spectral_page.restore_session_state({'csa_ppm': -165.0, 'use_087': False, 'n_samples': 30})
    assert spectral_page.csa_spin.value() == -165.0
    assert spectral_page.jwh_radio.isChecked()
    assert spectral_page.mc_spin.value() == 30
    integrated_page.restore_session_state({'csa_ppm': -165.0, 'analysis_method': 'single_jwh',
                                           'monte_carlo_iterations': 30})
    assert integrated_page.csa_spin.value() == -165.0
    assert _checked_method(integrated_page) == 'single_jwh'
    assert integrated_page.mc_spin.value() == 30


def test_a_page_falls_back_to_the_defaults_for_what_a_run_did_not_record(pages):
    spectral_page, integrated_page = pages
    for page, state in ((spectral_page, {'csa_ppm': -165.0, 'use_087': False, 'n_samples': 30}),
                        (integrated_page, {'csa_ppm': -165.0, 'analysis_method': 'single_jwh',
                                           'monte_carlo_iterations': 30})):
        page.restore_session_state(state)
        page.restore_session_state({'prefix': 'a run that recorded nothing else'})
    assert spectral_page.csa_spin.value() == CSA_PPM
    assert spectral_page.j087_radio.isChecked()
    assert spectral_page.mc_spin.value() == DENSITY_SAMPLES
    assert integrated_page.csa_spin.value() == CSA_PPM
    assert _checked_method(integrated_page).endswith("_087")
    assert integrated_page.mc_spin.value() == MODELFREE_SAMPLES
