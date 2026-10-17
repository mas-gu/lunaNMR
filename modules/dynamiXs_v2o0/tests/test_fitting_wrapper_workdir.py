# ABOUTME: A relaxation fit given no output prefix writes to a folder of its own, per call.
# ABOUTME: A pipeline writes into its output folder, whatever the working directory is.

import sys
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))
sys.path.insert(0, str(_DIR / "dynamiXs_T1_T2"))

MATRIX = (Path(__file__).resolve().parents[4] / "data_test" / "data_JA1_relaxation"
          / "run_20260825_1405" / "series" / "600_T1" / "peak_intensity_matrix.csv")


@pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")
def test_two_unprefixed_fits_do_not_share_files():
    from dynamiXs_integrated.fitting_wrapper import fit_relaxation_data
    _, first = fit_relaxation_data(str(MATRIX), "T1", field_freq=600.0)
    _, second = fit_relaxation_data(str(MATRIX), "T1", field_freq=600.0)
    assert Path(first).parent != Path(second).parent
    assert not Path(first).resolve().is_relative_to(_DIR)


RUN = MATRIX.parents[2]


@pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")
def test_a_pipeline_writes_into_its_output_folder_not_the_working_directory(tmp_path, monkeypatch):
    """The GUI runs the pipeline in place, without the CLI's change into the output
    folder: every file it writes must still land in that folder."""
    import os
    from dynamiXs_integrated.integrated_analysis import (IntegratedAnalysisPipeline,
                                                         IntegratedAnalysisParameters)
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    noe = RUN / "results" / "hetnoe"
    p = IntegratedAnalysisParameters()
    p.field1_t1_file = str(RUN / "series" / "600_T1" / "peak_intensity_matrix.csv")
    p.field1_t2_file = str(RUN / "series" / "600_T2" / "peak_intensity_matrix.csv")
    p.field1_noe_sat_file = str(noe / "f600_noe_sat.csv")
    p.field1_noe_unsat_file = str(noe / "f600_noe_unsat.csv")
    p.field1_freq_mhz = 600.0
    p.enable_dual_field = False
    p.analysis_method = "single_087"
    p.output_dir = str(tmp_path / "out")
    p.output_prefix = "run"
    IntegratedAnalysisPipeline(p).run_complete_analysis()
    assert os.listdir(cwd) == []
    assert (tmp_path / "out" / "run_field1_T1_fit_results.txt").exists()
