# ABOUTME: The relaxation and density pipelines put their module folder on sys.path once.
# ABOUTME: A long-lived process running model-free repeatedly otherwise grows sys.path per call.

import sys
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))
sys.path.insert(0, str(_DIR / "dynamiXs_T1_T2"))

MATRIX = (Path(__file__).resolve().parents[4] / "data_test" / "data_JA1_relaxation"
          / "run_20260825_1405" / "series" / "600_T1" / "peak_intensity_matrix.csv")


@pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")
def test_a_relaxation_fit(tmp_path):
    from dynamiXs_integrated.fitting_wrapper import fit_relaxation_data
    folder = str(_DIR / "dynamiXs_T1_T2")
    fit_relaxation_data(str(MATRIX), "T1", field_freq=600.0,
                        output_prefix=str(tmp_path / "a"))
    count = sys.path.count(folder)
    fit_relaxation_data(str(MATRIX), "T1", field_freq=600.0,
                        output_prefix=str(tmp_path / "b"))
    assert sys.path.count(folder) == count


def test_a_density_analysis(tmp_path):
    from dynamiXs_integrated.integrated_analysis import (IntegratedAnalysisPipeline,
                                                         IntegratedAnalysisParameters)
    folder = str(_DIR / "dynamiXs_density_functions")
    p = IntegratedAnalysisParameters()
    p.output_dir = str(tmp_path)
    p.analysis_method = "unknown"
    pipeline = IntegratedAnalysisPipeline(p)
    counts = []
    for _ in range(2):
        with pytest.raises(ValueError, match="Unknown analysis method"):
            pipeline._run_spectral_density_analysis({})
        counts.append(sys.path.count(folder))
    assert counts[0] == counts[1] >= 1
