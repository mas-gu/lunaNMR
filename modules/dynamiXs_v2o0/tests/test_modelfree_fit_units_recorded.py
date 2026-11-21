# ABOUTME: The model-free pipeline's T1/T2 fit JSONs record the delay units stated for each
# ABOUTME: series, which the Fit Viewer and the fit figures label their time axis with.

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))
sys.path.insert(0, str(_DIR / "dynamiXs_T1_T2"))

RESIDUES = ("3LysH", "8TyrH", "10ValH", "11LeuH")


def _matrix(path, delays, time_constant):
    rows = [{"Assignment": name,
             **{f"{d:g}": 1.0e6 * n * np.exp(-d / (time_constant * (1 + 0.02 * n)))
                for d in delays}}
            for n, name in enumerate(RESIDUES, start=1)]
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def _plane(path, scale):
    pd.DataFrame([[name, scale * 1.0e6] for name in RESIDUES]).to_csv(
        path, index=False, header=False)
    return str(path)


def test_each_fit_json_records_the_units_of_its_own_series(tmp_path):
    from dynamiXs_integrated.integrated_analysis import (IntegratedAnalysisPipeline,
                                                         IntegratedAnalysisParameters)
    params = IntegratedAnalysisParameters()
    params.output_dir = str(tmp_path / "out")
    params.output_prefix = "run"
    params.json_folder = str(tmp_path / "out" / "json")
    pipeline = IntegratedAnalysisPipeline(params, progress_callback=lambda *_: None)

    field = pipeline._process_field_data(
        field_freq=600.0,
        t1_file=_matrix(tmp_path / "T1.csv", [0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2], 0.7),
        t2_file=_matrix(tmp_path / "T2.csv", [0, 10, 20, 40, 80, 120, 160], 90.0),
        noe_sat_file=_plane(tmp_path / "sat.csv", 0.8),
        noe_unsat_file=_plane(tmp_path / "unsat.csv", 1.0),
        field_label="field1", t1_units="s", t2_units="ms")

    assert len(field["merged_data"]) == len(RESIDUES)
    json_folder = Path(params.json_folder)
    t1 = json.loads((json_folder / "run_field1_T1_fit_data.json").read_text())
    t2 = json.loads((json_folder / "run_field1_T2_fit_data.json").read_text())
    assert t1["metadata"]["time_units"] == "s"
    assert t2["metadata"]["time_units"] == "ms"
    # the rates are per second whatever unit each series was fitted in
    assert all(1.0 < r["R1"] < 2.0 and 8.0 < r["R2"] < 14.0 for r in field["merged_data"].values())
