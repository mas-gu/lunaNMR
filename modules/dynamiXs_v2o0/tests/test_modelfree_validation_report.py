# ABOUTME: The model-free pipeline's validation report covers both fields of a dual-field
# ABOUTME: run, and stays the one-field report for a single-field run.

import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))
sys.path.insert(0, str(_DIR / "dynamiXs_T1_T2"))


def _field(r2_of=None):
    """Merged rates of twelve residues; `r2_of` overrides R2 for the residues it names."""
    merged = {f"{n}AlaH": {"R1": 1.4, "R1_err": 0.03, "R2": 11.0, "R2_err": 0.2,
                           "hetNOE": 0.8, "hetNOE_err": 0.02} for n in range(3, 15)}
    for residue, r2 in (r2_of or {}).items():
        merged[residue]["R2"] = r2
    return merged


def _pipeline(tmp_path, dual):
    from dynamiXs_integrated.integrated_analysis import (IntegratedAnalysisPipeline,
                                                         IntegratedAnalysisParameters)
    params = IntegratedAnalysisParameters()
    params.output_dir = str(tmp_path)
    params.enable_dual_field = dual
    params.field1_freq_mhz, params.field2_freq_mhz = 600.0, 700.0
    pipeline = IntegratedAnalysisPipeline(params, progress_callback=lambda *_: None)
    pipeline.field1_results = {"merged_data": _field(), "excluded_residues": ["20SerH"]}
    pipeline.field2_results = {"merged_data": _field({"7AlaH": 0.9}),
                               "excluded_residues": ["31ThrH"]}
    return pipeline


def test_a_dual_field_report_names_what_only_field_2_shows(tmp_path):
    report = _pipeline(tmp_path, dual=True)._validate_datasets()

    assert "20SerH" in report and "31ThrH" in report
    assert "R2 < R1 (1 residues - WARNING)" in report and "7AlaH" in report
    field1, field2 = report.split("Field 2 (700.0 MHz)")
    assert field1.startswith("Field 1 (600.0 MHz)")
    assert "20SerH" in field1 and "NO VALIDATION ISSUES DETECTED" in field1
    assert "31ThrH" in field2 and "7AlaH" in field2 and "20SerH" not in field2


def test_a_single_field_report_is_the_field_1_report(tmp_path):
    from dynamiXs_integrated.data_validation import (generate_validation_report,
                                                     validate_relaxation_rates)
    pipeline = _pipeline(tmp_path, dual=False)

    report = pipeline._validate_datasets()

    assert report == generate_validation_report(
        validate_relaxation_rates(pipeline.field1_results["merged_data"]), ["20SerH"])
