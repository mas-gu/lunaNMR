# ABOUTME: A delay column that does not parse must be dropped, not left as None in the array.
# ABOUTME: One None makes x an object array and the first np.exp(-x/T) raises TypeError.

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_DIR))


def _matrix(tmp_path, headers):
    """A LunaNMR-format relaxation matrix whose delay columns carry `headers`."""
    decay = lambda t: 1000.0 * np.exp(-np.asarray(t, float) / 120.0) + 5.0
    times = [8.0, 51.0, 102.0, 204.0]
    rows = [["Peak_Number", "Assignment", "Reference_X", "Reference_Y"] + list(headers)]
    for i, name in enumerate(("K3", "E4", "V5"), start=1):
        rows.append([i, name, 8.0, 120.0] + [float(v) * (1 + 0.01 * i) for v in decay(times)])
    path = tmp_path / "peak_intensity_matrix.csv"
    pd.DataFrame(rows).to_csv(path, header=False, index=False)
    return path


class TestAnUnparseableColumnIsDropped:
    """`series` mixes label kinds in one matrix: a spectrum whose filename carried no
    parseable delay gets a stem-named column beside the delay-named ones. That column
    has no time, so it cannot be fitted — but it must not poison the ones that can.
    """

    def test_a_stem_named_column_between_delays_does_not_crash(self, tmp_path):
        """np.array([...None...]) is an object array, so the first reduction over it
        raises `TypeError: '<=' not supported between instances of 'float' and
        'NoneType'` from inside numpy — reported verbatim to the user, with nothing
        naming the column responsible."""
        from fit_Tx_NMRRE import run_analysis_with_params
        csv = _matrix(tmp_path, ["8", "51", "hetnoe_600_saturated", "204"])
        result = run_analysis_with_params({
            'input_csv_file': str(csv),
            'output_prefix': str(tmp_path / 'f'),
            'results_txt_file': str(tmp_path / 'f.txt'),
            'experiment_type': 'T2', 'time_units': 'ms',
            'error_method': 'analytical', 'n_bootstrap': 0,
            'field_name': 'field1', 'field_freq': 600.0, 'json_folder': None,
        })
        assert result['n_fitted'] >= 1

    def test_the_dropped_column_is_reported(self, tmp_path):
        """Silently fitting 3 of 4 points is how a series loses a spectrum unnoticed."""
        from fit_Tx_NMRRE import run_analysis_with_params
        csv = _matrix(tmp_path, ["8", "51", "hetnoe_600_saturated", "204"])
        result = run_analysis_with_params({
            'input_csv_file': str(csv),
            'output_prefix': str(tmp_path / 'f'),
            'results_txt_file': str(tmp_path / 'f.txt'),
            'experiment_type': 'T2', 'time_units': 'ms',
            'error_method': 'analytical', 'n_bootstrap': 0,
            'field_name': 'field1', 'field_freq': 600.0, 'json_folder': None,
        })
        assert result['dropped_columns'] == ['hetnoe_600_saturated']

    def test_every_column_parsing_drops_nothing(self, tmp_path):
        from fit_Tx_NMRRE import run_analysis_with_params
        csv = _matrix(tmp_path, ["8", "51", "102", "204"])
        result = run_analysis_with_params({
            'input_csv_file': str(csv),
            'output_prefix': str(tmp_path / 'f'),
            'results_txt_file': str(tmp_path / 'f.txt'),
            'experiment_type': 'T2', 'time_units': 'ms',
            'error_method': 'analytical', 'n_bootstrap': 0,
            'field_name': 'field1', 'field_freq': 600.0, 'json_folder': None,
        })
        assert result['dropped_columns'] == []
        assert result['n_fitted'] == 3


class TestAMatrixWithNoDelaysAtAll:
    """The limiting case, and the one real datasets hit: every spectrum in
    a real methyl dataset is named 03_2D_sample_ref_00N, so no column
    carries a delay. It used to fit them as 1, 2 and 3 ms and report success.
    """

    def test_the_error_names_the_real_problem(self, tmp_path):
        """fit_Tx defaulted delay_start_idx to 0 when nothing parsed, so it tried to
        read the Assignment column as intensities and failed with 'could not convert
        string to float: <a residue label>' — naming neither the file nor the cause."""
        from fit_Tx_NMRRE import run_analysis_with_params
        csv = _matrix(tmp_path, ["ref_noCa_001", "ref_noCa_002", "ref_noCa_003", "ref_noCa_004"])
        with pytest.raises(ValueError, match="delay"):
            run_analysis_with_params({
                'input_csv_file': str(csv),
                'output_prefix': str(tmp_path / 'f'),
                'results_txt_file': str(tmp_path / 'f.txt'),
                'experiment_type': 'T2', 'time_units': 'ms',
                'error_method': 'analytical', 'n_bootstrap': 0,
                'field_name': 'field1', 'field_freq': 600.0, 'json_folder': None,
            })

    def test_the_error_names_the_file_and_the_headers(self, tmp_path):
        from fit_Tx_NMRRE import run_analysis_with_params
        csv = _matrix(tmp_path, ["ref_noCa_001", "ref_noCa_002", "ref_noCa_003", "ref_noCa_004"])
        with pytest.raises(ValueError) as exc:
            run_analysis_with_params({
                'input_csv_file': str(csv),
                'output_prefix': str(tmp_path / 'f'),
                'results_txt_file': str(tmp_path / 'f.txt'),
                'experiment_type': 'T2', 'time_units': 'ms',
                'error_method': 'analytical', 'n_bootstrap': 0,
                'field_name': 'field1', 'field_freq': 600.0, 'json_folder': None,
            })
        assert "peak_intensity_matrix.csv" in str(exc.value)
        assert "ref_noCa_001" in str(exc.value)
