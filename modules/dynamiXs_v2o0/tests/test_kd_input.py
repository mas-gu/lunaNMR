# ABOUTME: Tests the titration input adapter that reads the LunaNMR series tidy CSV.
# ABOUTME: Per-residue observable arrays (positions + intensities) aligned to sorted points.

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_KD_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_Kd"
sys.path.insert(0, str(_KD_DIR))


def _tidy(tmp_path):
    rows = [
        # spectrum_name (point), assignment, ppm_x, ppm_y, height, volume
        ("0",   "R1", 8.00, 120.0, 1000.0, 2000.0),
        ("1",   "R1", 8.05, 120.3,  600.0, 1200.0),
        ("0.5", "R1", 8.02, 120.1,  800.0, 1600.0),
        ("0",   "R2", 7.00, 118.0, 1000.0, 2000.0),
        ("1",   "R2", 7.00, 118.0, 1000.0, 2000.0),   # non-mover
        ("0.5", "R2", 7.00, 118.0, 1000.0, 2000.0),
    ]
    df = pd.DataFrame(rows, columns=["spectrum_name", "assignment",
                                     "ppm_x", "ppm_y", "height", "volume"])
    p = tmp_path / "series_analysis_tidy.csv"
    df.to_csv(p, index=False)
    return str(p)


class TestLoadTidy:
    def test_points_sorted_ascending(self, tmp_path):
        from kd_input import load_titration_tidy
        points, _ = load_titration_tidy(_tidy(tmp_path))
        assert points == [0.0, 0.5, 1.0]

    def test_per_residue_arrays_aligned_to_points(self, tmp_path):
        from kd_input import load_titration_tidy
        _, residues = load_titration_tidy(_tidy(tmp_path))
        assert set(residues) == {"R1", "R2"}
        r1 = residues["R1"]
        # arrays follow the sorted point order [0, 0.5, 1]
        assert r1["ppm_x"] == [8.00, 8.02, 8.05]
        assert r1["height"] == [1000.0, 800.0, 600.0]

    def test_missing_point_is_nan(self, tmp_path):
        from kd_input import load_titration_tidy
        # drop R1 at point 1 (spectrum_name reads back as float 1.0)
        df = pd.read_csv(_tidy(tmp_path))
        df = df[~((df.assignment == "R1") & (df.spectrum_name.astype(float) == 1.0))]
        p = tmp_path / "t2.csv"
        df.to_csv(p, index=False)
        _, residues = load_titration_tidy(str(p))
        assert np.isnan(residues["R1"]["height"][-1])  # point 1 missing -> nan


class TestLoadTracking:
    def test_wide_tracking_csv(self, tmp_path):
        from kd_input import load_titration_tracking, load_titration
        # wide format: Assignment + {label}_Position_X/Y/Height/Volume
        df = pd.DataFrame([{
            'Assignment': 'R1',
            '0_Position_X': 8.0, '0_Position_Y': 120.0, '0_Height': 1000.0, '0_Volume': 2000.0,
            '0.5_Position_X': 8.02, '0.5_Position_Y': 120.1, '0.5_Height': 800.0, '0.5_Volume': 1600.0,
            '1_Position_X': 8.05, '1_Position_Y': 120.3, '1_Height': 600.0, '1_Volume': 1200.0,
        }])
        p = tmp_path / "comprehensive_peak_tracking.csv"
        df.to_csv(p, index=False)
        points, residues = load_titration_tracking(str(p))
        assert points == [0.0, 0.5, 1.0]
        assert residues['R1']['ppm_x'] == [8.0, 8.02, 8.05]
        # dispatcher routes a no-spectrum_name file to the tracking loader
        pts2, _ = load_titration(str(p))
        assert pts2 == [0.0, 0.5, 1.0]


class TestLoadIntensityMatrix:
    def test_matrix_rows_residues_cols_points(self, tmp_path):
        from kd_input import load_titration
        # peak_intensity matrix: Peak_Number, Assignment, Reference_X/Y, then point cols
        df = pd.DataFrame([
            {'Peak_Number': 1, 'Assignment': 'R1', 'Reference_X': 8.0, 'Reference_Y': 120.0,
             '0': 1000.0, '0.5': 700.0, '1': 400.0},
            {'Peak_Number': 2, 'Assignment': 'R2', 'Reference_X': 7.0, 'Reference_Y': 118.0,
             '0': 2000.0, '0.5': 2000.0, '1': 2000.0},
        ])
        p = tmp_path / "peak_intensity_matrix.csv"
        df.to_csv(p, index=False)
        points, residues = load_titration(str(p))
        assert points == [0.0, 0.5, 1.0]
        assert residues['R1']['height'] == [1000.0, 700.0, 400.0]
        # no positions in a matrix -> NaN (CSP unavailable)
        assert all(np.isnan(v) for v in residues['R1']['ppm_x'])


class TestCspSeries:
    def test_csp_series_relative_to_first_point(self, tmp_path):
        from kd_input import load_titration_tidy, csp_series
        _, residues = load_titration_tidy(_tidy(tmp_path))
        csp = csp_series(residues["R1"], alpha=0.14)
        assert csp[0] == pytest.approx(0.0)  # reference point
        assert csp[-1] > csp[0]              # mover grows
        # R2 (non-mover) -> all ~0
        csp2 = csp_series(residues["R2"], alpha=0.14)
        assert max(csp2) == pytest.approx(0.0)

    def test_intensity_ratio_bad_reference_is_nan(self):
        from kd_input import intensity_ratio_series
        import math
        for ref in (0.0, float('nan'), -5.0):           # zero / missing / negative ref
            r = intensity_ratio_series({'height': [ref, 100.0, 50.0]})
            assert all(math.isnan(x) for x in r)        # excluded, not raw garbage
        # valid reference -> normal ratio
        r = intensity_ratio_series({'height': [200.0, 100.0, 50.0]})
        assert r == [1.0, 0.5, 0.25]

    def test_undetected_position_zero_is_excluded(self):
        from kd_input import csp_series
        # point 2 undetected -> position 0.0 sentinel; must be NaN, not a huge CSP
        residue = {'ppm_x': [8.0, 8.05, 0.0], 'ppm_y': [120.0, 120.3, 0.0]}
        csp = csp_series(residue, alpha=0.14)
        assert csp[0] == pytest.approx(0.0)
        assert csp[1] > 0.0
        assert np.isnan(csp[2])              # excluded, not ~7 ppm spurious CSP
