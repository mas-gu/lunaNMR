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
