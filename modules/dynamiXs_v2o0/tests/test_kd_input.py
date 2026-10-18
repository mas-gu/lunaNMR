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


def test_tidy_collapses_duplicate_point_rows(tmp_path):
    # A residue listed TWICE at the same titration point must not crash (it used to
    # raise "float() ... not 'Series'"); the duplicates are averaged to one value.
    from kd_input import load_titration_tidy
    rows = [
        ("0", "R1", 8.00, 120.0, 1000.0, 2000.0),
        ("1", "R1", 8.10, 120.0,  500.0, 1000.0),
        ("1", "R1", 8.20, 120.0,  700.0, 1000.0),   # duplicate point for R1
    ]
    df = pd.DataFrame(rows, columns=["spectrum_name", "assignment",
                                     "ppm_x", "ppm_y", "height", "volume"])
    csv = tmp_path / "series_analysis_tidy.csv"
    df.to_csv(csv, index=False)
    points, res = load_titration_tidy(str(csv))
    assert points == [0.0, 1.0]
    r1 = res["R1"]
    assert all(isinstance(v, float) for v in r1["ppm_x"])   # scalars, no Series
    assert r1["ppm_x"][1] == pytest.approx(8.15)            # mean of 8.10 & 8.20


class TestDescriptiveLabels:
    def test_to_point_parses_descriptive_spectrum_labels(self):
        from kd_input import _to_point
        assert _to_point("0.0") == 0.0
        assert _to_point("10o5") == 10.5
        # Real spectrum names: trailing sequence number is the point index.
        assert _to_point("03_2D_series_ref_001") == 1.0
        assert _to_point("03_2D_series_ref_002") == 2.0
        assert _to_point("no_number_here") is None

    def test_load_tracking_reads_descriptive_labels(self, tmp_path):
        from kd_input import load_titration
        labels = ["03_ATP_001", "03_ATP_002", "03_ATP_003"]
        cols = {"Assignment": ["R1", "R2"]}
        for i, lbl in enumerate(labels):
            cols[f"{lbl}_Position_X"] = [8.0 + 0.01 * i, 9.0]
            cols[f"{lbl}_Position_Y"] = [120.0, 125.0]
            cols[f"{lbl}_Height"] = [1000.0 - 100 * i, 500.0]
            cols[f"{lbl}_Volume"] = [2000.0, 1000.0]
        csv = tmp_path / "comprehensive_peak_tracking.csv"
        pd.DataFrame(cols).to_csv(csv, index=False)
        points, residues = load_titration(str(csv))
        assert points == [1.0, 2.0, 3.0]
        assert set(residues) == {"R1", "R2"}


_DATA_TEST = Path(__file__).resolve().parents[4] / "data_test"


def _runs_with_failed_points():
    """Real point-labelled series runs that wrote both CSVs and have a Failed point."""
    from kd_input import _to_point
    runs = []
    for tracking in sorted(_DATA_TEST.glob("**/series_results_*/comprehensive_peak_tracking.csv")):
        tidy = tracking.with_name("series_analysis_tidy.csv")
        if not tidy.exists():
            continue
        rows = pd.read_csv(tidy)
        if ((rows["quality"] == "Failed").any()
                and all(_to_point(n) is not None for n in rows["spectrum_name"].unique())):
            runs.append((tidy, tracking))
    return runs


@pytest.mark.skipif(not _runs_with_failed_points(), reason="data_test not present")
@pytest.mark.parametrize("tidy,tracking", _runs_with_failed_points())
def test_tracking_blanks_the_points_the_tidy_file_of_the_same_run_marks_failed(tidy, tracking):
    # The same series must give the same Kd input whichever of its two CSVs is loaded.
    from kd_input import _to_point, load_titration_tidy, load_titration_tracking
    points, by_tidy = load_titration_tidy(str(tidy))
    _, by_tracking = load_titration_tracking(str(tracking))
    rows = pd.read_csv(tidy)
    failed = rows[rows["quality"] == "Failed"]
    assert len(failed)
    for _, row in failed.iterrows():
        rec = by_tracking[str(row["assignment"])]
        k = points.index(_to_point(row["spectrum_name"]))
        assert all(np.isnan(rec[key][k]) for key in ("ppm_x", "ppm_y", "height", "volume")), (
            row["assignment"], row["spectrum_name"])
    for name, rec in by_tidy.items():
        if name in by_tracking:
            assert by_tracking[name]["n_failed_points"] == rec["n_failed_points"], name


_BINDING_TIDY = _DATA_TEST / "data_binding" / "tidy" / "DNAJA1_HSPA8_tidy.csv"


class TestRepeatAndUnreadableLabels:
    """A second acquisition at a titration point is labelled '<point>_2' by the series
    writer. float() accepts underscores, so '1_2' read as point 12."""

    def test_repeat_suffix_names_the_same_point(self):
        from kd_input import _to_point
        assert _to_point("1_2") == 1.0
        assert _to_point("0.125_2") == 0.125
        assert _to_point("1o0_2") == 1.0

    @pytest.mark.skipif(not _BINDING_TIDY.exists(), reason="data_test not present")
    def test_a_repeat_acquisition_joins_its_point(self, tmp_path):
        from kd_input import load_titration_tidy
        rows = pd.read_csv(_BINDING_TIDY, dtype={"spectrum_name": str})
        repeat = rows[rows["spectrum_name"] == "2"].assign(spectrum_name="2_2")
        csv = tmp_path / "series_analysis_tidy.csv"
        pd.concat([rows, repeat]).to_csv(csv, index=False)
        points, _ = load_titration_tidy(str(csv))
        assert points == load_titration_tidy(str(_BINDING_TIDY))[0]

    @pytest.mark.skipif(not _BINDING_TIDY.exists(), reason="data_test not present")
    def test_a_spectrum_name_with_no_point_is_refused_not_dropped(self, tmp_path):
        # Dropping the point shifts every concentration after it by one place.
        from kd_input import load_titration_tidy
        rows = pd.read_csv(_BINDING_TIDY, dtype={"spectrum_name": str})
        rows.loc[rows["spectrum_name"] == "1", "spectrum_name"] = "2D_sample_1o0b"
        csv = tmp_path / "series_analysis_tidy.csv"
        rows.to_csv(csv, index=False)
        with pytest.raises(ValueError, match="2D_sample_1o0b"):
            load_titration_tidy(str(csv))

    def test_a_tracking_label_with_no_point_is_refused_not_dropped(self, tmp_path):
        from kd_input import load_titration_tracking
        cols = {"Assignment": ["R1"]}
        for lbl in ("0", "0o5", "2D_sample_1o0b"):
            for suffix in ("Position_X", "Position_Y", "Height", "Volume"):
                cols[f"{lbl}_{suffix}"] = [1.0]
        csv = tmp_path / "comprehensive_peak_tracking.csv"
        pd.DataFrame(cols).to_csv(csv, index=False)
        with pytest.raises(ValueError, match="2D_sample_1o0b"):
            load_titration_tracking(str(csv))


def _series_folder_with(*names):
    """A point-labelled real series run holding all of `names`."""
    for tidy, _ in _runs_with_failed_points():
        if all((tidy.parent / n).exists() for n in names):
            return tidy.parent
    return None


_MATRIX_FOLDER = _series_folder_with("peak_intensity_matrix.csv", "peak_detection_matrix.csv",
                                     "processing_summary.csv")


@pytest.mark.skipif(_MATRIX_FOLDER is None, reason="data_test not present")
class TestOnlyAnIntensityMatrixLoadsAsOne:
    """Any CSV that is neither tidy nor tracking reaches the matrix loader."""

    def test_an_intensity_matrix_loads(self):
        from kd_input import load_titration
        points, residues = load_titration(str(_MATRIX_FOLDER / "peak_intensity_matrix.csv"))
        assert points and residues

    def test_a_detection_flag_matrix_is_refused(self):
        # Every cell is 0 or 1; fitted as intensities, every residue 'fits' one Kd.
        from kd_input import load_titration
        with pytest.raises(ValueError, match="0 or 1"):
            load_titration(str(_MATRIX_FOLDER / "peak_detection_matrix.csv"))

    def test_a_csv_with_no_point_columns_is_refused(self):
        from kd_input import load_titration
        with pytest.raises(ValueError, match="no titration point"):
            load_titration(str(_MATRIX_FOLDER / "processing_summary.csv"))


def test_numeric_residue_names_read_the_same_from_every_input(tmp_path):
    # pandas reads an all-numeric Assignment column with a blank in it as float, so a
    # residue named 10 came back as '10.0' from the tracking and matrix files.
    from kd_input import load_titration
    names = ["10", "11", None]
    tidy = tmp_path / "series_analysis_tidy.csv"
    pd.DataFrame({"spectrum_name": ["0", "1"] * 3,
                  "assignment": [n for n in names for _ in (0, 1)],
                  "ppm_x": 8.0, "ppm_y": 120.0, "height": 1000.0,
                  "volume": 2000.0}).to_csv(tidy, index=False)
    tracking = tmp_path / "comprehensive_peak_tracking.csv"
    matrix = tmp_path / "peak_intensity_matrix.csv"
    wide = {"Assignment": names}
    for lbl in ("0", "1"):
        for suffix in ("Position_X", "Position_Y", "Height", "Volume"):
            wide[f"{lbl}_{suffix}"] = [1000.0] * 3
    pd.DataFrame(wide).to_csv(tracking, index=False)
    pd.DataFrame({"Assignment": names, "0": [1000.0] * 3, "1": [900.0] * 3}).to_csv(
        matrix, index=False)
    for csv in (tidy, tracking, matrix):
        assert set(load_titration(str(csv))[1]) == {"10", "11"}, csv.name
