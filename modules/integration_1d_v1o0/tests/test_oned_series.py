# ABOUTME: Tests for series intensity extraction across an ordered set of 1D spectra.
# ABOUTME: Covers peak selection from a zoom box, per-point tracking and the results table.

import sys
from pathlib import Path

import numpy as np
import pytest

_MODULE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MODULE_DIR))

REAL_DIR = _MODULE_DIR.parent.parent / "data_example" / "1D"
needs_series = pytest.mark.skipif(not list(REAL_DIR.glob("*.ft1")),
                                  reason="real 1D series not present")

PEAK_A, PEAK_B = 8.1847, 8.1752


def _synthetic_series(n_points=5, centres=(5.0, 5.01), start=(1.0, 2.0), end=(2.0, 1.0)):
    """A conversion series: peak 0 grows, peak 1 decays, sum conserved."""
    from oned_loader import from_arrays
    from oned_voigt import voigt_1d
    ppm = np.linspace(5.5, 4.5, 20001)
    spectra = []
    for k in range(n_points):
        f = k / max(1, n_points - 1)
        y = np.zeros_like(ppm)
        for c, a0, a1 in zip(centres, start, end):
            y = y + voigt_1d(ppm, a0 + f * (a1 - a0), c, 0.0008, 0.0004, 0.0)
        spectra.append(from_arrays(y, ppm))
    return spectra


class TestPeakSelection:
    def test_detects_peaks_inside_a_zoom_box(self):
        """Box select: the user drags a rectangle, peaks inside it are picked up."""
        from oned_series import peaks_in_range
        spec = _synthetic_series(1)[0]
        found = peaks_in_range(spec, 5.02, 4.99)
        assert len(found) == 2
        assert found[0]['ppm'] == pytest.approx(5.01, abs=1e-3)

    def test_box_accepts_bounds_in_either_order(self):
        from oned_series import peaks_in_range
        spec = _synthetic_series(1)[0]
        assert len(peaks_in_range(spec, 4.99, 5.02)) == 2

    def test_empty_box_finds_nothing(self):
        from oned_series import peaks_in_range
        spec = _synthetic_series(1)[0]
        assert peaks_in_range(spec, 4.7, 4.6) == []

    def test_manual_position_snaps_to_the_local_maximum(self):
        """Click select: the clicked ppm need not be exactly on the apex."""
        from oned_series import peak_at_position
        spec = _synthetic_series(1)[0]
        peak = peak_at_position(spec, 5.0008)
        assert peak['ppm'] == pytest.approx(5.0, abs=1e-3)

    def test_manual_position_far_from_any_peak_reports_the_click(self):
        from oned_series import peak_at_position
        spec = _synthetic_series(1)[0]
        peak = peak_at_position(spec, 4.7)
        assert peak['ppm'] == pytest.approx(4.7, abs=1e-2)
        assert peak['detected'] is False


class TestSeriesExtraction:
    def test_returns_one_row_per_spectrum_and_peak(self):
        from oned_series import integrate_series
        spectra = _synthetic_series(5)
        table = integrate_series(spectra, [{'assignment': 'A', 'position': 5.0},
                                           {'assignment': 'B', 'position': 5.01}])
        assert len(table) == 5 * 2

    def test_tracks_the_expected_trend(self):
        from oned_series import integrate_series
        spectra = _synthetic_series(5, start=(1.0, 2.0), end=(2.0, 1.0))
        table = integrate_series(spectra, [{'assignment': 'A', 'position': 5.0},
                                           {'assignment': 'B', 'position': 5.01}])
        a = [r['height'] for r in table if r['assignment'] == 'A']
        b = [r['height'] for r in table if r['assignment'] == 'B']
        assert a[-1] > a[0]
        assert b[-1] < b[0]

    def test_two_peaks_never_collapse_onto_one(self):
        """The failure that makes a ratio read 1.0; guarded at every series point."""
        from oned_series import integrate_series
        spectra = _synthetic_series(3)
        table = integrate_series(spectra, [{'assignment': 'A', 'position': 5.0},
                                           {'assignment': 'B', 'position': 5.01}])
        for k in range(3):
            row = [r for r in table if r['point'] == k]
            centres = [r['ppm'] for r in row if r['ppm'] is not None]
            assert len(set(centres)) == len(centres)

    def test_carries_point_index_and_spectrum_name(self):
        from oned_series import integrate_series
        spectra = _synthetic_series(3)
        table = integrate_series(spectra, [{'assignment': 'A', 'position': 5.0}],
                                 names=['first', 'second', 'third'])
        assert [r['spectrum'] for r in table] == ['first', 'second', 'third']
        assert [r['point'] for r in table] == [0, 1, 2]

    def test_reports_area_alongside_intensity(self):
        from oned_series import integrate_series
        table = integrate_series(_synthetic_series(2), [{'assignment': 'A', 'position': 5.0}])
        assert all(r['area'] is not None for r in table)

    def test_applies_per_spectrum_intensity_scale(self):
        from oned_series import integrate_series
        spectra = _synthetic_series(2)
        plain = integrate_series(spectra, [{'assignment': 'A', 'position': 5.0}])
        spectra[1].intensity_scale = 2.0
        scaled = integrate_series(spectra, [{'assignment': 'A', 'position': 5.0}])
        assert scaled[1]['height'] == pytest.approx(2.0 * plain[1]['height'], rel=1e-6)

    def test_series_to_matrix_pivots_peaks_against_points(self):
        from oned_series import integrate_series, series_matrix
        table = integrate_series(_synthetic_series(4),
                                 [{'assignment': 'A', 'position': 5.0},
                                  {'assignment': 'B', 'position': 5.01}])
        rows, matrix = series_matrix(table, value='height')
        assert rows == ['A', 'B']
        assert matrix.shape == (2, 4)


class TestCsvExport:
    def test_writes_one_row_per_spectrum_and_one_column_per_peak(self, tmp_path):
        from oned_series import integrate_series, write_series_csv
        table = integrate_series(_synthetic_series(4),
                                 [{'assignment': 'A', 'position': 5.0},
                                  {'assignment': 'B', 'position': 5.01}],
                                 names=['s1', 's2', 's3', 's4'])
        out = tmp_path / "series.csv"
        write_series_csv(table, out)

        lines = out.read_text().strip().splitlines()
        assert lines[0].split(',') == ['spectrum', 'A', 'B']
        assert len(lines) == 5                        # header + 4 spectra
        assert lines[1].split(',')[0] == 's1'

    def test_values_match_the_table(self, tmp_path):
        from oned_series import integrate_series, write_series_csv
        peaks = [{'assignment': 'A', 'position': 5.0}]
        table = integrate_series(_synthetic_series(3), peaks, names=['a', 'b', 'c'])
        out = tmp_path / "series.csv"
        write_series_csv(table, out)

        first = float(out.read_text().strip().splitlines()[1].split(',')[1])
        assert first == pytest.approx(table[0]['height'], rel=1e-6)

    def test_can_export_area_instead_of_intensity(self, tmp_path):
        from oned_series import integrate_series, write_series_csv
        table = integrate_series(_synthetic_series(2), [{'assignment': 'A', 'position': 5.0}])
        out = tmp_path / "area.csv"
        write_series_csv(table, out, value='area')
        assert float(out.read_text().strip().splitlines()[1].split(',')[1]) == pytest.approx(
            table[0]['area'], rel=1e-6)

    def test_unmeasured_peak_writes_an_empty_cell(self, tmp_path):
        from oned_series import write_series_csv
        table = [{'point': 0, 'spectrum': 's1', 'assignment': 'A', 'height': 1.0},
                 {'point': 0, 'spectrum': 's1', 'assignment': 'B', 'height': None}]
        out = tmp_path / "gap.csv"
        write_series_csv(table, out)
        assert out.read_text().strip().splitlines()[1] == 's1,1.0,'


@needs_series
class TestRealSeries:
    def _paths(self):
        return sorted(REAL_DIR.glob("*.ft1"))

    def test_extracts_the_hydrolysis_trend(self):
        from oned_series import integrate_series, load_series
        spectra = load_series(self._paths())
        table = integrate_series(spectra, [{'assignment': 'GDP', 'position': PEAK_A},
                                           {'assignment': 'GTP', 'position': PEAK_B}])
        a = [r['height'] for r in table if r['assignment'] == 'GDP']
        b = [r['height'] for r in table if r['assignment'] == 'GTP']
        assert len(a) == len(b) == 53
        assert a[-1] / a[0] > 5.0        # product grows strongly
        assert b[-1] / b[0] < 0.5        # substrate decays

    def test_the_pair_is_resolved_at_every_point(self):
        from oned_series import integrate_series, load_series
        spectra = load_series(self._paths())
        table = integrate_series(spectra, [{'assignment': 'GDP', 'position': PEAK_A},
                                           {'assignment': 'GTP', 'position': PEAK_B}])
        for k in range(53):
            centres = [r['ppm'] for r in table if r['point'] == k and r['ppm'] is not None]
            assert len(set(centres)) == 2, f"peaks merged at point {k}"
