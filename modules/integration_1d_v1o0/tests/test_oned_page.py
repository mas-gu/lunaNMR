# ABOUTME: Tests the interactive 1D page - peak picking by box and click, then series integration.
# ABOUTME: Runs headless via QT_QPA_PLATFORM=offscreen; drives the real widgets, not mocks.

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

_MODULE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MODULE_DIR))

pytest.importorskip('PySide6')

REAL_DIR = _MODULE_DIR.parent.parent / "data_example" / "1D"
needs_series = pytest.mark.skipif(not list(REAL_DIR.glob("*.ft1")),
                                  reason="real 1D series not present")

PEAK_A, PEAK_B = 8.1847, 8.1752


@pytest.fixture(scope='module')
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture
def page(app):
    from oned_page import OneDIntegrationPage
    return OneDIntegrationPage()


class TestFileSelection:
    def test_filter_offers_the_nmrpipe_1d_extension(self, page):
        """The picker must show .ft1 files, not only .ft/.ft2."""
        from oned_page import SPECTRUM_FILTER
        assert '*.ft1' in SPECTRUM_FILTER

    @needs_series
    def test_loads_an_explicit_list_of_ft1_files(self, page):
        chosen = [str(p) for p in sorted(REAL_DIR.glob('*.ft1'))[:4]]
        page.load_spectra(paths=chosen)
        assert len(page.paths) == 4
        assert page.spectrum is not None

    def test_unreadable_files_are_skipped_with_a_message(self, page, tmp_path):
        bogus = tmp_path / "broken.ft1"
        bogus.write_bytes(b"not NMR data")
        page.load_spectra(paths=[str(bogus)])
        assert page.paths == []
        assert 'could not' in page.status.text().lower()

    @needs_series
    def test_a_mixed_selection_keeps_the_readable_files(self, page, tmp_path):
        bogus = tmp_path / "broken.ft1"
        bogus.write_bytes(b"not NMR data")
        good = [str(p) for p in sorted(REAL_DIR.glob('*.ft1'))[:2]]
        page.load_spectra(paths=good + [str(bogus)])
        assert len(page.paths) == 2


class TestPageConstruction:
    def test_starts_with_no_peaks_and_integration_disabled(self, page):
        assert page.peaks == []
        assert page.integrate_button.isEnabled() is False

    def test_reports_an_empty_folder(self, page, tmp_path):
        page.load_folder(folder=str(tmp_path))
        assert "No 1D spectra" in page.status.text()


@needs_series
class TestRealWorkflow:
    def test_loads_the_series(self, page):
        page.load_folder(folder=str(REAL_DIR))
        assert len(page.paths) == 53
        assert page.spectrum is not None

    def test_box_select_finds_both_peaks(self, page):
        """Drag a rectangle over the region: both resonances are picked up."""
        from oned_series import peaks_in_range
        page.load_folder(folder=str(REAL_DIR))
        page.add_peaks(peaks_in_range(page.spectrum, 8.20, 8.16))
        positions = sorted(p['position'] for p in page.peaks)
        assert len(positions) == 2
        assert positions[1] == pytest.approx(PEAK_A, abs=2e-3)
        assert positions[0] == pytest.approx(PEAK_B, abs=2e-3)

    def test_click_select_snaps_to_the_peak(self, page):
        page.load_folder(folder=str(REAL_DIR))
        page.add_peak_at(8.1840)                       # deliberately off-apex
        assert page.peaks[0]['position'] == pytest.approx(PEAK_A, abs=2e-3)

    def test_duplicate_pick_is_ignored(self, page):
        page.load_folder(folder=str(REAL_DIR))
        page.add_peak_at(8.1845)
        page.add_peak_at(8.1845)
        assert len(page.peaks) == 1

    def test_removing_and_clearing_peaks(self, page):
        page.load_folder(folder=str(REAL_DIR))
        page.add_peak_at(PEAK_A)
        page.add_peak_at(PEAK_B)
        page.peak_table.setCurrentCell(0, 0)
        page.remove_selected_peak()
        assert len(page.peaks) == 1
        page.clear_peaks()
        assert page.peaks == []
        assert page.integrate_button.isEnabled() is False

    def test_integration_writes_the_expected_csv(self, page, tmp_path, monkeypatch):
        from PySide6.QtWidgets import QMessageBox
        monkeypatch.setattr(QMessageBox, 'information', lambda *a, **k: None)

        page.load_folder(folder=str(REAL_DIR))
        page.add_peak_at(PEAK_A)
        page.add_peak_at(PEAK_B)
        assert page.integrate_button.isEnabled() is True

        out = tmp_path / "series.csv"
        page.integrate(out_path=str(out))

        lines = out.read_text().strip().splitlines()
        assert lines[0] == 'spectrum,peak_1,peak_2'
        assert len(lines) == 54                          # header + 53 spectra
        assert lines[1].startswith('1D_KB_GTP_001')

        first = [float(v) for v in lines[1].split(',')[1:]]
        last = [float(v) for v in lines[-1].split(',')[1:]]
        assert last[0] / first[0] > 5.0                  # peak at 8.1847 grows
        assert last[1] / first[1] < 0.5                  # peak at 8.1752 decays

    def test_area_observable_writes_a_different_table(self, page, tmp_path, monkeypatch):
        from PySide6.QtWidgets import QMessageBox
        monkeypatch.setattr(QMessageBox, 'information', lambda *a, **k: None)

        page.load_folder(folder=str(REAL_DIR))
        page.add_peak_at(PEAK_A)

        intensity = tmp_path / "i.csv"
        page.integrate(out_path=str(intensity))
        page.observable.setCurrentIndex(1)
        area = tmp_path / "a.csv"
        page.integrate(out_path=str(area))

        assert intensity.read_text() != area.read_text()


@needs_series
class TestPlotterInteraction:
    def test_box_drag_emits_the_peaks_inside_it(self, app):
        from oned_loader import load_spectrum
        from oned_plotter import OneDPlotter

        plotter = OneDPlotter()
        plotter.set_spectrum(load_spectrum(sorted(REAL_DIR.glob('*.ft1'))[50]))

        captured = []
        plotter.peaks_selected.connect(captured.append)

        class _Event:
            def __init__(self, x):
                self.xdata = x
        plotter._on_box(_Event(8.20), _Event(8.16))

        # the box returns everything detectable inside it, which here includes a weak
        # third line at 8.195 (S/N 12) alongside the two strong ones
        assert len(captured) == 1
        found = sorted(p['ppm'] for p in captured[0])
        assert len(found) >= 2
        assert any(abs(p - PEAK_A) < 2e-3 for p in found)
        assert any(abs(p - PEAK_B) < 2e-3 for p in found)

    def test_y_scale_changes_the_axis_limits(self, app):
        from oned_loader import load_spectrum
        from oned_plotter import OneDPlotter

        plotter = OneDPlotter()
        plotter.set_spectrum(load_spectrum(sorted(REAL_DIR.glob('*.ft1'))[50]))
        top_before = plotter.axes.get_ylim()[1]
        plotter.y_scale.setValue(10.0)
        assert plotter.axes.get_ylim()[1] < top_before
