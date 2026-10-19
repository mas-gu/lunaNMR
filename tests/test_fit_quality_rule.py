# ABOUTME: Every path that labels a fit labels the same (R-squared, success, near-noise, S/N) identically.
# ABOUTME: Real DNAJA1 fits: the series pipeline, the integrator called directly, the overlay refit, the browser.

"""
The fit-quality label decides what the Kd loader and the matrices keep: they drop a
Failed point. The series wrote its own label from R-squared and S/N, while the fitters
that produced the fit labelled it by another rule, so the same fit could carry Fair in
batch_results.json and Poor in the tidy CSV, and a refit or a single-spectrum fit got a
label the series would not have given it.
"""

import contextlib
import io
import json
import shutil
import sys
from pathlib import Path

import pytest

from lunaNMR.utils.quality_categories import fit_quality_label

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
WEAKEST_POINT = "2D_DNAJA1_HSPA8_3o0.ft"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


def _member(fit, assignment):
    """The fitted peak of `assignment` among a 2D fit's members."""
    return next(p for p in fit['all_peaks'] if str(p.get('assignment')) == str(assignment))


def _expected(fit, assignment):
    member = _member(fit, assignment)
    return fit_quality_label(fit['r_squared'], True, bool(member.get('near_noise')),
                             fit.get('detection_sigma'))


# (point, assignment, R-squared, S/N) of real series fits with R-squared in [0.5, 0.6]:
# a DNAJA1 titration and a 700 MHz T1 series.
BELOW_MEASURED_SNR = [
    ("DNAJA1 1", "32LysH", 0.5948, 3.5096),
    ("DNAJA1 2", "68GlyH", 0.5307, 3.3196),
    ("DNAJA1 3", "18ThrH", 0.5708, 4.4859),
    ("DNAJA1 3", "20GluH", 0.5912, 4.2854),
    ("DNAJA1 3", "55LeuH", 0.5792, 4.7823),
    ("DNAJA1 3", "58AlaH", 0.5579, 4.5334),
    ("T1 1800", "67GlyH", 0.5751, 4.5093),
    ("T1 2400", "45PheH", 0.5202, 4.9253),
    ("T1 2400", "47GlnH", 0.5625, 4.0572),
    ("T1 2400", "57AspH", 0.5657, 4.2408),
]
ABOVE_MEASURED_SNR = [
    ("DNAJA1 1", "30AlaH", 0.5476, 5.4073),
    ("DNAJA1 2", "19GlnH", 0.5118, 5.5555),
    ("DNAJA1 3", "3LysH", 0.5755, 5.8168),
    ("T1 1800", "11LeuH", 0.5936, 5.2520),
    ("T1 2400", "25AlaH", 0.5082, 6.6606),
]


class TestFairNeedsRSquaredAbove06:
    @pytest.mark.parametrize("point, assignment, r2, snr", BELOW_MEASURED_SNR)
    def test_a_fit_at_or_below_06_is_poor_near_the_noise(self, point, assignment, r2, snr):
        assert fit_quality_label(r2, True, False, snr) == 'Poor'

    @pytest.mark.parametrize("point, assignment, r2, snr", ABOVE_MEASURED_SNR)
    def test_a_fit_at_or_below_06_is_poor_clear_of_the_noise(self, point, assignment, r2, snr):
        assert fit_quality_label(r2, True, False, snr) == 'Poor'

    @pytest.mark.parametrize("point, assignment, r2, snr", BELOW_MEASURED_SNR + ABOVE_MEASURED_SNR)
    def test_the_label_does_not_depend_on_snr(self, point, assignment, r2, snr):
        assert {fit_quality_label(r2, True, False, s) for s in (None, 1.0, snr, 50.0)} == {'Poor'}


@pytest.fixture(scope="module")
def stored_series(dnaja1_titration_series):
    return json.loads((dnaja1_titration_series / "batch_results.json").read_text())


class TestTheSeriesPipeline:
    def test_the_fit_and_the_series_give_a_fit_the_same_label(self, stored_series):
        disagreements = []
        for name, result in stored_series['results'].items():
            series_label = {str(p['assignment']): p['quality'] for p in result['integration_results']}
            for fit in result['fitted_results']:
                if series_label[str(fit['assignment'])] != fit['fitting_quality']:
                    disagreements.append((name, fit['assignment'], fit['r_squared'],
                                          fit['fitting_quality'], series_label[str(fit['assignment'])]))
        assert not disagreements, disagreements

    def test_every_label_is_the_shared_rule(self, stored_series):
        wrong = []
        for name, result in stored_series['results'].items():
            for fit in result['fitted_results']:
                if not fit.get('success'):
                    expected = 'Failed'
                else:
                    expected = _expected(fit, fit['assignment'])
                if fit['fitting_quality'] != expected:
                    wrong.append((name, fit['assignment'], fit['r_squared'], fit['fitting_quality'], expected))
        assert not wrong, wrong

    def test_the_near_noise_and_failed_cases_are_both_present(self, stored_series):
        """The two special cases this rule exists for occur in the run the tests read."""
        fits = [f for r in stored_series['results'].values() for f in r['fitted_results']
                if f.get('success')]
        low = [f for f in fits if f['r_squared'] < 0.2]
        assert any(_member(f, f['assignment']).get('near_noise') for f in low)
        assert any(f['fitting_quality'] == 'Failed' for f in low)


@pytest.fixture(scope="module")
def integrator():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.utils.file_manager import NMRFileManager
    set_ps2d_config('15N')
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(DNAJA1 / WEAKEST_POINT))
    it.noise_level = it._estimate_noise_level_advanced()
    it.peak_list = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
    it.gui_params = {'use_ps2d_multi_peak': True, 'use_adaptive_optimization': False}
    return it


def _context(integrator):
    rows = integrator.peak_list
    return [{'x_ppm': float(r['Position_X']), 'y_ppm': float(r['Position_Y']),
             'pos_x': float(r['Position_X']), 'pos_y': float(r['Position_Y']),
             'assignment': str(r['Assignment'])} for _, r in rows.iterrows()]


class TestTheIntegratorCalledDirectly:
    def test_a_2d_fit_is_labelled_by_the_shared_rule(self, integrator):
        context = _context(integrator)
        wrong, labels = [], set()
        for peak in context:
            with contextlib.redirect_stdout(io.StringIO()):
                fit = integrator.fit_peak_voigt_2d(peak['x_ppm'], peak['y_ppm'], peak['assignment'],
                                                   all_peaks_context=context)
            if fit is None:
                continue
            fit = dict(fit, r_squared=fit['avg_r_squared'])
            labels.add(fit['fitting_quality'])
            if fit['fitting_quality'] != _expected(fit, peak['assignment']):
                wrong.append((peak['assignment'], fit['r_squared'], fit['fitting_quality'],
                              _expected(fit, peak['assignment'])))
        assert not wrong, wrong
        assert {'Fair', 'Poor', 'Failed'} <= labels


class TestTheOverlayRefit:
    @pytest.fixture
    def viewer(self, dnaja1_titration_series, tmp_path, open_overlay_viewer):
        from PySide6.QtWidgets import QApplication
        QApplication.instance() or QApplication(sys.argv)
        folder = tmp_path / "out"
        shutil.copytree(dnaja1_titration_series, folder)
        viewer = open_overlay_viewer(folder)
        yield viewer
        viewer.deleteLater()

    def test_a_refit_is_labelled_by_the_shared_rule(self, viewer):
        idx = [s['name'] for s in viewer.spectra].index(WEAKEST_POINT)
        checked = 0
        for assignment in ('28LysH', '48IleH', '3LysH'):
            peak = viewer._find_peak_obj(idx, assignment)
            cx, cy = viewer._peak_center(peak)
            r2 = viewer._refit_peak_in_spectrum(idx, assignment, cx, cy)
            assert r2 is not None
            peak = viewer._find_peak_obj(idx, assignment)
            by_noise = {fit_quality_label(r2, True, near_noise, peak['detection_sigma'])
                        for near_noise in (False, True)}
            if len(by_noise) == 1:
                assert peak['quality'] == by_noise.pop()
                checked += 1
        assert checked


class TestTheResultsBrowser:
    def test_the_status_column_shows_the_label_the_fit_carries(self, stored_series):
        from PySide6.QtWidgets import QApplication
        from lunaNMR.gui.dialogs.results_browser_dialog import ResultsBrowserDialog
        from lunaNMR.processors.multi_spectrum_processor import BatchResults
        QApplication.instance() or QApplication(sys.argv)
        batch = BatchResults()
        for name, result in stored_series['results'].items():
            batch.add_result(name, result)
        dialog = ResultsBrowserDialog(batch_results=batch)
        try:
            carried = {}
            for name, result in stored_series['results'].items():
                for fit in result['fitted_results']:
                    carried.setdefault(str(fit['assignment']), {})[name] = fit['fitting_quality']
            shown = 0
            for assignment in ('28LysH', '48IleH', '3LysH'):
                dialog.peak_combo.setCurrentText(assignment)
                dialog.on_peak_selected(dialog.peak_combo.currentIndex())
                table = dialog.tracking_table
                for row in range(table.rowCount()):
                    name = table.item(row, 0).text()
                    assert table.item(row, 6).text() == carried[assignment][name]
                    shown += 1
            assert shown
        finally:
            dialog.deleteLater()
