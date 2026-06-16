# ABOUTME: Tests the titration-mode tracking policy (cascade + drift-off + wider window).
# ABOUTME: Titration peaks move; the policy stops the relaxation defaults suppressing that.

import pandas as pd
import pytest

from lunaNMR.core.ps2d_config import PS2DConfig
from lunaNMR.processors.multi_spectrum_processor import (
    MultiSpectrumProcessor,
    resolve_titration_tracking,
)


class TestTitrationPositionMargins:
    """Titration mode widens the fit's per-step position freedom (pos_margin)."""

    @pytest.mark.parametrize("nucleus", ["15N", "13C"])
    def test_titration_pos_margin_is_wider(self, nucleus):
        cfg = PS2DConfig(nucleus)
        # Must exceed the relaxation pos_margin so the fit can follow a moving peak.
        assert cfg.titration_pos_margin_f1 > cfg.pos_margin_f1
        assert cfg.titration_pos_margin_f2 > cfg.pos_margin_f2

    def test_integrator_defaults_to_no_override(self):
        from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
        integ = EnhancedVoigtIntegrator()
        # None means the fitter falls back to the relaxation config defaults.
        assert integ.titration_pos_margin_f1 is None
        assert integ.titration_pos_margin_f2 is None


class TestResolveTitrationTracking:
    """The policy: titration forces cascade propagation and disables the
    absolute (free-state) drift clamp so peaks can accumulate position changes.

    Note: the per-step movement is bounded by the fit's pos_margin, not by a
    search window (cascade spectrum 2+ skips detection), so this policy
    intentionally does NOT touch the search window.
    """

    def test_titration_forces_cascade_and_drift_off(self):
        mode, drift, overridden = resolve_titration_tracking(
            series_mode="titration",
            peak_source_mode="independent",   # user picked something else
            enable_drift_limit=True,          # user left it on
        )
        assert mode == "cascade"
        assert drift is False
        assert overridden is True

    def test_time_mode_is_passthrough(self):
        result = resolve_titration_tracking(
            series_mode="time",
            peak_source_mode="cascade",
            enable_drift_limit=True,
        )
        assert result == ("cascade", True, False)

    def test_time_mode_preserves_user_choices(self):
        result = resolve_titration_tracking("time", "independent", False)
        assert result == ("independent", False, False)


class TestTitrationWindow:
    """Per-step tracking window scales with the reference spectrum's average LW."""

    @staticmethod
    def _proc(ref_stats):
        proc = MultiSpectrumProcessor.__new__(MultiSpectrumProcessor)
        proc.reference_spectrum_statistics = ref_stats
        return proc

    def test_window_is_factor_times_reference_lw(self):
        proc = self._proc({'lw_f1_median': 0.16, 'lw_f2_median': 0.0227})
        wf1, wf2 = proc._titration_window(0.3, 0.05)
        assert round(wf1, 4) == 0.64    # 4 x 0.16
        assert round(wf2, 4) == 0.0908  # 4 x 0.0227

    def test_window_capped(self):
        proc = self._proc({'lw_f1_median': 0.5, 'lw_f2_median': 0.5})
        assert proc._titration_window(0.3, 0.05) == (
            MultiSpectrumProcessor.TITRATION_WINDOW_CAP_F1,
            MultiSpectrumProcessor.TITRATION_WINDOW_CAP_F2)

    def test_window_falls_back_without_reference_stats(self):
        assert self._proc(None)._titration_window(0.3, 0.05) == (0.3, 0.05)


class TestTidyResultsFile:
    """series_analysis_tidy.csv must be created (it was silently dropped by a
    duplicate _create_comprehensive_output_files that shadowed the tidy call)."""

    def test_tidy_file_created_with_expected_columns_and_precision(self, tmp_path):
        p = MultiSpectrumProcessor.__new__(MultiSpectrumProcessor)
        p.output_folder = str(tmp_path)
        batch = {'results': {'spec_1o0.ft': {'success': True, 'integration_results': [
            {'assignment': 'R1', 'peak_number': 1, 'ppm_x': 8.1234, 'ppm_y': 120.456,
             'height': 1000.0, 'volume': 2000.0, 'snr': 8.0, 'quality': 'Good', 'r_squared': 0.99},
        ]}}}
        p._create_tidy_results_file(batch)
        f = tmp_path / 'series_analysis_tidy.csv'
        assert f.exists()
        df = pd.read_csv(f)
        assert {'spectrum_name', 'assignment', 'ppm_x', 'ppm_y', 'height', 'volume'} <= set(df.columns)
        assert df['ppm_x'].iloc[0] == 8.1234        # 4-decimal precision (CSP needs it)


class TestCarryForwardSeed:
    """A failed peak keeps its last-good position as the seed for the next spectrum."""

    @staticmethod
    def _proc():
        return MultiSpectrumProcessor.__new__(MultiSpectrumProcessor)

    def test_failed_peak_carries_prior_position(self):
        prior = [{'assignment': 'A', 'success': True, 'ppm_x': 8.0}]
        current = [{'assignment': 'A', 'success': False, 'ppm_x': 0.0}]
        seed = self._proc()._carry_forward_seed(current, prior)
        assert seed[0]['ppm_x'] == 8.0
        assert seed[0]['success'] is True

    def test_successful_peak_keeps_current(self):
        prior = [{'assignment': 'A', 'success': True, 'ppm_x': 8.0}]
        current = [{'assignment': 'A', 'success': True, 'ppm_x': 8.2}]
        seed = self._proc()._carry_forward_seed(current, prior)
        assert seed[0]['ppm_x'] == 8.2

    def test_no_prior_returns_current(self):
        current = [{'assignment': 'A', 'success': False, 'ppm_x': 0.0}]
        assert self._proc()._carry_forward_seed(current, None) is current

    def test_failed_with_mismatched_assignment_kept_as_is(self):
        # Different assignments at the same index -> lists not aligned, don't carry.
        prior = [{'assignment': 'B', 'success': True, 'ppm_x': 7.0}]
        current = [{'assignment': 'A', 'success': False, 'ppm_x': 0.0}]
        seed = self._proc()._carry_forward_seed(current, prior)
        assert seed[0]['assignment'] == 'A' and seed[0]['success'] is False

    def test_assignmentless_failure_carries_by_index(self):
        # Real failures are {'success': False, 'error': ...} with NO assignment;
        # they must still carry the prior position (matched by list position).
        prior = [{'assignment': 'A', 'success': True, 'peak_x': 8.0, 'peak_y': 120.0}]
        current = [{'success': False, 'error': 'fit blew up'}]
        seed = self._proc()._carry_forward_seed(current, prior)
        assert seed[0]['peak_x'] == 8.0 and seed[0]['success'] is True

    def test_chained_failure_keeps_last_good(self):
        # Good at n, fails at n+1 and n+2 -> seed stays at n's position, not reset.
        proc = self._proc()
        n = [{'assignment': 'A', 'success': True, 'peak_x': 8.0}]
        seed1 = proc._carry_forward_seed([{'success': False}], n)       # n+1 fails
        assert seed1[0]['peak_x'] == 8.0
        seed2 = proc._carry_forward_seed([{'success': False}], seed1)   # n+2 fails
        assert seed2[0]['peak_x'] == 8.0

    def test_length_mismatch_falls_back(self):
        prior = [{'success': True, 'peak_x': 8.0}, {'success': True, 'peak_x': 7.0}]
        current = [{'success': False}]
        assert self._proc()._carry_forward_seed(current, prior) is current


class TestFailedFitFallback:
    """A failed fit in cascade/titration mode reports the n-1 seed position and
    the intensity sampled there, instead of collapsing to the reference position
    with zero intensity. A weak-but-visible peak the optimizer can't localize
    keeps the carried position and the value measured at it."""

    @staticmethod
    def _proc(mode):
        from types import SimpleNamespace
        p = MultiSpectrumProcessor.__new__(MultiSpectrumProcessor)
        p.peak_source_mode = mode
        # Reference (spectrum 1) position differs from the carried n-1 position.
        p.reference_peaks = pd.DataFrame({'Assignment': ['A'],
                                          'Position_X': [8.0], 'Position_Y': [120.0]})
        # integrator.peak_list holds the n-1 seed position and the intensity
        # sampled at that position in the current spectrum (index-aligned).
        p.integrator = SimpleNamespace(peak_list=pd.DataFrame({
            'Assignment': ['A'], 'Position_X': [8.3], 'Position_Y': [121.0],
            'Height': [5000.0], 'Intensity': [5000.0]}))
        return p

    def test_cascade_assignmentless_failure_uses_nminus1(self):
        # Real failures are {'success': False, 'error': ...} with NO assignment,
        # so the seed must be found by list position, not assignment.
        p = self._proc('cascade')
        out = p._convert_to_integration_format([{'success': False, 'error': 'weak'}])
        assert out[0]['ppm_x'] == 8.3        # n-1 seed, not reference 8.0
        assert out[0]['ppm_y'] == 121.0
        assert out[0]['height'] == 5000.0    # sampled value, not 0
        assert out[0]['volume'] == 0.0       # no integral computed for a failed fit
        assert out[0]['Position_X'] == 8.3   # legacy field too
        assert out[0]['Height'] == 5000.0
        # Assignment must come from the seed, not the synthetic 'Peak_N' — else the
        # tidy/tracking output files the value away from residue A's series.
        assert out[0]['assignment'] == 'A'
        assert out[0]['Assignment'] == 'A'

    def test_non_cascade_failure_keeps_reference_and_zero(self):
        p = self._proc('reference')
        out = p._convert_to_integration_format([{'assignment': 'A', 'success': False}])
        assert out[0]['ppm_x'] == 8.0        # reference, unchanged behavior
        assert out[0]['height'] == 0.0


class TestProcessorWiring:
    """process_nmr_series applies the titration policy at setup."""

    @staticmethod
    def _run(tmp_path, series_mode, requested_mode):
        proc = MultiSpectrumProcessor(
            {'detection_params': {}, 'gui_params': {}, 'fitting_params': {},
             'processing_options': {}})
        ref = pd.DataFrame({'Assignment': ['R1'],
                            'Position_X': [8.0], 'Position_Y': [120.0]})
        proc.process_nmr_series(
            nmr_files=[], reference_peaks=ref, output_folder=str(tmp_path),
            peak_source_mode=requested_mode, series_mode=series_mode)
        return proc

    def test_titration_forces_cascade(self, tmp_path):
        proc = self._run(tmp_path, "titration", "independent")
        assert proc.peak_source_mode == "cascade"
        assert proc.titration_tracking is True

    def test_time_mode_keeps_requested_mode(self, tmp_path):
        proc = self._run(tmp_path, "time", "independent")
        assert proc.peak_source_mode == "independent"
        assert proc.titration_tracking is False
