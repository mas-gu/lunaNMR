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
