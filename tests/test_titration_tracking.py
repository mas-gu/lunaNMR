# ABOUTME: Tests the titration-mode tracking policy (cascade + drift-off + wider window).
# ABOUTME: Titration peaks move; the policy stops the relaxation defaults suppressing that.

import pandas as pd
import pytest

from lunaNMR.core.ps2d_config import PS2DConfig
from lunaNMR.processors.multi_spectrum_processor import (
    MultiSpectrumProcessor,
    resolve_titration_tracking,
)


class TestTitrationConfigWindows:
    """Each nucleus carries an enlarged per-step search window for titration."""

    @pytest.mark.parametrize("nucleus", ["15N", "13C"])
    def test_titration_windows_exist_and_are_wider(self, nucleus):
        cfg = PS2DConfig(nucleus)
        # The titration per-step window must be wider than the relaxation window,
        # so cascade detection can follow a moving peak across a titration step.
        assert cfg.titration_search_window_f1 > cfg.search_window_f1
        assert cfg.titration_search_window_f2 > cfg.search_window_f2


class TestResolveTitrationTracking:
    """The policy: titration forces cascade + drift-off + enlarged windows."""

    def test_titration_forces_cascade_and_drift_off(self):
        cfg = PS2DConfig("15N")
        mode, drift, win_x, win_y, overridden = resolve_titration_tracking(
            series_mode="titration",
            peak_source_mode="independent",   # user picked something else
            enable_drift_limit=True,          # user left it on
            config=cfg,
        )
        assert mode == "cascade"
        assert drift is False
        assert overridden is True

    def test_titration_window_axes_mapping(self):
        # x = 1H (F2), y = 15N/13C (F1). Getting this backwards silently breaks tracking.
        cfg = PS2DConfig("15N")
        _, _, win_x, win_y, _ = resolve_titration_tracking(
            "titration", "cascade", False, cfg)
        assert win_x == cfg.titration_search_window_f2
        assert win_y == cfg.titration_search_window_f1

    def test_time_mode_is_passthrough(self):
        cfg = PS2DConfig("15N")
        result = resolve_titration_tracking(
            series_mode="time",
            peak_source_mode="cascade",
            enable_drift_limit=True,
            config=cfg,
        )
        assert result == ("cascade", True, None, None, False)

    def test_time_mode_preserves_user_choices(self):
        cfg = PS2DConfig("13C")
        mode, drift, win_x, win_y, overridden = resolve_titration_tracking(
            "time", "independent", False, cfg)
        assert (mode, drift, win_x, win_y, overridden) == (
            "independent", False, None, None, False)


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
