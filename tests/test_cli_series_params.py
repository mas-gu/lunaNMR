# ABOUTME: `series --params` exposes the tuning knobs docs/CLI.md tells you to tune.
# ABOUTME: Unknown keys must be refused: a typo that silently does nothing is the whole bug class.

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))


def _inputs(tmp_path):
    (tmp_path / "s_8ms.ft").write_bytes(b"")
    peaks = tmp_path / "p.txt"
    peaks.write_text("Assignment\tPosition_X\tPosition_Y\nK3\t8.0\t120.0\n")
    return peaks


def _captured_params(monkeypatch, tmp_path, extra):
    """The params dict the CLI hands MultiSpectrumProcessor."""
    import lunaNMR.processors.multi_spectrum_processor as msp
    captured = {}

    class _Spy:
        def __init__(self, params):
            captured.update(params)

        def process_nmr_series(self, *a, **kw):
            raise SystemExit(0)

    monkeypatch.setattr(msp, 'MultiSpectrumProcessor', _Spy)
    peaks = _inputs(tmp_path)
    from lunaNMR.cli import main
    with pytest.raises(SystemExit):
        main(["series", "--spectra", str(tmp_path), "--peaks", str(peaks),
              "--out", str(tmp_path / "o"), *extra])
    return captured


class TestTuningKnobsAreReachable:
    """docs/CLI.md advises tuning search_window/min_r_squared, and the CLI froze them at
    the GUI defaults. A series that fit badly had no headless remediation at all.
    """

    def test_a_nested_value_is_applied(self, monkeypatch, tmp_path):
        cfg = tmp_path / "tune.json"
        cfg.write_text(json.dumps({"detection_params": {"search_window_x": 0.05}}))
        params = _captured_params(monkeypatch, tmp_path, ["--params", str(cfg)])
        assert params["detection_params"]["search_window_x"] == 0.05

    def test_untouched_keys_keep_their_defaults(self, monkeypatch, tmp_path):
        """A merge, not a replacement: naming one knob must not blank the section."""
        cfg = tmp_path / "tune.json"
        cfg.write_text(json.dumps({"detection_params": {"search_window_x": 0.05}}))
        params = _captured_params(monkeypatch, tmp_path, ["--params", str(cfg)])
        assert params["detection_params"]["search_window_y"] == 0.070
        assert params["fitting_params"]["min_r_squared"] == 0.5

    def test_several_sections_at_once(self, monkeypatch, tmp_path):
        cfg = tmp_path / "tune.json"
        cfg.write_text(json.dumps({
            "detection_params": {"search_window_x": 0.04, "search_window_y": 0.04},
            "fitting_params": {"min_r_squared": 0.8},
            "gui_params": {"max_peaks_fit": 120}}))
        params = _captured_params(monkeypatch, tmp_path, ["--params", str(cfg)])
        assert params["detection_params"]["search_window_y"] == 0.04
        assert params["fitting_params"]["min_r_squared"] == 0.8
        assert params["gui_params"]["max_peaks_fit"] == 120


class TestAnUnknownKeyIsRefused:
    """The failure this whole surface exists to prevent: a setting that is accepted,
    does nothing, and reports success. A typo'd knob must be an error, not a no-op.
    """

    def _err(self, tmp_path, payload):
        cfg = tmp_path / "tune.json"
        cfg.write_text(json.dumps(payload))
        peaks = _inputs(tmp_path)
        from lunaNMR.cli import main
        return main(["series", "--spectra", str(tmp_path), "--peaks", str(peaks),
                     "--out", str(tmp_path / "o"), "--params", str(cfg)])

    def test_a_misspelled_knob_is_an_error(self, tmp_path, capsys):
        assert self._err(tmp_path, {"detection_params": {"search_window": 0.05}}) == 1
        assert "search_window" in capsys.readouterr().err

    def test_the_error_suggests_the_real_key(self, tmp_path, capsys):
        self._err(tmp_path, {"detection_params": {"search_window": 0.05}})
        assert "search_window_x" in capsys.readouterr().err

    def test_an_unknown_section_is_an_error(self, tmp_path):
        assert self._err(tmp_path, {"detection": {"search_window_x": 0.05}}) == 1

    def test_a_non_object_section_is_an_error(self, tmp_path):
        assert self._err(tmp_path, {"detection_params": 0.05}) == 1

    def test_malformed_json_is_an_error(self, tmp_path):
        cfg = tmp_path / "tune.json"
        cfg.write_text("{not json")
        peaks = _inputs(tmp_path)
        from lunaNMR.cli import main
        assert main(["series", "--spectra", str(tmp_path), "--peaks", str(peaks),
                     "--out", str(tmp_path / "o"), "--params", str(cfg)]) == 1


class TestCommandLineFlagsWin:
    """--params is a file; --parallel and --peak-source are what the caller just typed.
    The typed thing wins, or a stale file silently overrides the command line."""

    def test_parallel_flag_beats_the_file(self, monkeypatch, tmp_path):
        cfg = tmp_path / "tune.json"
        cfg.write_text(json.dumps({"processing_options": {"use_parallel_processing": False}}))
        params = _captured_params(monkeypatch, tmp_path, ["--params", str(cfg), "--parallel"])
        assert params["processing_options"]["use_parallel_processing"] is True

    def test_independent_mode_beats_the_file(self, monkeypatch, tmp_path):
        cfg = tmp_path / "tune.json"
        cfg.write_text(json.dumps(
            {"processing_options": {"rerun_adaptive_per_spectrum": False}}))
        params = _captured_params(monkeypatch, tmp_path,
                                  ["--params", str(cfg), "--peak-source", "independent"])
        assert params["processing_options"]["rerun_adaptive_per_spectrum"] is True
