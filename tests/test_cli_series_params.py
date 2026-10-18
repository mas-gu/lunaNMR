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
    """docs/CLI.md advises tuning the detection and fitting knobs, and the CLI froze them at
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
            "fitting_params": {"fitting_window_x": 0.3},
            "gui_params": {"fix_linewidths": True}}))
        params = _captured_params(monkeypatch, tmp_path, ["--params", str(cfg)])
        assert params["detection_params"]["search_window_y"] == 0.04
        assert params["fitting_params"]["fitting_window_x"] == 0.3
        assert params["gui_params"]["fix_linewidths"] is True


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


class TestAWrongTypeIsRefused:
    """A knob is checked for its type as well as its name. "false" is a non-empty string,
    so the engine read it as True: position fixing switched ON, every Stage 2 skipped,
    while the summary echoed the value the caller believed they had set.
    """

    @pytest.mark.parametrize("payload", [
        {"gui_params": {"fix_positions": "false"}},
        {"gui_params": {"use_ps2d_linewidth_reuse": 0}},
        {"detection_params": {"search_window_x": "0.05"}},
        {"fitting_params": {"max_iterations": None}},
        {"fitting_params": {"max_iterations": True}},
    ])
    def test_a_value_of_the_wrong_type_is_an_error(self, tmp_path, capsys, payload):
        assert TestAnUnknownKeyIsRefused()._err(tmp_path, payload) == 1
        (section, values), = payload.items()
        (key, _), = values.items()
        assert key in capsys.readouterr().err

    def test_the_dry_run_refuses_it_too(self, tmp_path):
        cfg = tmp_path / "tune.json"
        cfg.write_text(json.dumps({"gui_params": {"fix_positions": "false"}}))
        peaks = _inputs(tmp_path)
        from lunaNMR.cli import main
        assert main(["series", "--spectra", str(tmp_path), "--peaks", str(peaks),
                     "--out", str(tmp_path / "o"), "--params", str(cfg), "--dry-run"]) == 1

    def test_an_integer_is_accepted_for_a_real_valued_knob(self):
        from lunaNMR.cli import _default_series_params, _merge_series_params
        merged = _merge_series_params(_default_series_params(),
                                      {"fitting_params": {"fitting_window_y": 2}}, "tune.json")
        assert merged["fitting_params"]["fitting_window_y"] == 2


class TestCommandLineFlagsWin:
    """--params is a file; --parallel and --peak-source are what the caller just typed.
    The typed thing wins, or a stale file silently overrides the command line."""

    def test_parallel_is_the_flags_alone(self, monkeypatch, tmp_path):
        params = _captured_params(monkeypatch, tmp_path, ["--parallel"])
        assert params["processing_options"]["use_parallel_processing"] is True
        assert params["gui_params"]["use_parallel_processing"] is True

    def test_independent_mode_beats_the_file(self, monkeypatch, tmp_path):
        cfg = tmp_path / "tune.json"
        cfg.write_text(json.dumps(
            {"processing_options": {"rerun_adaptive_per_spectrum": False}}))
        params = _captured_params(monkeypatch, tmp_path,
                                  ["--params", str(cfg), "--peak-source", "independent"])
        assert params["processing_options"]["rerun_adaptive_per_spectrum"] is True


class TestAKnobWithNoEffectIsRefused:
    """A key the series engine never reads is accepted, echoed in the summary as the
    setting the run used, and changes nothing. Each one below was traced to no reader on
    the series path, and measured: a 2-point DNAJA1 run wrote a byte-identical tidy CSV
    with it changed, in reference, cascade, independent and titration modes.
    """

    @pytest.mark.parametrize("section,key,value", [
        ("gui_params", "use_centroid_refinement", False),
        ("gui_params", "centroid_noise_multiplier", 4.5),
        ("gui_params", "use_ps2d_multi_peak", False),
        ("gui_params", "height_threshold", 0.5),
        ("gui_params", "distance_factor", 50.0),
        ("gui_params", "prominence_threshold", 0.3),
        ("gui_params", "smoothing_sigma", 3.0),
        ("gui_params", "max_peaks_fit", 2),
        ("gui_params", "max_optimization_iterations", 2),
        ("fitting_params", "min_r_squared", 0.95),
        ("processing_options", "use_global_optimization", True),
    ])
    def test_it_is_refused_as_having_no_effect(self, tmp_path, capsys, section, key, value):
        assert TestAnUnknownKeyIsRefused()._err(tmp_path, {section: {key: value}}) == 1
        err = capsys.readouterr().err
        assert key in err and "no effect" in err

    @pytest.mark.parametrize("section", ["gui_params", "processing_options"])
    def test_parallel_processing_is_refused_in_favour_of_the_flag(self, tmp_path, capsys,
                                                                  section):
        payload = {section: {"use_parallel_processing": True}}
        assert TestAnUnknownKeyIsRefused()._err(tmp_path, payload) == 1
        err = capsys.readouterr().err
        assert "no effect" in err and "--parallel" in err

    def test_the_dry_run_refuses_it_too(self, tmp_path):
        cfg = tmp_path / "tune.json"
        cfg.write_text(json.dumps({"fitting_params": {"min_r_squared": 0.9}}))
        peaks = _inputs(tmp_path)
        from lunaNMR.cli import main
        assert main(["series", "--spectra", str(tmp_path), "--peaks", str(peaks),
                     "--out", str(tmp_path / "o"), "--params", str(cfg), "--dry-run"]) == 1

    def test_the_known_keys_it_lists_are_the_accepted_ones(self, tmp_path, capsys):
        TestAnUnknownKeyIsRefused()._err(tmp_path, {"gui_params": {"no_such_knob": 1}})
        listed = capsys.readouterr().err.split("Known keys:")[1]
        assert "fix_positions" in listed
        assert "use_centroid_refinement" not in listed

    def test_noise_threshold_is_still_accepted(self, monkeypatch, tmp_path):
        cfg = tmp_path / "tune.json"
        cfg.write_text(json.dumps({"detection_params": {"noise_threshold": 5.0}}))
        params = _captured_params(monkeypatch, tmp_path, ["--params", str(cfg)])
        assert params["detection_params"]["noise_threshold"] == 5.0


class TestTitrationRefusesTheFreeStateAnchor:
    """use_original_reference_for_detection is what --peak-source independent sets; in a
    titration it detects every point against the free-state list, the combination the
    CLI already refuses as --peak-source independent."""

    def test_it_is_refused_with_titration(self, tmp_path, capsys):
        cfg = tmp_path / "tune.json"
        cfg.write_text(json.dumps(
            {"processing_options": {"use_original_reference_for_detection": True}}))
        peaks = _inputs(tmp_path)
        from lunaNMR.cli import main
        assert main(["series", "--spectra", str(tmp_path), "--peaks", str(peaks),
                     "--out", str(tmp_path / "o"), "--params", str(cfg),
                     "--mode", "titration", "--dry-run"]) == 1
        assert "cascade" in capsys.readouterr().err

    def test_it_is_accepted_in_time_mode(self, monkeypatch, tmp_path):
        cfg = tmp_path / "tune.json"
        cfg.write_text(json.dumps(
            {"processing_options": {"use_original_reference_for_detection": True}}))
        params = _captured_params(monkeypatch, tmp_path,
                                  ["--params", str(cfg), "--peak-source", "cascade"])
        assert params["processing_options"]["use_original_reference_for_detection"] is True
