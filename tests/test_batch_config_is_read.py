# ABOUTME: Every setting `lunaNMR batch` accepts changes the run: a config's extensions and progress
# ABOUTME: interval are used, settings nothing reads are refused, and each default is declared once.

"""
A config file's `file_handling.extensions` was replaced by the --extensions default on
every run; `logging.progress_interval` was stored where nothing looked; nine more keys
(retry counts, optimiser limits, a log file name, summary switches) were accepted and
read by nothing, and the `high_throughput` preset was made of them. BatchProcessor also
declared its own quality threshold and recursion default beside ConfigManager's.
"""

import contextlib
import io
import json
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

needs_spectra = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


def _batch(*argv):
    from lunaNMR.cli import main
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        code = main(["batch", *argv, "--format", "json"])
    return code, json.loads(out.getvalue())


def _config(tmp_path, content):
    path = tmp_path / "config.json"
    path.write_text(json.dumps(content))
    return str(path)


@pytest.fixture
def folder(tmp_path):
    """One real spectrum, named `.ft`."""
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    (spectra / "2D_DNAJA1_0o0.ft").symlink_to(DNAJA1 / "2D_DNAJA1_0o0.ft")
    return spectra


@needs_spectra
def test_the_extensions_of_a_config_file_are_the_ones_searched(folder, tmp_path):
    config = _config(tmp_path, {"file_handling": {"extensions": [".ft2"]}})
    code, envelope = _batch(str(folder), "--dry-run", "--output", str(tmp_path / "out"),
                            "--config", config)
    assert envelope["spectra_found"] == 0
    assert code == 1


@needs_spectra
def test_extensions_on_the_command_line_replace_the_config_files(folder, tmp_path):
    config = _config(tmp_path, {"file_handling": {"extensions": [".ft2"]}})
    code, envelope = _batch(str(folder), "--dry-run", "--output", str(tmp_path / "out"),
                            "--config", config, "--extensions", ".ft")
    assert envelope["spectra_found"] == 1
    assert code == 0


@pytest.mark.parametrize("content, named", [
    ({"processing": {"max_fitting_attempts": 1}}, "processing.max_fitting_attempts"),
    ({"processing": {"max_optimization_iterations": 10}}, "processing.max_optimization_iterations"),
    ({"processing": {"baseline_max_iter": 10}}, "processing.baseline_max_iter"),
    ({"processing": {"convergence_tolerance": 1e-6}}, "processing.convergence_tolerance"),
    ({"processing": {"max_iterations": 500}}, "processing.max_iterations"),
    ({"processing": {"max_iter": 500}}, "processing.max_iter"),
    ({"processing_parameters": {"fitting_window_x": 0.05}}, "processing_parameters.fitting_window_x"),
    ({"processing_parameters": {"fitting_window_y": 1.0}}, "processing_parameters.fitting_window_y"),
    ({"processing_parameters": {"serch_window_x": 0.05}}, "processing_parameters.serch_window_x"),
    ({"logging": {"detailed_logging": False}}, "logging.detailed_logging"),
    ({"logging": {"log_file": "run.log"}}, "logging.log_file"),
    ({"output": {"generate_summary": False}}, "output"),
    ({"output": {"summary_format": "json"}}, "output"),
    ({"output": {"backup_results": True}}, "output"),
])
def test_a_setting_nothing_reads_is_refused(tmp_path, content, named):
    code, envelope = _batch(str(tmp_path), "--dry-run", "--output", str(tmp_path / "out"),
                            "--config", _config(tmp_path, content))
    assert code == 2
    assert envelope["error"]["type"] == "UsageError"
    assert named in envelope["error"]["message"]


def test_the_search_windows_of_a_config_are_read(tmp_path):
    from lunaNMR.batch_processing.batch_processor import BatchProcessor
    from lunaNMR.batch_processing.config_manager import ConfigManager
    config = ConfigManager().load_config(_config(tmp_path, {
        "processing_parameters": {"search_window_x": 0.05, "search_window_y": 0.2}}))
    with contextlib.redirect_stdout(io.StringIO()):
        processor = BatchProcessor(config, output_dir=tmp_path / "out")
    assert (processor.integrator.search_window_x, processor.integrator.search_window_y) == (0.05, 0.2)
    manager = processor._fit_all_parameter_manager('15N')
    assert (manager.current_params['search_window_x'],
            manager.current_params['search_window_y']) == (0.05, 0.2)


def test_a_processor_holds_no_iteration_limit_or_fitting_window(tmp_path):
    """The fit's iteration limit and its region are the PS2D configuration's."""
    from lunaNMR.batch_processing.batch_processor import BatchProcessor
    out = tmp_path / "out"
    with contextlib.redirect_stdout(io.StringIO()):
        processor = BatchProcessor(output_dir=out)
    assert 'max_iterations' not in processor.config
    log = next(out.glob("batch_processing_*.log")).read_text()
    assert "max_iterations" not in log and "fitting window" not in log


def test_no_configuration_batch_builds_carries_a_setting_nothing_reads():
    from lunaNMR.batch_processing.config_manager import ConfigManager
    manager = ConfigManager()
    configs = {'default': manager.default_config} | {
        name: manager.get_preset_config(name) for name in manager.presets}
    for name, config in configs.items():
        assert sorted(set(config) - {'description'}) == [
            'file_handling', 'logging', 'processing'], name
        assert sorted(config['logging']) == ['level', 'progress_interval'], name
        assert sorted(config['processing']) == [
            'expected_peaks', 'quality_threshold', 'skip_on_error', 'sn_thresholds'], name


def test_there_is_no_preset_that_changes_no_run(tmp_path):
    from lunaNMR.batch_processing.config_manager import ConfigManager
    assert 'high_throughput' not in ConfigManager().presets
    code, envelope = _batch(str(tmp_path), "--preset", "high_throughput", "--dry-run")
    assert code == 2
    assert envelope["error"]["type"] == "UsageError"


@needs_spectra
def test_the_progress_interval_of_a_config_sets_how_often_progress_is_logged(tmp_path):
    """Three spectra no peak is detected in (a threshold nothing reaches): progress is
    logged before the first and the third with an interval of 2."""
    from lunaNMR.batch_processing.batch_processor import BatchProcessor
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    for name in ("2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o5.ft", "2D_DNAJA1_HSPA8_1o0.ft"):
        (spectra / name).symlink_to(DNAJA1 / name)
    config = {"processing": {"sn_thresholds": {"15N1H": 1000.0}},
              "logging": {"progress_interval": 2}}
    out = tmp_path / "out"
    with contextlib.redirect_stdout(io.StringIO()):
        BatchProcessor(config, output_dir=out).process_folder(spectra)
    log = next(out.glob("batch_processing_*.log")).read_text()
    assert log.count("Progress:") == 2


def test_a_processor_without_a_config_runs_with_the_command_lines_defaults(tmp_path):
    from lunaNMR.batch_processing.batch_processor import BatchProcessor
    from lunaNMR.batch_processing.config_manager import ConfigManager
    defaults = ConfigManager().default_config
    with contextlib.redirect_stdout(io.StringIO()):
        processor = BatchProcessor(output_dir=tmp_path / "out")
    assert processor.config['quality_threshold'] == defaults['processing']['quality_threshold'] == 0.8
    assert processor.config['file_handling'] == defaults['file_handling']
    assert defaults['file_handling']['recursive_search'] is False


def test_a_processor_without_a_config_does_not_search_subfolders(tmp_path):
    from lunaNMR.batch_processing.batch_processor import BatchProcessor
    nested = tmp_path / "spectra" / "nested"
    nested.mkdir(parents=True)
    (nested / "spectrum.ft").write_bytes(b"")
    with contextlib.redirect_stdout(io.StringIO()):
        processor = BatchProcessor(output_dir=tmp_path / "out")
        assert processor.find_spectrum_files(tmp_path / "spectra") == []
