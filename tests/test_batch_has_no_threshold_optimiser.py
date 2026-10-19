# ABOUTME: batch has no threshold optimiser: no --optimize flag, no optimisation config keys.
# ABOUTME: While the detected list is capped every candidate threshold yields the same peaks.

"""
The optimiser tried several S/N thresholds per spectrum and kept the one whose peak
count fell in a target range. batch caps the detected list at the expected peak count,
so every candidate gave the same peaks (DNAJA1 1.8-50, methyl 1.5-10): it cost six
detections per spectrum and chose nothing. `--optimize` is an unknown argument, a
usage error, not a flag accepted and ignored.
"""

import contextlib
import importlib.util
import io
import json


def test_optimize_is_a_usage_error(tmp_path):
    from lunaNMR.cli import main
    stdout, stderr = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        code = main(["batch", str(tmp_path), "--optimize", "--dry-run", "--format", "json"])
    assert code == 2
    assert "unrecognized arguments: --optimize" in stderr.getvalue()
    envelope = json.loads(stdout.getvalue())
    assert envelope["ok"] is False
    assert envelope["error"]["type"] == "UsageError"


def test_no_configuration_carries_optimiser_settings():
    from lunaNMR.batch_processing.config_manager import ConfigManager
    manager = ConfigManager()
    configs = {'default': manager.default_config} | {
        name: manager.get_preset_config(name) for name in manager.presets}
    for name, config in configs.items():
        assert 'optimization' not in config, name
        assert 'auto_optimize' not in config['processing'], name


def test_the_optimiser_module_is_gone():
    assert importlib.util.find_spec("lunaNMR.batch_processing.parameter_optimizer") is None
