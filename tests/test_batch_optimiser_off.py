# ABOUTME: batch's threshold optimiser runs only when asked for (--optimize or a config that enables it).
# ABOUTME: While the detected list is capped every candidate threshold yields the same peaks, so by default it is off.


def test_the_default_config_does_not_optimise():
    from lunaNMR.batch_processing.config_manager import ConfigManager
    assert ConfigManager().default_config['optimization']['enable_auto_optimization'] is False


def test_every_preset_leaves_it_off():
    from lunaNMR.batch_processing.config_manager import ConfigManager
    manager = ConfigManager()
    for name in manager.presets:
        assert manager.get_preset_config(name)['optimization']['enable_auto_optimization'] is False, name
