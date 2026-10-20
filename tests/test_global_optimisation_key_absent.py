# ABOUTME: No parameter set or processing-options dict carries use_global_optimization: nothing reads it.
# ABOUTME: An old project or --params file that has it is ignored or refused, never silently obeyed.

import json
import sys

import pytest

KEY = 'use_global_optimization'


def test_the_parameter_manager_carries_none():
    from lunaNMR.utils.parameter_manager import NMRParameterManager
    manager = NMRParameterManager()
    assert KEY not in manager.parameter_definitions
    assert KEY not in manager.get_integrator_parameters()['processing_options']
    assert KEY not in str(manager.get_parameter_summary())


def test_the_cli_series_defaults_carry_none():
    from lunaNMR.cli import _default_series_params
    assert KEY not in _default_series_params()['processing_options']


def test_a_project_does_not_save_it():
    from lunaNMR.utils.project_manager import ProjectManager
    assert KEY not in ProjectManager.GUI_PARAMS


def test_an_old_series_params_file_with_it_is_refused(tmp_path, capsys):
    from lunaNMR.cli import main
    params = tmp_path / "p.json"
    params.write_text(json.dumps({'processing_options': {KEY: True}}))
    code = main(["series", "--dry-run", "--spectra", str(tmp_path), "--peaks", str(params),
                 "--out", str(tmp_path / "o"), "--params", str(params)])
    assert code != 0
    assert KEY in capsys.readouterr().err


@pytest.fixture
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def test_the_main_window_holds_none(app):
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    window = LunaNMRMainWindow()
    try:
        assert not hasattr(window, KEY)
    finally:
        window.close()
        window.deleteLater()
