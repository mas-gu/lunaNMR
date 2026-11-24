# ABOUTME: A CLI series and a GUI series start with the same cluster-locking choice.
# ABOUTME: Compares the CLI's series defaults with a real main window's setting.

import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_a_cli_series_locks_clusters_as_the_main_window_does():
    from PySide6.QtWidgets import QApplication
    from lunaNMR.cli import _default_series_params
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    app = QApplication.instance() or QApplication(sys.argv)
    set_ps2d_config('15N')
    window = LunaNMRMainWindow()
    try:
        cli = _default_series_params()['processing_options']['lock_cluster_assignments']
        assert cli is window.lock_cluster_assignments is True
    finally:
        window.close()
        window.deleteLater()
        app.processEvents()
        set_ps2d_config('15N')
