# ABOUTME: Standalone Kd / titration module dialog, launched from the lunaNMR Modules menu.
# ABOUTME: Hosts the KdTitrationPage as an independent module (sibling to DynamiXs).

import os
import sys
from pathlib import Path

from PySide6.QtWidgets import QVBoxLayout

from lunaNMR.gui.base.base_dialog import BaseDialog

# Reuse the dynamiXs module's page/worker/fitter; add it to the path for imports.
_dynamixs_path = Path(__file__).parent.parent.parent.parent / "modules" / "dynamiXs_v2o0"
if str(_dynamixs_path) not in sys.path:
    sys.path.insert(0, str(_dynamixs_path))

from kd_titration_page import KdTitrationPage


class KdTitrationDialog(BaseDialog):
    """Independent Kd / titration analysis module embedded in lunaNMR.

    Hosts the single KdTitrationPage. The page's 'back' button (BasePage calls
    main_window.show_main_menu) closes the module, returning to lunaNMR.
    """

    def __init__(self, parent=None, main_window=None):
        super().__init__(
            parent=parent,
            title="Kd / Titration Analysis",
            default_size=(960, 760),
            min_size=(820, 640),
            modal=False,
        )
        self.main_window = main_window
        # current_dir is read by the page's file pickers.
        self.current_dir = getattr(main_window, 'current_nmr_folder', None) or os.getcwd()

        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        self.kd_page = KdTitrationPage(self)
        layout.addWidget(self.kd_page)
        self.setLayout(layout)

        # Restore state saved from a previous session (set when the dialog closed
        # or when a project was loaded).
        if main_window is not None:
            if getattr(main_window, 'kd_state', None):
                self.set_state(main_window.kd_state)
            if getattr(main_window, 'kd_file_refs', None):
                self.set_file_refs(main_window.kd_file_refs)

    def show_main_menu(self):
        """BasePage back-button target: this module has one page, so close it."""
        self.close()

    # -------------------------------------------------------------------------
    # State management (for project save/load)
    # -------------------------------------------------------------------------

    def get_state(self) -> dict:
        """Collect Kd page state for project save."""
        return self.kd_page.get_session_state()

    def set_state(self, state: dict):
        """Restore Kd page state from project load."""
        if state:
            self.kd_page.restore_session_state(state)

    def get_file_refs(self) -> dict:
        """Collect Kd input file reference."""
        refs = {}
        if getattr(self.kd_page, 'input_file', None):
            refs['input_file'] = self.kd_page.input_file
        return refs

    def set_file_refs(self, refs: dict):
        """Restore Kd input file reference."""
        if refs and 'input_file' in refs:
            self.kd_page.input_file = refs['input_file']
            if hasattr(self.kd_page, 'file_drop'):
                self.kd_page.file_drop.setText(os.path.basename(refs['input_file']))

    def closeEvent(self, event):
        """Stop a running fit before the page is destroyed, so the worker thread
        can't emit signals into already-deleted widgets (use-after-free)."""
        worker = getattr(self.kd_page, 'worker', None)
        if worker is not None and worker.isRunning():
            worker.cancel()
            worker.wait(5000)
        # Transfer state to the main window so a later project save can capture it.
        if self.main_window is not None:
            self.main_window.kd_state = self.get_state()
            self.main_window.kd_file_refs = self.get_file_refs()
        super().closeEvent(event)
