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

    def show_main_menu(self):
        """BasePage back-button target: this module has one page, so close it."""
        self.close()
