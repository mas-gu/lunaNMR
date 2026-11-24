# ABOUTME: Standalone 1D integration module dialog, launched from the lunaNMR Modules menu.
# ABOUTME: Hosts the OneDIntegrationPage - pick peaks on a 1D spectrum, integrate the whole series.

from PySide6.QtWidgets import QVBoxLayout

from lunaNMR.analysis_modules import add_to_import_path
from lunaNMR.gui.base.base_dialog import BaseDialog

# The 1D module is self-contained; add it to the path for its bare-name imports.
add_to_import_path("integration_1d_v1o0")

from oned_page import OneDIntegrationPage


class Integration1DDialog(BaseDialog):
    """Independent 1D peak integration module embedded in lunaNMR."""

    def __init__(self, parent=None, main_window=None):
        super().__init__(
            parent=parent,
            title="1D Integration",
            default_size=(1180, 760),
            min_size=(900, 600),
            modal=False,
        )
        self.main_window = main_window

        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)

        self.page = OneDIntegrationPage(parent=self)
        layout.addWidget(self.page)

        self.setLayout(layout)

    def show_main_menu(self):
        """The page's back action closes the module, returning to lunaNMR."""
        self.close()
