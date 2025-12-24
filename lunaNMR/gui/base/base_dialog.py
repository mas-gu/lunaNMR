"""ABOUTME: Base dialog class for lunaNMR Qt interface providing common dialog functionality
ABOUTME: Automatically loads stylesheets and provides helpers for OK/Cancel buttons
"""

from pathlib import Path
from typing import Optional
import logging

from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton
)
from PySide6.QtCore import Qt, QSize
from PySide6.QtGui import QIcon

logger = logging.getLogger(__name__)


class BaseDialog(QDialog):
    """Base class for dialogs in lunaNMR.

    Provides automatic stylesheet loading, common dialog setup, and helper methods
    for standard OK/Cancel button boxes and dialog configuration.

    Features:
        - Automatic QSS stylesheet loading from lunaNMR/gui/styles/main.qss
        - Default modal behavior
        - Standard button box creation (OK/Cancel, Yes/No, etc.)
        - Default sizing with sensible defaults
        - Window icon handling
        - Error handling for missing resources

    Example:
        ```python
        from lunaNMR.gui.base import BaseDialog
        from PySide6.QtWidgets import QLabel, QVBoxLayout

        class SettingsDialog(BaseDialog):
            def __init__(self, parent=None):
                super().__init__(
                    parent=parent,
                    title="Settings",
                    default_size=(400, 300),
                    modal=True
                )

                # Create your dialog content
                layout = QVBoxLayout()
                layout.addWidget(QLabel("Settings content here"))

                # Add standard OK/Cancel buttons
                button_box = self.create_button_box(
                    buttons=QDialogButtonBox.Ok | QDialogButtonBox.Cancel
                )
                layout.addWidget(button_box)

                self.setLayout(layout)

        # Usage
        dialog = SettingsDialog(parent=main_window)
        if dialog.exec():
            # User clicked OK
            print("Settings accepted")
        ```
    """

    def __init__(
        self,
        parent=None,
        title: str = "Dialog",
        default_size: Optional[tuple[int, int]] = None,
        min_size: Optional[tuple[int, int]] = None,
        modal: bool = True,
        icon_path: Optional[str] = None
    ):
        """Initialize the base dialog with common setup.

        Args:
            parent: Parent widget, or None for independent dialog
            title: Dialog title text
            default_size: Default dialog size as (width, height), or None for auto-size
            min_size: Minimum dialog size as (width, height), or None for no minimum
            modal: Whether the dialog should be modal (blocking parent)
            icon_path: Path to dialog icon file, or None for no icon
        """
        super().__init__(parent)

        # Store configuration
        self._title = title
        self._modal = modal

        # Setup dialog properties
        self.setWindowTitle(title)
        self.setModal(modal)

        if default_size:
            self.resize(*default_size)

        if min_size:
            self.setMinimumSize(QSize(*min_size))

        # Setup icon if provided
        if icon_path:
            self._setup_icon(icon_path)

        # Load stylesheet
        self._load_stylesheet()

    def _load_stylesheet(self) -> None:
        """Load the main QSS stylesheet from lunaNMR/gui/styles/main.qss.

        Handles missing stylesheet gracefully by logging a warning and continuing
        without styling.
        """
        try:
            # Determine stylesheet path relative to this file
            base_dir = Path(__file__).parent.parent  # lunaNMR/gui/
            stylesheet_path = base_dir / "styles" / "main.qss"

            if not stylesheet_path.exists():
                logger.warning(
                    f"Stylesheet not found at {stylesheet_path}. "
                    "Dialog will use default Qt styling."
                )
                return

            # Read and apply stylesheet
            with open(stylesheet_path, 'r', encoding='utf-8') as f:
                stylesheet = f.read()

            self.setStyleSheet(stylesheet)
            logger.debug(f"Loaded stylesheet from {stylesheet_path}")

        except Exception as e:
            logger.warning(
                f"Failed to load stylesheet: {e}. "
                "Dialog will use default Qt styling."
            )

    def _setup_icon(self, icon_path: str) -> None:
        """Setup dialog icon from file path.

        Args:
            icon_path: Path to icon file
        """
        try:
            icon = QIcon(icon_path)
            if not icon.isNull():
                self.setWindowIcon(icon)
                logger.debug(f"Set dialog icon from {icon_path}")
            else:
                logger.warning(f"Icon file not found or invalid: {icon_path}")
        except Exception as e:
            logger.warning(f"Failed to set dialog icon: {e}")

    def create_button_box(
        self,
        buttons: QDialogButtonBox.StandardButton = (
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        ),
        connect_default: bool = True
    ) -> QDialogButtonBox:
        """Create a standard dialog button box.

        Args:
            buttons: Standard buttons to include (e.g., Ok | Cancel)
            connect_default: Whether to connect buttons to accept()/reject() slots

        Returns:
            QDialogButtonBox instance ready to add to layout

        Example:
            ```python
            # OK and Cancel buttons (default)
            button_box = self.create_button_box()

            # Yes and No buttons
            button_box = self.create_button_box(
                buttons=QDialogButtonBox.Yes | QDialogButtonBox.No
            )

            # Custom handling
            button_box = self.create_button_box(connect_default=False)
            button_box.accepted.connect(self.custom_accept_handler)
            button_box.rejected.connect(self.custom_reject_handler)
            ```
        """
        button_box = QDialogButtonBox(buttons)

        if connect_default:
            button_box.accepted.connect(self.accept)
            button_box.rejected.connect(self.reject)

        # Apply primary button style to accept buttons
        accept_button = button_box.button(QDialogButtonBox.Ok)
        if accept_button:
            accept_button.setProperty("class", "primary")

        yes_button = button_box.button(QDialogButtonBox.Yes)
        if yes_button:
            yes_button.setProperty("class", "primary")

        apply_button = button_box.button(QDialogButtonBox.Apply)
        if apply_button:
            apply_button.setProperty("class", "primary")

        return button_box

    def create_custom_buttons(
        self,
        buttons: list[tuple[str, callable]],
        add_cancel: bool = True
    ) -> QHBoxLayout:
        """Create a custom button layout with specified buttons.

        Args:
            buttons: List of (label, callback) tuples for custom buttons
            add_cancel: Whether to automatically add a Cancel button

        Returns:
            QHBoxLayout containing the buttons

        Example:
            ```python
            buttons = [
                ("Apply", self.apply_settings),
                ("Reset", self.reset_to_defaults),
            ]
            button_layout = self.create_custom_buttons(buttons, add_cancel=True)
            main_layout.addLayout(button_layout)
            ```
        """
        layout = QHBoxLayout()
        layout.addStretch()

        # Add custom buttons
        for label, callback in buttons:
            button = QPushButton(label)
            button.clicked.connect(callback)
            layout.addWidget(button)

        # Add cancel button if requested
        if add_cancel:
            cancel_button = QPushButton("Cancel")
            cancel_button.clicked.connect(self.reject)
            layout.addWidget(cancel_button)

        return layout

    def center_on_parent(self) -> None:
        """Center the dialog on its parent window.

        If no parent exists, centers on screen instead.

        Example:
            ```python
            dialog = BaseDialog(parent=main_window)
            dialog.center_on_parent()
            dialog.show()
            ```
        """
        if self.parent():
            parent_geometry = self.parent().frameGeometry()
            dialog_geometry = self.frameGeometry()
            center_point = parent_geometry.center()
            dialog_geometry.moveCenter(center_point)
            self.move(dialog_geometry.topLeft())
        else:
            self.center_on_screen()

    def center_on_screen(self) -> None:
        """Center the dialog on the screen.

        Example:
            ```python
            dialog = BaseDialog()
            dialog.center_on_screen()
            dialog.show()
            ```
        """
        frame_geometry = self.frameGeometry()
        screen_center = self.screen().availableGeometry().center()
        frame_geometry.moveCenter(screen_center)
        self.move(frame_geometry.topLeft())

    def get_stylesheet_path(self) -> Path:
        """Get the path to the main stylesheet file.

        Returns:
            Path object pointing to main.qss

        Note:
            This is useful for subclasses that need to load additional stylesheets
            or verify stylesheet location.
        """
        base_dir = Path(__file__).parent.parent
        return base_dir / "styles" / "main.qss"

    @staticmethod
    def create_simple_layout(spacing: int = 10) -> QVBoxLayout:
        """Create a simple vertical layout with standard spacing.

        Args:
            spacing: Spacing between widgets in pixels

        Returns:
            QVBoxLayout with configured spacing

        Example:
            ```python
            layout = BaseDialog.create_simple_layout()
            layout.addWidget(QLabel("Content"))
            button_box = self.create_button_box()
            layout.addWidget(button_box)
            self.setLayout(layout)
            ```
        """
        layout = QVBoxLayout()
        layout.setSpacing(spacing)
        return layout
