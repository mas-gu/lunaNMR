"""ABOUTME: Base window class for lunaNMR Qt interface providing common functionality
ABOUTME: Automatically loads stylesheets and provides helper methods for window setup
"""

from pathlib import Path
from typing import Optional
import logging

from PySide6.QtWidgets import QMainWindow, QStatusBar, QMenuBar
from PySide6.QtCore import QSize
from PySide6.QtGui import QIcon

logger = logging.getLogger(__name__)


class BaseWindow(QMainWindow):
    """Base class for main application windows in lunaNMR.

    Provides automatic stylesheet loading, common window setup, and helper methods
    for status bars, menu bars, and window configuration.

    Features:
        - Automatic QSS stylesheet loading from lunaNMR/gui/styles/main.qss
        - Default window sizing with sensible defaults
        - Status bar setup with helper methods
        - Menu bar creation helpers
        - Window icon handling
        - Error handling for missing resources

    Example:
        ```python
        from lunaNMR.gui.base import BaseWindow
        from PySide6.QtWidgets import QLabel

        class MainWindow(BaseWindow):
            def __init__(self):
                super().__init__(
                    title="lunaNMR - Main Window",
                    default_size=(1200, 800),
                    enable_status_bar=True
                )

                # Your custom initialization
                central_widget = QLabel("Hello lunaNMR")
                self.setCentralWidget(central_widget)

                # Use helper methods
                self.update_status("Application ready")
        ```
    """

    def __init__(
        self,
        title: str = "lunaNMR",
        default_size: tuple[int, int] = (800, 600),
        min_size: Optional[tuple[int, int]] = None,
        enable_status_bar: bool = False,
        icon_path: Optional[str] = None,
        parent=None
    ):
        """Initialize the base window with common setup.

        Args:
            title: Window title text
            default_size: Default window size as (width, height)
            min_size: Minimum window size as (width, height), or None for no minimum
            enable_status_bar: Whether to create and enable a status bar
            icon_path: Path to window icon file, or None for no icon
            parent: Parent widget, or None for top-level window
        """
        super().__init__(parent)

        # Store configuration
        self._title = title
        self._status_bar_enabled = enable_status_bar

        # Setup window properties
        self.setWindowTitle(title)
        self.resize(*default_size)

        if min_size:
            self.setMinimumSize(QSize(*min_size))

        # Setup icon if provided
        if icon_path:
            self._setup_icon(icon_path)

        # Setup status bar if requested
        if enable_status_bar:
            self._setup_status_bar()

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
                    "Window will use default Qt styling."
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
                "Window will use default Qt styling."
            )

    def _setup_icon(self, icon_path: str) -> None:
        """Setup window icon from file path.

        Args:
            icon_path: Path to icon file
        """
        try:
            icon = QIcon(icon_path)
            if not icon.isNull():
                self.setWindowIcon(icon)
                logger.debug(f"Set window icon from {icon_path}")
            else:
                logger.warning(f"Icon file not found or invalid: {icon_path}")
        except Exception as e:
            logger.warning(f"Failed to set window icon: {e}")

    def _setup_status_bar(self) -> None:
        """Create and configure the status bar."""
        status_bar = QStatusBar()
        self.setStatusBar(status_bar)
        logger.debug("Status bar enabled")

    def update_status(self, message: str, timeout: int = 0) -> None:
        """Update status bar message.

        Args:
            message: Status message to display
            timeout: Message timeout in milliseconds (0 = permanent)

        Note:
            If status bar is not enabled, this method does nothing.
        """
        if self._status_bar_enabled and self.statusBar():
            self.statusBar().showMessage(message, timeout)

    def create_menu_bar(self) -> QMenuBar:
        """Create and return the window's menu bar.

        Returns:
            The menu bar instance

        Example:
            ```python
            menu_bar = self.create_menu_bar()
            file_menu = menu_bar.addMenu("&File")
            file_menu.addAction("&Open", self.open_file)
            ```
        """
        menu_bar = self.menuBar()
        return menu_bar

    def center_on_screen(self) -> None:
        """Center the window on the screen.

        Example:
            ```python
            window = BaseWindow()
            window.center_on_screen()
            window.show()
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
