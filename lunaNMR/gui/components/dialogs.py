"""ABOUTME: Dialog components for lunaNMR Qt interface
ABOUTME: Provides the error message dialog
"""

from typing import Optional
import logging

from PySide6.QtWidgets import (
    QMessageBox
)

logger = logging.getLogger(__name__)


# SIMPLE MESSAGE DIALOGS (QMessageBox wrappers)

def show_error(
    parent=None,
    title: str = "Error",
    message: str = "",
    detailed_text: Optional[str] = None
) -> None:
    """Show an error message dialog.

    This is a direct replacement for tkinter.messagebox.showerror().

    Args:
        parent: Parent widget or None
        title: Dialog title
        message: Error message to display
        detailed_text: Optional detailed error information (collapsible)

    Example:
        ```python
        try:
            risky_operation()
        except Exception as e:
            show_error(
                self,
                "Operation Failed",
                "Could not complete operation",
                detailed_text=str(e)
            )
        ```
    """
    msg_box = QMessageBox(parent)
    msg_box.setWindowTitle(title)
    msg_box.setText(message)
    msg_box.setIcon(QMessageBox.Critical)
    msg_box.setStandardButtons(QMessageBox.Ok)

    if detailed_text:
        msg_box.setDetailedText(detailed_text)

    msg_box.exec()
