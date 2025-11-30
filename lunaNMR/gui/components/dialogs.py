"""ABOUTME: Dialog components for lunaNMR Qt interface
ABOUTME: Provides message dialogs (error, warning, info, question) and custom confirmation dialogs
"""

from typing import Optional
import logging

from PySide6.QtWidgets import (
    QMessageBox,
    QVBoxLayout,
    QLabel,
    QCheckBox,
    QDialogButtonBox
)
from PySide6.QtCore import Qt

from lunaNMR.gui.base import BaseDialog
from lunaNMR.gui.styles.design_system import (
    SPACING_MD,
    SPACING_SM,
    FONT_SIZE_BODY
)

logger = logging.getLogger(__name__)


# ============================================================================
# SIMPLE MESSAGE DIALOGS (QMessageBox wrappers)
# ============================================================================

def show_question(
    parent=None,
    title: str = "Question",
    message: str = "",
    default_yes: bool = True
) -> bool:
    """Show a yes/no question dialog.

    This is a direct replacement for tkinter.messagebox.askyesno().

    Args:
        parent: Parent widget or None
        title: Dialog title
        message: Question message to display
        default_yes: If True, Yes button is default; if False, No is default

    Returns:
        True if user clicked Yes, False if user clicked No

    Example:
        ```python
        if show_question(self, "Confirm", "Delete all data?", default_yes=False):
            # User clicked Yes
            delete_data()
        ```
    """
    msg_box = QMessageBox(parent)
    msg_box.setWindowTitle(title)
    msg_box.setText(message)
    msg_box.setIcon(QMessageBox.Question)
    msg_box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)

    if default_yes:
        msg_box.setDefaultButton(QMessageBox.Yes)
    else:
        msg_box.setDefaultButton(QMessageBox.No)

    result = msg_box.exec()
    return result == QMessageBox.Yes


def show_confirmation(
    parent=None,
    title: str = "Confirm",
    message: str = "",
    ok_text: str = "OK",
    cancel_text: str = "Cancel"
) -> bool:
    """Show an OK/Cancel confirmation dialog.

    This is a direct replacement for tkinter.messagebox.askokcancel().

    Args:
        parent: Parent widget or None
        title: Dialog title
        message: Confirmation message to display
        ok_text: Text for OK button (default: "OK")
        cancel_text: Text for Cancel button (default: "Cancel")

    Returns:
        True if user clicked OK, False if user clicked Cancel

    Example:
        ```python
        if show_confirmation(self, "Save", "Save changes before closing?"):
            # User clicked OK
            save_changes()
        ```
    """
    msg_box = QMessageBox(parent)
    msg_box.setWindowTitle(title)
    msg_box.setText(message)
    msg_box.setIcon(QMessageBox.Question)

    ok_button = msg_box.addButton(ok_text, QMessageBox.AcceptRole)
    cancel_button = msg_box.addButton(cancel_text, QMessageBox.RejectRole)
    msg_box.setDefaultButton(ok_button)

    # Apply primary button style to OK button
    ok_button.setProperty("class", "primary")

    msg_box.exec()
    return msg_box.clickedButton() == ok_button


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


def show_warning(
    parent=None,
    title: str = "Warning",
    message: str = "",
    detailed_text: Optional[str] = None
) -> None:
    """Show a warning message dialog.

    This is a direct replacement for tkinter.messagebox.showwarning().

    Args:
        parent: Parent widget or None
        title: Dialog title
        message: Warning message to display
        detailed_text: Optional detailed warning information (collapsible)

    Example:
        ```python
        if file_count > 100:
            show_warning(
                self,
                "Large Dataset",
                "Processing 100+ files may take several minutes"
            )
        ```
    """
    msg_box = QMessageBox(parent)
    msg_box.setWindowTitle(title)
    msg_box.setText(message)
    msg_box.setIcon(QMessageBox.Warning)
    msg_box.setStandardButtons(QMessageBox.Ok)

    if detailed_text:
        msg_box.setDetailedText(detailed_text)

    msg_box.exec()


def show_info(
    parent=None,
    title: str = "Information",
    message: str = "",
    detailed_text: Optional[str] = None
) -> None:
    """Show an information message dialog.

    This is a direct replacement for tkinter.messagebox.showinfo().

    Args:
        parent: Parent widget or None
        title: Dialog title
        message: Information message to display
        detailed_text: Optional detailed information (collapsible)

    Example:
        ```python
        show_info(
            self,
            "Processing Complete",
            f"Successfully processed {count} spectra in {elapsed} seconds"
        )
        ```
    """
    msg_box = QMessageBox(parent)
    msg_box.setWindowTitle(title)
    msg_box.setText(message)
    msg_box.setIcon(QMessageBox.Information)
    msg_box.setStandardButtons(QMessageBox.Ok)

    if detailed_text:
        msg_box.setDetailedText(detailed_text)

    msg_box.exec()


# ============================================================================
# CUSTOM DIALOGS (BaseDialog subclasses)
# ============================================================================

class ConfirmationDialog(BaseDialog):
    """Custom confirmation dialog with optional checkbox.

    This provides more flexibility than the simple show_confirmation() function,
    allowing for additional options like "Don't ask again" checkboxes.

    Features:
        - Custom message and icon
        - Optional checkbox (e.g., "Don't show this again")
        - Customizable button text
        - Styled with design system

    Example:
        ```python
        dialog = ConfirmationDialog(
            parent=self,
            title="Clear All Data",
            message="This will delete all peaks and results. Continue?",
            checkbox_text="Don't ask me again",
            ok_text="Delete",
            cancel_text="Cancel",
            destructive=True
        )

        if dialog.exec():
            # User clicked OK
            clear_all_data()

            if dialog.checkbox_checked():
                # User checked "Don't ask again"
                save_preference("skip_clear_warning", True)
        ```
    """

    def __init__(
        self,
        parent=None,
        title: str = "Confirm",
        message: str = "",
        checkbox_text: Optional[str] = None,
        ok_text: str = "OK",
        cancel_text: str = "Cancel",
        destructive: bool = False,
        icon: Optional[QMessageBox.Icon] = None
    ):
        """Initialize the confirmation dialog.

        Args:
            parent: Parent widget or None
            title: Dialog title
            message: Confirmation message to display
            checkbox_text: Optional checkbox label (e.g., "Don't ask again")
            ok_text: Text for OK/accept button
            cancel_text: Text for Cancel/reject button
            destructive: If True, OK button gets destructive styling (red)
            icon: Optional QMessageBox icon (Question, Warning, etc.)
        """
        super().__init__(
            parent=parent,
            title=title,
            default_size=(400, 200),
            modal=True
        )

        self._checkbox = None
        self._destructive = destructive

        # Create layout
        layout = self.create_simple_layout(spacing=SPACING_MD)

        # Icon and message
        message_label = QLabel(message)
        message_label.setWordWrap(True)
        message_label.setAlignment(Qt.AlignLeft | Qt.AlignTop)

        # Set font size from design system
        font = message_label.font()
        font.setPointSize(FONT_SIZE_BODY)
        message_label.setFont(font)

        layout.addWidget(message_label)

        # Optional checkbox
        if checkbox_text:
            self._checkbox = QCheckBox(checkbox_text)
            layout.addWidget(self._checkbox)

        layout.addStretch()

        # Custom buttons
        button_layout = self.create_custom_buttons(
            buttons=[(ok_text, self._on_ok_clicked)],
            add_cancel=True
        )

        # Apply destructive style to OK button if requested
        if destructive:
            ok_button = button_layout.itemAt(1).widget()  # First button after stretch
            if ok_button:
                ok_button.setProperty("class", "destructive")
                ok_button.style().polish(ok_button)  # Reapply stylesheet
        else:
            # Apply primary style to OK button
            ok_button = button_layout.itemAt(1).widget()
            if ok_button:
                ok_button.setProperty("class", "primary")
                ok_button.style().polish(ok_button)

        layout.addLayout(button_layout)

        self.setLayout(layout)
        self.center_on_parent()

    def _on_ok_clicked(self):
        """Handle OK button click."""
        self.accept()

    def checkbox_checked(self) -> bool:
        """Check if the optional checkbox was checked.

        Returns:
            True if checkbox exists and is checked, False otherwise
        """
        if self._checkbox:
            return self._checkbox.isChecked()
        return False


class InputDialog(BaseDialog):
    """Custom input dialog for getting text from user.

    This is a more flexible alternative to QInputDialog with design system styling.

    Example:
        ```python
        dialog = InputDialog(
            parent=self,
            title="Peak Label",
            message="Enter peak assignment:",
            default_value="H1",
            placeholder="e.g., H1, H2, etc."
        )

        if dialog.exec():
            assignment = dialog.get_value()
            update_peak_assignment(assignment)
        ```
    """

    def __init__(
        self,
        parent=None,
        title: str = "Input",
        message: str = "",
        default_value: str = "",
        placeholder: str = "",
        ok_text: str = "OK",
        cancel_text: str = "Cancel"
    ):
        """Initialize the input dialog.

        Args:
            parent: Parent widget or None
            title: Dialog title
            message: Prompt message
            default_value: Default text in input field
            placeholder: Placeholder text when field is empty
            ok_text: Text for OK button
            cancel_text: Text for Cancel button
        """
        from PySide6.QtWidgets import QLineEdit

        super().__init__(
            parent=parent,
            title=title,
            default_size=(400, 150),
            modal=True
        )

        # Create layout
        layout = self.create_simple_layout(spacing=SPACING_MD)

        # Message label
        if message:
            message_label = QLabel(message)
            message_label.setWordWrap(True)
            font = message_label.font()
            font.setPointSize(FONT_SIZE_BODY)
            message_label.setFont(font)
            layout.addWidget(message_label)

        # Input field
        self._input = QLineEdit()
        self._input.setText(default_value)
        self._input.setPlaceholderText(placeholder)
        layout.addWidget(self._input)

        layout.addStretch()

        # Buttons
        button_box = self.create_button_box()
        layout.addWidget(button_box)

        self.setLayout(layout)
        self.center_on_parent()

        # Focus input field
        self._input.setFocus()
        self._input.selectAll()

    def get_value(self) -> str:
        """Get the entered text value.

        Returns:
            The text entered in the input field
        """
        return self._input.text()


# ============================================================================
# CONVENIENCE FUNCTIONS
# ============================================================================

def ask_yes_no(
    parent=None,
    title: str = "Question",
    message: str = "",
    default_yes: bool = True
) -> bool:
    """Alias for show_question() to match CustomTkinter naming.

    Args:
        parent: Parent widget or None
        title: Dialog title
        message: Question message
        default_yes: Default button selection

    Returns:
        True if Yes was clicked, False if No was clicked
    """
    return show_question(parent, title, message, default_yes)


def ask_ok_cancel(
    parent=None,
    title: str = "Confirm",
    message: str = ""
) -> bool:
    """Alias for show_confirmation() to match tkinter naming.

    Args:
        parent: Parent widget or None
        title: Dialog title
        message: Confirmation message

    Returns:
        True if OK was clicked, False if Cancel was clicked
    """
    return show_confirmation(parent, title, message)
