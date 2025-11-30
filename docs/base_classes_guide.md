# lunaNMR Base Classes Guide

## Overview

The `lunaNMR.gui.base` module provides foundational classes for building Qt-based GUI components in lunaNMR v1.0. These classes handle common functionality like stylesheet loading, window setup, and standard dialog patterns.

## Installation

Ensure PySide6 is installed:

```bash
pip install -r requirements.txt
```

## Base Classes

### BaseWindow

Main application window base class that provides:

- Automatic QSS stylesheet loading from `lunaNMR/gui/styles/main.qss`
- Default window sizing and positioning
- Status bar management
- Menu bar creation helpers
- Window icon handling
- Error handling for missing resources

#### Basic Usage

```python
from lunaNMR.gui.base import BaseWindow
from PySide6.QtWidgets import QLabel

class MainWindow(BaseWindow):
    def __init__(self):
        super().__init__(
            title="lunaNMR - Main Window",
            default_size=(1200, 800),
            min_size=(800, 600),
            enable_status_bar=True
        )

        # Your custom initialization
        central_widget = QLabel("Hello lunaNMR")
        self.setCentralWidget(central_widget)

        # Use helper methods
        self.update_status("Application ready")
        self.center_on_screen()
```

#### Constructor Parameters

- `title` (str): Window title text (default: "lunaNMR")
- `default_size` (tuple): Default window size as (width, height) (default: (800, 600))
- `min_size` (tuple | None): Minimum window size or None
- `enable_status_bar` (bool): Whether to create a status bar (default: False)
- `icon_path` (str | None): Path to window icon file
- `parent` (QWidget | None): Parent widget

#### Helper Methods

**Status Bar**
```python
self.update_status("Processing...", timeout=5000)  # 5 second timeout
self.clear_status()
```

**Menu Bar**
```python
menu_bar = self.create_menu_bar()
file_menu = menu_bar.addMenu("&File")
file_menu.addAction("&Open", self.open_file)
```

**Window Positioning**
```python
self.center_on_screen()  # Center on screen
self.set_central_widget_margins(10, 10, 10, 10)  # Set margins
```

**Stylesheet Access**
```python
stylesheet_path = self.get_stylesheet_path()  # Get path to main.qss
```

---

### BaseDialog

Dialog base class that provides:

- Automatic QSS stylesheet loading
- Modal dialog configuration
- Standard button box creation (OK/Cancel, Yes/No, etc.)
- Custom button layout helpers
- Dialog positioning

#### Basic Usage

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

        # Create layout
        layout = self.create_simple_layout()

        # Add content
        layout.addWidget(QLabel("Settings content here"))

        # Add standard OK/Cancel buttons
        button_box = self.create_button_box()
        layout.addWidget(button_box)

        self.setLayout(layout)
        self.set_margins(15, 15, 15, 15)

# Usage
dialog = SettingsDialog(parent=main_window)
if dialog.exec():
    print("Settings accepted")
else:
    print("Settings cancelled")
```

#### Constructor Parameters

- `parent` (QWidget | None): Parent widget
- `title` (str): Dialog title (default: "Dialog")
- `default_size` (tuple | None): Default size as (width, height) or None for auto-size
- `min_size` (tuple | None): Minimum size or None
- `modal` (bool): Whether dialog should be modal (default: True)
- `icon_path` (str | None): Path to icon file

#### Helper Methods

**Button Boxes**
```python
# Standard OK/Cancel buttons
button_box = self.create_button_box()

# Yes/No buttons
button_box = self.create_button_box(
    buttons=QDialogButtonBox.Yes | QDialogButtonBox.No
)

# Custom handling
button_box = self.create_button_box(connect_default=False)
button_box.accepted.connect(self.custom_accept_handler)
```

**Custom Buttons**
```python
buttons = [
    ("Apply", self.apply_settings),
    ("Reset", self.reset_to_defaults),
]
button_layout = self.create_custom_buttons(buttons, add_cancel=True)
main_layout.addLayout(button_layout)
```

**Dialog Positioning**
```python
self.center_on_parent()  # Center on parent window
self.center_on_screen()  # Center on screen
```

**Layout Helpers**
```python
layout = BaseDialog.create_simple_layout(spacing=10)
self.set_margins(10, 10, 10, 10)
```

---

## Design Patterns

### Pattern 1: Main Application Window

```python
from lunaNMR.gui.base import BaseWindow
from PySide6.QtWidgets import QApplication, QVBoxLayout, QWidget
import sys

class LunaNMRMainWindow(BaseWindow):
    def __init__(self):
        super().__init__(
            title="lunaNMR v1.0",
            default_size=(1400, 900),
            min_size=(800, 600),
            enable_status_bar=True
        )

        self._setup_ui()
        self._create_menus()
        self.center_on_screen()
        self.update_status("Ready")

    def _setup_ui(self):
        """Setup the main user interface."""
        central_widget = QWidget()
        layout = QVBoxLayout(central_widget)
        # Add your widgets here
        self.setCentralWidget(central_widget)

    def _create_menus(self):
        """Create application menus."""
        menu_bar = self.create_menu_bar()

        file_menu = menu_bar.addMenu("&File")
        file_menu.addAction("&Open", self.open_file)
        file_menu.addSeparator()
        file_menu.addAction("E&xit", self.close)

def main():
    app = QApplication(sys.argv)
    window = LunaNMRMainWindow()
    window.show()
    sys.exit(app.exec())
```

### Pattern 2: Settings Dialog

```python
from lunaNMR.gui.base import BaseDialog
from PySide6.QtWidgets import QVBoxLayout, QLabel, QSpinBox, QGroupBox

class SettingsDialog(BaseDialog):
    def __init__(self, parent=None):
        super().__init__(
            parent=parent,
            title="Application Settings",
            default_size=(500, 400),
            modal=True
        )

        layout = self.create_simple_layout(spacing=15)

        # Add settings groups
        layout.addWidget(self._create_general_group())
        layout.addWidget(self._create_processing_group())

        # Add buttons
        button_box = self.create_button_box()
        layout.addWidget(button_box)

        self.setLayout(layout)
        self.set_margins(20, 20, 20, 20)
        self.center_on_parent()

    def _create_general_group(self):
        """Create general settings group."""
        group = QGroupBox("General Settings")
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Configuration options here"))
        group.setLayout(layout)
        return group

    def _create_processing_group(self):
        """Create processing settings group."""
        group = QGroupBox("Processing Settings")
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Processing options here"))
        group.setLayout(layout)
        return group
```

### Pattern 3: Confirmation Dialog

```python
from lunaNMR.gui.base import BaseDialog
from PySide6.QtWidgets import QDialogButtonBox, QLabel

class ConfirmDialog(BaseDialog):
    def __init__(self, message, title="Confirm", parent=None):
        super().__init__(
            parent=parent,
            title=title,
            default_size=(400, 150),
            modal=True
        )

        layout = self.create_simple_layout()
        layout.addWidget(QLabel(message))

        button_box = self.create_button_box(
            buttons=QDialogButtonBox.Yes | QDialogButtonBox.No
        )
        layout.addWidget(button_box)

        self.setLayout(layout)
        self.center_on_parent()

# Usage
dialog = ConfirmDialog("Are you sure you want to delete this?", parent=main_window)
if dialog.exec():
    # User clicked Yes
    perform_deletion()
```

### Pattern 4: Custom Button Dialog

```python
from lunaNMR.gui.base import BaseDialog
from PySide6.QtWidgets import QLabel

class ProcessDialog(BaseDialog):
    def __init__(self, parent=None):
        super().__init__(
            parent=parent,
            title="Process Spectra",
            default_size=(450, 250),
            modal=True
        )

        layout = self.create_simple_layout()
        layout.addWidget(QLabel("Choose processing action:"))

        # Custom buttons
        buttons = [
            ("Process All", self.process_all),
            ("Process Selected", self.process_selected),
        ]
        button_layout = self.create_custom_buttons(buttons, add_cancel=True)
        layout.addLayout(button_layout)

        self.setLayout(layout)
        self.center_on_parent()

    def process_all(self):
        """Process all spectra."""
        # Implementation here
        self.accept()

    def process_selected(self):
        """Process selected spectra."""
        # Implementation here
        self.accept()
```

---

## Stylesheet Customization

The base classes automatically load `lunaNMR/gui/styles/main.qss`. This stylesheet provides:

- Dark theme with professional color scheme
- Consistent styling across all Qt widgets
- Hover and focus states
- Disabled state styling
- Special "primary" button class for OK/Accept buttons

### Custom Widget Styling

To apply the primary button style:

```python
button = QPushButton("OK")
button.setProperty("class", "primary")
```

This is automatically applied to OK, Yes, and Apply buttons in `create_button_box()`.

---

## Error Handling

Both base classes handle errors gracefully:

**Missing Stylesheet**
- Warning logged to console
- Application continues with default Qt styling
- No crash or exception

**Missing Icon**
- Warning logged to console
- Window/dialog displays without icon
- No crash or exception

**Example Log Output:**
```
WARNING - Stylesheet not found at /path/to/main.qss. Window will use default Qt styling.
WARNING - Icon file not found or invalid: /path/to/icon.png
```

---

## Testing

A test script is provided at `test_base_classes.py`:

```bash
cd lunaNMR_v1o0
python3 test_base_classes.py
```

This demonstrates:
- BaseWindow with menu bar and status bar
- BaseDialog with standard buttons
- Parent-child dialog relationships
- Centering and positioning
- Stylesheet application

---

## Implementation Notes

### Design Decisions

1. **Automatic Stylesheet Loading**: Reduces boilerplate code in every window/dialog
2. **Sensible Defaults**: Common configurations (modal dialogs, status bars) are opt-in
3. **Flexible Button Creation**: Both standard and custom button patterns supported
4. **Error Tolerance**: Missing resources don't crash the application
5. **Comprehensive Docstrings**: Every method includes usage examples

### Best Practices

1. **Always call super().__init__()** with appropriate parameters
2. **Use helper methods** instead of reimplementing common patterns
3. **Center dialogs on parent** for better UX: `dialog.center_on_parent()`
4. **Set margins** for visual breathing room: `self.set_margins(15, 15, 15, 15)`
5. **Use create_simple_layout()** for consistent spacing

### Future Extensions

The base classes are designed to be extended. Potential additions:

- BaseToolWindow for floating tool palettes
- BaseSplashScreen for startup screens
- BaseProgressDialog for long-running operations
- Additional helper methods as patterns emerge

---

## File Structure

```
lunaNMR_v1o0/
├── lunaNMR/
│   └── gui/
│       ├── base/
│       │   ├── __init__.py           # Exports BaseWindow, BaseDialog
│       │   ├── base_window.py        # BaseWindow implementation
│       │   └── base_dialog.py        # BaseDialog implementation
│       └── styles/
│           └── main.qss              # Main stylesheet
└── test_base_classes.py              # Test/demo script
```

---

## Quick Reference

### Import
```python
from lunaNMR.gui.base import BaseWindow, BaseDialog
```

### Create Window
```python
window = BaseWindow(
    title="My Window",
    default_size=(800, 600),
    enable_status_bar=True
)
```

### Create Dialog
```python
dialog = BaseDialog(
    parent=parent_window,
    title="My Dialog",
    default_size=(400, 300),
    modal=True
)
```

### Standard Patterns
```python
# Status updates
self.update_status("Processing...", timeout=3000)

# Menu creation
menu_bar = self.create_menu_bar()

# Button boxes
button_box = self.create_button_box()

# Positioning
self.center_on_screen()
self.center_on_parent()
```
