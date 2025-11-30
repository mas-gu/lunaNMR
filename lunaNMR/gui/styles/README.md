# LunaNMR Qt Design System

This directory contains the complete design system for the Qt/PySide6 version of LunaNMR v1.0, ported from the CustomTkinter implementation.

## Files

- **`design_system.py`** - All design constants (colors, spacing, typography, sizing)
- **`main.qss`** - Qt stylesheet (QSS) implementing the Apple-style design
- **`__init__.py`** - Module initialization with convenience functions
- **`example_usage.py`** - Demo application showing styled widgets

## Quick Start

### 1. Load Stylesheet in Application

```python
from PySide6.QtWidgets import QApplication
from lunaNMR.gui.styles import load_stylesheet

app = QApplication(sys.argv)
app.setStyleSheet(load_stylesheet())
```

### 2. Use Design Constants

```python
from lunaNMR.gui.styles import (
    PRIMARY_BUTTON_BG,
    SPACING_MD,
    BUTTON_CORNER_RADIUS
)

# Use in code
button.setStyleSheet(f"background-color: {PRIMARY_BUTTON_BG};")
layout.setSpacing(SPACING_MD)
```

### 3. Apply Widget Classes

The stylesheet defines classes for different button types:

```python
from PySide6.QtWidgets import QPushButton

# Primary action button (blue)
primary_btn = QPushButton("Fit All Peaks")
primary_btn.setProperty("class", "primary")

# Secondary button (grey, default)
secondary_btn = QPushButton("Browse")
# No class needed - default styling

# Destructive button (red)
delete_btn = QPushButton("Clear All")
delete_btn.setProperty("class", "destructive")

# Success button (green)
success_btn = QPushButton("Detect Peaks")
success_btn.setProperty("class", "success")

# Icon button (compact)
icon_btn = QPushButton("🔄")
icon_btn.setProperty("class", "icon")
```

## Design System Overview

### Color Palette

**Backgrounds:**
- `BG_COLOR = "#FAFAFA"` - Main window background
- `PANEL_BG_COLOR = "#F5F5F7"` - Secondary panels
- `FRAME_BG_COLOR = "#FFFFFF"` - Cards and containers

**Text:**
- `PRIMARY_TEXT = "#1C1C1E"` - Main text
- `SECONDARY_TEXT = "#8E8E93"` - Help text
- `DISABLED_TEXT = "#C7C7CC"` - Disabled state

**Buttons:**
- Primary: `#5B9EE5` (soft blue)
- Secondary: `#E5E5EA` (light grey)
- Destructive: `#E8554E` (soft red)
- Success: `#34C759` (green)

**Accents:**
- Success: `#34C759`
- Warning: `#F0A04B`
- Error: `#E8554E`
- Info: `#5B9EE5`

### Spacing System (8pt Grid)

```python
SPACING_XS = 4   # Tight spacing
SPACING_SM = 8   # Between elements
SPACING_MD = 16  # Between sections
SPACING_LG = 24  # Major breaks
SPACING_XL = 32  # Window padding
```

### Corner Radius

```python
BUTTON_CORNER_RADIUS = 10  # All buttons
FRAME_CORNER_RADIUS = 12   # Panels and group boxes
DIALOG_CORNER_RADIUS = 14  # Modal dialogs
CARD_CORNER_RADIUS = 8     # Small UI cards
```

### Typography

**Font Sizes:**
- Large Header: 14pt (bold)
- Medium Header: 12pt (bold)
- Section Label: 11pt (bold)
- Body Text: 10pt
- Small Text: 9pt
- Tiny Text: 8pt

### Button Sizing

```python
BUTTON_WIDTH_ICON = 36      # Icon buttons
BUTTON_WIDTH_SMALL = 60     # Small buttons
BUTTON_WIDTH_STANDARD = 80  # Standard buttons
BUTTON_WIDTH_MEDIUM = 100   # Important actions
BUTTON_WIDTH_LARGE = 180    # Main actions
BUTTON_WIDTH_DIALOG = 120   # Dialog buttons

BUTTON_HEIGHT_STANDARD = 28 # Standard height
BUTTON_HEIGHT_DIALOG = 40   # Dialog height
```

## QSS Widget Styling

### Buttons (QPushButton)

Classes available:
- `.primary` - Blue primary action buttons
- `.secondary` - Grey secondary buttons (default)
- `.destructive` - Red destructive action buttons
- `.success` - Green success buttons
- `.icon` - Compact icon buttons (36×36px)
- `.dialog` - Dialog buttons (120×40px)

### Frames and Containers (QFrame, QGroupBox)

Classes available:
- `.panel` - Secondary panel background
- `.card` - White card background with border

QGroupBox automatically styled with rounded corners and proper title positioning.

### Labels (QLabel)

Classes available:
- `.header-large` - 14pt bold
- `.header-medium` - 12pt bold
- `.section-label` - 11pt bold
- `.secondary` - Grey secondary text
- `.disabled` - Disabled state

### Input Fields

All automatically styled:
- QLineEdit - White background, rounded corners, blue focus
- QTextEdit, QPlainTextEdit - Same styling as QLineEdit
- QComboBox - Dropdown with custom arrow
- QSpinBox, QDoubleSpinBox - Numeric inputs with styled buttons

### Other Widgets

All have Apple-style rounded, soft styling:
- QCheckBox - Rounded checkboxes with blue checked state
- QRadioButton - Circular radio buttons
- QSlider - Rounded groove with circular handle
- QProgressBar - Rounded with green fill (or warning/error colors)
- QTabWidget - Rounded tabs with blue selected state
- QMenuBar, QMenu - Soft rounded menus
- QTableView, QListView, QTreeView - Consistent selection colors
- QScrollBar - Thin, rounded scrollbars

## Mapping from CustomTkinter to Qt

| CustomTkinter | Qt Widget | Notes |
|--------------|-----------|-------|
| `CTkButton` | `QPushButton` | Use `.primary`, `.secondary` classes |
| `CTkFrame` | `QFrame` or `QGroupBox` | QGroupBox for labeled frames |
| `CTkLabel` | `QLabel` | Use classes for styling |
| `CTkEntry` | `QLineEdit` | Auto-styled |
| `CTkTextbox` | `QTextEdit` or `QPlainTextEdit` | Auto-styled |
| `CTkComboBox` | `QComboBox` | Auto-styled |
| `CTkCheckBox` | `QCheckBox` | Auto-styled |
| `CTkRadioButton` | `QRadioButton` | Auto-styled |
| `CTkSlider` | `QSlider` | Auto-styled |
| `CTkProgressBar` | `QProgressBar` | Auto-styled |
| `CTkTabview` | `QTabWidget` | Auto-styled |

## Important Notes

### Secondary Button Borders

Secondary buttons MUST have a 1px border to prevent grey-on-grey edge artifacts. This is automatically applied in the QSS:

```css
QPushButton.secondary, QPushButton {
    border: 1px solid #C8C8CD;  /* REQUIRED */
}
```

### Widget Classes in Qt

To apply a class to a widget, use `setProperty()`:

```python
button.setProperty("class", "primary")
```

After setting properties, you may need to refresh the stylesheet:

```python
button.style().unpolish(button)
button.style().polish(button)
```

### Color Consistency

All colors are centralized in `design_system.py`. When adding custom styling, always import and use these constants rather than hardcoding colors:

```python
# Good
from lunaNMR.gui.styles import PRIMARY_BUTTON_BG
widget.setStyleSheet(f"background-color: {PRIMARY_BUTTON_BG};")

# Bad
widget.setStyleSheet("background-color: #5B9EE5;")
```

## Testing

Run the validation:

```bash
python3 -m lunaNMR.gui.styles.design_system
```

Run the demo application (requires PySide6):

```bash
cd lunaNMR/gui/styles
python3 example_usage.py
```

## Design Philosophy

The design follows Apple's modern design language:

1. **Clean & Minimal** - Reduce visual noise
2. **Soft Colors** - Pleasant for extended use
3. **High Contrast** - Readable without harsh colors
4. **Consistent Spacing** - 8pt grid system
5. **Rounded Corners** - Modern, friendly aesthetic
6. **Contextual Colors** - Color communicates meaning

## Future Enhancements

- Dark mode support (alternative QSS)
- High contrast mode for accessibility
- Custom icon set matching design system
- Animation support for state transitions

## References

- Source: `LUNАНMR_UX_STYLE_GUIDE.md`
- CustomTkinter constants: `lunaNMR/gui/gui_components.py`
- Qt Style Sheets: https://doc.qt.io/qt-6/stylesheet.html
