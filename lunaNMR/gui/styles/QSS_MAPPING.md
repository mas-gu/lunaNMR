# CustomTkinter to Qt QSS Mapping Reference

This document shows how CustomTkinter widgets and styling map to Qt QSS selectors.

## Button Styling

### Primary Action Button (Blue)

**CustomTkinter:**
```python
ctk.CTkButton(
    text="Fit All Peaks",
    corner_radius=10,
    fg_color="#5B9EE5",
    hover_color="#4A8DD4",
    text_color="#FFFFFF",
    font=("TkDefaultFont", 14, "bold")
)
```

**Qt/PySide6:**
```python
button = QPushButton("Fit All Peaks")
button.setProperty("class", "primary")
```

**QSS Selector:**
```css
QPushButton.primary {
    background-color: #5B9EE5;
    color: #FFFFFF;
    border: none;
    border-radius: 10px;
    padding: 6px 16px;
    font-weight: bold;
    font-size: 14pt;
}

QPushButton.primary:hover {
    background-color: #4A8DD4;
}
```

---

### Secondary Button (Grey)

**CustomTkinter:**
```python
ctk.CTkButton(
    text="Browse",
    corner_radius=10,
    fg_color="#E5E5EA",
    hover_color="#D1D1D6",
    text_color="#1C1C1E",
    border_width=1,
    border_color="#C8C8CD"
)
```

**Qt/PySide6:**
```python
button = QPushButton("Browse")
# No class needed - default QPushButton styling
```

**QSS Selector:**
```css
QPushButton.secondary, QPushButton {
    background-color: #E5E5EA;
    color: #1C1C1E;
    border: 1px solid #C8C8CD;  /* CRITICAL for grey buttons */
    border-radius: 10px;
    padding: 6px 16px;
}

QPushButton.secondary:hover, QPushButton:hover {
    background-color: #D1D1D6;
}
```

---

### Destructive Button (Red)

**CustomTkinter:**
```python
ctk.CTkButton(
    text="Clear All",
    corner_radius=10,
    fg_color="#E8554E",
    hover_color="#D44943",
    text_color="#FFFFFF"
)
```

**Qt/PySide6:**
```python
button = QPushButton("Clear All")
button.setProperty("class", "destructive")
```

**QSS Selector:**
```css
QPushButton.destructive {
    background-color: #E8554E;
    color: #FFFFFF;
    border: none;
    border-radius: 10px;
    padding: 6px 16px;
}

QPushButton.destructive:hover {
    background-color: #D44943;
}
```

---

## Frame/Container Styling

### Label Frame (Panel)

**CustomTkinter:**
```python
frame = ctk.CTkFrame(
    parent,
    corner_radius=12,
    fg_color="#F5F5F7"
)

label = ctk.CTkLabel(frame, text="Section Title", font=("TkDefaultFont", 10, "bold"))
label.pack(anchor=tk.W, padx=10, pady=(10, 5))
```

**Qt/PySide6:**
```python
group = QGroupBox("Section Title")
# Automatically styled by QSS
```

**QSS Selector:**
```css
QGroupBox {
    background-color: #F5F5F7;
    border: 1px solid #D1D1D6;
    border-radius: 12px;
    margin-top: 12px;
    padding: 16px;
    font-weight: bold;
}

QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 10px;
    padding: 0 8px;
    color: #1C1C1E;
    background-color: #F5F5F7;
}
```

---

### Card/White Frame

**CustomTkinter:**
```python
frame = ctk.CTkFrame(
    parent,
    corner_radius=8,
    fg_color="#FFFFFF",
    border_width=1,
    border_color="#D1D1D6"
)
```

**Qt/PySide6:**
```python
frame = QFrame()
frame.setProperty("class", "card")
```

**QSS Selector:**
```css
QFrame.card {
    background-color: #FFFFFF;
    border: 1px solid #D1D1D6;
    border-radius: 8px;
}
```

---

## Text Input Styling

### Entry/Line Edit

**CustomTkinter:**
```python
entry = ctk.CTkEntry(
    parent,
    corner_radius=8,
    border_width=1,
    border_color="#D1D1D6"
)
```

**Qt/PySide6:**
```python
entry = QLineEdit()
# Automatically styled by QSS
```

**QSS Selector:**
```css
QLineEdit {
    background-color: #FFFFFF;
    color: #1C1C1E;
    border: 1px solid #D1D1D6;
    border-radius: 8px;
    padding: 6px 8px;
}

QLineEdit:focus {
    border-color: #5B9EE5;  /* Blue focus indicator */
}
```

---

## Label Styling

### Headers

**CustomTkinter:**
```python
# Large header
label = ctk.CTkLabel(parent, text="Main Title", font=("TkDefaultFont", 14, "bold"))

# Medium header
label = ctk.CTkLabel(parent, text="Section", font=("TkDefaultFont", 12, "bold"))

# Section label
label = ctk.CTkLabel(parent, text="Control", font=("TkDefaultFont", 11, "bold"))
```

**Qt/PySide6:**
```python
# Large header
label = QLabel("Main Title")
label.setProperty("class", "header-large")

# Medium header
label = QLabel("Section")
label.setProperty("class", "header-medium")

# Section label
label = QLabel("Control")
label.setProperty("class", "section-label")
```

**QSS Selectors:**
```css
QLabel.header-large {
    font-size: 14pt;
    font-weight: bold;
}

QLabel.header-medium {
    font-size: 12pt;
    font-weight: bold;
}

QLabel.section-label {
    font-size: 11pt;
    font-weight: bold;
}
```

---

### Secondary/Help Text

**CustomTkinter:**
```python
label = ctk.CTkLabel(
    parent,
    text="Help text",
    text_color="#8E8E93",
    font=("TkDefaultFont", 9)
)
```

**Qt/PySide6:**
```python
label = QLabel("Help text")
label.setProperty("class", "secondary")
```

**QSS Selector:**
```css
QLabel.secondary {
    color: #8E8E93;
    font-size: 9pt;
}
```

---

## Checkbox Styling

**CustomTkinter:**
```python
checkbox = ctk.CTkCheckBox(
    parent,
    text="Enable Option",
    corner_radius=4
)
```

**Qt/PySide6:**
```python
checkbox = QCheckBox("Enable Option")
# Automatically styled by QSS
```

**QSS Selector:**
```css
QCheckBox {
    color: #1C1C1E;
    spacing: 8px;
}

QCheckBox::indicator {
    width: 18px;
    height: 18px;
    border: 1px solid #D1D1D6;
    border-radius: 4px;
    background-color: #FFFFFF;
}

QCheckBox::indicator:checked {
    background-color: #5B9EE5;
    border-color: #5B9EE5;
}
```

---

## Tab Widget Styling

**CustomTkinter:**
```python
tabview = ctk.CTkTabview(parent, corner_radius=12)
tab1 = tabview.add("Tab 1")
tab2 = tabview.add("Tab 2")
```

**Qt/PySide6:**
```python
tabs = QTabWidget()
tabs.addTab(widget1, "Tab 1")
tabs.addTab(widget2, "Tab 2")
# Automatically styled by QSS
```

**QSS Selectors:**
```css
QTabWidget::pane {
    background-color: #FFFFFF;
    border: 1px solid #D1D1D6;
    border-radius: 12px;
}

QTabBar::tab {
    background-color: #E5E5EA;
    color: #1C1C1E;
    border: 1px solid #D1D1D6;
    border-bottom: none;
    border-top-left-radius: 8px;
    border-top-right-radius: 8px;
    padding: 8px 16px;
}

QTabBar::tab:selected {
    background-color: #FFFFFF;
    color: #5B9EE5;
    font-weight: bold;
}
```

---

## Progress Bar Styling

**CustomTkinter:**
```python
progress = ctk.CTkProgressBar(
    parent,
    corner_radius=10,
    progress_color="#34C759"
)
```

**Qt/PySide6:**
```python
progress = QProgressBar()
# Automatically styled by QSS
```

**QSS Selector:**
```css
QProgressBar {
    background-color: #F5F5F7;
    border: none;
    border-radius: 10px;
    height: 20px;
}

QProgressBar::chunk {
    background-color: #34C759;  /* Success green */
    border-radius: 10px;
}
```

---

## State Pseudo-Classes

Qt QSS supports powerful state-based styling:

| State | QSS Pseudo-Class | Description |
|-------|-----------------|-------------|
| Hover | `:hover` | Mouse over widget |
| Pressed | `:pressed` | Widget being clicked |
| Checked | `:checked` | Checkbox/radio checked |
| Disabled | `:disabled` | Widget disabled |
| Focus | `:focus` | Widget has keyboard focus |
| Selected | `:selected` | Item selected in list/tree |

**Example:**
```css
QPushButton:hover {
    background-color: #4A8DD4;
}

QPushButton:pressed {
    background-color: #3A7DC4;
}

QPushButton:disabled {
    background-color: #C7C7CC;
}
```

---

## Key Differences: CustomTkinter vs Qt QSS

### 1. Property Names

| CustomTkinter | Qt QSS |
|--------------|--------|
| `fg_color` | `background-color` |
| `text_color` | `color` |
| `hover_color` | `:hover { background-color }` |
| `border_width` | `border-width` or `border: 1px solid` |
| `corner_radius` | `border-radius` |

### 2. Class Application

**CustomTkinter:**
```python
# Properties set directly on widget
button = ctk.CTkButton(fg_color="#5B9EE5", corner_radius=10)
```

**Qt:**
```python
# Properties set via CSS classes
button = QPushButton()
button.setProperty("class", "primary")
# Styling comes from QSS
```

### 3. Inheritance

In Qt QSS, styling cascades:
- `QPushButton` styles ALL buttons
- `QPushButton.primary` styles only primary-class buttons
- More specific selectors override general ones

### 4. State Management

CustomTkinter uses separate properties (`hover_color`, `fg_color`), while Qt uses pseudo-classes (`:hover`, `:pressed`).

---

## Best Practices

1. **Use Classes**: Apply widget classes with `setProperty("class", "classname")` rather than inline styles
2. **Centralize Colors**: Import from `design_system.py` rather than hardcoding
3. **Follow Spacing**: Use `SPACING_*` constants for consistent layout
4. **Test States**: Always test hover, pressed, focus, and disabled states
5. **Refresh After Property Changes**: Call `widget.style().unpolish(widget)` then `widget.style().polish(widget)` if properties change dynamically

---

## Complete Widget Coverage

All these widgets are styled in `main.qss`:

- ✓ QPushButton (with 5 classes)
- ✓ QLabel (with 4 classes)
- ✓ QFrame (with 2 classes)
- ✓ QGroupBox
- ✓ QLineEdit
- ✓ QTextEdit, QPlainTextEdit
- ✓ QComboBox
- ✓ QSpinBox, QDoubleSpinBox
- ✓ QCheckBox
- ✓ QRadioButton
- ✓ QSlider
- ✓ QProgressBar
- ✓ QScrollBar
- ✓ QTabWidget, QTabBar
- ✓ QMenuBar, QMenu
- ✓ QToolTip
- ✓ QTableView
- ✓ QListView
- ✓ QTreeView
- ✓ QStatusBar
- ✓ QSplitter

All maintain the Apple-style design with soft colors, rounded corners, and consistent spacing.
