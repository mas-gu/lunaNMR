# MatplotlibWidget for Qt Integration

## Overview

The `MatplotlibWidget` provides a standardized way to embed matplotlib plots in PySide6/Qt applications. It replaces the TkAgg backend integration from CustomTkinter with proper Qt support using `FigureCanvasQTAgg`.

## Key Features

1. **Qt Backend Integration**
   - Uses `FigureCanvasQTAgg` (not `FigureCanvasTkAgg`)
   - Uses `NavigationToolbar2QT` (not `NavigationToolbar2Tk`)
   - Full Qt event system integration

2. **Design System Integration**
   - Automatic application of LunaNMR design system colors
   - Panel background: `#F5F5F7` (PANEL_BG_COLOR)
   - Figure background: `#FFFFFF` (FRAME_BG_COLOR)
   - Text colors match design system

3. **Helper Methods**
   - `clear()` - Clear plot data
   - `clear_full()` - Complete reset
   - `refresh()` - Redraw canvas
   - `set_title()` - Set plot title with styling
   - `set_labels()` - Set axis labels with styling
   - `configure_grid()` - Configure grid appearance
   - `save_figure()` - Export to file

4. **Multi-Axes Support**
   - `MatplotlibMultiAxesWidget` for subplot grids
   - Automatic layout management
   - Shared axes options

## Usage

### Basic Plot Widget

```python
from lunaNMR.gui.components.matplotlib_widget import MatplotlibWidget
import numpy as np

# Create widget
plot_widget = MatplotlibWidget(parent=self, toolbar=True, figsize=(8, 6))

# Add to layout
layout.addWidget(plot_widget)

# Plot data
x = np.linspace(0, 10, 100)
y = np.sin(x)
plot_widget.axes.plot(x, y)
plot_widget.set_title("Sine Wave")
plot_widget.set_labels("X", "Y")
plot_widget.configure_grid(True)
plot_widget.refresh()
```

### Multi-Axes Widget

```python
from lunaNMR.gui.components.matplotlib_widget import MatplotlibMultiAxesWidget

# Create 2x2 subplot grid
multi_plot = MatplotlibMultiAxesWidget(
    parent=self,
    nrows=2,
    ncols=2,
    toolbar=True
)

# Access individual subplots
multi_plot.axes_list[0].plot(x, y1)
multi_plot.axes_list[1].plot(x, y2)

# Or use grid indexing
ax = multi_plot.get_axes_at(row=0, col=1)
ax.plot(x, y3)

# Refresh all
multi_plot.refresh()
```

## Qt vs Tk Differences

### Import Changes

**Old (TkAgg):**
```python
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
```

**New (QtAgg):**
```python
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
```

### Backend Setup

**Critical: Import Order Matters**

The matplotlib backend must be set AFTER importing PySide6:

```python
# Correct order (as in launch_lunaNMR.py)
from PySide6.QtWidgets import QApplication
import matplotlib
matplotlib.use('QtAgg', force=True)

from lunaNMR.gui.components.matplotlib_widget import MatplotlibWidget
```

**Do NOT** call `matplotlib.use()` before importing PySide6, as this can cause Qt plugin loading errors.

### Canvas Integration

**Old (Tk):**
```python
canvas = FigureCanvasTkAgg(figure, parent)
canvas.get_tk_widget().pack()
```

**New (Qt):**
```python
canvas = FigureCanvasQTAgg(figure)
layout.addWidget(canvas)  # Qt uses layouts, not pack()
```

### Toolbar Integration

**Old (Tk):**
```python
toolbar = NavigationToolbar2Tk(canvas, parent)
toolbar.pack()
```

**New (Qt):**
```python
toolbar = NavigationToolbar2QT(canvas, parent)
layout.addWidget(toolbar)
```

## Design System Integration

The widget automatically applies LunaNMR design system colors:

### Figure Colors

```python
# Applied automatically in __init__
figure.patch.set_facecolor(PANEL_BG_COLOR)  # #F5F5F7
axes.set_facecolor(FRAME_BG_COLOR)          # #FFFFFF
```

### Text Colors

```python
# Applied via matplotlib rcParams
rcParams['text.color'] = PRIMARY_TEXT        # #1C1C1E
rcParams['axes.labelcolor'] = PRIMARY_TEXT
rcParams['xtick.color'] = PRIMARY_TEXT
rcParams['ytick.color'] = PRIMARY_TEXT
```

### Custom Colors

Use the color conversion utility:

```python
mpl_color = widget._hex_to_mpl_color("#5B9EE5")
axes.plot(x, y, color=mpl_color)
```

## API Reference

### MatplotlibWidget

#### Constructor

```python
MatplotlibWidget(
    parent=None,        # Parent Qt widget
    toolbar=True,       # Include navigation toolbar
    figsize=(8, 6),     # Figure size in inches
    dpi=100,           # Resolution
    tight_layout=True   # Use tight_layout
)
```

#### Attributes

- `figure` - matplotlib Figure instance
- `canvas` - FigureCanvasQTAgg instance
- `axes` - Primary Axes object
- `toolbar` - NavigationToolbar2QT or None

#### Methods

**clear()**
Clear axes without removing labels/settings.

**clear_full()**
Complete reset of axes.

**refresh()**
Redraw canvas with latest changes.

**set_title(title, **kwargs)**
Set plot title with design system styling.

**set_labels(xlabel=None, ylabel=None, **kwargs)**
Set axis labels with design system styling.

**configure_grid(visible=True, **kwargs)**
Configure grid appearance.

**set_axes_background(color=None)**
Set axes background color (defaults to FRAME_BG_COLOR).

**save_figure(filename, **kwargs)**
Export figure to file.

**enable_toolbar()** / **disable_toolbar()**
Show/hide navigation toolbar.

### MatplotlibMultiAxesWidget

Extends MatplotlibWidget with multi-subplot support.

#### Constructor

```python
MatplotlibMultiAxesWidget(
    parent=None,
    nrows=1,            # Number of rows
    ncols=1,            # Number of columns
    toolbar=True,
    figsize=(10, 8),
    dpi=100,
    sharex=False,       # Share x-axis
    sharey=False        # Share y-axis
)
```

#### Attributes

- `nrows` - Number of subplot rows
- `ncols` - Number of subplot columns
- `axes_list` - List of all Axes objects

#### Methods

**clear_all()**
Clear all subplots.

**get_axes_at(row, col)**
Get axes at specific grid position (0-indexed).

## Migration Guide

### Step 1: Replace Canvas Creation

**Before:**
```python
self.canvas = FigureCanvasTkAgg(self.figure, parent_frame)
self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
```

**After:**
```python
self.plot_widget = MatplotlibWidget(parent=parent_frame, toolbar=True)
layout.addWidget(self.plot_widget)
```

### Step 2: Update Plot References

**Before:**
```python
self.fig = Figure()
self.ax = self.fig.add_subplot(111)
self.ax.plot(x, y)
self.canvas.draw()
```

**After:**
```python
self.plot_widget.axes.plot(x, y)
self.plot_widget.refresh()
```

### Step 3: Update Toolbar

**Before:**
```python
toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
```

**After:**
```python
# Toolbar is already included if toolbar=True
# Access via self.plot_widget.toolbar if needed
```

### Step 4: Update Clear/Refresh Logic

**Before:**
```python
self.ax.clear()
self.canvas.draw()
```

**After:**
```python
self.plot_widget.clear()
self.plot_widget.refresh()
```

## Testing

### Run Unit Tests

```bash
# Simple import and instantiation test
python3 test_mpl_widget_simple.py

# Interactive demo with multiple tabs
python3 test_matplotlib_widget.py
```

### Verify Integration

1. Import order is correct (PySide6 first)
2. Backend is QtAgg (not TkAgg)
3. Canvas type is FigureCanvasQTAgg
4. Toolbar type is NavigationToolbar2QT
5. Design system colors applied
6. All helper methods work

## Common Issues

### Issue: "could not import module 'PySide6.QtGui'"

**Cause:** matplotlib.use() called before PySide6 import

**Solution:** Always import PySide6 before calling matplotlib.use()

```python
# Correct
from PySide6.QtWidgets import QApplication
import matplotlib
matplotlib.use('QtAgg', force=True)

# Wrong
import matplotlib
matplotlib.use('QtAgg', force=True)
from PySide6.QtWidgets import QApplication  # Too late!
```

### Issue: Plot not updating after changes

**Cause:** Forgot to call refresh()

**Solution:** Always call refresh() after plot modifications

```python
widget.axes.plot(x, y)
widget.set_title("New Title")
widget.refresh()  # Required!
```

### Issue: Design system colors not applied

**Cause:** Colors set before widget styling

**Solution:** Let widget apply colors, or use set_axes_background()

```python
# Automatic
widget = MatplotlibWidget()  # Colors applied

# Manual
widget.set_axes_background("#FFFFFF")
widget.refresh()
```

## Performance Tips

1. **Batch Updates:** Make all plot changes before calling refresh()
2. **Clear Efficiently:** Use clear() instead of clear_full() when possible
3. **Toolbar:** Disable toolbar if not needed to save memory
4. **Figure Size:** Choose appropriate figsize for your use case
5. **DPI:** Lower DPI (75-100) for faster rendering, higher (150-300) for export

## Examples

See the following files for complete examples:
- `test_mpl_widget_simple.py` - Basic usage
- `test_matplotlib_widget.py` - Interactive demo with all features
- `lunaNMR/gui/visualization.py` - Original TkAgg implementation (for comparison)

## Next Steps

1. Port existing plotters to use MatplotlibWidget
2. Update main_gui.py to use Qt canvas
3. Test with real NMR data
4. Verify toolbar functionality
5. Check design system color consistency

---

Author: Guillaume Mas
Date: 2025-01-24
Module: lunaNMR.gui.components.matplotlib_widget
