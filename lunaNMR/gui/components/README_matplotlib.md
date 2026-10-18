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
   - Figure background: `#F5F5F7` (PANEL_BG_COLOR)
   - Axes background: `#FFFFFF` (FRAME_BG_COLOR), set by `clear_full()`, `get_axes()` and
     on every `MatplotlibMultiAxesWidget` subplot
   - Text colors match design system

3. **Helper Methods**
   - `clear()` - Clear the axes
   - `clear_full()` - Clear the axes and reset their background
   - `refresh()` - Redraw canvas
   - `set_title()` - Set plot title with styling
   - `add_subplot()` / `get_axes()` - Add a subplot / get the primary axes
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
plot_widget.axes.set_xlabel("X")
plot_widget.axes.set_ylabel("Y")
plot_widget.axes.grid(True)
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

# axes_list is in row-major order: row 0, col 1 of a 2x2 grid is index 1
ax = multi_plot.axes_list[0 * multi_plot.ncols + 1]
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

**Critical: Set the Backend First**

`launch_lunaNMR.py` forces the QtAgg backend at its very top, before any other matplotlib
import, so that `NavigationToolbar2QT` renders properly:

```python
# As in launch_lunaNMR.py
import matplotlib
matplotlib.use('QtAgg', force=True)

from lunaNMR.gui.components.matplotlib_widget import MatplotlibWidget
```

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

# Applied by clear_full(), get_axes() and each MatplotlibMultiAxesWidget subplot
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
    toolbar=False,      # Include navigation toolbar
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
Clear the axes (`Axes.clear()`: data, title and axis labels).

**clear_full()**
`clear()`, then reset the axes background to FRAME_BG_COLOR.

**refresh()**
Redraw canvas with latest changes.

**set_title(title, **kwargs)**
Set plot title with design system styling.

**add_subplot(nrows, ncols, index)**
Add a subplot to the figure (1-indexed) and return its Axes.

**get_axes(create_if_missing=True)**
Return the primary axes, creating one (FRAME_BG_COLOR background) if the figure has none.

**save_figure(filename, **kwargs)**
Export figure to file (defaults: `dpi=300`, `bbox_inches='tight'`, white face color).

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

### Verify Integration

1. The QtAgg backend is set before any other matplotlib import
2. Backend is QtAgg (not TkAgg)
3. Canvas type is FigureCanvasQTAgg
4. Toolbar type is NavigationToolbar2QT
5. Design system colors applied
6. All helper methods work

## Common Issues

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

**Solution:** Let widget apply colors, or set the axes background yourself

```python
# Automatic
widget = MatplotlibWidget()  # Colors applied

# Manual
widget.axes.set_facecolor("#FFFFFF")
widget.refresh()
```

## Performance Tips

1. **Batch Updates:** Make all plot changes before calling refresh()
2. **Clear Efficiently:** Use clear() instead of clear_full() when possible
3. **Toolbar:** Leave `toolbar=False` (the default) if not needed to save memory
4. **Figure Size:** Choose appropriate figsize for your use case
5. **DPI:** Lower DPI (75-100) for faster rendering, higher (150-300) for export

## Examples

See the following files for complete examples:
- `lunaNMR/gui/components/spectrum_plotter.py` - `SpectrumPlotter(MatplotlibWidget)`
- `lunaNMR/gui/components/voigt_analysis_plotter.py` - `VoigtAnalysisPlotter(MatplotlibMultiAxesWidget)`, a 2x2 grid
- `lunaNMR/gui/dialogs/results_browser_dialog.py` - plain `MatplotlibWidget` charts

---

Author: Guillaume Mas
Date: 2025-01-24
Module: lunaNMR.gui.components.matplotlib_widget
