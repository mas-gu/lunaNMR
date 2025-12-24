#!/usr/bin/env python3
# ABOUTME: Mouse and keyboard navigation for 2D NMR spectrum widgets
# ABOUTME: Provides pan, zoom, and keyboard shortcuts with NMR axis convention support

"""
NMRNavigationHandler - Interactive Navigation for 2D NMR Spectrum Widgets

This module provides intuitive mouse and keyboard navigation for matplotlib-based
2D NMR spectrum widgets. It replaces the standard matplotlib toolbar with direct
manipulation controls optimized for NMR workflow.

Features:
- Left-click drag: Pan/move around the spectrum
- Scroll wheel: Zoom in/out centered on cursor position
- Arrow keys: Pan in increments
- +/- keys: Zoom in/out centered on view
- Home key: Reset to full spectrum view
- Shift+click: Delegate to peak addition
- Ctrl+click: Delegate to peak deletion

NMR Axis Convention:
- X-axis (1H): High ppm on LEFT, low ppm on RIGHT (inverted)
- Y-axis (15N): High ppm at BOTTOM, low ppm at TOP (inverted)

Author: Guillaume Mas
Date: 2025
"""

from typing import Tuple, Optional, Callable
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication


def calculate_zoom_limits(
    x_lim: Tuple[float, float],
    y_lim: Tuple[float, float],
    x_cursor: float,
    y_cursor: float,
    factor: float
) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """
    Calculate new axis limits for zoom operation, keeping cursor position fixed.

    This algorithm handles inverted NMR axes where x_min > x_max and y_min > y_max.

    Args:
        x_lim: Current X axis limits (min, max) - may be inverted
        y_lim: Current Y axis limits (min, max) - may be inverted
        x_cursor: Cursor X position in data coordinates
        y_cursor: Cursor Y position in data coordinates
        factor: Zoom factor (>1 = zoom in, <1 = zoom out)

    Returns:
        Tuple of (new_x_lim, new_y_lim)
    """
    x_min, x_max = x_lim
    y_min, y_max = y_lim

    # Current ranges (use signed difference to preserve direction)
    x_range = x_min - x_max  # Positive if inverted (high ppm on left)
    y_range = y_min - y_max  # Positive if inverted (high ppm on bottom)

    # New ranges after zoom (divide by factor to zoom in)
    new_x_range = x_range / factor
    new_y_range = y_range / factor

    # Calculate cursor's relative position in current view (0-1)
    # For inverted axes: rel = (cursor - max) / (min - max)
    if abs(x_range) > 1e-10:
        rel_x = (x_cursor - x_max) / x_range
    else:
        rel_x = 0.5

    if abs(y_range) > 1e-10:
        rel_y = (y_cursor - y_max) / y_range
    else:
        rel_y = 0.5

    # Calculate new limits keeping cursor at same relative position
    # new_max = cursor - rel * new_range
    # new_min = new_max + new_range
    new_x_max = x_cursor - rel_x * new_x_range
    new_x_min = new_x_max + new_x_range

    new_y_max = y_cursor - rel_y * new_y_range
    new_y_min = new_y_max + new_y_range

    return (new_x_min, new_x_max), (new_y_min, new_y_max)


def calculate_pan_limits(
    x_lim: Tuple[float, float],
    y_lim: Tuple[float, float],
    start_x: float,
    start_y: float,
    end_x: float,
    end_y: float
) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """
    Calculate new axis limits for pan operation.

    Uses "grab mode" - like grabbing and dragging the data/paper.
    Dragging right moves the data right (reveals higher ppm on left in NMR).

    Args:
        x_lim: Current X axis limits (min, max)
        y_lim: Current Y axis limits (min, max)
        start_x: Pan start X position in data coordinates
        start_y: Pan start Y position in data coordinates
        end_x: Pan end X position in data coordinates
        end_y: Pan end Y position in data coordinates

    Returns:
        Tuple of (new_x_lim, new_y_lim)
    """
    # Calculate delta in data coordinates
    delta_x = end_x - start_x
    delta_y = end_y - start_y

    # Grab mode: subtract delta to move data in drag direction
    new_x_min = x_lim[0] - delta_x
    new_x_max = x_lim[1] - delta_x
    new_y_min = y_lim[0] - delta_y
    new_y_max = y_lim[1] - delta_y

    return (new_x_min, new_x_max), (new_y_min, new_y_max)


class NMRNavigationHandler:
    """
    Navigation handler providing mouse and keyboard interaction for NMR spectrum widgets.

    This class attaches to a MatplotlibWidget and provides:
    - Left-click drag to pan
    - Scroll wheel to zoom at cursor
    - Arrow keys to pan, +/- to zoom, Home to reset
    - Modifier keys (Shift/Ctrl) to delegate to peak editing

    Usage:
        handler = NMRNavigationHandler()
        handler.attach(spectrum_plotter)  # Attach to widget
        handler.on_peak_edit = my_peak_edit_callback  # Optional callback

    Attributes:
        on_peak_edit: Optional callback for peak edit requests (x, y, modifiers)
        on_reset_zoom: Optional callback for reset zoom requests
    """

    # Zoom factor per scroll step
    ZOOM_FACTOR = 1.15

    # Pan step as fraction of view
    PAN_STEP = 0.1

    def __init__(self):
        """Initialize the navigation handler."""
        self._widget = None
        self._canvas = None
        self._ax = None
        self._connection_ids = []

        # Pan state - store both data and pixel positions
        self._pan_active = False
        self._pan_start_x_data = None
        self._pan_start_y_data = None
        self._pan_start_xlim = None
        self._pan_start_ylim = None

        # Callbacks
        self.on_peak_edit: Optional[Callable[[float, float, Qt.KeyboardModifiers], None]] = None
        self.on_reset_zoom: Optional[Callable[[], None]] = None

    def attach(self, widget):
        """
        Attach the navigation handler to a matplotlib widget.

        Args:
            widget: A MatplotlibWidget or similar with canvas and axes attributes
        """
        self.detach()  # Disconnect any previous widget

        self._widget = widget
        self._canvas = widget.canvas
        self._ax = widget.ax if hasattr(widget, 'ax') else widget.axes

        # Set focus policy to receive keyboard events
        self._canvas.setFocusPolicy(Qt.StrongFocus)

        # Connect matplotlib events
        self._connection_ids = [
            self._canvas.mpl_connect('button_press_event', self._on_mouse_press),
            self._canvas.mpl_connect('button_release_event', self._on_mouse_release),
            self._canvas.mpl_connect('motion_notify_event', self._on_mouse_motion),
            self._canvas.mpl_connect('scroll_event', self._on_scroll),
            self._canvas.mpl_connect('key_press_event', self._on_key_press),
        ]

        # Set default cursor
        self._canvas.setCursor(Qt.OpenHandCursor)

    def detach(self):
        """Disconnect the navigation handler from the current widget."""
        if self._canvas is not None:
            for cid in self._connection_ids:
                self._canvas.mpl_disconnect(cid)
            self._connection_ids = []

            # Restore default cursor
            self._canvas.setCursor(Qt.ArrowCursor)

        self._widget = None
        self._canvas = None
        self._ax = None
        self._pan_active = False

    def _on_mouse_press(self, event):
        """Handle mouse button press event."""
        if event.inaxes != self._ax:
            return

        if event.button != 1:  # Only handle left click
            return

        # Check for modifier keys (delegate to peak editing)
        modifiers = QApplication.keyboardModifiers()

        if modifiers & (Qt.ShiftModifier | Qt.ControlModifier):
            # Modifier held: delegate to peak editing callback
            if self.on_peak_edit is not None and event.xdata is not None:
                self.on_peak_edit(event.xdata, event.ydata, modifiers)
            return

        # No modifier: start pan operation
        if event.xdata is not None and event.ydata is not None:
            self._pan_active = True
            self._pan_start_x_data = event.xdata
            self._pan_start_y_data = event.ydata
            # Store original limits to avoid drift
            self._pan_start_xlim = self._ax.get_xlim()
            self._pan_start_ylim = self._ax.get_ylim()
            self._canvas.setCursor(Qt.ClosedHandCursor)

    def _on_mouse_release(self, event):
        """Handle mouse button release event."""
        if event.button == 1 and self._pan_active:
            self._pan_active = False
            self._canvas.setCursor(Qt.OpenHandCursor)

    def _on_mouse_motion(self, event):
        """Handle mouse motion event for panning."""
        if not self._pan_active:
            return

        if event.inaxes != self._ax:
            return

        if event.x is None or event.y is None:
            return

        # Convert current mouse position to data coordinates using ORIGINAL limits
        # This avoids drift that would occur if we used event.xdata directly
        # (which is calculated from current, possibly already-panned limits)
        inv = self._ax.transData.inverted()

        # Get the mouse position in screen (display) coordinates
        # and convert to data using current transform
        # But we want to know what the data position would be at the ORIGINAL limits

        # Simple approach: calculate offset in pixels, convert to data units
        bbox = self._ax.get_window_extent()

        # Calculate fractional position of mouse in axes (0-1)
        frac_x = (event.x - bbox.x0) / bbox.width
        frac_y = (event.y - bbox.y0) / bbox.height

        # Convert to data coordinates using ORIGINAL limits
        # For inverted NMR axes: xlim[0] is LEFT edge, xlim[1] is RIGHT edge
        orig_x_range = self._pan_start_xlim[0] - self._pan_start_xlim[1]
        orig_y_range = self._pan_start_ylim[0] - self._pan_start_ylim[1]

        # Current data position using original limits
        # At frac=0 (left/bottom edge), data = xlim[0]/ylim[0]
        # At frac=1 (right/top edge), data = xlim[1]/ylim[1]
        current_x_data = self._pan_start_xlim[0] - frac_x * orig_x_range
        current_y_data = self._pan_start_ylim[0] - frac_y * orig_y_range

        # Calculate pan delta from original click position
        new_x_lim, new_y_lim = calculate_pan_limits(
            self._pan_start_xlim, self._pan_start_ylim,
            self._pan_start_x_data, self._pan_start_y_data,
            current_x_data, current_y_data
        )

        # Apply new limits
        self._ax.set_xlim(new_x_lim)
        self._ax.set_ylim(new_y_lim)

        self._canvas.draw_idle()

    def _on_scroll(self, event):
        """Handle scroll wheel event for zooming."""
        if event.inaxes != self._ax:
            return

        if event.xdata is None or event.ydata is None:
            return

        # Determine zoom direction
        if event.button == 'up':
            factor = self.ZOOM_FACTOR
        elif event.button == 'down':
            factor = 1.0 / self.ZOOM_FACTOR
        else:
            return

        # Calculate new limits
        x_lim = self._ax.get_xlim()
        y_lim = self._ax.get_ylim()

        new_x_lim, new_y_lim = calculate_zoom_limits(
            x_lim, y_lim,
            event.xdata, event.ydata,
            factor
        )

        # Apply new limits
        self._ax.set_xlim(new_x_lim)
        self._ax.set_ylim(new_y_lim)

        self._canvas.draw_idle()

    def _on_key_press(self, event):
        """Handle keyboard events for navigation."""
        # Get current view
        x_lim = self._ax.get_xlim()
        y_lim = self._ax.get_ylim()

        x_range = abs(x_lim[0] - x_lim[1])
        y_range = abs(y_lim[0] - y_lim[1])

        # Calculate pan step
        pan_step_x = x_range * self.PAN_STEP
        pan_step_y = y_range * self.PAN_STEP

        if event.key == 'left':
            # Pan left (increase ppm values, since axis inverted)
            self._apply_pan(pan_step_x, 0)
        elif event.key == 'right':
            # Pan right (decrease ppm values)
            self._apply_pan(-pan_step_x, 0)
        elif event.key == 'up':
            # Pan up (increase ppm values in Y)
            self._apply_pan(0, pan_step_y)
        elif event.key == 'down':
            # Pan down (decrease ppm values in Y)
            self._apply_pan(0, -pan_step_y)
        elif event.key in ['+', '=']:
            # Zoom in (center of view)
            self._zoom_center(self.ZOOM_FACTOR)
        elif event.key == '-':
            # Zoom out
            self._zoom_center(1.0 / self.ZOOM_FACTOR)
        elif event.key == 'home':
            # Reset view
            if self.on_reset_zoom is not None:
                self.on_reset_zoom()
            elif hasattr(self._widget, 'reset_zoom'):
                self._widget.reset_zoom()

    def _apply_pan(self, delta_x: float, delta_y: float):
        """Apply a pan offset to the view."""
        x_lim = self._ax.get_xlim()
        y_lim = self._ax.get_ylim()

        self._ax.set_xlim(x_lim[0] + delta_x, x_lim[1] + delta_x)
        self._ax.set_ylim(y_lim[0] + delta_y, y_lim[1] + delta_y)

        self._canvas.draw_idle()

    def _zoom_center(self, factor: float):
        """Zoom in/out centered on the current view."""
        x_lim = self._ax.get_xlim()
        y_lim = self._ax.get_ylim()

        center_x = (x_lim[0] + x_lim[1]) / 2
        center_y = (y_lim[0] + y_lim[1]) / 2

        new_x_lim, new_y_lim = calculate_zoom_limits(
            x_lim, y_lim,
            center_x, center_y,
            factor
        )

        self._ax.set_xlim(new_x_lim)
        self._ax.set_ylim(new_y_lim)

        self._canvas.draw_idle()
