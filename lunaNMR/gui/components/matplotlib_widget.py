#!/usr/bin/env python3
# ABOUTME: Base Qt matplotlib widget for embedding plots in LunaNMR Qt interface
# ABOUTME: Provides FigureCanvasQTAgg integration with navigation toolbar and design system styling

"""
MatplotlibWidget - Base Qt Widget for Matplotlib Integration

This module provides a reusable Qt widget that embeds matplotlib figures with the Qt backend.
It replaces the TkAgg integration from the CustomTkinter version with proper Qt support.

Key Features:
- FigureCanvasQTAgg integration (not TkAgg)
- Optional NavigationToolbar2QT for pan/zoom
- Design system color integration
- Auto-configured layout management
- Helper methods for common operations

Usage:
    from lunaNMR.gui.components.matplotlib_widget import MatplotlibWidget

    # Create widget with default settings
    plot_widget = MatplotlibWidget(parent=self, toolbar=True)

    # Plot data
    plot_widget.axes.plot(x_data, y_data)
    plot_widget.canvas.draw()

    # Or use convenience methods
    plot_widget.clear()
    plot_widget.refresh()

Author: Guillaume Mas
Date: 2025-01-24
"""

from PySide6.QtWidgets import QWidget, QVBoxLayout, QSizePolicy
from PySide6.QtCore import Qt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
import numpy as np

from lunaNMR.gui.styles.design_system import (
    PANEL_BG_COLOR,
    FRAME_BG_COLOR,
    PRIMARY_TEXT,
    SPACING_SM,
)


class MatplotlibWidget(QWidget):
    """
    Base Qt widget for matplotlib figure embedding.

    This widget provides a standardized way to embed matplotlib plots in Qt
    interfaces. It handles figure creation, canvas setup, optional toolbar,
    and design system integration.

    Attributes:
        figure (Figure): Matplotlib figure instance
        canvas (FigureCanvasQTAgg): Qt canvas for rendering the figure
        axes (Axes): Primary axes object for plotting
        toolbar (NavigationToolbar2QT or None): Navigation toolbar if enabled

    Example:
        widget = MatplotlibWidget(parent=self, figsize=(8, 6))
        widget.axes.plot([1, 2, 3], [1, 4, 9])
        widget.refresh()
    """

    def __init__(
        self,
        parent=None,
        toolbar=False,
        figsize=(8, 6),
        dpi=100,
        tight_layout=True
    ):
        """
        Initialize the matplotlib widget.

        Args:
            parent: Parent Qt widget
            toolbar: If True, include navigation toolbar (default False)
            figsize: Figure size in inches (width, height)
            dpi: Dots per inch for figure resolution
            tight_layout: If True, use tight_layout for better spacing
        """
        super().__init__(parent)

        self.tight_layout = tight_layout

        # Create matplotlib figure with design system colors
        self.figure = Figure(figsize=figsize, dpi=dpi)
        self._apply_design_system_colors()

        # Create Qt canvas
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        # Create primary axes
        self.axes = self.figure.add_subplot(111)

        # Setup layout
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(SPACING_SM)

        # Add toolbar if requested
        self.toolbar = None
        if toolbar:
            self.toolbar = NavigationToolbar2QT(self.canvas, self)
            self.layout.addWidget(self.toolbar)

        # Add canvas
        self.layout.addWidget(self.canvas)

        # Apply tight layout if requested
        if self.tight_layout:
            self.figure.tight_layout()

    def _apply_design_system_colors(self):
        """
        Apply design system colors to the matplotlib figure.

        Sets appropriate background colors for the figure and plot area
        to match the LunaNMR design system.
        """
        # Convert hex colors to matplotlib format
        panel_bg = self._hex_to_mpl_color(PANEL_BG_COLOR)
        frame_bg = self._hex_to_mpl_color(FRAME_BG_COLOR)
        text_color = self._hex_to_mpl_color(PRIMARY_TEXT)

        # Set figure and axes backgrounds
        self.figure.patch.set_facecolor(panel_bg)
        self.figure.patch.set_alpha(1.0)

        # Set text color
        from matplotlib import rcParams
        rcParams['text.color'] = text_color
        rcParams['axes.labelcolor'] = text_color
        rcParams['xtick.color'] = text_color
        rcParams['ytick.color'] = text_color

    @staticmethod
    def _hex_to_mpl_color(hex_color: str) -> tuple:
        """
        Convert hex color to matplotlib RGB tuple.

        Args:
            hex_color: Hex color string (e.g., "#F5F5F7")

        Returns:
            Tuple of (r, g, b) values normalized to 0-1 range
        """
        hex_color = hex_color.lstrip('#')
        r, g, b = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
        return (r/255.0, g/255.0, b/255.0)

    def clear(self):
        """
        Clear the axes without removing labels or settings.

        This removes all plotted data but preserves axis labels, limits,
        and other settings. Use clear_full() to completely reset.
        """
        self.axes.clear()

    def clear_full(self):
        """
        Completely clear and reset the axes.

        This removes all data, labels, and settings, returning the axes
        to a fresh state.
        """
        self.axes.clear()
        self.axes.set_facecolor(self._hex_to_mpl_color(FRAME_BG_COLOR))

    def refresh(self):
        """
        Refresh the canvas to display latest changes.

        Call this after updating plot data to ensure changes are rendered.
        """
        if self.tight_layout:
            self.figure.tight_layout()
        self.canvas.draw()

    def set_title(self, title: str, **kwargs):
        """
        Set the plot title with design system styling.

        Args:
            title: Title text
            **kwargs: Additional arguments passed to axes.set_title()
        """
        default_kwargs = {
            'fontsize': 12,
            'color': PRIMARY_TEXT,
            'pad': SPACING_SM
        }
        default_kwargs.update(kwargs)
        self.axes.set_title(title, **default_kwargs)

    def add_subplot(self, nrows: int, ncols: int, index: int):
        """
        Add a subplot to the figure.

        Args:
            nrows: Number of subplot rows
            ncols: Number of subplot columns
            index: Subplot index (1-indexed)

        Returns:
            Axes object for the new subplot
        """
        return self.figure.add_subplot(nrows, ncols, index)

    def get_axes(self, create_if_missing: bool = True):
        """
        Get the primary axes object.

        Args:
            create_if_missing: If True and no axes exist, create one

        Returns:
            Primary axes object or None
        """
        if not self.figure.axes and create_if_missing:
            self.axes = self.figure.add_subplot(111)
            self.axes.set_facecolor(self._hex_to_mpl_color(FRAME_BG_COLOR))

        return self.axes if self.figure.axes else None

    def save_figure(self, filename: str, **kwargs):
        """
        Save the figure to a file.

        Args:
            filename: Output filename (extension determines format)
            **kwargs: Additional arguments passed to figure.savefig()
        """
        default_kwargs = {
            'dpi': 300,
            'bbox_inches': 'tight',
            'facecolor': 'white'
        }
        default_kwargs.update(kwargs)
        self.figure.savefig(filename, **default_kwargs)


class MatplotlibMultiAxesWidget(MatplotlibWidget):
    """
    Extended matplotlib widget supporting multiple axes/subplots.

    This class extends MatplotlibWidget to support multiple subplots
    in a grid layout. Useful for complex visualizations with multiple
    related plots.

    Example:
        widget = MatplotlibMultiAxesWidget(
            parent=self,
            nrows=2,
            ncols=2,
            toolbar=True
        )
        widget.axes_list[0].plot(data1)
        widget.axes_list[1].plot(data2)
        widget.refresh()
    """

    def __init__(
        self,
        parent=None,
        nrows=1,
        ncols=1,
        toolbar=True,
        figsize=(10, 8),
        dpi=100,
        sharex=False,
        sharey=False
    ):
        """
        Initialize multi-axes widget.

        Args:
            parent: Parent Qt widget
            nrows: Number of subplot rows
            ncols: Number of subplot columns
            toolbar: If True, include navigation toolbar
            figsize: Figure size in inches (width, height)
            dpi: Dots per inch for figure resolution
            sharex: Share x-axis across subplots
            sharey: Share y-axis across subplots
        """
        # Initialize base widget without creating default axes
        super().__init__(
            parent=parent,
            toolbar=toolbar,
            figsize=figsize,
            dpi=dpi,
            tight_layout=True
        )

        # Clear the default single axes
        self.figure.clear()

        # Create subplot grid
        self.nrows = nrows
        self.ncols = ncols
        self.axes_list = []

        for i in range(nrows * ncols):
            ax = self.figure.add_subplot(nrows, ncols, i + 1)
            ax.set_facecolor(self._hex_to_mpl_color(FRAME_BG_COLOR))
            self.axes_list.append(ax)

        # Set primary axes to first subplot
        self.axes = self.axes_list[0] if self.axes_list else None

        # Apply sharing if requested
        if sharex or sharey:
            for ax in self.axes_list[1:]:
                if sharex:
                    ax.sharex(self.axes_list[0])
                if sharey:
                    ax.sharey(self.axes_list[0])

        # Apply tight layout
        self.figure.tight_layout()

    def clear_all(self):
        """Clear all subplots."""
        for ax in self.axes_list:
            ax.clear()
