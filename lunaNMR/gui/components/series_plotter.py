#!/usr/bin/env python3
# ABOUTME: SeriesPlotter Qt component for time series analysis and visualization
# ABOUTME: Displays intensity vs time plots for multiple peaks with error bars and legend

"""
SeriesPlotter - Qt-based time series visualization component

This widget displays time series data (intensity vs time) for multiple peaks
with features including:
- Multiple series overlay with distinct colors
- Error bars for uncertainty visualization
- Professional legend with peak assignments
- Grid lines for readability
- Proper axis formatting

Author: Guillaume Mas
Date: 2025
"""

from typing import Optional, Dict
import numpy as np

from PySide6.QtWidgets import QWidget, QVBoxLayout

import matplotlib
matplotlib.use('QtAgg')  # Force Qt backend
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from lunaNMR.gui.styles.design_system import (
    BG_COLOR, PANEL_BG_COLOR, PRIMARY_TEXT, SECONDARY_TEXT,
    PRIMARY_BUTTON_BG, SUCCESS_GREEN, WARNING_ORANGE, ERROR_RED,
    INFO_BLUE
)


class SeriesPlotter(QWidget):
    """
    Time series visualization widget for NMR peak intensity analysis.

    Features:
    - Multiple time series overlay
    - Color cycling for different peaks
    - Error bars with lighter shading
    - Professional legend
    - Grid for readability
    - Scientific notation for y-axis if needed

    Usage:
        plotter = SeriesPlotter()
        plotter.plot_series(data_dict, labels_dict)
        plotter.add_error_bars(errors_dict)
    """

    # Color palette for multiple series (using design system colors + extras)
    SERIES_COLORS = [
        PRIMARY_BUTTON_BG,  # Blue
        SUCCESS_GREEN,      # Green
        WARNING_ORANGE,     # Orange
        ERROR_RED,          # Red
        INFO_BLUE,          # Alt blue
        "#9B59B6",          # Purple
        "#1ABC9C",          # Teal
        "#E67E22",          # Dark orange
        "#2ECC71",          # Bright green
        "#3498DB",          # Sky blue
    ]

    # Marker styles for series
    MARKER_STYLES = ['o', 's', '^', 'D', 'v', 'p', '*', 'h', 'x', '+']

    def __init__(self, parent=None):
        """
        Initialize SeriesPlotter widget.

        Args:
            parent: Parent QWidget (optional)
        """
        super().__init__(parent)

        # Create matplotlib figure and canvas
        self.figure = Figure(figsize=(8, 5), dpi=100, facecolor=BG_COLOR)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot(111)

        # Setup layout
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self.canvas)

        # Data storage
        self.series_data = {}  # {label: (x_data, y_data)}
        self.error_data = {}   # {label: error_values}
        self.plot_lines = {}   # {label: line_object}
        self.error_bars = {}   # {label: error_bar_object}

        # Color cycling
        self.color_index = 0
        self.marker_index = 0

        # Initialize plot appearance
        self._setup_plot_appearance()

    def _setup_plot_appearance(self):
        """Configure plot appearance with design system colors."""
        self.ax.set_facecolor(PANEL_BG_COLOR)
        self.ax.grid(True, alpha=0.3, linestyle='--', color=SECONDARY_TEXT)

        # Set axis label colors
        self.ax.tick_params(colors=PRIMARY_TEXT)
        self.ax.spines['bottom'].set_color(SECONDARY_TEXT)
        self.ax.spines['top'].set_color(SECONDARY_TEXT)
        self.ax.spines['left'].set_color(SECONDARY_TEXT)
        self.ax.spines['right'].set_color(SECONDARY_TEXT)

    def _get_next_color(self) -> str:
        """Get next color from palette with cycling."""
        color = self.SERIES_COLORS[self.color_index % len(self.SERIES_COLORS)]
        self.color_index += 1
        return color

    def _get_next_marker(self) -> str:
        """Get next marker style with cycling."""
        marker = self.MARKER_STYLES[self.marker_index % len(self.MARKER_STYLES)]
        self.marker_index += 1
        return marker

    def _color_to_lighter_shade(self, hex_color: str, alpha: float = 0.3) -> tuple:
        """
        Convert hex color to lighter RGBA tuple for error bars.

        Args:
            hex_color: Hex color string (e.g., "#5B9EE5")
            alpha: Alpha value for transparency (0.0 to 1.0)

        Returns:
            RGBA tuple (r, g, b, a) with values 0-1
        """
        hex_color = hex_color.lstrip('#')
        r, g, b = tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))
        return (r, g, b, alpha)

    def plot_series(self, data: Dict[str, tuple], labels: Optional[Dict[str, str]] = None):
        """
        Plot multiple time series.

        Args:
            data: Dictionary mapping series_id to (x_data, y_data) tuples
            labels: Optional dictionary mapping series_id to display labels
                   If None, uses series_id as label

        Example:
            data = {
                'peak_1': (time_points, intensities_1),
                'peak_2': (time_points, intensities_2),
            }
            labels = {'peak_1': 'Peak A (8.2 ppm)', 'peak_2': 'Peak B (7.8 ppm)'}
            plotter.plot_series(data, labels)
        """
        # Clear existing data
        self.clear()

        # Reset color/marker cycling
        self.color_index = 0
        self.marker_index = 0

        # Plot each series
        for series_id, (x_data, y_data) in data.items():
            # Get label
            label = labels.get(series_id, series_id) if labels else series_id

            # Get color and marker
            color = self._get_next_color()
            marker = self._get_next_marker()

            # Plot line with markers
            line = self.ax.plot(
                x_data, y_data,
                color=color,
                marker=marker,
                markersize=6,
                linewidth=2,
                label=label,
                alpha=0.8
            )[0]

            # Store data and line
            self.series_data[series_id] = (x_data, y_data, color)
            self.plot_lines[series_id] = line

        # Update plot
        self._update_plot()

    def add_error_bars(self, errors: Dict[str, np.ndarray]):
        """
        Add error bars to existing series.

        Args:
            errors: Dictionary mapping series_id to error values
                   Error values should match the length of y_data

        Example:
            errors = {
                'peak_1': np.array([0.1, 0.12, 0.08, ...]),
                'peak_2': np.array([0.15, 0.13, 0.11, ...]),
            }
            plotter.add_error_bars(errors)
        """
        for series_id, error_vals in errors.items():
            if series_id not in self.series_data:
                continue

            x_data, y_data, color = self.series_data[series_id]

            # Remove existing error bars if present
            if series_id in self.error_bars:
                self.error_bars[series_id].remove()

            # Get lighter shade of line color for error bars
            error_color = self._color_to_lighter_shade(color, alpha=0.3)

            # Plot error bars
            error_bar = self.ax.errorbar(
                x_data, y_data,
                yerr=error_vals,
                fmt='none',
                ecolor=error_color,
                elinewidth=1.5,
                capsize=3,
                capthick=1.5,
                alpha=0.5
            )

            self.error_bars[series_id] = error_bar
            self.error_data[series_id] = error_vals

        # Redraw
        self.canvas.draw()

    def set_axis_labels(self, xlabel: str = "Time / Spectrum Index",
                       ylabel: str = "Intensity (a.u.)"):
        """
        Set axis labels.

        Args:
            xlabel: X-axis label
            ylabel: Y-axis label
        """
        self.ax.set_xlabel(xlabel, color=PRIMARY_TEXT, fontsize=11)
        self.ax.set_ylabel(ylabel, color=PRIMARY_TEXT, fontsize=11)
        self.canvas.draw()

    def set_title(self, title: str):
        """
        Set plot title.

        Args:
            title: Plot title text
        """
        self.ax.set_title(title, color=PRIMARY_TEXT, fontsize=12, pad=10)
        self.canvas.draw()

    def update_legend(self, location: str = 'best'):
        """
        Update legend with current series.

        Args:
            location: Legend location ('best', 'upper right', 'upper left', etc.)
        """
        if self.plot_lines:
            legend = self.ax.legend(
                loc=location,
                frameon=True,
                fancybox=True,
                shadow=False,
                framealpha=0.9,
                fontsize=9
            )
            legend.get_frame().set_facecolor('white')
            legend.get_frame().set_edgecolor(SECONDARY_TEXT)
            self.canvas.draw()

    def clear(self):
        """Clear all series and reset plot."""
        self.ax.clear()
        self.series_data.clear()
        self.error_data.clear()
        self.plot_lines.clear()
        self.error_bars.clear()
        self._setup_plot_appearance()
        self.canvas.draw()

    def _update_plot(self):
        """Update plot with current settings."""
        # Set default labels if not set
        if not self.ax.get_xlabel():
            self.set_axis_labels()

        # Update legend
        self.update_legend()

        # Format y-axis with scientific notation if values are very large/small
        y_data_all = [y for _, y, _ in self.series_data.values()]
        if y_data_all:
            y_flat = np.concatenate(y_data_all)
            y_max = np.max(np.abs(y_flat))
            if y_max > 1e4 or (y_max < 1e-2 and y_max > 0):
                self.ax.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))

        # Tight layout
        self.figure.tight_layout()

        # Draw
        self.canvas.draw()

    def _plot_no_series_data(self):
        """Plot message when no series data is available."""
        self.clear()
        self.ax.text(0.5, 0.5,
                    'No series results available\nRun series integration to view analysis',
                    transform=self.ax.transAxes, ha='center', va='center',
                    fontsize=10, color=SECONDARY_TEXT,
                    bbox=dict(boxstyle="round,pad=0.5", facecolor=PANEL_BG_COLOR))
        self.ax.set_title('Series Analysis - No Data', color=PRIMARY_TEXT, fontsize=12)
        self.ax.axis('off')
        self.canvas.draw()
