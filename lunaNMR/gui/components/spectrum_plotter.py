#!/usr/bin/env python3
# ABOUTME: Qt-based 2D NMR spectrum visualization with contour plotting and peak overlays
# ABOUTME: Ported from CustomTkinter visualization.py SpectrumPlotter class for PySide6

"""
SpectrumPlotter - 2D NMR Spectrum Visualization Component

This module provides a Qt widget for displaying 2D NMR spectra with contour plots,
peak annotations, and interactive features. It uses matplotlib's Qt backend for
high-quality scientific visualization.

Key Features:
- 2D contour plots with customizable levels and colormaps
- Peak position overlays with quality-based color coding
- Reference peak display
- PS2D ellipse visualization (debug mode)
- Interactive zoom and pan via matplotlib toolbar
- Crosshair and peak highlighting
- Configurable colorbar display

Author: Guillaume Mas
Date: 2025
"""

import numpy as np
from PySide6.QtCore import Signal

from lunaNMR.gui.components.matplotlib_widget import MatplotlibWidget
from lunaNMR.gui.styles.design_system import (
    PRIMARY_BUTTON_BG,
    SUCCESS_GREEN,
    ERROR_RED,
    WARNING_ORANGE,
    FRAME_BG_COLOR,
    SECONDARY_TEXT,
    PRIMARY_TEXT
)


class SpectrumPlotter(MatplotlibWidget):
    """
    Main spectrum visualization widget with peak overlays.

    Displays 2D NMR spectra as contour plots with support for:
    - Peak position markers (detected and reference)
    - Quality-based color coding (excellent/good/fair/poor)
    - PS2D elliptical integration windows
    - Interactive zoom and pan
    - Crosshair positioning
    - Peak highlighting

    Inherits from MatplotlibWidget for consistent matplotlib integration.

    Signals:
        peak_clicked: Emitted when a peak is clicked (x_ppm, y_ppm)
    """

    peak_clicked = Signal(float, float)

    def __init__(self, parent=None, toolbar=True, figsize=(8, 6)):
        """
        Initialize the spectrum plotter widget.

        Args:
            parent: Parent Qt widget
            toolbar: Show matplotlib navigation toolbar
            figsize: Figure size in inches (width, height)
        """
        super().__init__(parent=parent, toolbar=toolbar, figsize=figsize, tight_layout=False)

        # Store references for compatibility with original API
        self.fig = self.figure
        self.ax = self.axes

        # Plotting state
        self.contour_plot = None
        self.colorbar = None
        self.colorbar_ax = None
        self.peak_annotations = []
        self.crosshair_lines = []

        # Store initial axes position to prevent shrinking
        self.initial_axes_position = self.ax.get_position()

        # Calculate fixed positioning for axes and colorbar
        self._setup_fixed_layout()

        # Plot settings
        self.settings = {
            'contour_levels': 15,
            'colormap': 'viridis',
            'show_colorbar': False,
            'line_width': 0.5,
            'alpha': 0.7
        }

        # Connect click events
        self.canvas.mpl_connect('button_press_event', self._on_click)

    def plot_spectrum(self, integrator, **kwargs):
        """
        Plot 2D NMR spectrum with contours.

        Args:
            integrator: EnhancedVoigtIntegrator object with nmr_data
            **kwargs: Additional plotting options:
                - contour_min_level: Minimum contour level (fraction of max)
                - contour_levels: Number of contour levels
                - contour_increment: Multiplicative increment between levels
                - colormap: Matplotlib colormap name
                - show_colorbar: Show colorbar
        """
        if not hasattr(integrator, 'nmr_data') or integrator.nmr_data is None:
            self.ax.clear()
            self.ax.text(0.5, 0.5, 'No spectrum data loaded\nPlease load NMR data first',
                        transform=self.ax.transAxes, ha='center', va='center',
                        fontsize='small', bbox=dict(boxstyle="round,pad=0.5",
                                                    facecolor="lightgray"))
            self.canvas.draw()
            return

        # Clear previous plot
        self.ax.clear()
        self.peak_annotations = []

        # Update settings
        self.settings.update(kwargs)

        # Create coordinate grids
        X, Y = np.meshgrid(integrator.ppm_x_axis, integrator.ppm_y_axis)
        data_abs = np.abs(integrator.nmr_data)

        # Determine contour levels from GUI controls or defaults
        levels = self._calculate_contour_levels(integrator, data_abs, kwargs)

        # Create contour plot
        self.contour_plot = self.ax.contour(
            X, Y, data_abs,
            levels=levels,
            cmap=self.settings['colormap'],
            linewidths=self.settings['line_width'],
            alpha=self.settings['alpha']
        )

        # Set up axes
        self.ax.set_xlabel('¹H Chemical Shift (ppm)', fontsize='small')
        self.ax.set_ylabel('¹⁵N/¹³C Chemical Shift (ppm)', fontsize='small')
        self.ax.set_title('lunaNMR v1.0', fontsize='small', fontweight='bold')
        # NMR convention: high ppm on left (X) and bottom (Y)
        self._ensure_nmr_axis_orientation()
        self.ax.grid(True, alpha=0.3)

        # Add colorbar if requested
        if self.settings['show_colorbar'] and self.contour_plot:
            self._create_colorbar_safely()

        self.canvas.draw()

    def _calculate_contour_levels(self, integrator, data_abs, kwargs):
        """
        Calculate contour levels from parameters or defaults.

        Args:
            integrator: Integrator object with threshold
            data_abs: Absolute value of NMR data
            kwargs: User-specified parameters

        Returns:
            np.ndarray: Contour levels
        """
        if 'contour_min_level' in kwargs and 'contour_levels' in kwargs and 'contour_increment' in kwargs:
            # Use GUI-specified parameters
            min_level = kwargs['contour_min_level'] * np.max(data_abs)
            num_levels = kwargs['contour_levels']
            increment = kwargs['contour_increment']

            # Generate levels using increment method
            levels = []
            current_level = min_level

            # Ensure increment is valid for generating increasing levels
            if increment <= 1.0:
                # Convert to additive increment if multiplicative increment is too small
                max_level = np.max(data_abs)
                total_range = max_level - min_level
                additive_increment = total_range / (num_levels - 1) if num_levels > 1 else total_range

                for i in range(num_levels):
                    level = min_level + i * additive_increment
                    if level > max_level:
                        break
                    levels.append(level)
            else:
                # Use multiplicative increment for increment > 1.0
                for i in range(num_levels):
                    if current_level > np.max(data_abs):
                        break
                    levels.append(current_level)
                    current_level *= increment

            levels = np.array(levels)

            # Ensure levels are strictly increasing and remove duplicates
            levels = np.unique(levels)
            levels = np.sort(levels)
        else:
            # Default behavior
            if hasattr(integrator, 'threshold') and integrator.threshold > 0:
                min_level = integrator.threshold
            else:
                min_level = np.max(data_abs) * 0.05

            max_level = np.max(data_abs)
            levels = np.linspace(min_level, max_level, self.settings['contour_levels'])

        # Final validation: ensure levels are strictly increasing and non-empty
        if len(levels) == 0:
            # Fallback: create basic levels
            min_val = np.max(data_abs) * 0.05
            max_val = np.max(data_abs)
            levels = np.linspace(min_val, max_val, self.settings['contour_levels'])
        elif len(levels) == 1:
            # Need at least 2 levels
            levels = np.array([levels[0], levels[0] * 1.5])
        elif not np.all(levels[1:] > levels[:-1]):
            # Not strictly increasing - force it
            levels = np.linspace(np.min(levels), np.max(levels), len(levels))

        return levels

    def _setup_fixed_layout(self):
        """Setup fixed layout positions for axes and colorbar."""
        # Define fixed positions (left, bottom, width, height) in figure coordinates
        main_axes_pos = [0.12, 0.15, 0.75, 0.75]  # Leave space for colorbar
        colorbar_pos = [0.89, 0.15, 0.03, 0.75]   # Fixed colorbar position

        # Store positions
        self.main_axes_position = main_axes_pos
        self.colorbar_position = colorbar_pos

    def _create_colorbar_safely(self):
        """Safely create or update colorbar with fixed positioning."""
        try:
            # Only create colorbar if contour plot exists and has valid collections
            if not (self.contour_plot is not None and
                   hasattr(self.contour_plot, 'collections') and
                   len(self.contour_plot.collections) > 0):
                return

            # Remove existing colorbar if it exists
            if self.colorbar is not None:
                try:
                    if hasattr(self.colorbar, 'ax') and self.colorbar.ax in self.fig.axes:
                        self.fig.delaxes(self.colorbar.ax)
                except:
                    pass
                self.colorbar = None
                self.colorbar_ax = None

            # Set main axes to fixed position
            self.ax.set_position(self.main_axes_position)

            # Create dedicated colorbar axes with fixed position
            self.colorbar_ax = self.fig.add_axes(self.colorbar_position)

            # Create colorbar in the dedicated axes
            self.colorbar = self.fig.colorbar(
                self.contour_plot,
                cax=self.colorbar_ax,
                label='Intensity'
            )

            # Ensure main axes position stays fixed
            self.ax.set_position(self.main_axes_position)

        except Exception as e:
            # Silently handle colorbar creation failures
            self.colorbar = None
            self.colorbar_ax = None
            # Ensure main axes position is maintained
            try:
                self.ax.set_position(self.main_axes_position)
            except:
                pass

    def plot_peaks(self, integrator, show_detected=True, show_assigned=True, **kwargs):
        """
        Overlay peak positions and assignments using a standardized data format.

        Args:
            integrator: Integrator object with fitted_peaks and/or peak_list
            show_detected: Show detected/fitted peaks
            show_assigned: Show reference peaks
            **kwargs: Additional options:
                - show_ellipses: Show PS2D elliptical windows
                - nucleus_type: '15N' or '13C' for PS2D config
        """
        # Clear previous annotations
        for annotation in self.peak_annotations:
            try:
                annotation.remove()
            except:
                pass
        self.peak_annotations.clear()

        # The integrator should have a list of standardized peak dictionaries
        peak_list = getattr(integrator, 'fitted_peaks', [])
        if not peak_list:
            peak_list = getattr(integrator, 'integration_results', [])

        # Plot detected/fitted peaks with quality-based color coding
        if show_detected and peak_list:
            self._plot_detected_peaks(peak_list)

        # Plot reference peaks (blue crosses)
        if show_assigned and hasattr(integrator, 'peak_list') and integrator.peak_list is not None:
            self._plot_reference_peaks(integrator.peak_list)

        # Plot PS2D ellipses if requested (debug tool)
        show_ellipses = kwargs.get('show_ellipses', False)
        nucleus_type = kwargs.get('nucleus_type', '15N')
        if show_ellipses:
            self._plot_ps2d_ellipses(integrator, nucleus_type)

        # Add legend
        handles, labels = self.ax.get_legend_handles_labels()
        if handles:
            self.ax.legend(loc='upper right', fontsize='small')

        self.canvas.draw()

    def _plot_detected_peaks(self, peak_list):
        """
        Plot detected/fitted peaks with quality-based color coding.

        Args:
            peak_list: List of peak dictionaries
        """
        detected_x = []
        detected_y = []
        colors = []

        for p in peak_list:
            # Try multiple possible field names for X coordinate
            x = p.get('ppm_x') or p.get('x_ppm') or p.get('peak_x') or p.get('center_x')
            # Try multiple possible field names for Y coordinate
            y = p.get('ppm_y') or p.get('y_ppm') or p.get('peak_y') or p.get('center_y')

            if x is not None and y is not None:
                detected_x.append(float(x))
                detected_y.append(float(y))

                # Extract quality from multiple possible field names
                quality = p.get('quality') or p.get('Quality') or p.get('fitting_quality') or 'Unknown'

                # Map quality to color
                if quality in ['Excellent', 'excellent']:
                    colors.append('lime')
                elif quality in ['Good', 'good']:
                    colors.append(SUCCESS_GREEN)
                elif quality in ['Fair', 'fair']:
                    colors.append(WARNING_ORANGE)
                elif quality in ['Poor', 'poor']:
                    colors.append(ERROR_RED)
                else:
                    colors.append(ERROR_RED)  # No quality info = detected but not fitted

        if detected_x:
            # Plot with quality-based colors
            self.ax.scatter(detected_x, detected_y, c=colors, marker='o', s=60,
                          alpha=0.8, edgecolors='black', linewidth=1, zorder=5,
                          label=f'Detected ({len(detected_x)})')

            # Add annotations with matching colors
            for peak, color in zip(peak_list, colors):
                # Use flexible field name lookup for assignment
                assignment_label = (peak.get('assignment') or
                                  peak.get('Assignment') or
                                  str(peak.get('peak_number', '')))
                # Use flexible field name lookup for coordinates
                x = peak.get('ppm_x') or peak.get('x_ppm') or peak.get('peak_x') or peak.get('center_x')
                y = peak.get('ppm_y') or peak.get('y_ppm') or peak.get('peak_y') or peak.get('center_y')

                if x is not None and y is not None:
                    annotation = self.ax.annotate(
                        assignment_label,
                        (float(x), float(y)),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize='small', color=color, fontweight='bold',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8,
                                 edgecolor=color)
                    )
                    self.peak_annotations.append(annotation)

    def _plot_reference_peaks(self, peak_df):
        """
        Plot reference peaks (blue crosses).

        Args:
            peak_df: DataFrame with reference peak positions
        """
        try:
            # Find X coordinate column (1H dimension)
            x_col = None
            for col_name in ['Position_X', 'position_x', 'pos_x', 'x_ppm', 'chemical_shift_x', 'H1', '1H']:
                if col_name in peak_df.columns:
                    x_col = col_name
                    break

            # Find Y coordinate column (15N/13C dimension)
            y_col = None
            for col_name in ['Position_Y', 'position_y', 'pos_y', 'y_ppm', 'chemical_shift_y', 'N15', '15N', 'C13', '13C']:
                if col_name in peak_df.columns:
                    y_col = col_name
                    break

            if x_col and y_col:
                ref_x = [float(row[x_col]) for _, row in peak_df.iterrows()]
                ref_y = [float(row[y_col]) for _, row in peak_df.iterrows()]

                if ref_x and ref_y:
                    self.ax.scatter(ref_x, ref_y, c='blue', marker='x', s=40,
                                  alpha=0.6, linewidth=2, zorder=4,
                                  label=f'Reference ({len(ref_x)})')

                    # Add peak numbers as annotations
                    for idx, (x, y) in enumerate(zip(ref_x, ref_y)):
                        # Get assignment if available, otherwise use peak number
                        if 'Assignment' in peak_df.columns:
                            assignment = str(peak_df.iloc[idx]['Assignment'])
                        else:
                            assignment = str(idx + 1)

                        annotation = self.ax.annotate(
                            assignment,
                            (x, y),
                            xytext=(8, 8), textcoords='offset points',
                            fontsize='small', color='blue', fontweight='bold',
                            bbox=dict(boxstyle="round,pad=0.2", facecolor="lightblue", alpha=0.7),
                            zorder=5
                        )
                        self.peak_annotations.append(annotation)

                    print(f"✅ Plotted {len(ref_x)} reference peaks using columns '{x_col}' and '{y_col}'")
            else:
                available_cols = list(peak_df.columns) if hasattr(peak_df, 'columns') else []
                print(f"⚠️  Could not find position columns in peak list.")
                print(f"    Available columns: {available_cols}")
                print(f"    Expected: Position_X/Position_Y or similar variants")

        except Exception as e:
            print(f"❌ Error plotting reference peaks: {e}")
            import traceback
            traceback.print_exc()

    def _plot_ps2d_ellipses(self, integrator, nucleus_type='15N'):
        """
        Plot PS2D elliptical integration windows for all peaks (debug tool).

        Args:
            integrator: EnhancedVoigtIntegrator with peak_list
            nucleus_type: '15N' or '13C' for config lookup
        """
        from matplotlib.patches import Ellipse as EllipsePatch

        # Get PS2D configuration
        try:
            from lunaNMR.core.ps2d_config import get_ps2d_config
            config = get_ps2d_config()

            radF1 = config.radF1
            radF2 = config.radF2
            radF1_selector = config.radF1_selector
            radF2_selector = config.radF2_selector

        except Exception as e:
            print(f"⚠️  Could not load PS2D config: {e}")
            return

        # Get peak positions - try multiple sources
        peaks_to_plot = []

        # Try detected/fitted peaks (list of dicts) - HIGHEST PRIORITY
        if hasattr(integrator, 'fitted_peaks') and integrator.fitted_peaks:
            try:
                for peak in integrator.fitted_peaks:
                    peaks_to_plot.append({
                        'x': peak.get('ppm_x', peak.get('x_ppm', 0)),
                        'y': peak.get('ppm_y', peak.get('y_ppm', 0)),
                        'label': peak.get('assignment', peak.get('Assignment', 'Unknown'))
                    })
                print(f"✅ Using {len(peaks_to_plot)} detected peak positions for ellipses")
            except Exception as e:
                print(f"⚠️  Could not extract peaks from fitted_peaks: {e}")

        # Try reference peak list (DataFrame) - FALLBACK ONLY
        if not peaks_to_plot and hasattr(integrator, 'peak_list') and integrator.peak_list is not None:
            try:
                peak_df = integrator.peak_list

                # Find position columns
                x_col, y_col = None, None
                for col_name in ['Position_X', 'position_x', 'pos_x', 'x_ppm', '1H']:
                    if col_name in peak_df.columns:
                        x_col = col_name
                        break

                for col_name in ['Position_Y', 'position_y', 'pos_y', 'y_ppm', '15N', '13C']:
                    if col_name in peak_df.columns:
                        y_col = col_name
                        break

                if x_col and y_col:
                    for idx, row in peak_df.iterrows():
                        assignment = str(row.get('Assignment', idx + 1))
                        peaks_to_plot.append({
                            'x': float(row[x_col]),
                            'y': float(row[y_col]),
                            'label': assignment
                        })
                    print(f"✅ Using {len(peaks_to_plot)} reference peak positions for ellipses")
            except Exception as e:
                print(f"⚠️  Could not extract peaks from peak_list: {e}")

        # Plot ellipses
        if peaks_to_plot:
            for peak in peaks_to_plot:
                # Fitting region ellipse (red, dashed)
                fit_ellipse = EllipsePatch(
                    (peak['x'], peak['y']),
                    width=2 * radF2,
                    height=2 * radF1,
                    edgecolor='red',
                    facecolor='none',
                    linestyle='--',
                    linewidth=1.5,
                    alpha=0.7,
                    zorder=3
                )
                self.ax.add_patch(fit_ellipse)

                # Data selector ellipse (blue, dotted)
                selector_ellipse = EllipsePatch(
                    (peak['x'], peak['y']),
                    width=2 * radF2_selector,
                    height=2 * radF1_selector,
                    edgecolor='blue',
                    facecolor='none',
                    linestyle=':',
                    linewidth=1.0,
                    alpha=0.5,
                    zorder=2
                )
                self.ax.add_patch(selector_ellipse)

    def highlight_peak(self, peak_x, peak_y, highlight_color='yellow', radius=0.1):
        """
        Highlight a specific peak with a colored circle.

        Args:
            peak_x: X coordinate (ppm)
            peak_y: Y coordinate (ppm)
            highlight_color: Color for highlight circle
            radius: Radius of highlight circle (ppm units)
        """
        from matplotlib.patches import Circle

        circle = Circle((peak_x, peak_y), radius,
                       color=highlight_color, alpha=0.3, zorder=10)
        self.ax.add_patch(circle)
        self.canvas.draw()

    def add_crosshair(self, x, y):
        """
        Add crosshair lines at specified position.

        Args:
            x: X coordinate (ppm)
            y: Y coordinate (ppm)
        """
        # Remove old crosshairs
        for line in self.crosshair_lines:
            try:
                line.remove()
            except:
                pass
        self.crosshair_lines.clear()

        # Add new crosshairs
        vline = self.ax.axvline(x, color='red', linestyle='--', alpha=0.5, linewidth=1)
        hline = self.ax.axhline(y, color='red', linestyle='--', alpha=0.5, linewidth=1)
        self.crosshair_lines.extend([vline, hline])

        self.canvas.draw()

    def set_zoom(self, x_center, y_center, x_range, y_range):
        """
        Set zoom to specific region.

        Args:
            x_center: X coordinate center (ppm)
            y_center: Y coordinate center (ppm)
            x_range: X range width (ppm)
            y_range: Y range height (ppm)
        """
        self.ax.set_xlim(x_center - x_range/2, x_center + x_range/2)
        self.ax.set_ylim(y_center - y_range/2, y_center + y_range/2)
        # Ensure NMR convention: high ppm on left/bottom
        self._ensure_nmr_axis_orientation()
        self.canvas.draw()

    def reset_zoom(self):
        """Reset zoom to show full spectrum."""
        self.ax.autoscale()
        # Ensure NMR convention: high ppm on left/bottom
        self._ensure_nmr_axis_orientation()
        self.canvas.draw()

    def _ensure_nmr_axis_orientation(self):
        """Ensure axes follow NMR convention (high ppm on left/bottom).

        NMR spectra are conventionally displayed with:
        - X-axis (1H): high ppm on LEFT, low ppm on RIGHT
        - Y-axis (15N): high ppm at BOTTOM, low ppm at TOP
        """
        if not self.ax.xaxis_inverted():
            self.ax.invert_xaxis()
        if not self.ax.yaxis_inverted():
            self.ax.invert_yaxis()

    def clear(self):
        """Clear the plot and all annotations."""
        super().clear()
        self.peak_annotations.clear()
        self.crosshair_lines.clear()
        self.contour_plot = None
        self.colorbar = None
        self.colorbar_ax = None
        self.refresh()

    def update_view(self):
        """Refresh the display."""
        self.refresh()

    def _on_click(self, event):
        """
        Handle mouse click events.

        Args:
            event: Matplotlib mouse event
        """
        if event.inaxes == self.ax and event.xdata and event.ydata:
            # Emit signal with clicked coordinates
            self.peak_clicked.emit(event.xdata, event.ydata)
