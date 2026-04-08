#!/usr/bin/env python3
# ABOUTME: Qt VoigtAnalysisPlotter for 2D/3D Voigt fitting visualization in LunaNMR Qt interface
# ABOUTME: Supports layer toggling, color presets, 3D surface plots, cross-sections, and residual analysis

"""
VoigtAnalysisPlotter - Advanced Voigt Analysis Visualization for Qt

This module provides a comprehensive visualization widget for 2D Voigt fitting results,
supporting both 2D contour plots and 3D surface visualizations.

Key Features:
- Layer toggling (experimental, fitted, residuals, individual peaks)
- Color scheme presets (Classic, Clean, Dark, Warm)
- 3D surface plots with peak labels
- Cross-section views (F1/F2)
- Interactive click events
- Residual overlay vs separate panel modes
- Intensity scaling
- Peak clipping to fit region

Port Details:
- CustomTkinter → PySide6
- TkAgg → FigureCanvasQTAgg
- 2×2 axes grid using MatplotlibMultiAxesWidget
- All features from original implementation preserved

Author: Guillaume Mas
Date: 2025-01-24
"""

from PySide6.QtCore import Signal
import numpy as np

from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, ListedColormap, Normalize
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from scipy.special import wofz

from lunaNMR.gui.components.matplotlib_widget import MatplotlibMultiAxesWidget


def _to_scalar(value):
    """Convert value to scalar, handling lists from JSON deserialization."""
    if value is None:
        return 0.0
    if isinstance(value, (list, tuple)):
        return float(value[0]) if len(value) > 0 else 0.0
    return float(value)


class VoigtAnalysisPlotter(MatplotlibMultiAxesWidget):
    """
    Voigt fitting analysis visualization with 2×2 axes grid.

    Grid Layout:
        [[ax_x, ax_y], [ax_2d, ax_residuals]]
        But accessed via axes_list[0-3] mapping to positions:
        [0, 1, 2, 3] → [[0,1], [2,3]]

    Features:
        - Layer toggling (experimental/fitted/residuals/individual peaks)
        - Color presets for 3D plots
        - Cross-sections (F1/F2)
        - Residual modes (overlay/separate)
        - Intensity scaling
        - Peak clipping control
        - Interactive click events (via Qt signal)

    Signals:
        plot_clicked(float, float): Emitted when user clicks on plot with (x, y) coordinates
    """

    # Qt signal for click events
    plot_clicked = Signal(float, float)

    def __init__(self, parent=None, toolbar=True):
        """
        Initialize VoigtAnalysisPlotter with 2×2 grid.

        Args:
            parent: Parent Qt widget
            toolbar: Include navigation toolbar if True
        """
        # Initialize with 2×2 grid
        super().__init__(
            parent=parent,
            nrows=2,
            ncols=2,
            toolbar=toolbar,
            figsize=(12, 10),
            dpi=100
        )

        # Map axes_list indices to logical positions
        # axes_list[0] = top-left    (ax_x)
        # axes_list[1] = top-right   (ax_y)
        # axes_list[2] = bottom-left (ax_2d)
        # axes_list[3] = bottom-right (ax_residuals)

        # Feature 1: Layer toggling
        self.show_experimental = True
        self.show_fitted = True
        self.show_residuals = False  # Hidden by default
        self.show_individual_peaks = False  # Hidden by default
        self.show_peak_labels = True  # Toggle 3D peak assignment labels

        # Feature 2: Cross-sections (disabled by default, code kept for future use)
        self.show_cross_sections = False
        self.ax_f1_cross = None
        self.ax_f2_cross = None
        self.click_position = None

        # Feature 3: Residual mode (separate panel vs overlay)
        self.residual_mode = 'overlay'  # Default to overlay mode

        # Feature 4: Intensity scaling
        self.intensity_scale_factor = 1.0  # 100% by default
        self.auto_z_min = None
        self.auto_z_max = None

        # Store current result for refresh operations
        self.current_result = None
        self.ax_data = None
        self.ax_resid = None

        # Store event connection ID to prevent accumulation
        self.click_event_cid = None

        # Peak color mapping - ensures same peak gets same color across spectra
        self._peak_color_map = {}  # assignment -> color_index

        # Feature 5: Color scheme presets for 3D plots
        self.color_presets = {
            'Classic': {
                'background': 'black',
                'experimental': {'color': 'silver', 'alpha': 0.4, 'linewidth': 2.0},
                'total_fit': {'color': 'lightsteelblue', 'alpha': 0.4},
                'grid': {'color': '#333333', 'alpha': 0.4},
                'pane': {'facecolor': 'black', 'alpha': 0.9}
            },
            'Clean': {
                'background': 'white',
                'experimental': {'color': 'black', 'alpha': 0.4, 'linewidth': 2.0},
                'total_fit': {'color': 'dodgerblue', 'alpha': 0.4},
                'grid': {'color': 'lightgray', 'alpha': 0.5},
                'pane': {'facecolor': 'white', 'alpha': 0.95}
            },
            'Dark': {
                'background': '#2C3E50',
                'experimental': {'color': 'silver', 'alpha': 0.4, 'linewidth': 2.0},
                'total_fit': {'color': 'cyan', 'alpha': 0.4},
                'grid': {'color': '#34495E', 'alpha': 0.6},
                'pane': {'facecolor': '#2C3E50', 'alpha': 0.9}
            },
            'Warm': {
                'background': '#404040',
                'experimental': {'color': 'wheat', 'alpha': 0.4, 'linewidth': 2.0},
                'total_fit': {'color': 'skyblue', 'alpha': 0.4},
                'grid': {'color': '#505050', 'alpha': 0.5},
                'pane': {'facecolor': '#404040', 'alpha': 0.9}
            }
        }
        self.current_color_scheme = 'Clean'  # Default scheme

        # Feature 6: Individual peak clipping to fit region
        self.clip_individual_peaks = True  # ON by default
        self.peak_display_multiplier = 2.2  # Show peaks up to 2.2× fit ellipse radius

        # Feature 7: Fixed axis limits for Peak Mode comparison
        self.use_fixed_z_limits = False
        self.fixed_z_min = None
        self.fixed_z_max = None
        self.use_fixed_xy_limits = False
        self.fixed_x_min = None  # F2 (1H)
        self.fixed_x_max = None
        self.fixed_y_min = None  # F1 (15N)
        self.fixed_y_max = None

        # Feature 8: Preserve 3D view orientation between plots
        self.preserved_elev = 25  # Default elevation
        self.preserved_azim = 45  # Default azimuth
        self.preserved_roll = 0   # Default roll

        # Initialize with placeholder message
        self._plot_no_data()

    # ===== Feature Control Methods =====

    def toggle_experimental(self, visible):
        """Toggle experimental data layer visibility"""
        self.show_experimental = visible
        self.refresh_plot()

    def toggle_fitted(self, visible):
        """Toggle fitted peaks layer visibility"""
        self.show_fitted = visible
        self.refresh_plot()

    def toggle_individual_peaks(self, visible):
        """Toggle individual peak surfaces visibility"""
        self.show_individual_peaks = visible
        self.refresh_plot()

    def toggle_peak_labels(self, visible):
        """Toggle 3D peak assignment labels visibility"""
        self.show_peak_labels = visible
        self.refresh_plot()

    def toggle_residuals(self, visible):
        """Toggle residuals layer visibility"""
        self.show_residuals = visible
        self.refresh_plot()

    def toggle_peak_clipping(self, enabled):
        """Toggle individual peak clipping to fit region

        Args:
            enabled: True to limit peaks to 2.5× fit ellipse radius, False for full extent
        """
        self.clip_individual_peaks = enabled
        self.refresh_plot()

    def _get_peak_color(self, assignment: str, colors: list) -> str:
        """Get consistent color for a peak based on its assignment.

        Assigns colors in order of first appearance. Same assignment always
        gets the same color across different spectra within this viewer.

        Args:
            assignment: Peak assignment string (e.g., 'A.142.ALA.H')
            colors: List of available colors

        Returns:
            Color string from the colors list
        """
        if assignment not in self._peak_color_map:
            # Assign next available color index
            self._peak_color_map[assignment] = len(self._peak_color_map)
        color_idx = self._peak_color_map[assignment] % len(colors)
        return colors[color_idx]

    def set_residual_mode(self, mode):
        """Set residual visualization mode

        Args:
            mode: 'separate' or 'overlay'
        """
        if mode not in ['separate', 'overlay']:
            raise ValueError(f"Invalid residual mode: {mode}. Must be 'separate' or 'overlay'")
        self.residual_mode = mode
        self.refresh_plot()

    def set_intensity_scale(self, scale_percent):
        """Set intensity scale factor from percentage (50-200%)

        Args:
            scale_percent: Scale as percentage (100 = normal, 50 = half, 200 = double)
        """
        self.intensity_scale_factor = scale_percent / 100.0
        self.apply_intensity_scale(self.intensity_scale_factor)

    def set_color_scheme(self, scheme_name):
        """Set color scheme for 3D plots

        Args:
            scheme_name: One of 'Classic', 'Clean', 'Dark', 'Warm'
        """
        if scheme_name in self.color_presets:
            self.current_color_scheme = scheme_name
            self.refresh_plot()

    def apply_intensity_scale(self, factor):
        """Apply intensity scaling to Z-axis limits"""
        if not hasattr(self, 'ax_data') or self.ax_data is None:
            return  # No plot yet

        if self.auto_z_min is None or self.auto_z_max is None:
            return  # No auto-limits stored yet

        # Calculate scaled limits
        z_center = (self.auto_z_max + self.auto_z_min) / 2
        z_half_range = (self.auto_z_max - self.auto_z_min) / 2

        scaled_half_range = z_half_range * factor
        new_z_min = z_center - scaled_half_range
        new_z_max = z_center + scaled_half_range

        # Update data panel
        self.ax_data.set_zlim(new_z_min, new_z_max)

        # Update residuals panel if in separate mode
        if self.residual_mode == 'separate' and self.ax_resid is not None:
            self.ax_resid.set_zlim(new_z_min, new_z_max)

        # Redraw canvas
        self.canvas.draw_idle()

    def set_fixed_z_limits(self, z_min: float, z_max: float):
        """Set fixed z-axis limits for Peak Mode comparison.

        When set, the z-axis will use these limits instead of auto-scaling
        to the current data. This allows visual comparison of intensity
        changes across spectra.

        Args:
            z_min: Minimum z-axis value
            z_max: Maximum z-axis value
        """
        self.fixed_z_min = z_min
        self.fixed_z_max = z_max
        self.use_fixed_z_limits = True

    def clear_fixed_z_limits(self):
        """Clear fixed z-limits and revert to auto-scaling."""
        self.use_fixed_z_limits = False
        self.fixed_z_min = None
        self.fixed_z_max = None

    def set_fixed_xy_limits(self, x_min: float, x_max: float, y_min: float, y_max: float):
        """Set fixed x/y axis limits for Peak Mode comparison.

        Args:
            x_min: Minimum x-axis (F2/1H) value in ppm
            x_max: Maximum x-axis (F2/1H) value in ppm
            y_min: Minimum y-axis (F1/15N) value in ppm
            y_max: Maximum y-axis (F1/15N) value in ppm
        """
        self.fixed_x_min = x_min
        self.fixed_x_max = x_max
        self.fixed_y_min = y_min
        self.fixed_y_max = y_max
        self.use_fixed_xy_limits = True

    def clear_fixed_xy_limits(self):
        """Clear fixed x/y limits and revert to auto-scaling."""
        self.use_fixed_xy_limits = False
        self.fixed_x_min = None
        self.fixed_x_max = None
        self.fixed_y_min = None
        self.fixed_y_max = None

    def reset_view_orientation(self):
        """Reset 3D view orientation to defaults."""
        self.preserved_elev = 25
        self.preserved_azim = 45
        self.preserved_roll = 0

    def refresh_plot(self):
        """Redraw current plot with new settings"""
        if self.current_result is not None:
            self.plot_voigt_analysis_3d(self.current_result)
            self.canvas.draw()

    def show_placeholder(self, message: str = "Select a peak to view Voigt analysis"):
        """Show placeholder message in all axes.

        Args:
            message: Text message to display
        """
        # Clear all axes
        for ax in self.axes_list:
            ax.clear()
            ax.text(0.5, 0.5, message,
                   ha='center', va='center',
                   transform=ax.transAxes,
                   fontsize=11, color='gray',
                   style='italic')
            ax.set_title('')
            ax.axis('off')

        self.current_result = None
        self.canvas.draw_idle()

    def update_cross_sections(self, x_ppm, y_ppm, result):
        """Update F1 and F2 cross-section plots

        Args:
            x_ppm: F2 (1H) position in ppm
            y_ppm: F1 (15N/13C) position in ppm
            result: Voigt analysis result dictionary
        """
        if not self.show_cross_sections:
            return

        if self.ax_f1_cross is None or self.ax_f2_cross is None:
            return  # Cross-section axes not created yet

        region_2d = result.get('region_2d')
        if region_2d is None:
            return

        # Convert arrays from lists to numpy arrays (may be lists from JSON serialization)
        for key in ['f1_ppm', 'f2_ppm', 'f1_grid', 'f2_grid', 'intensity']:
            if key in region_2d and isinstance(region_2d[key], list):
                region_2d[key] = np.array(region_2d[key])

        f1_ppm = region_2d['f1_ppm']  # (M,)
        f2_ppm = region_2d['f2_ppm']  # (N,)
        experimental = region_2d['intensity']  # (M, N)
        fitted = result.get('fitted_2d_surface')  # (M, N)

        if fitted is not None and isinstance(fitted, list):
            fitted = np.array(fitted)

        if fitted is None:
            return

        # Find nearest indices
        f2_idx = np.argmin(np.abs(f2_ppm - x_ppm))
        f1_idx = np.argmin(np.abs(f1_ppm - y_ppm))

        # Extract 1D slices
        f1_exp_slice = experimental[:, f2_idx]  # Vertical slice at x_ppm
        f1_fit_slice = fitted[:, f2_idx]

        f2_exp_slice = experimental[f1_idx, :]  # Horizontal slice at y_ppm
        f2_fit_slice = fitted[f1_idx, :]

        # --- Plot F1 cross-section ---
        self.ax_f1_cross.clear()
        self.ax_f1_cross.plot(f1_ppm, f1_exp_slice, 'o',
                              color='gray', markersize=3, label='Exp')
        self.ax_f1_cross.plot(f1_ppm, f1_fit_slice, '-',
                              color='red', linewidth=2, label='Fit')

        # Add individual peak components if available
        individual_surfaces = result.get('individual_surfaces', [])
        all_peaks = result.get('all_peaks', [])
        if individual_surfaces:
            colors = ['orange', 'purple', 'brown', 'pink', 'olive']
            for i, surf in enumerate(individual_surfaces):
                component_slice = surf[:, f2_idx]
                # Use assignment-based color for consistency across spectra
                peak_assignment = all_peaks[i].get('assignment', f'Peak {i+1}') if i < len(all_peaks) else f'Peak {i+1}'
                color = self._get_peak_color(peak_assignment, colors)
                self.ax_f1_cross.plot(f1_ppm, component_slice, '--',
                                      color=color, linewidth=1, alpha=0.7,
                                      label=peak_assignment)

        self.ax_f1_cross.set_xlabel('F1 (ppm)', fontsize=8)
        self.ax_f1_cross.set_ylabel('Intensity', fontsize=8)
        self.ax_f1_cross.set_title(f'F1 Cross-Section at F2={x_ppm:.3f} ppm', fontsize=9)
        self.ax_f1_cross.invert_xaxis()  # NMR convention
        self.ax_f1_cross.legend(fontsize=6, loc='best')
        self.ax_f1_cross.grid(True, alpha=0.3)

        # --- Plot F2 cross-section ---
        self.ax_f2_cross.clear()
        self.ax_f2_cross.plot(f2_ppm, f2_exp_slice, 'o',
                              color='gray', markersize=3, label='Exp')
        self.ax_f2_cross.plot(f2_ppm, f2_fit_slice, '-',
                              color='blue', linewidth=2, label='Fit')

        if individual_surfaces:
            colors = ['orange', 'purple', 'brown', 'pink', 'olive']
            for i, surf in enumerate(individual_surfaces):
                component_slice = surf[f1_idx, :]
                # Use assignment-based color for consistency across spectra
                peak_assignment = all_peaks[i].get('assignment', f'Peak {i+1}') if i < len(all_peaks) else f'Peak {i+1}'
                color = self._get_peak_color(peak_assignment, colors)
                self.ax_f2_cross.plot(f2_ppm, component_slice, '--',
                                      color=color, linewidth=1, alpha=0.7,
                                      label=peak_assignment)

        self.ax_f2_cross.set_xlabel('F2 (ppm)', fontsize=8)
        self.ax_f2_cross.set_ylabel('Intensity', fontsize=8)
        self.ax_f2_cross.set_title(f'F2 Cross-Section at F1={y_ppm:.2f} ppm', fontsize=9)
        self.ax_f2_cross.invert_xaxis()
        self.ax_f2_cross.legend(fontsize=6, loc='best')
        self.ax_f2_cross.grid(True, alpha=0.3)

        # Redraw canvas
        self.canvas.draw_idle()

    def plot_voigt_analysis(self, voigt_result):
        """Plot comprehensive Voigt analysis results (delegates to specialized methods)"""
        if not voigt_result:
            self._plot_no_data()
            return

        # Extract data
        method = voigt_result.get('method', '')

        # Check if this is 2D simultaneous fitting
        if '2d_simultaneous' in method and 'region_2d' in voigt_result:
            # Use 2D contour visualization
            self._plot_2d_contour(voigt_result)
        else:
            # Use traditional 1D cross-section visualization
            self._plot_1d_fits(voigt_result)

    def plot_voigt_analysis_3d(self, voigt_result):
        """
        Plot 3D surface visualization of Voigt analysis results

        This is the main entry point for 3D visualization
        """
        # Clear figure at start
        self.figure.clear()

        if not voigt_result:
            self._plot_no_data_3d()
            return

        # Extract data
        method = voigt_result.get('method', '')

        # Check if this is 2D simultaneous fitting
        if '2d_simultaneous' in method and 'region_2d' in voigt_result:
            # Use 3D surface visualization for 2D simultaneous fits
            self._plot_2d_surface_3d(voigt_result)
        else:
            # For 1D cross-sections, show message that 3D view is only for 2D fits
            self._plot_1d_not_applicable_3d(voigt_result)

    def _plot_1d_fits(self, voigt_result):
        """Plot 1D cross-section fitting results"""
        # Clear and recreate 2×1 layout for 1D fits
        self.figure.clear()

        ax_x = self.figure.add_subplot(2, 1, 1)
        ax_y = self.figure.add_subplot(2, 1, 2)

        x_fit = voigt_result.get('x_fit', {})
        y_fit = voigt_result.get('y_fit', {})
        assignment = voigt_result.get('assignment', 'Unknown')
        quality = voigt_result.get('fitting_quality', 'Unknown')

        # Plot X-dimension fit (top)
        self._plot_1d_fit(ax_x, x_fit, '¹H Chemical Shift (ppm)',
                         f'X-Dimension Fit - {assignment}')

        # Plot Y-dimension fit (bottom)
        self._plot_1d_fit(ax_y, y_fit, '¹⁵N/¹³C Chemical Shift (ppm)',
                         f'Y-Dimension Fit - {assignment}')

        # Set overall title
        self.figure.suptitle(f'Voigt Profile Analysis: {assignment} (Quality: {quality})',
                            fontsize=10, fontweight='bold')

        self.figure.tight_layout(rect=[0, 0, 1, 0.96])
        self.canvas.draw()

    def _plot_1d_fit(self, ax, fit_data, xlabel, title):
        """Plot 1D fitting results"""
        if not fit_data or 'ppm_scale' not in fit_data:
            ax.text(0.5, 0.5, 'No fitting data available',
                   transform=ax.transAxes, ha='center', va='center')
            ax.set_title(title)
            return

        ppm_scale = fit_data['ppm_scale']
        cross_section = fit_data['cross_section']
        fitted_curve = fit_data.get('fitted_curve', None)
        r_squared = fit_data.get('r_squared_local', fit_data.get('r_squared', 0))

        # Plot experimental data
        ax.plot(ppm_scale, cross_section, 'b-', linewidth=2, alpha=0.7, label='Experimental')

        # Plot fitted curve
        if fitted_curve is not None:
            ax.plot(ppm_scale, fitted_curve, 'r-', linewidth=2, label='Voigt Fit')

            # Calculate and plot residuals
            residuals = cross_section - fitted_curve
            ax.plot(ppm_scale, residuals + np.min(cross_section)*0.9, 'g--',
                   alpha=0.6, label='Residuals')

        ax.set_xlabel(xlabel)
        ax.set_ylabel('Intensity')

        # Show window size information in title if available
        window_info = ""
        if 'window_size' in fit_data:
            window_size = fit_data['window_size']
            gui_based = fit_data.get('gui_based', False)
            window_source = "GUI" if gui_based else "Auto"
            window_info = f", Window: ±{window_size/2:.3f} ppm ({window_source})"

        ax.set_title(f'{title} (R² = {r_squared:.3f}{window_info})')
        ax.legend()
        ax.grid(True, alpha=0.3)

        if '¹H' in xlabel:
            ax.invert_xaxis()

    def _plot_2d_contour(self, voigt_result):
        """
        Plot 2D contour visualization for simultaneous multi-peak fitting

        Layout (2x2 grid):
        [0,0] Experimental 2D Data    [0,1] Fitted 2D Surface
        [1,0] 2D Residuals            [1,1] Parameters Table
        """
        # Extract data
        region_2d = voigt_result.get('region_2d')
        fitted_surface = voigt_result.get('fitted_2d_surface')
        all_peaks = voigt_result.get('all_peaks', [])
        r_squared = voigt_result.get('avg_r_squared', 0)
        assignment = voigt_result.get('assignment', 'Unknown')
        quality = voigt_result.get('fitting_quality', 'Unknown')

        # Convert arrays from lists to numpy arrays (may be lists from JSON serialization)
        if region_2d is not None:
            for key in ['f1_ppm', 'f2_ppm', 'f1_grid', 'f2_grid', 'intensity']:
                if key in region_2d and isinstance(region_2d[key], list):
                    region_2d[key] = np.array(region_2d[key])

        if fitted_surface is not None and isinstance(fitted_surface, list):
            fitted_surface = np.array(fitted_surface)

        if region_2d is None or fitted_surface is None:
            self._plot_no_data()
            return

        f1_ppm = region_2d['f1_ppm']
        f2_ppm = region_2d['f2_ppm']
        experimental = region_2d['intensity']

        # Calculate residuals
        residuals = experimental - fitted_surface

        # Determine contour levels
        vmin = min(np.min(experimental), np.min(fitted_surface))
        vmax = max(np.max(experimental), np.max(fitted_surface))
        levels = np.linspace(vmin + (vmax - vmin) * 0.2, vmax, 20)

        # Clear and create 2x2 grid
        self.figure.clear()
        gs = self.figure.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
        ax_exp = self.figure.add_subplot(gs[0, 0])
        ax_fit = self.figure.add_subplot(gs[0, 1])
        ax_res = self.figure.add_subplot(gs[1, 0])
        ax_table = self.figure.add_subplot(gs[1, 1])

        # Plot 1: Experimental Data (top left) - BLUE filled contours
        ax_exp.contourf(f2_ppm, f1_ppm, experimental, levels=levels,
                        cmap='Blues', alpha=0.8)
        ax_exp.contour(f2_ppm, f1_ppm, experimental, levels=levels,
                       colors='blue', linewidths=0.5, alpha=0.5)
        ax_exp.set_xlabel('¹H Chemical Shift (ppm)', fontsize=9)
        ax_exp.set_ylabel('¹⁵N/¹³C Chemical Shift (ppm)', fontsize=9)
        ax_exp.set_title(f'Experimental - {assignment}', fontsize=9, fontweight='bold')
        ax_exp.invert_xaxis()
        ax_exp.invert_yaxis()
        ax_exp.grid(True, alpha=0.3)

        # Mark peak positions
        for i, peak in enumerate(all_peaks):
            pos_f2 = _to_scalar(peak['pos_f2'])
            pos_f1 = _to_scalar(peak['pos_f1'])
            ax_exp.plot(pos_f2, pos_f1, 'r+', markersize=12, markeredgewidth=2)
            peak_label = peak.get('assignment', str(i+1))
            ax_exp.text(pos_f2, pos_f1, f"  {peak_label}",
                       fontsize=8, color='red', fontweight='bold', va='center')

        # Plot 2: Individual Fitted Peaks (top right) - Multi-color visualization
        individual_surfaces = voigt_result.get('individual_surfaces', None)
        colors = ['red', 'orange', 'purple', 'brown', 'pink', 'olive', 'cyan', 'magenta']

        # Convert individual surfaces from lists to numpy arrays (may be lists from JSON)
        if individual_surfaces is not None:
            individual_surfaces = [
                np.array(s) if isinstance(s, list) else s
                for s in individual_surfaces
            ]

        if individual_surfaces is not None and len(individual_surfaces) > 0:
            # Calculate global contour levels
            baseline = region_2d['intensity'].min()
            all_peaks_with_baseline = [surf + baseline for surf in individual_surfaces]

            global_min = min(np.min(surf) for surf in all_peaks_with_baseline)
            global_max = max(np.max(surf) for surf in all_peaks_with_baseline)
            global_levels = np.linspace(global_min + (global_max - global_min) * 0.4, global_max, 15)

            # Plot each peak with same global levels
            for i, (peak_with_baseline, peak) in enumerate(zip(all_peaks_with_baseline, all_peaks)):
                # Use assignment-based color for consistency across spectra
                peak_assignment = peak.get('assignment', f'Peak {i+1}')
                color = self._get_peak_color(peak_assignment, colors)
                ax_fit.contour(f2_ppm, f1_ppm, peak_with_baseline, levels=global_levels,
                             colors=color, linewidths=1.2, alpha=0.7)

            ax_fit.set_title(f'Individual Fitted Peaks (R²={r_squared:.3f})', fontsize=9, fontweight='bold')
        else:
            # Fallback: plot summed surface
            ax_fit.contourf(f2_ppm, f1_ppm, fitted_surface, levels=levels,
                           cmap='Reds', alpha=0.8)
            ax_fit.contour(f2_ppm, f1_ppm, fitted_surface, levels=levels,
                          colors='red', linewidths=0.5, alpha=0.5)
            ax_fit.set_title(f'Fitted (R²={r_squared:.3f})', fontsize=9, fontweight='bold')

        ax_fit.set_xlabel('¹H Chemical Shift (ppm)', fontsize=9)
        ax_fit.set_ylabel('¹⁵N/¹³C Chemical Shift (ppm)', fontsize=9)
        ax_fit.invert_xaxis()
        ax_fit.invert_yaxis()
        ax_fit.grid(True, alpha=0.3)

        # Mark peak positions with matching colors
        for i, peak in enumerate(all_peaks):
            peak_assignment = peak.get('assignment', f'Peak {i+1}')
            if individual_surfaces is not None and len(individual_surfaces) > 0:
                # Use assignment-based color for consistency across spectra
                color = self._get_peak_color(peak_assignment, colors)
            else:
                color = 'blue'
            pos_f2 = _to_scalar(peak['pos_f2'])
            pos_f1 = _to_scalar(peak['pos_f1'])
            ax_fit.plot(pos_f2, pos_f1, '+', color=color,
                       markersize=12, markeredgewidth=2)
            peak_label = peak_assignment
            ax_fit.text(pos_f2, pos_f1, f"  {peak_label}",
                       fontsize=8, color=color, fontweight='bold', va='center')

        # Plot 3: Residuals (bottom left) - Color-coded quality regions
        relative_residuals = 100 * residuals / np.max(experimental)
        quality_map = np.abs(relative_residuals)

        # Define custom colormap
        colors_quality = ['#2ecc71', '#f1c40f', '#e74c3c']  # green, yellow, red
        boundaries = [0, 10, 20, 100]
        cmap_quality = ListedColormap(colors_quality)
        norm = BoundaryNorm(boundaries, cmap_quality.N)

        ax_res.contourf(f2_ppm, f1_ppm, quality_map,
                       levels=boundaries, cmap=cmap_quality, norm=norm, alpha=0.7)

        # Add peak positions for reference
        for i, peak in enumerate(all_peaks):
            ax_res.plot(_to_scalar(peak['pos_f2']), _to_scalar(peak['pos_f1']), 'k+', markersize=8, markeredgewidth=1.5)

        ax_res.set_xlabel('¹H Chemical Shift (ppm)', fontsize=9)
        ax_res.set_ylabel('¹⁵N/¹³C Chemical Shift (ppm)', fontsize=9)
        ax_res.set_title('Fit Quality (Residual %)', fontsize=9, fontweight='bold')
        ax_res.invert_xaxis()
        ax_res.invert_yaxis()
        ax_res.grid(True, alpha=0.3)

        # Add legend for quality regions
        legend_elements = [
            Patch(facecolor='#2ecc71', alpha=0.7, label='Good (<10%)'),
            Patch(facecolor='#f1c40f', alpha=0.7, label='Moderate (10-20%)'),
            Patch(facecolor='#e74c3c', alpha=0.7, label='Poor (>20%)')
        ]
        ax_res.legend(handles=legend_elements, loc='upper right', fontsize=7, framealpha=0.9)

        # Plot 4: Parameters Table (bottom right)
        ax_table.axis('off')

        # Create table text
        table_text = "Fitted Parameters\n" + "="*30 + "\n\n"

        for i, peak in enumerate(all_peaks):
            lw_f2 = _to_scalar(peak['lw_gau_f2']) + _to_scalar(peak['lw_lor_f2'])
            lw_f1 = _to_scalar(peak['lw_gau_f1']) + _to_scalar(peak['lw_lor_f1'])

            colors_list = ['red', 'orange', 'purple', 'brown', 'pink', 'olive', 'cyan', 'magenta']
            color = colors_list[i % len(colors_list)]

            peak_label = peak.get('assignment', f'Peak {i+1}')
            volume = _to_scalar(peak.get('volume', peak.get('intensity', 0.0)))
            height = _to_scalar(peak.get('height', peak.get('amplitude', 0.0)))

            table_text += f"{peak_label}\n"
            table_text += f"LW ¹H: {lw_f2:.4f}  LW ¹⁵N: {lw_f1:.3f}\n"
            table_text += f"V: {volume:.2e}  I: {height:.2e}\n"
            if i < len(all_peaks) - 1:
                table_text += "\n"

        ax_table.text(0.1, 0.85, table_text, transform=ax_table.transAxes,
                     fontsize=8, fontfamily='monospace', verticalalignment='top',
                     bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))

        # Add overall statistics at top
        stats_text = f"Quality: {quality}\nR² = {r_squared:.4f}\n{len(all_peaks)} peaks (simultaneous fit)"
        ax_table.text(0.5, 0.98, stats_text, ha='center', va='top',
                     fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Set overall title
        self.figure.suptitle(f'2D Simultaneous Voigt Fitting: {assignment}',
                            fontsize=10, fontweight='bold')

        self.figure.tight_layout(rect=[0, 0, 1, 0.96])
        self.canvas.draw()

    def _plot_2d_surface_3d(self, voigt_result):
        """
        Plot 3D surface visualization for 2D simultaneous multi-peak fitting

        Features:
        - Layer toggling (experimental/fitted/residuals)
        - Cross-sections (F1 and F2 slices)
        - Residual mode (separate panel or overlay)
        - Intensity scaling
        """
        # Store result for refresh operations
        self.current_result = voigt_result

        # Extract 2D data
        region_2d = voigt_result.get('region_2d')
        fitted_surface_original = voigt_result.get('fitted_2d_surface')
        individual_surfaces_original = voigt_result.get('individual_surfaces', None)
        all_peaks = voigt_result.get('all_peaks', [])
        r_squared = voigt_result.get('avg_r_squared', 0)
        baseline = voigt_result.get('baseline', 0.0)
        assignment = voigt_result.get('assignment', 'Unknown')

        # Convert arrays from lists to numpy arrays (may be lists from JSON serialization)
        if region_2d is not None:
            for key in ['f1_ppm', 'f2_ppm', 'f1_grid', 'f2_grid', 'intensity']:
                if key in region_2d and isinstance(region_2d[key], list):
                    region_2d[key] = np.array(region_2d[key])

        if fitted_surface_original is not None and isinstance(fitted_surface_original, list):
            fitted_surface_original = np.array(fitted_surface_original)

        if individual_surfaces_original is not None:
            individual_surfaces_original = [
                np.array(s) if isinstance(s, list) else s
                for s in individual_surfaces_original
            ]

        if region_2d is None or fitted_surface_original is None:
            self._plot_no_data_3d()
            return

        fitted_surface = fitted_surface_original

        # Re-reconstruct individual surfaces with clipping if enabled
        if self.clip_individual_peaks and individual_surfaces_original is not None and len(all_peaks) > 0:
            from lunaNMR.core.ps2d_config import get_ps2d_config

            config = get_ps2d_config()
            display_multiplier = self.peak_display_multiplier
            radF1_display = config.radF1 * display_multiplier
            radF2_display = config.radF2 * display_multiplier

            # Reconstruct with clipping
            individual_surfaces = []
            for peak in all_peaks:
                fitted_intensity = _to_scalar(peak.get('intensity', peak.get('amplitude', 1000)))

                SQRT_8LN2 = np.sqrt(8 * np.log(2))
                sigma_f1 = _to_scalar(peak['lw_gau_f1']) / SQRT_8LN2
                sigma_f2 = _to_scalar(peak['lw_gau_f2']) / SQRT_8LN2
                gamma_f1 = _to_scalar(peak['lw_lor_f1']) / 2.0
                gamma_f2 = _to_scalar(peak['lw_lor_f2']) / 2.0

                sigma_f1 = max(sigma_f1, 1e-10)
                sigma_f2 = max(sigma_f2, 1e-10)

                pos_f1 = _to_scalar(peak['pos_f1'])
                pos_f2 = _to_scalar(peak['pos_f2'])

                SQRT_2 = np.sqrt(2)
                z_f1 = ((pos_f1 - region_2d['f1_grid']) + 1j * gamma_f1) / (sigma_f1 * SQRT_2)
                z_f2 = ((pos_f2 - region_2d['f2_grid']) + 1j * gamma_f2) / (sigma_f2 * SQRT_2)

                fade_f1 = np.real(wofz(z_f1))
                fade_f2 = np.real(wofz(z_f2))

                peak_surface = fitted_intensity * fade_f1 * fade_f2 / (sigma_f1 * sigma_f2 * 2.0 * np.pi)

                # Apply elliptical clipping
                ellipse_distance_sq = ((region_2d['f1_grid'] - pos_f1) / radF1_display) ** 2 + \
                                      ((region_2d['f2_grid'] - pos_f2) / radF2_display) ** 2
                ellipse_mask = ellipse_distance_sq <= 1.0

                peak_surface = np.where(ellipse_mask, peak_surface, 0.0)
                individual_surfaces.append(peak_surface)
        else:
            individual_surfaces = individual_surfaces_original

        f1_ppm = region_2d['f1_ppm']
        f2_ppm = region_2d['f2_ppm']
        experimental = region_2d['intensity']
        residuals = experimental - (fitted_surface + baseline)

        # Create 2D meshgrids for 3D plotting
        F2_ppm, F1_ppm = np.meshgrid(f2_ppm, f1_ppm)

        # Save current view orientation before clearing
        if self.ax_data is not None:
            self.preserved_elev = self.ax_data.elev
            self.preserved_azim = self.ax_data.azim
            if hasattr(self.ax_data, 'roll'):
                self.preserved_roll = self.ax_data.roll

        # Clear figure completely
        self.figure.clear()
        self.canvas.draw_idle()

        # Determine layout based on features
        if self.show_cross_sections:
            if self.residual_mode == 'separate':
                gs = self.figure.add_gridspec(3, 2, height_ratios=[2, 2, 1.5], hspace=0.3, wspace=0.2)
            else:
                gs = self.figure.add_gridspec(2, 2, height_ratios=[2, 1.5], hspace=0.3, wspace=0.2)
        else:
            if self.residual_mode == 'separate':
                gs = self.figure.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.3)
            else:
                gs = self.figure.add_gridspec(1, 1)

        # Create axes
        ax_data = self.figure.add_subplot(
            gs[0, 0] if self.show_cross_sections or self.residual_mode == 'separate' else gs[0],
            projection='3d'
        )
        self.ax_data = ax_data

        # Get current color scheme
        scheme = self.color_presets[self.current_color_scheme]

        # Set 3D background pane colors from preset
        ax_data.xaxis.pane.set_facecolor(scheme['pane']['facecolor'])
        ax_data.yaxis.pane.set_facecolor(scheme['pane']['facecolor'])
        ax_data.zaxis.pane.set_facecolor(scheme['pane']['facecolor'])
        ax_data.xaxis.pane.set_alpha(scheme['pane']['alpha'])
        ax_data.yaxis.pane.set_alpha(scheme['pane']['alpha'])
        ax_data.zaxis.pane.set_alpha(scheme['pane']['alpha'])

        self.ax_resid = None

        if self.residual_mode == 'separate' and self.show_residuals:
            ax_resid = self.figure.add_subplot(gs[1, 0], projection='3d')
            self.ax_resid = ax_resid

            # Set 3D background pane colors for residuals panel
            ax_resid.xaxis.pane.set_facecolor(scheme['pane']['facecolor'])
            ax_resid.yaxis.pane.set_facecolor(scheme['pane']['facecolor'])
            ax_resid.zaxis.pane.set_facecolor(scheme['pane']['facecolor'])
            ax_resid.xaxis.pane.set_alpha(scheme['pane']['alpha'])
            ax_resid.yaxis.pane.set_alpha(scheme['pane']['alpha'])
            ax_resid.zaxis.pane.set_alpha(scheme['pane']['alpha'])

        if self.show_cross_sections:
            if self.residual_mode == 'separate':
                self.ax_f1_cross = self.figure.add_subplot(gs[2, 0])
                self.ax_f2_cross = self.figure.add_subplot(gs[2, 1])
            else:
                self.ax_f1_cross = self.figure.add_subplot(gs[1, 0])
                self.ax_f2_cross = self.figure.add_subplot(gs[1, 1])

        # Panel 1: Experimental data + Individual fitted peaks
        colors = ['red', 'orange', 'purple', 'brown', 'pink', 'olive', 'cyan', 'magenta']

        # Adaptive alpha scaling for individual peaks
        n_peaks = len(individual_surfaces) if individual_surfaces is not None else 0
        if n_peaks <= 2:
            individual_alpha = 0.9
        elif n_peaks <= 5:
            individual_alpha = 0.8
        else:
            individual_alpha = 0.7

        # Plot experimental data wireframe if enabled
        if self.show_experimental:
            ax_data.plot_wireframe(F2_ppm, F1_ppm, experimental,
                                   color=scheme['experimental']['color'],
                                   alpha=scheme['experimental']['alpha'],
                                   linewidth=scheme['experimental']['linewidth'],
                                   label='Experimental')

        # Plot total fitted surface as wireframe if enabled
        if self.show_fitted and fitted_surface is not None:
            ax_data.plot_wireframe(F2_ppm, F1_ppm, fitted_surface,
                                   color=scheme['total_fit']['color'],
                                   alpha=scheme['experimental']['alpha'],
                                   linewidth=scheme['experimental']['linewidth'],
                                   label='Total Fit')

        # Overlay individual fitted peaks if enabled
        if self.show_individual_peaks and individual_surfaces is not None and len(individual_surfaces) > 0:
            for i, (surf, peak) in enumerate(zip(individual_surfaces, all_peaks)):
                # Use assignment-based color for consistency across spectra
                peak_assignment = peak.get('assignment', f'Peak {i+1}')
                color = self._get_peak_color(peak_assignment, colors)

                # If clipping is enabled, mask where surface is zero
                if self.clip_individual_peaks:
                    threshold = np.max(surf) * 1e-6
                    surf_masked = np.where(surf > threshold, surf, np.nan)
                else:
                    surf_masked = surf

                ax_data.plot_wireframe(F2_ppm, F1_ppm, surf_masked,
                                       color=color, alpha=individual_alpha, linewidth=1.5,
                                       label=peak_assignment)

        # Add elliptical contour overlay showing PS2D fit regions
        if self.show_individual_peaks and all_peaks is not None and len(all_peaks) > 0:
            from lunaNMR.core.ps2d_config import get_ps2d_config

            config = get_ps2d_config()
            radF1 = config.radF1
            radF2 = config.radF2

            theta = np.linspace(0, 2*np.pi, 100)

            for i, peak in enumerate(all_peaks):
                # Use assignment-based color for consistency across spectra
                peak_assignment = peak.get('assignment', f'Peak {i+1}')
                color = self._get_peak_color(peak_assignment, colors)

                pos_f2 = _to_scalar(peak['pos_f2'])
                pos_f1 = _to_scalar(peak['pos_f1'])

                ellipse_f2 = pos_f2 + radF2 * np.cos(theta)
                ellipse_f1 = pos_f1 + radF1 * np.sin(theta)
                ellipse_z = np.zeros_like(theta)

                ax_data.plot(ellipse_f2, ellipse_f1, ellipse_z,
                            color=color, linewidth=2.5, alpha=0.7,
                            linestyle='--', zorder=10)

        # Add 3D text labels for peak assignments if enabled
        if self.show_peak_labels and individual_surfaces is not None and len(individual_surfaces) > 0:
            label_positions = []

            for i, (surf, peak) in enumerate(zip(individual_surfaces, all_peaks)):
                # Use assignment-based color for consistency across spectra
                peak_assignment = peak.get('assignment', f'Peak {i+1}')
                color = self._get_peak_color(peak_assignment, colors)

                pos_f2 = _to_scalar(peak['pos_f2'])
                pos_f1 = _to_scalar(peak['pos_f1'])

                # Find Z-height at peak maximum
                f2_idx = np.argmin(np.abs(f2_ppm - pos_f2))
                f1_idx = np.argmin(np.abs(f1_ppm - pos_f1))
                max_z = fitted_surface[f1_idx, f2_idx]

                # Position label 15% above peak maximum
                label_z = max_z * 1.15
                label_x = pos_f2
                label_y = pos_f1

                # Simple collision avoidance
                for prev_x, prev_y, prev_z in label_positions:
                    f2_dist = abs(label_x - prev_x) / 0.05
                    f1_dist = abs(label_y - prev_y) / 0.5

                    if f2_dist < 1.5 and f1_dist < 1.5:
                        label_z *= 1.1
                        label_x += 0.01 * (1 if i % 2 == 0 else -1)

                label_positions.append((label_x, label_y, label_z))

                # Draw arrow from label to peak
                ax_data.plot([label_x, pos_f2], [label_y, pos_f1], [label_z, max_z],
                            color=color, linewidth=2, alpha=0.8, linestyle='-')

                # Add 3D text label
                ax_data.text(label_x, label_y, label_z,
                            peak_assignment,
                            color=color,
                            fontsize=18,
                            fontweight='bold',
                            ha='center',
                            va='bottom',
                            bbox=dict(boxstyle='round,pad=0.3',
                                     facecolor='white',
                                     edgecolor=color,
                                     alpha=0.85))

        # Overlay residuals if in overlay mode
        if self.residual_mode == 'overlay' and self.show_residuals:
            max_intensity = np.max(experimental)
            residual_baseline = 0 - (0.2 * max_intensity)

            z_max_resid = np.max(np.abs(residuals))
            residual_range = 0.1 * max_intensity

            if z_max_resid > 0:
                normalized_residuals = residuals / z_max_resid
            else:
                normalized_residuals = residuals

            residual_heights = residual_baseline + (normalized_residuals + 1) * 0.5 * residual_range

            # Use custom colormap with LunaNMR colors
            colors_custom = [
                (0.357, 0.620, 0.898),  # Bright blue
                (1.0, 1.0, 1.0),        # White
                (0.910, 0.333, 0.306)   # Bright red
            ]
            cmap_rdbu = LinearSegmentedColormap.from_list('lunaNMR_residuals', colors_custom)

            norm = Normalize(vmin=-z_max_resid, vmax=z_max_resid)
            colors_mapped = cmap_rdbu(norm(residuals))

            residual_surf = ax_data.plot_surface(
                F2_ppm, F1_ppm, residual_heights,
                facecolors=colors_mapped,
                alpha=0.6,
                edgecolor='none',
                antialiased=True,
                shade=False
            )

        # Configure axes
        ax_data.set_xlabel('¹H Chemical Shift (ppm)', fontsize=9)
        ax_data.set_ylabel('¹⁵N/¹³C Chemical Shift (ppm)', fontsize=9)
        ax_data.set_zlabel('Intensity', fontsize=9)
        ax_data.invert_xaxis()
        ax_data.invert_yaxis()
        ax_data.view_init(elev=self.preserved_elev, azim=self.preserved_azim, roll=self.preserved_roll)

        # Limit tick labels
        ax_data.xaxis.set_major_locator(MaxNLocator(nbins=2))
        ax_data.yaxis.set_major_locator(MaxNLocator(nbins=2))

        # Z-axis height reduced by 30%
        ax_data.set_box_aspect([1, 1, 1.05])

        # Add legend (limit to 6 peaks)
        if len(all_peaks) <= 6 and (self.show_experimental or self.show_fitted or self.show_individual_peaks):
            ax_data.legend(loc='upper right', fontsize=7, framealpha=0.9)

        # Calculate and store auto Z-axis scale
        z_max_data = np.max(experimental)
        z_min_data = np.min(experimental)
        self.auto_z_min = z_min_data
        self.auto_z_max = z_max_data

        # Apply z-axis limits: fixed (for Peak Mode comparison) or auto-scaled
        if self.use_fixed_z_limits and self.fixed_z_min is not None:
            ax_data.set_zlim(self.fixed_z_min, self.fixed_z_max)
            ax_data.set_autoscalez_on(False)
        else:
            self.apply_intensity_scale(self.intensity_scale_factor)

        # Apply x/y axis limits: fixed (for Peak Mode comparison) or auto-scaled
        if self.use_fixed_xy_limits and self.fixed_x_min is not None:
            ax_data.set_xlim(self.fixed_x_min, self.fixed_x_max)
            ax_data.set_ylim(self.fixed_y_min, self.fixed_y_max)

        # Panel 2: Residuals in separate panel (if enabled)
        if self.residual_mode == 'separate' and self.show_residuals and self.ax_resid is not None:
            z_max_resid = np.max(np.abs(residuals))

            # Create custom colormap
            colors_custom = [
                (0.357, 0.620, 0.898),  # Bright blue
                (1.0, 1.0, 1.0),        # White
                (0.910, 0.333, 0.306)   # Bright red
            ]
            cmap_lunaNMR = LinearSegmentedColormap.from_list('lunaNMR_residuals', colors_custom)

            surf_resid = self.ax_resid.plot_surface(F2_ppm, F1_ppm, residuals,
                                                    cmap=cmap_lunaNMR, alpha=0.8,
                                                    vmin=-z_max_resid, vmax=z_max_resid,
                                                    antialiased=True)

            # Add colorbar
            cbar = self.figure.colorbar(surf_resid, ax=self.ax_resid, shrink=0.6, aspect=10)
            cbar.set_label('Residuals', fontsize=8)

            # Configure axes
            self.ax_resid.set_xlabel('¹H Chemical Shift (ppm)', fontsize=9)
            self.ax_resid.set_ylabel('¹⁵N/¹³C Chemical Shift (ppm)', fontsize=9)
            self.ax_resid.set_zlabel('Residuals', fontsize=9)
            self.ax_resid.invert_xaxis()
            self.ax_resid.invert_yaxis()
            self.ax_resid.view_init(elev=self.preserved_elev, azim=self.preserved_azim, roll=self.preserved_roll)

            self.ax_resid.xaxis.set_major_locator(MaxNLocator(nbins=2))
            self.ax_resid.yaxis.set_major_locator(MaxNLocator(nbins=2))

            # Match z-limits to data panel (fixed or auto)
            if self.use_fixed_z_limits and self.fixed_z_min is not None:
                self.ax_resid.set_zlim(self.fixed_z_min, self.fixed_z_max)
            else:
                self.ax_resid.set_zlim(z_min_data, z_max_data)

            # Match x/y limits to data panel if fixed
            if self.use_fixed_xy_limits and self.fixed_x_min is not None:
                self.ax_resid.set_xlim(self.fixed_x_min, self.fixed_x_max)
                self.ax_resid.set_ylim(self.fixed_y_min, self.fixed_y_max)

            self.ax_resid.set_title(f'Residuals (R²={r_squared:.3f})', fontsize=10, fontweight='bold')

        # Panel 3: Cross-sections (if enabled)
        if self.show_cross_sections and self.ax_f1_cross is not None:
            peak_x = voigt_result.get('peak_x', f2_ppm[len(f2_ppm)//2])
            peak_y = voigt_result.get('peak_y', f1_ppm[len(f1_ppm)//2])
            self.update_cross_sections(peak_x, peak_y, voigt_result)

            # Connect click event
            if self.click_event_cid is not None:
                self.canvas.mpl_disconnect(self.click_event_cid)

            def on_click(event):
                if event.inaxes == ax_data or (self.ax_resid and event.inaxes == self.ax_resid):
                    if event.xdata and event.ydata:
                        self.update_cross_sections(event.xdata, event.ydata, voigt_result)
                        # Emit Qt signal for external handlers
                        self.plot_clicked.emit(event.xdata, event.ydata)

            self.click_event_cid = self.canvas.mpl_connect('button_press_event', on_click)
        else:
            if self.click_event_cid is not None:
                self.canvas.mpl_disconnect(self.click_event_cid)
                self.click_event_cid = None

        # Set overall figure title
        self.figure.suptitle(f'3D Voigt Analysis: {assignment}',
                            fontsize=12, fontweight='bold')

        self.figure.tight_layout(rect=[0, 0, 1, 0.96])
        self.canvas.draw()

    def _plot_1d_not_applicable_3d(self, voigt_result):
        """Show message that 3D view is only applicable for 2D simultaneous fits"""
        self.figure.clear()

        ax = self.figure.add_subplot(111)
        ax.axis('off')

        assignment = voigt_result.get('assignment', 'Unknown')
        quality = voigt_result.get('fitting_quality', 'Unknown')

        message_text = f"""
3D Surface View

This view is only available for 2D simultaneous
multi-peak fitting results.

The current peak ({assignment}) was fitted using
1D cross-section fitting, which is best visualized
in the standard 2D Contour view.

Quality: {quality}

Please use the "2D Voigt Analysis" tab for
detailed 1D cross-section visualization.
"""

        ax.text(0.5, 0.5, message_text, transform=ax.transAxes,
               ha='center', va='center', fontsize=11,
               bbox=dict(boxstyle='round,pad=1.0', facecolor='lightyellow',
                        alpha=0.9, edgecolor='gray'))

        self.figure.suptitle('3D Voigt Analysis - Not Applicable',
                            fontsize=12, fontweight='bold')

        self.figure.tight_layout(rect=[0, 0, 1, 0.96])
        self.canvas.draw()

    def _plot_no_data(self):
        """Plot message when no data is available"""
        self.clear_all()

        # Use first axes for message
        ax = self.axes_list[0]
        ax.text(0.5, 0.5, 'No Voigt analysis data\nFit a peak to view results',
               transform=ax.transAxes, ha='center', va='center',
               fontsize=10, bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray"))
        ax.axis('off')

        # Hide other axes
        for ax in self.axes_list[1:]:
            ax.axis('off')

        self.canvas.draw()

    def _plot_no_data_3d(self):
        """Plot message when no data is available for 3D view"""
        self.figure.clear()

        ax = self.figure.add_subplot(111)
        ax.axis('off')

        ax.text(0.5, 0.5, 'No Voigt analysis data\nFit a peak to view 3D results',
               transform=ax.transAxes, ha='center', va='center',
               fontsize=12, bbox=dict(boxstyle="round,pad=1.0",
                                     facecolor="lightgray", alpha=0.9))

        self.canvas.draw()
