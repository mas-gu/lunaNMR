#!/usr/bin/env python3
"""
Visualization Module

This module handles all plotting and visualization functionality for the
NMR Peak Series application, including spectrum plots, Voigt fitting results,
and series analysis visualizations.

Classes:
- SpectrumPlotter: Main spectrum visualization
- VoigtAnalysisPlotter: Voigt fitting results visualization
- SeriesPlotter: Series analysis and statistics plots
- PlotManager: Coordinate plotting across multiple tabs

Author: Guillaume Mas
Date: 2025
"""

# Force TkAgg backend for GUI consistency
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle, Circle
import numpy as np
import pandas as pd
from datetime import datetime

class SpectrumPlotter:
    """Main spectrum visualization with peak overlays"""

    def __init__(self, figure, axis):
        self.fig = figure
        self.ax = axis
        self.contour_plot = None
        self.colorbar = None
        self.colorbar_ax = None  # Store colorbar axes reference
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

    def plot_spectrum(self, integrator, **kwargs):
        """Plot 2D NMR spectrum with contours"""
        if not hasattr(integrator, 'nmr_data') or integrator.nmr_data is None:
            self.ax.text(0.5, 0.5, 'No spectrum data loaded\nPlease load NMR data first',
                        transform=self.ax.transAxes, ha='center', va='center',
                        fontsize='small', bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray"))
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
        if 'contour_min_level' in kwargs and 'contour_levels' in kwargs and 'contour_increment' in kwargs:
            # Use GUI-specified parameters
            min_level = kwargs['contour_min_level'] * np.max(data_abs)
            num_levels = kwargs['contour_levels']
            increment = kwargs['contour_increment']

            # Generate levels using increment method (like original implementation)
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
        self.ax.set_title('lunaNMR v0.9', fontsize='small', fontweight='bold')
        self.ax.invert_xaxis()
        self.ax.invert_yaxis()
        self.ax.grid(True, alpha=0.3)

        # Add colorbar if requested
        if self.settings['show_colorbar'] and self.contour_plot:
            self._create_colorbar_safely()

    def _setup_fixed_layout(self):
        """Setup fixed layout positions for axes and colorbar"""
        # Get current figure size
        fig_width, fig_height = self.fig.get_size_inches()

        # Define fixed positions (left, bottom, width, height) in figure coordinates
        main_axes_pos = [0.12, 0.15, 0.75, 0.75]  # Leave space for colorbar
        colorbar_pos = [0.89, 0.15, 0.03, 0.75]   # Fixed colorbar position

        # Store positions
        self.main_axes_position = main_axes_pos
        self.colorbar_position = colorbar_pos

    def _create_colorbar_safely(self):
        """Safely create or update colorbar with fixed positioning"""
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
        """Overlay peak positions and assignments using a standardized data format."""
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

        if not peak_list:
            # No detected peaks, but we might still have reference peaks to show
            pass
        else:
            # Plot detected/fitted peaks (existing code for red circles)
            if show_detected:
                # Show ALL peaks in fitted_peaks list as "detected" (red circles)
                # This includes both matched peaks and reference-retained peaks
                detected_x = [p['ppm_x'] for p in peak_list]
                detected_y = [p['ppm_y'] for p in peak_list]
                if detected_x:
                    self.ax.scatter(detected_x, detected_y, c='red', marker='o', s=60,
                                  alpha=0.8, edgecolors='white', linewidth=1, zorder=5,
                                  label=f'Detected ({len(detected_x)})')
                    for peak in peak_list:
                        # Use flexible field name lookup for assignment
                        assignment_label = (peak.get('assignment') or
                                          peak.get('Assignment') or
                                          str(peak.get('peak_number', '')))
                        annotation = self.ax.annotate(
                            assignment_label,
                            (peak['ppm_x'], peak['ppm_y']),
                            xytext=(5, 5), textcoords='offset points',
                            fontsize='small', color='red', fontweight='bold',
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8)
                        )
                        self.peak_annotations.append(annotation)

        # ENHANCED: Plot reference peaks (blue crosses) - ROBUST VERSION
        if show_assigned and hasattr(integrator, 'peak_list') and integrator.peak_list is not None:
            try:
                # Handle multiple possible column name formats
                peak_df = integrator.peak_list

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

        # Plot PS2D ellipses if requested (debug tool)
        show_ellipses = kwargs.get('show_ellipses', False)
        nucleus_type = kwargs.get('nucleus_type', '15N')
        if show_ellipses:
            self._plot_ps2d_ellipses(integrator, nucleus_type)

        # Add legend
        handles, labels = self.ax.get_legend_handles_labels()
        if handles:
            self.ax.legend(loc='upper right', fontsize='small')

    def _plot_ps2d_ellipses(self, integrator, nucleus_type='15N'):
        """
        Plot PS2D elliptical integration windows for all peaks (debug tool).

        This visualization shows the elliptical windows that PS2D uses for:
        1. Data selection (radF1_selector × radF2_selector)
        2. Fitting region (radF1 × radF2)

        Args:
            integrator: EnhancedVoigtIntegrator with peak_list
            nucleus_type: '15N' or '13C' for config lookup
        """
        from matplotlib.patches import Ellipse as EllipsePatch

        # Get PS2D configuration
        try:
            from lunaNMR.core.ps2d_config import get_ps2d_config
            # DO NOT call set_ps2d_config() here - it would override user's manual changes!
            # Just read the current config which may have been modified by the user
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

        # Try 1: Reference peak list (DataFrame)
        if hasattr(integrator, 'peak_list') and integrator.peak_list is not None:
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
                        peaks_to_plot.append({
                            'x': float(row[x_col]),
                            'y': float(row[y_col]),
                            'label': row.get('Assignment', f'Peak_{idx+1}')
                        })
            except Exception as e:
                print(f"⚠️  Could not extract peaks from peak_list: {e}")

        # Try 2: Detected/fitted peaks (list of dicts)
        if not peaks_to_plot and hasattr(integrator, 'fitted_peaks') and integrator.fitted_peaks:
            try:
                for peak in integrator.fitted_peaks:
                    peaks_to_plot.append({
                        'x': peak.get('ppm_x', peak.get('x_ppm', 0)),
                        'y': peak.get('ppm_y', peak.get('y_ppm', 0)),
                        'label': peak.get('assignment', peak.get('Assignment', 'Unknown'))
                    })
            except Exception as e:
                print(f"⚠️  Could not extract peaks from fitted_peaks: {e}")

        if not peaks_to_plot:
            print("⚠️  No peaks available for ellipse visualization")
            return

        try:
            # Plot ellipses for each peak
            for idx, peak in enumerate(peaks_to_plot):
                peak_x = peak['x']
                peak_y = peak['y']

                # Data selector ellipse (outer, dashed)
                selector_ellipse = EllipsePatch(
                    (peak_x, peak_y),
                    width=2 * radF2_selector,  # diameter in 1H dimension
                    height=2 * radF1_selector,  # diameter in 15N/13C dimension
                    fill=False,
                    edgecolor='cyan',
                    linestyle='--',
                    linewidth=1.5,
                    alpha=0.6,
                    zorder=6,
                    label='Data Selector' if idx == 0 else ''
                )
                self.ax.add_patch(selector_ellipse)

                # Fitting region ellipse (inner, solid)
                fitting_ellipse = EllipsePatch(
                    (peak_x, peak_y),
                    width=2 * radF2,  # diameter in 1H dimension
                    height=2 * radF1,  # diameter in 15N/13C dimension
                    fill=False,
                    edgecolor='magenta',
                    linestyle='-',
                    linewidth=2,
                    alpha=0.8,
                    zorder=7,
                    label='Fitting Region' if idx == 0 else ''
                )
                self.ax.add_patch(fitting_ellipse)

            print(f"✅ Plotted PS2D ellipses for {len(peaks_to_plot)} peaks")
            print(f"   Fitting region: {radF1:.3f} × {radF2:.3f} ppm")
            print(f"   Data selector: {radF1_selector:.3f} × {radF2_selector:.3f} ppm")

        except Exception as e:
            print(f"❌ Error plotting PS2D ellipses: {e}")
            import traceback
            traceback.print_exc()

    def _plot_fitted_curves(self, voigt_fits):
        """Overlay fitted Voigt curves on spectrum"""
        for fit in voigt_fits[-5:]:  # Show only last 5 fits to avoid clutter
            if 'x_fit' in fit and 'y_fit' in fit:
                x_center = fit['x_fit'].get('center', 0)
                y_center = fit['y_fit'].get('center', 0)
                x_sigma = fit['x_fit'].get('sigma', 0.05)
                y_sigma = fit['y_fit'].get('sigma', 2.0)

                # Draw fitted peak region as ellipse
                ellipse = plt.Circle((x_center, y_center),
                                   radius=max(x_sigma*2, y_sigma*0.5),
                                   fill=False, color='orange', linewidth=2,
                                   alpha=0.8, zorder=6)
                self.ax.add_patch(ellipse)

    def _plot_detection_quality(self, integrator):
        """Overlay detection quality visualization on spectrum"""
        try:
            # Check if integrated detection-fitting results are available
            if hasattr(integrator, 'fitted_peaks') and integrator.fitted_peaks is not None:
                quality_x = []
                quality_y = []
                quality_scores = []
                quality_colors = []

                for peak in integrator.fitted_peaks:
                    if peak.get('detected', False):
                        x_pos = peak['ppm_x']
                        y_pos = peak['ppm_y']

                        # Get detection quality metrics
                        detection_confidence = peak.get('detection_confidence', 0.5)
                        fit_quality = peak.get('fit_r_squared', 0.0)
                        aic_score = peak.get('aic_score', 0.0)

                        # Calculate composite quality score (0-1 scale)
                        composite_quality = (detection_confidence * 0.4 +
                                           min(fit_quality, 1.0) * 0.4 +
                                           max(0, 1.0 - aic_score/10.0) * 0.2)

                        quality_x.append(x_pos)
                        quality_y.append(y_pos)
                        quality_scores.append(composite_quality)

                        # Color mapping: red (poor) -> yellow (moderate) -> green (high quality)
                        if composite_quality < 0.3:
                            color = 'red'
                        elif composite_quality < 0.7:
                            color = 'orange'
                        else:
                            color = 'green'
                        quality_colors.append(color)

                # Plot quality indicators as colored rings around peaks
                if quality_x:
                    for i, (x, y, score, color) in enumerate(zip(quality_x, quality_y, quality_scores, quality_colors)):
                        # Ring size based on quality score
                        ring_size = 0.05 + (score * 0.10)  # Size range: 0.05 to 0.15
                        ring_width = 0.02 + (score * 0.03)  # Width range: 0.02 to 0.05

                        # Outer ring (quality indicator)
                        outer_ring = Circle((x, y), radius=ring_size,
                                          fill=False, color=color, linewidth=3,
                                          alpha=0.8, zorder=8)
                        self.ax.add_patch(outer_ring)

                        # Quality score text overlay (small, unobtrusive)
                        score_text = self.ax.annotate(
                            f'{score:.2f}',
                            (x, y), xytext=(10, -25), textcoords='offset points',
                            fontsize='small', color=color, fontweight='bold',
                            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                                    edgecolor=color, alpha=0.9)
                        )
                        self.peak_annotations.append(score_text)

                    # Add quality legend
                    from matplotlib.patches import Patch
                    quality_legend_elements = [
                        Patch(facecolor='green', alpha=0.8, label='High Quality (≥0.7)'),
                        Patch(facecolor='orange', alpha=0.8, label='Moderate Quality (0.3-0.7)'),
                        Patch(facecolor='red', alpha=0.8, label='Low Quality (<0.3)')
                    ]

                    # Add quality legend to upper left (separate from main legend)
                    quality_legend = self.ax.legend(handles=quality_legend_elements,
                                                  loc='upper left', fontsize='small',
                                                  title='Detection Quality', title_fontsize=9,
                                                  framealpha=0.9, fancybox=True)
                    quality_legend.set_zorder(10)

                    print(f"✅ Detection quality visualization added for {len(quality_x)} peaks")

        except Exception as e:
            print(f"⚠️  Could not add detection quality visualization: {e}")

    def highlight_peak(self, peak_x, peak_y, highlight_color='yellow', radius=0.1):
        """Highlight a specific peak position"""
        # Remove previous highlights
        for annotation in self.peak_annotations:
            if hasattr(annotation, '_highlight'):
                annotation.remove()

        # Add highlight circle
        highlight = Circle((peak_x, peak_y), radius=radius,
                         facecolor=highlight_color, alpha=0.3,
                         edgecolor=highlight_color, linewidth=3, zorder=7)
        highlight._highlight = True
        self.ax.add_patch(highlight)

        # Add crosshair
        self.add_crosshair(peak_x, peak_y)

    def add_crosshair(self, x, y):
        """Add crosshair lines at specified position"""
        # Remove previous crosshairs
        for line in self.crosshair_lines:
            line.remove()
        self.crosshair_lines = []

        # Add new crosshair
        xlims = self.ax.get_xlim()
        ylims = self.ax.get_ylim()

        hline = self.ax.axhline(y=y, color='yellow', linestyle='--', alpha=0.7, zorder=6)
        vline = self.ax.axvline(x=x, color='yellow', linestyle='--', alpha=0.7, zorder=6)

        self.crosshair_lines = [hline, vline]

    def set_zoom(self, x_center, y_center, x_range, y_range):
        """Set zoom level around specified center"""
        self.ax.set_xlim(x_center + x_range/2, x_center - x_range/2)
        self.ax.set_ylim(y_center + y_range/2, y_center - y_range/2)

    def reset_zoom(self):
        """Reset zoom to show full spectrum"""
        self.ax.autoscale()

class VoigtAnalysisPlotter:
    """Voigt fitting analysis visualization"""

    def __init__(self, figure, axes):
        self.fig = figure
        self.axes = axes  # Should be 2x2 grid: [[ax_x, ax_y], [ax_2d, ax_residuals]]

    def plot_voigt_analysis(self, voigt_result):
        """Plot comprehensive Voigt analysis results"""
        if not voigt_result:
            self._plot_no_data()
            return

        # Extract data
        assignment = voigt_result.get('assignment', 'Unknown')
        quality = voigt_result.get('fitting_quality', 'Unknown')
        method = voigt_result.get('method', '')

        # Check if this is 2D simultaneous fitting
        if '2d_simultaneous' in method and 'region_2d' in voigt_result:
            # Use 2D contour visualization (clears figure and creates 2×2 grid internally)
            self._plot_2d_contour(voigt_result, assignment, quality)
        else:
            # Use traditional 1D cross-section visualization
            # CRITICAL: Clear figure and recreate 2×1 axes grid (in case we're switching from 2D mode)
            self.fig.clear()

            # Recreate original 2×1 axes layout
            self.axes = [
                [self.fig.add_subplot(2, 1, 1)],  # Top: X-dimension
                [self.fig.add_subplot(2, 1, 2)]   # Bottom: Y-dimension
            ]

            x_fit = voigt_result.get('x_fit', {})
            y_fit = voigt_result.get('y_fit', {})

            # Plot X-dimension fit (top)
            self._plot_1d_fit(self.axes[0][0], x_fit, '¹H Chemical Shift (ppm)',
                             f'X-Dimension Fit - {assignment}')

            # Plot Y-dimension fit (bottom)
            self._plot_1d_fit(self.axes[1][0], y_fit, '¹⁵N/¹³C Chemical Shift (ppm)',
                             f'Y-Dimension Fit - {assignment}')

            # Set overall title
            self.fig.suptitle(f'Voigt Profile Analysis: {assignment} (Quality: {quality})',
                             fontsize='small', fontweight='bold')

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
        # Use local R-squared if available, otherwise fall back to global
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
        # ENHANCEMENT: Show window size information in title if available
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

    def _plot_1d_residuals(self, ax, fit_data, xlabel, title):
        """Plot 1D fitting residuals"""
        if not fit_data or 'ppm_scale' not in fit_data:
            ax.text(0.5, 0.5, 'No residual data available',
                   transform=ax.transAxes, ha='center', va='center')
            ax.set_title(title)
            return

        ppm_scale = fit_data['ppm_scale']
        cross_section = fit_data['cross_section']
        fitted_curve = fit_data.get('fitted_curve', None)
        # Use local R-squared if available, otherwise fall back to global
        r_squared = fit_data.get('r_squared_local', fit_data.get('r_squared', 0))

        if fitted_curve is not None:
            # Calculate residuals
            residuals = cross_section - fitted_curve

            # Plot residuals
            ax.plot(ppm_scale, residuals, 'g-', linewidth=1.5, label='Residuals')

            # Add zero line for reference
            ax.axhline(y=0, color='k', linestyle='-', alpha=0.5, linewidth=1)

            # Calculate RMS residual for display
            rms_residual = np.sqrt(np.mean(residuals**2))
            ax.text(0.05, 0.95, f'RMS: {rms_residual:.1f}',
                   transform=ax.transAxes, bbox=dict(boxstyle="round,pad=0.3",
                   facecolor='white', alpha=0.8))
        else:
            ax.text(0.5, 0.5, 'No fitted curve available for residuals',
                   transform=ax.transAxes, ha='center', va='center')

        ax.set_xlabel(xlabel)
        ax.set_ylabel('Residual Intensity')
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)

        if '¹H' in xlabel:
            ax.invert_xaxis()

    def _plot_2d_contour(self, voigt_result, assignment, quality):
        """
        Plot 2D contour visualization for simultaneous multi-peak fitting

        Layout (dynamically creates 2x2 grid):
        [0,0] Experimental 2D Data    [0,1] Fitted 2D Surface
        [1,0] 2D Residuals            [1,1] Parameters Table
        """
        import numpy as np

        # Extract 2D data
        region_2d = voigt_result.get('region_2d')
        fitted_surface = voigt_result.get('fitted_2d_surface')
        all_peaks = voigt_result.get('all_peaks', [])
        r_squared = voigt_result.get('avg_r_squared', 0)

        if region_2d is None or fitted_surface is None:
            self._plot_no_data()
            return

        f1_ppm = region_2d['f1_ppm']
        f2_ppm = region_2d['f2_ppm']
        experimental = region_2d['intensity']

        # Calculate residuals
        residuals = experimental - fitted_surface

        # Determine contour levels (adaptive based on data range)
        # Focus on 20-100% of peak height to reveal saddle points between overlapping peaks
        # Use COMBINED range of experimental and fitted to avoid truncation
        vmin = min(np.min(experimental), np.min(fitted_surface))
        vmax = max(np.max(experimental), np.max(fitted_surface))
        levels = np.linspace(vmin + (vmax - vmin) * 0.2, vmax, 20)

        # CRITICAL: Clear figure and create 2x2 grid for 2D visualization
        # The default layout is 2x1, but we need 2x2 for 2D simultaneous fits
        self.fig.clear()
        gs = self.fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
        ax_exp = self.fig.add_subplot(gs[0, 0])
        ax_fit = self.fig.add_subplot(gs[0, 1])
        ax_res = self.fig.add_subplot(gs[1, 0])
        ax_table = self.fig.add_subplot(gs[1, 1])

        # Plot 1: Experimental Data (top left) - BLUE filled contours
        contour_exp = ax_exp.contourf(f2_ppm, f1_ppm, experimental, levels=levels,
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
            ax_exp.plot(peak['pos_f2'], peak['pos_f1'], 'r+', markersize=12, markeredgewidth=2)
            # Use assignment if available, otherwise show peak number
            peak_label = peak.get('assignment', str(i+1))
            ax_exp.text(peak['pos_f2'], peak['pos_f1'], f"  {peak_label}",
                       fontsize=8, color='red', fontweight='bold', va='center')

        # Plot 2: Individual Fitted Peaks (top right) - Multi-color visualization
        individual_surfaces = voigt_result.get('individual_surfaces', None)

        # Define color palette for different peaks (works for n peaks)
        colors = ['red', 'orange', 'purple', 'brown', 'pink', 'olive', 'cyan', 'magenta']

        if individual_surfaces is not None and len(individual_surfaces) > 0:
            # Plot each individual peak with a different color
            for i, (peak_surface, peak) in enumerate(zip(individual_surfaces, all_peaks)):
                color = colors[i % len(colors)]  # Cycle through colors if more than 8 peaks

                # Add baseline to individual surface for proper visualization
                peak_with_baseline = peak_surface + region_2d['intensity'].min()

                # Calculate adaptive levels for this peak (40-100% of its max)
                peak_min = np.min(peak_with_baseline)
                peak_max = np.max(peak_with_baseline)
                peak_levels = np.linspace(peak_min + (peak_max - peak_min) * 0.4, peak_max, 15)

                # Plot contour lines for this peak
                ax_fit.contour(f2_ppm, f1_ppm, peak_with_baseline, levels=peak_levels,
                             colors=color, linewidths=1.2, alpha=0.7)

            ax_fit.set_title(f'Individual Fitted Peaks (R²={r_squared:.3f})', fontsize=9, fontweight='bold')
        else:
            # Fallback: plot summed surface (original behavior)
            contour_fit = ax_fit.contourf(f2_ppm, f1_ppm, fitted_surface, levels=levels,
                                          cmap='Reds', alpha=0.8)
            ax_fit.contour(f2_ppm, f1_ppm, fitted_surface, levels=levels,
                           colors='red', linewidths=0.5, alpha=0.5)
            ax_fit.set_title(f'Fitted (R²={r_squared:.3f}) - SAME levels', fontsize=9, fontweight='bold')

        ax_fit.set_xlabel('¹H Chemical Shift (ppm)', fontsize=9)
        ax_fit.set_ylabel('¹⁵N/¹³C Chemical Shift (ppm)', fontsize=9)
        ax_fit.invert_xaxis()
        ax_fit.invert_yaxis()
        ax_fit.grid(True, alpha=0.3)

        # Mark peak positions with matching colors
        for i, peak in enumerate(all_peaks):
            if individual_surfaces is not None and len(individual_surfaces) > 0:
                color = colors[i % len(colors)]
            else:
                color = 'blue'
            ax_fit.plot(peak['pos_f2'], peak['pos_f1'], '+', color=color,
                       markersize=12, markeredgewidth=2)
            # Use assignment if available, otherwise show peak number
            peak_label = peak.get('assignment', str(i+1))
            ax_fit.text(peak['pos_f2'], peak['pos_f1'], f"  {peak_label}",
                       fontsize=8, color=color, fontweight='bold', va='center')

        # Plot 3: Residuals (bottom left) - Option D: Color-coded quality regions
        # Calculate relative residuals as percentage
        relative_residuals = 100 * residuals / np.max(experimental)

        # Create discrete quality regions with strict thresholds
        # good (<10%), moderate (10-20%), poor (>20%)
        quality_map = np.abs(relative_residuals)

        # Define custom colormap: green -> yellow -> red
        from matplotlib.colors import BoundaryNorm, ListedColormap
        colors_quality = ['#2ecc71', '#f1c40f', '#e74c3c']  # green, yellow, red
        boundaries = [0, 10, 20, 100]  # percentage boundaries (strict for well-resolved fits)
        cmap_quality = ListedColormap(colors_quality)
        norm = BoundaryNorm(boundaries, cmap_quality.N)

        # Plot filled contours with discrete colors
        contourf_res = ax_res.contourf(f2_ppm, f1_ppm, quality_map,
                                       levels=boundaries, cmap=cmap_quality, norm=norm, alpha=0.7)

        # Add peak positions for reference
        for i, peak in enumerate(all_peaks):
            ax_res.plot(peak['pos_f2'], peak['pos_f1'], 'k+', markersize=8, markeredgewidth=1.5)

        ax_res.set_xlabel('¹H Chemical Shift (ppm)', fontsize=9)
        ax_res.set_ylabel('¹⁵N/¹³C Chemical Shift (ppm)', fontsize=9)
        ax_res.set_title('Fit Quality (Residual %)', fontsize=9, fontweight='bold')
        ax_res.invert_xaxis()
        ax_res.invert_yaxis()
        ax_res.grid(True, alpha=0.3)

        # Add legend for quality regions (strict thresholds)
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#2ecc71', alpha=0.7, label='Good (<10%)'),
            Patch(facecolor='#f1c40f', alpha=0.7, label='Moderate (10-20%)'),
            Patch(facecolor='#e74c3c', alpha=0.7, label='Poor (>20%)')
        ]
        ax_res.legend(handles=legend_elements, loc='upper right', fontsize=7, framealpha=0.9)

        # Plot 4: Parameters Table (bottom right)
        ax_table.axis('off')

        # Create cleaner table with vertical layout
        table_text = "Fitted Parameters\n" + "="*30 + "\n\n"

        for i, peak in enumerate(all_peaks):
            # Calculate total linewidth (Gaussian + Lorentzian)
            lw_f2 = peak['lw_gau_f2'] + peak['lw_lor_f2']
            lw_f1 = peak['lw_gau_f1'] + peak['lw_lor_f1']

            # Use matching color for peak number
            colors_list = ['red', 'orange', 'purple', 'brown', 'pink', 'olive', 'cyan', 'magenta']
            color = colors_list[i % len(colors_list)]

            # Use assignment if available, otherwise fall back to "Peak N"
            peak_label = peak.get('assignment', f'Peak {i+1}')

            # Get volume (use intensity as fallback for backward compatibility)
            volume = peak.get('volume', peak.get('intensity', 0.0))

            # Line 1: Peak assignment
            table_text += f"{peak_label}\n"
            # Line 2: LW 1H and LW 15N
            table_text += f"LW ¹H: {lw_f2:.4f}  LW ¹⁵N: {lw_f1:.3f}\n"
            # Line 3: Volume
            table_text += f"Volume: {volume:.2e}\n"
            if i < len(all_peaks) - 1:
                table_text += "\n"

        # Display as monospace text
        ax_table.text(0.1, 0.95, table_text, transform=ax_table.transAxes,
                     fontsize=8, fontfamily='monospace', verticalalignment='top',
                     bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))

        # Add overall statistics
        stats_text = f"Quality: {quality}\nR² = {r_squared:.4f}\n{len(all_peaks)} peaks (simultaneous fit)"
        ax_table.text(0.5, 0.1, stats_text, ha='center', va='center',
                     fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Set overall title
        self.fig.suptitle(f'2D Simultaneous Voigt Fitting: {assignment}',
                         fontsize=10, fontweight='bold')

    def _plot_2d_summary(self, ax, voigt_result):
        """Plot 2D fit summary information"""
        ax.axis('off')

        # Extract key information
        assignment = voigt_result.get('assignment', 'Unknown')
        peak_pos = voigt_result.get('peak_position', (0, 0))
        quality = voigt_result.get('fitting_quality', 'Unknown')
        # Use local average R-squared if available, otherwise fall back to global
        avg_r2 = voigt_result.get('avg_r_squared_local', voigt_result.get('avg_r_squared', 0))
        timestamp = voigt_result.get('timestamp', 'Unknown')

        x_fit = voigt_result.get('x_fit', {})
        y_fit = voigt_result.get('y_fit', {})

        # Create summary text
        summary_text = f"""
2D Voigt Fit Results

Peak: {assignment}
Position: ({peak_pos[0]:.3f}, {peak_pos[1]:.1f}) ppm
Quality: {quality}
Average R²: {avg_r2:.3f}

X-Dimension Parameters:
  Center: {x_fit.get('center', 0):.3f} ppm
  σ (Gaussian): {x_fit.get('sigma', 0):.3f}
  γ (Lorentzian): {x_fit.get('gamma', 0):.3f}
  Amplitude: {x_fit.get('amplitude', 0):.0f}

Y-Dimension Parameters:
  Center: {y_fit.get('center', 0):.1f} ppm
  σ (Gaussian): {y_fit.get('sigma', 0):.1f}
  γ (Lorentzian): {y_fit.get('gamma', 0):.1f}
  Amplitude: {y_fit.get('amplitude', 0):.0f}

Fit Timestamp: {str(timestamp)[:19] if timestamp != 'Unknown' else 'Unknown'}
        """

        # Choose background color based on quality
        if quality == 'Excellent':
            bgcolor = 'lightgreen'
        elif quality == 'Good':
            bgcolor = 'lightblue'
        elif quality == 'Fair':
            bgcolor = 'lightyellow'
        else:
            bgcolor = 'lightcoral'

        ax.text(0.05, 0.95, summary_text.strip(), transform=ax.transAxes,
               fontsize='small', fontfamily='monospace', verticalalignment='top',
               bbox=dict(boxstyle="round,pad=0.5", facecolor=bgcolor, alpha=0.8))

    def _plot_parameters(self, ax, voigt_result):
        """Plot fitting parameters visualization"""
        ax.axis('off')

        x_fit = voigt_result.get('x_fit', {})
        y_fit = voigt_result.get('y_fit', {})

        if not x_fit or not y_fit:
            ax.text(0.5, 0.5, 'No parameter data available',
                   transform=ax.transAxes, ha='center', va='center')
            return

        # Create parameter comparison
        params = ['sigma', 'gamma', 'amplitude', 'r_squared']
        x_values = [x_fit.get(p, 0) for p in params]
        y_values = [y_fit.get(p, 0) for p in params]

        # Normalize values for comparison
        max_vals = [max(x_values[i], y_values[i]) for i in range(len(params))]
        x_norm = [x_values[i] / max_vals[i] if max_vals[i] > 0 else 0 for i in range(len(params))]
        y_norm = [y_values[i] / max_vals[i] if max_vals[i] > 0 else 0 for i in range(len(params))]

        # Create bar chart
        x_pos = np.arange(len(params))
        width = 0.35

        bars1 = ax.bar(x_pos - width/2, x_norm, width, label='X-dim', color='skyblue', alpha=0.8)
        bars2 = ax.bar(x_pos + width/2, y_norm, width, label='Y-dim', color='lightcoral', alpha=0.8)

        ax.set_xlabel('Parameters')
        ax.set_ylabel('Normalized Values')
        ax.set_title('Parameter Comparison (Normalized)')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(['σ', 'γ', 'Amplitude', 'R²'])
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Add value labels on bars
        for i, (bar1, bar2) in enumerate(zip(bars1, bars2)):
            ax.text(bar1.get_x() + bar1.get_width()/2, bar1.get_height() + 0.01,
                   f'{x_values[i]:.3f}', ha='center', va='bottom', fontsize='small')
            ax.text(bar2.get_x() + bar2.get_width()/2, bar2.get_height() + 0.01,
                   f'{y_values[i]:.3f}', ha='center', va='bottom', fontsize='small')

    def _plot_no_data(self):
        """Plot message when no data is available"""
        # Clear figure and recreate 2×1 axes layout (in case we're coming from 2D mode)
        self.fig.clear()

        # Recreate original 2×1 axes layout
        self.axes = [
            [self.fig.add_subplot(2, 1, 1)],  # Top
            [self.fig.add_subplot(2, 1, 2)]   # Bottom
        ]

        for ax_row in self.axes:
            for ax in ax_row:
                ax.text(0.5, 0.5, 'No Voigt analysis data\nFit a peak to view results',
                    transform=ax.transAxes, ha='center', va='center',
                    fontsize='small', bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray"))
                ax.axis('off')

class SeriesPlotter:
    """Series analysis and statistics visualization"""

    def __init__(self, figure, axis):
        self.fig = figure
        self.ax = axis

    def plot_series_overview(self, batch_results):
        """Plot series processing overview"""
        if not batch_results or not batch_results.results:
            self._plot_no_series_data()
            return

        self.ax.clear()

        # Extract data
        spectra = list(batch_results.results.keys())
        detection_rates = []
        statuses = []

        for spectrum_name, result in batch_results.results.items():
            detection_rates.append(result.get('detection_rate', 0.0))
            statuses.append(result.get('status', 'failed'))

        # Create color map based on status
        colors = ['green' if status == 'success' else 'red' for status in statuses]

        # Create bar plot
        x_pos = np.arange(len(spectra))
        bars = self.ax.bar(x_pos, detection_rates, color=colors, alpha=0.7,
                          edgecolor='black', linewidth=0.5)

        # Customize plot
        self.ax.set_xlabel('Spectrum Index')
        self.ax.set_ylabel('Detection Rate (%)')
        self.ax.set_title(f'Series Integration Results ({len(spectra)} spectra)')
        self.ax.set_ylim(0, 100)

        # Add reference lines
        self.ax.axhline(y=50, color='orange', linestyle='--', alpha=0.7, label='50% threshold')
        self.ax.axhline(y=80, color='green', linestyle='--', alpha=0.7, label='80% threshold')

        # Customize x-axis
        if len(spectra) <= 20:
            self.ax.set_xticks(x_pos)
            self.ax.set_xticklabels([s.replace('.ft', '') for s in spectra],
                                   rotation=45, ha='right')
        else:
            self.ax.set_xlabel('Spectrum Index (too many to label individually)')

        # Add statistics text
        summary = batch_results.get_summary()
        stats_text = (f"Success Rate: {summary['success_rate']:.1f}%\n"
                      f"Avg Detection: {np.mean(detection_rates):.1f}%\n"
                      f"Processing Time: {summary['duration']}")

        self.ax.text(0.02, 0.98, stats_text, transform=self.ax.transAxes,
                    verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3",
                    facecolor="white", alpha=0.8))

        self.ax.grid(True, alpha=0.3)
        self.ax.legend()

        # Add hover information (simplified)
        def on_hover(event):
            if event.inaxes == self.ax:
                for i, bar in enumerate(bars):
                    if bar.contains(event)[0]:
                        spectrum = spectra[i]
                        rate = detection_rates[i]
                        status = statuses[i]
                        self.ax.set_title(f'{spectrum}: {rate:.1f}% detection ({status})')
                        self.fig.canvas.draw_idle()
                        break

        self.fig.canvas.mpl_connect('motion_notify_event', on_hover)

    def plot_detection_statistics(self, batch_results):
        """Plot detailed detection statistics"""
        if not batch_results or not batch_results.statistics:
            return

        stats = batch_results.statistics

        # Create subplots within the main axis (simplified)
        self.ax.clear()

        if 'detection_rate' in stats:
            det_stats = stats['detection_rate']

            # Create histogram of detection rates
            detection_rates = []
            for result in batch_results.results.values():
                if result['status'] == 'success':
                    detection_rates.append(result.get('detection_rate', 0.0))

            if detection_rates:
                self.ax.hist(detection_rates, bins=20, alpha=0.7, color='skyblue',
                           edgecolor='black', linewidth=0.5)
                self.ax.axvline(det_stats['mean'], color='red', linestyle='--',
                              linewidth=2, label=f"Mean: {det_stats['mean']:.1f}%")
                self.ax.axvline(det_stats['median'], color='green', linestyle='--',
                              linewidth=2, label=f"Median: {det_stats['median']:.1f}%")

                self.ax.set_xlabel('Detection Rate (%)')
                self.ax.set_ylabel('Number of Spectra')
                self.ax.set_title('Detection Rate Distribution')
                self.ax.legend()
                self.ax.grid(True, alpha=0.3)

    def _plot_no_series_data(self):
        """Plot message when no series data is available"""
        self.ax.clear()
        self.ax.text(0.5, 0.5, 'No series results available\nRun series integration to view analysis',
                    transform=self.ax.transAxes, ha='center', va='center',
                    fontsize='small', bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray"))
        self.ax.set_title('Series Analysis - No Data')
        self.ax.axis('off')

class PlotManager:
    """Coordinate plotting across multiple tabs and figures"""

    def __init__(self):
        self.spectrum_plotter = None
        self.voigt_plotter = None
        self.series_plotter = None

        # Plot settings
        self.global_settings = {
            'dpi': 100,
            'style': 'default',
            'font_size': 10
        }

    def register_plotters(self, spectrum_plotter=None, voigt_plotter=None, series_plotter=None):
        """Register plotter instances"""
        if spectrum_plotter:
            self.spectrum_plotter = spectrum_plotter
        if voigt_plotter:
            self.voigt_plotter = voigt_plotter
        if series_plotter:
            self.series_plotter = series_plotter

    def update_all_plots(self, integrator=None, voigt_result=None, batch_results=None, **kwargs):
        """Update all plots with current data"""
        if self.spectrum_plotter and integrator:
            self.spectrum_plotter.plot_spectrum(integrator, **kwargs)
            self.spectrum_plotter.plot_peaks(integrator, **kwargs)
            self.spectrum_plotter.fig.canvas.draw()

        if self.voigt_plotter and voigt_result:
            self.voigt_plotter.plot_voigt_analysis(voigt_result)
            self.voigt_plotter.fig.canvas.draw()

        if self.series_plotter and batch_results:
            self.series_plotter.plot_series_overview(batch_results)
            self.series_plotter.fig.canvas.draw()

    def apply_global_settings(self):
        """Apply global plot settings"""
        plt.rcParams['figure.dpi'] = self.global_settings['dpi']
        plt.rcParams['font.size'] = self.global_settings['font_size']

        if self.global_settings['style'] in plt.style.available:
            plt.style.use(self.global_settings['style'])

    def save_all_plots(self, output_folder):
        """Save all current plots to files"""
        saved_files = []

        try:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

            if self.spectrum_plotter:
                spectrum_file = f"{output_folder}/spectrum_plot_{timestamp}.png"
                self.spectrum_plotter.fig.savefig(spectrum_file, dpi=300, bbox_inches='tight')
                saved_files.append(spectrum_file)

            if self.voigt_plotter:
                voigt_file = f"{output_folder}/voigt_analysis_{timestamp}.png"
                self.voigt_plotter.fig.savefig(voigt_file, dpi=300, bbox_inches='tight')
                saved_files.append(voigt_file)

            if self.series_plotter:
                series_file = f"{output_folder}/series_overview_{timestamp}.png"
                self.series_plotter.fig.savefig(series_file, dpi=300, bbox_inches='tight')
                saved_files.append(series_file)

        except Exception as e:
            print(f"Warning: Failed to save some plots: {e}")

        return saved_files
