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
                detected_x = [p['ppm_x'] for p in peak_list if p.get('detected')]
                detected_y = [p['ppm_y'] for p in peak_list if p.get('detected')]
                if detected_x:
                    self.ax.scatter(detected_x, detected_y, c='red', marker='o', s=60,
                                  alpha=0.8, edgecolors='white', linewidth=1, zorder=5,
                                  label=f'Detected ({len(detected_x)})')
                    for peak in peak_list:
                        if peak.get('detected'):
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

        # Add legend
        handles, labels = self.ax.get_legend_handles_labels()
        if handles:
            self.ax.legend(loc='upper right', fontsize='small')

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

        # Clear all axes
        for ax_row in self.axes:
            for ax in ax_row:
                ax.clear()

        # Extract data
        x_fit = voigt_result.get('x_fit', {})
        y_fit = voigt_result.get('y_fit', {})
        assignment = voigt_result.get('assignment', 'Unknown')
        quality = voigt_result.get('fitting_quality', 'Unknown')

        # Plot X-dimension fit (top left)
        self._plot_1d_fit(self.axes[0][0], x_fit, '¹H Chemical Shift (ppm)',
                         f'X-Dimension Fit - {assignment}')

        # Plot Y-dimension fit (top right)
        self._plot_1d_fit(self.axes[1][0], y_fit, '¹⁵N/¹³C Chemical Shift (ppm)',
                         f'Y-Dimension Fit - {assignment}')

        # Plot X-dimension residuals (bottom left)
        #self._plot_1d_residuals(self.axes[1][0], x_fit, '¹H Chemical Shift (ppm)',
        #                       f'X-Dimension Residuals - {assignment}')

        # Plot Y-dimension residuals (bottom right)
        #self._plot_1d_residuals(self.axes[1][1], y_fit, '¹⁵N/¹³C Chemical Shift (ppm)',
        #                       f'Y-Dimension Residuals - {assignment}')

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
        ax.set_title(f'{title} (R² = {r_squared:.3f})')
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
        for ax_row in self.axes:
            for ax in ax_row:
                ax.clear()
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
