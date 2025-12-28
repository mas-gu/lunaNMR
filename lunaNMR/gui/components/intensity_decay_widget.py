# ABOUTME: Widget for displaying intensity vs delay plots for NMR relaxation series
# ABOUTME: Shows volume decay across spectra with highlighted current spectrum point

"""
IntensityDecayWidget - Relaxation Decay Visualization

Displays peak volume vs delay time for T1/T2 relaxation analysis.
Shows data points across all spectra in a series with the current
spectrum highlighted for easy comparison with the 3D Voigt surface.

Features:
- Extracts delay values from spectrum names (numeric or filename patterns)
- Plots volume vs delay with highlighted current point
- Updates highlight as user navigates through spectra
- Minimal design for side-by-side display with 3D Voigt plot
"""

import re
from typing import List, Tuple, Optional, Dict, Any

from PySide6.QtWidgets import QWidget, QVBoxLayout, QSizePolicy, QLabel
from PySide6.QtCore import Signal
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from lunaNMR.gui.styles.design_system import (
    PANEL_BG_COLOR,
    PRIMARY_TEXT,
    SPACING_SM,
)


def extract_delay_from_spectrum_name(spectrum_name: str) -> Optional[float]:
    """
    Extract delay value in milliseconds from a spectrum name.

    Handles multiple formats:
    - Pure numeric: "50" -> 50.0
    - Numeric with suffix: "50_2" -> 50.0 (duplicate handling)
    - Filename with ms: "T1_50ms.ft" -> 50.0
    - Filename with s: "T1_1s.ft" -> 1000.0

    Args:
        spectrum_name: Name of the spectrum (may be filename or delay value)

    Returns:
        Delay in milliseconds, or None if unparseable
    """
    if not spectrum_name:
        return None

    # Remove file extension if present
    name = spectrum_name.rsplit('.', 1)[0] if '.' in spectrum_name else spectrum_name

    # Try numeric with duplicate suffix FIRST (e.g., "50_2")
    # Must check before pure numeric to avoid "50_2" -> 502.0
    match = re.match(r'^(\d+(?:\.\d+)?)_\d+$', name)
    if match:
        return float(match.group(1))

    # Try pure numeric (e.g., "50", "100")
    try:
        return float(name)
    except ValueError:
        pass

    # Try _XXms pattern (milliseconds)
    match = re.search(r'_(\d+(?:\.\d+)?)ms(?:$|\.)', spectrum_name, re.IGNORECASE)
    if match:
        return float(match.group(1))

    # Try _Xs pattern (seconds -> milliseconds)
    match = re.search(r'_(\d+(?:\.\d+)?)s(?:$|\.)', spectrum_name, re.IGNORECASE)
    if match:
        return float(match.group(1)) * 1000.0

    return None


def collect_decay_data(
    spectra: List[Dict[str, Any]],
    assignment: str
) -> Tuple[List[float], List[float], List[int]]:
    """
    Collect decay data (delay vs volume) for a peak across all spectra.

    Args:
        spectra: List of spectrum dicts with 'name' and 'fitted_peaks'
        assignment: Peak assignment to find

    Returns:
        Tuple of (delays_ms, volumes, spectrum_indices)
        Only includes spectra where peak exists and delay is parseable
    """
    delays = []
    volumes = []
    indices = []

    for idx, spec in enumerate(spectra):
        # Parse delay from spectrum name
        delay_ms = extract_delay_from_spectrum_name(spec.get('name', ''))
        if delay_ms is None:
            continue

        # Find peak in this spectrum
        fitted_peaks = spec.get('fitted_peaks', [])
        peak_volume = None

        for peak in fitted_peaks:
            peak_assignment = peak.get('assignment', peak.get('Assignment', ''))
            if peak_assignment == assignment:
                # Try multiple keys for volume
                peak_volume = (
                    peak.get('volume') or
                    peak.get('Volume') or
                    peak.get('intensity') or
                    peak.get('Intensity') or
                    0.0
                )
                break

        if peak_volume is not None:
            delays.append(delay_ms)
            volumes.append(peak_volume)
            indices.append(idx)

    return delays, volumes, indices


class IntensityDecayWidget(QWidget):
    """
    Widget displaying intensity (volume) vs delay plot.

    Shows decay curve for a selected peak across all spectra in a series,
    with the currently selected spectrum highlighted.

    Signals:
        point_clicked(int): Emitted when user clicks a data point, with spectrum index
    """

    point_clicked = Signal(int)  # Emits spectrum index when point is clicked

    def __init__(self, parent=None):
        """Initialize the intensity decay widget."""
        super().__init__(parent)

        self._delays = []
        self._volumes = []
        self._indices = []
        self._current_index = None
        self._assignment = ""

        # Exponential fit data from DynamiXs T1/T2 fitting
        self._fit_curve_time = None
        self._fit_curve_intensity = None
        self._t_value = None
        self._t_error = None
        self._time_units = 'ms'
        self._experiment_type = ''

        self._setup_ui()

    def _setup_ui(self):
        """Set up the widget UI."""
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(SPACING_SM)

        # Title label
        self.title_label = QLabel("Intensity vs Delay")
        self.title_label.setStyleSheet(f"""
            QLabel {{
                font-weight: bold;
                font-size: 13px;
                color: {PRIMARY_TEXT};
                padding: 4px;
            }}
        """)
        layout.addWidget(self.title_label)

        # Create matplotlib figure
        self.figure = Figure(figsize=(4, 3), dpi=100)
        self.figure.set_facecolor(PANEL_BG_COLOR)
        self.axes = self.figure.add_subplot(111)
        self._style_axes()

        # Create canvas
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.mpl_connect('button_press_event', self._on_click)

        layout.addWidget(self.canvas, stretch=1)

        self.setLayout(layout)

    def _style_axes(self):
        """Apply styling to axes."""
        self.axes.set_facecolor(PANEL_BG_COLOR)
        self.axes.tick_params(colors=PRIMARY_TEXT, labelsize=9)
        for spine in self.axes.spines.values():
            spine.set_color(PRIMARY_TEXT)
        self.axes.xaxis.label.set_color(PRIMARY_TEXT)
        self.axes.yaxis.label.set_color(PRIMARY_TEXT)

    def set_data(
        self,
        spectra: List[Dict[str, Any]],
        assignment: str,
        current_spectrum_index: Optional[int] = None
    ):
        """
        Set decay data from spectra list.

        Args:
            spectra: List of spectrum dicts with 'name' and 'fitted_peaks'
            assignment: Peak assignment to display
            current_spectrum_index: Index of currently selected spectrum (highlighted)
        """
        self._assignment = assignment
        self._delays, self._volumes, self._indices = collect_decay_data(spectra, assignment)
        self._current_index = current_spectrum_index

        self._update_plot()

    def set_highlight(self, spectrum_index: int):
        """
        Update which spectrum is highlighted.

        Args:
            spectrum_index: Index of spectrum to highlight
        """
        self._current_index = spectrum_index
        self._update_plot()

    def set_exponential_fit(self, fit_data: dict):
        """Set exponential fit data from DynamiXs T1/T2 fitting.

        Args:
            fit_data: Dict containing:
                - fit_curve: {'time': [...], 'intensity': [...]} for smooth curve
                - t_value: T1 or T2 value
                - t_error: uncertainty
                - time_units: time units string (e.g., 'ms')
                - experiment_type: 'T1' or 'T2'
        """
        if fit_data is None:
            self._fit_curve_time = None
            self._fit_curve_intensity = None
            self._t_value = None
            self._t_error = None
            self._time_units = 'ms'
            self._experiment_type = ''
        else:
            fit_curve = fit_data.get('fit_curve', {})
            self._fit_curve_time = fit_curve.get('time', [])
            self._fit_curve_intensity = fit_curve.get('intensity', [])
            self._t_value = fit_data.get('t_value')
            self._t_error = fit_data.get('t_error')
            self._time_units = fit_data.get('time_units', 'ms')
            self._experiment_type = fit_data.get('experiment_type', '')

        self._update_plot()

    def clear_exponential_fit(self):
        """Clear exponential fit data."""
        self.set_exponential_fit(None)

    def _update_plot(self):
        """Redraw the plot with current data."""
        self.axes.clear()
        self._style_axes()

        if not self._delays:
            self.axes.text(
                0.5, 0.5, "No decay data available",
                ha='center', va='center',
                transform=self.axes.transAxes,
                fontsize=10, color='gray'
            )
            self.canvas.draw()
            return

        # Plot all points
        self.axes.scatter(
            self._delays, self._volumes,
            s=50, c='#1f77b4', alpha=0.7, zorder=2,
            label='Data points'
        )

        # Connect with line (only if no exponential fit, otherwise fit curve replaces this)
        if not self._fit_curve_time:
            sorted_data = sorted(zip(self._delays, self._volumes, self._indices))
            sorted_delays = [d[0] for d in sorted_data]
            sorted_volumes = [d[1] for d in sorted_data]

            self.axes.plot(
                sorted_delays, sorted_volumes,
                '-', color='#1f77b4', alpha=0.3, linewidth=1, zorder=1
            )

        # Plot exponential fit curve from DynamiXs (if available)
        if self._fit_curve_time and self._fit_curve_intensity:
            self.axes.plot(
                self._fit_curve_time, self._fit_curve_intensity,
                '-', color='#e63946', linewidth=2, zorder=1.5,
                label='Exponential fit'
            )

        # Highlight current spectrum
        if self._current_index is not None and self._current_index in self._indices:
            idx_pos = self._indices.index(self._current_index)
            self.axes.scatter(
                [self._delays[idx_pos]], [self._volumes[idx_pos]],
                s=150, c='red', marker='o', zorder=3,
                edgecolors='darkred', linewidths=2,
                label='Current'
            )

        # Add T1/T2 annotation if available
        if self._t_value is not None:
            exp_type = self._experiment_type or 'T'
            if self._t_error is not None:
                annotation = f"{exp_type} = {self._t_value:.2f} ± {self._t_error:.2f} {self._time_units}"
            else:
                annotation = f"{exp_type} = {self._t_value:.2f} {self._time_units}"
            self.axes.text(
                0.95, 0.95, annotation,
                transform=self.axes.transAxes,
                fontsize=9, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='#cccccc')
            )

        # Labels (use time_units if available)
        time_label = f"Delay ({self._time_units})" if self._time_units else "Delay (ms)"
        self.axes.set_xlabel(time_label, fontsize=10)
        self.axes.set_ylabel("Volume", fontsize=10)

        # Title with assignment
        if self._assignment:
            self.title_label.setText(f"Intensity vs Delay: {self._assignment}")

        # Format y-axis with scientific notation if needed
        self.axes.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))

        self.figure.tight_layout()
        self.canvas.draw()

    def _on_click(self, event):
        """Handle click events on the plot."""
        if event.inaxes != self.axes:
            return

        if not self._delays:
            return

        # Find closest point
        click_x = event.xdata
        click_y = event.ydata

        if click_x is None or click_y is None:
            return

        # Normalize distances (x and y have different scales)
        x_range = max(self._delays) - min(self._delays) if len(self._delays) > 1 else 1
        y_range = max(self._volumes) - min(self._volumes) if len(self._volumes) > 1 else 1

        min_dist = float('inf')
        closest_idx = None

        for i, (d, v, spec_idx) in enumerate(zip(self._delays, self._volumes, self._indices)):
            dx = (d - click_x) / x_range if x_range > 0 else 0
            dy = (v - click_y) / y_range if y_range > 0 else 0
            dist = dx * dx + dy * dy

            if dist < min_dist:
                min_dist = dist
                closest_idx = spec_idx

        # Only emit if click is reasonably close to a point
        if closest_idx is not None and min_dist < 0.05:  # 5% threshold
            self.point_clicked.emit(closest_idx)

    def clear(self):
        """Clear the plot."""
        self._delays = []
        self._volumes = []
        self._indices = []
        self._current_index = None
        self._assignment = ""
        # Clear exponential fit data
        self._fit_curve_time = None
        self._fit_curve_intensity = None
        self._t_value = None
        self._t_error = None
        self._time_units = 'ms'
        self._experiment_type = ''
        self._update_plot()
