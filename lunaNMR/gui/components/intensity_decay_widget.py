# ABOUTME: Widget for displaying intensity vs delay plots for NMR relaxation series
# ABOUTME: Shows volume/height/detected intensity decay with selectable value type

"""
IntensityDecayWidget - Relaxation Decay Visualization

Displays peak values vs delay time for T1/T2 relaxation analysis.
Shows data points across all spectra in a series with the current
spectrum highlighted for easy comparison with the 3D Voigt surface.

Features:
- Selectable value type: Volume, Height, or Detected Intensity
- Extracts delay values from spectrum names (numeric or filename patterns)
- Plots selected value vs delay with highlighted current point
- Updates highlight as user navigates through spectra
- Minimal design for side-by-side display with 3D Voigt plot
"""

import re
from typing import List, Tuple, Optional, Dict, Any

from PySide6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QSizePolicy, QLabel, QComboBox
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


def extract_titration_from_spectrum_name(spectrum_name: str) -> Optional[float]:
    """
    Extract a dimensionless titration point from a spectrum name.

    Reads the trailing _<value> suffix using the filesystem-safe 'o'-for-'.'
    convention (e.g. "sample_1o0.ft" -> 1.0, "sample_1.0.ft" -> 1.0,
    "titr_0o5" -> 0.5). A decimal separator ('o' or '.') is required, so plain
    integer indices (e.g. "_001", "_2") are not mistaken for titration points.

    Returns:
        Titration point, or None if the name has no decimal-separated suffix.
    """
    if not spectrum_name:
        return None

    name = spectrum_name.rsplit('.', 1)[0] if '.' in spectrum_name else spectrum_name

    match = re.search(r'_(\d+[o.]\d+)$', name)
    if match:
        return float(match.group(1).replace('o', '.'))

    return None


def extract_index_from_spectrum_name(spectrum_name: str) -> Optional[int]:
    """
    Extract spectrum index from trailing 3+ digits in filename.

    Fallback for when no delay pattern is found. Extracts the last
    group of 3 or more digits before the file extension.

    Examples:
        "03_2D_series_ref_001.ft" -> 1
        "experiment_042.ucsf" -> 42
        "scan_0123.ft" -> 123

    Args:
        spectrum_name: Filename to parse

    Returns:
        Index as integer, or None if no 3+ digit pattern found
    """
    if not spectrum_name:
        return None

    # Match 3+ digits immediately before extension (e.g., _001.ft, _042.ucsf)
    match = re.search(r'_(\d{3,})(?:\.[^.]+)?$', spectrum_name)
    if match:
        return int(match.group(1))

    return None


def _extract_peak_value(peak: Dict[str, Any], value_type: str) -> Optional[float]:
    """Extract value from peak dict based on value_type."""
    if value_type == 'volume':
        return (
            peak.get('volume') or
            peak.get('Volume') or
            peak.get('intensity') or
            peak.get('Intensity') or
            0.0
        )
    elif value_type == 'height':
        return (
            peak.get('height') or
            peak.get('Height') or
            0.0
        )
    elif value_type == 'detected_intensity':
        return (
            peak.get('detected_intensity') or
            peak.get('Detected_Intensity') or
            0.0
        )
    else:
        return peak.get('volume', 0.0)


def collect_decay_data(
    spectra: List[Dict[str, Any]],
    assignment: str,
    value_type: str = 'volume'
) -> Tuple[List[float], List[float], List[int], str]:
    """
    Collect decay data (x-axis vs value) for a peak across all spectra.

    X-axis priority:
    1. Titration mode: o-decimal suffix (_1o0, _0o5) -> dimensionless point
    2. Delay mode: Extract delay from filename patterns (_50ms, _1s, etc.)
    3. Index mode (fallback): Extract 3+ digit index from filename (_001.ft)

    Args:
        spectra: List of spectrum dicts with 'name' and 'fitted_peaks'
        assignment: Peak assignment to find
        value_type: Type of value to extract ('volume', 'height', 'detected_intensity')

    Returns:
        Tuple of (x_values, values, spectrum_indices, mode)
        - x_values: titration points, delays in ms, or spectrum indices
        - values: peak values
        - spectrum_indices: original spectrum indices
        - mode: 'titration', 'delay', or 'index'
    """
    # First pass: try titration extraction (unambiguous o-decimal suffix)
    titration_x = []
    values = []
    indices = []

    for idx, spec in enumerate(spectra):
        point = extract_titration_from_spectrum_name(spec.get('name', ''))
        if point is None:
            continue

        # Find peak in this spectrum
        fitted_peaks = spec.get('fitted_peaks', [])
        for peak in fitted_peaks:
            peak_assignment = peak.get('assignment', peak.get('Assignment', ''))
            if peak_assignment == assignment:
                peak_value = _extract_peak_value(peak, value_type)
                if peak_value is not None:
                    titration_x.append(point)
                    values.append(peak_value)
                    indices.append(idx)
                break

    if titration_x:
        return titration_x, values, indices, 'titration'

    # Second pass: try delay extraction
    delays = []
    values = []
    indices = []

    for idx, spec in enumerate(spectra):
        delay_ms = extract_delay_from_spectrum_name(spec.get('name', ''))
        if delay_ms is None:
            continue

        # Find peak in this spectrum
        fitted_peaks = spec.get('fitted_peaks', [])
        for peak in fitted_peaks:
            peak_assignment = peak.get('assignment', peak.get('Assignment', ''))
            if peak_assignment == assignment:
                peak_value = _extract_peak_value(peak, value_type)
                if peak_value is not None:
                    delays.append(delay_ms)
                    values.append(peak_value)
                    indices.append(idx)
                break

    # If any delays found, use delay mode
    if delays:
        return delays, values, indices, 'delay'

    # Fallback: try index extraction from 3+ trailing digits
    x_values = []
    values = []
    indices = []

    for idx, spec in enumerate(spectra):
        spec_index = extract_index_from_spectrum_name(spec.get('name', ''))
        if spec_index is None:
            continue

        # Find peak in this spectrum
        fitted_peaks = spec.get('fitted_peaks', [])
        for peak in fitted_peaks:
            peak_assignment = peak.get('assignment', peak.get('Assignment', ''))
            if peak_assignment == assignment:
                peak_value = _extract_peak_value(peak, value_type)
                if peak_value is not None:
                    x_values.append(float(spec_index))
                    values.append(peak_value)
                    indices.append(idx)
                break

    return x_values, values, indices, 'index'


class IntensityDecayWidget(QWidget):
    """
    Widget displaying intensity (volume) vs delay plot.

    Shows decay curve for a selected peak across all spectra in a series,
    with the currently selected spectrum highlighted.

    Signals:
        point_clicked(int): Emitted when user clicks a data point, with spectrum index
    """

    point_clicked = Signal(int)  # Emits spectrum index when point is clicked

    # Value type options for the combo box
    VALUE_TYPES = [
        ('Volume', 'volume'),
        ('Height', 'height'),
        ('Detected Intensity', 'detected_intensity'),
    ]

    def __init__(self, parent=None):
        """Initialize the intensity decay widget."""
        super().__init__(parent)

        self._x_values = []  # Delay (ms) or spectrum index depending on mode
        self._values = []
        self._indices = []
        self._current_index = None
        self._assignment = ""
        self._value_type = 'volume'  # Default to volume
        self._x_mode = 'delay'  # 'delay' or 'index' - determines X-axis label
        self._spectra = []  # Store spectra reference for re-collection on mode change

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

        # Header row with title and value type selector
        header_layout = QHBoxLayout()
        header_layout.setContentsMargins(0, 0, 0, 0)
        header_layout.setSpacing(SPACING_SM)

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
        header_layout.addWidget(self.title_label)

        header_layout.addStretch()

        # Value type selector combo box
        self.value_type_combo = QComboBox()
        self.value_type_combo.setToolTip("Select which value to plot vs delay")
        for display_name, _ in self.VALUE_TYPES:
            self.value_type_combo.addItem(display_name)
        self.value_type_combo.setCurrentIndex(0)  # Default to Volume
        self.value_type_combo.currentIndexChanged.connect(self._on_value_type_changed)
        self.value_type_combo.setStyleSheet(f"""
            QComboBox {{
                font-size: 11px;
                padding: 2px 6px;
                min-width: 100px;
            }}
        """)
        header_layout.addWidget(self.value_type_combo)

        layout.addLayout(header_layout)

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
        self._spectra = spectra  # Store for re-collection on mode change
        self._assignment = assignment
        self._x_values, self._values, self._indices, self._x_mode = collect_decay_data(
            spectra, assignment, self._value_type
        )
        self._current_index = current_spectrum_index

        self._update_plot()

    def _on_value_type_changed(self, index: int):
        """Handle value type combo box selection change."""
        if 0 <= index < len(self.VALUE_TYPES):
            _, self._value_type = self.VALUE_TYPES[index]
            # Re-collect data with new value type
            if self._spectra and self._assignment:
                self._x_values, self._values, self._indices, self._x_mode = collect_decay_data(
                    self._spectra, self._assignment, self._value_type
                )
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

    def _get_y_axis_label(self) -> str:
        """Get Y-axis label based on current value type."""
        labels = {
            'volume': 'Volume',
            'height': 'Height',
            'detected_intensity': 'Detected Intensity',
        }
        return labels.get(self._value_type, 'Volume')

    def _get_x_axis_label(self) -> str:
        """Get X-axis label based on current mode."""
        if self._x_mode == 'titration':
            return "Titration point"
        if self._x_mode == 'delay':
            return f"Delay ({self._time_units})" if self._time_units else "Delay (ms)"
        else:
            return "Spectrum #"

    def _update_plot(self):
        """Redraw the plot with current data."""
        self.axes.clear()
        self._style_axes()

        if not self._x_values:
            self.axes.text(
                0.5, 0.5, "No data available",
                ha='center', va='center',
                transform=self.axes.transAxes,
                fontsize=10, color='gray'
            )
            self.canvas.draw()
            return

        # Plot all points
        self.axes.scatter(
            self._x_values, self._values,
            s=50, c='#1f77b4', alpha=0.7, zorder=2,
            label='Data points'
        )

        # Connect with line (only if no exponential fit, otherwise fit curve replaces this)
        if not self._fit_curve_time:
            sorted_data = sorted(zip(self._x_values, self._values, self._indices))
            sorted_x = [d[0] for d in sorted_data]
            sorted_values = [d[1] for d in sorted_data]

            self.axes.plot(
                sorted_x, sorted_values,
                '-', color='#1f77b4', alpha=0.3, linewidth=1, zorder=1
            )

        # Plot exponential fit curve from DynamiXs (if available, only in delay mode)
        if self._x_mode == 'delay' and self._fit_curve_time and self._fit_curve_intensity:
            self.axes.plot(
                self._fit_curve_time, self._fit_curve_intensity,
                '-', color='#e63946', linewidth=2, zorder=1.5,
                label='Exponential fit'
            )

        # Highlight current spectrum
        if self._current_index is not None and self._current_index in self._indices:
            idx_pos = self._indices.index(self._current_index)
            self.axes.scatter(
                [self._x_values[idx_pos]], [self._values[idx_pos]],
                s=150, c='red', marker='o', zorder=3,
                edgecolors='darkred', linewidths=2,
                label='Current'
            )

        # Add T1/T2 annotation if available (only in delay mode)
        if self._x_mode == 'delay' and self._t_value is not None:
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

        # Labels
        self.axes.set_xlabel(self._get_x_axis_label(), fontsize=10)
        self.axes.set_ylabel(self._get_y_axis_label(), fontsize=10)

        # Title with assignment and value type
        if self._assignment:
            x_label = {'titration': 'Titration', 'delay': 'Delay'}.get(self._x_mode, 'Spectrum #')
            self.title_label.setText(f"{self._get_y_axis_label()} vs {x_label}: {self._assignment}")

        # Format y-axis with scientific notation if needed
        self.axes.ticklabel_format(axis='y', style='scientific', scilimits=(0, 3))

        # Force y-axis to start at 0
        self.axes.set_ylim(bottom=0)

        self.figure.tight_layout()
        self.canvas.draw()

    def _on_click(self, event):
        """Handle click events on the plot."""
        if event.inaxes != self.axes:
            return

        if not self._x_values:
            return

        # Find closest point
        click_x = event.xdata
        click_y = event.ydata

        if click_x is None or click_y is None:
            return

        # Normalize distances (x and y have different scales)
        x_range = max(self._x_values) - min(self._x_values) if len(self._x_values) > 1 else 1
        y_range = max(self._values) - min(self._values) if len(self._values) > 1 else 1

        min_dist = float('inf')
        closest_idx = None

        for i, (x, v, spec_idx) in enumerate(zip(self._x_values, self._values, self._indices)):
            dx = (x - click_x) / x_range if x_range > 0 else 0
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
        self._x_values = []
        self._values = []
        self._indices = []
        self._current_index = None
        self._assignment = ""
        self._x_mode = 'delay'
        self._spectra = []
        # Clear exponential fit data
        self._fit_curve_time = None
        self._fit_curve_intensity = None
        self._t_value = None
        self._t_error = None
        self._time_units = 'ms'
        self._experiment_type = ''
        self._update_plot()
