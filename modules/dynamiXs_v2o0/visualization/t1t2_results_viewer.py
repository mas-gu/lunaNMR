# ABOUTME: Interactive viewer for T1/T2 fitting results from DynamiXs.
# ABOUTME: Displays relaxation time values vs residue with error bars and outlier handling.

"""
T1/T2 Results Viewer Module - Interactive visualization of T1/T2 fitting results
PySide6 Qt Implementation

Displays T1 or T2 relaxation times vs residue sequence with error bars.
Includes outlier detection and Y-axis control.
"""

import pandas as pd
import numpy as np
import os
import re
from pathlib import Path

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFrame,
    QLabel, QPushButton, QCheckBox, QDoubleSpinBox, QFileDialog,
    QScrollArea, QSizePolicy
)
from PySide6.QtCore import Qt

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

# Import GUI components
import sys
sys.path.append(str(Path(__file__).parent.parent))
from constants import (
    BG_COLOR, PANEL_BG_COLOR, FRAME_BG_COLOR,
    PRIMARY_TEXT, SECONDARY_TEXT,
    SPACING_XS, SPACING_SM, SPACING_MD,
    FONT_SIZE_SECTION_LABEL, FONT_SIZE_BODY, FONT_SIZE_SMALL
)
from gui_components import (
    create_primary_button, create_secondary_button, create_label,
    get_font, show_info, show_error, show_warning
)


class T1T2ResultsViewer(QMainWindow):
    """
    Interactive viewer for T1/T2 fitting results.

    Features:
    - Residue-based scatter plot with error bars
    - IQR-based outlier detection with visual distinction
    - Option to exclude outliers from plot
    - User-controlled Y-axis limits
    - Error bar display toggle
    - CSV/TSV file loading
    """

    def __init__(self, parent=None, results_file=None, experiment_type="T1", session_results=None):
        """
        Initialize the T1/T2 Results Viewer.

        Parameters:
        -----------
        parent : QWidget, optional
            Parent widget
        results_file : str, optional
            Path to the results file (TSV format from fitting)
        experiment_type : str
            Type of experiment: "T1" or "T2"
        session_results : dict, optional
            Session results from T1T2FittingPage containing multiple experiments.
            Format: {
                'field1_t1': {'results_file': '...', 'frequency_mhz': 600, 'experiment_type': 'T1'},
                'field1_t2': {'results_file': '...', 'frequency_mhz': 600, 'experiment_type': 'T2'},
                ...
            }
        """
        super().__init__(parent)

        self.setWindowTitle("T1/T2 Results Viewer")
        self.setMinimumSize(1200, 700)
        self.setStyleSheet(f"background-color: {BG_COLOR}; color: {PRIMARY_TEXT};")

        # Store parameters
        self.results_file = results_file
        self.experiment_type = experiment_type
        self.session_results = session_results

        # Multi-experiment data storage {key: DataFrame}
        self.experiment_data = {}

        # Data storage (legacy single file mode)
        self.df = None
        self.show_errors = True
        self.exclude_outliers = False

        # Per-experiment-type Y-axis settings
        self.y_axis_settings = {
            'T1': {'auto': True, 'min': 0.0, 'max': 2000.0},
            'T2': {'auto': True, 'min': 0.0, 'max': 500.0}
        }

        # Outlier detection results (cached)
        self._outlier_mask = None
        self._iqr_bounds = None

        # Build UI
        self._create_ui()

        # Load data
        if session_results:
            self._load_session_results(session_results)
        elif results_file and os.path.exists(results_file):
            self._load_results_file(results_file)
        else:
            self._show_blank_state()

    def _load_session_results(self, session_results):
        """Load data from multiple experiments in session."""
        self.experiment_data.clear()
        loaded_count = 0

        for key, entry in session_results.items():
            if entry is None:
                continue

            results_file = entry.get('results_file')
            if not results_file or not os.path.exists(results_file):
                continue

            try:
                # Load file
                try:
                    df = pd.read_csv(results_file, sep='\t')
                except:
                    df = pd.read_csv(results_file)

                # Store with metadata
                self.experiment_data[key] = {
                    'df': df,
                    'frequency_mhz': entry.get('frequency_mhz', 0),
                    'experiment_type': entry.get('experiment_type', 'T1'),
                    'file_path': results_file
                }
                loaded_count += 1

            except Exception as e:
                print(f"Warning: Could not load {results_file}: {e}")

        if loaded_count > 0:
            self._update_statistics()
            self._update_plot()
            self.file_label.setText(f"Loaded {loaded_count} experiment(s)")
        else:
            self._show_blank_state()

    def _create_ui(self):
        """Create the user interface."""
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        main_layout = QHBoxLayout(central_widget)
        main_layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        main_layout.setSpacing(SPACING_SM)

        # Left panel - Plot (75% width)
        left_panel = QFrame()
        left_panel.setStyleSheet(f"""
            QFrame {{
                background-color: {PANEL_BG_COLOR};
                border-radius: 8px;
            }}
        """)
        self._create_plot_panel(left_panel)
        main_layout.addWidget(left_panel, 3)

        # Right panel - Controls (25% width)
        right_panel = QFrame()
        right_panel.setStyleSheet(f"""
            QFrame {{
                background-color: {PANEL_BG_COLOR};
                border-radius: 8px;
            }}
        """)
        right_panel.setFixedWidth(300)
        self._create_control_panel(right_panel)
        main_layout.addWidget(right_panel, 0)

    def _create_plot_panel(self, parent):
        """Create the plot display panel."""
        layout = QVBoxLayout(parent)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        layout.setSpacing(SPACING_SM)

        # Header
        title_label = create_label("Relaxation Time Results")
        title_label.setFont(get_font(18, bold=True))
        layout.addWidget(title_label)

        # Plot area
        plot_frame = QFrame()
        plot_frame.setStyleSheet(f"background-color: {FRAME_BG_COLOR}; border-radius: 8px;")
        plot_layout = QVBoxLayout(plot_frame)
        plot_layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)

        # Create matplotlib figure
        self.figure = Figure(figsize=(10, 6), facecolor=PANEL_BG_COLOR)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        plot_layout.addWidget(self.canvas)

        # Toolbar
        self.toolbar = NavigationToolbar2QT(self.canvas, plot_frame)
        plot_layout.addWidget(self.toolbar)

        layout.addWidget(plot_frame, 1)

    def _create_control_panel(self, parent):
        """Create the control panel."""
        layout = QVBoxLayout(parent)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        layout.setSpacing(SPACING_MD)

        # Header
        title_label = create_label("Controls")
        title_label.setFont(get_font(16, bold=True))
        layout.addWidget(title_label)

        # File selection section
        file_section = self._create_section_frame("Results File:")
        file_layout = file_section.layout()

        self.file_label = QLabel("No file loaded")
        self.file_label.setFont(get_font(FONT_SIZE_SMALL))
        self.file_label.setStyleSheet(f"""
            QLabel {{
                background-color: {BG_COLOR};
                color: {SECONDARY_TEXT};
                padding: 6px;
                border-radius: 4px;
            }}
        """)
        self.file_label.setWordWrap(True)
        file_layout.addWidget(self.file_label)

        browse_btn = create_primary_button("Load Results File", clicked=self._browse_file)
        browse_btn.setFixedWidth(250)
        file_layout.addWidget(browse_btn)

        layout.addWidget(file_section)

        # Display options section
        options_section = self._create_section_frame("Display Options:")
        options_layout = options_section.layout()

        self.error_checkbox = QCheckBox("Show Error Bars")
        self.error_checkbox.setFont(get_font(FONT_SIZE_BODY))
        self.error_checkbox.setChecked(True)
        self.error_checkbox.stateChanged.connect(self._on_error_toggle)
        options_layout.addWidget(self.error_checkbox)

        self.exclude_outliers_checkbox = QCheckBox("Exclude Outliers (IQR)")
        self.exclude_outliers_checkbox.setFont(get_font(FONT_SIZE_BODY))
        self.exclude_outliers_checkbox.setChecked(False)
        self.exclude_outliers_checkbox.stateChanged.connect(self._on_exclude_outliers_toggle)
        options_layout.addWidget(self.exclude_outliers_checkbox)

        layout.addWidget(options_section)

        # T1 Y-Axis limits section
        t1_yaxis_section = self._create_section_frame("T1 Y-Axis Limits:")
        t1_yaxis_layout = t1_yaxis_section.layout()

        self.t1_auto_y_checkbox = QCheckBox("Auto")
        self.t1_auto_y_checkbox.setFont(get_font(FONT_SIZE_BODY))
        self.t1_auto_y_checkbox.setChecked(True)
        self.t1_auto_y_checkbox.stateChanged.connect(lambda s: self._on_auto_y_toggle('T1'))
        t1_yaxis_layout.addWidget(self.t1_auto_y_checkbox)

        # T1 Y-min/max row
        t1_minmax_row = QFrame()
        t1_minmax_layout = QHBoxLayout(t1_minmax_row)
        t1_minmax_layout.setContentsMargins(0, 0, 0, 0)
        t1_minmax_layout.setSpacing(SPACING_XS)

        t1_min_label = QLabel("Min:")
        t1_min_label.setFont(get_font(FONT_SIZE_SMALL))
        t1_minmax_layout.addWidget(t1_min_label)

        self.t1_ymin_spin = QDoubleSpinBox()
        self.t1_ymin_spin.setRange(0, 100000)
        self.t1_ymin_spin.setValue(0)
        self.t1_ymin_spin.setSuffix(" ms")
        self.t1_ymin_spin.setEnabled(False)
        self.t1_ymin_spin.valueChanged.connect(lambda v: self._on_yaxis_changed('T1'))
        t1_minmax_layout.addWidget(self.t1_ymin_spin)

        t1_max_label = QLabel("Max:")
        t1_max_label.setFont(get_font(FONT_SIZE_SMALL))
        t1_minmax_layout.addWidget(t1_max_label)

        self.t1_ymax_spin = QDoubleSpinBox()
        self.t1_ymax_spin.setRange(0, 100000)
        self.t1_ymax_spin.setValue(2000)
        self.t1_ymax_spin.setSuffix(" ms")
        self.t1_ymax_spin.setEnabled(False)
        self.t1_ymax_spin.valueChanged.connect(lambda v: self._on_yaxis_changed('T1'))
        t1_minmax_layout.addWidget(self.t1_ymax_spin)

        t1_yaxis_layout.addWidget(t1_minmax_row)
        layout.addWidget(t1_yaxis_section)

        # T2 Y-Axis limits section
        t2_yaxis_section = self._create_section_frame("T2 Y-Axis Limits:")
        t2_yaxis_layout = t2_yaxis_section.layout()

        self.t2_auto_y_checkbox = QCheckBox("Auto")
        self.t2_auto_y_checkbox.setFont(get_font(FONT_SIZE_BODY))
        self.t2_auto_y_checkbox.setChecked(True)
        self.t2_auto_y_checkbox.stateChanged.connect(lambda s: self._on_auto_y_toggle('T2'))
        t2_yaxis_layout.addWidget(self.t2_auto_y_checkbox)

        # T2 Y-min/max row
        t2_minmax_row = QFrame()
        t2_minmax_layout = QHBoxLayout(t2_minmax_row)
        t2_minmax_layout.setContentsMargins(0, 0, 0, 0)
        t2_minmax_layout.setSpacing(SPACING_XS)

        t2_min_label = QLabel("Min:")
        t2_min_label.setFont(get_font(FONT_SIZE_SMALL))
        t2_minmax_layout.addWidget(t2_min_label)

        self.t2_ymin_spin = QDoubleSpinBox()
        self.t2_ymin_spin.setRange(0, 100000)
        self.t2_ymin_spin.setValue(0)
        self.t2_ymin_spin.setSuffix(" ms")
        self.t2_ymin_spin.setEnabled(False)
        self.t2_ymin_spin.valueChanged.connect(lambda v: self._on_yaxis_changed('T2'))
        t2_minmax_layout.addWidget(self.t2_ymin_spin)

        t2_max_label = QLabel("Max:")
        t2_max_label.setFont(get_font(FONT_SIZE_SMALL))
        t2_minmax_layout.addWidget(t2_max_label)

        self.t2_ymax_spin = QDoubleSpinBox()
        self.t2_ymax_spin.setRange(0, 100000)
        self.t2_ymax_spin.setValue(500)
        self.t2_ymax_spin.setSuffix(" ms")
        self.t2_ymax_spin.setEnabled(False)
        self.t2_ymax_spin.valueChanged.connect(lambda v: self._on_yaxis_changed('T2'))
        t2_minmax_layout.addWidget(self.t2_ymax_spin)

        t2_yaxis_layout.addWidget(t2_minmax_row)
        layout.addWidget(t2_yaxis_section)

        # Statistics section
        stats_section = self._create_section_frame("Statistics:")
        stats_layout = stats_section.layout()

        self.stats_label = QLabel("No data loaded")
        self.stats_label.setFont(get_font(FONT_SIZE_SMALL))
        self.stats_label.setStyleSheet(f"color: {SECONDARY_TEXT};")
        self.stats_label.setWordWrap(True)
        stats_layout.addWidget(self.stats_label)

        layout.addWidget(stats_section)

        # Action buttons
        btn_frame = QFrame()
        btn_layout = QVBoxLayout(btn_frame)
        btn_layout.setContentsMargins(0, 0, 0, 0)
        btn_layout.setSpacing(SPACING_SM)

        refresh_btn = create_secondary_button("Refresh Plot", clicked=self._update_plot)
        refresh_btn.setFixedWidth(250)
        btn_layout.addWidget(refresh_btn)

        export_btn = create_secondary_button("Export Plot", clicked=self._export_plot)
        export_btn.setFixedWidth(250)
        btn_layout.addWidget(export_btn)

        layout.addWidget(btn_frame)

        layout.addStretch()

    def _create_section_frame(self, title):
        """Create a styled section frame with title."""
        frame = QFrame()
        frame.setStyleSheet(f"""
            QFrame {{
                background-color: {FRAME_BG_COLOR};
                border-radius: 8px;
            }}
        """)

        layout = QVBoxLayout(frame)
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_XS)

        title_label = create_label(title)
        title_label.setFont(get_font(FONT_SIZE_SECTION_LABEL, bold=True))
        layout.addWidget(title_label)

        return frame

    def _browse_file(self):
        """Browse for a results file."""
        initial_dir = os.path.dirname(self.results_file) if self.results_file else os.getcwd()
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Results File",
            initial_dir,
            "Text files (*.txt);;CSV files (*.csv);;All files (*.*)"
        )

        if file_path:
            self._load_results_file(file_path)

    def _load_results_file(self, file_path):
        """Load results from a TSV/CSV file."""
        try:
            # Try tab-delimited first (standard output from fitting)
            try:
                self.df = pd.read_csv(file_path, sep='\t')
            except:
                self.df = pd.read_csv(file_path)

            self.results_file = file_path
            self.file_label.setText(os.path.basename(file_path))

            # Detect experiment type from columns
            if 'T1' in self.df.columns:
                self.experiment_type = 'T1'
            elif 'T2' in self.df.columns:
                self.experiment_type = 'T2'

            # Reset outlier cache
            self._outlier_mask = None
            self._iqr_bounds = None

            # Update statistics
            self._update_statistics()

            # Update plot
            self._update_plot()

            show_info(self, "File Loaded",
                     f"Loaded {len(self.df)} residues from:\n{os.path.basename(file_path)}")

        except Exception as e:
            show_error(self, "Load Error", f"Error loading results file:\n{str(e)}")
            self.df = None
            self._show_blank_state()

    def _detect_outliers(self, values):
        """
        Detect outliers using IQR method.

        Parameters:
        -----------
        values : array-like
            Values to check for outliers

        Returns:
        --------
        tuple : (outlier_mask, lower_bound, upper_bound)
            - outlier_mask: boolean array, True for outliers
            - lower_bound: lower IQR bound
            - upper_bound: upper IQR bound
        """
        values = np.array(values)
        valid_values = values[~np.isnan(values)]

        if len(valid_values) < 4:
            # Not enough data for IQR
            return np.zeros(len(values), dtype=bool), np.nanmin(values), np.nanmax(values)

        Q1 = np.percentile(valid_values, 25)
        Q3 = np.percentile(valid_values, 75)
        IQR = Q3 - Q1

        # Use 1.5 * IQR rule
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR

        # Also apply physical constraints (T1/T2 should be positive and reasonable)
        lower_bound = max(lower_bound, 0)  # Can't have negative relaxation times

        outlier_mask = (values < lower_bound) | (values > upper_bound) | np.isnan(values)

        return outlier_mask, lower_bound, upper_bound

    def _update_statistics(self):
        """Update statistics display."""
        # Multi-experiment mode
        if self.experiment_data:
            self._update_statistics_multi()
            return

        # Legacy single-file mode
        if self.df is None:
            self.stats_label.setText("No data loaded")
            return

        # Get valid data (where Success == "Yes" or the column doesn't exist)
        if 'Success' in self.df.columns:
            valid_df = self.df[self.df['Success'] == 'Yes']
        else:
            valid_df = self.df

        exp_type = self.experiment_type

        # Find value column - handle variations
        value_col = None
        for col in valid_df.columns:
            if col == exp_type or col.startswith(f"{exp_type} ") or col.startswith(f"{exp_type}("):
                value_col = col
                break

        if value_col and value_col in valid_df.columns:
            values = valid_df[value_col].dropna()
            if len(values) > 0:
                # Detect outliers
                outlier_mask, lower_bound, upper_bound = self._detect_outliers(values.values)
                n_outliers = np.sum(outlier_mask)

                # Stats excluding outliers
                clean_values = values.values[~outlier_mask]

                if len(clean_values) > 0:
                    stats_text = (
                        f"Experiment: {exp_type}\n"
                        f"Total residues: {len(values)}\n"
                        f"Outliers: {n_outliers}\n"
                        f"─────────────\n"
                        f"Mean: {np.mean(clean_values):.2f} ms\n"
                        f"Std: {np.std(clean_values):.2f} ms\n"
                        f"Min: {np.min(clean_values):.2f} ms\n"
                        f"Max: {np.max(clean_values):.2f} ms\n"
                        f"─────────────\n"
                        f"IQR range: [{lower_bound:.1f}, {upper_bound:.1f}]"
                    )
                else:
                    stats_text = f"Experiment: {exp_type}\nAll {len(values)} values are outliers"

                self.stats_label.setText(stats_text)
                return

        self.stats_label.setText("No valid data")

    def _update_statistics_multi(self):
        """Update statistics for multi-experiment mode."""
        if not self.experiment_data:
            self.stats_label.setText("No data loaded")
            return

        stats_lines = []
        for key, data in sorted(self.experiment_data.items()):
            df = data['df']
            exp_type = data['experiment_type']

            # Filter for successful fits
            if 'Success' in df.columns:
                valid_df = df[df['Success'] == 'Yes']
            else:
                valid_df = df

            # Find value column
            value_col = None
            for col in valid_df.columns:
                if col == exp_type or col.startswith(f"{exp_type} ") or col.startswith(f"{exp_type}("):
                    value_col = col
                    break

            if value_col and len(valid_df) > 0:
                values = valid_df[value_col].dropna().values
                if len(values) > 0:
                    stats_lines.append(f"{key}: n={len(values)}, mean={np.mean(values):.1f} ms")
                else:
                    stats_lines.append(f"{key}: no valid values")
            else:
                stats_lines.append(f"{key}: column not found")

        self.stats_label.setText("\n".join(stats_lines) if stats_lines else "No data")

    def _extract_residue_numbers(self, residue_labels):
        """
        Extract numeric residue numbers from residue labels.

        Handles formats like: "75.GLY" -> 75, "GLY75" -> 75, "75" -> 75
        """
        numeric_residues = []

        for label in residue_labels:
            if isinstance(label, (int, float, np.integer, np.floating)):
                numeric_residues.append(int(label))
                continue

            label_str = str(label)
            numbers = re.findall(r'\d+', label_str)

            if numbers:
                numeric_residues.append(int(numbers[0]))
            else:
                numeric_residues.append(len(numeric_residues) + 1)

        return np.array(numeric_residues)

    def _add_missing_residue_bars(self, ax, residues, values):
        """Add grey bars for missing residues to show gaps in sequence."""
        if len(residues) == 0:
            return

        min_residue = 1
        max_residue = int(np.max(residues))

        present_residues = set(residues)
        all_residues = set(range(min_residue, max_residue + 1))
        missing_residues = sorted(all_residues - present_residues)

        if not missing_residues:
            return

        for res in missing_residues:
            ax.axvspan(res - 0.5, res + 0.5,
                      facecolor='lightgrey', alpha=0.3, zorder=0)

    def _format_xaxis(self, ax, residues):
        """Format x-axis with ticks every 10 residues."""
        if len(residues) == 0:
            return

        min_residue = 1
        max_residue = int(np.max(residues))

        ax.set_xlim(min_residue - 0.5, max_residue + 0.5)

        if max_residue >= 10:
            first_tick = 10
            ticks = list(range(first_tick, max_residue + 1, 10))
            if 1 not in ticks:
                ticks = [1] + ticks
            if max_residue not in ticks and max_residue - ticks[-1] > 2:
                ticks.append(max_residue)
        else:
            ticks = list(range(1, max_residue + 1))

        ax.set_xticks(ticks)
        ax.set_xticklabels([str(t) for t in ticks])

    def _update_plot(self):
        """Update the plot with current data (single or multi-experiment)."""
        # Use multi-experiment mode if data available
        if self.experiment_data:
            self._update_plot_multi()
            return

        # Legacy single-file mode
        if self.df is None:
            self._show_blank_state()
            return

        self.figure.clear()
        ax = self.figure.add_subplot(111)

        # Style the plot like Model-free viewer
        ax.set_facecolor(FRAME_BG_COLOR)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        exp_type = self.experiment_type

        # Filter for successful fits
        if 'Success' in self.df.columns:
            valid_df = self.df[self.df['Success'] == 'Yes'].copy()
        else:
            valid_df = self.df.copy()

        if len(valid_df) == 0:
            ax.text(0.5, 0.5, "No successful fits to display",
                   ha='center', va='center', fontsize=14, color=SECONDARY_TEXT,
                   transform=ax.transAxes)
            self.canvas.draw()
            return

        # Extract residue numbers for x-axis (actual residue numbers, not indices)
        residues = self._extract_residue_numbers(valid_df['Residue'].values)
        valid_df['ResNum'] = residues

        # Sort by residue number
        valid_df = valid_df.sort_values('ResNum')
        residues = valid_df['ResNum'].values

        # Get values
        y = valid_df[exp_type].values

        # Detect outliers
        outlier_mask, lower_bound, upper_bound = self._detect_outliers(y)
        self._outlier_mask = outlier_mask
        self._iqr_bounds = (lower_bound, upper_bound)

        # Add grey bars for missing residues
        self._add_missing_residue_bars(ax, residues, y)

        # Get errors if available
        err_col = f"{exp_type}_err"
        has_errors = err_col in valid_df.columns and self.show_errors
        if has_errors:
            yerr = valid_df[err_col].values
        else:
            yerr = None

        # Separate normal and outlier data
        normal_mask = ~outlier_mask

        if self.exclude_outliers:
            # Only plot normal points
            x_normal = residues[normal_mask]
            y_normal = y[normal_mask]
            if has_errors:
                yerr_normal = yerr[normal_mask]
                ax.errorbar(x_normal, y_normal, yerr=yerr_normal, fmt='o', capsize=3,
                           color='#2196F3', ecolor=SECONDARY_TEXT,
                           markersize=6, label=f'{exp_type} values', zorder=3)
            else:
                ax.plot(x_normal, y_normal, 'o', color='#2196F3', markersize=6,
                       label=f'{exp_type} values', zorder=3)
        else:
            # Plot normal points
            x_normal = residues[normal_mask]
            y_normal = y[normal_mask]
            if has_errors:
                yerr_normal = yerr[normal_mask]
                ax.errorbar(x_normal, y_normal, yerr=yerr_normal, fmt='o', capsize=3,
                           color='#2196F3', ecolor=SECONDARY_TEXT,
                           markersize=6, label=f'{exp_type} values', zorder=3)
            else:
                ax.plot(x_normal, y_normal, 'o', color='#2196F3', markersize=6,
                       label=f'{exp_type} values', zorder=3)

            # Plot outliers with distinct style (gray, smaller, different marker)
            if np.any(outlier_mask):
                x_outlier = residues[outlier_mask]
                y_outlier = y[outlier_mask]
                if has_errors:
                    yerr_outlier = yerr[outlier_mask]
                    ax.errorbar(x_outlier, y_outlier, yerr=yerr_outlier, fmt='x', capsize=2,
                               color='#9E9E9E', ecolor='#BDBDBD',
                               markersize=5, label=f'Outliers ({np.sum(outlier_mask)})', zorder=3)
                else:
                    ax.plot(x_outlier, y_outlier, 'x', color='#9E9E9E', markersize=5,
                           label=f'Outliers ({np.sum(outlier_mask)})', zorder=3)

        # Set Y-axis limits - always start at 0
        # Use per-experiment-type settings
        auto_checkbox = self.t1_auto_y_checkbox if exp_type == 'T1' else self.t2_auto_y_checkbox
        ymin_spin = self.t1_ymin_spin if exp_type == 'T1' else self.t2_ymin_spin
        ymax_spin = self.t1_ymax_spin if exp_type == 'T1' else self.t2_ymax_spin

        if auto_checkbox.isChecked():
            # Auto mode: always start at 0, use IQR upper bound with padding
            y_range = upper_bound - lower_bound
            padding = y_range * 0.1 if y_range > 0 else 50
            ax.set_ylim(0, upper_bound + padding)

            # Update spinboxes to show current auto values
            ymin_spin.blockSignals(True)
            ymax_spin.blockSignals(True)
            ymin_spin.setValue(0)
            ymax_spin.setValue(upper_bound + padding)
            ymin_spin.blockSignals(False)
            ymax_spin.blockSignals(False)
        else:
            # Manual mode - still enforce min >= 0
            ax.set_ylim(max(0, ymin_spin.value()), ymax_spin.value())

        # Labels and title
        ax.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
        ax.set_ylabel(f"{exp_type} (ms)", fontsize=12, fontweight='bold')
        ax.set_title(f"{exp_type} Relaxation Times vs Residue Sequence", fontsize=14, fontweight='bold', pad=20)

        # Format X-axis with ticks every 10 residues
        self._format_xaxis(ax, residues)

        # Grid styling like Model-free viewer
        ax.grid(True, alpha=0.3, linestyle='--', zorder=1)
        ax.legend(loc='upper right')

        self.figure.tight_layout()
        self.canvas.draw()

    def _update_plot_multi(self):
        """Update plot with stacked subplots for multiple experiments."""
        n_experiments = len(self.experiment_data)
        if n_experiments == 0:
            self._show_blank_state()
            return

        self.figure.clear()

        # Colors for different fields
        colors = {
            'field1': '#2196F3',  # Blue
            'field2': '#FF9800',  # Orange
        }

        # Create stacked subplots
        axes = self.figure.subplots(n_experiments, 1, sharex=True)
        if n_experiments == 1:
            axes = [axes]  # Ensure iterable

        # Collect all residue numbers to determine global X-axis range
        all_residues = set()
        for key, data in self.experiment_data.items():
            df = data['df']
            if 'Success' in df.columns:
                valid_df = df[df['Success'] == 'Yes']
            else:
                valid_df = df
            if len(valid_df) > 0:
                residues = self._extract_residue_numbers(valid_df['Residue'].values)
                all_residues.update(residues)

        if not all_residues:
            self._show_blank_state()
            return

        max_residue = max(all_residues)

        for idx, (key, data) in enumerate(sorted(self.experiment_data.items())):
            ax = axes[idx]

            # Style the plot like Model-free viewer
            ax.set_facecolor(FRAME_BG_COLOR)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            df = data['df']
            exp_type = data['experiment_type']
            freq_mhz = data.get('frequency_mhz', 0)

            # Parse field from key (e.g., 'field1_t1' -> 'field1')
            field_name = key.rsplit('_', 1)[0]
            color = colors.get(field_name, '#2196F3')

            # Filter for successful fits
            if 'Success' in df.columns:
                valid_df = df[df['Success'] == 'Yes'].copy()
            else:
                valid_df = df.copy()

            if len(valid_df) == 0:
                ax.text(0.5, 0.5, f"No successful fits for {key}",
                       ha='center', va='center', fontsize=10, color=SECONDARY_TEXT,
                       transform=ax.transAxes)
                continue

            # Extract residue numbers for x-axis (actual residue numbers)
            residues = self._extract_residue_numbers(valid_df['Residue'].values)
            valid_df['ResNum'] = residues
            valid_df = valid_df.sort_values('ResNum')
            residues = valid_df['ResNum'].values

            # Find the value column - handle variations like 'T1', 'T2', 'T1 (ms)', etc.
            value_col = None
            for col in valid_df.columns:
                if col == exp_type or col.startswith(f"{exp_type} ") or col.startswith(f"{exp_type}("):
                    value_col = col
                    break

            if value_col is None:
                ax.text(0.5, 0.5, f"Column '{exp_type}' not found in data",
                       ha='center', va='center', fontsize=10, color=SECONDARY_TEXT,
                       transform=ax.transAxes)
                continue

            # Get values
            y = valid_df[value_col].values

            # Detect outliers
            outlier_mask, lower_bound, upper_bound = self._detect_outliers(y)
            normal_mask = ~outlier_mask

            # Add grey bars for missing residues
            self._add_missing_residue_bars(ax, residues, y)

            # Get errors - try variations
            err_col = None
            for col in valid_df.columns:
                if col == f"{exp_type}_err" or col.startswith(f"{exp_type}_err"):
                    err_col = col
                    break
            has_errors = err_col is not None and self.show_errors
            yerr = valid_df[err_col].values if has_errors else None

            # Plot based on exclude_outliers setting
            if self.exclude_outliers:
                x_plot = residues[normal_mask]
                y_plot = y[normal_mask]
                yerr_plot = yerr[normal_mask] if has_errors else None
            else:
                x_plot = residues
                y_plot = y
                yerr_plot = yerr

            if has_errors:
                ax.errorbar(x_plot, y_plot, yerr=yerr_plot, fmt='o', capsize=2,
                           color=color, ecolor=SECONDARY_TEXT, markersize=5, zorder=3)
            else:
                ax.plot(x_plot, y_plot, 'o', color=color, markersize=5, zorder=3)

            # Plot outliers if not excluded
            if not self.exclude_outliers and np.any(outlier_mask):
                ax.scatter(residues[outlier_mask], y[outlier_mask],
                          c='#9E9E9E', s=15, marker='x', alpha=0.5, zorder=3)

            # Y-axis - always start at 0, use per-experiment-type settings
            auto_checkbox = self.t1_auto_y_checkbox if exp_type == 'T1' else self.t2_auto_y_checkbox
            ymin_spin = self.t1_ymin_spin if exp_type == 'T1' else self.t2_ymin_spin
            ymax_spin = self.t1_ymax_spin if exp_type == 'T1' else self.t2_ymax_spin

            if auto_checkbox.isChecked():
                y_range = upper_bound - lower_bound
                padding = y_range * 0.1 if y_range > 0 else 50
                ax.set_ylim(0, upper_bound + padding)
            else:
                # Manual mode - use the spinbox values
                ax.set_ylim(max(0, ymin_spin.value()), ymax_spin.value())

            # Labels
            field_label = field_name.replace('field', 'Field ')
            freq_label = f" ({freq_mhz} MHz)" if freq_mhz else ""
            ax.set_ylabel(f"{exp_type} (ms)", fontsize=10, fontweight='bold')
            ax.set_title(f"{exp_type} - {field_label}{freq_label}", fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='--', zorder=1)

            # X-axis formatting - use residue numbers with ticks every 10
            self._format_xaxis(ax, residues)

            # X-axis label on bottom subplot only
            if idx == n_experiments - 1:
                ax.set_xlabel('Residue Number', fontsize=10, fontweight='bold')

        self.figure.tight_layout()
        self.canvas.draw()

    def _show_blank_state(self):
        """Show blank plot with instruction message."""
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.text(0.5, 0.5,
               "Load a results file to display relaxation time data",
               ha='center', va='center',
               fontsize=16, color=SECONDARY_TEXT,
               transform=ax.transAxes)
        ax.axis('off')
        self.canvas.draw()

    def _on_error_toggle(self, state):
        """Handle error bar toggle."""
        self.show_errors = self.error_checkbox.isChecked()
        self._update_plot()

    def _on_exclude_outliers_toggle(self, state):
        """Handle exclude outliers toggle."""
        self.exclude_outliers = self.exclude_outliers_checkbox.isChecked()
        self._update_plot()

    def _on_auto_y_toggle(self, exp_type):
        """Handle auto Y-axis toggle for specific experiment type."""
        if exp_type == 'T1':
            auto = self.t1_auto_y_checkbox.isChecked()
            self.y_axis_settings['T1']['auto'] = auto
            self.t1_ymin_spin.setEnabled(not auto)
            self.t1_ymax_spin.setEnabled(not auto)
        else:  # T2
            auto = self.t2_auto_y_checkbox.isChecked()
            self.y_axis_settings['T2']['auto'] = auto
            self.t2_ymin_spin.setEnabled(not auto)
            self.t2_ymax_spin.setEnabled(not auto)
        self._update_plot()

    def _on_yaxis_changed(self, exp_type):
        """Handle Y-axis limit change for specific experiment type."""
        if exp_type == 'T1':
            if not self.t1_auto_y_checkbox.isChecked():
                self.y_axis_settings['T1']['min'] = self.t1_ymin_spin.value()
                self.y_axis_settings['T1']['max'] = self.t1_ymax_spin.value()
                self._update_plot()
        else:  # T2
            if not self.t2_auto_y_checkbox.isChecked():
                self.y_axis_settings['T2']['min'] = self.t2_ymin_spin.value()
                self.y_axis_settings['T2']['max'] = self.t2_ymax_spin.value()
                self._update_plot()

    def _export_plot(self):
        """Export the current plot to a file."""
        if self.df is None:
            show_warning(self, "No Data", "No data to export.")
            return

        initial_dir = os.path.dirname(self.results_file) if self.results_file else os.getcwd()
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Plot",
            os.path.join(initial_dir, f"{self.experiment_type}_results.pdf"),
            "PDF files (*.pdf);;PNG files (*.png);;SVG files (*.svg)"
        )

        if file_path:
            try:
                self.figure.savefig(file_path, dpi=300, bbox_inches='tight')
                show_info(self, "Export Complete", f"Plot exported to:\n{file_path}")
            except Exception as e:
                show_error(self, "Export Error", f"Error exporting plot:\n{str(e)}")


def open_t1t2_results_viewer(parent=None, results_file=None, experiment_type="T1"):
    """
    Convenience function to open the T1/T2 results viewer.

    Parameters:
    -----------
    parent : QWidget, optional
        Parent widget
    results_file : str, optional
        Path to results file
    experiment_type : str
        "T1" or "T2"

    Returns:
    --------
    T1T2ResultsViewer : The viewer window instance
    """
    viewer = T1T2ResultsViewer(parent, results_file, experiment_type)
    viewer.show()
    return viewer
