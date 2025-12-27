"""
Results Viewer Module - Interactive visualization of model-free analysis results
PySide6 Qt Implementation

Displays fitted parameters (R1, R2, hetNOE, J values, S², τe, Rex) vs residue sequence
"""

import pandas as pd
import numpy as np
import os
import re
from pathlib import Path

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFrame,
    QLabel, QPushButton, QCheckBox, QRadioButton, QButtonGroup,
    QFileDialog, QMessageBox, QScrollArea, QSizePolicy
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
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER,
    SUCCESS_GREEN, WARNING_ORANGE, ERROR_RED,
    SPACING_XS, SPACING_SM, SPACING_MD,
    FONT_SIZE_SECTION_LABEL, FONT_SIZE_BODY, FONT_SIZE_SMALL
)
from gui_components import (
    create_primary_button, create_secondary_button, create_label,
    get_font, show_info, show_error, show_warning
)


class ResultsViewer(QMainWindow):
    """
    Interactive viewer for model-free analysis results

    Features:
    - Parameter selection (R1, R2, hetNOE, J values, S², τe, Rex)
    - Field selection (Field 1, Field 2, Overlay)
    - Error bar display
    - Export functionality
    """

    def __init__(self, parent=None, results_file=None, field1_freq=600.0, field2_freq=None, is_dual_field=False):
        """
        Initialize the Results Viewer

        Parameters:
        -----------
        parent : QWidget, optional
            Parent widget
        results_file : str, optional
            Path to the CSV results file (can be None for blank state)
        field1_freq : float
            Field 1 frequency in MHz (default: 600.0)
        field2_freq : float, optional
            Field 2 frequency in MHz (for dual-field analysis)
        is_dual_field : bool
            Whether dual-field analysis was performed
        """
        super().__init__(parent)

        self.setWindowTitle("Model-Free Results Viewer")
        self.setMinimumSize(1400, 900)
        self.setStyleSheet(f"background-color: {BG_COLOR};")

        # Store parameters
        self.results_file = results_file
        self.field1_freq = field1_freq
        self.field2_freq = field2_freq
        self.is_dual_field = is_dual_field

        # State variables
        self.df = None
        self.current_parameter = None
        self.current_field = "field1"
        self.show_errors = True
        self.selected_params = []

        # UI references
        self.field_section = None
        self.field_button_group = None
        self.parameter_buttons = {}
        self.show_all_btn = None
        self.error_checkbox = None

        # Load results if file provided
        if results_file:
            self._load_results()
            self._detect_dual_field()

        # Build UI
        self._create_ui()

        # Initial blank state
        self._show_blank_state()

    def _load_results(self):
        """Load results from CSV file"""
        try:
            if os.path.exists(self.results_file):
                self.df = pd.read_csv(self.results_file)
                print(f"Loaded results: {len(self.df)} residues")
                print(f"Columns: {list(self.df.columns)}")
            else:
                show_warning(
                    self, "File Not Found",
                    f"Results file not found:\n{self.results_file}\n\nYou can load a CSV file using the 'Load CSV File' button."
                )
                self.df = None
        except Exception as e:
            show_error(self, "Load Error", f"Error loading results:\n{str(e)}")
            self.df = None

    def _load_csv_file(self):
        """Open file dialog to load a CSV file"""
        initial_dir = os.path.dirname(self.results_file) if self.results_file else os.getcwd()
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Results CSV File",
            initial_dir,
            "CSV files (*.csv);;All files (*.*)"
        )

        if not file_path:
            return  # User cancelled

        try:
            self.df = pd.read_csv(file_path)
            self.results_file = file_path
            print(f"Loaded results: {len(self.df)} residues")
            print(f"Columns: {list(self.df.columns)}")

            # Detect if dual-field data is present
            self._detect_dual_field()

            # Update field controls visibility
            self._update_field_controls()

            # Update parameter buttons availability
            self._update_parameter_buttons()

            # Show blank state (waiting for parameter selection)
            self._show_blank_state()

            show_info(
                self, "File Loaded",
                f"Successfully loaded {len(self.df)} residues from:\n{os.path.basename(file_path)}"
            )
        except Exception as e:
            show_error(self, "Load Error", f"Error loading CSV file:\n{str(e)}")

    def _detect_dual_field(self):
        """Detect if loaded CSV contains dual-field data"""
        if self.df is None:
            return

        has_field2 = any('_f2' in col or '_field2' in col or '_2' in col
                        for col in self.df.columns)

        if has_field2:
            self.is_dual_field = True
            if self.field2_freq is None:
                self.field2_freq = 700.0
        else:
            self.is_dual_field = False

    def _update_field_controls(self):
        """Show or hide field controls based on dual-field detection"""
        if self.field_section is None:
            return

        # Update visibility
        self.field_section.setVisible(self.is_dual_field)

        # Update field labels if visible
        if self.is_dual_field and hasattr(self, 'field1_radio'):
            self.field1_radio.setText(f"Field 1 ({self.field1_freq} MHz)")
            self.field2_radio.setText(f"Field 2 ({self.field2_freq} MHz)")

    def _update_parameter_buttons(self):
        """Update parameter button availability based on loaded data"""
        if self.df is None:
            return

        for param, btn in self.parameter_buttons.items():
            if self._find_column_with_field(param):
                btn.setEnabled(True)
            else:
                btn.setEnabled(False)

    def _create_ui(self):
        """Create the user interface"""
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        main_layout = QVBoxLayout(central_widget)
        main_layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        main_layout.setSpacing(SPACING_MD)

        # Header
        self._create_header(main_layout)

        # Control panel
        self._create_control_panel(main_layout)

        # Plot area
        self._create_plot_area(main_layout)

    def _create_header(self, parent_layout):
        """Create header with title and buttons"""
        header_frame = QFrame()
        header_layout = QHBoxLayout(header_frame)
        header_layout.setContentsMargins(0, 0, 0, 0)

        title_label = create_label("Model-Free Results Viewer")
        title_label.setFont(get_font(20, bold=True))
        header_layout.addWidget(title_label)

        header_layout.addStretch()

        # Load CSV button
        load_csv_btn = create_primary_button("Load CSV File", clicked=self._load_csv_file)
        load_csv_btn.setFixedWidth(120)
        header_layout.addWidget(load_csv_btn)

        # Close button
        close_btn = create_secondary_button("Close", clicked=self.close)
        close_btn.setFixedWidth(80)
        header_layout.addWidget(close_btn)

        parent_layout.addWidget(header_frame)

    def _create_control_panel(self, parent_layout):
        """Create the control panel with parameter/field selectors"""
        control_frame = QFrame()
        control_frame.setStyleSheet(f"""
            QFrame {{
                background-color: {PANEL_BG_COLOR};
                border-radius: 8px;
            }}
        """)

        control_layout = QVBoxLayout(control_frame)
        control_layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        control_layout.setSpacing(SPACING_SM)

        # Parameter selection section
        param_section = QFrame()
        param_layout = QHBoxLayout(param_section)
        param_layout.setContentsMargins(0, 0, 0, 0)
        param_layout.setSpacing(SPACING_MD)

        param_label = create_label("Display Parameter:")
        param_label.setFont(get_font(FONT_SIZE_SECTION_LABEL, bold=True))
        param_layout.addWidget(param_label)

        # Show All button
        self.show_all_btn = create_primary_button("Show All", clicked=self._show_all_parameters)
        self.show_all_btn.setFixedWidth(100)
        param_layout.addWidget(self.show_all_btn)

        # Parameter groups
        self.param_groups = {
            "Relaxation": ["R1", "R2", "hetNOE"],
            "Spectral Density": ["J(0)", "J(wN)", "J(0.87wH)"],
            "Model-Free": ["S2", "te", "Rex"]
        }

        for group_name, params in self.param_groups.items():
            group_frame = QFrame()
            group_layout = QVBoxLayout(group_frame)
            group_layout.setContentsMargins(0, 0, 0, 0)
            group_layout.setSpacing(2)

            group_label = create_label(f"{group_name}:")
            group_label.setFont(get_font(FONT_SIZE_SMALL))
            group_layout.addWidget(group_label)

            button_row = QFrame()
            button_layout = QHBoxLayout(button_row)
            button_layout.setContentsMargins(0, 0, 0, 0)
            button_layout.setSpacing(SPACING_XS)

            for param in params:
                btn = create_secondary_button(param)
                btn.setFixedWidth(100)
                btn.setCheckable(True)

                # Check if parameter exists in dataframe
                if self._find_column_with_field(param):
                    btn.setEnabled(True)
                else:
                    btn.setEnabled(False)

                btn.clicked.connect(lambda checked, p=param: self._toggle_parameter(p))
                button_layout.addWidget(btn)
                self.parameter_buttons[param] = btn

            button_layout.addStretch()
            group_layout.addWidget(button_row)
            param_layout.addWidget(group_frame)

        param_layout.addStretch()
        control_layout.addWidget(param_section)

        # Field selection section (only visible if dual-field)
        self.field_section = QFrame()
        field_layout = QHBoxLayout(self.field_section)
        field_layout.setContentsMargins(0, 0, 0, 0)
        field_layout.setSpacing(SPACING_SM)

        field_label = create_label("Field:")
        field_label.setFont(get_font(FONT_SIZE_SECTION_LABEL, bold=True))
        field_layout.addWidget(field_label)

        self.field_button_group = QButtonGroup(self)
        self.field_button_group.buttonClicked.connect(self._on_field_change)

        self.field1_radio = QRadioButton(f"Field 1 ({self.field1_freq} MHz)")
        self.field1_radio.setFont(get_font(FONT_SIZE_BODY))
        self.field1_radio.setChecked(True)
        self.field_button_group.addButton(self.field1_radio, 1)
        field_layout.addWidget(self.field1_radio)

        self.field2_radio = QRadioButton(f"Field 2 ({self.field2_freq or 700.0} MHz)")
        self.field2_radio.setFont(get_font(FONT_SIZE_BODY))
        self.field_button_group.addButton(self.field2_radio, 2)
        field_layout.addWidget(self.field2_radio)

        self.overlay_radio = QRadioButton("Overlay Both")
        self.overlay_radio.setFont(get_font(FONT_SIZE_BODY))
        self.field_button_group.addButton(self.overlay_radio, 3)
        field_layout.addWidget(self.overlay_radio)

        field_layout.addStretch()
        self.field_section.setVisible(self.is_dual_field)
        control_layout.addWidget(self.field_section)

        # Options section
        options_section = QFrame()
        options_layout = QHBoxLayout(options_section)
        options_layout.setContentsMargins(0, 0, 0, 0)

        self.error_checkbox = QCheckBox("Show error bars")
        self.error_checkbox.setFont(get_font(FONT_SIZE_BODY))
        self.error_checkbox.setChecked(True)
        self.error_checkbox.stateChanged.connect(self._on_error_checkbox_change)
        options_layout.addWidget(self.error_checkbox)

        options_layout.addStretch()
        control_layout.addWidget(options_section)

        parent_layout.addWidget(control_frame)

    def _create_plot_area(self, parent_layout):
        """Create the matplotlib plot area"""
        plot_frame = QFrame()
        plot_frame.setStyleSheet(f"""
            QFrame {{
                background-color: {PANEL_BG_COLOR};
                border-radius: 8px;
            }}
        """)

        plot_layout = QVBoxLayout(plot_frame)
        plot_layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        # Create matplotlib figure
        self.figure = Figure(figsize=(14, 7), facecolor=PANEL_BG_COLOR)
        self.ax = self.figure.add_subplot(111)

        # Style the plot
        self.ax.set_facecolor(FRAME_BG_COLOR)
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)

        # Create canvas
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        plot_layout.addWidget(self.canvas)

        # Add toolbar (stored as instance variable to prevent garbage collection)
        self.toolbar = NavigationToolbar2QT(self.canvas, plot_frame)
        plot_layout.addWidget(self.toolbar)

        parent_layout.addWidget(plot_frame, 1)  # stretch factor 1

    def _find_column_with_field(self, param_name, field=None):
        """
        Find column name in dataframe that matches parameter and field

        Parameters:
        -----------
        param_name : str
            Parameter name (e.g., 'R1', 'S2')
        field : str, optional
            Field specification ('field1', 'field2', or None for any)

        Returns:
        --------
        str or None : Column name if found
        """
        if self.df is None:
            return None

        variations = {
            "R1": ["R1"],
            "R2": ["R2"],
            "hetNOE": ["hetNOE", "NOE", "noe"],
            "J(0)": ["J0", "J_0", "J(0)"],
            "J(wN)": ["JwN", "J_wN", "J(wN)", "JN"],
            "J(wH)": ["JwH", "J_wH", "J(wH)", "JH"],
            "J(0.87wH)": ["JwH_087", "J0.87wH", "J_0.87wH", "J(0.87wH)", "JwH0.87"],
            "S2": ["S2", "S²", "order_parameter"],
            "te": ["te", "tau_e", "taue"],
            "Rex": ["Rex", "R_ex"],
            "tc": ["tc", "tau_c", "tauc"]
        }

        base_names = [param_name]
        if param_name in variations:
            base_names.extend(variations[param_name])

        if field == "field1":
            field_suffixes = ["_f1", "_field1", "_1", ""]
        elif field == "field2":
            field_suffixes = ["_f2", "_field2", "_2", ""]
        else:
            field_suffixes = ["", "_f1", "_f2", "_field1", "_field2", "_1", "_2"]

        for base in base_names:
            for suffix in field_suffixes:
                col_name = f"{base}{suffix}"
                if col_name in self.df.columns:
                    return col_name

        return None

    def _find_column(self, param_name):
        """Find column name based on current field selection"""
        global_params = ["S2", "tc", "te"]

        if param_name in global_params:
            return self._find_column_with_field(param_name, None)

        if self.current_field == "field1":
            return self._find_column_with_field(param_name, "field1")
        elif self.current_field == "field2":
            return self._find_column_with_field(param_name, "field2")
        else:
            return self._find_column_with_field(param_name, "field1")

    def _show_blank_state(self):
        """Show blank plot with instruction message"""
        self.ax.clear()

        if self.df is None:
            message = "No data loaded.\n\nClick 'Load CSV File' to load results."
        else:
            message = "Select a parameter above to display results"

        self.ax.text(0.5, 0.5, message,
                    ha='center', va='center',
                    fontsize=16, color=SECONDARY_TEXT,
                    transform=self.ax.transAxes)
        self.ax.axis('off')
        self.canvas.draw()

    def _toggle_parameter(self, param):
        """Toggle parameter selection (multi-select within same group)"""
        if param in self.selected_params:
            self.selected_params.remove(param)
        else:
            param_group = self._get_param_group(param)
            # Remove params from different groups
            self.selected_params = [p for p in self.selected_params
                                   if self._get_param_group(p) == param_group]
            self.selected_params.append(param)

        # Update button visual states
        for p, btn in self.parameter_buttons.items():
            btn.setChecked(p in self.selected_params)

        # Reset Show All button
        self.show_all_btn.setProperty("class", "primary")
        self.show_all_btn.style().unpolish(self.show_all_btn)
        self.show_all_btn.style().polish(self.show_all_btn)

        # Update plot
        if self.selected_params:
            self.current_parameter = self.selected_params[0]
            self._update_plot_stacked()
        else:
            self._show_blank_state()

    def _get_param_group(self, param):
        """Get the group name for a parameter"""
        for group_name, params in self.param_groups.items():
            if param in params:
                return group_name
        return None

    def _show_all_parameters(self):
        """Display all available parameters in subplots"""
        self.current_parameter = "ALL"
        self.selected_params = []

        # Reset parameter button states
        for btn in self.parameter_buttons.values():
            btn.setChecked(False)

        # Highlight Show All button
        self.show_all_btn.setStyleSheet(f"background-color: {SUCCESS_GREEN}; color: white;")

        self._update_plot_all()

    def _on_field_change(self, button):
        """Handle field selection change"""
        btn_id = self.field_button_group.id(button)
        if btn_id == 1:
            self.current_field = "field1"
        elif btn_id == 2:
            self.current_field = "field2"
        else:
            self.current_field = "overlay"

        if self.selected_params:
            self._update_plot_stacked()
        elif self.current_parameter == "ALL":
            self._update_plot_all()
        else:
            self._update_plot()

    def _on_error_checkbox_change(self, state):
        """Handle error checkbox state change"""
        self.show_errors = self.error_checkbox.isChecked()
        if self.selected_params:
            self._update_plot_stacked()
        elif self.current_parameter == "ALL":
            self._update_plot_all()
        elif self.current_parameter:
            self._update_plot()

    def _update_plot(self):
        """Update the plot with current selection"""
        if self.current_parameter is None or self.current_parameter == "ALL":
            self._show_blank_state()
            return

        self.ax.clear()

        residue_col = None
        for col in ['Residue', 'residue', 'Res', 'res', 'ResidueNum']:
            if col in self.df.columns:
                residue_col = col
                break

        if residue_col is None:
            residues = self.df.index.values
        else:
            residues = self._extract_residue_numbers(self.df[residue_col].values)

        if self.is_dual_field and self.current_field == "overlay":
            self._plot_overlay(residues)
        else:
            self._plot_single_field(residues)

        self.figure.tight_layout()
        self.canvas.draw()

    def _update_plot_stacked(self):
        """Display multiple selected parameters in stacked vertical plots"""
        if not self.selected_params:
            self._show_blank_state()
            return

        self.figure.clear()
        n_params = len(self.selected_params)

        residue_col = None
        for col in ['Residue', 'residue', 'Res', 'res', 'ResidueNum']:
            if col in self.df.columns:
                residue_col = col
                break

        if residue_col is None:
            residues = self.df.index.values
        else:
            residues = self._extract_residue_numbers(self.df[residue_col].values)

        for idx, param in enumerate(self.selected_params):
            ax = self.figure.add_subplot(n_params, 1, idx + 1)

            if self.is_dual_field and self.current_field == "overlay":
                self._plot_param_overlay(ax, param, residues)
            else:
                self._plot_param_single(ax, param, residues)

            if idx == n_params - 1:
                ax.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
                self._format_xaxis_for_ax(ax, residues)
                ax.tick_params(labelbottom=True)
            else:
                ax.set_xlabel('')
                self._format_xaxis_for_ax(ax, residues)
                ax.tick_params(labelbottom=False)

            ax.set_ylabel(self._get_ylabel(param), fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='--', zorder=1)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            ylim = ax.get_ylim()
            ax.set_ylim(0, ylim[1] * 1.1)

        group_name = self._get_param_group(self.selected_params[0])
        title = f"{group_name} Parameters vs Residue Sequence"
        if self.is_dual_field and self.current_field == "overlay":
            title += " (Overlay)"
        self.figure.suptitle(title, fontsize=14, fontweight='bold')

        self.figure.tight_layout(rect=[0, 0, 1, 0.97])
        self.canvas.draw()

    def _plot_param_single(self, ax, param, residues):
        """Plot a single parameter on given axis (single field)"""
        col_name = self._find_column_with_field(
            param,
            "field1" if self.current_field == "field1" else
            "field2" if self.current_field == "field2" else None
        )

        if col_name is None:
            ax.text(0.5, 0.5, f"'{param}' not found",
                   ha='center', va='center', transform=ax.transAxes)
            ax.axis('off')
            return

        error_col = None
        for suffix in ['err', '_err', '_error', 'Error']:
            potential_err = f"{col_name}{suffix}"
            if potential_err in self.df.columns:
                error_col = potential_err
                break

        values = self.df[col_name].values
        errors = self.df[error_col].values if error_col and self.show_errors else None

        self._add_missing_residue_bars_for_ax(ax, residues, values)

        if errors is not None and self.show_errors:
            ax.errorbar(residues, values, yerr=errors,
                       fmt='o', markersize=5, capsize=3,
                       color=PRIMARY_BUTTON_BG, ecolor=SECONDARY_TEXT, zorder=3)
        else:
            ax.plot(residues, values, 'o', markersize=5,
                   color=PRIMARY_BUTTON_BG, zorder=3)

    def _plot_param_overlay(self, ax, param, residues):
        """Plot a single parameter on given axis (overlay both fields)"""
        col_field1 = self._find_column_with_field(param, "field1")
        col_field2 = self._find_column_with_field(param, "field2")

        all_values = []
        if col_field1:
            all_values.extend(self.df[col_field1].dropna().values)
        if col_field2:
            all_values.extend(self.df[col_field2].dropna().values)

        if all_values:
            self._add_missing_residue_bars_for_ax(ax, residues, np.array(all_values))

        colors = ['#5B9EE5', '#E8554E']

        for col_name, field_label, color in [
            (col_field1, f"F1 ({self.field1_freq} MHz)", colors[0]),
            (col_field2, f"F2 ({self.field2_freq} MHz)", colors[1])
        ]:
            if col_name is None:
                continue

            error_col = None
            for suffix in ['err', '_err', '_error', 'Error']:
                potential_err = f"{col_name}{suffix}"
                if potential_err in self.df.columns:
                    error_col = potential_err
                    break

            values = self.df[col_name].values
            errors = self.df[error_col].values if error_col and self.show_errors else None

            if errors is not None and self.show_errors:
                ax.errorbar(residues, values, yerr=errors,
                           fmt='o', markersize=4, capsize=2,
                           color=color, ecolor=color, alpha=0.7,
                           label=field_label, zorder=3)
            else:
                ax.plot(residues, values, 'o', markersize=4,
                       color=color, alpha=0.7,
                       label=field_label, zorder=3)

        ax.legend(loc='best', fontsize=8)

    def _plot_single_field(self, residues):
        """Plot data for a single field"""
        col_name = self._find_column(self.current_parameter)
        if col_name is None:
            self.ax.text(0.5, 0.5,
                        f"Parameter '{self.current_parameter}' not found in results",
                        ha='center', va='center',
                        fontsize=14, color=ERROR_RED,
                        transform=self.ax.transAxes)
            self.ax.axis('off')
            return

        error_col = None
        for suffix in ['err', '_err', '_error', 'Error']:
            potential_err = f"{col_name}{suffix}"
            if potential_err in self.df.columns:
                error_col = potential_err
                break

        values = self.df[col_name].values
        errors = self.df[error_col].values if error_col and self.show_errors else None

        self._add_missing_residue_bars(residues, values)

        field_label = f"{self.current_parameter}"
        if self.is_dual_field:
            field_name = f"Field {self.current_field[-1]}" if "field" in self.current_field else self.current_field
            field_label = f"{self.current_parameter} ({field_name})"

        if errors is not None and self.show_errors:
            self.ax.errorbar(residues, values, yerr=errors,
                           fmt='o', markersize=6, capsize=3,
                           color=PRIMARY_BUTTON_BG, ecolor=SECONDARY_TEXT,
                           label=field_label, zorder=3)
        else:
            self.ax.plot(residues, values, 'o', markersize=6,
                       color=PRIMARY_BUTTON_BG,
                       label=field_label, zorder=3)

        self.ax.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
        self.ax.set_ylabel(self._get_ylabel(self.current_parameter),
                          fontsize=12, fontweight='bold')
        self.ax.set_title(f"{self.current_parameter} vs Residue Sequence",
                         fontsize=14, fontweight='bold', pad=20)
        self.ax.grid(True, alpha=0.3, linestyle='--', zorder=1)
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)

        self._format_xaxis(residues)

        ylim = self.ax.get_ylim()
        self.ax.set_ylim(0, ylim[1] * 1.1)

    def _plot_overlay(self, residues):
        """Plot overlay of both fields"""
        col_field1 = self._find_column_with_field(self.current_parameter, "field1")
        col_field2 = self._find_column_with_field(self.current_parameter, "field2")

        all_values = []
        if col_field1:
            all_values.extend(self.df[col_field1].dropna().values)
        if col_field2:
            all_values.extend(self.df[col_field2].dropna().values)

        if all_values:
            self._add_missing_residue_bars(residues, np.array(all_values))

        colors = ['#5B9EE5', '#E8554E']

        for col_name, field_label, color in [
            (col_field1, f"Field 1 ({self.field1_freq} MHz)", colors[0]),
            (col_field2, f"Field 2 ({self.field2_freq} MHz)", colors[1])
        ]:
            if col_name is None:
                continue

            error_col = None
            for suffix in ['err', '_err', '_error', 'Error']:
                potential_err = f"{col_name}{suffix}"
                if potential_err in self.df.columns:
                    error_col = potential_err
                    break

            values = self.df[col_name].values
            errors = self.df[error_col].values if error_col and self.show_errors else None

            if errors is not None and self.show_errors:
                self.ax.errorbar(residues, values, yerr=errors,
                               fmt='o', markersize=6, capsize=3,
                               color=color, ecolor=color, alpha=0.7,
                               label=field_label, zorder=3)
            else:
                self.ax.plot(residues, values, 'o', markersize=6,
                           color=color, alpha=0.7,
                           label=field_label, zorder=3)

        self.ax.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
        self.ax.set_ylabel(self._get_ylabel(self.current_parameter),
                          fontsize=12, fontweight='bold')
        self.ax.set_title(f"{self.current_parameter} vs Residue Sequence (Overlay)",
                         fontsize=14, fontweight='bold', pad=20)
        self.ax.grid(True, alpha=0.3, linestyle='--', zorder=1)
        self.ax.legend(loc='best')
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)

        self._format_xaxis(residues)

        ylim = self.ax.get_ylim()
        self.ax.set_ylim(0, ylim[1] * 1.1)

    def _update_plot_all(self):
        """Display all parameters in a grid of subplots"""
        self.figure.clear()

        available_params = list(self.parameter_buttons.keys())
        n_params = len(available_params)

        if n_params == 0:
            self._show_blank_state()
            return

        n_cols = min(3, n_params)
        n_rows = (n_params + n_cols - 1) // n_cols

        residue_col = None
        for col in ['Residue', 'residue', 'Res', 'res', 'ResidueNum']:
            if col in self.df.columns:
                residue_col = col
                break

        if residue_col is None:
            residues = self.df.index.values
        else:
            residues = self._extract_residue_numbers(self.df[residue_col].values)

        for idx, param in enumerate(available_params):
            ax = self.figure.add_subplot(n_rows, n_cols, idx + 1)

            if self.is_dual_field and self.current_field == "overlay":
                self._plot_param_overlay_small(ax, param, residues)
            else:
                self._plot_param_single_small(ax, param, residues)

            min_residue = 1
            max_residue = int(np.max(residues))
            ax.set_xlim(min_residue - 0.5, max_residue + 0.5)

            if max_residue >= 20:
                tick_step = 20 if max_residue > 50 else 10
                ticks = [1] + list(range(tick_step, max_residue + 1, tick_step))
                if max_residue not in ticks:
                    ticks.append(max_residue)
                ax.set_xticks(ticks)
                ax.set_xticklabels([str(t) for t in ticks])
            else:
                ticks = [1] + list(range(5, max_residue + 1, 5))
                if max_residue not in ticks:
                    ticks.append(max_residue)
                ax.set_xticks(ticks)

            ax.set_title(param, fontsize=10, fontweight='bold')
            ax.set_xlabel('Residue', fontsize=8)
            ax.set_ylabel(self._get_ylabel(param), fontsize=8)
            ax.grid(True, alpha=0.2, linestyle='--', zorder=1)
            ax.tick_params(labelsize=7)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            ylim = ax.get_ylim()
            ax.set_ylim(0, ylim[1] * 1.1)

        self.figure.tight_layout()
        self.canvas.draw()

    def _plot_param_single_small(self, ax, param, residues):
        """Plot a single parameter on small subplot (single field)"""
        col_name = self._find_column_with_field(
            param,
            "field1" if self.current_field == "field1" else
            "field2" if self.current_field == "field2" else None
        )

        if col_name is None:
            ax.text(0.5, 0.5, f"'{param}' not found",
                   ha='center', va='center', transform=ax.transAxes, fontsize=8)
            ax.axis('off')
            return

        error_col = None
        for suffix in ['err', '_err', '_error', 'Error']:
            potential_err = f"{col_name}{suffix}"
            if potential_err in self.df.columns:
                error_col = potential_err
                break

        values = self.df[col_name].values
        errors = self.df[error_col].values if error_col and self.show_errors else None

        self._add_missing_residue_bars_for_ax(ax, residues, values)

        if errors is not None and self.show_errors:
            ax.errorbar(residues, values, yerr=errors,
                       fmt='o', markersize=3, capsize=2,
                       color=PRIMARY_BUTTON_BG, ecolor=SECONDARY_TEXT, zorder=3)
        else:
            ax.plot(residues, values, 'o', markersize=3,
                   color=PRIMARY_BUTTON_BG, zorder=3)

    def _plot_param_overlay_small(self, ax, param, residues):
        """Plot a single parameter on small subplot (overlay both fields)"""
        col_field1 = self._find_column_with_field(param, "field1")
        col_field2 = self._find_column_with_field(param, "field2")

        all_values = []
        if col_field1:
            all_values.extend(self.df[col_field1].dropna().values)
        if col_field2:
            all_values.extend(self.df[col_field2].dropna().values)

        if all_values:
            self._add_missing_residue_bars_for_ax(ax, residues, np.array(all_values))

        colors = ['#5B9EE5', '#E8554E']

        for col_name, color in [(col_field1, colors[0]), (col_field2, colors[1])]:
            if col_name is None:
                continue

            error_col = None
            for suffix in ['err', '_err', '_error', 'Error']:
                potential_err = f"{col_name}{suffix}"
                if potential_err in self.df.columns:
                    error_col = potential_err
                    break

            values = self.df[col_name].values
            errors = self.df[error_col].values if error_col and self.show_errors else None

            if errors is not None and self.show_errors:
                ax.errorbar(residues, values, yerr=errors,
                           fmt='o', markersize=2, capsize=1,
                           color=color, ecolor=color, alpha=0.7, zorder=3)
            else:
                ax.plot(residues, values, 'o', markersize=2,
                       color=color, alpha=0.7, zorder=3)

    def _extract_residue_numbers(self, residue_labels):
        """
        Extract numeric residue numbers from residue labels

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

    def _add_missing_residue_bars_for_ax(self, ax, residues, values):
        """Add grey bars for missing residues to a specific axis"""
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

    def _format_xaxis_for_ax(self, ax, residues):
        """Format x-axis with ticks every 10 residues for a specific axis"""
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

    def _add_missing_residue_bars(self, residues, values):
        """Add grey bars for missing residues"""
        min_residue = 1
        max_residue = int(np.max(residues))

        present_residues = set(residues)
        all_residues = set(range(min_residue, max_residue + 1))
        missing_residues = sorted(all_residues - present_residues)

        if not missing_residues:
            return

        valid_values = values[~np.isnan(values)]
        if len(valid_values) == 0:
            return

        for res in missing_residues:
            self.ax.axvspan(res - 0.5, res + 0.5,
                          facecolor='lightgrey', alpha=0.3, zorder=0)

    def _format_xaxis(self, residues):
        """Format x-axis with ticks every 10 residues"""
        min_residue = 1
        max_residue = int(np.max(residues))

        self.ax.set_xlim(min_residue - 0.5, max_residue + 0.5)

        if max_residue >= 10:
            first_tick = 10
            ticks = list(range(first_tick, max_residue + 1, 10))
            if 1 not in ticks:
                ticks = [1] + ticks
            if max_residue not in ticks and max_residue - ticks[-1] > 2:
                ticks.append(max_residue)
        else:
            ticks = list(range(1, max_residue + 1))

        self.ax.set_xticks(ticks)
        self.ax.set_xticklabels([str(t) for t in ticks])

    def _get_ylabel(self, param):
        """Get appropriate y-axis label for parameter"""
        labels = {
            "R1": "R₁ (s⁻¹)",
            "R2": "R₂ (s⁻¹)",
            "hetNOE": "hetNOE",
            "J(0)": "J(0) (ns/rad)",
            "J(wN)": "J(ωₙ) (ns/rad)",
            "J(wH)": "J(ωₕ) (ns/rad)",
            "J(0.87wH)": "J(0.87ωₕ) (ns/rad)",
            "S2": "S²",
            "tc": "τc (ns)",
            "te": "τₑ (ps)",
            "Rex": "Rₑₓ (s⁻¹)"
        }
        return labels.get(param, param)


def open_results_viewer(parent=None, results_file=None, field1_freq=600.0, field2_freq=None, is_dual_field=False):
    """
    Convenience function to open the results viewer

    Parameters:
    -----------
    parent : QWidget, optional
        Parent widget
    results_file : str
        Path to results CSV file
    field1_freq : float
        Field 1 frequency in MHz
    field2_freq : float, optional
        Field 2 frequency in MHz
    is_dual_field : bool
        Whether dual-field analysis was performed

    Returns:
    --------
    ResultsViewer : The viewer window instance
    """
    viewer = ResultsViewer(parent, results_file, field1_freq, field2_freq, is_dual_field)
    viewer.show()
    return viewer
