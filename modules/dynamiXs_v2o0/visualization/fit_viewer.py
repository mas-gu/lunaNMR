"""
T1/T2 Fit Viewer Module - Interactive visualization of T1/T2 fitting results
PySide6 Qt Implementation

Displays exponential decay fits with measured data points and fitted curves
"""

import json
import os
import re
import numpy as np
from pathlib import Path

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFrame,
    QLabel, QPushButton, QCheckBox, QRadioButton, QButtonGroup,
    QFileDialog, QMessageBox, QScrollArea, QSizePolicy, QSplitter,
    QListWidget, QListWidgetItem, QAbstractItemView
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
    SPACING_XS, SPACING_SM, SPACING_MD,
    FONT_SIZE_SECTION_LABEL, FONT_SIZE_BODY, FONT_SIZE_SMALL
)
from gui_components import (
    create_primary_button, create_secondary_button, create_label,
    get_font, show_info, show_error, show_warning
)


class FitViewer(QMainWindow):
    """
    Interactive viewer for T1/T2 fitting results

    Features:
    - JSON folder selection
    - Multi-residue selection (max 4 plots)
    - Field selection (F1, F2, Overlay)
    - Measurement type selection (T1, T2)
    - 4x1 vertical subplot layout
    """

    def __init__(self, parent=None, json_folder=None, json_folders=None):
        """
        Initialize the Fit Viewer

        Parameters:
        -----------
        parent : QWidget, optional
            Parent widget
        json_folder : str, optional
            Path to JSON data folder (legacy single folder)
        json_folders : list, optional
            List of JSON data folder paths (for multi-experiment support)
        """
        super().__init__(parent)

        self.setWindowTitle("T1/T2 Fit Viewer")
        self.setMinimumSize(1400, 900)
        self.setStyleSheet(f"background-color: {BG_COLOR};")

        # Data storage
        self.json_folders = json_folders or ([json_folder] if json_folder else [])
        self.json_folder = self.json_folders[0] if self.json_folders else None  # For UI display
        self.data = {}  # {field_name_type: json_data}
        self.available_residues = []
        self.selected_residues = []

        # UI state
        self.field_mode = "field1"
        self.show_t1 = True
        self.show_t2 = True

        # Build UI
        self._create_ui()

        # Load data from all folders
        if self.json_folders:
            self._load_json_folders(self.json_folders)

    def _create_ui(self):
        """Create the user interface"""
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        main_layout = QHBoxLayout(central_widget)
        main_layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        main_layout.setSpacing(SPACING_SM)

        # Left panel - Plots (75% width)
        left_panel = QFrame()
        left_panel.setStyleSheet(f"""
            QFrame {{
                background-color: {PANEL_BG_COLOR};
                border-radius: 8px;
            }}
        """)
        self._create_plot_panel(left_panel)
        main_layout.addWidget(left_panel, 3)  # stretch factor 3

        # Right panel - Navigator (25% width)
        right_panel = QFrame()
        right_panel.setStyleSheet(f"""
            QFrame {{
                background-color: {PANEL_BG_COLOR};
                border-radius: 8px;
            }}
        """)
        right_panel.setFixedWidth(350)
        self._create_navigator_panel(right_panel)
        main_layout.addWidget(right_panel, 0)

    def _create_plot_panel(self, parent):
        """Create the plot display panel"""
        layout = QVBoxLayout(parent)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        layout.setSpacing(SPACING_SM)

        # Header
        title_label = create_label("Fit Visualization")
        title_label.setFont(get_font(18, bold=True))
        layout.addWidget(title_label)

        # Plot area
        plot_frame = QFrame()
        plot_frame.setStyleSheet(f"background-color: {FRAME_BG_COLOR}; border-radius: 8px;")
        plot_layout = QVBoxLayout(plot_frame)
        plot_layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)

        # Create matplotlib figure (4x1 layout)
        self.figure = Figure(figsize=(10, 12), facecolor=PANEL_BG_COLOR)
        self.axes = []

        # Create canvas (stored as instance variable to prevent garbage collection)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        plot_layout.addWidget(self.canvas)

        # Toolbar (stored as instance variable to prevent garbage collection)
        self.toolbar = NavigationToolbar2QT(self.canvas, plot_frame)
        plot_layout.addWidget(self.toolbar)

        layout.addWidget(plot_frame, 1)

        # Initial blank state
        self._show_blank_state()

    def _create_navigator_panel(self, parent):
        """Create the navigator control panel"""
        layout = QVBoxLayout(parent)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        layout.setSpacing(SPACING_SM)

        # Header
        title_label = create_label("Peak Navigator")
        title_label.setFont(get_font(16, bold=True))
        layout.addWidget(title_label)

        # Scrollable content
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        scroll_area.setStyleSheet("background-color: transparent;")

        scroll_content = QWidget()
        scroll_layout = QVBoxLayout(scroll_content)
        scroll_layout.setContentsMargins(0, 0, 0, 0)
        scroll_layout.setSpacing(SPACING_MD)

        # === JSON Folder Selection ===
        folder_section = self._create_section_frame("JSON Data Folder:")
        folder_layout = folder_section.layout()

        self.folder_label = QLabel(self.json_folder if self.json_folder else "No folder selected")
        self.folder_label.setFont(get_font(FONT_SIZE_SMALL))
        self.folder_label.setStyleSheet(f"""
            QLabel {{
                background-color: {BG_COLOR};
                color: {SECONDARY_TEXT};
                padding: 6px;
                border-radius: 4px;
            }}
        """)
        self.folder_label.setWordWrap(True)
        folder_layout.addWidget(self.folder_label)

        browse_btn = create_primary_button("Browse Folder", clicked=self._browse_json_folder)
        browse_btn.setFixedWidth(280)
        folder_layout.addWidget(browse_btn)

        scroll_layout.addWidget(folder_section)

        # === Field Selection ===
        field_section = self._create_section_frame("Field:")
        field_layout = field_section.layout()

        self.field_button_group = QButtonGroup(self)
        self.field_button_group.buttonClicked.connect(self._on_field_change)

        self.field1_radio = QRadioButton("Field 1")
        self.field1_radio.setFont(get_font(FONT_SIZE_BODY))
        self.field1_radio.setChecked(True)
        self.field_button_group.addButton(self.field1_radio, 1)
        field_layout.addWidget(self.field1_radio)

        self.field2_radio = QRadioButton("Field 2")
        self.field2_radio.setFont(get_font(FONT_SIZE_BODY))
        self.field_button_group.addButton(self.field2_radio, 2)
        field_layout.addWidget(self.field2_radio)

        self.overlay_radio = QRadioButton("Overlay Both")
        self.overlay_radio.setFont(get_font(FONT_SIZE_BODY))
        self.field_button_group.addButton(self.overlay_radio, 3)
        field_layout.addWidget(self.overlay_radio)

        scroll_layout.addWidget(field_section)

        # === Measurement Type ===
        type_section = self._create_section_frame("Measurement Type:")
        type_layout = type_section.layout()

        self.t1_checkbox = QCheckBox("T1")
        self.t1_checkbox.setFont(get_font(FONT_SIZE_BODY))
        self.t1_checkbox.setChecked(True)
        self.t1_checkbox.stateChanged.connect(self._update_max_residues)
        type_layout.addWidget(self.t1_checkbox)

        self.t2_checkbox = QCheckBox("T2")
        self.t2_checkbox.setFont(get_font(FONT_SIZE_BODY))
        self.t2_checkbox.setChecked(True)
        self.t2_checkbox.stateChanged.connect(self._update_max_residues)
        type_layout.addWidget(self.t2_checkbox)

        scroll_layout.addWidget(type_section)

        # === Residue Selection ===
        residue_section = self._create_section_frame("Select Residue:")
        residue_layout = residue_section.layout()

        # Instruction label
        hint_label = create_label("Click or use ↑↓ arrow keys to navigate")
        hint_label.setFont(get_font(FONT_SIZE_SMALL))
        hint_label.setStyleSheet(f"color: {SECONDARY_TEXT};")
        residue_layout.addWidget(hint_label)

        # Residue list widget (replaces checkboxes)
        self.residue_list = QListWidget()
        self.residue_list.setMinimumHeight(300)
        self.residue_list.setStyleSheet(f"""
            QListWidget {{
                background-color: {BG_COLOR};
                border: none;
                border-radius: 4px;
            }}
            QListWidget::item {{
                padding: 6px 8px;
                border-radius: 4px;
            }}
            QListWidget::item:selected {{
                background-color: {PRIMARY_BUTTON_BG};
                color: white;
            }}
            QListWidget::item:hover:!selected {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
        """)
        self.residue_list.setFont(get_font(FONT_SIZE_BODY))
        # Single selection mode - simple navigation
        self.residue_list.setSelectionMode(QAbstractItemView.SingleSelection)
        # Connect both click and keyboard navigation (arrow keys)
        self.residue_list.itemClicked.connect(self._on_residue_clicked)
        self.residue_list.currentRowChanged.connect(self._on_residue_row_changed)
        residue_layout.addWidget(self.residue_list)

        scroll_layout.addWidget(residue_section)

        # === Action Buttons ===
        button_frame = QFrame()
        button_layout = QHBoxLayout(button_frame)
        button_layout.setContentsMargins(0, 0, 0, 0)
        button_layout.setSpacing(SPACING_XS)

        clear_btn = create_secondary_button("Clear Selection", clicked=self._clear_selection)
        clear_btn.setFixedWidth(280)
        button_layout.addWidget(clear_btn)

        button_layout.addStretch()
        scroll_layout.addWidget(button_frame)

        scroll_layout.addStretch()
        scroll_area.setWidget(scroll_content)
        layout.addWidget(scroll_area, 1)

    def _create_section_frame(self, title):
        """Create a styled section frame with title"""
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

    def _browse_json_folder(self):
        """Browse for JSON data folder"""
        initial_dir = self.json_folder if self.json_folder else os.getcwd()
        folder_path = QFileDialog.getExistingDirectory(
            self,
            "Select JSON Data Folder",
            initial_dir
        )

        if folder_path:
            self._load_json_folder(folder_path)

    def _load_json_folders(self, folder_paths):
        """Load JSON files from multiple folders."""
        try:
            self.data = {}
            self.available_residues = []
            residues_set = set()

            # Update folder label with count
            if len(folder_paths) == 1:
                self.folder_label.setText(folder_paths[0])
            else:
                self.folder_label.setText(f"{len(folder_paths)} folders")

            # Load from each folder
            for folder_path in folder_paths:
                if not os.path.exists(folder_path):
                    continue

                # Scan for JSON files
                json_files = list(Path(folder_path).glob("*_fit_data.json"))

                for json_file in json_files:
                    try:
                        with open(json_file, 'r') as f:
                            data = json.load(f)

                        # Extract file key (e.g., "field1_T1")
                        filename = json_file.stem
                        file_key = filename.replace("_fit_data", "")

                        self.data[file_key] = data

                        # Collect residues
                        for fit in data['fits']:
                            residues_set.add(fit['residue'])

                    except Exception as e:
                        print(f"Error loading {json_file}: {e}")

            if not self.data:
                show_warning(self, "No Data", "No fit data JSON files found in the folders.")
                self._populate_residue_list()
                return

            # Sort residues
            self.available_residues = sorted(list(residues_set), key=self._residue_sort_key)

            # Update UI
            self._update_field_controls()
            self._populate_residue_list()

            # Show status
            n_files = len(self.data)
            n_residues = len(self.available_residues)
            show_info(self, "Data Loaded",
                     f"Loaded {n_files} dataset(s) with {n_residues} residue(s)")

        except Exception as e:
            show_error(self, "Load Error", f"Error loading JSON folders:\n{str(e)}")

    def _load_json_folder(self, folder_path):
        """Load JSON files from a single folder."""
        # Use the multi-folder method with a single folder
        self._load_json_folders([folder_path])

    def _residue_sort_key(self, residue_str):
        """Extract numeric part of residue for sorting"""
        numbers = re.findall(r'\d+', residue_str)
        return int(numbers[0]) if numbers else 0

    def _update_field_controls(self):
        """Enable/disable field controls based on available data"""
        has_field1 = any('field1' in key for key in self.data.keys())
        has_field2 = any('field2' in key for key in self.data.keys())

        # Enable/disable radio buttons
        self.field1_radio.setEnabled(has_field1)
        self.field2_radio.setEnabled(has_field2)
        self.overlay_radio.setEnabled(has_field1 and has_field2)

        # Set default selection
        if has_field1 and has_field2:
            self.overlay_radio.setChecked(True)
            self.field_mode = "overlay"
        elif has_field1:
            self.field1_radio.setChecked(True)
            self.field_mode = "field1"
        elif has_field2:
            self.field2_radio.setChecked(True)
            self.field_mode = "field2"

    def _populate_residue_list(self):
        """Populate the residue selection list"""
        self.residue_list.clear()

        if not self.available_residues:
            # Show placeholder item
            placeholder = QListWidgetItem("No data available")
            placeholder.setFlags(Qt.NoItemFlags)  # Not selectable
            self.residue_list.addItem(placeholder)
            return

        # Add residue items
        for residue in self.available_residues:
            item = QListWidgetItem(residue)
            item.setData(Qt.UserRole, residue)
            self.residue_list.addItem(item)

        # Select first residue by default
        if self.available_residues:
            self.residue_list.setCurrentRow(0)

    def _on_residue_clicked(self, item: QListWidgetItem):
        """Handle residue selection via click - updates plot display"""
        residue = item.data(Qt.UserRole)
        if residue:
            self._display_residue(residue)

    def _on_residue_row_changed(self, row: int):
        """Handle residue selection via keyboard navigation (arrow keys)"""
        if row < 0:
            return

        item = self.residue_list.item(row)
        if item is None:
            return

        residue = item.data(Qt.UserRole)
        if residue:
            self._display_residue(residue)

    def _display_residue(self, residue: str):
        """Display the selected residue's fit data"""
        if not self.t1_checkbox.isChecked() and not self.t2_checkbox.isChecked():
            return

        self._generate_plots([residue])

    def _update_max_residues(self):
        """Update measurement types and refresh current plot"""
        self.show_t1 = self.t1_checkbox.isChecked()
        self.show_t2 = self.t2_checkbox.isChecked()

        # Refresh current selection
        current_item = self.residue_list.currentItem()
        if current_item:
            residue = current_item.data(Qt.UserRole)
            if residue:
                self._display_residue(residue)

    def _on_field_change(self, button):
        """Handle field selection change and refresh plots"""
        btn_id = self.field_button_group.id(button)
        if btn_id == 1:
            self.field_mode = "field1"
        elif btn_id == 2:
            self.field_mode = "field2"
        else:
            self.field_mode = "overlay"

        # Auto-update current plot
        current_item = self.residue_list.currentItem()
        if current_item and (self.t1_checkbox.isChecked() or self.t2_checkbox.isChecked()):
            residue = current_item.data(Qt.UserRole)
            if residue:
                self._display_residue(residue)

    def _clear_selection(self):
        """Clear residue selection"""
        self.residue_list.clearSelection()
        self._show_blank_state()

    def _update_plots(self):
        """Update plots with current selection (legacy method kept for compatibility)"""
        current_item = self.residue_list.currentItem()
        if current_item:
            residue = current_item.data(Qt.UserRole)
            if residue:
                self._display_residue(residue)

    def _generate_plots(self, residues):
        """Generate plots for selected residues"""
        self.figure.clear()
        self.axes = []

        # Determine plot count
        plots = []
        for residue in sorted(residues, key=self._residue_sort_key):
            if self.t1_checkbox.isChecked():
                plots.append((residue, "T1"))
            if self.t2_checkbox.isChecked():
                plots.append((residue, "T2"))

        n_plots = len(plots)
        if n_plots == 0:
            self._show_blank_state()
            return

        # Create subplots (4x1 layout)
        for idx, (residue, meas_type) in enumerate(plots):
            ax = self.figure.add_subplot(n_plots, 1, idx + 1)
            self.axes.append(ax)

            # Plot data
            self._plot_single_fit(ax, residue, meas_type)

        self.figure.tight_layout()
        self.canvas.draw()

    def _plot_single_fit(self, ax, residue, meas_type):
        """Plot single fit on given axes"""
        if self.field_mode == "overlay":
            self._plot_overlay(ax, residue, meas_type)
        else:
            self._plot_single_field(ax, residue, meas_type, self.field_mode)

    def _plot_single_field(self, ax, residue, meas_type, field):
        """Plot single field data"""
        # Get data key
        data_key = f"{field}_{meas_type}"

        if data_key not in self.data:
            ax.text(0.5, 0.5, f"No {meas_type} data for {field}",
                   ha='center', va='center', transform=ax.transAxes)
            ax.axis('off')
            return

        # Find residue in data
        fit_data = None
        for fit in self.data[data_key]['fits']:
            if fit['residue'] == residue:
                fit_data = fit
                break

        if not fit_data:
            ax.text(0.5, 0.5, f"No data for residue {residue}",
                   ha='center', va='center', transform=ax.transAxes)
            ax.axis('off')
            return

        # Extract data
        metadata = self.data[data_key]['metadata']
        time_points = metadata['time_points']
        intensities = fit_data['intensities']
        fit_time = fit_data['fit_curve']['time']
        fit_intensity = fit_data['fit_curve']['intensity']
        t_value = fit_data['t2']
        t_error = fit_data['t2_err']

        # Plot
        ax.plot(time_points, intensities, 'bo', markersize=8, label='Data')
        ax.plot(fit_time, fit_intensity, 'b-', linewidth=2, label='Fit')

        # Annotation
        field_freq = metadata['field_freq']
        time_units = metadata['time_units']
        textstr = f"{meas_type} ({field_freq} MHz) = {t_value:.2f} ± {t_error:.2f} {time_units}"
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Labels
        ax.set_xlabel(f"Time ({time_units})")
        ax.set_ylabel("Signal Intensity")
        ax.set_title(f"Residue {residue} - {meas_type}")
        ax.legend()
        ax.grid(True, alpha=0.3)

    def _plot_overlay(self, ax, residue, meas_type):
        """Plot overlay of both fields"""
        # Get data for both fields
        field1_key = f"field1_{meas_type}"
        field2_key = f"field2_{meas_type}"

        has_field1 = field1_key in self.data
        has_field2 = field2_key in self.data

        if not has_field1 and not has_field2:
            ax.text(0.5, 0.5, f"No {meas_type} data available",
                   ha='center', va='center', transform=ax.transAxes)
            ax.axis('off')
            return

        # Find residue data
        fit1 = None
        fit2 = None

        if has_field1:
            for fit in self.data[field1_key]['fits']:
                if fit['residue'] == residue:
                    fit1 = fit
                    break

        if has_field2:
            for fit in self.data[field2_key]['fits']:
                if fit['residue'] == residue:
                    fit2 = fit
                    break

        if not fit1 and not fit2:
            ax.text(0.5, 0.5, f"No data for residue {residue}",
                   ha='center', va='center', transform=ax.transAxes)
            ax.axis('off')
            return

        # Plot Field 1 (blue)
        annotations = []
        time_units = "ms"

        if fit1:
            metadata1 = self.data[field1_key]['metadata']
            time_points1 = metadata1['time_points']
            intensities1 = fit1['intensities']
            fit_time1 = fit1['fit_curve']['time']
            fit_intensity1 = fit1['fit_curve']['intensity']

            ax.plot(time_points1, intensities1, 'bo', markersize=8, label=f"{metadata1['field_freq']} MHz data")
            ax.plot(fit_time1, fit_intensity1, 'b-', linewidth=2, label=f"{metadata1['field_freq']} MHz fit")

            annotations.append(f"{meas_type} ({metadata1['field_freq']} MHz) = {fit1['t2']:.2f} ± {fit1['t2_err']:.2f} {metadata1['time_units']}")
            time_units = metadata1['time_units']

        # Plot Field 2 (red)
        if fit2:
            metadata2 = self.data[field2_key]['metadata']
            time_points2 = metadata2['time_points']
            intensities2 = fit2['intensities']
            fit_time2 = fit2['fit_curve']['time']
            fit_intensity2 = fit2['fit_curve']['intensity']

            ax.plot(time_points2, intensities2, 'ro', markersize=8, label=f"{metadata2['field_freq']} MHz data")
            ax.plot(fit_time2, fit_intensity2, 'r-', linewidth=2, label=f"{metadata2['field_freq']} MHz fit")

            annotations.append(f"{meas_type} ({metadata2['field_freq']} MHz) = {fit2['t2']:.2f} ± {fit2['t2_err']:.2f} {metadata2['time_units']}")
            time_units = metadata2['time_units']

        # Annotation
        textstr = '\n'.join(annotations)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Labels
        ax.set_xlabel(f"Time ({time_units})")
        ax.set_ylabel("Signal Intensity")
        ax.set_title(f"Residue {residue} - {meas_type}")
        ax.legend()
        ax.grid(True, alpha=0.3)

    def _show_blank_state(self):
        """Show blank plot with instruction message"""
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.text(0.5, 0.5,
               "Select residues to display fits",
               ha='center', va='center',
               fontsize=16, color=SECONDARY_TEXT,
               transform=ax.transAxes)
        ax.axis('off')
        self.canvas.draw()


def open_fit_viewer(parent=None, json_folder=None, json_folders=None):
    """
    Convenience function to open the fit viewer

    Parameters:
    -----------
    parent : QWidget, optional
        Parent widget
    json_folder : str, optional
        Path to JSON data folder (legacy single folder)
    json_folders : list, optional
        List of JSON data folder paths (for multi-experiment support)

    Returns:
    --------
    FitViewer : The viewer window instance
    """
    viewer = FitViewer(parent, json_folder=json_folder, json_folders=json_folders)
    viewer.show()
    return viewer
