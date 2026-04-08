# ABOUTME: GUI dialog for Series QC analysis and comparison.
# ABOUTME: Provides drag-drop interface and visual QC reports for series integration results.

"""
Series QC Dialog for analyzing and comparing series integration results.

Provides:
- Drag-drop or browse interface for loading series results
- Summary dashboard with key QC metrics
- Flagged peaks table for problematic peaks
- Comparison table for multi-series analysis
"""

import os
import json
from pathlib import Path
from typing import Optional, List, Dict

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QListWidget, QListWidgetItem,
    QTableWidget, QTableWidgetItem, QHeaderView,
    QProgressBar, QTabWidget, QWidget, QGroupBox,
    QFileDialog, QMessageBox, QAbstractItemView,
    QFrame, QSplitter, QSizePolicy, QComboBox
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont, QColor, QBrush

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.utils.series_qc_analyzer import (
    SeriesQCAnalyzer,
    SeriesQCReport,
    SeriesComparator,
    SeriesType
)


class SeriesListWidget(QListWidget):
    """List widget for displaying loaded series with remove buttons."""

    series_removed = Signal(str)  # Emits series name when removed

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.setMinimumHeight(100)
        # Ensure proper focus and mouse handling on Linux
        self.setFocusPolicy(Qt.StrongFocus)
        self.setMouseTracking(True)
        self.setStyleSheet("""
            QListWidget {
                border: 2px dashed #999;
                border-radius: 8px;
                padding: 8px;
            }
            QListWidget::item {
                padding: 8px;
                border-radius: 4px;
                margin: 2px 0px;
            }
            QListWidget::item:selected {
                background-color: #cde;
            }
            QListWidget::item:hover {
                background-color: #eee;
            }
        """)


class SeriesQCDialog(BaseDialog):
    """Main QC dialog with browse, analysis, and comparison functionality."""

    def __init__(self, parent=None, main_window=None):
        super().__init__(
            parent=parent,
            title="Series QC Analysis",
            default_size=(950, 750),
            min_size=(800, 600),
            modal=False
        )

        self.main_window = main_window
        self.series_type = SeriesType.REPLICATE
        self.comparator = SeriesComparator(series_type=self.series_type)
        self.series_paths: Dict[str, str] = {}  # name -> path mapping

        self._setup_ui()

    def _setup_ui(self):
        """Build the main UI layout."""
        layout = QVBoxLayout(self)
        layout.setSpacing(12)

        # Section 0: Series type selector
        self._create_series_type_selector(layout)

        # Section 1: Series selection area
        self._create_series_section(layout)

        # Section 2: Action buttons
        self._create_action_buttons(layout)

        # Section 3: Results tabs
        self._create_results_tabs(layout)

    def _create_series_type_selector(self, parent_layout):
        """Create the series type selector row."""
        row_layout = QHBoxLayout()

        label = QLabel("Series Type:")
        label.setStyleSheet("font-weight: bold;")
        row_layout.addWidget(label)

        self.series_type_combo = QComboBox()
        self.series_type_combo.addItem("Replicate (volumes stable)", SeriesType.REPLICATE)
        self.series_type_combo.addItem("Relaxation (T1/T2 decay)", SeriesType.RELAXATION)
        self.series_type_combo.addItem("Titration (intensity changes)", SeriesType.TITRATION)
        self.series_type_combo.addItem("Kinetic (volume changes)", SeriesType.KINETIC)
        self.series_type_combo.setMinimumWidth(200)
        self.series_type_combo.currentIndexChanged.connect(self._on_series_type_changed)
        row_layout.addWidget(self.series_type_combo)

        # Help text
        help_label = QLabel(
            "Replicate: flags volume/height CV. Others: skip intensity flags, keep LW/position/R² flags."
        )
        help_label.setStyleSheet("color: #666; font-style: italic;")
        row_layout.addWidget(help_label)

        row_layout.addStretch()
        parent_layout.addLayout(row_layout)

    def _on_series_type_changed(self, index: int):
        """Handle series type combo box change."""
        self.series_type = self.series_type_combo.currentData()

        # Clear loaded series and recreate comparator with new type
        self.series_paths.clear()
        self.series_list.clear()
        self.comparator = SeriesComparator(series_type=self.series_type)
        self._update_button_states()

        # Update count label
        self.series_count_label.setText("0 series loaded")

        # Reset summary display to placeholder
        self._reset_summary_display()

    def _reset_summary_display(self):
        """Reset summary tab to placeholder state."""
        self.summary_label.setVisible(True)
        self.summary_label.setText(
            "Load series results and click 'Analyze Single' to see summary."
        )
        self.summary_content.setVisible(False)

        # Reset flagged tab
        self.flagged_table.setRowCount(0)
        self.results_tabs.setTabText(1, "Flagged Peaks")

        # Reset comparison tab
        self.comparison_label.setVisible(True)
        self.comparison_table.setVisible(False)
        self.comparison_table.setRowCount(0)
        self.comparison_summary.setVisible(False)

    def _create_series_section(self, parent_layout):
        """Create the series selection area with available and selected lists."""
        # Main horizontal layout with two panels
        h_layout = QHBoxLayout()

        # Left panel: Available Series
        left_group = QGroupBox("Available Series")
        left_layout = QVBoxLayout(left_group)

        # Available series list
        self.available_list = QListWidget()
        self.available_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.available_list.setMinimumHeight(100)
        self.available_list.setFocusPolicy(Qt.StrongFocus)
        self.available_list.setMouseTracking(True)
        # Allow both single-click selection and double-click to add
        self.available_list.itemClicked.connect(self._on_available_item_clicked)
        self.available_list.itemDoubleClicked.connect(self._on_add_selected)
        left_layout.addWidget(self.available_list)

        # Refresh button
        refresh_btn = QPushButton("Refresh")
        refresh_btn.clicked.connect(self._refresh_available_series)
        left_layout.addWidget(refresh_btn)

        h_layout.addWidget(left_group)

        # Center: Transfer buttons
        center_layout = QVBoxLayout()
        center_layout.addStretch()

        add_btn = QPushButton("→")
        add_btn.setFixedWidth(40)
        add_btn.setToolTip("Add selected series to analysis")
        add_btn.clicked.connect(self._on_add_selected)
        center_layout.addWidget(add_btn)

        remove_btn = QPushButton("←")
        remove_btn.setFixedWidth(40)
        remove_btn.setToolTip("Remove from analysis")
        remove_btn.clicked.connect(self._on_remove_selected)
        center_layout.addWidget(remove_btn)

        center_layout.addStretch()
        h_layout.addLayout(center_layout)

        # Right panel: Selected for Analysis
        right_group = QGroupBox("Selected for Analysis")
        right_layout = QVBoxLayout(right_group)

        # Selected series list
        self.series_list = SeriesListWidget()
        self.series_list.setMinimumHeight(100)
        right_layout.addWidget(self.series_list)

        # Browse and count
        btn_layout = QHBoxLayout()

        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.clicked.connect(self._on_browse)
        btn_layout.addWidget(self.browse_btn)

        btn_layout.addStretch()

        self.series_count_label = QLabel("0 series loaded")
        self.series_count_label.setStyleSheet("color: #666;")
        btn_layout.addWidget(self.series_count_label)

        right_layout.addLayout(btn_layout)
        h_layout.addWidget(right_group)

        parent_layout.addLayout(h_layout)

        # Connect selection change
        self.series_list.itemSelectionChanged.connect(self._on_selection_changed)

        # Load available series on startup
        self._refresh_available_series()

    def _create_action_buttons(self, parent_layout):
        """Create the Analyze and Compare action buttons."""
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()

        self.analyze_btn = QPushButton("Analyze Single")
        self.analyze_btn.setEnabled(False)
        self.analyze_btn.clicked.connect(self._on_analyze_single)
        self.analyze_btn.setMinimumWidth(120)
        btn_layout.addWidget(self.analyze_btn)

        self.compare_btn = QPushButton("Compare Selected")
        self.compare_btn.setEnabled(False)
        self.compare_btn.clicked.connect(self._on_compare)
        self.compare_btn.setMinimumWidth(120)
        btn_layout.addWidget(self.compare_btn)

        btn_layout.addStretch()
        parent_layout.addLayout(btn_layout)

    def _create_results_tabs(self, parent_layout):
        """Create the tabbed results area."""
        self.results_tabs = QTabWidget()

        # Tab 1: Summary Dashboard
        self.summary_widget = self._create_summary_tab()
        self.results_tabs.addTab(self.summary_widget, "Summary")

        # Tab 2: Flagged Peaks
        self.flagged_widget = self._create_flagged_tab()
        self.results_tabs.addTab(self.flagged_widget, "Flagged Peaks")

        # Tab 3: Comparison
        self.comparison_widget = self._create_comparison_tab()
        self.results_tabs.addTab(self.comparison_widget, "Comparison")

        parent_layout.addWidget(self.results_tabs)

    def _create_summary_tab(self) -> QWidget:
        """Create the summary dashboard tab."""
        widget = QWidget()
        layout = QVBoxLayout(widget)

        # Placeholder - will be populated after analysis
        self.summary_label = QLabel("Load series results and click 'Analyze Single' to see summary.")
        self.summary_label.setAlignment(Qt.AlignCenter)
        self.summary_label.setStyleSheet("color: #666; font-size: 14px;")
        layout.addWidget(self.summary_label)

        # Summary content (hidden until analysis)
        self.summary_content = QWidget()
        self.summary_content.setVisible(False)
        summary_layout = QGridLayout(self.summary_content)
        summary_layout.setColumnStretch(1, 1)
        summary_layout.setColumnStretch(3, 2)

        row = 0

        # Row 0: Basic info
        summary_layout.addWidget(QLabel("Series:"), row, 0)
        self.info_series_name = QLabel("-")
        self.info_series_name.setStyleSheet("font-weight: bold;")
        summary_layout.addWidget(self.info_series_name, row, 1)

        summary_layout.addWidget(QLabel("Peaks:"), row, 2)
        self.info_n_peaks = QLabel("-")
        summary_layout.addWidget(self.info_n_peaks, row, 3)

        summary_layout.addWidget(QLabel("Spectra:"), row, 4)
        self.info_n_spectra = QLabel("-")
        summary_layout.addWidget(self.info_n_spectra, row, 5)
        row += 1

        # Row 1: Series Type
        summary_layout.addWidget(QLabel("Analysis Mode:"), row, 0)
        self.info_series_type = QLabel("-")
        self.info_series_type.setStyleSheet("font-weight: bold; color: #444;")
        summary_layout.addWidget(self.info_series_type, row, 1, 1, 5)
        row += 1

        # Row 2: Volume CV (only for Replicate)
        self.label_vol_cv = QLabel("Volume CV:")
        summary_layout.addWidget(self.label_vol_cv, row, 0)
        self.info_vol_cv = QLabel("-")
        summary_layout.addWidget(self.info_vol_cv, row, 1)
        self.bar_vol_cv = QProgressBar()
        self.bar_vol_cv.setMaximum(100)
        self.bar_vol_cv.setTextVisible(False)
        summary_layout.addWidget(self.bar_vol_cv, row, 2, 1, 2)
        self.status_vol_cv = QLabel("-")
        summary_layout.addWidget(self.status_vol_cv, row, 4, 1, 2)
        # Store widgets for show/hide
        self.vol_cv_widgets = [self.label_vol_cv, self.info_vol_cv, self.bar_vol_cv, self.status_vol_cv]
        row += 1

        # Row 3: LW Stability F1 (15N/13C)
        summary_layout.addWidget(QLabel("LW Stability F1:"), row, 0)
        self.info_lw_f1 = QLabel("-")
        summary_layout.addWidget(self.info_lw_f1, row, 1)
        self.bar_lw_f1 = QProgressBar()
        self.bar_lw_f1.setMaximum(100)
        self.bar_lw_f1.setTextVisible(False)
        summary_layout.addWidget(self.bar_lw_f1, row, 2, 1, 2)
        self.status_lw_f1 = QLabel("-")
        summary_layout.addWidget(self.status_lw_f1, row, 4, 1, 2)
        row += 1

        # Row 4: LW Stability F2 (1H)
        summary_layout.addWidget(QLabel("LW Stability F2:"), row, 0)
        self.info_lw_f2 = QLabel("-")
        summary_layout.addWidget(self.info_lw_f2, row, 1)
        self.bar_lw_f2 = QProgressBar()
        self.bar_lw_f2.setMaximum(100)
        self.bar_lw_f2.setTextVisible(False)
        summary_layout.addWidget(self.bar_lw_f2, row, 2, 1, 2)
        self.status_lw_f2 = QLabel("-")
        summary_layout.addWidget(self.status_lw_f2, row, 4, 1, 2)
        row += 1

        # Row 5: Position Drift F1
        summary_layout.addWidget(QLabel("Position Drift F1:"), row, 0)
        self.info_pos_f1 = QLabel("-")
        summary_layout.addWidget(self.info_pos_f1, row, 1)
        self.bar_pos_f1 = QProgressBar()
        self.bar_pos_f1.setMaximum(100)
        self.bar_pos_f1.setTextVisible(False)
        summary_layout.addWidget(self.bar_pos_f1, row, 2, 1, 2)
        self.status_pos_f1 = QLabel("-")
        summary_layout.addWidget(self.status_pos_f1, row, 4, 1, 2)
        row += 1

        # Row 6: Position Drift F2
        summary_layout.addWidget(QLabel("Position Drift F2:"), row, 0)
        self.info_pos_f2 = QLabel("-")
        summary_layout.addWidget(self.info_pos_f2, row, 1)
        self.bar_pos_f2 = QProgressBar()
        self.bar_pos_f2.setMaximum(100)
        self.bar_pos_f2.setTextVisible(False)
        summary_layout.addWidget(self.bar_pos_f2, row, 2, 1, 2)
        self.status_pos_f2 = QLabel("-")
        summary_layout.addWidget(self.status_pos_f2, row, 4, 1, 2)
        row += 1

        # Row 7: R² Quality
        summary_layout.addWidget(QLabel("R² Quality:"), row, 0)
        self.info_r2 = QLabel("-")
        summary_layout.addWidget(self.info_r2, row, 1)
        self.bar_r2 = QProgressBar()
        self.bar_r2.setMaximum(100)
        self.bar_r2.setTextVisible(False)
        summary_layout.addWidget(self.bar_r2, row, 2, 1, 2)
        self.status_r2 = QLabel("-")
        summary_layout.addWidget(self.status_r2, row, 4, 1, 2)
        row += 1

        # Row 8: Cluster Stable
        summary_layout.addWidget(QLabel("Cluster Stable:"), row, 0)
        self.info_cluster = QLabel("-")
        summary_layout.addWidget(self.info_cluster, row, 1)
        self.bar_cluster = QProgressBar()
        self.bar_cluster.setMaximum(100)
        self.bar_cluster.setTextVisible(False)
        summary_layout.addWidget(self.bar_cluster, row, 2, 1, 2)
        self.status_cluster = QLabel("-")
        summary_layout.addWidget(self.status_cluster, row, 4, 1, 2)
        row += 1

        # Row 9: Flagged peaks
        summary_layout.addWidget(QLabel("Flagged Peaks:"), row, 0)
        self.info_flagged = QLabel("-")
        self.info_flagged.setStyleSheet("color: #b86800; font-weight: bold;")
        summary_layout.addWidget(self.info_flagged, row, 1, 1, 5)

        layout.addWidget(self.summary_content)
        layout.addStretch()

        return widget

    def _create_flagged_tab(self) -> QWidget:
        """Create the flagged peaks table tab."""
        widget = QWidget()
        layout = QVBoxLayout(widget)

        # Table for flagged peaks
        self.flagged_table = QTableWidget()
        self.flagged_table.setColumnCount(6)
        self.flagged_table.setHorizontalHeaderLabels([
            "Peak ID", "Assignment", "Issue", "Vol CV", "R² mean", "Details"
        ])

        # Configure table
        header = self.flagged_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.Stretch)
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(4, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(5, QHeaderView.Stretch)

        self.flagged_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.flagged_table.setAlternatingRowColors(True)
        self.flagged_table.setSortingEnabled(True)

        layout.addWidget(self.flagged_table)

        return widget

    def _create_comparison_tab(self) -> QWidget:
        """Create the comparison table tab."""
        widget = QWidget()
        layout = QVBoxLayout(widget)

        # Placeholder label
        self.comparison_label = QLabel(
            "Add 2 or more series and click 'Compare Selected' to see comparison."
        )
        self.comparison_label.setAlignment(Qt.AlignCenter)
        self.comparison_label.setStyleSheet("color: #666; font-size: 14px;")
        layout.addWidget(self.comparison_label)

        # Comparison table (hidden until comparison)
        self.comparison_table = QTableWidget()
        self.comparison_table.setVisible(False)
        self.comparison_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.comparison_table.setAlternatingRowColors(True)
        layout.addWidget(self.comparison_table)

        # Summary info (hidden until comparison)
        self.comparison_summary = QLabel()
        self.comparison_summary.setVisible(False)
        self.comparison_summary.setStyleSheet("color: #666; margin-top: 8px;")
        layout.addWidget(self.comparison_summary)

        layout.addStretch()

        return widget

    # =========================================================================
    # Event handlers
    # =========================================================================

    def _on_browse(self):
        """Handle Browse button click - open folder selection dialog."""
        # Use project folder if available, otherwise home
        initial_dir = str(Path.home())
        if self.main_window and hasattr(self.main_window, 'get_preferred_data_directory'):
            initial_dir = self.main_window.get_preferred_data_directory()

        folder = QFileDialog.getExistingDirectory(
            self,
            "Select Series Results Folder",
            initial_dir,
            QFileDialog.ShowDirsOnly
        )

        if not folder:
            return

        # Search for training_data.json
        training_data_path = self._find_training_data(folder)

        if not training_data_path:
            QMessageBox.warning(
                self,
                "Not Found",
                f"Could not find training_data.json in:\n{folder}\n\n"
                "Checked: folder itself, data_collection/ subfolder, and all subfolders."
            )
            return

        # Generate series name from folder
        series_name = Path(folder).name

        # Check for duplicate
        if series_name in self.series_paths:
            QMessageBox.warning(
                self,
                "Duplicate",
                f"Series '{series_name}' is already loaded."
            )
            return

        # Validate by loading
        if not self.comparator.add_series(series_name, training_data_path):
            QMessageBox.warning(
                self,
                "Invalid Data",
                f"Could not load training data from:\n{training_data_path}\n\n"
                "The file may be corrupted or in an unsupported format."
            )
            return

        # Add to list
        self.series_paths[series_name] = training_data_path
        self._add_series_to_list(series_name, training_data_path)
        self._update_button_states()

    def _find_training_data(self, folder: str) -> Optional[str]:
        """Search for training_data.json in folder and subfolders."""
        folder_path = Path(folder)

        # Check folder itself
        direct = folder_path / "training_data.json"
        if direct.exists():
            return str(direct)

        # Check data_collection subfolder
        data_coll = folder_path / "data_collection" / "training_data.json"
        if data_coll.exists():
            return str(data_coll)

        # Check any subfolder
        for subdir in folder_path.iterdir():
            if subdir.is_dir():
                candidate = subdir / "training_data.json"
                if candidate.exists():
                    return str(candidate)

        return None

    def _refresh_available_series(self):
        """Discover and display available series from main_window.saved_series."""
        self.available_list.clear()
        self._available_series_paths = {}  # name -> path mapping

        # Get in-session series from main_window.saved_series
        if self.main_window and hasattr(self.main_window, 'saved_series'):
            saved_series = getattr(self.main_window, 'saved_series', {}) or {}
            for series_name, batch_results in saved_series.items():
                training_path = None

                if hasattr(batch_results, 'metadata'):
                    metadata = batch_results.metadata or {}

                    # Look for training_data.json in output_folder/data_collection/
                    output_folder = metadata.get('output_folder')
                    if output_folder:
                        candidate = Path(output_folder) / "data_collection" / "training_data.json"
                        if candidate.exists():
                            training_path = str(candidate)

                if training_path:
                    self._available_series_paths[series_name] = training_path
                    # Show spectrum count if available
                    n_spectra = metadata.get('total_spectra', metadata.get('spectrum_count', '?'))
                    item = QListWidgetItem(f"● {series_name}  ({n_spectra} spectra)")
                    item.setData(Qt.UserRole, series_name)
                    item.setData(Qt.UserRole + 1, training_path)
                    item.setToolTip(f"Series: {series_name}\nPath: {training_path}")
                    # Explicitly set item flags for Linux compatibility
                    item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
                    self.available_list.addItem(item)
                else:
                    # Series exists but no training_data.json found
                    item = QListWidgetItem(f"○ {series_name}  (no QC data)")
                    item.setData(Qt.UserRole, series_name)
                    item.setData(Qt.UserRole + 1, None)
                    item.setToolTip(f"Series: {series_name}\nNo training_data.json found")
                    item.setFlags(item.flags() & ~Qt.ItemIsEnabled)  # Disable item
                    self.available_list.addItem(item)

    def _on_available_item_clicked(self, item):
        """Handle single-click on available series item (for selection feedback).

        This ensures the item is properly selected even on Linux systems
        where selection behavior can be inconsistent.
        """
        # Item is already selected by Qt, but we log for debugging
        import logging
        logger = logging.getLogger(__name__)
        series_name = item.data(Qt.UserRole)
        logger.debug(f"Available series item clicked: {series_name}")

    def _on_add_selected(self):
        """Add selected series from available list to analysis list."""
        selected = self.available_list.selectedItems()
        if not selected:
            return

        for item in selected:
            series_name = item.data(Qt.UserRole)
            training_path = item.data(Qt.UserRole + 1)

            # Skip if already in analysis list
            if series_name in self.series_paths:
                continue

            # Validate and add to comparator
            if self.comparator.add_series(series_name, training_path):
                self.series_paths[series_name] = training_path
                self._add_series_to_list(series_name, training_path)

        self._update_button_states()

    def _add_series_to_list(self, name: str, path: str):
        """Add a series entry to the list widget."""
        report = self.comparator.series.get(name)
        if report:
            n_spectra = report.n_spectra
            n_peaks = report.n_peaks
            text = f"{name}  ({n_peaks} peaks, {n_spectra} spectra)"
        else:
            text = name

        item = QListWidgetItem(text)
        item.setData(Qt.UserRole, name)
        item.setData(Qt.UserRole + 1, path)
        # Explicitly set item flags for Linux compatibility
        item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
        self.series_list.addItem(item)

        # Update count
        count = self.series_list.count()
        self.series_count_label.setText(f"{count} series loaded")

    def _on_remove_selected(self):
        """Remove selected series from list."""
        selected = self.series_list.selectedItems()
        for item in selected:
            name = item.data(Qt.UserRole)
            if name in self.series_paths:
                del self.series_paths[name]
            if name in self.comparator.series:
                del self.comparator.series[name]
            self.series_list.takeItem(self.series_list.row(item))

        count = self.series_list.count()
        self.series_count_label.setText(f"{count} series loaded")
        self._update_button_states()

    def _on_selection_changed(self):
        """Handle selection change in series list."""
        self._update_button_states()

    def _update_button_states(self):
        """Update enabled state of action buttons."""
        count = self.series_list.count()
        selected = len(self.series_list.selectedItems())

        self.analyze_btn.setEnabled(count >= 1)
        self.compare_btn.setEnabled(count >= 2)

    def _on_analyze_single(self):
        """Analyze the first loaded series and display results."""
        if not self.comparator.series:
            return

        # Get first series (or selected if one is selected)
        selected = self.series_list.selectedItems()
        if selected:
            name = selected[0].data(Qt.UserRole)
        else:
            name = list(self.comparator.series.keys())[0]

        report = self.comparator.series.get(name)
        if not report:
            return

        self._display_summary(name, report)
        self._display_flagged_peaks(report)

        # Save report
        output_dir = Path(self.series_paths[name]).parent
        analyzer = SeriesQCAnalyzer()
        analyzer.report = report
        save_path = analyzer.save_report(str(output_dir))

        QMessageBox.information(
            self,
            "Analysis Complete",
            f"QC report saved to:\n{save_path}"
        )

    def _display_summary(self, name: str, report: SeriesQCReport):
        """Display summary dashboard for a report."""
        self.summary_label.setVisible(False)
        self.summary_content.setVisible(True)

        # Basic info
        self.info_series_name.setText(name)
        self.info_n_peaks.setText(str(report.n_peaks))
        self.info_n_spectra.setText(str(report.n_spectra))

        # Series type display
        type_labels = {
            SeriesType.REPLICATE: "Replicate (volume stability checked)",
            SeriesType.RELAXATION: "Relaxation (volume changes expected)",
            SeriesType.TITRATION: "Titration (intensity changes expected)",
            SeriesType.KINETIC: "Kinetic (volume changes expected)",
        }
        self.info_series_type.setText(type_labels.get(report.series_type, str(report.series_type)))

        # Show/hide Volume CV row based on series type
        show_volume_cv = (report.series_type == SeriesType.REPLICATE)
        for w in self.vol_cv_widgets:
            w.setVisible(show_volume_cv)

        # Volume CV (lower is better, scale 0-0.5 to 0-100)
        if show_volume_cv:
            vol_cv = report.volume_cv_median
            self.info_vol_cv.setText(f"{vol_cv:.3f}")
            vol_pct = min(100, int((1 - vol_cv / 0.5) * 100))
            self.bar_vol_cv.setValue(max(0, vol_pct))
            self._set_status(self.status_vol_cv, vol_cv, thresholds=(0.15, 0.3), lower_better=True)

        # LW Stability F1 (lower CV is better)
        n_peaks = report.n_peaks or 1
        lw_f1_cv = report.lw_f1_cv_median
        self.info_lw_f1.setText(f"CV={lw_f1_cv:.3f} ({report.n_lw_f1_stable}/{n_peaks} stable)")
        lw_f1_pct = min(100, int((1 - lw_f1_cv / 0.5) * 100))
        self.bar_lw_f1.setValue(max(0, lw_f1_pct))
        self._set_status(self.status_lw_f1, lw_f1_cv, thresholds=(0.15, 0.3), lower_better=True)

        # LW Stability F2 (lower CV is better)
        lw_f2_cv = report.lw_f2_cv_median
        self.info_lw_f2.setText(f"CV={lw_f2_cv:.3f} ({report.n_lw_f2_stable}/{n_peaks} stable)")
        lw_f2_pct = min(100, int((1 - lw_f2_cv / 0.5) * 100))
        self.bar_lw_f2.setValue(max(0, lw_f2_pct))
        self._set_status(self.status_lw_f2, lw_f2_cv, thresholds=(0.15, 0.3), lower_better=True)

        # Position Drift F1 (lower is better, threshold 0.05 ppm)
        pos_f1 = report.pos_f1_drift_median
        self.info_pos_f1.setText(f"{pos_f1:.4f} ppm ({report.n_pos_f1_stable}/{n_peaks} stable)")
        pos_f1_pct = min(100, int((1 - pos_f1 / 0.1) * 100))
        self.bar_pos_f1.setValue(max(0, pos_f1_pct))
        self._set_status(self.status_pos_f1, pos_f1, thresholds=(0.03, 0.05), lower_better=True)

        # Position Drift F2 (lower is better, threshold 0.01 ppm)
        pos_f2 = report.pos_f2_drift_median
        self.info_pos_f2.setText(f"{pos_f2:.4f} ppm ({report.n_pos_f2_stable}/{n_peaks} stable)")
        pos_f2_pct = min(100, int((1 - pos_f2 / 0.02) * 100))
        self.bar_pos_f2.setValue(max(0, pos_f2_pct))
        self._set_status(self.status_pos_f2, pos_f2, thresholds=(0.005, 0.01), lower_better=True)

        # R² Quality (higher is better, scale 0.5-1.0 to 0-100)
        r2 = report.r2_mean
        self.info_r2.setText(f"{r2:.3f}")
        r2_pct = min(100, int((r2 - 0.5) / 0.5 * 100))
        self.bar_r2.setValue(max(0, r2_pct))
        self._set_status(self.status_r2, r2, thresholds=(0.7, 0.9), lower_better=False)

        # Cluster Stable
        stable_pct = int(report.n_cluster_stable / n_peaks * 100)
        self.info_cluster.setText(f"{report.n_cluster_stable}/{report.n_peaks}")
        self.bar_cluster.setValue(stable_pct)
        self._set_status(self.status_cluster, stable_pct / 100.0, thresholds=(0.7, 0.9), lower_better=False)

        # Flagged peaks
        n_flagged = len(report.flagged_peaks)
        self.info_flagged.setText(f"{n_flagged} peaks flagged for review")

        # Switch to summary tab
        self.results_tabs.setCurrentIndex(0)

    def _set_status(self, label: QLabel, value: float, thresholds: tuple, lower_better: bool):
        """Set status label color based on value and thresholds."""
        t1, t2 = thresholds
        if lower_better:
            if value <= t1:
                label.setText("Good")
                label.setStyleSheet("color: #2a8; font-weight: bold;")
            elif value <= t2:
                label.setText("OK")
                label.setStyleSheet("color: #b86800; font-weight: bold;")
            else:
                label.setText("Warning")
                label.setStyleSheet("color: #c44; font-weight: bold;")
        else:
            if value >= t2:
                label.setText("Good")
                label.setStyleSheet("color: #2a8; font-weight: bold;")
            elif value >= t1:
                label.setText("OK")
                label.setStyleSheet("color: #b86800; font-weight: bold;")
            else:
                label.setText("Poor")
                label.setStyleSheet("color: #c44; font-weight: bold;")

    def _display_flagged_peaks(self, report: SeriesQCReport):
        """Populate the flagged peaks table."""
        self.flagged_table.setRowCount(0)

        flagged_peaks = report.flagged_peaks
        self.flagged_table.setRowCount(len(flagged_peaks))

        for row, pid in enumerate(flagged_peaks):
            stats = report.peak_stats.get(pid)
            if not stats:
                continue

            # Peak ID
            self.flagged_table.setItem(row, 0, QTableWidgetItem(pid))

            # Assignment
            self.flagged_table.setItem(row, 1, QTableWidgetItem(stats.assignment))

            # Issue (primary flag)
            issue = stats.flags[0] if stats.flags else "-"
            self.flagged_table.setItem(row, 2, QTableWidgetItem(issue))

            # Volume CV
            vol_item = QTableWidgetItem(f"{stats.volume_cv:.3f}")
            if stats.volume_cv > 0.3:
                vol_item.setForeground(QBrush(QColor("#c44")))
            self.flagged_table.setItem(row, 3, vol_item)

            # R² mean
            r2_item = QTableWidgetItem(f"{stats.r2_mean:.3f}")
            if stats.r2_mean < 0.7:
                r2_item.setForeground(QBrush(QColor("#c44")))
            self.flagged_table.setItem(row, 4, r2_item)

            # Details (all flags)
            details = ", ".join(stats.flags)
            self.flagged_table.setItem(row, 5, QTableWidgetItem(details))

        # Update tab label
        self.results_tabs.setTabText(1, f"Flagged ({len(flagged_peaks)})")

    def _on_compare(self):
        """Compare all loaded series and display results."""
        if len(self.comparator.series) < 2:
            QMessageBox.warning(
                self,
                "Need More Series",
                "Please load at least 2 series to compare."
            )
            return

        result = self.comparator.compare()
        self._display_comparison(result)

        # Switch to comparison tab
        self.results_tabs.setCurrentIndex(2)

    def _display_comparison(self, result: Dict):
        """Display the comparison results."""
        self.comparison_label.setVisible(False)
        self.comparison_table.setVisible(True)
        self.comparison_summary.setVisible(True)

        series_names = list(self.comparator.series.keys())
        n_series = len(series_names)

        # Setup table: Metric | Series A | Series B | ... | Winner
        self.comparison_table.setColumnCount(n_series + 2)
        headers = ["Metric"] + series_names + ["Winner"]
        self.comparison_table.setHorizontalHeaderLabels(headers)

        # Metrics to display
        metrics = ["volume_cv_median", "r2_mean", "n_flagged", "n_cluster_stable"]
        metric_labels = {
            "volume_cv_median": "Volume CV (median)",
            "r2_mean": "R² (mean)",
            "n_flagged": "Flagged Peaks",
            "n_cluster_stable": "Cluster Stable",
        }

        self.comparison_table.setRowCount(len(metrics))

        for row, metric in enumerate(metrics):
            metric_data = result.get(metric, {})

            # Metric name
            self.comparison_table.setItem(row, 0, QTableWidgetItem(metric_labels.get(metric, metric)))

            # Series values
            values = metric_data.get("values", {})
            for col, name in enumerate(series_names, start=1):
                val = values.get(name, "-")
                if isinstance(val, float):
                    text = f"{val:.3f}"
                else:
                    text = str(val)
                self.comparison_table.setItem(row, col, QTableWidgetItem(text))

            # Winner
            winner = metric_data.get("winner", "-")
            winner_item = QTableWidgetItem(winner)
            winner_item.setForeground(QBrush(QColor("#2a8")))
            winner_item.setFont(QFont("", -1, QFont.Bold))
            self.comparison_table.setItem(row, n_series + 1, winner_item)

        # Configure table
        header = self.comparison_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        for i in range(1, n_series + 2):
            header.setSectionResizeMode(i, QHeaderView.ResizeToContents)

        # Summary text
        n_common = result.get("n_common_peaks", 0)
        n_unique = result.get("n_unique_peaks", {})
        unique_text = ", ".join([f"{name}: {count}" for name, count in n_unique.items()])
        self.comparison_summary.setText(
            f"Common peaks: {n_common}  |  Unique peaks: {unique_text}"
        )
