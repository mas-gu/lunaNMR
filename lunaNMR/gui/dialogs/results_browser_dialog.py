# ABOUTME: Results browser dialog for viewing batch/series integration results
# ABOUTME: Provides tabbed interface with overview stats, peak tracking, and evolution plots

import logging
from typing import Dict, List

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QTabWidget, QWidget, QLabel,
    QTableWidget, QTableWidgetItem, QHeaderView, QPushButton,
    QComboBox, QGroupBox, QSplitter, QMessageBox
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.gui.components.matplotlib_widget import MatplotlibWidget
from lunaNMR.gui.styles.design_system import (
    SPACING_SM, SPACING_MD, SPACING_LG,
    FONT_SIZE_BODY, FONT_SIZE_SECTION_LABEL, FONT_SIZE_SMALL,
    PRIMARY_TEXT, SECONDARY_TEXT, SUCCESS_GREEN, ERROR_RED, WARNING_ORANGE,
    INFO_BLUE,
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER, BUTTON_CORNER_RADIUS, BUTTON_HEIGHT_DIALOG,
    FRAME_BG_COLOR
)

logger = logging.getLogger(__name__)


class ResultsBrowserDialog(BaseDialog):
    """Dialog for browsing batch/series integration results.

    Provides multi-tab interface for:
        - Overview statistics and summary
        - Peak tracking across spectra
        - Evolution plots over series
        - Detailed per-spectrum analysis

    Based on v0.9 results browser (main_gui.py:6785-7080)

    Example:
        ```python
        dialog = ResultsBrowserDialog(parent, batch_results)
        dialog.exec()
        ```
    """

    def __init__(self, parent=None, batch_results=None, main_window=None):
        """Initialize the results browser dialog.

        Args:
            parent: Parent widget
            batch_results: BatchResults object with series integration results
            main_window: Reference to MainWindow for accessing spectra data
        """
        super().__init__(
            parent=parent,
            title="Results Browser - lunaNMR v1.0",
            default_size=(1000, 700),
            modal=False  # Non-modal so user can work while browsing
        )

        self.batch_results = batch_results
        self.main_window = main_window

        # Setup UI
        self.setup_ui()

        # Populate with data
        self.populate_data()

        # Center on parent
        if parent:
            self.center_on_parent()
        else:
            self.center_on_screen()

        logger.debug("ResultsBrowserDialog initialized")

    def setup_ui(self):
        """Setup the dialog user interface."""
        layout = QVBoxLayout()
        layout.setSpacing(SPACING_MD)
        layout.setContentsMargins(SPACING_LG, SPACING_LG, SPACING_LG, SPACING_LG)

        # Title
        title_label = QLabel("Series Integration Results")
        title_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SECTION_LABEL}pt;
                font-weight: bold;
                color: {PRIMARY_TEXT};
                padding-bottom: {SPACING_SM}px;
            }}
        """)
        layout.addWidget(title_label)

        # Tab widget
        self.tab_widget = QTabWidget()
        self.tab_widget.setStyleSheet(f"""
            QTabWidget::pane {{
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                background-color: {FRAME_BG_COLOR};
            }}
            QTabBar::tab {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_TEXT};
                padding: {SPACING_SM}px {SPACING_MD}px;
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-bottom: none;
                border-top-left-radius: {BUTTON_CORNER_RADIUS}px;
                border-top-right-radius: {BUTTON_CORNER_RADIUS}px;
                margin-right: 2px;
            }}
            QTabBar::tab:selected {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
            }}
        """)

        # Create tabs
        self.overview_tab = self.create_overview_tab()
        self.tracking_tab = self.create_tracking_tab()
        self.evolution_tab = self.create_evolution_tab()
        self.details_tab = self.create_details_tab()

        self.tab_widget.addTab(self.overview_tab, "Overview")
        self.tab_widget.addTab(self.tracking_tab, "Peak Tracking")
        self.tab_widget.addTab(self.evolution_tab, "Evolution Plots")
        self.tab_widget.addTab(self.details_tab, "Details")

        layout.addWidget(self.tab_widget, stretch=1)

        # Button row
        button_layout = QHBoxLayout()
        button_layout.addStretch()

        export_btn = QPushButton("Export Data...")
        export_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        export_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}pt;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
        """)
        export_btn.clicked.connect(self.export_results)
        button_layout.addWidget(export_btn)

        close_btn = QPushButton("Close")
        close_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        close_btn.setMinimumWidth(100)
        close_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}pt;
                font-weight: bold;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
        """)
        close_btn.clicked.connect(self.accept)
        button_layout.addWidget(close_btn)

        layout.addLayout(button_layout)

        self.setLayout(layout)

    def create_overview_tab(self) -> QWidget:
        """Create overview statistics tab."""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setSpacing(SPACING_MD)

        # Statistics group
        stats_group = QGroupBox("Summary Statistics")
        stats_group.setStyleSheet(self._group_style())
        stats_layout = QVBoxLayout(stats_group)

        # Statistics labels
        self.total_spectra_label = QLabel("Total Spectra: -")
        self.successful_spectra_label = QLabel("Successful: -")
        self.failed_spectra_label = QLabel("Failed: -")
        self.total_peaks_label = QLabel("Total Peaks Tracked: -")
        self.avg_detection_label = QLabel("Avg Detection Rate: -")
        self.avg_r_squared_label = QLabel("Avg R²: -")

        for lbl in [self.total_spectra_label, self.successful_spectra_label,
                    self.failed_spectra_label, self.total_peaks_label,
                    self.avg_detection_label, self.avg_r_squared_label]:
            lbl.setStyleSheet(f"font-size: {FONT_SIZE_BODY}pt; padding: 2px;")
            stats_layout.addWidget(lbl)

        layout.addWidget(stats_group)

        # Quality distribution chart
        quality_group = QGroupBox("Quality Distribution")
        quality_group.setStyleSheet(self._group_style())
        quality_layout = QVBoxLayout(quality_group)

        self.quality_chart = MatplotlibWidget(figsize=(6, 3))
        quality_layout.addWidget(self.quality_chart)

        layout.addWidget(quality_group)

        layout.addStretch()

        return widget

    def create_tracking_tab(self) -> QWidget:
        """Create peak tracking table tab."""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setSpacing(SPACING_MD)

        # Peak selector
        selector_layout = QHBoxLayout()
        selector_layout.addWidget(QLabel("Select Peak:"))

        self.peak_combo = QComboBox()
        self.peak_combo.setMinimumWidth(200)
        self.peak_combo.currentIndexChanged.connect(self.on_peak_selected)
        selector_layout.addWidget(self.peak_combo)

        selector_layout.addStretch()
        layout.addLayout(selector_layout)

        # Tracking table
        self.tracking_table = QTableWidget()
        self.tracking_table.setColumnCount(7)
        self.tracking_table.setHorizontalHeaderLabels([
            "Spectrum", "Pos F2 (ppm)", "Pos F1 (ppm)",
            "Intensity", "Volume", "R²", "Status"
        ])
        self.tracking_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tracking_table.setAlternatingRowColors(True)
        self.tracking_table.setStyleSheet(f"""
            QTableWidget {{
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                font-size: {FONT_SIZE_SMALL}pt;
            }}
            QTableWidget::item {{
                padding: 4px;
            }}
        """)
        layout.addWidget(self.tracking_table)

        return widget

    def create_evolution_tab(self) -> QWidget:
        """Create evolution plots tab."""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setSpacing(SPACING_MD)

        # Peak selector for evolution
        selector_layout = QHBoxLayout()
        selector_layout.addWidget(QLabel("Select Peak:"))

        self.evolution_peak_combo = QComboBox()
        self.evolution_peak_combo.setMinimumWidth(200)
        self.evolution_peak_combo.currentIndexChanged.connect(self.on_evolution_peak_selected)
        selector_layout.addWidget(self.evolution_peak_combo)

        # Property selector
        selector_layout.addWidget(QLabel("Property:"))
        self.property_combo = QComboBox()
        self.property_combo.addItems(["Intensity", "Volume", "R²", "Position F2", "Position F1"])
        self.property_combo.currentIndexChanged.connect(self.on_evolution_peak_selected)
        selector_layout.addWidget(self.property_combo)

        selector_layout.addStretch()
        layout.addLayout(selector_layout)

        # Evolution chart
        self.evolution_chart = MatplotlibWidget(figsize=(8, 5))
        layout.addWidget(self.evolution_chart)

        return widget

    def create_details_tab(self) -> QWidget:
        """Create detailed per-spectrum analysis tab."""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setSpacing(SPACING_MD)

        # Spectrum selector
        selector_layout = QHBoxLayout()
        selector_layout.addWidget(QLabel("Select Spectrum:"))

        self.spectrum_combo = QComboBox()
        self.spectrum_combo.setMinimumWidth(300)
        self.spectrum_combo.currentIndexChanged.connect(self.on_spectrum_selected)
        selector_layout.addWidget(self.spectrum_combo)

        selector_layout.addStretch()
        layout.addLayout(selector_layout)

        # Splitter for info and table
        splitter = QSplitter(Qt.Vertical)

        # Spectrum info
        info_widget = QWidget()
        info_layout = QVBoxLayout(info_widget)
        info_layout.setContentsMargins(0, 0, 0, 0)

        info_group = QGroupBox("Spectrum Information")
        info_group.setStyleSheet(self._group_style())
        info_inner = QVBoxLayout(info_group)

        self.spectrum_status_label = QLabel("Status: -")
        self.spectrum_peaks_label = QLabel("Peaks: -")
        self.spectrum_detection_label = QLabel("Detection Rate: -")
        self.spectrum_avg_r2_label = QLabel("Avg R²: -")

        for lbl in [self.spectrum_status_label, self.spectrum_peaks_label,
                    self.spectrum_detection_label, self.spectrum_avg_r2_label]:
            lbl.setStyleSheet(f"font-size: {FONT_SIZE_BODY}pt; padding: 2px;")
            info_inner.addWidget(lbl)

        info_layout.addWidget(info_group)
        splitter.addWidget(info_widget)

        # Peak details table
        self.details_table = QTableWidget()
        self.details_table.setColumnCount(8)
        self.details_table.setHorizontalHeaderLabels([
            "Assignment", "Pos F2", "Pos F1", "Intensity",
            "Volume", "LW F2", "LW F1", "R²"
        ])
        self.details_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.details_table.setAlternatingRowColors(True)
        self.details_table.setStyleSheet(f"""
            QTableWidget {{
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                font-size: {FONT_SIZE_SMALL}pt;
            }}
        """)
        splitter.addWidget(self.details_table)

        splitter.setSizes([150, 400])
        layout.addWidget(splitter)

        return widget

    def _group_style(self) -> str:
        """Return standard QGroupBox stylesheet."""
        return f"""
            QGroupBox {{
                font-size: {FONT_SIZE_BODY}pt;
                font-weight: bold;
                color: {PRIMARY_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                margin-top: {SPACING_SM}px;
                padding-top: {SPACING_MD}px;
                background-color: {FRAME_BG_COLOR};
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                left: {SPACING_SM}px;
                padding: 0 {SPACING_SM}px;
            }}
        """

    def populate_data(self):
        """Populate all tabs with batch results data."""
        if not self.batch_results or not hasattr(self.batch_results, 'results'):
            return

        results = self.batch_results.results

        # Populate overview tab
        self.populate_overview(results)

        # Populate peak combos
        self.populate_peak_combos(results)

        # Populate spectrum combo
        self.populate_spectrum_combo(results)

    def populate_overview(self, results: Dict):
        """Populate overview statistics."""
        total_spectra = len(results)
        successful = sum(1 for d in results.values() if d.get('status') == 'success')
        failed = total_spectra - successful

        # Collect stats
        all_r_squared = []
        all_detection_rates = []
        total_peaks = 0

        for spectrum_name, data in results.items():
            fitted = data.get('fitted_results', [])
            for peak in fitted:
                if isinstance(peak, dict) and 'r_squared' in peak:
                    all_r_squared.append(peak['r_squared'])
            total_peaks += data.get('total_peaks', 0)
            if 'detection_rate' in data:
                all_detection_rates.append(data['detection_rate'])

        avg_r2 = sum(all_r_squared) / len(all_r_squared) if all_r_squared else 0
        avg_detection = sum(all_detection_rates) / len(all_detection_rates) if all_detection_rates else 0

        # Update labels
        self.total_spectra_label.setText(f"Total Spectra: {total_spectra}")
        self.successful_spectra_label.setText(f"Successful: {successful}")
        self.failed_spectra_label.setText(f"Failed: {failed}")
        self.total_peaks_label.setText(f"Total Peak Fits: {len(all_r_squared)}")
        self.avg_detection_label.setText(f"Avg Detection Rate: {avg_detection:.1%}")
        self.avg_r_squared_label.setText(f"Avg R²: {avg_r2:.4f}")

        # Color code success/failure
        self.successful_spectra_label.setStyleSheet(
            f"font-size: {FONT_SIZE_BODY}pt; color: {SUCCESS_GREEN}; padding: 2px;")
        if failed > 0:
            self.failed_spectra_label.setStyleSheet(
                f"font-size: {FONT_SIZE_BODY}pt; color: {ERROR_RED}; padding: 2px;")

        # Draw quality distribution chart
        self.draw_quality_chart(all_r_squared)

    def draw_quality_chart(self, r_squared_values: List[float]):
        """Draw R² distribution histogram."""
        ax = self.quality_chart.figure.add_subplot(111)
        ax.clear()

        if r_squared_values:
            ax.hist(r_squared_values, bins=20, color=INFO_BLUE, alpha=0.7, edgecolor='white')
            ax.axvline(x=0.8, color=ERROR_RED, linestyle='--', label='R²=0.8 threshold')
            ax.set_xlabel('R² Value')
            ax.set_ylabel('Count')
            ax.set_title('Fitting Quality Distribution')
            ax.legend()
        else:
            ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)

        self.quality_chart.figure.tight_layout()
        self.quality_chart.canvas.draw()

    def populate_peak_combos(self, results: Dict):
        """Populate peak selection combo boxes."""
        # Get unique assignments across all spectra
        assignments = set()

        for spectrum_name, data in results.items():
            fitted = data.get('fitted_results', [])
            for peak in fitted:
                if isinstance(peak, dict):
                    assignment = peak.get('assignment') or peak.get('Assignment', '')
                    if assignment:
                        assignments.add(assignment)

        # Sort assignments
        sorted_assignments = sorted(assignments)

        # Update both combos
        self.peak_combo.clear()
        self.peak_combo.addItems(sorted_assignments)

        self.evolution_peak_combo.clear()
        self.evolution_peak_combo.addItems(sorted_assignments)

    def populate_spectrum_combo(self, results: Dict):
        """Populate spectrum selection combo."""
        self.spectrum_combo.clear()
        self.spectrum_combo.addItems(sorted(results.keys()))

    def on_peak_selected(self, index: int):
        """Handle peak selection in tracking tab."""
        if index < 0:
            return

        assignment = self.peak_combo.currentText()
        if not assignment:
            return

        results = self.batch_results.results

        # Clear table
        self.tracking_table.setRowCount(0)

        # Find this peak across all spectra
        for spectrum_name in sorted(results.keys()):
            data = results[spectrum_name]
            fitted = data.get('fitted_results', [])

            for peak in fitted:
                if isinstance(peak, dict):
                    peak_assignment = peak.get('assignment') or peak.get('Assignment', '')
                    if peak_assignment == assignment:
                        row = self.tracking_table.rowCount()
                        self.tracking_table.insertRow(row)

                        # Spectrum name
                        self.tracking_table.setItem(row, 0, QTableWidgetItem(spectrum_name))

                        # Position F2
                        pos_f2 = peak.get('pos_f2') or peak.get('Position_X', 0)
                        self.tracking_table.setItem(row, 1, QTableWidgetItem(f"{pos_f2:.4f}"))

                        # Position F1
                        pos_f1 = peak.get('pos_f1') or peak.get('Position_Y', 0)
                        self.tracking_table.setItem(row, 2, QTableWidgetItem(f"{pos_f1:.3f}"))

                        # Intensity
                        intensity = peak.get('intensity', 0)
                        self.tracking_table.setItem(row, 3, QTableWidgetItem(f"{intensity:.2e}"))

                        # Volume
                        volume = peak.get('volume', 0)
                        self.tracking_table.setItem(row, 4, QTableWidgetItem(f"{volume:.2e}"))

                        # R²
                        r2 = peak.get('r_squared', 0)
                        r2_item = QTableWidgetItem(f"{r2:.4f}")
                        if r2 < 0.7:
                            r2_item.setForeground(QColor(ERROR_RED))
                        elif r2 < 0.85:
                            r2_item.setForeground(QColor(WARNING_ORANGE))
                        else:
                            r2_item.setForeground(QColor(SUCCESS_GREEN))
                        self.tracking_table.setItem(row, 5, r2_item)

                        # Status
                        status = "Good" if r2 >= 0.8 else ("Fair" if r2 >= 0.7 else "Poor")
                        self.tracking_table.setItem(row, 6, QTableWidgetItem(status))

    def on_evolution_peak_selected(self, index: int = None):
        """Handle peak/property selection for evolution plot."""
        assignment = self.evolution_peak_combo.currentText()
        property_name = self.property_combo.currentText()

        if not assignment:
            return

        results = self.batch_results.results

        # Collect data series
        spectra_names = []
        values = []

        property_map = {
            "Intensity": "intensity",
            "Volume": "volume",
            "R²": "r_squared",
            "Position F2": "pos_f2",
            "Position F1": "pos_f1"
        }
        key = property_map.get(property_name, "intensity")

        for spectrum_name in sorted(results.keys()):
            data = results[spectrum_name]
            fitted = data.get('fitted_results', [])

            for peak in fitted:
                if isinstance(peak, dict):
                    peak_assignment = peak.get('assignment') or peak.get('Assignment', '')
                    if peak_assignment == assignment:
                        spectra_names.append(spectrum_name)
                        # Handle alternate key names
                        if key == 'pos_f2':
                            val = peak.get('pos_f2') or peak.get('Position_X', 0)
                        elif key == 'pos_f1':
                            val = peak.get('pos_f1') or peak.get('Position_Y', 0)
                        else:
                            val = peak.get(key, 0)
                        values.append(val)
                        break

        # Draw evolution plot
        self.draw_evolution_plot(spectra_names, values, property_name, assignment)

    def draw_evolution_plot(self, spectra: List[str], values: List[float],
                            property_name: str, assignment: str):
        """Draw evolution plot for selected peak."""
        ax = self.evolution_chart.figure.clear()
        ax = self.evolution_chart.figure.add_subplot(111)

        if values:
            x = range(len(values))
            ax.plot(x, values, 'o-', color=INFO_BLUE, linewidth=2, markersize=6)
            ax.set_xticks(x)
            ax.set_xticklabels([s[:15] + '...' if len(s) > 15 else s for s in spectra],
                              rotation=45, ha='right')
            ax.set_ylabel(property_name)
            ax.set_title(f'{property_name} Evolution: {assignment}')
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)

        self.evolution_chart.figure.tight_layout()
        self.evolution_chart.canvas.draw()

    def on_spectrum_selected(self, index: int):
        """Handle spectrum selection in details tab."""
        if index < 0:
            return

        spectrum_name = self.spectrum_combo.currentText()
        if not spectrum_name:
            return

        results = self.batch_results.results
        data = results.get(spectrum_name, {})

        # Update info labels
        status = data.get('status', 'unknown')
        status_color = SUCCESS_GREEN if status == 'success' else ERROR_RED
        self.spectrum_status_label.setText(f"Status: {status}")
        self.spectrum_status_label.setStyleSheet(
            f"font-size: {FONT_SIZE_BODY}pt; color: {status_color}; padding: 2px;")

        total_peaks = data.get('total_peaks', 0)
        detected = data.get('detected_peaks', 0)
        self.spectrum_peaks_label.setText(f"Peaks: {detected}/{total_peaks}")

        detection_rate = data.get('detection_rate', 0)
        self.spectrum_detection_label.setText(f"Detection Rate: {detection_rate:.1%}")

        # Calculate avg R²
        fitted = data.get('fitted_results', [])
        r2_values = [p.get('r_squared', 0) for p in fitted if isinstance(p, dict) and 'r_squared' in p]
        avg_r2 = sum(r2_values) / len(r2_values) if r2_values else 0
        self.spectrum_avg_r2_label.setText(f"Avg R²: {avg_r2:.4f}")

        # Populate details table
        self.details_table.setRowCount(0)

        for peak in fitted:
            if not isinstance(peak, dict):
                continue

            row = self.details_table.rowCount()
            self.details_table.insertRow(row)

            # Assignment
            assignment = peak.get('assignment') or peak.get('Assignment', '-')
            self.details_table.setItem(row, 0, QTableWidgetItem(str(assignment)))

            # Position F2
            pos_f2 = peak.get('pos_f2') or peak.get('Position_X', 0)
            self.details_table.setItem(row, 1, QTableWidgetItem(f"{pos_f2:.4f}"))

            # Position F1
            pos_f1 = peak.get('pos_f1') or peak.get('Position_Y', 0)
            self.details_table.setItem(row, 2, QTableWidgetItem(f"{pos_f1:.3f}"))

            # Intensity
            intensity = peak.get('intensity', 0)
            self.details_table.setItem(row, 3, QTableWidgetItem(f"{intensity:.2e}"))

            # Volume
            volume = peak.get('volume', 0)
            self.details_table.setItem(row, 4, QTableWidgetItem(f"{volume:.2e}"))

            # Linewidth F2
            lw_f2 = peak.get('lw_f2') or peak.get('linewidth_x', 0)
            self.details_table.setItem(row, 5, QTableWidgetItem(f"{lw_f2:.4f}"))

            # Linewidth F1
            lw_f1 = peak.get('lw_f1') or peak.get('linewidth_y', 0)
            self.details_table.setItem(row, 6, QTableWidgetItem(f"{lw_f1:.3f}"))

            # R²
            r2 = peak.get('r_squared', 0)
            r2_item = QTableWidgetItem(f"{r2:.4f}")
            if r2 < 0.7:
                r2_item.setForeground(QColor(ERROR_RED))
            elif r2 < 0.85:
                r2_item.setForeground(QColor(WARNING_ORANGE))
            else:
                r2_item.setForeground(QColor(SUCCESS_GREEN))
            self.details_table.setItem(row, 7, r2_item)

    def export_results(self):
        """Export current results view."""
        if self.main_window and hasattr(self.main_window, 'batch_export_results'):
            self.main_window.batch_export_results()
        else:
            QMessageBox.information(self, "Export",
                "Use File > Export menu from main window for export options.")
