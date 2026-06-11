# ABOUTME: Viewer for LunaNMR series integration results within DynamiXs
# ABOUTME: Shows fit quality, intensity matrix heatmap, and export options

"""
Series Results Viewer - Display LunaNMR series integration results in DynamiXs

Shows processing results from MultiSpectrumProcessor including:
- Summary statistics
- Per-residue fit quality (R² plot)
- Intensity matrix heatmap
- Quality filtering and export
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, List, Any

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFrame, QTabWidget,
    QLabel, QPushButton, QTableWidget, QTableWidgetItem,
    QHeaderView, QSplitter, QScrollArea, QWidget, QComboBox,
    QCheckBox, QSpinBox, QGroupBox, QGridLayout, QFileDialog,
    QMessageBox
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

# Import DynamiXs style components
import sys
sys.path.append(str(Path(__file__).parent.parent))
from constants import (
    BG_COLOR, PANEL_BG_COLOR, FRAME_BG_COLOR,
    PRIMARY_TEXT, SECONDARY_TEXT,
    SUCCESS_GREEN, WARNING_ORANGE, ERROR_RED,
    SPACING_XS, SPACING_SM, SPACING_MD,
    FONT_SIZE_SECTION_LABEL, FONT_SIZE_BODY, FONT_SIZE_SMALL
)
from gui_components import (
    create_primary_button, create_secondary_button, create_label,
    create_secondary_label, get_font, show_info, show_error,
    StyledGroupBox
)


class SeriesResultsViewer(QDialog):
    """
    Dialog for viewing LunaNMR series integration results.

    Displays processing summary, fit quality metrics, and intensity data
    from BatchResults or tidy CSV output.
    """

    # Quality thresholds for color coding
    R2_EXCELLENT = 0.95
    R2_GOOD = 0.85
    R2_POOR = 0.70

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        batch_results: Optional[Any] = None,
        tidy_csv_path: Optional[str] = None,
        experiment_type: str = "T1"
    ):
        """
        Initialize the Series Results Viewer.

        Args:
            parent: Parent widget
            batch_results: BatchResults object from MultiSpectrumProcessor
            tidy_csv_path: Path to tidy CSV output file (alternative to batch_results)
            experiment_type: Type of experiment (T1, T2, hetNOE)
        """
        super().__init__(parent)

        self.batch_results = batch_results
        self.tidy_csv_path = tidy_csv_path
        self.experiment_type = experiment_type

        self._df = None  # Loaded data frame
        self._quality_data = None  # Per-residue quality metrics

        self._setup_ui()
        self._load_data()

    def _setup_ui(self):
        """Set up the dialog UI."""
        self.setWindowTitle(f"Series Integration Results - {self.experiment_type}")
        self.setMinimumSize(1200, 800)
        self.setStyleSheet(f"background-color: {BG_COLOR}; color: {PRIMARY_TEXT};")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        layout.setSpacing(SPACING_MD)

        # Tab widget for different views
        self.tab_widget = QTabWidget()
        self.tab_widget.setStyleSheet(f"""
            QTabWidget::pane {{
                background-color: {PANEL_BG_COLOR};
                border: 1px solid {FRAME_BG_COLOR};
                border-radius: 4px;
            }}
            QTabBar::tab {{
                background-color: {FRAME_BG_COLOR};
                color: {PRIMARY_TEXT};
                padding: 8px 16px;
                margin-right: 2px;
                border-top-left-radius: 4px;
                border-top-right-radius: 4px;
            }}
            QTabBar::tab:selected {{
                background-color: {PANEL_BG_COLOR};
            }}
        """)

        # Summary tab
        self._create_summary_tab()

        # Quality plot tab
        self._create_quality_tab()

        # Intensity heatmap tab
        self._create_heatmap_tab()

        # Data table tab
        self._create_table_tab()

        layout.addWidget(self.tab_widget, stretch=1)

        # Bottom button bar
        button_bar = QFrame()
        button_layout = QHBoxLayout(button_bar)
        button_layout.setContentsMargins(0, 0, 0, 0)

        self.export_btn = create_primary_button("Export Results", self._on_export)
        button_layout.addWidget(self.export_btn)

        button_layout.addStretch()

        self.close_btn = create_secondary_button("Close", self.accept)
        button_layout.addWidget(self.close_btn)

        layout.addWidget(button_bar)

    def _create_summary_tab(self):
        """Create the summary statistics tab."""
        summary_widget = QWidget()
        layout = QVBoxLayout(summary_widget)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        # Summary group
        summary_group = StyledGroupBox("Processing Summary")
        summary_layout = QGridLayout()
        summary_group.setLayout(summary_layout)

        # Summary labels (will be populated in _load_data)
        self.total_spectra_label = create_label("Total Spectra: --")
        self.successful_label = create_label("Successful: --")
        self.failed_label = create_label("Failed: --")
        self.total_peaks_label = create_label("Total Peaks: --")
        self.avg_r2_label = create_label("Avg R²: --")

        summary_layout.addWidget(self.total_spectra_label, 0, 0)
        summary_layout.addWidget(self.successful_label, 0, 1)
        summary_layout.addWidget(self.failed_label, 0, 2)
        summary_layout.addWidget(self.total_peaks_label, 1, 0)
        summary_layout.addWidget(self.avg_r2_label, 1, 1)

        layout.addWidget(summary_group)

        # Quality distribution group
        quality_group = StyledGroupBox("Fit Quality Distribution")
        quality_layout = QGridLayout()
        quality_group.setLayout(quality_layout)

        self.excellent_label = create_label(f"Excellent (R² ≥ {self.R2_EXCELLENT}): --")
        self.good_label = create_label(f"Good ({self.R2_GOOD} ≤ R² < {self.R2_EXCELLENT}): --")
        self.poor_label = create_label(f"Poor ({self.R2_POOR} ≤ R² < {self.R2_GOOD}): --")
        self.failed_fits_label = create_label(f"Failed (R² < {self.R2_POOR}): --")

        quality_layout.addWidget(self.excellent_label, 0, 0)
        quality_layout.addWidget(self.good_label, 0, 1)
        quality_layout.addWidget(self.poor_label, 1, 0)
        quality_layout.addWidget(self.failed_fits_label, 1, 1)

        layout.addWidget(quality_group)

        layout.addStretch()

        self.tab_widget.addTab(summary_widget, "Summary")

    def _create_quality_tab(self):
        """Create the fit quality plot tab."""
        quality_widget = QWidget()
        layout = QVBoxLayout(quality_widget)
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)

        # Matplotlib figure for R² plot
        self.quality_figure = Figure(figsize=(12, 6), dpi=100)
        self.quality_figure.patch.set_facecolor('#1e1e1e')
        self.quality_canvas = FigureCanvasQTAgg(self.quality_figure)
        self.quality_toolbar = NavigationToolbar2QT(self.quality_canvas, quality_widget)

        layout.addWidget(self.quality_toolbar)
        layout.addWidget(self.quality_canvas, stretch=1)

        self.tab_widget.addTab(quality_widget, "Fit Quality")

    def _create_heatmap_tab(self):
        """Create the intensity heatmap tab."""
        heatmap_widget = QWidget()
        layout = QVBoxLayout(heatmap_widget)
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)

        # Matplotlib figure for heatmap
        self.heatmap_figure = Figure(figsize=(12, 8), dpi=100)
        self.heatmap_figure.patch.set_facecolor('#1e1e1e')
        self.heatmap_canvas = FigureCanvasQTAgg(self.heatmap_figure)
        self.heatmap_toolbar = NavigationToolbar2QT(self.heatmap_canvas, heatmap_widget)

        layout.addWidget(self.heatmap_toolbar)
        layout.addWidget(self.heatmap_canvas, stretch=1)

        self.tab_widget.addTab(heatmap_widget, "Intensity Heatmap")

    def _create_table_tab(self):
        """Create the data table tab."""
        table_widget = QWidget()
        layout = QVBoxLayout(table_widget)
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)

        # Filter controls
        filter_frame = QFrame()
        filter_layout = QHBoxLayout(filter_frame)
        filter_layout.setContentsMargins(0, 0, 0, SPACING_SM)

        filter_layout.addWidget(create_label("Quality Filter:"))

        self.quality_filter_combo = QComboBox()
        self.quality_filter_combo.addItems([
            "All",
            f"Excellent (R² ≥ {self.R2_EXCELLENT})",
            f"Good+ (R² ≥ {self.R2_GOOD})",
            f"Poor+ (R² ≥ {self.R2_POOR})"
        ])
        self.quality_filter_combo.currentIndexChanged.connect(self._apply_table_filter)
        filter_layout.addWidget(self.quality_filter_combo)

        filter_layout.addStretch()

        layout.addWidget(filter_frame)

        # Data table
        self.data_table = QTableWidget()
        self.data_table.setStyleSheet(f"""
            QTableWidget {{
                background-color: {PANEL_BG_COLOR};
                color: {PRIMARY_TEXT};
                gridline-color: {FRAME_BG_COLOR};
            }}
            QHeaderView::section {{
                background-color: {FRAME_BG_COLOR};
                color: {PRIMARY_TEXT};
                padding: 4px;
                border: none;
            }}
        """)
        self.data_table.setAlternatingRowColors(True)
        self.data_table.horizontalHeader().setStretchLastSection(True)

        layout.addWidget(self.data_table, stretch=1)

        self.tab_widget.addTab(table_widget, "Data Table")

    def _load_data(self):
        """Load data from batch_results or tidy CSV."""
        if self.batch_results is not None:
            self._load_from_batch_results()
        elif self.tidy_csv_path:
            self._load_from_tidy_csv()
        else:
            self._show_empty_state()

    def _load_from_batch_results(self):
        """Load data from BatchResults object."""
        if not self.batch_results or not self.batch_results.results:
            self._show_empty_state()
            return

        # Extract summary
        summary = self.batch_results.get_summary()
        self.total_spectra_label.setText(f"Total Spectra: {summary.get('total_spectra', 0)}")
        self.successful_label.setText(f"Successful: {summary.get('successful', 0)}")
        self.failed_label.setText(f"Failed: {summary.get('failed', 0)}")

        # Build dataframe from results
        all_data = []
        for spectrum_name, result in self.batch_results.results.items():
            if result.get('integration_results'):
                for integration in result['integration_results']:
                    row = {
                        'spectrum_name': spectrum_name,
                        'assignment': integration.get('assignment', ''),
                        'height': integration.get('height', 0),
                        'volume': integration.get('volume', 0),
                        'r_squared': integration.get('r_squared', 0),
                        'ppm_x': integration.get('ppm_x', 0),
                        'ppm_y': integration.get('ppm_y', 0),
                    }
                    all_data.append(row)

        if all_data:
            self._df = pd.DataFrame(all_data)
            self._calculate_quality_metrics()
            self._update_summary_stats()
            self._plot_quality()
            self._plot_heatmap()
            self._populate_table()
        else:
            self._show_empty_state()

    def _load_from_tidy_csv(self):
        """Load data from tidy CSV file."""
        try:
            self._df = pd.read_csv(self.tidy_csv_path)

            # Map column names to expected format
            column_mapping = {
                'spectrum': 'spectrum_name',
                'Spectrum': 'spectrum_name',
            }
            self._df = self._df.rename(columns=column_mapping)

            # Count unique spectra
            n_spectra = self._df['spectrum_name'].nunique() if 'spectrum_name' in self._df.columns else 0
            self.total_spectra_label.setText(f"Total Spectra: {n_spectra}")
            self.successful_label.setText(f"Successful: {n_spectra}")
            self.failed_label.setText("Failed: 0")

            self._calculate_quality_metrics()
            self._update_summary_stats()
            self._plot_quality()
            self._plot_heatmap()
            self._populate_table()

        except Exception as e:
            show_error(self, "Error", f"Failed to load CSV: {e}")
            self._show_empty_state()

    def _calculate_quality_metrics(self):
        """Calculate per-residue quality metrics from the dataframe."""
        if self._df is None or self._df.empty:
            return

        # Group by assignment and calculate mean R²
        if 'r_squared' in self._df.columns and 'assignment' in self._df.columns:
            self._quality_data = self._df.groupby('assignment').agg({
                'r_squared': 'mean',
                'height': 'mean',
                'volume': 'mean'
            }).reset_index()
            self._quality_data = self._quality_data.rename(columns={'r_squared': 'mean_r2'})
        else:
            self._quality_data = None

    def _update_summary_stats(self):
        """Update summary statistics labels."""
        if self._df is None:
            return

        # Total peaks
        n_peaks = self._df['assignment'].nunique() if 'assignment' in self._df.columns else 0
        self.total_peaks_label.setText(f"Total Peaks: {n_peaks}")

        # Average R²
        if 'r_squared' in self._df.columns:
            avg_r2 = self._df['r_squared'].mean()
            self.avg_r2_label.setText(f"Avg R²: {avg_r2:.3f}")

            # Quality distribution
            r2_values = self._df['r_squared'].dropna()
            n_excellent = (r2_values >= self.R2_EXCELLENT).sum()
            n_good = ((r2_values >= self.R2_GOOD) & (r2_values < self.R2_EXCELLENT)).sum()
            n_poor = ((r2_values >= self.R2_POOR) & (r2_values < self.R2_GOOD)).sum()
            n_failed = (r2_values < self.R2_POOR).sum()

            self.excellent_label.setText(f"Excellent (R² ≥ {self.R2_EXCELLENT}): {n_excellent}")
            self.good_label.setText(f"Good ({self.R2_GOOD} ≤ R² < {self.R2_EXCELLENT}): {n_good}")
            self.poor_label.setText(f"Poor ({self.R2_POOR} ≤ R² < {self.R2_GOOD}): {n_poor}")
            self.failed_fits_label.setText(f"Failed (R² < {self.R2_POOR}): {n_failed}")

    def _plot_quality(self):
        """Plot the fit quality (R² vs residue) chart."""
        self.quality_figure.clear()

        if self._quality_data is None or self._quality_data.empty:
            ax = self.quality_figure.add_subplot(111)
            ax.text(0.5, 0.5, "No quality data available",
                   ha='center', va='center', color='white', fontsize=12)
            ax.set_facecolor('#1e1e1e')
            self.quality_canvas.draw()
            return

        ax = self.quality_figure.add_subplot(111)
        ax.set_facecolor('#2d2d2d')

        # Sort by residue number if possible
        data = self._quality_data.copy()
        try:
            # Extract residue number from assignment (e.g., "142.ALA" -> 142)
            data['residue_num'] = data['assignment'].str.extract(r'(\d+)').astype(float)
            data = data.sort_values('residue_num')
        except Exception:
            pass

        # Plot R² values with color coding
        colors = []
        for r2 in data['mean_r2']:
            if r2 >= self.R2_EXCELLENT:
                colors.append('#4CAF50')  # Green
            elif r2 >= self.R2_GOOD:
                colors.append('#2196F3')  # Blue
            elif r2 >= self.R2_POOR:
                colors.append('#FF9800')  # Orange
            else:
                colors.append('#F44336')  # Red

        x_positions = range(len(data))
        bars = ax.bar(x_positions, data['mean_r2'], color=colors, edgecolor='none', alpha=0.8)

        # Add threshold lines
        ax.axhline(y=self.R2_EXCELLENT, color='#4CAF50', linestyle='--', alpha=0.5,
                   label=f'Excellent ({self.R2_EXCELLENT})')
        ax.axhline(y=self.R2_GOOD, color='#2196F3', linestyle='--', alpha=0.5,
                   label=f'Good ({self.R2_GOOD})')
        ax.axhline(y=self.R2_POOR, color='#FF9800', linestyle='--', alpha=0.5,
                   label=f'Poor ({self.R2_POOR})')

        # Formatting
        ax.set_xlabel('Residue', color='white', fontsize=10)
        ax.set_ylabel('R²', color='white', fontsize=10)
        ax.set_title(f'{self.experiment_type} Fit Quality by Residue', color='white', fontsize=12)
        ax.tick_params(colors='white', labelsize=8)
        ax.set_ylim(0, 1.05)

        # X-axis labels (show every Nth if too many)
        labels = data['assignment'].tolist()
        if len(labels) > 30:
            step = len(labels) // 20
            ax.set_xticks(x_positions[::step])
            ax.set_xticklabels(labels[::step], rotation=45, ha='right')
        else:
            ax.set_xticks(x_positions)
            ax.set_xticklabels(labels, rotation=45, ha='right')

        ax.legend(loc='lower right', fontsize=8, facecolor='#2d2d2d', edgecolor='gray',
                 labelcolor='white')

        self.quality_figure.tight_layout()
        self.quality_canvas.draw()

    def _plot_heatmap(self):
        """Plot intensity heatmap (residues vs spectra)."""
        self.heatmap_figure.clear()

        if self._df is None or self._df.empty:
            ax = self.heatmap_figure.add_subplot(111)
            ax.text(0.5, 0.5, "No data available for heatmap",
                   ha='center', va='center', color='white', fontsize=12)
            ax.set_facecolor('#1e1e1e')
            self.heatmap_canvas.draw()
            return

        # Pivot to create matrix (assignment x spectrum_name)
        try:
            intensity_col = 'height' if 'height' in self._df.columns else 'volume'
            pivot_df = self._df.pivot_table(
                index='assignment',
                columns='spectrum_name',
                values=intensity_col,
                aggfunc='first'
            )

            # Sort residues by number
            try:
                residue_order = sorted(pivot_df.index,
                                       key=lambda x: float(x.split('.')[0]) if '.' in str(x) else float(x))
                pivot_df = pivot_df.reindex(residue_order)
            except Exception:
                pass

            ax = self.heatmap_figure.add_subplot(111)

            # Plot heatmap
            im = ax.imshow(pivot_df.values, aspect='auto', cmap='viridis')

            # Colorbar
            cbar = self.heatmap_figure.colorbar(im, ax=ax)
            cbar.set_label(intensity_col.capitalize(), color='white')
            cbar.ax.yaxis.set_tick_params(color='white')
            plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

            # Labels
            ax.set_xlabel('Spectrum', color='white', fontsize=10)
            ax.set_ylabel('Residue', color='white', fontsize=10)
            ax.set_title(f'{self.experiment_type} Intensity Matrix', color='white', fontsize=12)

            # Axis ticks
            ax.set_xticks(range(len(pivot_df.columns)))
            ax.set_xticklabels(pivot_df.columns, rotation=45, ha='right', fontsize=7, color='white')

            if len(pivot_df.index) <= 50:
                ax.set_yticks(range(len(pivot_df.index)))
                ax.set_yticklabels(pivot_df.index, fontsize=7, color='white')
            else:
                # Show subset
                step = len(pivot_df.index) // 30
                ax.set_yticks(range(0, len(pivot_df.index), step))
                ax.set_yticklabels(pivot_df.index[::step], fontsize=7, color='white')

            ax.tick_params(colors='white')

        except Exception as e:
            ax = self.heatmap_figure.add_subplot(111)
            ax.text(0.5, 0.5, f"Error creating heatmap: {e}",
                   ha='center', va='center', color='white', fontsize=10)
            ax.set_facecolor('#1e1e1e')

        self.heatmap_figure.tight_layout()
        self.heatmap_canvas.draw()

    def _populate_table(self):
        """Populate the data table with results."""
        if self._df is None or self._df.empty:
            return

        self.data_table.clear()

        # Set columns
        columns = list(self._df.columns)
        self.data_table.setColumnCount(len(columns))
        self.data_table.setHorizontalHeaderLabels(columns)

        # Populate rows
        self.data_table.setRowCount(len(self._df))
        for i, (_, row) in enumerate(self._df.iterrows()):
            for j, col in enumerate(columns):
                value = row[col]
                if isinstance(value, float):
                    item = QTableWidgetItem(f"{value:.4f}")
                else:
                    item = QTableWidgetItem(str(value))

                # Color code R² values
                if col == 'r_squared' and isinstance(value, (int, float)):
                    if value >= self.R2_EXCELLENT:
                        item.setBackground(QColor(76, 175, 80, 100))  # Green
                    elif value >= self.R2_GOOD:
                        item.setBackground(QColor(33, 150, 243, 100))  # Blue
                    elif value >= self.R2_POOR:
                        item.setBackground(QColor(255, 152, 0, 100))  # Orange
                    else:
                        item.setBackground(QColor(244, 67, 54, 100))  # Red

                self.data_table.setItem(i, j, item)

        self.data_table.resizeColumnsToContents()

    def _apply_table_filter(self):
        """Apply quality filter to the table."""
        if self._df is None:
            return

        filter_index = self.quality_filter_combo.currentIndex()
        thresholds = [0, self.R2_EXCELLENT, self.R2_GOOD, self.R2_POOR]

        if filter_index == 0 or 'r_squared' not in self._df.columns:
            # Show all
            filtered_df = self._df
        else:
            threshold = thresholds[filter_index]
            filtered_df = self._df[self._df['r_squared'] >= threshold]

        # Repopulate table with filtered data
        columns = list(filtered_df.columns)
        self.data_table.setRowCount(len(filtered_df))
        for i, (_, row) in enumerate(filtered_df.iterrows()):
            for j, col in enumerate(columns):
                value = row[col]
                if isinstance(value, float):
                    item = QTableWidgetItem(f"{value:.4f}")
                else:
                    item = QTableWidgetItem(str(value))
                self.data_table.setItem(i, j, item)

    def _show_empty_state(self):
        """Show empty state when no data is available."""
        self.total_spectra_label.setText("Total Spectra: 0")
        self.successful_label.setText("Successful: 0")
        self.failed_label.setText("Failed: 0")
        self.total_peaks_label.setText("Total Peaks: 0")
        self.avg_r2_label.setText("Avg R²: --")

    def _on_export(self):
        """Export results to CSV."""
        if self._df is None or self._df.empty:
            show_error(self, "Error", "No data to export")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Results",
            f"{self.experiment_type}_series_results.csv",
            "CSV Files (*.csv)"
        )

        if file_path:
            try:
                self._df.to_csv(file_path, index=False)
                show_info(self, "Success", f"Results exported to:\n{file_path}")
            except Exception as e:
                show_error(self, "Error", f"Failed to export: {e}")
