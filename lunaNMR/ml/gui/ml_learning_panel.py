# ABOUTME: ML Learning Center dialog for viewing and controlling ML settings
# ABOUTME: Provides status, predictions, history, and configuration controls

import logging
from typing import Optional

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QTextEdit, QTableWidget, QTableWidgetItem,
    QDoubleSpinBox, QCheckBox, QMessageBox, QHeaderView
)
from PySide6.QtCore import Qt

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.gui.styles.design_system import (
    SPACING_SM, SPACING_MD, SPACING_LG,
    FONT_SIZE_BODY, FONT_SIZE_SECTION_LABEL, FONT_SIZE_SMALL,
    PRIMARY_TEXT, PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER, BUTTON_CORNER_RADIUS, BUTTON_HEIGHT_DIALOG,
    FRAME_BG_COLOR, SUCCESS_GREEN
)

logger = logging.getLogger(__name__)


class MLLearningPanel(BaseDialog):
    """Dialog for viewing and controlling ML Learning Center settings.

    Provides:
        - ML status display (enabled/disabled, active predictions)
        - Model information (version, sample counts, blend weights)
        - Current predictions table (ML vs Stats vs Default)
        - Learning history (sample counts per nucleus)
        - Configuration controls (confidence threshold, fallback options)
        - Action buttons (retrain, reset, export)

    Example:
        ```python
        from lunaNMR.ml.gui import MLLearningPanel
        from lunaNMR.ml import ModelManager

        manager = ModelManager()
        panel = MLLearningPanel(parent, model_manager=manager)
        panel.exec()
        ```
    """

    def __init__(self, parent=None, model_manager=None, main_window=None):
        """Initialize the ML Learning Center panel.

        Args:
            parent: Parent widget
            model_manager: Optional ModelManager instance for ML operations
            main_window: Reference to MainWindow for accessing config
        """
        super().__init__(
            parent=parent,
            title="ML Learning Center",
            default_size=(650, 600),
            modal=True
        )

        self.model_manager = model_manager
        self.main_window = main_window

        # Build UI
        self.setup_ui()

        # Update display
        self.update_status()

        # Center on parent
        if parent:
            self.center_on_parent()
        else:
            self.center_on_screen()

        logger.debug("MLLearningPanel initialized")

    def setup_ui(self):
        """Setup the dialog user interface."""
        layout = QVBoxLayout()
        layout.setSpacing(SPACING_MD)
        layout.setContentsMargins(SPACING_LG, SPACING_LG, SPACING_LG, SPACING_LG)

        # Title and status row
        header_layout = QHBoxLayout()

        title_label = QLabel("ML Learning Center")
        title_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
            }}
        """)
        header_layout.addWidget(title_label)

        header_layout.addStretch()

        # Status indicator
        self.status_label = QLabel("Status: Checking...")
        self.status_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_BODY}px;
                color: {PRIMARY_TEXT};
                padding: {SPACING_SM}px;
            }}
        """)
        header_layout.addWidget(self.status_label)

        layout.addLayout(header_layout)

        # Model Information section
        model_info_group = self._create_model_info_section()
        layout.addWidget(model_info_group)

        # Current Predictions section
        predictions_group = self._create_predictions_section()
        layout.addWidget(predictions_group, stretch=1)

        # Learning History section
        history_group = self._create_history_section()
        layout.addWidget(history_group)

        # Advanced Configuration section
        config_group = self._create_config_section()
        layout.addWidget(config_group)

        # Action buttons
        buttons_layout = self._create_buttons_layout()
        layout.addLayout(buttons_layout)

        self.setLayout(layout)

    def _create_model_info_section(self) -> QGroupBox:
        """Create the model information section."""
        group = QGroupBox("Model Information")
        group.setStyleSheet(self._group_box_style())

        layout = QVBoxLayout()
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        self.model_info_text = QTextEdit()
        self.model_info_text.setReadOnly(True)
        self.model_info_text.setMaximumHeight(100)
        self.model_info_text.setStyleSheet(f"""
            QTextEdit {{
                font-family: 'Courier New', monospace;
                font-size: {FONT_SIZE_SMALL}px;
                background-color: white;
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px;
            }}
        """)
        layout.addWidget(self.model_info_text)

        group.setLayout(layout)
        return group

    def _create_predictions_section(self) -> QGroupBox:
        """Create the current predictions table section."""
        group = QGroupBox("Current Predictions")
        group.setStyleSheet(self._group_box_style())

        layout = QVBoxLayout()
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        # Predictions table
        self.predictions_table = QTableWidget()
        self.predictions_table.setColumnCount(5)
        self.predictions_table.setHorizontalHeaderLabels([
            "Parameter", "ML Pred", "Stats", "Default", "Using"
        ])
        self.predictions_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.predictions_table.setAlternatingRowColors(True)
        self.predictions_table.setStyleSheet(f"""
            QTableWidget {{
                font-size: {FONT_SIZE_SMALL}px;
                background-color: white;
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
            }}
            QTableWidget::item {{
                padding: {SPACING_SM}px;
            }}
            QHeaderView::section {{
                background-color: {SECONDARY_BUTTON_BG};
                padding: {SPACING_SM}px;
                border: none;
                font-weight: bold;
            }}
        """)
        layout.addWidget(self.predictions_table)

        # Info label
        info_label = QLabel("Note: Predictions shown for currently loaded spectrum")
        info_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SMALL}px;
                color: {SECONDARY_BUTTON_TEXT};
                font-style: italic;
            }}
        """)
        layout.addWidget(info_label)

        group.setLayout(layout)
        return group

    def _create_history_section(self) -> QGroupBox:
        """Create the learning history section."""
        group = QGroupBox("Learning History")
        group.setStyleSheet(self._group_box_style())

        layout = QHBoxLayout()
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        self.history_label = QLabel("No training data collected yet")
        self.history_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_BODY}px;
                color: {PRIMARY_TEXT};
            }}
        """)
        layout.addWidget(self.history_label)

        layout.addStretch()

        group.setLayout(layout)
        return group

    def _create_config_section(self) -> QGroupBox:
        """Create the advanced configuration section."""
        group = QGroupBox("Advanced Configuration")
        group.setStyleSheet(self._group_box_style())

        layout = QVBoxLayout()
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        # Confidence threshold row
        threshold_layout = QHBoxLayout()

        threshold_label = QLabel("Confidence Threshold:")
        threshold_label.setStyleSheet(f"font-size: {FONT_SIZE_BODY}px;")
        threshold_layout.addWidget(threshold_label)

        self.confidence_spin = QDoubleSpinBox()
        self.confidence_spin.setRange(0.0, 1.0)
        self.confidence_spin.setSingleStep(0.1)
        self.confidence_spin.setValue(0.6)
        self.confidence_spin.setDecimals(2)
        self.confidence_spin.setStyleSheet(f"""
            QDoubleSpinBox {{
                font-size: {FONT_SIZE_BODY}px;
                padding: {SPACING_SM}px;
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
            }}
        """)
        threshold_layout.addWidget(self.confidence_spin)

        threshold_layout.addStretch()
        layout.addLayout(threshold_layout)

        # Statistical fallback checkbox
        self.stats_fallback_check = QCheckBox("Prefer statistical fallback when available")
        self.stats_fallback_check.setStyleSheet(f"font-size: {FONT_SIZE_BODY}px;")
        layout.addWidget(self.stats_fallback_check)

        # Training data collection checkbox
        self.collect_data_check = QCheckBox("Collect training data from fitting results")
        self.collect_data_check.setStyleSheet(f"font-size: {FONT_SIZE_BODY}px;")
        self.collect_data_check.setToolTip(
            "When enabled, successful fitting results are saved to train the ML model.\n"
            "Data is stored locally in lunaNMR_v1o0/ml_training_data/"
        )
        self.collect_data_check.stateChanged.connect(self._on_collect_data_changed)
        layout.addWidget(self.collect_data_check)

        group.setLayout(layout)
        return group

    def _create_buttons_layout(self) -> QHBoxLayout:
        """Create the action buttons layout."""
        layout = QHBoxLayout()

        # Retrain button
        self.retrain_button = QPushButton("Force Retrain Now")
        self.retrain_button.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.retrain_button.setStyleSheet(self._secondary_button_style())
        self.retrain_button.clicked.connect(self._on_retrain_clicked)
        layout.addWidget(self.retrain_button)

        # Reset button
        self.reset_button = QPushButton("Reset to Pretrained")
        self.reset_button.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.reset_button.setStyleSheet(f"""
            QPushButton {{
                background-color: #FF3B30;
                color: white;
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
            }}
            QPushButton:hover {{
                background-color: #FF2D20;
            }}
        """)
        self.reset_button.clicked.connect(self._on_reset_clicked)
        layout.addWidget(self.reset_button)

        layout.addStretch()

        # Close button
        close_button = QPushButton("Close")
        close_button.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        close_button.setMinimumWidth(100)
        close_button.setStyleSheet(self._primary_button_style())
        close_button.clicked.connect(self.accept)
        layout.addWidget(close_button)

        return layout

    def _group_box_style(self) -> str:
        """Return common QGroupBox styling."""
        return f"""
            QGroupBox {{
                font-size: {FONT_SIZE_BODY}px;
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

    def _primary_button_style(self) -> str:
        """Return primary button styling."""
        return f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
        """

    def _secondary_button_style(self) -> str:
        """Return secondary button styling."""
        return f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
        """

    def update_status(self):
        """Update all status displays with current ML state."""
        self._update_status_indicator()
        self._update_model_info()
        self._update_predictions_table()
        self._update_history_display()
        self._update_config_controls()

    def _update_status_indicator(self):
        """Update the main status indicator."""
        if self.model_manager is None:
            self.status_label.setText("Status: No Model Manager")
            self.status_label.setStyleSheet(f"""
                QLabel {{
                    font-size: {FONT_SIZE_BODY}px;
                    color: #FF3B30;
                    padding: {SPACING_SM}px;
                }}
            """)
            return

        status = self.model_manager.get_status()
        ml_available = status.get("ml_available", False)
        stats_available = status.get("stats_available", False)

        if ml_available:
            self.status_label.setText("Status: Active (ML Ready)")
            color = SUCCESS_GREEN
        elif stats_available:
            self.status_label.setText("Status: Active (Stats Only)")
            color = "#FF9500"  # Orange
        else:
            self.status_label.setText("Status: Inactive")
            color = "#FF3B30"  # Red

        self.status_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_BODY}px;
                color: {color};
                padding: {SPACING_SM}px;
                font-weight: bold;
            }}
        """)

    def _update_model_info(self):
        """Update the model information display."""
        if self.model_manager is None:
            self.model_info_text.setPlainText("No model manager available")
            return

        status = self.model_manager.get_status()

        info_lines = []

        # Pre-trained model info
        if status.get("ml_available"):
            info_lines.append("Pre-trained model: Available")
        else:
            info_lines.append("Pre-trained model: Not loaded")

        # Stats predictor info
        stats_15n = status.get("stats_samples_15N", 0)
        stats_13c = status.get("stats_samples_13C", 0)
        info_lines.append(f"Statistical samples: 15N={stats_15n}, 13C={stats_13c}")

        # Local adaptation info
        adaptation_status = status.get("adaptation", {})
        if adaptation_status:
            local_15n = adaptation_status.get("15N", {}).get("has_local_model", False)
            local_13c = adaptation_status.get("13C", {}).get("has_local_model", False)
            info_lines.append(f"Local adaptation: 15N={'Yes' if local_15n else 'No'}, "
                            f"13C={'Yes' if local_13c else 'No'}")

        self.model_info_text.setPlainText("\n".join(info_lines))

    def _update_predictions_table(self):
        """Update the predictions table with current values."""
        # Define parameters to display
        params = [
            ("Linewidth F1 (ppm)", "lw_f1_median", 0.40),
            ("Linewidth F2 (ppm)", "lw_f2_median", 0.020),
            ("Window rad F1", "rad_f1", 0.40),
            ("Window rad F2", "rad_f2", 0.04),
            ("Overlap F1", "overlap_threshold_f1", 0.40),
            ("Overlap F2", "overlap_threshold_f2", 0.04),
            ("Expected R²", "achievable_r2", 0.85),
        ]

        self.predictions_table.setRowCount(len(params))

        for row, (name, key, default) in enumerate(params):
            # Parameter name
            self.predictions_table.setItem(row, 0, QTableWidgetItem(name))

            # ML prediction
            ml_val = self._get_ml_prediction(key)
            self.predictions_table.setItem(row, 1, QTableWidgetItem(
                f"{ml_val:.3f}" if ml_val is not None else "-"
            ))

            # Stats prediction
            stats_val = self._get_stats_prediction(key)
            self.predictions_table.setItem(row, 2, QTableWidgetItem(
                f"{stats_val:.3f}" if stats_val is not None else "-"
            ))

            # Default
            self.predictions_table.setItem(row, 3, QTableWidgetItem(f"{default:.3f}"))

            # Currently using
            using = self._get_active_source(key, ml_val, stats_val, default)
            self.predictions_table.setItem(row, 4, QTableWidgetItem(using))

    def _get_ml_prediction(self, key: str) -> Optional[float]:
        """Get ML prediction for a parameter."""
        if self.model_manager is None:
            return None
        # Would get from current prediction if available
        return None

    def _get_stats_prediction(self, key: str) -> Optional[float]:
        """Get stats prediction for a parameter."""
        if self.model_manager is None:
            return None
        # Would get from stats predictor if available
        return None

    def _get_active_source(self, key: str, ml_val, stats_val, default) -> str:
        """Determine which source is currently active."""
        if ml_val is not None:
            return "ML"
        elif stats_val is not None:
            return "Stats"
        else:
            return "Default"

    def _update_history_display(self):
        """Update the learning history display."""
        if self.model_manager is None:
            self.history_label.setText("No training data collected yet")
            return

        status = self.model_manager.get_status()
        stats_15n = status.get("stats_samples_15N", 0)
        stats_13c = status.get("stats_samples_13C", 0)

        if stats_15n == 0 and stats_13c == 0:
            self.history_label.setText("No training data collected yet")
        else:
            self.history_label.setText(
                f"Collected samples: 15N = {stats_15n} | 13C = {stats_13c}"
            )

    def _update_config_controls(self):
        """Update configuration controls from config manager."""
        if self.main_window and hasattr(self.main_window, 'config_manager'):
            config = self.main_window.config_manager.config
            ml_config = config.get('ml_learning', {})

            self.confidence_spin.setValue(
                ml_config.get('ml_confidence_threshold', 0.6)
            )
            self.stats_fallback_check.setChecked(
                ml_config.get('use_statistical_fallback', True)
            )
            self.collect_data_check.setChecked(
                ml_config.get('collect_training_data', False)
            )

    def _on_collect_data_changed(self, state):
        """Handle collect training data checkbox change."""
        if self.main_window and hasattr(self.main_window, 'config_manager'):
            config = self.main_window.config_manager.config
            if 'ml_learning' not in config:
                config['ml_learning'] = {}

            enabled = self.collect_data_check.isChecked()
            config['ml_learning']['collect_training_data'] = enabled
            self.main_window.config_manager.save_config()

            status = "enabled" if enabled else "disabled"
            logger.info(f"Training data collection {status}")

    def _on_retrain_clicked(self):
        """Handle Force Retrain button click."""
        reply = QMessageBox.question(
            self,
            "Force Retrain",
            "Retrain local models from collected data?\n\n"
            "This will update the local adaptation based on your data.",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            if self.model_manager:
                # Trigger retraining
                logger.info("Force retrain requested")
                QMessageBox.information(
                    self,
                    "Retrain",
                    "Retraining not yet implemented.\n\n"
                    "Local models will automatically retrain\n"
                    "when enough samples are collected."
                )
            else:
                QMessageBox.warning(
                    self,
                    "Warning",
                    "No model manager available"
                )

    def _on_reset_clicked(self):
        """Handle Reset to Pretrained button click."""
        reply = QMessageBox.question(
            self,
            "Reset to Pretrained",
            "Reset all local adaptations and return to pretrained models?\n\n"
            "This will clear all collected training data.",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            if self.model_manager:
                # Clear local adaptations
                self.model_manager.stats_predictor.clear()
                logger.info("Reset to pretrained models")
                QMessageBox.information(
                    self,
                    "Reset Complete",
                    "Local adaptations cleared.\n"
                    "Using pretrained models only."
                )
                self.update_status()
            else:
                QMessageBox.warning(
                    self,
                    "Warning",
                    "No model manager available"
                )

    def apply_config(self):
        """Apply configuration changes to config manager."""
        if self.main_window and hasattr(self.main_window, 'config_manager'):
            config = self.main_window.config_manager.config
            if 'ml_learning' not in config:
                config['ml_learning'] = {}

            config['ml_learning']['ml_confidence_threshold'] = self.confidence_spin.value()
            config['ml_learning']['use_statistical_fallback'] = self.stats_fallback_check.isChecked()
            config['ml_learning']['collect_training_data'] = self.collect_data_check.isChecked()

            self.main_window.config_manager.save_config()
            logger.info("ML configuration saved")
