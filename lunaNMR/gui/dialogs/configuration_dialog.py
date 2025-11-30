# ABOUTME: Configuration dialog for viewing and managing lunaNMR settings
# ABOUTME: Provides save/load/reset functionality and displays current parameters

import os
import logging
from typing import Optional

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QTextEdit, QFileDialog, QMessageBox, QScrollArea, QWidget
)
from PySide6.QtCore import Qt

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.gui.styles.design_system import (
    SPACING_SM, SPACING_MD, SPACING_LG,
    FONT_SIZE_BODY, FONT_SIZE_SECTION_LABEL, FONT_SIZE_SMALL,
    PRIMARY_TEXT, SECONDARY_TEXT,
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER, BUTTON_CORNER_RADIUS, BUTTON_HEIGHT_DIALOG,
    FRAME_BG_COLOR
)

logger = logging.getLogger(__name__)


class ConfigurationDialog(BaseDialog):
    """Dialog for viewing and managing lunaNMR configuration settings.

    This dialog provides:
        - Save configuration to file
        - Load configuration from file
        - Reset configuration to defaults
        - Display of all current settings

    Based on v0.9 setup_config_tab() (main_gui.py:2017-2058)

    Example:
        ```python
        dialog = ConfigurationDialog(parent, main_window)
        dialog.exec()
        ```
    """

    def __init__(self, parent=None, main_window=None):
        """Initialize the configuration dialog.

        Args:
            parent: Parent widget
            main_window: Reference to MainWindow for accessing config_manager and parameters
        """
        super().__init__(
            parent=parent,
            title="Configuration - lunaNMR v1.0",
            default_size=(600, 500),
            modal=True
        )

        self.main_window = main_window

        # Build UI
        self.setup_ui()

        # Update display with current settings
        self.update_config_display()

        # Center on parent
        if parent:
            self.center_on_parent()
        else:
            self.center_on_screen()

        logger.debug("ConfigurationDialog initialized")

    def setup_ui(self):
        """Setup the dialog user interface."""
        layout = QVBoxLayout()
        layout.setSpacing(SPACING_MD)
        layout.setContentsMargins(SPACING_LG, SPACING_LG, SPACING_LG, SPACING_LG)

        # Title label
        title_label = QLabel("Configuration Management")
        title_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
                padding-bottom: {SPACING_SM}px;
            }}
        """)
        layout.addWidget(title_label)

        # Configuration Files section (buttons)
        config_buttons_group = self.create_buttons_section()
        layout.addWidget(config_buttons_group)

        # Current Settings section (text display)
        settings_group = self.create_settings_section()
        layout.addWidget(settings_group, stretch=1)

        # Close button
        button_layout = QHBoxLayout()
        button_layout.addStretch()

        close_button = QPushButton("Close")
        close_button.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        close_button.setMinimumWidth(120)
        close_button.setStyleSheet(f"""
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
        """)
        close_button.clicked.connect(self.accept)
        button_layout.addWidget(close_button)

        button_layout.addStretch()
        layout.addLayout(button_layout)

        self.setLayout(layout)

    def create_buttons_section(self) -> QGroupBox:
        """Create the configuration file buttons section.

        Returns:
            QGroupBox containing Save/Load/Reset buttons
        """
        group = QGroupBox("Configuration Files")
        group.setStyleSheet(f"""
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
        """)

        layout = QHBoxLayout()
        layout.setSpacing(SPACING_SM)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        # Save Configuration button (primary)
        save_button = QPushButton("Save Configuration")
        save_button.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        save_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
        """)
        save_button.clicked.connect(self.save_config)
        layout.addWidget(save_button)

        # Load Configuration button (secondary)
        load_button = QPushButton("Load Configuration")
        load_button.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        load_button.setStyleSheet(f"""
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
        """)
        load_button.clicked.connect(self.load_config)
        layout.addWidget(load_button)

        # Reset to Defaults button (destructive)
        reset_button = QPushButton("Reset to Defaults")
        reset_button.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        reset_button.setStyleSheet(f"""
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
        reset_button.clicked.connect(self.reset_config)
        layout.addWidget(reset_button)

        layout.addStretch()

        group.setLayout(layout)
        return group

    def create_settings_section(self) -> QGroupBox:
        """Create the current settings display section.

        Returns:
            QGroupBox containing settings text display
        """
        group = QGroupBox("Current Settings")
        group.setStyleSheet(f"""
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
        """)

        layout = QVBoxLayout()
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        # Settings text display (read-only)
        self.config_text = QTextEdit()
        self.config_text.setReadOnly(True)
        self.config_text.setStyleSheet(f"""
            QTextEdit {{
                font-family: 'Courier New', monospace;
                font-size: {FONT_SIZE_SMALL}px;
                background-color: white;
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px;
            }}
        """)
        layout.addWidget(self.config_text)

        group.setLayout(layout)
        return group

    def update_config_display(self):
        """Update configuration display with current settings.

        Based on v0.9 update_config_display() (main_gui.py:6166-6205)
        """
        if not self.main_window:
            self.config_text.setPlainText("No main window reference - cannot display settings")
            return

        mw = self.main_window

        # Build configuration info string (matching v0.9 format)
        config_info = "Configuration Status:\n\n"

        # Detection Parameters
        config_info += "Detection Parameters:\n"
        config_info += f"  Noise Threshold: {getattr(mw, 'noise_threshold', 3.0):.2f}\n"
        config_info += f"  Search Window X: ±{getattr(mw, 'search_window_x', 0.05):.3f} ppm\n"
        config_info += f"  Search Window Y: ±{getattr(mw, 'search_window_y', 0.5):.2f} ppm\n"
        config_info += f"  Use Reference Detection: {getattr(mw, 'use_reference_detection', True)}\n"
        config_info += "\n"

        # Fitting Parameters
        config_info += "Fitting Parameters:\n"
        config_info += f"  Fitting Window X: {getattr(mw, 'fitting_window_x', 0.2):.1f} ppm\n"
        config_info += f"  Fitting Window Y: {getattr(mw, 'fitting_window_y', 2.0):.1f} ppm\n"
        config_info += f"  Min R²: {getattr(mw, 'min_r_squared', 0.5):.2f}\n"
        config_info += f"  Max Iterations: {getattr(mw, 'max_iterations', 1000)}\n"
        config_info += f"  Fix Positions: {getattr(mw, 'fix_positions', False)}\n"
        config_info += f"  Fix Linewidths: {getattr(mw, 'fix_linewidths', False)}\n"
        config_info += "\n"

        # Display Options
        config_info += "Display Options:\n"
        config_info += f"  Show Detected: {getattr(mw, 'show_detected', True)}\n"
        config_info += f"  Show Assigned: {getattr(mw, 'show_assigned', True)}\n"
        config_info += f"  Show Fitted Curves: {getattr(mw, 'show_fitted_curves', True)}\n"
        config_info += "\n"

        # Parallel Processing
        config_info += "Parallel Processing:\n"
        config_info += f"  Use Parallel: {getattr(mw, 'use_parallel', False)}\n"
        config_info += "\n"

        # File Paths
        config_info += "File Paths:\n"
        nmr_file = getattr(mw, 'current_nmr_file', None)
        peak_file = getattr(mw, 'current_peak_file', None)
        config_info += f"  Current NMR File: {os.path.basename(nmr_file) if nmr_file else 'None'}\n"
        config_info += f"  Current Peak File: {os.path.basename(peak_file) if peak_file else 'None'}\n"
        config_info += "\n"

        # Configuration File Info
        if hasattr(mw, 'config_manager'):
            config_info += f"Configuration File: {mw.config_manager.config_file}\n"
            last_updated = mw.config_manager.config.get('last_updated', 'Never')
            config_info += f"Last Updated: {last_updated}\n"

        self.config_text.setPlainText(config_info)

    def save_config(self):
        """Save configuration to file.

        Based on v0.9 save_config() (main_gui.py:6128-6140)
        """
        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save Configuration",
            "",
            "JSON files (*.json);;All files (*.*)"
        )

        if filename:
            if self.main_window and hasattr(self.main_window, 'config_manager'):
                if self.main_window.config_manager.export_config(filename):
                    QMessageBox.information(
                        self,
                        "Success",
                        f"Configuration saved to:\n{filename}"
                    )
                    self.update_config_display()
                else:
                    QMessageBox.critical(
                        self,
                        "Error",
                        "Failed to save configuration"
                    )
            else:
                QMessageBox.warning(
                    self,
                    "Warning",
                    "Config manager not available"
                )

    def load_config(self):
        """Load configuration from file.

        Based on v0.9 load_config() (main_gui.py:6142-6156)
        """
        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Load Configuration",
            "",
            "JSON files (*.json);;All files (*.*)"
        )

        if filename:
            if self.main_window and hasattr(self.main_window, 'config_manager'):
                if self.main_window.config_manager.import_config(filename):
                    self._apply_config_to_gui()
                    QMessageBox.information(
                        self,
                        "Success",
                        f"Configuration loaded from:\n{filename}"
                    )
                    self.update_config_display()
                else:
                    QMessageBox.critical(
                        self,
                        "Error",
                        "Failed to load configuration"
                    )
            else:
                QMessageBox.warning(
                    self,
                    "Warning",
                    "Config manager not available"
                )

    def reset_config(self):
        """Reset configuration to defaults.

        Based on v0.9 reset_config() (main_gui.py:6158-6165)
        """
        reply = QMessageBox.question(
            self,
            "Reset Configuration",
            "Reset all settings to default values?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            if self.main_window and hasattr(self.main_window, 'config_manager'):
                self.main_window.config_manager.reset_to_defaults()
                self._apply_config_to_gui()
                QMessageBox.information(
                    self,
                    "Success",
                    "Configuration reset to defaults"
                )
                self.update_config_display()
            else:
                QMessageBox.warning(
                    self,
                    "Warning",
                    "Config manager not available"
                )

    def _apply_config_to_gui(self):
        """Apply loaded config values to main window and GUI controls.

        Maps config_manager values to main_window attributes and updates
        the corresponding GUI spinboxes and checkboxes.
        """
        if not self.main_window:
            return

        mw = self.main_window
        config = mw.config_manager.config

        # Get processing parameters from config
        proc_params = config.get('processing_parameters', {})
        series_opts = config.get('series_options', {})

        # Update main window attributes from config
        mw.noise_threshold = proc_params.get('noise_threshold', 3.0)
        mw.search_window_x = proc_params.get('search_window_x', 0.08)
        mw.search_window_y = proc_params.get('search_window_y', 0.8)
        mw.fix_positions = proc_params.get('fix_positions', False)
        mw.use_parallel_processing = series_opts.get('parallel_processing', False)

        # Update GUI controls if they exist
        try:
            # Expert panel controls
            if hasattr(mw, 'expert_noise_spin'):
                mw.expert_noise_spin.setValue(mw.noise_threshold)
            if hasattr(mw, 'expert_search_x_spin'):
                mw.expert_search_x_spin.setValue(mw.search_window_x)
            if hasattr(mw, 'expert_search_y_spin'):
                mw.expert_search_y_spin.setValue(mw.search_window_y)
            if hasattr(mw, 'expert_fix_pos_check'):
                mw.expert_fix_pos_check.setChecked(mw.fix_positions)
            if hasattr(mw, 'expert_parallel_check'):
                mw.expert_parallel_check.setChecked(mw.use_parallel_processing)

            # Update integrator if available
            if hasattr(mw, 'integrator') and mw.integrator:
                mw.integrator.set_search_window(mw.search_window_x, mw.search_window_y)
                mw.integrator.set_threshold_multiplier(mw.noise_threshold)

            logger.info(f"Config applied to GUI: noise={mw.noise_threshold}, "
                       f"search_x={mw.search_window_x}, search_y={mw.search_window_y}")

        except Exception as e:
            logger.warning(f"Error applying config to GUI controls: {e}")