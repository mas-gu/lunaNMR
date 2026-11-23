#!/usr/bin/env python3
"""
Configuration Management Module

This module handles saving and loading of application settings,
user preferences, and processing parameters.

Classes:
- ConfigurationManager: Main configuration management
- UserPreferences: User interface preferences
- ProcessingParameters: NMR processing parameters
"""

import os
import json
from pathlib import Path
from lunaNMR.utils.output_manager import log_warning, log_error

class ConfigurationManager:
    """Main configuration management system"""

    def __init__(self, app_name="NMRPeakSeries"):
        self.app_name = app_name
        self.config_folder = self._get_config_folder()
        self.config_file = os.path.join(self.config_folder, "config.json")

        os.makedirs(self.config_folder, exist_ok=True)

        self.default_config = {
            "version": "1.1.0",
            "last_updated": None,
            "user_preferences": {
                "window_geometry": "1400x900",
                "window_state": "normal",
                "theme": "default",
                "font_size": 10,
                "auto_save": True,
                "show_tooltips": True,
                "confirm_exit": True
            },
            "processing_parameters": {
                "search_window_x": 0.2,
                "search_window_y": 3.0,
                "fitting_window_x": 0.3,
                "fitting_window_y": 8.0,
                "max_iterations": 1000,
                "use_reference_detection": True,
                "processing_mode": "full_detection"
            },
            "display_options": {
                "show_detected": True,
                "show_assigned": True,
                "show_fitted_curves": False,
                "contour_levels": 15,
                "colormap": "viridis",
                "plot_dpi": 100
            },
            "series_options": {
                "auto_process_series": True,
                "save_individual_results": True,
                "create_summary_plots": True,
                "parallel_processing": False,
                "max_workers": 4
            },
            "file_paths": {
                "last_nmr_folder": "",
                "last_peak_folder": "",
                "last_output_folder": ""
            }
        }

        self.config = self.load_config()

    def _get_config_folder(self):
        """Get platform-appropriate configuration folder"""
        home = Path.home()

        if os.name == 'nt':  # Windows
            config_folder = home / "AppData" / "Local" / self.app_name
        elif os.name == 'posix':
            if os.uname().sysname == 'Darwin':  # macOS
                config_folder = home / "Library" / "Application Support" / self.app_name
            else:  # Linux
                config_folder = home / ".config" / self.app_name
        else:
            config_folder = home / f".{self.app_name.lower()}"

        return str(config_folder)

    def load_config(self):
        """Load configuration from file"""
        try:
            if os.path.exists(self.config_file):
                with open(self.config_file, 'r', encoding='utf-8') as f:
                    loaded_config = json.load(f)

                # Merge with defaults (handles new keys)
                config = self._merge_configs(self.default_config, loaded_config)
                return config
            else:
                return self.default_config.copy()

        except Exception as e:
            log_warning(f"Failed to load config, using defaults: {e}")
            return self.default_config.copy()

    def _merge_configs(self, default, loaded):
        """Recursively merge configurations, keeping defaults for missing keys"""
        merged = default.copy()

        for key, value in loaded.items():
            if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
                merged[key] = self._merge_configs(merged[key], value)
            else:
                merged[key] = value

        return merged

    def get(self, key_path, default=None):
        """Get configuration value using dot notation (e.g., 'user_preferences.font_size')"""
        keys = key_path.split('.')
        value = self.config

        try:
            for key in keys:
                value = value[key]
            return value
        except (KeyError, TypeError):
            return default

    def reset_to_defaults(self, section=None):
        """Reset configuration to defaults"""
        if section:
            if section in self.default_config:
                self.config[section] = self.default_config[section].copy()
        else:
            self.config = self.default_config.copy()

    def export_config(self, file_path):
        """Export configuration to file"""
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(self.config, f, indent=2)
            return True
        except Exception as e:
            log_error(f"Failed to export config: {e}")
            return False

    def import_config(self, file_path):
        """Import configuration from file"""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                imported_config = json.load(f)

            self.config = self._merge_configs(self.config, imported_config)
            return True

        except Exception as e:
            log_error(f"Failed to import config: {e}")
            return False

class UserPreferences:
    """User interface preferences manager"""

    def __init__(self, config_manager):
        self.config_manager = config_manager

class ProcessingParameters:
    """NMR processing parameters manager"""

    def __init__(self, config_manager):
        self.config_manager = config_manager

    def get_detection_params(self):
        """Get peak detection parameters"""
        return {
            'search_window_x': self.config_manager.get("processing_parameters.search_window_x"),
            'search_window_y': self.config_manager.get("processing_parameters.search_window_y"),
            'use_reference_detection': self.config_manager.get("processing_parameters.use_reference_detection", True)
        }

    def get_fitting_params(self):
        """Get Voigt fitting parameters"""
        return {
            'fitting_window_x': self.config_manager.get("processing_parameters.fitting_window_x", 0.3),
            'fitting_window_y': self.config_manager.get("processing_parameters.fitting_window_y", 8.0),
            'max_iterations': self.config_manager.get("processing_parameters.max_iterations", 1000)
        }
