#!/usr/bin/env python3
"""
Configuration Management Module

This module handles saving and loading of application settings,
user preferences, and processing parameters.

Classes:
- ConfigurationManager: Main configuration management
- UserPreferences: User interface preferences
- ProcessingParameters: NMR processing parameters

Author: Guillaume Mas
Date: 2025
"""

import os
import json
from pathlib import Path
from datetime import datetime
#import tkinter as tk
from lunaNMR.utils.output_manager import log_info, log_warning, log_error

class ConfigurationManager:
    """Main configuration management system"""

    def __init__(self, app_name="NMRPeakSeries"):
        self.app_name = app_name
        self.config_folder = self._get_config_folder()
        self.config_file = os.path.join(self.config_folder, "config.json")
        self.recent_files_file = os.path.join(self.config_folder, "recent_files.json")

        # Ensure config folder exists
        os.makedirs(self.config_folder, exist_ok=True)

        # Default configuration
        self.default_config = {
            "version": "1.0.0",
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
                "noise_threshold": 3.0,
                "search_window_x": 0.2,
                "search_window_y": 3.0,
                "fitting_window_x": 0.3,
                "fitting_window_y": 8.0,
                "min_r_squared": 0.7,
                "max_iterations": 1000,
                "use_reference_detection": True,
                "processing_mode": "full_detection",
                # Multi-peak detection parameters
                "multi_peak_r2_threshold": 0.7,
                "multi_peak_improvement_threshold": 0.1,
                "peak_detection_sensitivity": 1.5,
                "overlap_detection_factor": 0.8,
                "residual_analysis_threshold": 1.5
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
            },
            "overlap_resolution": {
                "enabled": False,  # Default OFF for backward compatibility
                "auto_detect": True,
                "detection_threshold": 0.5,  # ppm separation threshold
                "staged_fitting": True,
                "jackknife_validation": True,
                "correlation_analysis": True,
                "n_resamples": 50,
                "max_peaks": 10,
                "model_selection_criterion": "BIC",
                "cv_threshold": 0.2,
                "amplitude_correlation_threshold": 0.7
            }
        }

        # Load existing configuration
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

    def save_config(self):
        """Save configuration to file"""
        try:
            self.config["last_updated"] = datetime.now().isoformat()

            with open(self.config_file, 'w', encoding='utf-8') as f:
                json.dump(self.config, f, indent=2)

            return True

        except Exception as e:
            log_warning(f"Failed to save config: {e}")
            return False

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

    def set(self, key_path, value):
        """Set configuration value using dot notation"""
        keys = key_path.split('.')
        config_section = self.config

        # Navigate to the parent of the target key
        for key in keys[:-1]:
            if key not in config_section:
                config_section[key] = {}
            config_section = config_section[key]

        # Set the value
        config_section[keys[-1]] = value

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

            # Merge with current config
            self.config = self._merge_configs(self.config, imported_config)
            return True

        except Exception as e:
            log_error(f"Failed to import config: {e}")
            return False

    # Recent Files Management
    def load_recent_files(self):
        """Load recent files list"""
        try:
            if os.path.exists(self.recent_files_file):
                with open(self.recent_files_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            else:
                return {"nmr_files": [], "peak_files": []}
        except Exception:
            return {"nmr_files": [], "peak_files": []}

    def save_recent_files(self, recent_files):
        """Save recent files list"""
        try:
            with open(self.recent_files_file, 'w', encoding='utf-8') as f:
                json.dump(recent_files, f, indent=2)
            return True
        except Exception:
            return False

    def add_recent_file(self, file_path, file_type):
        """Add file to recent files list"""
        recent = self.load_recent_files()

        key = f"{file_type}_files"
        if key not in recent:
            recent[key] = []

        file_path = os.path.abspath(file_path)

        # Remove if already exists
        if file_path in recent[key]:
            recent[key].remove(file_path)

        # Add to beginning
        recent[key].insert(0, file_path)

        # Limit size
        max_recent = 10
        if len(recent[key]) > max_recent:
            recent[key] = recent[key][:max_recent]

        self.save_recent_files(recent)

    def load_overlap_config(self):
        """
        Load overlap resolution configuration

        Returns:
            Dictionary with overlap resolution settings
        """
        return self.config.get('overlap_resolution', {
            'enabled': False,
            'auto_detect': True,
            'detection_threshold': 0.5,
            'staged_fitting': True,
            'jackknife_validation': True,
            'correlation_analysis': True,
            'n_resamples': 50,
            'max_peaks': 10,
            'model_selection_criterion': 'BIC',
            'cv_threshold': 0.2,
            'amplitude_correlation_threshold': 0.7
        })

    def enable_overlap_resolution(self, enable=True):
        """
        Quick enable/disable overlap resolution

        Args:
            enable: True to enable, False to disable
        """
        if 'overlap_resolution' not in self.config:
            self.config['overlap_resolution'] = self.load_overlap_config()

        self.config['overlap_resolution']['enabled'] = enable
        self.save_config()
        status = 'enabled' if enable else 'disabled'
        log_info(f"Overlap resolution {status}")

class UserPreferences:
    """User interface preferences manager"""

    def __init__(self, config_manager):
        self.config_manager = config_manager

    def get_window_geometry(self):
        """Get saved window geometry"""
        return self.config_manager.get("user_preferences.window_geometry", "1400x900")

    def set_window_geometry(self, geometry):
        """Save window geometry"""
        self.config_manager.set("user_preferences.window_geometry", geometry)

    def get_window_state(self):
        """Get saved window state"""
        return self.config_manager.get("user_preferences.window_state", "normal")

    def set_window_state(self, state):
        """Save window state"""
        self.config_manager.set("user_preferences.window_state", state)


class ProcessingParameters:
    """NMR processing parameters manager"""

    def __init__(self, config_manager):
        self.config_manager = config_manager

    def get_detection_params(self):
        """Get peak detection parameters"""
        return {
            'noise_threshold': self.config_manager.get("processing_parameters.noise_threshold", 3.0),
            'search_window_x': self.config_manager.get("processing_parameters.search_window_x", 0.2),
            'search_window_y': self.config_manager.get("processing_parameters.search_window_y", 3.0),
            'use_reference_detection': self.config_manager.get("processing_parameters.use_reference_detection", True)
        }

    def set_detection_params(self, **params):
        """Set peak detection parameters"""
        for key, value in params.items():
            self.config_manager.set(f"processing_parameters.{key}", value)

    def get_fitting_params(self):
        """Get Voigt fitting parameters"""
        return {
            'fitting_window_x': self.config_manager.get("processing_parameters.fitting_window_x", 0.3),
            'fitting_window_y': self.config_manager.get("processing_parameters.fitting_window_y", 8.0),
            'min_r_squared': self.config_manager.get("processing_parameters.min_r_squared", 0.7),
            'max_iterations': self.config_manager.get("processing_parameters.max_iterations", 1000),
            # Multi-peak detection parameters
            'multi_peak_r2_threshold': self.config_manager.get("processing_parameters.multi_peak_r2_threshold", 0.7),
            'multi_peak_improvement_threshold': self.config_manager.get("processing_parameters.multi_peak_improvement_threshold", 0.1),
            'peak_detection_sensitivity': self.config_manager.get("processing_parameters.peak_detection_sensitivity", 1.5),
            'overlap_detection_factor': self.config_manager.get("processing_parameters.overlap_detection_factor", 0.8),
            'residual_analysis_threshold': self.config_manager.get("processing_parameters.residual_analysis_threshold", 1.5)
        }

    def set_fitting_params(self, **params):
        """Set Voigt fitting parameters"""
        for key, value in params.items():
            self.config_manager.set(f"processing_parameters.{key}", value)

    def get_processing_mode(self):
        """Get processing mode"""
        return self.config_manager.get("processing_parameters.processing_mode", "full_detection")

    def set_processing_mode(self, mode):
        """Set processing mode"""
        self.config_manager.set("processing_parameters.processing_mode", mode)

    def get_display_options(self):
        """Get display options"""
        return {
            'show_detected': self.config_manager.get("display_options.show_detected", True),
            'show_assigned': self.config_manager.get("display_options.show_assigned", True),
            'show_fitted_curves': self.config_manager.get("display_options.show_fitted_curves", False),
            'contour_levels': self.config_manager.get("display_options.contour_levels", 15),
            'colormap': self.config_manager.get("display_options.colormap", "viridis"),
            'plot_dpi': self.config_manager.get("display_options.plot_dpi", 100)
        }

    def set_display_options(self, **options):
        """Set display options"""
        for key, value in options.items():
            self.config_manager.set(f"display_options.{key}", value)

    def get_series_options(self):
        """Get series processing options"""
        return {
            'auto_process_series': self.config_manager.get("series_options.auto_process_series", True),
            'save_individual_results': self.config_manager.get("series_options.save_individual_results", True),
            'create_summary_plots': self.config_manager.get("series_options.create_summary_plots", True),
            'parallel_processing': self.config_manager.get("series_options.parallel_processing", False),
            'max_workers': self.config_manager.get("series_options.max_workers", 4)
        }

    def set_series_options(self, **options):
        """Set series processing options"""
        for key, value in options.items():
            self.config_manager.set(f"series_options.{key}", value)
