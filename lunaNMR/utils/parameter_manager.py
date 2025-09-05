"""
Centralized Parameter Management for NMR Processing
Handles all parameter validation, conversion, and synchronization

Author: Guillaume Mas
Date: 2025
"""

import tkinter as tk
from typing import Dict, Any, List, Optional

class NMRParameterManager:
    """
    Centralized parameter management for NMR processing

    This class eliminates parameter coupling by providing a single source
    of truth for all processing parameters used by both single-spectrum
    and multi-spectrum workflows.
    """

    def __init__(self):
        """Initialize with default parameters"""

        # Define all default parameters with validation ranges
        self.parameter_definitions = {
            # Detection Parameters
            'search_window_x': {'default': 0.1, 'min': 0.01, 'max': 1.0, 'type': float},
            'search_window_y': {'default': 0.5, 'min': 0.01, 'max': 2.0, 'type': float},
            'noise_threshold': {'default': 0.1, 'min': 0.01, 'max': 10.0, 'type': float},

            # Fitting Parameters
            'fitting_window_x': {'default': 0.15, 'min': 0.01, 'max': 0.5, 'type': float},  # Increased from 0.05 to ensure adequate data points
            'fitting_window_y': {'default': 2, 'min': 0.01, 'max': 10.0, 'type': float},
            'min_r_squared': {'default': 0.7, 'min': 0.0, 'max': 1.0, 'type': float},  # Lowered from 0.85 to be less restrictive
            'max_iterations': {'default': 1000, 'min': 10, 'max': 1000, 'type': int},  # Increased from 100 to ensure convergence

            # Peak Detection Parameters
            'height_threshold': {'default': 0.1, 'min': 0.01, 'max': 1.0, 'type': float},
            'distance_factor': {'default': 2.0, 'min': 1.0, 'max': 100.0, 'type': float},
            'prominence_threshold': {'default': 0.05, 'min': 0.01, 'max': 0.5, 'type': float},
            'smoothing_sigma': {'default': 1.0, 'min': 0.1, 'max': 5.0, 'type': float},
            'max_peaks_fit': {'default': 50, 'min': 1, 'max': 200, 'type': int},
            'max_optimization_iterations': {'default': 50, 'min': 1, 'max': 100, 'type': int},  # Increased from 10 to 50, max raised to 100

            # Processing Options (boolean parameters)
            'use_parallel_processing': {'default': True, 'type': bool},
            'use_global_optimization': {'default': False, 'type': bool},
            'use_centroid_refinement': {'default': True, 'type': bool},

            # Advanced Parameters
            'centroid_window_x_ppm': {'default': 0.02, 'min': 0.005, 'max': 0.1, 'type': float},
            'centroid_window_y_ppm': {'default': 1.0, 'min': 0.01, 'max': 5.0, 'type': float},
            'centroid_noise_multiplier': {'default': 2.0, 'min': 1.0, 'max': 5.0, 'type': float}
        }

        # Initialize current parameters with defaults
        self.current_params = {}
        for param_name, definition in self.parameter_definitions.items():
            self.current_params[param_name] = definition['default']

    def update_from_gui_variables(self, gui_object) -> Dict[str, Any]:
        """
        Update parameters from GUI tkinter variables

        Args:
            gui_object: GUI object containing tkinter variables

        Returns:
            Dictionary of updated parameters
        """
        # Define mapping from GUI variable names to our parameter names
        gui_variable_mapping = {
            'search_window_x': 'search_window_x',
            'search_window_y': 'search_window_y',
            'noise_threshold': 'noise_threshold',
            'fitting_window_x': 'fitting_window_x',
            'fitting_window_y': 'fitting_window_y',
            'min_r_squared': 'min_r_squared',
            'max_iterations': 'max_iterations',
            'peak_height_threshold': 'height_threshold',
            'peak_distance_factor': 'distance_factor',
            'peak_prominence_threshold': 'prominence_threshold',
            'smoothing_sigma': 'smoothing_sigma',
            'max_peaks_fit': 'max_peaks_fit',
            'max_optimization_iterations': 'max_optimization_iterations',
            'use_parallel_processing': 'use_parallel_processing',
            'use_global_optimization': 'use_global_optimization',
            'use_centroid_refinement': 'use_centroid_refinement',
            'centroid_window_x_ppm': 'centroid_window_x_ppm',
            'centroid_window_y_ppm': 'centroid_window_y_ppm',
            'centroid_noise_multiplier': 'centroid_noise_multiplier'
        }

        updated_params = {}

        for gui_var_name, param_name in gui_variable_mapping.items():
            if hasattr(gui_object, gui_var_name):
                gui_var = getattr(gui_object, gui_var_name)
                if hasattr(gui_var, 'get'):
                    try:
                        value = gui_var.get()
                        # Validate the parameter
                        validated_value = self._validate_parameter(param_name, value)
                        self.current_params[param_name] = validated_value
                        updated_params[param_name] = validated_value
                    except Exception as e:
                        print(f"⚠️ Error updating parameter {param_name}: {e}")

        print(f"✅ Parameter Manager updated {len(updated_params)} parameters from GUI")
        return updated_params.copy()

    def _validate_parameter(self, param_name: str, value: Any) -> Any:
        """Validate a parameter value against its definition"""

        if param_name not in self.parameter_definitions:
            print(f"⚠️ Unknown parameter: {param_name}")
            return value

        definition = self.parameter_definitions[param_name]

        # Type conversion
        try:
            if definition['type'] == bool:
                validated_value = bool(value)
            elif definition['type'] == int:
                validated_value = int(value)
            elif definition['type'] == float:
                validated_value = float(value)
            else:
                validated_value = value
        except (ValueError, TypeError):
            print(f"⚠️ Type conversion failed for {param_name}: {value}")
            return definition['default']

        # Range validation (for numeric types)
        if definition['type'] in [int, float]:
            if 'min' in definition and validated_value < definition['min']:
                print(f"⚠️ {param_name} value {validated_value} below minimum {definition['min']}")
                validated_value = definition['min']
            elif 'max' in definition and validated_value > definition['max']:
                print(f"⚠️ {param_name} value {validated_value} above maximum {definition['max']}")
                validated_value = definition['max']

        return validated_value

    def get_integrator_parameters(self) -> Dict[str, Dict[str, Any]]:
        """Get parameters formatted for core integrator"""

        return {
            'detection_params': {
                'search_window_x': self.current_params['search_window_x'],
                'search_window_y': self.current_params['search_window_y'],
                'noise_threshold': self.current_params['noise_threshold']
            },
            'fitting_params': {
                'fitting_window_x': self.current_params['fitting_window_x'],
                'fitting_window_y': self.current_params['fitting_window_y'],
                'min_r_squared': self.current_params['min_r_squared'],
                'max_iterations': self.current_params['max_iterations']
            },
            'gui_params': {
                'height_threshold': self.current_params['height_threshold'],
                'distance_factor': self.current_params['distance_factor'],
                'prominence_threshold': self.current_params['prominence_threshold'],
                'smoothing_sigma': self.current_params['smoothing_sigma'],
                'max_peaks_fit': self.current_params['max_peaks_fit'],
                'max_optimization_iterations': self.current_params['max_optimization_iterations'],
                'use_parallel_processing': self.current_params['use_parallel_processing'],
                'use_centroid_refinement': self.current_params['use_centroid_refinement'],
                'centroid_window_x_ppm': self.current_params['centroid_window_x_ppm'],
                'centroid_window_y_ppm': self.current_params['centroid_window_y_ppm'],
                'centroid_noise_multiplier': self.current_params['centroid_noise_multiplier']
            },
            'processing_options': {
                'use_parallel_processing': self.current_params['use_parallel_processing'],
                'use_global_optimization': self.current_params['use_global_optimization']
            }
        }

    def validate_all_parameters(self) -> List[str]:
        """Validate all current parameters and return list of errors"""

        errors = []

        for param_name, current_value in self.current_params.items():
            if param_name in self.parameter_definitions:
                definition = self.parameter_definitions[param_name]

                # Check type
                expected_type = definition['type']
                if not isinstance(current_value, expected_type):
                    errors.append(f"{param_name}: expected {expected_type.__name__}, got {type(current_value).__name__}")

                # Check range for numeric types
                if expected_type in [int, float]:
                    if 'min' in definition and current_value < definition['min']:
                        errors.append(f"{param_name}: value {current_value} below minimum {definition['min']}")
                    if 'max' in definition and current_value > definition['max']:
                        errors.append(f"{param_name}: value {current_value} above maximum {definition['max']}")

        return errors

    def get_parameter_summary(self) -> str:
        """Get human-readable summary of current parameters"""

        lines = ["Parameter Summary:"]
        lines.append("-" * 40)

        categories = {
            'Detection': ['search_window_x', 'search_window_y', 'noise_threshold'],
            'Fitting': ['fitting_window_x', 'fitting_window_y', 'min_r_squared', 'max_iterations'],
            'Peak Detection': ['height_threshold', 'distance_factor', 'prominence_threshold', 'smoothing_sigma'],
            'Processing': ['use_parallel_processing', 'use_global_optimization', 'max_peaks_fit']
        }

        for category, param_names in categories.items():
            lines.append(f"\n{category}:")
            for param_name in param_names:
                if param_name in self.current_params:
                    value = self.current_params[param_name]
                    lines.append(f"  {param_name}: {value}")

        return "\n".join(lines)

    def reset_to_defaults(self):
        """Reset all parameters to default values"""
        for param_name, definition in self.parameter_definitions.items():
            self.current_params[param_name] = definition['default']
        print("✅ All parameters reset to defaults")

    def export_parameters(self, filename: str):
        """Export current parameters to a file"""
        try:
            import json
            with open(filename, 'w') as f:
                json.dump(self.current_params, f, indent=2)
            print(f"✅ Parameters exported to {filename}")
        except Exception as e:
            print(f"❌ Failed to export parameters: {e}")

    def import_parameters(self, filename: str):
        """Import parameters from a file"""
        try:
            import json
            with open(filename, 'r') as f:
                imported_params = json.load(f)

            # Validate and update parameters
            for param_name, value in imported_params.items():
                if param_name in self.parameter_definitions:
                    validated_value = self._validate_parameter(param_name, value)
                    self.current_params[param_name] = validated_value

            print(f"✅ Parameters imported from {filename}")
        except Exception as e:
            print(f"❌ Failed to import parameters: {e}")
