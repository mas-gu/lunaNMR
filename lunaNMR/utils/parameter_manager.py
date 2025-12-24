# ABOUTME: GUI parameter synchronization and validation with simplified/legacy mode support.
# ABOUTME: Overrides calculated defaults with GUI spinbox values (critical for centroid windows in simplified mode).
"""
Centralized Parameter Management for NMR Processing
Handles all parameter validation, conversion, and synchronization

Author: Guillaume Mas
Date: 2025
"""

#import tkinter as tk
from typing import Dict, Any, List, Optional
from lunaNMR.utils.output_manager import log_info, log_warning

class NMRParameterManager:
    """
    Centralized parameter management for NMR processing

    This class eliminates parameter coupling by providing a single source
    of truth for all processing parameters used by both single-spectrum
    and multi-spectrum workflows.
    """

    def __init__(self):
        """Initialize with default parameters"""

        # Initialize simplified parameter manager for automated fitting
        try:
            from .simplified_parameter_manager import SimplifiedParameterManager, ParameterAdapter
            self.simplified_manager = SimplifiedParameterManager()
            self.parameter_adapter = ParameterAdapter(self.simplified_manager)
            self.use_simplified_mode = False  # Can be toggled by user
        except ImportError:
            self.simplified_manager = None
            self.parameter_adapter = None
            self.use_simplified_mode = False

        # Define all default parameters with validation ranges
        self.parameter_definitions = {
            # Detection Parameters
            'search_window_x': {'default': 0.1, 'min': 0.01, 'max': 1.0, 'type': float},
            'search_window_y': {'default': 0.5, 'min': 0.01, 'max': 2.0, 'type': float},
            'noise_threshold': {'default': 0.1, 'min': 0.01, 'max': 10.0, 'type': float},

            # Fitting Parameters
            'fitting_window_x': {'default': 0.15, 'min': 0.01, 'max': 0.5, 'type': float},  # Increased from 0.05 to ensure adequate data points
            'fitting_window_y': {'default': 2, 'min': 0.01, 'max': 10.0, 'type': float},
            'min_r_squared': {'default': 0.1, 'min': 0.0, 'max': 1.0, 'type': float},  # Lowered to 0.1 for debugging/QC
            'max_iterations': {'default': 1000, 'min': 10, 'max': 1000, 'type': int},  # Increased from 100 to ensure convergence

            # Peak Detection Parameters
            'height_threshold': {'default': 0.1, 'min': 0.001, 'max': 1.0, 'type': float},
            'distance_factor': {'default': 2.0, 'min': 1.0, 'max': 300.0, 'type': float},
            'prominence_threshold': {'default': 0.05, 'min': 0.001, 'max': 0.5, 'type': float},
            'smoothing_sigma': {'default': 1.0, 'min': 0.1, 'max': 5.0, 'type': float},
            'max_peaks_fit': {'default': 50, 'min': 1, 'max': 200, 'type': int},
            'max_optimization_iterations': {'default': 50, 'min': 1, 'max': 100, 'type': int},  # Increased from 10 to 50, max raised to 100

            # Processing Options (boolean parameters)
            'use_parallel_processing': {'default': False, 'type': bool},
            'use_global_optimization': {'default': False, 'type': bool},
            'use_centroid_refinement': {'default': True, 'type': bool},

            # Multi-Peak Detection Parameters
            'multi_peak_r2_threshold': {'default': 0.7, 'min': 0.1, 'max': 1.0, 'type': float},  # R² threshold to trigger multi-peak detection
            'multi_peak_improvement_threshold': {'default': 0.1, 'min': 0.05, 'max': 0.5, 'type': float},  # Minimum improvement required from multi-peak fit
            'peak_detection_sensitivity': {'default': 1.5, 'min': 0.5, 'max': 5.0, 'type': float},  # Height threshold multiplier (lower = more sensitive)
            'overlap_detection_factor': {'default': 0.8, 'min': 0.3, 'max': 1.5, 'type': float},  # Separation factor for overlap detection
            'residual_analysis_threshold': {'default': 1.5, 'min': 0.5, 'max': 3.0, 'type': float},  # Threshold for residual peak detection

            # Advanced Parameters
            'centroid_window_x_ppm': {'default': 0.02, 'min': 0.005, 'max': 0.1, 'type': float},
            'centroid_window_y_ppm': {'default': 1.0, 'min': 0.01, 'max': 5.0, 'type': float},
            'centroid_noise_multiplier': {'default': 2.0, 'min': 1.0, 'max': 5.0, 'type': float},

            # PS2D Multi-Peak Fitting Parameters
            'use_ps2d_multi_peak': {'default': True, 'type': bool},
            'fix_linewidths': {'default': False, 'type': bool},
            'fix_positions': {'default': False, 'type': bool},
            'lw_lorentz_1h': {'default': None, 'type': (type(None), float)},
            'lw_gauss_1h': {'default': None, 'type': (type(None), float)},
            'lw_lorentz_15n': {'default': None, 'type': (type(None), float)},
            'lw_gauss_15n': {'default': None, 'type': (type(None), float)},

            # Series Integration / PS2D Linewidth Reuse Parameters
            'use_ps2d_linewidth_reuse': {'default': False, 'type': bool}  # Enable PS2D linewidth reuse from reference spectrum (C++ peakfit.cpp:586-607)
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
            'centroid_noise_multiplier': 'centroid_noise_multiplier',
            # Multi-peak detection parameters
            'multi_peak_r2_threshold': 'multi_peak_r2_threshold',
            'multi_peak_improvement_threshold': 'multi_peak_improvement_threshold',
            'peak_detection_sensitivity': 'peak_detection_sensitivity',
            'overlap_detection_factor': 'overlap_detection_factor',
            'residual_analysis_threshold': 'residual_analysis_threshold',
            # PS2D multi-peak fitting parameters (CRITICAL for PS2D activation)
            'use_ps2d_multi_peak': 'use_ps2d_multi_peak',
            'fix_linewidths': 'fix_linewidths',
            'fix_positions': 'fix_positions',
            'lw_lorentz_1h': 'lw_lorentz_1h',
            'lw_gauss_1h': 'lw_gauss_1h',
            'lw_lorentz_15n': 'lw_lorentz_15n',
            'lw_gauss_15n': 'lw_gauss_15n',
            # Series Integration / PS2D Linewidth Reuse
            'use_ps2d_linewidth_reuse': 'use_ps2d_linewidth_reuse'
        }

        updated_params = {}

        for gui_var_name, param_name in gui_variable_mapping.items():
            if hasattr(gui_object, gui_var_name):
                gui_var = getattr(gui_object, gui_var_name)
                try:
                    # Handle both Qt/tkinter widgets (.get() method) and plain attributes
                    if hasattr(gui_var, 'get') and callable(gui_var.get):
                        # Widget with .get() method (tkinter variable or similar)
                        value = gui_var.get()
                    elif hasattr(gui_var, 'value') and callable(gui_var.value):
                        # Qt widget with .value() method (QSpinBox, QDoubleSpinBox)
                        value = gui_var.value()
                    elif hasattr(gui_var, 'isChecked') and callable(gui_var.isChecked):
                        # Qt checkbox with .isChecked() method
                        value = gui_var.isChecked()
                    elif isinstance(gui_var, (int, float, bool, str, type(None))):
                        # Plain Python attribute (int, float, bool, str, None)
                        value = gui_var
                    else:
                        # Unknown type - skip
                        continue

                    # Validate the parameter
                    validated_value = self._validate_parameter(param_name, value)
                    self.current_params[param_name] = validated_value
                    updated_params[param_name] = validated_value
                except Exception as e:
                    log_warning(f"Error updating parameter {param_name}: {e}")

        log_info(f"Parameter Manager updated {len(updated_params)} parameters from GUI")
        return updated_params.copy()

    def _validate_parameter(self, param_name: str, value: Any) -> Any:
        """Validate a parameter value against its definition"""

        if param_name not in self.parameter_definitions:
            log_warning(f"Unknown parameter: {param_name}")
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
            log_warning(f"Type conversion failed for {param_name}: {value}")
            return definition['default']

        # Range validation (for numeric types)
        if definition['type'] in [int, float]:
            if 'min' in definition and validated_value < definition['min']:
                log_warning(f"{param_name} value {validated_value} below minimum {definition['min']}")
                validated_value = definition['min']
            elif 'max' in definition and validated_value > definition['max']:
                log_warning(f"{param_name} value {validated_value} above maximum {definition['max']}")
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
                'centroid_noise_multiplier': self.current_params['centroid_noise_multiplier'],
                # PS2D multi-peak fitting parameters (CRITICAL for PS2D activation)
                'use_ps2d_multi_peak': self.current_params.get('use_ps2d_multi_peak', True),
                'fix_linewidths': self.current_params.get('fix_linewidths', False),
                'fix_positions': self.current_params.get('fix_positions', False),
                'lw_lorentz_1h': self.current_params.get('lw_lorentz_1h', None),
                'lw_gauss_1h': self.current_params.get('lw_gauss_1h', None),
                'lw_lorentz_15n': self.current_params.get('lw_lorentz_15n', None),
                'lw_gauss_15n': self.current_params.get('lw_gauss_15n', None),
                # Series Integration / PS2D Linewidth Reuse
                'use_ps2d_linewidth_reuse': self.current_params.get('use_ps2d_linewidth_reuse', False),
                # Adaptive parameter optimization
                'use_adaptive_optimization': self.current_params.get('use_adaptive_optimization', True)
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

                # Check type (handle both single types and tuples for Optional types)
                expected_type = definition['type']
                if isinstance(expected_type, tuple):
                    # Multiple acceptable types (e.g., Optional[float] = (NoneType, float))
                    if not isinstance(current_value, expected_type):
                        type_names = " or ".join(t.__name__ for t in expected_type)
                        errors.append(f"{param_name}: expected {type_names}, got {type(current_value).__name__}")
                else:
                    # Single type
                    if not isinstance(current_value, expected_type):
                        errors.append(f"{param_name}: expected {expected_type.__name__}, got {type(current_value).__name__}")

                # Check range for numeric types
                if expected_type in [int, float] or (isinstance(expected_type, tuple) and float in expected_type):
                    if current_value is not None:  # Skip range check for None values
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
        log_info("All parameters reset to defaults")

    def update_simplified_parameters(self, **kwargs):
        """Update simplified parameters when in simplified mode"""
        if not self.use_simplified_mode or not self.simplified_manager:
            log_warning("Not in simplified mode or manager not available")
            return

        self.simplified_manager.update_simplified_parameters(**kwargs)

    def get_effective_parameters(self, nucleus_type='default'):
        """
        Get effective parameters based on current mode.

        Returns either simplified-derived parameters or legacy parameters.
        """
        if self.use_simplified_mode and self.parameter_adapter:
            # Use simplified parameters with legacy compatibility
            params = self.parameter_adapter.get_integrator_parameters(nucleus_type)

            # CRITICAL FIX: Detection parameters (search windows) must ALWAYS use GUI values
            # Simplified mode only affects FITTING parameters, NOT detection parameters
            params['detection_params']['search_window_x'] = self.current_params['search_window_x']
            params['detection_params']['search_window_y'] = self.current_params['search_window_y']
            params['detection_params']['noise_threshold'] = self.current_params['noise_threshold']

            # CRITICAL FIX: GUI control parameters must ALWAYS use GUI values
            # Simplified mode should NOT override fix_positions, fix_linewidths, use_parallel_processing
            params['gui_params']['fix_positions'] = self.current_params.get('fix_positions', False)
            params['gui_params']['fix_linewidths'] = self.current_params.get('fix_linewidths', False)
            params['gui_params']['use_parallel_processing'] = self.current_params.get('use_parallel_processing', False)

            # CRITICAL FIX: Centroid parameters must ALWAYS use GUI values
            # Simplified mode should NOT override user-specified centroid window sizes
            params['gui_params']['centroid_window_x_ppm'] = self.current_params.get('centroid_window_x_ppm', 0.01)
            params['gui_params']['centroid_window_y_ppm'] = self.current_params.get('centroid_window_y_ppm', 0.1)
            params['gui_params']['centroid_noise_multiplier'] = self.current_params.get('centroid_noise_multiplier', 2.0)

            # Adaptive optimization parameter
            params['gui_params']['use_adaptive_optimization'] = self.current_params.get('use_adaptive_optimization', True)

            return params
        else:
            # Use legacy parameters
            return self.get_integrator_parameters()

    def validate_simplified_parameters(self):
        """Validate simplified parameters"""
        if self.simplified_manager:
            is_valid, errors = self.simplified_manager.validate_simplified_parameters()
            return is_valid, errors
        else:
            return True, []
