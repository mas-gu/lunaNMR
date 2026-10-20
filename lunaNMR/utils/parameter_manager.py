# ABOUTME: GUI parameter synchronization and validation.
# ABOUTME: Overrides calculated defaults with GUI spinbox values and groups them for the integrator.
"""
Centralized Parameter Management for NMR Processing
Handles all parameter validation, conversion, and synchronization
"""

from typing import Dict, Any, List
from lunaNMR.core.ps2d_config import SEARCH_WINDOW_PPM, get_ps2d_config
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

        self.parameter_definitions = {
            # Detection Parameters
            'search_window_x': {'default': SEARCH_WINDOW_PPM, 'min': 0.001, 'max': 1.0, 'type': float},  # F2/1H
            'search_window_y': {'default': get_ps2d_config().search_window_f1, 'min': 0.001, 'max': 2.0, 'type': float},  # F1/15N or 13C

            # Fitting Parameters
            'fitting_window_x': {'default': 0.15, 'min': 0.01, 'max': 0.5, 'type': float},  # Increased from 0.05 to ensure adequate data points
            'fitting_window_y': {'default': 2.0, 'min': 0.01, 'max': 10.0, 'type': float},  # Fixed: float not int
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
            'use_centroid_refinement': {'default': True, 'type': bool},

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
            'use_ps2d_linewidth_reuse': {'default': False, 'type': bool}  
        }

        self.current_params = {}
        for param_name, definition in self.parameter_definitions.items():
            self.current_params[param_name] = definition['default']

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
                'search_window_y': self.current_params['search_window_y']
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
                'use_adaptive_optimization': self.current_params.get('use_adaptive_optimization', False)
            },
            'processing_options': {
                'use_parallel_processing': self.current_params['use_parallel_processing']
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
            'Detection': ['search_window_x', 'search_window_y'],
            'Fitting': ['fitting_window_x', 'fitting_window_y', 'min_r_squared', 'max_iterations'],
            'Peak Detection': ['height_threshold', 'distance_factor', 'prominence_threshold', 'smoothing_sigma'],
            'Processing': ['use_parallel_processing', 'max_peaks_fit']
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
