# ABOUTME: GUI parameter synchronization and validation.
# ABOUTME: Overrides calculated defaults with GUI spinbox values and groups them for the integrator.
"""
Centralized Parameter Management for NMR Processing
Handles all parameter validation, conversion, and synchronization
"""

from typing import Dict, Any, List
from lunaNMR.core.ps2d_config import SEARCH_WINDOW_PPM, get_ps2d_config

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
            'search_window_y': {'default': get_ps2d_config().search_window_f1, 'min': 0.001, 'max': 4.0, 'type': float},  # F1/15N or 13C

            # Processing Options (boolean parameters)
            'use_parallel_processing': {'default': False, 'type': bool},

            # Advanced Parameters
            'centroid_window_x_ppm': {'default': 0.02, 'min': 0.005, 'max': 0.2, 'type': float},
            'centroid_window_y_ppm': {'default': 1.0, 'min': 0.01, 'max': 5.0, 'type': float},

            # PS2D Multi-Peak Fitting Parameters
            'fix_linewidths': {'default': False, 'type': bool},
            'fix_positions': {'default': False, 'type': bool},

            # Series Integration / PS2D Linewidth Reuse Parameters
            'use_ps2d_linewidth_reuse': {'default': False, 'type': bool}  
        }

        self.current_params = {}
        for param_name, definition in self.parameter_definitions.items():
            self.current_params[param_name] = definition['default']

    def get_integrator_parameters(self) -> Dict[str, Dict[str, Any]]:
        """Get parameters formatted for core integrator"""

        return {
            'detection_params': {
                'search_window_x': self.current_params['search_window_x'],
                'search_window_y': self.current_params['search_window_y']
            },
            'gui_params': {
                'use_parallel_processing': self.current_params['use_parallel_processing'],
                'centroid_window_x_ppm': self.current_params['centroid_window_x_ppm'],
                'centroid_window_y_ppm': self.current_params['centroid_window_y_ppm'],
                'fix_linewidths': self.current_params.get('fix_linewidths', False),
                'fix_positions': self.current_params.get('fix_positions', False),
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
            'Processing': ['use_parallel_processing']
        }

        for category, param_names in categories.items():
            lines.append(f"\n{category}:")
            for param_name in param_names:
                if param_name in self.current_params:
                    value = self.current_params[param_name]
                    lines.append(f"  {param_name}: {value}")

        return "\n".join(lines)
