"""
Simplified Parameter Manager for Automated NMR Fitting
=====================================================

This module implements the Priority 1 improvements by reducing parameter complexity
from 25+ parameters to 3-5 core parameters with natural parameter space and
adaptive quality thresholds.

Author: Guillaume Mas
Date: 2025
"""

import numpy as np
import warnings
from typing import Dict, Any, Tuple, Optional
from dataclasses import dataclass

@dataclass
class SimplifiedFittingParameters:
    """Core simplified parameters for automated Voigt fitting ONLY (does not affect peak detection)"""

    # Core Parameters (3-5 main parameters) - VOIGT FITTING ONLY
    sensitivity: float = 0.5          # Legacy parameter (no longer affects detection - kept for compatibility)
    window_scale: float = 1.0         # Fitting window scale factor (affects fitting windows only)
    quality_target: float = 0.85      # Target fitting quality (affects fitting R² thresholds only)

    # Advanced Parameters (optional, for fine-tuning)
    noise_estimation_method: str = 'auto'  # 'auto', 'percentile', 'robust'
    baseline_method: str = 'auto'          # 'auto', 'polynomial', 'iterative'

class SimplifiedParameterManager:
    """
    Simplified parameter management for Voigt fitting ONLY (does not affect peak detection).

    Reduces Voigt fitting complexity from 25+ to 3-5 core parameters while keeping
    peak detection parameters at their complex defaults for maximum quality.

    This class implements natural parameter space transformations and adaptive
    thresholds while maintaining full backward compatibility with existing workflows.
    """

    def __init__(self):
        """Initialize with simplified parameter structure"""

        # Core simplified parameters
        self.simplified_params = SimplifiedFittingParameters()

        # Natural parameter space mappings (log-transformed to avoid bounds violations)
        self.natural_space_mappings = {
            'sensitivity': {
                'transform': lambda x: np.clip(x, 0.01, 0.99),  # Keep in bounds
                'inverse': lambda x: x,
                'default_range': (0.1, 0.9)
            },
            'window_scale': {
                'transform': lambda x: np.exp(np.clip(np.log(np.maximum(x, 0.01)), -3, 3)),
                'inverse': lambda x: np.log(np.maximum(x, 0.01)),
                'default_range': (0.1, 10.0)
            },
            'quality_target': {
                'transform': lambda x: np.clip(x, 0.3, 0.99),
                'inverse': lambda x: x,
                'default_range': (0.5, 0.95)
            }
        }

        # Nucleus-specific adaptive parameters
        self.nucleus_parameters = {
            '1H': {
                'typical_linewidth': 0.01,    # ppm
                'snr_threshold': 5.0,
                'window_base': {'x': 0.15, 'y': 2.0}
            },
            '15N': {
                'typical_linewidth': 0.5,     # ppm
                'snr_threshold': 3.0,
                'window_base': {'x': 0.3, 'y': 4.0}
            },
            '13C': {
                'typical_linewidth': 0.6,     # ppm
                'snr_threshold': 3.0,
                'window_base': {'x': 0.25, 'y': 3.0}
            },
            'default': {
                'typical_linewidth': 0.1,     # ppm
                'snr_threshold': 4.0,
                'window_base': {'x': 0.2, 'y': 2.0}
            }
        }

        # Cache for derived parameters to avoid recalculation
        self._derived_params_cache = {}
        self._cache_valid = False

    def update_simplified_parameters(self, **kwargs) -> None:
        """Update simplified parameters and invalidate cache"""

        for param_name, value in kwargs.items():
            if hasattr(self.simplified_params, param_name):
                # Apply natural space transformation
                if param_name in self.natural_space_mappings:
                    value = self.natural_space_mappings[param_name]['transform'](value)

                setattr(self.simplified_params, param_name, value)
                print(f"✅ Updated simplified parameter {param_name} = {value}")
            else:
                warnings.warn(f"Unknown simplified parameter: {param_name}")

        # Invalidate cache when parameters change
        self._cache_valid = False

    def calculate_adaptive_quality_threshold(self, spectrum_data: np.ndarray,
                                           nucleus_type: str = 'default') -> float:
        """
        Calculate adaptive quality threshold based on SNR instead of fixed thresholds.

        Parameters:
        -----------
        spectrum_data : np.ndarray
            The spectrum intensity data
        nucleus_type : str
            Type of nucleus ('1H', '15N', '13C', 'default')

        Returns:
        --------
        float : Adaptive quality threshold (R²)
        """

        # Estimate SNR
        snr = self._estimate_snr(spectrum_data)

        # Get nucleus-specific baseline threshold
        nucleus_params = self.nucleus_parameters.get(nucleus_type, self.nucleus_parameters['default'])
        base_snr = nucleus_params['snr_threshold']

        # Adaptive threshold calculation
        # High SNR -> higher quality expected
        # Low SNR -> more tolerant thresholds
        if snr >= base_snr * 2:
            # High SNR: expect excellent fitting
            adaptive_threshold = 0.9
        elif snr >= base_snr:
            # Good SNR: expect good fitting
            adaptive_threshold = 0.8
        elif snr >= base_snr * 0.5:
            # Medium SNR: be more tolerant
            adaptive_threshold = 0.6
        else:
            # Low SNR: very tolerant
            adaptive_threshold = 0.4

        # Apply user's quality target as a scaling factor
        target_factor = self.simplified_params.quality_target / 0.85  # normalize to 0.85 default
        adaptive_threshold *= target_factor

        # Ensure reasonable bounds
        adaptive_threshold = np.clip(adaptive_threshold, 0.3, 0.99)

        return adaptive_threshold

    def _estimate_snr(self, spectrum_data: np.ndarray) -> float:
        """Estimate signal-to-noise ratio of spectrum"""

        if len(spectrum_data) == 0:
            return 1.0

        # Method 1: Robust noise estimation using MAD (Median Absolute Deviation)
        if self.simplified_params.noise_estimation_method in ['auto', 'robust']:
            # Use outer 20% of spectrum for noise estimation
            n_points = len(spectrum_data)
            edge_points = int(0.2 * n_points)

            noise_regions = np.concatenate([
                spectrum_data[:edge_points],
                spectrum_data[-edge_points:]
            ])

            # MAD-based noise estimation (more robust than std)
            median_noise = np.median(noise_regions)
            noise_mad = np.median(np.abs(noise_regions - median_noise))
            noise_level = 1.4826 * noise_mad  # Convert MAD to equivalent std

        # Method 2: Percentile-based noise estimation
        elif self.simplified_params.noise_estimation_method == 'percentile':
            noise_level = np.percentile(spectrum_data, 25) - np.percentile(spectrum_data, 10)
            noise_level = max(noise_level, np.std(spectrum_data) * 0.1)

        else:
            # Fallback: simple standard deviation
            noise_level = np.std(spectrum_data)

        # Signal estimation
        signal_level = np.max(spectrum_data) - np.median(spectrum_data)

        # Calculate SNR
        snr = signal_level / max(noise_level, 1e-10)  # Avoid division by zero

        return max(snr, 0.1)  # Minimum SNR to avoid issues

    def derive_legacy_parameters(self, nucleus_type: str = 'default',
                               spectrum_info: Optional[Dict] = None) -> Dict[str, Any]:
        """
        Derive legacy parameter set from simplified parameters for backward compatibility.

        This ensures all existing workflows continue to work unchanged.

        Parameters:
        -----------
        nucleus_type : str
            Type of nucleus for adaptive parameter calculation
        spectrum_info : dict, optional
            Additional spectrum information for adaptive calculation

        Returns:
        --------
        dict : Complete legacy parameter set
        """

        # Check cache validity
        cache_key = f"{nucleus_type}_{hash(str(spectrum_info))}"
        if self._cache_valid and cache_key in self._derived_params_cache:
            return self._derived_params_cache[cache_key].copy()

        # Get nucleus-specific parameters
        nucleus_params = self.nucleus_parameters.get(nucleus_type, self.nucleus_parameters['default'])

        # Calculate adaptive windows based on simplified parameters
        base_window = nucleus_params['window_base']
        scale = self.simplified_params.window_scale

        # Sensitivity affects detection thresholds
        sensitivity = self.simplified_params.sensitivity

        # Derive complete legacy parameter set
        # IMPORTANT: Simplified parameters ONLY affect Voigt fitting, NOT peak detection
        legacy_params = {
            # Detection Parameters - KEEP COMPLEX DEFAULTS (not affected by simplified mode)
            'search_window_x': base_window['x'],  # Use base window, not scaled
            'search_window_y': base_window['y'],  # Use base window, not scaled
            'noise_threshold': 0.1,  # Fixed default, not sensitivity-based

            # Fitting Parameters - AFFECTED BY SIMPLIFIED MODE
            'fitting_window_x': base_window['x'] * scale,
            'fitting_window_y': base_window['y'] * scale,
            'min_r_squared': self.simplified_params.quality_target,
            'max_iterations': 1000,  # Keep high for robustness

            # Peak Detection Parameters - KEEP COMPLEX DEFAULTS (not affected by simplified mode)
            'height_threshold': 0.1,  # Fixed default, not sensitivity-based
            'distance_factor': 2.0,  # Fixed default, not sensitivity-based
            'prominence_threshold': 0.05,  # Fixed default, not sensitivity-based
            'smoothing_sigma': 1.0,  # Keep constant for stability
            'max_peaks_fit': 25,  # Fixed default, not sensitivity-based
            'max_optimization_iterations': 50,

            # Processing Options - MIXED (fitting-related ones affected by simplified mode)
            'use_parallel_processing': True,
            'use_global_optimization': self.simplified_params.quality_target > 0.8,  # Fitting-related
            'use_centroid_refinement': True,

            # Multi-Peak Detection Parameters - KEEP COMPLEX DEFAULTS (not affected by simplified mode)
            'multi_peak_r2_threshold': 0.7,  # Fixed default, not quality_target-based
            'multi_peak_improvement_threshold': 0.1,  # Fixed default, not sensitivity-based
            'peak_detection_sensitivity': 1.5,  # Fixed default, not sensitivity-based
            'overlap_detection_factor': 0.8,  # Keep constant
            'residual_analysis_threshold': 1.5,  # Keep constant

            # Advanced Parameters (keep defaults for stability)
            'centroid_window_x_ppm': nucleus_params['typical_linewidth'] * 2,
            'centroid_window_y_ppm': nucleus_params['typical_linewidth'] * 10,
            'centroid_noise_multiplier': 2.0
        }

        # Apply spectrum-specific adaptations if available
        if spectrum_info and 'snr' in spectrum_info:
            snr = spectrum_info['snr']
            # Adapt thresholds based on actual SNR
            if snr < 3.0:
                legacy_params['min_r_squared'] *= 0.8  # Be more tolerant for low SNR
                legacy_params['noise_threshold'] *= 1.2
            elif snr > 10.0:
                legacy_params['min_r_squared'] = min(0.95, legacy_params['min_r_squared'] * 1.1)

        # Cache the results
        self._derived_params_cache[cache_key] = legacy_params.copy()
        self._cache_valid = True

        return legacy_params

    def get_parameter_summary(self) -> str:
        """Get human-readable summary of simplified parameters"""

        summary = ["🔧 Simplified Parameter Summary:"]
        summary.append("-" * 40)
        summary.append(f"Sensitivity: {self.simplified_params.sensitivity:.2f}")
        summary.append(f"Window Scale: {self.simplified_params.window_scale:.2f}")
        summary.append(f"Quality Target: {self.simplified_params.quality_target:.2f}")
        summary.append(f"Noise Method: {self.simplified_params.noise_estimation_method}")
        summary.append(f"Baseline Method: {self.simplified_params.baseline_method}")

        return "\n".join(summary)

    def validate_simplified_parameters(self) -> Tuple[bool, list]:
        """Validate simplified parameters"""

        errors = []

        # Check sensitivity range
        if not (0.0 <= self.simplified_params.sensitivity <= 1.0):
            errors.append("Sensitivity must be between 0.0 and 1.0")

        # Check window scale range
        if not (0.1 <= self.simplified_params.window_scale <= 10.0):
            errors.append("Window scale must be between 0.1 and 10.0")

        # Check quality target range
        if not (0.3 <= self.simplified_params.quality_target <= 1.0):
            errors.append("Quality target must be between 0.3 and 1.0")

        # Check method validity
        valid_noise_methods = ['auto', 'percentile', 'robust']
        if self.simplified_params.noise_estimation_method not in valid_noise_methods:
            errors.append(f"Noise method must be one of {valid_noise_methods}")

        valid_baseline_methods = ['auto', 'polynomial', 'iterative']
        if self.simplified_params.baseline_method not in valid_baseline_methods:
            errors.append(f"Baseline method must be one of {valid_baseline_methods}")

        return len(errors) == 0, errors

    def reset_to_defaults(self):
        """Reset to default simplified parameters"""

        self.simplified_params = SimplifiedFittingParameters()
        self._cache_valid = False
        print("✅ Simplified parameters reset to defaults")

# Backward compatibility adapter
class ParameterAdapter:
    """Adapter to maintain backward compatibility with existing parameter manager"""

    def __init__(self, simplified_manager: SimplifiedParameterManager):
        self.simplified_manager = simplified_manager
        self.legacy_manager = None  # Will be set if needed

    def get_integrator_parameters(self, nucleus_type: str = 'default') -> Dict[str, Dict[str, Any]]:
        """Get parameters formatted for core integrator (backward compatible)"""

        legacy_params = self.simplified_manager.derive_legacy_parameters(nucleus_type)

        return {
            'detection_params': {
                'search_window_x': legacy_params['search_window_x'],
                'search_window_y': legacy_params['search_window_y'],
                'noise_threshold': legacy_params['noise_threshold']
            },
            'fitting_params': {
                'fitting_window_x': legacy_params['fitting_window_x'],
                'fitting_window_y': legacy_params['fitting_window_y'],
                'min_r_squared': legacy_params['min_r_squared'],
                'max_iterations': legacy_params['max_iterations']
            },
            'gui_params': {
                'height_threshold': legacy_params['height_threshold'],
                'distance_factor': legacy_params['distance_factor'],
                'prominence_threshold': legacy_params['prominence_threshold'],
                'smoothing_sigma': legacy_params['smoothing_sigma'],
                'max_peaks_fit': legacy_params['max_peaks_fit'],
                'max_optimization_iterations': legacy_params['max_optimization_iterations'],
                'use_parallel_processing': legacy_params['use_parallel_processing'],
                'use_centroid_refinement': legacy_params['use_centroid_refinement'],
                'centroid_window_x_ppm': legacy_params['centroid_window_x_ppm'],
                'centroid_window_y_ppm': legacy_params['centroid_window_y_ppm'],
                'centroid_noise_multiplier': legacy_params['centroid_noise_multiplier']
            },
            'processing_options': {
                'use_parallel_processing': legacy_params['use_parallel_processing'],
                'use_global_optimization': legacy_params['use_global_optimization']
            }
        }