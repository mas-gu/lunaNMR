"""
Configuration module for overlap resolution in LunaNMR.

This module provides a flexible configuration system for the overlap resolution
engine with validation, presets, and JSON serialization support.

Classes:
    OverlapResolutionConfig: Main configuration container with validation

Functions:
    create_default_config: Factory for default configuration
    create_fast_config: Factory for speed-optimized configuration
    create_thorough_config: Factory for accuracy-optimized configuration

Author: LunaNMR Development Team
Version: 0.9
"""

import json
import copy
from typing import Dict, Any, List, Optional, Union
import warnings
from pathlib import Path


class OverlapResolutionConfig:
    """
    Configuration container for overlap resolution.

    Provides centralized configuration management with validation,
    dot notation access, and JSON serialization support.

    Attributes:
        config (Dict): Current configuration dictionary

    Example:
        >>> config = OverlapResolutionConfig()
        >>> config.get('jackknife.n_resamples')
        50
        >>> config.set('jackknife.n_resamples', 100)
        >>> config.to_json('my_config.json')
    """

    DEFAULT_CONFIG = {
        # Staged Fitting Configuration
        'staged_fitting': {
            'enabled': True,
            'stage1_position_lock': True,
            'stage2_position_tolerance': 0.05,  # ±5% of position
            'stage1_max_iterations': 1000,
            'stage2_max_iterations': 1500,
            'stage3_max_iterations': 2000,
            'stage1_ftol': 1e-6,
            'stage2_ftol': 1e-7,
            'stage3_ftol': 1e-8
        },

        # Jackknife Validation Configuration
        'jackknife': {
            'enabled': True,
            'n_resamples': 50,
            'delete_fraction': 'sqrt_n',  # Delete √N points, can also be float (0-1)
            'cv_threshold': 0.2,  # Coefficient of variation threshold
            'min_successful_refits': 25,
            'max_refit_iterations': 500
        },

        # Correlation Analysis Configuration
        'correlation': {
            'enabled': True,
            'amplitude_correlation_threshold': 0.7,
            'sigma_gamma_correlation_threshold': 0.85,
            'auto_merge_threshold': 0.9
        },

        # Model Selection Configuration
        'model_selection': {
            'criterion': 'BIC',  # 'AIC' or 'BIC'
            'max_peaks': 10,
            'min_peak_separation_ppm': 0.01,
            'aic_penalty_factor': 2.0,
            'try_single_peak_first': True
        },

        # Quality Thresholds
        'quality': {
            'min_r_squared': 0.7,
            'excellent_r_squared': 0.95,
            'good_r_squared': 0.85,
            'fair_r_squared': 0.70
        },

        # Fallback Behavior
        'fallback': {
            'enabled': True,
            'single_peak_on_failure': True,
            'max_retry_attempts': 2,
            'retry_with_fewer_peaks': True
        },

        # Performance Settings
        'performance': {
            'enable_parallel': False,  # Future: parallel jackknife
            'verbose': True,
            'log_level': 'INFO'  # 'DEBUG', 'INFO', 'WARNING', 'ERROR'
        }
    }

    def __init__(self, user_config: Optional[Dict] = None):
        """
        Initialize configuration with optional user overrides.

        Args:
            user_config: Dictionary of user configuration overrides

        Raises:
            ValueError: If configuration validation fails
        """
        self.config = self._deep_copy_dict(self.DEFAULT_CONFIG)
        if user_config:
            self._update_config(user_config)

        # Validate and collect warnings
        warnings_list = self._validate_config()
        for warning_msg in warnings_list:
            warnings.warn(warning_msg, UserWarning)

    def get(self, key_path: str, default: Any = None) -> Any:
        """
        Get configuration value using dot notation.

        Args:
            key_path: Dot-separated path (e.g., 'staged_fitting.enabled')
            default: Default value if key not found

        Returns:
            Configuration value or default

        Example:
            >>> config.get('jackknife.n_resamples')
            50
            >>> config.get('nonexistent.key', 'default_value')
            'default_value'
        """
        keys = key_path.split('.')
        value = self.config

        try:
            for key in keys:
                value = value[key]
            return value
        except (KeyError, TypeError):
            return default

    def set(self, key_path: str, value: Any):
        """
        Set configuration value using dot notation.

        Args:
            key_path: Dot-separated path (e.g., 'jackknife.n_resamples')
            value: Value to set

        Raises:
            ValueError: If validation fails for the new value

        Example:
            >>> config.set('jackknife.n_resamples', 100)
            >>> config.set('correlation.enabled', False)
        """
        keys = key_path.split('.')
        d = self.config

        # Navigate to parent dictionary
        for key in keys[:-1]:
            if key not in d:
                d[key] = {}
            d = d[key]

        # Set the value
        d[keys[-1]] = value

        # Validate the change
        warnings_list = self._validate_config()
        if warnings_list:
            # Validation warnings - inform user but don't fail
            for warning_msg in warnings_list:
                warnings.warn(warning_msg, UserWarning)

    def update(self, updates: Dict):
        """
        Update configuration with dictionary of changes.

        Args:
            updates: Dictionary containing configuration updates

        Example:
            >>> config.update({
            ...     'jackknife': {'n_resamples': 100},
            ...     'correlation': {'enabled': False}
            ... })
        """
        self._update_config(updates)

        # Validate after updates
        warnings_list = self._validate_config()
        for warning_msg in warnings_list:
            warnings.warn(warning_msg, UserWarning)

    def _update_config(self, user_config: Dict):
        """
        Recursively update configuration dictionary.

        Args:
            user_config: User configuration dictionary
        """
        def recursive_update(base: Dict, updates: Dict):
            """Recursively merge updates into base dictionary."""
            for key, value in updates.items():
                if key in base and isinstance(base[key], dict) and isinstance(value, dict):
                    recursive_update(base[key], value)
                else:
                    base[key] = value

        recursive_update(self.config, user_config)

    def _validate_config(self) -> List[str]:
        """
        Validate configuration and return list of warnings.

        Performs comprehensive validation of all configuration parameters
        including range checks, type validation, and logical consistency.

        Returns:
            List of warning messages (empty if valid)
        """
        warnings_list = []

        # Validate jackknife configuration
        n_resamples = self.get('jackknife.n_resamples')
        if n_resamples is not None and n_resamples <= 0:
            warnings_list.append(
                f"jackknife.n_resamples must be > 0, got {n_resamples}"
            )

        cv_threshold = self.get('jackknife.cv_threshold')
        if cv_threshold is not None and not (0 <= cv_threshold <= 1):
            warnings_list.append(
                f"jackknife.cv_threshold must be between 0 and 1, got {cv_threshold}"
            )

        min_successful = self.get('jackknife.min_successful_refits')
        if min_successful is not None and n_resamples is not None:
            if min_successful > n_resamples:
                warnings_list.append(
                    f"jackknife.min_successful_refits ({min_successful}) > "
                    f"n_resamples ({n_resamples})"
                )

        delete_fraction = self.get('jackknife.delete_fraction')
        if delete_fraction not in ['sqrt_n', None]:
            if isinstance(delete_fraction, (int, float)):
                if not (0 < delete_fraction < 1):
                    warnings_list.append(
                        f"jackknife.delete_fraction must be 'sqrt_n' or float in (0,1), "
                        f"got {delete_fraction}"
                    )

        # Validate correlation thresholds
        for param in ['amplitude_correlation_threshold',
                      'sigma_gamma_correlation_threshold',
                      'auto_merge_threshold']:
            value = self.get(f'correlation.{param}')
            if value is not None and not (0 <= value <= 1):
                warnings_list.append(
                    f"correlation.{param} must be between 0 and 1, got {value}"
                )

        # Validate model selection
        max_peaks = self.get('model_selection.max_peaks')
        if max_peaks is not None and max_peaks < 1:
            warnings_list.append(
                f"model_selection.max_peaks must be >= 1, got {max_peaks}"
            )

        criterion = self.get('model_selection.criterion')
        if criterion not in ['AIC', 'BIC', None]:
            warnings_list.append(
                f"model_selection.criterion must be 'AIC' or 'BIC', got {criterion}"
            )

        min_sep = self.get('model_selection.min_peak_separation_ppm')
        if min_sep is not None and min_sep <= 0:
            warnings_list.append(
                f"model_selection.min_peak_separation_ppm must be > 0, got {min_sep}"
            )

        # Validate quality thresholds
        r_squared_params = ['min_r_squared', 'excellent_r_squared',
                           'good_r_squared', 'fair_r_squared']
        for param in r_squared_params:
            value = self.get(f'quality.{param}')
            if value is not None and not (0 <= value <= 1):
                warnings_list.append(
                    f"quality.{param} must be between 0 and 1, got {value}"
                )

        # Check logical ordering of R² thresholds
        min_r2 = self.get('quality.min_r_squared', 0)
        fair_r2 = self.get('quality.fair_r_squared', 0)
        good_r2 = self.get('quality.good_r_squared', 0)
        excellent_r2 = self.get('quality.excellent_r_squared', 0)

        if not (min_r2 <= fair_r2 <= good_r2 <= excellent_r2):
            warnings_list.append(
                "quality.r_squared thresholds should be ordered: "
                "min <= fair <= good <= excellent"
            )

        # Validate staged fitting parameters
        for stage in [1, 2, 3]:
            max_iter = self.get(f'staged_fitting.stage{stage}_max_iterations')
            if max_iter is not None and max_iter <= 0:
                warnings_list.append(
                    f"staged_fitting.stage{stage}_max_iterations must be > 0, "
                    f"got {max_iter}"
                )

            if stage <= 3:
                ftol = self.get(f'staged_fitting.stage{stage}_ftol')
                if ftol is not None and ftol <= 0:
                    warnings_list.append(
                        f"staged_fitting.stage{stage}_ftol must be > 0, got {ftol}"
                    )

        pos_tol = self.get('staged_fitting.stage2_position_tolerance')
        if pos_tol is not None and not (0 < pos_tol < 1):
            warnings_list.append(
                f"staged_fitting.stage2_position_tolerance should be between 0 and 1, "
                f"got {pos_tol}"
            )

        # Validate fallback configuration
        max_retries = self.get('fallback.max_retry_attempts')
        if max_retries is not None and max_retries < 0:
            warnings_list.append(
                f"fallback.max_retry_attempts must be >= 0, got {max_retries}"
            )

        # Validate performance settings
        log_level = self.get('performance.log_level')
        valid_levels = ['DEBUG', 'INFO', 'WARNING', 'ERROR']
        if log_level not in valid_levels and log_level is not None:
            warnings_list.append(
                f"performance.log_level must be one of {valid_levels}, got {log_level}"
            )

        return warnings_list

    def _deep_copy_dict(self, d: Dict) -> Dict:
        """
        Create deep copy of nested dictionary.

        Args:
            d: Dictionary to copy

        Returns:
            Deep copy of input dictionary
        """
        return copy.deepcopy(d)

    def to_dict(self) -> Dict:
        """
        Export configuration as dictionary.

        Returns:
            Deep copy of configuration dictionary
        """
        return self._deep_copy_dict(self.config)

    def to_json(self, filepath: Union[str, Path], indent: int = 2):
        """
        Save configuration to JSON file.

        Args:
            filepath: Path to output JSON file
            indent: Number of spaces for indentation (default: 2)

        Example:
            >>> config.to_json('my_config.json')
        """
        filepath = Path(filepath)

        # Ensure directory exists
        filepath.parent.mkdir(parents=True, exist_ok=True)

        with open(filepath, 'w') as f:
            json.dump(self.config, f, indent=indent)

    @classmethod
    def from_json(cls, filepath: Union[str, Path]) -> 'OverlapResolutionConfig':
        """
        Load configuration from JSON file.

        Args:
            filepath: Path to JSON configuration file

        Returns:
            OverlapResolutionConfig instance

        Raises:
            FileNotFoundError: If file does not exist
            json.JSONDecodeError: If file is not valid JSON

        Example:
            >>> config = OverlapResolutionConfig.from_json('my_config.json')
        """
        filepath = Path(filepath)

        with open(filepath, 'r') as f:
            user_config = json.load(f)

        return cls(user_config=user_config)

    def __repr__(self) -> str:
        """String representation of configuration."""
        return f"OverlapResolutionConfig({self.config})"

    def __str__(self) -> str:
        """Human-readable string representation."""
        lines = ["OverlapResolutionConfig:"]
        lines.append(f"  Staged Fitting: {self.get('staged_fitting.enabled')}")
        lines.append(f"  Jackknife: {self.get('jackknife.enabled')} "
                    f"({self.get('jackknife.n_resamples')} resamples)")
        lines.append(f"  Correlation Analysis: {self.get('correlation.enabled')}")
        lines.append(f"  Model Selection: {self.get('model_selection.criterion')} "
                    f"(max {self.get('model_selection.max_peaks')} peaks)")
        lines.append(f"  Min R²: {self.get('quality.min_r_squared')}")
        return '\n'.join(lines)


# Factory functions for common configurations

def create_default_config() -> OverlapResolutionConfig:
    """
    Factory function for default configuration.

    Provides balanced settings suitable for most use cases:
    - 50 jackknife resamples
    - Moderate iteration limits
    - Standard quality thresholds

    Returns:
        OverlapResolutionConfig with default settings

    Example:
        >>> config = create_default_config()
    """
    return OverlapResolutionConfig()


def create_fast_config() -> OverlapResolutionConfig:
    """
    Configuration optimized for speed.

    Reduces computational overhead while maintaining reasonable accuracy:
    - 25 jackknife resamples (reduced from 50)
    - Lower iteration limits for faster convergence
    - Same quality thresholds

    Returns:
        OverlapResolutionConfig optimized for speed

    Example:
        >>> config = create_fast_config()
    """
    return OverlapResolutionConfig({
        'jackknife': {
            'enabled': True,
            'n_resamples': 25,
            'min_successful_refits': 15
        },
        'staged_fitting': {
            'stage1_max_iterations': 500,
            'stage2_max_iterations': 750,
            'stage3_max_iterations': 1000
        }
    })


def create_thorough_config() -> OverlapResolutionConfig:
    """
    Configuration optimized for accuracy.

    Maximizes accuracy at the cost of computational time:
    - 100 jackknife resamples (increased from 50)
    - Higher iteration limits for better convergence
    - More peaks allowed in model selection
    - Stricter convergence tolerances

    Returns:
        OverlapResolutionConfig optimized for accuracy

    Example:
        >>> config = create_thorough_config()
    """
    return OverlapResolutionConfig({
        'jackknife': {
            'enabled': True,
            'n_resamples': 100,
            'min_successful_refits': 60,
            'cv_threshold': 0.15  # Stricter uncertainty requirement
        },
        'staged_fitting': {
            'stage1_max_iterations': 2000,
            'stage2_max_iterations': 2500,
            'stage3_max_iterations': 3000,
            'stage1_ftol': 1e-7,
            'stage2_ftol': 1e-8,
            'stage3_ftol': 1e-9
        },
        'model_selection': {
            'max_peaks': 15
        },
        'quality': {
            'min_r_squared': 0.75  # Higher minimum threshold
        }
    })


def create_minimal_config() -> OverlapResolutionConfig:
    """
    Minimal configuration for quick testing.

    Disables advanced features for fastest possible execution:
    - Jackknife disabled
    - Correlation analysis disabled
    - Minimal iterations
    - Single peak preference

    Returns:
        OverlapResolutionConfig with minimal settings

    Example:
        >>> config = create_minimal_config()
    """
    return OverlapResolutionConfig({
        'jackknife': {
            'enabled': False
        },
        'correlation': {
            'enabled': False
        },
        'staged_fitting': {
            'stage1_max_iterations': 300,
            'stage2_max_iterations': 500,
            'stage3_max_iterations': 700
        },
        'model_selection': {
            'try_single_peak_first': True,
            'max_peaks': 5
        }
    })


# Module-level convenience function
def load_config(filepath: Optional[Union[str, Path]] = None,
                preset: str = 'default') -> OverlapResolutionConfig:
    """
    Load configuration from file or create preset.

    Args:
        filepath: Path to JSON config file (optional)
        preset: Preset name if no filepath given ('default', 'fast', 'thorough', 'minimal')

    Returns:
        OverlapResolutionConfig instance

    Example:
        >>> config = load_config('custom.json')
        >>> config = load_config(preset='fast')
    """
    if filepath:
        return OverlapResolutionConfig.from_json(filepath)

    preset_factories = {
        'default': create_default_config,
        'fast': create_fast_config,
        'thorough': create_thorough_config,
        'minimal': create_minimal_config
    }

    if preset not in preset_factories:
        raise ValueError(
            f"Unknown preset '{preset}'. Choose from: {list(preset_factories.keys())}"
        )

    return preset_factories[preset]()
