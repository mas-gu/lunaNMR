#!/usr/bin/env python3
"""
Configuration Management for Batch Processing

This module handles configuration file management, parameter validation,
and preset configurations for different types of NMR batch processing tasks.

Features:
- JSON configuration file support
- Parameter validation and defaults
- Preset configurations for common nucleus types
- Configuration merging and overrides
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Any, Optional, Union
import logging

class ConfigManager:
    """
    Configuration management for batch processing operations.

    This class handles loading, validation, and management of configuration
    settings for batch NMR processing tasks.
    """

    def __init__(self):
        """Initialize the configuration manager."""
        self.logger = logging.getLogger(__name__)

        # Default configuration template
        self.default_config = {
            'processing': {
                'sn_thresholds': {
                    '15N1H': 2.5,
                    '13C1H': 2.0,
                    '1H': 2.0,
                    'default': 2.2
                },
                'expected_peaks': {
                    '15N1H': 100,
                    '13C1H': 80,
                    '1H': 60,
                    'default': 150
                },
                'quality_threshold': 0.8,
                'max_fitting_attempts': 3,
                'skip_on_error': True,
                'auto_optimize': False
            },
            'file_handling': {
                'extensions': ['.ft2', '.ft', '.pipe', '.ucsf', '.nmrpipe'],
                'recursive_search': False,
                'pattern_filters': [],
                'exclude_patterns': ['*backup*', '*temp*']
            },
            'optimization': {
                'enable_auto_optimization': True,
                'optimization_strategy': 'balanced',  # 'conservative', 'balanced', 'aggressive'
                'sn_test_range': {
                    '15N1H': [1.8, 2.0, 2.5, 3.0, 3.5],
                    '13C1H': [1.5, 1.8, 2.0, 2.5, 3.0],
                    '1H': [1.5, 1.8, 2.0, 2.5],
                    'default': [1.8, 2.0, 2.5, 3.0]
                },
                'optimal_peak_ranges': {
                    '15N1H': [30, 80],
                    '13C1H': [20, 60],
                    '1H': [15, 40],
                    'default': [20, 70]
                }
            },
            'logging': {
                'level': 'INFO',
                'detailed_logging': True,
                'progress_interval': 5,
                'log_file': None  # If None, auto-generate
            },
            'output': {
                'generate_summary': True,
                'summary_format': 'text',  # 'text', 'json', 'both'
                'ml_data_path': None,  # If None, use default ML data directory
                'backup_results': False
            }
        }

        # Preset configurations for common use cases
        self.presets = {
            '15N1H': {
                'description': 'Optimized settings for 15N-1H HSQC spectra',
                'processing': {
                    'sn_thresholds': {'default': 2.5},
                    'expected_peaks': {'default': 60}
                },
                'optimization': {
                    'optimal_peak_ranges': {'default': [30, 80]}
                }
            },
            '13C1H': {
                'description': 'Optimized settings for 13C-1H HMQC/HSQC spectra',
                'processing': {
                    'sn_thresholds': {'default': 2.0},
                    'expected_peaks': {'default': 40}
                },
                'optimization': {
                    'optimal_peak_ranges': {'default': [20, 60]}
                }
            },
            '1H': {
                'description': 'Optimized settings for 1H spectra',
                'processing': {
                    'sn_thresholds': {'default': 2.0},
                    'expected_peaks': {'default': 30}
                },
                'optimization': {
                    'optimal_peak_ranges': {'default': [15, 40]}
                }
            },
            'conservative': {
                'description': 'Conservative settings for high-quality data only',
                'processing': {
                    'sn_thresholds': {
                        '15N1H': 3.0,
                        '13C1H': 2.5,
                        '1H': 2.5,
                        'default': 2.8
                    },
                    'expected_peaks': {
                        '15N1H': 30,
                        '13C1H': 25,
                        '1H': 20,
                        'default': 20
                    },
                    'quality_threshold': 0.9
                },
                'optimization': {
                    'optimization_strategy': 'conservative'
                }
            },
            'aggressive': {
                'description': 'Aggressive settings for maximum data collection',
                'processing': {
                    'sn_thresholds': {
                        '15N1H': 1.8,
                        '13C1H': 1.5,
                        '1H': 1.5,
                        'default': 1.6
                    },
                    'expected_peaks': {
                        '15N1H': 200,
                        '13C1H': 150,
                        '1H': 100,
                        'default': 300
                    },
                    'quality_threshold': 0.6
                },
                'optimization': {
                    'optimization_strategy': 'aggressive'
                }
            },
            'high_throughput': {
                'description': 'Settings optimized for processing speed',
                'processing': {
                    'max_fitting_attempts': 1,
                    'skip_on_error': True
                },
                'optimization': {
                    'enable_auto_optimization': False
                },
                'logging': {
                    'progress_interval': 10,
                    'detailed_logging': False
                }
            }
        }

    def load_config(self, config_path: Union[str, Path]) -> Dict[str, Any]:
        """
        Load configuration from a JSON file.

        Args:
            config_path: Path to the configuration file

        Returns:
            Dictionary with loaded configuration

        Raises:
            FileNotFoundError: If config file doesn't exist
            ValueError: If config file is invalid JSON or has invalid parameters
        """
        config_path = Path(config_path)

        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        try:
            with open(config_path, 'r') as f:
                user_config = json.load(f)

            self.logger.info(f"Loaded configuration from: {config_path}")

            # Merge with defaults
            merged_config = self._merge_configs(self.default_config, user_config)

            # Validate configuration
            self._validate_config(merged_config)

            return merged_config

        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in configuration file {config_path}: {e}")
        except Exception as e:
            raise ValueError(f"Error loading configuration from {config_path}: {e}")

    def save_config(self, config: Dict[str, Any], config_path: Union[str, Path]):
        """
        Save configuration to a JSON file.

        Args:
            config: Configuration dictionary to save
            config_path: Path where to save the configuration

        Raises:
            ValueError: If configuration is invalid
        """
        # Validate before saving
        self._validate_config(config)

        config_path = Path(config_path)

        try:
            # Ensure directory exists
            config_path.parent.mkdir(parents=True, exist_ok=True)

            # Save with pretty formatting
            with open(config_path, 'w') as f:
                json.dump(config, f, indent=2, sort_keys=True)

            self.logger.info(f"Configuration saved to: {config_path}")

        except Exception as e:
            raise ValueError(f"Error saving configuration to {config_path}: {e}")

    def get_preset_config(self, preset_name: str) -> Dict[str, Any]:
        """
        Get a preset configuration.

        Args:
            preset_name: Name of the preset ('15N1H', '13C1H', etc.)

        Returns:
            Dictionary with preset configuration

        Raises:
            ValueError: If preset doesn't exist
        """
        if preset_name not in self.presets:
            available = ', '.join(self.presets.keys())
            raise ValueError(f"Unknown preset '{preset_name}'. Available presets: {available}")

        preset = self.presets[preset_name]

        # Merge preset with defaults
        config = self._merge_configs(self.default_config, preset)

        self.logger.info(f"Using preset configuration: {preset_name}")
        self.logger.info(f"Description: {preset.get('description', 'No description')}")

        return config

    def create_preset_config_files(self, output_dir: Union[str, Path]):
        """
        Create configuration files for all available presets.

        Args:
            output_dir: Directory where to create the configuration files
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        created_files = []

        for preset_name in self.presets.keys():
            config = self.get_preset_config(preset_name)
            config_file = output_dir / f"batch_config_{preset_name}.json"

            self.save_config(config, config_file)
            created_files.append(config_file)

        self.logger.info(f"Created {len(created_files)} preset configuration files in {output_dir}")
        return created_files

    def _merge_configs(self, base_config: Dict[str, Any], override_config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Recursively merge two configuration dictionaries.

        Args:
            base_config: Base configuration dictionary
            override_config: Configuration to merge/override with

        Returns:
            Merged configuration dictionary
        """
        merged = base_config.copy()

        for key, value in override_config.items():
            if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
                # Recursively merge nested dictionaries
                merged[key] = self._merge_configs(merged[key], value)
            else:
                # Override value
                merged[key] = value

        return merged

    def _validate_config(self, config: Dict[str, Any]):
        """
        Validate configuration parameters.

        Args:
            config: Configuration dictionary to validate

        Raises:
            ValueError: If configuration contains invalid parameters
        """
        errors = []

        # Validate processing parameters
        if 'processing' in config:
            proc_config = config['processing']

            # Validate S/N thresholds
            if 'sn_thresholds' in proc_config:
                for nucleus, threshold in proc_config['sn_thresholds'].items():
                    if not isinstance(threshold, (int, float)) or threshold <= 0:
                        errors.append(f"Invalid S/N threshold for {nucleus}: {threshold}")

            # Validate expected peaks
            if 'expected_peaks' in proc_config:
                for nucleus, peaks in proc_config['expected_peaks'].items():
                    if not isinstance(peaks, int) or peaks <= 0:
                        errors.append(f"Invalid expected peaks for {nucleus}: {peaks}")

            # Validate quality threshold
            if 'quality_threshold' in proc_config:
                qt = proc_config['quality_threshold']
                if not isinstance(qt, (int, float)) or not (0 < qt <= 1):
                    errors.append(f"Quality threshold must be between 0 and 1: {qt}")

        # Validate file handling parameters
        if 'file_handling' in config:
            fh_config = config['file_handling']

            if 'extensions' in fh_config:
                if not isinstance(fh_config['extensions'], list):
                    errors.append("File extensions must be a list")
                else:
                    for ext in fh_config['extensions']:
                        if not isinstance(ext, str) or not ext.startswith('.'):
                            errors.append(f"Invalid file extension: {ext}")

        # Validate optimization parameters
        if 'optimization' in config:
            opt_config = config['optimization']

            if 'optimization_strategy' in opt_config:
                valid_strategies = ['conservative', 'balanced', 'aggressive']
                if opt_config['optimization_strategy'] not in valid_strategies:
                    errors.append(f"Invalid optimization strategy. Must be one of: {valid_strategies}")

        # Validate logging parameters
        if 'logging' in config:
            log_config = config['logging']

            if 'level' in log_config:
                valid_levels = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
                if log_config['level'] not in valid_levels:
                    errors.append(f"Invalid log level. Must be one of: {valid_levels}")

        if errors:
            error_msg = "Configuration validation failed:\n" + "\n".join(f"  - {error}" for error in errors)
            raise ValueError(error_msg)

    def get_available_presets(self) -> Dict[str, str]:
        """
        Get available preset configurations with descriptions.

        Returns:
            Dictionary mapping preset names to descriptions
        """
        return {
            name: preset.get('description', 'No description')
            for name, preset in self.presets.items()
        }

    def create_example_config(self, output_path: Union[str, Path]):
        """
        Create an example configuration file with comments.

        Args:
            output_path: Path where to create the example configuration
        """
        output_path = Path(output_path)

        # Create a comprehensive example config
        example_config = {
            "_comment": "Example configuration for lunaNMR batch processing",
            "_description": {
                "processing": "Core processing parameters",
                "file_handling": "File discovery and filtering settings",
                "optimization": "Automatic parameter optimization settings",
                "logging": "Logging and progress reporting settings",
                "output": "Output and result handling settings"
            },
            **self.default_config
        }

        # Add comments to sections
        example_config["processing"]["_comment"] = "S/N thresholds and peak count expectations by nucleus type"
        example_config["optimization"]["_comment"] = "Parameters for automatic optimization of processing settings"
        example_config["logging"]["_comment"] = "Logging level and progress reporting settings"

        with open(output_path, 'w') as f:
            json.dump(example_config, f, indent=2, sort_keys=True)

        self.logger.info(f"Example configuration created: {output_path}")

if __name__ == "__main__":
    # Example usage
    cm = ConfigManager()

    # Create example configurations
    print("Available presets:")
    for name, desc in cm.get_available_presets().items():
        print(f"  {name}: {desc}")

    # Create preset configuration files
    cm.create_preset_config_files("config_presets")

    # Create example configuration
    cm.create_example_config("example_config.json")

    print("\nConfiguration files created in current directory")