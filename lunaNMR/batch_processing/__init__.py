"""
Batch Processing Module for LunaNMR

This module provides batch processing capabilities for automated NMR spectrum
analysis and ML training data generation. It operates completely independently
from the GUI interface.

Features:
- Batch processing of multiple NMR spectra
- Automatic parameter optimization
- Command-line interface for easy usage
- Integration with ML training data collection
- Robust error handling and progress tracking

Components:
- BatchProcessor: Core batch processing engine
- ParameterOptimizer: Automatic parameter tuning
- CLIInterface: Command-line interface
- ConfigManager: Configuration file handling

Usage:
    from lunaNMR.batch_processing import BatchProcessor

    processor = BatchProcessor()
    processor.process_folder("/path/to/spectra")
"""

__version__ = "1.0.0"
__author__ = "lunaNMR Development Team"

# Import main components when available
try:
    from .batch_processor import BatchProcessor
    from .parameter_optimizer import ParameterOptimizer
    from .config_manager import ConfigManager
    BATCH_PROCESSING_AVAILABLE = True
except ImportError:
    BATCH_PROCESSING_AVAILABLE = False

# CLI interface imported separately to avoid dependency issues
try:
    from .cli_interface import CLIInterface
    CLI_AVAILABLE = True
except ImportError:
    CLI_AVAILABLE = False

__all__ = [
    'BatchProcessor',
    'ParameterOptimizer',
    'ConfigManager',
    'CLIInterface',
    'BATCH_PROCESSING_AVAILABLE',
    'CLI_AVAILABLE'
]