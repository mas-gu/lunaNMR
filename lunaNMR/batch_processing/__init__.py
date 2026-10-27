"""
Batch Processing Module for LunaNMR

This module provides batch processing capabilities for automated NMR spectrum
analysis. It operates completely independently from the GUI interface.

Features:
- Batch processing of multiple NMR spectra
- Command-line interface for easy usage
- Robust error handling and progress tracking

Components:
- BatchProcessor: Core batch processing engine
- CLIInterface: Command-line interface
- ConfigManager: Configuration file handling

Usage:
    from lunaNMR.batch_processing import BatchProcessor

    processor = BatchProcessor()
    processor.process_folder("/path/to/spectra")
"""

__version__ = "1.1.0"
__author__ = "lunaNMR Development Team"

from .batch_processor import BatchProcessor
from .config_manager import ConfigManager
from .cli_interface import CLIInterface

__all__ = [
    'BatchProcessor',
    'ConfigManager',
    'CLIInterface',
]