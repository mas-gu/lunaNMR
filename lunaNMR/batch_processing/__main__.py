#!/usr/bin/env python3
"""
Main entry point for lunaNMR batch processing module.

This allows the batch processing system to be run as a module using:
    python -m lunaNMR.batch_processing

The module provides a complete command-line interface for batch processing
NMR spectra with automatic parameter optimization and configuration management.
"""

from .cli_interface import main

if __name__ == "__main__":
    main()