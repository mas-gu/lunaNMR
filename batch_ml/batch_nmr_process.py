#!/usr/bin/env python3
"""
Standalone Batch Processing Script for LunaNMR

This is a standalone script that can be run independently to process
multiple NMR spectra for ML training data generation. It provides a
simple interface to the comprehensive batch processing system.

Usage:
    python batch_nmr_process.py /path/to/spectra
    python batch_nmr_process.py /path/to/spectra --nucleus 15N1H
    python batch_nmr_process.py /path/to/spectra --optimize --preset conservative

This script is completely independent of the lunaNMR GUI interface.
"""

import sys
import os
from pathlib import Path

# Add the lunaNMR path to Python path (adjusted for new location)
script_dir = Path(__file__).parent
lunanmr_path = script_dir.parent  # Go up to lunaNMR_v0o9 folder
sys.path.insert(0, str(lunanmr_path))

try:
    from lunaNMR.batch_processing.cli_interface import CLIInterface

    def main():
        """Main entry point for standalone script."""
        print("LunaNMR Standalone Batch Processor")
        print("=" * 50)

        # Create and run CLI interface
        cli = CLIInterface()
        return cli.main()

    if __name__ == "__main__":
        sys.exit(main())

except ImportError as e:
    print("Error: Could not import lunaNMR batch processing components")
    print(f"Details: {e}")
    print("\nPlease ensure:")
    print("1. lunaNMR is properly installed")
    print("2. The batch processing module is available")
    print("3. All dependencies are installed")
    sys.exit(1)