#!/usr/bin/env python3
"""
Launch script for DynamiXs v2.0 GUI

This script launches the DynamiXs NMR relaxation analysis GUI using PySide6.
Simply run this script to start the graphical interface.

Usage:
    python run_dynamixs_gui.py

Requirements:
    - PySide6 (pip install PySide6)
    - numpy
    - matplotlib
    - pandas
    - scipy
    - lmfit
"""

import sys
import os
from pathlib import Path

# Add the current directory to the Python path
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))


def main():
    """Launch the DynamiXs GUI application."""
    try:
        # Import Qt modules
        from PySide6.QtWidgets import QApplication
        from PySide6.QtCore import Qt

        # Import our modules
        from gui_components import load_stylesheet
        from dynamiXs_gui import DynamiXsMainWindow

        # Create the application
        app = QApplication(sys.argv)

        # Set application metadata
        app.setApplicationName("DynamiXs")
        app.setApplicationVersion("2.0")
        app.setOrganizationName("LunaNMR")

        # Load the stylesheet
        load_stylesheet(app)

        # Create and show the main window
        window = DynamiXsMainWindow()
        window.show()

        # Start the event loop
        sys.exit(app.exec())

    except ImportError as e:
        print(f"Error importing required modules: {e}")
        print("\nPlease ensure all required packages are installed:")
        print("  pip install PySide6 numpy matplotlib pandas scipy lmfit")
        print("\nNote: DynamiXs v2.0 uses PySide6 (Qt) instead of CustomTkinter.")
        sys.exit(1)

    except Exception as e:
        print(f"Error starting GUI: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    print("Starting DynamiXs v2.0 GUI...")
    print("GUI should open in a new window.")
    main()
