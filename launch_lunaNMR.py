#!/usr/bin/env python3
# ABOUTME: Starts the LunaNMR Qt application: sets up paths and the matplotlib backend,
# ABOUTME: checks dependencies, and reports start-up failures in a tkinter error dialog.
"""
LunaNMR Suite - Launch Script

Launches LunaNMR; the optional modules (DynamiXs, Kd, 1D integration) open from its
Modules menu. This script handles environment setup and graceful error handling.
"""

import logging
import sys
import os

# CRITICAL: Set matplotlib backend BEFORE any matplotlib imports anywhere in the application
# This ensures NavigationToolbar2Qt renders properly on all systems
import matplotlib

matplotlib.use("QtAgg", force=True)

# Keep text EDITABLE in exported vector figures (so they can be tweaked in Illustrator
# etc.): SVG keeps labels as <text> rather than outlined paths; PDF embeds TrueType so
# glyphs stay selectable. Applies to every figure export app-wide, including the
# matplotlib toolbar's save button in each viewer.
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42

import tkinter as tk
from tkinter import messagebox


def check_dependencies():
    """Check if all required dependencies are available"""
    required_modules = ["pandas", "numpy", "matplotlib", "scipy", "nmrglue"]

    missing_modules = []
    for module in required_modules:
        try:
            __import__(module)
        except ImportError:
            missing_modules.append(module)

    return missing_modules


def setup_paths():
    """Setup Python paths for module imports"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)

    sys.path.insert(0, current_dir)
    sys.path.insert(0, parent_dir)

    return current_dir, parent_dir


def launch_lunaNMR():
    """Launch LunaNMR application (Qt/PySide6 version)"""
    try:
        from lunaNMR.gui.main_window import main as gui_main

        gui_main()
        return True
    except ImportError as e:
        error_msg = f"Failed to import LunaNMR modules:\n{str(e)}\n\nPlease ensure the lunaNMR package is properly installed."
        print(f"Error: {error_msg}")
        messagebox.showerror("Import Error", error_msg)
        return False


def main():
    """Main launcher function - directly launches LunaNMR"""
    print("LunaNMR Suite v1.0 - Starting...")

    # The terminal the app is started from shows lunaNMR's INFO records
    logging.basicConfig(level=logging.INFO)

    current_dir, parent_dir = setup_paths()

    missing = check_dependencies()
    if missing:
        error_msg = f"Missing required modules: {', '.join(missing)}\n\nPlease install with:\npip install {' '.join(missing)}"

        try:
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror("Missing Dependencies", error_msg)
            root.destroy()
        except:
            pass

        return False

    try:
        print("Launching LunaNMR...")
        success = launch_lunaNMR()

        if success:
            print("Application closed normally")
        else:
            print("Application closed with errors")

        return success

    except Exception as e:
        error_msg = f"Application error:\n{str(e)}"
        print(f"Error: {error_msg}")

        try:
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror("Application Error", error_msg)
            root.destroy()
        except:
            pass

        import traceback

        traceback.print_exc()
        return False


if __name__ == "__main__":
    try:
        success = main()
        print(f"Launch completed with {'success' if success else 'errors'}")
    except KeyboardInterrupt:
        print("\nLaunch interrupted by user")
        success = True  # User interrupt is not an error
    except Exception as e:
        print(f"Critical launcher error: {e}")
        import traceback

        traceback.print_exc()
        success = False
    finally:
        # Ensure clean exit on Ubuntu
        import platform

        success = locals().get("success", True)  # Default to success if undefined
        if platform.system() == "Linux":
            # Force clean exit
            os._exit(0 if success else 1)
        elif platform.system() == "Darwin":  # Add macOS support
            # Force clean exit on macOS too
            os._exit(0 if success else 1)
        else:
            sys.exit(0 if success else 1)
