#!/usr/bin/env python3
"""
LunaNMR Suite - Unified Launch Script

Unified launcher for the LunaNMR suite including:
- LunaNMR: Advanced NMR Peak Analysis and Integration
- DynamiXs: Dynamic exchange analysis (optional module)

This script handles environment setup, application selection, and graceful error handling.
"""

import sys
import os
import tkinter as tk
from tkinter import messagebox, ttk

def check_dependencies():
    """Check if all required dependencies are available"""
    required_modules = [
        'pandas', 'numpy', 'matplotlib', 'scipy', 'sklearn', 'nmrglue'
    ]

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

    # Add directories to Python path
    sys.path.insert(0, current_dir)
    sys.path.insert(0, parent_dir)

    return current_dir, parent_dir

def check_dynamixs_availability():
    """Check if DynamiXs module is available"""
    try:
        from modules.dynamiXs import DynamiXsGUI, run_dynamixs
        return True
    except ImportError:
        return False

def show_application_selector():
    """Show application selector GUI"""
    root = tk.Tk()
    root.title("LunaNMR Suite v0.9")
    root.geometry("400x300")
    root.resizable(False, False)

    # Center the window
    root.update_idletasks()
    x = (root.winfo_screenwidth() // 2) - (400 // 2)
    y = (root.winfo_screenheight() // 2) - (300 // 2)
    root.geometry(f"400x300+{x}+{y}")

    selected_app = tk.StringVar(value="lunaNMR")

    # Main frame
    main_frame = ttk.Frame(root, padding="20")
    main_frame.pack(fill=tk.BOTH, expand=True)

    # Title
    title_label = ttk.Label(main_frame, text="LunaNMR Suite v0.9",
                           font=("Arial", 16, "bold"))
    title_label.pack(pady=(0, 20))

    # Subtitle
    subtitle_label = ttk.Label(main_frame, text="Select Application to Launch:")
    subtitle_label.pack(pady=(0, 15))

    # Application selection frame
    app_frame = ttk.LabelFrame(main_frame, text="Available Applications", padding="10")
    app_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 20))

    # LunaNMR option
    lunaNMR_frame = ttk.Frame(app_frame)
    lunaNMR_frame.pack(fill=tk.X, pady=(0, 10))

    ttk.Radiobutton(lunaNMR_frame, text="LunaNMR", variable=selected_app,
                   value="lunaNMR").pack(anchor=tk.W)
    ttk.Label(lunaNMR_frame, text="Advanced NMR Peak Analysis and Integration",
             foreground="gray").pack(anchor=tk.W, padx=(20, 0))

    # DynamiXs option
    dynamixs_frame = ttk.Frame(app_frame)
    dynamixs_frame.pack(fill=tk.X)

    dynamixs_available = check_dynamixs_availability()
    dynamixs_radio = ttk.Radiobutton(dynamixs_frame, text="DynamiXs",
                                    variable=selected_app, value="dynamixs")

    if dynamixs_available:
        dynamixs_radio.pack(anchor=tk.W)
        ttk.Label(dynamixs_frame, text="Dynamic Analysis", 
                 foreground="gray").pack(anchor=tk.W, padx=(20, 0))
    else:
        dynamixs_radio.configure(state=tk.DISABLED)
        dynamixs_radio.pack(anchor=tk.W)
        ttk.Label(dynamixs_frame, text="Dynamic Exchange Analysis (Not Available)",
                 foreground="red").pack(anchor=tk.W, padx=(20, 0))

    # Button frame
    button_frame = ttk.Frame(main_frame)
    button_frame.pack(fill=tk.X)

    result = {"app": None}

    def launch_selected():
        if selected_app.get() == "dynamixs" and not dynamixs_available:
            messagebox.showerror("Module Not Available",
                               "DynamiXs module is not available.\nPlease check the modules/dynamiXs directory.")
            return
        result["app"] = selected_app.get()
        root.quit()

    def cancel_launch():
        result["app"] = None
        root.quit()

    ttk.Button(button_frame, text="Launch", command=launch_selected).pack(side=tk.RIGHT, padx=(10, 0))
    ttk.Button(button_frame, text="Cancel", command=cancel_launch).pack(side=tk.RIGHT)

    root.mainloop()
    root.destroy()

    return result["app"]

def launch_lunaNMR():
    """Launch LunaNMR application"""
    try:
        from lunaNMR.gui.main_gui import main as gui_main
        gui_main()
        return True
    except ImportError as e:
        error_msg = f"Failed to import LunaNMR modules:\n{str(e)}\n\nPlease ensure the lunaNMR package is properly installed."
        print(f"Error: {error_msg}")
        messagebox.showerror("Import Error", error_msg)
        return False

def launch_dynamixs():
    """Launch DynamiXs application"""
    try:
        from modules.dynamiXs import run_dynamixs
        run_dynamixs()
        return True
    except ImportError as e:
        error_msg = f"Failed to import DynamiXs modules:\n{str(e)}\n\nPlease ensure the DynamiXs module is properly installed."
        print(f"Error: {error_msg}")
        messagebox.showerror("Import Error", error_msg)
        return False

def main():
    """Main launcher function"""
    print("LunaNMR Suite v0.9 - Starting...")

    # Setup paths
    current_dir, parent_dir = setup_paths()

    # Check dependencies
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

    # Show application selector
    selected_app = show_application_selector()

    if selected_app is None:
        print("Launch cancelled by user")
        return True

    # Launch selected application
    try:
        if selected_app == "lunaNMR":
            print("Launching LunaNMR...")
            success = launch_lunaNMR()
        elif selected_app == "dynamixs":
            print("Launching DynamiXs...")
            success = launch_dynamixs()
        else:
            print(f"Unknown application: {selected_app}")
            return False

        if success:
            print("Application closed normally")

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
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        #print("\n⏹️ Launch interrupted by user")
        sys.exit(0)
    except Exception as e:
        #print(f"💥 Critical launcher error: {e}")
        sys.exit(1)
