# DynamiXs Folder Integration Guide

## Overview

The `dynamiXs` folder contains a complete submodule for dynamic NMR analysis with its own GUI and specialized functions. This guide explains how to properly integrate it into the new professional package structure.

## Current dynamiXs Structure

```
dynamiXs/
├── run_dynamixs_gui.py              # Main launcher for dynamiXs GUI
├── dynamiXs_gui.py                  # Main dynamiXs GUI
├── dynamiXs_cpmg/                   # CPMG relaxation dispersion
│   ├── cpmg_RD_bootstrap.py
│   └── cpmg_RD.py
├── dynamiXs_density_functions/      # Spectral density functions
│   ├── ZZ_multi_2fields_density.py
│   ├── ZZ_density.py
│   ├── ZZ_2fields_density_claude_rex_mcmc_error.py
│   ├── ZZ_2fields_density_087.py
│   ├── ZZ_multi_density_087.py
│   ├── ZZ_2fields_density.py
│   ├── ZZ_multi_density.py
│   └── ZZ_multi_2fields_density_087.py
├── dynamiXs_format/                 # Data formatting utilities
│   ├── NH_rename_list_iter.py
│   └── NH_compile_list_horizontal.py
├── dynamiXs_plot/                   # Plotting and visualization
│   ├── compare_datasets_087.py
│   ├── compare_datasets.py
│   ├── plot_dual_field_datasets.py
│   └── plot_dual_field_datasets_087.py
└── dynamiXs_T1_T2/                  # T1/T2 relaxation analysis
    ├── fit_Tx_NMRRE.py
    └── fitmulti__Tx_NMRRE.py
```

## **Option 1: Keep DynamiXs as Independent Submodule (Recommended)**

### **Step 1: Move DynamiXs to Dedicated Module Location**

```bash
# Create dedicated modules directory for independent tools
mkdir -p modules

# Move dynamiXs as a complete submodule
mv dynamiXs modules/dynamiXs

# Create init file for modules directory
touch modules/__init__.py
```

### **Step 2: Create DynamiXs Package Init**

**modules/dynamiXs/__init__.py:**
```python
"""
DynamiXs: Dynamic NMR Analysis Toolkit

A comprehensive toolkit for NMR dynamics analysis including:
- CPMG relaxation dispersion analysis
- Spectral density function modeling
- T1/T2 relaxation analysis
- Multi-field data comparison
- Advanced plotting and visualization
"""

__version__ = "1.0.0"
__author__ = "Guillaume"
__description__ = "Dynamic NMR Analysis Toolkit"

# Main GUI access
try:
    from .dynamiXs_gui import DynamiXsGUI
except ImportError:
    pass

# Main launcher
try:
    from .run_dynamixs_gui import main as run_dynamixs
except ImportError:
    pass

__all__ = [
    'DynamiXsGUI',
    'run_dynamixs'
]
```

### **Step 3: Update Main lunaNMR Package to Include DynamiXs**

**lunaNMR/__init__.py** (add to existing content):
```python
# DynamiXs Integration (optional submodule)
try:
    import sys
    import os
    # Add modules directory to path
    modules_path = os.path.join(os.path.dirname(__file__), '..', 'modules')
    if modules_path not in sys.path:
        sys.path.append(modules_path)
    
    from dynamiXs import DynamiXsGUI, run_dynamixs
    __all__.extend(['DynamiXsGUI', 'run_dynamixs'])
except ImportError:
    pass
```

### **Step 4: Create Unified Launcher**

**launch_lunaNMR.py** (add to existing main() function):
```python
def show_application_selector():
    """Show dialog to select which application to launch"""
    root = tk.Tk()
    root.withdraw()
    
    choice = messagebox.askyesnocancel(
        "LunaNMR Application Selector",
        "Select Application to Launch:\n\n"
        "YES: LunaNMR Peak Analysis (Main)\n"
        "NO: DynamiXs Dynamics Analysis\n"
        "CANCEL: Exit"
    )
    
    root.destroy()
    
    if choice is True:
        return "lunaNMR"
    elif choice is False:
        return "dynamiXs"
    else:
        return None

# Add to main() function before launching GUI:
def main():
    """Main entry point for LunaNMR applications"""
    try:
        # Setup paths
        setup_paths()
        
        # Check dependencies
        missing = check_dependencies()
        if missing:
            # ... existing error handling ...
        
        # Show application selector
        app_choice = show_application_selector()
        
        if app_choice == "lunaNMR":
            # Launch main LunaNMR GUI (existing code)
            from lunaNMR.gui.main_gui import NMRPeaksSeriesGUI
            root = tk.Tk()
            app = NMRPeaksSeriesGUI(root)
            root.title("LunaNMR v0.9 - Advanced NMR Peak Analysis")
            root.mainloop()
            
        elif app_choice == "dynamiXs":
            # Launch DynamiXs GUI
            try:
                from modules.dynamiXs import run_dynamixs
                run_dynamixs()
            except ImportError:
                messagebox.showerror(
                    "DynamiXs Not Available", 
                    "DynamiXs module not found or not properly installed."
                )
        
        # If None (cancelled), just exit
            
    except Exception as e:
        # ... existing error handling ...
```

## **Option 2: Full Integration into LunaNMR Package**

### **Step 1: Move DynamiXs Components into Main Package**

```bash
# Create dynamics module within main package
mkdir -p lunaNMR/dynamics

# Move and organize dynamiXs components
mv dynamiXs/dynamiXs_gui.py lunaNMR/dynamics/dynamics_gui.py
mv dynamiXs/run_dynamixs_gui.py lunaNMR/dynamics/launcher.py

# Create submodules for different analysis types
mkdir -p lunaNMR/dynamics/{cpmg,density_functions,formatting,plotting,relaxation}

# Move CPMG files
mv dynamiXs/dynamiXs_cpmg/* lunaNMR/dynamics/cpmg/

# Move density functions
mv dynamiXs/dynamiXs_density_functions/* lunaNMR/dynamics/density_functions/

# Move formatting utilities
mv dynamiXs/dynamiXs_format/* lunaNMR/dynamics/formatting/

# Move plotting functions
mv dynamiXs/dynamiXs_plot/* lunaNMR/dynamics/plotting/

# Move T1/T2 analysis
mv dynamiXs/dynamiXs_T1_T2/* lunaNMR/dynamics/relaxation/

# Create __init__.py files for all subdirectories
touch lunaNMR/dynamics/__init__.py
touch lunaNMR/dynamics/{cpmg,density_functions,formatting,plotting,relaxation}/__init__.py

# Remove empty dynamiXs folder
rmdir dynamiXs/dynamiXs_*
rmdir dynamiXs
```

### **Step 2: Create Dynamics Package Structure**

**lunaNMR/dynamics/__init__.py:**
```python
"""
Dynamics Analysis Module for LunaNMR

Advanced NMR dynamics analysis including relaxation dispersion,
spectral density modeling, and multi-field analysis.
"""

try:
    from .dynamics_gui import DynamiXsGUI
except ImportError:
    pass

try:
    from .launcher import main as launch_dynamics
except ImportError:
    pass

# Submodule imports
try:
    from . import cpmg
    from . import density_functions
    from . import formatting
    from . import plotting
    from . import relaxation
except ImportError:
    pass

__all__ = [
    'DynamiXsGUI',
    'launch_dynamics',
    'cpmg',
    'density_functions', 
    'formatting',
    'plotting',
    'relaxation'
]
```

## **Updated Final Directory Structure**

### **Option 1 Structure (Recommended):**
```
lunaNMR_v0o9/
├── README.md
├── requirements.txt  
├── launch_lunaNMR.py                  # Unified launcher with app selector
├── 
├── lunaNMR/                          # Main NMR analysis package
│   ├── __init__.py                   # Includes dynamiXs integration
│   ├── gui/
│   ├── core/
│   ├── processors/
│   ├── integrators/
│   ├── utils/
│   └── validation/
│
├── modules/                          # Independent modules ✅ NEW
│   ├── __init__.py                   # ✅ NEW
│   └── dynamiXs/                     # ✅ MOVED from root
│       ├── __init__.py               # ✅ NEW
│       ├── dynamiXs_gui.py           # ✅ MOVED
│       ├── run_dynamixs_gui.py       # ✅ MOVED
│       ├── dynamiXs_cpmg/            # ✅ MOVED
│       ├── dynamiXs_density_functions/ # ✅ MOVED
│       ├── dynamiXs_format/          # ✅ MOVED
│       ├── dynamiXs_plot/            # ✅ MOVED
│       └── dynamiXs_T1_T2/           # ✅ MOVED
│
├── docs/
└── md/
```

### **Option 2 Structure:**
```
lunaNMR_v0o9/
├── README.md
├── requirements.txt
├── launch_lunaNMR.py
├── 
├── lunaNMR/                          # Unified package
│   ├── __init__.py
│   ├── gui/
│   ├── core/
│   ├── processors/
│   ├── integrators/
│   ├── utils/
│   ├── validation/
│   └── dynamics/                     # ✅ NEW - integrated dynamics
│       ├── __init__.py               # ✅ NEW
│       ├── dynamics_gui.py           # ✅ RENAMED from dynamiXs_gui.py
│       ├── launcher.py               # ✅ RENAMED from run_dynamixs_gui.py
│       ├── cpmg/                     # ✅ REORGANIZED
│       ├── density_functions/        # ✅ REORGANIZED
│       ├── formatting/               # ✅ REORGANIZED
│       ├── plotting/                 # ✅ REORGANIZED
│       └── relaxation/               # ✅ REORGANIZED
│
├── docs/
└── md/
```

## **Recommendation**

**Use Option 1 (Independent Submodule)** because:

✅ **Preserves modularity** - DynamiXs remains a distinct tool  
✅ **Easier maintenance** - Separate development paths  
✅ **Cleaner architecture** - Clear separation of concerns  
✅ **User choice** - Users can install just what they need  
✅ **Minimal disruption** - Existing DynamiXs code needs fewer changes  

## **Implementation Commands for Option 1**

```bash
# Execute these commands after completing the main reorganization:

# Create modules directory
mkdir -p modules
touch modules/__init__.py

# Move dynamiXs as complete submodule
mv dynamiXs modules/dynamiXs

# Create dynamiXs package init
cat > modules/dynamiXs/__init__.py << 'EOF'
"""
DynamiXs: Dynamic NMR Analysis Toolkit
"""
__version__ = "1.0.0"

try:
    from .dynamiXs_gui import DynamiXsGUI
    from .run_dynamixs_gui import main as run_dynamixs
except ImportError:
    pass

__all__ = ['DynamiXsGUI', 'run_dynamixs']
EOF

echo "✅ DynamiXs integration complete!"
echo "📁 DynamiXs moved to: modules/dynamiXs/"
echo "🚀 Launch with unified launcher: python3 launch_lunaNMR.py"
```

This approach maintains the professional package structure while preserving DynamiXs as a powerful, independent analysis toolkit.