# Package Structure Guide - lunaNMR v0.9

Complete guide to the new modular package structure of lunaNMR v0.9, designed for professional Python development, maintainability, and extensibility.

## 📋 Table of Contents

1. [Overview](#overview)
2. [Directory Structure](#directory-structure)
3. [Module Descriptions](#module-descriptions)
4. [Import Patterns](#import-patterns)
5. [Development Workflow](#development-workflow)
6. [Extension Guidelines](#extension-guidelines)

---

## 🎯 Overview

lunaNMR v0.9 has been completely restructured from a flat file organization to a professional Python package architecture. This change provides:

**Benefits of New Structure:**
- ✅ **Modular Design**: Clear separation of concerns
- ✅ **Professional Organization**: Standard Python package layout  
- ✅ **Easy Extension**: Simple to add new components
- ✅ **Better Testing**: Isolated modules for unit testing
- ✅ **Import Clarity**: Explicit import paths
- ✅ **Optional Components**: DynamiXs as optional submodule
- ✅ **Documentation**: Self-documenting code organization

**Migration from v0.8:**
- Flat file structure → Hierarchical package structure
- Single GUI entry point → Unified launcher with app selector
- Mixed responsibilities → Separated concerns
- Implicit dependencies → Explicit import statements

---

## 📁 Directory Structure

### **Complete Package Layout**

```
lunaNMR_v0o9/                     # Root directory
│
├── 🚀 launch_lunaNMR.py          # Unified launcher with app selector
├── 📄 README.md                  # Project documentation
├── ⚙️ requirements.txt           # Dependencies
│
├── 📦 lunaNMR/                   # MAIN PACKAGE
│   ├── __init__.py              # Package initialization
│   │
│   ├── 🖥️ gui/                   # Graphical User Interface
│   │   ├── __init__.py
│   │   ├── main_gui.py          # Primary application interface
│   │   ├── gui_components.py    # Reusable GUI widgets
│   │   └── spectrum_browser.py  # Interactive spectrum viewer
│   │
│   ├── ⚙️ core/                  # Core Processing Engines
│   │   ├── __init__.py
│   │   ├── core_integrator.py           # Base integration engine
│   │   ├── enhanced_voigt_fitter.py     # Advanced Voigt fitting
│   │   ├── enhanced_peak_picker.py      # Peak detection algorithms
│   │   └── integrated_detection_fitter.py # Unified detection+fitting
│   │
│   ├── 📊 processors/            # Spectrum Processors
│   │   ├── __init__.py
│   │   ├── series_processor.py         # Time-series analysis (T₁, T₂)
│   │   ├── multi_spectrum_processor.py # Multi-condition analysis
│   │   ├── single_spectrum_processor.py # Individual spectrum analysis
│   │   └── parallel_fitting.py         # Parallel processing support
│   │
│   ├── 🔧 integrators/           # Specialized Integrators
│   │   ├── __init__.py
│   │   ├── inplace_advanced_nmr_integrator.py # Advanced in-memory analysis
│   │   ├── inplace_series_nmr_integrator.py   # Series with memory optimization
│   │   └── simple_pattern_matcher.py          # Pattern recognition
│   │
│   ├── 🛠️ utils/                 # Utility Classes
│   │   ├── __init__.py
│   │   ├── config_manager.py            # Configuration management
│   │   ├── file_manager.py              # File I/O operations
│   │   ├── parameter_manager.py         # Parameter handling
│   │   └── global_optimization_manager.py # System-wide optimization
│   │
│   └── ✅ validation/            # Validation & Testing
│       ├── __init__.py
│       └── verify_installation.py       # Installation verification
│
├── 🧩 modules/                   # Optional Modules
│   ├── __init__.py
│   └── dynamiXs/                # DynamiXs Integration (Optional)
│       ├── __init__.py
│       ├── dynamiXs_GUI.py      # DynamiXs graphical interface
│       ├── main_dynamiXs.py     # DynamiXs main application
│       └── [...other DynamiXs files...]
│
└── 📚 docs/                      # Documentation
    ├── INSTALLATION.md           # Installation guide
    ├── QUICKSTART.md            # Quick start tutorial
    ├── ARCHITECTURE.md          # Architecture overview
    ├── API_REFERENCE.md         # API documentation
    ├── PACKAGE_STRUCTURE.md    # This document
    ├── EXAMPLES.md              # Usage examples
    ├── TROUBLESHOOTING.md       # Common issues
    ├── FILE_FORMATS.md          # Supported file formats
    └── CHANGELOG.md             # Version history
```

---

## 📋 Module Descriptions

### **🚀 Launch System**

**File**: `launch_lunaNMR.py`
**Purpose**: Unified application launcher with professional GUI selector

**Features**:
- Application selection interface (LunaNMR vs DynamiXs)
- Dependency checking and validation  
- Graceful error handling with user feedback
- Cross-platform compatibility

**Usage**:
```bash
python3 launch_lunaNMR.py
```

### ** lunaNMR/ - Main Package**

#### **lunaNMR/__init__.py**
**Purpose**: Package initialization and public API exports

**Key Exports**:
```python
# Core functionality
from .core import CoreIntegrator, EnhancedVoigtFitter, EnhancedPeakPicker
from .utils import ConfigManager, FileManager

# Version information
__version__ = "0.9.0"
__all__ = [...]
```

#### **🖥️ lunaNMR/gui/ - User Interface**

**Components**:
- `main_gui.py`: Primary application with complete workflow
- `gui_components.py`: Reusable widgets and dialogs
- `spectrum_browser.py`: Interactive spectrum visualization

**Responsibilities**:
- User interaction and event handling
- Data visualization and plotting
- Configuration interfaces
- Results presentation and export

**Architecture**: Model-View-Controller pattern
- **Model**: Analysis results and configurations
- **View**: GUI components and visualizations
- **Controller**: Event handlers and business logic

#### **⚙️ lunaNMR/core/ - Processing Engines**

**Components**:
- `core_integrator.py`: Base integration class with fundamental operations
- `enhanced_voigt_fitter.py`: Advanced Voigt fitting with multi-peak support
- `enhanced_peak_picker.py`: Professional peak detection with clustering
- `integrated_detection_fitter.py`: Unified detection and fitting pipeline

**Key Features**:
- Mathematical algorithms for NMR analysis
- Baseline correction methods (ArPLS, polynomial, iterative)
- Multi-peak deconvolution
- Quality assessment and uncertainty estimation

#### **📊 lunaNMR/processors/ - Analysis Workflows**

**Components**:
- `series_processor.py`: T₁, T₂, NOE relaxation analysis
- `multi_spectrum_processor.py`: Temperature, pH, titration studies
- `single_spectrum_processor.py`: Individual spectrum analysis
- `parallel_fitting.py`: Multi-core processing support

**Design Pattern**: Strategy pattern for different analysis types

#### **🔧 lunaNMR/integrators/ - Specialized Analysis**

**Components**:
- `inplace_advanced_nmr_integrator.py`: Memory-efficient advanced analysis
- `inplace_series_nmr_integrator.py`: Series analysis with optimization
- `simple_pattern_matcher.py`: Automated pattern recognition

**Purpose**: Domain-specific integration methods for specialized experiments

#### **🛠️ lunaNMR/utils/ - Cross-Cutting Concerns**

**Components**:
- `config_manager.py`: Centralized configuration with validation
- `file_manager.py`: Multi-format file I/O (Bruker, Varian, etc.)
- `parameter_manager.py`: Parameter validation and optimization
- `global_optimization_manager.py`: System-wide performance optimization

**Design Pattern**: Singleton pattern for managers, Factory pattern for file handlers

#### **✅ lunaNMR/validation/ - Quality Assurance**

**Components**:
- `verify_installation.py`: Comprehensive installation verification

**Features**:
- Dependency checking
- System compatibility validation
- Performance benchmarking
- Installation troubleshooting

### **🧩 modules/ - Optional Extensions**

#### **modules/dynamiXs/ - Dynamic Exchange Analysis**

**Purpose**: Optional submodule for advanced relaxation analysis
**Integration**: Plugin architecture with independent development
**Availability**: Detected and loaded dynamically

**Benefits of Optional Design**:
- Reduces core package dependencies
- Enables independent development cycles
- Allows specialized licensing
- Maintains modularity

---

##  Import Patterns

### **Recommended Import Styles**

#### **Core Components**
```python
# Preferred: Explicit imports
from lunaNMR.core import CoreIntegrator, EnhancedVoigtFitter
from lunaNMR.utils import ConfigManager, FileManager
from lunaNMR.processors import SeriesProcessor

# Alternative: Module imports
from lunaNMR import core, utils, processors
integrator = core.CoreIntegrator()
config = utils.ConfigManager()
```

#### **GUI Components**
```python
# For GUI development
from lunaNMR.gui import main_gui, gui_components, spectrum_browser

# For launching GUI
from lunaNMR.gui.main_gui import main
main()  # Launch GUI directly
```

#### **Optional Modules**
```python
# Safe optional imports
try:
    from modules.dynamiXs import DynamiXsGUI, run_dynamixs
    DYNAMIXS_AVAILABLE = True
except ImportError:
    DYNAMIXS_AVAILABLE = False
```

### **Import Hierarchy**

```
Application Level
    ↑ imports from
GUI Level (lunaNMR.gui)
    ↑ imports from  
Processor Level (lunaNMR.processors)
    ↑ imports from
Core Level (lunaNMR.core)
    ↑ imports from
Utility Level (lunaNMR.utils)
```

**Rules**:
- ✅ Higher levels can import from lower levels
- ❌ Lower levels should not import from higher levels
- ✅ Same-level imports are allowed with caution
- ✅ Utils can be imported by any level

---

## 🛠️ Development Workflow

### **Adding New Components**

#### **1. New Core Algorithm**
```python
# lunaNMR/core/custom_fitter.py
from .core_integrator import BaseIntegrator

class CustomFitter(BaseIntegrator):
    def __init__(self, config=None):
        super().__init__(config)

    def fit_custom_model(self, data):
        # Implementation
        pass

# lunaNMR/core/__init__.py - Add export
from .custom_fitter import CustomFitter
__all__.append('CustomFitter')
```

#### **2. New Processor**
```python
# lunaNMR/processors/custom_processor.py
from ..core import CoreIntegrator
from ..utils import ConfigManager

class CustomProcessor:
    def __init__(self, config=None):
        self.config = config or ConfigManager()
        self.integrator = CoreIntegrator(config)

    def process_custom_experiment(self, data):
        # Implementation
        pass

# lunaNMR/processors/__init__.py - Add export
from .custom_processor import CustomProcessor
__all__.append('CustomProcessor')
```

#### **3. New GUI Component**
```python
# lunaNMR/gui/custom_dialog.py
import tkinter as tk
from tkinter import ttk

class CustomDialog:
    def __init__(self, parent):
        self.dialog = tk.Toplevel(parent)
        self.setup_ui()

    def setup_ui(self):
        # GUI implementation
        pass

# lunaNMR/gui/__init__.py - Add export
from .custom_dialog import CustomDialog
__all__.append('CustomDialog')
```


This new package structure provides a solid foundation for current and future development while maintaining backward compatibility where possible. The modular design ensures that lunaNMR can grow and adapt to new requirements while maintaining code quality and ease of use.
