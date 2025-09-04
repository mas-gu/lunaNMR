# Architecture Guide - lunaNMR v0.9

Complete architectural overview of lunaNMR v0.9.

## Table of Contents

1. [Package Structure](#package-structure)
2. [Core Architecture](#core-architecture)
3. [Module Descriptions](#module-descriptions)
4. [Data Flow](#data-flow)
5. [Extension Points](#extension-points)
6. [Integration Patterns](#integration-patterns)
7. [Development Guidelines](#development-guidelines)

---

## Package Structure

### ** Modular Architecture**

lunaNMR v0.9 implements a professional Python package structure with clear separation of concerns:

```
lunaNMR_v0o9/
├── lunaNMR/                    # Main package
│   ├── __init__.py            # Package initialization & exports
│   ├── gui/                   # Graphical user interface
│   │   ├── __init__.py
│   │   ├── main_gui.py        # Primary GUI application
│   │   ├── gui_components.py  # Reusable GUI components
│   │   └── spectrum_browser.py # Interactive spectrum viewer
│   ├── core/                  # Core processing engines
│   │   ├── __init__.py
│   │   ├── core_integrator.py           # Base integration engine
│   │   ├── enhanced_voigt_fitter.py     # Advanced Voigt fitting
│   │   ├── enhanced_peak_picker.py      # Peak detection algorithms
│   │   └── integrated_detection_fitter.py # Unified detection+fitting
│   ├── processors/            # Spectrum processors
│   │   ├── __init__.py
│   │   ├── series_processor.py         # Time-series analysis
│   │   ├── multi_spectrum_processor.py # Multi-condition analysis
│   │   ├── single_spectrum_processor.py # Single spectrum analysis
│   │   └── parallel_fitting.py         # Parallel processing
│   ├── integrators/           # Specialized integrators
│   │   ├── __init__.py
│   │   ├── inplace_advanced_nmr_integrator.py
│   │   ├── inplace_series_nmr_integrator.py
│   │   └── simple_pattern_matcher.py
│   ├── utils/                 # Utility classes
│   │   ├── __init__.py
│   │   ├── config_manager.py          # Configuration management
│   │   ├── file_manager.py            # File I/O operations
│   │   ├── parameter_manager.py       # Parameter handling
│   │   └── global_optimization_manager.py
│   └── validation/            # Validation & testing tools
│       ├── __init__.py
│       └── verify_installation.py
├── modules/                   # Optional modules
│   ├── __init__.py
│   └── dynamiXs/             # DynamiXs integration (optional)
│       ├── __init__.py
│       └── [DynamiXs files...]
├── docs/                     # Documentation
├── launch_lunaNMR.py         # Unified launcher with app selector
└── README.md
```

---

## Core Architecture

### **Layered Architecture Pattern**

lunaNMR follows a layered architecture with clear responsibilities:

```
┌─────────────────────────────────────────┐
│             Presentation Layer          │
│          lunaNMR.gui.main_gui           │
│     User Interface & Visualization      │
├─────────────────────────────────────────┤
│            Application Layer            │
│        lunaNMR.processors.*             │
│      Business Logic & Workflows         │
├─────────────────────────────────────────┤
│              Service Layer              │
│         lunaNMR.integrators.*           │
│      Specialized Analysis Services      │
├─────────────────────────────────────────┤
│               Core Layer                │
│           lunaNMR.core.*                │
│     Fundamental Processing Engines      │
├─────────────────────────────────────────┤
│            Infrastructure Layer         │
│           lunaNMR.utils.*               │
│    Configuration, I/O, Validation       │
└─────────────────────────────────────────┘
```

### **Component Interaction Model**

```python
# Example of how components interact
from lunaNMR.core import EnhancedVoigtFitter
from lunaNMR.processors import SingleSpectrumProcessor
from lunaNMR.utils import ConfigManager, FileManager
from lunaNMR.gui import main_gui

# Configuration flows down through all layers
config = ConfigManager()

# Core components provide fundamental algorithms
fitter = EnhancedVoigtFitter(config)

# Processors orchestrate workflows
processor = SingleSpectrumProcessor(fitter, config)

# GUI coordinates user interaction
gui = main_gui.NMRPeaksSeriesGUI(processor)
```

---

## 📋 Module Descriptions

### **lunaNMR.core - Core Processing Engines**

**Purpose**: Fundamental NMR analysis algorithms and mathematical operations.

**Key Components**:
- `CoreIntegrator`: Base class for all integration operations
- `EnhancedVoigtFitter`: Advanced Voigt profile fitting with multi-peak support
- `EnhancedPeakPicker`: Professional peak detection with clustering algorithms
- `IntegratedDetectionFitter`: Unified peak detection and fitting pipeline

**Responsibilities**:
- Mathematical peak fitting algorithms
- Signal processing operations  
- Baseline correction methods
- Quality assessment metrics

**Import Pattern**:
```python
from lunaNMR.core import EnhancedVoigtFitter, EnhancedPeakPicker
from lunaNMR.core.core_integrator import VoigtIntegrator
```

### **lunaNMR.processors - Spectrum Processors**

**Purpose**: High-level analysis workflows for different experimental scenarios.

**Key Components**:
- `SeriesProcessor`: Time-series relaxation analysis (T₁, T₂, NOE)
- `MultiSpectrumProcessor`: Multi-condition experiments (titrations, temperature)
- `SingleSpectrumProcessor`: Individual spectrum analysis
- `ParallelFitting`: Parallel processing for large datasets

**Responsibilities**:
- Orchestrating multi-step analysis workflows
- Managing data transformations
- Coordinating between core components
- Batch processing capabilities

**Import Pattern**:
```python
from lunaNMR.processors import SeriesProcessor, MultiSpectrumProcessor
```

### **lunaNMR.integrators - Specialized Integrators**

**Purpose**: Domain-specific integration methods for specialized NMR experiments.

**Key Components**:
- `InplaceAdvancedNMRIntegrator`: In-memory advanced analysis
- `InplaceSeriesNMRIntegrator`: Series analysis with memory optimization
- `SimplePatternMatcher`: Pattern recognition for automated assignments

**Responsibilities**:
- Specialized analysis algorithms
- Memory-efficient processing
- Pattern recognition and classification

### **lunaNMR.gui - Graphical User Interface**

**Purpose**: User interface components and visualization tools.

**Key Components**:
- `main_gui.py`: Primary application interface
- `gui_components.py`: Reusable UI widgets
- `spectrum_browser.py`: Interactive spectrum visualization

**Responsibilities**:
- User interaction management
- Data visualization
- Configuration interfaces
- Results presentation

**Architecture Pattern**: Model-View-Controller (MVC)
- **Model**: Data and analysis results
- **View**: GUI components and plots  
- **Controller**: Event handlers and workflow coordination

### **lunaNMR.utils - Utility Classes**

**Purpose**: Cross-cutting concerns and infrastructure services.

**Key Components**:
- `ConfigManager`: Configuration management
- `FileManager`: File I/O operations with format support
- `ParameterManager`: Parameter validation and optimization
- `GlobalOptimizationManager`: System-wide optimization

**Responsibilities**:
- Configuration management
- File format handling
- Parameter validation
- System utilities

### **lunaNMR.validation - Validation & Testing**

**Purpose**: Installation verification and system testing.

**Key Components**:
- `verify_installation.py`: Installation verification script

**Responsibilities**:
- Dependency checking
- System compatibility validation
- Installation troubleshooting

### **modules/ - Optional Modules**

**Purpose**: Optional functionality that extends lunaNMR capabilities.

**Structure**:
- `modules.dynamiXs`: Dynamic exchange analysis (T₁, T₂, NOE fitting)

**Design Pattern**: Plugin architecture allowing independent development and deployment.

---

## Data Flow

### **Typical Analysis Workflow**

```python
# 1. Configuration & Setup
config = ConfigManager()
file_manager = FileManager()

# 2. Data Loading
spectrum_data = file_manager.load_spectrum('data.ft2')
peak_list = file_manager.load_peak_list('peaks.csv')

# 3. Core Processing
peak_picker = EnhancedPeakPicker(config)
detected_peaks = peak_picker.detect_peaks(spectrum_data)

voigt_fitter = EnhancedVoigtFitter(config)
fitting_results = voigt_fitter.fit_peaks(detected_peaks, spectrum_data)

# 4. High-Level Processing
processor = SingleSpectrumProcessor(voigt_fitter, config)
analysis_results = processor.process_spectrum(spectrum_data, peak_list)

# 5. Results Export
file_manager.export_results(analysis_results, 'results.csv')
```

### **Data Flow Diagram**

```
User Input → ConfigManager → FileManager → Core Engines → Processors → GUI Display
    ↓              ↓             ↓            ↓             ↓           ↓
Configuration → File I/O → Peak Detection → Analysis → Visualization → Export
```

### **Memory Management**

- **Lazy Loading**: Large datasets loaded only when needed
- **Streaming Processing**: Large series processed in chunks
- **Memory Monitoring**: Automatic memory usage optimization
- **Garbage Collection**: Explicit cleanup of temporary objects

---

## 🔧 Extension Points

### **Adding Custom Fitting Algorithms**

```python
# 1. Create custom fitter inheriting from base classes
from lunaNMR.core.core_integrator import BaseIntegrator

class CustomFitter(BaseIntegrator):
    def fit_peak(self, x_data, y_data, initial_params):
        # Implement custom fitting logic
        pass

# 2. Register with the system
from lunaNMR.utils import ParameterManager
ParameterManager.register_fitter('custom', CustomFitter)
```

### **Creating New Processors**

```python
# 1. Inherit from base processor
from lunaNMR.processors.base_processor import BaseProcessor

class TemperatureSeriesProcessor(BaseProcessor):
    def process_temperature_series(self, temperature_data):
        # Implement temperature-dependent analysis
        pass

# 2. Integration with GUI
# Add menu item in main_gui.py to access new processor
```

### **Plugin Development**

```python
# modules/custom_plugin/__init__.py
class CustomPlugin:
    def __init__(self):
        self.name = "Custom Analysis Plugin"

    def get_analysis_methods(self):
        return {'custom_analysis': self.custom_analysis}

    def custom_analysis(self, data):
        # Custom analysis implementation
        pass
```

---

## Integration Patterns

### **Dependency Injection**

Components receive dependencies through constructors:

```python
# Dependencies injected at initialization
class EnhancedVoigtFitter:
    def __init__(self, config_manager, parameter_manager):
        self.config = config_manager
        self.params = parameter_manager
```

### **Observer Pattern**

Progress updates and notifications:

```python
class AnalysisProgress:
    def __init__(self):
        self.observers = []

    def add_observer(self, observer):
        self.observers.append(observer)

    def notify_progress(self, current, total):
        for observer in self.observers:
            observer.update_progress(current, total)
```

### **Factory Pattern**

Dynamic component creation:

```python
class ProcessorFactory:
    @staticmethod
    def create_processor(processor_type, config):
        if processor_type == 'series':
            return SeriesProcessor(config)
        elif processor_type == 'multi_spectrum':
            return MultiSpectrumProcessor(config)
        # etc.
```

---

## Development Guidelines

### **Import Structure**

**Absolute Imports**: Always use absolute imports from the package root:
```python
# Correct
from lunaNMR.core import EnhancedVoigtFitter
from lunaNMR.utils.config_manager import ConfigManager

# Avoid relative imports in most cases
```

### **Configuration Management**

**Centralized Configuration**: All components should accept configuration:
```python
class NewComponent:
    def __init__(self, config=None):
        self.config = config or ConfigManager.get_default()
```

### **Error Handling**

**Graceful Degradation**: Components should handle errors gracefully:
```python
try:
    result = enhanced_analysis(data)
except AnalysisError:
    # Fall back to basic analysis
    result = basic_analysis(data)
```

### **Testing Strategy**

**Unit Tests**: Each module should have comprehensive unit tests
**Integration Tests**: Test component interactions
**Validation Tests**: Verify against known standards

### **Documentation Standards**

**Docstring Format**: Use NumPy-style docstrings
**Type Hints**: Include type annotations for all public APIs
**Examples**: Provide usage examples in docstrings

---

## 🚀 Performance Considerations

### **Optimization Strategies**

1. **Vectorized Operations**: Use NumPy for mathematical operations
2. **Parallel Processing**: Leverage multiprocessing for independent tasks
3. **Memory Management**: Implement memory-efficient algorithms
4. **Caching**: Cache expensive computations
5. **Lazy Loading**: Load data only when needed

### **Scalability**

The modular architecture enables:
- **Horizontal Scaling**: Multiple instances for different datasets
- **Vertical Scaling**: Optimized algorithms for larger datasets
- **Distributed Processing**: Framework ready for cluster deployment

---

This architecture provides a solid foundation for current NMR analysis needs while maintaining flexibility for future enhancements and extensions. The modular design ensures maintainability, testability, and extensibility as the software continues to evolve.
