# 🌙 LunaNMR v0.9

**Advanced NMR Peak Analysis and Integration Suite**

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey.svg)](https://github.com/your-username/lunaNMR)

LunaNMR is a professional-grade software suite for Nuclear Magnetic Resonance (NMR) peak detection, Voigt profile fitting, and comprehensive analysis.

---

## ✨ Key Features

###  **Advanced Peak Analysis**
- **Advanced Voigt Fitting**: Multi-peak deconvolution with Gaussian and Lorentzian components
- **ArPLS Baseline Correction**: Asymmetrically Reweighted Penalized Least Squares algorithm
- **Multi-Peak Detection**: Automatic clustering and overlapping peak resolution
- **Quality Assessment**: R² values, uncertainties, and confidence intervals

###  **Comprehensive Analysis Workflows**
- **Single Spectrum Analysis**: Individual spectrum processing with full uncertainty analysis
- **Series Analysis**: T₁, T₂, and hetNOE relaxation studies with exponential fitting
- **Multi-Spectrum Processing**: Temperature, pH, and titration series analysis
- **Batch Processing**: Automated analysis of large datasets with parallel processing

###  **User-Friendly Interface**
- **Unified Launcher**: Application selector for LunaNMR and DynamiXs modules
- **Interactive GUI**: Point-and-click interface with real-time visualization
- **Spectrum Browser**: Interactive peak picking and annotation tools
- **Publication-Ready Plots**: High-quality figures with customizable styling

### **Advanced Capabilities**
- **Multiple File Formats**: Bruker, Varian, NMRPipe, SPARKY, and more
- **Parallel Processing**: Multi-core support for large-scale analysis
- **Extensible Architecture**: Plugin system for custom analysis methods
- **Professional Configuration**: Centralized parameter management and validation

---

## 🚀 Quick Start

### **Installation**

1. **Clone the repository:**
```bash
git clone https://github.com/your-username/lunaNMR.git
cd lunaNMR/lunaNMR_v0o9
```

2. **Install dependencies:**
```bash
pip install -r requirements.txt
```

3. **Verify installation:**
```bash
python3 lunaNMR/validation/verify_installation.py
```

4. **Launch the application:**
```bash
python3 launch_lunaNMR.py
```

### **5-Minute Quick Start**

```python
from lunaNMR.core import CoreIntegrator, EnhancedVoigtFitter
from lunaNMR.utils import ConfigManager

# Initialize components
config = ConfigManager()
integrator = CoreIntegrator(config=config)
fitter = EnhancedVoigtFitter(config=config)

# Load your data
integrator.load_data('peak_list.csv', 'spectrum.ft2')

# Detect and fit peaks
results = fitter.fit_peaks_enhanced(
    integrator.ppm_1h,
    integrator.intensities,
    fitting_method='voigt'
)

# View results
for i, result in enumerate(results):
    if result['success']:
        print(f"Peak {i+1}: {result['center']:.3f} ppm (R² = {result['r_squared']:.3f})")
```

---

## 📦 Package Structure

```
lunaNMR_v0o9/
├──  launch_lunaNMR.py          # Unified application launcher
├──  lunaNMR/                   # Main package
│   ├── ️ gui/                   # Graphical user interface
│   ├──  core/                  # Core processing engines
│   ├──  processors/            # Analysis workflows  
│   ├──  integrators/           # Specialized integrators
│   ├── ️ utils/                 # Configuration and utilities
│   └──  validation/            # Installation verification
├──  modules/                   # Optional modules
│   └── dynamiXs/                # Dynamic exchange analysis
└──  docs/                      # Comprehensive documentation
```

### **Core Components**

| Component | Purpose | Key Classes |
|-----------|---------|-------------|
| **lunaNMR.core** | Fundamental algorithms | `CoreIntegrator`, `EnhancedVoigtFitter`, `EnhancedPeakPicker` |
| **lunaNMR.processors** | Analysis workflows | `SeriesProcessor`, `MultiSpectrumProcessor` |
| **lunaNMR.gui** | User interface | `main_gui`, `spectrum_browser` |
| **lunaNMR.utils** | Configuration & I/O | `ConfigManager`, `FileManager` |
| **modules.dynamiXs** | Relaxation analysis | `DynamiXsGUI` (optional) |

---

##  Scientific Background

### **Voigt Profile Fitting**

LunaNMR implements professional-grade Voigt profile fitting, combining Gaussian and Lorentzian components to accurately model NMR lineshapes:

**Mathematical Foundation:**
- **Voigt Function**: Convolution of Gaussian (instrumental broadening) and Lorentzian (natural linewidth) functions
- **ArPLS Baseline Correction**: Asymmetrically Reweighted Penalized Least Squares with automatic λ selection
- **Multi-Peak Deconvolution**: Simultaneous fitting of overlapping peaks with constraints

**Key Algorithms:**
- **Peak Detection**: Network-based clustering with prominence analysis
- **Fitting Optimization**: Levenberg-Marquardt with robust parameter estimation
- **Uncertainty Analysis**: Bootstrap and covariance-based error estimation

---

## 🖥️ Graphical Interface

### **Application Selector**

Launch `python3 launch_lunaNMR.py` to access the unified launcher:

```
┌─────────────────────────────────────────┐
│           LunaNMR Suite v0.9            │
│                                         │
│  Select Application to Launch:          │
│                                         │
│  ○ LunaNMR                             │
│    Advanced NMR Peak Analysis          │
│                                         │
│  ○ DynamiXs                            │
│    Dynamic Exchange Analysis           │
│                                         │
│         [Launch]    [Cancel]            │
└─────────────────────────────────────────┘
```

---

---

## 📚 Documentation

### **Complete Documentation**

- **📋 [Installation Guide](docs/INSTALLATION.md)** - Comprehensive setup instructions
- **🚀 [Quick Start](docs/QUICKSTART.md)** - 5-minute tutorial
- **🏗️ [Architecture](docs/ARCHITECTURE.md)** - Technical architecture overview
- **📦 [Package Structure](docs/PACKAGE_STRUCTURE.md)** - Detailed package organization
- **📖 [API Reference](docs/API_REFERENCE.md)** - Complete API documentation
- **📂 [File Formats](docs/FILE_FORMATS.md)** - Supported data formats

## 🧪 Supported File Formats

### **Input Formats**

| Format | Extension | Vendor | Support |
|--------|-----------|---------|---------|
| **Bruker TopSpin** | `2ii`, `.2rr` | Bruker | ✅ Full |
| **Varian/Agilent** | `.fid`, `.ft` | Agilent | ✅ Full |
| **NMRPipe** | `.ft`, `.pipe` | NIH | ✅ Full |
| **SPARKY** | `.ucsf` | UCSF | ✅ Full |

### **Export Formats**

- **Data**: CSV, Excel, JSON, HDF5
- **Figures**: PNG, PDF, SVG, EPS
- **Reports**: HTML, LaTeX, Markdown

---

## 🤝 Contributing

We welcome contributions from the NMR community! Please see our contributing guidelines:

### **Development Setup**

```bash
# Clone repository
git clone https://github.com/your-username/lunaNMR.git
cd lunaNMR/lunaNMR_v0o9

# Create development environment
python3 -m venv dev_env
source dev_env/bin/activate  # or dev_env\Scripts\activate on Windows

# Install development dependencies
pip install -r requirements.txt
pip install -r requirements-dev.txt

# Run tests
python -m pytest tests/
```

### **Contributing Guidelines**

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
3. **Add tests** for new functionality
4. **Update documentation** as needed
5. **Commit** changes (`git commit -m 'Add amazing feature'`)
6. **Push** to branch (`git push origin feature/amazing-feature`)
7. **Create** a Pull Request

### **Code Standards**

- **Python**: PEP 8 style guide
- **Documentation**: NumPy-style docstrings
- **Testing**: pytest with >90% coverage
- **Type Hints**: Required for all public APIs

---

## 📈 Performance

### **Benchmarks**

| Operation | Dataset Size | Time | Memory |
|-----------|--------------|------|--------|
| **Peak Detection** | 1000 points | <1s | 50MB |
| **Voigt Fitting** | 10 peaks | 2-5s | 100MB |
| **Series Analysis** | 100 spectra | 1-2min | 500MB |
| **Batch Processing** | 1000 spectra | 10-20min | 2GB |

### **Optimization Tips**

- **Parallel Processing**: Enable for large datasets (`parallel=True`)
- **Memory Management**: Use streaming for very large series
- **Parameter Tuning**: Optimize fitting windows for your data
- **Caching**: Enable result caching for repeated analyses

---

## 🏆 Citation

If you use LunaNMR in your research, please cite:

```bibtex
@software{lunaNMR2025,
  title = {LunaNMR: Advanced NMR Peak Analysis and Integration Suite},
  author = {LunaNMR Contributors},
  year = {2025},
  version = {0.9},
  url = {https://github.com/mas-gu/lunaNMR},
  license = {MIT}
}
```

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
