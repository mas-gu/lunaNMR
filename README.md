# LunaNMR v1.0

> ⚠️ **Warning**
> This project is experimental and under active development. Expect instability, breaking API changes, incomplete features, and limited documentation. No support or warranty is provided. Use at your own risk.


**Advanced NMR Peak Analysis and Integration Suite**

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey.svg)](https://github.com/your-username/lunaNMR)

LunaNMR is a professional-grade software suite for Nuclear Magnetic Resonance (NMR) peak detection, Voigt profile fitting, and comprehensive analysis, now extended with binding-affinity (Kd / titration) fitting and protein dynamics analysis (T₁/T₂, methyl T₂, spectral density, and CPMG). A major focus of this release is **automation** — from automatic peak detection and overlap-aware multi-peak fitting through to end-to-end batch processing of large multi-spectrum datasets with minimal manual intervention.

---

## Screenshots

<div align="center">
 LunaNMR interface <br> 
  <img width="776" height="494" alt="luna_interface" src="https://github.com/user-attachments/assets/21e991e3-e622-4e60-8d32-c3b6ccd6d365" />
</div> 

<div align="center">
 LunaNMR serie integration <br> 
<img width="405" height="494" alt="serie_integration" src="https://github.com/user-attachments/assets/1b71dcea-e2a6-46d8-ad68-e29f8049ef71" /> 
</div> 

<div align="center">
 Individual-peak fitting vs volume integration <br> 
<img width="776" height="494" alt="individualpeak" src="https://github.com/user-attachments/assets/899ef2f6-3f1c-42aa-87c4-3b4b3702de66" />
</div> 

<div align="center">
 Multi-peaks fitting vs volume integration  <br> 
<img width="776" height="494" alt="multipeaks" src="https://github.com/user-attachments/assets/163800b4-081c-42e3-988f-ad015242c306" />
</div> 

---

## Key Features

###  **Advanced Peak Analysis**
- **Advanced 2D Voigt Fitting**: Multi-peak deconvolution with Gaussian and Lorentzian components
- **Multi-Peak Detection**: Automatic clustering and overlapping peak resolution
- **Quality Assessment**: R² values, uncertainties, and confidence intervals

###  **Comprehensive Analysis Workflows**
- **Single Spectrum Analysis**: Individual spectrum processing with full uncertainty analysis
- **Series Analysis**: T₁, T₂, and hetNOE relaxation studies with exponential fitting
- **Multi-Spectrum Processing**: Temperature, pH, and titration series analysis
- **Batch Processing**: Automated analysis of large datasets with parallel processing

###  **Automation**
- **Automatic Peak Detection**: Network-based clustering and prominence analysis with no manual seeding required
- **Overlap-Aware Fitting**: Automatic routing of closely-spaced peaks to 2D simultaneous multi-peak deconvolution
- **End-to-End Batch Pipelines**: CLI-driven processing of full multi-spectrum series, from detection to fitted results
- **Hands-Off Series Workflows**: Automated peak tracking, fitting, and output generation across titration, relaxation, and dynamics series

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

## Quick Start

### **Installation**

1. **Clone the repository:**
```bash
git clone https://github.com/your-username/lunaNMR.git
cd lunaNMR/lunaNMR_v1o0
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
---

## Package Structure

```
lunaNMR_v1o0/
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

---

##  Scientific Background

### **Voigt Profile Fitting**

LunaNMR implements professional-grade Voigt profile fitting, combining Gaussian and Lorentzian components to accurately model NMR lineshapes:

**Mathematical Foundation:**
- **Voigt Function**: Convolution of Gaussian (instrumental broadening) and Lorentzian (natural linewidth) functions
- **Multi-Peak Deconvolution**: Simultaneous fitting of overlapping peaks with constraints

**Key Algorithms:**
- **Peak Detection**: Network-based clustering with prominence analysis
- **Fitting Optimization**: Levenberg-Marquardt with robust parameter estimation
- **Uncertainty Analysis**: Bootstrap and covariance-based error estimation

---

##  Graphical Interface

### **Application Selector**

Launch `python3 launch_lunaNMR.py` to access the unified launcher:

---

---

## Documentation

### **Complete Documentation**

## Supported File Formats

### **Input Formats**

| Format | Extension | Support |
|--------|-----------|---------|
| **Bruker TopSpin** | `2ii`, `.2rr`  | ✅ Full |
| **Varian/Agilent** | `.fid`, `.ft`  | ✅ Full |
| **NMRPipe** | `.ft`, `.pipe` | ✅ Full |
| **SPARKY** | `.ucsf` | ✅ Full |

### **Export Formats**

- **Data**: CSV, Excel, JSON, 
- **Figures**: PNG, PDF, SVG, EPS

---

## Contributing

We welcome contributions from the NMR community! Please see our contributing guidelines:

### **Development Setup**

```bash
# Clone repository
git clone https://github.com/your-username/lunaNMR.git
cd lunaNMR/lunaNMR_v1o0

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

### **Optimization Tips**

- **Parallel Processing**: Enable for large datasets (`parallel=True`)
- **Memory Management**: Use streaming for very large series
- **Parameter Tuning**: Optimize fitting windows for your data
- **Caching**: Enable result caching for repeated analyses

---

## Citation

If you use LunaNMR in your research, please cite:

```bibtex
@software{lunaNMR2025,
  title = {LunaNMR: Advanced NMR Peak Analysis and Integration Suite},
  author = {LunaNMR Contributors},
  year = {2025},
  version = {1.0},
  url = {https://github.com/mas-gu/lunaNMR},
  license = {MIT}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
