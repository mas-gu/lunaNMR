# DynamiXs GUI v2.0 - NMR Relaxation Analysis Suite

A comprehensive graphical user interface for NMR relaxation data analysis built with Python and PySide6 (Qt).

## Features

### 1. T1/T2 Fitting Analysis
- **Single and multi-core processing** for exponential decay fitting
- **Bootstrap error estimation** for reliable parameter uncertainties
- **Interactive parameter adjustment** (initial amplitude, time constant, bootstrap iterations)
- **Automatic plot generation** with customizable layouts
- **Real-time progress monitoring** and results display

### 2. Data Plotting Tools
- **Dataset comparison** with difference plots and statistical analysis
- **Dual field plotting** for multi-field NMR experiments  
- **Interactive parameter selection** for flexible data visualization
- **PDF output** with publication-ready formatting

### 3. Cross-Platform Compatibility
- **macOS, Windows, and Linux** support through PySide6 (Qt)
- **Modern Qt-based interface** with LunaNMR styling
- **User-friendly interface** requiring no scripting knowledge

## Quick Start

1. **Launch the GUI:**
   ```bash
   cd dynamiXs/
   python run_dynamixs_gui.py
   ```

2. **Select your analysis type** from the main menu

3. **Import your CSV data files** using the browse buttons

4. **Configure analysis parameters** through the interface

5. **Run analysis** and view results in real-time

## File Format Requirements

### T1/T2 Data Files
- **CSV format** with comma separation
- **First row:** Time delays (headers: t1, t2, t3, ...)
- **First column:** Residue names (rows 2+)
- **Data matrix:** Signal intensities for each residue and time point

Example:
```
,0.05,0.1,0.3,0.6,0.9,1.2,1.8,2.4
A5,125.3,98.2,67.4,45.1,32.8,24.6,15.2,9.8
L6,156.7,123.4,89.2,61.5,42.3,28.9,17.1,10.2
```

## Analysis Modules

### T1/T2 Fitting
- **Input:** CSV file with relaxation decay data
- **Processing:** Single-core or multi-core exponential fitting
- **Output:** 
  - Results text file with fitted parameters and errors
  - Multi-page PDF plots showing data and fits
  - Statistical summaries

### Dataset Comparison
- **Input:** Two CSV files with analysis results
- **Processing:** Statistical comparison and difference calculation
- **Output:**
  - Comparison plots (original data + differences)
  - Statistical summaries (mean, std, range)

## Mathematical Framework

### Model-Free Analysis Equations

The dynamiXs package implements dual-field model-free analysis using the Lipari-Szabo approach with chemical exchange (Rex) fitting. The core equations are:

#### NMR Relaxation Rate Equations
**Longitudinal relaxation (R₁):**
```
R₁ = A × (3J(ωN) + 6J(ωH + ωN) + J(|ωH - ωN|)) + C × J(ωN)
```

**Transverse relaxation (R₂):**
```
R₂ = 0.5 × A × (4J(0) + 3J(ωN) + 6J(ωH + ωN) + 6J(0.87ωH) + J(|ωH - ωN|))
   + C × (2J(0)/3 + 0.5J(ωN)) + Rex
```

**Heteronuclear NOE:**
```
NOE = 1 + (A × γH/γN × T₁ × (6J(ωH + ωN) - J(|ωH - ωN|)))
```

#### Physical Constants
```
A = d² = (1/4) × [μ₀ℏγNγH/(4πr³NH)]²    (Dipolar coupling constant)
C = c² = (1/3) × (ωN × Δσ)²              (CSA constant)
```

**Key Parameters:**
- `γH = 2.675 × 10⁸ rad s⁻¹ T⁻¹` (¹H gyromagnetic ratio)
- `γN = -2.713 × 10⁷ rad s⁻¹ T⁻¹` (¹⁵N gyromagnetic ratio)
- `rNH = 1.015 × 10⁻¹⁰ m` (N-H bond length)
- `Δσ = -160 ppm` (¹⁵N CSA)
- `μ₀ = 1.257 × 10⁻⁶ H/m` (Permeability of vacuum)
- `ℏ = 1.055 × 10⁻³⁴ J·s` (Reduced Planck constant)

#### Model-Free Spectral Density Functions
**Lipari-Szabo Model:**
```
J(ω) = (2/5) × [S² × τc/(1 + (ωτc)²) + (1-S²) × τeff/(1 + (ωτeff)²)]
```

Where:
- `τeff = 1/(1/τc + 1/τe)` (Effective correlation time)
- `S²` = Generalized order parameter (0-1)
- `τc` = Overall tumbling correlation time
- `τe` = Internal motion correlation time
- `Rex` = Chemical exchange contribution to R₂

#### Enhanced J(0.87ωH) Implementation
For improved accuracy, the analysis uses `J(0.87ωH)` instead of `J(ωH)` to better account for cross-correlation effects and dipolar-CSA interference.

#### Dual-Field Analysis
The dual-field approach provides better parameter separation by using data from two magnetic fields:
- **Field 1** (e.g., 600 MHz): Lower chemical exchange sensitivity
- **Field 2** (e.g., 700-800 MHz): Higher chemical exchange sensitivity

**Expected Rex scaling:** `Rex₂/Rex₁ ≈ (B₂/B₁)²`

#### Reduced Spectral Density Mapping
Alternative analysis using direct calculation (Farrow et al. 1995):
```
J(0) = (3/[2((3d²/4) + c²)]) × (-0.5R₁ + R₂ - (3/5)σNOE)
J(ωN) = (1/((3d²/4) + c²)) × (R₁ - (7/5)σNOE)
J(0.87ωH) = [4/(5d²)] × σNOE
```

Where:
- `σNOE = (NOE - 1) × R₁ × (γN/γH)`
- `d² = [μ₀ħγNγH/(8π²r³NH)]²` (dipolar coupling constant)
- `c² = (ωN × Δσ)²/3` (CSA constant)

**Reference:** Farrow et al. J. Biomol. NMR, 6 (1995) 153-162, Equations 5-7

## Technical Details

### Dependencies
```
- Python 3.9+
- PySide6
- numpy
- matplotlib
- pandas
- lmfit
- scipy
```

Install with: `pip install PySide6 numpy matplotlib pandas lmfit scipy`

### Architecture
- **Main GUI:** `dynamiXs_gui.py` - Central Qt interface controller
- **GUI Components:** `gui_components.py` - Reusable Qt widgets and factories
- **Design Constants:** `constants.py` - LunaNMR styling constants
- **Stylesheets:** `styles/main.qss` - Qt stylesheet for consistent theming
- **Workers:** `workers.py` - QThread workers for background analysis
- **Visualization:** `visualization/` - Results and fit viewers
- **Launcher:** `run_dynamixs_gui.py` - Application entry point

### Performance
- **Multi-core support** for T1/T2 fitting using all available CPU cores
- **Threading** for non-blocking GUI during long calculations
- **Memory efficient** processing for large datasets
- **Progress tracking** for user feedback

## File Organization

```
dynamiXs_v2o0/
├── run_dynamixs_gui.py       # Application entry point
├── dynamiXs_gui.py           # Main Qt GUI application
├── gui_components.py         # Reusable Qt widgets
├── constants.py              # Design tokens and colors
├── workers.py                # QThread background workers
├── README.md                 # This file
├── styles/
│   └── main.qss              # Qt stylesheet
├── visualization/
│   ├── results_viewer.py     # Model-free results viewer
│   └── fit_viewer.py         # T1/T2 fit viewer
├── dynamiXs_T1_T2/           # T1/T2 fitting modules
│   ├── fit_Tx_NMRRE.py       # Single-core fitting
│   └── fitmulti__Tx_NMRRE.py # Multi-core fitting
├── dynamiXs_integrated/      # Integrated analysis pipeline
├── dynamiXs_density_functions/ # Spectral density analysis
├── dynamiXs_cpmg/            # CPMG analysis modules
├── dynamiXs_plot/            # Plotting utilities
└── dynamiXs_format/          # Data formatting tools
```

## Current Implementation Status

✅ **Completed:**
- Main GUI with 5-option menu
- T1/T2 fitting interface (single and multi-core)
- Dataset comparison interface
- Dual field plotting interface
- File handling and parameter configuration
- Real-time progress monitoring

🚧 **In Development:**
- Spectral density analysis module
- CPMG analysis integration
- Data formatting utilities

## Usage Notes

- **Working Directory:** The GUI operates in your current working directory for file I/O
- **File Validation:** Input files are checked for existence and format
- **Error Handling:** Comprehensive error messages and user feedback
- **Threading:** Long calculations run in background threads to keep GUI responsive

## Troubleshooting

### GUI Won't Start
- Ensure PySide6 is installed: `pip install PySide6`
- Check Python version: requires 3.9+
- Verify all dependencies are installed

### Analysis Fails
- Check input file format matches requirements
- Ensure file paths are accessible
- Review error messages in results panel

### Performance Issues
- Use multi-core option for large datasets
- Close other applications to free memory
- Consider reducing bootstrap iterations for faster results

### Threading Warnings
- Matplotlib threading warnings have been resolved by using the Agg backend
- All analysis scripts now use non-interactive plotting for thread safety
- Generated plots are saved as PDF files without requiring GUI interaction

## Future Enhancements

- **CPMG relaxation dispersion** fitting
- **Batch processing** for multiple datasets
- **Parameter optimization** suggestions
- **Export to common formats** (Excel, Origin, etc.)
