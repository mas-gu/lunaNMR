# File Formats Guide - lunaNMR v0.9

This guide provides comprehensive information about all file formats supported by lunaNMR for input, output, and configuration. Understanding these formats is essential for preparing your data and interpreting results.

## 📋 Table of Contents

1. [Input File Formats](#input-file-formats)
2. [Output File Formats](#output-file-formats)
3. [Configuration Files](#configuration-files)
4. [Data Validation](#data-validation)
5. [Format Conversion](#format-conversion)
6. [Best Practices](#best-practices)

---

## Input File Formats

### **1. NMR Data Files**

lunaNMR uses [nmrglue](https://github.com/jjhelmus/nmrglue) to support a wide range of NMR data formats:

#### **Supported Formats**
- **Bruker TopSpin:** processed files (`1r`, `2rr`)
- **JEOL:** Delta format files
- **NMRPipe:** `.ft` files
- **Sparky:** `.ucsf` files

#### **Bruker Format Example**
```
your_experiment/
├── acqus          # Acquisition parameters
├── acqu2s         # 2D acquisition parameters (if applicable)
├── procs          # Processing parameters
├── proc2s         # 2D processing parameters (if applicable)
├── fid            # Raw FID data
├── 1r             # Processed 1D spectrum
└── 2rr            # Processed 2D spectrum
```

#### **File Requirements**
- **2D spectra:** Both dimensions must be properly calibrated
- **Processing:** Processed data preferred for peak detection

### **2. Peak List Files**

Peak lists define the positions of interest for analysis. Multiple formats are supported:

#### **Format : Simple Text Format**
```
Assignment, Position_X, Position_Y
1, 8, 110
2, 9, 112
3, 10, 111
```


#### **Field Specifications**
- **Assignment** Unique identifier for each peak
- **Position_X** ¹H chemical shift in ppm (float, typically 6.0-12.0)
- **Position_Y** ¹⁵N chemical shift in ppm (float, typically 100-140) [optional]


### **3. dynamiXs Data Files**

For time-series analysis (T₁, T₂, relaxation studies):

#### **Format: Multi-Column CSV**
```csv
# T1 Relaxation Series Data
# Time delays: 0.01, 0.05, 0.1, 0.3, 0.6, 0.9, 1.2, 1.8, 2.4 seconds
Residue,0.01,0.05,0.1,0.3,0.6,0.9,1.2,1.8,2.4
A5,1250.3,1098.2,967.4,745.1,532.8,424.6,315.2,198.7,125.4
L6,1456.7,1323.4,1189.2,961.5,742.3,628.9,417.1,302.8,201.2
V7,1123.8,1087.1,1034.5,823.2,651.7,498.3,387.9,245.1,156.7
I8,2234.1,2156.8,2087.3,1823.7,1456.2,1198.4,897.6,623.4,412.8
F9,1876.4,1798.2,1723.9,1456.7,1234.5,987.3,756.8,534.2,387.1
```

#### **Requirements**
- **First column:** Residue identifiers matching peak list
- **Header row:** Time points or experimental parameters
- **Data matrix:** Intensity values for each residue/time point
- **Missing data:** Use `NaN` or empty cells for missing values

---

## 📤 Output File Formats

### **1. Analysis Results**

#### **Standard CSV Output**
```csv
# lunaNMR Analysis Results - Generated on 2025-09-03
# Software: lunaNMR v0.9, Method: Enhanced Voigt Fitting
Index,Assignment,1H_center,1H_uncertainty,15N_center,15N_uncertainty,Amplitude,Amp_uncertainty,Gaussian_width,Gauss_uncertainty,Lorentzian_width,Lorentz_uncertainty,R_squared,Fit_quality,SNR,Processing_time
1,A5-HN,8.247,0.001,120.45,0.02,1523.4,45.2,0.018,0.002,0.012,0.001,0.987,Excellent,15.2,0.034
2,L6-HN,7.893,0.001,118.23,0.02,2187.6,52.8,0.016,0.001,0.014,0.001,0.992,Excellent,22.1,0.028
3,V7-HN,8.420,0.002,125.67,0.03,1834.7,48.9,0.020,0.002,0.011,0.001,0.989,Excellent,18.7,0.031
```


*Understanding these file formats will help you get the most out of lunaNMR's analysis capabilities. Always validate your input data and review output quality metrics for reliable results.*
