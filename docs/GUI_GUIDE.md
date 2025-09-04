# GUI User Guide - lunaNMR v0.9

Complete user guide for the lunaNMR graphical interface. This guide covers all features, menus, and workflows available in the user-friendly GUI application for NMR peak analysis.

## 📋 Table of Contents

1. [Application Launch](#application-launch)
2. [Main Interface](#main-interface)  
3. [Single Spectrum Analysis](#single-spectrum-analysis)
4. [Series Analysis](#series-analysis)
5. [Spectrum Browser](#spectrum-browser)
6. [Configuration Settings](#configuration-settings)
7. [Export Options](#export-options)
8. [Troubleshooting](#troubleshooting)

---

## 🚀 Application Launch

### **Starting lunaNMR**

1. **Open terminal/command prompt** and navigate to the lunaNMR folder:
```bash
cd lunaNMR_v0o9
```

2. **Launch the application:**
```bash
python3 launch_lunaNMR.py
```

3. **Application Selector Window:**
```
┌─────────────────────────────────────────┐
│           LunaNMR Suite v0.9            │
│                                         │
│  Select Application to Launch:          │
│                                         │
│  ● LunaNMR                              │
│    Advanced NMR Peak Analysis           │
│                                         │
│  ○ DynamiXs                             │
│    Dynamic Exchange Analysis            │
│                                         │
│         [Launch]    [Cancel]            │
└─────────────────────────────────────────┘
```

4. **Select "LunaNMR"** and click **[Launch]**

---

## 🖥️ Main Interface

### **Main Window Layout**

After launching, you'll see the main lunaNMR interface:

```
┌─────────────────────────────────────────────────────────────┐
│ lunaNMR v0.9                                      [X][□][-] │
├─────────────────────────────────────────────────────────────┤
│ File   Edit   Analysis   View   Tools   Help                │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  ┌─────────────────────┐  ┌─────────────────────┐           │
│  │   Single Spectrum   │  │   Series Analysis   │           │
│  │      Analysis       │  │                     │           │
│  └─────────────────────┘  └─────────────────────┘           │
│                                                             │  
│  ┌─────────────────────┐  ┌─────────────────────┐           │
│  │     Spectrum        │  │  ⚙️ Configuration    │           │
│  │      Browser        │  │                     │           │
│  └─────────────────────┘  └─────────────────────┘           │
│                                                             │
│  ┌─────────────────────┐                                    │
│  │   DynamiXs          │  [Recent Files...]                 │
│  │   Relaxation        │                                    │
│  └─────────────────────┘                                    │
└─────────────────────────────────────────────────────────────┘
```

### **Menu Bar Options**

#### **File Menu**
- **Open Spectrum...** - Load NMR data files
- **Open Peak List...** - Load peak list files  
- **Recent Files** - Access recently used files
- **Export Results...** - Save analysis results
- **Exit** - Close the application

#### **Edit Menu**
- **Undo/Redo** - Undo/redo recent actions
- **Copy Results** - Copy results to clipboard
- **Preferences...** - Open settings dialog

#### **Analysis Menu**
- **Single Spectrum** - Analyze individual spectra
- **Series Analysis** - Process time-series data
- **Batch Processing** - Analyze multiple files
- **Reset Analysis** - Clear current results

#### **View Menu**
- **Zoom In/Out** - Adjust plot zoom
- **Fit to Window** - Auto-scale plots
- **Show/Hide Panels** - Toggle interface elements
- **Full Screen** - Maximize analysis view

#### **Tools Menu**
- **Configuration** - Advanced settings
- **Spectrum Browser** - Interactive viewer
- **Data Validation** - Check data integrity
- **Performance Monitor** - System resource usage

#### **Help Menu**
- **Quick Start Guide** - Basic tutorial
- **User Manual** - Complete documentation
- **About lunaNMR** - Version and credits

---

## 📊 Single Spectrum Analysis

### **Opening Single Spectrum Analysis**

1. **Click "Single Spectrum Analysis"** from the main interface
2. **New window opens** with the single spectrum workflow

### **Analysis Window Layout**

```
┌───────────────────────────────────────────────────────────────────┐
│ Single Spectrum Analysis                                    [X][-] │
├───────────────────────────────────────────────────────────────────┤
│ ┌─ Data Loading ─────────┐ ┌─ Analysis Settings ──────────────┐   │
│ │ Peak List: [Browse...] │ │ Mode: ● Enhanced  ○ Professional │   │
│ │ NMR Data:  [Browse...] │ │                                  │   │
│ └────────────────────────┘ │ Prominence:  [0.1    ] ▲▼       │   │
│                            │ Window X:    [0.2    ] ▲▼       │   │
│ ┌─ Progress ─────────────┐ │ Window Y:    [1.5    ] ▲▼       │   │
│ │ Status: Ready          │ │ Baseline:    [Auto ▼]           │   │
│ │ [████████████████████] │ └──────────────────────────────────┘   │
│ └────────────────────────┘                                       │
│                            [Detect Peaks] [Reset] [Export]       │
├───────────────────────────────────────────────────────────────────┤
│ ┌─ Results Table ───────────────────────────────────────────────┐ │
│ │ Peak | 1H (ppm) | 15N (ppm) | Intensity |  R²   | Quality   │ │
│ │   1  │  8.247   │   120.5   │   1523    │ 0.987 │ Excellent │ │
│ │   2  │  7.891   │   118.2   │   2205    │ 0.952 │ Excellent │ │
│ │   3  │  8.415   │   125.7   │   1805    │ 0.889 │ Good      │ │
│ └───────────────────────────────────────────────────────────────┘ │
├───────────────────────────────────────────────────────────────────┤
│ ┌─ Spectrum Plot ───────────────────────────────────────────────┐ │
│ │     ^                                                         │ │
│ │  I  │  ∩     ∩                    ∩                          │ │
│ │  n  │ ╱ ╲   ╱ ╲                  ╱ ╲                         │ │
│ │  t  │╱   ╲ ╱   ╲                ╱   ╲                        │ │
│ │  e  │     ╲╱     ╲              ╱     ╲                       │ │
│ │  n  │              ╲──────────╱──      ╲───────────────────  │ │
│ │  s  │                                                        │ │
│ │  i  └────────────────────────────────────────────────────▶   │ │
│ │  t                    1H Chemical Shift (ppm)                │ │
│ │  y                9     8     7     6     5                  │ │
│ └───────────────────────────────────────────────────────────────┘ │
└───────────────────────────────────────────────────────────────────┘
```

### **Step-by-Step Workflow**

#### **Step 1: Load Data**
1. **Peak List File:**
   - Click **[Browse...]** next to "Peak List"
   - Select your peak list file (`.txt`, `.csv`)
   - Supported formats: Tab-delimited, comma-separated

2. **NMR Data File:**
   - Click **[Browse...]** next to "NMR Data"
   - Select your spectrum file (`.ft2`, `.fid`, `.ft`, `.ucsf`)
   - File format is detected automatically

#### **Step 2: Configure Analysis**

**Enhanced Mode (Recommended):**
- Automatic parameter selection
- Best for most users
- Click **Enhanced** radio button

**Professional Mode:**
- Manual parameter control
- For experienced users
- Adjustable parameters:
  - **Prominence:** Peak detection sensitivity (0.05-0.5)
  - **Window X:** Fitting window in 1H dimension (ppm)
  - **Window Y:** Fitting window in 15N dimension (ppm)  
  - **Baseline:** Correction method (Auto, Manual, None)

#### **Step 3: Run Analysis**
1. **Click [Detect Peaks]** to start analysis
2. **Monitor progress** in the progress bar
3. **Wait for completion** (typically 30 seconds to 2 minutes)

#### **Step 4: Review Results**

**Results Table:**
- **Peak:** Sequential peak number
- **1H (ppm):** Fitted 1H chemical shift
- **15N (ppm):** Fitted 15N chemical shift (if 2D)
- **Intensity:** Peak amplitude
- **R²:** Goodness of fit (0-1 scale)
- **Quality:** Color-coded quality assessment

**Quality Color Coding:**
- 🟢 **Excellent (R² > 0.95):** Green row
- 🟡 **Good (R² > 0.90):** Yellow row
- 🟠 **Fair (R² > 0.80):** Orange row  
- 🔴 **Poor (R² < 0.80):** Red row

**Interactive Plot:**
- **Blue line:** Original spectrum
- **Red markers:** Peak centers
- **Green curves:** Fitted Voigt profiles
- **Click peaks** for detailed information
- **Zoom/pan** with mouse controls

---

## 📈 Series Analysis

### **Opening Series Analysis**

1. **Click "Series Analysis"** from main interface
2. **Select experiment type:**
   - T₁ Relaxation
   - T₂ Relaxation
   - Heteronuclear NOE

### **Series Analysis Interface**

```
┌───────────────────────────────────────────────────────────────────┐
│ Series Analysis - T₁ Relaxation                            [X][-] │
├───────────────────────────────────────────────────────────────────┤
│ ┌─ Experiment Setup ────────────┐ ┌─ Fitting Options ───────────┐ │
│ │ Type: [T1 Relaxation ▼]       │ │ Model: [Monoexponential ▼] │ │
│ │ Data: [Browse...]             │ │ Method: [Robust ▼]         │ │
│ │ Format: [CSV ▼]               │ │ Bootstrap: ☑               │ │
│ └───────────────────────────────┘ └─────────────────────────────┘ │
│                                                                   │
│ ┌─ Time Points ─────────────────────────────────────────────────┐ │
│ │ Delay (s): 0.01  0.1  0.3  0.6  1.0  1.5  2.0  3.0  5.0     │ │
│ │ Points:     ☑    ☑    ☑    ☑    ☑    ☑    ☑    ☑    ☑       │ │
│ └───────────────────────────────────────────────────────────────┘ │
│                                                                   │
│         [Process Series] [Reset] [Export Results]                │
├───────────────────────────────────────────────────────────────────┤
│ ┌─ Relaxation Results ──────────────────────────────────────────┐ │
│ │ Residue | T₁ (s) | Uncertainty | R² | Quality | Notes        │ │
│ │  G10    │  1.23  │    0.05     │0.98│Excellent│ ✅           │ │
│ │  A11    │  1.45  │    0.08     │0.95│Excellent│ ✅           │ │
│ │  L12    │  0.89  │    0.12     │0.87│Good     │ ✅           │ │
│ └───────────────────────────────────────────────────────────────┘ │
├───────────────────────────────────────────────────────────────────┤
│ ┌─ Relaxation Curves ───────────────────────────────────────────┐ │
│ │   ^                                                           │ │
│ │   │ ●                                                         │ │
│ │ I │  ╲●                                                       │ │
│ │ n │   ╲                                                       │ │
│ │ t │    ●╲                                                     │ │
│ │ e │     ╲●                                                    │ │
│ │ n │      ╲                                                    │ │
│ │ s │       ●──●───●                                            │ │
│ │ i │                                                           │ │
│ │ t │                                                           │ │
│ │ y └─────────────────────────────────────────────────────────▶ │ │
│ │             Time Delay (s)                                    │ │
│ │     0     1     2     3     4     5                           │ │
│ └───────────────────────────────────────────────────────────────┘ │
└───────────────────────────────────────────────────────────────────┘
```

### **Series Analysis Workflow**

#### **Step 1: Experiment Setup**
1. **Select Type:** Choose T₁, T₂, or hetNOE from dropdown
2. **Load Data:** Browse for your series data file
3. **Format:** Specify file format (CSV, Excel, JSON)

#### **Step 2: Configure Fitting**
1. **Model Selection:**
   - **Monoexponential:** Single exponential decay
   - **Biexponential:** Two-component decay
   - **Stretched:** Non-exponential relaxation

2. **Fitting Method:**
   - **Robust:** Outlier-resistant fitting
   - **Standard:** Least squares fitting
   - **Bootstrap:** Error estimation via resampling

#### **Step 3: Time Point Selection**
- **Enable/disable** individual time points
- **Review delay values** for accuracy
- **Quality control** for problematic points

#### **Step 4: Process and Review**
1. **Click [Process Series]** to start fitting
2. **Monitor progress** for each residue
3. **Review results table** for T₁/T₂ values and uncertainties
4. **Inspect curves** in the plot panel

---

## 🔍 Spectrum Browser

### **Opening Spectrum Browser**

1. **Click "Spectrum Browser"** from main interface
2. **Interactive spectrum viewer** opens

### **Browser Interface Features**

```
┌───────────────────────────────────────────────────────────────────┐
│ Spectrum Browser                                            [X][-] │
├───────────────────────────────────────────────────────────────────┤
│ File: spectrum.ft2 │ [Open...] [Save...] │ Zoom: [100%▼] [Fit]    │
├───────────────────────────────────────────────────────────────────┤
│ ┌─ Navigation ──┐ ┌─ Peak Tools ───────┐ ┌─ Display ─────────┐   │
│ │ ◀◀ ◀ ■ ▶ ▶▶ │ │ + Pick  🎯 Center  │ │ Grid: ☑  Labels: ☑│   │
│ │              │ │ - Delete ✏ Edit    │ │ Scale: [Auto ▼]   │   │
│ └──────────────┘ └────────────────────┘ └───────────────────┘   │
├───────────────────────────────────────────────────────────────────┤
│ ┌─ Main Spectrum View ──────────────────────────────────────────┐ │
│ │  15N                                                          │ │
│ │   ^                                                           │ │
│ │   │   ●                                                       │ │
│ │120│      ●                                                    │ │
│ │   │         ●                                                 │ │
│ │   │            ●                                              │ │
│ │110│               ●                                           │ │
│ │   │                                                           │ │
│ │   └─────────────────────────────────────────────────────────▶ │ │
│ │     9.0    8.5    8.0    7.5    7.0    6.5         1H       │ │
│ │                                                               │ │
│ │  [Crosshair: 8.25, 120.5 ppm] [Intensity: 1523]            │ │
│ └───────────────────────────────────────────────────────────────┘ │
├───────────────────────────────────────────────────────────────────┤
│ ┌─ Peak List ─────────────────────────────────────────────────┐   │
│ │ # │ 1H(ppm) │15N(ppm)│ Int  │ Vol  │ Note                  │   │
│ │ 1 │  8.25   │ 120.5  │ 1523 │ 2.1e6│ G10-HN                │   │
│ │ 2 │  7.89   │ 118.2  │ 2205 │ 3.2e6│ A11-HN                │   │
│ │ 3 │  8.42   │ 125.7  │ 1805 │ 2.8e6│ L12-HN                │   │
│ └─────────────────────────────────────────────────────────────────┘ │
└───────────────────────────────────────────────────────────────────┘
```

### **Interactive Features**

#### **Navigation Controls**
- **◀◀ ◀ ■ ▶ ▶▶:** Navigate through spectral regions
- **Mouse wheel:** Zoom in/out
- **Click and drag:** Pan across spectrum
- **Double-click:** Center on position

#### **Peak Tools**
- **+ Pick:** Click to add new peak at cursor position
- **- Delete:** Remove selected peak
- **🎯 Center:** Auto-center peak at maximum
- **✏ Edit:** Modify peak parameters manually

#### **Display Options**
- **Grid:** Toggle coordinate grid
- **Labels:** Show/hide peak labels
- **Scale:** Adjust intensity scaling
- **Contours:** Modify contour levels (2D spectra)

#### **Peak List Management**
- **Click rows** to highlight peaks in spectrum
- **Edit values** directly in the table
- **Add notes** and assignments
- **Export peak list** to various formats

---

## ⚙️ Configuration Settings

### **Opening Configuration**

1. **Click "Configuration"** from main interface or Tools menu
2. **Settings dialog** opens with multiple tabs

### **Configuration Tabs**

```
┌─────────────────────────────────────────────────────────────┐
│ Configuration Settings                                [X]   │
├─────────────────────────────────────────────────────────────┤
│ [General] [Analysis] [Display] [Export] [Advanced]         │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Analysis Parameters ─────────────────────────────────────┐│
│ │                                                           ││
│ │ Peak Detection:                                           ││
│ │   Prominence Threshold: [0.10] (0.01-0.5)               ││
│ │   Minimum Distance:     [20  ] points                    ││
│ │   SNR Threshold:        [3.0 ] ratio                     ││
│ │                                                           ││
│ │ Fitting Parameters:                                       ││
│ │   Max Iterations:       [1000] iterations                ││
│ │   Tolerance:            [1e-6] convergence               ││
│ │   Window 1H:            [0.2 ] ppm                       ││
│ │   Window 15N:           [1.5 ] ppm                       ││
│ │                                                           ││
│ │ Baseline Correction:                                      ││
│ │   Method:               [ArPLS ▼]                        ││
│ │   Lambda:               [Auto ▼]                         ││
│ │   Iterations:           [10  ]                           ││
│ │                                                           ││
│ └───────────────────────────────────────────────────────────┘│
├─────────────────────────────────────────────────────────────┤
│                           [Apply] [Reset] [OK] [Cancel]     │
└─────────────────────────────────────────────────────────────┘
```

### **Configuration Sections**

#### **General Tab**
- **Default Directories:** Set default paths for data and results
- **Recent Files:** Configure recent files list size
- **Backup Settings:** Enable automatic backups
- **Update Preferences:** Check for software updates

#### **Analysis Tab**  
- **Peak Detection:** Prominence threshold, distance, SNR
- **Fitting Parameters:** Iterations, tolerance, windows
- **Baseline Correction:** Method selection and parameters
- **Quality Thresholds:** R² limits for quality classification

#### **Display Tab**
- **Plot Settings:** Colors, fonts, line styles
- **Table Format:** Decimal places, column visibility
- **Interface Theme:** Light, dark, or system theme
- **Language:** Interface language selection

#### **Export Tab**
- **Default Format:** CSV, Excel, JSON preferences
- **Figure Settings:** DPI, size, file format for plots
- **Data Precision:** Number of decimal places in exports
- **Metadata:** Include analysis parameters in exports

#### **Advanced Tab**
- **Performance:** Parallel processing, memory limits
- **Debugging:** Enable debug logs, verbose output
- **Experimental:** Beta features and advanced options
- **Plugin Settings:** Configure installed plugins

### **Saving Configurations**

1. **Project Configurations:** Save/load settings for specific projects
2. **User Defaults:** Set global default preferences
3. **Configuration Files:** Export/import settings as JSON files
4. **Reset Options:** Restore factory defaults

---

## 📁 Export Options

### **Accessing Export Functions**

Export options are available from:
- **File → Export Results...** (main menu)
- **[Export]** button in analysis windows
- **Right-click → Export** in results tables

### **Export Dialog**

```
┌─────────────────────────────────────────────────────────────┐
│ Export Results                                        [X]   │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Export Type ─────────────────────────────────────────────┐│
│ │ ○ Data Only          ○ Figures Only     ● Both           ││
│ └───────────────────────────────────────────────────────────┘│
│                                                             │
│ ┌─ Data Export ─────────────────────────────────────────────┐│
│ │ Format:    [CSV ▼]                                       ││
│ │ Include:   ☑ Peak positions    ☑ Intensities            ││
│ │            ☑ Uncertainties     ☑ Quality metrics        ││
│ │            ☑ Analysis params   ☑ Metadata               ││
│ │                                                           ││
│ │ Filename:  [results_2024.csv                        ]    ││
│ │ Location:  [/Users/name/nmr_data/                  ...]   ││
│ └───────────────────────────────────────────────────────────┘│
│                                                             │
│ ┌─ Figure Export ───────────────────────────────────────────┐│
│ │ Format:    [PNG ▼]     DPI: [300 ▼]                     ││
│ │ Size:      [12x8] inches                                 ││
│ │ Include:   ☑ Spectrum plot     ☑ Fitting curves         ││
│ │            ☑ Results table     ☑ Quality summary        ││
│ │                                                           ││
│ │ Filename:  [analysis_figure.png                     ]    ││
│ │ Location:  [/Users/name/figures/                   ...]   ││
│ └───────────────────────────────────────────────────────────┘│
├─────────────────────────────────────────────────────────────┤
│                     [Preview] [Export] [Cancel]            │
└─────────────────────────────────────────────────────────────┘
```

### **Export Format Options**

#### **Data Formats**
- **CSV:** Comma-separated values (Excel compatible)
- **Excel:** Native Excel format with multiple sheets
- **JSON:** Structured data with metadata
- **TSV:** Tab-separated values (origin/Igor compatible)
- **HDF5:** High-performance format for large datasets

#### **Figure Formats**
- **PNG:** Web/presentation graphics (300+ DPI recommended)
- **PDF:** Vector format for publications
- **SVG:** Scalable vector graphics  
- **EPS:** Encapsulated PostScript for journals
- **TIFF:** High-quality raster format

### **Export Content Options**

#### **Data Content**
- **Peak Positions:** Chemical shifts in all dimensions
- **Intensities:** Peak amplitudes and volumes
- **Uncertainties:** Bootstrap or covariance-based errors
- **Quality Metrics:** R² values, fit quality classifications
- **Analysis Parameters:** Settings used for analysis
- **Metadata:** Timestamps, software version, file info

#### **Figure Content**
- **Spectrum Plot:** Original data with fits overlayed
- **Fitting Curves:** Individual peak fits with residuals
- **Results Table:** Formatted table of peak parameters
- **Quality Summary:** Statistics and quality distributions

### **Batch Export**

For multiple analyses:

1. **Select multiple results** in the results browser
2. **Choose "Batch Export"** from the export menu
3. **Configure naming** pattern for files
4. **Set common export** parameters
5. **Process all selected** analyses automatically

---

## 🔧 Troubleshooting

### **Common GUI Issues**

#### **Application Won't Launch**
**Symptoms:** Double-clicking does nothing, or error messages appear

**Solutions:**
1. **Check Python Installation:**
```bash
python3 --version
# Should show Python 3.8 or higher
```

2. **Verify Dependencies:**
```bash
cd lunaNMR_v0o9
python3 lunaNMR/validation/verify_installation.py
```

3. **Launch from Command Line:**
```bash
python3 launch_lunaNMR.py
# Check for error messages
```

#### **GUI Appears Corrupted or Blank**
**Symptoms:** Missing buttons, blank windows, corrupted display

**Solutions:**
1. **Check Display Settings:** Ensure adequate screen resolution (1024x768+)
2. **Update Graphics Drivers:** Especially on Windows systems
3. **Try Different Themes:** Go to Configuration → Display → Theme
4. **Reset GUI Settings:** Delete user configuration files

#### **Files Won't Load**
**Symptoms:** "File not found" or "Format not supported" errors

**Solutions:**
1. **Check File Path:** Ensure no special characters or spaces in path
2. **Verify File Format:** Use supported formats (see File Formats guide)
3. **Check File Permissions:** Ensure files are readable
4. **Try Different Location:** Copy files to desktop and try again

#### **Analysis Hangs or Crashes**
**Symptoms:** Progress bar stops, application becomes unresponsive

**Solutions:**
1. **Check System Resources:** Close other applications
2. **Reduce Dataset Size:** Try with smaller test dataset
3. **Adjust Parameters:** Use less aggressive fitting settings
4. **Enable Debug Mode:** Check logs for specific errors

### **Analysis Quality Issues**

#### **Poor Peak Detection**
**Symptoms:** Missing peaks, false positive peaks

**Solutions:**
1. **Adjust Prominence:** Lower for more peaks, raise for fewer
2. **Check Baseline:** Poor baseline affects peak detection
3. **Review SNR:** Low signal-to-noise requires parameter adjustment
4. **Manual Peak Picking:** Use Spectrum Browser for difficult cases

#### **Bad Fitting Results**
**Symptoms:** Low R² values, unrealistic parameters

**Solutions:**
1. **Check Peak List:** Ensure positions match spectrum
2. **Adjust Fitting Windows:** Wider windows for overlapped peaks
3. **Review Baseline:** Poor baseline correction affects fitting
4. **Try Different Model:** Switch between fitting methods
5. **Manual Inspection:** Use Spectrum Browser to verify issues

#### **Inconsistent Results**
**Symptoms:** Results vary between runs, unreproducible outcomes

**Solutions:**
1. **Fix Random Seeds:** Enable deterministic fitting in Advanced settings
2. **Increase Iterations:** Use more iterations for convergence
3. **Check Data Quality:** Verify input data integrity
4. **Document Parameters:** Keep consistent analysis settings

### **Performance Issues**

#### **Slow Analysis**
**Symptoms:** Long processing times, unresponsive interface

**Solutions:**
1. **Enable Parallel Processing:** Use multiple CPU cores
2. **Optimize Parameters:** Reduce unnecessary iterations
3. **Check System Load:** Close other CPU-intensive applications
4. **Increase RAM:** Consider system memory upgrade for large datasets

#### **Memory Errors**
**Symptoms:** "Out of memory" messages, system slowdown

**Solutions:**
1. **Close Other Applications:** Free up system memory
2. **Process Smaller Batches:** Break large analyses into chunks
3. **Enable Streaming:** Use memory-efficient processing modes
4. **System Upgrade:** Consider adding more RAM

### **Getting Additional Help**

#### **Log Files**
Check log files for detailed error information:
- **Location:** `lunaNMR_v0o9/logs/`
- **Files:** `debug.log`, `error.log`, `analysis.log`

#### **Debug Mode**
Enable debug mode for detailed troubleshooting:
1. **Configuration → Advanced → Debugging**
2. **Enable "Verbose Output"**
3. **Check console output** for detailed messages

#### **Community Support**
- **Documentation:** Check complete docs in `docs/` folder
- **GitHub Issues:** Report bugs and feature requests
- **User Forums:** Community discussions and solutions
- **Email Support:** Contact developers for critical issues

---

This comprehensive GUI guide covers all aspects of using lunaNMR through its graphical interface. The application is designed to be intuitive for users without programming experience while providing powerful analysis capabilities for advanced NMR research.
