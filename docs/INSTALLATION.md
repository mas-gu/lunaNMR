# Installation Guide - lunaNMR v0.9

This guide will walk you through installing the lunaNMR graphical application on your system. lunaNMR provides a user-friendly interface for NMR peak analysis without requiring programming knowledge.

## 📋 Table of Contents

1. [System Requirements](#system-requirements)
2. [Pre-Installation Checklist](#pre-installation-checklist)
3. [Installation Methods](#installation-methods)
4. [Dependency Installation](#dependency-installation)
5. [Verification](#verification)
6. [Troubleshooting](#troubleshooting)
7. [Platform-Specific Notes](#platform-specific-notes)

---

## 🖥️ System Requirements

### **Minimum Requirements**
- **Python:** 3.6 or higher (Python 3.8+ strongly recommended)
- **Architecture:** 64-bit (required for large datasets)
- **Memory:** 4 GB RAM minimum (8 GB+ recommended for series analysis)
- **Storage:** 500 MB free space (plus space for your NMR data)
- **Display:** 1024x768 minimum resolution for GUI

### **Recommended Configuration**
- **Python:** 3.9 or 3.10 (optimal compatibility)
- **Memory:** 16 GB RAM for large-scale analysis
- **CPU:** Multi-core processor (4+ cores) for parallel processing

### **Operating System Support**

| Platform | Status | Notes |
|----------|--------|-------|
| **macOS** | ✅ Fully Supported | Tested on macOS 10.14+ |
| **Linux** | ✅ Fully Supported | Ubuntu 18.04+, CentOS 7+ |
| **Windows** | ✅ Fully Supported | Windows 10+, some GUI limitations |

---

## ✅ Pre-Installation Checklist

### **1. Python Installation Verification**
```bash
# Check Python version
python3 --version
# Should show Python 3.6.0 or higher

# Check if pip is available
python3 -m pip --version
```

If Python is not installed or the version is too old:

#### **macOS (using Homebrew)**
```bash
# Install Homebrew if not present
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Python
brew install python3
```

#### **Linux (Ubuntu/Debian)**
```bash
# Update package list
sudo apt update

# Install Python and pip
sudo apt install python3 python3-pip python3-venv python3-tk
```

#### **Linux (CentOS/RHEL)**
```bash
# Install Python and development tools
sudo yum install python3 python3-pip python3-devel python3-tkinter
```

#### **Windows**
1. Download Python from [python.org](https://www.python.org/downloads/windows/)
2. **Important:** Check "Add Python to PATH" during installation
3. Select "Install for all users" if you have admin privileges

### **2. Tkinter Verification**
lunaNMR uses Tkinter for its GUI. Test if it's available:

```bash
python3 -c "import tkinter; print('Tkinter is available')"
```

If Tkinter is missing:
- **Linux:** Install `python3-tk` package
- **macOS/Windows:** Usually included with Python

---

## 🚀 Installation Methods

### **Method 1: Direct Download (Recommended)**

1. **Download the lunaNMR package**
   ```bash
   # If using Git
   git clone https://github.com/yourusername/lunaNMR.git
   cd lunaNMR/lunaNMR_v0o9

   # Or download and extract the ZIP file
   # Navigate to the lunaNMR_v0o9 directory
   ```

2. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

3. **Verify installation**
   ```bash
   python3 verify_installation.py
   ```

### **Method 2: Alternative Installation**

```bash
# Navigate to the lunaNMR directory
cd lunaNMR/lunaNMR_v0o9

# Install dependencies directly (system-wide)
pip install -r requirements.txt

# Test the installation
python3 launch_lunaNMR.py
```

---

## Dependency Installation

### **Core Dependencies**

lunaNMR requires several scientific Python packages. These are automatically installed with `pip install -r requirements.txt`:

```bash
# Core scientific computing
numpy>=1.20.0        # Numerical arrays and mathematical functions
pandas>=1.3.0        # Data manipulation and analysis
matplotlib>=3.5.0    # Plotting and visualization
scipy>=1.7.0         # Scientific computing and optimization

# Specialized libraries
networkx>=2.5        # Graph algorithms for peak clustering
nmrglue>=0.9.0       # NMR file format reading
psutil>=5.8.0        # System resource monitoring
scikit-learn>=1.0.0  # Machine learning algorithms

# GUI framework (usually included with Python)
# tkinter - Part of Python standard library
```

### **Optional Dependencies**

For enhanced functionality, you may also install:

```bash
# Jupyter notebook support
pip install jupyter

# Code formatting and linting
pip install black flake8

# Testing framework
pip install pytest

# Performance profiling
pip install line_profiler memory_profiler
```

### **Manual Dependency Installation**

If `requirements.txt` installation fails, install packages individually:

```bash
pip install numpy pandas matplotlib scipy
pip install networkx nmrglue psutil scikit-learn
```

---

##  Verification

### **1. Run the Installation Verification Script**

```bash
cd lunaNMR_v0o9
python3 verify_installation.py
```

**Expected Output:**
```
LunaNMR v0.9 - Installation Verification
=====================================

✅ Python version: 3.9.7 (compatible)
✅ All required packages are installed:
   - numpy: 1.26.4
   - pandas: 2.2.3
   - matplotlib: 3.9.2
   - scipy: 1.13.1
   - networkx: 3.2.1
   - nmrglue: 0.11
   - psutil: 5.9.0
   - scikit-learn: 1.5.2

✅ GUI framework (tkinter) is available
✅ All core modules can be imported successfully

🎉 Installation successful! You can now launch lunaNMR.
```

### **2. Test GUI Launch**

```bash
python3 launch_lunaNMR.py
```

This should open the lunaNMR Suite selector, allowing you to choose between:
- **LunaNMR**: Advanced NMR Peak Analysis and Integration
- **DynamiXs**: Dynamic Exchange Analysis (if available)

If the application selector appears, your installation is successful!

### **3. Test the GUI Application**

After installation, test that the GUI opens properly:

```bash
python3 launch_lunaNMR.py
```

You should see the lunaNMR Suite application selector. If the window opens and you can select "LunaNMR", your installation is complete and ready to use!

---

## 🔧 Troubleshooting

### **Common Installation Issues**

#### **Issue: "No module named tkinter"**
**Solution:**
```bash
# Linux
sudo apt install python3-tk

# macOS (if using Homebrew Python)
brew install python-tk

# Windows - Reinstall Python with Tkinter option checked
```

#### **Issue: "Permission denied" errors**
**Solution:**
```bash
# Use user installation instead of system-wide
pip install --user -r requirements.txt

# Or use virtual environment (recommended)
python3 -m venv lunaNMR_env
source lunaNMR_env/bin/activate
pip install -r requirements.txt
```

#### **Issue: nmrglue installation fails**
**Solution:**
```bash
# Install build dependencies first
# Linux:
sudo apt install python3-dev build-essential

# macOS:
xcode-select --install

# Then retry nmrglue installation
pip install nmrglue
```

#### **Issue: Slow performance or memory errors**
**Solution:**
- Ensure you have sufficient RAM (8GB+ recommended)
- Close other memory-intensive applications
- Use virtual environment to avoid package conflicts
- Consider upgrading to Python 3.9+ for better memory management

---

## 🖥️ Platform-Specific Notes

### **macOS**

#### **Apple Silicon (M1/M2) Macs**
```bash
# Use native Python (recommended)
# Download Python from python.org for Apple Silicon

# Or use Homebrew
arch -arm64 brew install python3

# Install dependencies
pip3 install -r requirements.txt
```

#### **Intel Macs**
```bash
# Standard installation works
brew install python3
pip3 install -r requirements.txt
```

#### **macOS Catalina and earlier**
- Python 2.7 is default; ensure you use `python3`
- May need to install Xcode Command Line Tools: `xcode-select --install`

### **Linux**

#### **Ubuntu/Debian**
```bash
# Complete development environment
sudo apt update
sudo apt install python3 python3-pip python3-venv python3-dev
sudo apt install python3-tk build-essential
sudo apt install libblas-dev liblapack-dev  # For NumPy optimization
```

#### **CentOS/RHEL**
```bash
# Enable EPEL repository for additional packages
sudo yum install epel-release
sudo yum install python3 python3-pip python3-devel
sudo yum install python3-tkinter gcc gcc-c++
```

#### **Arch Linux**
```bash
# Install Python and dependencies
sudo pacman -S python python-pip python-virtualenv tk
```

### **Windows**

#### **Windows 10/11**
1. **Install Python from python.org** (not Microsoft Store version)
2. **Check "Add Python to PATH"** during installation
3. **Install Microsoft C++ Build Tools** if compilation is needed
4. **Use Command Prompt or PowerShell** (not Git Bash for Python)

```cmd
# Verify installation
python --version
pip --version

# Install lunaNMR
pip install -r requirements.txt
```

#### **Windows-Specific Issues**
- **Long path names:** Enable long path support in Windows settings
- **Antivirus software:** May interfere with Python execution
- **User permissions:** Run Command Prompt as Administrator if needed

---

## 🔍 Advanced Configuration

### **Environment Variables**

You can set these environment variables for optimal performance:

```bash
# Linux/macOS
export OMP_NUM_THREADS=4           # Limit OpenMP threads
export MKL_NUM_THREADS=4           # Limit Intel MKL threads
export OPENBLAS_NUM_THREADS=4      # Limit OpenBLAS threads

# Add to ~/.bashrc or ~/.zshrc to make permanent
```

### **Memory Settings**

For large datasets:

```bash
# Increase Python recursion limit if needed
python3 -c "
import sys
print(f'Current recursion limit: {sys.getrecursionlimit()}')
sys.setrecursionlimit(10000)  # Increase if needed
"
```

---

**Next Steps:** Once installation is complete, proceed to the [Quick Start Guide](QUICKSTART.md) to learn how to use the lunaNMR graphical interface!
