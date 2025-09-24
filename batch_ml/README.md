# LunaNMR Batch Processing Scripts

This folder contains two scripts for batch processing NMR spectra to generate ML training data. Both scripts are completely independent of the lunaNMR GUI interface and process multiple spectra automatically.

## 📁 Available Scripts

### 1. **quick_batch.py** - Simple & Fast

**Best for:** Quick processing with minimal setup

```bash
python quick_batch.py /path/to/spectra
python quick_batch.py /path/to/spectra 15N1H    # specify nucleus type
python quick_batch.py /path/to/spectra auto     # enable optimization
```

#### ✅ Advantages:
- **Minimal arguments** - Just folder path required
- **Fast to run** - No complex setup needed
- **Simple output** - Clean, easy-to-read results
- **Perfect for beginners** - Straightforward usage
- **Sensible defaults** - Good settings out-of-the-box

#### ❌ Disadvantages:
- **Limited customization** - Few configuration options
- **Basic error reporting** - Less detailed diagnostics
- **No advanced features** - No presets or config files

---

### 2. **batch_nmr_process.py** - Full-Featured

**Best for:** Advanced users who need maximum control

```bash
python batch_nmr_process.py /path/to/spectra
python batch_nmr_process.py /path/to/spectra --nucleus 15N1H --optimize
python batch_nmr_process.py /path/to/spectra --preset conservative --config custom.json
```

#### ✅ Advantages:
- **Full CLI interface** - Complete command-line options
- **Advanced configuration** - Presets, config files, detailed settings
- **Comprehensive logging** - Detailed progress and error reporting
- **Maximum flexibility** - Fine-tune every parameter
- **Professional features** - Dry-run, validation, extensive options

#### ❌ Disadvantages:
- **More complex** - Many command-line options to learn
- **Slower setup** - Requires understanding of various parameters
- **Can be overwhelming** - Many options for simple use cases

---

## 🚀 Quick Start - Recommended Usage

### For Your 100 Spectra Folder:

**Option 1 - Simplest (Recommended for first run):**
```bash
cd lunaNMR_v0o9/batch_ml
python quick_batch.py /path/to/your/100/spectra auto
```

**Option 2 - Advanced (If you need more control):**
```bash
cd lunaNMR_v0o9/batch_ml
python batch_nmr_process.py /path/to/your/100/spectra --optimize --preset aggressive
# Use --preset aggressive for maximum ML training data collection
# See "Available Presets" section below for all 6 preset options
```

---

## 📊 Expected Results

Both scripts will:
- ✅ **Auto-detect nucleus types** (15N1H, 13C1H) from filenames
- ✅ **Optimize S/N thresholds** automatically for maximum peak detection
- ✅ **Process all spectra** in the folder without manual intervention
- ✅ **Collect ML training data** from successfully fitted peaks
- ✅ **Continue on errors** - won't stop if individual spectra fail
- ✅ **Provide progress updates** and final summary

**Sample Output:**
```
🚀 Starting batch processing...
   Folder: /data/spectra
   Nucleus: auto-detect
   Optimization: enabled

Progress: 20.0% (20/100)
Progress: 40.0% (40/100)
...

✅ BATCH PROCESSING COMPLETED
📁 Total files: 100
✅ Successfully processed: 95
❌ Failed: 5
🧬 Peaks detected: 2,847
📊 ML samples collected: 1,923
📈 Success rate: 95.0%

🎉 SUCCESS: 1,923 ML training samples collected!
```

---

## 🛠️ Script Comparison

| Feature | quick_batch.py | batch_nmr_process.py |
|---------|---------------|---------------------|
| **Ease of use** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
| **Setup time** | ⭐⭐⭐⭐⭐ | ⭐⭐ |
| **Customization** | ⭐⭐ | ⭐⭐⭐⭐⭐ |
| **Error reporting** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **Advanced features** | ⭐ | ⭐⭐⭐⭐⭐ |
| **Best for beginners** | ✅ Yes | ❌ No |
| **Best for experts** | ❌ No | ✅ Yes |

---

## 🎯 Recommendation

**Start with `quick_batch.py`** for your first run - it will handle your 100 spectra perfectly with minimal fuss. If you need more control or advanced features later, switch to `batch_nmr_process.py`.

Both scripts are completely independent and safe to use alongside the lunaNMR GUI without any interference.

---

## ⚙️ **Available Presets - Complete Guide**

The `batch_nmr_process.py` script supports **6 powerful presets** for different NMR analysis scenarios. Each preset optimizes parameters for specific conditions and data types.

### **🧬 Nucleus-Specific Presets**

#### **1. `15N1H` - 15N-1H HSQC Optimized**
```bash
python batch_nmr_process.py /path/to/spectra --preset 15N1H
```
**Best for:** 15N-1H HSQC/HMQC protein spectra
- **S/N Threshold:** 2.5 (optimized for HSQC sensitivity)
- **Expected Peaks:** ~60 (typical protein backbone)
- **Peak Range:** 30-80 peaks (conservative protein estimate)
- **Quality Threshold:** 0.8 (standard high-quality)
- **Use Case:** Protein backbone assignment, structure studies

#### **2. `13C1H` - 13C-1H HSQC Optimized**
```bash
python batch_nmr_process.py /path/to/spectra --preset 13C1H
```
**Best for:** 13C-1H HSQC/HMQC spectra
- **S/N Threshold:** 2.0 (carbon detection sensitivity)
- **Expected Peaks:** ~40 (typical carbon count)
- **Peak Range:** 20-60 peaks (carbon-optimized)
- **Quality Threshold:** 0.8 (standard high-quality)
- **Use Case:** Carbon assignment, metabolomics, small molecules

#### **3. `1H` - 1H NMR Optimized**
```bash
python batch_nmr_process.py /path/to/spectra --preset 1H
```
**Best for:** 1D 1H NMR spectra
- **S/N Threshold:** 2.0 (proton sensitivity)
- **Expected Peaks:** ~30 (typical 1H spectrum)
- **Peak Range:** 15-40 peaks (1H-optimized)
- **Quality Threshold:** 0.8 (standard high-quality)
- **Use Case:** 1D proton NMR, mixture analysis, quantification

---

### **🎯 Strategy-Based Presets**

#### **4. `conservative` - High-Quality Data Only**
```bash
python batch_nmr_process.py /path/to/spectra --preset conservative
```
**Best for:** High-quality datasets, publication-ready results
- **S/N Thresholds:** Higher across all nuclei (15N: 3.0, 13C: 2.5, 1H: 2.5)
- **Quality Threshold:** 0.9 (very high quality requirement)
- **Strategy:** Conservative optimization approach
- **Result:** Fewer peaks but higher confidence
- **Use Case:** Final analysis, publication figures, validation studies

**Advantages:**
- ✅ **Highest Data Quality** - Only excellent fits collected
- ✅ **Reduced False Positives** - Stringent S/N filtering
- ✅ **Publication Ready** - Conservative, reliable results
- ✅ **Clean Datasets** - Minimal noise contamination

**Trade-offs:**
- ❌ **Lower Sensitivity** - May miss weak but real peaks
- ❌ **Reduced Data** - Fewer ML training samples
- ❌ **Strict Requirements** - May exclude borderline cases

#### **5. `aggressive` - Maximum Data Collection**
```bash
python batch_nmr_process.py /path/to/spectra --preset aggressive
```
**Best for:** ML training data collection, exploratory analysis
- **S/N Thresholds:** Lower across all nuclei (15N: 1.8, 13C: 1.5, 1H: 1.5)
- **Quality Threshold:** 0.6 (accept moderate quality)
- **Strategy:** Aggressive optimization approach
- **Result:** Maximum peak detection and data collection
- **Use Case:** ML training datasets, initial screening, comprehensive analysis

**Advantages:**
- ✅ **Maximum Sensitivity** - Detects weak peaks
- ✅ **Rich ML Datasets** - More training samples collected
- ✅ **Comprehensive Coverage** - Captures edge cases
- ✅ **Exploratory Power** - Find unexpected peaks

**Trade-offs:**
- ❌ **More Noise** - Increased false positive rate
- ❌ **Quality Variation** - Mixed quality training data
- ❌ **Processing Time** - More peaks to fit

#### **6. `high_throughput` - Speed Optimized**
```bash
python batch_nmr_process.py /path/to/spectra --preset high_throughput
```
**Best for:** Large datasets, screening applications, time-critical analysis
- **Fitting Attempts:** 1 (single attempt per peak)
- **Auto-Optimization:** Disabled (use fixed parameters)
- **Error Handling:** Skip errors immediately
- **Logging:** Minimal progress reporting
- **Result:** Fastest processing with acceptable quality

**Advantages:**
- ✅ **Maximum Speed** - Optimized for throughput
- ✅ **Scalable** - Handle hundreds of spectra
- ✅ **Resource Efficient** - Minimal computational overhead
- ✅ **Robust** - Continues on errors

**Trade-offs:**
- ❌ **Lower Quality** - No optimization refinement
- ❌ **Less Detailed** - Minimal diagnostic information
- ❌ **Fixed Parameters** - No adaptive optimization

---

### **📊 Preset Comparison Table**

| Preset | S/N Thresh. | Quality Req. | Speed | Data Quality | ML Samples | Best Use Case |
|--------|-------------|--------------|-------|--------------|------------|---------------|
| **`15N1H`** | 2.5 | 0.8 | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | Protein HSQC |
| **`13C1H`** | 2.0 | 0.8 | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | Carbon HSQC |
| **`1H`** | 2.0 | 0.8 | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | 1D Proton |
| **`conservative`** | 2.5-3.0 | 0.9 | ⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐ | Publication |
| **`aggressive`** | 1.5-1.8 | 0.6 | ⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ML Training |
| **`high_throughput`** | Default | 0.8 | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐ | Large Scale |

---

### **🚀 Preset Usage Examples**

#### **For Different Research Scenarios:**

```bash
# Protein Structure Study (high quality required)
python batch_nmr_process.py /data/protein_hsqc/ --preset conservative

# Drug Discovery Screening (speed critical)
python batch_nmr_process.py /data/compound_library/ --preset high_throughput

# ML Model Training (need lots of data)
python batch_nmr_process.py /data/training_spectra/ --preset aggressive

# Metabolomics Analysis (carbon focus)
python batch_nmr_process.py /data/metabolites/ --preset 13C1H --optimize

# Mixed Dataset (let auto-detection work)
python batch_nmr_process.py /data/mixed_spectra/ --preset balanced
```

#### **Combining Presets with Other Options:**

```bash
# Conservative preset + auto-optimization
python batch_nmr_process.py /path/to/spectra --preset conservative --optimize

# Aggressive preset + recursive search
python batch_nmr_process.py /path/to/spectra --preset aggressive --recursive

# High-throughput + specific nucleus
python batch_nmr_process.py /path/to/spectra --preset high_throughput --nucleus 15N1H

# Custom extensions + conservative preset
python batch_nmr_process.py /path/to/spectra --preset conservative --extensions .ft2 .pipe
```

---

### **📋 Advanced Preset Usage**

#### **Creating Custom Configurations:**

```bash
# Generate example config files to modify
python batch_nmr_process.py --create-examples /path/to/config_dir

# List all available presets
python batch_nmr_process.py --list-presets

# Use custom config based on preset
python batch_nmr_process.py /path/to/spectra --config my_custom_config.json
```

#### **Preset Override Examples:**

```json
{
  "// Custom config based on conservative preset": "",
  "processing": {
    "sn_thresholds": {
      "15N1H": 2.8,
      "default": 2.5
    },
    "quality_threshold": 0.85
  },
  "optimization": {
    "enable_auto_optimization": true,
    "optimization_strategy": "conservative"
  }
}
```

---

### **🎯 Recommendation Guide**

#### **Choose Your Preset Based On:**

**📊 Data Quality Priority:**
- **High Quality Needed:** `conservative` preset
- **Balanced Quality/Quantity:** Nucleus-specific presets (`15N1H`, `13C1H`, `1H`)
- **Maximum Data Collection:** `aggressive` preset

**⚡ Processing Speed Priority:**
- **Speed Critical:** `high_throughput` preset
- **Balanced Speed/Quality:** Nucleus-specific presets
- **Quality over Speed:** `conservative` preset

**🧬 Application Type:**
- **Protein Studies:** `15N1H` or `conservative`
- **Metabolomics:** `13C1H` or `aggressive`
- **Drug Screening:** `high_throughput` or `aggressive`
- **ML Training:** `aggressive` preset
- **Publication:** `conservative` preset

**💾 Dataset Size:**
- **Small Dataset (<50 spectra):** `conservative` or nucleus-specific
- **Medium Dataset (50-200 spectra):** Nucleus-specific presets
- **Large Dataset (>200 spectra):** `high_throughput` or `aggressive`

### **✅ Preset Selection Made Simple**

1. **Start with nucleus-specific presets** (`15N1H`, `13C1H`, `1H`) for most cases
2. **Use `conservative`** when quality is paramount
3. **Use `aggressive`** when collecting ML training data
4. **Use `high_throughput`** for large-scale screening
5. **Combine with `--optimize`** for best results (except high_throughput)

**Most Common Commands:**
```bash
# Recommended for most users (auto-detects nucleus type)
python batch_nmr_process.py /path/to/spectra --optimize

# For protein spectra
python batch_nmr_process.py /path/to/spectra --preset 15N1H --optimize

# For ML training data collection
python batch_nmr_process.py /path/to/spectra --preset aggressive --optimize
```