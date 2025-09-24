# Enhanced ML Training Data Collection v2.0 📊

This directory contains ML training data collected by the **Enhanced ML Training Data Collection System v2.0**, designed for advanced NMR peak fitting and analysis model development.

## 🎯 **Overview**

The Enhanced ML Data Collection v2.0 system transparently captures comprehensive training data from NMR peak fitting operations, providing rich datasets for developing sophisticated machine learning models for:

- **Peak Detection & Classification**
- **Automated Peak Fitting**
- **Quality Assessment & Validation**
- **Hard Case Learning & Robustness**
- **Cross-Method Consensus Analysis**

## 📁 **Data Structure**

### **Batch Files**
- `enhanced_training_batch_YYYYMMDD_HHMMSS.pkl` - Training data batches (Python pickle format)
- `enhanced_training_batch_YYYYMMDD_HHMMSS.json` - Metadata and collection statistics

### **Sample Structure**
Each training sample contains **8 comprehensive feature categories** with **92+ total features**:

```json
{
  "timestamp": "2025-09-16T10:30:45.123456",
  "collection_version": "2.0",
  "sample_id": "uuid-string",

  "spectral_features": {
    "peak_height": 120000.0,
    "baseline_level": 30000.0,
    "snr": 45.2,
    "peak_width_fwhm": 0.008,
    "peak_width_moment": 0.0095,
    "peak_asymmetry": 0.12,
    "noise_variance": 1250.0,
    "spectral_resolution": 0.004,
    // ... 27 total spectral features
  },

  "chemical_features": {
    "nucleus_type": "1H",
    "chemical_shift_center": 7.95,
    "chemical_region": "aromatic",
    "coupling_environment": "doublet",
    "relaxation_estimate": 1.5,
    // ... 11 total chemical features
  },

  "optimization_features": {
    "initial_parameters": {...},
    "convergence_iterations": 5,
    "optimization_difficulty": 0.15,
    "parameter_correlation_matrix": {...},
    // ... 15 total optimization features
  },

  "quality_metrics": {
    "r_squared": 0.985,
    "label": "good",  // "good", "reject", or filtered
    "rmse_normalized": 0.023,
    "aic": 145.2,
    "bic": 162.8
  }
}
```

## 🏷️ **Quality Labeling System**

### **Three-Tier Quality Classification:**

| Label | R² Range | Purpose | Collection |
|-------|----------|---------|------------|
| **`good`** | ≥ 0.8 | High-quality training samples | ✅ Always collected |
| **`reject`** | 0.3 ≤ R² < 0.8 | Hard cases for adversarial learning | ⚙️ Configurable (`collect_rejected_samples`) |
| **Filtered** | < 0.3 | Poor quality, unusable | ❌ Never collected |

## 📊 **Feature Categories Breakdown**

### **1. Spectral Features (27 features)**
- **Peak Characteristics**: Height, width measurements (FWHM, moment, percentiles)
- **Baseline Analysis**: Level, drift, slope, curvature
- **Noise Assessment**: Variance, correlation, robust estimates
- **Signal Quality**: SNR, dynamic range, intensity distribution

### **2. Chemical Features (11 features)**
- **Nuclear Properties**: Nucleus type, chemical shifts, regions
- **Coupling Information**: J-coupling estimates, multiplicity predictions
- **Relaxation**: T1/T2 estimates, chemical shift anisotropy

### **3. Context Features (11 features)**
- **Detection Confidence**: Algorithm confidence scores
- **Peak Environment**: Nearby peaks, overlap severity, isolation
- **Data Quality**: Baseline stability, interference assessment

### **4. Optimization Features (15 features)**
- **Convergence Analysis**: Iteration counts, difficulty metrics
- **Parameter Evolution**: Initial vs final parameters, correlation matrices
- **Algorithm Diagnostics**: Gradient norms, Hessian conditions

### **5. Multi-Method Results (8 features)**
- **Method Comparison**: Alternative fitting results, consensus analysis
- **Stability Assessment**: Cross-method parameter agreement
- **Robustness Metrics**: R² ranges, parameter stability scores

### **6. Physics Validation (5 features)**
- **Parameter Validity**: Physical realism checks
- **Spectroscopic Constraints**: Linewidth consistency
- **Shape Analysis**: Peak shape realism scores

## 🎨 **Enhanced v2.0 Capabilities**

### **New in Version 2.0:**
- ✅ **Optimization Trajectory Tracking** - Full parameter evolution history
- ✅ **Multi-Method Consensus** - Results from alternative fitting approaches
- ✅ **Physics-Based Validation** - Spectroscopically realistic parameter constraints
- ✅ **Hard Case Collection** - Rejected samples for adversarial training
- ✅ **Advanced Quality Metrics** - AIC/BIC model selection, residual analysis
- ✅ **Uncertainty Quantification** - Parameter confidence intervals
- ✅ **Cross-Validation Support** - Hold-out validation metrics

## ⚙️ **Configuration Options**

```python
from lunaNMR.ml.training_data_collector import MLTrainingDataCollector

collector = MLTrainingDataCollector()

# Quality thresholds
collector.min_r_squared = 0.8          # Good vs reject threshold
collector.min_r_squared_absolute = 0.3  # Absolute minimum

# Collection behavior
collector.collect_rejected_samples = True  # Collect hard cases
collector.batch_size = 100                # Auto-save frequency
collector.max_samples_per_session = 1000  # Memory protection
```

## 🔬 **Integration Points**

The Enhanced ML Collection system integrates seamlessly with:

### **Core Components:**
- ✅ **EnhancedVoigtFitter** - Peak fitting with optimization tracking
- ✅ **CoreIntegrator** - S/N detection with quality assessment
- ✅ **Level 2 Parameter Estimation** - Robust multi-method consensus

### **Backward Compatibility:**
- ✅ **v1.0 Interface Support** - Legacy method compatibility
- ✅ **Existing Workflows** - No disruption to current operations
- ✅ **Progressive Enhancement** - Opt-in advanced features

## 📈 **Data Statistics Example**

```json
{
  "collection_stats": {
    "total_attempts": 1523,
    "successfully_collected": 1445,
    "quality_filtered": 78,
    "hard_cases_collected": 234,
    "optimization_trajectories_saved": 1445,
    "multi_method_comparisons": 892,
    "enhanced_features_extracted": 1445,
    "physics_violations": 12,
    "batch_saves": 15
  },
  "success_rate": "94.9%",
  "quality_distribution": {
    "good": 1211,
    "reject": 234
  }
}
```

## 🚀 **Usage for ML Model Development**

### **Loading Training Data:**
```python
import pickle
import pandas as pd

# Load batch data
with open('enhanced_training_batch_20250916_143022.pkl', 'rb') as f:
    training_samples = pickle.load(f)

# Convert to DataFrame for analysis
features_data = []
labels = []

for sample in training_samples:
    # Extract all features
    row = {}
    row.update(sample['spectral_features'])
    row.update(sample['chemical_features'])
    row.update(sample['context_features'])

    features_data.append(row)
    labels.append(sample['quality_metrics']['label'])

df = pd.DataFrame(features_data)
quality_labels = pd.Series(labels)
```

### **Model Development Applications:**
1. **Binary Classification**: Good vs reject quality prediction
2. **Parameter Regression**: Direct parameter prediction from spectra
3. **Anomaly Detection**: Unusual peak behavior identification
4. **Method Selection**: Optimal fitting algorithm recommendation
5. **Uncertainty Estimation**: Confidence interval prediction

## 🛡️ **Data Quality Assurance**

- **Comprehensive Validation**: Physics-based parameter checking
- **Outlier Detection**: Statistical anomaly identification
- **Feature Schema Alignment**: Metadata matches actual data structure
- **Version Tracking**: Full provenance and reproducibility
- **Error Handling**: Graceful degradation, no workflow disruption

## 📋 **Maintenance & Monitoring**

### **Health Monitoring:**
- Collection success rates tracked automatically
- Error statistics logged for debugging
- Feature extraction completeness verified
- Integration point status monitored

### **Data Management:**
- Automatic batch saving prevents memory issues
- Configurable storage paths and retention policies
- Metadata preservation for batch traceability
- Compression for efficient storage

---

## 🎓 **For Researchers & Developers**

This Enhanced ML Training Data Collection v2.0 system provides unprecedented insight into NMR peak fitting processes, enabling development of:

- **Next-generation automated NMR analysis tools**
- **Robust peak detection algorithms**
- **Quality-aware fitting methods**
- **Physics-informed machine learning models**
- **Uncertainty-quantified spectral analysis**

The rich, multi-dimensional dataset supports both **supervised learning** (using quality labels) and **unsupervised discovery** of spectral patterns and relationships.

## 🔗 **ML Parameter Estimator Compatibility**

The Enhanced ML Data Collection v2.0 system is **100% compatible** with ML Parameter Estimator development:

### **✅ Complete Interface Support:**
- **Input Data**: Raw spectral data (x_data, y_data) preserved
- **Target Parameters**: All 5 Voigt parameters + uncertainties
- **Context Features**: Detection confidence, peak environment (11 features)
- **Quality Labels**: Good/reject classification for robust training

### **✅ Advanced Training Capabilities:**
- **92+ Features**: Comprehensive spectral, chemical, optimization analysis
- **Physics-Informed**: Spectroscopic constraints and validation
- **Multi-Method Consensus**: Alternative fitting approaches captured
- **Uncertainty Quantification**: Parameter confidence intervals included

### **✅ Production-Ready Pipeline:**
```python
# Direct compatibility with ML Parameter Estimator interface
def estimate_initial_parameters(x_data, y_data, peak_center, nucleus_type, context):
    # v2.0 data provides ALL required inputs and enhanced training targets
    return trained_model.predict(extract_v2_features(x_data, y_data, context))
```

**See [ML_PARAMETER_ESTIMATOR_COMPATIBILITY.md](ML_PARAMETER_ESTIMATOR_COMPATIBILITY.md) for complete implementation guide.**

---

**Happy Model Building!** 🚀🧬📊