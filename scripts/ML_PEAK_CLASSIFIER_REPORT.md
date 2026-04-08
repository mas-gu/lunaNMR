# ML Peak Classifier Implementation Report

**Date:** 2026-01-04
**Session Summary:** Implementation and evaluation of ML-based peak classifiers for NMR spectra

---

## Next Steps & Recommendations

### CNN Architecture Tuning (Priority)

The CNN achieved 80.93% F1 vs RF's 83.75%. To improve:

1. **Increase Patch Size** (currently 21×21)
   ```bash
   python scripts/train_cnn_classifier.py \
       --training-dir ml_training_data/adaptativev2 \
       --spectra-dirs /path/to/spectra \
       --output lunaNMR/ml/pretrained/cnn_patch31.pt \
       --patch-size 31
   ```

2. **Increase Dropout** (currently 0.3)
   - Edit `lunaNMR/ml/cnn_peak_classifier.py` line ~45: `dropout: float = 0.5`

3. **Train Longer**
   ```bash
   --epochs 100 --patience 20
   ```

4. **Multi-scale Architecture** - Feed 15×15, 21×21, 31×31 patches simultaneously

5. **Add Attention Mechanism** - Self-attention to focus on peak center

### Data Collection (High Impact)

- **13C spectra:** Only 8 samples in training data (vs 126 for 15N)
- Collecting more 13C data would enable nucleus-specific models
- Current models are trained on combined 15N+13C data

### Hybrid Approach Refinement

- RF hybrid outperforms CNN hybrid (83.75% vs 80.93%)
- Consider: CNN features + RF classifier (best of both worlds)

---

## Results Summary

### Hybrid Mode (S/N Detection → ML Filter)

| Model | F1 Score | Precision | Recall | Filter Rate | vs S/N Baseline |
|-------|----------|-----------|--------|-------------|-----------------|
| **RF (16 features)** | **83.75%** | 83.17% | 84.33% | 22.6% | **+4.0%** |
| LightGBM | 82.86% | 84.51% | 81.26% | 26.4% | +3.3% |
| XGBoost | 82.78% | 87.55% | 78.50% | 31.6% | +3.0% |
| CNN (21×21 patches) | 80.93% | 79.43% | 82.49% | 7.6% | +0.6% |
| S/N Baseline | 79.75% | 75.87% | 85.24% | - | - |

### Key Observations

1. **RF is the best performer** with balanced precision/recall
2. **XGBoost has highest precision** (87.55%) but filters too aggressively
3. **CNN underperforms** - only filters 7.6% vs RF's 22.6%
4. **All ML methods improve over S/N baseline** in F1 score

---

## Feature Engineering

### Current Features (16 total)

| # | Feature | Description | Importance (RF) |
|---|---------|-------------|-----------------|
| 1 | intensity | Peak intensity value | 17.05% |
| 2 | gaussian_fit_r2 | R² of 2D Gaussian fit | 12.20% |
| 3 | gradient_magnitude | Gradient at peak apex | 12.09% |
| 4 | prominence_f2 | Prominence in F2 (1H) | 11.43% |
| 5 | prominence_f1 | Prominence in F1 (15N/13C) | 11.01% |
| 6 | local_snr | Local signal-to-noise ratio | 10.14% |
| 7 | template_correlation | Correlation with Gaussian template | 9.04% |
| 8 | symmetry_f1 | Peak symmetry in F1 | 2.50% |
| 9 | fwhm_f1 | Full-width half-max in F1 | 2.08% |
| 10 | neighborhood_peaks | Density of nearby peaks | 1.98% |
| 11 | fwhm_f2 | Full-width half-max in F2 | 1.92% |
| 12 | curvature | 2D curvature at apex | 1.84% |
| 13 | local_max_count | Count of local maxima nearby | 1.77% |
| 14 | peak_sharpness | 2nd derivative magnitude | 1.73% |
| 15 | symmetry_f2 | Peak symmetry in F2 | 1.68% |
| 16 | aspect_ratio | FWHM_f1 / FWHM_f2 ratio | 1.53% |

### Features Added This Session

1. **template_correlation** - Normalized cross-correlation with 2D Gaussian template
2. **gaussian_fit_r2** - R² from fitting 2D Gaussian to local region
3. **peak_sharpness** - 5-point stencil 2nd derivative normalized by intensity
4. **neighborhood_peaks** - Count of significant peaks in 15×15 region

---

## CNN Implementation

### Architecture

```
Input: (batch, 1, 21, 21)
│
├─ Conv2D(1→16, 3×3) → BatchNorm → ReLU → MaxPool(2×2)
│  Output: (batch, 16, 10, 10)
│
├─ Conv2D(16→32, 3×3) → BatchNorm → ReLU → MaxPool(2×2)
│  Output: (batch, 32, 5, 5)
│
├─ Conv2D(32→64, 3×3) → BatchNorm → ReLU → AdaptiveAvgPool(1×1)
│  Output: (batch, 64, 1, 1)
│
├─ Flatten → Linear(64→32) → ReLU → Dropout(0.3)
│
└─ Linear(32→1) → Sigmoid
   Output: probability
```

**Parameters:** 25,633
**Inference time:** ~0.6ms per candidate (CPU)

### Training Results

```
Training samples: 35,829
Validation samples: 7,637
Test samples: 7,080
Best validation loss: 0.2577
Best validation accuracy: 90.95%
Epochs trained: 49
Training time: 793.8s
```

### Files Created

| File | Description |
|------|-------------|
| `lunaNMR/ml/cnn_peak_classifier.py` | CNN model definition and inference wrapper |
| `scripts/train_cnn_classifier.py` | Training script with augmentation |
| `scripts/evaluate_cnn_classifier.py` | Evaluation script (hybrid mode) |
| `lunaNMR/ml/pretrained/cnn_peak_classifier.pt` | Trained model |

---

## Training Data

### Dataset Statistics

- **Total spectra:** 134
- **Total peaks:** 22,201
- **Good quality peaks (R² ≥ 0.85):** 19,897 (89.6%)

### Nucleus Distribution

| Nucleus | Spectra | Notes |
|---------|---------|-------|
| 15N-HSQC | 126 | Sufficient for training |
| 13C-HSQC | 8 | **Needs more data** |

### Train/Val/Test Split

- **Spectrum-level split** (entire spectra held out, not individual peaks)
- Train: ~70% of spectra
- Validation: ~15%
- Test: ~15%

---

## Commands Reference

### Training

**RF/XGBoost/LightGBM (16 features):**
```bash
python scripts/train_peak_classifier.py \
    --training-dir ml_training_data/adaptativev2 \
    --spectra-dirs /path/to/spectra \
    --output lunaNMR/ml/pretrained/model.joblib \
    --model-type rf  # or xgb, lgb
```

**CNN:**
```bash
# Requires separate conda environment with Python 3.12
conda activate cnn_training
python scripts/train_cnn_classifier.py \
    --training-dir ml_training_data/adaptativev2 \
    --spectra-dirs /path/to/spectra \
    --output lunaNMR/ml/pretrained/cnn_peak_classifier.pt \
    --epochs 50 --batch-size 64
```

### Evaluation

**RF/XGBoost/LightGBM:**
```bash
python scripts/evaluate_peak_classifier.py \
    --model lunaNMR/ml/pretrained/model.joblib \
    --spectra-dirs /path/to/spectra \
    --hybrid --test-only --per-nucleus
```

**CNN:**
```bash
KMP_DUPLICATE_LIB_OK=TRUE python scripts/evaluate_cnn_classifier.py \
    --model lunaNMR/ml/pretrained/cnn_peak_classifier.pt \
    --spectra-dirs /path/to/spectra \
    --hybrid --test-only --hybrid-threshold 0.3
```

### Nucleus-Specific Models

```bash
# Train 15N-specific model
python scripts/train_peak_classifier.py \
    --training-dir ml_training_data/adaptativev2 \
    --spectra-dirs /path/to/spectra \
    --output lunaNMR/ml/pretrained/peak_classifier_15N_rf.joblib \
    --model-type rf --nucleus-type 15N

# Evaluate with nucleus-specific models
python scripts/evaluate_peak_classifier.py \
    --model-15N lunaNMR/ml/pretrained/peak_classifier_15N_rf.joblib \
    --model-13C lunaNMR/ml/pretrained/peak_classifier_13C_rf.joblib \
    --spectra-dirs /path/to/spectra \
    --hybrid --per-nucleus --test-only
```

---

## Environment Setup

### Main Environment (Python 3.13)

Works for RF/XGBoost/LightGBM evaluation.

### CNN Training Environment

PyTorch doesn't support Python 3.13 yet. Create separate environment:

```bash
conda create -n cnn_training python=3.12 -y
conda activate cnn_training
pip install torch "numpy<2.0" scipy pandas nmrglue scikit-learn
```

**Known Issues:**
- OpenMP conflict: Set `KMP_DUPLICATE_LIB_OK=TRUE`
- Numpy 2.x incompatibility: Use `numpy<2.0`

---

## Issues Encountered & Solutions

### 1. Feature Count Mismatch

**Problem:** Model trained with 14 features, extractor produces 16 features
**Solution:** Retrain models after adding new features

### 2. 13C Spectra Not Evaluated

**Problem:** 13C test spectra not found during evaluation
**Solution:** Ensure `--spectra-dirs` includes path to 13C data

### 3. PyTorch + Numpy Compatibility

**Problem:** `RuntimeError: Numpy is not available`
**Solution:** Use `torch.tensor()` instead of `torch.from_numpy()`, and `probs.cpu().tolist()` instead of `.numpy()`

### 4. OpenMP Conflict on macOS

**Problem:** `libiomp5.dylib already initialized`
**Solution:** `export KMP_DUPLICATE_LIB_OK=TRUE`

---

## Future Work

1. **Collect more 13C training data** - Only 8 spectra currently
2. **Hybrid CNN+RF** - Use CNN embeddings as features for RF
3. **Architecture search** - Try different patch sizes, depths
4. **Transfer learning** - Pre-train on synthetic peaks
5. **Active learning** - Identify uncertain predictions for manual labeling
6. **Integration into lunaNMR GUI** - Add ML detection option

---

## File Locations

```
lunaNMR_v1o0/
├── lunaNMR/ml/
│   ├── cnn_peak_classifier.py      # CNN model
│   ├── peak_feature_extractor.py   # 16-feature extractor
│   ├── peak_classifier.py          # RF/XGB/LGB wrapper
│   └── pretrained/
│       ├── cnn_peak_classifier.pt  # Trained CNN
│       ├── rf_modelv2.joblib       # RF with 16 features
│       ├── xgb_modelv2.joblib      # XGBoost with 16 features
│       └── lgb_modelv2.joblib      # LightGBM with 16 features
├── scripts/
│   ├── train_peak_classifier.py    # Train RF/XGB/LGB
│   ├── train_cnn_classifier.py     # Train CNN
│   ├── evaluate_peak_classifier.py # Evaluate RF/XGB/LGB
│   └── evaluate_cnn_classifier.py  # Evaluate CNN
└── ml_training_data/
    └── adaptativev2/
        └── training_data.json      # 134 spectra, 22K peaks
```
