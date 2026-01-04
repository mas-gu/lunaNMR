#!/usr/bin/env python3
# ABOUTME: Pre-training script for lunaNMR ML parameter prediction models
# ABOUTME: Trains RandomForest models from JSON data collected by ComprehensiveTrainingCollector

"""
Pre-training script for lunaNMR ML models.

Reads training samples from JSON files (collected by ComprehensiveTrainingCollector),
extracts features and targets, trains separate models for 15N and 13C spectra,
and exports .joblib model files.

Usage:
    # Basic usage (uses default paths)
    python scripts/pretrain_models.py

    # Single file
    python scripts/pretrain_models.py --data-file my_data.json

    # Multiple files (merged for training)
    python scripts/pretrain_models.py --data-file lab1_data.json lab2_data.json idp_data.json

    # Glob pattern
    python scripts/pretrain_models.py --data-file "collected_data/*.json"

    # Full example with options
    python scripts/pretrain_models.py \\
        --data-file data/protein_*.json \\
        --output-dir lunaNMR/ml/pretrained \\
        --min-r2 0.80 \\
        --min-samples 20

    # Dry run (show stats without training)
    python scripts/pretrain_models.py --dry-run
"""

import argparse
import json
import logging
from datetime import datetime
from glob import glob
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union

import numpy as np

# Check for scikit-learn
try:
    import joblib
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import cross_val_score, train_test_split
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

# Check for tqdm (optional, for progress bars)
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    # Fallback: no-op wrapper
    def tqdm(iterable, **kwargs):
        return iterable

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


# =============================================================================
# CONSTANTS
# =============================================================================

FEATURE_NAMES = [
    # Original 11 features
    'nucleus_type',        # Binary: 15N=1, 13C=0
    'snr_estimate',
    'noise_level',
    'dynamic_range',
    'peak_count',
    'peak_density',
    'shift_range_f1_min',
    'shift_range_f1_max',
    'shift_range_f2_min',
    'shift_range_f2_max',
    'field_strength_mhz',
    # Pre-fitting features from spectrum/peak list (9 fields)
    'baseline_level',
    'baseline_std',
    'mean_peak_separation_f1',
    'mean_peak_separation_f2',
    'min_peak_separation_f1',
    'min_peak_separation_f2',
    'f1_dispersion',
    'f2_dispersion',
    'is_idp_like',         # Binary: True=1, False=0
    # Detected intensity stats (from peak detection, before fitting) - 6 fields
    'detected_intensity_mean',
    'detected_intensity_std',
    'detected_intensity_cv',
    'detected_intensity_dynamic_range',
    'detected_intensity_skewness',
    'detected_intensity_kurtosis',
    # Additional pre-fitting features (11 fields)
    'dispersion_ratio',    # f1_dispersion / f2_dispersion
    'n_close_pairs_f1',    # Peak pairs within 2*typical_lw in F1
    'n_close_pairs_f2',    # Peak pairs within 2*typical_lw in F2
    'n_close_pairs_2d',    # Peak pairs within elliptical distance
    'fraction_potentially_overlapping',
    'peaklist_intensity_mean',
    'peaklist_intensity_std',
    'peaklist_intensity_cv',
    'peaklist_intensity_dynamic_range',
    'max_local_density',
    'crowding_hotspot_fraction',
]

TARGET_NAMES = [
    'lw_f1_median',
    'lw_f2_median',
    'rad_f1',
    'rad_f2',
    'overlap_threshold_f1',
    'overlap_threshold_f2',
    'achievable_r2',
]

DEFAULT_MODEL_PARAMS = {
    'n_estimators': 100,
    'max_depth': 10,
    'min_samples_leaf': 5,
    'n_jobs': -1,
    'random_state': 42,
}

# =============================================================================
# PER-PEAK TRAINING CONSTANTS
# =============================================================================

# Peak-specific features (from PeakTrainingData)
PEAK_FEATURE_NAMES = [
    # Context flags
    'is_isolated',           # 1 if isolated, 0 if clustered
    'cluster_size',          # Size of cluster (1 for isolated)
    'tooclose_flag',         # Has peaks within tooclose threshold
    'heavy_overlap_flag',    # Has heavy overlap
    # Position context (normalized to spectrum range)
    'relative_pos_f1',       # (pos_f1 - f1_min) / f1_dispersion
    'relative_pos_f2',       # (pos_f2 - f2_min) / f2_dispersion
    # Intensity context (normalized to spectrum)
    'relative_intensity',    # detected_intensity / detected_intensity_mean
    'log_intensity',         # log10(detected_intensity) for scale invariance
    # Fitting constraints applied
    'fix_positions_applied',
    'fix_linewidths_applied',
    'fitting_mode_2d',       # 1 if 2D, 0 if 1D
    # Spectrum context (from parent sample)
    'nucleus_type',          # 15N=1, 13C=0
    'snr_estimate',
    'noise_level_normalized', # noise_level / detected_intensity_mean
    'peak_density',
    'fraction_clustered',
    'is_idp_like',
    'field_strength_mhz',
    'pass1_lw_f1_median',    # Spectrum-level median linewidth
    'pass1_lw_f2_median',
    'fraction_potentially_overlapping',
    'f1_dispersion',
    'f2_dispersion',
]

# Per-peak linewidth targets
PEAK_LW_TARGET_NAMES = [
    'lw_total_f1',           # Total F1 linewidth (ppm)
    'lw_total_f2',           # Total F2 linewidth (ppm)
]

# Per-peak quality targets
PEAK_QUALITY_TARGET_NAMES = [
    'r_squared',             # Fit quality (0-1)
]


# =============================================================================
# DATA LOADING
# =============================================================================

def load_training_data(data_source: Union[Path, str, List[Path]]) -> List[Dict]:
    """
    Load training samples from one or more JSON files.

    Parameters
    ----------
    data_source : Path, str, or list of Paths
        - Single Path: load from that file
        - String with '*': glob pattern to match multiple files
        - List of Paths: load and merge all files

    Returns
    -------
    list
        List of sample dictionaries (merged if multiple files)
    """
    # Determine list of files to load
    if isinstance(data_source, list):
        # List of paths provided
        files = [Path(f) for f in data_source]
    elif isinstance(data_source, str) and '*' in data_source:
        # Glob pattern
        matched = glob(data_source)
        if not matched:
            raise FileNotFoundError(f"No files match pattern: {data_source}")
        files = [Path(f) for f in sorted(matched)]
        logger.info(f"Glob pattern matched {len(files)} files")
    else:
        # Single file
        files = [Path(data_source)]

    # Validate all files exist
    for f in files:
        if not f.exists():
            raise FileNotFoundError(f"Training data file not found: {f}")

    # Load and merge samples from all files
    all_samples = []
    for data_file in files:
        logger.info(f"Loading training data from {data_file}...")
        with open(data_file) as f:
            data = json.load(f)
        samples = data.get('samples', [])
        version = data.get('version', 'unknown')
        logger.info(f"  Loaded {len(samples)} samples (version {version})")
        all_samples.extend(samples)

    if len(files) > 1:
        logger.info(f"Total: {len(all_samples)} samples from {len(files)} files")
    else:
        logger.info(f"Loaded {len(all_samples)} training samples")

    return all_samples


# =============================================================================
# QUALITY FILTERING
# =============================================================================

def filter_quality_samples(
    samples: List[Dict],
    min_r2: float = 0.80,
    min_peaks: int = 20
) -> List[Dict]:
    """
    Filter samples by quality criteria.

    Parameters
    ----------
    samples : list
        List of sample dictionaries
    min_r2 : float
        Minimum R² threshold (default 0.80)
    min_peaks : int
        Minimum peak count (default 20)

    Returns
    -------
    list
        Filtered samples meeting quality criteria
    """
    logger.info(f"Filtering by quality (R² >= {min_r2}, adaptive_success=True, peaks >= {min_peaks})...")

    good_samples = []
    for sample in samples:
        # Check R² threshold
        if sample.get('overall_r2_pass2', 0) < min_r2:
            continue

        # Check adaptive optimization succeeded
        if not sample.get('adaptive_success', False):
            continue

        # Check minimum peak count
        if sample.get('peak_count', 0) < min_peaks:
            continue

        # Check required fields exist and non-zero
        if sample.get('pass1_lw_f1_median', 0) <= 0:
            continue
        if sample.get('pass1_lw_f2_median', 0) <= 0:
            continue

        good_samples.append(sample)

    logger.info(f"After filtering: {len(good_samples)} quality samples")
    return good_samples


# =============================================================================
# NUCLEUS SPLIT
# =============================================================================

def split_by_nucleus(samples: List[Dict]) -> Tuple[List[Dict], List[Dict]]:
    """
    Split samples by nucleus type.

    Parameters
    ----------
    samples : list
        List of sample dictionaries

    Returns
    -------
    tuple
        (samples_15N, samples_13C)
    """
    samples_15N = [s for s in samples if s.get('nucleus_type') == '15N']
    samples_13C = [s for s in samples if s.get('nucleus_type') == '13C']

    logger.info(f"Split by nucleus: 15N={len(samples_15N)}, 13C={len(samples_13C)}")
    return samples_15N, samples_13C


# =============================================================================
# FEATURE EXTRACTION
# =============================================================================

def extract_features(sample: Dict) -> List[float]:
    """
    Extract 37-element feature vector from a single sample.

    Parameters
    ----------
    sample : dict
        Sample dictionary from training data

    Returns
    -------
    list
        37-element feature vector (all pre-fitting features)
    """
    return [
        # Original 11 features
        1.0 if sample.get('nucleus_type') == '15N' else 0.0,
        sample.get('snr_estimate', 0),
        sample.get('noise_level', 0),
        sample.get('dynamic_range', 0),
        sample.get('peak_count', 0),
        sample.get('peak_density', 0),
        sample.get('shift_range_f1_min', 0),
        sample.get('shift_range_f1_max', 0),
        sample.get('shift_range_f2_min', 0),
        sample.get('shift_range_f2_max', 0),
        sample.get('field_strength_mhz', 600),
        # Pre-fitting features from spectrum/peak list (9 fields)
        sample.get('baseline_level', 0),
        sample.get('baseline_std', 0),
        sample.get('mean_peak_separation_f1', 0),
        sample.get('mean_peak_separation_f2', 0),
        sample.get('min_peak_separation_f1', 0),
        sample.get('min_peak_separation_f2', 0),
        sample.get('f1_dispersion', 0),
        sample.get('f2_dispersion', 0),
        1.0 if sample.get('is_idp_like', False) else 0.0,
        # Detected intensity stats (from peak detection, before fitting) - 6 fields
        sample.get('detected_intensity_mean', 0),
        sample.get('detected_intensity_std', 0),
        sample.get('detected_intensity_cv', 0),
        sample.get('detected_intensity_dynamic_range', 0),
        sample.get('detected_intensity_skewness', 0),
        sample.get('detected_intensity_kurtosis', 0),
        # Additional pre-fitting features (11 fields)
        sample.get('dispersion_ratio', 0),
        sample.get('n_close_pairs_f1', 0),
        sample.get('n_close_pairs_f2', 0),
        sample.get('n_close_pairs_2d', 0),
        sample.get('fraction_potentially_overlapping', 0),
        sample.get('peaklist_intensity_mean', 0),
        sample.get('peaklist_intensity_std', 0),
        sample.get('peaklist_intensity_cv', 0),
        sample.get('peaklist_intensity_dynamic_range', 0),
        sample.get('max_local_density', 0),
        sample.get('crowding_hotspot_fraction', 0),
    ]


def extract_targets(sample: Dict) -> List[float]:
    """
    Extract 7-element target vector from a single sample.

    Parameters
    ----------
    sample : dict
        Sample dictionary from training data

    Returns
    -------
    list
        7-element target vector
    """
    return [
        sample.get('pass1_lw_f1_median', 0),
        sample.get('pass1_lw_f2_median', 0),
        sample.get('adaptive_radF1', 0),
        sample.get('adaptive_radF2', 0),
        sample.get('adaptive_overlap_threshold_y', 0),  # F1 threshold
        sample.get('adaptive_overlap_threshold_x', 0),  # F2 threshold
        sample.get('overall_r2_pass2', 0),
    ]


# =============================================================================
# PER-PEAK FEATURE/TARGET EXTRACTION
# =============================================================================

def extract_spectrum_context(sample: Dict) -> Dict[str, float]:
    """
    Extract spectrum-level context for per-peak features.

    Parameters
    ----------
    sample : dict
        Spectrum-level sample dictionary

    Returns
    -------
    dict
        Context values needed for per-peak feature extraction
    """
    nucleus = sample.get('nucleus_type', '15N')
    nucleus_val = 1.0 if nucleus == '15N' else 0.0

    # Get detected intensity mean for normalization
    det_int_mean = sample.get('detected_intensity_mean', 1.0)
    if det_int_mean == 0:
        det_int_mean = 1.0  # Avoid division by zero

    return {
        'nucleus_type': nucleus_val,
        'snr_estimate': sample.get('snr_estimate', 0),
        'noise_level': sample.get('noise_level', 0),
        'noise_level_normalized': sample.get('noise_level', 0) / det_int_mean,
        'peak_density': sample.get('peak_density', 0),
        'fraction_clustered': sample.get('fraction_clustered', 0),
        'is_idp_like': 1.0 if sample.get('is_idp_like', False) else 0.0,
        'field_strength_mhz': sample.get('field_strength_mhz', 600.0),
        'pass1_lw_f1_median': sample.get('pass1_lw_f1_median', 0),
        'pass1_lw_f2_median': sample.get('pass1_lw_f2_median', 0),
        'fraction_potentially_overlapping': sample.get('fraction_potentially_overlapping', 0),
        'f1_dispersion': sample.get('f1_dispersion', 1.0),
        'f2_dispersion': sample.get('f2_dispersion', 1.0),
        'shift_range_f1_min': sample.get('shift_range_f1_min', 0),
        'shift_range_f2_min': sample.get('shift_range_f2_min', 0),
        'detected_intensity_mean': det_int_mean,
    }


def extract_peak_features(peak: Dict, context: Dict) -> List[float]:
    """
    Extract feature vector for a single peak.

    Parameters
    ----------
    peak : dict
        Peak data dictionary (from sample['peaks'])
    context : dict
        Spectrum context from extract_spectrum_context()

    Returns
    -------
    list
        Feature vector matching PEAK_FEATURE_NAMES
    """
    # Position normalization
    f1_disp = context.get('f1_dispersion', 1.0)
    f2_disp = context.get('f2_dispersion', 1.0)
    f1_min = context.get('shift_range_f1_min', 0)
    f2_min = context.get('shift_range_f2_min', 0)

    pos_f1 = peak.get('pos_f1', 0)
    pos_f2 = peak.get('pos_f2', 0)
    relative_pos_f1 = (pos_f1 - f1_min) / f1_disp if f1_disp > 0 else 0
    relative_pos_f2 = (pos_f2 - f2_min) / f2_disp if f2_disp > 0 else 0

    # Intensity normalization
    det_int = peak.get('detected_intensity', 0)
    det_int_mean = context.get('detected_intensity_mean', 1.0)
    relative_intensity = det_int / det_int_mean if det_int_mean > 0 else 0
    log_intensity = np.log10(max(det_int, 1.0))  # Avoid log(0)

    # Fitting mode
    fitting_mode = peak.get('fitting_mode', '1D')
    fitting_mode_2d = 1.0 if fitting_mode == '2D' else 0.0

    return [
        # Context flags
        1.0 if peak.get('is_isolated', False) else 0.0,
        float(peak.get('cluster_size', 1)),
        1.0 if peak.get('tooclose_flag', False) else 0.0,
        1.0 if peak.get('heavy_overlap_flag', False) else 0.0,
        # Position context
        relative_pos_f1,
        relative_pos_f2,
        # Intensity context
        relative_intensity,
        log_intensity,
        # Fitting constraints
        1.0 if peak.get('fix_positions_applied', False) else 0.0,
        1.0 if peak.get('fix_linewidths_applied', False) else 0.0,
        fitting_mode_2d,
        # Spectrum context
        context.get('nucleus_type', 1.0),
        context.get('snr_estimate', 0),
        context.get('noise_level_normalized', 0),
        context.get('peak_density', 0),
        context.get('fraction_clustered', 0),
        context.get('is_idp_like', 0),
        context.get('field_strength_mhz', 600.0),
        context.get('pass1_lw_f1_median', 0),
        context.get('pass1_lw_f2_median', 0),
        context.get('fraction_potentially_overlapping', 0),
        context.get('f1_dispersion', 1.0),
        context.get('f2_dispersion', 1.0),
    ]


def extract_peak_linewidth_targets(peak: Dict) -> List[float]:
    """
    Extract linewidth targets for a single peak.

    Parameters
    ----------
    peak : dict
        Peak data dictionary

    Returns
    -------
    list
        [lw_total_f1, lw_total_f2]
    """
    return [
        peak.get('lw_total_f1', 0),
        peak.get('lw_total_f2', 0),
    ]


def extract_peak_quality_targets(peak: Dict) -> List[float]:
    """
    Extract quality targets for a single peak.

    Parameters
    ----------
    peak : dict
        Peak data dictionary

    Returns
    -------
    list
        [r_squared]
    """
    return [
        peak.get('r_squared', 0),
    ]


def flatten_peak_data(
    samples: List[Dict],
    min_r2: float = 0.0,
    require_success: bool = True
) -> Tuple[List[Dict], List[Dict]]:
    """
    Flatten nested samples to individual peak records, split by nucleus type.

    Parameters
    ----------
    samples : list
        List of spectrum-level sample dictionaries
    min_r2 : float
        Minimum R-squared to include peak (default 0.0 = include all)
    require_success : bool
        If True, only include peaks with success=True

    Returns
    -------
    tuple
        (peaks_15N, peaks_13C) - lists of peak records with spectrum context
    """
    peaks_15N = []
    peaks_13C = []

    for sample in tqdm(samples, desc="Flattening peaks", unit="spectrum", disable=not TQDM_AVAILABLE):
        nucleus = sample.get('nucleus_type', '15N')
        context = extract_spectrum_context(sample)
        peaks = sample.get('peaks', [])

        for peak in peaks:
            # Filter by success
            if require_success and not peak.get('success', False):
                continue

            # Filter by R-squared
            r2 = peak.get('r_squared', 0)
            if r2 < min_r2:
                continue

            # Skip peaks with invalid linewidths
            lw_f1 = peak.get('lw_total_f1', 0)
            lw_f2 = peak.get('lw_total_f2', 0)
            if lw_f1 <= 0 or lw_f2 <= 0:
                continue

            # Create peak record with context
            peak_record = {
                **peak,
                '_context': context,
            }

            if nucleus == '15N':
                peaks_15N.append(peak_record)
            else:
                peaks_13C.append(peak_record)

    return peaks_15N, peaks_13C


def build_peak_matrices(
    peaks: List[Dict],
    target_type: str = 'linewidth'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build feature matrix X and target matrix y for per-peak training.

    Parameters
    ----------
    peaks : list
        List of peak records from flatten_peak_data()
    target_type : str
        'linewidth' or 'quality'

    Returns
    -------
    tuple
        (X, y) numpy arrays
    """
    if not peaks:
        return np.array([]), np.array([])

    X = []
    y = []

    desc = f"Building {target_type} matrix"
    for peak in tqdm(peaks, desc=desc, unit="peak", disable=not TQDM_AVAILABLE):
        context = peak.get('_context', {})
        features = extract_peak_features(peak, context)
        X.append(features)

        if target_type == 'linewidth':
            targets = extract_peak_linewidth_targets(peak)
        else:
            targets = extract_peak_quality_targets(peak)
        y.append(targets)

    X = np.array(X, dtype=np.float64)
    y = np.array(y, dtype=np.float64)

    # For single-target models, sklearn expects 1D array
    if y.shape[1] == 1:
        y = y.ravel()

    return X, y


# =============================================================================
# MATRIX BUILDING
# =============================================================================

def build_matrices(samples: List[Dict]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build feature matrix X and target matrix y from samples.

    Parameters
    ----------
    samples : list
        List of sample dictionaries

    Returns
    -------
    tuple
        (X, y) numpy arrays with shapes (n_samples, 37) and (n_samples, 7)
    """
    X = np.array([
        extract_features(s) for s in tqdm(
            samples, desc="Building spectrum matrix", unit="sample", disable=not TQDM_AVAILABLE
        )
    ], dtype=np.float64)
    y = np.array([extract_targets(s) for s in samples], dtype=np.float64)

    logger.info(f"  Feature matrix shape: {X.shape}")
    logger.info(f"  Target matrix shape: {y.shape}")

    return X, y


# =============================================================================
# MODEL TRAINING
# =============================================================================

def train_model(
    X: np.ndarray,
    y: np.ndarray,
    nucleus: str,
    cv_folds: int = 5
) -> Tuple[RandomForestRegressor, float, float]:
    """
    Train RandomForest model with cross-validation.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix (n_samples, 11)
    y : np.ndarray
        Target matrix (n_samples, 7)
    nucleus : str
        Nucleus type for logging ("15N" or "13C")
    cv_folds : int
        Number of cross-validation folds

    Returns
    -------
    tuple
        (model, cv_mean, cv_std)
    """
    if not SKLEARN_AVAILABLE:
        raise ImportError("scikit-learn is required for training. Install with: pip install scikit-learn")

    logger.info(f"Training {nucleus} model ({len(X)} samples)...")

    model = RandomForestRegressor(**DEFAULT_MODEL_PARAMS)

    # Cross-validation (use fewer folds if not enough samples)
    actual_cv = min(cv_folds, len(X))
    if actual_cv < 2:
        logger.warning(f"  Not enough samples for cross-validation, skipping CV")
        model.fit(X, y)
        return model, 0.0, 0.0

    cv_scores = cross_val_score(model, X, y, cv=actual_cv, scoring='r2')
    cv_mean = float(cv_scores.mean())
    cv_std = float(cv_scores.std())

    logger.info(f"  Cross-validation R²: {cv_mean:.3f} ± {cv_std:.3f}")

    # Train final model on all data
    model.fit(X, y)

    return model, cv_mean, cv_std


def evaluate_on_test_set(
    model: RandomForestRegressor,
    X_test: np.ndarray,
    y_test: np.ndarray,
    model_name: str
) -> float:
    """
    Evaluate a trained model on held-out test data.

    Parameters
    ----------
    model : RandomForestRegressor
        Trained model
    X_test : np.ndarray
        Test features
    y_test : np.ndarray
        Test targets
    model_name : str
        Name for logging

    Returns
    -------
    float
        Test R² score
    """
    from sklearn.metrics import r2_score

    y_pred = model.predict(X_test)
    test_r2 = r2_score(y_test, y_pred)
    logger.info(f"  Test set R²: {test_r2:.3f} ({len(X_test)} samples)")
    return test_r2


def display_feature_importance(
    model: RandomForestRegressor,
    feature_names: List[str],
    model_name: str,
    top_n: int = 10
) -> List[Tuple[str, float]]:
    """
    Display feature importance rankings for a trained model.

    Parameters
    ----------
    model : RandomForestRegressor
        Trained model
    feature_names : list
        Names of features
    model_name : str
        Name for display
    top_n : int
        Number of top features to show

    Returns
    -------
    list
        List of (feature_name, importance) tuples sorted by importance
    """
    importances = model.feature_importances_
    indices = np.argsort(importances)[::-1]

    logger.info(f"\n  Feature Importance ({model_name}):")
    logger.info(f"  {'─' * 45}")

    ranked = []
    for i in range(min(top_n, len(feature_names))):
        idx = indices[i]
        name = feature_names[idx] if idx < len(feature_names) else f"feature_{idx}"
        importance = importances[idx]
        ranked.append((name, importance))

        # Visual bar
        bar_len = int(importance * 40)
        bar = '█' * bar_len
        logger.info(f"  {i+1:2d}. {name:30s} {importance:.3f} {bar}")

    return ranked


def train_peak_models(
    samples: List[Dict],
    output_dir: Path,
    min_r2_linewidth: float = 0.85,
    min_r2_quality: float = 0.0,
    cv_folds: int = 5,
    min_samples: int = 100,
    test_split: float = 0.0,
    show_feature_importance: bool = False,
) -> Dict[str, any]:
    """
    Train per-peak models for linewidth and quality prediction.

    Parameters
    ----------
    samples : list
        List of spectrum-level sample dictionaries (with 'peaks' field)
    output_dir : Path
        Directory to save models
    min_r2_linewidth : float
        Minimum R² for peaks used in linewidth training
    min_r2_quality : float
        Minimum R² for peaks used in quality training (0 = include all)
    cv_folds : int
        Number of cross-validation folds
    min_samples : int
        Minimum peaks required to train a model
    test_split : float
        Fraction of data to hold out for testing (0 = no split)
    show_feature_importance : bool
        If True, display feature importance rankings

    Returns
    -------
    dict
        Training results and statistics
    """
    if not SKLEARN_AVAILABLE:
        raise ImportError("scikit-learn is required. Install with: pip install scikit-learn")

    results = {
        'linewidth_models': {},
        'quality_models': {},
        'statistics': {},
    }

    logger.info("\n" + "=" * 60)
    logger.info("PER-PEAK MODEL TRAINING")
    logger.info("=" * 60)

    # === LINEWIDTH MODELS ===
    logger.info(f"\n--- Linewidth Models (min R² = {min_r2_linewidth}) ---")

    # Flatten peaks with quality filter for linewidth training
    peaks_15N_lw, peaks_13C_lw = flatten_peak_data(
        samples, min_r2=min_r2_linewidth, require_success=True
    )

    logger.info(f"  15N peaks for linewidth training: {len(peaks_15N_lw)}")
    logger.info(f"  13C peaks for linewidth training: {len(peaks_13C_lw)}")

    # Train 15N linewidth model
    if len(peaks_15N_lw) >= min_samples:
        X, y = build_peak_matrices(peaks_15N_lw, target_type='linewidth')
        logger.info(f"  15N linewidth matrix: X={X.shape}, y={y.shape}")

        # Split into train/test if requested
        if test_split > 0:
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=test_split, random_state=42
            )
            logger.info(f"  Train/test split: {len(X_train)}/{len(X_test)} peaks")
        else:
            X_train, y_train = X, y
            X_test, y_test = None, None

        model, cv_mean, cv_std = train_model(X_train, y_train, "15N peak-linewidth", cv_folds)

        # Evaluate on test set
        if X_test is not None:
            test_r2 = evaluate_on_test_set(model, X_test, y_test, "15N linewidth")
        else:
            test_r2 = None

        # Show feature importance
        if show_feature_importance:
            display_feature_importance(model, PEAK_FEATURE_NAMES, "15N linewidth")

        model_path = output_dir / "peak_lw_model_15N.joblib"
        joblib.dump(model, model_path)
        logger.info(f"  Saved: {model_path}")

        results['linewidth_models']['15N'] = {
            'model_path': str(model_path),
            'n_samples': len(peaks_15N_lw),
            'cv_r2_mean': cv_mean,
            'cv_r2_std': cv_std,
            'test_r2': test_r2,
            'n_features': X.shape[1],
            'n_targets': y.shape[1] if y.ndim > 1 else 1,
        }
    else:
        logger.warning(f"  Skipping 15N linewidth model: {len(peaks_15N_lw)} < {min_samples} peaks")

    # Train 13C linewidth model
    if len(peaks_13C_lw) >= min_samples:
        X, y = build_peak_matrices(peaks_13C_lw, target_type='linewidth')
        logger.info(f"  13C linewidth matrix: X={X.shape}, y={y.shape}")

        # Split into train/test if requested
        if test_split > 0:
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=test_split, random_state=42
            )
            logger.info(f"  Train/test split: {len(X_train)}/{len(X_test)} peaks")
        else:
            X_train, y_train = X, y
            X_test, y_test = None, None

        model, cv_mean, cv_std = train_model(X_train, y_train, "13C peak-linewidth", cv_folds)

        # Evaluate on test set
        if X_test is not None:
            test_r2 = evaluate_on_test_set(model, X_test, y_test, "13C linewidth")
        else:
            test_r2 = None

        # Show feature importance
        if show_feature_importance:
            display_feature_importance(model, PEAK_FEATURE_NAMES, "13C linewidth")

        model_path = output_dir / "peak_lw_model_13C.joblib"
        joblib.dump(model, model_path)
        logger.info(f"  Saved: {model_path}")

        results['linewidth_models']['13C'] = {
            'model_path': str(model_path),
            'n_samples': len(peaks_13C_lw),
            'cv_r2_mean': cv_mean,
            'cv_r2_std': cv_std,
            'test_r2': test_r2,
            'n_features': X.shape[1],
            'n_targets': y.shape[1] if y.ndim > 1 else 1,
        }
    else:
        logger.warning(f"  Skipping 13C linewidth model: {len(peaks_13C_lw)} < {min_samples} peaks")

    # === QUALITY MODELS ===
    logger.info(f"\n--- Quality Models (min R² = {min_r2_quality}) ---")

    # Flatten peaks for quality training (include all R² values as targets)
    peaks_15N_q, peaks_13C_q = flatten_peak_data(
        samples, min_r2=min_r2_quality, require_success=True
    )

    logger.info(f"  15N peaks for quality training: {len(peaks_15N_q)}")
    logger.info(f"  13C peaks for quality training: {len(peaks_13C_q)}")

    # Train 15N quality model
    if len(peaks_15N_q) >= min_samples:
        X, y = build_peak_matrices(peaks_15N_q, target_type='quality')
        logger.info(f"  15N quality matrix: X={X.shape}, y={y.shape}")

        # Split into train/test if requested
        if test_split > 0:
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=test_split, random_state=42
            )
            logger.info(f"  Train/test split: {len(X_train)}/{len(X_test)} peaks")
        else:
            X_train, y_train = X, y
            X_test, y_test = None, None

        model, cv_mean, cv_std = train_model(X_train, y_train, "15N peak-quality", cv_folds)

        # Evaluate on test set
        if X_test is not None:
            test_r2 = evaluate_on_test_set(model, X_test, y_test, "15N quality")
        else:
            test_r2 = None

        # Show feature importance
        if show_feature_importance:
            display_feature_importance(model, PEAK_FEATURE_NAMES, "15N quality")

        model_path = output_dir / "peak_quality_model_15N.joblib"
        joblib.dump(model, model_path)
        logger.info(f"  Saved: {model_path}")

        results['quality_models']['15N'] = {
            'model_path': str(model_path),
            'n_samples': len(peaks_15N_q),
            'cv_r2_mean': cv_mean,
            'cv_r2_std': cv_std,
            'test_r2': test_r2,
            'n_features': X.shape[1],
            'n_targets': y.shape[1] if y.ndim > 1 else 1,
        }
    else:
        logger.warning(f"  Skipping 15N quality model: {len(peaks_15N_q)} < {min_samples} peaks")

    # Train 13C quality model
    if len(peaks_13C_q) >= min_samples:
        X, y = build_peak_matrices(peaks_13C_q, target_type='quality')
        logger.info(f"  13C quality matrix: X={X.shape}, y={y.shape}")

        # Split into train/test if requested
        if test_split > 0:
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=test_split, random_state=42
            )
            logger.info(f"  Train/test split: {len(X_train)}/{len(X_test)} peaks")
        else:
            X_train, y_train = X, y
            X_test, y_test = None, None

        model, cv_mean, cv_std = train_model(X_train, y_train, "13C peak-quality", cv_folds)

        # Evaluate on test set
        if X_test is not None:
            test_r2 = evaluate_on_test_set(model, X_test, y_test, "13C quality")
        else:
            test_r2 = None

        # Show feature importance
        if show_feature_importance:
            display_feature_importance(model, PEAK_FEATURE_NAMES, "13C quality")

        model_path = output_dir / "peak_quality_model_13C.joblib"
        joblib.dump(model, model_path)
        logger.info(f"  Saved: {model_path}")

        results['quality_models']['13C'] = {
            'model_path': str(model_path),
            'n_samples': len(peaks_13C_q),
            'cv_r2_mean': cv_mean,
            'cv_r2_std': cv_std,
            'test_r2': test_r2,
            'n_features': X.shape[1],
            'n_targets': y.shape[1] if y.ndim > 1 else 1,
        }
    else:
        logger.warning(f"  Skipping 13C quality model: {len(peaks_13C_q)} < {min_samples} peaks")

    # Store statistics
    results['statistics'] = {
        'total_peaks_15N': len(peaks_15N_q),
        'total_peaks_13C': len(peaks_13C_q),
        'high_quality_peaks_15N': len(peaks_15N_lw),
        'high_quality_peaks_13C': len(peaks_13C_lw),
        'min_r2_linewidth': min_r2_linewidth,
        'min_r2_quality': min_r2_quality,
        'feature_names': PEAK_FEATURE_NAMES,
        'lw_target_names': PEAK_LW_TARGET_NAMES,
        'quality_target_names': PEAK_QUALITY_TARGET_NAMES,
    }

    return results


# =============================================================================
# MODEL EXPORT
# =============================================================================

def export_model(
    model: RandomForestRegressor,
    output_dir: Path,
    nucleus: str
) -> Path:
    """
    Save trained model to joblib file.

    Parameters
    ----------
    model : RandomForestRegressor
        Trained model
    output_dir : Path
        Output directory
    nucleus : str
        Nucleus type ("15N" or "13C")

    Returns
    -------
    Path
        Path to saved model file
    """
    if not SKLEARN_AVAILABLE:
        raise ImportError("joblib is required for export. Install with: pip install scikit-learn")

    model_file = output_dir / f'base_model_{nucleus}.joblib'
    joblib.dump(model, model_file)
    logger.info(f"  Model saved: {model_file}")
    return model_file


def export_metadata(
    output_dir: Path,
    samples_15N: int,
    samples_13C: int,
    cv_score_15N: Optional[float],
    cv_score_13C: Optional[float],
    min_r2: float
) -> Path:
    """
    Save training metadata to JSON file.

    Parameters
    ----------
    output_dir : Path
        Output directory
    samples_15N : int
        Number of 15N training samples
    samples_13C : int
        Number of 13C training samples
    cv_score_15N : float or None
        Cross-validation score for 15N model
    cv_score_13C : float or None
        Cross-validation score for 13C model
    min_r2 : float
        R² threshold used for filtering

    Returns
    -------
    Path
        Path to saved metadata file
    """
    metadata = {
        'version': '1.0.0',
        'created': datetime.now().isoformat(),
        'samples_15N': samples_15N,
        'samples_13C': samples_13C,
        'cv_score_15N': cv_score_15N,
        'cv_score_13C': cv_score_13C,
        'min_r2_threshold': min_r2,
        'feature_names': FEATURE_NAMES,
        'target_names': TARGET_NAMES,
        'model_params': DEFAULT_MODEL_PARAMS,
    }

    metadata_file = output_dir / 'training_metadata.json'
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)

    logger.info(f"Metadata saved: {metadata_file}")
    return metadata_file


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Main entry point for CLI."""
    parser = argparse.ArgumentParser(
        description='Train ML models for lunaNMR parameter prediction',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--data-file',
        nargs='+',
        default=['ml_training_data/training_data.json'],
        help='Training data JSON file(s). Accepts: single file, multiple files, '
             'or glob pattern (e.g., "data/*.json"). Default: ml_training_data/training_data.json'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('lunaNMR/ml/pretrained'),
        help='Directory to save trained models (default: lunaNMR/ml/pretrained)'
    )
    parser.add_argument(
        '--min-r2',
        type=float,
        default=0.80,
        help='Minimum R² threshold for quality filtering (default: 0.80)'
    )
    parser.add_argument(
        '--min-samples',
        type=int,
        default=20,
        help='Minimum samples required to train a model (default: 20)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show statistics without training'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )

    # Per-peak training arguments
    parser.add_argument(
        '--train-peak-models',
        action='store_true',
        help='Enable per-peak model training (linewidth + quality)'
    )
    parser.add_argument(
        '--peak-lw-min-r2',
        type=float,
        default=0.85,
        help='Minimum R² for peaks used in linewidth training (default: 0.85)'
    )
    parser.add_argument(
        '--peak-min-samples',
        type=int,
        default=100,
        help='Minimum peaks required to train a peak model (default: 100)'
    )
    parser.add_argument(
        '--skip-spectrum-models',
        action='store_true',
        help='Skip spectrum-level models, train only peak models'
    )

    # Evaluation arguments
    parser.add_argument(
        '--test-split',
        type=float,
        default=0.0,
        help='Fraction of data to hold out for final evaluation (e.g., 0.2 for 20%%)'
    )
    parser.add_argument(
        '--show-feature-importance',
        action='store_true',
        help='Display feature importance rankings after training'
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if not SKLEARN_AVAILABLE:
        logger.error("ERROR: scikit-learn is required for training.")
        logger.error("Install with: pip install scikit-learn")
        return 1

    try:
        # Load training data
        # Handle single file, multiple files, or glob pattern
        data_source = args.data_file[0] if len(args.data_file) == 1 else args.data_file
        samples = load_training_data(data_source)

        # Filter for spectrum-level training
        good_samples = filter_quality_samples(samples, args.min_r2)
        samples_15N, samples_13C = split_by_nucleus(good_samples)

        # Dry run: show statistics
        if args.dry_run:
            logger.info("\n[Dry run - no models trained]")
            if not args.skip_spectrum_models:
                logger.info(f"  Would train 15N spectrum model: {len(samples_15N)} samples")
                logger.info(f"  Would train 13C spectrum model: {len(samples_13C)} samples")
            if args.train_peak_models:
                peaks_15N_lw, peaks_13C_lw = flatten_peak_data(
                    samples, min_r2=args.peak_lw_min_r2, require_success=True
                )
                peaks_15N_q, peaks_13C_q = flatten_peak_data(
                    samples, min_r2=0.0, require_success=True
                )
                logger.info(f"  Would train 15N linewidth model: {len(peaks_15N_lw)} peaks")
                logger.info(f"  Would train 13C linewidth model: {len(peaks_13C_lw)} peaks")
                logger.info(f"  Would train 15N quality model: {len(peaks_15N_q)} peaks")
                logger.info(f"  Would train 13C quality model: {len(peaks_13C_q)} peaks")
            return 0

        # Create output directory
        args.output_dir.mkdir(parents=True, exist_ok=True)

        # Track results for metadata
        cv_score_15N = None
        cv_score_13C = None
        peak_results = None

        # === SPECTRUM-LEVEL MODELS ===
        if not args.skip_spectrum_models:
            logger.info("\n" + "=" * 60)
            logger.info("SPECTRUM-LEVEL MODEL TRAINING")
            logger.info("=" * 60)

            if len(samples_15N) >= args.min_samples:
                X, y = build_matrices(samples_15N)

                # Split into train/test if requested
                if args.test_split > 0:
                    X_train, X_test, y_train, y_test = train_test_split(
                        X, y, test_size=args.test_split, random_state=42
                    )
                    logger.info(f"  Train/test split: {len(X_train)}/{len(X_test)} samples")
                else:
                    X_train, y_train = X, y
                    X_test, y_test = None, None

                model, cv_score_15N, _ = train_model(X_train, y_train, '15N')

                # Evaluate on test set
                if X_test is not None:
                    evaluate_on_test_set(model, X_test, y_test, '15N spectrum')

                # Show feature importance
                if args.show_feature_importance:
                    display_feature_importance(model, FEATURE_NAMES, '15N spectrum')

                export_model(model, args.output_dir, '15N')
            else:
                logger.info(f"Skipping 15N model ({len(samples_15N)} < {args.min_samples} samples)")

            if len(samples_13C) >= args.min_samples:
                X, y = build_matrices(samples_13C)

                # Split into train/test if requested
                if args.test_split > 0:
                    X_train, X_test, y_train, y_test = train_test_split(
                        X, y, test_size=args.test_split, random_state=42
                    )
                    logger.info(f"  Train/test split: {len(X_train)}/{len(X_test)} samples")
                else:
                    X_train, y_train = X, y
                    X_test, y_test = None, None

                model, cv_score_13C, _ = train_model(X_train, y_train, '13C')

                # Evaluate on test set
                if X_test is not None:
                    evaluate_on_test_set(model, X_test, y_test, '13C spectrum')

                # Show feature importance
                if args.show_feature_importance:
                    display_feature_importance(model, FEATURE_NAMES, '13C spectrum')

                export_model(model, args.output_dir, '13C')
            else:
                logger.info(f"Skipping 13C model ({len(samples_13C)} < {args.min_samples} samples)")

        # === PER-PEAK MODELS ===
        if args.train_peak_models:
            peak_results = train_peak_models(
                samples,
                args.output_dir,
                min_r2_linewidth=args.peak_lw_min_r2,
                min_r2_quality=0.0,
                cv_folds=5,
                min_samples=args.peak_min_samples,
                test_split=args.test_split,
                show_feature_importance=args.show_feature_importance,
            )

        # Export metadata (with peak model info if trained)
        export_metadata(
            args.output_dir,
            len(samples_15N) if not args.skip_spectrum_models else 0,
            len(samples_13C) if not args.skip_spectrum_models else 0,
            cv_score_15N,
            cv_score_13C,
            args.min_r2
        )

        # Append peak model stats to metadata if trained
        if peak_results:
            metadata_file = args.output_dir / 'training_metadata.json'
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
            metadata['peak_models'] = {
                'linewidth': peak_results['linewidth_models'],
                'quality': peak_results['quality_models'],
                'statistics': peak_results['statistics'],
            }
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)

        logger.info("\nTraining complete!")
        logger.info(f"  Models saved to: {args.output_dir}")

        return 0

    except FileNotFoundError as e:
        logger.error(f"ERROR: {e}")
        return 1
    except Exception as e:
        logger.error(f"ERROR: Training failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    exit(main())
