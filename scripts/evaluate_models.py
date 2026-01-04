#!/usr/bin/env python3
# ABOUTME: Post-hoc evaluation script for trained ML models
# ABOUTME: Analyzes model performance, feature importance, and prediction quality

"""
Model evaluation script for lunaNMR ML models.

Loads trained models and evaluates them on test data, providing:
- Performance metrics (R², MAE, RMSE)
- Feature importance rankings
- Prediction vs actual analysis
- Error distribution

Usage:
    # Evaluate a spectrum-level model
    python scripts/evaluate_models.py lunaNMR/ml/pretrained/base_model_15N.joblib

    # Evaluate on specific test data
    python scripts/evaluate_models.py model.joblib --data test_data.json

    # Show top 15 features
    python scripts/evaluate_models.py model.joblib --top-features 15

    # Evaluate peak-level model
    python scripts/evaluate_models.py peak_lw_model_15N.joblib --model-type peak-linewidth
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import joblib
    from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

# Try matplotlib for plots
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


# =============================================================================
# CONSTANTS (from pretrain_models.py)
# =============================================================================

SPECTRUM_FEATURE_NAMES = [
    'nucleus_type', 'snr_estimate', 'noise_level', 'dynamic_range',
    'peak_count', 'peak_density', 'shift_range_f1_min', 'shift_range_f1_max',
    'shift_range_f2_min', 'shift_range_f2_max', 'field_strength_mhz',
    'baseline_level', 'baseline_std', 'mean_peak_separation_f1',
    'mean_peak_separation_f2', 'min_peak_separation_f1', 'min_peak_separation_f2',
    'f1_dispersion', 'f2_dispersion', 'is_idp_like',
    'detected_intensity_mean', 'detected_intensity_std', 'detected_intensity_cv',
    'detected_intensity_dynamic_range', 'detected_intensity_skewness',
    'detected_intensity_kurtosis', 'dispersion_ratio',
    'n_close_pairs_f1', 'n_close_pairs_f2', 'n_close_pairs_2d',
    'fraction_potentially_overlapping', 'peaklist_intensity_mean',
    'peaklist_intensity_std', 'peaklist_intensity_cv',
    'peaklist_intensity_dynamic_range', 'max_local_density', 'crowding_hotspot_fraction',
]

SPECTRUM_TARGET_NAMES = [
    'lw_f1_median', 'lw_f2_median', 'rad_f1', 'rad_f2',
    'overlap_threshold_f1', 'overlap_threshold_f2', 'achievable_r2',
]

PEAK_FEATURE_NAMES = [
    'is_isolated', 'cluster_size', 'tooclose_flag', 'heavy_overlap_flag',
    'relative_pos_f1', 'relative_pos_f2', 'relative_intensity', 'log_intensity',
    'fix_positions_applied', 'fix_linewidths_applied', 'fitting_mode_2d',
    'nucleus_type', 'snr_estimate', 'noise_level_normalized', 'peak_density',
    'fraction_clustered', 'is_idp_like', 'field_strength_mhz',
    'pass1_lw_f1_median', 'pass1_lw_f2_median', 'fraction_potentially_overlapping',
    'f1_dispersion', 'f2_dispersion',
]

PEAK_LW_TARGET_NAMES = ['lw_total_f1', 'lw_total_f2']
PEAK_QUALITY_TARGET_NAMES = ['r_squared']


# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def load_model(model_path: Path):
    """Load a trained model from joblib file."""
    if not SKLEARN_AVAILABLE:
        raise ImportError("scikit-learn required. Install with: pip install scikit-learn")
    return joblib.load(model_path)


def get_feature_importance(model, feature_names: List[str], top_n: int = 10) -> List[Tuple[str, float]]:
    """Extract and rank feature importances."""
    importances = model.feature_importances_
    indices = np.argsort(importances)[::-1]

    ranked = []
    for i in range(min(top_n, len(feature_names))):
        idx = indices[i]
        name = feature_names[idx] if idx < len(feature_names) else f"feature_{idx}"
        ranked.append((name, importances[idx]))

    return ranked


def display_feature_importance(ranked: List[Tuple[str, float]], title: str = "Feature Importance"):
    """Display feature importance rankings."""
    print(f"\n{title}")
    print("─" * 55)

    for i, (name, importance) in enumerate(ranked, 1):
        bar_len = int(importance * 50)
        bar = '█' * bar_len
        print(f"{i:2d}. {name:35s} {importance:.4f} {bar}")


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    """Compute evaluation metrics."""
    metrics = {
        'r2': r2_score(y_true, y_pred),
        'mae': mean_absolute_error(y_true, y_pred),
        'rmse': np.sqrt(mean_squared_error(y_true, y_pred)),
    }

    # Per-target metrics if multi-output
    if y_true.ndim > 1 and y_true.shape[1] > 1:
        for i in range(y_true.shape[1]):
            metrics[f'r2_target_{i}'] = r2_score(y_true[:, i], y_pred[:, i])
            metrics[f'mae_target_{i}'] = mean_absolute_error(y_true[:, i], y_pred[:, i])

    return metrics


def display_metrics(metrics: Dict[str, float], title: str = "Performance Metrics"):
    """Display evaluation metrics."""
    print(f"\n{title}")
    print("─" * 40)

    for key, value in metrics.items():
        print(f"  {key:20s}: {value:.4f}")


def analyze_errors(y_true: np.ndarray, y_pred: np.ndarray, percentiles: List[int] = [50, 90, 95, 99]):
    """Analyze prediction errors."""
    errors = np.abs(y_true - y_pred)

    if errors.ndim > 1:
        errors = errors.mean(axis=1)  # Average across targets

    print("\nError Distribution")
    print("─" * 40)
    print(f"  Mean absolute error: {errors.mean():.4f}")
    print(f"  Std of errors:       {errors.std():.4f}")
    print(f"  Min error:           {errors.min():.4f}")
    print(f"  Max error:           {errors.max():.4f}")

    print("\n  Percentiles:")
    for p in percentiles:
        val = np.percentile(errors, p)
        print(f"    {p}th percentile:    {val:.4f}")


def plot_predictions(y_true: np.ndarray, y_pred: np.ndarray, output_path: Optional[Path] = None):
    """Create prediction vs actual scatter plot."""
    if not MATPLOTLIB_AVAILABLE:
        print("  [matplotlib not available, skipping plot]")
        return

    if y_true.ndim > 1:
        n_targets = y_true.shape[1]
        fig, axes = plt.subplots(1, n_targets, figsize=(5 * n_targets, 5))
        if n_targets == 1:
            axes = [axes]

        for i, ax in enumerate(axes):
            ax.scatter(y_true[:, i], y_pred[:, i], alpha=0.3, s=10)
            lims = [
                min(y_true[:, i].min(), y_pred[:, i].min()),
                max(y_true[:, i].max(), y_pred[:, i].max())
            ]
            ax.plot(lims, lims, 'r--', alpha=0.8)
            ax.set_xlabel(f'Actual (target {i})')
            ax.set_ylabel(f'Predicted (target {i})')
            ax.set_title(f'Target {i}: R² = {r2_score(y_true[:, i], y_pred[:, i]):.3f}')
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(y_true, y_pred, alpha=0.3, s=10)
        lims = [min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())]
        ax.plot(lims, lims, 'r--', alpha=0.8)
        ax.set_xlabel('Actual')
        ax.set_ylabel('Predicted')
        ax.set_title(f'R² = {r2_score(y_true, y_pred):.3f}')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"\n  Plot saved to: {output_path}")
    else:
        plt.show()


def load_test_data(data_path: Path, model_type: str) -> Tuple[np.ndarray, np.ndarray]:
    """Load test data and extract features/targets based on model type."""
    # Import extraction functions from pretrain_models
    sys.path.insert(0, str(Path(__file__).parent))
    from pretrain_models import (
        load_training_data, filter_quality_samples, split_by_nucleus,
        build_matrices, flatten_peak_data, build_peak_matrices
    )

    samples = load_training_data(str(data_path))

    if model_type == 'spectrum':
        good_samples = filter_quality_samples(samples, min_r2=0.8)
        samples_15N, samples_13C = split_by_nucleus(good_samples)
        # Use 15N by default
        if samples_15N:
            X, y = build_matrices(samples_15N)
        else:
            X, y = build_matrices(samples_13C)
        return X, y

    elif model_type == 'peak-linewidth':
        peaks_15N, peaks_13C = flatten_peak_data(samples, min_r2=0.85, require_success=True)
        if peaks_15N:
            X, y = build_peak_matrices(peaks_15N, target_type='linewidth')
        else:
            X, y = build_peak_matrices(peaks_13C, target_type='linewidth')
        return X, y

    elif model_type == 'peak-quality':
        peaks_15N, peaks_13C = flatten_peak_data(samples, min_r2=0.0, require_success=True)
        if peaks_15N:
            X, y = build_peak_matrices(peaks_15N, target_type='quality')
        else:
            X, y = build_peak_matrices(peaks_13C, target_type='quality')
        return X, y

    else:
        raise ValueError(f"Unknown model type: {model_type}")


def infer_model_type(model_path: Path) -> Tuple[str, List[str]]:
    """Infer model type from filename."""
    name = model_path.stem.lower()

    if 'peak_lw' in name or 'peak_linewidth' in name:
        return 'peak-linewidth', PEAK_FEATURE_NAMES
    elif 'peak_quality' in name:
        return 'peak-quality', PEAK_FEATURE_NAMES
    elif 'base_model' in name:
        return 'spectrum', SPECTRUM_FEATURE_NAMES
    else:
        # Default to spectrum
        return 'spectrum', SPECTRUM_FEATURE_NAMES


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Evaluate trained ML models for lunaNMR',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        'model',
        type=Path,
        help='Path to trained model (.joblib file)'
    )
    parser.add_argument(
        '--data', '-d',
        type=Path,
        help='Test data file (JSON). If not provided, only shows model info.'
    )
    parser.add_argument(
        '--model-type', '-t',
        choices=['spectrum', 'peak-linewidth', 'peak-quality'],
        help='Model type (auto-detected from filename if not specified)'
    )
    parser.add_argument(
        '--top-features', '-n',
        type=int,
        default=10,
        help='Number of top features to display (default: 10)'
    )
    parser.add_argument(
        '--plot', '-p',
        type=Path,
        help='Save prediction plot to file (requires matplotlib)'
    )
    parser.add_argument(
        '--json', '-j',
        action='store_true',
        help='Output results as JSON'
    )

    args = parser.parse_args()

    if not SKLEARN_AVAILABLE:
        print("ERROR: scikit-learn required. Install with: pip install scikit-learn")
        return 1

    if not args.model.exists():
        print(f"ERROR: Model file not found: {args.model}")
        return 1

    # Load model
    print(f"Loading model: {args.model}")
    model = load_model(args.model)

    # Infer model type if not specified
    if args.model_type:
        model_type = args.model_type
        feature_names = PEAK_FEATURE_NAMES if 'peak' in model_type else SPECTRUM_FEATURE_NAMES
    else:
        model_type, feature_names = infer_model_type(args.model)

    print(f"Model type: {model_type}")
    print(f"Number of features: {len(feature_names)}")
    print(f"Number of trees: {model.n_estimators}")

    # Feature importance
    ranked = get_feature_importance(model, feature_names, args.top_features)

    if args.json:
        results = {
            'model_path': str(args.model),
            'model_type': model_type,
            'n_features': len(feature_names),
            'n_estimators': model.n_estimators,
            'feature_importance': [{'name': n, 'importance': float(i)} for n, i in ranked],
        }
    else:
        display_feature_importance(ranked, f"Top {args.top_features} Features ({model_type})")

    # Evaluate on test data if provided
    if args.data:
        if not args.data.exists():
            print(f"\nERROR: Test data file not found: {args.data}")
            return 1

        print(f"\nLoading test data: {args.data}")
        X, y = load_test_data(args.data, model_type)
        print(f"Test samples: {len(X)}")

        # Make predictions
        y_pred = model.predict(X)

        # Compute metrics
        metrics = compute_metrics(y, y_pred)

        if args.json:
            results['test_data'] = str(args.data)
            results['n_test_samples'] = len(X)
            results['metrics'] = metrics
        else:
            display_metrics(metrics)
            analyze_errors(y, y_pred)

        # Plot if requested
        if args.plot:
            plot_predictions(y, y_pred, args.plot)

    if args.json:
        print(json.dumps(results, indent=2))

    return 0


if __name__ == '__main__':
    sys.exit(main())
