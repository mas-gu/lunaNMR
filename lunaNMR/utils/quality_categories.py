#!/usr/bin/env python3
"""
Quality Categorization System for LunaNMR
==========================================

Centralized system for R² quality thresholds with support for:
- Scientific quality categories (strict, traditional)
- ML-friendly quality categories (relaxed for training diversity)
- Full retrocompatibility with existing code
- Configurable thresholds via JSON config

Author: LunaNMR Team
Version: 1.0
"""

import numpy as np
from typing import Dict, List, Union, Optional
from dataclasses import dataclass

@dataclass
class QualityThresholds:
    """Quality threshold configuration"""
    excellent: float
    good: float
    fair: float
    name: str
    description: str

class QualityCategorizer:
    """
    Centralized quality categorization system with multiple threshold sets.

    Features:
    - Multiple predefined threshold sets (scientific, ml_friendly)
    - Custom threshold support
    - Backward compatibility with existing code
    - Batch analysis with statistics
    """

    # Predefined threshold sets
    THRESHOLD_SETS = {
        'scientific': QualityThresholds(
            excellent=0.95,
            good=0.85,
            fair=0.70,
            name="Scientific",
            description="Traditional scientific quality thresholds (strict)"
        ),
        'ml_friendly': QualityThresholds(
            excellent=0.90,
            good=0.75,
            fair=0.60,
            name="ML-Friendly",
            description="Relaxed thresholds optimized for ML training diversity"
        ),
        'permissive': QualityThresholds(
            excellent=0.85,
            good=0.70,
            fair=0.55,
            name="Permissive",
            description="Very relaxed thresholds for maximum data inclusion"
        )
    }

    def __init__(self, threshold_set: str = 'scientific', custom_thresholds: Optional[Dict] = None):
        """
        Initialize quality categorizer.

        Args:
            threshold_set: Name of predefined threshold set or 'custom'
            custom_thresholds: Dict with 'excellent', 'good', 'fair' keys (if threshold_set='custom')
        """
        if threshold_set == 'custom' and custom_thresholds:
            self.thresholds = QualityThresholds(
                excellent=custom_thresholds['excellent'],
                good=custom_thresholds['good'],
                fair=custom_thresholds['fair'],
                name="Custom",
                description="User-defined custom thresholds"
            )
        elif threshold_set in self.THRESHOLD_SETS:
            self.thresholds = self.THRESHOLD_SETS[threshold_set]
        else:
            # Fallback to scientific for backward compatibility
            self.thresholds = self.THRESHOLD_SETS['scientific']

        self.threshold_set_name = threshold_set

    def categorize_r_squared(self, r_squared: float) -> str:
        """
        Categorize a single R² value.

        Args:
            r_squared: R² value to categorize

        Returns:
            Category string: 'excellent', 'good', 'fair', or 'poor'
        """
        if r_squared >= self.thresholds.excellent:
            return 'excellent'
        elif r_squared >= self.thresholds.good:
            return 'good'
        elif r_squared >= self.thresholds.fair:
            return 'fair'
        else:
            return 'poor'

    def categorize_batch(self, r_squared_values: Union[List[float], np.ndarray]) -> Dict:
        """
        Categorize a batch of R² values with full statistics.

        Args:
            r_squared_values: Array or list of R² values

        Returns:
            Dict with categories, counts, percentages, and statistics
        """
        r_squared_values = np.array(r_squared_values)

        # Calculate category counts
        excellent_mask = r_squared_values >= self.thresholds.excellent
        good_mask = (r_squared_values >= self.thresholds.good) & (r_squared_values < self.thresholds.excellent)
        fair_mask = (r_squared_values >= self.thresholds.fair) & (r_squared_values < self.thresholds.good)
        poor_mask = r_squared_values < self.thresholds.fair

        excellent = np.sum(excellent_mask)
        good = np.sum(good_mask)
        fair = np.sum(fair_mask)
        poor = np.sum(poor_mask)
        total = len(r_squared_values)

        return {
            'threshold_set': self.threshold_set_name,
            'threshold_info': {
                'excellent': f"≥{self.thresholds.excellent:.2f}",
                'good': f"{self.thresholds.good:.2f}-{self.thresholds.excellent:.2f}",
                'fair': f"{self.thresholds.fair:.2f}-{self.thresholds.good:.2f}",
                'poor': f"<{self.thresholds.fair:.2f}"
            },
            'statistics': {
                'total_samples': total,
                'average_r_squared': float(np.mean(r_squared_values)),
                'median_r_squared': float(np.median(r_squared_values)),
                'std_r_squared': float(np.std(r_squared_values)),
                'min_r_squared': float(np.min(r_squared_values)),
                'max_r_squared': float(np.max(r_squared_values))
            },
            'categories': {
                'excellent': {
                    'count': int(excellent),
                    'percentage': float(excellent / total * 100) if total > 0 else 0.0,
                    'threshold': f"≥{self.thresholds.excellent:.2f}"
                },
                'good': {
                    'count': int(good),
                    'percentage': float(good / total * 100) if total > 0 else 0.0,
                    'threshold': f"{self.thresholds.good:.2f}-{self.thresholds.excellent:.2f}"
                },
                'fair': {
                    'count': int(fair),
                    'percentage': float(fair / total * 100) if total > 0 else 0.0,
                    'threshold': f"{self.thresholds.fair:.2f}-{self.thresholds.good:.2f}"
                },
                'poor': {
                    'count': int(poor),
                    'percentage': float(poor / total * 100) if total > 0 else 0.0,
                    'threshold': f"<{self.thresholds.fair:.2f}"
                }
            },
            'ml_suitability': self._assess_ml_suitability(r_squared_values)
        }

    def _assess_ml_suitability(self, r_squared_values: np.ndarray) -> Dict:
        """Assess suitability of data for ML training"""
        suitable = np.sum(r_squared_values >= self.thresholds.fair)  # fair+ quality
        marginal = np.sum((r_squared_values >= self.thresholds.fair * 0.8) & (r_squared_values < self.thresholds.fair))
        unsuitable = np.sum(r_squared_values < self.thresholds.fair * 0.8)
        total = len(r_squared_values)

        return {
            'suitable_for_training': {
                'count': int(suitable),
                'percentage': float(suitable / total * 100) if total > 0 else 0.0,
                'description': f"R² ≥ {self.thresholds.fair:.2f}"
            },
            'marginal_for_training': {
                'count': int(marginal),
                'percentage': float(marginal / total * 100) if total > 0 else 0.0,
                'description': f"R² {self.thresholds.fair*0.8:.2f}-{self.thresholds.fair:.2f}"
            },
            'unsuitable_for_training': {
                'count': int(unsuitable),
                'percentage': float(unsuitable / total * 100) if total > 0 else 0.0,
                'description': f"R² < {self.thresholds.fair*0.8:.2f}"
            }
        }


# Convenience functions for backward compatibility
def get_quality_categorizer(mode: str = 'scientific', config: Optional[Dict] = None) -> QualityCategorizer:
    """
    Factory function to create quality categorizer with different modes.

    Args:
        mode: 'scientific', 'ml_friendly', 'permissive', or 'auto'
        config: Optional configuration dict with quality_threshold_mode

    Returns:
        QualityCategorizer instance
    """
    # Check config for mode override
    if config and 'quality_threshold_mode' in config:
        mode = config['quality_threshold_mode']

    # Auto mode: choose based on context
    if mode == 'auto':
        if config and config.get('ml_data_collection_enabled', False):
            mode = 'ml_friendly'
        else:
            mode = 'scientific'

    return QualityCategorizer(threshold_set=mode)

def categorize_quality_legacy(r_squared: float) -> str:
    """Legacy function for backward compatibility (uses scientific thresholds)"""
    categorizer = QualityCategorizer('scientific')
    return categorizer.categorize_r_squared(r_squared)

# Export key functions for easy importing
__all__ = [
    'QualityCategorizer',
    'QualityThresholds',
    'get_quality_categorizer',
    'categorize_quality_legacy'
]