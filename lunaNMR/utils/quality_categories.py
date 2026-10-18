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
    - Multiple predefined threshold sets (scientific, permissive)
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
            }
        }


# Convenience functions for backward compatibility
def get_quality_categorizer(mode: str = 'scientific', config: Optional[Dict] = None) -> QualityCategorizer:
    """
    Factory function to create quality categorizer with different modes.

    Args:
        mode: 'scientific', 'permissive', or 'auto'
        config: Optional configuration dict with quality_threshold_mode

    Returns:
        QualityCategorizer instance
    """
    # Check config for mode override
    if config and 'quality_threshold_mode' in config:
        mode = config['quality_threshold_mode']

    if mode == 'auto':
        mode = 'scientific'

    return QualityCategorizer(threshold_set=mode)

def categorize_quality_legacy(r_squared: float) -> str:
    """Legacy function for backward compatibility (uses scientific thresholds)"""
    categorizer = QualityCategorizer('scientific')
    return categorizer.categorize_r_squared(r_squared)

# Export key functions for easy importing
# ---------------------------------------------------------------------------
# Which points of a series carry a measurement
# ---------------------------------------------------------------------------
# Below this many sigma a peak is not detected. The same floor holdPositionsAtNoise
# uses; 3 sigma is the usual bar for calling a signal present at all.
DETECTION_SIGMA_FLOOR = 3.0

# Observables that need a POSITION. A peak that has decayed has none, so it cannot
# contribute a chemical-shift perturbation - but it still contributes an intensity.
POSITION_OBSERVABLES = ('csp', 'shift', 'position')


def usable_points(frame, observable):
    """The rows of a series that carry a measurement of `observable`.

    Two states look identical in a series and need opposite handling, which is the
    whole reason this exists.

    A residue whose peak has DECAYED reports a near-zero intensity, and that zero IS
    the measurement: dropping it from an I/I0 series invents a plateau that is not
    there. It has no position, though, so it cannot contribute a CSP.

    A residue whose fit stopped on its INTENSITY floor reports a number that came from
    the bound rather than from the data. The floor is cluster-relative (0.1% of the
    brightest member's seed), so that number is a product of its NEIGHBOURS' seeds and
    its own seeded linewidths, carrying nothing from the data at its own position. It
    is not a measurement of anything and is dropped from every observable.

    A residue whose POSITION stopped on its margin has no measured position, so it is
    dropped from the position observables. Its intensity is kept: on a 600 MHz hetNOE pair
    most such peaks sit on the 0.020 ppm 15N margin, a quarter of one 15N point, and
    widening it to 0.05 ppm moved the freed peaks' heights by a median 0.0%.

    Neither state is visible in `success` or in R-squared. `success` is an OR whose
    first arm is the optimiser's "chi2 stopped changing", true immediately on flat
    noise; R-squared is per CLUSTER, so a member with nothing at its position rides the
    one its healthy neighbours earn - measured, 0.58 sigma under a cluster R-squared of
    0.9884.

    A row from before these flags existed carries NaN. It is kept: refusing every
    legacy series would be worse than the risk it represents, and the flags are absent
    rather than negative. Re-run the series to get the filter.

    :param frame: a DataFrame of series rows, one per (spectrum, peak)
    :param observable: what is being measured; see POSITION_OBSERVABLES
    :return: the subset of `frame` that can be used, index preserved
    """
    keep = frame.index == frame.index                 # all True, right length
    positional = str(observable).lower() in POSITION_OBSERVABLES
    columns = ('position_pinned', 'intensity_pinned') if positional else ('intensity_pinned',)
    for column in columns:
        if column in frame:
            # `== True` and not truthiness: NaN must not read as pinned, and None
            # must not read as "not pinned" by accident either - both are unknown,
            # and unknown is kept.
            keep &= ~(frame[column].fillna(False).astype(bool))

    if positional and 'detection_sigma' in frame:
        detected = frame['detection_sigma'] >= DETECTION_SIGMA_FLOOR
        keep &= detected | frame['detection_sigma'].isna()

    return frame[keep]


__all__ = [
    'usable_points',
    'DETECTION_SIGMA_FLOOR',
    'QualityCategorizer',
    'QualityThresholds',
    'get_quality_categorizer',
    'categorize_quality_legacy'
]