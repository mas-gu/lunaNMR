"""
PS2D-Style Multi-Peak Voigt Fitter for LunaNMR
==============================================

Date: 2025-10-04
Version: 1.0 - Guillaume Mas
"""

import numpy as np

from typing import Dict, List, Tuple


# ============================================================================
# MATHEMATICAL CONSTANTS 
# ============================================================================

SQRT_2 = 1.4142135624                    # √2
SQRT_2PI = np.sqrt(2.0 * np.pi)         # √(2π)
SQRT_8LN2 = np.sqrt(8.0 * np.log(2.0))  # √(8·ln(2)) = 2.35482...

# Derivative step size 
DERIV_STEP = 1.000001  # Relative 1e-6 step (0.0001%)

# Levenberg-Marquardt optimizer parameters 
LAMBDA_INIT = 0.001
LAMBDA_UP = 10.0
LAMBDA_DOWN = 0.1
MAX_ITER = 1000
CONV_TOL = 1e-6

# NPAR_VOIGT = 6 (2D), but we use 1D simplification
NPAR_VOIGT_1D = 4  # pos, lw_lor, lw_gauss, intensity (per peak) - NO baseline
NFIX_1D = 3        # pos, lw_lor, lw_gauss (fixed shape parameters)


# ============================================================================
# ============================================================================

class Ps2dLinewidthEstimator:
    """
    Apriori linewidth estimation following  logic with enhancements.

    Implements three-tier hierarchy matching 
    1. User overrides (explicit -LINEWIDTH assignments)
    2. Reference peak templates (fitted peaks become templates for subsequent fits)
    3. Automatic estimation (spatial analysis of well-resolved peaks, enhancement)
    4. Nucleus defaults 

    This is completely self-contained with NO dependencies on enhanced_voigt_fitter.py

    """

    def __init__(self):
        self.user_overrides = {}  # {assignment: {'dimension': {'lw_lorentz': ..., 'lw_gauss': ...}}}

        # Reference peaks (fitted peaks become templates)
        self.reference_peaks = {}  # {dimension: {'lw_lorentz': ..., 'lw_gauss': ...}}

        # Spatial analysis cache (automatic estimation from well-resolved peaks)
        self.spatial_estimates = {}  # {dimension: {'lw_lorentz': ..., 'lw_gauss': ..., 'quality': ...}}

        # Nucleus defaults (CORRECTED - 1H is narrow, 15N is broad)
        self.defaults = {
            '1H': {'lw_lorentz': 0.001, 'lw_gauss': 0.03},     # 1H: narrow peaks
            '15N': {'lw_lorentz': 0.0001, 'lw_gauss': 0.3},    # 15N: broad peaks
            '13C': {'lw_lorentz': 0.0001, 'lw_gauss': 0.3}     # 13C: similar to 15N
        }

        # Nucleus detection ranges
        self.nucleus_ranges = {
            '1H': (5.5, 12.0),
            '15N': (100.0, 140.0),
            '13C': (5.0, 50.0)
        }

    def detect_nucleus_from_ppm_range(self, ppm_min: float, ppm_max: float) -> str:
        """Detect nucleus type from PPM range (self-contained)"""
        ppm_center = (ppm_min + ppm_max) / 2

        if 5.5 <= ppm_center <= 12.0:
            return '1H'
        elif 100.0 <= ppm_center <= 140.0:
            return '15N'
        elif 5.0 <= ppm_center <= 50.0:
            return '13C'
        else:
            # Default to 1H if ambiguous
            return '1H'

    def register_fitted_peak(self, dimension: str, lw_lorentz: float, lw_gauss: float,
                            assignment: str = None):
        """
        Register a fitted peak as reference template 

        from reference peak 'p' parameter.

        Parameters
        ----------
        dimension : str
            Dimension identifier
        lw_lorentz : float
            Fitted Lorentzian FWHM (ppm)
        lw_gauss : float
            Fitted Gaussian FWHM (ppm)
        assignment : str, optional
            Peak assignment for tracking
        """
        self.reference_peaks[dimension] = {
            'lw_lorentz': lw_lorentz,
            'lw_gauss': lw_gauss,
            'assignment': assignment
        }

    def estimate_from_spatial_analysis(self, x_data: np.ndarray, y_data: np.ndarray,
                                      all_peak_positions: List[float],
                                      dimension: str = None) -> Dict:
        """
        Estimate linewidths from well-resolved peaks (enhancement, self-contained)

        Algorithm:
        1. Detect nucleus type from x_data range
        2. Identify isolated peaks (≥ 6× typical width separation)
        3. Estimate linewidth using FWHM method (robust for NMR)
        4. Use median of estimates (robust against outliers)
        5. Quality check: coefficient of variation < 50%

        Parameters
        ----------
        x_data : np.ndarray
            PPM values
        y_data : np.ndarray
            Intensity values
        all_peak_positions : List[float]
            All detected peak positions for context
        dimension : str, optional
            Dimension identifier (auto-detect if None)

        Returns
        -------
        dict : Estimation results with quality indicator
        """
        # Detect nucleus if not specified
        if dimension is None:
            dimension = self.detect_nucleus_from_ppm_range(x_data.min(), x_data.max())

        # Map x/y to nucleus type
        nucleus_map = {'x': '1H', 'y': '15N'}
        nucleus_type = nucleus_map.get(dimension, dimension)

        # Get typical width for this nucleus
        typical_width = self.defaults[nucleus_type]['lw_gauss']

        # Find isolated peaks (≥ 6× typical width separation from ALL neighbors)
        min_separation = typical_width * 6.0
        isolated_peaks = []

        for i, pos in enumerate(all_peak_positions):
            # Find all neighbors within threshold
            neighbors = [p for j, p in enumerate(all_peak_positions)
                        if j != i and abs(p - pos) < min_separation]

            if len(neighbors) == 0:  # Completely isolated
                isolated_peaks.append(pos)

        result = {
            'lw_lorentz': self.defaults[nucleus_type]['lw_lorentz'],
            'lw_gauss': self.defaults[nucleus_type]['lw_gauss'],
            'isolated_count': len(isolated_peaks),
            'quality': 'default',
            'method': 'nucleus_default'
        }

        # Need at least 2 isolated peaks for robust estimation
        if len(isolated_peaks) >= 2:
            linewidth_estimates = []

            # Analyze each isolated peak (up to 5 to avoid over-computation)
            for peak_pos in isolated_peaks[:5]:
                try:
                    # Extract local region around peak (±4× typical width)
                    window_half = typical_width * 4.0
                    mask = (x_data >= peak_pos - window_half) & (x_data <= peak_pos + window_half)

                    if np.sum(mask) < 10:  # Need sufficient points
                        continue

                    x_local = x_data[mask]
                    y_local = y_data[mask]

                    # FWHM-based estimation (most reliable for NMR)
                    baseline = np.median(y_local)  # Robust baseline
                    peak_height = y_local.max() - baseline

                    if peak_height > 0:
                        half_max = baseline + peak_height / 2.0
                        above_half = y_local > half_max

                        if np.sum(above_half) >= 3:
                            indices = np.where(above_half)[0]
                            fwhm = x_local[indices[-1]] - x_local[indices[0]]

                            # Sanity check: FWHM should be within reasonable range
                            if typical_width * 0.1 < fwhm < typical_width * 10:
                                linewidth_estimates.append(fwhm)

                except Exception:
                    continue

            # Calculate robust statistics if we got enough estimates
            if len(linewidth_estimates) >= 2:
                median_lw = np.median(linewidth_estimates)
                std_lw = np.std(linewidth_estimates)
                cv = std_lw / median_lw if median_lw > 0 else 1.0

                # Quality check: CV < 50% indicates consistent estimates
                if cv < 0.5:
                    result['lw_gauss'] = median_lw
                    result['lw_lorentz'] = median_lw / 300.0
                    result['quality'] = 'good'
                    result['method'] = 'spatial_analysis'
                    result['cv'] = cv
                    result['n_peaks_analyzed'] = len(linewidth_estimates)
                else:
                    result['quality'] = 'inconsistent'
                    result['cv'] = cv

        # Cache result
        self.spatial_estimates[dimension] = result

        return result

    def get_linewidth(self, dimension: str, assignment: str = None,
                     x_data: np.ndarray = None, y_data: np.ndarray = None,
                     all_peak_positions: List[float] = None,
                     nucleus_type_override: str = None) -> Tuple[float, float]:
        """
        Get linewidth following REVISED hierarchy (spatial analysis prioritized)

        Priority (REVISED for better automation):
        1. Assignment-specific user override (if specified)
        2. Spatial analysis from well-resolved peaks (if data provided and quality good)
        3. Reference peak template / GUI defaults (fallback)
        4. Nucleus defaults (final fallback)

        Parameters
        ----------
        dimension : str
            Dimension identifier ('1H', '15N', '13C', 'x', 'y')
        assignment : str, optional
            Peak assignment for user override lookup
        x_data, y_data : np.ndarray, optional
            Data for spatial analysis (if not in cache)
        all_peak_positions : List[float], optional
            Peak context for spatial analysis
        nucleus_type_override : str, optional
            Explicit nucleus type ('1H', '15N', '13C') to override dimension mapping

        Returns
        -------
        tuple : (lw_lorentz, lw_gauss) in ppm
        """
        # Map x/y to nucleus type (with override support for 13C)
        if nucleus_type_override:
            nucleus_type = nucleus_type_override
        else:
            nucleus_map = {'x': '1H', 'y': '15N'}
            nucleus_type = nucleus_map.get(dimension, dimension)

        if assignment and assignment in self.user_overrides:
            if dimension in self.user_overrides[assignment]:
                override = self.user_overrides[assignment][dimension]
                return override['lw_lorentz'], override['lw_gauss']

        # Tier 2: Spatial analysis (PRIORITIZED - automatic from well-resolved peaks)
        if (x_data is not None and y_data is not None and
            all_peak_positions is not None and len(all_peak_positions) > 1):

            # Check cache first
            if dimension not in self.spatial_estimates:
                self.estimate_from_spatial_analysis(x_data, y_data,
                                                    all_peak_positions, dimension)

            estimate = self.spatial_estimates.get(dimension, {})
            if estimate.get('quality') == 'good':
                # Spatial analysis succeeded - use it!
                return estimate['lw_lorentz'], estimate['lw_gauss']

        # Tier 3: Reference peak template / GUI defaults (fallback if spatial fails)
        if dimension in self.reference_peaks:
            ref = self.reference_peaks[dimension]
            return ref['lw_lorentz'], ref['lw_gauss']

        defaults = self.defaults.get(nucleus_type, self.defaults['1H'])
        return defaults['lw_lorentz'], defaults['lw_gauss']


# ============================================================================
# ============================================================================
# This global instance enables linewidth learning: fitted peaks become templates
_global_linewidth_estimator = Ps2dLinewidthEstimator()


def get_global_estimator() -> Ps2dLinewidthEstimator:
    """
    Get the global linewidth estimator instance.

    This enables  behavior where fitted peaks become reference templates
    for subsequent fits in the same session.

    Returns
    -------
    Ps2dLinewidthEstimator : Global estimator with accumulated knowledge
    """
    return _global_linewidth_estimator
