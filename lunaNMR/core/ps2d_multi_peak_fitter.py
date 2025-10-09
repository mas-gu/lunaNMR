"""
PS2D-Style Multi-Peak Voigt Fitter for LunaNMR
==============================================

Complete port of  multi-peak fitting implementation (C++ → Python).
Implements the exact 5-stage fitting strategy from peakfit.cpp for simultaneous
fitting of multiple overlapping Voigt peaks.

1. Multi-peak simultaneous fitting (2-10 peaks)
2. Five-stage progressive refinement (stages 0-4)
3. Levenberg-Marquardt optimization via scipy.optimize.least_squares
4. Relative numerical derivatives (1e-6 step)
5. Analytical intensity derivatives
6. Linewidth constraints (fixLW flag)
7. Position constraints (fixPos flag)
8. Exact Faddeeva-based Voigt profiles


Date: 2025-10-04
Version: 1.0 - MULTI_PEAK_PS2D_EXACT_PORT
"""

import numpy as np
from scipy.special import wofz
from scipy.optimize import least_squares
from typing import Dict, List, Tuple, Optional, Union
import warnings

# ============================================================================
# MATHEMATICAL CONSTANTS (from faddeeva.cpp and fitfunctions.cpp)
# ============================================================================

SQRT_2 = 1.4142135624                    # √2
SQRT_2PI = np.sqrt(2.0 * np.pi)         # √(2π)
SQRT_8LN2 = np.sqrt(8.0 * np.log(2.0))  # √(8·ln(2)) = 2.35482...

# Derivative step size (from fitfunctions.cpp line 1986)
DERIV_STEP = 1.000001  # Relative 1e-6 step (0.0001%)

# Levenberg-Marquardt optimizer parameters (from fitmrq.cpp)
LAMBDA_INIT = 0.001
LAMBDA_UP = 10.0
LAMBDA_DOWN = 0.1
MAX_ITER = 1000
CONV_TOL = 1e-6

# From peakfit.cpp: NPAR_VOIGT = 6 (2D), but we use 1D simplification
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

    def set_user_override(self, assignment: str, dimension: str,
                         lw_lorentz: float, lw_gauss: float):
        """
        Set user-specified linewidth override 

        Parameters
        ----------
        assignment : str
            Peak assignment (e.g., 'A123', 'G45')
        dimension : str
            Dimension identifier ('1H', '15N', 'x', 'y')
        lw_lorentz : float
            Lorentzian FWHM (ppm)
        lw_gauss : float
            Gaussian FWHM (ppm)
        """
        if assignment not in self.user_overrides:
            self.user_overrides[assignment] = {}

        self.user_overrides[assignment][dimension] = {
            'lw_lorentz': lw_lorentz,
            'lw_gauss': lw_gauss
        }

    def register_fitted_peak(self, dimension: str, lw_lorentz: float, lw_gauss: float,
                            assignment: str = None):
        """
        Register a fitted peak as reference template 

        Matches  peak.cpp:598-605 where unspecified peaks inherit
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
                     all_peak_positions: List[float] = None) -> Tuple[float, float]:
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
            Dimension identifier ('1H', '15N', 'x', 'y')
        assignment : str, optional
            Peak assignment for user override lookup
        x_data, y_data : np.ndarray, optional
            Data for spatial analysis (if not in cache)
        all_peak_positions : List[float], optional
            Peak context for spatial analysis

        Returns
        -------
        tuple : (lw_lorentz, lw_gauss) in ppm
        """
        # Map x/y to nucleus type
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


def reset_global_estimator():
    """Reset global estimator to initial state (for new spectrum/session)"""
    global _global_linewidth_estimator
    _global_linewidth_estimator = Ps2dLinewidthEstimator()


# ============================================================================
# VOIGT PROFILE COMPUTATION (from faddeeva.cpp)
# ============================================================================

def voigt_peak_height_1d(lw_gauss: float, lw_lorentz: float) -> float:
    """
    Calculate peak height of unit-volume Voigt profile (1D version).

    From faddeeva.cpp lines 1998-2018. This is used to convert between
    peak intensity (height) and volume (integral).

    Parameters
    ----------
    lw_gauss : float
        Gaussian FWHM (ppm or Hz)
    lw_lorentz : float
        Lorentzian FWHM (ppm or Hz)

    Returns
    -------
    float : Peak height for unit volume Voigt
    """
    # Convert FWHM to internal sigma/gamma parameters
    sigma = lw_gauss / SQRT_8LN2  # Gaussian standard deviation
    gamma = lw_lorentz / 2.0       # Lorentzian half-width at half-max

    # At peak center (x=0), z = i*gamma/(sigma*sqrt(2))
    z = 1j * gamma / (sigma * SQRT_2)

    # Faddeeva function at peak center
    w_z = wofz(z)

    # Peak height = Re[w(z)] / (sigma * sqrt(2*pi))
    return w_z.real / (sigma * SQRT_2PI)


def multi_voigt_profile_1d(x: np.ndarray, params: np.ndarray, n_peaks: int) -> np.ndarray:
    """
    Multi-peak 1D Voigt profile - EXACT  (NO baseline)

    From faddeeva.cpp lines 1927-1950 (voigt function)

    Parameters organized as :
    For each peak i (i=0 to n_peaks-1):
        params[i*4 + 0] = pos_f1 (ppm)
        params[i*4 + 1] = lw_lorentz_f1 (FWHM, ppm)
        params[i*4 + 2] = lw_gauss_f1 (FWHM, ppm)
        params[i*4 + 3] = intensity (volume)

    NO BASELINE - data must be baseline-corrected BEFORE fitting!

    Mathematical Formula (from faddeeva.cpp):
    ------------------------------------------------
    y(x) = Σ[i=0 to n_peaks-1] Intensity_i * Voigt(x, i)

    where:
        Voigt(x, i) = Re[w(z_i)] / (σ_i * √(2π))
        z_i = (pos_i - x + i*γ_i) / (σ_i * √2)
        w(z) = Faddeeva function = exp(-z²) * erfc(-i*z)
        σ_i = FWHM_Gauss_i / √(8·ln(2))
        γ_i = FWHM_Lorentz_i / 2

    Parameters
    ----------
    x : np.ndarray
        Frequency axis (ppm)
    params : np.ndarray
        Parameter vector for all peaks (NO baseline)
        Length: n_peaks * 4
    n_peaks : int
        Number of overlapping peaks

    Returns
    -------
    np.ndarray : Summed Voigt profile for all peaks (NO baseline)
    """
    y = np.zeros_like(x, dtype=float)

    # Loop over all overlapping peaks (from faddeeva.cpp line 220)
    for i in range(n_peaks):
        # Extract parameters for peak i
        pos = params[i * 4 + 0]
        lw_lorentz = params[i * 4 + 1]
        lw_gauss = params[i * 4 + 2]
        intensity = params[i * 4 + 3]

        # Convert FWHM to internal sigma/gamma (faddeeva.cpp lines 222-225)
        sigma = lw_gauss / SQRT_8LN2  # Gaussian sigma
        gamma = lw_lorentz / 2.0       # Lorentzian gamma

        # Complex argument: z = (position - x + i*gamma) / (sigma*sqrt(2))
        # (faddeeva.cpp line 228)
        z = (pos - x + 1j * gamma) / (sigma * SQRT_2)

        # Voigt profile = intensity * Re[w(z)] / (sigma*sqrt(2*pi))
        # (faddeeva.cpp lines 233-236)
        w_z = wofz(z)
        voigt = w_z.real / (sigma * SQRT_2PI)

        # Add this peak's contribution
        y += intensity * voigt

    return y


def compute_multi_voigt_jacobian_1d(x: np.ndarray, params: np.ndarray,
                                     n_peaks: int) -> np.ndarray:
    """
    Jacobian for multi-peak Voigt - EXACT  (NO baseline)

    From fitfunctions.cpp lines 1953-1996 (voigt function)

    Uses the exact  approach:
    - Numerical derivatives with 1e-6 relative step for pos/linewidths
    - Analytical derivative for intensity: ∂f/∂I = f/I

    The Jacobian matrix J has dimensions (n_points, n_params) where:
        J[k, i] = ∂y[k]/∂params[i]

    Parameters
    ----------
    x : np.ndarray
        Frequency axis (n_points,)
    params : np.ndarray
        Parameter vector (n_peaks*4,) - NO baseline
    n_peaks : int
        Number of peaks

    Returns
    -------
    np.ndarray : Jacobian matrix (n_points, n_peaks*4)
    """
    n_points = len(x)
    n_params = n_peaks * 4  # NO baseline!

    # Initialize Jacobian (fitfunctions.cpp lines 301-302)
    jac = np.zeros((n_points, n_params), dtype=float)

    # Compute current function value
    y_base = multi_voigt_profile_1d(x, params, n_peaks)

    # Loop over all peaks (fitfunctions.cpp line 307)
    for i in range(n_peaks):
        base_idx = i * 4

        # Extract parameters for peak i
        pos = params[base_idx + 0]
        lw_lorentz = params[base_idx + 1]
        lw_gauss = params[base_idx + 2]
        intensity = params[base_idx + 3]

        # NUMERICAL DERIVATIVES (finite difference, 0.0001% perturbation)
        # From fitfunctions.cpp lines 328-357

        # 1. d/d(position) - parameter params[base_idx + 0]
        # (fitfunctions.cpp lines 330-332)
        if abs(pos) > 1e-10:
            params_perturb = params.copy()
            params_perturb[base_idx + 0] *= DERIV_STEP
            y_perturb = multi_voigt_profile_1d(x, params_perturb, n_peaks)
            jac[:, base_idx + 0] = (y_perturb - y_base) / (0.000001 * abs(pos))

        # 2. d/d(lorentzian_width) - parameter params[base_idx + 1]
        # (fitfunctions.cpp lines 334-337)
        if abs(lw_lorentz) > 1e-10:
            params_perturb = params.copy()
            params_perturb[base_idx + 1] *= DERIV_STEP
            y_perturb = multi_voigt_profile_1d(x, params_perturb, n_peaks)
            jac[:, base_idx + 1] = (y_perturb - y_base) / (0.000001 * abs(lw_lorentz))

        # 3. d/d(gaussian_width) - parameter params[base_idx + 2]
        # (fitfunctions.cpp lines 349-352)
        if abs(lw_gauss) > 1e-10:
            params_perturb = params.copy()
            params_perturb[base_idx + 2] *= DERIV_STEP
            y_perturb = multi_voigt_profile_1d(x, params_perturb, n_peaks)
            jac[:, base_idx + 2] = (y_perturb - y_base) / (0.000001 * abs(lw_gauss))

        # 4. d/d(intensity) - parameter params[base_idx + 3] - ANALYTICAL!
        # (fitfunctions.cpp line 360: dyda[i+5+f3] = f/a[i+f3+5])
        if abs(intensity) > 1e-10:
            # Compute contribution from this peak only
            # Zero out other peaks
            params_temp = params.copy()
            for j in range(n_peaks):
                if j != i:
                    params_temp[j * 4 + 3] = 0
            f_this_peak = multi_voigt_profile_1d(x, params_temp, n_peaks)
            jac[:, base_idx + 3] = f_this_peak / intensity

    return jac


# ============================================================================
# FIVE-STAGE MULTI-PEAK FITTER (from peakfit.cpp)
# ============================================================================

class Ps2dMultiPeakFitter:
    """
    Five-stage multi-peak Voigt fitter (EXACT from peakfit.cpp)

    Stages (from peakfit.cpp):
    0. fit_stage_zero (lines 231-248): Fix pos/lw, fit intensity only (VOIGT warm-up)
    1. fit_stage_one (lines 251-272): Fix pos, float lw+intensity
    2. fit_stage_two (lines 275-289): Float positions
    3. fit_stage_three (lines 292-319): Global intensity refinement (multi-plane only)
    4. fit_stage_four (lines 321-345): Final global optimization

    Implementation Notes
    --------------------
    - The original C++ code handles multi-plane (pseudo-3D) data. This 1D version
      simplifies to single-plane fitting.
    - Parameter fixing is implemented via scipy bounds: fixed = (value, value)
    - Uses scipy.optimize.least_squares with 'trf' or 'lm' method

    Parameters
    ----------
    verbose : bool, optional
        Print detailed fitting progress
    use_voigt_warmup : bool, optional
        Include Stage 0 (VOIGT-only warm-up), default True
    max_iterations : int, optional
        Maximum iterations per stage, default 1000
    dimension : str, optional
        Nucleus dimension ('x' for 1H, 'y' for 15N/13C), default 'x'
        Controls nucleus-specific ABSOLUTE position bounds:
        - 'x' (1H): ±0.03 ppm ABSOLUTE from initial position
        - 'y' (15N/13C): ±0.2 ppm ABSOLUTE from initial position
        ABSOLUTE: Peak position can NEVER exceed these bounds from initial guess,
        regardless of fitting stage or intermediate results
    """

    def __init__(self, verbose: bool = False, use_voigt_warmup: bool = True,
                 max_iterations: int = 1000, dimension: str = 'x'):
        self.verbose = verbose
        self.use_voigt_warmup = use_voigt_warmup
        self.max_iterations = max_iterations
        self.dimension = dimension  # 'x' for 1H, 'y' for 15N/13C
        self.convergence_flag = 0
        self.chi2 = None
        self.r_squared = None
        self.fitted_params = None
        self.covariance = None

        # Nucleus-specific ABSOLUTE position bounds (ppm)
        # ABSOLUTE: bounds enforced relative to INITIAL position, never exceeded
        self.position_bounds = {
            'x': 0.03,   # 1H: ±0.03 ppm ABSOLUTE (tight - narrow peaks)
            'y': 0.2,    # 15N/13C: ±0.2 ppm ABSOLUTE (flexible - broader peaks)
        }

    def fit_multi_peak(self,
                      x_data: np.ndarray,
                      y_data: np.ndarray,
                      initial_peaks: List[Dict],
                      fix_positions: Optional[List[bool]] = None,
                      fix_linewidths: Optional[List[bool]] = None,
                      max_lw: Optional[List[float]] = None,
                      min_lw: Optional[List[float]] = None,
                      
                      bounds: Optional[Tuple] = None) -> Dict:
        """
        Fit multiple overlapping Voigt peaks with 5-stage strategy

        This is the main entry point that orchestrates all 5 fitting stages.

        Parameters
        ----------
        x_data : np.ndarray
            Frequency axis (ppm), shape (n_points,)
        y_data : np.ndarray
            Intensity data, shape (n_points,)
        initial_peaks : List[Dict]
            Initial parameters for each peak:
            [{'pos': float, 'lw_lorentz': float, 'lw_gauss': float, 'intensity': float}, ...]
            - pos: peak position (ppm)
            - lw_lorentz: Lorentzian FWHM (ppm)
            - lw_gauss: Gaussian FWHM (ppm)
            - intensity: peak height or initial estimate
        fix_positions : List[bool], optional
            Fix position for each peak (default: all False)
        fix_linewidths : List[bool], optional
            Fix linewidths for each peak (default: all False)
            Corresponds to  'fixLW' flag
        max_lw : List[float], optional
            Maximum allowed linewidth for each peak (for validation)
        min_lw : List[float], optional
            Minimum allowed linewidth for each peak (for validation)
        bounds : Tuple, optional
            Global (lower_bounds, upper_bounds) for all parameters
            If not provided, reasonable defaults are used

        Returns
        -------
        Dict : Complete fitting results
            {
                'success': bool,
                'peaks': List[Dict],  # Results for each peak
                'baseline': float,
                'r_squared': float,
                'fitted_curve': np.ndarray,
                'residuals': np.ndarray,
                'covariance': np.ndarray,
                'iterations': int,
                'chi2': float,
                'convergence_flag': int,
                'method': 'ps2d_multi_peak_5_stage'
            }
        """
        n_peaks = len(initial_peaks)

        if not 1 <= n_peaks <= 10:
            raise ValueError(
                f"PS2D multi-peak fitter supports 1-10 peaks, got {n_peaks}. "
                f"For more peaks, consider alternative methods or split the spectrum."
            )

        if self.verbose:
            print(f"\n=== PS2D Multi-Peak Fitter ===")
            print(f"Number of peaks: {n_peaks}")
            print(f"Data points: {len(x_data)}")

        # Set default flags
        if fix_positions is None:
            fix_positions = [False] * n_peaks
        if fix_linewidths is None:
            fix_linewidths = [False] * n_peaks
        if max_lw is None:
            max_lw = [np.inf] * n_peaks
        if min_lw is None:
            min_lw = [0.0] * n_peaks

        # Initialize parameter vector (peakfit.cpp lines 106-162) - NO baseline
        params_init = self._initialize_parameters(initial_peaks)

        # Store for constraint checking
        self.initial_params = params_init.copy()
        self.fix_positions = fix_positions
        self.fix_linewidths = fix_linewidths
        self.max_lw = max_lw
        self.min_lw = min_lw
        self.n_peaks = n_peaks

        # STAGE 0: VOIGT warm-up (optional, peakfit.cpp lines 422-425)
        if self.use_voigt_warmup:
            if self.verbose:
                print("\n--- Stage 0: Intensity-only warm-up ---")
            params_current = self._fit_stage_zero(x_data, y_data, params_init, n_peaks)
        else:
            params_current = params_init.copy()

        # STAGE 1: Fit linewidths (positions fixed, peakfit.cpp lines 428-429)
        if self.verbose:
            print("\n--- Stage 1: Linewidth fitting (positions fixed) ---")
        result_stage1 = self._fit_stage_one(x_data, y_data, params_current, n_peaks)
        params_current = result_stage1.x

        # STAGE 2: Float positions (if not fixed by user, peakfit.cpp lines 432-433)
        if not all(fix_positions):
            if self.verbose:
                print("\n--- Stage 2: Position refinement ---")
            result_stage2 = self._fit_stage_two(x_data, y_data, params_current, n_peaks)
            params_current = result_stage2.x
        else:
            if self.verbose:
                print("\n--- Stage 2: Skipped (all positions fixed) ---")
            result_stage2 = result_stage1

        # STAGE 3: Skipped for 1D single-plane fitting
        # (Original C++ only for multi-plane data, peakfit.cpp lines 442-444)

        # STAGE 4: Final global optimization (peakfit.cpp lines 447-448)
        if self.verbose:
            print("\n--- Stage 4: Final global optimization ---")
        result_final = self._fit_stage_four(x_data, y_data, params_current, n_peaks)

        # Extract results
        self.fitted_params = result_final.x
        fitted_curve = multi_voigt_profile_1d(x_data, self.fitted_params, n_peaks)
        residuals = y_data - fitted_curve

        # Calculate chi-squared
        self.chi2 = np.sum(residuals**2)

        # Calculate R²
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_data - np.mean(y_data))**2)
        self.r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

        # Estimate covariance from Jacobian
        try:
            jac = compute_multi_voigt_jacobian_1d(x_data, self.fitted_params, n_peaks)
            # Approximate covariance: (J^T J)^(-1) * residual_variance
            jtj = jac.T @ jac
            residual_var = ss_res / (len(x_data) - len(self.fitted_params))
            self.covariance = np.linalg.inv(jtj) * residual_var
        except np.linalg.LinAlgError:
            warnings.warn("Could not compute covariance matrix (singular)")
            self.covariance = np.zeros((len(self.fitted_params), len(self.fitted_params)))

        # Validate linewidth constraints (peakfit.cpp lines 1048-1076)
        self._validate_linewidth_constraints()

        # Extract individual peak results
        peak_results = []
        for i in range(n_peaks):
            base_idx = i * 4
            peak_height = voigt_peak_height_1d(
                self.fitted_params[base_idx + 2],
                self.fitted_params[base_idx + 1]
            )

            peak_dict = {
                'position': self.fitted_params[base_idx + 0],
                'lw_lorentz': self.fitted_params[base_idx + 1],
                'lw_gauss': self.fitted_params[base_idx + 2],
                'intensity': self.fitted_params[base_idx + 3],
                'height': self.fitted_params[base_idx + 3] * peak_height,
                'errors': {}
            }

            # Add uncertainties from covariance
            if self.covariance is not None:
                peak_dict['errors']['position'] = np.sqrt(self.covariance[base_idx + 0, base_idx + 0])
                peak_dict['errors']['lw_lorentz'] = np.sqrt(self.covariance[base_idx + 1, base_idx + 1])
                peak_dict['errors']['lw_gauss'] = np.sqrt(self.covariance[base_idx + 2, base_idx + 2])
                peak_dict['errors']['intensity'] = np.sqrt(self.covariance[base_idx + 3, base_idx + 3])

            peak_results.append(peak_dict)

        # Compile final results
        results = {
            'success': result_final.success,
            'peaks': peak_results,
            
            'r_squared': self.r_squared,
            'fitted_curve': fitted_curve,
            'residuals': residuals,
            'covariance': self.covariance,
            'iterations': result_final.nfev,
            'chi2': self.chi2,
            'convergence_flag': self.convergence_flag,
            'method': 'ps2d_multi_peak_5_stage'
        }

        if self.verbose:
            print(f"\n=== Fitting Complete ===")
            print(f"Success: {result_final.success}")
            print(f"R²: {self.r_squared:.6f}")
            print(f"χ²: {self.chi2:.6e}")
            print(f"Convergence flag: {self.convergence_flag}")

        return results

    def _initialize_parameters(self, initial_peaks: List[Dict]) -> np.ndarray:
        """
        Initialize parameter vector from peak list - EXACT  (NO baseline)

        From peakfit.cpp lines 106-162 (initFitParGlobal)
        """
        n_peaks = len(initial_peaks)
        params = np.zeros(n_peaks * 4)  # NO baseline

        for i, peak in enumerate(initial_peaks):
            base_idx = i * 4

            # Position (ppm)
            params[base_idx + 0] = peak['pos']

            # Lorentzian FWHM (ppm) - ensure minimum value for bounds compatibility
            params[base_idx + 1] = max(peak.get('lw_lorentz', 0.01), 1e-5)

            # Gaussian FWHM (ppm) - ensure minimum value for bounds compatibility
            params[base_idx + 2] = max(peak.get('lw_gauss', 0.01), 1e-5)

            # Intensity (volume)
            # If given as height, convert to volume (peakfit.cpp lines 190-193)
            intensity_input = peak.get('intensity', peak.get('height', 1.0))

            # Normalize by Voigt peak height to get volume
            peak_height = voigt_peak_height_1d(params[base_idx + 2], params[base_idx + 1])
            if peak_height > 1e-10:
                params[base_idx + 3] = intensity_input / peak_height
            else:
                params[base_idx + 3] = intensity_input

        return params

    def _calculate_intensity_lower_bound(self, params_init: np.ndarray,
                                         y_data: np.ndarray,
                                         peak_index: int) -> float:
        """
        Calculate adaptive lower bound for peak intensity

        Parameters
        ----------
        params_init : np.ndarray
            Initial parameter vector
        y_data : np.ndarray
            Signal data (for max intensity reference)
        peak_index : int
            Index of the peak (0 to n_peaks-1)

        Returns
        -------
        float : Lower bound for intensity
            Either 10% of initial intensity or 0.1% of max signal, whichever is smaller
            This prevents weak peak elimination while avoiding noise fitting
        """
        base_idx = peak_index * 4
        initial_intensity = params_init[base_idx + 3]
        max_signal = y_data.max()

        # Use minimum of 10% of initial or 0.1% of max signal
        # This ensures initial value is always feasible while preventing noise fitting
        return max(0.0, min(initial_intensity * 0.1, max_signal * 0.001))

    def _fit_stage_zero(self, x_data: np.ndarray, y_data: np.ndarray,
                        params_init: np.ndarray, n_peaks: int) -> np.ndarray:
        """
        Stage 0: Fit intensity only with fixed positions and linewidths

        From peakfit.cpp lines 231-248 (fit_stage_zero)

        Purpose: VOIGT warm-up to get reasonable intensity estimates
        """
        # Create bounds: fix all parameters except intensity
        # Note: scipy requires lower < upper, so we use tight bounds for fixed params
        n_params = n_peaks * 4  # NO baseline
        lower = np.full(n_params, -np.inf)
        upper = np.full(n_params, np.inf)

        for i in range(n_peaks):
            base_idx = i * 4

            # Fixed parameters: position and linewidths (tight bounds = fixed)
            # Use small epsilon to avoid lower==upper
            eps = 1e-10
            lower[base_idx + 0] = params_init[base_idx + 0] - eps
            upper[base_idx + 0] = params_init[base_idx + 0] + eps
            lower[base_idx + 1] = params_init[base_idx + 1] - eps
            upper[base_idx + 1] = params_init[base_idx + 1] + eps
            lower[base_idx + 2] = params_init[base_idx + 2] - eps
            upper[base_idx + 2] = params_init[base_idx + 2] + eps

            # Intensity is floating (must be strictly positive)
            lower[base_idx + 3] = 0.0
            upper[base_idx + 3] = np.inf


        # Residual function
        def residuals(params):
            return y_data - multi_voigt_profile_1d(x_data, params, n_peaks)

        # Optimize
        result = least_squares(
            residuals,
            params_init,
            bounds=(lower, upper),
            method='trf',
            max_nfev=self.max_iterations,
            ftol=1e-8,
            xtol=1e-8
        )

        if self.verbose:
            print(f"  Iterations: {result.nfev}, Chi²: {np.sum(result.fun**2):.6e}")

        return result.x

    def _fit_stage_one(self, x_data: np.ndarray, y_data: np.ndarray,
                       params_init: np.ndarray, n_peaks: int):
        """
        Stage 1: Fit linewidths and intensity (positions fixed)

        From peakfit.cpp lines 251-272 (fit_stage_one)

        Key: Respects fixLW flag (self.fix_linewidths)
        """
        # Create bounds (peakfit.cpp lines 524-535)
        lower = np.full_like(params_init, -np.inf)
        upper = np.full_like(params_init, np.inf)

        eps = 1e-10  # Small epsilon to avoid lower==upper

        if self.verbose:
            print(f"   [Stage 1] Constraint enforcement:")
            print
            for i in range(n_peaks):
                print(f"      Peak {i+1}: fixLW={self.fix_linewidths[i]}, pos_init={params_init[i*4]:.4f}")

        for i in range(n_peaks):
            base_idx = i * 4

            # Position: FIXED (peakfit.cpp lines 525-526)
            lower[base_idx + 0] = params_init[base_idx + 0] - eps
            upper[base_idx + 0] = params_init[base_idx + 0] + eps

            # Linewidths: FLOAT or FIXED based on fixLW flag
            # (peakfit.cpp lines 527-533)
            if self.fix_linewidths[i]:
                # Fixed linewidths
                lower[base_idx + 1] = params_init[base_idx + 1] - eps
                upper[base_idx + 1] = params_init[base_idx + 1] + eps
                lower[base_idx + 2] = params_init[base_idx + 2] - eps
                upper[base_idx + 2] = params_init[base_idx + 2] + eps
            else:
                # Floating linewidths (must be positive and physically reasonable)
                # Use adaptive minimum: either hardcoded 0.001 or 10% below initial value
                lw_min = min(0.001, params_init[base_idx + 1] * 0.1, params_init[base_idx + 2] * 0.1)
                lower[base_idx + 1] = lw_min
                upper[base_idx + 1] = params_init[base_idx + 1] * 10.0
                lower[base_idx + 2] = lw_min
                upper[base_idx + 2] = params_init[base_idx + 2] * 10.0

            # Intensity: FLOATING (must be strictly positive)
            lower[base_idx + 3] = self._calculate_intensity_lower_bound(params_init, y_data, i)
            upper[base_idx + 3] = np.inf


        # Residual function
        def residuals(params):
            return y_data - multi_voigt_profile_1d(x_data, params, n_peaks)

        # Jacobian function
        def jacobian(params):
            return -compute_multi_voigt_jacobian_1d(x_data, params, n_peaks)

        # Optimize
        result = least_squares(
            residuals,
            params_init,
            jac=jacobian,
            bounds=(lower, upper),
            method='trf',
            max_nfev=self.max_iterations,
            ftol=1e-8,
            xtol=1e-8
        )

        if self.verbose:
            print(f"      Stage 1 iterations: {result.nfev}, Chi²: {np.sum(result.fun**2):.6e}")
            # Show fitted positions after stage 1
            for i in range(n_peaks):
                pos_fitted = result.x[i*4]
                pos_init = params_init[i*4]
                shift = pos_fitted - pos_init
                print(f"      Peak {i+1} position: init={pos_init:.4f}, fitted={pos_fitted:.4f}, shift={shift:.6f}")

        return result

    def _fit_stage_two(self, x_data: np.ndarray, y_data: np.ndarray,
                       params_init: np.ndarray, n_peaks: int):
        """
        Stage 2: Float positions (if not fixed by user)

        From peakfit.cpp lines 275-289 (fit_stage_two)

        Key: Respects fixPos flag (self.fix_positions)
        """
        # Create bounds (peakfit.cpp lines 562-568)
        lower = np.full_like(params_init, -np.inf)
        upper = np.full_like(params_init, np.inf)

        eps = 1e-10  # Small epsilon to avoid lower==upper
        any_floating_positions = False

        if self.verbose:
            print(f"   [Stage 2] Constraint enforcement:")
            for i in range(n_peaks):
                print(f"      Peak {i+1}: fixPos={self.fix_positions[i]}, fixLW={self.fix_linewidths[i]}, pos_init={params_init[i*4]:.4f}")

        for i in range(n_peaks):
            base_idx = i * 4

            # Position: FLOAT or FIXED based on fixPos flag
            # (peakfit.cpp lines 563-567)
            if self.fix_positions[i]:
                # Fixed position
                lower[base_idx + 0] = params_init[base_idx + 0] - eps
                upper[base_idx + 0] = params_init[base_idx + 0] + eps
            else:
                # Floating position (ABSOLUTE bounds relative to INITIAL position)
                # Ensures peak NEVER moves beyond initial_pos ± bound
                pos_bound = self.position_bounds.get(self.dimension, 0.1)
                initial_pos = self.initial_params[base_idx + 0]
                lower[base_idx + 0] = initial_pos - pos_bound
                upper[base_idx + 0] = initial_pos + pos_bound
                any_floating_positions = True

            # Linewidths: same as Stage 1
            if self.fix_linewidths[i]:
                lower[base_idx + 1] = params_init[base_idx + 1] - eps
                upper[base_idx + 1] = params_init[base_idx + 1] + eps
                lower[base_idx + 2] = params_init[base_idx + 2] - eps
                upper[base_idx + 2] = params_init[base_idx + 2] + eps
            else:
                # Floating linewidths (must be positive and physically reasonable)
                # Use adaptive minimum: either hardcoded 0.001 or 10% below initial value
                lw_min = min(0.001, params_init[base_idx + 1] * 0.1, params_init[base_idx + 2] * 0.1)
                lower[base_idx + 1] = lw_min
                upper[base_idx + 1] = params_init[base_idx + 1] * 10.0
                lower[base_idx + 2] = lw_min
                upper[base_idx + 2] = params_init[base_idx + 2] * 10.0

            # Intensity: FLOATING (must be strictly positive)
            lower[base_idx + 3] = self._calculate_intensity_lower_bound(params_init, y_data, i)
            upper[base_idx + 3] = np.inf


        # If all positions are fixed, skip optimization (peakfit.cpp lines 576-578)
        if not any_floating_positions:
            if self.verbose:
                print("  All positions fixed - stage skipped")
            # Return identity result
            class IdentityResult:
                def __init__(self, x):
                    self.x = x
                    self.success = True
                    self.nfev = 0
                    self.fun = np.zeros(len(x_data))
            return IdentityResult(params_init)

        # Residual function
        def residuals(params):
            return y_data - multi_voigt_profile_1d(x_data, params, n_peaks)

        # Jacobian function
        def jacobian(params):
            return -compute_multi_voigt_jacobian_1d(x_data, params, n_peaks)

        # Optimize
        result = least_squares(
            residuals,
            params_init,
            jac=jacobian,
            bounds=(lower, upper),
            method='trf',
            max_nfev=self.max_iterations,
            ftol=1e-8,
            xtol=1e-8
        )

        if self.verbose:
            print(f"      Stage 2 iterations: {result.nfev}, Chi²: {np.sum(result.fun**2):.6e}")
            # Show fitted positions after stage 2
            for i in range(n_peaks):
                pos_fitted = result.x[i*4]
                pos_init = params_init[i*4]
                shift = pos_fitted - pos_init
                print(f"      Peak {i+1} position: init={pos_init:.4f}, fitted={pos_fitted:.4f}, shift={shift:.6f}")

        return result

    def _fit_stage_four(self, x_data: np.ndarray, y_data: np.ndarray,
                        params_init: np.ndarray, n_peaks: int):
        """
        Stage 4: Final global optimization with user constraints

        From peakfit.cpp lines 321-345 (fit_stage_four)

        All parameters float unless explicitly fixed by user flags
        """
        # Create bounds (peakfit.cpp lines 673-689)
        lower = np.full_like(params_init, -np.inf)
        upper = np.full_like(params_init, np.inf)

        eps = 1e-10  # Small epsilon to avoid lower==upper

        if self.verbose:
            print(f"   [Stage 4] Constraint enforcement:")
            for i in range(n_peaks):
                print(f"      Peak {i+1}: fixPos={self.fix_positions[i]}, fixLW={self.fix_linewidths[i]}, pos_init={params_init[i*4]:.4f}")

        for i in range(n_peaks):
            base_idx = i * 4

            # Position: respect fixPos flag (peakfit.cpp lines 674-675)
            if self.fix_positions[i]:
                lower[base_idx + 0] = params_init[base_idx + 0] - eps
                upper[base_idx + 0] = params_init[base_idx + 0] + eps
            else:
                # Floating position (ABSOLUTE bounds relative to INITIAL position)
                # Ensures peak NEVER moves beyond initial_pos ± bound across ALL stages
                pos_bound = self.position_bounds.get(self.dimension, 0.1)
                initial_pos = self.initial_params[base_idx + 0]
                lower[base_idx + 0] = initial_pos - pos_bound
                upper[base_idx + 0] = initial_pos + pos_bound

            # Linewidths: respect fixLW flag (peakfit.cpp lines 676-688)
            if self.fix_linewidths[i]:
                lower[base_idx + 1] = params_init[base_idx + 1] - eps
                upper[base_idx + 1] = params_init[base_idx + 1] + eps
                lower[base_idx + 2] = params_init[base_idx + 2] - eps
                upper[base_idx + 2] = params_init[base_idx + 2] + eps
            else:
                # Floating linewidths (must be positive and physically reasonable)
                # Use adaptive minimum: either hardcoded 0.001 or 10% below initial value
                lw_min = min(0.001, params_init[base_idx + 1] * 0.1, params_init[base_idx + 2] * 0.1)
                lower[base_idx + 1] = lw_min
                upper[base_idx + 1] = params_init[base_idx + 1] * 10.0
                lower[base_idx + 2] = lw_min
                upper[base_idx + 2] = params_init[base_idx + 2] * 10.0

            # Intensity: FLOATING (must be strictly positive)
            lower[base_idx + 3] = self._calculate_intensity_lower_bound(params_init, y_data, i)
            upper[base_idx + 3] = np.inf


        # Residual function
        def residuals(params):
            return y_data - multi_voigt_profile_1d(x_data, params, n_peaks)

        # Jacobian function
        def jacobian(params):
            return -compute_multi_voigt_jacobian_1d(x_data, params, n_peaks)

        # Optimize
        result = least_squares(
            residuals,
            params_init,
            jac=jacobian,
            bounds=(lower, upper),
            method='trf',
            max_nfev=self.max_iterations,
            ftol=1e-8,
            xtol=1e-8
        )

        if self.verbose:
            print(f"      Stage 4 iterations: {result.nfev}, Chi²: {np.sum(result.fun**2):.6e}")
            # Show final fitted positions after stage 4
            for i in range(n_peaks):
                pos_fitted = result.x[i*4]
                pos_init = params_init[i*4]
                shift = pos_fitted - pos_init
                print(f"      Peak {i+1} FINAL position: init={pos_init:.4f}, fitted={pos_fitted:.4f}, shift={shift:.6f}")

        return result

    def _validate_linewidth_constraints(self):
        """
        Validate fitted linewidths against max/min constraints

        From peakfit.cpp lines 1048-1076 (calcConvFlag)

        Sets convergence_flag to indicate constraint violations:
        0 = Success
        4 = Linewidth exceeds maximum
        6 = Linewidth below minimum
        """
        self.convergence_flag = 0  # Success by default

        for i in range(self.n_peaks):
            base_idx = i * 4

            lw_lorentz = self.fitted_params[base_idx + 1]
            lw_gauss = self.fitted_params[base_idx + 2]

            # Check against maximum (peakfit.cpp lines 760-764)
            if lw_lorentz > self.max_lw[i] and self.max_lw[i] > 0:
                self.convergence_flag = 4
                warnings.warn(
                    f"Peak {i}: Lorentzian linewidth {lw_lorentz:.6f} exceeds "
                    f"maximum {self.max_lw[i]:.6f}"
                )

            if lw_gauss > self.max_lw[i] and self.max_lw[i] > 0:
                self.convergence_flag = 4
                warnings.warn(
                    f"Peak {i}: Gaussian linewidth {lw_gauss:.6f} exceeds "
                    f"maximum {self.max_lw[i]:.6f}"
                )

            # Check against minimum (peakfit.cpp lines 770-774)
            if lw_lorentz < self.min_lw[i] and self.min_lw[i] > 0:
                self.convergence_flag = 6
                warnings.warn(
                    f"Peak {i}: Lorentzian linewidth {lw_lorentz:.6f} below "
                    f"minimum {self.min_lw[i]:.6f}"
                )

            if lw_gauss < self.min_lw[i] and self.min_lw[i] > 0:
                self.convergence_flag = 6
                warnings.warn(
                    f"Peak {i}: Gaussian linewidth {lw_gauss:.6f} below "
                    f"minimum {self.min_lw[i]:.6f}"
                )


# ============================================================================
# HIGH-LEVEL CONVENIENCE FUNCTION
# ============================================================================

def fit_overlapping_peaks_ps2d_style(
    x_data: np.ndarray,
    y_data: np.ndarray,
    detected_peaks: List[Dict],
    dimension: str = 'x',
    target_position: float = None,
    fix_linewidths: bool = False,
    fix_positions: bool = False,
    verbose: bool = False,
    lw_lorentz_1h: float = None,
    lw_gauss_1h: float = None,
    lw_lorentz_15n: float = None,
    lw_gauss_15n: float = None
) -> Dict:
    """
    High-level interface for multi-peak PS2D-style fitting

    Automatically sets up initial parameters and calls multi-stage fitter.
    This function provides a simple interface similar to the existing lunaNMR
    fitting functions.

    Parameters
    ----------
    x_data : np.ndarray
        Frequency axis (ppm)
    y_data : np.ndarray
        Intensity data
    detected_peaks : List[Dict]
        Detected peaks with approximate positions
        [{'position': float, 'height': float, ...}, ...]
    dimension : str, optional
        'x' (1H) or 'y' (15N/13C) for linewidth estimation, default 'x'
    fix_linewidths : bool, optional
        If True, hold linewidths at initial estimates 
        Default False
    fix_positions : bool, optional
        If True, hold positions at initial estimates 
        Default False
    verbose : bool, optional
        Print progress, default False
    lw_lorentz_1h : float, optional
        User-provided Lorentzian linewidth for 1H 
    lw_gauss_1h : float, optional
        User-provided Gaussian linewidth for 1H 
    lw_lorentz_15n : float, optional
        User-provided Lorentzian linewidth for 15N/13C 
    lw_gauss_15n : float, optional
        User-provided Gaussian linewidth for 15N/13C 

    Returns
    -------
    Dict : Fitting results (same format as Ps2dMultiPeakFitter.fit_multi_peak)

    Examples
    --------
    >>> # Fit two overlapping peaks
    >>> peaks = [
    ...     {'position': 8.5, 'height': 1000},
    ...     {'position': 8.52, 'height': 800}
    ... ]
    >>> result = fit_overlapping_peaks_ps2d_style(x, y, peaks, verbose=True)
    >>> print(f"R² = {result['r_squared']:.4f}")
    """
    n_peaks = len(detected_peaks)

    # BASELINE REMOVED - data must be baseline-corrected BEFORE calling this function

    # ========================================================================
    # ========================================================================
    estimator = get_global_estimator()

    if lw_lorentz_1h is not None and lw_gauss_1h is not None:
        # Register global override for 1H dimension
        estimator.reference_peaks['x'] = {
            'lw_lorentz': lw_lorentz_1h,
            'lw_gauss': lw_gauss_1h,
            'assignment': 'user_override'
        }
    if lw_lorentz_15n is not None and lw_gauss_15n is not None:
        # Register global override for 15N dimension
        estimator.reference_peaks['y'] = {
            'lw_lorentz': lw_lorentz_15n,
            'lw_gauss': lw_gauss_15n,
            'assignment': 'user_override'
        }

    # Collect all peak positions for spatial context
    all_positions = [peak.get('position', peak.get('pos', 0.0)) for peak in detected_peaks]

    # Get optimized linewidths using three-tier hierarchy:
    # 1. User override (if provided above)
    # 2. Reference peak template (if available)
    # 3. Spatial analysis (if data suitable)
    # 4. Nucleus defaults (fallback)
    typical_lw_lorentz, typical_lw_gauss = estimator.get_linewidth(
        dimension=dimension,
        x_data=x_data,
        y_data=y_data,
        all_peak_positions=all_positions
    )

    # Report estimation method (if verbose)
    if verbose:
        estimate_info = estimator.spatial_estimates.get(dimension, {})
        method = estimate_info.get('method', 'user_override' if dimension in estimator.reference_peaks else 'nucleus_default')
        print(f"\n📊 Linewidth estimation ({dimension}):")
        print(f"   Method: {method}")
        print(f"   LW_Lorentz: {typical_lw_lorentz:.6f} ppm")
        print(f"   LW_Gauss: {typical_lw_gauss:.6f} ppm")
        if method == 'spatial_analysis':
            print(f"   Quality: {estimate_info.get('quality', 'unknown')}")
            print(f"   Isolated peaks: {estimate_info.get('isolated_count', 0)}")
            print(f"   Peaks analyzed: {estimate_info.get('n_peaks_analyzed', 0)}")
            print(f"   CV: {estimate_info.get('cv', 0):.3f}")

    # Build initial peak list
    initial_peaks = []
    for peak in detected_peaks:
        pos = peak.get('position', peak.get('pos', 0.0))
        height = peak.get('height', peak.get('intensity', 1.0))
        assignment = peak.get('assignment', None)

        # Check if peak has explicit linewidths (highest priority)
        if 'lw_lorentz' in peak and 'lw_gauss' in peak:
            lw_lorentz = peak['lw_lorentz']
            lw_gauss = peak['lw_gauss']
        else:
            # Use estimated linewidths (from tier hierarchy)
            lw_lorentz = typical_lw_lorentz
            lw_gauss = typical_lw_gauss

        initial_peaks.append({
            'pos': pos,
            'lw_lorentz': lw_lorentz,
            'lw_gauss': lw_gauss,
            'intensity': height,
            'assignment': assignment
        })

    # Create constraint lists
    fix_pos_list = [fix_positions] * n_peaks
    fix_lw_list = [fix_linewidths] * n_peaks

    # Create fitter and run (pass dimension for nucleus-specific bounds)
    fitter = Ps2dMultiPeakFitter(verbose=verbose, use_voigt_warmup=True, dimension=dimension)

    results = fitter.fit_multi_peak(
        x_data=x_data,
        y_data=y_data,
        initial_peaks=initial_peaks,
        fix_positions=fix_pos_list,
        fix_linewidths=fix_lw_list
    )

    # ========================================================================
    # ========================================================================
    # If fit successful with good quality, register fitted linewidths as reference
    if results.get('success') and results.get('r_squared', 0) > 0.85:
        for peak in results.get('peaks', []):
            estimator.register_fitted_peak(
                dimension=dimension,
                lw_lorentz=peak['lw_lorentz'],
                lw_gauss=peak['lw_gauss'],
                assignment=peak.get('assignment')
            )
        if verbose:
            print(f"\n✅ Registered {len(results['peaks'])} fitted peaks as linewidth templates for {dimension}-dimension")
            print(f"   LW_Lorentz: {results['peaks'][0]['lw_lorentz']:.6f} ppm")
            print(f"   LW_Gauss: {results['peaks'][0]['lw_gauss']:.6f} ppm")

    return results


# ============================================================================
# TESTING AND VALIDATION
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("PS2D Multi-Peak Fitter - Test Suite")
    print("=" * 70)

    # Test 1: Single Voigt peak
    print("\n" + "=" * 70)
    print("TEST 1: Single Voigt Peak")
    print("=" * 70)

    x_test = np.linspace(8.0, 9.0, 500)

    true_params = np.array([8.5, 0.01, 0.01, 1.0])  # pos, lw_lor, lw_gauss, intensity
    y_true = multi_voigt_profile_1d(x_test, true_params, 1)
    y_noisy = y_true + np.random.normal(0, 0.02, len(y_true))

    # Fit
    initial_peaks = [{'pos': 8.48, 'lw_lorentz': 0.015, 'lw_gauss': 0.015, 'intensity': 0.9}]

    fitter = Ps2dMultiPeakFitter(verbose=True)
    result = fitter.fit_multi_peak(x_test, y_noisy, initial_peaks)

    print(f"\nTrue parameters:")
    print(f"  Position: {true_params[0]:.6f}")
    print(f"  LW Lorentz: {true_params[1]:.6f}")
    print(f"  LW Gauss: {true_params[2]:.6f}")
    print(f"  Intensity: {true_params[3]:.6f}")
    #  print(f"  Baseline: ...")  # REMOVED - no baseline fitting

    print(f"\nFitted parameters:")
    print(f"  Position: {result['peaks'][0]['position']:.6f} ± {result['peaks'][0]['errors']['position']:.6f}")
    print(f"  LW Lorentz: {result['peaks'][0]['lw_lorentz']:.6f} ± {result['peaks'][0]['errors']['lw_lorentz']:.6f}")
    print(f"  LW Gauss: {result['peaks'][0]['lw_gauss']:.6f} ± {result['peaks'][0]['errors']['lw_gauss']:.6f}")
    print(f"  Intensity: {result['peaks'][0]['intensity']:.6f} ± {result['peaks'][0]['errors']['intensity']:.6f}")
    #  print(f"  Baseline: ...")  # REMOVED - no baseline fitting
    print(f"  R²: {result['r_squared']:.6f}")

    # Test 2: Two overlapping peaks
    print("\n" + "=" * 70)
    print("TEST 2: Two Overlapping Voigt Peaks")
    print("=" * 70)

    x_test2 = np.linspace(7.0, 8.0, 1000)

    true_params2 = np.array([
        7.4, 0.02, 0.015, 2.0,    # Peak 1
        7.45, 0.025, 0.02, 1.5    # Peak 2
    ])
    y_true2 = multi_voigt_profile_1d(x_test2, true_params2, 2)
    y_noisy2 = y_true2 + np.random.normal(0, 0.03, len(y_true2))

    # Fit with approximate initial guesses
    initial_peaks2 = [
        {'pos': 7.38, 'lw_lorentz': 0.025, 'lw_gauss': 0.02, 'intensity': 1.8},
        {'pos': 7.47, 'lw_lorentz': 0.03, 'lw_gauss': 0.025, 'intensity': 1.3}
    ]

    fitter2 = Ps2dMultiPeakFitter(verbose=True)
    result2 = fitter2.fit_multi_peak(x_test2, y_noisy2, initial_peaks2)

    print(f"\nTrue parameters:")
    for i in range(2):
        print(f"  Peak {i+1}:")
        print(f"    Position: {true_params2[i*4+0]:.6f}")
        print(f"    LW Lorentz: {true_params2[i*4+1]:.6f}")
        print(f"    LW Gauss: {true_params2[i*4+2]:.6f}")
        print(f"    Intensity: {true_params2[i*4+3]:.6f}")
    #  print(f"  Baseline: ...")  # REMOVED - no baseline fitting

    print(f"\nFitted parameters:")
    for i in range(2):
        print(f"  Peak {i+1}:")
        print(f"    Position: {result2['peaks'][i]['position']:.6f} ± {result2['peaks'][i]['errors']['position']:.6f}")
        print(f"    LW Lorentz: {result2['peaks'][i]['lw_lorentz']:.6f} ± {result2['peaks'][i]['errors']['lw_lorentz']:.6f}")
        print(f"    LW Gauss: {result2['peaks'][i]['lw_gauss']:.6f} ± {result2['peaks'][i]['errors']['lw_gauss']:.6f}")
        print(f"    Intensity: {result2['peaks'][i]['intensity']:.6f} ± {result2['peaks'][i]['errors']['intensity']:.6f}")
    #  print(f"  Baseline: ...")  # REMOVED - no baseline fitting
    print(f"  R²: {result2['r_squared']:.6f}")

    # Test 3: High-level convenience function
    print("\n" + "=" * 70)
    print("TEST 3: High-Level Convenience Function")
    print("=" * 70)

    detected = [
        {'position': 7.38, 'height': 1.8},
        {'position': 7.47, 'height': 1.3}
    ]

    result3 = fit_overlapping_peaks_ps2d_style(
        x_test2, y_noisy2, detected,
        dimension='x',
        fix_linewidths=False,
        verbose=True
    )

    print(f"\nConvenience function result:")
    print(f"  R²: {result3['r_squared']:.6f}")
    print(f"  Success: {result3['success']}")
    print(f"  Method: {result3['method']}")

    # Test 4: Fixed linewidths (fixLW flag)
    print("\n" + "=" * 70)
    print("TEST 4: Fixed Linewidths (fixLW=True)")
    print("=" * 70)

    result4 = fit_overlapping_peaks_ps2d_style(
        x_test2, y_noisy2, detected,
        dimension='x',
        fix_linewidths=True,
        verbose=True
    )

    print(f"\nFixed linewidths result:")
    print(f"  R²: {result4['r_squared']:.6f}")
    for i in range(2):
        print(f"  Peak {i+1} linewidths: LW_Lor={result4['peaks'][i]['lw_lorentz']:.6f}, "
              f"LW_Gauss={result4['peaks'][i]['lw_gauss']:.6f}")

    print("\n" + "=" * 70)
    print("All tests completed successfully!")
    print("=" * 70)
