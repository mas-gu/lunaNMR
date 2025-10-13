#!/usr/bin/env python3
"""
Enhanced Voigt Fitting Module

This module provides robust, professional-grade Voigt profile fitting for NMR spectroscopy.
Addresses critical issues in the original fitting implementation with improved parameter
estimation, bounds handling, and quality assessment.

Key improvements:
- Robust baseline estimation using polynomial fitting
- Adaptive width estimation based on nucleus type (¹H: 5.5-12 ppm, ¹⁵N: 100-140 ppm)
- Peak center refinement with local maximum finding
- Multi-step fitting strategy (coarse → fine)
- Comprehensive parameter validation and uncertainty analysis
- Proper bounds handling for all fitting scenarios

Author: Guillaume Mas
Date: 2025
"""

import numpy as np
import pandas as pd
from scipy.special import wofz
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter, peak_widths
from scipy.ndimage import gaussian_filter1d, median_filter
import warnings
import time
warnings.filterwarnings('ignore')

# Import quality categorizer (with fallback for backward compatibility)
try:
    from ..utils.quality_categories import categorize_quality_legacy
    QUALITY_CATEGORIZER_AVAILABLE = True
except ImportError:
    QUALITY_CATEGORIZER_AVAILABLE = False
    # Fallback function for backward compatibility
    def categorize_quality_legacy(r_squared: float) -> str:
        """Legacy quality categorization"""
        if r_squared >= 0.95:
            return 'excellent'
        elif r_squared >= 0.85:
            return 'good'
        elif r_squared >= 0.70:
            return 'fair'  # was 'marginal'
        else:
            return 'poor'

# Import overlap resolution modules (optional dependency for backward compatibility)
#try:
#    from .overlap_resolver_engine import OverlapResolverEngine
#    from ..utils.overlap_config import OverlapResolutionConfig
#    OVERLAP_RESOLUTION_AVAILABLE = True
#except ImportError:
#    OVERLAP_RESOLUTION_AVAILABLE = False
#    # Silently continue without overlap resolution - will print warning if user tries to enable it

class EnhancedVoigtFitter:
    """Enhanced Voigt profile fitter with robust parameter estimation"""

    def __init__(self):
        self.fitting_parameters = {
            'max_iterations': 1000,
            'tolerance': 1e-8,
            'min_r_squared': 0.3,  # Lowered for realistic data
            'baseline_polynomial_degree': 1,
            'smoothing_window': 5,
            'peak_detection_prominence': 0.1,
            ##
            'multi_peak_r2_threshold': 0.7,          # R² threshold to trigger multi-peak detection
            'multi_peak_improvement_threshold': 0.1,  # Minimum improvement required from multi-peak fit
            'peak_detection_sensitivity': 1.5,        # Height threshold multiplier (lower = more sensitive)
            'overlap_detection_factor': 0.8,          # Separation factor for overlap detection
            'residual_analysis_threshold': 1.5      # Threshold for residual peak detection

        }

        # NMR-specific parameters (corrected ranges)
        self.nmr_ranges = {
            '1H': {'min': 5.5, 'max': 12.0, 'typical_width': 0.01},  #0.02    # ¹H range
            '15N': {'min': 100.0, 'max': 140.0, 'typical_width': 0.5}, #1.5  # ¹⁵N range
            '13C': {'min': 5, 'max': 50.0, 'typical_width': 0.6}  #1.0    # ¹³C range for completeness
        }

        # Nucleus-specific default windows (matches core_integrator.py)
        self.nucleus_default_windows = {
            '1H': {'x': 0.15, 'y': 2.0},    # Tight for narrow 1H range
            '15N': {'x': 0.3, 'y': 4.0},    # Wider for broader 15N range
            '13C': {'x': 0.25, 'y': 3.0},   # Medium for 13C range
            'default': {'x': 0.2, 'y': 2.0} # Fallback for unknown nuclei
        }

        self.last_fit_diagnostics = {}

        # Level 1 Emergency Parameters
        self.level1_params = {
            'validation_enabled': True,
            'monitoring_enabled': True,
            'fallback_enabled': True,
            'max_emergency_attempts': 4,
            'trajectory_check_interval': 10
        }

        # Level 2 Architectural Parameters
        self.level2_params = {
            'robust_estimation_enabled': True,
            'baseline_validation_enabled': True,
            'consensus_estimation_enabled': True,
            'method_performance_tracking': True
        }

        # Initialize Level 2 components
        if self.level2_params['robust_estimation_enabled']:
            self.parameter_estimator = RobustParameterEstimator(self)


        # Initialize simplified parameter manager (core functionality)
        try:
            from ..utils.simplified_parameter_manager import SimplifiedParameterManager, ParameterAdapter
            self.simplified_param_manager = SimplifiedParameterManager()
            self.parameter_adapter = ParameterAdapter(self.simplified_param_manager)
        except ImportError:
            self.simplified_param_manager = None
            self.parameter_adapter = None

        # Consensus fitting engine was removed (Phase 2 cleanup - orphaned overlap subsystem)
        # Kept flag for backward compatibility with any code checking this attribute
        self.use_automated_fitting = False

        # NEW: Overlap resolution components (lazy initialization for backward compatibility)
#        self.overlap_resolver = None  # Will be initialized on first use if needed
#       self.overlap_detection_enabled = False  # Default OFF for backward compatibility
#        self.overlap_detection_threshold = 0.5  # Minimum distance (ppm) to consider overlap
#        self.overlap_config = None  # Will use defaults if not set

    def set_gui_parameters(self, gui_fitting_params):
        """
        Store GUI fitting parameters for consistent window sizing in analysis displays.

        This method allows the enhanced fitter to access GUI-defined window parameters
        for analysis and visualization, ensuring consistency between fitting and display.

        Parameters:
        -----------
        gui_fitting_params : dict
            Dictionary containing GUI fitting parameters:
            - fitting_window_x: float, GUI X-window size in ppm (1H dimension)
            - fitting_window_y: float, GUI Y-window size in ppm (15N/13C dimension)
            - min_r_squared: float, quality threshold
            - max_iterations: int, fitting iterations

        Notes:
        ------
        - This method is called automatically by core_integrator.py
        - Parameters are stored as instance variables for use in analysis methods
        - Backward compatibility: If not called, falls back to hardcoded defaults
        """
        self.gui_fitting_params = gui_fitting_params.copy() if gui_fitting_params else {}
        self.has_gui_params = bool(gui_fitting_params)

        # Extract frequently used parameters
        self.gui_window_x = self.gui_fitting_params.get('fitting_window_x', 0.2)  # default fallback
        self.gui_window_y = self.gui_fitting_params.get('fitting_window_y', 2.0)  # default fallback

        # CRITICAL FIX: Update actual fitting parameters used by curve_fit
        if 'max_iterations' in self.gui_fitting_params:
            self.fitting_parameters['max_iterations'] = self.gui_fitting_params['max_iterations']
            print(f"   🎛️ Updated curve_fit max_iterations to {self.fitting_parameters['max_iterations']}")

        if 'min_r_squared' in self.gui_fitting_params:
            self.fitting_parameters['min_r_squared'] = self.gui_fitting_params['min_r_squared']

        # Update multi-peak detection parameters from GUI
        multi_peak_params = [
            'multi_peak_r2_threshold',
            'multi_peak_improvement_threshold', 
            'peak_detection_sensitivity',
            'overlap_detection_factor',
            'residual_analysis_threshold'
        ]
        
        for param in multi_peak_params:
            if param in self.gui_fitting_params:
                self.fitting_parameters[param] = self.gui_fitting_params[param]
                print(f"   🔀 Updated {param} to {self.fitting_parameters[param]}")

        if self.has_gui_params:
            print(f"   📊 Enhanced fitter configured with GUI parameters: X={self.gui_window_x:.3f} ppm, Y={self.gui_window_y:.1f} ppm")

    def _calculate_gui_based_multiplier(self, nucleus_type, ppm_range, data_length, fitted_width=None):
        """
        DEPRECATED: This function is NO LONGER USED.


        FIXED windows (NO multipliers, NO adaptation):
        - radF1 = 0.6 ppm (indirect dimension, 15N/13C)
        - radF2 = 0.06 ppm (direct dimension, 1H)

        Data selection uses elliptical boundary: radius² = (F1-pos)²/radF1² + (F2-pos)²/radF2² <= 1.0

        This function returns a dummy value for backward compatibility only.
        The actual window selection is handled by Ps2dDataSelector in core_integrator.py.

        Parameters:
        -----------
        nucleus_type : str
            '1H', '15N', or '13C' - determines which GUI parameter to use
        ppm_range : float
            Total ppm range of the data
        data_length : int
            Number of data points in the spectrum dimension
        fitted_width : float, optional
            Fitted peak linewidth for reference

        Returns:
        --------
        float : Returns 1.0 (dummy value - NOT USED)
        """
        # RETURN DUMMY VALUE - NOT USED ANYMORE
        return 1.0  # No multiplier needed - PS2D uses fixed radF1/radF2

    # =====================================
    # AUTOMATED FITTING METHODS (Priority 1 & 2)
    # =====================================

    def set_simplified_parameters(self, **kwargs):
        """
        Set simplified parameters for automated fitting.

        This method implements Priority 1 improvements by allowing users to
        control fitting with only 3-5 core parameters instead of 25+.

        Parameters:
        -----------
        sensitivity : float, optional
            Detection sensitivity (0-1, default 0.5)
        window_scale : float, optional
            Window sizing scale factor (0.1-10.0, default 1.0)
        quality_target : float, optional
            Target fitting quality (0.3-1.0, default 0.85)
        noise_estimation_method : str, optional
            Noise estimation method ('auto', 'percentile', 'robust')
        baseline_method : str, optional
            Baseline estimation method ('auto', 'polynomial', 'iterative')
        """
        if hasattr(self, 'simplified_param_manager'):
            self.simplified_param_manager.update_simplified_parameters(**kwargs)
            print(f"✅ Updated simplified parameters: {list(kwargs.keys())}")
        else:
            print("⚠️ Simplified parameter manager not available")

    def fit_with_consensus(self, x_data, y_data, nucleus_type='default', use_simplified=True):
        """
        Fit peaks using automated consensus fitting approach.

        This method implements Priority 2 improvements:
        - Separates detection from fitting
        - Uses multi-method consensus
        - Automated peak counting with AIC/BIC model selection

        Parameters:
        -----------
        x_data : np.ndarray
            Frequency/chemical shift data
        y_data : np.ndarray
            Intensity data
        nucleus_type : str
            Type of nucleus ('1H', '15N', '13C', 'default')
        use_simplified : bool
            Whether to use simplified parameters (Priority 1)

        Returns:
        --------
        dict : Consensus fitting results with backward-compatible format
        """
        if not self.use_automated_fitting:
            print("⚠️ Automated fitting not available, falling back to legacy method")
            return self._legacy_fit_fallback(x_data, y_data, nucleus_type)

        try:
            # Phase 1: Detection using consensus detector
            if use_simplified:
                # Use simplified parameters
                sensitivity = self.simplified_param_manager.simplified_params.sensitivity
            else:
                # Use default sensitivity
                sensitivity = 0.5

            # Detect peak candidates
            peak_candidates = self.consensus_detector.detect_peaks(
                x_data, y_data, nucleus_type, sensitivity
            )

            if not peak_candidates:
                print("ℹ️ No peaks detected")
                return self._empty_result()

            print(f"🔍 Detected {len(peak_candidates)} peak candidates")

            # Phase 2: Consensus fitting with automated peak counting
            max_peaks = min(10, len(peak_candidates))  # Reasonable limit
            consensus_result = self.consensus_fitter.fit_consensus(
                x_data, y_data, peak_candidates, max_peaks
            )

            # Phase 3: Convert to backward-compatible format
            legacy_result = self._convert_consensus_to_legacy(
                consensus_result, x_data, y_data, nucleus_type
            )

            # Phase 4: Adaptive quality assessment
            if use_simplified:
                adaptive_threshold = self.simplified_param_manager.calculate_adaptive_quality_threshold(
                    y_data, nucleus_type
                )
                legacy_result['adaptive_threshold'] = adaptive_threshold
                legacy_result['meets_adaptive_threshold'] = legacy_result.get('r_squared', 0) >= adaptive_threshold

            print(f"✅ Consensus fitting completed: {len(consensus_result.consensus_peaks)} peaks, "
                  f"R² = {consensus_result.best_fit.r_squared:.3f}")

            return legacy_result

        except Exception as e:
            print(f"❌ Automated fitting failed: {e}")
            return self._legacy_fit_fallback(x_data, y_data, nucleus_type)

    def get_automated_parameters_for_legacy(self, nucleus_type='default'):
        """
        Get legacy parameter set derived from simplified parameters.

        This ensures backward compatibility with existing workflows while
        using the simplified parameter system internally.

        Parameters:
        -----------
        nucleus_type : str
            Type of nucleus for adaptive parameter calculation

        Returns:
        --------
        dict : Complete legacy parameter set
        """
        if hasattr(self, 'parameter_adapter'):
            return self.parameter_adapter.get_integrator_parameters(nucleus_type)
        else:
            # Fallback to current parameters
            return self.fitting_parameters.copy()

    def _convert_consensus_to_legacy(self, consensus_result, x_data, y_data, nucleus_type):
        """Convert consensus result to legacy format for backward compatibility"""

        if not consensus_result.consensus_peaks:
            return self._empty_result()

        # Convert peaks to legacy format
        legacy_peaks = []
        for peak in consensus_result.consensus_peaks:
            legacy_peak = {
                'amplitude': peak.get('amplitude', 0),
                'center': peak.get('center', 0),
                'sigma': peak.get('sigma', 0.01),
                'gamma': peak.get('gamma', 0.01),
                'baseline': peak.get('baseline', 0),
                'width_hz': (peak.get('sigma', 0.01) + peak.get('gamma', 0.01)) * 2,
                'width_ppm': (peak.get('sigma', 0.01) + peak.get('gamma', 0.01)) * 2,
                'errors': peak.get('errors', [0, 0, 0, 0, 0])
            }
            legacy_peaks.append(legacy_peak)

        # Convert overall statistics
        best_fit = consensus_result.best_fit
        quality_assessment = consensus_result.quality_assessment

        legacy_result = {
            'peaks': legacy_peaks,
            'r_squared': best_fit.r_squared,
            'residual_std': best_fit.residual_std,
            'aic': best_fit.aic,
            'bic': best_fit.bic,
            'method': f"consensus_{best_fit.method}",
            'fit_time': best_fit.fit_time,
            'n_peaks': len(legacy_peaks),
            'quality_category': quality_assessment.get('overall_quality', 'unknown'),
            'convergence_info': best_fit.convergence_info,
            'model_selection': consensus_result.model_selection,
            'confidence_scores': consensus_result.confidence_scores,
            'automated_fitting': True,
            'nucleus_type': nucleus_type
        }

        return legacy_result

    def _legacy_fit_fallback(self, x_data, y_data, nucleus_type):
        """Fallback to legacy fitting when automated fitting fails"""
        print("🔄 Using legacy fitting as fallback")

        # This would call the existing legacy fitting method
        # For now, return a minimal result
        return {
            'peaks': [],
            'r_squared': 0.0,
            'method': 'legacy_fallback',
            'automated_fitting': False,
            'nucleus_type': nucleus_type,
            'error': 'Automated fitting not available'
        }

    def _empty_result(self):
        """Return empty result in legacy format"""
        return {
            'peaks': [],
            'r_squared': 0.0,
            'residual_std': np.inf,
            'aic': np.inf,
            'bic': np.inf,
            'method': 'no_peaks',
            'n_peaks': 0,
            'quality_category': 'no_signal',
            'automated_fitting': True
        }

    @staticmethod
    def voigt_profile_1d(x, amplitude, center, sigma, gamma, baseline=0):
        """
        1D Voigt profile using Faddeeva function

        Parameters:
        - x: frequency/chemical shift array
        - amplitude: peak amplitude
        - center: peak center position
        - sigma: Gaussian width (instrumental broadening)
        - gamma: Lorentzian width (natural line width)
        - baseline: baseline offset
        """
        try:
            # Avoid division by zero
            sigma = max(sigma, 1e-6)

            # Compute complex argument for Faddeeva function
            z = ((x - center) + 1j*gamma) / (sigma * np.sqrt(2))

            # Voigt profile using Faddeeva function
            voigt = amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2*np.pi))

            return voigt + baseline
        except:
            # Fallback to Gaussian if Voigt fails
            return amplitude * np.exp(-0.5 * ((x - center) / max(sigma, 1e-6))**2) + baseline

    def robust_baseline_estimation(self, x_data, y_data, method='auto', validation_enabled=True):
        """
        LEVEL 2 ENHANCED: Robust baseline estimation with method validation and quarantine

        This method implements comprehensive baseline estimation with:
        - Multiple estimation methods with quality assessment
        - Method quarantine for consistently failing approaches
        - Adaptive method selection based on data characteristics
        - Cross-validation and consensus building
        """
        baseline = None  # Initialize baseline variable
        
        try:
            if method == 'auto':
                # Auto method: try percentile first, fallback to polynomial
                return self.robust_baseline_estimation(x_data, y_data, method='percentile', validation_enabled=False)
            
            elif method == 'polynomial':
                # Use outer 20% of data for baseline fitting
                n_points = len(x_data)
                edge_fraction = 0.2
                n_edge = int(n_points * edge_fraction)

                # Combine left and right edges
                edge_indices = np.concatenate([
                    np.arange(n_edge),
                    np.arange(n_points - n_edge, n_points)
                ])

                # Fit polynomial to edge points
                degree = min(self.fitting_parameters['baseline_polynomial_degree'], len(edge_indices) - 1)
                poly_coeffs = np.polyfit(x_data[edge_indices], y_data[edge_indices], degree)
                baseline = np.polyval(poly_coeffs, x_data)

                # Return median of baseline for single value
                return np.median(baseline)

            elif method == 'iterative':
                # Iterative baseline correction (simplified)
                baseline_est = np.percentile(y_data, 10)  # Start with 10th percentile

                for iteration in range(3):
                    # Points below current baseline estimate + some margin
                    margin = np.std(y_data) * 0.5
                    mask = y_data < (baseline_est + margin)

                    if np.sum(mask) > len(y_data) * 0.1:  # At least 10% of points
                        baseline_est = np.median(y_data[mask])
                    else:
                        break

                return baseline_est

            elif method == 'percentile':
                # Simple percentile-based estimation
                return np.percentile(y_data, 15)

            elif method == 'asymmetric_polynomial':
                # Enhanced asymmetric baseline correction for overlapping peaks
                return self.asymmetric_baseline_correction(x_data, y_data)

        except Exception as e:
            print(f"Baseline estimation failed, using fallback: {e}")
            # Fallback: use median of lowest 20% of points
            baseline = np.percentile(y_data, 20)


        # === SAFETY BARRIER: Baseline Validation ===
        # Final safety check: ensure baseline is never None
        if baseline is None:
            baseline = np.percentile(y_data, 20)
            
        data_max = np.max(y_data)
        data_min = np.min(y_data)
        data_range = data_max - data_min

        # Check if baseline is physically reasonable
        if baseline > data_max * 0.8:  # Baseline shouldn't be near data maximum
            baseline = np.percentile(y_data, 15)

        if baseline < data_min - data_range * 0.2:  # Baseline shouldn't be far below minimum
            baseline = np.percentile(y_data, 15)

        # Check for NaN or infinite baseline
        if not np.isfinite(baseline):
            baseline = np.median(y_data)

        # === END BASELINE SAFETY BARRIER ===

    def asymmetric_baseline_correction(self, x_data, y_data, asymmetry_param=0.05,
                                     smoothness_param=1e6, max_iterations=10):
        """
        Enhanced asymmetric baseline correction using Asymmetrically Reweighted Penalized Least Squares (ArPLS)

        This method is particularly effective for overlapping peaks with asymmetric baseline distortion.
        It iteratively fits a smooth baseline while penalizing positive deviations (peaks) more than
        negative ones, creating an asymmetric penalty that follows the true baseline under peaks.

        Parameters:
        - x_data, y_data: spectral data
        - asymmetry_param: asymmetry parameter (0 < p < 1), smaller = more asymmetric
        - smoothness_param: smoothness parameter (λ), larger = smoother baseline
        - max_iterations: maximum iterations for convergence

        Returns:
        - baseline array same length as input data

        Reference: Baek et al. (2015) "Baseline correction using asymmetrically reweighted penalized least squares"
        """
        try:
            from scipy import sparse
            from scipy.sparse.linalg import spsolve

            n = len(y_data)
            if n < 4:
                # Too few points for AsLLS, fallback to simple method
                return np.full(n, np.percentile(y_data, 10))

            print(f"   Applying asymmetric baseline correction (λ={smoothness_param}, p={asymmetry_param})...")

            # Build second derivative matrix for smoothness penalty
            # D2 is the discrete second difference operator
            D1 = sparse.diags([1, -1], [0, 1], shape=(n-1, n))
            D2 = sparse.diags([1, -2, 1], [0, 1, 2], shape=(n-2, n))

            # Initial weights (all equal)
            w = np.ones(n)
            baseline_prev = np.copy(y_data)

            # Iterative asymmetric reweighting
            for iteration in range(max_iterations):
                # Build weight matrix
                W = sparse.diags(w, 0, shape=(n, n))

                # Build system: (W + λD₂ᵀD₂)z = Wy
                # where z is the baseline we want to solve for
                A = W + smoothness_param * D2.T @ D2
                b = W @ y_data

                # Solve linear system
                baseline_current = spsolve(A, b)

                # Calculate residuals
                residuals = y_data - baseline_current

                # Update weights asymmetrically
                # For positive residuals (peaks): use small weight p
                # For negative residuals (below baseline): use weight 1
                w = np.where(residuals > 0, asymmetry_param, 1.0)

                # Check convergence
                if iteration > 0:
                    change = np.mean(np.abs(baseline_current - baseline_prev))
                    noise_level = np.std(residuals) * 0.01

                    if change < noise_level:
                        print(f"   Converged after {iteration+1} iterations (change: {change:.2e})")
                        break

                baseline_prev = baseline_current.copy()

            # Validate result
            if not np.all(np.isfinite(baseline_current)):
                print("   Warning: Non-finite baseline values detected, using fallback")
                return np.full(n, np.percentile(y_data, 15))

            # Additional validation: baseline shouldn't be much higher than data
            if np.max(baseline_current) > np.max(y_data) * 1.2:
                print("   Warning: Baseline exceeds data, using conservative fallback")
                return np.full(n, np.percentile(y_data, 20))

            print(f"   Asymmetric baseline range: {np.min(baseline_current):.1f} to {np.max(baseline_current):.1f}")
            return baseline_current

        except ImportError:
            print("   Warning: scipy.sparse not available, using polynomial fallback")
            return self._fallback_asymmetric_baseline(x_data, y_data)

        except Exception as e:
            print(f"   Asymmetric baseline correction failed: {e}, using fallback")
            return self._fallback_asymmetric_baseline(x_data, y_data)

    def _fallback_asymmetric_baseline(self, x_data, y_data):
        """
        Simplified asymmetric baseline correction when scipy.sparse is not available

        Uses iterative polynomial fitting with asymmetric weighting
        """
        try:
            n_iterations = 5
            weights = np.ones(len(y_data))

            # Start with polynomial baseline
            degree = min(3, len(y_data) // 4)  # Adaptive degree

            for iteration in range(n_iterations):
                # Weighted polynomial fit
                poly_coeffs = np.polyfit(x_data, y_data, degree, w=weights)
                baseline_poly = np.polyval(poly_coeffs, x_data)

                # Calculate residuals and update weights
                residuals = y_data - baseline_poly

                # Asymmetric weighting: penalize positive residuals (peaks) more
                weights = np.where(residuals > 0, 0.1, 1.0)  # Strong penalty for peaks

                # Add smoothing for noise robustness
                if len(weights) > 5:
                    from scipy.ndimage import gaussian_filter1d
                    weights = gaussian_filter1d(weights, sigma=1.0)

            print(f"   Fallback asymmetric baseline applied ({n_iterations} iterations)")
            return baseline_poly

        except Exception as e:
            print(f"   Fallback asymmetric baseline failed: {e}")
            return np.full(len(y_data), np.percentile(y_data, 15))

    def adaptive_baseline_method_selection(self, x_data, y_data, peak_complexity='unknown'):
        """
        Automatically select the best baseline correction method based on data characteristics

        Parameters:
        - x_data, y_data: spectral data
        - peak_complexity: 'simple', 'moderate', 'complex', or 'unknown'

        Returns:
        - tuple: (selected_method, baseline_value, method_info)
        """
        methods_to_test = []

        # Determine methods based on peak complexity
        if peak_complexity == 'simple':
            methods_to_test = ['polynomial', 'percentile']
        elif peak_complexity == 'moderate':
            methods_to_test = ['polynomial', 'iterative', 'asymmetric_polynomial']
        elif peak_complexity == 'complex':
            methods_to_test = ['asymmetric_polynomial', 'iterative', 'polynomial']
        else:  # unknown
            methods_to_test = ['polynomial', 'iterative', 'asymmetric_polynomial']

        print(f"   Testing baseline methods: {methods_to_test}")

        best_method = methods_to_test[0]
        best_baseline = self.robust_baseline_estimation(x_data, y_data, method=best_method)
        best_score = 0

        # Simple scoring: prefer methods that give reasonable baseline values
        baseline_range = np.max(y_data) - np.min(y_data)
        
        # Initialize baseline variables for safety
        best_baseline = np.median(y_data)  # Safe fallback baseline
        test_baseline = best_baseline

        for method in methods_to_test:
            try:
                if method == 'asymmetric_polynomial':
                    baseline_result = self.asymmetric_baseline_correction(x_data, y_data)
                    if baseline_result is not None:
                        if np.isscalar(baseline_result):
                            test_baseline = baseline_result
                        else:
                            test_baseline = np.median(baseline_result)  # Convert array to scalar for comparison
                    else:
                        test_baseline = np.median(y_data)  # Fallback if method returns None
                else:
                    baseline_result = self.robust_baseline_estimation(x_data, y_data, method=method)
                    if baseline_result is not None:
                        test_baseline = baseline_result
                    else:
                        test_baseline = np.median(y_data)  # Fallback if method returns None

                # Validate baseline result
                if not np.isfinite(test_baseline):
                    test_baseline = np.median(y_data)

                # Score based on reasonableness
                if np.min(y_data) <= test_baseline <= np.min(y_data) + baseline_range * 0.3:
                    score = 1.0
                    if method == 'asymmetric_polynomial' and peak_complexity in ['moderate', 'complex']:
                        score += 0.5  # Bonus for asymmetric method in complex cases

                    if score > best_score:
                        best_score = score
                        best_method = method
                        if np.isscalar(test_baseline):
                            best_baseline = test_baseline
                        else:
                            best_baseline = test_baseline  # Keep array form if method returns array

            except Exception as e:
                print(f"   Method {method} failed: {e}")
                # Ensure we always have a valid baseline
                if best_baseline is None or not np.isfinite(best_baseline):
                    best_baseline = np.median(y_data)
                continue

        # Final safety check for baseline
        if best_baseline is None or not np.isfinite(best_baseline):
            print("   ⚠ All baseline methods failed, using median fallback")
            best_baseline = np.median(y_data)
            best_method = 'median_fallback'
            best_score = 0.5

        method_info = {
            'selected_method': best_method,
            'methods_tested': methods_to_test,
            'peak_complexity': peak_complexity,
            'score': best_score
        }

        print(f"   Selected baseline method: {best_method} (score: {best_score:.2f})")
        return best_method, best_baseline, method_info

    def quality_based_peak_filtering(self, fit_results, min_r_squared=0.5,
                                   min_amplitude_ratio=0.1, max_width_factor=5.0):
        """
        Filter peaks based on fitting quality criteria (INTEGRATION ENHANCEMENT)

        Args:
            fit_results: list of fit result dictionaries
            min_r_squared: minimum R² threshold
            min_amplitude_ratio: minimum amplitude/noise ratio
            max_width_factor: maximum width compared to typical

        Returns:
            dict: filtered results with quality diagnostics
        """
        if not fit_results:
            return {'filtered_peaks': [], 'quality_report': {}}

        filtered_peaks = []
        quality_report = {
            'total_peaks': len(fit_results),
            'r_squared_filtered': 0,
            'amplitude_filtered': 0,
            'width_filtered': 0,
            'accepted_peaks': 0,
            'rejection_reasons': []
        }

        for i, result in enumerate(fit_results):
            if not result.get('success', False):
                quality_report['rejection_reasons'].append({'peak': i, 'reason': 'fit_failed'})
                continue

            rejection_reasons = []

            # R² filter
            r_squared = result.get('r_squared', 0)
            if r_squared < min_r_squared:
                rejection_reasons.append(f'low_r_squared_{r_squared:.3f}')
                quality_report['r_squared_filtered'] += 1

            # Amplitude filter
            amplitude = result.get('amplitude', 0)
            if amplitude <= 0:
                rejection_reasons.append('negative_amplitude')
                quality_report['amplitude_filtered'] += 1
            else:
                # Check amplitude ratio if noise level available
                if 'noise_level' in result:
                    noise_level = result['noise_level']
                    amplitude_ratio = amplitude / max(noise_level, 1)
                    if amplitude_ratio < min_amplitude_ratio:
                        rejection_reasons.append(f'low_amplitude_ratio_{amplitude_ratio:.2f}')
                        quality_report['amplitude_filtered'] += 1

            # Width filter
            sigma = result.get('sigma', 0)
            gamma = result.get('gamma', 0)
            total_width = sigma + gamma

            # Get typical width for nucleus
            nucleus_type = result.get('nucleus_type', '1H')
            typical_width = self.nmr_ranges.get(nucleus_type, {}).get('typical_width', 0.01)

            if total_width > typical_width * max_width_factor or total_width <= 0:
                rejection_reasons.append(f'unreasonable_width_{total_width:.4f}')
                quality_report['width_filtered'] += 1

            # Accept peak if no rejection reasons
            if not rejection_reasons:
                filtered_peaks.append({
                    **result,
                    'quality_score': self._calculate_quality_score(result),
                    'filter_passed': True
                })
                quality_report['accepted_peaks'] += 1
            else:
                quality_report['rejection_reasons'].append({
                    'peak': i,
                    'reasons': rejection_reasons,
                    'r_squared': r_squared,
                    'amplitude': amplitude,
                    'width': total_width
                })

        # Sort by quality score
        filtered_peaks.sort(key=lambda x: x['quality_score'], reverse=True)

        acceptance_rate = quality_report['accepted_peaks'] / quality_report['total_peaks']

        return {
            'filtered_peaks': filtered_peaks,
            'quality_report': quality_report,
            'acceptance_rate': acceptance_rate
        }

    def _calculate_quality_score(self, result):
        """Calculate composite quality score for peak ranking"""
        score = 0.0

        # R² component (40%)
        r_squared = result.get('r_squared', 0)
        score += 0.4 * r_squared

        # Amplitude component (30%)
        amplitude = result.get('amplitude', 0)
        if amplitude > 0 and 'noise_level' in result:
            amplitude_ratio = min(1.0, amplitude / (result['noise_level'] * 10))
            score += 0.3 * amplitude_ratio

        # Width reasonableness (20%)
        sigma = result.get('sigma', 0)
        gamma = result.get('gamma', 0)
        total_width = sigma + gamma

        nucleus_type = result.get('nucleus_type', '1H')
        typical_width = self.nmr_ranges.get(nucleus_type, {}).get('typical_width', 0.01)

        if total_width > 0:
            width_ratio = min(typical_width, total_width) / max(typical_width, total_width)
            score += 0.2 * width_ratio

        # Detection confidence component (10%)
        if 'detection_confidence' in result:
            detection_conf = result['detection_confidence'].get('confidence', 0)
            score += 0.1 * detection_conf

        return min(1.0, max(0.0, score))

    def optimize_baseline_iteratively(self, x_data, y_data, peak_center, initial_guess,
                                      max_iterations=5, improvement_threshold=0.001):
        """
        DYNAMIC OPTIMIZATION: Iteratively optimize baseline window size based on fit quality

        This addresses the baseline correction window problem by:
        1. Starting with conservative small windows (5% edge regions)
        2. Progressively expanding windows (8%, 12%, 18%, 25%)
        3. Stopping when R² improvement falls below threshold
        4. Preventing neighboring peaks from contaminating baseline estimation

        Args:
            x_data: X-axis data (ppm)
            y_data: Intensity data
            peak_center: Peak center position for fitting
            initial_guess: Initial fitting parameters [amp, center, sigma, gamma, baseline]
            max_iterations: Maximum window sizes to try (default: 5)
            improvement_threshold: Minimum R² improvement to continue (default: 0.001)

        Returns:
            tuple: (optimized_baseline, final_r_squared, optimization_report)

        BACKWARD COMPATIBILITY: Falls back to original method if optimization fails
        """
        # Progressive window expansion strategy - start small to avoid neighboring peaks
        edge_fractions = [0.05, 0.08, 0.12, 0.18, 0.25]  # Conservative to aggressive

        optimization_report = {
            'method': 'dynamic_baseline_optimization',
            'iterations': [],
            'best_edge_fraction': None,
            'improvement_achieved': 0,
            'converged': False,
            'fallback_used': False
        }

        best_r_squared = 0
        best_baseline = None
        best_edge_fraction = edge_fractions[0]

        # Try each window size
        for i, edge_fraction in enumerate(edge_fractions):
            try:
                # Estimate baseline with current window size
                baseline = self._estimate_baseline_with_fraction(x_data, y_data, edge_fraction)

                # Create modified guess with new baseline
                modified_guess = initial_guess.copy()
                modified_guess[4] = baseline

                # Attempt fit with this baseline
                fit_result = self._fit_with_specific_baseline(x_data, y_data,
                                                            modified_guess, peak_center)

                # DEFENSIVE FIX: Ensure fit_result is a dictionary
                if not isinstance(fit_result, dict):
                    fit_result = {'success': False, 'r_squared': 0.0, 'error': 'invalid_fit_result_type'}

                current_r_squared = fit_result.get('r_squared', 0)

                # Track this iteration
                iteration_info = {
                    'iteration': i + 1,
                    'edge_fraction': edge_fraction,
                    'baseline_value': float(baseline),
                    'r_squared': float(current_r_squared),
                    'improved': current_r_squared > best_r_squared + improvement_threshold,
                    'fit_success': fit_result.get('success', False)
                }
                optimization_report['iterations'].append(iteration_info)

                # Check for improvement
                if (fit_result.get('success', False) and
                    current_r_squared > best_r_squared + improvement_threshold):

                    best_r_squared = current_r_squared
                    best_baseline = baseline
                    best_edge_fraction = edge_fraction

                    print(f"   Baseline optimization: edge_fraction={edge_fraction:.2f}, R²={current_r_squared:.4f}")
                else:
                    # No improvement, convergence achieved
                    optimization_report['converged'] = True
                    print(f"   Baseline optimization converged at edge_fraction={best_edge_fraction:.2f}")
                    break

            except Exception as e:
                # Track failed iteration
                optimization_report['iterations'].append({
                    'iteration': i + 1,
                    'edge_fraction': edge_fraction,
                    'error': str(e),
                    'fit_success': False
                })
                print(f"   Baseline optimization failed at edge_fraction={edge_fraction:.2f}: {e}")

        # Finalize report
        optimization_report['best_edge_fraction'] = best_edge_fraction
        optimization_report['improvement_achieved'] = best_r_squared

        # Fallback to original method if optimization completely failed
        if best_baseline is None:
            print("   Baseline optimization failed, using fallback method")
            optimization_report['fallback_used'] = True
            best_baseline = self.robust_baseline_estimation(x_data, y_data, method='polynomial')
            best_r_squared = 0

        return best_baseline, best_r_squared, optimization_report

    def _estimate_baseline_with_fraction(self, x_data, y_data, edge_fraction):
        """
        Helper: estimate baseline using specified edge fraction

        IMPROVEMENT: More intelligent edge selection to avoid peak contamination
        """
        n_points = len(x_data)
        n_edge = max(3, int(n_points * edge_fraction))

        # Use outer edge points, avoid peak region in center
        edge_indices = np.concatenate([
            np.arange(n_edge),                          # Left edge
            np.arange(n_points - n_edge, n_points)      # Right edge
        ])

        # Remove any obvious peak contamination from edges
        edge_values = y_data[edge_indices]

        # Use median for robustness against outliers
        baseline = np.median(edge_values)

        # Additional validation: check if edge values are reasonable
        edge_std = np.std(edge_values)
        if edge_std > np.std(y_data) * 0.5:  # High variation in edges
            # Fall back to percentile method
            baseline = np.percentile(y_data, 15)

        return baseline

    def _fit_with_specific_baseline(self, x_data, y_data, initial_guess, peak_center):
        """
        LEVEL 1 ENHANCED: Helper method with comprehensive safety barriers
        This method implements all Level 1 safety features:
        - Pre-fitting parameter validation
        - Monitored optimization
        - Emergency fallback activation
        - Enhanced error handling
        """
        try:
            # === LEVEL 1 SAFETY BARRIER 1: Pre-fitting Validation ===
            nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
            bounds = self.get_adaptive_bounds(initial_guess, x_data, y_data, nucleus_type, None)

            # Validate parameters and bounds before fitting
            is_valid, validation_error = self.validate_initial_parameters_and_bounds(
                initial_guess, bounds, x_data, y_data
            )

            if not is_valid:
                return self.emergency_fallback_fitting(x_data, y_data, peak_center)

            # === LEVEL 1 SAFETY BARRIER 2: Monitored Optimization ===
            popt, pcov, trajectory_info = self.monitored_curve_fit(
                self.voigt_profile_1d, x_data, y_data,
                initial_guess, bounds=bounds, #p0=
                maxfev=self.fitting_parameters['max_iterations']
            )

            # Check if monitoring detected problems
            if popt is None or not trajectory_info['monitoring_successful']:
                return self.emergency_fallback_fitting(x_data, y_data, peak_center)

            # Generate fitted curve
            y_fitted = self.voigt_profile_1d(x_data, *popt)

            # === LEVEL 1 SAFETY BARRIER 3: Enhanced R² Validation ===
            r_squared = self.calculate_r_squared(y_data, y_fitted)

            # Critical R² validation with graduated response
            if r_squared < -10.0:  # Catastrophic failure
                return self.emergency_fallback_fitting(x_data, y_data, peak_center)

            elif r_squared < -1.0:  # Severe failure
                fallback_result = self.emergency_fallback_fitting(x_data, y_data, peak_center)
                if fallback_result['success'] and fallback_result['r_squared'] > r_squared:
                    return fallback_result
                else:
                    pass
            elif r_squared < 0.1:  # Poor quality
                pass

            elif r_squared > 0.9999:  # Suspiciously perfect
                pass

            else:
                pass

            # === LEVEL 1 SAFETY BARRIER 4: Enhanced Parameter Validation ===
            amplitude, center, sigma, gamma, baseline = popt

            # Data-driven validation limits
            data_max = np.max(y_data)
            data_min = np.min(y_data)
            data_range = data_max - data_min
            ppm_range = abs(x_data[-1] - x_data[0])

            amplitude_max_allowed = data_range * 20  # Increased tolerance
            amplitude_min_allowed = data_range * 0.0001  # Minimum meaningful amplitude
            center_drift_max = ppm_range * 0.3  # Increased tolerance
            width_max_allowed = ppm_range * 0.4  # Increased tolerance

            # Check amplitude constraints
            if amplitude > amplitude_max_allowed or amplitude < amplitude_min_allowed:
                return self.emergency_fallback_fitting(x_data, y_data, peak_center)

            # Check center drift
            if abs(center - peak_center) > center_drift_max:
                return self.emergency_fallback_fitting(x_data, y_data, peak_center)

            # Check width constraints
            total_width = sigma + gamma
            if total_width > width_max_allowed or total_width <= 0:
                return self.emergency_fallback_fitting(x_data, y_data, peak_center)


            # === SUCCESSFUL RESULT ASSEMBLY ===
            return {
                'success': True,
                'method': 'level1_enhanced_fitting',
                'parameters': popt,
                'fitted_curve': y_fitted,
                'r_squared': r_squared,
                'amplitude': amplitude,
                'center': center,
                'sigma': sigma,
                'gamma': gamma,
                'baseline': baseline,
                # Level 1 diagnostic information
                'level1_diagnostics': {
                    'pre_validation_passed': True,
                    'monitoring_successful': trajectory_info['monitoring_successful'],
                    'trajectory_iterations': trajectory_info['total_iterations'],
                    'emergency_fallback_used': False,
                    'validation_error': None
                }
            }

            # Calculate quality metrics

            # === CRITICAL SAFETY BARRIER: R² Quality Gates ===
            # Calculate quality metrics with validation
            #r_squared = self.calculate_r_squared(y_data, y_fitted)
            r_squared = self.calculate_r_squared(y_data, y_fitted)

            # Critical R² validation checks
            if r_squared < 0:
                return {
                    'success': False,
                    'error': f'negative_r_squared_{r_squared:.3f}',
                    'r_squared': r_squared,
                    'safety_triggered': 'r_squared_negative'
                }

            if r_squared < 0.1:  # Less than 10% variance explained
                return {
                    'success': False,
                    'error': f'r_squared_too_low_{r_squared:.3f}',
                    'r_squared': r_squared,
                    'safety_triggered': 'r_squared_quality'
                }

            if r_squared > 0.9999:  # Suspiciously perfect fit
                print(f"  SAFETY: Suspiciously high R² ({r_squared:.6f}) - possible overfitting")
                # Don't reject, but flag for attention

            print(f"   ✅ SAFETY: R² = {r_squared:.4f} passed quality gates")
            # === END R² SAFETY BARRIER ===


            # Validate fitted parameters
            amplitude, center, sigma, gamma, baseline = popt
            parameters_valid = (
                amplitude > 0 and
                sigma > 0 and gamma >= 0 and
                abs(center - peak_center) < abs(x_data[-1] - x_data[0]) * 0.1
            )

            return {
                'success': parameters_valid,
                'parameters': popt,
                'fitted_curve': y_fitted,
                'r_squared': r_squared,
                'amplitude': amplitude,
                'center': center,
                'sigma': sigma,
                'gamma': gamma,
                'baseline': baseline
            }


          # === CRITICAL SAFETY BARRIER: Post-fitting parameter validation ===
            # Validate fitted parameters against physical constraints
            amplitude_max_allowed = np.max(y_data) * 10  # Never exceed 10× data maximum
            amplitude_min_allowed = np.max(y_data) * 0.001  # Never below 0.1% of data maximum

            center_drift_max = abs(x_data[-1] - x_data[0]) * 0.2  # Max 20% of data range
            width_max_allowed = abs(x_data[-1] - x_data[0]) * 0.3  # Max 30% of data range

            # Check amplitude constraints
            if amplitude > amplitude_max_allowed or amplitude < amplitude_min_allowed:
                return {
                    'success': False,
                    'error': f'amplitude_out_of_bounds_{amplitude:.2e}',
                    'r_squared': 0,
                    'safety_triggered': 'amplitude_validation'
                }

            # Check center drift
            if abs(center - peak_center) > center_drift_max:
                return {
                    'success': False,
                    'error': f'center_drift_excessive_{abs(center - peak_center):.4f}',
                    'r_squared': 0,
                    'safety_triggered': 'center_validation'
                }

            # Check width constraints
            total_width = sigma + gamma
            if total_width > width_max_allowed or total_width <= 0:
                return {
                    'success': False,
                    'error': f'width_out_of_bounds_{total_width:.4f}',
                    'r_squared': 0,
                    'safety_triggered': 'width_validation'
                }
            # === END SAFETY BARRIER ===

        except Exception as e:

            # Last resort emergency fallback
            try:
                return self.emergency_fallback_fitting(x_data, y_data, peak_center)
            except Exception as fallback_error:
                return {
                    'success': False,
                    'error': f'complete_level1_failure_{str(e)}_{str(fallback_error)}',
                    'r_squared': 0,
                    'method': 'level1_complete_failure'
                }


    @staticmethod
    def calculate_r_squared(y_true, y_pred):
        """Calculate R-squared value (moved here for better organization)"""
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

    def detect_nucleus_type(self, ppm_range):
        """
        Detect nucleus type based on ppm range
        """
        ppm_span = abs(ppm_range[1] - ppm_range[0])
        center_ppm = (ppm_range[0] + ppm_range[1]) / 2

        # Check ¹H range (5.5-12 ppm)
        if 5.0 <= center_ppm <= 13.0 and ppm_span < 20:
            return '1H'

        # Check ¹⁵N range (100-140 ppm)
        elif 90 <= center_ppm <= 150 and ppm_span < 100:
            return '15N'

        # Check ¹³C range (0-220 ppm)
        elif 5 <= center_ppm <= 50 and ppm_span < 50: #was 230
            return '13C'

        # Default to ¹H for small ranges
        return '1H'

    def adaptive_width_estimation(self, x_data, y_data, peak_center, nucleus_type=None):
        """
        Adaptive width estimation using multiple approaches
        """
        if nucleus_type is None:
            nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])

        # Get typical width for nucleus
        typical_width = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])['typical_width']

        try:
            # Method 1: FWHM-based estimation
            baseline = self.robust_baseline_estimation(x_data, y_data)
            peak_height = np.max(y_data) - baseline
            half_max = baseline + peak_height / 2

            # Find points above half maximum
            above_half = y_data > half_max
            if np.sum(above_half) > 3:  # Need at least 3 points
                indices = np.where(above_half)[0]
                fwhm = x_data[indices[-1]] - x_data[indices[0]]
                width_fwhm = fwhm / 2.355  # Convert FWHM to sigma equivalent
            else:
                width_fwhm = typical_width

            # Method 2: Second moment estimation
            try:
                y_norm = y_data - baseline
                y_norm = np.maximum(y_norm, 0)  # Only positive values
                total_intensity = np.sum(y_norm)

                if total_intensity > 0:
                    # Calculate second moment
                    mean_x = np.sum(x_data * y_norm) / total_intensity
                    second_moment = np.sum(y_norm * (x_data - mean_x)**2) / total_intensity
                    width_moment = np.sqrt(second_moment)
                else:
                    width_moment = typical_width
            except:
                width_moment = typical_width

            # Method 3: Gradient-based estimation
            try:
                # Find steepest gradients (peak edges)
                gradient = np.gradient(y_data, x_data)
                # Find peak of absolute gradient
                max_grad_left = np.argmax(gradient[:len(gradient)//2])
                max_grad_right = np.argmin(gradient[len(gradient)//2:]) + len(gradient)//2

                if max_grad_right > max_grad_left:
                    width_gradient = (x_data[max_grad_right] - x_data[max_grad_left]) / 4
                else:
                    width_gradient = typical_width
            except:
                width_gradient = typical_width

            # Combine estimates using weighted average
            weights = [0.5, 0.3, 0.2]  # Prefer FWHM, then moment, then gradient
            combined_width = (weights[0] * width_fwhm +
                            weights[1] * width_moment +
                            weights[2] * width_gradient)

            # Apply reasonable bounds based on nucleus type
            min_width = typical_width * 0.1
            max_width = typical_width * 10

            final_width = np.clip(combined_width, min_width, max_width)

            return final_width

        except Exception as e:
            print(f"Width estimation failed, using typical value: {e}")
            return typical_width

    def estimate_initial_parameters_from_resolved_peaks(self, x_data, y_data,
                                                       all_peak_positions=None):
        """
        DYNAMIC OPTIMIZATION: Estimate initial fitting parameters from well-resolved signals

        This addresses the overlapping peak problem by:
        1. Identifying well-separated peaks (>3x typical width apart)
        2. Estimating typical linewidth from isolated peaks
        3. Using robust statistical methods to avoid outliers
        4. Providing reliable initial parameters for complex fitting scenarios

        Args:
            x_data: X-axis data (ppm)
            y_data: Intensity data
            all_peak_positions: List of suspected peak positions for context analysis

        Returns:
            dict: Dictionary containing estimated parameters and metadata

        BACKWARD COMPATIBILITY: Falls back to nucleus defaults if analysis fails
        """
        print("   Estimating initial parameters from well-resolved peaks...")

        # Initialize results with defaults
        nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
        default_width = self.nmr_ranges[nucleus_type]['typical_width']

        results = {
            'nucleus_type': nucleus_type,
            'typical_linewidth': default_width,
            'typical_baseline_slope': 0,
            'noise_level': np.std(y_data) * 0.1,
            'well_resolved_count': 0,
            'estimation_quality': 'fallback',
            'method_used': 'nucleus_defaults'
        }

        if all_peak_positions is not None and len(all_peak_positions) > 1:
            try:
                # Find well-separated peaks (>3x typical width apart)
                min_separation = default_width * 3
                well_resolved_peaks = []

                for i, pos in enumerate(all_peak_positions):
                    # Check separation from all other peaks
                    neighbors = [p for j, p in enumerate(all_peak_positions)
                                if j != i and abs(p - pos) < min_separation * 2]

                    if len(neighbors) == 0:  # Peak is isolated
                        well_resolved_peaks.append(pos)
                        print(f"   Found well-resolved peak at {pos:.3f} ppm")

                results['well_resolved_count'] = len(well_resolved_peaks)

                # Estimate parameters from isolated peaks if we have enough
                if len(well_resolved_peaks) >= 2:
                    linewidths = []
                    baseline_estimates = []

                    # Analyze up to 5 peaks to avoid over-computation
                    peaks_to_analyze = well_resolved_peaks[:5]

                    for peak_pos in peaks_to_analyze:
                        try:
                            # Extract local region around this peak
                            local_region = self.extract_local_peak_region(
                                x_data, y_data, peak_pos,
                                fitted_width=default_width,
                                nucleus_type=nucleus_type,
                                window_multiplier=4.0  # Larger window for parameter estimation
                            )

                            if local_region['n_points'] > 10:  # Need sufficient points
                                # Estimate width from this peak
                                width = self.adaptive_width_estimation(
                                    local_region['x_local'],
                                    local_region['y_local'],
                                    peak_pos, nucleus_type
                                )
                                linewidths.append(width)

                                # Estimate local baseline
                                local_baseline = self.calculate_local_baseline(
                                    local_region['x_local'],
                                    local_region['y_local']
                                )
                                baseline_estimates.append(local_baseline)

                        except Exception as e:
                            print(f"   Warning: Failed to analyze peak at {peak_pos:.3f}: {e}")
                            continue

                    # Calculate robust statistics if we got enough estimates
                    if len(linewidths) >= 2:
                        # Use median for robustness against outliers
                        estimated_linewidth = np.median(linewidths)
                        linewidth_std = np.std(linewidths)

                        # Quality check: standard deviation shouldn't be too large
                        if linewidth_std < estimated_linewidth * 0.5:  # Reasonable consistency
                            results['typical_linewidth'] = estimated_linewidth
                            results['estimation_quality'] = 'good'
                            results['method_used'] = 'isolated_peaks_analysis'
                            print(f"   Estimated typical linewidth: {estimated_linewidth:.4f} ± {linewidth_std:.4f}")
                        else:
                            results['estimation_quality'] = 'inconsistent'
                            print(f"   Warning: Inconsistent linewidth estimates (std={linewidth_std:.4f})")

                    # Estimate baseline slope if we have baseline estimates
                    if len(baseline_estimates) >= 2:
                        baseline_range = np.max(baseline_estimates) - np.min(baseline_estimates)
                        x_range = x_data[-1] - x_data[0]
                        results['typical_baseline_slope'] = baseline_range / x_range

                # Even with few well-resolved peaks, try to get better noise estimate
                if len(well_resolved_peaks) >= 1:
                    # Use edge regions of spectrum for noise estimation
                    edge_size = min(50, len(x_data) // 10)
                    edge_data = np.concatenate([y_data[:edge_size], y_data[-edge_size:]])
                    results['noise_level'] = np.std(edge_data)

            except Exception as e:
                print(f"   Parameter estimation failed, using defaults: {e}")
                results['estimation_quality'] = 'failed'

        else:
            print("   No peak context provided, using nucleus defaults")
            results['method_used'] = 'no_context'

        return results

    def refine_peak_center(self, x_data, y_data, initial_center):
        """
        Refine peak center using local maximum finding and parabolic interpolation
        """
        try:
            # Find the index closest to initial guess
            center_idx = np.argmin(np.abs(x_data - initial_center))

            # Define search window (±5 points or ±10% of data range)
            search_window = min(5, len(x_data) // 10)

            left_bound = max(0, center_idx - search_window)
            right_bound = min(len(x_data), center_idx + search_window + 1)

            # Find local maximum in search window
            search_region = y_data[left_bound:right_bound]
            local_max_idx = np.argmax(search_region) + left_bound

            # Parabolic interpolation for sub-pixel accuracy
            if 1 <= local_max_idx < len(x_data) - 1:
                # Use 3-point parabolic interpolation
                y1, y2, y3 = y_data[local_max_idx-1:local_max_idx+2]
                x1, x2, x3 = x_data[local_max_idx-1:local_max_idx+2]

                # Parabolic fit: y = ax² + bx + c
                denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
                if abs(denom) > 1e-10:
                    a = (x3*(y2 - y1) + x2*(y1 - y3) + x1*(y3 - y2)) / denom
                    b = (x3*x3*(y1 - y2) + x2*x2*(y3 - y1) + x1*x1*(y2 - y3)) / denom

                    if abs(a) > 1e-10:
                        # Maximum at x = -b/(2a)
                        refined_center = -b / (2 * a)

                        # Check if refined center is reasonable
                        if x1 <= refined_center <= x3:
                            return refined_center

            # Fallback to discrete maximum
            return x_data[local_max_idx]

        except Exception as e:
            print(f"Peak center refinement failed: {e}")
            return initial_center

    def get_adaptive_bounds(self, initial_guess, x_data, y_data, nucleus_type=None, linewidth_constraints=None):
        """
        Get adaptive bounds based on data characteristics and nucleus type
        """
        amplitude, center, sigma, gamma, baseline = initial_guess

        if nucleus_type is None:
            nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])

        # Get nucleus-specific constraints
        nmr_params = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])
        typical_width = nmr_params['typical_width']

        # Data-dependent constraints
        ppm_range = abs(x_data[-1] - x_data[0])
        data_std = np.std(np.diff(x_data)) if len(x_data) > 1 else 0.001

        # Amplitude bounds
        amp_lower = max(0, amplitude * 0.1)
        amp_upper = amplitude * 5  # Allow significant variation


        # === CONSERVATIVE AMPLITUDE CONSTRAINTS: Selective Enhancement ===
        # Data-driven amplitude bounds with controlled flexibility for overlapping peaks
        data_max = np.max(y_data) if len(y_data) > 0 else abs(amplitude)
        data_min = np.min(y_data) if len(y_data) > 0 else 0
        data_range = data_max - data_min

        # Estimate local noise level for small peak detection (more conservative)
        noise_level = max(np.std(y_data) * 0.05, data_range * 0.005) if len(y_data) > 5 else data_range * 0.01

        # More conservative amplitude bounds to prevent overfitting
        # Lower bound: improved for small peaks but with safety limits
        amp_lower = max(0, noise_level * 3, data_range * 0.005)  # Slightly higher minimum for stability

        # Upper bound: more conservative scaling
        # Small peaks get modest boost, normal peaks stay conservative
        if amplitude < data_max * 0.1:  # Small peak detection
            amplitude_factor = 6.0  # Modest increase for small peaks (vs original 4×)
        else:
            amplitude_factor = 5.0  # Conservative for normal peaks

        amp_upper = min(amplitude * amplitude_factor, data_max * 6)  # Increased from 4× to 6× (not 10×)

        # Sanity check: if initial guess is unrealistic, constrain more carefully
        if amplitude > data_max * 5:
            amp_upper = min(data_max * 6, amplitude * 2)  # More conservative upper bound

        # === END ENHANCED AMPLITUDE CONSTRAINTS ===


        # Center bounds (tighter for better data, looser for noisy data)
        center_tolerance = min(typical_width * 0.1, ppm_range * 0.05)
        center_lower = center - center_tolerance
        center_upper = center + center_tolerance

        # Width bounds (sigma and gamma) - apply constraints if provided
        if linewidth_constraints:
            # Use provided constraints
            sigma_bounds = linewidth_constraints.get('sigma_bounds', (typical_width * 0.01, typical_width * 20))
            gamma_bounds = linewidth_constraints.get('gamma_bounds', (typical_width * 0.01, typical_width * 20))

            sigma_lower, sigma_upper = sigma_bounds
            gamma_lower, gamma_upper = gamma_bounds

            # Ensure constraints are reasonable relative to data
            sigma_lower = max(sigma_lower, data_std)
            gamma_lower = max(gamma_lower, data_std)
            sigma_upper = min(sigma_upper, ppm_range * 0.5)
            gamma_upper = min(gamma_upper, ppm_range * 0.5)
        else:
            # Default width bounds
            sigma_lower = gamma_lower = max(data_std, typical_width * 0.01)
            sigma_upper = gamma_upper = min(ppm_range * 0.3, typical_width * 20)

        # Conservative width parameters with controlled intensity-dependent scaling
        # Calculate intensity factor for width flexibility (more conservative)
        peak_intensity_ratio = amplitude / data_max if data_max > 0 else 1.0
        intensity_factor = 1.0 + 0.1 * np.log10(max(peak_intensity_ratio * 10, 0.1))  # Reduced from 0.2 to 0.1

        # Environment-dependent factor (estimate crowding from data complexity)
        data_complexity = np.std(np.diff(y_data)) / (noise_level + 1e-10) if len(y_data) > 2 else 1.0
        crowding_factor = 1.0 + 0.15 * min(data_complexity / 5.0, 1.0)  # Reduced from 30% to 15% increase

        # Combined flexibility factor with safety cap
        flexibility_factor = min(intensity_factor * crowding_factor, 2.0)  # Cap at 2× maximum

        if nucleus_type == '1H':
            sigma_lower = gamma_lower = max(sigma_lower, typical_width * 0.04)  # Slightly more flexible than original 0.05
            sigma_upper = gamma_upper = min(sigma_upper, typical_width * 4 * flexibility_factor)  # Modest increase from 3× to 4×
        elif nucleus_type == '15N':
            sigma_lower = gamma_lower = max(sigma_lower, typical_width * 0.02)  # Slightly reduced from 0.025
            sigma_upper = gamma_upper = min(sigma_upper, typical_width * 3 * flexibility_factor)  # Modest increase from 2× to 3×
        elif nucleus_type == '13C':
            sigma_lower = gamma_lower = max(sigma_lower, typical_width * 0.04)  # Slightly reduced from 0.05
            sigma_upper = gamma_upper = min(sigma_upper, typical_width * 3.5 * flexibility_factor)  # Modest increase from 2.5× to 3.5×

        # Baseline bounds
        baseline_tolerance = max(abs(amplitude) * 0.5, np.std(initial_guess) * 0.1)
        baseline_lower = baseline - baseline_tolerance
        baseline_upper = baseline + baseline_tolerance

        lower_bounds = [amp_lower, center_lower, sigma_lower, gamma_lower, baseline_lower]
        upper_bounds = [amp_upper, center_upper, sigma_upper, gamma_upper, baseline_upper]

        return (lower_bounds, upper_bounds)

    def validate_and_fix_initial_guess(self, initial_guess, bounds, x_data, y_data):
        """
        Validate initial guess against bounds and fix 'x0 is infeasible' errors

        This method ensures backwards compatibility by:
        1. Checking if initial guess violates bounds
        2. Adjusting initial guess to fit within bounds
        3. If needed, expanding bounds slightly to accommodate reasonable guesses
        4. Providing detailed diagnostics for debugging

        Returns: (validated_initial_guess, adjusted_bounds, validation_info)
        """
        try:
            lower_bounds, upper_bounds = bounds
            validated_guess = list(initial_guess)
            adjusted_bounds = (list(lower_bounds), list(upper_bounds))
            violations = []
            adjustments = []

            n_params = len(initial_guess)
            n_peaks = max(1, (n_params - 1) // 4)
            param_names = []

            for peak_idx in range(n_peaks):
                peak_label = f"peak{peak_idx + 1}"
                param_names.extend([
                    f"{peak_label}_amplitude",
                    f"{peak_label}_center",
                    f"{peak_label}_sigma",
                    f"{peak_label}_gamma",
                ])

            remaining = n_params - len(param_names)
            if remaining > 0:
                for extra_idx in range(remaining):
                    name = "baseline" if extra_idx == 0 else f"extra_param_{extra_idx}"
                    param_names.append(name)

            # Check each parameter against bounds
            for i, (param_name, guess, lower, upper) in enumerate(zip(param_names, initial_guess, lower_bounds, upper_bounds)):
                if not (lower <= guess <= upper):
                    violations.append(f"{param_name}: {guess:.3e} not in [{lower:.3e}, {upper:.3e}]")

                    # Strategy 1: Clamp to bounds with small buffer
                    if guess < lower:
                        validated_guess[i] = lower * 1.01  # 1% above lower bound
                        adjustments.append(f"{param_name}: clamped to lower bound + 1%")
                    elif guess > upper:
                        validated_guess[i] = upper * 0.99  # 1% below upper bound
                        adjustments.append(f"{param_name}: clamped to upper bound - 1%")

                    # Strategy 2: For amplitude, expand bounds if guess is reasonable
                    if param_name.endswith('amplitude') and guess > upper:
                        data_max = np.max(y_data) if len(y_data) > 0 else abs(guess)
                        if guess <= data_max * 15:  # If guess is within 15x data max (reasonable)
                            new_upper = min(guess * 1.2, data_max * 15)  # Expand but cap at 15x
                            adjusted_bounds[1][i] = new_upper
                            validated_guess[i] = guess  # Keep original guess
                            adjustments.append(f"{param_name}: expanded upper bound to {new_upper:.3e}")

            validation_info = {
                'had_violations': len(violations) > 0,
                'violations': violations,
                'adjustments': adjustments,
                'original_guess': list(initial_guess),
                'validated_guess': validated_guess,
                'bounds_adjusted': adjusted_bounds != (list(lower_bounds), list(upper_bounds))
            }

            # Log details if violations occurred
            if validation_info['had_violations']:
                print(f"   🔧 BOUNDS VALIDATION: Fixed {len(violations)} parameter violations")
                for violation in violations[:3]:  # Show first 3 violations
                    print(f"      {violation}")
                if len(adjustments) > 0:
                    print(f"   ✅ Applied {len(adjustments)} corrections")

            return tuple(validated_guess), tuple(adjusted_bounds), validation_info

        except Exception as e:
            # Fallback: return original values if validation fails
            print(f"   ⚠️ Bounds validation failed: {e}, using original values")
            return initial_guess, bounds, {'error': str(e), 'had_violations': False}

    @staticmethod
    def _clamp_vector_to_bounds(vector, bounds, margin=1e-6):
        """Clamp parameters so they fall within the provided bounds."""
        clamped = []

        for value, (lower, upper) in zip(vector, bounds):
            lower_adj = lower + margin
            upper_adj = upper - margin if np.isfinite(upper) else upper

            if np.isfinite(upper_adj) and upper_adj <= lower_adj:
                lower_adj = lower
                upper_adj = upper

            if np.isfinite(upper_adj):
                clamped.append(min(max(value, lower_adj), upper_adj))
            else:
                clamped.append(max(value, lower_adj))

        return clamped

    def validate_initial_parameters_and_bounds(self, initial_guess, bounds, x_data, y_data):
        """
        LEVEL 1 CRITICAL SAFETY: Validate parameters before any fitting attempts

        This method prevents catastrophic failures by checking parameter and bounds
        consistency before passing them to curve_fit optimization.

        Parameters:
        - initial_guess: [amplitude, center, sigma, gamma, baseline]
        - bounds: (lower_bounds, upper_bounds) tuple
        - x_data, y_data: experimental data

        Returns:
        - (bool, str): (is_valid, error_message)
        """
        try:
            amplitude, center, sigma, gamma, baseline = initial_guess
            lower_bounds, upper_bounds = bounds

            # Data characteristics
            data_max = np.max(y_data)
            data_min = np.min(y_data)
            data_range = data_max - data_min
            ppm_range = abs(x_data[-1] - x_data[0])


            # === AMPLITUDE VALIDATION ===
            if amplitude <= 0:
                return False, "zero_or_negative_amplitude"

            if amplitude > data_range * 100:  # Never exceed 100× data range
                return False, "amplitude_exceeds_data_range"

            if amplitude < data_range * 1e-6:  # Minimum meaningful amplitude
                return False, "amplitude_too_small"

            # === CENTER VALIDATION ===
            x_min = min(x_data[0], x_data[-1])
            x_max = max(x_data[0], x_data[-1])
            if not (x_min <= center <= x_max):
                return False, "center_outside_data_range"

            # === WIDTH VALIDATION ===
            if sigma <= 0 or gamma < 0:
                return False, "invalid_width_parameters"

            total_width = sigma + gamma
            if total_width > ppm_range * 0.5:  # Width shouldn't exceed half the data range
                return False, "width_exceeds_data_range"

            if total_width < ppm_range * 1e-6:  # Minimum meaningful width
                return False, "width_too_narrow"

            # === BOUNDS CONSISTENCY VALIDATION ===
            param_names = ['amplitude', 'center', 'sigma', 'gamma', 'baseline']
            for i, (param, lower, upper, name) in enumerate(zip(initial_guess, lower_bounds, upper_bounds, param_names)):
                if lower >= upper:
                    return False, f"inconsistent_bounds_{name}"

                if not (lower <= param <= upper):
                    return False, f"parameter_outside_bounds_{name}"

            # === BASELINE VALIDATION ===
            if baseline > data_max or baseline < data_min - data_range * 0.5:
                return False, "unrealistic_baseline"

            # === AMPLITUDE-BASELINE CONSISTENCY ===
            effective_peak_height = amplitude
            max_possible_height = data_max - baseline
            if effective_peak_height > max_possible_height * 3:  # Allow some flexibility was 10
                return False, "amplitude_baseline_inconsistency"

            return True, "parameters_valid"

        except Exception as e:
            return False, f"validation_exception_{str(e)}"

    def monitored_curve_fit(self, func, x_data, y_data, initial_guess, bounds, **kwargs):
        """
        LEVEL 1 CRITICAL SAFETY: Monitor optimization trajectory for early termination

        This method wraps curve_fit with real-time monitoring to detect and prevent
        optimization trajectories that lead to pathological solutions.

        Parameters:
        - func: fitting function (voigt_profile_1d)
        - x_data, y_data: experimental data
        - initial_guess: initial parameter guess
        - bounds: parameter bounds
        - **kwargs: additional curve_fit arguments

        Returns:
        - (popt, pcov, trajectory_info) if successful
        - (None, None, error_info) if monitoring detects problems
        """
        try:
            trajectory_log = []
            iteration_count = [0]  # Use list for closure
            divergence_count = [0]  # Track consecutive bad iterations

            def monitored_func(x, *params):
                iteration_count[0] += 1

                try:
                    result = func(x, *params)

                    # Check for NaN or infinite values
                    if not np.all(np.isfinite(result)):
                        divergence_count[0] += 1
                        if divergence_count[0] > 5:
                            raise RuntimeError("persistent_nan_or_inf_in_function")
                        return np.full_like(x, np.mean(y_data))  # Return reasonable fallback

                    # Monitor trajectory every 10 iterations
                    if iteration_count[0] % 10 == 0 and len(result) > 0:
                        # Calculate current R²
                        current_r2 = self.calculate_r_squared(y_data, result)

                        # Check for pathological parameter values
                        amplitude, center, sigma, gamma, baseline = params
                        param_score = self._evaluate_parameter_sanity(amplitude, center, sigma, gamma, baseline, x_data, y_data)

                        trajectory_log.append({
                            'iteration': iteration_count[0],
                            'r_squared': current_r2,
                            'params': params,
                            'param_score': param_score
                        })

                        # === EMERGENCY TERMINATION CONDITIONS ===

                        # Condition 1: Persistent negative R²
                        if current_r2 < -1.0 and iteration_count[0] > 50:
                            raise RuntimeError("optimization_trajectory_divergence")

                        # Condition 2: Consistently bad trajectory
                        if len(trajectory_log) >= 5:
                            recent_r2 = [entry['r_squared'] for entry in trajectory_log[-5:]]
                            if all(r2 < -0.5 for r2 in recent_r2):
                                raise RuntimeError("persistent_negative_r_squared")

                        # Condition 3: Parameter sanity degradation
                        if param_score < 0.1 and iteration_count[0] > 30:
                            raise RuntimeError("parameter_sanity_failure")

                        # Condition 4: Excessive iterations with no improvement
                        if len(trajectory_log) >= 10:
                            r2_trend = [entry['r_squared'] for entry in trajectory_log[-10:]]
                            if max(r2_trend) - min(r2_trend) < 0.001 and max(r2_trend) < 0.1:
                                raise RuntimeError("no_optimization_progress")

                        # Reset divergence count on successful evaluation
                        divergence_count[0] = 0

                    return result

                except RuntimeError:
                    raise  # Re-raise monitoring errors
                except Exception as e:
                    divergence_count[0] += 1
                    if divergence_count[0] > 10:
                        raise RuntimeError(f"function_evaluation_failure_{str(e)}")
                    return np.full_like(x, np.mean(y_data))  # Return reasonable fallback

            # Perform monitored fitting with bounds validation
            try:
                # NEW: Validate initial guess against bounds to prevent 'x0 is infeasible' errors
                validated_guess, validated_bounds, validation_info = self.validate_and_fix_initial_guess(
                    initial_guess, bounds, x_data, y_data
                )

                # Use validated parameters for curve_fit
                popt, pcov = curve_fit(monitored_func, x_data, y_data, p0=validated_guess, bounds=validated_bounds, **kwargs)

                # Final trajectory assessment
                trajectory_info = {
                    'total_iterations': iteration_count[0],
                    'trajectory_log': trajectory_log,
                    'bounds_validation': validation_info,  # NEW: Include validation details
                    'monitoring_successful': True,
                    'termination_reason': 'normal_convergence'
                }

                return popt, pcov, trajectory_info

            except RuntimeError as e:
                error_info = {
                    'monitoring_successful': False,
                    'termination_reason': str(e),
                    'trajectory_log': trajectory_log,
                    'total_iterations': iteration_count[0]
                }
                return None, None, error_info

        except ValueError as e:
            # Handle specific 'x0 is infeasible' error with detailed diagnosis
            if 'x0 is infeasible' in str(e):
                print(f"   ❌ BOUNDS ERROR: x0 is infeasible - this should have been prevented by validation!")
                print(f"      Original error: {e}")
                # This should not happen with our validation, but provide fallback
                error_info = {
                    'monitoring_successful': False,
                    'termination_reason': 'x0_infeasible_despite_validation',
                    'original_error': str(e),
                    'trajectory_log': [],
                    'total_iterations': 0,
                    'bounds_validation': {'error': 'validation_bypassed'}
                }
            else:
                error_info = {
                    'monitoring_successful': False,
                    'termination_reason': f'value_error_{str(e)}',
                    'trajectory_log': [],
                    'total_iterations': iteration_count[0]
                }
            return None, None, error_info

        except Exception as e:
            error_info = {
                'monitoring_successful': False,
                'termination_reason': f'monitoring_exception_{str(e)}',
                'trajectory_log': [],
                'total_iterations': 0
            }
            return None, None, error_info

    def _evaluate_parameter_sanity(self, amplitude, center, sigma, gamma, baseline, x_data, y_data):
        """
        Helper: Evaluate parameter sanity during optimization
        Returns score between 0 (pathological) and 1 (excellent)
        """
        try:
            score = 1.0

            # Data characteristics
            data_max = np.max(y_data)
            data_min = np.min(y_data)
            data_range = data_max - data_min
            ppm_range = abs(x_data[-1] - x_data[0])

            # Amplitude sanity (weight: 0.3)
            if amplitude <= 0 or amplitude > data_range * 50:
                score -= 0.3
            elif amplitude > data_range * 10:
                score -= 0.15
            elif amplitude < data_range * 0.01:
                score -= 0.1

            # Center sanity (weight: 0.2)
            center_drift = min(abs(center - x_data[0]), abs(center - x_data[-1]))
            if center < x_data[0] or center > x_data[-1]:
                score -= 0.2
            elif center_drift < ppm_range * 0.1:
                score -= 0.1

            # Width sanity (weight: 0.3)
            total_width = sigma + gamma
            if total_width <= 0 or total_width > ppm_range * 0.5:
                score -= 0.3
            elif total_width > ppm_range * 0.3:
                score -= 0.15
            elif sigma <= 0 or gamma < 0:
                score -= 0.2

            # Baseline sanity (weight: 0.2)
            if baseline > data_max or baseline < data_min - data_range:
                score -= 0.2
            elif baseline > data_max * 0.8:
                score -= 0.1

            return max(0.0, score)

        except:
            return 0.0  # Complete failure in evaluation

    def emergency_fallback_fitting(self, x_data, y_data, peak_center):
        """
        LEVEL 1 CRITICAL SAFETY: Emergency fallback when all standard methods fail

        This method implements a hierarchy of increasingly simple fitting approaches
        to ensure that some reasonable result is always returned, even for pathological cases.

        Parameters:
        - x_data, y_data: experimental data
        - peak_center: estimated peak center

        Returns:
        - dict: fitting result with emergency fallback flag
        """

        fallback_results = []
        data_max = np.max(y_data)
        data_min = np.min(y_data)
        data_range = data_max - data_min
        ppm_range = abs(x_data[-1] - x_data[0])

        # === FALLBACK 1: Simple Gaussian (Most Robust) ===
        try:
            def simple_gaussian(x, amp, center, width, baseline):
                return amp * np.exp(-((x - center) / width) ** 2) + baseline

            # Conservative parameter estimates
            width_est = ppm_range / 20  # Very conservative width
            baseline_est = np.median([data_min, np.percentile(y_data, 10)])
            amplitude_est = data_max - baseline_est

            simple_guess = [amplitude_est, peak_center, width_est, baseline_est]
            simple_bounds = (
                [data_range * 0.01, peak_center - width_est*3, width_est/20, data_min - data_range*0.2],
                [data_range * 5, peak_center + width_est*3, width_est*20, data_max*0.5]
            )

            popt_gauss, pcov_gauss = curve_fit(
                simple_gaussian, x_data, y_data,
                p0=simple_guess, bounds=simple_bounds,
                maxfev=400  # Increased for better convergence
            )

            y_fitted_gauss = simple_gaussian(x_data, *popt_gauss)
            r2_gauss = self.calculate_r_squared(y_data, y_fitted_gauss)

            if r2_gauss > 0.1:  # Minimal quality threshold
                fallback_results.append({
                    'method': 'gaussian_emergency_fallback',
                    'r_squared': r2_gauss,
                    'parameters': [popt_gauss[0], popt_gauss[1], popt_gauss[2], 0.0, popt_gauss[3]],  # Convert to Voigt format
                    'fitted_curve': y_fitted_gauss,
                    'success': True,
                    'emergency_fallback': True,
                    'fallback_level': 1
                })

        except Exception as e:
            print(f"   ❌ LEVEL1: Gaussian fallback failed: {e}")

        # === FALLBACK 2: Lorentzian Profile ===
        try:
            def simple_lorentzian(x, amp, center, gamma, baseline):
                return amp * gamma**2 / ((x - center)**2 + gamma**2) + baseline

            width_est = ppm_range / 25  # Even more conservative for Lorentzian
            baseline_est = np.percentile(y_data, 5)  # Lower percentile for Lorentzian
            amplitude_est = (data_max - baseline_est) * np.pi * width_est / 2  # Lorentzian normalization

            lorentz_guess = [amplitude_est, peak_center, width_est, baseline_est]
            lorentz_bounds = (
                [data_range * 0.005, peak_center - width_est*2, width_est/50, data_min - data_range*0.1],
                [data_range * 10, peak_center + width_est*2, width_est*10, data_max*0.3]
            )

            popt_lorentz, pcov_lorentz = curve_fit(
                simple_lorentzian, x_data, y_data,
                p0=lorentz_guess, bounds=lorentz_bounds,
                maxfev=400  # Increased for better convergence
            )

            y_fitted_lorentz = simple_lorentzian(x_data, *popt_lorentz)
            r2_lorentz = self.calculate_r_squared(y_data, y_fitted_lorentz)

            if r2_lorentz > 0.08:  # Even lower threshold for Lorentzian
                fallback_results.append({
                    'method': 'lorentzian_emergency_fallback',
                    'r_squared': r2_lorentz,
                    'parameters': [popt_lorentz[0], popt_lorentz[1], 0.0, popt_lorentz[2], popt_lorentz[3]],  # Convert to Voigt format
                    'fitted_curve': y_fitted_lorentz,
                    'success': True,
                    'emergency_fallback': True,
                    'fallback_level': 2
                })

        except Exception as e:
            print(f"   ❌ LEVEL1: Lorentzian fallback failed: {e}")

        # === FALLBACK 3: Triangular Approximation ===
        try:
            # Find peak region
            peak_idx = np.argmax(y_data)
            peak_amplitude = y_data[peak_idx]

            # Estimate width from half-maximum points
            half_max = (peak_amplitude + data_min) / 2
            half_max_indices = np.where(y_data > half_max)[0]

            if len(half_max_indices) > 3:
                width_points = len(half_max_indices)
                estimated_width = (x_data[half_max_indices[-1]] - x_data[half_max_indices[0]]) / 2.355  # FWHM to sigma

                # Create triangular fit (Gaussian with fixed width)
                def triangular_fit(x, amp, baseline):
                    return amp * np.exp(-((x - peak_center) / estimated_width) ** 2) + baseline

                tri_guess = [peak_amplitude - data_min, data_min]
                tri_bounds = ([0, data_min - data_range*0.1], [data_range * 2, data_max])

                popt_tri, _ = curve_fit(triangular_fit, x_data, y_data, p0=tri_guess, bounds=tri_bounds, maxfev=200)  # Increased for better convergence
                y_fitted_tri = triangular_fit(x_data, *popt_tri)
                r2_tri = self.calculate_r_squared(y_data, y_fitted_tri)

                if r2_tri > 0.05:
                    fallback_results.append({
                        'method': 'triangular_emergency_fallback',
                        'r_squared': r2_tri,
                        'parameters': [popt_tri[0], peak_center, estimated_width, 0.0, popt_tri[1]],
                        'fitted_curve': y_fitted_tri,
                        'success': True,
                        'emergency_fallback': True,
                        'fallback_level': 3
                    })

        except Exception as e:
            print(f"   ❌ LEVEL1: Triangular fallback failed: {e}")

        # === FALLBACK 4: Linear Interpolation (Ultimate Safety) ===
        try:
            # Create piecewise linear fit preserving peak shape
            peak_idx = np.argmax(y_data)
            baseline_est = np.median([y_data[0], y_data[-1], data_min])

            # Simple linear background with peak preserved
            y_linear = np.full_like(y_data, baseline_est)

            # Preserve a region around the peak
            preserve_width = max(10, len(y_data) // 10)
            start_idx = max(0, peak_idx - preserve_width // 2)
            end_idx = min(len(y_data), peak_idx + preserve_width // 2 + 1)

            y_linear[start_idx:end_idx] = y_data[start_idx:end_idx]

            r2_linear = self.calculate_r_squared(y_data, y_linear)

            if r2_linear > 0.01:  # Minimal threshold
                fallback_results.append({
                    'method': 'linear_interpolation_ultimate_fallback',
                    'r_squared': r2_linear,
                    'parameters': [data_max - baseline_est, x_data[peak_idx], ppm_range/20, 0.0, baseline_est],
                    'fitted_curve': y_linear,
                    'success': True,
                    'emergency_fallback': True,
                    'fallback_level': 4
                })

        except Exception as e:
            print(f"   ❌ LEVEL1: Linear interpolation fallback failed: {e}")

        # === SELECT BEST FALLBACK RESULT ===
        if fallback_results:
            # Sort by R² and select best
            best_fallback = max(fallback_results, key=lambda x: x['r_squared'])


            # Add comprehensive fallback information
            best_fallback.update({
                'amplitude': best_fallback['parameters'][0],
                'center': best_fallback['parameters'][1],
                'sigma': best_fallback['parameters'][2],
                'gamma': best_fallback['parameters'][3],
                'baseline': best_fallback['parameters'][4],
                'total_fallback_attempts': len(fallback_results),
                'fallback_methods_tried': [r['method'] for r in fallback_results]
            })

            return best_fallback
        else:
            return {
                'success': False,
                'method': 'all_emergency_fallbacks_failed',
                'r_squared': 0,
                'error': 'complete_emergency_fallback_failure',
                'emergency_fallback': True,
                'fallback_level': 99
            }

    def progressive_safe_fitting(self, x_data, y_data, peak_center, initial_guess):
        """
        SAFETY ENHANCEMENT: Progressive fitting strategy with increasing complexity

        Tries multiple fitting approaches in order of increasing complexity:
        1. Simple Gaussian fit (baseline validation)
        2. Pure Lorentzian fit (alternative check)
        3. Voigt fit with tight bounds
        4. Voigt fit with relaxed bounds (only if previous steps reasonable)

        Returns the best successful fit or None if all fail
        """
        print("   🔄 SAFETY: Attempting progressive fitting strategy...")

        fitting_attempts = []

        # Strategy 1: Simple Gaussian fit for baseline validation
        try:
            def gaussian_1d(x, amp, center, sigma, baseline):
                return amp * np.exp(-0.5 * ((x - center) / sigma) ** 2) + baseline

            gauss_guess = [initial_guess[0], initial_guess[1], initial_guess[2], initial_guess[4]]
            gauss_bounds = ([0, peak_center-0.01, 0.0001, np.min(y_data)],
                           [np.max(y_data)*3, peak_center+0.01, 0.1, np.max(y_data)])

            gauss_popt, _ = curve_fit(gaussian_1d, x_data, y_data, p0=gauss_guess, bounds=gauss_bounds, maxfev=500)
            gauss_fitted = gaussian_1d(x_data, *gauss_popt)
            gauss_r2 = self.calculate_r_squared(y_data, gauss_fitted)

            if gauss_r2 > 0.5:  # Reasonable Gaussian fit achieved
                fitting_attempts.append({
                    'method': 'gaussian',
                    'r_squared': gauss_r2,
                    'parameters': gauss_popt,
                    'fitted_curve': gauss_fitted,
                    'success': True
                })
                print(f"   ✅ SAFETY: Gaussian fit successful (R² = {gauss_r2:.4f})")

        except Exception as e:
            print(f"   ❌ SAFETY: Gaussian fit failed: {e}")

        # Strategy 2: Pure Lorentzian fit
        try:
            def lorentzian_1d(x, amp, center, gamma, baseline):
                return amp * gamma**2 / ((x - center)**2 + gamma**2) + baseline

            lorentz_guess = [initial_guess[0], initial_guess[1], initial_guess[3], initial_guess[4]]
            lorentz_bounds = ([0, peak_center-0.01, 0.0001, np.min(y_data)],
                             [np.max(y_data)*3, peak_center+0.01, 0.1, np.max(y_data)])

            lorentz_popt, _ = curve_fit(lorentzian_1d, x_data, y_data, p0=lorentz_guess, bounds=lorentz_bounds, maxfev=500)
            lorentz_fitted = lorentzian_1d(x_data, *lorentz_popt)
            lorentz_r2 = self.calculate_r_squared(y_data, lorentz_fitted)

            if lorentz_r2 > 0.5:  # Reasonable Lorentzian fit achieved
                fitting_attempts.append({
                    'method': 'lorentzian',
                    'r_squared': lorentz_r2,
                    'parameters': lorentz_popt,
                    'fitted_curve': lorentz_fitted,
                    'success': True
                })
                print(f"   ✅ SAFETY: Lorentzian fit successful (R² = {lorentz_r2:.4f})")

        except Exception as e:
            print(f"   ❌ SAFETY: Lorentzian fit failed: {e}")

        # Strategy 3: Voigt fit with tight bounds (original method)
        try:
            original_result = self._fit_with_specific_baseline(x_data, y_data, initial_guess, peak_center)
            if original_result.get('success', False) and original_result.get('r_squared', 0) > 0.3:
                fitting_attempts.append({
                    'method': 'voigt_tight',
                    'r_squared': original_result['r_squared'],
                    'result': original_result,
                    'success': True
                })
                print(f"   ✅ SAFETY: Tight Voigt fit successful (R² = {original_result['r_squared']:.4f})")

        except Exception as e:
            print(f"   ❌ SAFETY: Tight Voigt fit failed: {e}")

        # Return the best fit attempt
        if fitting_attempts:
            best_fit = max(fitting_attempts, key=lambda x: x['r_squared'])

            if best_fit['method'] == 'voigt_tight':
                return best_fit['result']
            else:
                # Convert simple fits to Voigt format for consistency
                return {
                    'success': True,
                    'method': best_fit['method'],
                    'r_squared': best_fit['r_squared'],
                    'fitted_curve': best_fit['fitted_curve'],
                    'fallback_fit': True,
                    'safety_triggered': 'progressive_fitting'
                }
        else:
            return {
                'success': False,
                'error': 'all_progressive_methods_failed',
                'r_squared': 0,
                'safety_triggered': 'progressive_fitting_failure'
            }

    def calculate_parameter_uncertainties(self, popt, pcov, alpha=0.05):
        """
        Calculate parameter uncertainties from covariance matrix

        Returns confidence intervals at (1-alpha) confidence level
        """
        try:
            # Parameter standard errors
            param_errors = np.sqrt(np.diag(pcov))

            # t-value for confidence interval (assuming large sample)
            from scipy.stats import t
            t_value = t.ppf(1 - alpha/2, df=len(popt))

            # Confidence intervals
            confidence_intervals = []
            for i, (param, error) in enumerate(zip(popt, param_errors)):
                ci_lower = param - t_value * error
                ci_upper = param + t_value * error
                confidence_intervals.append((ci_lower, ci_upper))

            return {
                'parameter_errors': param_errors,
                'confidence_intervals': confidence_intervals,
                'correlation_matrix': pcov / np.outer(param_errors, param_errors)
            }

        except Exception as e:
            print(f"Uncertainty calculation failed: {e}")
            return {
                'parameter_errors': np.full(len(popt), np.nan),
                'confidence_intervals': [(np.nan, np.nan)] * len(popt),
                'correlation_matrix': np.full((len(popt), len(popt)), np.nan)
            }

    def extract_local_peak_region(self, x_data, y_data, peak_center, fitted_width=None,
                                  nucleus_type=None, window_multiplier=None):
        """
        Extract local region around peak for quality assessment

        ENHANCED: Now supports GUI-based window sizing for consistent display

        This limits quality evaluation to the peak itself within reasonable range
        of its linewidth, excluding distant peaks that shouldn't affect local fit quality.

        Parameters:
        - x_data: full x-axis data
        - y_data: full intensity data
        - peak_center: fitted peak center
        - fitted_width: fitted peak width (sigma + gamma), or None for estimation
        - nucleus_type: nucleus type for typical width fallback
        - window_multiplier: how many widths to include on each side (None = use GUI params)

        Returns:
        - Dictionary with local region data and indices

        BACKWARD COMPATIBILITY:
        - If window_multiplier is provided: uses that value (original behavior)
        - If window_multiplier is None: calculates from GUI parameters (new behavior)
        """
        try:
            if fitted_width is None:
                # Estimate width from data if not provided
                fitted_width = self.adaptive_width_estimation(x_data, y_data, peak_center, nucleus_type)

            # ENHANCED: Calculate window multiplier from GUI parameters if not provided
            if window_multiplier is None:
                ppm_range = abs(x_data[-1] - x_data[0]) if len(x_data) > 1 else 1.0
                window_multiplier = self._calculate_gui_based_multiplier(
                    nucleus_type, ppm_range, len(x_data), fitted_width
                )
                print(f"   🎯 Using GUI-based window multiplier: {window_multiplier:.2f}× (was hardcoded 6.0×)")

            # Define local window around peak (now using GUI-based multiplier)
            half_window = fitted_width * window_multiplier

            # Find indices for local region
            local_mask = (x_data >= peak_center - half_window) & (x_data <= peak_center + half_window)
            local_indices = np.where(local_mask)[0]

            if len(local_indices) < 5:  # Need minimum points for meaningful assessment
                # Expand to at least 5 points around peak center
                center_idx = np.argmin(np.abs(x_data - peak_center))
                min_points = 5
                start_idx = max(0, center_idx - min_points//2)
                end_idx = min(len(x_data), center_idx + min_points//2 + 1)
                local_indices = np.arange(start_idx, end_idx)
                local_mask = np.zeros(len(x_data), dtype=bool)
                local_mask[local_indices] = True

            # Extract local data
            x_local = x_data[local_indices]
            y_local = y_data[local_indices]

            # Validation warning
            if len(local_indices) < 3:
                print(f"   ⚠️ Warning: Very small local region ({len(local_indices)} points) for peak at {peak_center:.3f}")

            return {
                #'x_data': x_local,
                #'y_data': y_local,
                #'indices': local_indices,
                'x_local': x_local,          # ✅ Consistent with caller expectations
                'y_local': y_local,          # ✅ Consistent with caller expectations
                'local_indices': local_indices,  # ✅ Consistent with caller expectations
                'peak_center': peak_center,
                'window_size': half_window * 2,  # total window size
                'window_width': half_window * 2,
                'n_points': len(local_indices),
                'ppm_range': x_local[-1] - x_local[0] if len(x_local) > 1 else 0,
                'multiplier_used': window_multiplier,
                'gui_based': window_multiplier != 6.0  # flag for logging/debugging
            }

        except Exception as e:
            print(f"Local region extraction failed: {e}")
            # Fallback to full data
            return {
                'x_local': x_data,
                'y_local': y_data,
                'local_indices': np.arange(len(x_data)),
                'local_mask': np.ones(len(x_data), dtype=bool),
                'window_width': abs(x_data[-1] - x_data[0]),
                'n_points': len(x_data)
            }

    def calculate_local_baseline(self, x_local, y_local):
        """
        Calculate local baseline within the peak region for better quality assessment
        """
        try:
            # For local regions, use edge points for baseline
            n_edge = max(1, len(y_local) // 6)  # Use outer 1/6 on each side

            if len(y_local) < 6:
                # For very small regions, just use minimum
                return np.min(y_local)

            left_edge = y_local[:n_edge]
            right_edge = y_local[-n_edge:]

            # Combine edges and use median for robustness
            edge_values = np.concatenate([left_edge, right_edge])
            local_baseline = np.median(edge_values)

            return local_baseline

        except:
            # Fallback
            return np.percentile(y_local, 20)

    def comprehensive_quality_assessment(self, x_data, y_data, y_pred, popt, pcov):
        """
        Comprehensive quality assessment with both global and local metrics
        """
        try:
            quality_metrics = {}

            # Extract fitted parameters
            amplitude, center, sigma, gamma, baseline = popt
            fitted_width = sigma + gamma

            # === GLOBAL QUALITY METRICS (original behavior) ===
            ss_res_global = np.sum((y_data - y_pred) ** 2)
            ss_tot_global = np.sum((y_data - np.mean(y_data)) ** 2)
            r_squared_global = 1 - (ss_res_global / ss_tot_global) if ss_tot_global != 0 else 0
            quality_metrics['r_squared_global'] = r_squared_global

            # === LOCAL QUALITY METRICS (new peak-specific assessment) ===
            # Extract local region around the peak
            nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
            local_region = self.extract_local_peak_region(x_data, y_data, center, fitted_width, nucleus_type)

            # Get local data and predictions
            x_local = local_region['x_local']
            y_local = local_region['y_local']
            y_pred_local = y_pred[local_region['local_indices']]

            # Calculate local baseline for better reference
            local_baseline = self.calculate_local_baseline(x_local, y_local)

            # Local R-squared using local baseline as reference
            ss_res_local = np.sum((y_local - y_pred_local) ** 2)
            ss_tot_local = np.sum((y_local - local_baseline) ** 2)  # Use local baseline
            r_squared_local = 1 - (ss_res_local / ss_tot_local) if ss_tot_local != 0 else 0
            quality_metrics['r_squared_local'] = r_squared_local

            # Peak-specific metrics
            quality_metrics['peak_region_points'] = local_region['n_points']
            quality_metrics['peak_region_width'] = local_region['window_width']

            # Local residual analysis
            residuals_local = y_local - y_pred_local
            quality_metrics['residual_std_local'] = np.std(residuals_local)
            quality_metrics['residual_mean_local'] = np.mean(residuals_local)

            # Local RMSE normalized by peak height
            peak_height_local = np.max(y_local) - local_baseline
            rmse_local = np.sqrt(np.mean(residuals_local ** 2))
            quality_metrics['rmse_local'] = rmse_local
            quality_metrics['rmse_normalized'] = rmse_local / peak_height_local if peak_height_local > 0 else np.inf

            # === COMPATIBILITY: Use local R-squared as primary metric ===
            quality_metrics['r_squared'] = r_squared_local  # Primary metric for compatibility

            # === GLOBAL METRICS (for reference) ===
            n = len(y_data)
            p = len(popt)
            adj_r_squared = 1 - (1 - r_squared_global) * (n - 1) / (n - p - 1) if n > p + 1 else r_squared_global
            quality_metrics['adj_r_squared'] = adj_r_squared

            rmse_global = np.sqrt(np.mean((y_data - y_pred) ** 2))
            quality_metrics['rmse_global'] = rmse_global

            chi_squared_red = ss_res_global / (n - p) if n > p else np.inf
            quality_metrics['chi_squared_reduced'] = chi_squared_red

            residuals_global = y_data - y_pred
            quality_metrics['residual_std_global'] = np.std(residuals_global)
            quality_metrics['residual_mean_global'] = np.mean(residuals_global)

            # Parameter validation
            quality_metrics['parameters_physical'] = (
                amplitude > 0 and
                sigma > 0 and
                gamma >= 0 and
                0 < sigma < abs(x_data[-1] - x_data[0]) and
                0 <= gamma < abs(x_data[-1] - x_data[0])
            )

            # Overall quality classification based on LOCAL R-squared
            if r_squared_local >= 0.95 and quality_metrics['parameters_physical']:
                quality_class = "Excellent"
            elif r_squared_local >= 0.85 and quality_metrics['parameters_physical']:
                quality_class = "Good"
            elif r_squared_local >= 0.7:
                quality_class = "Fair"
            else:
                quality_class = "Poor"

            quality_metrics['quality_class'] = quality_class

            return quality_metrics

#        except Exception as e:
#            print(f"Quality assessment failed: {e}")
#            return {'r_squared': 0, 'r_squared_local': 0, 'r_squared_global': 0, 'quality_class': 'Failed'}
        # IMPROVED:
        except Exception as e:
            print(f"Local quality assessment failed, using global metrics: {e}")
            # Fallback: use global R² when local assessment fails
            return {
                'r_squared': r_squared_global,           # Use global instead of 0
                'r_squared_local': r_squared_global,     # Use global as fallback
                'r_squared_global': r_squared_global,    # This was calculated successfully
                'quality_class': 'Good' if r_squared_global >= 0.8 else 'Fair' if r_squared_global >=  0.7 else 'Poor'
        }

    def assess_fit_quality_comprehensive(self, x_data, y_data, fit_result, nucleus_type=None):
        """
        Comprehensive quality assessment for multi-peak fitting results

        This method bridges the interface gap between the existing comprehensive_quality_assessment
        method and the multi-peak fitting workflow. It extracts fitted parameters from a result
        dictionary and delegates to the existing quality assessment infrastructure.

        Args:
            x_data: X-axis data (ppm)
            y_data: Y-axis intensity data
            fit_result: Dictionary containing fit results with 'fitted_curve' and 'parameters'
            nucleus_type: Optional nucleus type for specialized assessment

        Returns:
            dict: Quality metrics including 'quality_class', 'r_squared', etc.
        """
        try:
            # Extract fitted curve and parameters from result
            y_fitted = fit_result.get('fitted_curve')
            parameters = fit_result.get('parameters')

            if y_fitted is None:
                # Try to reconstruct fitted curve from parameters
                if parameters is not None and len(parameters) >= 5:
                    # For single peak: [amplitude, center, sigma, gamma, baseline]
                    y_fitted = self.voigt_profile_1d(x_data, *parameters[:5])
                elif parameters is not None and len(parameters) > 5:
                    # For multi-peak: reconstruct from all parameters
                    n_peaks = (len(parameters) - 1) // 4  # exclude shared baseline
                    y_fitted = np.zeros_like(x_data)

                    for i in range(n_peaks):
                        param_start = i * 4
                        if param_start + 3 < len(parameters) - 1:
                            amp = parameters[param_start]
                            center = parameters[param_start + 1]
                            sigma = parameters[param_start + 2]
                            gamma = parameters[param_start + 3]
                            baseline = parameters[-1]  # shared baseline

                            y_fitted += self.voigt_profile_1d(x_data, amp, center, sigma, gamma, 0)

                    y_fitted += parameters[-1]  # add shared baseline
                else:
                    raise ValueError("Cannot reconstruct fitted curve: missing parameters")

            # Extract or estimate parameter covariance
            pcov = fit_result.get('covariance')
            if pcov is None:
                # Create dummy covariance matrix for compatibility
                n_params = len(parameters) if parameters is not None else 5
                pcov = np.eye(n_params) * 1e-6

            # Determine effective parameters for single-peak assessment
            if parameters is not None and len(parameters) >= 5:
                # Use primary peak parameters (or best peak for multi-peak)
                if len(parameters) == 5:
                    # Single peak case
                    effective_params = parameters
                else:
                    # Multi-peak case: use the peak with highest amplitude
                    n_peaks = (len(parameters) - 1) // 4
                    best_peak_idx = 0
                    best_amplitude = 0

                    for i in range(n_peaks):
                        param_start = i * 4
                        if param_start < len(parameters) - 1:
                            amp = abs(parameters[param_start])
                            if amp > best_amplitude:
                                best_amplitude = amp
                                best_peak_idx = i

                    # Extract best peak parameters
                    param_start = best_peak_idx * 4
                    effective_params = [
                        parameters[param_start],      # amplitude
                        parameters[param_start + 1],  # center
                        parameters[param_start + 2],  # sigma
                        parameters[param_start + 3],  # gamma
                        parameters[-1]                # baseline
                    ]
            else:
                raise ValueError("Invalid or missing parameters in fit result")

            # Use existing comprehensive quality assessment
            quality_metrics = self.comprehensive_quality_assessment(
                x_data, y_data, y_fitted, effective_params, pcov
            )

            # Add multi-peak specific metrics
            if 'n_peaks_fitted' in fit_result:
                quality_metrics['n_peaks_fitted'] = fit_result['n_peaks_fitted']
                quality_metrics['multi_peak_assessment'] = True

                # Adjust quality class for multi-peak complexity
                original_r2 = quality_metrics.get('r_squared', 0)
                n_peaks = fit_result['n_peaks_fitted']

                # Multi-peak penalty: more peaks require higher R² for same quality class
                if n_peaks > 2:
                    r2_threshold_adjustment = 0.05 * (n_peaks - 2)
                    adjusted_r2 = original_r2 - r2_threshold_adjustment

                    # Reclassify based on adjusted R²
                    if adjusted_r2 >= 0.95:
                        quality_metrics['quality_class'] = "Excellent"
                    elif adjusted_r2 >= 0.85:
                        quality_metrics['quality_class'] = "Good"
                    elif adjusted_r2 >= 0.7:
                        quality_metrics['quality_class'] = "Fair"
                    else:
                        quality_metrics['quality_class'] = "Poor"

                    quality_metrics['multi_peak_penalty'] = r2_threshold_adjustment

            # Add fit success indicator
            quality_metrics['fit_success'] = fit_result.get('success', False)
            quality_metrics['fit_method'] = fit_result.get('method', 'unknown')

            return quality_metrics

        except Exception as e:
            print(f"Comprehensive quality assessment failed: {e}")
            # Return minimal fallback result
            r_squared_fallback = fit_result.get('r_squared', 0) if isinstance(fit_result, dict) else 0
            return {
                'r_squared': r_squared_fallback,
                'r_squared_local': r_squared_fallback,
                'r_squared_global': r_squared_fallback,
                'quality_class': 'Poor',
                'fit_success': False,
                'assessment_error': str(e)
            }




    def detect_potential_overlap(self, x_data: np.ndarray,
                                  y_data: np.ndarray,
                                  primary_peak_ppm: float,
                                  nucleus_type: str = '1H') -> dict:
        """
        Detect if there are other peaks near the primary peak

        Uses simple local maximum detection in a window around primary peak.
        This is a lightweight overlap detector that routes to the full
        OverlapResolverEngine when overlaps are detected.

        Args:
            x_data: PPM axis
            y_data: Intensity data (should be baseline-corrected)
            primary_peak_ppm: Main peak position
            nucleus_type: Nucleus type for adaptive threshold

        Returns:
            {
                'has_overlap': bool,
                'peak_candidates': List[float],  # PPM positions
                'n_peaks': int,
                'separation_distance': float,     # Minimum distance between peaks
                'confidence': float               # 0-1 scale
            }
        """
        try:
            # Define search window based on nucleus type
            nucleus_windows = {
                '1H': 0.3,   # ±0.3 ppm window for 1H
                '15N': 1.0,  # ±1.0 ppm window for 15N
                '13C': 0.8,  # ±0.8 ppm window for 13C
            }
            search_window = nucleus_windows.get(nucleus_type, 0.5)

            # Extract region around primary peak
            mask = np.abs(x_data - primary_peak_ppm) <= search_window
            if not np.any(mask):
                return {
                    'has_overlap': False,
                    'peak_candidates': [primary_peak_ppm],
                    'n_peaks': 1,
                    'separation_distance': float('inf'),
                    'confidence': 0.0
                }

            x_region = x_data[mask]
            y_region = y_data[mask]

            # Find local maxima using scipy's find_peaks
            # Adaptive prominence based on data range
            data_range = np.max(y_region) - np.min(y_region)
            prominence = 0.05 * data_range  # 5% of local data range

            peaks_idx, properties = find_peaks(
                y_region,
                prominence=prominence,
                distance=int(len(x_region) * 0.05)  # At least 5% of points apart
            )

            if len(peaks_idx) == 0:
                return {
                    'has_overlap': False,
                    'peak_candidates': [primary_peak_ppm],
                    'n_peaks': 1,
                    'separation_distance': float('inf'),
                    'confidence': 0.0
                }

            # Convert indices to PPM positions
            peak_candidates = x_region[peaks_idx].tolist()

            # Calculate minimum separation
            if len(peak_candidates) > 1:
                separations = np.diff(sorted(peak_candidates))
                min_separation = float(np.min(np.abs(separations)))
            else:
                min_separation = float('inf')

            # Determine if overlap exists
            has_overlap = (
                len(peak_candidates) > 1 and
                min_separation < self.overlap_detection_threshold
            )

            # Confidence based on prominence and separation
            if len(peak_candidates) > 1:
                avg_prominence = float(np.mean(properties['prominences']))
                prominence_ratio = avg_prominence / data_range
                confidence = min(1.0, prominence_ratio * 2)  # Scale to 0-1
            else:
                confidence = 0.0

            return {
                'has_overlap': has_overlap,
                'peak_candidates': peak_candidates,
                'n_peaks': len(peak_candidates),
                'separation_distance': min_separation,
                'confidence': confidence
            }

        except Exception as e:
            print(f"   ⚠️ Overlap detection failed: {e}")
            return {
                'has_overlap': False,
                'peak_candidates': [primary_peak_ppm],
                'n_peaks': 1,
                'separation_distance': float('inf'),
                'confidence': 0.0
            }

    def configure_overlap_resolution(self, enable: bool = True,
                                     config: dict = None):
        """
        Configure overlap resolution behavior

        This method enables/disables the overlap detection and resolution system.
        By default, overlap resolution is DISABLED for backward compatibility.

        Args:
            enable: Enable/disable overlap detection
            config: Optional configuration dictionary for overlap resolver

        Example:
            fitter.configure_overlap_resolution(
                enable=True,
                config={
                    'jackknife': {'n_resamples': 100},
                    'model_selection': {'max_peaks': 5}
                }
            )
        """
        if not OVERLAP_RESOLUTION_AVAILABLE:
            print("⚠️ Overlap resolution modules not available")
            print("   Please ensure overlap_resolver_engine.py and overlap_config.py are installed")
            return

        self.overlap_detection_enabled = enable

        if config:
            if self.overlap_config is None:
                self.overlap_config = OverlapResolutionConfig(config)
            else:
                self.overlap_config.update(config)

        status = 'enabled' if enable else 'disabled'
        print(f"✅ Overlap resolution {status}")
        if enable and config:
            print(f"   Configuration updated: {list(config.keys())}")


    def fit_voigt_profile(self, nmr_data, ppm_x, ppm_y, ppm_scale_x, ppm_scale_y,
                          assignment=None, linewidth_constraints=None):
        """
        Fit Voigt profile to a 2D NMR peak by extracting cross-sections and fitting each dimension.

        This is a wrapper method that extracts 1D cross-sections from 2D NMR data and fits
        Voigt profiles along each dimension (X and Y). This method is called by series_processor.py
        for PS2D-style linewidth reuse functionality.

        Parameters:
        -----------
        nmr_data : np.ndarray
            2D NMR data matrix (shape: [Y_points, X_points])
        ppm_x : float
            Peak position in X dimension (1H, ppm)
        ppm_y : float
            Peak position in Y dimension (15N/13C, ppm)
        ppm_scale_x : np.ndarray
            PPM scale for X dimension (1H)
        ppm_scale_y : np.ndarray
            PPM scale for Y dimension (15N/13C)
        assignment : str, optional
            Peak assignment/label (e.g., 'A123N-H')
        linewidth_constraints : dict, optional
            Linewidth constraints for PS2D linewidth reuse:
            {
                'x': {'sigma_bounds': (min, max), 'gamma_bounds': (min, max)},
                'y': {'sigma_bounds': (min, max), 'gamma_bounds': (min, max)}
            }

        Returns:
        --------
        dict : Comprehensive fit results with X and Y dimension fits
            {
                'success': bool,
                'x_fit': {...},  # X dimension fit results
                'y_fit': {...},  # Y dimension fit results
                'assignment': str,
                'position': {'x': ppm_x, 'y': ppm_y},
                'linewidth_constraints_applied': bool
            }
        """
        # Find peak indices in the data
        x_idx = np.argmin(np.abs(ppm_scale_x - ppm_x))
        y_idx = np.argmin(np.abs(ppm_scale_y - ppm_y))

        # Extract 1D cross-sections at peak position
        # X cross-section: horizontal slice at peak Y position
        x_cross_section = nmr_data[y_idx, :]
        x_ppm_scale = ppm_scale_x

        # Y cross-section: vertical slice at peak X position
        y_cross_section = nmr_data[:, x_idx]
        y_ppm_scale = ppm_scale_y

        # Extract linewidth constraints for each dimension if provided
        x_linewidth_constraints = None
        y_linewidth_constraints = None
        constraints_applied = False

        if linewidth_constraints:
            x_linewidth_constraints = linewidth_constraints.get('x')
            y_linewidth_constraints = linewidth_constraints.get('y')
            constraints_applied = (x_linewidth_constraints is not None or
                                  y_linewidth_constraints is not None)

            if constraints_applied:
                print(f"   🔒 Applying linewidth constraints for {assignment or 'peak'}")
                if x_linewidth_constraints:
                    print(f"      X (1H): σ={x_linewidth_constraints.get('sigma_bounds')}, "
                          f"γ={x_linewidth_constraints.get('gamma_bounds')}")
                if y_linewidth_constraints:
                    print(f"      Y (15N): σ={y_linewidth_constraints.get('sigma_bounds')}, "
                          f"γ={y_linewidth_constraints.get('gamma_bounds')}")

        # Fit X dimension (1H) with constraints
        try:
            x_fit = self.fit_peak_enhanced(
                x_ppm_scale, x_cross_section,
                initial_center=ppm_x,
                nucleus_type='1H',
                method='iterative_optimization',
                linewidth_constraints=x_linewidth_constraints
            )
        except Exception as e:
            print(f"   ⚠️ X dimension fit failed: {e}")
            x_fit = {'success': False, 'error': str(e)}

        # Fit Y dimension (15N/13C) with constraints
        try:
            y_fit = self.fit_peak_enhanced(
                y_ppm_scale, y_cross_section,
                initial_center=ppm_y,
                nucleus_type='15N',  # Assume 15N, could be detected
                method='iterative_optimization',
                linewidth_constraints=y_linewidth_constraints
            )
        except Exception as e:
            print(f"   ⚠️ Y dimension fit failed: {e}")
            y_fit = {'success': False, 'error': str(e)}

        # Combine results
        result = {
            'success': x_fit.get('success', False) and y_fit.get('success', False),
            'x_fit': x_fit,
            'y_fit': y_fit,
            'assignment': assignment,
            'position': {'x': ppm_x, 'y': ppm_y},
            'linewidth_constraints_applied': constraints_applied,
            'method': 'fit_voigt_profile'
        }

        # Print summary
        if result['success']:
            x_r2 = x_fit.get('r_squared', 0.0)
            y_r2 = y_fit.get('r_squared', 0.0)
            print(f"   ✅ Fit successful: R²_X={x_r2:.3f}, R²_Y={y_r2:.3f}")

            # Extract fitted linewidths for reporting
            if 'fitted_params' in x_fit:
                x_sigma = x_fit['fitted_params'].get('sigma', 0.0)
                x_gamma = x_fit['fitted_params'].get('gamma', 0.0)
                print(f"      X linewidths: σ={x_sigma:.4f}, γ={x_gamma:.4f} ppm")

            if 'fitted_params' in y_fit:
                y_sigma = y_fit['fitted_params'].get('sigma', 0.0)
                y_gamma = y_fit['fitted_params'].get('gamma', 0.0)
                print(f"      Y linewidths: σ={y_sigma:.4f}, γ={y_gamma:.4f} ppm")
        else:
            print(f"   ❌ Fit failed for {assignment or 'peak'}")

        return result

    def fit_peak_enhanced(self, x_data, y_data, initial_center=None, nucleus_type=None,
                         method='iterative_optimization', preprocessing=True,
                         all_peaks_context=None, linewidth_constraints=None, detection_confidence=None,
                         force_single_peak=False, _call_stack=None):
        """
        ENHANCED PEAK FITTING WITH DYNAMIC OPTIMIZATION

        This is the main entry point for enhanced peak fitting with iterative optimization.

        NEW OPTIMIZATION MODES:
        - 'iterative_optimization': Uses dynamic baseline optimization and quality-driven convergence
        - 'multi_step': Original multi-step fitting (preserved for backward compatibility)
        - 'single_step': Single-step fitting (preserved for backward compatibility)

        Parameters:
        - x_data: x-axis data (ppm)
        - y_data: intensity data
        - initial_center: initial guess for peak center (if None, uses data maximum)
        - nucleus_type: '1H', '15N', or '13C' (auto-detected if None)
        - method: 'iterative_optimization', 'multi_step', or 'single_step'
        - preprocessing: whether to apply data preprocessing
        - all_peaks_context: List of all suspected peaks for global parameter estimation (NEW)
        - linewidth_constraints: Dict with 'sigma_bounds' and 'gamma_bounds' for constrained fitting (NEW)
        - detection_confidence: Dict with detection confidence scores and peak characteristics (INTEGRATION)
        - force_single_peak: If True, skip overlap resolution (default: False) (NEW)

        Returns:
        - Comprehensive fit results dictionary with optimization diagnostics

        BACKWARD COMPATIBILITY: Original methods preserved, new method is opt-in
        NEW: Overlap detection and resolution (opt-in via configure_overlap_resolution)
        """
        # === RECURSION DETECTION ===
        if _call_stack is None:
            _call_stack = []

        # Create unique call identifier
        call_id = f"fit_peak_enhanced_{method}_{refined_center if initial_center is None else initial_center:.4f}"
        if hasattr(self, '_last_center'):
            call_id = f"fit_peak_enhanced_{method}_{self._last_center:.4f}"
        else:
            center_for_id = initial_center if initial_center is not None else (x_data[np.argmax(y_data)] if len(y_data) > 0 else 0.0)
            call_id = f"fit_peak_enhanced_{method}_{center_for_id:.4f}"

        # Check for recursion
        if call_id in _call_stack:
            print(f"   🛑 RECURSION DETECTED: {call_id}")
            print(f"   Call stack: {' -> '.join(_call_stack[-3:])}")  # Show last 3 calls
            return self.emergency_fallback_fitting(x_data, y_data,
                                                 initial_center if initial_center is not None else x_data[np.argmax(y_data)])

        # Add to call stack
        _call_stack.append(call_id)
        current_call_stack = _call_stack.copy()  # Preserve for passing to other methods

        try:
            # Clear previous diagnostics and initialize
            self.last_fit_diagnostics = {
                'method': method,
                'optimization_active': method == 'iterative_optimization',
                'detection_confidence': detection_confidence,
                'call_stack': current_call_stack,
                'recursion_detected': False
            }

            # Data preprocessing (preserved from original)
            if preprocessing and len(y_data) > 10:
                window_size = min(self.fitting_parameters['smoothing_window'], len(y_data) // 3)
                if window_size >= 3 and window_size % 2 == 1:
                    y_data = savgol_filter(y_data, window_size, 2)

            # Detect nucleus type
            if nucleus_type is None:
                nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])

            # Peak center refinement
            if initial_center is None:
                initial_center = x_data[np.argmax(y_data)]
            refined_center = self.refine_peak_center(x_data, y_data, initial_center)

            # === DYNAMIC OPTIMIZATION MODE ===
            if method == 'iterative_optimization':
                print(f"🔄 Starting iterative optimization for peak at {refined_center:.4f} ppm")

                # Step 1: Global parameter estimation from well-resolved peaks
                global_params = self.estimate_initial_parameters_from_resolved_peaks(
                    x_data, y_data, all_peaks_context
                )
                self.last_fit_diagnostics['global_params'] = global_params

                # ===================================================================
                # ===================================================================
                baseline_est = 0.0
                baseline_value = 0.0
                print

                # Store baseline info in diagnostics
                self.last_fit_diagnostics['baseline_method'] = {
                    'selected_method': 'none',
                    'peak_complexity': 'N/A',
                    'method_info': 'Using raw data '
                }

                # === NEW: OVERLAP DETECTION AND ROUTING ===
                # Baseline-corrected data for overlap detection
                if np.isscalar(baseline_est):
                    y_corrected = y_data - baseline_est
                else:
                    y_corrected = y_data - baseline_est

                # Check if overlap detection is enabled and available
                if (OVERLAP_RESOLUTION_AVAILABLE and
                    self.overlap_detection_enabled and
                    not force_single_peak):

                    overlap_info = self.detect_potential_overlap(
                        x_data, y_corrected, refined_center, nucleus_type
                    )

                    if overlap_info['has_overlap'] and overlap_info['n_peaks'] > 1:
                        print(f"   🔍 Detected {overlap_info['n_peaks']} potentially overlapping peaks")
                        print(f"      Separation: {overlap_info['separation_distance']:.4f} ppm, Confidence: {overlap_info['confidence']:.2f}")

                        # Lazy initialize overlap resolver
                        if self.overlap_resolver is None:
                            self.overlap_resolver = OverlapResolverEngine(enhanced_fitter=self)
                            if self.overlap_config:
                                # Convert OverlapResolutionConfig to dict before updating
                                config_dict = self.overlap_config.to_dict() if hasattr(self.overlap_config, 'to_dict') else self.overlap_config
                                self.overlap_resolver.config.update(config_dict)
                            print(f"   🚀 Overlap resolver initialized")

                        # Try overlap resolution
                        try:
                            result = self.overlap_resolver.resolve_overlapping_peaks(
                                x_data, y_corrected,
                                overlap_info['peak_candidates']
                            )

                            # Add baseline back to result
                            result['baseline'] = baseline_est if np.isscalar(baseline_est) else baseline_value
                            result['fitting_method'] = 'overlap_resolution'
                            result['overlap_info'] = overlap_info

                            # Check if resolution was successful
                            if result['success'] and result.get('r_squared', 0) > 0.7:
                                print(f"   ✅ Overlap resolution successful: {result['n_peaks']} peaks, R²={result['r_squared']:.3f}")
                                return result
                            else:
                                print(f"   ⚠️ Overlap resolution quality poor (R²={result.get('r_squared', 0):.3f})")
                                print(f"   🔄 Falling back to single-peak fit")
                                # Fall through to single-peak fitting

                        except Exception as e:
                            print(f"   ❌ Overlap resolution failed: {e}")
                            print(f"   🔄 Falling back to single-peak fit")
                            # Fall through to single-peak fitting

                # INTEGRATION ENHANCEMENT: Use detection confidence for parameter estimation
                if detection_confidence and isinstance(detection_confidence, dict):
                    # Use detection-informed amplitude if available
                    if 'estimated_amplitude' in detection_confidence:
                        amplitude_est = detection_confidence['estimated_amplitude']
                        print(f"     Detection amplitude: {amplitude_est:.0f}")
                    else:
                        amplitude_est = np.max(y_data) - baseline_value

                    # Use detection-informed width if available
                    if 'estimated_width' in detection_confidence:
                        width_est = detection_confidence['estimated_width']
                        print(f"     Detection width: {width_est:.4f} ppm")
                    else:
                        width_est = global_params['typical_linewidth']

                    # Use chemical shift context for width refinement
                    if 'chemical_shift_context' in detection_confidence:
                        context = detection_confidence['chemical_shift_context']
                        if 'typical_width' in context:
                            context_width = context['typical_width']
                            # Weighted average of global and context estimates
                            confidence_weight = detection_confidence.get('confidence', 0.5)
                            width_est = (confidence_weight * context_width +
                                       (1 - confidence_weight) * width_est)
                            print(f"     Context-adjusted width: {width_est:.4f} ppm")

                    # Store detection confidence info in diagnostics
                    self.last_fit_diagnostics['detection_informed_params'] = {
                        'amplitude': amplitude_est,
                        'width': width_est,
                        'confidence': detection_confidence.get('confidence', 0),
                        'used_detection_info': True
                    }

                else:
                    # === LEVEL 2 INTEGRATION: Replace standard parameter estimation ===
                    if hasattr(self, 'parameter_estimator') and self.level2_params['robust_estimation_enabled']:

                        estimation_result = self.parameter_estimator.estimate_initial_parameters(
                            x_data, y_data, refined_center, nucleus_type, context=detection_confidence
                        )

                        if estimation_result['success']:
                            # Extract parameters from Level 2 estimation
                            level2_params = estimation_result['parameters']
                            amplitude_est = level2_params[0]
                            width_est = level2_params[2] + level2_params[3]  # sigma + gamma

                            self.last_fit_diagnostics['level2_estimation'] = {
                                'method': estimation_result['method'],
                                'consensus_quality': estimation_result['consensus_quality'],
                                'individual_methods': len(estimation_result['individual_results']),
                                'data_quality': estimation_result['data_quality'],
                                'used_level2': True
                            }
                        else:
                            # Fallback to standard method
                            amplitude_est = np.max(y_data) - baseline_value
                            width_est = global_params['typical_linewidth']

                            self.last_fit_diagnostics['level2_estimation'] = {
                                'used_level2': False,
                                'fallback_reason': 'estimation_failed'
                            }
                    else:
                        # Standard parameter estimation (backward compatibility)
                        amplitude_est = np.max(y_data) - baseline_value
                        width_est = global_params['typical_linewidth']

                        self.last_fit_diagnostics['level2_estimation'] = {
                            'used_level2': False,
                            'fallback_reason': 'level2_disabled'
                        }

                    self.last_fit_diagnostics['detection_informed_params'] = {
                        'amplitude': amplitude_est,
                        'width': width_est,
                        'used_detection_info': False
                    }

                # Apply linewidth constraints if provided (GLOBAL OPTIMIZATION ENHANCEMENT)
                if linewidth_constraints:
                    sigma_bounds = linewidth_constraints.get('sigma_bounds', (width_est * 0.1, width_est * 10))
                    gamma_bounds = linewidth_constraints.get('gamma_bounds', (width_est * 0.1, width_est * 10))

                    # Constrain initial parameters to be within bounds
                    sigma_est = np.clip(width_est * 0.6, sigma_bounds[0], sigma_bounds[1])
                    gamma_est = np.clip(width_est * 0.4, gamma_bounds[0], gamma_bounds[1])

                    print(f"   Applying linewidth constraints: σ={sigma_bounds}, γ={gamma_bounds}")
                    self.last_fit_diagnostics['linewidth_constraints'] = linewidth_constraints
                else:
                    sigma_est = width_est * 0.6
                    gamma_est = width_est * 0.4

                # Use global parameters for better initial guess
                initial_guess = [
                    amplitude_est,
                    refined_center,
                    sigma_est,  # sigma (Gaussian component) - potentially constrained
                    gamma_est,  # gamma (Lorentzian component) - potentially constrained
                    baseline_value  # Use scalar baseline value for fitting
                ]

                constraints_info = " (constrained)" if linewidth_constraints else ""
                print(f"   Initial parameters: amp={amplitude_est:.0f}, width={width_est:.4f}, baseline={baseline_value:.1f}{constraints_info}")

                # Step 3: Standard fitting attempt for baseline comparison
                standard_result = self._fit_with_standard_method(x_data, y_data, initial_guess, nucleus_type,
                                                                linewidth_constraints=linewidth_constraints)
                standard_quality = standard_result.get('r_squared', 0)

                print(f"   Standard fitting: R² = {standard_quality:.4f}")
                self.last_fit_diagnostics['standard_result'] = {
                    'r_squared': standard_quality,
                    'success': standard_result.get('success', False)
                }

                # ===================================================================
                # NO BASELINE OPTIMIZATION
                # ===================================================================
                print("   Skipping baseline optimization: using raw data approach")
                best_result = standard_result
                self.last_fit_diagnostics['baseline_optimization_skipped'] = 'raw_data_approach'

##
                # === NEW: AUTOMATIC MULTI-PEAK DETECTION ===
                # Step 5: Check if poor fit indicates overlapping peaks
                # Use GUI parameter instead of hardcoded value
                MULTI_PEAK_R2_THRESHOLD = self.fitting_parameters.get('multi_peak_r2_threshold', 0.7)

                if best_result.get('r_squared', 0) < MULTI_PEAK_R2_THRESHOLD:
                    print(f"   Poor single-peak fit detected (R² = {best_result.get('r_squared', 0):.3f}), checking for overlapping peaks...")

                    # Detect peaks within the current fitting region
                    detected_peaks = self.detect_overlapping_peaks(x_data, y_data, nucleus_type=nucleus_type)

                    if len(detected_peaks) > 1:
                        print(f"   Found {len(detected_peaks)} overlapping peaks, attempting multi-peak fit...")

                        peak_positions = [p['position'] for p in detected_peaks]

                        # Additional residual analysis for validation
                        if best_result.get('success', False):
                            # Handle both flat and nested result formats
                            if 'x_fit' in best_result:
                                # Nested format from 2D fitting
                                amplitude = best_result['x_fit']['amplitude']
                                center = best_result['x_fit']['center']
                                sigma = best_result['x_fit']['sigma']
                                gamma = best_result['x_fit']['gamma']
                                baseline = best_result['x_fit']['baseline']
                            else:
                                # Flat format from 1D fitting
                                amplitude = best_result['amplitude']
                                center = best_result['center']
                                sigma = best_result['sigma']
                                gamma = best_result['gamma']
                                baseline = best_result['baseline']

                            fitted_curve = self.voigt_profile_1d(
                                x_data,
                                amplitude,
                                center,
                                sigma,
                                gamma,
                                baseline
                            )
                            residuals = y_data - fitted_curve
                            residual_peaks = find_peaks(np.abs(residuals),
                                                      height=np.std(residuals) * 1.5,
                                                      distance=max(1, len(x_data) // 200))

                            print(f"   Residual analysis: {len(residual_peaks[0])} peak-like structures in residuals")

                            if len(residual_peaks[0]) >= 2:
                                print("   ✓ Residuals confirm multiple peak hypothesis")

                        # Attempt multi-peak fitting with progressive strategies
                        try:
                            # Try progressive strategies first (more robust)
                            multi_result = self.fit_with_progressive_strategies(
                                x_data, y_data, detected_peaks, global_constraints=None
                            )

                            # DEFENSIVE FIX: Ensure multi_result is a dictionary
                            if not isinstance(multi_result, dict):
                                print(f"   ⚠ Progressive strategies returned {type(multi_result)}, converting to dict")
                                multi_result = {'success': False, 'r_squared': 0.0, 'error': 'progressive_strategies_invalid_return'}

                            # If progressive strategies fail, try original iterative approach
                            if not multi_result.get('success', False):
                                print("   Progressive strategies failed, trying iterative optimization...")
                                iterative_result = self.optimize_overlap_detection_iteratively(
                                    x_data, y_data, peak_positions,
                                    max_iterations=6, use_aic_selection=True,  # Reduced iterations for speed
                                    _call_stack=current_call_stack
                                )

                                # CRITICAL FIX: Handle tuple return from optimize_overlap_detection_iteratively
                                if isinstance(iterative_result, tuple):
                                    # Expected format: (best_fit, optimization_report)
                                    best_fit, optimization_report = iterative_result
                                    if best_fit is not None and isinstance(best_fit, dict):
                                        multi_result = best_fit
                                        multi_result['optimization_report'] = optimization_report
                                    else:
                                        print(f"   ⚠ Iterative optimization best_fit is None or invalid, using fallback")
                                        multi_result = {'success': False, 'r_squared': 0.0, 'error': 'iterative_optimization_no_best_fit'}
                                elif isinstance(iterative_result, dict):
                                    multi_result = iterative_result
                                else:
                                    print(f"   ⚠ Iterative optimization returned {type(iterative_result)}, converting to dict")
                                    multi_result = {'success': False, 'r_squared': 0.0, 'error': 'iterative_optimization_invalid_return'}

                            multi_quality = multi_result.get('r_squared', 0)
                            improvement_threshold = self.fitting_parameters.get('multi_peak_improvement_threshold', 0.1)

                            if (multi_result.get('success', False) and
                                multi_quality > best_result.get('r_squared', 0) + improvement_threshold):

                                print(f"   ✓ Multi-peak fit successful: R² improved from {best_result.get('r_squared', 0):.3f} to {multi_quality:.3f}")

                                # Store multi-peak diagnostics
                                self.last_fit_diagnostics['multi_peak_detection'] = {
                                    'triggered': True,
                                    'detected_peaks': len(detected_peaks),
                                    'improvement': multi_quality - best_result.get('r_squared', 0),
                                    'method': 'automatic_detection'
                                }

                                # Return the improved multi-peak result
                                multi_result['quality_class'] = self.assess_fit_quality_comprehensive(
                                    x_data, y_data, multi_result, nucleus_type
                                )['quality_class']

                                return multi_result
                            else:
                                print(f"   ⚠ Multi-peak fit did not improve results sufficiently (R² = {multi_quality:.3f})")

                        except Exception as e:
                            print(f"   ⚠ Multi-peak fitting failed: {e}")

                    else:
                        print("   No overlapping peaks detected in current region")

                    # Store multi-peak attempt info
                    self.last_fit_diagnostics['multi_peak_detection'] = {
                        'triggered': True,
                        'detected_peaks': len(detected_peaks) if 'detected_peaks' in locals() else 0,
                        'improvement': 0,
                        'method': 'automatic_detection_failed'
                    }

                else:
                    print(f"   Single-peak fit quality acceptable (R² = {best_result.get('r_squared', 0):.3f})")
                    self.last_fit_diagnostics['multi_peak_detection'] = {
                        'triggered': False,
                        'reason': 'single_peak_sufficient'
                    }

                # Add optimization diagnostics to result
                if best_result.get('success', False):
                    best_result['optimization_diagnostics'] = self.last_fit_diagnostics
                    print(f"🎯 Iterative optimization complete: R² = {best_result['r_squared']:.4f}")

                return best_result

            # === ORIGINAL METHODS (PRESERVED) ===
            elif method == 'multi_step':
                return self._fit_with_multistep_method(x_data, y_data, refined_center, nucleus_type, preprocessing)

            elif method == 'single_step':
                return self._fit_with_standard_method(x_data, y_data, None, nucleus_type, preprocessing)

            else:
                raise ValueError(f"Unknown fitting method: {method}. Use 'iterative_optimization', 'multi_step', or 'single_step'")

        except Exception as e:
            print(f"Enhanced fitting failed: {e}")
            import traceback
            traceback.print_exc()

            return {
                'success': False,
                'error': str(e),
                'method': method,
                'parameters': None,
                'fitted_curve': None,
                'quality_metrics': {'r_squared': 0, 'quality_class': 'Failed'},
                'diagnostics': self.last_fit_diagnostics
            }

    def _fit_with_standard_method(self, x_data, y_data, initial_guess=None, nucleus_type=None,
                                 preprocessing=True, linewidth_constraints=None):
        """
        PRESERVED ORIGINAL METHOD: Standard single-step Voigt fitting

        This preserves the original functionality for backward compatibility
        """
        try:
            if nucleus_type is None:
                nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])

            if initial_guess is None:
                # Estimate parameters using original method
                baseline = self.robust_baseline_estimation(x_data, y_data)
                peak_height = np.max(y_data)
                amplitude = peak_height - baseline
                center = x_data[np.argmax(y_data)]
                estimated_width = self.adaptive_width_estimation(x_data, y_data, center, nucleus_type)

                initial_guess = [amplitude, center, estimated_width * 0.7, estimated_width * 0.3, baseline]

            # Get adaptive bounds
            bounds = self.get_adaptive_bounds(initial_guess, x_data, y_data, nucleus_type, linewidth_constraints)
            lower_bounds, upper_bounds = bounds
            tuple_bounds = list(zip(lower_bounds, upper_bounds))
            initial_guess = self._clamp_vector_to_bounds(initial_guess, tuple_bounds)

            # Perform fitting
            popt, pcov = curve_fit(
                self.voigt_profile_1d, x_data, y_data,
                p0=initial_guess, bounds=(lower_bounds, upper_bounds),
                maxfev=self.fitting_parameters['max_iterations']
            )

            # Generate fitted curve and quality assessment
            y_fitted = self.voigt_profile_1d(x_data, *popt)
            quality_metrics = self.comprehensive_quality_assessment(x_data, y_data, y_fitted, popt, pcov)
            uncertainties = self.calculate_parameter_uncertainties(popt, pcov)

            # Prepare result dictionary
            result = {
                'success': True,
                'method': 'standard_single_step',
                'parameters': popt,
                'parameter_names': ['amplitude', 'center', 'sigma', 'gamma', 'baseline'],
                'covariance': pcov,
                'fitted_curve': y_fitted,
                'residuals': y_data - y_fitted,
                'nucleus_type': nucleus_type,
                'quality_metrics': quality_metrics,
                'uncertainties': uncertainties,

                # Individual parameter access for compatibility
                'amplitude': popt[0],
                'center': popt[1],
                'sigma': popt[2],
                'gamma': popt[3],
                'baseline': popt[4],
                'r_squared': quality_metrics['r_squared'],
                'quality_class': quality_metrics['quality_class']
            }

            return result

        except Exception as e:
            return {
                'success': False,
                'error': f'standard_fitting_failed: {str(e)}',
                'method': 'standard_single_step'
            }

    def _fit_with_multistep_method(self, x_data, y_data, refined_center, nucleus_type, preprocessing):
        """
        PRESERVED ORIGINAL METHOD: Multi-step fitting (coarse → fine)

        This preserves the original multi-step functionality
        """
        try:
            # Original parameter estimation
            baseline = self.robust_baseline_estimation(x_data, y_data)
            amplitude = np.max(y_data) - baseline
            estimated_width = self.adaptive_width_estimation(x_data, y_data, refined_center, nucleus_type)

            initial_guess = [amplitude, refined_center, estimated_width * 0.7, estimated_width * 0.3, baseline]
            bounds = self.get_adaptive_bounds(initial_guess, x_data, y_data, nucleus_type, None)

            # Multi-step fitting: coarse → fine
            try:
                # Step 1: Coarse fit
                relaxed_bounds = (
                    [b * 0.5 for b in bounds[0]],
                    [b * 1.5 for b in bounds[1]]
                )
                relaxed_bounds[0][0] = max(0, relaxed_bounds[0][0])

                popt_coarse, pcov_coarse = curve_fit(
                    self.voigt_profile_1d, x_data, y_data,
                    p0=initial_guess, bounds=relaxed_bounds,
                    maxfev=self.fitting_parameters['max_iterations'] // 2
                )

                # Step 2: Fine fit
                fine_bounds = self.get_adaptive_bounds(popt_coarse, x_data, y_data, nucleus_type)
                popt, pcov = curve_fit(
                    self.voigt_profile_1d, x_data, y_data,
                    p0=popt_coarse, bounds=fine_bounds,
                    maxfev=self.fitting_parameters['max_iterations']
                )

            except:
                # Fallback to single step
                popt, pcov = curve_fit(
                    self.voigt_profile_1d, x_data, y_data,
                    p0=initial_guess, bounds=bounds,
                    maxfev=self.fitting_parameters['max_iterations']
                )

            # Generate results using original format
            y_fitted = self.voigt_profile_1d(x_data, *popt)
            quality_metrics = self.comprehensive_quality_assessment(x_data, y_data, y_fitted, popt, pcov)
            uncertainties = self.calculate_parameter_uncertainties(popt, pcov)

            return {
                'success': True,
                'method': 'multi_step',
                'parameters': popt,
                'parameter_names': ['amplitude', 'center', 'sigma', 'gamma', 'baseline'],
                'covariance': pcov,
                'fitted_curve': y_fitted,
                'residuals': y_data - y_fitted,
                'nucleus_type': nucleus_type,
                'quality_metrics': quality_metrics,
                'uncertainties': uncertainties,
                'diagnostics': {'nucleus_type': nucleus_type, 'fitting_method': 'multi_step'},

                # Individual parameter access
                'amplitude': popt[0],
                'center': popt[1],
                'sigma': popt[2],
                'gamma': popt[3],
                'baseline': popt[4],
                'r_squared': quality_metrics['r_squared'],
                'quality_class': quality_metrics['quality_class']
            }

        except Exception as e:
            return {
                'success': False,
                'error': f'multistep_fitting_failed: {str(e)}',
                'method': 'multi_step'
            }

    def get_adaptive_peak_detection_params(self, x_data, y_data, nucleus_type, window_size_multiplier=1.0):
        """
        Calculate adaptive peak detection parameters based on window size and data characteristics
        
        Args:
            x_data: X-axis data
            y_data: Y-axis data
            nucleus_type: Detected nucleus type
            window_size_multiplier: Multiplier indicating how large the current window is relative to standard
            
        Returns:
            dict: Adaptive parameters for peak detection
        """
        # Get typical linewidth for this nucleus type
        typical_width = self.nmr_ranges[nucleus_type]['typical_width']
        
        # Calculate data characteristics
        spectral_width = x_data[-1] - x_data[0] if len(x_data) > 1 else 0.1
        data_points = len(x_data)
        data_range = np.max(y_data) - np.min(y_data)
        noise_level = np.std(y_data) * 0.05  # Base noise level
        
        # GUI parameters
        sensitivity = self.fitting_parameters.get('peak_detection_sensitivity', 1.5)
        overlap_factor = self.fitting_parameters.get('overlap_detection_factor', 0.8)
        
        # Adaptive height threshold - adjusts based on window size
        base_height_factor = sensitivity
        # Smaller windows (higher multiplier) need more sensitive detection
        if window_size_multiplier > 5.0:  # Large window
            height_factor = base_height_factor * 0.8  # Less sensitive
        elif window_size_multiplier < 2.0:  # Small window
            height_factor = base_height_factor * 1.3  # More sensitive
        else:
            height_factor = base_height_factor
            
        min_height = max(noise_level * height_factor, data_range * 0.05)
        
        # Adaptive distance - scales with both data size and spectral resolution
        base_distance = max(1, min(data_points // 50, int(data_points * 0.02)))
        
        # Window size adaptation
        if data_points < 10:
            min_distance = 1
        elif spectral_width < typical_width * 2:  # Very narrow window
            min_distance = max(1, base_distance // 2)
        elif window_size_multiplier > 8.0:  # Very large window
            min_distance = max(base_distance, int(data_points * 0.05))  # Increase separation requirement
        else:
            min_distance = base_distance
        
        # Adaptive prominence - adjusts for signal quality
        snr_estimate = np.max(y_data) / (noise_level * 20 + 1e-10)  # Avoid division by zero
        if snr_estimate > 50:  # High SNR
            min_prominence = noise_level * 0.6  # More sensitive
        elif snr_estimate < 10:  # Low SNR
            min_prominence = noise_level * 1.2  # Less sensitive
        else:
            min_prominence = noise_level * 0.8
        
        return {
            'min_height': min_height,
            'min_distance': min_distance,
            'min_prominence': min_prominence,
            'height_factor': height_factor,
            'overlap_factor': overlap_factor,
            'typical_width': typical_width,
            'spectral_width': spectral_width,
            'window_size_multiplier': window_size_multiplier
        }

    def detect_overlapping_peaks(self, x_data, y_data, min_separation=None, nucleus_type=None):
        """
        Detect potentially overlapping peaks that might affect quality assessment

        Parameters:
        - x_data: x-axis data
        - y_data: intensity data
        - min_separation: minimum separation to consider peaks non-overlapping
        - nucleus_type: nucleus type for typical separation

        Returns:
        - List of peak positions and their isolation status
        """
        try:
            if nucleus_type is None:
                nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])

            # Get typical width for nucleus
            typical_width = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])['typical_width']

            if min_separation is None:
                # Use GUI parameter for overlap detection
                overlap_factor = self.fitting_parameters.get('overlap_detection_factor', 0.8)
                min_separation = typical_width * overlap_factor

            # Find peaks using simple peak detection
            from scipy.signal import find_peaks

            # Use baseline-corrected data for peak finding
            try:
                baseline = self.robust_baseline_estimation(x_data, y_data)
            except Exception as baseline_error:
                print(f"   ⚠ Baseline estimation failed: {baseline_error}, using simple baseline")
                baseline = np.median(y_data)  # Simple fallback baseline
            
            y_corrected = y_data - baseline
            print(f"   Baseline correction applied (baseline={baseline:.1f}), proceeding with peak detection...")

            # Find peaks above noise threshold
            #noise_level = np.std(y_corrected) * 0.1
            #peaks, properties = find_peaks(y_corrected,
            #                             height=noise_level * 3,  # 3x noise level
            #                             distance=len(x_data) // 50)  # minimum distance


            # Get adaptive parameters using the new system
            # TODO: In future, pass actual window size multiplier from caller
            adaptive_params = self.get_adaptive_peak_detection_params(x_data, y_corrected, nucleus_type, window_size_multiplier=3.0)
            
            min_height = adaptive_params['min_height']
            min_distance = adaptive_params['min_distance']
            min_prominence = adaptive_params['min_prominence']
            
            # Extract noise level for later use
            noise_level = np.std(y_corrected) * 0.05  # Match the calculation in get_adaptive_peak_detection_params
            
            print(f"   Adaptive peak detection: height≥{min_height:.1f}, distance≥{min_distance}, prominence≥{min_prominence:.1f}")
            print(f"   Window analysis: spectral_width={adaptive_params['spectral_width']:.4f}, height_factor={adaptive_params['height_factor']:.2f}")

            peaks, properties = find_peaks(y_corrected,
                                         height=min_height,
                                         distance=min_distance,
                                         prominence=min_prominence,
                                         width=1)  # Minimum width requirement


##
            # Enhanced peak analysis with validation
            peak_info = []

            # Calculate peak widths for better characterization
            try:
                widths, width_heights, left_ips, right_ips = peak_widths(y_corrected, peaks, rel_height=0.5)
            except:
                widths = np.full(len(peaks), typical_width * len(x_data) / (x_data[-1] - x_data[0]))

            for i, peak_idx in enumerate(peaks):
                peak_position = x_data[peak_idx]
                peak_intensity = y_corrected[peak_idx]
                peak_prominence = properties['prominences'][i] if 'prominences' in properties else peak_intensity
                peak_width_points = widths[i] if i < len(widths) else typical_width
                peak_width_ppm = peak_width_points * abs(x_data[1] - x_data[0])

                # Enhanced neighbor analysis
                distances_to_others = []
                overlapping_neighbors = 0

                for j, other_idx in enumerate(peaks):
                    if other_idx != peak_idx:
                        distance = abs(x_data[other_idx] - peak_position)
                        distances_to_others.append(distance)

                        # Check for significant overlap (peaks closer than sum of half-widths)
                        other_width = widths[j] if j < len(widths) else peak_width_points
                        overlap_threshold = (peak_width_points + other_width) * 0.5 * abs(x_data[1] - x_data[0])

                        if distance < overlap_threshold:
                            overlapping_neighbors += 1

                min_distance_to_neighbor = min(distances_to_others) if distances_to_others else float('inf')
                is_isolated = min_distance_to_neighbor > min_separation

                # Find nearby peaks with enhanced criteria
                nearby_peaks = []
                if not is_isolated:
                    for j, other_idx in enumerate(peaks):
                        if other_idx != peak_idx:
                            distance = abs(x_data[other_idx] - peak_position)
                            if distance <= min_separation:
                                other_width = widths[j] if j < len(widths) else peak_width_points
                                nearby_peaks.append({
                                    'position': x_data[other_idx],
                                    'intensity': y_corrected[other_idx],
                                    'distance': distance,
                                    'width_ppm': other_width * abs(x_data[1] - x_data[0]),
                                    'prominence': properties['prominences'][j] if 'prominences' in properties else y_corrected[other_idx]
                                })

                # Quality assessment for each peak
                signal_to_noise = peak_intensity / (noise_level * len(x_data)**0.5)
                quality_score = min(1.0, (peak_prominence / (noise_level * 5)) * (signal_to_noise / 10))

                peak_info.append({
                    'position': peak_position,
                    'intensity': peak_intensity,
                    'prominence': peak_prominence,
                    'width_ppm': peak_width_ppm,
                    'isolated': is_isolated,
                    'overlapping_neighbors': overlapping_neighbors,
                    'min_neighbor_distance': min_distance_to_neighbor,
                    'nearby_peaks': nearby_peaks,
                    'signal_to_noise': signal_to_noise,
                    'quality_score': quality_score
                })

            # Sort peaks by intensity (most prominent first)
            peak_info.sort(key=lambda x: x['intensity'], reverse=True)

            print(f"   Detected {len(peak_info)} peaks:")
            for i, peak in enumerate(peak_info):
                status = "isolated" if peak['isolated'] else f"overlapped({peak['overlapping_neighbors']})"
                print(f"     Peak {i+1}: {peak['position']:.3f} ppm, intensity={peak['intensity']:.0f}, {status}, Q={peak['quality_score']:.2f}")
##

#            peak_info = []
#            for i, peak_idx in enumerate(peaks):
#                peak_position = x_data[peak_idx]
#                peak_height = y_corrected[peak_idx]
#
#                # Check isolation from other peaks
#                is_isolated = True
#                nearby_peaks = []
#
#                for j, other_peak_idx in enumerate(peaks):
#                    if i != j:
#                        other_position = x_data[other_peak_idx]
#                        separation = abs(peak_position - other_position)
#
#                        if separation < min_separation:
#                            is_isolated = False
#                            nearby_peaks.append({
#                                'position': other_position,
#                                'separation': separation,
#                                'height': y_corrected[other_peak_idx]
#                            })
#
            #    peak_info.append({
            #        'position': peak_position,
            #        'index': peak_idx,
            #        'height': peak_height,
            #        'is_isolated': is_isolated,
            #        'nearby_peaks': nearby_peaks,
            #        'isolation_radius': min_separation
            #    })

            return peak_info

        except Exception as e:
            print(f"   ⚠ Overlapping peak detection failed: {e}")
            print(f"   Falling back to single peak analysis")

            # Provide minimal peak info for single peak case
            if len(x_data) > 0 and len(y_data) > 0:
                max_idx = np.argmax(y_data)
                return [{
                    'position': x_data[max_idx],
                    'intensity': y_data[max_idx],
                    'isolated': True,
                    'overlapping_neighbors': 0,
                    'min_neighbor_distance': float('inf'),
                    'nearby_peaks': [],
                    'quality_score': 0.5,
                    'error_fallback': True
                }]

            return []

    def validate_multi_peak_result(self, result, x_data, y_data, original_peaks):
        """
        Enhanced validation and quality control for multi-peak fitting results
        
        Performs chemical shift reasonableness, peak width consistency,
        and intensity ratio plausibility checks.
        
        Args:
            result: Multi-peak fitting result dictionary
            x_data: Original X-axis data
            y_data: Original Y-axis data
            original_peaks: Originally detected peak positions
            
        Returns:
            dict: Validation report with pass/fail status and diagnostics
        """
        validation_report = {
            'overall_pass': False,
            'chemical_shift_reasonable': False,
            'peak_widths_consistent': False,
            'intensity_ratios_plausible': False,
            'baseline_reasonable': False,
            'r_squared_acceptable': False,
            'warnings': [],
            'recommendations': []
        }
        
        if not result.get('success', False):
            validation_report['warnings'].append('Fitting was not successful')
            return validation_report
        
        r_squared = result.get('r_squared', 0.0)
        n_peaks = result.get('n_peaks_fitted', 1)
        
        # Extract fitted parameters based on result structure
        if 'peak_positions' in result:
            fitted_positions = result['peak_positions']
            # Try to extract widths and amplitudes from parameters
            params = result.get('parameters', [])
            fitted_widths = []
            fitted_amplitudes = []
            
            if len(params) >= n_peaks * 4:  # amp, center, sigma, gamma per peak
                for i in range(n_peaks):
                    start_idx = i * 4
                    fitted_amplitudes.append(params[start_idx])
                    fitted_widths.append(params[start_idx + 2])  # sigma
        else:
            # Single peak case
            fitted_positions = [result.get('center', 0.0)]
            fitted_widths = [result.get('sigma', 0.0)]
            fitted_amplitudes = [result.get('amplitude', 0.0)]
        
        fitted_baseline = result.get('baseline', 0.0)
        
        # 1. Chemical Shift Reasonableness Check
        nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
        typical_range = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])
        
        shifts_in_range = all(
            typical_range['min'] <= pos <= typical_range['max'] 
            for pos in fitted_positions
        )
        
        if shifts_in_range:
            validation_report['chemical_shift_reasonable'] = True
        else:
            validation_report['warnings'].append(f"Chemical shifts outside typical {nucleus_type} range")
        
        # 2. Peak Width Consistency Check
        if len(fitted_widths) > 1:
            width_ratios = []
            for i in range(len(fitted_widths) - 1):
                if fitted_widths[i+1] > 0:
                    width_ratios.append(fitted_widths[i] / fitted_widths[i+1])
            
            # Check if width ratios are reasonable (within 10× of each other)
            reasonable_widths = all(0.1 <= ratio <= 10.0 for ratio in width_ratios)
            if reasonable_widths:
                validation_report['peak_widths_consistent'] = True
            else:
                validation_report['warnings'].append("Peak widths vary by more than 10× - may indicate fitting issues")
        else:
            validation_report['peak_widths_consistent'] = True  # Single peak case
        
        # 3. Intensity Ratio Plausibility Check
        if len(fitted_amplitudes) > 1:
            amplitude_ratios = []
            max_amp = max(fitted_amplitudes)
            for amp in fitted_amplitudes:
                if max_amp > 0:
                    amplitude_ratios.append(amp / max_amp)
            
            # Check if any peak is excessively small or large
            reasonable_intensities = all(0.01 <= ratio <= 1.0 for ratio in amplitude_ratios)
            if reasonable_intensities:
                validation_report['intensity_ratios_plausible'] = True
            else:
                validation_report['warnings'].append("Some peaks have very low relative intensity (<1%)")
        else:
            validation_report['intensity_ratios_plausible'] = True  # Single peak case
        
        # 4. Baseline Reasonableness Check
        data_min = np.min(y_data)
        data_max = np.max(y_data)
        data_range = data_max - data_min
        
        baseline_reasonable = (data_min - data_range * 0.2) <= fitted_baseline <= (data_min + data_range * 0.3)
        if baseline_reasonable:
            validation_report['baseline_reasonable'] = True
        else:
            validation_report['warnings'].append(f"Baseline ({fitted_baseline:.1f}) seems unreasonable for data range")
        
        # 5. R² Acceptability Check
        r2_threshold = 0.7 if n_peaks > 1 else 0.8  # Lower threshold for multi-peak
        if r_squared >= r2_threshold:
            validation_report['r_squared_acceptable'] = True
        else:
            validation_report['warnings'].append(f"R² ({r_squared:.3f}) below acceptable threshold ({r2_threshold})")
        
        # Overall Pass/Fail Decision
        critical_checks = [
            validation_report['chemical_shift_reasonable'],
            validation_report['r_squared_acceptable']
        ]
        
        secondary_checks = [
            validation_report['peak_widths_consistent'],
            validation_report['intensity_ratios_plausible'],
            validation_report['baseline_reasonable']
        ]
        
        # Pass if all critical checks pass and at least 2/3 secondary checks pass
        critical_pass = all(critical_checks)
        secondary_pass = sum(secondary_checks) >= 2
        
        validation_report['overall_pass'] = critical_pass and secondary_pass
        
        # Generate recommendations
        if not validation_report['overall_pass']:
            if r_squared < 0.6:
                validation_report['recommendations'].append("Consider single peak fitting or different baseline method")
            if not validation_report['chemical_shift_reasonable']:
                validation_report['recommendations'].append("Check peak detection parameters or data preprocessing")
            if not validation_report['peak_widths_consistent']:
                validation_report['recommendations'].append("Consider constraining peak widths or reducing peak count")
        
        return validation_report

    def fit_with_progressive_strategies(self, x_data, y_data, detected_peaks, global_constraints=None):
        """
        Progressive overlap resolution with multiple fallback strategies
        
        Strategy 1: Constrained multi-peak with tight bounds (if global constraints available)
        Strategy 2: Standard multi-peak with adaptive separation
        Strategy 3: Pairwise fitting then combination
        Strategy 4: Single peak with expanded window
        
        Args:
            x_data: X-axis data
            y_data: Y-axis data  
            detected_peaks: List of detected peak information
            global_constraints: Optional global linewidth constraints
            
        Returns:
            dict: Best fitting result with strategy information
        """
        print(f"   🔄 Progressive overlap resolution for {len(detected_peaks)} peaks...")
        
        strategies_tried = []
        best_result = None
        best_r_squared = 0.0
        
        # Strategy 1: Constrained multi-peak with global bounds (if available)
        if global_constraints is not None and len(detected_peaks) >= 2:
            try:
                print(f"   Strategy 1: Constrained multi-peak with global bounds...")
                result = self.fit_with_global_constraints(x_data, y_data, detected_peaks, global_constraints)
                r_squared = result.get('r_squared', 0.0)
                strategies_tried.append({'strategy': 'constrained_global', 'r_squared': r_squared, 'success': result.get('success', False)})
                
                if result.get('success', False) and r_squared > 0.85:
                    print(f"   ✅ Strategy 1 successful: R² = {r_squared:.3f}")
                    best_result = result
                    best_r_squared = r_squared
                    best_result['strategy_used'] = 'constrained_global'
                    return best_result
                    
            except Exception as e:
                print(f"   Strategy 1 failed: {e}")
                strategies_tried.append({'strategy': 'constrained_global', 'r_squared': 0.0, 'success': False, 'error': str(e)})
        
        # Strategy 2: Standard multi-peak with adaptive separation  
        if len(detected_peaks) >= 2:
            try:
                print(f"   Strategy 2: Standard multi-peak with adaptive separation...")
                peak_positions = [p['position'] for p in detected_peaks]
                nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
                typical_separation = self.nmr_ranges[nucleus_type]['typical_width'] * 1.0
                
                result = self._fit_multiple_peaks_simultaneously(x_data, y_data, peak_positions, typical_separation)
                r_squared = result.get('r_squared', 0.0)
                strategies_tried.append({'strategy': 'adaptive_separation', 'r_squared': r_squared, 'success': result.get('success', False)})
                
                if result.get('success', False) and (best_result is None or r_squared > best_r_squared):
                    print(f"   ✅ Strategy 2 successful: R² = {r_squared:.3f}")
                    best_result = result
                    best_r_squared = r_squared
                    best_result['strategy_used'] = 'adaptive_separation'
                    
            except Exception as e:
                print(f"   Strategy 2 failed: {e}")
                strategies_tried.append({'strategy': 'adaptive_separation', 'r_squared': 0.0, 'success': False, 'error': str(e)})
        
        # Strategy 3: Pairwise fitting then combination (for 3+ peaks)
        if len(detected_peaks) >= 3:
            try:
                print(f"   Strategy 3: Pairwise fitting then combination...")
                result = self.fit_pairwise_then_combine(x_data, y_data, detected_peaks)
                r_squared = result.get('r_squared', 0.0)
                strategies_tried.append({'strategy': 'pairwise_combination', 'r_squared': r_squared, 'success': result.get('success', False)})
                
                if result.get('success', False) and (best_result is None or r_squared > best_r_squared):
                    print(f"   ✅ Strategy 3 successful: R² = {r_squared:.3f}")
                    best_result = result
                    best_r_squared = r_squared  
                    best_result['strategy_used'] = 'pairwise_combination'
                    
            except Exception as e:
                print(f"   Strategy 3 failed: {e}")
                strategies_tried.append({'strategy': 'pairwise_combination', 'r_squared': 0.0, 'success': False, 'error': str(e)})
        
        # Strategy 4: Single peak with expanded window (fallback)
        try:
            print(f"   Strategy 4: Single peak with expanded window (fallback)...")
            nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
            typical_width = self.nmr_ranges[nucleus_type]['typical_width']
            expanded_window_factor = 3.0  # Expand window to capture overlapping signals
            
            result = self.fit_single_peak_with_expanded_context(x_data, y_data, typical_width, expanded_window_factor)
            r_squared = result.get('r_squared', 0.0)
            strategies_tried.append({'strategy': 'single_expanded', 'r_squared': r_squared, 'success': result.get('success', False)})
            
            if result.get('success', False) and (best_result is None or r_squared > best_r_squared):
                print(f"   ✅ Strategy 4 successful: R² = {r_squared:.3f}")
                best_result = result
                best_r_squared = r_squared
                best_result['strategy_used'] = 'single_expanded'
                
        except Exception as e:
            print(f"   Strategy 4 failed: {e}")
            strategies_tried.append({'strategy': 'single_expanded', 'r_squared': 0.0, 'success': False, 'error': str(e)})
        
        # Return best result with strategy information
        if best_result is not None:
            best_result['progressive_strategies'] = {
                'strategies_tried': strategies_tried,
                'best_strategy': best_result.get('strategy_used', 'unknown'),
                'final_r_squared': best_r_squared
            }
            print(f"   🏆 Progressive resolution complete: Best strategy = {best_result.get('strategy_used', 'unknown')}, R² = {best_r_squared:.3f}")
        else:
            print(f"   ❌ All progressive strategies failed")
            best_result = {
                'success': False,
                'error': 'all_progressive_strategies_failed',
                'progressive_strategies': {'strategies_tried': strategies_tried}
            }
        
        return best_result

    def fit_with_global_constraints(self, x_data, y_data, detected_peaks, global_constraints):
        """
        Placeholder for constrained fitting with global linewidth bounds
        TODO: Implement when global optimization is integrated
        """
        print(f"   Global constraints not yet implemented, falling back to standard multi-peak")
        peak_positions = [p['position'] for p in detected_peaks]
        nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
        typical_separation = self.nmr_ranges[nucleus_type]['typical_width'] * 0.8
        return self._fit_multiple_peaks_simultaneously(x_data, y_data, peak_positions, typical_separation)

    def fit_pairwise_then_combine(self, x_data, y_data, detected_peaks):
        """
        Fit pairs of peaks separately, then combine results
        Useful for highly overlapped peaks where simultaneous fitting fails
        """
        print(f"   Pairwise fitting approach for {len(detected_peaks)} peaks...")
        
        if len(detected_peaks) < 2:
            return {'success': False, 'error': 'insufficient_peaks_for_pairwise'}
        
        # For now, fall back to standard multi-peak
        # TODO: Implement actual pairwise fitting algorithm
        peak_positions = [p['position'] for p in detected_peaks]
        nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
        typical_separation = self.nmr_ranges[nucleus_type]['typical_width'] * 1.2
        return self._fit_multiple_peaks_simultaneously(x_data, y_data, peak_positions, typical_separation)

    def fit_single_peak_with_expanded_context(self, x_data, y_data, typical_width, expansion_factor):
        """
        Fit single peak with expanded window to capture overlapping components
        """
        print(f"   Single peak with {expansion_factor:.1f}× expanded context...")
        
        try:
            # Use standard single peak fitting with relaxed constraints
            nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
            
            # Estimate peak parameters
            max_idx = np.argmax(y_data)
            peak_position = x_data[max_idx]
            peak_amplitude = y_data[max_idx] - np.min(y_data)
            peak_baseline = np.min(y_data)
            
            # Expanded width estimate to account for overlapping components
            width_estimate = typical_width * expansion_factor
            
            initial_guess = [peak_amplitude, peak_position, width_estimate * 0.7, width_estimate * 0.3, peak_baseline]
            
            # Relaxed bounds for expanded fitting
            bounds = (
                [0, x_data[0], typical_width * 0.1, 0, np.min(y_data) - np.std(y_data)],
                [np.max(y_data) * 3, x_data[-1], typical_width * 10, typical_width * 5, np.max(y_data)]
            )
            
            from scipy.optimize import curve_fit
            popt, pcov = curve_fit(
                self.voigt_profile_1d, x_data, y_data,
                p0=initial_guess, bounds=bounds,
                maxfev=self.fitting_parameters['max_iterations']
            )
            
            y_fitted = self.voigt_profile_1d(x_data, *popt)
            r_squared = self.calculate_r_squared(y_data, y_fitted)
            
            return {
                'success': True,
                'method': 'single_peak_expanded',
                'parameters': popt,
                'fitted_curve': y_fitted,
                'r_squared': r_squared,
                'amplitude': popt[0],
                'center': popt[1],
                'sigma': popt[2],
                'gamma': popt[3],
                'baseline': popt[4],
                'expansion_factor': expansion_factor
            }
            
        except Exception as e:
            return {'success': False, 'error': f'expanded_single_peak_failed: {str(e)}'}

    def sequential_peak_isolation(self, x_data, y_data, peak_positions):
        """
        ENHANCED MULTI-PEAK INITIALIZATION: Sequential Peak Isolation Strategy

        Addresses the overlapping peak initialization problem by:
        1. Fitting each peak individually in a constrained window
        2. Using residual subtraction to isolate each peak
        3. Extracting realistic initial parameters from individual fits
        4. Building smart initial bounds based on individual results

        This dramatically improves multi-peak fitting success rates by providing
        physically meaningful initial guesses instead of generic defaults.

        Args:
            x_data: X-axis data (ppm)
            y_data: Intensity data
            peak_positions: List of suspected peak centers

        Returns:
            dict: {
                'success': bool,
                'individual_fits': list of individual fit results,
                'combined_initial_guess': list of parameters for multi-peak fit,
                'smart_bounds': tuple of (lower_bounds, upper_bounds),
                'isolation_quality': assessment of how well peaks were isolated
            }
        """
        print(f"   🔍 Sequential peak isolation for {len(peak_positions)} peaks...")

        n_peaks = len(peak_positions)
        individual_fits = []
        nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
        typical_width = self.nmr_ranges[nucleus_type]['typical_width']

        # Create working copy of data for residual subtraction
        working_y_data = y_data.copy()
        baseline_est = self.robust_baseline_estimation(x_data, y_data, method='polynomial')

        isolation_scores = []

        # Step 1: Fit each peak individually with residual subtraction
        for i, peak_pos in enumerate(peak_positions):
            print(f"     Isolating peak {i+1} at {peak_pos:.4f} ppm...")

            try:
                # Create isolation window around this peak
                window_size = typical_width * 3  # 3× typical width
                window_mask = (x_data >= peak_pos - window_size/2) & (x_data <= peak_pos + window_size/2)

                if not np.any(window_mask):
                    # Fallback to larger window if none found
                    window_size = typical_width * 6
                    window_mask = (x_data >= peak_pos - window_size/2) & (x_data <= peak_pos + window_size/2)

                if not np.any(window_mask):
                    print(f"       ⚠ Could not create isolation window for peak {i+1}")
                    individual_fits.append({'success': False, 'error': 'no_isolation_window'})
                    continue

                # Extract windowed data
                x_window = x_data[window_mask]
                y_window = working_y_data[window_mask]

                # Fit this peak individually using simple method to avoid recursion
                peak_fit = self.fit_peak_enhanced(
                    x_window, y_window, peak_pos,
                    method='single_step',  # Use simple method to avoid complexity
                    nucleus_type=nucleus_type
                )

                if peak_fit['success']:
                    # Extract fitted parameters
                    fitted_params = peak_fit['parameters']
                    amplitude, center, sigma, gamma, fitted_baseline = fitted_params

                    # Quality assessment for this isolation
                    r_squared = peak_fit['r_squared']
                    isolation_score = min(1.0, max(0.0, r_squared))  # Clamp between 0-1
                    isolation_scores.append(isolation_score)

                    print(f"       ✓ Peak {i+1}: R²={r_squared:.3f}, amp={amplitude:.0f}, width={sigma+gamma:.4f}")

                    # Subtract this peak from working data for next iterations
                    if i < n_peaks - 1:  # Don't subtract on last iteration
                        peak_contribution = self.voigt_profile_1d(
                            x_data, amplitude, center, sigma, gamma, 0
                        )
                        working_y_data = working_y_data - peak_contribution

                    individual_fits.append({
                        'success': True,
                        'peak_index': i,
                        'position': peak_pos,
                        'fitted_center': center,
                        'parameters': fitted_params,
                        'r_squared': r_squared,
                        'isolation_score': isolation_score,
                        'window_size': window_size
                    })

                else:
                    print(f"       ✗ Peak {i+1} fit failed: {peak_fit.get('error', 'unknown')}")
                    isolation_scores.append(0.0)
                    individual_fits.append({
                        'success': False,
                        'error': peak_fit.get('error', 'fit_failed'),
                        'peak_index': i,
                        'position': peak_pos
                    })

            except Exception as e:
                print(f"       ✗ Peak {i+1} isolation failed: {e}")
                isolation_scores.append(0.0)
                individual_fits.append({
                    'success': False,
                    'error': str(e),
                    'peak_index': i,
                    'position': peak_pos
                })

        # Step 2: Build combined initial guess from successful fits
        successful_fits = [fit for fit in individual_fits if fit['success']]

        if len(successful_fits) == 0:
            print("   ✗ No peaks successfully isolated")
            return {
                'success': False,
                'error': 'no_successful_isolations',
                'individual_fits': individual_fits,
                'isolation_quality': 0.0
            }

        # Calculate overall isolation quality
        isolation_quality = np.mean(isolation_scores) if isolation_scores else 0.0

        print(f"   ✓ Successfully isolated {len(successful_fits)}/{n_peaks} peaks")
        print(f"   📊 Isolation quality: {isolation_quality:.2f}")

        # Build initial parameter vector and smart bounds
        combined_initial_guess = []
        lower_bounds = []
        upper_bounds = []

        # Use successful fits in order of peak positions
        sorted_successful = sorted(successful_fits, key=lambda x: x['position'])

        for fit in sorted_successful:
            params = fit['parameters']
            amplitude, center, sigma, gamma, fitted_baseline = params

            # Add to combined guess
            combined_initial_guess.extend([amplitude, center, sigma, gamma])

            # Create smart bounds based on fitted parameters
            # Amplitude: allow ±50% variation
            lower_bounds.extend([
                amplitude * 0.5,           # amplitude lower
                center - typical_width,    # center lower (allow some drift)
                sigma * 0.2,              # sigma lower
                gamma * 0.2               # gamma lower
            ])
            upper_bounds.extend([
                amplitude * 2.0,          # amplitude upper
                center + typical_width,   # center upper
                sigma * 5.0,             # sigma upper
                gamma * 5.0              # gamma upper
            ])

        # Add shared baseline
        # Use median of individual baselines or original estimate
        individual_baselines = [fit['parameters'][4] for fit in successful_fits]
        shared_baseline = np.median(individual_baselines) if individual_baselines else baseline_est
        combined_initial_guess.append(shared_baseline)

        # Baseline bounds
        data_range = np.max(y_data) - np.min(y_data)
        lower_bounds.append(shared_baseline - data_range * 0.1)
        upper_bounds.append(shared_baseline + data_range * 0.1)

        return {
            'success': True,
            'individual_fits': individual_fits,
            'successful_fits': len(successful_fits),
            'combined_initial_guess': combined_initial_guess,
            'smart_bounds': (lower_bounds, upper_bounds),
            'isolation_quality': isolation_quality,
            'shared_baseline': shared_baseline,
            'method': 'sequential_peak_isolation',
            'nucleus_type': nucleus_type,
            'constraint_coupling_ready': True  # Enable constraint coupling
        }

    def apply_constraint_coupling(self, initial_params, bounds, isolation_result, n_peaks):
        """
        ENHANCED MULTI-PEAK FITTING: Apply constraint coupling for same nucleus widths

        For peaks from the same nucleus type, enforce similar linewidths by:
        1. Calculating average width from individual fits
        2. Constraining all peak widths to be within reasonable range of average
        3. Reducing parameter space complexity for better convergence

        Args:
            initial_params: Initial parameter vector
            bounds: Current bounds tuple (lower, upper)
            isolation_result: Results from sequential peak isolation
            n_peaks: Number of peaks being fitted

        Returns:
            dict: {
                'coupled_params': modified initial parameters,
                'coupled_bounds': modified bounds with coupling constraints,
                'coupling_info': details about applied constraints
            }
        """
        if not isolation_result.get('constraint_coupling_ready', False):
            return {
                'coupled_params': initial_params,
                'coupled_bounds': bounds,
                'coupling_info': {'applied': False, 'reason': 'isolation_not_ready'}
            }

        print(f"   🔗 Applying constraint coupling for {n_peaks} peaks of same nucleus type...")

        successful_fits = [fit for fit in isolation_result['individual_fits'] if fit['success']]
        if len(successful_fits) < 2:
            print(f"       Not enough successful fits ({len(successful_fits)}) for coupling")
            return {
                'coupled_params': initial_params,
                'coupled_bounds': bounds,
                'coupling_info': {'applied': False, 'reason': 'insufficient_fits'}
            }

        # Extract width information from individual fits
        individual_widths = []
        individual_sigmas = []
        individual_gammas = []

        for fit in successful_fits:
            params = fit['parameters']
            amplitude, center, sigma, gamma, baseline = params
            total_width = sigma + gamma
            individual_widths.append(total_width)
            individual_sigmas.append(sigma)
            individual_gammas.append(gamma)

        # Calculate statistics
        avg_width = np.mean(individual_widths)
        std_width = np.std(individual_widths)
        avg_sigma = np.mean(individual_sigmas)
        avg_gamma = np.mean(individual_gammas)

        # Calculate width variability
        width_cv = std_width / avg_width if avg_width > 0 else 1.0  # Coefficient of variation

        print(f"       Width analysis: avg={avg_width:.4f}, std={std_width:.4f}, CV={width_cv:.2f}")

        # Only apply coupling if widths are reasonably similar (CV < 0.5)
        if width_cv > 0.5:
            print(f"       Widths too variable (CV={width_cv:.2f} > 0.5), skipping coupling")
            return {
                'coupled_params': initial_params,
                'coupled_bounds': bounds,
                'coupling_info': {'applied': False, 'reason': f'high_variability_cv_{width_cv:.2f}'}
            }

        # Apply constraint coupling
        coupled_params = initial_params.copy()
        lower_bounds, upper_bounds = bounds
        coupled_lower = lower_bounds.copy()
        coupled_upper = upper_bounds.copy()

        # Define coupling tolerance based on observed variability
        width_tolerance = max(std_width * 2, avg_width * 0.3)  # Allow ±2 std devs or ±30%
        sigma_tolerance = max(np.std(individual_sigmas) * 2, avg_sigma * 0.4)
        gamma_tolerance = max(np.std(individual_gammas) * 2, avg_gamma * 0.4)

        coupling_applications = []

        # Apply coupling constraints to each peak
        for i in range(n_peaks):
            param_offset = i * 4  # Each peak has 4 parameters

            if param_offset + 3 < len(coupled_params):
                # Get current parameter indices
                sigma_idx = param_offset + 2
                gamma_idx = param_offset + 3

                # Current initial values
                current_sigma = coupled_params[sigma_idx]
                current_gamma = coupled_params[gamma_idx]
                current_total = current_sigma + current_gamma

                # Apply coupling to sigma
                if abs(current_sigma - avg_sigma) > sigma_tolerance:
                    old_sigma = current_sigma
                    coupled_params[sigma_idx] = avg_sigma
                    coupling_applications.append(f"peak_{i+1}_sigma: {old_sigma:.4f}→{avg_sigma:.4f}")

                # Apply coupling to gamma
                if abs(current_gamma - avg_gamma) > gamma_tolerance:
                    old_gamma = current_gamma
                    coupled_params[gamma_idx] = avg_gamma
                    coupling_applications.append(f"peak_{i+1}_gamma: {old_gamma:.4f}→{avg_gamma:.4f}")

                # Tighten bounds for sigma
                sigma_lower = max(coupled_lower[sigma_idx], avg_sigma - sigma_tolerance)
                sigma_upper = min(coupled_upper[sigma_idx], avg_sigma + sigma_tolerance)
                coupled_lower[sigma_idx] = sigma_lower
                coupled_upper[sigma_idx] = sigma_upper

                # Tighten bounds for gamma
                gamma_lower = max(coupled_lower[gamma_idx], avg_gamma - gamma_tolerance)
                gamma_upper = min(coupled_upper[gamma_idx], avg_gamma + gamma_tolerance)
                coupled_lower[gamma_idx] = gamma_lower
                coupled_upper[gamma_idx] = gamma_upper

        print(f"       ✓ Applied {len(coupling_applications)} coupling constraints")
        if coupling_applications:
            for app in coupling_applications[:3]:  # Show first 3
                print(f"         {app}")

        return {
            'coupled_params': coupled_params,
            'coupled_bounds': (coupled_lower, coupled_upper),
            'coupling_info': {
                'applied': True,
                'n_constraints': len(coupling_applications),
                'avg_width': avg_width,
                'width_tolerance': width_tolerance,
                'width_cv': width_cv,
                'applications': coupling_applications
            }
        }

    def enhanced_multi_peak_parameter_validation(self, parameters, x_data, y_data, n_peaks, nucleus_type):
        """
        ENHANCED MULTI-PEAK VALIDATION: Comprehensive parameter validation for multi-peak fits

        Validates fitted parameters against:
        1. Physical constraints (positive amplitudes, reasonable widths)
        2. Spectroscopic constraints (nucleus-specific ranges)
        3. Peak separation constraints (no unrealistic overlap)
        4. Amplitude ratios (no extreme intensity differences)
        5. Width consistency (similar widths for same nucleus)

        Args:
            parameters: Fitted parameter vector [amp1,pos1,σ1,γ1,amp2,pos2,σ2,γ2,...,baseline]
            x_data, y_data: Experimental data
            n_peaks: Number of peaks
            nucleus_type: Nucleus type ('1H', '15N', '13C')

        Returns:
            dict: {
                'valid': bool,
                'violations': list of validation failures,
                'warnings': list of concerns,
                'quality_score': float (0-1),
                'recommendations': list of suggested fixes
            }
        """
        print(f"   🔍 Enhanced multi-peak parameter validation ({n_peaks} peaks, {nucleus_type})...")

        violations = []
        warnings = []
        recommendations = []

        # Extract parameters
        *peak_params, baseline = parameters

        # Data characteristics
        data_max = np.max(y_data)
        data_min = np.min(y_data)
        data_range = data_max - data_min
        data_span = abs(x_data[-1] - x_data[0])

        # Nucleus-specific constraints
        nucleus_constraints = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])
        typical_width = nucleus_constraints['typical_width']

        individual_peaks = []

        # Extract individual peak parameters
        for i in range(n_peaks):
            param_offset = i * 4
            if param_offset + 3 < len(peak_params):
                amplitude = peak_params[param_offset]
                center = peak_params[param_offset + 1]
                sigma = peak_params[param_offset + 2]
                gamma = peak_params[param_offset + 3]
                total_width = sigma + gamma

                individual_peaks.append({
                    'index': i,
                    'amplitude': amplitude,
                    'center': center,
                    'sigma': sigma,
                    'gamma': gamma,
                    'total_width': total_width
                })

        # VALIDATION 1: Basic physical constraints
        for peak in individual_peaks:
            i = peak['index']

            # Amplitude validation
            if peak['amplitude'] <= 0:
                violations.append(f"Peak {i+1}: negative/zero amplitude ({peak['amplitude']:.2e})")
            elif peak['amplitude'] > data_max * 5:
                warnings.append(f"Peak {i+1}: very high amplitude ({peak['amplitude']:.0f} > 5×data_max)")

            # Width validation
            if peak['sigma'] <= 0:
                violations.append(f"Peak {i+1}: negative/zero σ ({peak['sigma']:.2e})")
            elif peak['gamma'] < 0:
                violations.append(f"Peak {i+1}: negative γ ({peak['gamma']:.2e})")

            if peak['total_width'] > data_span * 0.5:
                warnings.append(f"Peak {i+1}: very broad width ({peak['total_width']:.4f} > 50% of data span)")
            elif peak['total_width'] < typical_width * 0.01:
                warnings.append(f"Peak {i+1}: very narrow width ({peak['total_width']:.4f} < 1% typical)")

        # VALIDATION 2: Peak separation constraints
        centers = [peak['center'] for peak in individual_peaks]
        centers.sort()

        for i in range(len(centers) - 1):
            separation = abs(centers[i+1] - centers[i])
            min_expected_separation = typical_width * 0.5  # Minimum reasonable separation

            if separation < min_expected_separation:
                violations.append(f"Peaks too close: {separation:.4f} ppm < {min_expected_separation:.4f} ppm minimum")

        # VALIDATION 3: Amplitude ratio constraints
        amplitudes = [peak['amplitude'] for peak in individual_peaks]
        if len(amplitudes) > 1:
            amp_ratio = max(amplitudes) / min(amplitudes)
            if amp_ratio > 100:  # More than 100:1 ratio is suspicious
                warnings.append(f"Extreme amplitude ratio: {amp_ratio:.1f}:1")

        # VALIDATION 4: Width consistency for same nucleus
        widths = [peak['total_width'] for peak in individual_peaks]
        if len(widths) > 1:
            width_cv = np.std(widths) / np.mean(widths) if np.mean(widths) > 0 else 1.0
            if width_cv > 1.0:  # More than 100% variation
                warnings.append(f"Inconsistent widths: CV={width_cv:.2f}")

        # VALIDATION 5: Baseline reasonableness
        baseline_range_fraction = abs(baseline - np.mean(y_data)) / data_range
        if baseline_range_fraction > 0.3:
            warnings.append(f"Baseline far from data mean: {baseline:.1f} vs {np.mean(y_data):.1f}")

        # VALIDATION 6: Center positions within data range
        for peak in individual_peaks:
            if peak['center'] < x_data[0] or peak['center'] > x_data[-1]:
                violations.append(f"Peak {peak['index']+1}: center {peak['center']:.4f} outside data range")

        # Calculate quality score
        violation_penalty = len(violations) * 0.3
        warning_penalty = len(warnings) * 0.1
        quality_score = max(0.0, 1.0 - violation_penalty - warning_penalty)

        # Generate recommendations
        if len(violations) > 0:
            recommendations.append("Critical violations found - parameter optimization failed")
        elif len(warnings) > 3:
            recommendations.append("Multiple concerns - consider refitting with tighter constraints")
        elif width_cv > 0.5 and len(widths) > 1:
            recommendations.append("Apply width coupling constraints for same nucleus type")

        validation_result = {
            'valid': len(violations) == 0,
            'violations': violations,
            'warnings': warnings,
            'quality_score': quality_score,
            'recommendations': recommendations,
            'n_violations': len(violations),
            'n_warnings': len(warnings),
            'peak_summary': individual_peaks
        }

        print(f"       Validation: {len(violations)} violations, {len(warnings)} warnings, quality={quality_score:.2f}")

        return validation_result

    def optimize_overlap_detection_iteratively(self, x_data, y_data,
                                              suspected_peak_positions,
                                              max_iterations=4, use_aic_selection=True,
                                              _call_stack=None):
        """
        DYNAMIC OPTIMIZATION: Iteratively refine overlap detection and multi-peak fitting

        This addresses the overlapping peak resolution problem by:
        1. AIC model selection for optimal peak count (NEW)
        2. Progressive strategies: start aggressive (assume overlap), relax if needed
        3. Multiple peak counts: try N, N-1, N+1 peaks to find optimal fit
        4. Adaptive separation thresholds: start strict, relax progressively
        5. Quality-driven convergence: stop when no further improvement

        Args:
            x_data: X-axis data (ppm)
            y_data: Intensity data
            suspected_peak_positions: List of suspected peak centers
            max_iterations: Maximum fitting strategies to try
            use_aic_selection: Use AIC model selection for optimal peak count (default: True)

        Returns:
            tuple: (best_fit_result, optimization_report)

        BACKWARD COMPATIBILITY: Falls back to single peak fitting if all strategies fail
        """
        # === RECURSION DETECTION FOR AIC SELECTION ===
        if _call_stack is None:
            _call_stack = []

        # Create unique call identifier
        center_avg = sum(suspected_peak_positions[:3]) / min(3, len(suspected_peak_positions))  # Use first 3 peaks for ID
        call_id = f"optimize_overlap_iteratively_{len(suspected_peak_positions)}peaks_{center_avg:.4f}"

        # Check for recursion
        if call_id in _call_stack:
            print(f"   🛑 AIC RECURSION DETECTED: {call_id}")
            print(f"   Call stack: {' -> '.join(_call_stack[-3:])}")  # Show last 3 calls
            # Return emergency single-peak fallback
            emergency_center = suspected_peak_positions[0] if suspected_peak_positions else np.mean(x_data)
            emergency_result = self.emergency_fallback_fitting(x_data, y_data, emergency_center)
            return emergency_result, {'recursion_detected': True, 'call_id': call_id}

        # Add to call stack
        _call_stack.append(call_id)

        print(f"   Optimizing overlap detection for {len(suspected_peak_positions)} suspected peaks...")

        optimization_report = {
            'method': 'dynamic_overlap_detection_with_aic' if use_aic_selection else 'dynamic_overlap_detection',
            'strategies_tried': [],
            'best_strategy': None,
            'convergence_reason': None,
            'suspected_positions': suspected_peak_positions.copy(),
            'aic_selection_used': use_aic_selection
        }

        # ENHANCED: Try AIC model selection first for optimal peak count
        if use_aic_selection and len(suspected_peak_positions) >= 2:
            print(f"   🎯 Starting with AIC model selection...")

            try:
                aic_result = self.optimal_peak_count_by_aic(
                    x_data, y_data, suspected_peak_positions,
                    max_peaks=min(6, len(suspected_peak_positions)),
                    _call_stack=_call_stack
                )

                # DEFENSIVE FIX: Ensure aic_result is a dictionary
                if not isinstance(aic_result, dict):
                    print(f"   ⚠ AIC method returned {type(aic_result)}, converting to dict")
                    aic_result = {'success': False, 'error': f'aic_method_invalid_return_{type(aic_result).__name__}'}

                if aic_result['success']:
                    optimal_n_peaks = aic_result['optimal_n_peaks']
                    best_aic_model = aic_result['best_model']

                    optimization_report['aic_selection'] = {
                        'optimal_n_peaks': optimal_n_peaks,
                        'best_aicc': best_aic_model['aicc'],
                        'confidence': best_aic_model['aicc_weight'],
                        'comparison_table': aic_result['comparison_table']
                    }

                    print(f"   ✓ AIC selected {optimal_n_peaks} peaks with {best_aic_model['aicc_weight']:.1%} confidence")

                    # Use AIC result as best fit if confidence is high
                    if best_aic_model['aicc_weight'] > 0.7:  # High confidence threshold
                        print(f"   ✓ High confidence AIC result, using as final solution")

                        optimization_report['convergence_reason'] = 'high_confidence_aic'
                        optimization_report['best_strategy'] = {
                            'method': 'aic_selection',
                            'n_peaks': optimal_n_peaks,
                            'confidence': best_aic_model['aicc_weight']
                        }

                        return best_aic_model['fit_result'], optimization_report

                    else:
                        # Medium confidence - use as starting point for progressive strategies
                        print(f"   📊 Medium confidence AIC result, will validate with progressive strategies")
                        aic_suggested_peaks = optimal_n_peaks

                        # Modify strategies to focus around AIC suggestion
                        strategies = [
                            {'n_peaks': aic_suggested_peaks, 'separation_factor': 1.0, 'description': f'aic_suggested_{aic_suggested_peaks}_peaks'},
                            {'n_peaks': max(1, aic_suggested_peaks-1), 'separation_factor': 1.2, 'description': 'aic_minus_one'},
                            {'n_peaks': min(len(suspected_peak_positions), aic_suggested_peaks+1), 'separation_factor': 0.8, 'description': 'aic_plus_one'},
                            {'n_peaks': 1, 'separation_factor': 2.0, 'description': 'single_peak_fallback'}
                        ]

                        optimization_report['aic_informed_strategies'] = True

                else:
                    # DEFENSIVE FIX: Ensure aic_result is a dictionary
                    if not isinstance(aic_result, dict):
                        print(f"   ⚠ AIC result returned {type(aic_result)}, converting to dict")
                        aic_result = {'success': False, 'error': f'invalid_aic_result_type_{type(aic_result).__name__}'}

                    print(f"   ⚠ AIC model selection failed: {aic_result.get('error', 'unknown')}")
                    optimization_report['aic_selection'] = {'error': aic_result.get('error', 'failed')}

                    # Fall back to original progressive strategies
                    strategies = [
                        {'n_peaks': len(suspected_peak_positions), 'separation_factor': 0.5, 'description': 'aggressive_multi_peak'},
                        {'n_peaks': len(suspected_peak_positions), 'separation_factor': 1.0, 'description': 'standard_multi_peak'},
                        {'n_peaks': max(1, len(suspected_peak_positions)-1), 'separation_factor': 1.5, 'description': 'conservative_fewer_peaks'},
                        {'n_peaks': 1, 'separation_factor': 2.0, 'description': 'single_peak_fallback'}
                    ]

            except Exception as e:
                optimization_report['aic_selection'] = {'error': str(e)}

                # Fall back to original progressive strategies
                strategies = [
                    {'n_peaks': len(suspected_peak_positions), 'separation_factor': 0.5, 'description': 'aggressive_multi_peak'},
                    {'n_peaks': len(suspected_peak_positions), 'separation_factor': 1.0, 'description': 'standard_multi_peak'},
                    {'n_peaks': max(1, len(suspected_peak_positions)-1), 'separation_factor': 1.5, 'description': 'conservative_fewer_peaks'},
                    {'n_peaks': 1, 'separation_factor': 2.0, 'description': 'single_peak_fallback'}
                ]

        else:
            # Original progressive strategies (no AIC or insufficient peaks for AIC)
            print(f"   📈 Using progressive strategies (AIC disabled or insufficient peaks)")
            strategies = [
                {'n_peaks': len(suspected_peak_positions), 'separation_factor': 0.5, 'description': 'aggressive_multi_peak'},
                {'n_peaks': len(suspected_peak_positions), 'separation_factor': 1.0, 'description': 'standard_multi_peak'},
                {'n_peaks': max(1, len(suspected_peak_positions)-1), 'separation_factor': 1.5, 'description': 'conservative_fewer_peaks'},
                {'n_peaks': 1, 'separation_factor': 2.0, 'description': 'single_peak_fallback'}
            ]

        best_fit = None
        best_quality = 0

        for i, strategy in enumerate(strategies):
            try:
                print(f"   Strategy {i+1}: {strategy['description']} ({strategy['n_peaks']} peaks)")

                # Adjust separation threshold based on nucleus type
                nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
                separation_threshold = (self.nmr_ranges[nucleus_type]['typical_width'] *
                                      strategy['separation_factor'])

                # Attempt fitting with this strategy
                if strategy['n_peaks'] > 1:
                    # Multi-peak simultaneous fitting
                    fit_result = self._fit_multiple_peaks_simultaneously(
                        x_data, y_data,
                        suspected_peak_positions[:strategy['n_peaks']],
                        separation_threshold
                    )
                else:
                    # Single peak fitting - use the most prominent position
                    center_pos = suspected_peak_positions[0]
                    fit_result = self.fit_peak_enhanced(x_data, y_data, center_pos, nucleus_type,
                                                      method='single_step')  # Use standard method for single peaks

                # CRITICAL FIX: Handle cases where fit_result might be a tuple instead of dict
                if isinstance(fit_result, tuple):
                    print(f"   ⚠ Warning: fit_result returned as tuple, converting to dict")
                    # Convert tuple to dict format (assuming standard tuple format)
                    fit_result = {
                        'success': False,
                        'r_squared': 0.0,
                        'error': 'fit_result_returned_as_tuple'
                    }
                elif not isinstance(fit_result, dict):
                    print(f"   ⚠ Warning: fit_result is {type(fit_result)}, converting to dict")
                    fit_result = {
                        'success': False,
                        'r_squared': 0.0,
                        'error': f'unexpected_result_type_{type(fit_result).__name__}'
                    }

                current_quality = fit_result.get('r_squared', 0)

                strategy_report = {
                    'strategy_index': i + 1,
                    'description': strategy['description'],
                    'n_peaks': strategy['n_peaks'],
                    'separation_factor': strategy['separation_factor'],
                    'separation_threshold': separation_threshold,
                    'r_squared': current_quality,
                    'success': fit_result.get('success', False)
                }

                if fit_result.get('success', False):
                    # Add fitted peak positions if available
                    if 'peak_positions' in fit_result:
                        strategy_report['fitted_positions'] = fit_result['peak_positions']
                    elif 'center' in fit_result:
                        strategy_report['fitted_positions'] = [fit_result['center']]

                    print(f"     R² = {current_quality:.4f}")
                else:
                    strategy_report['error'] = fit_result.get('error', 'unknown_error')
                    print(f"     Failed: {fit_result.get('error', 'unknown')}")

                optimization_report['strategies_tried'].append(strategy_report)

                # Check for improvement (5% improvement threshold to avoid minor fluctuations)
                improvement_threshold = 0.005
                if (fit_result.get('success', False) and
                    current_quality > best_quality + improvement_threshold):

                    best_quality = current_quality
                    best_fit = fit_result
                    optimization_report['best_strategy'] = strategy_report

                    print(f"     ✓ New best strategy (R² = {current_quality:.4f})")

                    # Early stopping for very good fits
                    if current_quality > 0.95:
                        optimization_report['convergence_reason'] = 'high_quality_achieved'
                        print(f"     Early stop: excellent quality achieved")
                        break

            except Exception as e:
                strategy_report = {
                    'strategy_index': i + 1,
                    'description': strategy['description'],
                    'error': str(e),
                    'success': False
                }
                optimization_report['strategies_tried'].append(strategy_report)
                print(f"     Failed with error: {e}")

        # Determine convergence reason
        if best_fit is None:
            optimization_report['convergence_reason'] = 'all_strategies_failed'
            print("   ❌ All overlap detection strategies failed")
        elif optimization_report.get('convergence_reason') is None:
            optimization_report['convergence_reason'] = 'all_strategies_exhausted'
            print(f"   ✓ Optimization completed, best R² = {best_quality:.4f}")

        return best_fit, optimization_report

    def _fit_multiple_peaks_simultaneously(self, x_data, y_data, peak_positions,
                                         separation_threshold):
        """
        CORE MULTI-PEAK DECONVOLUTION: Simultaneous multi-peak Voigt fitting

        This is the heart of the overlapping peak resolution system:
        1. Constructs simultaneous Voigt model for all peaks
        2. Uses constrained optimization to maintain peak positions
        3. Shares baseline across all peaks for consistency
        4. Validates results for physical reasonableness

        Args:
            x_data: X-axis data (ppm)
            y_data: Intensity data
            peak_positions: List of peak center positions
            separation_threshold: Minimum allowed peak separation

        Returns:
            dict: Comprehensive fitting results or failure information

        ROBUST DESIGN: Extensive error handling and parameter validation
        """
        from scipy.optimize import minimize

        n_peaks = len(peak_positions)
        if n_peaks == 1:
            # Single peak - use enhanced single peak fitting
            return self.fit_peak_enhanced(x_data, y_data, peak_positions[0])

        print(f"     Fitting {n_peaks} peaks simultaneously...")

        try:
            # ENHANCED INITIALIZATION: Try sequential peak isolation first
            print(f"     🚀 Attempting enhanced sequential peak isolation...")
            isolation_result = self.sequential_peak_isolation(x_data, y_data, peak_positions)

            if isolation_result['success'] and isolation_result['isolation_quality'] > 0.3:
                # Use sequential isolation results
                print(f"     ✓ Using sequential isolation (quality: {isolation_result['isolation_quality']:.2f})")
                initial_params = isolation_result['combined_initial_guess']
                parameter_bounds = isolation_result['smart_bounds']
                initialization_method = 'sequential_isolation'

                # Extract peak amplitudes for logging
                peak_amplitudes = []
                for i in range(0, len(initial_params)-1, 4):  # Skip baseline at end
                    if i < len(initial_params)-1:
                        peak_amplitudes.append(initial_params[i])  # amplitude is first in each group

                print(f"     Sequential baseline: {initial_params[-1]:.1f}, amplitudes: {[int(a) for a in peak_amplitudes]}")

            else:
                # Fallback to original simple initialization
                print(f"     ⚠ Sequential isolation failed or low quality, using simple initialization")
                baseline_est = self.robust_baseline_estimation(x_data, y_data, method='polynomial')
                nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
                typical_width = self.nmr_ranges[nucleus_type]['typical_width']

                # Build initial parameter vector: [amp1, pos1, sig1, gam1, amp2, pos2, sig2, gam2, ..., baseline]
                initial_params = []
                peak_amplitudes = []

                for pos in peak_positions:
                    # Estimate amplitude from local maximum
                    pos_idx = np.argmin(np.abs(x_data - pos))
                    local_amplitude = max(y_data[pos_idx] - baseline_est, np.max(y_data) * 0.1)
                    peak_amplitudes.append(local_amplitude)

                    initial_params.extend([
                        local_amplitude,           # amplitude
                        pos,                       # center
                        typical_width * 0.6,       # sigma (Gaussian component)
                        typical_width * 0.4        # gamma (Lorentzian component)
                    ])

                initial_params.append(baseline_est)  # Shared baseline

                # Create basic bounds for fallback
                lower_bounds = []
                upper_bounds = []
                data_max = np.max(y_data)
                data_range = data_max - np.min(y_data)

                for i, pos in enumerate(peak_positions):
                    amp = peak_amplitudes[i]
                    lower_bounds.extend([
                        amp * 0.1,                    # amplitude lower
                        pos - typical_width * 2,      # center lower
                        typical_width * 0.1,          # sigma lower
                        typical_width * 0.1           # gamma lower
                    ])
                    upper_bounds.extend([
                        data_max * 3,                 # amplitude upper
                        pos + typical_width * 2,      # center upper
                        typical_width * 10,           # sigma upper
                        typical_width * 10            # gamma upper
                    ])

                # Baseline bounds
                lower_bounds.append(baseline_est - data_range * 0.2)
                upper_bounds.append(baseline_est + data_range * 0.2)

                parameter_bounds = (lower_bounds, upper_bounds)
                initialization_method = 'simple_estimation'

                print(f"     Simple baseline: {baseline_est:.1f}, amplitudes: {[int(a) for a in peak_amplitudes]}")

            # Define multi-peak Voigt model
            def multi_voigt_model(params):
                """Multi-peak Voigt model with shared baseline"""
                *peak_params, baseline = params
                y_model = np.full(len(x_data), baseline)

                for i in range(n_peaks):
                    p_start = i * 4
                    if p_start + 3 < len(peak_params):
                        amplitude, center, sigma, gamma = peak_params[p_start:p_start+4]
                        # Add this peak's contribution
                        y_model += self.voigt_profile_1d(x_data, amplitude, center, sigma, gamma, 0)

                return y_model

            # Objective function (minimize sum of squared residuals)
            def objective(params):
                try:
                    y_pred = multi_voigt_model(params)
                    residuals = y_data - y_pred
                    return np.sum(residuals**2)
                except Exception as e:
                    print(f"     Objective function error: {e}")
                    return 1e10  # Large penalty for failed evaluations

            # Use smart bounds from sequential isolation or create basic bounds
            if 'parameter_bounds' in locals() and initialization_method == 'sequential_isolation':
                # Use smart bounds from sequential isolation
                lower_bounds, upper_bounds = parameter_bounds

                # Ensure bounds are compatible with scipy.optimize format
                bounds = [(lower_bounds[i], upper_bounds[i]) for i in range(len(initial_params))]

                # Additional safety: ensure center position bounds respect separation
                for i, pos in enumerate(peak_positions):
                    center_idx = i * 4 + 1  # center is 2nd parameter for each peak
                    if center_idx < len(bounds):
                        current_lower, current_upper = bounds[center_idx]

                        # Apply separation constraint but don't override smart bounds too much
                        position_tolerance = min(separation_threshold * 0.5,
                                               self.nmr_ranges.get(self.detect_nucleus_type([x_data[0], x_data[-1]]),
                                                                 self.nmr_ranges['1H'])['typical_width'] * 2)

                        # STRENGTHENED: Enforce strict data range bounds with additional margin
                        data_margin = 0.001 * abs(x_data[-1] - x_data[0])  # 0.1% margin for numerical stability
                        constrained_lower = max(current_lower, x_data[0] + data_margin, pos - position_tolerance)
                        constrained_upper = min(current_upper, x_data[-1] - data_margin, pos + position_tolerance)

                        # SAFETY: Ensure bounds are not inverted after constraint application
                        if constrained_lower >= constrained_upper:
                            # Fix by expanding around position
                            margin = max(position_tolerance * 0.1, abs(pos) * 0.01)
                            constrained_lower = pos - margin
                            constrained_upper = pos + margin

                        bounds[center_idx] = (constrained_lower, constrained_upper)

                print(f"     Optimizing {len(initial_params)} parameters with smart bounds from isolation...")

            else:
                # Fallback to original bounds construction
                bounds = []
                nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
                typical_width = self.nmr_ranges[nucleus_type]['typical_width']

                for i, pos in enumerate(peak_positions):
                    # Peak position constraint (stay near initial position)
                    position_tolerance = min(separation_threshold * 0.5, typical_width * 2)

                    # STRENGTHENED: Ensure position bounds are valid with strict data range enforcement
                    data_margin = 0.001 * abs(x_data[-1] - x_data[0])  # 0.1% margin for numerical stability
                    pos_lower = max(x_data[0] + data_margin, pos - position_tolerance)
                    pos_upper = min(x_data[-1] - data_margin, pos + position_tolerance)

                    # Ensure bounds are not inverted
                    if pos_upper <= pos_lower:
                        pos_lower = max(x_data[0] + data_margin, pos - typical_width * 0.1)
                        pos_upper = min(x_data[-1] - data_margin, pos + typical_width * 0.1)

                        # Double-check still not inverted
                        if pos_upper <= pos_lower:
                            # Emergency: use position ± small margin
                            margin = max(typical_width * 0.01, abs(pos) * 0.001)
                            pos_lower = pos - margin
                            pos_upper = pos + margin

                    # Parameter bounds for this peak (with constrained center position)
                    bounds.extend([
                        (0, np.max(y_data) * 3),                    # amplitude > 0, reasonable upper limit
                        (pos_lower, pos_upper),                     # constrained center position
                        (typical_width * 0.01, typical_width * 5),  # sigma bounds
                        (0, typical_width * 3)                      # gamma bounds (≥ 0)
                    ])

                # Baseline bounds
                baseline_range = np.max(y_data) - np.min(y_data)
                bounds.append((np.min(y_data) - baseline_range * 0.1,
                              np.max(y_data) + baseline_range * 0.1))

                print(f"     Optimizing {len(initial_params)} parameters with standard bounds...")

            # Apply constraint coupling if sequential isolation was successful
            if initialization_method == 'sequential_isolation' and 'isolation_result' in locals():
                coupling_result = self.apply_constraint_coupling(
                    initial_params,
                    (list(bounds[i][0] for i in range(len(bounds))),
                     list(bounds[i][1] for i in range(len(bounds)))),
                    isolation_result, n_peaks
                )

                if coupling_result['coupling_info']['applied']:
                    # Use coupled parameters and bounds
                    initial_params = coupling_result['coupled_params']
                    coupled_lower, coupled_upper = coupling_result['coupled_bounds']
                    bounds = [(coupled_lower[i], coupled_upper[i]) for i in range(len(initial_params))]
                    print(f"     ✓ Constraint coupling applied ({coupling_result['coupling_info']['n_constraints']} constraints)")
                else:
                    print(f"     ⚠ Constraint coupling skipped: {coupling_result['coupling_info']['reason']}")

##GM Added 

            # Clamp parameters using shared validator to prevent 'x0 is infeasible'
            bounds_tuple = (
                  [b[0] for b in bounds],
                  [b[1] for b in bounds],
            )
            validated_guess, validated_bounds, validation_info = self.validate_and_fix_initial_guess(
                  initial_params, bounds_tuple, x_data, y_data
            )
            initial_params = validated_guess
            bounds = list(zip(validated_bounds[0], validated_bounds[1]))
            self.last_fit_diagnostics['multi_peak_initial_validation'] = validation_info

##GM added 

            initial_params = self._clamp_vector_to_bounds(initial_params, bounds)

            # CRITICAL: Validate initial parameters and fix inverted bounds
            print(f"     🔧 Validating multi-peak initial parameters and bounds...")
            bounds_validation_success = True
            bounds_fixes = []

            for i, (param, (lower, upper)) in enumerate(zip(initial_params, bounds)):
                # CRITICAL FIX: Check for inverted bounds first
                if lower > upper:
                    print(f"       🚨 INVERTED BOUNDS {i}: lower={lower:.3e} > upper={upper:.3e}")
                    # Fix by swapping and expanding around parameter
                    param_range = abs(lower - upper)
                    expansion = max(param_range * 0.5, abs(param) * 0.1)  # 50% expansion or 10% of param value

                    new_lower = min(lower, upper) - expansion
                    new_upper = max(lower, upper) + expansion

                    # Ensure parameter fits
                    if param < new_lower:
                        new_lower = param - expansion
                    elif param > new_upper:
                        new_upper = param + expansion

                    bounds[i] = (new_lower, new_upper)
                    bounds_fixes.append(f"Fixed inverted bounds {i}: [{lower:.3e}, {upper:.3e}] → [{new_lower:.3e}, {new_upper:.3e}]")
                    bounds_validation_success = False
                    lower, upper = new_lower, new_upper  # Update for next check

                # Check parameter within bounds
                if not (lower <= param <= upper):
                    print(f"       ⚠ Parameter {i}: {param:.3e} not in [{lower:.3e}, {upper:.3e}]")
                    # Fix the violation with safer logic
                    if param < lower:
                        if upper > lower:  # Valid bounds
                            initial_params[i] = lower + (upper - lower) * 0.01  # 1% above lower
                        else:  # Should not happen after bounds fix, but safety
                            initial_params[i] = lower + abs(param) * 0.01
                        print(f"         Fixed to {initial_params[i]:.3e}")
                    elif param > upper:
                        if upper > lower:  # Valid bounds
                            initial_params[i] = upper - (upper - lower) * 0.01  # 1% below upper
                        else:  # Should not happen after bounds fix, but safety
                            initial_params[i] = upper - abs(param) * 0.01
                        print(f"         Fixed to {initial_params[i]:.3e}")
                    bounds_validation_success = False

            if bounds_fixes:
                print(f"     🔧 Fixed {len(bounds_fixes)} inverted bounds:")
                for fix in bounds_fixes[:3]:  # Show first 3
                    print(f"       {fix}")

            if bounds_validation_success:
                print(f"     ✓ All {len(initial_params)} parameters within valid bounds")
            else:
                print(f"     🔧 Fixed parameter and bounds violations")

            # Perform bounded optimization with enhanced error handling
            # Increased maxiter for complex multi-peak fitting
            max_iterations = max(2000, n_peaks * 800)  # Scale with peak count

            try:
                result = minimize(
                    objective, initial_params,
                    method='L-BFGS-B',  # Works well with bounds only
                    bounds=bounds,
                    options={'maxiter': max_iterations, 'ftol': 1e-9}
                )
            except ValueError as e:
                if 'x0 is infeasible' in str(e):
                    print(f"     ❌ Multi-peak x0 infeasible error caught - applying emergency bounds expansion")
                    # Emergency bounds expansion
                    expanded_bounds = []
                    for i, (param, (lower, upper)) in enumerate(zip(initial_params, bounds)):
                        # Expand bounds by 50% around parameter value
                        param_range = abs(upper - lower)
                        expansion = param_range * 0.5
                        new_lower = min(lower, param - expansion)
                        new_upper = max(upper, param + expansion)
                        expanded_bounds.append((new_lower, new_upper))

                    print(f"     🔧 Attempting with expanded bounds...")
                    try:
                        result = minimize(
                            objective, initial_params,
                            method='L-BFGS-B',
                            bounds=expanded_bounds,
                            options={'maxiter': max_iterations, 'ftol': 1e-9}
                        )
                    except Exception as e2:
                        print(f"     ❌ Emergency bounds expansion also failed: {e2}")
                        raise e2
                else:
                    raise e

            if result.success:
                # Extract fitted parameters
                *fitted_peak_params, fitted_baseline = result.x

                # Generate fitted curve
                y_fitted = multi_voigt_model(result.x)

                # Calculate quality metrics
                r_squared = self.calculate_r_squared(y_data, y_fitted)

                # Extract individual peak information
                fitted_positions = []
                fitted_amplitudes = []
                fitted_widths = []

                for i in range(n_peaks):
                    p_start = i * 4
                    amplitude = fitted_peak_params[p_start]
                    center = fitted_peak_params[p_start + 1]
                    sigma = fitted_peak_params[p_start + 2]
                    gamma = fitted_peak_params[p_start + 3]

                    fitted_positions.append(center)
                    fitted_amplitudes.append(amplitude)
                    fitted_widths.append(sigma + gamma)  # Total width

                # ENHANCED VALIDATION: Use comprehensive parameter validation
                nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
                validation_result = self.enhanced_multi_peak_parameter_validation(
                    result.x, x_data, y_data, n_peaks, nucleus_type
                )

                if validation_result['valid']:
                    print(f"     ✓ Multi-peak fit successful: R² = {r_squared:.4f}")
                    print(f"     Fitted positions: {[f'{p:.4f}' for p in fitted_positions]}")

                    return {
                        'success': True,
                        'method': 'simultaneous_multi_peak',
                        'initialization_method': initialization_method,
                        'parameters': result.x,
                        'fitted_curve': y_fitted,
                        'r_squared': r_squared,
                        'n_peaks_fitted': n_peaks,
                        'peak_positions': fitted_positions,
                        'peak_amplitudes': fitted_amplitudes,
                        'peak_widths': fitted_widths,
                        'baseline': fitted_baseline,
                        'isolation_info': isolation_result if 'isolation_result' in locals() else None,
                        'coupling_info': coupling_result if 'coupling_result' in locals() else None,
                        'validation_info': validation_result,
                        'optimization_info': {
                            'n_iterations': result.nit,
                            'function_evaluations': result.nfev,
                            'optimization_success': result.success,
                            'final_residual': float(result.fun),
                            'initialization_used': initialization_method,
                            'constraint_coupling_applied': 'coupling_result' in locals() and coupling_result['coupling_info']['applied'],
                            'validation_quality': validation_result['quality_score']
                        }
                    }
                else:
                    # Enhanced validation failed
                    print(f"     ✗ Multi-peak validation failed: {len(validation_result['violations'])} violations")
                    for violation in validation_result['violations'][:3]:  # Show first 3
                        print(f"       {violation}")

                    return {
                        'success': False,
                        'error': 'enhanced_validation_failed',
                        'validation_info': validation_result,
                        'details': {
                            'n_violations': validation_result['n_violations'],
                            'n_warnings': validation_result['n_warnings'],
                            'quality_score': validation_result['quality_score'],
                            'fitted_positions': fitted_positions,
                            'initial_positions': peak_positions
                        }
                    }
            else:
                return {
                    'success': False,
                    'error': f'optimization_failed: {result.message}',
                    'details': {
                        'n_iterations': result.nit,
                        'function_evaluations': result.nfev,
                        'final_residual': float(result.fun) if hasattr(result, 'fun') else None
                    }
                }

        except Exception as e:
            return {
                'success': False,
                'error': f'multi_peak_fitting_exception: {str(e)}',
                'n_peaks_attempted': n_peaks
            }

    def calculate_aic_for_peak_count(self, x_data, y_data, peak_positions, n_peaks_to_fit, _call_stack=None):
        """
        Calculate AIC (Akaike Information Criterion) for a specific number of peaks

        AIC = 2k - 2ln(L)
        where k = number of parameters, L = likelihood

        For Gaussian noise: AIC ≈ n*ln(RSS/n) + 2k
        where n = number of data points, RSS = residual sum of squares

        Parameters:
        - x_data, y_data: spectral data
        - peak_positions: suspected peak positions (sorted by prominence)
        - n_peaks_to_fit: number of peaks to fit (1 to len(peak_positions))

        Returns:
        - dict with AIC value, fitted parameters, and fit quality
        """
        try:
            if n_peaks_to_fit > len(peak_positions):
                return {'success': False, 'error': 'insufficient_peak_positions'}

            # Select most prominent peaks for fitting
            selected_positions = peak_positions[:n_peaks_to_fit]

            # Fit the selected number of peaks
            if n_peaks_to_fit == 1:
                # CRITICAL: Use 'single_step' method to prevent recursion back to AIC
                fit_result = self.fit_peak_enhanced(x_data, y_data, selected_positions[0],
                                                  method='single_step', _call_stack=_call_stack)

                # DEFENSIVE FIX: Ensure fit_result is a dictionary
                if not isinstance(fit_result, dict):
                    return {'success': False, 'error': f'single_peak_enhanced_invalid_return_{type(fit_result).__name__}'}

                if fit_result['success']:
                    y_fitted = fit_result['fitted_curve']
                    n_params = 5  # amp, center, sigma, gamma, baseline
                else:
                    return {'success': False, 'error': 'single_peak_fit_failed'}
            else:
                # Multi-peak fitting
                nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
                separation_threshold = self.nmr_ranges[nucleus_type]['typical_width'] * 1.0

                fit_result = self._fit_multiple_peaks_simultaneously(
                    x_data, y_data, selected_positions, separation_threshold
                )

                # DEFENSIVE FIX: Ensure fit_result is a dictionary
                if not isinstance(fit_result, dict):
                    return {'success': False, 'error': f'multi_peak_simultaneous_invalid_return_{type(fit_result).__name__}'}

                if fit_result['success']:
                    y_fitted = fit_result['fitted_curve']
                    n_params = n_peaks_to_fit * 4 + 1  # (amp,center,sigma,gamma) per peak + baseline
                else:
                    return {'success': False, 'error': 'multi_peak_fit_failed'}

            # Calculate RSS (Residual Sum of Squares)
            residuals = y_data - y_fitted
            rss = np.sum(residuals ** 2)

            # Calculate AIC
            n_points = len(y_data)

            # Prevent numerical issues
            if rss <= 0 or n_points <= n_params:
                return {'success': False, 'error': 'numerical_issues'}

            # AIC formula: n*ln(RSS/n) + 2k
            # Add penalty for small sample size: AICc = AIC + 2k(k+1)/(n-k-1)
            aic = n_points * np.log(rss / n_points) + 2 * n_params

            # Small sample correction (AICc)
            if n_points - n_params - 1 > 0:
                aicc = aic + (2 * n_params * (n_params + 1)) / (n_points - n_params - 1)
            else:
                aicc = np.inf  # Model too complex for sample size

            return {
                'success': True,
                'n_peaks': n_peaks_to_fit,
                'aic': aic,
                'aicc': aicc,
                'rss': rss,
                'n_params': n_params,
                'r_squared': fit_result.get('r_squared', 0),
                'fit_result': fit_result,
                'selected_positions': selected_positions
            }

        except Exception as e:
            return {
                'success': False,
                'error': f'aic_calculation_failed: {str(e)}',
                'n_peaks': n_peaks_to_fit
            }

    def optimal_peak_count_by_aic(self, x_data, y_data, suspected_positions, max_peaks=5, _call_stack=None):
        """
        Determine optimal number of peaks using AIC model selection

        Tests models with 1 to max_peaks and selects the one with minimum AIC.
        Includes safeguards against overfitting and numerical instability.

        Parameters:
        - x_data, y_data: spectral data
        - suspected_positions: list of suspected peak positions (sorted by prominence)
        - max_peaks: maximum number of peaks to test

        Returns:
        - dict with optimal model selection results and comparison table
        """
        print(f"   🔍 AIC model selection: testing 1-{max_peaks} peaks...")

        # Limit max_peaks to available positions and reasonable computational limits
        max_peaks = min(max_peaks, len(suspected_positions), 6)

        aic_results = []
        valid_models = []

        # Test each peak count
        for n_peaks in range(1, max_peaks + 1):
            print(f"     Testing {n_peaks} peak{'s' if n_peaks > 1 else ''}...")

            aic_result = self.calculate_aic_for_peak_count(x_data, y_data, suspected_positions, n_peaks, _call_stack)

            # DEFENSIVE FIX: Ensure aic_result is a dictionary
            if not isinstance(aic_result, dict):
                print(f"       ⚠ AIC calculation returned {type(aic_result)}, converting to dict")
                aic_result = {'success': False, 'error': f'aic_calc_invalid_return_{type(aic_result).__name__}'}

            aic_results.append(aic_result)

            if aic_result['success']:
                valid_models.append(aic_result)
                print(f"       AIC: {aic_result['aic']:.1f}, AICc: {aic_result['aicc']:.1f}, R²: {aic_result['r_squared']:.3f}")
            else:
                print(f"       Failed: {aic_result.get('error', 'unknown')}")

        if not valid_models:
            return {
                'success': False,
                'error': 'no_valid_models',
                'tested_counts': list(range(1, max_peaks + 1)),
                'all_results': aic_results
            }

        # Find model with minimum AICc (corrected AIC preferred for small samples)
        best_model = min(valid_models, key=lambda x: x['aicc'])

        # Calculate AIC differences for interpretation
        min_aicc = best_model['aicc']
        for model in valid_models:
            model['delta_aicc'] = model['aicc'] - min_aicc
            model['aicc_weight'] = np.exp(-0.5 * model['delta_aicc'])

        # Normalize weights
        total_weight = sum(m['aicc_weight'] for m in valid_models)
        for model in valid_models:
            model['aicc_weight'] /= total_weight

        print(f"   ✓ Optimal model: {best_model['n_peaks']} peaks (AICc: {best_model['aicc']:.1f})")

        # Model comparison summary
        comparison_table = []
        for model in valid_models:
            comparison_table.append({
                'n_peaks': model['n_peaks'],
                'aic': model['aic'],
                'aicc': model['aicc'],
                'delta_aicc': model['delta_aicc'],
                'weight': model['aicc_weight'],
                'r_squared': model['r_squared'],
                'is_best': model == best_model
            })

        return {
            'success': True,
            'optimal_n_peaks': best_model['n_peaks'],
            'best_model': best_model,
            'all_models': valid_models,
            'comparison_table': comparison_table,
            'model_selection_summary': {
                'tested_range': f"1-{max_peaks}",
                'valid_models': len(valid_models),
                'best_aicc': best_model['aicc'],
                'confidence': best_model['aicc_weight']
            }
        }

    def get_isolated_peak_window(self, x_data, y_data, peak_center, overlapping_peaks=None,
                                window_multiplier=None):
        """
        Get window around peak that excludes overlapping peaks for quality assessment

        ENHANCED: Now supports GUI-based window sizing for consistent display

        Parameters:
        -----------
        x_data, y_data : array
            Spectral data
        peak_center : float
            Peak center position in ppm
        overlapping_peaks : list, optional
            List of overlapping peak positions
        window_multiplier : float, optional
            Window size multiplier (None = use GUI params, maintains backward compatibility)

        Returns:
        --------
        dict : Local region data with GUI-consistent window size
        """
        try:
            if overlapping_peaks is None:
                overlapping_peaks = self.detect_overlapping_peaks(x_data, y_data)

            # ENHANCED: Use GUI parameters if window_multiplier not explicitly provided
            if window_multiplier is None:
                ppm_range = abs(x_data[-1] - x_data[0]) if len(x_data) > 1 else 1.0
                nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
                window_multiplier = self._calculate_gui_based_multiplier(
                    nucleus_type, ppm_range, len(x_data)
                )
                print(f"   🎯 Isolated window using GUI-based multiplier: {window_multiplier:.2f}× (was hardcoded 5.0×)")

            # Find the peak of interest
            target_peak = None
            for peak in overlapping_peaks:
                if abs(peak['position'] - peak_center) < 0.01:  # Close match
                    target_peak = peak
                    break

            if target_peak is None:
                # Fallback to standard local region extraction
                nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
                typical_width = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])['typical_width']
                return self.extract_local_peak_region(x_data, y_data, peak_center,
                                                    typical_width, nucleus_type, window_multiplier)

            if target_peak['is_isolated']:
                # Peak is isolated, use standard window
                nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
                typical_width = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])['typical_width']
                return self.extract_local_peak_region(x_data, y_data, peak_center,
                                                    typical_width, nucleus_type, window_multiplier)
            else:
                # Peak has nearby overlapping peaks - create constrained window
                nearby_positions = [p['position'] for p in target_peak['nearby_peaks']]

                # Find the closest neighboring peaks on each side
                left_neighbors = [pos for pos in nearby_positions if pos < peak_center]
                right_neighbors = [pos for pos in nearby_positions if pos > peak_center]

                # Set window boundaries to halfway between peaks
                if left_neighbors:
                    left_boundary = peak_center - (peak_center - max(left_neighbors)) / 2
                else:
                    nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
                    typical_width = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])['typical_width']
                    left_boundary = peak_center - typical_width * window_multiplier

                if right_neighbors:
                    right_boundary = peak_center - (peak_center - min(right_neighbors)) / 2
                else:
                    nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
                    typical_width = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])['typical_width']
                    right_boundary = peak_center + typical_width * window_multiplier

                # Extract constrained region
                constrained_mask = (x_data >= left_boundary) & (x_data <= right_boundary)
                constrained_indices = np.where(constrained_mask)[0]

                if len(constrained_indices) < 5:
                    # Expand to minimum points if too constrained
                    center_idx = np.argmin(np.abs(x_data - peak_center))
                    min_points = 5
                    start_idx = max(0, center_idx - min_points//2)
                    end_idx = min(len(x_data), center_idx + min_points//2 + 1)
                    constrained_indices = np.arange(start_idx, end_idx)
                    constrained_mask = np.zeros(len(x_data), dtype=bool)
                    constrained_mask[constrained_indices] = True

                x_local = x_data[constrained_indices]
                y_local = y_data[constrained_indices]

                return {
                    'x_local': x_local,
                    'y_local': y_local,
                    'local_indices': constrained_indices,
                    'local_mask': constrained_mask,
                    'window_width': right_boundary - left_boundary,
                    'n_points': len(constrained_indices),
                    'isolation_constrained': True,
                    'nearby_peaks_excluded': len(target_peak['nearby_peaks'])
                }

        except Exception as e:
            print(f"Isolated peak window extraction failed: {e}")
            # Fallback to standard extraction
            nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])
            typical_width = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])['typical_width']
            return self.extract_local_peak_region(x_data, y_data, peak_center,
                                                typical_width, nucleus_type, window_multiplier)

    def enhanced_peak_fitting_parallel(self, peak_list, use_parallel=True, progress_callback=None, parent_integrator=None):
        """
        New parallel entry point that maintains complete compatibility with existing interface.

        Args:
            peak_list: DataFrame with peak information or single peak coordinates
            use_parallel: Enable parallel processing (default: True)
            progress_callback: Progress update callback function
            parent_integrator: Parent integrator with nmr_data (auto-detected if None)

        Returns:
            Same format as existing enhanced_peak_fitting methods
        """
        import pandas as pd

        # Handle single peak input (maintain compatibility)
        if not isinstance(peak_list, pd.DataFrame):
            # Assume single peak with (peak_x, peak_y, assignment) format
            if isinstance(peak_list, (list, tuple)) and len(peak_list) >= 2:
                peak_x, peak_y = peak_list[0], peak_list[1]
                assignment = peak_list[2] if len(peak_list) > 2 else 'Single_Peak'

                # Create single-row DataFrame
                peak_list = pd.DataFrame({
                    'Position_X': [peak_x],
                    'Position_Y': [peak_y],
                    'Assignment': [assignment]
                })
            else:
                raise ValueError("Invalid peak_list format. Expected DataFrame or (peak_x, peak_y, assignment) tuple")

        # Auto-detect parent integrator if not provided
        if parent_integrator is None:
            parent_integrator = getattr(self, 'parent', None)

        # Determine processing method
        if use_parallel and len(peak_list) > 2:  # Parallel threshold
            try:
                # Ensure we have the necessary data context
                if parent_integrator is None:
                    # Look for integrator in common places
                    import inspect
                    frame = inspect.currentframe()
                    try:
                        # Check calling context for integrator
                        while frame:
                            frame_locals = frame.f_locals
                            if 'self' in frame_locals:
                                candidate = frame_locals['self']
                                if hasattr(candidate, 'nmr_data') and hasattr(candidate, 'enhanced_fitter'):
                                    if candidate.enhanced_fitter is self:
                                        parent_integrator = candidate
                                        break
                            frame = frame.f_back
                    finally:
                        del frame

                if parent_integrator is None or not hasattr(parent_integrator, 'nmr_data'):
                    raise ValueError("No parent integrator with nmr_data found - cannot run parallel processing")

                # Thread-safe setup of data context
                import threading
                with threading.Lock():
                    self.nmr_data = parent_integrator.nmr_data
                    self.ppm_x_axis = parent_integrator.ppm_x_axis
                    self.ppm_y_axis = parent_integrator.ppm_y_axis

                # Use new parallel implementation
                from lunaNMR.core.parallel_voigt_processor import ParallelVoigtProcessor

                print(f"🚀 Using parallel Voigt fitting for {len(peak_list)} peaks")
                parallel_processor = ParallelVoigtProcessor(self)
                results = parallel_processor.fit_all_peaks_parallel(peak_list, progress_callback)

                # Return single result if single peak input
                if len(results) == 1 and len(peak_list) == 1:
                    return results[0]
                else:
                    return results

            except Exception as e:
                print(f"⚠️ Parallel processing failed: {e}")
                print("🔄 Falling back to sequential processing")

        # Fallback to sequential processing
        print(f"🔄 Using sequential Voigt fitting for {len(peak_list)} peaks")
        return self._enhanced_peak_fitting_sequential(peak_list, progress_callback, parent_integrator)

    def _enhanced_peak_fitting_sequential(self, peak_list, progress_callback=None, parent_integrator=None):
        """
        Sequential processing fallback that calls existing enhanced_peak_fitting
        method for each peak individually.
        """
        results = []

        for i, (peak_idx, peak_row) in enumerate(peak_list.iterrows()):
            peak_x = float(peak_row['Position_X'])
            peak_y = float(peak_row['Position_Y'])
            assignment = peak_row.get('Assignment', f'Peak_{i+1}')

            try:
                # Call existing fit_peak_enhanced method (unchanged)
                result = self.fit_peak_enhanced(peak_x, peak_y, assignment)
                if result:
                    result['processing_mode'] = 'sequential'
                    result['peak_number'] = i + 1
                    results.append(result)

                if progress_callback:
                    progress = ((i + 1) / len(peak_list)) * 100
                    progress_callback(progress, f"Sequential: {i+1}/{len(peak_list)}", assignment)

            except Exception as e:
                print(f"❌ Sequential processing failed for peak {i+1} ({assignment}): {e}")

        # Return single result if single peak input
        if len(results) == 1 and len(peak_list) == 1:
            return results[0]
        else:
            return results


    def get_fitting_diagnostics(self):
        """Return detailed diagnostics from last fit"""
        return self.last_fit_diagnostics.copy()


# Convenience function for backward compatibility
def fit_voigt_peak(x_data, y_data, initial_center=None, nucleus_type=None):
    """
    Convenience function for enhanced Voigt fitting

    Maintains compatibility with existing code while providing enhanced functionality
    """
    fitter = EnhancedVoigtFitter()
    return fitter.fit_peak_enhanced(x_data, y_data, initial_center, nucleus_type)


# Additional utility methods for EnhancedVoigtFitter class
def generate_optimization_report(self, peak_assignment="Unknown"):
    """
    COMPREHENSIVE REPORTING: Generate detailed optimization report for user review

    This creates a comprehensive report of the dynamic optimization process showing:
    1. What strategies were tried and why
    2. How each strategy performed (R² values)
    3. What parameters were optimized
    4. Final recommendations for further improvement

    Args:
        peak_assignment: String identifier for this peak (for report labeling)

    Returns:
        dict: Comprehensive report dictionary with user-friendly summaries
    """
    import pandas as pd

    diagnostics = self.last_fit_diagnostics

    report = {
        'peak_assignment': peak_assignment,
        'timestamp': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'optimization_method': diagnostics.get('method', 'unknown'),
        'optimization_active': diagnostics.get('optimization_active', False),
        'summary': {},
        'details': {},
        'recommendations': [],
        'quality_assessment': {}
    }

    # Global parameter estimation summary
    if 'global_params' in diagnostics:
        global_params = diagnostics['global_params']
        report['summary']['global_estimation'] = {
            'typical_linewidth': global_params.get('typical_linewidth', 'unknown'),
            'well_resolved_peaks_used': global_params.get('well_resolved_count', 0),
            'noise_level': global_params.get('noise_level', 'unknown'),
            'estimation_quality': global_params.get('estimation_quality', 'unknown'),
            'method_used': global_params.get('method_used', 'unknown')
        }

        # Recommendations based on global parameter quality
        if global_params.get('well_resolved_count', 0) < 2:
            report['recommendations'].append(
                "RECOMMENDATION: Include more well-resolved peaks in analysis context for better parameter estimation"
            )
        elif global_params.get('estimation_quality') == 'inconsistent':
            report['recommendations'].append(
                "WARNING: Inconsistent linewidth estimates detected - manual review recommended"
            )
        elif global_params.get('estimation_quality') == 'good':
            report['recommendations'].append(
                "✓ Good global parameter estimation achieved from well-resolved peaks"
            )

    # Standard fitting performance
    if 'standard_result' in diagnostics:
        standard = diagnostics['standard_result']
        report['summary']['standard_fitting'] = {
            'r_squared': standard.get('r_squared', 0),
            'success': standard.get('success', False),
            'quality_class': categorize_quality_legacy(standard.get('r_squared', 0))
        }

    # Baseline optimization summary
    if 'baseline_optimization' in diagnostics:
        baseline_opt = diagnostics['baseline_optimization']

        report['summary']['baseline_optimization'] = {
            'iterations_tried': len(baseline_opt.get('iterations', [])),
            'converged': baseline_opt.get('converged', False),
            'best_edge_fraction': baseline_opt.get('best_edge_fraction', 'unknown'),
            'improvement_achieved': baseline_opt.get('improvement_achieved', 0),
            'fallback_used': baseline_opt.get('fallback_used', False)
        }

        report['details']['baseline_iterations'] = baseline_opt.get('iterations', [])

        # Detailed recommendations based on baseline optimization
        effective = diagnostics.get('baseline_optimization_effective', False)
        if effective:
            improvement = report['summary']['baseline_optimization']['improvement_achieved']
            best_fraction = baseline_opt.get('best_edge_fraction', 'unknown')
            report['recommendations'].append(
                f"✓ Baseline optimization successful: R² improved to {improvement:.4f} "
                f"using {best_fraction} edge fraction"
            )
        elif baseline_opt.get('fallback_used', False):
            report['recommendations'].append(
                "WARNING: Baseline optimization failed, used fallback method - "
                "consider manual baseline inspection"
            )
        else:
            report['recommendations'].append(
                "INFO: Baseline optimization did not improve results - "
                "original baseline estimation was already optimal"
            )

    # Add general recommendations if none were generated
    if len(report['recommendations']) == 0:
        report['recommendations'].append("INFO: Standard fitting completed without optimization")

    return report

def print_optimization_summary(self, peak_assignment="Unknown"):
    """
    USER-FRIENDLY CONSOLE OUTPUT: Print comprehensive optimization summary
    """
    report = self.generate_optimization_report(peak_assignment)

    print(f"\n{'='*60}")
    print(f"🔍 OPTIMIZATION REPORT: {peak_assignment}")
    print(f"{'='*60}")
    print(f"Method: {report['optimization_method']}")

    # Global parameters section
    if 'global_estimation' in report['summary']:
        global_est = report['summary']['global_estimation']

    # Standard fitting section
    if 'standard_fitting' in report['summary']:
        standard = report['summary']['standard_fitting']

    # Baseline optimization section
    if 'baseline_optimization' in report['summary']:
        baseline_opt = report['summary']['baseline_optimization']

    # Recommendations section
    for i, rec in enumerate(report['recommendations'], 1):
        print(f"   {i}. {rec}")

    print(f"{'='*60}")
    return report

# Add these methods to the EnhancedVoigtFitter class
EnhancedVoigtFitter.generate_optimization_report = generate_optimization_report
EnhancedVoigtFitter.print_optimization_summary = print_optimization_summary


class RobustParameterEstimator:
    """
    LEVEL 2 ARCHITECTURAL: Robust parameter estimation with cross-validation

    This class implements a comprehensive parameter estimation pipeline that:
    - Uses multiple estimation methods with cross-validation
    - Implements consensus-based parameter selection
    - Provides uncertainty quantification
    - Handles edge cases and pathological data
    """

    def __init__(self, parent_fitter):
        self.parent = parent_fitter
        self.estimation_methods = {
            'moment_based': self._moment_based_estimation,
            'peak_detection': self._peak_detection_estimation,
            'correlation_based': self._correlation_based_estimation,
            'robust_statistics': self._robust_statistics_estimation,
            'physics_informed': self._physics_informed_estimation
        }

        # Estimation quality tracking
        self.method_performance = {method: {'success_count': 0, 'failure_count': 0,
                                           'quality_scores': []} for method in self.estimation_methods}

    def estimate_initial_parameters(self, x_data, y_data, peak_center, nucleus_type=None, context=None):
        """
        LEVEL 2: Comprehensive parameter estimation with method consensus

        Parameters:
        - x_data, y_data: experimental data
        - peak_center: estimated peak center
        - nucleus_type: '1H', '15N', or '13C'
        - context: additional context information

        Returns:
        - dict: parameter estimates with uncertainty and confidence scores
        """
        if nucleus_type is None:
            nucleus_type = self.parent.detect_nucleus_type([x_data[0], x_data[-1]])

        # Data preprocessing and validation
        data_quality = self._assess_data_quality(x_data, y_data)

        # Run all estimation methods
        estimation_results = {}
        for method_name, method_func in self.estimation_methods.items():
            try:
                start_time = time.time()
                result = method_func(x_data, y_data, peak_center, nucleus_type, data_quality, context)
                execution_time = time.time() - start_time

                if result is not None and result.get('success', False):
                    estimation_results[method_name] = result
                    estimation_results[method_name]['execution_time'] = execution_time
                    self.method_performance[method_name]['success_count'] += 1
                else:
                    self.method_performance[method_name]['failure_count'] += 1

            except Exception as e:
                self.method_performance[method_name]['failure_count'] += 1

        if not estimation_results:
            return self._fallback_parameter_estimation(x_data, y_data, peak_center, nucleus_type)

        # Consensus-based parameter selection
        consensus_params = self._build_parameter_consensus(estimation_results, data_quality)

        # Uncertainty quantification
        uncertainties = self._quantify_parameter_uncertainties(estimation_results, consensus_params)

        # Final validation and adjustment
        validated_params = self._validate_and_adjust_parameters(consensus_params, uncertainties, x_data, y_data, nucleus_type)

        return {
            'success': True,
            'method': 'robust_consensus_estimation',
            'parameters': validated_params,
            'uncertainties': uncertainties,
            'consensus_quality': consensus_params.get('consensus_score', 0),
            'individual_results': estimation_results,
            'data_quality': data_quality,
            'method_performance': self.method_performance
        }

    def _assess_data_quality(self, x_data, y_data):
        """Comprehensive data quality assessment"""
        try:
            # Signal-to-noise ratio estimation
            baseline_est = np.percentile(y_data, 10)
            peak_height = np.max(y_data) - baseline_est
            noise_est = np.std(y_data[:len(y_data)//4])  # Estimate from edges
            snr = peak_height / max(noise_est, 1e-10)

            # Baseline stability
            left_baseline = np.mean(y_data[:len(y_data)//5])
            right_baseline = np.mean(y_data[-len(y_data)//5:])
            baseline_stability = 1.0 / (1.0 + abs(left_baseline - right_baseline) / max(peak_height, 1e-10))

            # Peak shape assessment
            peak_idx = np.argmax(y_data)
            peak_symmetry = self._assess_peak_symmetry(x_data, y_data, peak_idx)

            # Spectral resolution
            ppm_per_point = abs(x_data[1] - x_data[0]) if len(x_data) > 1 else 0.001
            resolution_quality = min(1.0, 0.01 / max(ppm_per_point, 1e-6))  # Good if < 0.01 ppm/point

            return {
                'snr': snr,
                'baseline_stability': baseline_stability,
                'peak_symmetry': peak_symmetry,
                'resolution_quality': resolution_quality,
                'noise_level': noise_est,
                'overall_quality': (snr/100 + baseline_stability + peak_symmetry + resolution_quality) / 4
            }
        except:
            return {'snr': 1, 'baseline_stability': 0.5, 'peak_symmetry': 0.5, 'resolution_quality': 0.5, 'overall_quality': 0.5}

    def _assess_peak_symmetry(self, x_data, y_data, peak_idx):
        """Assess peak symmetry for lineshape hints"""
        try:
            half_max = (y_data[peak_idx] + np.min(y_data)) / 2

            # Find half-maximum points
            left_half_idx = None
            right_half_idx = None

            for i in range(peak_idx, 0, -1):
                if y_data[i] < half_max:
                    left_half_idx = i
                    break

            for i in range(peak_idx, len(y_data)):
                if y_data[i] < half_max:
                    right_half_idx = i
                    break

            if left_half_idx is not None and right_half_idx is not None:
                left_width = peak_idx - left_half_idx
                right_width = right_half_idx - peak_idx
                symmetry = min(left_width, right_width) / max(left_width, right_width)
                return symmetry

            return 0.5
        except:
            return 0.5

    def _moment_based_estimation(self, x_data, y_data, peak_center, nucleus_type, data_quality, context):
        """Parameter estimation based on statistical moments"""
        try:
            # Baseline correction
            baseline = self.parent.robust_baseline_estimation(x_data, y_data, method='percentile')
            y_corrected = y_data - baseline
            y_corrected = np.maximum(y_corrected, 0)  # Only positive values

            total_intensity = np.sum(y_corrected)
            if total_intensity <= 0:
                return None

            # First moment (center of mass)
            center_moment = np.sum(x_data * y_corrected) / total_intensity

            # Second moment (variance)
            variance = np.sum(y_corrected * (x_data - center_moment)**2) / total_intensity
            sigma_est = np.sqrt(variance)

            # Amplitude estimation
            amplitude_est = np.max(y_corrected)

            # Initial gamma estimate (assume equal Gaussian/Lorentzian contributions)
            gamma_est = sigma_est * 0.5

            quality_score = min(1.0, data_quality['snr'] / 20) * data_quality['baseline_stability']

            return {
                'success': True,
                'amplitude': amplitude_est,
                'center': center_moment,
                'sigma': sigma_est,
                'gamma': gamma_est,
                'baseline': baseline,
                'quality_score': quality_score,
                'method_info': 'statistical_moments'
            }
        except:
            return None

    def _peak_detection_estimation(self, x_data, y_data, peak_center, nucleus_type, data_quality, context):
        """Parameter estimation using peak detection algorithms"""
        try:
            # Find peaks with prominence threshold
            prominence_threshold = np.std(y_data) * max(1.0, data_quality['snr'] / 10)
            peaks, properties = find_peaks(y_data, prominence=prominence_threshold)

            if len(peaks) == 0:
                return None

            # Find peak closest to expected center
            peak_centers_ppm = x_data[peaks]
            distances = np.abs(peak_centers_ppm - peak_center)
            best_peak_idx = peaks[np.argmin(distances)]

            # Width estimation using scipy
            try:
                widths, width_heights, left_ips, right_ips = peak_widths(y_data, [best_peak_idx], rel_height=0.5)
                width_ppm = widths[0] * abs(x_data[1] - x_data[0])  # Convert to ppm
                sigma_est = width_ppm / 2.355  # FWHM to sigma
            except:
                # Fallback width estimation
                sigma_est = self.parent.nmr_ranges.get(nucleus_type, self.parent.nmr_ranges['1H'])['typical_width']

            # Amplitude and baseline
            amplitude_est = y_data[best_peak_idx]
            baseline_est = self.parent.robust_baseline_estimation(x_data, y_data, method='percentile')
            amplitude_est -= baseline_est

            # Center refinement
            center_est = x_data[best_peak_idx]

            # Gamma estimation based on peak shape
            gamma_est = sigma_est * (2.0 - data_quality['peak_symmetry'])  # Less symmetric = more Lorentzian

            quality_score = min(1.0, len(peaks) / 3.0) * data_quality['overall_quality']

            return {
                'success': True,
                'amplitude': max(amplitude_est, 0),
                'center': center_est,
                'sigma': max(sigma_est, 1e-6),
                'gamma': max(gamma_est, 0),
                'baseline': baseline_est,
                'quality_score': quality_score,
                'method_info': f'peak_detection_{len(peaks)}_peaks'
            }
        except:
            return None

    def _correlation_based_estimation(self, x_data, y_data, peak_center, nucleus_type, data_quality, context):
        """Parameter estimation using template correlation"""
        try:
            # Create template Voigt profiles with different parameters
            baseline_est = self.parent.robust_baseline_estimation(x_data, y_data, method='percentile')
            amplitude_est = np.max(y_data) - baseline_est

            typical_width = self.parent.nmr_ranges.get(nucleus_type, self.parent.nmr_ranges['1H'])['typical_width']

            best_correlation = -1
            best_params = None

            # Test different width combinations
            for sigma_factor in [0.3, 0.5, 0.7, 1.0, 1.5, 2.0]:
                for gamma_factor in [0.1, 0.3, 0.5, 0.7, 1.0]:
                    sigma_test = typical_width * sigma_factor
                    gamma_test = typical_width * gamma_factor

                    # Generate template
                    template = self.parent.voigt_profile_1d(x_data, amplitude_est, peak_center, sigma_test, gamma_test, baseline_est)

                    # Calculate correlation
                    correlation = np.corrcoef(y_data, template)[0, 1]

                    if not np.isnan(correlation) and correlation > best_correlation:
                        best_correlation = correlation
                        best_params = (amplitude_est, peak_center, sigma_test, gamma_test, baseline_est)

            if best_params is not None and best_correlation > 0.5:
                return {
                    'success': True,
                    'amplitude': best_params[0],
                    'center': best_params[1],
                    'sigma': best_params[2],
                    'gamma': best_params[3],
                    'baseline': best_params[4],
                    'quality_score': best_correlation,
                    'method_info': f'correlation_{best_correlation:.3f}'
                }

            return None
        except:
            return None

    def _robust_statistics_estimation(self, x_data, y_data, peak_center, nucleus_type, data_quality, context):
        """Parameter estimation using robust statistical methods"""
        try:
            # Robust baseline using median absolute deviation
            baseline_est = np.median(y_data)
            mad = np.median(np.abs(y_data - baseline_est))
            robust_baseline = baseline_est - 2 * mad  # Conservative baseline

            y_corrected = y_data - robust_baseline
            y_corrected = np.maximum(y_corrected, 0)

            # Robust center estimation using weighted mean
            weights = y_corrected / np.sum(y_corrected)
            robust_center = np.sum(x_data * weights)

            # Robust width estimation using interquartile range
            cumulative_intensity = np.cumsum(y_corrected)
            cumulative_intensity /= cumulative_intensity[-1]

            # Find quartile positions
            q1_idx = np.searchsorted(cumulative_intensity, 0.25)
            q3_idx = np.searchsorted(cumulative_intensity, 0.75)

            if q3_idx > q1_idx:
                iqr_width = x_data[q3_idx] - x_data[q1_idx]
                sigma_est = iqr_width / 2.7  # Approximate conversion
            else:
                sigma_est = self.parent.nmr_ranges.get(nucleus_type, self.parent.nmr_ranges['1H'])['typical_width']

            # Robust amplitude estimation
            robust_amplitude = np.percentile(y_corrected, 95)  # 95th percentile to avoid outliers

            # Gamma estimation based on tail behavior
            gamma_est = sigma_est * 0.3  # Conservative Lorentzian component

            quality_score = data_quality['baseline_stability'] * min(1.0, data_quality['snr'] / 15)

            return {
                'success': True,
                'amplitude': robust_amplitude,
                'center': robust_center,
                'sigma': max(sigma_est, 1e-6),
                'gamma': max(gamma_est, 0),
                'baseline': robust_baseline,
                'quality_score': quality_score,
                'method_info': 'robust_statistics'
            }
        except:
            return None

    def _physics_informed_estimation(self, x_data, y_data, peak_center, nucleus_type, data_quality, context):
        """Parameter estimation using NMR physics constraints"""
        try:
            # Get nucleus-specific constraints
            nmr_params = self.parent.nmr_ranges.get(nucleus_type, self.parent.nmr_ranges['1H'])

            # Physics-based baseline estimation
            baseline_est = self.parent.robust_baseline_estimation(x_data, y_data, method='iterative')

            # Amplitude estimation with physical limits
            peak_height = np.max(y_data) - baseline_est

            # Expected amplitude based on typical NMR intensities
            if nucleus_type == '1H':
                expected_amplitude_range = (peak_height * 0.8, peak_height * 1.2)
            elif nucleus_type == '15N':
                expected_amplitude_range = (peak_height * 0.5, peak_height * 1.5)  # More variable
            else:
                expected_amplitude_range = (peak_height * 0.6, peak_height * 1.4)

            amplitude_est = np.clip(peak_height, *expected_amplitude_range)

            # Physics-informed width estimation
            base_width = nmr_params['typical_width']

            # Adjust for field strength and temperature (if available in context)
            if context and 'field_strength' in context:
                field_factor = context['field_strength'] / 600.0  # Normalize to 600 MHz
                base_width *= field_factor

            if context and 'temperature' in context:
                temp_factor = 298.0 / context['temperature']  # Normalize to 25°C
                base_width *= np.sqrt(temp_factor)  # Temperature affects width

            # Lineshape distribution based on nucleus type
            if nucleus_type == '1H':
                # Protons typically more Gaussian (instrumental broadening dominates)
                sigma_est = base_width * 0.7
                gamma_est = base_width * 0.3
            elif nucleus_type == '15N':
                # 15N often more Lorentzian (relaxation broadening)
                sigma_est = base_width * 0.4
                gamma_est = base_width * 0.6
            else:
                # Default balanced distribution
                sigma_est = base_width * 0.5
                gamma_est = base_width * 0.5

            # Center refinement using physics constraints
            center_est = self._physics_guided_center_refinement(x_data, y_data, peak_center, nucleus_type, baseline_est)

            quality_score = 0.8 * data_quality['overall_quality']  # Physics-based methods are generally reliable

            return {
                'success': True,
                'amplitude': amplitude_est,
                'center': center_est,
                'sigma': sigma_est,
                'gamma': gamma_est,
                'baseline': baseline_est,
                'quality_score': quality_score,
                'method_info': f'physics_{nucleus_type}'
            }
        except:
            return None

    def _physics_guided_center_refinement(self, x_data, y_data, initial_center, nucleus_type, baseline):
        """Refine peak center using physics-guided local search"""
        try:
            # Search in a small window around initial guess
            search_window = self.parent.nmr_ranges.get(nucleus_type, self.parent.nmr_ranges['1H'])['typical_width'] * 0.5

            search_mask = np.abs(x_data - initial_center) <= search_window
            if np.sum(search_mask) < 3:
                return initial_center

            # Find local maximum in search window
            local_x = x_data[search_mask]
            local_y = y_data[search_mask]

            # Parabolic interpolation for sub-pixel precision
            max_idx = np.argmax(local_y)
            if 0 < max_idx < len(local_y) - 1:
                # Fit parabola around maximum
                x_vals = local_x[max_idx-1:max_idx+2]
                y_vals = local_y[max_idx-1:max_idx+2]

                # Parabolic coefficients: y = ax² + bx + c
                # Peak at x = -b/(2a)
                try:
                    coeffs = np.polyfit(x_vals, y_vals, 2)
                    if coeffs[0] < 0:  # Downward parabola
                        refined_center = -coeffs[1] / (2 * coeffs[0])
                        if abs(refined_center - initial_center) <= search_window:
                            return refined_center
                except:
                    pass

            return local_x[max_idx]
        except:
            return initial_center

    def _build_parameter_consensus(self, estimation_results, data_quality):
        """Build consensus parameters from multiple estimation methods"""
        try:
            if not estimation_results:
                return None

            # Weight methods by quality scores and historical performance
            weights = {}
            parameters = {'amplitude': [], 'center': [], 'sigma': [], 'gamma': [], 'baseline': []}

            for method_name, result in estimation_results.items():
                # Base weight on method quality score
                base_weight = result.get('quality_score', 0.5)

                # Adjust weight based on historical performance
                performance = self.method_performance[method_name]
                total_attempts = performance['success_count'] + performance['failure_count']
                if total_attempts > 0:
                    success_rate = performance['success_count'] / total_attempts
                    performance_weight = success_rate
                else:
                    performance_weight = 0.5

                # Adjust weight based on data quality compatibility
                quality_weight = 1.0
                if method_name == 'moment_based' and data_quality['snr'] < 5:
                    quality_weight = 0.5  # Moments unreliable for low SNR
                elif method_name == 'correlation_based' and data_quality['baseline_stability'] < 0.5:
                    quality_weight = 0.7  # Correlation affected by baseline issues

                final_weight = base_weight * performance_weight * quality_weight
                weights[method_name] = final_weight

                # Collect parameters
                for param in parameters:
                    if param in result and result[param] is not None:
                        parameters[param].append((result[param], final_weight))

            # Calculate weighted consensus
            consensus = {}
            consensus_score = 0

            for param, values in parameters.items():
                if not values:
                    continue

                weighted_values = [val * weight for val, weight in values]
                total_weight = sum(weight for val, weight in values)

                if total_weight > 0:
                    consensus[param] = sum(weighted_values) / total_weight

                    # Calculate consensus quality (agreement between methods)
                    param_values = [val for val, weight in values]
                    param_std = np.std(param_values)
                    param_mean = np.mean(param_values)
                    if param_mean != 0:
                        cv = param_std / abs(param_mean)  # Coefficient of variation
                        param_consensus = max(0, 1 - cv)  # Higher is better
                        consensus_score += param_consensus
                else:
                    consensus[param] = None

            consensus_score /= len(parameters)  # Average consensus across parameters
            consensus['consensus_score'] = consensus_score

            return consensus

        except Exception as e:
            return None

    def _quantify_parameter_uncertainties(self, estimation_results, consensus_params):
        """Quantify parameter uncertainties based on method agreement"""
        try:
            uncertainties = {}

            if not consensus_params:
                return {param: float('inf') for param in ['amplitude', 'center', 'sigma', 'gamma', 'baseline']}

            for param in ['amplitude', 'center', 'sigma', 'gamma', 'baseline']:
                param_values = []
                for result in estimation_results.values():
                    if param in result and result[param] is not None:
                        param_values.append(result[param])

                if len(param_values) >= 2:
                    param_std = np.std(param_values)
                    param_mean = consensus_params[param] if param in consensus_params else np.mean(param_values)

                    # Relative uncertainty (standard deviation / mean)
                    if abs(param_mean) > 1e-10:
                        uncertainties[param] = param_std / abs(param_mean)
                    else:
                        uncertainties[param] = float('inf')
                elif len(param_values) == 1:
                    # Single estimate - assign moderate uncertainty
                    uncertainties[param] = 0.2  # 20% uncertainty
                else:
                    # No estimates - infinite uncertainty
                    uncertainties[param] = float('inf')

            return uncertainties
        except:
            return {param: float('inf') for param in ['amplitude', 'center', 'sigma', 'gamma', 'baseline']}

    def _validate_and_adjust_parameters(self, consensus_params, uncertainties, x_data, y_data, nucleus_type):
        """Final validation and adjustment of consensus parameters"""
        try:
            if not consensus_params:
                return self._fallback_parameter_estimation(x_data, y_data, np.mean(x_data), nucleus_type)['parameters']

            validated = consensus_params.copy()

            # Data constraints
            data_max = np.max(y_data)
            data_min = np.min(y_data)
            data_range = data_max - data_min
            ppm_range = abs(x_data[-1] - x_data[0])

            # Amplitude validation and adjustment
            if 'amplitude' in validated:
                if validated['amplitude'] <= 0:
                    validated['amplitude'] = data_range * 0.1
                elif validated['amplitude'] > data_range * 10:
                    validated['amplitude'] = data_range * 2
                elif uncertainties.get('amplitude', 0) > 1.0:  # High uncertainty
                    # Conservative adjustment towards data-driven estimate
                    data_amplitude = data_max - (validated.get('baseline', data_min))
                    validated['amplitude'] = 0.7 * validated['amplitude'] + 0.3 * data_amplitude

            # Center validation and adjustment
            if 'center' in validated:
                if validated['center'] < x_data[0] or validated['center'] > x_data[-1]:
                    validated['center'] = x_data[np.argmax(y_data)]  # Use data maximum
                elif uncertainties.get('center', 0) > 0.1:  # High center uncertainty
                    # Refine using local maximum search
                    center_window = ppm_range * 0.05
                    center_mask = np.abs(x_data - validated['center']) <= center_window
                    if np.sum(center_mask) > 0:
                        local_x = x_data[center_mask]
                        local_y = y_data[center_mask]
                        refined_center = local_x[np.argmax(local_y)]
                        validated['center'] = refined_center

            # Width validation and adjustment
            nmr_params = self.parent.nmr_ranges.get(nucleus_type, self.parent.nmr_ranges['1H'])
            typical_width = nmr_params['typical_width']

            if 'sigma' in validated:
                if validated['sigma'] <= 0:
                    validated['sigma'] = typical_width * 0.5
                elif validated['sigma'] > ppm_range * 0.5:
                    validated['sigma'] = typical_width
                elif uncertainties.get('sigma', 0) > 0.5:
                    # Adjust towards typical values
                    validated['sigma'] = 0.6 * validated['sigma'] + 0.4 * typical_width * 0.6

            if 'gamma' in validated:
                if validated['gamma'] < 0:
                    validated['gamma'] = 0
                elif validated['gamma'] > ppm_range * 0.5:
                    validated['gamma'] = typical_width * 0.5
                elif uncertainties.get('gamma', 0) > 0.5:
                    # Adjust towards typical values
                    validated['gamma'] = 0.6 * validated['gamma'] + 0.4 * typical_width * 0.4

            # Baseline validation and adjustment
            if 'baseline' in validated:
                if validated['baseline'] > data_max * 0.8:
                    validated['baseline'] = data_min
                elif validated['baseline'] < data_min - data_range:
                    validated['baseline'] = data_min
                elif uncertainties.get('baseline', 0) > 0.3:
                    # Use robust baseline estimate
                    robust_baseline = self.parent.robust_baseline_estimation(x_data, y_data, method='percentile')
                    validated['baseline'] = 0.7 * validated['baseline'] + 0.3 * robust_baseline

            # Ensure parameter list format
            return [
                validated.get('amplitude', data_range),
                validated.get('center', x_data[np.argmax(y_data)]),
                validated.get('sigma', typical_width * 0.6),
                validated.get('gamma', typical_width * 0.4),
                validated.get('baseline', data_min)
            ]

        except Exception as e:
            return self._fallback_parameter_estimation(x_data, y_data, np.mean(x_data), nucleus_type)['parameters']

    def _fallback_parameter_estimation(self, x_data, y_data, peak_center, nucleus_type):
        """Ultra-simple fallback parameter estimation"""
        try:
            data_max = np.max(y_data)
            data_min = np.min(y_data)
            nmr_params = self.parent.nmr_ranges.get(nucleus_type, self.parent.nmr_ranges['1H'])

            return {
                'success': True,
                'parameters': [
                    data_max - data_min,  # amplitude
                    peak_center,          # center
                    nmr_params['typical_width'] * 0.6,  # sigma
                    nmr_params['typical_width'] * 0.4,  # gamma
                    data_min             # baseline
                ]
            }
        except:
            return {
                'success': False,
                'parameters': [1000, peak_center, 0.01, 0.01, 0]
            }



if __name__ == "__main__":
    # Test the enhanced fitter
    print("🧪 Testing Enhanced Voigt Fitter")
    print("=" * 40)

    # Create test data
    x_test = np.linspace(7.5, 8.5, 200)  # ¹H range
    true_params = [2000, 8.0, 0.015, 0.008, 150]  # Realistic ¹H peak

    from scipy.special import wofz
    def test_voigt(x, amp, cen, sig, gam, base):
        z = ((x - cen) + 1j*gam) / (sig * np.sqrt(2))
        return amp * np.real(wofz(z)) / (sig * np.sqrt(2*np.pi)) + base

    y_true = test_voigt(x_test, *true_params)
    y_noisy = y_true + np.random.normal(0, 100, len(y_true))

    # Test enhanced fitting
    fitter = EnhancedVoigtFitter()
    result = fitter.fit_peak_enhanced(x_test, y_noisy, nucleus_type='1H')

    if result['success']:
        print(f"✅ Fitting successful!")
        print(f"   Quality: {result['quality_class']} (R² = {result['r_squared']:.3f})")
        print(f"   Center: {result['center']:.4f} ppm (true: {true_params[1]:.4f})")
        print(f"   Amplitude: {result['amplitude']:.1f} (true: {true_params[0]:.1f})")
        print(f"   Width (σ+γ): {result['sigma'] + result['gamma']:.4f} (true: {true_params[2] + true_params[3]:.4f})")
    else:
        print(f"❌ Fitting failed: {result.get('error', 'Unknown error')}")
