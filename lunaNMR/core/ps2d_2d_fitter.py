# ABOUTME: PS2D 2D simultaneous multi-peak Voigt fitter using 5-stage Levenberg-Marquardt optimization.
# ABOUTME: Handles overlapping peak clusters with elliptical data selection and absolute position/linewidth constraints.
"""
PS2D 2D Multi-Peak Fitter - Simultaneous 2D Voigt Fitting
==========================================================

for resolving closely-spaced overlapping peaks in 2D NMR spectra.

1. Simultaneous 2D fitting of multiple overlapping peaks
2. n-stage Levenberg-Marquardt optimization strategy
3. Union of elliptical windows for overlap group data selection
4. Relative numerical derivatives with analytical intensity derivatives

Critical Use Case:
------------------
Resolves peaks that overlap heavily in one dimension but are separated in another.
Example: Peaks 144-145 at (8.747, 112.9) and (8.728, 113.3)
- 1H separation: 0.019 ppm (heavily overlapping)
- 15N separation: 0.4 ppm (well separated)

Traditional 1D cross-section fitting fails because:
- X cross-sections see single blended peak → wrong position
- Y cross-sections contaminated by 1H overlap → wrong intensity

2D simultaneous fitting resolves both peaks correctly by using full 2D lineshape.

Date: 2025-10-06
Version: 1.0 
"""

import numpy as np
from typing import Dict, List
import sys

# Import building blocks
from .ps2d_style_fitter import (
    multi_voigt_profile_2d,
    compute_multi_voigt_jacobian_2d,
    Ps2dStyleLevenbergMarquardt,
    SQRT_2,
    SQRT_8LN2
)
from .ps2d_config import get_ps2d_config
from scipy.special import wofz


# Import training data collector (optional)
try:
    from ..ml import PS2DTrainingDataCollector
    ML_COLLECTOR_AVAILABLE = True
except ImportError:
    PS2DTrainingDataCollector = None
    ML_COLLECTOR_AVAILABLE = False


def calculate_peak_height(lw_lor_f1: float, lw_gau_f1: float,
                          lw_lor_f2: float, lw_gau_f2: float,
                          intensity: float) -> float:
    """
    Calculate peak height (maximum value at peak center) from Voigt parameters.

    For 2D Voigt profile at center (Δf1=0, Δf2=0):
    height = intensity × V_f1(0) × V_f2(0)
    where V(0) = Re[w(i·γ/σ√2)] / (σ√(2π))

    Parameters:
    -----------
    lw_lor_f1, lw_lor_f2 : float
        Lorentzian FWHM (ppm)
    lw_gau_f1, lw_gau_f2 : float
        Gaussian FWHM (ppm)
    intensity : float
        Fitted intensity parameter (= volume for normalized Voigt)

    Returns:
    --------
    float : Peak height at center position
    """
    # Convert FWHM to sigma/gamma
    sigma_f1 = lw_gau_f1 / SQRT_8LN2
    sigma_f2 = lw_gau_f2 / SQRT_8LN2
    gamma_f1 = lw_lor_f1 / 2.0
    gamma_f2 = lw_lor_f2 / 2.0

    # Prevent division by zero
    sigma_f1 = max(sigma_f1, 1e-10)
    sigma_f2 = max(sigma_f2, 1e-10)

    # At peak center (Δf=0), z = i·γ/(σ√2)
    z1 = 1j * gamma_f1 / (sigma_f1 * SQRT_2)
    z2 = 1j * gamma_f2 / (sigma_f2 * SQRT_2)

    # Compute Faddeeva at center
    fade_f1_center = np.real(wofz(z1))
    fade_f2_center = np.real(wofz(z2))

    # Height = intensity × normalized Voigt at center
    height = intensity * fade_f1_center * fade_f2_center / (sigma_f1 * sigma_f2 * 2.0 * np.pi)

    return height


class Ps2dMultiPeakFitter2D:
    """
    2D multi-peak fitter using EXACT  strategy

    Implements simultaneous Levenberg-Marquardt fitting of multiple 2D Voigt
    """

    def __init__(self, verbose: bool = False, training_collector=None, lg_penalty_weight: float = 0.0, constrain_intensity_ratios: bool = True, intensity_ratio_penalty_weight: float = 10.0):
        """
        Initialize 2D multi-peak fitter

        Parameters:
        -----------
        verbose : bool
            Print detailed fitting progress
        training_collector : PS2DTrainingDataCollector, optional
            Training data collector for ML model development
        lg_penalty_weight : float
            Penalty weight for enforcing similar L/G ratios across cluster.
            Default 10.0 provides gentle guidance (peaks encouraged to similar shapes).
            Set to 0.0 to disable constraint (unconstrained optimization).
            Only applied when heavy overlap detected.
        constrain_intensity_ratios : bool
            If True, add soft penalty to guide intensity ratios toward measured ratios.
            Default True (true soft constraint via penalty in objective function).
            Set to False to allow intensities to float freely without guidance.
            Only applied when heavy overlap detected.
        intensity_ratio_penalty_weight : float
            Penalty weight for intensity ratio guidance (λ parameter).
            Default 10.0 provides gentle guidance (balanced with L/G constraint).
            Recommended range: 10-100.
            Higher values = stronger guidance toward measured ratios.
            Only applied when heavy overlap detected.
        """
        self.verbose = verbose
        self.lg_penalty_weight = lg_penalty_weight
        self.constrain_intensity_ratios = constrain_intensity_ratios
        self.intensity_ratio_penalty_weight = intensity_ratio_penalty_weight

        # Get max_iterations from ps2d_config
        config = get_ps2d_config()
        max_iterations = config.max_iterations

        self.optimizer = Ps2dStyleLevenbergMarquardt(
            verbose=verbose,
            max_iter=max_iterations
        )
        self.training_collector = training_collector

    def _calculate_lg_ratio(self, lw_lor: float, lw_gau: float) -> float:
        """
        Calculate Lorentzian fraction (L/G ratio).

        Parameters:
        -----------
        lw_lor : float
            Lorentzian FWHM (ppm)
        lw_gau : float
            Gaussian FWHM (ppm)

        Returns:
        --------
        float : L/(L+G) ratio, range [0, 1]
                0 = pure Gaussian, 1 = pure Lorentzian
        """
        total_lw = lw_lor + lw_gau
        if total_lw < 1e-10:  # Avoid division by zero
            return 0.5  # Default to 50/50 if both components are zero
        return lw_lor / total_lw

    def _calculate_lg_penalty(self, params: np.ndarray, n_peaks: int) -> float:
        """
        Calculate penalty for L/G ratio variance across peaks in cluster.

        Enforces that all peaks have similar Lorentzian/Gaussian character
        (prevents artifacts like one peak being 100% Lorentzian, another 50%).

        Parameters:
        -----------
        params : np.ndarray
            Flattened parameter array for all peaks
        n_peaks : int
            Number of peaks in cluster

        Returns:
        --------
        float : Variance of L/G ratios (0 = all peaks identical, >0 = different)
        """
        NPAR_VOIGT = 8

        # VECTORIZED: Extract all L/G ratios at once using NumPy indexing
        # Parameter order: pos_f1(0), lw_lor_f1(1), pos_f2(2), lw_lor_f2(3),
        #                  lw_gau_f1(4), lw_gau_f2(5), intensity(6), baseline(7)
        offsets = np.arange(n_peaks) * NPAR_VOIGT

        # F1 dimension (15N/13C): lw_lor_f1 at index 1, lw_gau_f1 at index 4
        lw_lor_f1 = params[offsets + 1]
        lw_gau_f1 = params[offsets + 4]
        total_f1 = lw_lor_f1 + lw_gau_f1
        # Vectorized L/G ratio with safe division
        lg_ratios_f1 = np.where(total_f1 > 1e-10, lw_lor_f1 / total_f1, 0.5)

        # F2 dimension (1H): lw_lor_f2 at index 3, lw_gau_f2 at index 5
        lw_lor_f2 = params[offsets + 3]
        lw_gau_f2 = params[offsets + 5]
        total_f2 = lw_lor_f2 + lw_gau_f2
        lg_ratios_f2 = np.where(total_f2 > 1e-10, lw_lor_f2 / total_f2, 0.5)

        # Calculate variance for each dimension (np.var uses N denominator by default)
        variance_f1 = np.var(lg_ratios_f1)
        variance_f2 = np.var(lg_ratios_f2)

        # Total penalty: sum of variances in both dimensions
        total_penalty = variance_f1 + variance_f2

        return total_penalty

    def _fit_stage1_with_lg_penalty(self, model_function, jacobian_function, x, y, p0, bounds, fixed_params, n_peaks):
        """
        Custom Stage 1 fitting with L/G ratio penalty.

        Uses iterative approach: optimize with standard LM, then nudge linewidths
        toward cluster mean L/G ratio, repeat until convergence.

        This avoids the residual vector dimension mismatch issue.
        """
        NPAR_VOIGT = 8
        params = p0.copy()

        # Run standard optimization first
        params, cov, info = self.optimizer.fit(
            func=model_function,
            jacobian=jacobian_function,
            x=x,
            y=y,
            p0=params,
            bounds=bounds,
            fixed_params=fixed_params
        )

        # Now apply L/G ratio constraint via iterative refinement
        max_penalty_iterations = 10
        for penalty_iter in range(max_penalty_iterations):
            # Calculate current L/G ratios
            lg_ratios_f1 = []
            lg_ratios_f2 = []
            for i in range(n_peaks):
                offset = i * NPAR_VOIGT
                lg_f1 = self._calculate_lg_ratio(params[offset + 1], params[offset + 2])
                lg_f2 = self._calculate_lg_ratio(params[offset + 4], params[offset + 5])
                lg_ratios_f1.append(lg_f1)
                lg_ratios_f2.append(lg_f2)

            # Calculate cluster mean L/G ratios
            mean_lg_f1 = np.mean(lg_ratios_f1)
            mean_lg_f2 = np.mean(lg_ratios_f2)

            # Check convergence: if L/G ratios are already similar, stop
            variance = self._calculate_lg_penalty(params, n_peaks)
            if variance < 1e-4:  # Converged (L/G ratios within ~1%)
                break

            # Directly enforce cluster mean L/G ratio (no nudging - direct assignment)
            # This is a hard constraint: total linewidth conserved, L/G ratio forced to cluster mean
            for i in range(n_peaks):
                offset = i * NPAR_VOIGT

                # F1 dimension: force to cluster mean
                total_lw_f1 = params[offset + 1] + params[offset + 2]
                params[offset + 1] = total_lw_f1 * mean_lg_f1  # Lorentzian component
                params[offset + 2] = total_lw_f1 * (1.0 - mean_lg_f1)  # Gaussian component

                # F2 dimension: force to cluster mean
                total_lw_f2 = params[offset + 4] + params[offset + 5]
                params[offset + 4] = total_lw_f2 * mean_lg_f2  # Lorentzian component
                params[offset + 5] = total_lw_f2 * (1.0 - mean_lg_f2)  # Gaussian component

            # Clip to bounds
            params = np.clip(params, bounds[0], bounds[1])

            # Restore fixed parameters (positions and spare)
            for param_idx, fixed_value in fixed_params.items():
                params[param_idx] = fixed_value

        return params, cov, info

    def _fit_stage1_with_intensity_penalty(self, model_function, jacobian_function, x, y, p0, bounds,
                                          fixed_params, n_peaks, target_ratios, lambda_ratio):
        """
        Custom Stage 1 fitting with intensity ratio penalty PROPERLY integrated into objective.

        This is a TRUE soft constraint - the penalty is part of what the optimizer minimizes
        at every iteration, not post-processing. The optimizer naturally finds parameters
        that balance good data fit (low χ²) with staying near target ratios (low penalty).

        Mathematical formulation:
        ------------------------
        Standard LM minimizes: Σ(y - model)²

        With penalty, we augment the objective:
          y_augmented = [y₁, y₂, ..., yₙ, 0, 0, ..., 0]  ← n_peaks zeros for penalty targets
          model_augmented = [m₁, m₂, ..., mₙ, √λ·p₁, √λ·p₂, ..., √λ·pₙ]

          where pᵢ = (Iᵢ/Iₘₐₓ - target_ratioᵢ) = deviation from target ratio

        LM minimizes: Σ(y - model)² + Σ(√λ·pᵢ)²
                    = χ²_data + λ × Penalty_ratio ✅

        Parameters:
        -----------
        target_ratios : List[float]
            Target intensity ratios from raw data (normalized to max=1.0)
        lambda_ratio : float
            Penalty weight (e.g., 10-100). Higher = stronger constraint.
            Recommended: 10-50 (weak guidance, just inform the fit)
        """
        NPAR_VOIGT = 8

        # Create augmented model function that includes ratio penalty residuals
        def augmented_model_function(x_dummy, *params_tuple):
            """Model function with penalty residuals appended"""
            # Standard model evaluation (optimizer passes params as *args)
            params_flat = np.array(params_tuple)
            model_values = model_function(x_dummy, *params_tuple)  # Shape: (N_data,)

            # Calculate ratio penalty residuals (one per peak)
            if n_peaks < 2:  # Single peak, no ratio penalty
                # Still need to augment for consistency, but penalty is zero
                penalty_residuals = np.zeros(n_peaks)
                return np.concatenate([model_values, penalty_residuals])

            # Extract current intensities
            intensities = [params_flat[i*NPAR_VOIGT + 6] for i in range(n_peaks)]
            max_intensity = max(intensities)

            if max_intensity < 1e-10:  # Avoid division by zero
                penalty_residuals = np.zeros(n_peaks)
                return np.concatenate([model_values, penalty_residuals])

            # Calculate current ratios (normalized to max=1.0)
            current_ratios = [I / max_intensity for I in intensities]

            # Penalty residuals: √λ × (current_ratio - target_ratio)
            penalty_residuals = []
            for i in range(n_peaks):
                deviation = current_ratios[i] - target_ratios[i]
                penalty_residuals.append(np.sqrt(lambda_ratio) * deviation)

            # Concatenate: data residuals + penalty residuals
            return np.concatenate([model_values, penalty_residuals])

        # Create augmented Jacobian that includes penalty derivatives
        def augmented_jacobian_function(x_dummy, params_array):
            """Jacobian with penalty derivatives appended"""
            # Standard Jacobian evaluation
            J_data = jacobian_function(x_dummy, params_array)  # Shape: (N_data, N_params)

            # Create Jacobian for penalty residuals
            # For penalty p_i = √λ × (I_i/I_max - target_i), we need ∂p_i/∂I_j
            J_penalty = np.zeros((n_peaks, len(params_array)))

            if n_peaks >= 2:
                # Extract intensities
                intensities = np.array([params_array[i*NPAR_VOIGT + 6] for i in range(n_peaks)])
                max_intensity = np.max(intensities)

                if max_intensity > 1e-10:
                    # Index of peak with max intensity
                    max_idx = np.argmax(intensities)

                    # Derivatives: ∂(I_i/I_max)/∂I_j
                    for i in range(n_peaks):
                        for j in range(n_peaks):
                            param_idx = j * NPAR_VOIGT + 6  # Intensity parameter
                            if i == j == max_idx:
                                # ∂(I_max/I_max)/∂I_max = 0
                                deriv = 0.0
                            elif i == j:
                                # ∂(I_i/I_max)/∂I_i = 1/I_max
                                deriv = 1.0 / max_intensity
                            elif j == max_idx:
                                # ∂(I_i/I_max)/∂I_max = -I_i/I_max²
                                deriv = -intensities[i] / (max_intensity**2)
                            else:
                                # No dependence
                                deriv = 0.0

                            J_penalty[i, param_idx] = np.sqrt(lambda_ratio) * deriv

            # Concatenate: [J_data; J_penalty]
            return np.vstack([J_data, J_penalty])

        # Create augmented target data (original y + zeros for penalty targets)
        y_augmented = np.concatenate([y, np.zeros(n_peaks)])

        # Fit with augmented objective
        params, cov, info = self.optimizer.fit(
            func=augmented_model_function,
            jacobian=augmented_jacobian_function,
            x=x,
            y=y_augmented,
            p0=p0,
            bounds=bounds,
            fixed_params=fixed_params
        )

        return params, cov, info

    def _fit_stage4_with_intensity_penalty(self, model_function, jacobian_function, x, y, p0, bounds,
                                          fixed_params, n_peaks, target_ratios, lambda_ratio, optimizer):
        """
        Custom Stage 4 fitting with intensity ratio penalty PROPERLY integrated into objective.

        Same mathematical formulation as Stage 1, but used in global refinement stage.

        Mathematical formulation:
        ------------------------
        Standard LM minimizes: Σ(y - model)²

        With penalty, we augment the objective:
          y_augmented = [y₁, y₂, ..., yₙ, 0, 0, ..., 0]  ← n_peaks zeros for penalty targets
          model_augmented = [m₁, m₂, ..., mₙ, √λ·p₁, √λ·p₂, ..., √λ·pₙ]

          where pᵢ = (Iᵢ/Iₘₐₓ - target_ratioᵢ) = deviation from target ratio

        LM minimizes: Σ(y - model)² + Σ(√λ·pᵢ)²
                    = χ²_data + λ × Penalty_ratio ✅

        Parameters:
        -----------
        optimizer : Ps2dStyleLevenbergMarquardt
            Optimizer instance with max_iterations from ps2d_config

        Other parameters same as _fit_stage1_with_intensity_penalty()
        """
        NPAR_VOIGT = 8

        # Create augmented model function that includes ratio penalty residuals
        def augmented_model_function(x_dummy, *params_tuple):
            """Model function with penalty residuals appended"""
            # Standard model evaluation (optimizer passes params as *args)
            params_flat = np.array(params_tuple)
            model_values = model_function(x_dummy, *params_tuple)  # Shape: (N_data,)

            # Calculate ratio penalty residuals (one per peak)
            if n_peaks < 2:  # Single peak, no ratio penalty
                penalty_residuals = np.zeros(n_peaks)
                return np.concatenate([model_values, penalty_residuals])

            # Extract current intensities
            intensities = [params_flat[i*NPAR_VOIGT + 6] for i in range(n_peaks)]
            max_intensity = max(intensities)

            if max_intensity < 1e-10:  # Avoid division by zero
                penalty_residuals = np.zeros(n_peaks)
                return np.concatenate([model_values, penalty_residuals])

            # Calculate current ratios (normalized to max=1.0)
            current_ratios = [I / max_intensity for I in intensities]

            # Penalty residuals: √λ × (current_ratio - target_ratio)
            penalty_residuals = []
            for i in range(n_peaks):
                deviation = current_ratios[i] - target_ratios[i]
                penalty_residuals.append(np.sqrt(lambda_ratio) * deviation)

            # Concatenate: data residuals + penalty residuals
            return np.concatenate([model_values, penalty_residuals])

        # Create augmented Jacobian that includes penalty derivatives
        def augmented_jacobian_function(x_dummy, params_array):
            """Jacobian with penalty derivatives appended"""
            # Standard Jacobian evaluation
            J_data = jacobian_function(x_dummy, params_array)  # Shape: (N_data, N_params)

            # Create Jacobian for penalty residuals
            # For penalty p_i = √λ × (I_i/I_max - target_i), we need ∂p_i/∂I_j
            J_penalty = np.zeros((n_peaks, len(params_array)))

            if n_peaks >= 2:
                # Extract intensities
                intensities = np.array([params_array[i*NPAR_VOIGT + 6] for i in range(n_peaks)])
                max_intensity = np.max(intensities)

                if max_intensity > 1e-10:
                    # Index of peak with max intensity
                    max_idx = np.argmax(intensities)

                    # Derivatives: ∂(I_i/I_max)/∂I_j
                    for i in range(n_peaks):
                        for j in range(n_peaks):
                            param_idx = j * NPAR_VOIGT + 6  # Intensity parameter
                            if i == j == max_idx:
                                # ∂(I_max/I_max)/∂I_max = 0
                                deriv = 0.0
                            elif i == j:
                                # ∂(I_i/I_max)/∂I_i = 1/I_max
                                deriv = 1.0 / max_intensity
                            elif j == max_idx:
                                # ∂(I_i/I_max)/∂I_max = -I_i/I_max²
                                deriv = -intensities[i] / (max_intensity**2)
                            else:
                                # No dependence
                                deriv = 0.0

                            J_penalty[i, param_idx] = np.sqrt(lambda_ratio) * deriv

            # Concatenate: [J_data; J_penalty]
            return np.vstack([J_data, J_penalty])

        # Create augmented target data (original y + zeros for penalty targets)
        y_augmented = np.concatenate([y, np.zeros(n_peaks)])

        # Fit with augmented objective using provided optimizer
        params, cov, info = optimizer.fit(
            func=augmented_model_function,
            jacobian=augmented_jacobian_function,
            x=x,
            y=y_augmented,
            p0=p0,
            bounds=bounds,
            fixed_params=fixed_params
        )

        return params, cov, info

    def fit_multi_peak_2d(self,
                          f1_grid: np.ndarray,
                          f2_grid: np.ndarray,
                          intensity: np.ndarray,
                          initial_peaks: List[Dict],
                          fix_positions: bool = False,
                          fix_linewidths: bool = False,
                          data_mask: np.ndarray = None,
                          spectrum_statistics: Dict = None,
                          reference_positions: Dict = None,
                          peak_assignments: List[str] = None) -> Dict:
        """
        Fit multiple overlapping 2D Voigt peaks simultaneously

        This is the CORE function that implements  2D multi-peak fitting.
        Uses 5-stage Levenberg-Marquardt strategy for stable convergence.

        Parameters:
        -----------
        f1_grid, f2_grid : np.ndarray
            2D frequency grids (from meshgrid), same shape
        intensity : np.ndarray
            2D intensity data, same shape as grids
        initial_peaks : list of dict
            Initial parameter guesses for each peak:
            [{'pos_f1', 'lw_lor_f1', 'lw_gau_f1',
              'pos_f2', 'lw_lor_f2', 'lw_gau_f2', 'intensity'}, ...]
        fix_positions : bool
            If True, peak positions are fixed during fitting
        fix_linewidths : bool
            If True, linewidths are fixed during fitting
        data_mask : np.ndarray, optional
            Boolean mask indicating which data points to use for fitting.
        spectrum_statistics : dict, optional
            Learned linewidth statistics from PASS1 isolated peak fitting.
            Used to compute adaptive bounds: min(config_max, median + 5×MAD).
            Keys: 'lw_f1_median', 'lw_f1_mad', 'lw_f2_median', 'lw_f2_mad'.
            PS2D uses union of elliptical windows
            If None, all data points are used (NOT recommended for overlapping peaks).
        reference_positions : dict, optional
            Mapping of assignment -> (x_ppm, y_ppm) for cascade mode absolute drift limiting.
            If provided, position bounds are clipped to stay within max_drift of original reference.
        peak_assignments : list of str, optional
            Assignment names for each peak, needed for reference_positions lookup.

        Returns:
        --------
        dict : Fitting results with keys:
            'success': bool - convergence flag
            'n_peaks': int - number of peaks fitted
            'peaks': list of dict - fitted parameters for each peak
            'r_squared': float - goodness of fit
            'chi2': float - final chi-squared value
            'iterations': int - total iterations across all stages
            'fitted_2d': np.ndarray - fitted 2D surface

        Notes:
        ------
        5-stage fitting strategy
        - Stage 0: Fix positions/widths, fit intensities only
        - Stage 1: Fix positions, float widths + intensities
        - Stage 2: Float positions (if allowed)
        - Stage 3: For future implementations
        - Stage 4: Final global refinement (all parameters float)
        """

        n_peaks = len(initial_peaks)

        #if self.verbose:
        #    print("=" * 70)
        #    print(f"PS2D 2D Multi-Peak Fitting: {n_peaks} peaks")
        #    print("=" * 70)
        #    sys.stdout.flush()

        # Convert initial peaks to flat parameter array
        # Layout: [pos_f1, lw_lor_f1, lw_gau_f1, pos_f2, lw_lor_f2, lw_gau_f2, intensity, spare] × n_peaks
        params = []
        for peak in initial_peaks:
            params.extend([
                peak['pos_f1'], peak['lw_lor_f1'], peak['lw_gau_f1'],
                peak['pos_f2'], peak['lw_lor_f2'], peak['lw_gau_f2'],
                peak['intensity'], 0.0  # spare parameter unused
            ])
        params = np.array(params)

        # Check if cluster has peaks that are too close (apply L/G and intensity constraints)
        # DECOUPLED from heavy_overlap (which only affects linewidth estimation)
        # tooclose flag indicates spectral ambiguity requiring constraints
        has_tooclose = any(peak.get('tooclose', False) for peak in initial_peaks)

        # Store initial intensity ratios from raw data (only used if peaks are too close)
        # CRITICAL: Use original_height (measured peak height) NOT intensity (volume estimate)
        # These ratios represent the measured peak heights and should be preserved
        if has_tooclose:
            initial_intensities = [peak.get('original_height', peak['intensity']) for peak in initial_peaks]
            max_initial_intensity = max(initial_intensities)
            if max_initial_intensity > 0:
                initial_intensity_ratios = [intensity / max_initial_intensity for intensity in initial_intensities]
            else:
                initial_intensity_ratios = [1.0 / n_peaks] * n_peaks  # Equal if all zero
        else:
            # Peaks not too close - disable intensity ratio constraint for this cluster
            initial_intensity_ratios = []

        # Flatten grids and intensity for fitting
        f1_flat = f1_grid.ravel()
        f2_flat = f2_grid.ravel()
        y_flat = intensity.ravel()

        # Apply data mask if provided (PS2D elliptical window filtering)
        if data_mask is not None:
            mask_flat = data_mask.ravel()
            y_flat_masked = y_flat[mask_flat]
            n_masked = np.sum(mask_flat)
            if self.verbose:
                n_total = len(y_flat)
                sys.stdout.flush()
        else:
            mask_flat = np.ones(len(y_flat), dtype=bool)
            y_flat_masked = y_flat

        # Create wrapper functions for optimizer
        def model_function(_f1_f2_dummy, *p):
            """Model function for optimizer (ignores x, uses stored grids)"""
            y_pred_full = multi_voigt_profile_2d(f1_grid, f2_grid, np.array(p), n_peaks).ravel()
            return y_pred_full[mask_flat]  # Return only masked points

        def jacobian_function(_f1_f2_dummy, p):
            """Jacobian function for optimizer"""
            jac_full = compute_multi_voigt_jacobian_2d(f1_grid, f2_grid, np.array(p), n_peaks)
            return jac_full[mask_flat, :]  # Return only masked rows

        # L/G penalty will be applied via custom optimizer wrapper (not by extending residuals)

        # Set parameter bounds
        NPAR_VOIGT = 8

        # Calculate median linewidths for diagnostics only (not used for bounds)
        initial_lw_lor_f1 = [peak['lw_lor_f1'] for peak in initial_peaks]
        initial_lw_gau_f1 = [peak['lw_gau_f1'] for peak in initial_peaks]
        initial_lw_lor_f2 = [peak['lw_lor_f2'] for peak in initial_peaks]
        initial_lw_gau_f2 = [peak['lw_gau_f2'] for peak in initial_peaks]

        median_lw_lor_f1 = np.median(initial_lw_lor_f1)
        median_lw_gau_f1 = np.median(initial_lw_gau_f1)
        median_lw_lor_f2 = np.median(initial_lw_lor_f2)
        median_lw_gau_f2 = np.median(initial_lw_gau_f2)

        # Get nucleus-adaptive parameters from centralized config
        config = get_ps2d_config()

        # ====================================================================
        # ADAPTIVE LINEWIDTH BOUNDS (PASS1 → PASS2 Learning)
        # ====================================================================
        # Two-pass fitting strategy:
        #   PASS1: Fit isolated peaks, collect linewidth statistics (median, MAD)
        #   PASS2: Fit clusters using learned bounds from PASS1
        #
        # Bound formula: min(config_hard_cap, learned_median + 5×MAD)
        #   - 5×MAD ≈ 7.4σ for Gaussian, captures virtually all real variation
        #   - MAD (Median Absolute Deviation) is robust to outliers
        #   - Config max acts as hard safety cap (prevents runaway optimization)
        #
        # Benefits:
        #   - Data-driven: adapts to actual spectrum linewidth distribution
        #   - Tighter bounds prevent optimizer wandering into unrealistic values
        #   - Works automatically - no manual parameter tuning needed
        # ====================================================================

        # Minimum bounds (allow nearly pure Gaussian peaks)
        absolute_min_lw_f1 = 0.0001   # Minimum realistic 15N/13C linewidth (ppm)
        absolute_min_lw_f2 = 0.000005 # Minimum realistic 1H linewidth (ppm)

        # Maximum bounds: use learned statistics if available, else config defaults
        if spectrum_statistics is not None:
            # ================================================================
            # LEARNED BOUNDS WITH MINIMUM MAD FLOOR
            # ================================================================
            # Formula: median + α×MAD, capped by config hard limit
            # - α is sample-size-dependent (7 for n<10, 5 for n≥10)
            # - MAD has minimum floor (5% of typical) to prevent overly tight bounds
            # ================================================================

            # Get sample-size-dependent alpha from statistics (default 5.0)
            alpha = spectrum_statistics.get('alpha', 5.0)

            # Minimum MAD floor: 5% of typical linewidth
            # Prevents zero or near-zero MAD from creating overly tight bounds
            MIN_MAD_F1 = config.typical_linewidth_f1 * 0.05
            MIN_MAD_F2 = config.typical_linewidth_f2 * 0.05

            # Apply MAD floor
            mad_f1 = max(spectrum_statistics.get('lw_f1_mad', MIN_MAD_F1), MIN_MAD_F1)
            mad_f2 = max(spectrum_statistics.get('lw_f2_mad', MIN_MAD_F2), MIN_MAD_F2)

            # Calculate learned bounds: median + α×MAD
            median_f1 = spectrum_statistics.get('lw_f1_median', config.typical_linewidth_f1)
            median_f2 = spectrum_statistics.get('lw_f2_median', config.typical_linewidth_f2)

            learned_max_f1 = median_f1 + alpha * mad_f1
            learned_max_f2 = median_f2 + alpha * mad_f2

            # Cap by config maximum (safety limit)
            absolute_max_lw_f1 = min(config.max_linewidth_f1, learned_max_f1)
            absolute_max_lw_f2 = min(config.max_linewidth_f2, learned_max_f2)

        else:
            # Fallback: use config defaults (nucleus-adaptive, 2-3× typical)
            absolute_max_lw_f1 = config.max_linewidth_f1
            absolute_max_lw_f2 = config.max_linewidth_f2

        # Calculate cluster-wide intensity bounds to prevent artificial constraints
        # Root cause: Per-peak bounds from contaminated 1D cross-sections can trap peaks
        # Solution: Use max initial intensity across cluster to set bounds for all peaks
        max_initial_intensity = max(peak['intensity'] for peak in initial_peaks)


        lower_bounds = []
        upper_bounds = []
        for i in range(n_peaks):
            offset = i * NPAR_VOIGT
            peak = initial_peaks[i]

            # F1 position bounds - constrain to ~1 linewidth of movement
            # FIXED 2025-10-13: lw_gau_f1 IS the Gaussian FWHM (not half-width)
            #fwhm_f1 = peak['lw_gau_f1']  # NEW: lw_gau IS the FWHM (no compensation needed)
            fwhm_f1 = 2.0 * peak['lw_gau_f1']  # OLD: Compensated for FWHM/2 storage bug
            pos_f1_margin = max(0.04 * fwhm_f1, config.pos_margin_f1)  # Nucleus-adaptive minimum #GM max(1.5
            lower_bounds.append(peak['pos_f1'] - pos_f1_margin)
            upper_bounds.append(peak['pos_f1'] + pos_f1_margin)

            # F1 linewidths - absolute bounds: give optimizer full freedom
            # Independent of initial guesses, based on physical limits
            # Applies same bounds to both Lorentzian and Gaussian components
            lower_bounds.extend([absolute_min_lw_f1, absolute_min_lw_f1])  # Lor, Gau
            upper_bounds.extend([absolute_max_lw_f1, absolute_max_lw_f1])  # Lor, Gau

            # F2 position bounds - constrain to ~1 linewidth of movement
            # FIXED 2025-10-13: lw_gau_f2 IS the Gaussian FWHM (not half-width)
            #fwhm_f2 = peak['lw_gau_f2']  # NEW: lw_gau IS the FWHM (no compensation needed)
            fwhm_f2 = 2.0 * peak['lw_gau_f2']  # OLD: Compensated for FWHM/2 storage bug
            pos_f2_margin = max(0.04 * fwhm_f2, config.pos_margin_f2)  # Nucleus-adaptive minimum #GM max(1.5
            lower_bounds.append(peak['pos_f2'] - pos_f2_margin)
            upper_bounds.append(peak['pos_f2'] + pos_f2_margin)

            # F2 linewidths - absolute bounds: give optimizer full freedom
            # Independent of initial guesses, based on physical limits
            # Applies same bounds to both Lorentzian and Gaussian components
            lower_bounds.extend([absolute_min_lw_f2, absolute_min_lw_f2])  # Lor, Gau
            upper_bounds.extend([absolute_max_lw_f2, absolute_max_lw_f2])  # Lor, Gau

            # Intensity bounds - cluster-relative to prevent trapping
            # All peaks can reach 0.1% to 500% of the brightest peak in cluster
            # This prevents per-peak bounds from artificially constraining overlapping peaks
            # where 1D cross-sections give contaminated (too low) initial guesses

            # Intensity bounds - cluster-relative (wide bounds)
            # Intensity ratio guidance is now via penalty in objective, not bounds
            lower_bounds.append(max_initial_intensity * 0.001)  # Min 0.1% of max
            upper_bounds.append(max_initial_intensity * 5.0)    # Max 5× max

            # Spare (always zero)
            lower_bounds.append(0.0)
            upper_bounds.append(0.0)

        bounds = (np.array(lower_bounds), np.array(upper_bounds))

        # ====================================================================
        # CASCADE MODE: Apply absolute drift limits from reference positions
        # ====================================================================
        # In cascade mode, peaks can drift from spectrum to spectrum.
        # Without absolute bounds, drift accumulates: N spectra × pos_margin.
        # This clipping enforces ABSOLUTE bounds relative to ORIGINAL reference.
        #
        # Bounds intersection: final_bound = max(per_spectrum_bound, absolute_bound)
        # This clips positions to stay within max_drift of original reference.
        # ====================================================================
        if reference_positions and peak_assignments:
            drift_clipped_count = 0
            for i in range(n_peaks):
                if i < len(peak_assignments):
                    assignment = peak_assignments[i]
                    if assignment in reference_positions:
                        ref_x, ref_y = reference_positions[assignment]  # (x=F2=1H, y=F1=15N/13C)

                        # Index of F1 position in bounds array (first param per peak)
                        f1_idx = i * NPAR_VOIGT + 0
                        # Index of F2 position in bounds array (fourth param per peak)
                        f2_idx = i * NPAR_VOIGT + 3

                        # Get current per-spectrum bounds
                        current_lower_f1 = bounds[0][f1_idx]
                        current_upper_f1 = bounds[1][f1_idx]
                        current_lower_f2 = bounds[0][f2_idx]
                        current_upper_f2 = bounds[1][f2_idx]

                        # Compute absolute bounds from reference
                        absolute_lower_f1 = ref_y - config.max_drift_f1
                        absolute_upper_f1 = ref_y + config.max_drift_f1
                        absolute_lower_f2 = ref_x - config.max_drift_f2
                        absolute_upper_f2 = ref_x + config.max_drift_f2

                        # Apply intersection (most restrictive bounds)
                        new_lower_f1 = max(current_lower_f1, absolute_lower_f1)
                        new_upper_f1 = min(current_upper_f1, absolute_upper_f1)
                        new_lower_f2 = max(current_lower_f2, absolute_lower_f2)
                        new_upper_f2 = min(current_upper_f2, absolute_upper_f2)

                        # Check if clipping occurred
                        if (new_lower_f1 != current_lower_f1 or new_upper_f1 != current_upper_f1 or
                            new_lower_f2 != current_lower_f2 or new_upper_f2 != current_upper_f2):
                            drift_clipped_count += 1

                        # Update bounds
                        bounds[0][f1_idx] = new_lower_f1
                        bounds[1][f1_idx] = new_upper_f1
                        bounds[0][f2_idx] = new_lower_f2
                        bounds[1][f2_idx] = new_upper_f2


        # Log position constraints for diagnostics (showing last peak's margins as example)
        #if self.verbose:
        #    print(f"   🔒 Position constraints (adaptive, ~1.5× FWHM):")
        #    print(f"      Example: F1 = ±{pos_f1_margin:.3f} ppm, F2 = ±{pos_f2_margin:.4f} ppm")
        #    sys.stdout.flush()

        # Validate and clip initial parameters
        params = np.clip(params, bounds[0], bounds[1])

        # DIAGNOSTIC: Log initial intensity values before any fitting
        #if self.verbose:
        #    print(f"\n   📊 Initial intensity estimates (normalized, before fitting):")
        #    for i in range(n_peaks):
        #        offset = i * NPAR_VOIGT
        #        intensity = params[offset + 6]
        #        lower = bounds[0][offset + 6]
        #        upper = bounds[1][offset + 6]
        #        print(f"      Peak {i+1}: {intensity:.6e} (bounds: [{lower:.6e}, {upper:.6e}])")
        #    sys.stdout.flush()

        total_iterations = 0

        # ====================================================================
        # STAGE 0: DISABLED - Causes intensity collapse with wrong initial linewidths
        # ====================================================================
        # Root cause: Fitting intensity-only with fixed (wrong) linewidths is unstable.
        # When initial linewidths don't match true peak shape (common for 13C peaks),
        # the optimizer drives intensity → 0 to minimize χ² (intensity hits lower bound).
        #
        # Why it fails for 13C but works for 15N:
        # - 15N peaks: FWHM ~0.10 ppm → 1D cross-section guess ~0.10 ppm (close enough)
        # - 13C peaks: FWHM ~0.05 ppm → 1D cross-section guess ~0.10 ppm (2× too broad)
        #
        # Solution: Skip Stage 0, start with Stage 1 (linewidths + intensity together)
        # This allows optimizer to adjust both simultaneously, preventing collapse.
        # ====================================================================

        # Compute initial χ² for convergence tracking (before any optimization)
        y_pred_initial = multi_voigt_profile_2d(f1_grid, f2_grid, params, n_peaks).ravel()
        y_pred_initial_masked = y_pred_initial[mask_flat]
        self._stage0_initial_chi2 = np.sum((y_flat_masked - y_pred_initial_masked)**2)


        # ====================================================================
        # STAGE 1: Fix positions, float linewidths + intensities
        # ====================================================================
        if not fix_linewidths:
            #if self.verbose:
            #    print("\nStage 1: Linewidth + intensity fit (positions fixed)")
            #    if self.lg_penalty_weight > 0.0 and has_tooclose:
            #        print(f"   🔗 L/G ratio constraint enabled (λ={self.lg_penalty_weight:.0f})")
            #        print(f"      Peaks too close: enforcing similar L/G character")
            #    elif self.lg_penalty_weight > 0.0 and not has_tooclose:
            #        print(f"   ℹ️  L/G ratio constraint DISABLED for this cluster")
            #        print(f"      Peaks sufficiently separated: independent L/G ratios allowed")
            #    if self.constrain_intensity_ratios and has_tooclose and len(initial_intensity_ratios) == n_peaks:
            #        print(f"   🔗 Intensity ratio penalty enabled (λ={self.intensity_ratio_penalty_weight:.0f})")
            #        print(f"      Peaks too close: applying ratio constraint")
            #    elif self.constrain_intensity_ratios and not has_tooclose:
            #        print(f"   ℹ️  Intensity ratio penalty DISABLED for this cluster")
            #        print(f"      Peaks sufficiently separated: independent intensities allowed")
            #    sys.stdout.flush()

            # Fix only positions
            fixed_stage1 = {}
            for i in range(n_peaks):
                offset = i * NPAR_VOIGT
                fixed_stage1[offset + 0] = params[offset + 0]  # Fix pos_f1
                fixed_stage1[offset + 3] = params[offset + 3]  # Fix pos_f2
                fixed_stage1[offset + 7] = 0.0  # Fix spare

            # Stage 1: Apply constraints via penalty functions (true soft constraints)
            # Both L/G ratio and intensity ratio penalties integrated into objective
            # CRITICAL: Both constraints ONLY applied when peaks are too close (spectral ambiguity)
            # DECOUPLED from heavy_overlap (which only affects linewidth estimation)

            # Check which penalties to apply
            apply_lg_penalty = (self.lg_penalty_weight > 0.0 and
                               has_tooclose)  # Only for peaks that are too close
            apply_intensity_penalty = (self.constrain_intensity_ratios and
                                      len(initial_intensity_ratios) == n_peaks and
                                      n_peaks >= 2 and
                                      has_tooclose)  # Only for peaks that are too close

            if apply_lg_penalty and apply_intensity_penalty:
                # Apply L/G penalty first
                params, cov, info = self._fit_stage1_with_lg_penalty(
                    model_function=model_function,
                    jacobian_function=jacobian_function,
                    x=f1_flat,
                    y=y_flat_masked,
                    p0=params,
                    bounds=bounds,
                    fixed_params=fixed_stage1,
                    n_peaks=n_peaks
                )
                # Then apply intensity ratio penalty (proper penalty in objective)
                params, cov, info = self._fit_stage1_with_intensity_penalty(
                    model_function=model_function,
                    jacobian_function=jacobian_function,
                    x=f1_flat,
                    y=y_flat_masked,
                    p0=params,
                    bounds=bounds,
                    fixed_params=fixed_stage1,
                    n_peaks=n_peaks,
                    target_ratios=initial_intensity_ratios,
                    lambda_ratio=self.intensity_ratio_penalty_weight
                )
            elif apply_lg_penalty:
                # Apply only L/G penalty
                params, cov, info = self._fit_stage1_with_lg_penalty(
                    model_function=model_function,
                    jacobian_function=jacobian_function,
                    x=f1_flat,
                    y=y_flat_masked,
                    p0=params,
                    bounds=bounds,
                    fixed_params=fixed_stage1,
                    n_peaks=n_peaks
                )
            elif apply_intensity_penalty:
                # Apply only intensity ratio penalty
                params, cov, info = self._fit_stage1_with_intensity_penalty(
                    model_function=model_function,
                    jacobian_function=jacobian_function,
                    x=f1_flat,
                    y=y_flat_masked,
                    p0=params,
                    bounds=bounds,
                    fixed_params=fixed_stage1,
                    n_peaks=n_peaks,
                    target_ratios=initial_intensity_ratios,
                    lambda_ratio=self.intensity_ratio_penalty_weight
                )
            else:
                # Standard optimization without penalties
                params, cov, info = self.optimizer.fit(
                    func=model_function,
                    jacobian=jacobian_function,
                    x=f1_flat,
                    y=y_flat_masked,
                    p0=params,
                    bounds=bounds,
                    fixed_params=fixed_stage1
                )

            total_iterations += info['iterations']

        # ====================================================================
        # STAGE 2: Float positions (if allowed)
        # ====================================================================
        if not fix_positions:

            # Fix spare parameters + respect fix_linewidths flag
            fixed_stage2 = {}
            for i in range(n_peaks):
                offset = i * NPAR_VOIGT
                fixed_stage2[offset + 7] = 0.0  # Fix spare

                # CRITICAL: Respect fix_linewidths flag during position refinement
                # Without this, linewidths would float in Stage 2 even when user wants them fixed
                if fix_linewidths:
                    fixed_stage2[offset + 1] = params[offset + 1]  # Fix lw_lor_f1
                    fixed_stage2[offset + 2] = params[offset + 2]  # Fix lw_gau_f1
                    fixed_stage2[offset + 4] = params[offset + 4]  # Fix lw_lor_f2
                    fixed_stage2[offset + 5] = params[offset + 5]  # Fix lw_gau_f2

            params, cov, info = self.optimizer.fit(
                func=model_function,
                jacobian=jacobian_function,
                x=f1_flat,
                y=y_flat_masked,  # Use masked data (union of elliptical windows)
                p0=params,
                bounds=bounds,
                fixed_params=fixed_stage2
            )

            total_iterations += info['iterations']

        # ====================================================================
        # STAGE 3: Final global refinement
        # ====================================================================
        # DIAGNOSTIC: Save positions before Stage 4 for comparison
        positions_before_stage4 = []
        for i in range(n_peaks):
            offset = i * NPAR_VOIGT
            positions_before_stage4.append({
                'f1': params[offset + 0],
                'f2': params[offset + 3]
            })

        # Fix spare parameters + respect fix_positions/fix_linewidths flags
        # USER EXPECTATION: fix_positions/fix_linewidths are ABSOLUTE constraints
        # Positions/linewidths must remain unchanged regardless of χ² improvement potential
        fixed_stage4 = {}
        for i in range(n_peaks):
            offset = i * NPAR_VOIGT
            fixed_stage4[offset + 7] = 0.0  # Always fix spare parameter

            # CRITICAL: Respect fix_positions flag (user's absolute constraint)
            if fix_positions:
                fixed_stage4[offset + 0] = params[offset + 0]  # Fix pos_f1
                fixed_stage4[offset + 3] = params[offset + 3]  # Fix pos_f2

            # CRITICAL: Respect fix_linewidths flag (user's absolute constraint)
            if fix_linewidths:
                fixed_stage4[offset + 1] = params[offset + 1]  # Fix lw_lor_f1
                fixed_stage4[offset + 2] = params[offset + 2]  # Fix lw_gau_f1
                fixed_stage4[offset + 4] = params[offset + 4]  # Fix lw_lor_f2
                fixed_stage4[offset + 5] = params[offset + 5]  # Fix lw_gau_f2

        # Always use default optimizer (max_iterations from ps2d_config)
        # Iterations are controlled by config, not by fix_positions flag
        stage4_optimizer = self.optimizer

        # Apply intensity ratio penalty if enabled (same as Stage 1)
        # This ensures the penalty continues to guide the fit during global refinement
        apply_intensity_penalty = (self.constrain_intensity_ratios and
                                  len(initial_intensity_ratios) == n_peaks and
                                  n_peaks >= 2)

        if apply_intensity_penalty:

            params, cov, info = self._fit_stage4_with_intensity_penalty(
                model_function=model_function,
                jacobian_function=jacobian_function,
                x=f1_flat,
                y=y_flat_masked,  # Use masked data (union of elliptical windows)
                p0=params,
                bounds=bounds,
                fixed_params=fixed_stage4,
                n_peaks=n_peaks,
                target_ratios=initial_intensity_ratios,
                lambda_ratio=self.intensity_ratio_penalty_weight,
                optimizer=stage4_optimizer
            )
        else:
            params, cov, info = stage4_optimizer.fit(
                func=model_function,
                jacobian=jacobian_function,
                x=f1_flat,
                y=y_flat_masked,  # Use masked data (union of elliptical windows)
                p0=params,
                bounds=bounds,
                fixed_params=fixed_stage4
            )

        total_iterations += info['iterations']

        # ====================================================================
        # L/G RATIO CONSTRAINT: Apply after Stage 4 if enabled
        # ====================================================================
        # CRITICAL: Stage 4 allows linewidths to float (unless fix_linewidths=True)
        # This can undo the L/G constraint applied in Stage 1
        # Solution: Re-apply constraint after Stage 4 optimization completes
        if self.lg_penalty_weight > 0.0 and not fix_linewidths:

            # Calculate current L/G ratios for F1 and F2 dimensions
            lg_ratios_f1 = []
            lg_ratios_f2 = []
            for i in range(n_peaks):
                offset = i * NPAR_VOIGT
                lw_lor_f1 = params[offset + 1]
                lw_gau_f1 = params[offset + 2]
                lw_lor_f2 = params[offset + 4]
                lw_gau_f2 = params[offset + 5]

                lg_f1 = self._calculate_lg_ratio(lw_lor_f1, lw_gau_f1)
                lg_f2 = self._calculate_lg_ratio(lw_lor_f2, lw_gau_f2)
                lg_ratios_f1.append(lg_f1)
                lg_ratios_f2.append(lg_f2)

            # Calculate cluster mean L/G ratios
            mean_lg_f1 = np.mean(lg_ratios_f1)
            mean_lg_f2 = np.mean(lg_ratios_f2)

            # Directly enforce cluster mean L/G ratio
            for i in range(n_peaks):
                offset = i * NPAR_VOIGT

                # F1: force to cluster mean
                total_lw_f1 = params[offset + 1] + params[offset + 2]
                params[offset + 1] = total_lw_f1 * mean_lg_f1  # lw_lor_f1
                params[offset + 2] = total_lw_f1 * (1.0 - mean_lg_f1)  # lw_gau_f1

                # F2: force to cluster mean
                total_lw_f2 = params[offset + 4] + params[offset + 5]
                params[offset + 4] = total_lw_f2 * mean_lg_f2  # lw_lor_f2
                params[offset + 5] = total_lw_f2 * (1.0 - mean_lg_f2)  # lw_gau_f2

            # Clip to bounds and restore fixed parameters
            params = np.clip(params, bounds[0], bounds[1])
            for param_idx, fixed_value in fixed_stage4.items():
                params[param_idx] = fixed_value

        # ====================================================================
        # INTENSITY RATIO CONSTRAINT: Apply after Stage 4 if enabled
        # ====================================================================
        # PHYSICAL PRINCIPLE: Measured peak height ratios provide guidance for volume ratios.
        # Use soft constraint: ratios can vary within tolerance to minimize χ², but
        # cannot deviate too far from measured values (prevents unphysical solutions).
        #
        # APPLIES REGARDLESS OF fix_positions FLAG:
        # Initial intensities come from detected peak positions (already refined by
        # parabolic interpolation), so ratios are valid whether positions float or not.
        # Intensity ratio guidance is now applied in Stage 1 via gentle nudging
        # No post-Stage-4 hard clipping needed (hard clipping destroyed fit quality)


        # ====================================================================
        # Compute final R² and prepare results
        # ====================================================================
        y_fit_flat = multi_voigt_profile_2d(f1_grid, f2_grid, params, n_peaks).ravel()
        y_fit_2d = y_fit_flat.reshape(f1_grid.shape)  # Use grid shape instead of intensity.shape

        # Compute R² using ONLY masked data (PS2D approach)
        # This prevents unfitted peaks outside ellipses from degrading R²
        y_fit_masked = y_fit_flat[mask_flat]
        ss_res = np.sum((y_flat_masked - y_fit_masked)**2)
        ss_tot = np.sum((y_flat_masked - np.mean(y_flat_masked))**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

        # Pragmatic success criterion for 2D multi-peak fitting:
        # Accept result if ANY of:
        # (1) Formal convergence achieved (optimizer converged within tolerance)
        # (2) R² > 0.2 (acceptable fit for overlapping NMR peaks)
        # (3) χ² reduced by > 100× from initial (proves convergence even if R² broken)
        #
        # Rationale: For very narrow peaks on baseline-heavy regions, R² can be
        # negative even when fitting is excellent (proven by massive χ² reduction).
        formal_convergence = info['success']
        pragmatic_r2_success = r_squared > 0.2

        # Check χ² reduction from Stage 0 initial χ²
        # We track Stage 0's first iteration as the baseline
        # Require 100× reduction as evidence of successful convergence
        chi2_reduction_success = False
        chi2_initial = 0.0
        chi2_reduction_factor = 0.0
        if hasattr(self, '_stage0_initial_chi2'):
            chi2_initial = self._stage0_initial_chi2
            chi2_reduction_factor = chi2_initial / info['final_chi2'] if info['final_chi2'] > 0 else 0.0
            chi2_reduction_success = chi2_reduction_factor > 100

        final_success = formal_convergence or pragmatic_r2_success or chi2_reduction_success

        # Extract fitted peaks and calculate derived quantities
        fitted_peaks = []
        for i in range(n_peaks):
            offset = i * NPAR_VOIGT

            # Extract fitted parameters
            lw_lor_f1 = params[offset + 1]
            lw_gau_f1 = params[offset + 2]
            lw_lor_f2 = params[offset + 4]
            lw_gau_f2 = params[offset + 5]
            intensity = params[offset + 6]

            # Calculate derived quantities
            # Volume = intensity (for normalized Voigt, this IS the volume)
            volume = intensity

            # Height = peak maximum at center position
            height = calculate_peak_height(lw_lor_f1, lw_gau_f1, lw_lor_f2, lw_gau_f2, intensity)

            # Amplitude = height (NMR convention)
            amplitude = height

            fitted_peaks.append({
                'pos_f1': params[offset + 0],
                'lw_lor_f1': lw_lor_f1,
                'lw_gau_f1': lw_gau_f1,
                'pos_f2': params[offset + 3],
                'lw_lor_f2': lw_lor_f2,
                'lw_gau_f2': lw_gau_f2,
                'intensity': intensity,
                'volume': volume,
                'height': height,
                'amplitude': amplitude
            })

        # Prepare result dictionary
        result = {
            'success': final_success,
            'formal_convergence': formal_convergence,
            'pragmatic_acceptance': pragmatic_r2_success and not formal_convergence,
            'n_peaks': n_peaks,
            'peaks': fitted_peaks,
            'r_squared': r_squared,
            'chi2': info['final_chi2'],
            'chi2_initial': chi2_initial,
            'chi2_reduction_factor': chi2_reduction_factor,
            'chi2_reduction_success': chi2_reduction_success,
            'iterations': total_iterations,
            'fitted_2d': y_fit_2d,
            'params': params,
            'covariance': cov
        }

        # Collect training data if collector is available and fit was successful
        if self.training_collector is not None and final_success:
            self.training_collector.collect_ps2d_fit(result)

        return result
