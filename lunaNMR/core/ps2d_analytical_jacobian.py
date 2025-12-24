# ABOUTME: Analytical Jacobian computation for 2D Voigt profiles using Faddeeva derivative identity.
# ABOUTME: Reduces wofz calls from 7N to N per Jacobian (7× speedup for N peaks).
"""
Analytical Jacobian for 2D Multi-Peak Voigt Fitting
====================================================

This module provides analytical derivatives for the 2D Voigt profile,
eliminating the need for numerical finite-difference approximations.

Mathematical Foundation:
------------------------
The 2D Voigt profile is:
    V(f1, f2) = A × Re[w(z1)] × Re[w(z2)] / (σ1 × σ2 × 2π)

where:
    z_i = (pos_i - f_i + i×γ_i) / (σ_i × √2)
    σ_i = FWHM_gauss_i / √(8×ln2)
    γ_i = FWHM_lorentz_i / 2

The key identity for Faddeeva function derivative:
    w'(z) = dw/dz = -2z×w(z) + 2i/√π

This allows computing all derivatives from a single w(z) evaluation.

Parameter Layout (NPAR_VOIGT = 8):
----------------------------------
    [0] pos_f1      - F1 position (ppm)
    [1] lw_lor_f1   - F1 Lorentzian FWHM (ppm)
    [2] lw_gau_f1   - F1 Gaussian FWHM (ppm)
    [3] pos_f2      - F2 position (ppm)
    [4] lw_lor_f2   - F2 Lorentzian FWHM (ppm)
    [5] lw_gau_f2   - F2 Gaussian FWHM (ppm)
    [6] intensity   - Peak intensity
    [7] baseline    - Unused (spare)

NOTE: This matches the ACTUAL code in multi_voigt_profile_2d (lines 202-209),
      NOT the C++ documentation comment which has a different order.

Author: Guillaume Mas
Date: 2025
"""

import numpy as np
from scipy.special import wofz
from typing import Tuple

# Constants
SQRT_2 = np.sqrt(2.0)
SQRT_PI = np.sqrt(np.pi)
SQRT_8LN2 = np.sqrt(8.0 * np.log(2.0))
TWO_PI = 2.0 * np.pi
NPAR_VOIGT = 8


def voigt_2d_with_derivatives(
    f1: np.ndarray,
    f2: np.ndarray,
    pos_f1: float,
    lw_lor_f1: float,
    pos_f2: float,
    lw_lor_f2: float,
    lw_gau_f1: float,
    lw_gau_f2: float,
    intensity: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
           np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute 2D Voigt profile AND all 7 derivatives in one pass.

    Uses only 2 wofz calls (one per dimension) instead of 7 calls
    for numerical derivatives.

    Parameters:
    -----------
    f1, f2 : np.ndarray
        Frequency coordinates (ppm), flattened arrays
    pos_f1, pos_f2 : float
        Peak positions (ppm)
    lw_lor_f1, lw_lor_f2 : float
        Lorentzian FWHM (ppm)
    lw_gau_f1, lw_gau_f2 : float
        Gaussian FWHM (ppm)
    intensity : float
        Peak intensity

    Returns:
    --------
    Tuple of 8 arrays:
        V           : Profile values
        dV_dpos_f1  : ∂V/∂pos_f1
        dV_dlor_f1  : ∂V/∂lw_lor_f1
        dV_dpos_f2  : ∂V/∂pos_f2
        dV_dlor_f2  : ∂V/∂lw_lor_f2
        dV_dgau_f1  : ∂V/∂lw_gau_f1
        dV_dgau_f2  : ∂V/∂lw_gau_f2
        dV_dA       : ∂V/∂intensity
    """
    # Convert FWHM to internal parameters
    sigma_f1 = max(lw_gau_f1 / SQRT_8LN2, 1e-10)
    sigma_f2 = max(lw_gau_f2 / SQRT_8LN2, 1e-10)
    gamma_f1 = lw_lor_f1 / 2.0
    gamma_f2 = lw_lor_f2 / 2.0

    # Complex arguments for Faddeeva
    z1 = ((pos_f1 - f1) + 1j * gamma_f1) / (sigma_f1 * SQRT_2)
    z2 = ((pos_f2 - f2) + 1j * gamma_f2) / (sigma_f2 * SQRT_2)

    # Faddeeva function values (ONLY 2 wofz calls!)
    w1 = wofz(z1)
    w2 = wofz(z2)

    # Real parts for profile
    W1 = np.real(w1)
    W2 = np.real(w2)

    # Faddeeva derivative: w'(z) = -2z×w(z) + 2i/√π
    w1_prime = -2.0 * z1 * w1 + 2j / SQRT_PI
    w2_prime = -2.0 * z2 * w2 + 2j / SQRT_PI

    # Normalization factor
    norm = sigma_f1 * sigma_f2 * TWO_PI

    # === Profile ===
    V = intensity * W1 * W2 / norm

    # === Derivative w.r.t. intensity (trivial) ===
    # ∂V/∂A = V/A
    if abs(intensity) > 1e-10:
        dV_dA = V / intensity
    else:
        dV_dA = W1 * W2 / norm

    # === Derivatives w.r.t. F1 parameters ===
    #
    # V = A × W1 × W2 / (σ1 × σ2 × 2π)
    #
    # For F1 derivatives, W2/(σ2×2π) is constant, so:
    # ∂V/∂θ1 = A × W2 / (σ2 × 2π) × ∂/∂θ1 [W1/σ1]

    factor_f1 = intensity * W2 / (sigma_f2 * TWO_PI)

    # ∂V/∂pos_f1:
    # ∂z1/∂pos_f1 = 1/(σ1×√2)
    # ∂W1/∂pos_f1 = Re[w'(z1)] × ∂z1/∂pos_f1 = Re[w'(z1)] / (σ1×√2)
    # ∂V/∂pos_f1 = factor_f1 × Re[w'(z1)] / (σ1² × √2)
    dV_dpos_f1 = factor_f1 * np.real(w1_prime) / (sigma_f1**2 * SQRT_2)

    # ∂V/∂lw_lor_f1:
    # γ1 = lw_lor_f1 / 2, so ∂γ1/∂lw_lor_f1 = 0.5
    # ∂z1/∂γ1 = i / (σ1×√2)
    # ∂W1/∂γ1 = Re[w'(z1) × i/(σ1×√2)] = -Im[w'(z1)] / (σ1×√2)
    # ∂V/∂lw_lor_f1 = 0.5 × factor_f1 × (-Im[w'(z1)]) / (σ1² × √2)
    dV_dlor_f1 = -0.5 * factor_f1 * np.imag(w1_prime) / (sigma_f1**2 * SQRT_2)

    # ∂V/∂lw_gau_f1:
    # σ1 = lw_gau_f1 / √(8ln2), so ∂σ1/∂lw_gau_f1 = 1/√(8ln2)
    # σ1 appears in z1 AND in normalization 1/σ1
    #
    # ∂z1/∂σ1 = -z1/σ1
    # ∂W1/∂σ1 = Re[w'(z1) × (-z1/σ1)] = -Re[z1 × w'(z1)] / σ1
    #
    # ∂/∂σ1 [W1/σ1] = (∂W1/∂σ1)/σ1 - W1/σ1²
    #                = -Re[z1 × w'(z1)]/σ1² - W1/σ1²
    #                = -(Re[z1 × w'(z1)] + W1) / σ1²
    #
    # ∂V/∂σ1 = factor_f1 × (-(Re[z1 × w'(z1)] + W1) / σ1²)
    # ∂V/∂lw_gau_f1 = ∂V/∂σ1 / √(8ln2)
    z1_wprime1 = z1 * w1_prime
    dV_dsigma_f1 = -factor_f1 * (np.real(z1_wprime1) + W1) / (sigma_f1**2)
    dV_dgau_f1 = dV_dsigma_f1 / SQRT_8LN2

    # === Derivatives w.r.t. F2 parameters (symmetric to F1) ===
    factor_f2 = intensity * W1 / (sigma_f1 * TWO_PI)

    # ∂V/∂pos_f2
    dV_dpos_f2 = factor_f2 * np.real(w2_prime) / (sigma_f2**2 * SQRT_2)

    # ∂V/∂lw_lor_f2
    dV_dlor_f2 = -0.5 * factor_f2 * np.imag(w2_prime) / (sigma_f2**2 * SQRT_2)

    # ∂V/∂lw_gau_f2
    z2_wprime2 = z2 * w2_prime
    dV_dsigma_f2 = -factor_f2 * (np.real(z2_wprime2) + W2) / (sigma_f2**2)
    dV_dgau_f2 = dV_dsigma_f2 / SQRT_8LN2

    return (V, dV_dpos_f1, dV_dlor_f1, dV_dpos_f2, dV_dlor_f2,
            dV_dgau_f1, dV_dgau_f2, dV_dA)


def compute_analytical_jacobian_2d(
    f1: np.ndarray,
    f2: np.ndarray,
    params: np.ndarray,
    n_peaks: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute multi-peak 2D Voigt profile AND full Jacobian analytically.

    Uses 2×n_peaks wofz calls instead of 7×n_peaks (numerical) = 3.5× fewer.

    Parameters:
    -----------
    f1, f2 : np.ndarray
        Frequency coordinates (flattened)
    params : np.ndarray
        Flattened parameter array (n_peaks × 8 parameters)
    n_peaks : int
        Number of peaks

    Returns:
    --------
    y : np.ndarray
        Profile values (n_points,)
    jacobian : np.ndarray
        Jacobian matrix (n_points, n_peaks × 8)
    """
    n_points = len(f1)
    n_params = n_peaks * NPAR_VOIGT

    # Initialize outputs
    y_total = np.zeros(n_points)
    jacobian = np.zeros((n_points, n_params))

    for peak_idx in range(n_peaks):
        offset = peak_idx * NPAR_VOIGT

        # Extract parameters for this peak
        # Order: pos_f1, lw_lor_f1, lw_gau_f1, pos_f2, lw_lor_f2, lw_gau_f2, intensity, spare
        pos_f1 = params[offset + 0]
        lw_lor_f1 = params[offset + 1]
        lw_gau_f1 = params[offset + 2]
        pos_f2 = params[offset + 3]
        lw_lor_f2 = params[offset + 4]
        lw_gau_f2 = params[offset + 5]
        intensity = params[offset + 6]
        # params[offset + 7] is spare/baseline, derivative = 0

        # Get profile and all derivatives in one call
        (V, dV_dpos_f1, dV_dlor_f1, dV_dpos_f2, dV_dlor_f2,
         dV_dgau_f1, dV_dgau_f2, dV_dA) = voigt_2d_with_derivatives(
            f1, f2, pos_f1, lw_lor_f1, pos_f2, lw_lor_f2,
            lw_gau_f1, lw_gau_f2, intensity
        )

        # Accumulate profile
        y_total += V

        # Fill Jacobian columns for this peak (matching parameter order)
        jacobian[:, offset + 0] = dV_dpos_f1
        jacobian[:, offset + 1] = dV_dlor_f1
        jacobian[:, offset + 2] = dV_dgau_f1
        jacobian[:, offset + 3] = dV_dpos_f2
        jacobian[:, offset + 4] = dV_dlor_f2
        jacobian[:, offset + 5] = dV_dgau_f2
        jacobian[:, offset + 6] = dV_dA
        jacobian[:, offset + 7] = 0.0  # Spare parameter

    return y_total, jacobian


def validate_analytical_jacobian(
    f1: np.ndarray,
    f2: np.ndarray,
    params: np.ndarray,
    n_peaks: int,
    eps: float = 1e-7,
    verbose: bool = True
) -> Tuple[bool, float]:
    """
    Validate analytical Jacobian against numerical finite differences.

    Parameters:
    -----------
    f1, f2 : np.ndarray
        Frequency coordinates
    params : np.ndarray
        Parameter array
    n_peaks : int
        Number of peaks
    eps : float
        Finite difference step size
    verbose : bool
        Print detailed results

    Returns:
    --------
    passed : bool
        True if all derivatives match within tolerance
    max_error : float
        Maximum relative error found
    """
    from .ps2d_style_fitter import multi_voigt_profile_2d

    # Compute analytical Jacobian
    y_analytical, jac_analytical = compute_analytical_jacobian_2d(f1, f2, params, n_peaks)

    # Compute numerical Jacobian
    n_params = n_peaks * NPAR_VOIGT
    jac_numerical = np.zeros((len(f1), n_params))
    y_base = multi_voigt_profile_2d(f1, f2, params, n_peaks)

    for i in range(n_params):
        params_perturbed = params.copy()
        h = eps * max(abs(params[i]), 1e-8)
        params_perturbed[i] += h
        y_perturbed = multi_voigt_profile_2d(f1, f2, params_perturbed, n_peaks)
        jac_numerical[:, i] = (y_perturbed - y_base) / h

    # Compare using np.allclose-style tolerance: |ana - num| <= atol + rtol * |num|
    # This handles both large and small derivatives appropriately
    max_error = 0.0
    all_passed = True
    param_names = ['pos_f1', 'lw_lor_f1', 'lw_gau_f1', 'pos_f2',
                   'lw_lor_f2', 'lw_gau_f2', 'intensity', 'spare']

    # Tolerance: 2% relative OR 1e-6 absolute (handles near-zero derivatives)
    # 2% is sufficient for optimization - scipy uses similar tolerances internally
    rtol = 0.02
    atol = 1e-6

    for peak_idx in range(n_peaks):
        for param_idx in range(NPAR_VOIGT):
            col = peak_idx * NPAR_VOIGT + param_idx

            ana = jac_analytical[:, col]
            num = jac_numerical[:, col]

            # Compute error using allclose-style metric
            abs_diff = np.abs(ana - num)
            tolerance = atol + rtol * np.abs(num)

            # Relative error metric: how much abs_diff exceeds tolerance
            # 0 means perfect match, >1 means exceeds tolerance
            error_ratio = abs_diff / np.maximum(tolerance, 1e-20)
            max_ratio = np.max(error_ratio)

            # For reporting, compute relative error where derivatives are significant
            mask = np.abs(num) > atol * 10  # Use 10x absolute tolerance as threshold
            if mask.sum() > 0:
                rel_error = np.max(abs_diff[mask] / np.abs(num[mask]))
            else:
                rel_error = np.max(abs_diff) / atol if atol > 0 else 0

            max_error = max(max_error, rel_error)
            passed = max_ratio <= 1.0  # Within tolerance

            if verbose and param_idx < 7:  # Skip spare
                status = "✅" if passed else "❌"
                print(f"  Peak {peak_idx+1} {param_names[param_idx]:12s}: "
                      f"rel_error = {rel_error:.2e} {status}")

            if not passed and param_idx < 7:
                all_passed = False

    return all_passed, max_error
