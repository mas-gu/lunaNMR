"""
PS2D Exact Overlap Integration Wrapper
========================================

This module provides integration functions that connect the PS2D exact overlap detector
(ellipsoid-based geometric overlap detection + transitive closure grouping) with the
lunaNMR workflow.

This is the missing piece that enables the EXACT C++ PS2D algorithm in lunaNMR GUI.

Usage:
------
Instead of calling fit_overlapping_peaks_ps2d_style() directly, use:
    fit_peaks_with_exact_overlap_detection()

This ensures peaks are grouped by geometric overlap BEFORE fitting, exactly matching
C++ PS2D behavior from peak.cpp:1693-1797 (autoOverLapPeaks).

Author: Guillaume Mas
Date: 2025-10-05
Version: 1.0 - INTEGRATION_WRAPPER
"""

import numpy as np
from typing import List, Dict, Optional

# Import the exact overlap detector components
try:
    from lunaNMR.core.ps2d_exact_overlap_detector import (
        EllipsoidOverlapDetector,
        TransitiveOverlapGrouper,
        Ps2dMultiPeakDeconvolver,
        detect_and_fit_overlapping_peaks
    )
    OVERLAP_DETECTOR_AVAILABLE = True
except ImportError as e:
    OVERLAP_DETECTOR_AVAILABLE = False
    print(f"⚠️ Cannot import ps2d_exact_overlap_detector: {e}")

# Import the multi-peak fitter (used for fitting each group)
try:
    from lunaNMR.core.ps2d_multi_peak_fitter import (
        Ps2dMultiPeakFitter,
        fit_overlapping_peaks_ps2d_style,
        multi_voigt_profile_1d
    )
    MULTI_PEAK_FITTER_AVAILABLE = True
except ImportError as e:
    MULTI_PEAK_FITTER_AVAILABLE = False
    print(f"⚠️ Cannot import ps2d_multi_peak_fitter: {e}")


def fit_peaks_with_exact_overlap_detection(
    x_data: np.ndarray,
    y_data: np.ndarray,
    detected_peaks: List[Dict],
    dimension: str = 'x',
    target_position: Optional[float] = None,
    fix_linewidths: bool = False,
    fix_positions: bool = False,
    verbose: bool = False,
    lw_lorentz_1h: Optional[float] = None,
    lw_gauss_1h: Optional[float] = None,
    lw_lorentz_15n: Optional[float] = None,
    lw_gauss_15n: Optional[float] = None,
    use_exact_overlap_detection: bool = True
) -> Dict:
    """
    Fit peaks using EXACT C++ PS2D overlap detection algorithm

    This function implements the complete C++ PS2D workflow:
    1. Ellipsoid overlap detection (geometric intersection test)
    2. Transitive closure grouping (if A overlaps B and B overlaps C → all grouped)
    3. Simultaneous fitting of each overlap group (5-stage deconvolution)
    4. Individual fitting of isolated peaks

    This is the EXACT replica of C++  workflow from:
    - peak.cpp:1693-1797 (autoOverLapPeaks function)
    - peak.cpp:1800-1819 (overlapLoop recursive function)
    - peakfit.cpp:347-455 (fitGlobal simultaneous fitting)

    Parameters
    ----------
    x_data : np.ndarray
        Frequency axis (ppm)
    y_data : np.ndarray
        Intensity data (must be baseline-corrected)
    detected_peaks : List[Dict]
        All detected peaks in the spectrum
        Each peak must have: 'position', 'height' (or 'intensity')
        Optional: 'lw_lorentz', 'lw_gauss', 'radius'
    dimension : str, optional
        'x' (1H) or 'y' (15N/13C) for nucleus-specific settings
    fix_linewidths : bool, optional
        Fix linewidths at initial estimates (C++ fixLW flag)
    fix_positions : bool, optional
        Fix positions at initial estimates (C++ fixPos flag)
    verbose : bool, optional
        Print detailed progress
    lw_lorentz_1h : float, optional
        User-provided Lorentzian linewidth for 1H (ppm)
    lw_gauss_1h : float, optional
        User-provided Gaussian linewidth for 1H (ppm)
    lw_lorentz_15n : float, optional
        User-provided Lorentzian linewidth for 15N/13C (ppm)
    lw_gauss_15n : float, optional
        User-provided Gaussian linewidth for 15N/13C (ppm)
    use_exact_overlap_detection : bool, optional
        Use ellipsoid overlap detection (True) or fit all peaks together (False)
        Default: True (EXACT C++ behavior)

    Returns
    -------
    Dict : Fitting results compatible with lunaNMR format
        {
            'success': bool,
            'peaks': List[Dict],  # All fitted peaks
            'r_squared': float,
            'fitted_curve': np.ndarray,
            'overlap_groups': List[List[int]],  # NEW: overlap group information
            'n_overlap_groups': int,  # NEW: number of groups found
            'n_isolated_peaks': int,  # NEW: number of isolated peaks
            'method': 'ps2d_exact_overlap' or 'ps2d_multi_peak_5_stage'
        }

    Notes
    -----
    - If use_exact_overlap_detection=False, falls back to fitting all peaks together
      (original behavior, NOT exact C++ PS2D)
    - If overlap detector not available, automatically falls back to original method
    - Peak radii for ellipsoid detection: estimated from linewidths if not provided
      (radius ≈ 1.5 × FWHM, matching C++ PS2D convention)
    """

    if verbose:
        print(f"\n{'='*70}")
        print(f"PS2D Exact Overlap Detection + Fitting")
        print(f"{'='*70}")
        print(f"Input: {len(detected_peaks)} peaks, dimension={dimension}")
        print(f"Exact overlap detection: {use_exact_overlap_detection}")

    # ========================================================================
    # FALLBACK: If exact overlap detection disabled or unavailable
    # ========================================================================
    if not use_exact_overlap_detection or not OVERLAP_DETECTOR_AVAILABLE:
        if verbose:
            reason = "disabled by user" if not use_exact_overlap_detection else "not available"
            print(f"⚠️ Exact overlap detection {reason} - using standard multi-peak fitting")

        # Fall back to standard multi-peak fitting (fits all peaks together)
        if MULTI_PEAK_FITTER_AVAILABLE:
            result = fit_overlapping_peaks_ps2d_style(
                x_data=x_data,
                y_data=y_data,
                detected_peaks=detected_peaks,
                dimension=dimension,
                fix_linewidths=fix_linewidths,
                fix_positions=fix_positions,
                verbose=verbose,
                lw_lorentz_1h=lw_lorentz_1h,
                lw_gauss_1h=lw_gauss_1h,
                lw_lorentz_15n=lw_lorentz_15n,
                lw_gauss_15n=lw_gauss_15n
            )
            # Add metadata
            result['overlap_groups'] = [list(range(len(detected_peaks)))]  # All peaks in one group
            result['n_overlap_groups'] = 1
            result['n_isolated_peaks'] = 0
            return result
        else:
            raise RuntimeError("Multi-peak fitter not available")

    # ========================================================================
    # STEP 1: PREPARE PEAKS WITH RADII FOR ELLIPSOID DETECTION
    # ========================================================================
    # C++ PS2D uses peak radii (radF1, radF2, radF3) for ellipsoid overlap detection
    # If not provided, estimate from linewidths: radius ≈ 1.5 × FWHM

    if verbose:
        print(f"\n--- Step 1: Preparing peaks for ellipsoid detection ---")

    # Estimate typical linewidths for radius calculation
    # Use user-provided values if available, otherwise use nucleus defaults
    if dimension == 'x' and lw_gauss_1h is not None:
        typical_lw = lw_gauss_1h
    elif dimension == 'y' and lw_gauss_15n is not None:
        typical_lw = lw_gauss_15n
    elif dimension == 'x':
        typical_lw = 0.03  # 1H default (from Ps2dLinewidthEstimator)
    else:  # dimension == 'y'
        typical_lw = 0.3   # 15N default (from Ps2dLinewidthEstimator)

    # Prepare peaks with radii
    peaks_with_radii = []
    for i, peak in enumerate(detected_peaks):
        # Extract position
        pos = peak.get('position', peak.get('pos', 0.0))

        # Get or estimate linewidth for radius calculation
        if 'lw_gauss' in peak:
            lw = peak['lw_gauss']
        elif 'linewidth' in peak:
            lw = peak['linewidth']
        else:
            lw = typical_lw

        # Calculate radius (C++ PS2D convention: radius ≈ 1.5 × FWHM)
        # This defines the ellipsoid extent for overlap detection
        radius = lw * 1.5

        # For 1D spectra: put radius in first dimension, zero in others
        # For 2D/3D spectra: would distribute across dimensions
        peak_with_radius = {
            'position': pos,  # 1D position
            'f1': pos,        # Also store as f1 (C++ PS2D convention)
            'f2': 0.0,        # Not used for 1D
            'f3': 0.0,        # Not used for 1D
            'radius': radius,
            'radF1': radius,
            'radF2': radius / 10.0,  # Approximate (not critical for 1D)
            'radF3': radius / 10.0,  # Approximate (not critical for 1D)
            'original_index': i,
            **peak  # Keep all original data
        }
        peaks_with_radii.append(peak_with_radius)

    if verbose:
        print(f"   Prepared {len(peaks_with_radii)} peaks with radii")
        print(f"   Typical radius: {typical_lw * 1.5:.4f} ppm (1.5 × {typical_lw:.4f} FWHM)")

    # ========================================================================
    # STEP 2: DETECT OVERLAP GROUPS USING EXACT C++ ALGORITHM
    # ========================================================================
    # Uses ellipsoid intersection test + transitive closure grouping
    # Exact replica of C++ peak.cpp:1693-1797

    if verbose:
        print(f"\n--- Step 2: Detecting overlap groups (exact C++ algorithm) ---")

    # Create overlap detector (ellipsoid geometric intersection)
    detector = EllipsoidOverlapDetector(
        verbose=verbose,
        use_preliminary_only=False  # Use full ellipsoid test (exact C++)
    )

    # Create grouper (recursive transitive closure)
    grouper = TransitiveOverlapGrouper(detector, verbose=verbose)

    # Run auto overlap detection (EXACT C++ algorithm)
    overlap_groups = grouper.auto_overlap_peaks(peaks_with_radii)

    # Identify isolated peaks (not in any group)
    peaks_in_groups = set()
    for group in overlap_groups:
        peaks_in_groups.update(group)

    isolated_peak_indices = [
        i for i in range(len(peaks_with_radii))
        if i not in peaks_in_groups
    ]

    if verbose:
        print(f"\n   Overlap detection results:")
        print(f"   - Total peaks: {len(peaks_with_radii)}")
        print(f"   - Overlap groups: {len(overlap_groups)}")
        print(f"   - Isolated peaks: {len(isolated_peak_indices)}")
        for i, group in enumerate(overlap_groups):
            positions = [peaks_with_radii[idx]['position'] for idx in group]
            print(f"      Group {i+1}: {len(group)} peaks at {[f'{p:.4f}' for p in positions]}")

    # ========================================================================
    # STEP 3: FIT EACH OVERLAP GROUP SIMULTANEOUSLY
    # ========================================================================
    # Each group is fitted as a SINGLE simultaneous optimization
    # Exact replica of C++ peakfit.cpp:347-455 (fitGlobal)

    if verbose:
        print(f"\n--- Step 3: Fitting overlap groups (simultaneous deconvolution) ---")

    all_fitted_peaks = []
    all_fit_metadata = []

    # Create deconvolver (wraps Ps2dMultiPeakFitter with exact C++ settings)
    deconvolver = Ps2dMultiPeakDeconvolver(verbose=verbose, dimension=dimension)

    # Fit each overlap group
    for group_idx, group_indices in enumerate(overlap_groups):
        if verbose:
            print(f"\n   Fitting overlap group {group_idx + 1}/{len(overlap_groups)} ({len(group_indices)} peaks)")

        # Extract peaks in this group
        group_peaks = [peaks_with_radii[i] for i in group_indices]

        # Prepare fix flags for this group
        fix_pos_group = [fix_positions] * len(group_indices)
        fix_lw_group = [fix_linewidths] * len(group_indices)

        # Fit group (ALL peaks SIMULTANEOUSLY)
        try:
            group_result = deconvolver.fit_overlapping_group(
                x_data=x_data,
                y_data=y_data,
                peak_group=group_peaks,
                fix_positions=fix_pos_group,
                fix_linewidths=fix_lw_group
            )

            if group_result.get('success'):
                # Extract fitted peaks
                for fitted_peak in group_result['peaks']:
                    all_fitted_peaks.append(fitted_peak)

                all_fit_metadata.append({
                    'group_size': len(group_indices),
                    'r_squared': group_result.get('r_squared', 0),
                    'method': 'simultaneous_deconvolution'
                })

                if verbose:
                    print(f"      ✅ Group {group_idx + 1} fitted: R² = {group_result.get('r_squared', 0):.6f}")
            else:
                if verbose:
                    print(f"      ⚠️ Group {group_idx + 1} fit failed")

        except Exception as e:
            if verbose:
                print(f"      ❌ Group {group_idx + 1} exception: {e}")

    # ========================================================================
    # STEP 4: FIT ISOLATED PEAKS (NOT IN ANY GROUP)
    # ========================================================================
    # Isolated peaks are fitted individually using the same 5-stage algorithm

    if len(isolated_peak_indices) > 0:
        if verbose:
            print(f"\n--- Step 4: Fitting isolated peaks ({len(isolated_peak_indices)} peaks) ---")

        for peak_idx in isolated_peak_indices:
            peak = peaks_with_radii[peak_idx]

            if verbose:
                print(f"   Fitting isolated peak at {peak['position']:.4f} ppm")

            # Fit single peak as a "group of 1"
            try:
                single_result = deconvolver.fit_overlapping_group(
                    x_data=x_data,
                    y_data=y_data,
                    peak_group=[peak],
                    fix_positions=[fix_positions],
                    fix_linewidths=[fix_linewidths]
                )

                if single_result.get('success'):
                    all_fitted_peaks.extend(single_result['peaks'])
                    all_fit_metadata.append({
                        'group_size': 1,
                        'r_squared': single_result.get('r_squared', 0),
                        'method': 'isolated_peak'
                    })

                    if verbose:
                        print(f"      ✅ Isolated peak fitted: R² = {single_result.get('r_squared', 0):.6f}")
                else:
                    if verbose:
                        print(f"      ⚠️ Isolated peak fit failed")

            except Exception as e:
                if verbose:
                    print(f"      ❌ Isolated peak exception: {e}")

    # ========================================================================
    # STEP 5: COMBINE RESULTS AND COMPUTE GLOBAL METRICS
    # ========================================================================

    if verbose:
        print(f"\n--- Step 5: Computing global fit quality ---")

    # Reconstruct full fitted curve by summing all peaks
    if len(all_fitted_peaks) > 0:
        # Build parameter vector for all fitted peaks
        all_params = []
        for peak in all_fitted_peaks:
            all_params.extend([
                peak['position'],
                peak['lw_lorentz'],
                peak['lw_gauss'],
                peak['intensity']
            ])
        all_params = np.array(all_params)

        # Compute fitted curve
        fitted_curve = multi_voigt_profile_1d(x_data, all_params, len(all_fitted_peaks))

        # Compute global R²
        residuals = y_data - fitted_curve
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_data - np.mean(y_data))**2)
        r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

        if verbose:
            print(f"   Total fitted peaks: {len(all_fitted_peaks)}")
            print(f"   Global R²: {r_squared:.6f}")

        success = True
    else:
        # No peaks fitted successfully
        fitted_curve = np.zeros_like(y_data)
        r_squared = 0.0
        success = False

        if verbose:
            print(f"   ⚠️ No peaks fitted successfully")

    # ========================================================================
    # RETURN RESULTS IN LUNANMR-COMPATIBLE FORMAT
    # ========================================================================

    result = {
        'success': success,
        'peaks': all_fitted_peaks,
        'r_squared': r_squared,
        'fitted_curve': fitted_curve,
        'residuals': y_data - fitted_curve if success else y_data,
        'overlap_groups': overlap_groups,
        'n_overlap_groups': len(overlap_groups),
        'n_isolated_peaks': len(isolated_peak_indices),
        'fit_metadata': all_fit_metadata,
        'method': 'ps2d_exact_overlap',
        'baseline': 0.0  # NO baseline (exact C++ PS2D)
    }

    if verbose:
        print(f"\n{'='*70}")
        print(f"PS2D Exact Overlap Detection Complete")
        print(f"{'='*70}\n")

    return result


# ============================================================================
# CONVENIENCE ALIASES
# ============================================================================

# Alias for backward compatibility
fit_peaks_with_ps2d_exact_overlap = fit_peaks_with_exact_overlap_detection

# Export all public functions
__all__ = [
    'fit_peaks_with_exact_overlap_detection',
    'fit_peaks_with_ps2d_exact_overlap',
    'OVERLAP_DETECTOR_AVAILABLE',
    'MULTI_PEAK_FITTER_AVAILABLE'
]
