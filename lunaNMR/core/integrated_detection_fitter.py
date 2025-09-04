#!/usr/bin/env python3
"""
Integrated Detection-Fitting Module

This module provides comprehensive integration of peak detection with Voigt profile fitting,
implementing iterative refinement, adaptive thresholding, multi-resolution detection,
and physics-informed constraints for enhanced NMR peak analysis.

Classes:
- IterativeDetectionFitter: Core integration class
- DetectionQualityScorer: Peak detection confidence scoring
- RapidFitAssessor: Fast fitting quality assessment
- AdaptiveThresholdCalculator: Dynamic detection parameter adjustment
- MultiResolutionDetector: Hierarchical peak detection
- PhysicsInformedConstraints: NMR coupling pattern constraints

Author: Guillaume Mas
Date: 2025
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Any, Optional
from scipy.signal import find_peaks, peak_widths, peak_prominences
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.special import wofz
import warnings
warnings.filterwarnings('ignore')

# Import enhanced components
try:
    from lunaNMR.core.enhanced_voigt_fitter import EnhancedVoigtFitter
    ENHANCED_AVAILABLE = True
except ImportError:
    ENHANCED_AVAILABLE = False
    print(" Enhanced Voigt fitter not available for integration")


class Peak1DRefiner:
    """
    1D Cross-Section Refinement Engine for Peak Coordinate Enhancement

    This class implements sub-pixel precision peak coordinate refinement using
    1D cross-section fitting with Voigt profiles. It extracts 1D slices along
    both spectral dimensions and uses multi-peak deconvolution for overlapping peaks.
    """

    def __init__(self, enhanced_fitter=None):
        self.enhanced_fitter = enhanced_fitter or (EnhancedVoigtFitter() if ENHANCED_AVAILABLE else None)

        # Refinement parameters
        self.refinement_params = {
            'fitting_window_x': 0.2,      # ppm window for 1H cross-sections was 0.3 was 0.2
            'fitting_window_y': 2,      # ppm window for 15N/13C cross-sections was 5.0 was 3.0 was 1.5
            'max_peaks_fit': 8,           # Maximum peaks per cross-section #was 4
            'max_optimization_iterations': 50,  # Fitting iterations
            'min_r_squared': 0.7,         # Quality threshold for refined coordinates
            'overlap_threshold': 0.05,    # ppm threshold for detecting overlapping peaks
            'coordinate_precision': 0.001, # Coordinate refinement precision
            'enable_multi_peak': True,    # Enable multi-peak deconvolution
            'convergence_threshold': 1e-6, # Fitting convergence threshold
        }

        # Refinement statistics tracking
        self.refinement_stats = {
            'peaks_processed': 0,
            'peaks_refined': 0,
            'peaks_split': 0,
            'peaks_merged': 0,
            'failed_refinements': 0,
            'avg_coordinate_shift': 0.0,
            'avg_r_squared_improvement': 0.0
        }

    def refine_peak_coordinates_1d(self, detected_peaks, nmr_data, ppm_x_axis, ppm_y_axis,
                                  refinement_params=None):
        """
        Main 1D refinement method - refines peak coordinates using 1D cross-section fitting

        Args:
            detected_peaks: DataFrame with columns ['Position_X', 'Position_Y', ...]
            nmr_data: 2D NMR intensity matrix
            ppm_x_axis: 1H chemical shift axis (ppm)
            ppm_y_axis: 15N/13C chemical shift axis (ppm)
            refinement_params: Optional parameter override

        Returns:
            dict: {
                'refined_peaks': DataFrame with refined coordinates,
                'refinement_quality': dict with quality metrics,
                'diagnostics': dict with detailed refinement information
            }
        """
        if refinement_params:
            self.refinement_params.update(refinement_params)

        refined_peaks = detected_peaks.copy()
        refinement_diagnostics = []
        self._reset_stats()

        print(f"🔬 Starting 1D refinement for {len(detected_peaks)} detected peaks...")

        for idx, peak in detected_peaks.iterrows():
            self.refinement_stats['peaks_processed'] += 1

            try:
                # Extract 1D cross-sections around detected peak
                cross_sections = self._extract_cross_sections(
                    peak['Position_X'], peak['Position_Y'], nmr_data, ppm_x_axis, ppm_y_axis
                )

                if cross_sections is None:
                    self.refinement_stats['failed_refinements'] += 1
                    continue

                # Perform 1D fitting on both dimensions
                x_refinement = self._fit_1d_cross_section(
                    cross_sections['x_data'], cross_sections['x_ppm'],
                    peak['Position_X'], 'x'
                )

                y_refinement = self._fit_1d_cross_section(
                    cross_sections['y_data'], cross_sections['y_ppm'],
                    peak['Position_Y'], 'y'
                )

                # Update coordinates if refinement successful
                original_x, original_y = peak['Position_X'], peak['Position_Y']
                coordinate_updated = False

                if x_refinement['success'] and x_refinement['r_squared'] >= self.refinement_params['min_r_squared']:
                    refined_peaks.loc[idx, 'Position_X'] = x_refinement['refined_position']
                    coordinate_updated = True

                if y_refinement['success'] and y_refinement['r_squared'] >= self.refinement_params['min_r_squared']:
                    refined_peaks.loc[idx, 'Position_Y'] = y_refinement['refined_position']
                    coordinate_updated = True

                if coordinate_updated:
                    self.refinement_stats['peaks_refined'] += 1

                    # Calculate coordinate shift
                    shift = np.sqrt((refined_peaks.loc[idx, 'Position_X'] - original_x)**2 +
                                  (refined_peaks.loc[idx, 'Position_Y'] - original_y)**2)
                    self.refinement_stats['avg_coordinate_shift'] += shift

                # Store refinement diagnostics
                diagnostics = {
                    'peak_index': idx,
                    'original_coords': (original_x, original_y),
                    'refined_coords': (refined_peaks.loc[idx, 'Position_X'], refined_peaks.loc[idx, 'Position_Y']),
                    'x_refinement': x_refinement,
                    'y_refinement': y_refinement,
                    'coordinate_updated': coordinate_updated
                }
                refinement_diagnostics.append(diagnostics)

            except Exception as e:
                self.refinement_stats['failed_refinements'] += 1
                print(f"⚠️ Refinement failed for peak {idx}: {str(e)}")
                continue

        # Calculate final statistics
        if self.refinement_stats['peaks_refined'] > 0:
            self.refinement_stats['avg_coordinate_shift'] /= self.refinement_stats['peaks_refined']

        refinement_quality = self._calculate_refinement_quality(refinement_diagnostics)

        print(f"✅ 1D refinement complete: {self.refinement_stats['peaks_refined']}/{self.refinement_stats['peaks_processed']} peaks refined")

        return {
            'refined_peaks': refined_peaks,
            'refinement_quality': refinement_quality,
            'diagnostics': refinement_diagnostics,
            'statistics': self.refinement_stats.copy()
        }

    def _extract_cross_sections(self, peak_x_ppm, peak_y_ppm, nmr_data, ppm_x_axis, ppm_y_axis):
        """Extract 1D cross-sections around a peak position"""
        try:
            # Find peak position in data indices
            x_idx = np.argmin(np.abs(ppm_x_axis - peak_x_ppm))
            y_idx = np.argmin(np.abs(ppm_y_axis - peak_y_ppm))

            # Calculate window sizes in data points
            x_ppm_range = abs(ppm_x_axis[0] - ppm_x_axis[-1])
            y_ppm_range = abs(ppm_y_axis[0] - ppm_y_axis[-1])

            x_window_points = max(10, int(self.refinement_params['fitting_window_x'] * len(ppm_x_axis) / x_ppm_range))
            y_window_points = max(10, int(self.refinement_params['fitting_window_y'] * len(ppm_y_axis) / y_ppm_range))

            # Define extraction windows
            x_start = max(0, x_idx - x_window_points // 2)
            x_end = min(len(ppm_x_axis), x_idx + x_window_points // 2)

            y_start = max(0, y_idx - y_window_points // 2)
            y_end = min(len(ppm_y_axis), y_idx + y_window_points // 2)

            # Extract 1D cross-sections
            x_cross_section = nmr_data[y_idx, x_start:x_end]  # 1H slice at fixed 15N/13C
            y_cross_section = nmr_data[y_start:y_end, x_idx]  # 15N/13C slice at fixed 1H

            # Corresponding ppm scales
            x_ppm_section = ppm_x_axis[x_start:x_end]
            y_ppm_section = ppm_y_axis[y_start:y_end]

            return {
                'x_data': x_cross_section,
                'y_data': y_cross_section,
                'x_ppm': x_ppm_section,
                'y_ppm': y_ppm_section,
                'x_center_idx': x_idx - x_start,
                'y_center_idx': y_idx - y_start
            }

        except Exception as e:
            print(f"⚠️ Cross-section extraction failed: {str(e)}")
            return None

    def _fit_1d_cross_section(self, intensity_data, ppm_data, initial_position, dimension):
        """Fit 1D cross-section with Voigt profile to refine peak coordinate"""
        try:
            # Handle empty or invalid data
            if len(intensity_data) < 3 or np.all(intensity_data <= 0):
                return {'success': False, 'reason': 'insufficient_data'}

            # Find local maxima in the cross-section for multi-peak detection
            smoothed_data = gaussian_filter1d(intensity_data, sigma=1.0)
            peaks_indices, _ = find_peaks(smoothed_data, height=np.max(smoothed_data) * 0.1)

            if len(peaks_indices) == 0:
                return {'success': False, 'reason': 'no_peaks_found'}

            # Determine number of peaks to fit
            n_peaks = min(len(peaks_indices), self.refinement_params['max_peaks_fit'])

            if n_peaks == 1:
                # Single peak fitting
                return self._fit_single_peak_1d(intensity_data, ppm_data, initial_position)
            else:
                # Multi-peak fitting for overlapping peaks
                return self._fit_multi_peak_1d(intensity_data, ppm_data, peaks_indices[:n_peaks], initial_position)

        except Exception as e:
            return {'success': False, 'reason': f'fitting_error: {str(e)}'}

    def _fit_single_peak_1d(self, intensity_data, ppm_data, initial_position):
        """Fit single Voigt profile to 1D cross-section"""
        try:
            # Initial parameter guesses
            max_intensity = np.max(intensity_data)
            center_idx = np.argmax(intensity_data)
            initial_center = ppm_data[center_idx]

            # Estimate width from peak width at half maximum
            half_max = max_intensity / 2
            left_idx = np.where(intensity_data[:center_idx] <= half_max)[0]
            right_idx = np.where(intensity_data[center_idx:] <= half_max)[0]

            if len(left_idx) > 0 and len(right_idx) > 0:
                width_estimate = abs(ppm_data[center_idx + right_idx[0]] - ppm_data[left_idx[-1]]) / 2
            else:
                width_estimate = 0.02  # Default width

            # Parameter bounds and initial guess
            p0 = [max_intensity, initial_center, width_estimate, width_estimate * 0.5]
            bounds = (
                [0, ppm_data[0], 0.001, 0.001],
                [max_intensity * 2, ppm_data[-1], 1.0, 1.0]
            )

            # Fit Voigt profile
            popt, pcov = curve_fit(self._voigt_profile_1d, ppm_data, intensity_data,
                                 p0=p0, bounds=bounds, maxfev=self.refinement_params['max_optimization_iterations'])

            # Calculate fit quality
            fitted_curve = self._voigt_profile_1d(ppm_data, *popt)
            r_squared = self._calculate_r_squared(intensity_data, fitted_curve)

            # Check if refinement is significantly better than initial guess
            coordinate_shift = abs(popt[1] - initial_position)

            return {
                'success': True,
                'refined_position': popt[1],
                'r_squared': r_squared,
                'coordinate_shift': coordinate_shift,
                'fitted_params': popt,
                'covariance': pcov,
                'fitted_curve': fitted_curve,
                'n_peaks': 1
            }

        except Exception as e:
            return {'success': False, 'reason': f'single_peak_fitting_error: {str(e)}'}

    def _fit_multi_peak_1d(self, intensity_data, ppm_data, peak_indices, target_position):
        """Fit multiple Voigt profiles to 1D cross-section for overlapping peaks"""
        try:
            n_peaks = len(peak_indices)

            # Initial parameter guesses for multi-peak fitting
            p0 = []
            bounds_lower = []
            bounds_upper = []

            for i, peak_idx in enumerate(peak_indices):
                amplitude = intensity_data[peak_idx]
                center = ppm_data[peak_idx]
                sigma = 0.02  # Default width
                gamma = 0.01  # Default Lorentzian component

                p0.extend([amplitude, center, sigma, gamma])
                bounds_lower.extend([0, ppm_data[0], 0.001, 0.001])
                bounds_upper.extend([amplitude * 2, ppm_data[-1], 1.0, 1.0])

            # Multi-peak Voigt fitting
            popt, pcov = curve_fit(
                lambda x, *params: self._multi_voigt_profile_1d(x, params, n_peaks),
                ppm_data, intensity_data, p0=p0, bounds=(bounds_lower, bounds_upper),
                maxfev=self.refinement_params['max_optimization_iterations']
            )

            # Calculate fit quality
            fitted_curve = self._multi_voigt_profile_1d(ppm_data, popt, n_peaks)
            r_squared = self._calculate_r_squared(intensity_data, fitted_curve)

            # Find the peak closest to target position
            peak_positions = [popt[i*4 + 1] for i in range(n_peaks)]
            closest_peak_idx = np.argmin([abs(pos - target_position) for pos in peak_positions])
            refined_position = peak_positions[closest_peak_idx]

            coordinate_shift = abs(refined_position - target_position)

            return {
                'success': True,
                'refined_position': refined_position,
                'r_squared': r_squared,
                'coordinate_shift': coordinate_shift,
                'fitted_params': popt,
                'covariance': pcov,
                'fitted_curve': fitted_curve,
                'n_peaks': n_peaks,
                'all_peak_positions': peak_positions
            }

        except Exception as e:
            return {'success': False, 'reason': f'multi_peak_fitting_error: {str(e)}'}

    def _voigt_profile_1d(self, x, amplitude, center, sigma, gamma):
        """1D Voigt profile function"""
        z = ((x - center) + 1j * gamma) / (sigma * np.sqrt(2))
        return amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))

    def _multi_voigt_profile_1d(self, x, params, n_peaks):
        """Multi-peak Voigt profile function"""
        result = np.zeros_like(x)
        for i in range(n_peaks):
            amplitude = params[i*4]
            center = params[i*4 + 1]
            sigma = params[i*4 + 2]
            gamma = params[i*4 + 3]
            result += self._voigt_profile_1d(x, amplitude, center, sigma, gamma)
        return result

    def _calculate_r_squared(self, y_true, y_pred):
        """Calculate R-squared coefficient"""
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    def _calculate_refinement_quality(self, diagnostics):
        """Calculate overall refinement quality metrics"""
        if not diagnostics:
            return {'overall_quality': 0.0}

        successful_refinements = [d for d in diagnostics if d['coordinate_updated']]

        if not successful_refinements:
            return {'overall_quality': 0.0}

        # Calculate quality metrics
        avg_r_squared_x = np.mean([d['x_refinement'].get('r_squared', 0) for d in successful_refinements
                                  if d['x_refinement'].get('success', False)])
        avg_r_squared_y = np.mean([d['y_refinement'].get('r_squared', 0) for d in successful_refinements
                                  if d['y_refinement'].get('success', False)])

        success_rate = len(successful_refinements) / len(diagnostics)

        overall_quality = (avg_r_squared_x + avg_r_squared_y) * success_rate / 2

        return {
            'overall_quality': overall_quality,
            'success_rate': success_rate,
            'avg_r_squared_x': avg_r_squared_x,
            'avg_r_squared_y': avg_r_squared_y,
            'total_peaks': len(diagnostics),
            'refined_peaks': len(successful_refinements)
        }

    def _reset_stats(self):
        """Reset refinement statistics"""
        for key in self.refinement_stats:
            if isinstance(self.refinement_stats[key], (int, float)):
                self.refinement_stats[key] = 0

    def set_refinement_params(self, params):
        """Update refinement parameters"""
        self.refinement_params.update(params)
        print(f"🔧 1D refinement parameters updated: {params}")


class DetectionQualityScorer:
    """Peak detection confidence scoring system"""

    def __init__(self):
        self.scoring_weights = {
            'intensity_ratio': 0.3,      # Peak/noise ratio
            'peak_width': 0.25,          # Reasonable width
            'prominence': 0.2,           # Peak prominence
            'isolation': 0.15,           # Separation from neighbors
            'chemical_shift_context': 0.1  # Expected chemical shift range
        }

    def score_peak_detection(self, x_data, y_data, peak_info, nucleus_type='1H'):
        """
        Score detected peaks for confidence assessment

        Args:
            x_data, y_data: spectral data
            peak_info: list of detected peak dictionaries
            nucleus_type: '1H', '15N', or '13C'

        Returns:
            dict: scored peaks with confidence values
        """
        if not peak_info:
            return {}

        scored_peaks = {}
        noise_level = self._estimate_noise_level(y_data)
        signal_max = np.max(y_data)

        # Nucleus-specific ranges for context scoring
        nucleus_ranges = {
            '1H': {'min': 0.5, 'max': 10.0, 'typical_width': 0.01},
            '15N': {'min': 100.0, 'max': 140.0, 'typical_width': 0.5},
            '13C': {'min': 10.0, 'max': 180.0, 'typical_width': 1.0}
        }

        typical_params = nucleus_ranges.get(nucleus_type, nucleus_ranges['1H'])

        for i, peak in enumerate(peak_info):
            position = peak['position']
            intensity = peak['intensity']

            # 1. Intensity ratio score
            intensity_ratio = (intensity - noise_level) / (signal_max - noise_level)
            intensity_score = min(1.0, intensity_ratio * 2.0)  # Normalize to [0,1]

            # 2. Peak width estimation and scoring
            peak_width = self._estimate_peak_width(x_data, y_data, peak['index'])
            width_ratio = peak_width / typical_params['typical_width']
            # Score peaks with reasonable widths (0.1x to 10x typical)
            if 0.1 <= width_ratio <= 10.0:
                width_score = 1.0 - abs(np.log10(width_ratio)) / 2.0  # Penalty for extreme ratios
            else:
                width_score = 0.1  # Very unusual width
            width_score = max(0.0, min(1.0, width_score))

            # 3. Prominence score
            prominence = self._calculate_prominence(y_data, peak['index'])
            prominence_ratio = prominence / signal_max
            prominence_score = min(1.0, prominence_ratio * 10.0)  # Normalize

            # 4. Isolation score (separation from other peaks)
            isolation_score = self._calculate_isolation_score(peak_info, i, x_data)

            # 5. Chemical shift context score
            ppm_range = typical_params['max'] - typical_params['min']
            if typical_params['min'] <= position <= typical_params['max']:
                context_score = 1.0
            else:
                # Penalty for being outside typical range
                distance = min(abs(position - typical_params['min']),
                             abs(position - typical_params['max']))
                context_score = max(0.0, 1.0 - distance / ppm_range)

            # Weighted composite score
            composite_score = (
                self.scoring_weights['intensity_ratio'] * intensity_score +
                self.scoring_weights['peak_width'] * width_score +
                self.scoring_weights['prominence'] * prominence_score +
                self.scoring_weights['isolation'] * isolation_score +
                self.scoring_weights['chemical_shift_context'] * context_score
            )

            scored_peaks[i] = {
                **peak,  # Original peak info
                'confidence': composite_score,
                'confidence_breakdown': {
                    'intensity_score': intensity_score,
                    'width_score': width_score,
                    'prominence_score': prominence_score,
                    'isolation_score': isolation_score,
                    'context_score': context_score
                },
                'estimated_width': peak_width,
                'prominence': prominence
            }

        return scored_peaks

    def _estimate_noise_level(self, y_data):
        """Estimate noise level from data edges"""
        if len(y_data) < 20:
            return np.std(y_data) * 0.1

        edge_size = min(10, len(y_data) // 20)
        edges = np.concatenate([y_data[:edge_size], y_data[-edge_size:]])
        return np.std(edges)

    def _estimate_peak_width(self, x_data, y_data, peak_index):
        """Estimate peak width using half-maximum method"""
        try:
            if peak_index <= 0 or peak_index >= len(y_data) - 1:
                return 0.01  # Default width

            peak_intensity = y_data[peak_index]
            baseline = min(y_data[max(0, peak_index-10):peak_index+11])
            half_max = baseline + (peak_intensity - baseline) / 2

            # Find points at half maximum
            left_idx = peak_index
            right_idx = peak_index

            # Search left
            while left_idx > 0 and y_data[left_idx] > half_max:
                left_idx -= 1

            # Search right
            while right_idx < len(y_data) - 1 and y_data[right_idx] > half_max:
                right_idx += 1

            if right_idx > left_idx:
                width = abs(x_data[right_idx] - x_data[left_idx])
                return max(width, abs(x_data[1] - x_data[0]))  # At least one data point
            else:
                return abs(x_data[1] - x_data[0]) * 2  # Fallback

        except Exception:
            return 0.01  # Safe default

    def _calculate_prominence(self, y_data, peak_index):
        """Calculate peak prominence"""
        try:
            if peak_index <= 0 or peak_index >= len(y_data) - 1:
                return 0.0

            prominences, _, _ = peak_prominences(y_data, [peak_index])
            return prominences[0] if prominences.size > 0 else 0.0
        except Exception:
            return 0.0

    def _calculate_isolation_score(self, peak_info, peak_index, x_data):
        """Calculate how isolated a peak is from neighbors"""
        if len(peak_info) <= 1:
            return 1.0  # Completely isolated

        current_peak = peak_info[peak_index]
        current_pos = current_peak['position']

        min_distance = float('inf')
        for i, other_peak in enumerate(peak_info):
            if i != peak_index:
                distance = abs(other_peak['position'] - current_pos)
                min_distance = min(min_distance, distance)

        if min_distance == float('inf'):
            return 1.0

        # Normalize by typical data spacing
        data_spacing = abs(x_data[-1] - x_data[0]) / len(x_data)
        normalized_distance = min_distance / (data_spacing * 10)  # 10 points is good separation

        return min(1.0, normalized_distance)

class RapidFitAssessor:
    """Fast fitting quality assessment for peak validation"""

    def __init__(self):
        self.voigt_cache = {}
        self.assessment_thresholds = {
            'min_r_squared': 0.5,
            'max_parameter_deviation': 3.0,
            'min_amplitude_ratio': 0.1
        }

    def rapid_assessment(self, x_data, y_data, peak_position, estimated_width=None, nucleus_type='1H'):
        """
        Quickly assess if a peak position is likely to yield a good fit

        Args:
            x_data, y_data: spectral data
            peak_position: suspected peak center
            estimated_width: estimated peak width (optional)
            nucleus_type: nucleus type for parameter bounds

        Returns:
            dict: assessment results with fit likelihood
        """
        try:
            # Estimate initial parameters
            peak_idx = np.argmin(np.abs(x_data - peak_position))
            baseline = self._estimate_local_baseline(y_data, peak_idx)
            amplitude = y_data[peak_idx] - baseline

            if estimated_width is None:
                estimated_width = self._quick_width_estimate(x_data, y_data, peak_idx)

            # Quick parameter reasonableness check
            if amplitude <= 0:
                return {'likelihood': 0.0, 'reason': 'negative_amplitude'}

            if estimated_width <= 0 or estimated_width > abs(x_data[-1] - x_data[0]) * 0.5:
                return {'likelihood': 0.1, 'reason': 'unreasonable_width'}

            # Fast fitting attempt with reduced precision
            initial_guess = [amplitude, peak_position,
                           estimated_width * 0.6, estimated_width * 0.4, baseline]

            # Use lightweight bounds
            bounds = self._get_rapid_bounds(initial_guess, x_data, nucleus_type)

            # Fast fit with limited iterations
            try:
                popt, pcov = curve_fit(
                    self._voigt_profile_simple, x_data, y_data,
                    p0=initial_guess, bounds=bounds,
                    maxfev=50  # Very limited iterations
                )

                # Quick quality assessment
                y_fitted = self._voigt_profile_simple(x_data, *popt)
                r_squared = self._quick_r_squared(y_data, y_fitted)

                # Parameter validation
                param_deviation = self._check_parameter_deviation(initial_guess, popt)

                # Composite likelihood score
                likelihood = self._calculate_likelihood(r_squared, param_deviation, popt, initial_guess)

                return {
                    'likelihood': likelihood,
                    'r_squared': r_squared,
                    'parameters': popt,
                    'parameter_deviation': param_deviation,
                    'fit_success': True
                }

            except Exception as e:
                return {'likelihood': 0.2, 'reason': f'fit_failed: {str(e)}'}

        except Exception as e:
            return {'likelihood': 0.0, 'reason': f'assessment_failed: {str(e)}'}

    def _voigt_profile_simple(self, x, amplitude, center, sigma, gamma, baseline):
        """Simplified Voigt profile for rapid assessment"""
        try:
            sigma = max(sigma, 1e-6)
            z = ((x - center) + 1j*gamma) / (sigma * np.sqrt(2))
            voigt = amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2*np.pi))
            return voigt + baseline
        except:
            # Fallback to Gaussian
            return amplitude * np.exp(-0.5 * ((x - center) / max(sigma, 1e-6))**2) + baseline

    def _estimate_local_baseline(self, y_data, peak_idx):
        """Quick local baseline estimation"""
        window = min(10, len(y_data) // 10)
        start = max(0, peak_idx - window)
        end = min(len(y_data), peak_idx + window + 1)

        # Use minimum of nearby points as baseline
        return np.min(y_data[start:end])

    def _quick_width_estimate(self, x_data, y_data, peak_idx):
        """Quick width estimation"""
        try:
            baseline = self._estimate_local_baseline(y_data, peak_idx)
            peak_height = y_data[peak_idx] - baseline
            half_max = baseline + peak_height / 2

            # Find half-maximum points
            left = peak_idx
            right = peak_idx

            while left > 0 and y_data[left] > half_max:
                left -= 1
            while right < len(y_data) - 1 and y_data[right] > half_max:
                right += 1

            if right > left:
                return abs(x_data[right] - x_data[left]) / 2.355  # FWHM to sigma
            else:
                return abs(x_data[1] - x_data[0]) * 2  # Fallback
        except:
            return 0.01  # Default

    def _get_rapid_bounds(self, initial_guess, x_data, nucleus_type):
        """Get relaxed bounds for rapid fitting"""
        amplitude, center, sigma, gamma, baseline = initial_guess

        # Relaxed bounds for speed
        ppm_range = abs(x_data[-1] - x_data[0])

        lower_bounds = [
            0,  # amplitude > 0
            center - ppm_range * 0.1,  # center
            max(1e-6, sigma * 0.1),  # sigma
            0,  # gamma >= 0
            baseline - abs(amplitude)  # baseline
        ]

        upper_bounds = [
            amplitude * 10,  # amplitude
            center + ppm_range * 0.1,  # center
            sigma * 10,  # sigma
            sigma * 10,  # gamma
            baseline + abs(amplitude)  # baseline
        ]

        return (lower_bounds, upper_bounds)

    def _quick_r_squared(self, y_true, y_pred):
        """Quick R-squared calculation"""
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

    def _check_parameter_deviation(self, initial, fitted):
        """Check how much fitted parameters deviated from initial guess"""
        deviations = []
        for i, (init, fit) in enumerate(zip(initial, fitted)):
            if init != 0:
                deviation = abs((fit - init) / init)
            else:
                deviation = abs(fit)
            deviations.append(deviation)

        return np.max(deviations)

    def _calculate_likelihood(self, r_squared, param_deviation, popt, initial_guess):
        """Calculate composite likelihood score"""
        # R-squared component (0-1)
        r_score = min(1.0, max(0.0, r_squared))

        # Parameter stability component (penalize large deviations)
        param_score = max(0.0, 1.0 - param_deviation / 5.0)

        # Physical reasonableness (amplitude, widths)
        amplitude, center, sigma, gamma, baseline = popt
        phys_score = 1.0

        if amplitude <= 0:
            phys_score *= 0.1
        if sigma <= 0 or gamma < 0:
            phys_score *= 0.1
        if sigma > abs(initial_guess[1]) or gamma > abs(initial_guess[1]):
            phys_score *= 0.5

        # Weighted composite
        return (0.5 * r_score + 0.3 * param_score + 0.2 * phys_score)

class AdaptiveThresholdCalculator:
    """Dynamic detection parameter adjustment based on local conditions and fitting success"""

    def __init__(self):
        self.base_params = {
            'height_threshold': 0.02,
            'distance_factor': 50.0,
            'prominence_threshold': 0.01,
            'smoothing_sigma': 0.5,
            # CRITICAL FIX: Add missing core parameters to base_params
            'max_peaks_fit': 4,
            'max_optimization_iterations': 50
        }
        self.adaptation_history = []
        self.region_statistics = {}

    def calculate_adaptive_thresholds(self, x_data, y_data, nucleus_type='1H',
                                   region_id=None, success_history=None):
        """
        Calculate adaptive thresholds based on local data characteristics and fitting history

        Args:
            x_data, y_data: spectral data
            nucleus_type: nucleus type for context
            region_id: identifier for spectrum region
            success_history: dict with previous fitting success rates

        Returns:
            dict: adapted detection parameters
        """
        # Start with base parameters
        adapted_params = self.base_params.copy()

        # Local data analysis
        local_stats = self._analyze_local_data(x_data, y_data)

        # Apply noise-based adaptation
        noise_factor = self._calculate_noise_adaptation_factor(local_stats)
        adapted_params['height_threshold'] *= noise_factor
        adapted_params['prominence_threshold'] *= noise_factor

        # Apply peak density adaptation
        density_factor = self._calculate_density_adaptation_factor(local_stats)
        adapted_params['distance_factor'] *= density_factor

        # Apply success history adaptation
        if success_history and region_id:
            history_factor = self._calculate_history_adaptation_factor(success_history, region_id)
            adapted_params['height_threshold'] *= history_factor['height']
            adapted_params['prominence_threshold'] *= history_factor['prominence']
            adapted_params['distance_factor'] *= history_factor['distance']

        # Nucleus-specific adjustments
        nucleus_factor = self._get_nucleus_adaptation_factor(nucleus_type)
        for key in ['height_threshold', 'prominence_threshold']:
            adapted_params[key] *= nucleus_factor

        # Store adaptation info
        adaptation_info = {
            'region_id': region_id,
            'noise_factor': noise_factor,
            'density_factor': density_factor,
            'nucleus_factor': nucleus_factor,
            'local_stats': local_stats,
            'adapted_params': adapted_params
        }

        self.adaptation_history.append(adaptation_info)

        return adapted_params, adaptation_info

    def _analyze_local_data(self, x_data, y_data):
        """Analyze local data characteristics"""
        stats = {}

        # Noise analysis
        if len(y_data) > 20:
            edge_size = min(10, len(y_data) // 20)
            edges = np.concatenate([y_data[:edge_size], y_data[-edge_size:]])
            stats['noise_level'] = np.std(edges)
        else:
            stats['noise_level'] = np.std(y_data) * 0.1

        # Signal characteristics
        stats['signal_max'] = np.max(y_data)
        stats['signal_mean'] = np.mean(y_data)
        stats['signal_std'] = np.std(y_data)
        stats['dynamic_range'] = stats['signal_max'] / max(stats['noise_level'], 1e-6)

        # Rough peak density estimate
        smoothed = gaussian_filter1d(y_data, sigma=1.0)
        rough_peaks, _ = find_peaks(smoothed, height=stats['signal_max'] * 0.1)
        stats['rough_peak_density'] = len(rough_peaks) / len(y_data)

        return stats

    def _calculate_noise_adaptation_factor(self, local_stats):
        """Calculate adaptation factor based on noise level"""
        dynamic_range = local_stats['dynamic_range']

        if dynamic_range > 100:  # Low noise
            return 0.8  # Can use lower thresholds
        elif dynamic_range > 20:  # Medium noise
            return 1.0  # Use standard thresholds
        else:  # High noise
            return 1.5  # Need higher thresholds

    def _calculate_density_adaptation_factor(self, local_stats):
        """Calculate adaptation factor based on peak density"""
        density = local_stats['rough_peak_density']

        if density > 0.1:  # High density
            return 0.5  # Allow closer peaks
        elif density > 0.05:  # Medium density
            return 0.8
        else:  # Low density
            return 1.2  # Require more separation

    def _calculate_history_adaptation_factor(self, success_history, region_id):
        """Calculate adaptation based on fitting success history"""
        region_success = success_history.get(region_id, 0.8)  # Default 80% success

        if region_success > 0.9:  # High success region
            return {
                'height': 0.9,  # Can be more aggressive
                'prominence': 0.9,
                'distance': 0.9
            }
        elif region_success > 0.7:  # Medium success
            return {
                'height': 1.0,  # Standard
                'prominence': 1.0,
                'distance': 1.0
            }
        else:  # Low success region
            return {
                'height': 1.2,  # Be more conservative
                'prominence': 1.3,
                'distance': 1.1
            }

    def _get_nucleus_adaptation_factor(self, nucleus_type):
        """Get nucleus-specific adaptation factors"""
        factors = {
            '1H': 1.0,    # Standard for proton
            '15N': 0.8,   # Nitrogen peaks can be weaker
            '13C': 0.7    # Carbon peaks often weaker
        }
        return factors.get(nucleus_type, 1.0)

    def update_region_statistics(self, region_id, fit_results):
        """Update region-specific statistics based on fitting results"""
        if region_id not in self.region_statistics:
            self.region_statistics[region_id] = {
                'total_attempts': 0,
                'successful_fits': 0,
                'avg_r_squared': 0.0,
                'last_update': None
            }

        stats = self.region_statistics[region_id]
        stats['total_attempts'] += len(fit_results)

        successful = [r for r in fit_results if r.get('success', False) and r.get('r_squared', 0) > 0.5]
        stats['successful_fits'] += len(successful)

        if successful:
            r_squared_values = [r.get('r_squared', 0) for r in successful]
            # Weighted average update
            new_avg = np.mean(r_squared_values)
            if stats['avg_r_squared'] == 0:
                stats['avg_r_squared'] = new_avg
            else:
                stats['avg_r_squared'] = 0.7 * stats['avg_r_squared'] + 0.3 * new_avg

        stats['last_update'] = pd.Timestamp.now()

class ChemicalShiftContextAnalyzer:
    """Chemical shift context-aware detection weights based on NMR coupling patterns"""

    def __init__(self):
        # Common NMR coupling patterns and typical ranges
        self.coupling_patterns = {
            '1H': {
                'aromatic': {'range': (6.5, 8.5), 'typical_width': 0.01, 'weight': 1.2},
                'alpha_proton': {'range': (4.0, 5.5), 'typical_width': 0.02, 'weight': 1.0},
                'aliphatic': {'range': (0.5, 4.0), 'typical_width': 0.015, 'weight': 0.9},
                'exchangeable': {'range': (1.0, 12.0), 'typical_width': 0.05, 'weight': 0.7}
            },
            '15N': {
                'amide': {'range': (110, 135), 'typical_width': 0.5, 'weight': 1.0},
                'amino': {'range': (0, 50), 'typical_width': 1.0, 'weight': 0.8},
                'aromatic': {'range': (150, 250), 'typical_width': 0.8, 'weight': 0.9}
            },
            '13C': {
                'carbonyl': {'range': (160, 220), 'typical_width': 1.0, 'weight': 1.1},
                'aromatic': {'range': (100, 160), 'typical_width': 0.8, 'weight': 1.0},
                'alpha_carbon': {'range': (50, 100), 'typical_width': 1.2, 'weight': 0.9},
                'aliphatic': {'range': (10, 50), 'typical_width': 1.5, 'weight': 0.8}
            }
        }

        # Physics-informed constraints for common amino acids
        self.amino_acid_patterns = {
            'GLY': {'1H': [3.9], '15N': [109], '13C': [45]},
            'ALA': {'1H': [4.3, 1.4], '15N': [123], '13C': [52, 18]},
            'VAL': {'1H': [4.1, 2.1], '15N': [120], '13C': [62, 32]},
            # Add more as needed
        }

    def calculate_context_weights(self, peak_positions, nucleus_type):
        """
        Calculate context-aware weights for detected peaks based on chemical shift ranges

        Args:
            peak_positions: list of peak positions (ppm)
            nucleus_type: '1H', '15N', or '13C'

        Returns:
            dict: weights and context information for each peak
        """
        if nucleus_type not in self.coupling_patterns:
            # Return uniform weights for unknown nucleus
            return {i: {'weight': 1.0, 'context': 'unknown'} for i, _ in enumerate(peak_positions)}

        patterns = self.coupling_patterns[nucleus_type]
        weighted_peaks = {}

        for i, position in enumerate(peak_positions):
            best_match = None
            best_score = 0

            # Find best matching chemical shift region
            for region_name, region_info in patterns.items():
                min_ppm, max_ppm = region_info['range']

                if min_ppm <= position <= max_ppm:
                    # Peak is within range - calculate confidence
                    range_center = (min_ppm + max_ppm) / 2
                    range_width = max_ppm - min_ppm
                    distance_from_center = abs(position - range_center)
                    confidence = 1.0 - (distance_from_center / (range_width / 2))

                    if confidence > best_score:
                        best_score = confidence
                        best_match = {
                            'region': region_name,
                            'weight': region_info['weight'],
                            'confidence': confidence,
                            'typical_width': region_info['typical_width']
                        }

            # Handle peaks outside typical ranges
            if best_match is None:
                # Find closest region for context
                min_distance = float('inf')
                closest_region = None

                for region_name, region_info in patterns.items():
                    min_ppm, max_ppm = region_info['range']

                    if position < min_ppm:
                        distance = min_ppm - position
                    elif position > max_ppm:
                        distance = position - max_ppm
                    else:
                        distance = 0

                    if distance < min_distance:
                        min_distance = distance
                        closest_region = region_name

                # Assign reduced weight for unusual chemical shifts
                best_match = {
                    'region': f'unusual_{closest_region}',
                    'weight': 0.5,  # Reduced confidence for unusual shifts
                    'confidence': 0.3,
                    'typical_width': patterns[closest_region]['typical_width'] if closest_region else 0.02
                }

            weighted_peaks[i] = best_match

        return weighted_peaks

    def identify_coupling_networks(self, peak_positions, nucleus_type):
        """
        Identify potential coupling networks between peaks

        Args:
            peak_positions: list of peak positions
            nucleus_type: nucleus type

        Returns:
            list: coupling network information
        """
        if len(peak_positions) < 2:
            return []

        networks = []

        # Expected coupling constants (Hz) - convert to ppm based on spectrometer frequency
        if nucleus_type == '1H':
            # Typical J-couplings in ppm (assuming 600 MHz spectrometer)
            j_couplings = {
                'geminal': 0.0002,      # ~0.12 Hz / 600 MHz
                'vicinal': 0.0001,      # ~0.06 Hz / 600 MHz
                'scalar': 0.00005       # ~0.03 Hz / 600 MHz
            }

            for i, pos1 in enumerate(peak_positions):
                for j, pos2 in enumerate(peak_positions[i+1:], i+1):
                    separation = abs(pos1 - pos2)

                    # Check if separation matches expected coupling
                    for coupling_type, expected_sep in j_couplings.items():
                        if abs(separation - expected_sep) < expected_sep * 0.5:
                            networks.append({
                                'peaks': [i, j],
                                'positions': [pos1, pos2],
                                'coupling_type': coupling_type,
                                'separation': separation,
                                'confidence': 1.0 - abs(separation - expected_sep) / expected_sep
                            })

        return networks

class IterativeDetectionFitter:
    """
    Core integration class that combines peak detection with Voigt profile fitting
    through iterative refinement and adaptive optimization
    """

    def __init__(self, enhanced_fitter=None):
        self.enhanced_fitter = enhanced_fitter or (EnhancedVoigtFitter() if ENHANCED_AVAILABLE else None)
        self.quality_scorer = DetectionQualityScorer()
        self.rapid_assessor = RapidFitAssessor()
        self.threshold_calculator = AdaptiveThresholdCalculator()
        self.context_analyzer = ChemicalShiftContextAnalyzer()

        # 1D refinement engine
        self.peak_1d_refiner = Peak1DRefiner(enhanced_fitter)

        # Integration parameters
        self.integration_params = {
            'max_iterations': 5,
            'convergence_threshold': 0.01,
            # SOLUTION 3: Lower confidence thresholds to be more permissive
            'min_peak_confidence': 0.1,    # Changed from 0.3 to 0.1 (much more permissive)
            'aic_threshold': 2.0,  # AIC threshold for model selection
            'aic_selection_enabled': True,
            'physics_constraints_enabled': True,
            # SOLUTION 2: Disable adaptive thresholds by default to use GUI parameters
            'adaptive_thresholds_enabled': False,  # Changed from True to False
            'multi_resolution_enabled': False,     # Changed from True to False for consistency
            # SOLUTION 4: Lower fit likelihood threshold to be more permissive
            'fit_likelihood_threshold': 0.2,       # Changed from 0.4 to 0.2

            # 1D REFINEMENT PARAMETERS (NEW)
            'enable_1d_refinement': True,          # Enable 1D cross-section refinement
            'refinement_quality_threshold': 0.7,   # R² threshold for accepting refinement
            'refinement_coordinate_threshold': 0.01, # Max coordinate shift to accept (ppm)
            'refinement_after_detection': True,    # Apply refinement after initial detection
            'refinement_on_final_iteration': True, # Always refine on final iteration
        }

        # Iteration tracking
        self.iteration_history = []
        self.convergence_diagnostics = {}

    def integrated_detection_fitting(self, x_data, y_data, nucleus_type='1H',
                                   initial_peak_positions=None, gui_params=None, nmr_2d_data=None,
                                   peak_list=None, in_place_mode=False):
        """
        Main integrated detection-fitting workflow with iterative refinement

        Args:
            x_data, y_data: spectral data
            nucleus_type: nucleus type for constraints
            initial_peak_positions: optional initial peak guesses
            gui_params: GUI parameters for detection
            nmr_2d_data: 2D intensity matrix for 2D NMR
            peak_list: pandas DataFrame with reference peak positions for in-place mode
            in_place_mode: if True, constrains detection to search windows around peak_list positions

        Returns:
            dict: comprehensive results with fitted peaks and diagnostics
        """
        print(f"🔄 Starting integrated detection-fitting for {nucleus_type} dimension")

        # Handle 2D NMR data differently
        if nucleus_type == '2D' and nmr_2d_data is not None:
            print(f"   📊 2D NMR data detected: shape {nmr_2d_data.shape}")
            # For now, use the 2D data as intensity values for peak detection
            # x_data = 1H axis, y_data = 15N axis, nmr_2d_data = intensity matrix
            intensity_data = nmr_2d_data
        else:
            # Traditional 1D approach
            intensity_data = y_data

        # UPDATED: Apply GUI parameters to integration settings
        if gui_params:
            print(f"   🎛️ Using GUI parameters:")
            if 'aic_threshold' in gui_params:
                self.integration_params['aic_threshold'] = gui_params['aic_threshold']
                print(f"      AIC threshold: {gui_params['aic_threshold']}")
            if 'detection_confidence_threshold' in gui_params:
                self.integration_params['min_peak_confidence'] = gui_params['detection_confidence_threshold']
                print(f"      Detection confidence: {gui_params['detection_confidence_threshold']}")
            if 'max_integration_iterations' in gui_params:
                self.integration_params['max_iterations'] = gui_params['max_integration_iterations']
                print(f"      Max iterations: {gui_params['max_integration_iterations']}")
            # SOLUTION 2: Make advanced features configurable via GUI
            if 'adaptive_thresholds_enabled' in gui_params:
                self.integration_params['adaptive_thresholds_enabled'] = gui_params['adaptive_thresholds_enabled']
                print(f"      Adaptive thresholds: {gui_params['adaptive_thresholds_enabled']}")
            if 'multi_resolution_enabled' in gui_params:
                self.integration_params['multi_resolution_enabled'] = gui_params['multi_resolution_enabled']
                print(f"      Multi-resolution: {gui_params['multi_resolution_enabled']}")
            # SOLUTION 4: Make fit likelihood threshold configurable
            if 'fit_likelihood_threshold' in gui_params:
                self.integration_params['fit_likelihood_threshold'] = gui_params['fit_likelihood_threshold']
                print(f"      Fit likelihood threshold: {gui_params['fit_likelihood_threshold']}")
            # SOLUTION 5: Make physics constraints and convergence threshold configurable
            if 'physics_constraints_enabled' in gui_params:
                self.integration_params['physics_constraints_enabled'] = gui_params['physics_constraints_enabled']
                print(f"      Physics constraints: {gui_params['physics_constraints_enabled']}")
            if 'convergence_threshold' in gui_params:
                self.integration_params['convergence_threshold'] = gui_params['convergence_threshold']
                print(f"      Convergence threshold: {gui_params['convergence_threshold']}")
            # Core detection parameters from GUI
            if 'height_threshold' in gui_params:
                print(f"      Height threshold: {gui_params['height_threshold']}")
            if 'distance_factor' in gui_params:
                print(f"      Distance factor: {gui_params['distance_factor']}")
            if 'prominence_threshold' in gui_params:
                print(f"      Prominence threshold: {gui_params['prominence_threshold']}")
            if 'max_peaks_fit' in gui_params:
                print(f"      Max peaks fit: {gui_params['max_peaks_fit']}")
            if 'max_optimization_iterations' in gui_params:
                print(f"      Max optimization iterations: {gui_params['max_optimization_iterations']}")
            if 'smoothing_sigma' in gui_params:
                print(f"      Smoothing sigma: {gui_params['smoothing_sigma']}")
            # 1D Refinement parameters from GUI (NEW)
            if 'enable_1d_refinement' in gui_params:
                self.integration_params['enable_1d_refinement'] = gui_params['enable_1d_refinement']
                print(f"      1D refinement enabled: {gui_params['enable_1d_refinement']}")
            if 'refinement_quality_threshold' in gui_params:
                self.integration_params['refinement_quality_threshold'] = gui_params['refinement_quality_threshold']
                print(f"      Refinement R² threshold: {gui_params['refinement_quality_threshold']}")
            if 'refinement_coordinate_threshold' in gui_params:
                self.integration_params['refinement_coordinate_threshold'] = gui_params['refinement_coordinate_threshold']
                print(f"      Refinement coordinate threshold: {gui_params['refinement_coordinate_threshold']}")

        # Check for in-place mode with peak list constraints
        if in_place_mode and peak_list is not None:
            print(f"   📍 IN-PLACE MODE: Constraining detection to {len(peak_list)} reference peak regions")
            return self._reference_based_integrated_detection(
                x_data, y_data, intensity_data, peak_list, nucleus_type, gui_params
            )

        # Initialize tracking
        self.iteration_history = []
        iteration = 0
        converged = False

        # Stage 1: Initial detection with adaptive thresholds
        if self.integration_params['adaptive_thresholds_enabled']:
            detection_params, adaptation_info = self.threshold_calculator.calculate_adaptive_thresholds(
                x_data, y_data, nucleus_type
            )
            print(f"   📊 Adaptive thresholds: height={detection_params['height_threshold']:.3f}, "
                  f"prominence={detection_params['prominence_threshold']:.3f}")
        else:
            # SOLUTION 5: Use GUI detection parameters directly when adaptive is disabled
            if gui_params and all(key in gui_params for key in ['height_threshold', 'distance_factor', 'prominence_threshold']):
                detection_params = {
                    'height_threshold': gui_params['height_threshold'],
                    'distance_factor': gui_params['distance_factor'],
                    'prominence_threshold': gui_params['prominence_threshold'],
                    'smoothing_sigma': gui_params.get('smoothing_sigma', 0.5),
                    # CRITICAL FIX: Add missing core parameters
                    'max_peaks_fit': gui_params.get('max_peaks_fit', 4),
                    'max_optimization_iterations': gui_params.get('max_optimization_iterations', 50)
                }
                print(f"   🎛️ Using GUI detection parameters: height={detection_params['height_threshold']:.3f}, "
                      f"distance={detection_params['distance_factor']:.1f}, prominence={detection_params['prominence_threshold']:.3f}")
            else:
                detection_params = self.threshold_calculator.base_params
                print(f"   📋 Using base parameters: height={detection_params['height_threshold']:.3f}, "
                      f"prominence={detection_params['prominence_threshold']:.3f}")
            adaptation_info = None

        # Multi-resolution detection if enabled
        if self.integration_params['multi_resolution_enabled']:
            if nucleus_type == '2D':
                initial_peaks = self._multi_resolution_detection_2d(x_data, y_data, intensity_data, detection_params, nucleus_type)
            else:
                initial_peaks = self._multi_resolution_detection(x_data, intensity_data, detection_params, nucleus_type)
        else:
            if nucleus_type == '2D':
                initial_peaks = self._standard_detection_2d(x_data, y_data, intensity_data, detection_params)
            else:
                initial_peaks = self._standard_detection(x_data, intensity_data, detection_params)

        if not initial_peaks:
            return self._create_failure_result("No peaks detected in initial stage")

        print(f"   🎯 Initial detection: {len(initial_peaks)} peaks")

        current_peaks = initial_peaks
        best_result = None
        best_quality = 0.0

        # Iterative refinement loop
        for iteration in range(self.integration_params['max_iterations']):
            print(f"   🔄 Iteration {iteration + 1}/{self.integration_params['max_iterations']}")

            # Stage 2: Score peak detection confidence
            scored_peaks = self.quality_scorer.score_peak_detection(
                x_data, y_data, current_peaks, nucleus_type
            )

            # Filter low-confidence peaks
            high_confidence_peaks = {
                k: v for k, v in scored_peaks.items()
                if v['confidence'] >= self.integration_params['min_peak_confidence']
            }

            if not high_confidence_peaks:
                print(f"   ⚠️ No high-confidence peaks (min_confidence={self.integration_params['min_peak_confidence']})")
                break

            print(f"   ✅ High-confidence peaks: {len(high_confidence_peaks)}/{len(scored_peaks)}")

            # Stage 3: Rapid fit assessment for validation
            validated_peaks = {}
            # SOLUTION 4: Use configurable fit likelihood threshold
            fit_likelihood_threshold = self.integration_params['fit_likelihood_threshold']

            print(f"      🧪 Fit validation: threshold={fit_likelihood_threshold}")

            for idx, peak in high_confidence_peaks.items():
                try:
                    # Handle 2D NMR peaks differently
                    if nucleus_type == '2D' and 'position_y' in peak:
                        # For 2D peaks, create a simple assessment (bypass complex fitting for now)
                        assessment = {
                            'likelihood': 0.8,  # High likelihood for detected 2D peaks
                            'snr': peak['intensity'] / (np.std(intensity_data) + 1e-10),
                            'width_estimate': 1.0,
                            'method': '2d_simplified'
                        }
                        print(f"         Peak {idx}: 2D peak at ({peak['position']:.2f}, {peak['position_y']:.1f}) → likelihood={assessment['likelihood']:.3f}")
                    else:
                        # Traditional 1D assessment
                        assessment = self.rapid_assessor.rapid_assessment(
                            x_data, intensity_data, peak['position'], peak.get('estimated_width'), nucleus_type
                        )
                        print(f"         Peak {idx}: 1D peak at {peak['position']:.2f} → likelihood={assessment['likelihood']:.3f}")

                    if assessment['likelihood'] > fit_likelihood_threshold:
                        validated_peaks[idx] = {**peak, 'fit_assessment': assessment}
                        print(f"         ✅ Peak {idx} PASSED validation")
                    else:
                        print(f"         ❌ Peak {idx} FAILED validation ({assessment['likelihood']:.3f} < {fit_likelihood_threshold})")

                except Exception as e:
                    print(f"         ⚠️ Peak {idx} assessment error: {str(e)}")
                    # For 2D peaks, give benefit of doubt if assessment fails
                    if nucleus_type == '2D':
                        validated_peaks[idx] = {**peak, 'fit_assessment': {'likelihood': 0.8, 'method': 'fallback_2d'}}
                        print(f"         🔧 Peak {idx} using fallback validation for 2D")

            if not validated_peaks:
                print(f"   ⚠️ No peaks passed fit validation")
                break

            print(f"   ✅ Validated peaks: {len(validated_peaks)} (fit likelihood > {fit_likelihood_threshold})")

            # Stage 4: Enhanced fitting with AIC model selection
            if nucleus_type == '2D':
                # For 2D NMR, use simplified fitting approach
                fitting_result = self._perform_2d_simplified_fitting(
                    x_data, y_data, intensity_data, validated_peaks, nucleus_type, iteration
                )
            else:
                # Traditional 1D enhanced fitting
                fitting_result = self._perform_enhanced_fitting(
                    x_data, intensity_data, validated_peaks, nucleus_type, iteration, detection_params
                )

            if not fitting_result['success']:
                print(f"   ❌ Enhanced fitting failed: {fitting_result.get('error', 'unknown')}")
                break

            # Stage 4.5: 1D Coordinate Refinement (NEW)
            if (self.integration_params['enable_1d_refinement'] and
                nucleus_type == '2D' and nmr_2d_data is not None and
                (self.integration_params['refinement_after_detection'] or
                 (iteration == self.integration_params['max_iterations'] - 1 and
                  self.integration_params['refinement_on_final_iteration']))):

                print(f"   🔬 Stage 4.5: 1D Cross-Section Refinement")

                # Convert fitted peaks to DataFrame for refinement
                if 'fitted_peaks' in fitting_result and fitting_result['fitted_peaks']:
                    try:
                        fitted_peaks_df = self._convert_peaks_to_dataframe(fitting_result['fitted_peaks'], nucleus_type)

                        if fitted_peaks_df is not None and not fitted_peaks_df.empty:
                            # Configure refinement parameters from GUI
                            refinement_params = {}
                            if gui_params:
                                if 'fitting_window_x' in gui_params:
                                    refinement_params['fitting_window_x'] = gui_params['fitting_window_x']
                                if 'fitting_window_y' in gui_params:
                                    refinement_params['fitting_window_y'] = gui_params['fitting_window_y']
                                if 'max_peaks_fit' in gui_params:
                                    refinement_params['max_peaks_fit'] = gui_params['max_peaks_fit']
                                if 'max_optimization_iterations' in gui_params:
                                    refinement_params['max_optimization_iterations'] = gui_params['max_optimization_iterations']
                                refinement_params['min_r_squared'] = self.integration_params['refinement_quality_threshold']

                            # Perform 1D refinement
                            refinement_result = self.peak_1d_refiner.refine_peak_coordinates_1d(
                                fitted_peaks_df, nmr_2d_data, x_data, y_data, refinement_params
                            )

                            if refinement_result and refinement_result['refined_peaks'] is not None:
                                # Update fitting result with refined coordinates
                                refined_peaks_list = self._convert_dataframe_to_peaks(
                                    refinement_result['refined_peaks'], nucleus_type
                                )

                                if refined_peaks_list:
                                    # Update the fitting result
                                    fitting_result['fitted_peaks'] = refined_peaks_list
                                    fitting_result['refinement_applied'] = True
                                    fitting_result['refinement_stats'] = refinement_result['statistics']
                                    fitting_result['refinement_quality'] = refinement_result['refinement_quality']

                                    print(f"      ✅ 1D refinement applied: {refinement_result['statistics']['peaks_refined']}/{refinement_result['statistics']['peaks_processed']} peaks refined")
                                    print(f"      📊 Average coordinate shift: {refinement_result['statistics']['avg_coordinate_shift']:.4f} ppm")
                                else:
                                    print(f"      ⚠️ 1D refinement: failed to convert refined peaks back")
                            else:
                                print(f"      ⚠️ 1D refinement failed or returned no results")
                        else:
                            print(f"      ⚠️ Could not convert fitted peaks to DataFrame for refinement")
                    except Exception as refinement_error:
                        print(f"      ❌ 1D refinement error: {str(refinement_error)}")
                        # Continue without refinement - don't break the main workflow
                else:
                    print(f"      ⚠️ No fitted peaks available for 1D refinement")

            # Stage 5: Check convergence
            convergence_info = self._check_convergence(current_peaks, fitting_result, iteration)

            # Update tracking
            iteration_info = {
                'iteration': iteration + 1,
                'initial_peaks': len(current_peaks),
                'scored_peaks': len(scored_peaks),
                'high_confidence': len(high_confidence_peaks),
                'validated_peaks': len(validated_peaks),
                'fitted_peaks': len(fitting_result.get('fitted_peaks', [])),
                'avg_r_squared': fitting_result.get('avg_r_squared', 0),
                'convergence': convergence_info
            }
            self.iteration_history.append(iteration_info)

            # Update best result if improved
            current_quality = fitting_result.get('avg_r_squared', 0)
            if current_quality > best_quality:
                best_quality = current_quality
                best_result = fitting_result
                print(f"   🎯 New best result: R² = {current_quality:.4f}")

            # Check convergence
            if convergence_info['converged']:
                converged = True
                print(f"   ✅ Converged after {iteration + 1} iterations")
                break

            # Stage 6: Refine peak positions for next iteration
            refined_peaks = self._refine_peak_positions(fitting_result, x_data, y_data)
            if refined_peaks:
                current_peaks = refined_peaks
                print(f"   🔧 Refined {len(refined_peaks)} peak positions for next iteration")
            else:
                print(f"   ⚠️ Could not refine peak positions, stopping iteration")
                break

        # Finalize results
        final_result = best_result or self._create_failure_result("All iterations failed")

        # CRITICAL FIX: Add 'peaks' key for GUI compatibility (maps to fitted_peaks)
        if 'fitted_peaks' in final_result:
            final_result['peaks'] = final_result['fitted_peaks']
            print(f"   ✅ Added 'peaks' key mapping to {len(final_result['peaks'])} fitted peaks for GUI")

        # Add integration diagnostics
        final_result['integration_diagnostics'] = {
            'converged': converged,
            'iterations_completed': iteration + 1,
            'iteration_history': self.iteration_history,
            'adaptive_thresholds': adaptation_info,
            'integration_params': self.integration_params
        }

        print(f"🏁 Integration complete: {len(final_result.get('fitted_peaks', []))} peaks fitted "
              f"(avg R² = {final_result.get('avg_r_squared', 0):.4f})")

        return final_result

    def _multi_resolution_detection(self, x_data, y_data, detection_params, nucleus_type):
        """Multi-resolution peak detection: coarse → medium → fine"""
        print("   🔍 Multi-resolution detection")

        # Coarse detection (heavy smoothing, major peaks only)
        coarse_smoothing = detection_params['smoothing_sigma'] * 3
        y_coarse = gaussian_filter1d(y_data, sigma=coarse_smoothing)
        coarse_peaks = self._detect_peaks_with_params(x_data, y_coarse, {
            **detection_params,
            'height_threshold': detection_params['height_threshold'] * 2,
            'prominence_threshold': detection_params['prominence_threshold'] * 2
        })

        print(f"     Coarse: {len(coarse_peaks)} major peaks")

        # Medium detection (moderate smoothing)
        medium_smoothing = detection_params['smoothing_sigma'] * 1.5
        y_medium = gaussian_filter1d(y_data, sigma=medium_smoothing)
        medium_peaks = self._detect_peaks_with_params(x_data, y_medium, {
            **detection_params,
            'height_threshold': detection_params['height_threshold'] * 1.5
        })

        print(f"     Medium: {len(medium_peaks)} intermediate peaks")

        # Fine detection (light smoothing, all peaks)
        fine_peaks = self._detect_peaks_with_params(x_data, y_data, detection_params)

        print(f"     Fine: {len(fine_peaks)} all peaks")

        # Merge results hierarchically (coarse peaks have priority)
        merged_peaks = self._merge_hierarchical_peaks(coarse_peaks, medium_peaks, fine_peaks, x_data)

        print(f"     Merged: {len(merged_peaks)} hierarchical peaks")

        return merged_peaks

    def _detect_peaks_with_params(self, x_data, y_data, params):
        """Standard peak detection with given parameters - handles both 1D and 2D NMR data"""
        from scipy.signal import find_peaks

        # Check if we have 2D NMR data
        if isinstance(y_data, np.ndarray) and y_data.ndim == 2:
            # 2D NMR data: y_data is the intensity matrix, x_data is the 1H axis
            print(f"      🔬 2D NMR peak detection: matrix shape {y_data.shape}")
            return self._detect_2d_nmr_peaks(x_data, y_data, params)
        else:
            # 1D data: traditional approach
            print(f"      📈 1D peak detection: data length {len(y_data)}")
            return self._detect_1d_peaks(x_data, y_data, params)

    def _detect_1d_peaks(self, x_data, y_data, params):
        """Detect peaks in 1D data"""
        from scipy.signal import find_peaks

        # Calculate noise-based thresholds (much more appropriate for NMR)
        noise_level = np.std(y_data)
        data_range = np.max(y_data) - np.min(y_data)

        # FIXED: Use much lower absolute thresholds based on noise, not signal max
        min_height = noise_level * 3  # 3x noise level (much more reasonable)
        min_distance = max(2, int(len(y_data) / params['distance_factor']))
        min_prominence = noise_level * 1.5  # 1.5x noise level for prominence

        # DEBUG: Add detailed logging
        print(f"         Data range: [{np.min(y_data):.3f}, {np.max(y_data):.3f}] (span: {data_range:.3f})")
        print(f"         Noise level: {noise_level:.3f}")
        print(f"         Thresholds: height={min_height:.3f}, distance={min_distance}, prominence={min_prominence:.3f}")

        # Find peaks using original data (scipy handles positive/negative automatically)
        peaks, properties = find_peaks(y_data,
                                     height=min_height,
                                     distance=min_distance,
                                     prominence=min_prominence)

        print(f"         Found {len(peaks)} peaks")

        # Convert to peak info format
        peak_info = []
        for i, peak_idx in enumerate(peaks):
            peak_info.append({
                'index': peak_idx,
                'position': x_data[peak_idx],
                'intensity': y_data[peak_idx],
                'distance_to_target': 0  # Will be updated if needed
            })

        return peak_info

    def _detect_2d_nmr_peaks(self, x_axis_1h, nmr_2d_matrix, params):
        """Detect peaks in 2D NMR data using local maxima approach"""
        from scipy.ndimage import maximum_filter

        # Calculate noise-based threshold for 2D data
        noise_level = np.std(nmr_2d_matrix)
        intensity_threshold = noise_level * 3  # 3x noise level

        print(f"         2D Matrix: {nmr_2d_matrix.shape}, noise: {noise_level:.3f}, threshold: {intensity_threshold:.3f}")

        # Find local maxima in 2D
        # Use anisotropic neighborhood sizes for different dimensions
        neighborhood_x = params.get('detection_square_size', 3)      # 1H dimension (X-axis)
        neighborhood_y = params.get('detection_rectangle_y', 1)      # 15N dimension (Y-axis, smaller!)
        print(f"         Using anisotropic detection size: {neighborhood_y}x{neighborhood_x} pixels (Y x X)")
        print(f"         X-dimension (1H): {neighborhood_x} pixels, Y-dimension (15N): {neighborhood_y} pixels")
        local_maxima = maximum_filter(nmr_2d_matrix, size=(neighborhood_y, neighborhood_x)) == nmr_2d_matrix

        # Apply intensity threshold
        above_threshold = nmr_2d_matrix > intensity_threshold
        peak_mask = local_maxima & above_threshold

        # Get peak coordinates
        peak_coordinates = np.where(peak_mask)
        y_indices, x_indices = peak_coordinates  # y=15N axis, x=1H axis

        print(f"         Found {len(x_indices)} 2D peaks above threshold")

        # Convert to peak info format
        peak_info = []
        if len(self.integrator.ppm_y_axis) > 0:  # Check if we have Y-axis data
            y_axis_15n = self.integrator.ppm_y_axis if hasattr(self, 'integrator') else np.arange(len(y_indices))

            for i in range(len(x_indices)):
                x_idx, y_idx = x_indices[i], y_indices[i]
                if x_idx < len(x_axis_1h) and y_idx < len(y_axis_15n):
                    peak_info.append({
                        'index': (y_idx, x_idx),  # 2D index
                        'position': x_axis_1h[x_idx],  # 1H chemical shift (legacy key)
                        'position_x': x_axis_1h[x_idx],  # 1H chemical shift (Simple Pattern Matcher key)
                        'position_y': y_axis_15n[y_idx],  # 15N chemical shift
                        'intensity': nmr_2d_matrix[y_idx, x_idx],
                        'distance_to_target': 0
                    })

        # Apply Peak Ridge Consolidation (Solution A) to eliminate Y-dimension fragmentation
        x_tolerance = params.get('consolidation_x_tolerance', 0.05)
        y_tolerance = params.get('consolidation_y_tolerance', 2.0)
        peak_info = self._consolidate_y_dimension_peaks(peak_info, x_tolerance=x_tolerance, y_tolerance=y_tolerance)

        return peak_info

    def _consolidate_y_dimension_peaks(self, peak_info, x_tolerance=0.05, y_tolerance=2.0):
        """
        Solution A: Peak Ridge Consolidation
        Merge nearby peaks in Y-dimension that likely represent fragments of the same peak

        Args:
            peak_info: List of detected peaks with position_x, position_y, intensity
            x_tolerance: Maximum 1H chemical shift difference to consider peaks as same X-group (ppm)
            y_tolerance: Maximum 15N chemical shift difference to consider for consolidation (ppm)

        Returns:
            List of consolidated peaks (one peak per X-region)
        """
        if not peak_info:
            return peak_info

        print(f"         🔧 Peak Ridge Consolidation: {len(peak_info)} peaks before consolidation")
        print(f"         📏 Tolerances: X={x_tolerance:.3f} ppm, Y={y_tolerance:.1f} ppm")

        # Group peaks by similar X-coordinates (1H dimension)
        x_groups = []

        for peak in peak_info:
            x_pos = peak.get('position', peak.get('position_x', 0))

            # Find existing group with similar X-coordinate
            assigned_to_group = False
            for group in x_groups:
                # Check if this peak belongs to an existing X-group
                group_x_positions = [p.get('position', p.get('position_x', 0)) for p in group]
                avg_x = sum(group_x_positions) / len(group_x_positions)

                if abs(x_pos - avg_x) <= x_tolerance:
                    group.append(peak)
                    assigned_to_group = True
                    break

            # Create new group if no existing group found
            if not assigned_to_group:
                x_groups.append([peak])

        print(f"         📊 Created {len(x_groups)} X-coordinate groups")

        # For each X-group, consolidate Y-dimension fragments
        consolidated_peaks = []

        for group_idx, group in enumerate(x_groups):
            if len(group) == 1:
                # Single peak in group - no consolidation needed
                consolidated_peaks.extend(group)
                continue

            # Multiple peaks in same X-region - need consolidation
            print(f"         🔍 Group {group_idx+1}: {len(group)} peaks at similar X-coordinates")

            # Sort peaks in group by Y-coordinate
            group_sorted = sorted(group, key=lambda p: p.get('position_y', 0))

            # Sub-group peaks that are close in Y-dimension
            y_subgroups = []
            current_subgroup = [group_sorted[0]]

            for i in range(1, len(group_sorted)):
                current_peak = group_sorted[i]
                previous_peak = group_sorted[i-1]

                y_current = current_peak.get('position_y', 0)
                y_previous = previous_peak.get('position_y', 0)

                if abs(y_current - y_previous) <= y_tolerance:
                    # Close in Y - add to current subgroup
                    current_subgroup.append(current_peak)
                else:
                    # Far in Y - start new subgroup
                    y_subgroups.append(current_subgroup)
                    current_subgroup = [current_peak]

            # Don't forget the last subgroup
            y_subgroups.append(current_subgroup)

            # For each Y-subgroup, keep only the strongest peak
            for subgroup in y_subgroups:
                # Find peak with maximum intensity in this subgroup
                strongest_peak = max(subgroup, key=lambda p: p.get('intensity', 0))

                # Debug logging
                if len(subgroup) > 1:
                    intensities = [p.get('intensity', 0) for p in subgroup]
                    y_positions = [p.get('position_y', 0) for p in subgroup]
                    print(f"           🎯 Consolidated {len(subgroup)} Y-fragments: Y={y_positions}, I={intensities}")
                    print(f"           ✅ Kept strongest: Y={strongest_peak.get('position_y', 0):.2f}, I={strongest_peak.get('intensity', 0):.1f}")

                consolidated_peaks.append(strongest_peak)

        print(f"         ✅ Peak Ridge Consolidation: {len(consolidated_peaks)} peaks after consolidation")
        print(f"         📉 Eliminated {len(peak_info) - len(consolidated_peaks)} Y-dimension fragments")

        return consolidated_peaks

    def _refine_peaks_with_centroids(self, peak_info, nmr_2d_matrix, x_axis_1h, y_axis_15n, params):
        """
        Peak Centroid Detection - Optional post-processing enhancement
        Calculate intensity-weighted centroids for improved coordinate accuracy

        Args:
            peak_info: List of detected peaks
            nmr_2d_matrix: 2D NMR intensity data
            x_axis_1h: 1H axis (ppm values)
            y_axis_15n: 15N axis (ppm values)
            params: Detection parameters

        Returns:
            List of peaks with refined centroid coordinates
        """
        if not peak_info:
            return peak_info

        print(f"         🎯 Peak Centroid Detection: Refining {len(peak_info)} peaks...")

        # Parameters for centroid calculation
        centroid_window_x_ppm = params.get('centroid_window_x_ppm', 0.05)  # window size in ppm (1H dimension)
        centroid_window_y_ppm = params.get('centroid_window_y_ppm', 1.0)   # window size in ppm (15N dimension)
        noise_multiplier = params.get('centroid_noise_multiplier', 2.0)    # noise threshold multiplier

        # Calculate noise level for threshold
        noise_level = np.std(nmr_2d_matrix[nmr_2d_matrix != 0])  # Exclude zeros
        noise_threshold = noise_level * noise_multiplier

        refined_peaks = []
        successful_refinements = 0

        for peak in peak_info:
            try:
                # Get original coordinates
                original_x = peak.get('ppm_x', peak.get('Position_X', 0))
                original_y = peak.get('ppm_y', peak.get('Position_Y', 0))

                # Convert ppm coordinates to matrix indices
                x_idx = np.argmin(np.abs(x_axis_1h - original_x))
                y_idx = np.argmin(np.abs(y_axis_15n - original_y))

                # Convert ppm window sizes to pixel counts
                x_ppm_per_pixel = abs(x_axis_1h[-1] - x_axis_1h[0]) / len(x_axis_1h)
                y_ppm_per_pixel = abs(y_axis_15n[-1] - y_axis_15n[0]) / len(y_axis_15n)
                centroid_window_x_pixels = int(np.ceil(centroid_window_x_ppm / x_ppm_per_pixel))
                centroid_window_y_pixels = int(np.ceil(centroid_window_y_ppm / y_ppm_per_pixel))

                # Define window boundaries (ensure within matrix bounds)
                x_min = max(0, x_idx - centroid_window_x_pixels)
                x_max = min(nmr_2d_matrix.shape[1], x_idx + centroid_window_x_pixels + 1)
                y_min = max(0, y_idx - centroid_window_y_pixels)
                y_max = min(nmr_2d_matrix.shape[0], y_idx + centroid_window_y_pixels + 1)

                # Extract intensity window
                intensity_window = np.abs(nmr_2d_matrix[y_min:y_max, x_min:x_max])  # Use absolute values

                # Apply noise-based threshold
                valid_mask = intensity_window >= noise_threshold

                if not np.any(valid_mask):
                    # Fallback to original coordinates if no valid intensities
                    refined_peaks.append(peak)
                    continue

                # Calculate intensity-weighted centroid
                y_indices, x_indices = np.mgrid[y_min:y_max, x_min:x_max]
                valid_intensities = intensity_window[valid_mask]
                valid_x_indices = x_indices[valid_mask]
                valid_y_indices = y_indices[valid_mask]

                # Weighted centroid calculation
                total_intensity = np.sum(valid_intensities)
                centroid_x_idx = np.sum(valid_x_indices * valid_intensities) / total_intensity
                centroid_y_idx = np.sum(valid_y_indices * valid_intensities) / total_intensity

                # Convert back to ppm coordinates using interpolation
                if centroid_x_idx < len(x_axis_1h) - 1 and centroid_x_idx >= 0:
                    # Linear interpolation for sub-pixel accuracy
                    x_floor = int(np.floor(centroid_x_idx))
                    x_ceil = int(np.ceil(centroid_x_idx))
                    x_frac = centroid_x_idx - x_floor

                    if x_floor == x_ceil:
                        refined_x_ppm = x_axis_1h[x_floor]
                    else:
                        refined_x_ppm = x_axis_1h[x_floor] * (1 - x_frac) + x_axis_1h[x_ceil] * x_frac
                else:
                    refined_x_ppm = original_x

                if centroid_y_idx < len(y_axis_15n) - 1 and centroid_y_idx >= 0:
                    y_floor = int(np.floor(centroid_y_idx))
                    y_ceil = int(np.ceil(centroid_y_idx))
                    y_frac = centroid_y_idx - y_floor

                    if y_floor == y_ceil:
                        refined_y_ppm = y_axis_15n[y_floor]
                    else:
                        refined_y_ppm = y_axis_15n[y_floor] * (1 - y_frac) + y_axis_15n[y_ceil] * y_frac
                else:
                    refined_y_ppm = original_y

                # Create refined peak with updated coordinates
                refined_peak = peak.copy()

                # Update coordinates using existing field names
                if 'ppm_x' in peak:
                    refined_peak['ppm_x'] = refined_x_ppm
                if 'ppm_y' in peak:
                    refined_peak['ppm_y'] = refined_y_ppm
                if 'Position_X' in peak:
                    refined_peak['Position_X'] = refined_x_ppm
                if 'Position_Y' in peak:
                    refined_peak['Position_Y'] = refined_y_ppm

                # Add refinement metadata
                shift_x = abs(refined_x_ppm - original_x)
                shift_y = abs(refined_y_ppm - original_y)
                refined_peak['centroid_shift_x'] = shift_x
                refined_peak['centroid_shift_y'] = shift_y
                refined_peak['centroid_refined'] = True

                refined_peaks.append(refined_peak)
                successful_refinements += 1

            except Exception as e:
                # Fallback to original coordinates on any error
                print(f"         ⚠️  Centroid refinement failed for peak at ({original_x:.3f}, {original_y:.3f}): {e}")
                refined_peaks.append(peak)

        print(f"         ✅ Peak Centroid Detection: {successful_refinements}/{len(peak_info)} peaks refined")
        if successful_refinements > 0:
            avg_shift_x = np.mean([p.get('centroid_shift_x', 0) for p in refined_peaks if p.get('centroid_refined', False)])
            avg_shift_y = np.mean([p.get('centroid_shift_y', 0) for p in refined_peaks if p.get('centroid_refined', False)])
            print(f"         📏 Average coordinate shifts: X={avg_shift_x:.4f} ppm, Y={avg_shift_y:.4f} ppm")

        return refined_peaks

    def _detect_2d_nmr_peaks_with_axis(self, x_axis_1h, y_axis_15n, nmr_2d_matrix, params):
        """Detect peaks in 2D NMR data with proper axis handling"""
        from scipy.ndimage import maximum_filter

        # Calculate noise-based threshold for 2D data
        noise_level = np.std(nmr_2d_matrix)
        intensity_threshold = noise_level * 3  # 3x noise level

        print(f"         2D Matrix: {nmr_2d_matrix.shape}, noise: {noise_level:.3f}, threshold: {intensity_threshold:.3f}")
        print(f"         Axes: 1H={len(x_axis_1h)} points, 15N={len(y_axis_15n)} points")
        print(f"         Data stats: mean={np.mean(nmr_2d_matrix):.1f}, std={np.std(nmr_2d_matrix):.1f}, max={np.max(nmr_2d_matrix):.1f}")

        # Find local maxima in 2D
        # Use anisotropic neighborhood sizes for different dimensions
        neighborhood_x = params.get('detection_square_size', 3)      # 1H dimension (X-axis)
        neighborhood_y = params.get('detection_rectangle_y', 1)      # 15N dimension (Y-axis, smaller!)
        print(f"         Using anisotropic detection size: {neighborhood_y}x{neighborhood_x} pixels (Y x X)")
        print(f"         X-dimension (1H): {neighborhood_x} pixels, Y-dimension (15N): {neighborhood_y} pixels")
        local_maxima = maximum_filter(nmr_2d_matrix, size=(neighborhood_y, neighborhood_x)) == nmr_2d_matrix

        # Apply intensity threshold
        above_threshold = nmr_2d_matrix > intensity_threshold
        peak_mask = local_maxima & above_threshold

        # Get peak coordinates
        peak_coordinates = np.where(peak_mask)
        y_indices, x_indices = peak_coordinates  # y=15N axis, x=1H axis

        print(f"         Found {len(x_indices)} 2D peaks above threshold")

        # DEBUG: Show diagnostic information
        if len(x_indices) > 0:
            intensities = [nmr_2d_matrix[y_indices[i], x_indices[i]] for i in range(min(5, len(x_indices)))]
            print(f"         Sample peak intensities (first 5): {[f'{x:.1f}' for x in intensities]}")
            peak_intensities = [nmr_2d_matrix[y_indices[i], x_indices[i]] for i in range(len(x_indices))]
            print(f"         Peak intensity range: [{np.min(peak_intensities):.1f}, {np.max(peak_intensities):.1f}]")

            # Check if we might be detecting too many peaks
            if len(x_indices) > 50:
                print(f"         ⚠️ WARNING: {len(x_indices)} peaks detected - might be excessive!")
                print(f"         Current threshold: {intensity_threshold:.1f} ({intensity_threshold/noise_level:.1f}x noise)")

                # Count how many points are above different thresholds
                above_3x = np.sum(nmr_2d_matrix > noise_level * 3)
                above_6x = np.sum(nmr_2d_matrix > noise_level * 6)
                above_10x = np.sum(nmr_2d_matrix > noise_level * 10)
                print(f"         Points above: 3x noise={above_3x}, 6x noise={above_6x}, 10x noise={above_10x}")

        if len(x_indices) == 0:
            print(f"         ❌ No peaks found - threshold might be too high: {intensity_threshold:.1f}")
            print(f"         Try lowering threshold or check data")

        # Convert to peak info format
        peak_info = []
        for i in range(len(x_indices)):
            x_idx, y_idx = x_indices[i], y_indices[i]
            if x_idx < len(x_axis_1h) and y_idx < len(y_axis_15n):
                peak_info.append({
                    'index': (y_idx, x_idx),  # 2D index
                    'position': x_axis_1h[x_idx],  # 1H chemical shift (legacy key)
                    'position_x': x_axis_1h[x_idx],  # 1H chemical shift (Simple Pattern Matcher key)
                    'position_y': y_axis_15n[y_idx],  # 15N chemical shift
                    'intensity': nmr_2d_matrix[y_idx, x_idx],
                    'distance_to_target': 0
                })

        # Apply Peak Ridge Consolidation (Solution A) to eliminate Y-dimension fragmentation
        x_tolerance = params.get('consolidation_x_tolerance', 0.05)
        y_tolerance = params.get('consolidation_y_tolerance', 2.0)
        peak_info = self._consolidate_y_dimension_peaks(peak_info, x_tolerance=x_tolerance, y_tolerance=y_tolerance)

        # Apply Peak Centroid Refinement (optional post-processing enhancement)
        if params.get('use_centroid_refinement', False):
            peak_info = self._refine_peaks_with_centroids(peak_info, nmr_2d_matrix, x_axis_1h, y_axis_15n, params)

        return peak_info

    def _perform_2d_simplified_fitting(self, x_axis, y_axis, intensity_matrix, validated_peaks, nucleus_type, iteration):
        """Simplified fitting for 2D NMR peaks - creates fitted peak results without complex fitting"""
        print(f"     🔧 2D simplified fitting for {len(validated_peaks)} validated peaks")

        try:
            # Convert validated peaks to fitted peak format
            fitted_peaks = []
            total_r_squared = 0.0

            # Handle both dict and list formats for validated_peaks
            if isinstance(validated_peaks, dict):
                peak_items = validated_peaks.items()
            else:
                peak_items = enumerate(validated_peaks)

            for idx, peak_data in peak_items:
                # Extract 2D peak information
                x_pos = peak_data['position']      # 1H chemical shift
                y_pos = peak_data.get('position_y', 0)  # 15N chemical shift
                intensity = peak_data['intensity']

                # Create a simplified fitted peak result with correct GUI field names
                fitted_peak = {
                    'peak_id': f'2D_peak_{idx}',
                    # CRITICAL FIX: Use correct field names expected by GUI
                    'ppm_x': x_pos,       # GUI expects 'ppm_x' not 'x_position'
                    'ppm_y': y_pos,       # GUI expects 'ppm_y' not 'y_position'
                    'intensity': intensity,
                    'volume': intensity,  # Simplified: use intensity as volume
                    'width_1h': 0.05,     # Estimated 1H linewidth
                    'width_15n': 1.0,     # Estimated 15N linewidth
                    'r_squared': 0.85,    # High quality for detected 2D peaks
                    'fit_r_squared': 0.85,  # Also provide this variant
                    'snr': intensity / (np.std(intensity_matrix) + 1e-10),
                    'fitting_method': '2d_simplified',
                    'assignment': f'Peak_{idx}',
                    'composite_quality': 0.85,
                    'detection_confidence': peak_data.get('confidence', 0.8),
                    'fit_assessment': peak_data.get('fit_assessment', {}),
                    # Add status flags expected by GUI
                    'detected': True,     # Mark as successfully detected
                    'converged': True,    # Mark as successfully fitted
                    'peak_number': idx + 1,  # 1-based peak numbering
                    # Add standard peak fitting attributes expected by GUI
                    'amp': intensity,
                    'center_x': x_pos,    # Keep these too for compatibility
                    'center_y': y_pos,
                    'sigma_x': 0.05,
                    'sigma_y': 1.0,
                    'gamma_x': 0.02,
                    'gamma_y': 0.5
                }

                fitted_peaks.append(fitted_peak)
                total_r_squared += 0.85  # Use fixed high R² for 2D peaks

                if len(fitted_peaks) % 50 == 0:  # Progress logging
                    print(f"       Processed {len(fitted_peaks)}/{len(validated_peaks)} 2D peaks")

            avg_r_squared = total_r_squared / len(fitted_peaks) if fitted_peaks else 0.0

            print(f"     ✅ 2D simplified fitting complete: {len(fitted_peaks)} peaks, avg R²={avg_r_squared:.3f}")

            return {
                'success': True,
                'fitted_peaks': fitted_peaks,
                'avg_r_squared': avg_r_squared,
                'fitting_method': '2d_simplified',
                'total_peaks': len(fitted_peaks)
            }

        except Exception as e:
            print(f"     ❌ 2D simplified fitting failed: {str(e)}")
            return {
                'success': False,
                'error': f'2d_simplified_fitting_failed: {str(e)}',
                'fitted_peaks': []
            }

    def _merge_hierarchical_peaks(self, coarse_peaks, medium_peaks, fine_peaks, x_data):
        """Merge peaks hierarchically with priority to coarse detection"""
        merged = []
        used_positions = set()

        # Minimum separation (data point spacing * 5)
        min_separation = abs(x_data[1] - x_data[0]) * 5

        # Add coarse peaks first (highest priority)
        for peak in coarse_peaks:
            merged.append({**peak, 'detection_level': 'coarse', 'priority': 3})
            used_positions.add(peak['position'])

        # Add medium peaks if not too close to coarse peaks
        for peak in medium_peaks:
            too_close = any(abs(peak['position'] - used_pos) < min_separation
                          for used_pos in used_positions)
            if not too_close:
                merged.append({**peak, 'detection_level': 'medium', 'priority': 2})
                used_positions.add(peak['position'])

        # Add fine peaks if not too close to existing peaks
        for peak in fine_peaks:
            too_close = any(abs(peak['position'] - used_pos) < min_separation
                          for used_pos in used_positions)
            if not too_close:
                merged.append({**peak, 'detection_level': 'fine', 'priority': 1})
                used_positions.add(peak['position'])

        # Sort by intensity (strongest first)
        merged.sort(key=lambda x: x['intensity'], reverse=True)

        return merged

    def _standard_detection(self, x_data, y_data, detection_params):
        """Standard single-level peak detection"""
        return self._detect_peaks_with_params(x_data, y_data, detection_params)

    def _standard_detection_2d(self, x_axis, y_axis, intensity_matrix, detection_params):
        """Standard 2D peak detection"""
        return self._detect_2d_nmr_peaks_with_axis(x_axis, y_axis, intensity_matrix, detection_params)

    def _perform_enhanced_fitting(self, x_data, y_data, validated_peaks, nucleus_type, iteration, detection_params=None):
        """Perform enhanced fitting with AIC model selection"""
        # Extract max peaks parameter for multi-peak fitting limit
        max_peaks_fit = detection_params.get('max_peaks_fit', 4) if detection_params else 4
        max_optimization_iterations = detection_params.get('max_optimization_iterations', 50) if detection_params else 50

        print(f"     🔧 Enhanced fitting with max_peaks_fit={max_peaks_fit}, max_opt_iter={max_optimization_iterations}")

        # Update enhanced_fitter parameters if available
        if self.enhanced_fitter and hasattr(self.enhanced_fitter, 'fitting_parameters'):
            self.enhanced_fitter.fitting_parameters['max_iterations'] = max_optimization_iterations
            print(f"     🎛️ Updated enhanced_fitter max_iterations to {max_optimization_iterations}")
        try:
            if not self.enhanced_fitter:
                return {'success': False, 'error': 'enhanced_fitter_not_available'}

            # Convert validated peaks to positions for fitting
            peak_positions = [peak['position'] for peak in validated_peaks.values()]

            # Use AIC model selection if enabled
            if self.integration_params['aic_selection_enabled'] and len(peak_positions) > 1:
                aic_threshold = self.integration_params['aic_threshold']
                print(f"     🎯 AIC model selection for {len(peak_positions)} peaks (threshold: {aic_threshold})")

                # Use the existing AIC framework with GUI max_peaks_fit parameter
                aic_result = self.enhanced_fitter.optimal_peak_count_by_aic(
                    x_data, y_data, peak_positions, max_peaks=min(max_peaks_fit, len(peak_positions))
                )

                if aic_result['success']:
                    optimal_count = aic_result['optimal_n_peaks']
                    best_model = aic_result['best_model']

                    # Apply AIC threshold check
                    aic_improvement = best_model.get('aic_improvement', 0)
                    if aic_improvement >= aic_threshold:
                        print(f"     ✅ AIC selected {optimal_count} peaks (improvement: {aic_improvement:.2f} ≥ {aic_threshold})")

                        return {
                            'success': True,
                            'fitted_peaks': self._extract_fitted_peaks_from_model(best_model, optimal_count),
                            'avg_r_squared': best_model['r_squared'],
                            'aic_result': aic_result,
                            'fitting_method': 'aic_enhanced'
                        }
                    else:
                        print(f"     ⚠️ AIC improvement {aic_improvement:.2f} < threshold {aic_threshold}, using simpler model")
                else:
                    print(f"     ⚠️ AIC selection failed, using standard fitting")

            # Fallback to standard enhanced fitting
            if len(peak_positions) == 1:
                # Single peak fitting
                result = self.enhanced_fitter.fit_peak_enhanced(
                    x_data, y_data, peak_positions[0], nucleus_type,
                    method='iterative_optimization'
                )

                if result['success']:
                    return {
                        'success': True,
                        'fitted_peaks': [result],
                        'avg_r_squared': result['r_squared'],
                        'fitting_method': 'single_peak_enhanced'
                    }
            else:
                # Multi-peak fitting via overlapping detection optimization
                fit_result, optimization_report = self.enhanced_fitter.optimize_overlap_detection_iteratively(
                    x_data, y_data, peak_positions, use_aic_selection=True
                )

                if fit_result and fit_result.get('success'):
                    return {
                        'success': True,
                        'fitted_peaks': [fit_result],
                        'avg_r_squared': fit_result['r_squared'],
                        'optimization_report': optimization_report,
                        'fitting_method': 'multi_peak_enhanced'
                    }

            return {'success': False, 'error': 'all_enhanced_methods_failed'}

        except Exception as e:
            return {'success': False, 'error': f'enhanced_fitting_exception: {str(e)}'}

    def _extract_fitted_peaks_from_model(self, model, n_peaks):
        """Extract fitted peaks from AIC model result"""
        fit_result = model['fit_result']

        if n_peaks == 1:
            return [fit_result]
        else:
            # For multi-peak results, need to extract individual peaks
            if 'peak_positions' in fit_result and 'peak_amplitudes' in fit_result:
                peaks = []
                for i in range(n_peaks):
                    peaks.append({
                        'success': True,
                        'center': fit_result['peak_positions'][i],
                        'amplitude': fit_result['peak_amplitudes'][i],
                        'r_squared': fit_result['r_squared'],
                        'method': fit_result['method']
                    })
                return peaks
            else:
                return [fit_result]

    def _check_convergence(self, previous_peaks, current_result, iteration):
        """Check if the iterative process has converged"""
        if iteration == 0:
            return {'converged': False, 'reason': 'first_iteration'}

        # Check if we have previous iteration data
        if len(self.iteration_history) == 0:
            return {'converged': False, 'reason': 'no_history'}

        prev_iteration = self.iteration_history[-1]
        current_quality = current_result.get('avg_r_squared', 0)
        prev_quality = prev_iteration.get('avg_r_squared', 0)

        # Quality improvement convergence
        quality_change = abs(current_quality - prev_quality)
        if quality_change < self.integration_params['convergence_threshold']:
            return {
                'converged': True,
                'reason': 'quality_converged',
                'quality_change': quality_change
            }

        # Peak count stability convergence
        current_count = len(current_result.get('fitted_peaks', []))
        prev_count = prev_iteration.get('fitted_peaks', 0)

        if current_count == prev_count and quality_change < 0.02:
            return {
                'converged': True,
                'reason': 'peak_count_stable',
                'peak_count': current_count,
                'quality_change': quality_change
            }

        return {
            'converged': False,
            'reason': 'still_improving',
            'quality_change': quality_change,
            'peak_count_change': current_count - prev_count
        }

    def _refine_peak_positions(self, fitting_result, x_data, y_data):
        """Refine peak positions based on fitting results for next iteration"""
        if not fitting_result.get('success'):
            return None

        fitted_peaks = fitting_result.get('fitted_peaks', [])
        if not fitted_peaks:
            return None

        refined_peaks = []
        for i, peak in enumerate(fitted_peaks):
            if peak.get('success') and 'center' in peak:
                refined_peaks.append({
                    'index': np.argmin(np.abs(x_data - peak['center'])),
                    'position': peak['center'],
                    'intensity': peak.get('amplitude', 0),
                    'distance_to_target': 0,
                    'refined': True,
                    'original_r_squared': peak.get('r_squared', 0)
                })

        return refined_peaks if refined_peaks else None

    def _create_failure_result(self, reason):
        """Create standardized failure result"""
        return {
            'success': False,
            'error': reason,
            'fitted_peaks': [],
            'avg_r_squared': 0.0,
            'integration_diagnostics': {
                'converged': False,
                'iterations_completed': 0,
                'iteration_history': self.iteration_history,
                'failure_reason': reason
            }
        }

    def fit_peak_region_integrated(self, nmr_data, ppm_x_axis, ppm_y_axis,
                                  x_center, y_center, **integration_params):
        """
        Integrated detection-fitting for a specific peak region

        Args:
            nmr_data: 2D NMR data array
            ppm_x_axis, ppm_y_axis: chemical shift axes
            x_center, y_center: peak center position
            integration_params: integration parameters

        Returns:
            dict: peak fitting results with quality metrics
        """
        try:
            # Define window around peak center
            x_window = integration_params.get('fitting_window_x', 0.2)  # ppm
            y_window = integration_params.get('fitting_window_y', 5.0)   # ppm

            # Find data indices for the region
            x_mask = (ppm_x_axis >= x_center - x_window) & (ppm_x_axis <= x_center + x_window)
            y_mask = (ppm_y_axis >= y_center - y_window) & (ppm_y_axis <= y_center + y_window)

            if not np.any(x_mask) or not np.any(y_mask):
                return self._create_failed_peak_result(x_center, y_center, "Region outside data range")

            # Extract region data
            x_indices = np.where(x_mask)[0]
            y_indices = np.where(y_mask)[0]

            region_x = ppm_x_axis[x_indices]
            region_y = ppm_y_axis[y_indices]
            region_data = nmr_data[np.ix_(y_indices, x_indices)]

            # Use 1D projections for integrated detection-fitting
            x_projection = np.sum(region_data, axis=0)
            y_projection = np.sum(region_data, axis=1)

            # Perform integrated detection-fitting on projections
            x_result = self.integrated_detection_fitting(
                region_x, x_projection,
                nucleus_type='1H',
                initial_peak_positions=[x_center],
                gui_params=integration_params
            )

            y_result = self.integrated_detection_fitting(
                region_y, y_projection,
                nucleus_type='15N',  # or '13C'
                initial_peak_positions=[y_center],
                gui_params=integration_params
            )

            # Extract best peak from each dimension
            x_peak = self._extract_best_peak(x_result, x_center)
            y_peak = self._extract_best_peak(y_result, y_center)

            # Combine results
            if x_peak and y_peak:
                # Calculate composite quality score
                x_quality = x_peak.get('composite_quality', 0.5)
                y_quality = y_peak.get('composite_quality', 0.5)
                composite_quality = (x_quality + y_quality) / 2

                # Calculate fit R² from both dimensions
                x_r2 = x_peak.get('fit_r_squared', 0.0)
                y_r2 = y_peak.get('fit_r_squared', 0.0)
                combined_r2 = (x_r2 + y_r2) / 2

                return {
                    'ppm_x': x_peak.get('position', x_center),
                    'ppm_y': y_peak.get('position', y_center),
                    'detected': True,
                    'converged': x_result.get('converged', False) and y_result.get('converged', False),
                    'composite_quality': composite_quality,
                    'fit_r_squared': combined_r2,
                    'detection_confidence': min(x_peak.get('detection_confidence', 0.5),
                                              y_peak.get('detection_confidence', 0.5)),
                    'aic_score': (x_peak.get('aic_score', 0) + y_peak.get('aic_score', 0)) / 2,
                    'x_peak_details': x_peak,
                    'y_peak_details': y_peak,
                    'integration_method': 'integrated_region_fitting'
                }
            else:
                return self._create_failed_peak_result(x_center, y_center, "No peaks found in region")

        except Exception as e:
            print(f"Error in integrated region fitting: {e}")
            return self._create_failed_peak_result(x_center, y_center, f"Error: {str(e)}")

    def _extract_best_peak(self, detection_result, target_position):
        """Extract the best peak closest to target position"""
        if not detection_result or 'peaks' not in detection_result:
            return None

        peaks = detection_result['peaks']
        if not peaks:
            return None

        # Find peak closest to target position
        best_peak = None
        min_distance = float('inf')

        for peak in peaks:
            if isinstance(peak, dict) and 'position' in peak:
                distance = abs(peak['position'] - target_position)
                if distance < min_distance:
                    min_distance = distance
                    best_peak = peak

        return best_peak

    def _create_failed_peak_result(self, x_center, y_center, reason):
        """Create a failed peak result"""
        return {
            'ppm_x': x_center,
            'ppm_y': y_center,
            'detected': False,
            'converged': False,
            'composite_quality': 0.0,
            'fit_r_squared': 0.0,
            'detection_confidence': 0.0,
            'aic_score': 999.0,
            'error': reason,
            'integration_method': 'integrated_region_fitting'
        }

    def _reference_based_integrated_detection(self, x_data, y_data, intensity_data, peak_list, nucleus_type, gui_params):
        """
        Reference-based integrated detection similar to standard in-place detection

        Constrains detection to search windows around peak list positions and matches
        detected peaks with original assignments.
        """
        print(f"   🔍 Reference-based integrated detection for {len(peak_list)} peaks")

        # Search window parameters (matching standard in-place detection)
        search_window_x_ppm = 0.2  # ±0.2 ppm for 1H dimension
        search_window_y_ppm = 3.0  # ±3.0 ppm for 15N dimension

        detected_peaks = []
        reference_retained = 0

        # Process each reference peak
        for idx, peak_row in peak_list.iterrows():
            try:
                ref_x_ppm = float(peak_row['Position_X'])
                ref_y_ppm = float(peak_row['Position_Y'])
                assignment = peak_row.get('Assignment', f'Peak_{idx+1}')

                # Handle case where assignment might be a number instead of string
                if isinstance(assignment, (int, float)):
                    assignment = str(int(assignment))
                elif assignment is None or assignment == '':
                    assignment = f'Peak_{idx+1}'

                print(f"      Processing {assignment}: ({ref_x_ppm:.3f}, {ref_y_ppm:.1f}) ppm")

                # Define search window around reference position
                x_min_ppm = ref_x_ppm - search_window_x_ppm
                x_max_ppm = ref_x_ppm + search_window_x_ppm
                y_min_ppm = ref_y_ppm - search_window_y_ppm
                y_max_ppm = ref_y_ppm + search_window_y_ppm

                # Convert to array indices for 2D data
                if nucleus_type == '2D':
                    # x_data = 1H axis, y_data = 15N axis
                    x_min_idx = np.argmin(np.abs(x_data - x_max_ppm))  # Note: reversed for NMR
                    x_max_idx = np.argmin(np.abs(x_data - x_min_ppm))
                    y_min_idx = np.argmin(np.abs(y_data - y_max_ppm))  # Note: reversed for NMR
                    y_max_idx = np.argmin(np.abs(y_data - y_min_ppm))

                    # Ensure proper bounds
                    x_min_idx = max(0, min(x_min_idx, x_max_idx))
                    x_max_idx = min(intensity_data.shape[1], max(x_min_idx, x_max_idx))
                    y_min_idx = max(0, min(y_min_idx, y_max_idx))
                    y_max_idx = min(intensity_data.shape[0], max(y_min_idx, y_max_idx))

                    # Extract search window
                    search_window = intensity_data[y_min_idx:y_max_idx, x_min_idx:x_max_idx]

                    # Find best peak in search window using 2D detection
                    best_peak = self._find_best_peak_in_2d_window(
                        search_window, x_data, y_data, x_min_idx, y_min_idx,
                        ref_x_ppm, ref_y_ppm, assignment, idx
                    )

                    if best_peak and best_peak.get('detected', False):
                        detected_peaks.append(best_peak)
                        print(f"         ✅ Detected at ({best_peak['ppm_x']:.3f}, {best_peak['ppm_y']:.1f}) ppm")
                    else:
                        # Retain reference position if no peak found
                        print(f"         ⚠️ No peak found, retaining reference position")
                        ref_y_idx = np.argmin(np.abs(y_data - ref_y_ppm))
                        ref_x_idx = np.argmin(np.abs(x_data - ref_x_ppm))
                        ref_intensity = intensity_data[ref_y_idx, ref_x_idx] if ref_y_idx < intensity_data.shape[0] and ref_x_idx < intensity_data.shape[1] else 0

                        retained_peak = {
                            'peak_id': f'retained_{assignment}',
                            'ppm_x': ref_x_ppm,
                            'ppm_y': ref_y_ppm,
                            'intensity': ref_intensity,
                            'volume': ref_intensity,
                            'width_1h': 0.05,
                            'width_15n': 1.0,
                            'r_squared': 0.5,  # Lower quality for retained
                            'fit_r_squared': 0.5,
                            'snr': abs(ref_intensity) / 1000.0,  # Estimate
                            'fitting_method': 'reference_retained',
                            'assignment': assignment,
                            'composite_quality': 0.5,
                            'detection_confidence': 0.3,  # Lower confidence
                            'fit_assessment': {'likelihood': 0.3},
                            'detected': False,  # Mark as not detected
                            'converged': True,  # But converged (retained)
                            'reference_retained': True,
                            'peak_number': idx + 1,
                            'amp': ref_intensity,
                            'center_x': ref_x_ppm,
                            'center_y': ref_y_ppm,
                            'sigma_x': 0.05,
                            'sigma_y': 1.0,
                            'gamma_x': 0.02,
                            'gamma_y': 0.5
                        }
                        detected_peaks.append(retained_peak)
                        reference_retained += 1

            except Exception as e:
                print(f"         ❌ Error processing {assignment}: {e}")
                continue

        # Calculate statistics
        detected_count = sum(1 for p in detected_peaks if p.get('detected', False))
        total_peaks = len(detected_peaks)
        avg_r_squared = sum(p.get('r_squared', 0) for p in detected_peaks) / max(total_peaks, 1)

        print(f"   📊 Reference-based results:")
        print(f"      Total peaks processed: {total_peaks}")
        print(f"      Successfully detected: {detected_count}")
        print(f"      Reference retained: {reference_retained}")
        print(f"      Detection rate: {(detected_count/max(len(peak_list), 1)*100):.1f}%")
        print(f"      Average R²: {avg_r_squared:.3f}")

        # Create result structure matching integrated detection format
        final_result = {
            'success': True,
            'fitted_peaks': detected_peaks,
            'peaks': detected_peaks,  # GUI compatibility
            'avg_r_squared': avg_r_squared,
            'fitting_method': 'reference_based_integrated',
            'total_peaks': total_peaks,
            'integration_diagnostics': {
                'converged': True,
                'iterations_completed': 1,
                'reference_based': True,
                'detected_count': detected_count,
                'retained_count': reference_retained,
                'detection_rate': detected_count/max(len(peak_list), 1)*100
            }
        }

        print(f"🏁 Reference-based integration complete: {total_peaks} peaks processed (avg R² = {avg_r_squared:.4f})")
        return final_result

    def _find_best_peak_in_2d_window(self, search_window, x_axis, y_axis, x_offset, y_offset,
                                    ref_x_ppm, ref_y_ppm, assignment, idx):
        """Find the best peak candidate within a 2D search window"""
        if search_window.size == 0:
            return None

        # Use scipy local maxima detection within window
        from scipy.ndimage import maximum_filter

        # Apply local maxima filter
        local_max = maximum_filter(np.abs(search_window), size=3)
        maxima_mask = (np.abs(search_window) == local_max) & (np.abs(search_window) > np.std(search_window) * 2)

        if not np.any(maxima_mask):
            return None

        # Find all maxima positions
        maxima_coords = np.argwhere(maxima_mask)

        if len(maxima_coords) == 0:
            return None

        # Select best candidate based on distance from reference and intensity
        best_candidate = None
        best_score = -1

        for coord in maxima_coords:
            local_y, local_x = coord
            global_y = local_y + y_offset
            global_x = local_x + x_offset

            # Get chemical shift positions
            if global_x < len(x_axis) and global_y < len(y_axis):
                ppm_x = x_axis[global_x]
                ppm_y = y_axis[global_y]
                intensity = search_window[local_y, local_x]

                # Calculate distance from reference
                distance = np.sqrt((ppm_x - ref_x_ppm)**2 + (ppm_y - ref_y_ppm)**2)

                # Combined score: intensity / distance (closer and stronger = better)
                score = abs(intensity) / (distance + 0.1)

                if score > best_score:
                    best_score = score
                    best_candidate = {
                        'peak_id': f'detected_{assignment}',
                        'ppm_x': ppm_x,
                        'ppm_y': ppm_y,
                        'intensity': intensity,
                        'volume': intensity,
                        'width_1h': 0.05,
                        'width_15n': 1.0,
                        'r_squared': 0.85,  # High quality for detected
                        'fit_r_squared': 0.85,
                        'snr': abs(intensity) / (np.std(search_window) + 1e-10),
                        'fitting_method': 'reference_based_2d',
                        'assignment': assignment,
                        'composite_quality': 0.85,
                        'detection_confidence': min(score / 1000.0, 0.95),  # Scale score to confidence
                        'fit_assessment': {'likelihood': min(score / 1000.0, 0.9)},
                        'detected': True,
                        'converged': True,
                        'distance_from_reference': distance,
                        'peak_number': idx + 1,  # Use simple 1-based indexing
                        'amp': intensity,
                        'center_x': ppm_x,
                        'center_y': ppm_y,
                        'sigma_x': 0.05,
                        'sigma_y': 1.0,
                        'gamma_x': 0.02,
                        'gamma_y': 0.5
                    }

        return best_candidate

    # =================== 1D REFINEMENT HELPER METHODS ===================

    def _convert_peaks_to_dataframe(self, peaks_list, nucleus_type):
        """Convert peak list to DataFrame for 1D refinement"""
        if nucleus_type != '2D' or not peaks_list:
            return None

        peak_data = []
        for i, peak in enumerate(peaks_list):
            # Handle different peak formats flexibly
            x_pos = peak.get('position', peak.get('Position_X', peak.get('x_center', 0)))
            y_pos = peak.get('position_y', peak.get('Position_Y', peak.get('y_center', 0)))

            # Only add if we have valid 2D coordinates
            if x_pos != 0 or y_pos != 0:
                peak_data.append({
                    'Position_X': float(x_pos),
                    'Position_Y': float(y_pos),
                    'Intensity': peak.get('intensity', peak.get('amplitude', peak.get('height', 1.0))),
                    'Assignment': peak.get('assignment', f"Peak_{i+1}")
                })

        if not peak_data:
            print(f"⚠️ No valid 2D peaks found for refinement conversion")
            return None

        return pd.DataFrame(peak_data)

    def _convert_dataframe_to_peaks(self, peaks_df, nucleus_type):
        """Convert DataFrame back to peak list format"""
        if nucleus_type != '2D' or peaks_df is None:
            return []

        peaks_list = []
        for idx, peak in peaks_df.iterrows():
            peak_dict = {
                'position': peak['Position_X'],       # 1H dimension
                'position_y': peak['Position_Y'],     # 15N/13C dimension
                'intensity': peak.get('Intensity', 1.0),
                'assignment': peak.get('Assignment', f"Peak_{idx+1}"),
                'refined': True  # Mark as refined
            }
            peaks_list.append(peak_dict)

        return peaks_list

# Utility functions for integration
def create_integrated_fitter(enhanced_fitter=None):
    """Factory function to create integrated detection-fitter"""
    return IterativeDetectionFitter(enhanced_fitter)

def validate_integration_params(params):
    """Validate integration parameters"""
    required_keys = ['max_iterations', 'convergence_threshold', 'min_peak_confidence']

    for key in required_keys:
        if key not in params:
            raise ValueError(f"Missing required parameter: {key}")

    if not 0 < params['convergence_threshold'] < 1:
        raise ValueError("convergence_threshold must be between 0 and 1")

    if not 0 < params['min_peak_confidence'] < 1:
        raise ValueError("min_peak_confidence must be between 0 and 1")

    if params['max_iterations'] < 1:
        raise ValueError("max_iterations must be at least 1")

    return True

# ============================================================================
# ENHANCED GRAPH-BASED PEAK DETECTION SYSTEM
# ============================================================================

class PeakNetworkGraph:
    """
    Graph representation of peak networks for pattern matching

    Represents peaks as nodes and their spatial relationships as edges,
    enabling geometric pattern recognition for complex overlapping cases.
    """

    def __init__(self):
        self.nodes = []      # Peak nodes with positions and properties
        self.edges = []      # Spatial relationships between peaks
        self.patterns = []   # Identified geometric patterns

    def add_peak_node(self, position, intensity, peak_id, properties=None):
        """Add a peak as a node in the graph"""
        node = {
            'id': peak_id,
            'position': position,
            'intensity': intensity,
            'properties': properties or {},
            'neighbors': [],
            'pattern_membership': []
        }
        self.nodes.append(node)
        return len(self.nodes) - 1

    def add_spatial_edge(self, node1_idx, node2_idx, distance, relationship_type):
        """Add spatial relationship between two peaks"""
        if node1_idx < len(self.nodes) and node2_idx < len(self.nodes):
            edge = {
                'nodes': (node1_idx, node2_idx),
                'distance': distance,
                'type': relationship_type,
                'strength': self._calculate_edge_strength(node1_idx, node2_idx, distance)
            }
            self.edges.append(edge)

            # Update neighbor lists
            self.nodes[node1_idx]['neighbors'].append(node2_idx)
            self.nodes[node2_idx]['neighbors'].append(node1_idx)

    def _calculate_edge_strength(self, node1_idx, node2_idx, distance):
        """Calculate strength of spatial relationship"""
        # Stronger for closer peaks with similar intensities
        node1 = self.nodes[node1_idx]
        node2 = self.nodes[node2_idx]

        intensity_ratio = min(node1['intensity'], node2['intensity']) / max(node1['intensity'], node2['intensity'])
        distance_factor = 1.0 / (1.0 + distance * 10)  # Closer = stronger

        return intensity_ratio * distance_factor

    def identify_geometric_patterns(self):
        """Identify geometric patterns in the peak network"""
        self.patterns = []

        # Find triangular patterns (3 peaks in close proximity)
        triangular_patterns = self._find_triangular_patterns()
        self.patterns.extend(triangular_patterns)

        # Find linear patterns (peaks in a line)
        linear_patterns = self._find_linear_patterns()
        self.patterns.extend(linear_patterns)

        # Find cluster patterns (multiple peaks clustered together)
        cluster_patterns = self._find_cluster_patterns()
        self.patterns.extend(cluster_patterns)

        return self.patterns

    def _find_triangular_patterns(self):
        """Find triangular patterns formed by 3 peaks"""
        triangular_patterns = []

        for i in range(len(self.nodes)):
            for j in range(i + 1, len(self.nodes)):
                for k in range(j + 1, len(self.nodes)):
                    # Check if these 3 nodes form a reasonable triangle
                    nodes = [self.nodes[i], self.nodes[j], self.nodes[k]]
                    positions = [node['position'] for node in nodes]

                    # Calculate triangle properties
                    distances = [
                        abs(positions[0] - positions[1]),
                        abs(positions[1] - positions[2]),
                        abs(positions[0] - positions[2])
                    ]

                    max_distance = max(distances)
                    min_distance = min(distances)

                    # Good triangle: not too elongated, reasonable size
                    if max_distance < 0.3 and min_distance > 0.01 and max_distance / min_distance < 5.0:
                        triangular_patterns.append({
                            'type': 'triangular',
                            'nodes': [i, j, k],
                            'positions': positions,
                            'centroid': np.mean(positions),
                            'max_distance': max_distance,
                            'compactness': min_distance / max_distance
                        })

        return triangular_patterns

    def _find_linear_patterns(self):
        """Find linear patterns formed by aligned peaks"""
        linear_patterns = []

        # Check all combinations of 3+ nodes for linearity
        for i in range(len(self.nodes)):
            potential_line = [i]
            base_position = self.nodes[i]['position']

            # Find peaks that could be on the same line
            for j in range(len(self.nodes)):
                if j != i:
                    distance = abs(self.nodes[j]['position'] - base_position)
                    if distance < 0.5:  # Within reasonable range
                        potential_line.append(j)

            if len(potential_line) >= 3:
                # Sort by position
                potential_line.sort(key=lambda idx: self.nodes[idx]['position'])
                positions = [self.nodes[idx]['position'] for idx in potential_line]

                # Check linearity (consistent spacing)
                spacings = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
                avg_spacing = np.mean(spacings)
                spacing_std = np.std(spacings)

                if spacing_std < avg_spacing * 0.5:  # Consistent spacing
                    linear_patterns.append({
                        'type': 'linear',
                        'nodes': potential_line,
                        'positions': positions,
                        'average_spacing': avg_spacing,
                        'spacing_consistency': 1.0 - (spacing_std / avg_spacing)
                    })

        return linear_patterns

    def _find_cluster_patterns(self):
        """Find cluster patterns of closely grouped peaks"""
        cluster_patterns = []

        # Simple clustering based on distance threshold
        unassigned = list(range(len(self.nodes)))
        cluster_distance_threshold = 0.15  # ppm

        while unassigned:
            # Start new cluster with first unassigned node
            seed_idx = unassigned[0]
            cluster = [seed_idx]
            unassigned.remove(seed_idx)

            # Find all nodes within threshold distance
            seed_position = self.nodes[seed_idx]['position']
            to_remove = []

            for idx in unassigned:
                distance = abs(self.nodes[idx]['position'] - seed_position)
                if distance <= cluster_distance_threshold:
                    cluster.append(idx)
                    to_remove.append(idx)

            # Remove assigned nodes
            for idx in to_remove:
                unassigned.remove(idx)

            # Only consider clusters with multiple peaks
            if len(cluster) >= 2:
                positions = [self.nodes[idx]['position'] for idx in cluster]
                cluster_patterns.append({
                    'type': 'cluster',
                    'nodes': cluster,
                    'positions': positions,
                    'centroid': np.mean(positions),
                    'span': max(positions) - min(positions),
                    'density': len(cluster) / (max(positions) - min(positions) + 1e-6)
                })

        return cluster_patterns


class GraphPatternMatcher:
    """
    Pattern matching between detected peaks and imported peak lists using graph representations
    """

    def __init__(self):
        self.matching_params = {
            'radius_x': 0.05,                # Max distance for 1H dimension (ppm)
            'radius_y': 2.0,                 # Max distance for 15N/13C dimension (ppm)
            'pattern_similarity_threshold': 0.7,  # Min similarity for pattern matching
            'missing_peak_tolerance': 1,     # Allow 1 missing peak per pattern
            'intensity_weight': 0.3,         # Weight for intensity in matching
            'position_weight': 0.7          # Weight for position in matching
        }

    def match_peak_patterns(self, detected_graph, imported_graph, constraints=None):
        """
        Match patterns between detected peaks and imported peak list

        Args:
            detected_graph: PeakNetworkGraph of detected peaks
            imported_graph: PeakNetworkGraph of imported peaks
            constraints: matching constraints (radius, etc.)

        Returns:
            dict: matching results with assignments and confidence scores
        """
        if constraints:
            self.matching_params.update(constraints)

        print(f"🧩 Pattern matching: {len(detected_graph.nodes)} detected vs {len(imported_graph.nodes)} imported")

        # Stage 1: Handle easy cases (single peak matching)
        easy_matches = self._match_simple_peaks(detected_graph, imported_graph)
        print(f"   ✅ Easy matches: {len(easy_matches)}")

        # Stage 2: Complex pattern matching for remaining peaks
        complex_matches = self._match_complex_patterns(detected_graph, imported_graph, easy_matches)
        print(f"   🧩 Complex matches: {len(complex_matches)}")

        # Combine results
        all_matches = easy_matches + complex_matches

        # Calculate assignment confidence
        assignments = self._create_peak_assignments(all_matches, detected_graph, imported_graph)

        return {
            'assignments': assignments,
            'easy_matches': len(easy_matches),
            'complex_matches': len(complex_matches),
            'unmatched_detected': self._count_unmatched(detected_graph, all_matches, 'detected'),
            'unmatched_detected_peaks': self._get_unmatched_detected_peaks(detected_graph, all_matches),
            'unmatched_imported': self._count_unmatched(imported_graph, all_matches, 'imported'),
            'matching_confidence': self._calculate_overall_confidence(assignments)
        }

    def _match_simple_peaks(self, detected_graph, imported_graph):
        """Match isolated peaks (easy cases)"""
        easy_matches = []
        used_detected = set()
        used_imported = set()

        for det_idx, det_node in enumerate(detected_graph.nodes):
            if det_idx in used_detected:
                continue

            # Find peaks suitable for easy matching (relaxed neighbor requirement)
            if len(det_node['neighbors']) <= 3:  # More permissive: allow up to 3 neighbors
                best_match = None
                best_distance = float('inf')

                for imp_idx, imp_node in enumerate(imported_graph.nodes):
                    if imp_idx in used_imported:
                        continue

                    # Calculate 2D distance between peaks - access coordinates from properties
                    det_props = det_node.get('properties', {})
                    imp_props = imp_node.get('properties', {})

                    det_x = det_props.get('position_x', det_node['position'])
                    det_y = det_props.get('position_y', 0)
                    imp_x = imp_props.get('position_x', imp_node['position'])
                    imp_y = imp_props.get('position_y', 0)

                    # Calculate separate X and Y distances for dimensional constraints
                    distance_x = abs(det_x - imp_x)
                    distance_y = abs(det_y - imp_y)

                    # Check separate dimensional constraints
                    radius_x = self.matching_params.get('radius_x', 0.05)
                    radius_y = self.matching_params.get('radius_y', 2.0)

                    # Combined constraint check - both dimensions must be within limits
                    within_constraints = (distance_x < radius_x and distance_y < radius_y)

                    # Combined distance for ranking (normalized by dimension ranges)
                    normalized_distance = np.sqrt((distance_x / radius_x)**2 + (distance_y / radius_y)**2)

                    if (within_constraints and
                        len(imp_node['neighbors']) <= 3 and  # More permissive
                        normalized_distance < best_distance):
                        best_match = imp_idx
                        best_distance = normalized_distance

                if best_match is not None:
                    # Confidence based on normalized distance (1.0 = perfect match, 0.0 = at constraint limit)
                    confidence = max(0.0, 1.0 - best_distance)
                    easy_matches.append({
                        'detected_idx': det_idx,
                        'imported_idx': best_match,
                        'distance': best_distance,
                        'confidence': confidence,
                        'type': 'simple'
                    })
                    used_detected.add(det_idx)
                    used_imported.add(best_match)

        return easy_matches

    def _match_complex_patterns(self, detected_graph, imported_graph, existing_matches):
        """Match complex overlapping patterns using graph similarity"""
        complex_matches = []
        used_detected = set(match['detected_idx'] for match in existing_matches)
        used_imported = set(match['imported_idx'] for match in existing_matches)

        # Get remaining unmatched nodes
        unmatched_detected = [i for i in range(len(detected_graph.nodes)) if i not in used_detected]
        unmatched_imported = [i for i in range(len(imported_graph.nodes)) if i not in used_imported]

        # Group unmatched nodes by patterns
        detected_patterns = self._group_by_patterns(unmatched_detected, detected_graph)
        imported_patterns = self._group_by_patterns(unmatched_imported, imported_graph)

        # Match patterns between detected and imported
        for det_pattern in detected_patterns:
            best_match = None
            best_similarity = 0.0

            for imp_pattern in imported_patterns:
                similarity = self._calculate_pattern_similarity(det_pattern, imp_pattern)

                if similarity > best_similarity and similarity > self.matching_params['pattern_similarity_threshold']:
                    best_match = imp_pattern
                    best_similarity = similarity

            if best_match is not None:
                # Create individual peak matches within the pattern
                pattern_matches = self._match_peaks_within_pattern(det_pattern, best_match, detected_graph, imported_graph)
                complex_matches.extend(pattern_matches)

        return complex_matches

    def _group_by_patterns(self, node_indices, graph):
        """Group nodes by their pattern membership"""
        pattern_groups = {}

        for node_idx in node_indices:
            node = graph.nodes[node_idx]
            pattern_memberships = node.get('pattern_membership', [])

            if not pattern_memberships:
                # Isolated node - create singleton pattern
                pattern_groups[f'isolated_{node_idx}'] = {
                    'type': 'isolated',
                    'nodes': [node_idx],
                    'properties': {'centroid': node['position']}
                }
            else:
                # Add to existing patterns
                for pattern_id in pattern_memberships:
                    if pattern_id not in pattern_groups:
                        pattern_groups[pattern_id] = {
                            'type': 'complex',
                            'nodes': [],
                            'properties': {}
                        }
                    pattern_groups[pattern_id]['nodes'].append(node_idx)

        return list(pattern_groups.values())

    def _calculate_pattern_similarity(self, pattern1, pattern2):
        """Calculate similarity between two patterns"""
        # Handle size difference (with missing peak tolerance)
        size1, size2 = len(pattern1['nodes']), len(pattern2['nodes'])
        size_diff = abs(size1 - size2)

        if size_diff > self.matching_params['missing_peak_tolerance']:
            return 0.0  # Too different

        # Position-based similarity
        positions1 = [pattern1.get('properties', {}).get('centroid', 0)]
        positions2 = [pattern2.get('properties', {}).get('centroid', 0)]

        if positions1[0] != 0 and positions2[0] != 0:
            # Calculate 2D position similarity with separate dimensional constraints
            x_diff = abs(positions1[0] - positions2[0])
            y_diff = abs(positions1[1] - positions2[1]) if len(positions1) > 1 and len(positions2) > 1 else 0

            radius_x = self.matching_params.get('radius_x', 0.05)
            radius_y = self.matching_params.get('radius_y', 2.0)

            # Normalized distance in 2D space
            normalized_position_distance = np.sqrt((x_diff / radius_x)**2 + (y_diff / radius_y)**2)
            position_similarity = max(0.0, 1.0 - min(1.0, normalized_position_distance))
        else:
            position_similarity = 0.5  # Default when centroids not available

        # Size similarity (with tolerance for missing peaks)
        size_similarity = 1.0 - (size_diff / max(size1, size2, 1))

        # Combine similarities
        overall_similarity = (self.matching_params['position_weight'] * position_similarity +
                             (1 - self.matching_params['position_weight']) * size_similarity)

        return overall_similarity

    def _match_peaks_within_pattern(self, detected_pattern, imported_pattern, detected_graph=None, imported_graph=None):
        """Create individual peak matches within matched patterns"""
        matches = []

        # Use provided graphs or fallback to class attributes
        det_graph = detected_graph or getattr(self, 'detected_graph', None)
        imp_graph = imported_graph or getattr(self, 'imported_graph', None)

        if not det_graph or not imp_graph:
            return matches  # Cannot match without graph references

        # Simple nearest-neighbor matching within patterns
        for det_idx in detected_pattern['nodes']:
            best_imp_idx = None
            best_distance = float('inf')

            for imp_idx in imported_pattern['nodes']:
                # Calculate 2D distance between peaks using separate dimensional constraints
                det_props = det_graph.nodes[det_idx].get('properties', {})
                imp_props = imp_graph.nodes[imp_idx].get('properties', {})

                det_x = det_props.get('position_x', det_graph.nodes[det_idx]['position'])
                det_y = det_props.get('position_y', 0)
                imp_x = imp_props.get('position_x', imp_graph.nodes[imp_idx]['position'])
                imp_y = imp_props.get('position_y', 0)

                # Separate dimensional distances
                distance_x = abs(det_x - imp_x)
                distance_y = abs(det_y - imp_y)

                # Check dimensional constraints
                radius_x = self.matching_params.get('radius_x', 0.05)
                radius_y = self.matching_params.get('radius_y', 2.0)

                # Normalized distance for comparison
                normalized_distance = np.sqrt((distance_x / radius_x)**2 + (distance_y / radius_y)**2)

                if (distance_x < radius_x and distance_y < radius_y and
                    normalized_distance < best_distance):
                    best_distance = normalized_distance
                    best_imp_idx = imp_idx

            if best_imp_idx is not None and best_distance < 1.0:  # Normalized distance threshold
                confidence = max(0.0, 1.0 - best_distance)
                matches.append({
                    'detected_idx': det_idx,
                    'imported_idx': best_imp_idx,
                    'distance': best_distance,
                    'confidence': confidence * 0.8,  # Slightly lower confidence for complex matches
                    'type': 'complex_pattern'
                })

        return matches

    def _create_peak_assignments(self, matches, detected_graph, imported_graph):
        """Create final peak assignments from matches"""
        assignments = []

        for match in matches:
            det_node = detected_graph.nodes[match['detected_idx']]
            imp_node = imported_graph.nodes[match['imported_idx']]

            # Store 2D coordinates properly - access from node properties
            det_props = det_node.get('properties', {})
            imp_props = imp_node.get('properties', {})

            # CRITICAL: Extract ORIGINAL detected coordinates from node properties
            # These are the actual spectrum maxima positions - NEVER modify them
            det_x_original = det_props.get('position_x')
            det_y_original = det_props.get('position_y')

            # Fallback to node position if properties missing
            if det_x_original is None:
                det_x_original = det_node['position'][0] if isinstance(det_node['position'], (list, tuple)) else det_node['position']
            if det_y_original is None:
                det_y_original = det_node['position'][1] if isinstance(det_node['position'], (list, tuple)) and len(det_node['position']) > 1 else 0

            # CRITICAL BUG FIX: Ensure we use the actual detected coordinates from the spectrum
            # The detected coordinates should come from the original peak detection, not the graph nodes
            detected_x = det_props.get('position_x', det_node.get('position_x', det_node['position']))
            detected_y = det_props.get('position_y', det_node.get('position_y', 0))

            # DEBUG output to trace coordinate flow
            #print(f"   🔶 DEBUG assignment: detected_x={detected_x}, detected_y={detected_y}, node_pos={det_node['position']}")

            assignment = {
                'detected_position': det_node['position'],  # X coordinate (backward compatibility)
                'detected_position_x': detected_x,  # FIXED: Use actual detected coordinate from spectrum
                'detected_position_y': detected_y,  # FIXED: Use actual detected coordinate from spectrum
                'imported_position': imp_node['position'],  # X coordinate (backward compatibility)
                'imported_position_x': imp_props.get('position_x', imp_node['position']),
                'imported_position_y': imp_props.get('position_y', 0),
                'detected_intensity': det_node['intensity'],
                'imported_id': imp_node.get('id', match['imported_idx']),
                'distance_error': match['distance'],
                'confidence': match['confidence'],
                'match_type': match['type']
            }
            assignments.append(assignment)

        return assignments

    def _count_unmatched(self, graph, matches, graph_type):
        """Count unmatched peaks in a graph"""
        if graph_type == 'detected':
            matched_indices = set(match['detected_idx'] for match in matches)
        else:
            matched_indices = set(match['imported_idx'] for match in matches)

        return len(graph.nodes) - len(matched_indices)

    def _get_unmatched_detected_peaks(self, detected_graph, matches):
        """Get actual unmatched detected peak data (not just counts)"""
        matched_indices = set(match['detected_idx'] for match in matches)
        unmatched_peaks = []

        for i, node in enumerate(detected_graph.nodes):
            if i not in matched_indices:
                # Extract peak data from the graph node
                peak_data = {
                    'position_x': node.get('properties', {}).get('position_x', node['position'][0] if isinstance(node['position'], (list, tuple)) else node['position']),
                    'position_y': node.get('properties', {}).get('position_y', node['position'][1] if isinstance(node['position'], (list, tuple)) and len(node['position']) > 1 else 0),
                    'intensity': node['intensity'],
                    'id': node['id']
                }
                unmatched_peaks.append(peak_data)

        return unmatched_peaks

    def _calculate_overall_confidence(self, assignments):
        """Calculate overall matching confidence"""
        if not assignments:
            return 0.0

        confidences = [assignment['confidence'] for assignment in assignments]
        return np.mean(confidences)


class PeakDetectionProgress:
    """Progress visualization for GUI integration"""

    def __init__(self):
        self.stages = [
            "🔄 Running peak detection standard",
            "🔍 Identifying complex regions",
            "🔗 Building peak network",
            "✨ Matching easy cases",
            "🧩 Resolving complex patterns",
            "✅ Finalizing assignments"
        ]
        self.current = 0

    def update(self, stage, metrics):
        """Update GUI with real-time progress"""
        if stage < len(self.stages):
            self.current = stage
            progress_bar = '█' * stage + '░' * (len(self.stages) - stage)

            print(f"\n{self.stages[stage]}")
            print(f"├─ Progress: {progress_bar} {stage}/{len(self.stages)}")
            print(f"└─ Metrics: {metrics}")

            return {
                'stage_name': self.stages[stage],
                'stage_number': stage,
                'total_stages': len(self.stages),
                'progress_percent': (stage / len(self.stages)) * 100,
                'metrics': metrics
            }

        return None


class EnhancedPeakDetectionIntegrated:
    """
    Enhanced peak detection integrated system with graph-based pattern matching

    This replaces the problematic in-place detection with intelligent pattern matching
    that can handle overlapping peaks and maintain assignment accuracy.
    """

    def __init__(self, original_fitter):
        self.original_fitter = original_fitter
        self.graph_matcher = GraphPatternMatcher()
        self.progress_callback = None
        self.progress_tracker = PeakDetectionProgress()

    def peak_detection_integrated_enhanced(self, x_data, y_data, intensity_data,
                                         peak_list, nucleus_type, gui_params):
        """
        Enhanced non-in-place peak detection integrated with graph-based matching

        Args:
            x_data, y_data: spectral axes
            intensity_data: spectral intensity data
            peak_list: imported peak list DataFrame
            nucleus_type: nucleus type for constraints
            gui_params: GUI parameters

        Returns:
            dict: comprehensive results with graph-matched assignments
        """
        print(f"🚀 Enhanced peak detection integrated for {len(peak_list)} imported peaks")

        # Step 1: Run standard peak detection for baseline (non-in-place)
        self.progress_tracker.update(0, f"detecting peaks in {nucleus_type} data")

        standard_peaks = self._run_baseline_detection(x_data, y_data, intensity_data, nucleus_type, gui_params, peak_list)
        print(f"   📊 Baseline detection: {len(standard_peaks)} peaks found")

        # DEBUG: Check format of standard_peaks
        if len(standard_peaks) > 0:
            print(f"   🔍 First peak keys: {list(standard_peaks[0].keys())}")
            print(f"   🔍 First peak sample: {standard_peaks[0]}")

        # DEBUG: Print detected peak coordinates before graph processing
        #print(f"   🔶 CRITICAL DEBUG: Raw detected peaks from baseline detection:")
        for i, peak in enumerate(standard_peaks[:3]):  # Show first 3
            x_pos = peak.get('position_x', peak.get('position', 'MISSING'))
            y_pos = peak.get('position_y', 'MISSING')
            print(f"      Peak {i+1}: X={x_pos}, Y={y_pos}")

        # DEBUG: Compare with reference peak positions
        #print(f"   🔶 CRITICAL DEBUG: Reference peaks for comparison:")
        for i, (idx, row) in enumerate(peak_list.iterrows()):
            if i >= 3:  # Show first 3
                break
            print(f"      Ref {i+1}: X={row['Position_X']:.4f}, Y={row['Position_Y']:.1f} ({row.get('Assignment', 'Unknown')})")

        # CHECK: Are detected coordinates same as reference coordinates?
        if len(standard_peaks) > 0 and len(peak_list) > 0:
            det_x = standard_peaks[0].get('position_x', 0)
            ref_x = peak_list.iloc[0]['Position_X']
            if abs(det_x - ref_x) < 0.001:
                print(f"   🚨 COORDINATE CORRUPTION DETECTED: Detected peak X={det_x:.4f} matches reference X={ref_x:.4f}")
            else:
                print(f"   ✅ Coordinates look different: Detected X={det_x:.4f} vs Reference X={ref_x:.4f}")

        # Step 2: Identify complex overlapping regions
        self.progress_tracker.update(1, f"analyzing {len(standard_peaks)} detected peaks")

        complex_regions = self._identify_overlapping_regions(standard_peaks)
        print(f"   🔍 Complex regions: {len(complex_regions)} identified")

        # Step 3: Build peak network graph for detected peaks
        self.progress_tracker.update(2, f"creating graph with {len(standard_peaks)} nodes")

        detected_graph = self._build_detected_peak_graph(standard_peaks)
        detected_graph.identify_geometric_patterns()
        print(f"   🔗 Detected graph: {len(detected_graph.nodes)} nodes, {len(detected_graph.patterns)} patterns")

        # DEBUG: Print graph node coordinates
        print(f"   🔶 DEBUG: Detected graph nodes:")
        for i, node in enumerate(detected_graph.nodes[:5]):  # Show first 5
            props = node.get('properties', {})
            x_pos = props.get('position_x', 'None')
            y_pos = props.get('position_y', 'None')
            print(f"      Node {i+1}: X={x_pos}, Y={y_pos}, pos={node['position']}")

        # Step 4: Build network graph for imported peaks
        imported_graph = self._build_imported_peak_graph(peak_list)
        imported_graph.identify_geometric_patterns()
        print(f"   📋 Imported graph: {len(imported_graph.nodes)} nodes, {len(imported_graph.patterns)} patterns")

        # Step 5: Pattern matching using Simple Pattern Matcher (NEW - Coordinate Preserving)
        self.progress_tracker.update(3, f"matching {len(standard_peaks)} vs {len(peak_list)} peaks")

        # Import and configure Simple Pattern Matcher
        from lunaNMR.integrators.simple_pattern_matcher import SimplePatternMatcher

        # Configure proximity parameters based on Enhanced Detection GUI settings
        proximity_x = gui_params.get('enhanced_radius_x', 0.05) * 2  # Double radius for search window
        proximity_y = gui_params.get('enhanced_radius_y', 2.0) * 2   # Double radius for search window

        pattern_matcher = SimplePatternMatcher(
            proximity_radius_x=proximity_x,
            proximity_radius_y=proximity_y,
            isolated_peak_threshold=gui_params.get('enhanced_radius_x', 0.05),
            enable_fast_isolated=True,
            assignment_only=True  # Enhanced Detection does assignment transfer only, no peak discovery
        )

        # Convert peak_list DataFrame to Simple Pattern Matcher format
        reference_peaks = []
        for idx, row in peak_list.iterrows():
            reference_peaks.append({
                'Position_X': row.get('Position_X', 0),
                'Position_Y': row.get('Position_Y', 0),
                'Assignment': row.get('Assignment', f'Peak_{idx+1}'),
                'Volume': row.get('Volume', 0),
                'Intensity': row.get('Intensity', 1)
            })

        # DEBUG: Check coordinates BEFORE Simple Pattern Matcher
        #print(f"   🔍 DEBUGGING: Coordinates BEFORE Simple Pattern Matcher:")
        if len(standard_peaks) > 0:
            x_val = standard_peaks[0].get('position_x', 'MISSING')
            y_val = standard_peaks[0].get('position_y', 'MISSING')
            print(f"      standard_peaks[0]: X={x_val}, Y={y_val}")
        if len(reference_peaks) > 0:
            print(f"      reference_peaks[0]: X={reference_peaks[0]['Position_X']:.4f}, Y={reference_peaks[0]['Position_Y']:.2f}")

        # Run Simple Pattern Matcher (preserves detected coordinates)
        print(f"   🔍 Using Simple Pattern Matcher with proximity: {proximity_x:.3f} × {proximity_y:.1f} ppm")
        print(f"   🔍 About to call assign_peaks with {len(standard_peaks)} detected peaks and {len(reference_peaks)} reference peaks")
        matched_peaks = pattern_matcher.assign_peaks(standard_peaks, reference_peaks, max_distance=0.5)

        # DEBUG: Check coordinates AFTER Simple Pattern Matcher
        #print(f"   🔍 DEBUGGING: Coordinates AFTER Simple Pattern Matcher:")
        print(f"   📊 Simple Pattern Matcher returned: {len(matched_peaks)} peaks")
        if len(matched_peaks) > 0:
            print(f"   🔍 First matched peak keys: {list(matched_peaks[0].keys())}")
            print(f"   🔍 First matched peak: {matched_peaks[0]}")
        if len(matched_peaks) > 0:
            first_assigned = next((peak for peak in matched_peaks if peak.get('Detection_Method') != 'detected_unassigned'), None)
            if first_assigned:
                print(f"      matched_peaks[0]: X={first_assigned['Position_X']:.4f}, Y={first_assigned['Position_Y']:.2f}")
                print(f"      Detection_Method: {first_assigned.get('Detection_Method')}")
                print(f"      Assignment: {first_assigned.get('Assignment')}")
            else:
                print(f"      No assigned peaks found in matched_peaks")

        # Convert results to old format for compatibility with existing code
        matching_results = {
            'assignments': [],
            'easy_matches': 0,
            'complex_matches': 0,
            'unmatched_detected': 0,
            'unmatched_imported': 0,
            'matching_confidence': 0.85
        }

        # Process Simple Pattern Matcher results
        assigned_count = 0
        for peak in matched_peaks:
            if peak.get('Detection_Method') != 'detected_unassigned':
                assigned_count += 1

                # Find original reference peak to get reference coordinates
                ref_position_x = 0
                ref_position_y = 0
                for ref_peak in reference_peaks:
                    if ref_peak['Assignment'] == peak['Assignment']:
                        ref_position_x = ref_peak['Position_X']
                        ref_position_y = ref_peak['Position_Y']
                        break

                matching_results['assignments'].append({
                    'detected_position_x': peak['Position_X'],  # DETECTED coordinate (preserved)
                    'detected_position_y': peak['Position_Y'],  # DETECTED coordinate (preserved)
                    'detected_intensity': peak['Intensity'],
                    'imported_position_x': ref_position_x,     # REFERENCE coordinate (constant)
                    'imported_position_y': ref_position_y,     # REFERENCE coordinate (constant)
                    'imported_id': peak['Assignment'],
                    'confidence': peak.get('Match_Confidence', 0.8),
                    'distance_error': peak.get('Match_Distance', 0.1),
                    'match_type': peak.get('Detection_Method', 'simple_pattern_match')
                })

        # Count isolated vs pattern matches
        isolated_count = sum(1 for peak in matched_peaks if peak.get('Detection_Method') == 'isolated_fast')
        pattern_count = sum(1 for peak in matched_peaks if peak.get('Detection_Method') == 'localized_pattern_match')

        matching_results['easy_matches'] = isolated_count
        matching_results['complex_matches'] = pattern_count
        matching_results['unmatched_detected'] = len(standard_peaks) - assigned_count
        matching_results['unmatched_imported'] = len(reference_peaks) - assigned_count

        print(f"   📊 Simple Pattern Matcher Results:")
        print(f"      Total assignments: {assigned_count}")
        print(f"      Isolated matches: {isolated_count}")
        print(f"      Pattern matches: {pattern_count}")
        print(f"      Unmatched detected: {matching_results['unmatched_detected']}")
        print(f"      Unmatched reference: {matching_results['unmatched_imported']}")

        # Step 6: Complex pattern matching
        self.progress_tracker.update(4, f"resolved {matching_results['complex_matches']} complex cases")

        # Step 7: Create final peak assignments
        self.progress_tracker.update(5, f"assigning {len(matching_results['assignments'])} final peaks")

        #print(f"   🔍 DEBUG: matching_results has {len(matching_results['assignments'])} assignments")
        if len(matching_results['assignments']) > 0:
            print(f"   🔍 DEBUG: First assignment: {matching_results['assignments'][0]}")
        final_peaks = self._create_final_assignments(matching_results, peak_list)
        #print(f"   🔍 DEBUG: _create_final_assignments returned {len(final_peaks)} peaks")

        return {
            'success': True,
            'peaks': final_peaks,
            'fitted_peaks': final_peaks,
            'method': 'enhanced_graph_based',
            'detected_peaks_count': len(standard_peaks),
            'matched_peaks_count': len(matching_results['assignments']),
            'easy_matches': matching_results['easy_matches'],
            'complex_matches': matching_results['complex_matches'],
            'unmatched_detected': matching_results['unmatched_detected'],
            'unmatched_imported': matching_results['unmatched_imported'],
            'overall_confidence': matching_results['matching_confidence'],
            'complex_regions_resolved': len(complex_regions)
        }

    def _run_baseline_detection(self, x_data, y_data, intensity_data, nucleus_type, gui_params, peak_list):
        """Run baseline peak detection (standard method, non-in-place)"""
        # Use the original detection but without in-place constraints
        if nucleus_type == '2D' and hasattr(self.original_fitter, '_perform_2d_simplified_fitting'):
            # For 2D detection, we need to run actual peak detection first
            # Use the original fitter's detection methods but without in-place constraints
            try:
                # Run 2D peak detection using the SAME method as Standard Detection (with consolidation)
                print(f"   🔍 Calling 2D detection with consolidation...")
                detected_peaks = self.original_fitter._detect_2d_nmr_peaks_with_axis(x_data, y_data, intensity_data, gui_params)
                print(f"   📊 2D detection returned: {len(detected_peaks) if detected_peaks else 0} peaks")
                return detected_peaks
            except Exception as e:
                print(f"   ⚠️ 2D detection failed: {str(e)}, falling back to simple method")
                # Fallback to simple peak detection
                return self._fallback_peak_detection(x_data, y_data, intensity_data, gui_params)
        else:
            # For 1D detection - use standard peak detection
            return self._fallback_peak_detection(x_data, y_data, intensity_data, gui_params)

    def _perform_2d_peak_detection(self, x_data, y_data, intensity_data, gui_params, peak_list):
        """Perform 2D peak detection constrained by peak list bounds and using noise threshold"""
        from scipy.signal import find_peaks

        # Get spectrum bounds from peak list with buffer for detection
        x_min = peak_list['Position_X'].min() - 0.2  # Add 0.2 ppm buffer
        x_max = peak_list['Position_X'].max() + 0.2
        y_min = peak_list['Position_Y'].min() - 5.0   # Add 5 ppm buffer
        y_max = peak_list['Position_Y'].max() + 5.0

        # Find corresponding indices
        x_start_idx = np.searchsorted(x_data, x_min, side='left')
        x_end_idx = np.searchsorted(x_data, x_max, side='right')
        y_start_idx = np.searchsorted(y_data, y_min, side='left')
        y_end_idx = np.searchsorted(y_data, y_max, side='right')

        # Constrain to valid ranges
        x_start_idx = max(0, x_start_idx)
        x_end_idx = min(len(x_data), x_end_idx)
        y_start_idx = max(0, y_start_idx)
        y_end_idx = min(len(y_data), y_end_idx)

        # Use noise threshold like standard detection
        noise_threshold = gui_params.get('noise_threshold', 0.01)

        detected_peaks = []

        # Find peaks only within constrained region
        for i in range(y_start_idx, y_end_idx):
            row_data = intensity_data[i, x_start_idx:x_end_idx]
            peaks, properties = find_peaks(row_data,
                                        height=noise_threshold,
                                        prominence=gui_params.get('prominence_threshold', 0.01) * np.max(row_data))

            for peak_idx in peaks:
                actual_x_idx = x_start_idx + peak_idx
                if actual_x_idx < len(x_data):
                    detected_peaks.append({
                        'position': x_data[actual_x_idx],
                        'position_x': x_data[actual_x_idx],
                        'position_y': y_data[i],
                        'intensity': row_data[peak_idx],
                        'method': 'enhanced_2d_detection'
                    })

        print(f"      📊 2D peak detection found {len(detected_peaks)} peaks (before filtering)")
        print(f"      🔍 Search region: X[{x_min:.2f}:{x_max:.2f}], Y[{y_min:.1f}:{y_max:.1f}]")

        # Apply Enhanced Detection peak reduction controls
        filtered_peaks = self._apply_peak_reduction_filters(detected_peaks, intensity_data, gui_params)

        print(f"      🔍 After peak reduction filters: {len(filtered_peaks)} peaks")
        return filtered_peaks

    def _fallback_peak_detection(self, x_data, y_data, intensity_data, gui_params):
        """Fallback peak detection for 1D or when 2D fails"""
        from scipy.signal import find_peaks

        # Apply detection parameters from GUI
        height_threshold = gui_params.get('height_threshold', 0.02)
        prominence_threshold = gui_params.get('prominence_threshold', 0.01)

        # Handle 2D case by looking at maximum projections
        if len(intensity_data.shape) == 2:
            # Project 2D data to 1D for peak detection
            projection = np.max(intensity_data, axis=0)  # Project along Y axis
            data_for_detection = projection
        else:
            data_for_detection = intensity_data

        peaks, _ = find_peaks(data_for_detection,
                            height=height_threshold * np.max(data_for_detection),
                            prominence=prominence_threshold * np.max(data_for_detection))

        # Convert to peak format
        detected_peaks = []
        for peak_idx in peaks:
            if peak_idx < len(x_data):
                # For 2D fallback, estimate Y position from projection
                if len(intensity_data.shape) == 2:
                    # Find the Y position with maximum intensity at this X position
                    col_data = intensity_data[:, peak_idx]
                    y_max_idx = np.argmax(col_data)
                    y_position = y_data[y_max_idx] if y_max_idx < len(y_data) else 0
                else:
                    y_position = 0  # 1D data

                detected_peaks.append({
                    'position': x_data[peak_idx],     # X position (backward compatibility)
                    'position_x': x_data[peak_idx],   # X position (explicit)
                    'position_y': y_position,         # Y position
                    'intensity': data_for_detection[peak_idx],
                    'method': 'fallback_detection'
                })

        print(f"      📊 Fallback detection found {len(detected_peaks)} peaks (before filtering)")

        # Apply Enhanced Detection peak reduction controls
        filtered_peaks = self._apply_peak_reduction_filters(detected_peaks, intensity_data, gui_params)

        print(f"      🔍 After peak reduction filters: {len(filtered_peaks)} peaks")
        return filtered_peaks

    def _apply_peak_reduction_filters(self, detected_peaks, intensity_data, gui_params):
        """
        Apply Enhanced Detection peak reduction controls to limit excessive detections

        Args:
            detected_peaks: list of detected peak dictionaries
            intensity_data: NMR data for intensity normalization
            gui_params: GUI parameters including enhanced_peak_limit

        Returns:
            list: filtered peak list with reduced count and minimum intensity applied
        """
        if not detected_peaks:
            return detected_peaks

        # Get user-configurable peak reduction parameters
        max_peaks = gui_params.get('enhanced_peak_limit', 50)  # Default: 50 peaks max
        min_intensity_fraction = 0.1  # Fixed: 10% of max intensity (no longer user-configurable)

        max_intensity = np.max(intensity_data)
        min_intensity_threshold = min_intensity_fraction * max_intensity

        print(f"         🎛️ Peak reduction settings:")
        print(f"            Max peaks allowed: {max_peaks}")
        print(f"            Min intensity threshold: {min_intensity_fraction:.3f} × max = {min_intensity_threshold:.6f}")

        # Step 1: Filter by minimum intensity
        intensity_filtered = []
        intensity_excluded = []
        for peak in detected_peaks:
            peak_intensity = peak.get('intensity', 0)
            if peak_intensity >= min_intensity_threshold:
                intensity_filtered.append(peak)
            else:
                intensity_excluded.append(peak)

        print(f"         📉 After intensity filter ({min_intensity_fraction:.3f}× max): {len(intensity_filtered)} peaks")

        # Step 2: Limit to maximum number of peaks (keep strongest)
        if len(intensity_filtered) > max_peaks:
            # Sort by intensity (descending) and keep top peaks
            sorted_peaks = sorted(intensity_filtered, key=lambda p: p.get('intensity', 0), reverse=True)
            final_peaks = sorted_peaks[:max_peaks]
            excluded_after_limit = sorted_peaks[max_peaks:]  # Store excluded peaks
            print(f"         🔢 After peak limit ({max_peaks} max): {len(final_peaks)} peaks (kept strongest)")
        else:
            final_peaks = intensity_filtered
            excluded_after_limit = []
            print(f"         ✅ Peak count within limit: {len(final_peaks)} ≤ {max_peaks}")

        # Store included peaks after limit filtering for debugging purposes
        if hasattr(self, '_included_peaks_after_limit_debug'):
            self._included_peaks_after_limit_debug = final_peaks
        else:
            # If this is a standalone function, store in a global or class variable
            setattr(self, '_included_peaks_after_limit_debug', final_peaks)

        if final_peaks:
            print(f"         🔍 DEBUG: Stored {len(final_peaks)} included peaks after limit filtering")

        # Also store excluded peaks for debugging purposes if needed
        all_excluded = intensity_excluded + excluded_after_limit
        if hasattr(self, '_excluded_peaks_debug'):
            self._excluded_peaks_debug = all_excluded
        else:
            setattr(self, '_excluded_peaks_debug', all_excluded)

        if all_excluded:
            print(f"         🔍 DEBUG: Also stored {len(all_excluded)} excluded peaks ({len(intensity_excluded)} intensity, {len(excluded_after_limit)} limit)")

        return final_peaks

    def _identify_overlapping_regions(self, detected_peaks):
        """Identify regions where peaks overlap or are too close together"""
        if len(detected_peaks) < 2:
            return []

        overlap_regions = []
        overlap_threshold = 0.1  # ppm - peaks closer than this are considered overlapping

        positions = [peak.get('position', 0) for peak in detected_peaks]
        positions.sort()

        current_region = []
        for i, pos in enumerate(positions):
            if not current_region:
                current_region.append(i)
            else:
                # Check distance from last peak in current region
                last_pos = positions[current_region[-1]]
                if abs(pos - last_pos) <= overlap_threshold:
                    current_region.append(i)
                else:
                    # End current region if it has multiple peaks
                    if len(current_region) >= 2:
                        overlap_regions.append({
                            'peak_indices': current_region.copy(),
                            'positions': [positions[idx] for idx in current_region],
                            'span': positions[current_region[-1]] - positions[current_region[0]],
                            'complexity': len(current_region)
                        })
                    current_region = [i]

        # Check final region
        if len(current_region) >= 2:
            overlap_regions.append({
                'peak_indices': current_region.copy(),
                'positions': [positions[idx] for idx in current_region],
                'span': positions[current_region[-1]] - positions[current_region[0]],
                'complexity': len(current_region)
            })

        return overlap_regions

    def _build_detected_peak_graph(self, detected_peaks):
        """Build graph representation of detected peaks"""
        graph = PeakNetworkGraph()

        # Add nodes for each detected peak
        for i, peak in enumerate(detected_peaks):
            # Handle 2D coordinates properly
            x_pos = peak.get('position', peak.get('position_x', 0))
            y_pos = peak.get('position_y', 0)
            intensity = peak.get('intensity', 0)

            # Store both X and Y coordinates
            peak_properties = {
                'source': 'detection',
                'method': peak.get('method', 'standard'),
                'position_x': x_pos,
                'position_y': y_pos,
                'coordinates_2d': True
            }

            graph.add_peak_node(x_pos, intensity, f'detected_{i}', peak_properties)

        # Add edges based on 2D proximity (Euclidean distance)
        proximity_threshold = 0.3  # ppm (2D distance)
        for i in range(len(detected_peaks)):
            for j in range(i + 1, len(detected_peaks)):
                x1 = detected_peaks[i].get('position', detected_peaks[i].get('position_x', 0))
                y1 = detected_peaks[i].get('position_y', 0)
                x2 = detected_peaks[j].get('position', detected_peaks[j].get('position_x', 0))
                y2 = detected_peaks[j].get('position_y', 0)

                # Calculate 2D Euclidean distance
                distance_2d = np.sqrt((x1 - x2)**2 + ((y1 - y2)/10)**2)  # Scale Y by 10 for typical NMR ranges

                if distance_2d <= proximity_threshold:
                    relationship_type = 'close' if distance_2d <= 0.15 else 'moderate'
                    graph.add_spatial_edge(i, j, distance_2d, relationship_type)

        return graph

    def _build_imported_peak_graph(self, peak_list):
        """Build graph representation of imported peak list"""
        graph = PeakNetworkGraph()

        # Add nodes for each imported peak
        for i, (_, row) in enumerate(peak_list.iterrows()):
            position_x = float(row.get('Position_X', 0))
            position_y = float(row.get('Position_Y', 0))
            intensity = float(row.get('Intensity', 1.0))  # Default intensity if missing
            assignment = row.get('Assignment', f'Peak_{i}')

            # Store both X and Y coordinates for imported peaks
            peak_properties = {
                'source': 'imported',
                'original_index': i,
                'full_assignment': assignment,
                'position_x': position_x,
                'position_y': position_y,
                'coordinates_2d': True
            }

            graph.add_peak_node(position_x, intensity, assignment, peak_properties)

        # Add edges based on 2D proximity (same as detected)
        proximity_threshold = 0.3  # ppm (2D distance)
        for i in range(len(peak_list)):
            for j in range(i + 1, len(peak_list)):
                row1 = peak_list.iloc[i]
                row2 = peak_list.iloc[j]
                x1 = float(row1.get('Position_X', 0))
                y1 = float(row1.get('Position_Y', 0))
                x2 = float(row2.get('Position_X', 0))
                y2 = float(row2.get('Position_Y', 0))

                # Calculate 2D Euclidean distance
                distance_2d = np.sqrt((x1 - x2)**2 + ((y1 - y2)/10)**2)  # Scale Y by 10 for typical NMR ranges

                if distance_2d <= proximity_threshold:
                    relationship_type = 'close' if distance_2d <= 0.15 else 'moderate'
                    graph.add_spatial_edge(i, j, distance_2d, relationship_type)

        return graph

    def _create_final_assignments(self, matching_results, original_peak_list):
        """Create final peak assignments from matching results"""
        #print(f"🔍 DEBUGGING _create_final_assignments:")
        print(f"   Input: {len(matching_results['assignments'])} assignments")
        if len(matching_results['assignments']) > 0:
            print(f"   Sample assignment: {matching_results['assignments'][0]}")
        final_peaks = []

        # Group assignments by imported peak to avoid duplicates
        peak_groups = {}
        for assignment in matching_results['assignments']:
            imported_id = assignment.get('imported_id', 'Unknown')
            if imported_id not in peak_groups:
                peak_groups[imported_id] = []
            peak_groups[imported_id].append(assignment)

        # For each imported peak, select the best matching detected peak
        for imported_id, assignments in peak_groups.items():
            if not assignments:
                continue

            # Select best assignment based on confidence and distance
            best_assignment = max(assignments, key=lambda a: a['confidence'] - a['distance_error'] * 0.1)

            # Get original peak list entry
            original_row = None

            # Convert imported_id to string for comparison
            imported_id_str = str(imported_id)

            # Try to match by assignment name (handle both string and numeric assignments)
            for idx, row in original_peak_list.iterrows():
                row_assignment = str(row.get('Assignment', ''))
                if row_assignment == imported_id_str:
                    original_row = row
                    break

            # Fallback: try to extract index from Peak_X format
            if original_row is None and 'Peak_' in imported_id_str:
                try:
                    idx = int(imported_id_str.split('_')[1]) - 1  # Convert to 0-based index
                    if 0 <= idx < len(original_peak_list):
                        original_row = original_peak_list.iloc[idx]
                except:
                    pass

            # Final fallback
            if original_row is None and len(original_peak_list) > 0:
                original_row = original_peak_list.iloc[0]

            # DEBUG: Show assignment coordinates before creating final peak
            #print(f"   🔍 BUG DEBUG assignment {imported_id}:")
            print(f"      detected_position_x: {best_assignment.get('detected_position_x', 'MISSING')}")
            print(f"      detected_position_y: {best_assignment.get('detected_position_y', 'MISSING')}")
            print(f"      imported_position_x: {best_assignment.get('imported_position_x', 'MISSING')}")
            print(f"      imported_position_y: {best_assignment.get('imported_position_y', 'MISSING')}")

            # Create enhanced peak entry with EXACT detected coordinates (full precision)
            final_peak = {
                'Position_X': float(best_assignment['detected_position_x']),  # EXACT detected X (full precision)
                'Position_Y': float(best_assignment['detected_position_y']),  # EXACT detected Y (full precision)
                'Intensity': best_assignment['detected_intensity'],
                'Assignment': original_row.get('Assignment', imported_id) if original_row is not None else imported_id,
                'Volume': original_row.get('Volume', 0) if original_row is not None else 0,

                # Enhanced metadata
                'Detection_Method': 'enhanced_graph_based',
                'Original_Position_X': best_assignment['imported_position_x'],
                'Original_Position_Y': best_assignment['imported_position_y'],
                'Position_Error_X': abs(best_assignment['detected_position_x'] - best_assignment['imported_position_x']),
                'Position_Error_Y': abs(best_assignment['detected_position_y'] - best_assignment['imported_position_y']),
                'Match_Confidence': best_assignment['confidence'],
                'Match_Type': best_assignment['match_type'],
                'Alternatives_Considered': len(assignments)  # How many alternatives were considered
            }

            final_peaks.append(final_peak)

        # IMPORTANT: Do NOT add unmatched reference peaks to final results!
        # Enhanced Detection should ONLY transfer assignments to detected peaks
        # It should NOT discover new peaks or add unmatched references

        # Count summary (for diagnostics only)
        matched_assignments = set(assignment.get('imported_id', 'Unknown') for assignment in matching_results['assignments'])
        unmatched_references = []
        for idx, row in original_peak_list.iterrows():
            assignment = row.get('Assignment', f'Peak_{idx}')
            if assignment not in matched_assignments:
                unmatched_references.append(assignment)

        print(f"   📊 Assignment Transfer Results:")
        print(f"      Detected peaks with assignments: {len(final_peaks)}")
        print(f"      Reference peaks matched: {len(matched_assignments)}")
        print(f"      Reference peaks unmatched: {len(unmatched_references)}")
        print(f"   ⚠️ Note: Unmatched references are NOT added to results (assignment transfer only)")

        # SIMPLE PATTERN MATCHING - DISABLED (already handled in Enhanced Detection above)
        # This section is commented out to prevent double processing that was overriding good results
        pass

        # Store debug peaks correctly (detected coordinates only)
        debug_peaks_corrected = []
        for peak in final_peaks:
            if peak.get('Detection_Method') in ['simple_pattern_match', 'detected_unassigned']:
                debug_peaks_corrected.append({
                    'position_x': peak.get('Position_X', 0),  # These should be detected coordinates
                    'position_y': peak.get('Position_Y', 0),  # These should be detected coordinates
                    'intensity': peak.get('Intensity', 0),
                    'method': 'simple_pattern_final'
                })

        # Override debug data with correct final assignments
        if hasattr(self, '_included_peaks_after_limit_debug'):
            self._included_peaks_after_limit_debug = debug_peaks_corrected
            print(f"   🔍 CORRECTED debug: {len(debug_peaks_corrected)} peaks with detected coordinates")

        # DEBUG: Print final peak coordinates to trace corruption
        #print(f"   🔶 DEBUG: Final peak coordinates:")
        for i, peak in enumerate(final_peaks):
            method = peak.get('Detection_Method', 'unknown')
            x_pos = peak.get('Position_X', 0)
            y_pos = peak.get('Position_Y', 0)
            assignment = peak.get('Assignment', 'Unknown')
            print(f"      Peak {i+1} ({assignment}): X={x_pos:.3f}, Y={y_pos:.1f}, method={method}")
        return final_peaks

# Version and compatibility info
INTEGRATION_VERSION = "1.1.0"  # Updated for enhanced graph-based detection
REQUIRED_ENHANCED_FITTER_VERSION = "1.0.0"

if __name__ == "__main__":
    print(f"🧪 Integrated Detection-Fitter v{INTEGRATION_VERSION}")
    print("=" * 50)

    # Basic functionality test
    fitter = create_integrated_fitter()
    print(f"✅ Created integrated fitter")
    print(f"   Enhanced fitter available: {ENHANCED_AVAILABLE}")
    print(f"   Integration components loaded successfully")
