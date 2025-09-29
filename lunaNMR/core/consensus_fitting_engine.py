"""
Consensus Fitting Engine for Automated NMR Peak Analysis
========================================================

This module implements Priority 2 improvements including multi-peak consensus fitting,
separation of detection from fitting, and automated peak counting using AIC/BIC model selection.

Implements consensus approach for robust peak fitting adapted for NMR spectroscopy.

Author: Guillaume Mas
Date: 2025
"""

import numpy as np
import warnings
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.special import wofz
import time

@dataclass
class PeakCandidate:
    """Data class for peak candidates from detection phase"""
    position: float              # Peak position (ppm)
    intensity: float            # Peak intensity
    width_estimate: float       # Initial width estimate
    prominence: float           # Peak prominence
    confidence: float           # Detection confidence (0-1)
    method: str                 # Detection method used

@dataclass
class FittingResult:
    """Data class for individual fitting results"""
    peaks: List[Dict]           # Fitted peak parameters
    r_squared: float            # Goodness of fit
    aic: float                  # Akaike Information Criterion
    bic: float                  # Bayesian Information Criterion
    residual_std: float         # Residual standard deviation
    method: str                 # Fitting method used
    convergence_info: Dict      # Convergence information
    fit_time: float             # Fitting time in seconds

@dataclass
class ConsensusResult:
    """Data class for consensus fitting results"""
    best_fit: FittingResult     # Best individual fit
    consensus_peaks: List[Dict] # Consensus peak parameters
    model_selection: Dict       # Model selection results
    confidence_scores: Dict     # Confidence scores for parameters
    quality_assessment: Dict    # Overall quality metrics

class ConsensusDetector:
    """
    Peak detection engine that separates detection from fitting.

    Uses multiple detection methods to identify peak candidates that are
    then passed to the consensus fitting engine.
    """

    def __init__(self):
        """Initialize detector with multiple methods"""

        self.detection_methods = {
            'prominence': self._detect_by_prominence,
            'derivative': self._detect_by_derivative,
            'template': self._detect_by_template,
            'adaptive': self._detect_adaptive
        }

        # Detection parameters (internally managed)
        self.detection_params = {
            'min_prominence_factor': 0.1,    # Fraction of max intensity
            'min_distance_factor': 0.5,      # Minimum separation in typical linewidths
            'noise_threshold_factor': 2.0,   # SNR threshold for detection
            'smoothing_window': 3             # Smoothing for derivative methods
        }

    def detect_peaks(self, x_data: np.ndarray, y_data: np.ndarray,
                    nucleus_type: str = 'default',
                    detection_sensitivity: float = 0.5) -> List[PeakCandidate]:
        """
        Detect peak candidates using multiple methods.

        Parameters:
        -----------
        x_data : np.ndarray
            Frequency/chemical shift data
        y_data : np.ndarray
            Intensity data
        nucleus_type : str
            Type of nucleus for parameter adaptation
        detection_sensitivity : float
            Detection sensitivity (0-1, higher = more sensitive)

        Returns:
        --------
        List[PeakCandidate] : List of detected peak candidates
        """

        if len(x_data) == 0 or len(y_data) == 0:
            return []

        # Estimate data characteristics
        data_info = self._analyze_spectrum(x_data, y_data, nucleus_type)

        # Apply all detection methods
        all_candidates = []

        for method_name, method_func in self.detection_methods.items():
            try:
                candidates = method_func(x_data, y_data, data_info, detection_sensitivity)
                for candidate in candidates:
                    candidate.method = method_name
                all_candidates.extend(candidates)
            except Exception as e:
                warnings.warn(f"Detection method {method_name} failed: {e}")

        # Merge overlapping candidates
        merged_candidates = self._merge_candidates(all_candidates, data_info)

        # Sort by confidence
        merged_candidates.sort(key=lambda x: x.confidence, reverse=True)

        return merged_candidates

    def _analyze_spectrum(self, x_data: np.ndarray, y_data: np.ndarray,
                         nucleus_type: str) -> Dict:
        """Analyze spectrum characteristics for adaptive detection"""

        # Estimate noise level
        n_points = len(y_data)
        edge_fraction = 0.2
        edge_points = int(n_points * edge_fraction)

        noise_regions = np.concatenate([
            y_data[:edge_points],
            y_data[-edge_points:]
        ])

        noise_level = np.std(noise_regions)
        signal_level = np.max(y_data) - np.median(y_data)
        snr = signal_level / max(noise_level, 1e-10)

        # Estimate typical linewidth based on nucleus
        linewidth_estimates = {
            '1H': 0.01,    # ppm
            '15N': 0.5,    # ppm
            '13C': 0.6,    # ppm
            'default': 0.1  # ppm
        }

        typical_linewidth = linewidth_estimates.get(nucleus_type, 0.1)

        # Convert to data points
        ppm_range = abs(x_data[-1] - x_data[0])
        points_per_ppm = len(x_data) / ppm_range if ppm_range > 0 else 1
        typical_width_points = typical_linewidth * points_per_ppm

        return {
            'noise_level': noise_level,
            'signal_level': signal_level,
            'snr': snr,
            'typical_linewidth_ppm': typical_linewidth,
            'typical_width_points': typical_width_points,
            'points_per_ppm': points_per_ppm,
            'x_range': ppm_range
        }

    def _detect_by_prominence(self, x_data: np.ndarray, y_data: np.ndarray,
                            data_info: Dict, sensitivity: float) -> List[PeakCandidate]:
        """Detect peaks using prominence analysis"""

        # Adaptive prominence threshold
        min_prominence = data_info['signal_level'] * self.detection_params['min_prominence_factor']
        min_prominence *= (1.0 - sensitivity * 0.8)  # More sensitive -> lower threshold

        # Adaptive distance
        min_distance = max(1, int(data_info['typical_width_points'] *
                                self.detection_params['min_distance_factor']))

        # Find peaks
        peak_indices, properties = find_peaks(
            y_data,
            prominence=min_prominence,
            distance=min_distance,
            height=data_info['noise_level'] * self.detection_params['noise_threshold_factor']
        )

        candidates = []
        for i, peak_idx in enumerate(peak_indices):
            if 0 <= peak_idx < len(x_data):
                # Estimate width from prominence
                width_estimate = self._estimate_peak_width(
                    x_data, y_data, peak_idx, data_info['typical_linewidth_ppm']
                )

                # Calculate confidence based on prominence and SNR
                prominence = properties['prominences'][i] if i < len(properties['prominences']) else 0
                confidence = min(1.0, prominence / data_info['signal_level'] * 2.0)

                candidate = PeakCandidate(
                    position=x_data[peak_idx],
                    intensity=y_data[peak_idx],
                    width_estimate=width_estimate,
                    prominence=prominence,
                    confidence=confidence,
                    method='prominence'
                )
                candidates.append(candidate)

        return candidates

    def _detect_by_derivative(self, x_data: np.ndarray, y_data: np.ndarray,
                            data_info: Dict, sensitivity: float) -> List[PeakCandidate]:
        """Detect peaks using derivative analysis"""

        if len(y_data) < 5:
            return []

        # Smooth data slightly
        from scipy.ndimage import gaussian_filter1d
        smoothed = gaussian_filter1d(y_data, sigma=self.detection_params['smoothing_window'])

        # First derivative
        dx = np.mean(np.diff(x_data)) if len(x_data) > 1 else 1.0
        first_deriv = np.gradient(smoothed, dx)

        # Find zero crossings in first derivative (peak positions)
        zero_crossings = []
        for i in range(1, len(first_deriv)):
            if first_deriv[i-1] > 0 and first_deriv[i] <= 0:  # Positive to negative = peak
                zero_crossings.append(i)

        candidates = []
        for peak_idx in zero_crossings:
            if 0 < peak_idx < len(x_data) - 1:
                # Check if it's a significant peak
                intensity = y_data[peak_idx]
                if intensity > data_info['noise_level'] * 3.0:  # Minimum SNR

                    width_estimate = self._estimate_peak_width(
                        x_data, y_data, peak_idx, data_info['typical_linewidth_ppm']
                    )

                    # Confidence based on derivative magnitude and intensity
                    deriv_mag = abs(first_deriv[peak_idx-1] - first_deriv[peak_idx+1])
                    confidence = min(1.0, (intensity / data_info['signal_level']) *
                                         (deriv_mag / np.max(np.abs(first_deriv))))

                    candidate = PeakCandidate(
                        position=x_data[peak_idx],
                        intensity=intensity,
                        width_estimate=width_estimate,
                        prominence=intensity - np.median(y_data),
                        confidence=confidence,
                        method='derivative'
                    )
                    candidates.append(candidate)

        return candidates

    def _detect_by_template(self, x_data: np.ndarray, y_data: np.ndarray,
                          data_info: Dict, sensitivity: float) -> List[PeakCandidate]:
        """Detect peaks using template matching"""

        # Create Voigt template with typical linewidth
        template_width = data_info['typical_width_points']
        template_size = int(template_width * 6)  # 6x width for template

        if template_size < 5 or template_size > len(y_data) // 4:
            return []  # Template too large or too small

        # Generate Voigt template
        template_x = np.linspace(-3*template_width, 3*template_width, template_size)
        template_y = self._voigt_template(template_x, amplitude=1.0, center=0.0,
                                        sigma=template_width/4, gamma=template_width/4)

        # Normalize template
        template_y = template_y / np.max(template_y)

        # Cross-correlation
        from scipy.signal import correlate
        correlation = correlate(y_data, template_y, mode='same')

        # Find correlation peaks
        corr_threshold = np.max(correlation) * (1.0 - sensitivity) * 0.5 + 0.3
        corr_peaks, _ = find_peaks(correlation, height=corr_threshold)

        candidates = []
        for peak_idx in corr_peaks:
            if 0 <= peak_idx < len(x_data):
                intensity = y_data[peak_idx]
                if intensity > data_info['noise_level'] * 2.0:

                    width_estimate = data_info['typical_linewidth_ppm']
                    confidence = min(1.0, correlation[peak_idx] / np.max(correlation))

                    candidate = PeakCandidate(
                        position=x_data[peak_idx],
                        intensity=intensity,
                        width_estimate=width_estimate,
                        prominence=intensity - np.median(y_data),
                        confidence=confidence,
                        method='template'
                    )
                    candidates.append(candidate)

        return candidates

    def _detect_adaptive(self, x_data: np.ndarray, y_data: np.ndarray,
                        data_info: Dict, sensitivity: float) -> List[PeakCandidate]:
        """Adaptive detection that adjusts to data characteristics"""

        # Choose detection strategy based on SNR
        if data_info['snr'] > 10:
            # High SNR: use strict prominence detection
            return self._detect_by_prominence(x_data, y_data, data_info, sensitivity * 0.8)
        elif data_info['snr'] > 3:
            # Medium SNR: use template matching
            return self._detect_by_template(x_data, y_data, data_info, sensitivity)
        else:
            # Low SNR: use derivative with higher sensitivity
            return self._detect_by_derivative(x_data, y_data, data_info, min(1.0, sensitivity * 1.2))

    def _estimate_peak_width(self, x_data: np.ndarray, y_data: np.ndarray,
                           peak_idx: int, typical_width: float) -> float:
        """Estimate peak width from data"""

        if peak_idx < 1 or peak_idx >= len(y_data) - 1:
            return typical_width

        peak_intensity = y_data[peak_idx]
        half_height = peak_intensity / 2.0

        # Find half-width points
        left_idx = peak_idx
        right_idx = peak_idx

        # Search left
        while left_idx > 0 and y_data[left_idx] > half_height:
            left_idx -= 1

        # Search right
        while right_idx < len(y_data) - 1 and y_data[right_idx] > half_height:
            right_idx += 1

        # Calculate FWHM
        if right_idx > left_idx:
            fwhm = abs(x_data[right_idx] - x_data[left_idx])
            return max(fwhm, typical_width * 0.1)  # Minimum width
        else:
            return typical_width

    def _merge_candidates(self, candidates: List[PeakCandidate],
                         data_info: Dict) -> List[PeakCandidate]:
        """Merge overlapping peak candidates from different methods"""

        if not candidates:
            return []

        # Sort by position
        candidates.sort(key=lambda x: x.position)

        merged = []
        merge_distance = data_info['typical_linewidth_ppm'] * 0.5  # Merge within half linewidth

        for candidate in candidates:
            # Check if we should merge with existing candidate
            merged_with_existing = False

            for existing in merged:
                if abs(candidate.position - existing.position) < merge_distance:
                    # Merge: keep candidate with higher confidence
                    if candidate.confidence > existing.confidence:
                        # Update existing with better candidate
                        existing.position = candidate.position
                        existing.intensity = candidate.intensity
                        existing.width_estimate = candidate.width_estimate
                        existing.prominence = max(existing.prominence, candidate.prominence)
                        existing.confidence = candidate.confidence
                        existing.method = f"{existing.method}+{candidate.method}"
                    merged_with_existing = True
                    break

            if not merged_with_existing:
                merged.append(candidate)

        return merged

    @staticmethod
    def _voigt_template(x, amplitude, center, sigma, gamma):
        """Generate Voigt profile template"""
        try:
            sigma = max(sigma, 1e-6)
            z = ((x - center) + 1j*gamma) / (sigma * np.sqrt(2))
            voigt = amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2*np.pi))
            return voigt
        except:
            # Fallback to Gaussian
            return amplitude * np.exp(-0.5 * ((x - center) / max(sigma, 1e-6))**2)

class ConsensusFittingEngine:
    """
    Multi-method consensus fitting engine for robust peak analysis.

    Implements multiple fitting approaches and uses model selection criteria
    to determine optimal peak count and parameters.
    """

    def __init__(self):
        """Initialize consensus fitting engine"""

        self.fitting_methods = {
            'sequential': self._fit_sequential,
            'simultaneous': self._fit_simultaneous,
            'iterative': self._fit_iterative
        }

        # Model selection criteria weights
        self.model_selection_weights = {
            'aic': 0.4,      # Akaike Information Criterion
            'bic': 0.4,      # Bayesian Information Criterion
            'r_squared': 0.2  # Goodness of fit
        }

    def fit_consensus(self, x_data: np.ndarray, y_data: np.ndarray,
                     peak_candidates: List[PeakCandidate],
                     max_peaks: int = 10) -> ConsensusResult:
        """
        Perform consensus fitting with automated peak counting.

        Parameters:
        -----------
        x_data : np.ndarray
            Frequency/chemical shift data
        y_data : np.ndarray
            Intensity data
        peak_candidates : List[PeakCandidate]
            Peak candidates from detection phase
        max_peaks : int
            Maximum number of peaks to consider

        Returns:
        --------
        ConsensusResult : Consensus fitting results
        """

        if not peak_candidates or len(x_data) == 0:
            # Return empty result
            empty_fit = FittingResult(
                peaks=[], r_squared=0.0, aic=np.inf, bic=np.inf,
                residual_std=np.inf, method='none', convergence_info={}, fit_time=0.0
            )
            return ConsensusResult(
                best_fit=empty_fit,
                consensus_peaks=[],
                model_selection={},
                confidence_scores={},
                quality_assessment={}
            )

        # Sort candidates by confidence and limit number
        peak_candidates.sort(key=lambda x: x.confidence, reverse=True)
        peak_candidates = peak_candidates[:max_peaks]

        # Try different numbers of peaks using model selection
        model_results = {}

        for n_peaks in range(1, min(len(peak_candidates) + 1, max_peaks + 1)):
            # Use top n_peaks candidates
            selected_candidates = peak_candidates[:n_peaks]

            # Try all fitting methods for this peak count
            method_results = {}

            for method_name, method_func in self.fitting_methods.items():
                try:
                    start_time = time.time()
                    result = method_func(x_data, y_data, selected_candidates)
                    result.fit_time = time.time() - start_time
                    result.method = method_name

                    method_results[method_name] = result
                except Exception as e:
                    warnings.warn(f"Fitting method {method_name} failed for {n_peaks} peaks: {e}")

            if method_results:
                # Select best method for this peak count
                best_method_result = min(method_results.values(), key=lambda x: x.aic)
                model_results[n_peaks] = best_method_result

        # Model selection: find optimal peak count
        if not model_results:
            # Fallback: return empty result
            empty_fit = FittingResult(
                peaks=[], r_squared=0.0, aic=np.inf, bic=np.inf,
                residual_std=np.inf, method='failed', convergence_info={}, fit_time=0.0
            )
            return ConsensusResult(
                best_fit=empty_fit,
                consensus_peaks=[],
                model_selection={'error': 'All fitting methods failed'},
                confidence_scores={},
                quality_assessment={'status': 'failed'}
            )

        # Find optimal model using weighted criteria
        best_n_peaks = self._select_optimal_model(model_results)
        best_fit = model_results[best_n_peaks]

        # Generate consensus parameters by averaging stable parameters
        consensus_peaks = self._generate_consensus_peaks(model_results, best_n_peaks)

        # Calculate confidence scores
        confidence_scores = self._calculate_confidence_scores(model_results, best_n_peaks)

        # Quality assessment
        quality_assessment = self._assess_quality(best_fit, len(y_data))

        return ConsensusResult(
            best_fit=best_fit,
            consensus_peaks=consensus_peaks,
            model_selection={'optimal_peaks': best_n_peaks, 'models_tested': list(model_results.keys())},
            confidence_scores=confidence_scores,
            quality_assessment=quality_assessment
        )

    def _fit_sequential(self, x_data: np.ndarray, y_data: np.ndarray,
                       candidates: List[PeakCandidate]) -> FittingResult:
        """Fit peaks sequentially (one at a time)"""

        fitted_peaks = []
        residual = y_data.copy()

        for candidate in candidates:
            # Fit single peak to residual
            try:
                # Initial guess
                p0 = [candidate.intensity, candidate.position,
                      candidate.width_estimate/2, candidate.width_estimate/2, 0]

                # Bounds
                bounds = (
                    [0, candidate.position - candidate.width_estimate*2, 1e-6, 1e-6, -np.inf],
                    [np.inf, candidate.position + candidate.width_estimate*2, np.inf, np.inf, np.inf]
                )

                # Fit
                popt, pcov = curve_fit(self._voigt_1d, x_data, residual, p0=p0, bounds=bounds,
                                     maxfev=1000)

                # Calculate fitted curve
                fitted_curve = self._voigt_1d(x_data, *popt)

                # Update residual
                residual = residual - fitted_curve

                # Store peak parameters
                peak_params = {
                    'amplitude': popt[0],
                    'center': popt[1],
                    'sigma': popt[2],
                    'gamma': popt[3],
                    'baseline': popt[4],
                    'errors': np.sqrt(np.diag(pcov)) if pcov is not None else [0]*5
                }
                fitted_peaks.append(peak_params)

            except Exception as e:
                warnings.warn(f"Sequential fitting failed for peak at {candidate.position}: {e}")

        # Calculate overall statistics
        if fitted_peaks:
            total_fit = np.sum([self._voigt_1d(x_data, p['amplitude'], p['center'],
                                              p['sigma'], p['gamma'], 0) for p in fitted_peaks], axis=0)
            if len(fitted_peaks) > 0:
                total_fit += fitted_peaks[0]['baseline']  # Add baseline once

            r_squared = self._calculate_r_squared(y_data, total_fit)
            residual_std = np.std(y_data - total_fit)
            n_params = len(fitted_peaks) * 4 + 1  # 4 params per peak + baseline
            aic = self._calculate_aic(y_data, total_fit, n_params)
            bic = self._calculate_bic(y_data, total_fit, n_params)
        else:
            r_squared = 0.0
            residual_std = np.std(y_data)
            aic = np.inf
            bic = np.inf

        return FittingResult(
            peaks=fitted_peaks,
            r_squared=r_squared,
            aic=aic,
            bic=bic,
            residual_std=residual_std,
            method='sequential',
            convergence_info={'n_peaks_fitted': len(fitted_peaks)},
            fit_time=0.0
        )

    def _fit_simultaneous(self, x_data: np.ndarray, y_data: np.ndarray,
                         candidates: List[PeakCandidate]) -> FittingResult:
        """Fit all peaks simultaneously"""

        if not candidates:
            return FittingResult([], 0.0, np.inf, np.inf, np.inf, 'simultaneous', {}, 0.0)

        n_peaks = len(candidates)

        # Build initial guess and bounds
        p0 = []
        lower_bounds = []
        upper_bounds = []

        for candidate in candidates:
            # [amplitude, center, sigma, gamma] for each peak
            p0.extend([candidate.intensity, candidate.position,
                      candidate.width_estimate/2, candidate.width_estimate/2])

            lower_bounds.extend([0, candidate.position - candidate.width_estimate*2, 1e-6, 1e-6])
            upper_bounds.extend([np.inf, candidate.position + candidate.width_estimate*2,
                               np.inf, np.inf])

        # Add baseline
        p0.append(0)
        lower_bounds.append(-np.inf)
        upper_bounds.append(np.inf)

        bounds = (lower_bounds, upper_bounds)

        try:
            # Multi-peak Voigt function
            def multi_voigt(x, *params):
                baseline = params[-1]
                result = np.full_like(x, baseline)

                for i in range(n_peaks):
                    amp = params[i*4]
                    center = params[i*4 + 1]
                    sigma = params[i*4 + 2]
                    gamma = params[i*4 + 3]
                    result += self._voigt_1d(x, amp, center, sigma, gamma, 0)

                return result

            # Fit
            popt, pcov = curve_fit(multi_voigt, x_data, y_data, p0=p0, bounds=bounds,
                                 maxfev=2000)

            # Extract peak parameters
            fitted_peaks = []
            for i in range(n_peaks):
                peak_params = {
                    'amplitude': popt[i*4],
                    'center': popt[i*4 + 1],
                    'sigma': popt[i*4 + 2],
                    'gamma': popt[i*4 + 3],
                    'baseline': popt[-1],
                    'errors': np.sqrt(np.diag(pcov))[i*4:(i+1)*4].tolist() if pcov is not None else [0]*4
                }
                fitted_peaks.append(peak_params)

            # Calculate fit statistics
            fitted_curve = multi_voigt(x_data, *popt)
            r_squared = self._calculate_r_squared(y_data, fitted_curve)
            residual_std = np.std(y_data - fitted_curve)
            n_params = len(popt)
            aic = self._calculate_aic(y_data, fitted_curve, n_params)
            bic = self._calculate_bic(y_data, fitted_curve, n_params)

            return FittingResult(
                peaks=fitted_peaks,
                r_squared=r_squared,
                aic=aic,
                bic=bic,
                residual_std=residual_std,
                method='simultaneous',
                convergence_info={'converged': True},
                fit_time=0.0
            )

        except Exception as e:
            warnings.warn(f"Simultaneous fitting failed: {e}")
            return FittingResult([], 0.0, np.inf, np.inf, np.inf, 'simultaneous',
                               {'converged': False, 'error': str(e)}, 0.0)

    def _fit_iterative(self, x_data: np.ndarray, y_data: np.ndarray,
                      candidates: List[PeakCandidate]) -> FittingResult:
        """Fit peaks using iterative refinement"""

        if not candidates:
            return FittingResult([], 0.0, np.inf, np.inf, np.inf, 'iterative', {}, 0.0)

        # Start with sequential fit
        current_result = self._fit_sequential(x_data, y_data, candidates)

        # Iteratively refine
        max_iterations = 3
        for iteration in range(max_iterations):
            if not current_result.peaks:
                break

            try:
                # Use current result as starting point for simultaneous fit
                refined_candidates = []
                for peak in current_result.peaks:
                    refined_candidate = PeakCandidate(
                        position=peak['center'],
                        intensity=peak['amplitude'],
                        width_estimate=(peak['sigma'] + peak['gamma']) * 2,
                        prominence=peak['amplitude'],
                        confidence=1.0,
                        method='refined'
                    )
                    refined_candidates.append(refined_candidate)

                # Try simultaneous fit with refined parameters
                refined_result = self._fit_simultaneous(x_data, y_data, refined_candidates)

                # Accept if improvement
                if refined_result.aic < current_result.aic:
                    current_result = refined_result
                else:
                    break  # No improvement, stop iterating

            except Exception:
                break  # Refinement failed, keep current result

        current_result.method = 'iterative'
        return current_result

    def _select_optimal_model(self, model_results: Dict[int, FittingResult]) -> int:
        """Select optimal number of peaks using model selection criteria"""

        if not model_results:
            return 1

        # Calculate weighted scores
        scores = {}

        for n_peaks, result in model_results.items():
            # Normalize criteria (lower is better for AIC/BIC, higher for R²)
            aic_score = 1.0 / (1.0 + result.aic)  # Invert AIC
            bic_score = 1.0 / (1.0 + result.bic)  # Invert BIC
            r2_score = result.r_squared

            # Weight and combine
            weighted_score = (
                self.model_selection_weights['aic'] * aic_score +
                self.model_selection_weights['bic'] * bic_score +
                self.model_selection_weights['r_squared'] * r2_score
            )

            scores[n_peaks] = weighted_score

        # Return n_peaks with highest score
        return max(scores, key=scores.get)

    def _generate_consensus_peaks(self, model_results: Dict[int, FittingResult],
                                optimal_n_peaks: int) -> List[Dict]:
        """Generate consensus peak parameters from multiple fits"""

        optimal_result = model_results[optimal_n_peaks]
        return optimal_result.peaks  # For now, just return optimal result

        # TODO: Implement true consensus by averaging compatible peaks across models

    def _calculate_confidence_scores(self, model_results: Dict[int, FittingResult],
                                   optimal_n_peaks: int) -> Dict:
        """Calculate confidence scores for fitted parameters"""

        optimal_result = model_results[optimal_n_peaks]

        # Simple confidence based on fit quality and parameter errors
        confidence_scores = {
            'overall': min(1.0, optimal_result.r_squared),
            'peak_confidence': []
        }

        for i, peak in enumerate(optimal_result.peaks):
            # Confidence based on parameter errors and amplitude
            if 'errors' in peak and len(peak['errors']) >= 4:
                amp_error = peak['errors'][0] if peak['errors'][0] > 0 else 0.1
                relative_error = amp_error / max(peak['amplitude'], 1e-10)
                peak_confidence = max(0.0, 1.0 - relative_error)
            else:
                peak_confidence = 0.8  # Default moderate confidence

            confidence_scores['peak_confidence'].append(peak_confidence)

        return confidence_scores

    def _assess_quality(self, fit_result: FittingResult, data_length: int) -> Dict:
        """Assess overall quality of the fit"""

        if fit_result.r_squared >= 0.9:
            quality = 'excellent'
        elif fit_result.r_squared >= 0.8:
            quality = 'good'
        elif fit_result.r_squared >= 0.6:
            quality = 'fair'
        else:
            quality = 'poor'

        return {
            'overall_quality': quality,
            'r_squared': fit_result.r_squared,
            'residual_std': fit_result.residual_std,
            'n_peaks': len(fit_result.peaks),
            'n_data_points': data_length,
            'degrees_of_freedom': data_length - len(fit_result.peaks) * 4 - 1
        }

    @staticmethod
    def _voigt_1d(x, amplitude, center, sigma, gamma, baseline=0):
        """1D Voigt profile"""
        try:
            sigma = max(sigma, 1e-6)
            z = ((x - center) + 1j*gamma) / (sigma * np.sqrt(2))
            voigt = amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2*np.pi))
            return voigt + baseline
        except:
            return amplitude * np.exp(-0.5 * ((x - center) / max(sigma, 1e-6))**2) + baseline

    @staticmethod
    def _calculate_r_squared(y_true, y_pred):
        """Calculate R² coefficient of determination"""
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / max(ss_tot, 1e-10))

    @staticmethod
    def _calculate_aic(y_true, y_pred, n_params):
        """Calculate Akaike Information Criterion"""
        n = len(y_true)
        mse = np.mean((y_true - y_pred) ** 2)
        aic = n * np.log(max(mse, 1e-10)) + 2 * n_params
        return aic

    @staticmethod
    def _calculate_bic(y_true, y_pred, n_params):
        """Calculate Bayesian Information Criterion"""
        n = len(y_true)
        mse = np.mean((y_true - y_pred) ** 2)
        bic = n * np.log(max(mse, 1e-10)) + n_params * np.log(n)
        return bic