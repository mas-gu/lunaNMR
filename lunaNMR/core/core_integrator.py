#!/usr/bin/env python3
"""
Core NMR Integration Module

This module contains the core technical functions for NMR peak detection,
integration, and Voigt profile fitting. It provides both standard and
enhanced integration capabilities.

Classes:
- VoigtIntegrator: Basic Voigt fitting integrator
- EnhancedVoigtIntegrator: Advanced integrator with enhanced features

Author: Guillaume Mas
Date: 2025
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Any
from scipy.special import wofz  # Faddeeva function for Voigt profiles
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

# Import enhanced Voigt fitter
try:
    from lunaNMR.core.enhanced_voigt_fitter import EnhancedVoigtFitter
    ENHANCED_FITTING_AVAILABLE = True
    print("✅ Enhanced Voigt fitting available")
except ImportError:
    ENHANCED_FITTING_AVAILABLE = False
    print("⚠️ Enhanced Voigt fitting not available, using standard fitting")

# Import integrated detection-fitter (NEW INTEGRATION)
try:
    from lunaNMR.core.integrated_detection_fitter import (
        IterativeDetectionFitter,
        DetectionQualityScorer,
        RapidFitAssessor,
        AdaptiveThresholdCalculator,
        ChemicalShiftContextAnalyzer,
        create_integrated_fitter
    )
    INTEGRATED_DETECTION_AVAILABLE = True
    print("🚀 Integrated detection-fitting available")
except ImportError as e:
    INTEGRATED_DETECTION_AVAILABLE = False
    print(f"⚠️ Integrated detection-fitting not available: {e}")

try:
    from lunaNMR.integrators.inplace_advanced_nmr_integrator import InPlaceAdvancedNMRIntegrator
    from lunaNMR.integrators.inplace_series_nmr_integrator import InPlaceSeriesNMRIntegrator
    BaseIntegrator = InPlaceAdvancedNMRIntegrator
except ImportError as e:
    print(f"Import error: {e}")
    print("Creating fallback base integrator class...")

    # Fallback base class with essential methods
    class BaseIntegrator:
        def __init__(self):
            self.nmr_data = None
            self.nmr_dict = None
            self.ppm_x_axis = None
            self.ppm_y_axis = None
            self.peak_list = None
            self.nmr_file_path = None
            self.peak_list_path = None

        def _calculate_ppm_axes(self):
            """Calculate PPM axes from NMR dictionary"""
            if self.nmr_dict is None or self.nmr_data is None:
                return
            try:
                import nmrglue as ng
                # Calculate ppm axes using nmrglue
                uc_x = ng.pipe.make_uc(self.nmr_dict, self.nmr_data, dim=1)
                uc_y = ng.pipe.make_uc(self.nmr_dict, self.nmr_data, dim=0)

                self.ppm_x_axis = uc_x.ppm_scale()
                self.ppm_y_axis = uc_y.ppm_scale()

                print(f"PPM axes calculated - X: {self.ppm_x_axis[0]:.2f} to {self.ppm_x_axis[-1]:.2f} ppm")
                print(f"                      Y: {self.ppm_y_axis[0]:.1f} to {self.ppm_y_axis[-1]:.1f} ppm")
            except Exception as e:
                print(f"Error calculating PPM axes: {e}")

        def _estimate_noise_level(self):
            """Estimate noise level from spectrum corners"""
            if self.nmr_data is None:
                return
            try:
                # Sample corners of the spectrum for noise estimation
                h, w = self.nmr_data.shape
                corner_size = min(20, h//10, w//10)

                corners = [
                    self.nmr_data[:corner_size, :corner_size],
                    self.nmr_data[:corner_size, -corner_size:],
                    self.nmr_data[-corner_size:, :corner_size],
                    self.nmr_data[-corner_size:, -corner_size:]
                ]

                noise_data = np.concatenate([corner.flatten() for corner in corners])
                self.noise_level = np.std(noise_data)
                print(f"Estimated noise level: {self.noise_level:.2e}")
            except Exception as e:
                print(f"Error estimating noise: {e}")
                self.noise_level = 1e6  # Default fallback

        def load_data(self, peak_list_file, nmr_file):
            """Load both peak list and NMR data"""
            success = True
            success &= self.load_peak_list_file(peak_list_file)
            success &= self._load_nmr_data_only(nmr_file)
            return success

        def detect_peaks_professional(self, **kwargs):
            """Fallback peak detection method"""
            print("Using fallback peak detection - limited functionality")
            return []

class VoigtIntegrator(BaseIntegrator):
    """Enhanced NMR integrator with Voigt profile fitting capabilities"""

    def __init__(self):
        super().__init__()
        self.voigt_fits = []  # Store detailed Voigt fitting results
        self.fitting_parameters = {
            'max_iterations': 1000,
            'tolerance': 1e-6,
            'min_r_squared': 0.7,
            'fitting_window_x': 0.2,  # ppm @GM was 0.4 was 0.4 was 0.15
            'fitting_window_y': 2   # ppm  @GM was 7 was 5 was 3 was 1.5
        }
        self.processing_mode = 'full_detection'  # 'full_detection' or 'in_place'

        # Initialize enhanced fitter if available
        if ENHANCED_FITTING_AVAILABLE:
            self.enhanced_fitter = EnhancedVoigtFitter()
            self.enhanced_fitter.parent = self  # Set parent reference for parallel processing
            print("✅ Enhanced Voigt fitter initialized")
        else:
            self.enhanced_fitter = None

    def load_nmr_data(self, data_2d, ppm_x_axis, ppm_y_axis):
        """
        TESTING SUPPORT: Load NMR data directly from arrays

        This method is specifically for testing and direct data input.
        For file-based loading, use load_data() or load_nmr_file().

        Parameters:
        - data_2d: 2D numpy array with NMR data
        - ppm_x_axis: 1D array with X-axis ppm values
        - ppm_y_axis: 1D array with Y-axis ppm values
        """
        self.nmr_data = data_2d
        self.ppm_x_axis = ppm_x_axis
        self.ppm_y_axis = ppm_y_axis

        print(f"✅ Loaded NMR data directly: {data_2d.shape}")
        print(f"   X-axis: {ppm_x_axis[0]:.2f} to {ppm_x_axis[-1]:.2f} ppm")
        print(f"   Y-axis: {ppm_y_axis[0]:.1f} to {ppm_y_axis[-1]:.1f} ppm")

        # Set dummy nmr_dict for compatibility
        self.nmr_dict = {'dummy': 'for_testing'}

        return True

    def set_processing_mode(self, mode):
        """Set processing mode: 'full_detection' or 'in_place'"""
        if mode in ['full_detection', 'in_place']:
            self.processing_mode = mode
        else:
            raise ValueError("Mode must be 'full_detection' or 'in_place'")

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
        # Compute complex argument for Faddeeva function
        z = ((x - center) + 1j*gamma) / (sigma * np.sqrt(2))

        # Voigt profile using Faddeeva function
        voigt = amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2*np.pi))

        return voigt + baseline

    @staticmethod
    def gaussian_profile_1d(x, amplitude, center, sigma, baseline=0):
        """1D Gaussian profile for comparison"""
        return amplitude * np.exp(-0.5 * ((x - center) / sigma)**2) + baseline

    @staticmethod
    def lorentzian_profile_1d(x, amplitude, center, gamma, baseline=0):
        """1D Lorentzian profile for comparison"""
        return amplitude * gamma / (np.pi * (gamma**2 + (x - center)**2)) + baseline

    def _estimate_peak_width(self, ppm_scale, cross_section, center_ppm):
        """Estimate peak width by finding FWHM"""
        try:
            # Find the peak index
            center_idx = np.argmin(np.abs(ppm_scale - center_ppm))
            peak_intensity = cross_section[center_idx]

            # Find baseline (average of edges)
            baseline = (np.mean(cross_section[:3]) + np.mean(cross_section[-3:])) / 2

            # Half maximum intensity
            half_max = baseline + (peak_intensity - baseline) / 2

            # Find FWHM indices
            left_idx = center_idx
            right_idx = center_idx

            # Search left side
            while left_idx > 0 and cross_section[left_idx] > half_max:
                left_idx -= 1

            # Search right side
            while right_idx < len(cross_section) - 1 and cross_section[right_idx] > half_max:
                right_idx += 1

            # Calculate FWHM in ppm
            if right_idx > left_idx:
                fwhm = abs(ppm_scale[right_idx] - ppm_scale[left_idx])
                return max(fwhm, 0.01)  # Minimum width constraint
            else:
                # Fallback to reasonable defaults
                return 0.05 if len(ppm_scale) < 50 else 2.0

        except:
            # Fallback width estimation
            return 0.05 if len(ppm_scale) < 50 else 2.0

    def _get_fitting_bounds(self, initial_guess, ppm_scale, dimension):
        """Get reasonable bounds for fitting parameters"""
        amplitude, center, sigma, gamma, baseline = initial_guess
        ppm_range = abs(ppm_scale[-1] - ppm_scale[0])

        if dimension == 'x':
            # X-dimension bounds (narrower peaks)
            lower_bounds = [
                0,                          # amplitude > 0
                center - ppm_range * 0.1,   # center can shift slightly
                0.001,                      # minimum sigma
                0.001,                      # minimum gamma
                baseline - abs(amplitude) * 0.5  # baseline can vary
            ]
            upper_bounds = [
                amplitude * 3,              # amplitude upper limit
                center + ppm_range * 0.1,   # center can shift slightly
                ppm_range * 0.2,           # maximum sigma
                ppm_range * 0.2,           # maximum gamma
                baseline + abs(amplitude) * 0.5  # baseline can vary
            ]
        else:
            # Y-dimension bounds (broader peaks)
            lower_bounds = [
                0,                          # amplitude > 0
                center - ppm_range * 0.2,   # center can shift more in Y
                0.1,                        # minimum sigma (larger than X)
                0.1,                        # minimum gamma (larger than X)
                baseline - abs(amplitude) * 0.8  # baseline can vary more
            ]
            upper_bounds = [
                amplitude * 3,              # amplitude upper limit
                center + ppm_range * 0.2,   # center can shift more in Y
                ppm_range * 0.5,           # maximum sigma (broader than X)
                ppm_range * 0.3,           # maximum gamma
                baseline + abs(amplitude) * 0.8  # baseline can vary more
            ]

        return (lower_bounds, upper_bounds)

    def fit_peak_1d(self, x_data, y_data, initial_guess, profile_type='voigt', bounds=None):
        """
        Fit 1D peak profile to data

        Parameters:
        - x_data: x-axis data (ppm)
        - y_data: intensity data
        - initial_guess: [amplitude, center, width1, width2, baseline]
        - profile_type: 'voigt', 'gaussian', or 'lorentzian'
        """
        try:
            if profile_type == 'voigt':
                # Voigt: amplitude, center, sigma, gamma, baseline
                if bounds is not None:
                    popt, pcov = curve_fit(self.voigt_profile_1d, x_data, y_data,
                                         p0=initial_guess, bounds=bounds,
                                         maxfev=self.fitting_parameters['max_iterations'])
                else:
                    popt, pcov = curve_fit(self.voigt_profile_1d, x_data, y_data, p0=initial_guess,
                                         maxfev=self.fitting_parameters['max_iterations'])

                # Calculate R-squared and quality metrics
                y_pred = self.voigt_profile_1d(x_data, *popt)
                r_squared = self.calculate_r_squared(y_data, y_pred)

                # Add local quality assessment using enhanced fitter if available
                quality_metrics = {'r_squared': r_squared, 'r_squared_global': r_squared}

                if self.enhanced_fitter is not None:
                    try:
                        # Extract parameters for local assessment
                        amplitude, center, sigma, gamma, baseline = popt
                        fitted_width = sigma + gamma

                        # Use enhanced fitter's local quality assessment
                        nucleus_type = self.enhanced_fitter.detect_nucleus_type([x_data[0], x_data[-1]])
                        local_region = self.enhanced_fitter.extract_local_peak_region(
                            x_data, y_data, center, fitted_width, nucleus_type
                        )

                        # Get local data and predictions
                        x_local = local_region['x_local']
                        y_local = local_region['y_local']
                        y_pred_local = y_pred[local_region['local_indices']]

                        # Calculate local baseline and R-squared
                        local_baseline = self.enhanced_fitter.calculate_local_baseline(x_local, y_local)
                        ss_res_local = np.sum((y_local - y_pred_local) ** 2)
                        ss_tot_local = np.sum((y_local - local_baseline) ** 2)
                        r_squared_local = 1 - (ss_res_local / ss_tot_local) if ss_tot_local != 0 else 0

                        # Update quality metrics
                        quality_metrics.update({
                            'r_squared_local': r_squared_local,
                            'r_squared': r_squared_local,  # Use local as primary
                            'peak_region_points': local_region['n_points'],
                            'peak_region_width': local_region['window_width']
                        })

                    except Exception as e:
                        print(f"Local quality assessment failed, using global: {e}")
                        quality_metrics['r_squared_local'] = r_squared
                else:
                    quality_metrics['r_squared_local'] = r_squared

                return {
                    'success': True,
                    'parameters': popt,
                    'covariance': pcov,
                    'r_squared': quality_metrics['r_squared'],
                    'fitted_curve': y_pred,
                    'profile_type': 'voigt',
                    'amplitude': popt[0],
                    'center': popt[1],
                    'sigma': popt[2],
                    'gamma': popt[3],
                    'baseline': popt[4],
                    'quality_metrics': quality_metrics
                }

            elif profile_type == 'gaussian':
                # Gaussian: amplitude, center, sigma, baseline
                popt, pcov = curve_fit(self.gaussian_profile_1d, x_data, y_data,
                                     p0=initial_guess[:4], maxfev=self.fitting_parameters['max_iterations'])
                y_pred = self.gaussian_profile_1d(x_data, *popt)
                r_squared = self.calculate_r_squared(y_data, y_pred)

                # Add local quality assessment for Gaussian fits too
                quality_metrics = {'r_squared': r_squared, 'r_squared_global': r_squared}

                if self.enhanced_fitter is not None:
                    try:
                        # Extract parameters for local assessment (Gaussian: amp, center, sigma, baseline)
                        amplitude, center, sigma, baseline = popt
                        fitted_width = sigma * 2.355  # Convert sigma to FWHM equivalent

                        # Use enhanced fitter's local quality assessment
                        nucleus_type = self.enhanced_fitter.detect_nucleus_type([x_data[0], x_data[-1]])
                        local_region = self.enhanced_fitter.extract_local_peak_region(
                            x_data, y_data, center, fitted_width, nucleus_type
                        )

                        # Get local data and predictions
                        x_local = local_region['x_local']
                        y_local = local_region['y_local']
                        y_pred_local = y_pred[local_region['local_indices']]

                        # Calculate local baseline and R-squared
                        local_baseline = self.enhanced_fitter.calculate_local_baseline(x_local, y_local)
                        ss_res_local = np.sum((y_local - y_pred_local) ** 2)
                        ss_tot_local = np.sum((y_local - local_baseline) ** 2)
                        r_squared_local = 1 - (ss_res_local / ss_tot_local) if ss_tot_local != 0 else 0

                        # Update quality metrics
                        quality_metrics.update({
                            'r_squared_local': r_squared_local,
                            'r_squared': r_squared_local,  # Use local as primary
                            'peak_region_points': local_region['n_points'],
                            'peak_region_width': local_region['window_width']
                        })

                    except Exception as e:
                        print(f"Local quality assessment failed, using global: {e}")
                        quality_metrics['r_squared_local'] = r_squared
                else:
                    quality_metrics['r_squared_local'] = r_squared

                return {
                    'success': True,
                    'parameters': popt,
                    'covariance': pcov,
                    'r_squared': quality_metrics['r_squared'],
                    'fitted_curve': y_pred,
                    'profile_type': 'gaussian',
                    'amplitude': popt[0],
                    'center': popt[1],
                    'sigma': popt[2],
                    'baseline': popt[3],
                    'quality_metrics': quality_metrics
                }

        except Exception as e:
            print(f"1D fitting failed: {e}")
            return {'success': False, 'error': str(e)}

    def fit_peak_1d_enhanced(self, x_data, y_data, initial_center=None, nucleus_type=None):
        """
        Enhanced 1D peak fitting using the improved fitter

        Falls back to standard fitting if enhanced fitter not available
        """
        if self.enhanced_fitter is not None:
            try:
                # Use enhanced fitter
                result = self.enhanced_fitter.fit_peak_enhanced(
                    x_data, y_data, initial_center, nucleus_type
                )
                return result
            except Exception as e:
                print(f"Enhanced fitting failed, falling back to standard: {e}")

        # Fallback to standard fitting
        if initial_center is None:
            initial_center = x_data[np.argmax(y_data)]

        # Standard parameter estimation
        baseline = (np.mean(y_data[:5]) + np.mean(y_data[-5:])) / 2
        amplitude = np.max(y_data) - baseline

        # Use nucleus-appropriate width estimates (CORRECTED VALUES)
        if nucleus_type == '15N':
            width_est = 1.5  # ¹⁵N typical width
        elif nucleus_type == '13C':
            width_est = 1.0  # ¹³C typical width
        else:  # Default to ¹H
            width_est = 0.02  # ¹H typical width (5.5-12 ppm range)

        initial_guess = [amplitude, initial_center, width_est * 0.7, width_est * 0.3, baseline]

        # Get bounds for standard fitting
        bounds = self._get_fitting_bounds(initial_guess, x_data, 'x' if nucleus_type == '1H' else 'y')

        return self.fit_peak_1d(x_data, y_data, initial_guess, 'voigt', bounds)

    def _calculate_initial_guess_1d(self, x_data, y_data, target_position):
        """Calculate initial guess for 1D Voigt fitting"""
        # Baseline estimation using edge points
        baseline = np.mean([np.mean(y_data[:3]), np.mean(y_data[-3:])])

        # Amplitude estimation
        peak_intensity = np.max(y_data)
        amplitude = abs(peak_intensity - baseline) * 1.2  # 20% buffer

        # Width estimation
        width_estimate = self._estimate_peak_width(x_data, y_data, target_position)

        # Voigt parameters: [amplitude, center, sigma, gamma, baseline]
        return [
            amplitude,                    # amplitude
            target_position,             # center
            width_estimate * 0.6,        # sigma (Gaussian component)
            width_estimate * 0.4,        # gamma (Lorentzian component)
            baseline                     # baseline
        ]

    # =================== ADAPTIVE MULTI-PEAK FITTING ===================

    def detect_peaks_1d(self, x_data, y_data, target_position=None, gui_params=None):
        """Detect peaks in 1D cross-section using scipy peak detection optimized for overlapping peaks"""
        from scipy.signal import find_peaks, peak_widths, peak_prominences

        # Get parameters from GUI if available, otherwise use defaults
        if gui_params:
            height_threshold = gui_params.get('height_threshold', 0.02)
            distance_factor = gui_params.get('distance_factor', 50.0)
            prominence_threshold = gui_params.get('prominence_threshold', 0.01)
            smoothing_sigma = gui_params.get('smoothing_sigma', 0.5)
        else:
            # Default values
            height_threshold = 0.02
            distance_factor = 50.0
            prominence_threshold = 0.01
            smoothing_sigma = 0.5

        # Light smoothing to reduce noise but preserve peak separation
        from scipy.ndimage import gaussian_filter1d
        y_smooth = gaussian_filter1d(y_data, sigma=smoothing_sigma)

        # Calculate noise level for adaptive thresholds
        noise_level = np.std(y_smooth[:10] if len(y_smooth) > 10 else y_smooth)
        signal_max = np.max(y_smooth)

        # Adaptive parameters for overlapping peak detection
        min_height = max(signal_max * height_threshold, noise_level * 3)
        min_distance = max(2, int(len(y_smooth) / distance_factor))
        min_prominence = max(signal_max * prominence_threshold, noise_level * 2)

        # Primary peak detection with relaxed criteria
        peaks, properties = find_peaks(y_smooth,
                                     height=min_height,
                                     distance=min_distance,
                                     prominence=min_prominence)

        print(f"   Initial detection: {len(peaks)} peaks (height≥{min_height:.0f}, dist≥{min_distance}, prom≥{min_prominence:.0f})")

        # Secondary detection for overlapping peaks (more sensitive)
        if len(peaks) <= 3:  # If few peaks detected, try more sensitive detection
            min_height_sensitive = max(signal_max * 0.01, noise_level * 2)  # Even lower threshold
            min_distance_sensitive = max(1, len(y_smooth) // 100)           # Very close peaks allowed

            peaks_sensitive, _ = find_peaks(y_smooth,
                                          height=min_height_sensitive,
                                          distance=min_distance_sensitive,
                                          prominence=noise_level)

            if len(peaks_sensitive) > len(peaks):
                peaks = peaks_sensitive
                print(f"   Sensitive detection: {len(peaks)} peaks (improved overlap detection)")

        # Filter out peaks too close to edges (boundary artifacts)
        edge_buffer = max(2, len(y_smooth) // 100)
        valid_peaks = peaks[(peaks >= edge_buffer) & (peaks < len(y_smooth) - edge_buffer)]

        if len(valid_peaks) != len(peaks):
            print(f"   Edge filtering: {len(valid_peaks)} peaks (removed {len(peaks) - len(valid_peaks)} edge artifacts)")

        if len(valid_peaks) == 0:
            return []

        # Get peak positions in ppm
        peak_positions = x_data[valid_peaks]
        peak_intensities = y_data[valid_peaks]

        # Sort by intensity (strongest first)
        sorted_indices = np.argsort(peak_intensities)[::-1]

        peak_info = []
        for idx in sorted_indices:
            peak_info.append({
                'index': valid_peaks[idx],
                'position': peak_positions[idx],
                'intensity': peak_intensities[idx],
                'distance_to_target': abs(peak_positions[idx] - target_position) if target_position else 0
            })

        return peak_info

    def fit_multi_peak_1d(self, x_data, y_data, peak_info, max_peaks=4):
        """Fit multiple Voigt profiles simultaneously"""
        if len(peak_info) == 0:
            return None

        # Limit number of peaks to avoid overfitting
        n_peaks = min(len(peak_info), max_peaks)
        peaks_to_fit = peak_info[:n_peaks]

        # Create multi-Voigt model
        def multi_voigt(x, *params):
            """Multi-peak Voigt model: amplitude, center, sigma, gamma for each peak + baseline"""
            n_peaks = (len(params) - 1) // 4  # -1 for baseline
            baseline = params[-1]

            result = np.full_like(x, baseline)
            for i in range(n_peaks):
                A, x0, sigma, gamma = params[i*4:(i+1)*4]
                # Add individual Voigt contribution
                z = ((x - x0) + 1j * gamma) / (sigma * np.sqrt(2))
                voigt = A * np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))
                result += voigt

            return result

        # Prepare initial guess
        baseline_est = np.mean([np.mean(y_data[:5]), np.mean(y_data[-5:])])
        initial_guess = []

        for peak in peaks_to_fit:
            # Estimate parameters for each peak
            amplitude = peak['intensity'] - baseline_est
            center = peak['position']
            # Estimate width from local data
            width_est = abs(x_data[1] - x_data[0]) * 5  # 5 data points width
            sigma = width_est * 0.6  # Gaussian component
            gamma = width_est * 0.4  # Lorentzian component

            initial_guess.extend([amplitude, center, sigma, gamma])

        initial_guess.append(baseline_est)

        # Set bounds
        bounds_lower = []
        bounds_upper = []

        for peak in peaks_to_fit:
            # Bounds for [amplitude, center, sigma, gamma]
            bounds_lower.extend([0, peak['position'] - abs(x_data[-1] - x_data[0]) * 0.1,
                               abs(x_data[1] - x_data[0]) * 0.5, abs(x_data[1] - x_data[0]) * 0.1])
            bounds_upper.extend([peak['intensity'] * 3, peak['position'] + abs(x_data[-1] - x_data[0]) * 0.1,
                               abs(x_data[-1] - x_data[0]) * 0.5, abs(x_data[-1] - x_data[0]) * 0.5])

        # Baseline bounds
        bounds_lower.append(baseline_est - abs(baseline_est))
        bounds_upper.append(baseline_est + abs(baseline_est))

        try:
            popt, pcov = curve_fit(multi_voigt, x_data, y_data,
                                 p0=initial_guess,
                                 bounds=(bounds_lower, bounds_upper),
                                 maxfev=1000)  # Reduced to prevent timeout

            # Calculate fit quality
            y_pred = multi_voigt(x_data, *popt)
            r_squared = self.calculate_r_squared(y_data, y_pred)

            # Extract parameters for each peak
            peaks_fitted = []
            for i in range(n_peaks):
                A, x0, sigma, gamma = popt[i*4:(i+1)*4]
                peaks_fitted.append({
                    'amplitude': A,
                    'center': x0,
                    'sigma': sigma,
                    'gamma': gamma,
                    'original_peak': peaks_to_fit[i]
                })

            return {
                'success': True,
                'n_peaks': n_peaks,
                'peaks': peaks_fitted,
                'baseline': popt[-1],
                'fitted_curve': y_pred,
                'r_squared': r_squared,
                'parameters': popt,
                'covariance': pcov
            }

        except Exception as e:
            if "maximum number of function evaluations" in str(e).lower():
                print(f"Multi-peak fitting timeout (maxfev exceeded) - continuing with best available fit")
            else:
                print(f"Multi-peak fitting failed: {e}")
            return None

    def select_target_peak(self, fit_result, target_position):
        """Select the peak closest to target position from multi-peak fit"""
        if not fit_result or 'peaks' not in fit_result:
            return None

        best_peak = None
        min_distance = float('inf')

        for peak in fit_result['peaks']:
            distance = abs(peak['center'] - target_position)
            if distance < min_distance:
                min_distance = distance
                best_peak = peak

        return best_peak

    def _try_gaussian_mixture_model(self, x_data, y_data, target_position, n_peaks_hint=None):
        """Try Gaussian Mixture Model for very close overlapping peaks"""
        try:
            from lunaNMR.processors.parallel_fitting import GaussianMixtureModel

            # Use GMM for close peak fitting
            gmm = GaussianMixtureModel(max_iter=100)
            gmm_result = gmm.fit_overlapping_peaks(x_data, y_data, n_peaks_hint)

            if gmm_result and gmm_result.get('success', False):
                # Find peak closest to target position
                best_peak = None
                min_distance = float('inf')

                for peak in gmm_result['peaks']:
                    distance = abs(peak['center'] - target_position)
                    if distance < min_distance:
                        min_distance = distance
                        best_peak = peak

                if best_peak:
                    return {
                        'success': True,
                        'r_squared': gmm_result['r_squared'],
                        'fitted_curve': gmm_result['fitted_curve'],
                        'amplitude': best_peak['amplitude'],
                        'center': best_peak['center'],
                        'sigma': best_peak['sigma'],
                        'gamma': best_peak['gamma'],
                        'baseline': 0,  # GMM doesn't separate baseline
                        'gmm_info': gmm_result
                    }

            return None

        except ImportError:
            print("   ⚠️ Gaussian Mixture Model not available (sklearn required)")
            return None
        except Exception as e:
            print(f"   GMM fitting failed: {e}")
            return None

    def iterative_parameter_optimization(self, x_data, y_data, target_position, dimension='x', gui_params=None):
        """Iterative parameter optimization for failed peak fits"""
        print(f"   🔄 Starting iterative parameter optimization for {dimension}-dimension")

        # Define parameter exploration space
        height_thresholds = [0.005, 0.01, 0.02, 0.05, 0.1]
        distance_factors = [20, 35, 50, 75, 100]
        prominence_thresholds = [0.005, 0.01, 0.02, 0.05]
        smoothing_sigmas = [0.3, 0.5, 0.8]

        best_fit = None
        best_r_squared = -float('inf')
        best_params = None
        attempts = 0
        max_attempts = gui_params.get('max_optimization_iterations', 50) if gui_params else 50

        # Start with more sensitive parameters (lower thresholds first)
        for height_thresh in height_thresholds:
            for distance_factor in distance_factors:
                for prominence_thresh in prominence_thresholds:
                    for smoothing_sigma in smoothing_sigmas:
                        if attempts >= max_attempts:
                            break

                        attempts += 1

                        # Create parameter set
                        test_params = {
                            'height_threshold': height_thresh,
                            'distance_factor': distance_factor,
                            'prominence_threshold': prominence_thresh,
                            'smoothing_sigma': smoothing_sigma,
                            'max_peaks_fit': gui_params.get('max_peaks_fit', 4) if gui_params else 4
                        }

                        try:
                            # Test peak detection
                            detected_peaks = self.detect_peaks_1d(x_data, y_data, target_position, test_params)
                            n_detected = len(detected_peaks)

                            if n_detected == 0:
                                continue  # No peaks detected, skip

                            if n_detected > 8:
                                continue  # Too many peaks, likely noise

                            # Perform fitting with detected peaks
                            if n_detected == 1:
                                # Single peak fitting
                                initial_guess = self._calculate_initial_guess_1d(x_data, y_data, target_position)
                                bounds = self._get_fitting_bounds(initial_guess, x_data, dimension)
                                fit_result = self.fit_peak_1d(x_data, y_data, initial_guess, 'voigt', bounds)
                            else:
                                # Multi-peak fitting
                                multi_result = self.fit_multi_peak_1d(x_data, y_data, detected_peaks, test_params['max_peaks_fit'])
                                if multi_result and multi_result['success']:
                                    target_peak = self.select_target_peak(multi_result, target_position)
                                    if target_peak:
                                        fit_result = {
                                            'success': True,
                                            'r_squared': multi_result['r_squared'],
                                            'fitted_curve': multi_result['fitted_curve'],
                                            'amplitude': target_peak['amplitude'],
                                            'center': target_peak['center'],
                                            'sigma': target_peak['sigma'],
                                            'gamma': target_peak['gamma'],
                                            'baseline': multi_result['baseline'],
                                            'multi_peak_info': {
                                                'n_peaks': multi_result['n_peaks'],
                                                'all_peaks': multi_result['peaks']
                                            }
                                        }
                                    else:
                                        fit_result = None
                                else:
                                    fit_result = None

                            # Evaluate fit quality
                            if fit_result and fit_result.get('success', False):
                                r_squared = fit_result.get('r_squared', -float('inf'))

                                if r_squared > best_r_squared:
                                    best_r_squared = r_squared
                                    best_fit = fit_result
                                    best_params = test_params.copy()

                                    print(f"   🎯 Improved fit: R²={r_squared:.3f}, peaks={n_detected}, "
                                          f"h={height_thresh:.3f}, d={distance_factor}, p={prominence_thresh:.3f}")

                                    # Early termination for excellent fits
                                    if r_squared > 0.9:
                                        print(f"   ✅ Excellent fit achieved, stopping optimization")
                                        break

                        except Exception as e:
                            continue  # Skip failed parameter combinations

                    if best_r_squared > 0.9:  # Break out of nested loops
                        break
                if best_r_squared > 0.9:
                    break
            if best_r_squared > 0.9:
                break

        if best_fit is not None:
            print(f"   ✅ Iterative optimization successful: R²={best_r_squared:.3f} after {attempts} attempts")
            print(f"   📊 Optimal parameters: h={best_params['height_threshold']:.3f}, "
                  f"d={best_params['distance_factor']}, p={best_params['prominence_threshold']:.3f}")
            return best_fit
        else:
            print(f"   ❌ Iterative optimization failed after {attempts} attempts")
            return None

    def adaptive_fit_1d(self, x_data, y_data, target_position, dimension='x', gui_params=None):
        """Adaptive fitting strategy with iterative parameter optimization fallback"""
        # INTEGRATION ENHANCEMENT: Check if integrated mode is enabled
        if (self.integration_parameters.get('enable_integrated_mode', False) and
            self.integrated_fitter is not None):

            print(f"🚀 Using integrated detection-fitting for {dimension.upper()}-dimension")

            # Use integrated detection-fitting workflow
            nucleus_type = '1H' if dimension.lower() == 'x' else '15N'  # Assume 15N for Y dimension

            integrated_result = self.integrated_fitter.integrated_detection_fitting(
                x_data, y_data, nucleus_type=nucleus_type,
                initial_peak_positions=[target_position] if target_position else None,
                gui_params=gui_params
            )

            if integrated_result['success'] and integrated_result.get('fitted_peaks'):
                # Convert integrated result to standard format
                best_peak = integrated_result['fitted_peaks'][0]  # Take first (best quality) peak

                fit_result = {
                    'success': True,
                    'r_squared': best_peak.get('r_squared', integrated_result.get('avg_r_squared', 0)),
                    'fitted_curve': best_peak.get('fitted_curve', []),
                    'amplitude': best_peak.get('amplitude', 0),
                    'center': best_peak.get('center', target_position),
                    'sigma': best_peak.get('sigma', 0.01),
                    'gamma': best_peak.get('gamma', 0.01),
                    'baseline': best_peak.get('baseline', 0),
                    'integration_diagnostics': integrated_result.get('integration_diagnostics', {}),
                    'method': 'integrated_detection_fitting'
                }

                print(f"   ✅ Integrated result: R²={fit_result['r_squared']:.4f}, "
                      f"converged={integrated_result['integration_diagnostics'].get('converged', False)}")

                return fit_result
            else:
                print(f"   ⚠️ Integrated detection-fitting failed, falling back to standard method")
                # Fall through to standard method

        # STANDARD METHOD (backward compatibility)
        # Get max peaks parameter
        max_peaks_fit = gui_params.get('max_peaks_fit', 4) if gui_params else 4

        # First attempt: use GUI parameters
        detected_peaks = self.detect_peaks_1d(x_data, y_data, target_position, gui_params)
        n_detected = len(detected_peaks)

        print(f"   {dimension.upper()}-dimension: {n_detected} peaks detected")

        # Attempt fitting with current parameters
        fit_result = None

        if n_detected <= 1:
            # Single peak or no peaks - use standard fitting
            initial_guess = self._calculate_initial_guess_1d(x_data, y_data, target_position)
            bounds = self._get_fitting_bounds(initial_guess, x_data, dimension)
            fit_result = self.fit_peak_1d(x_data, y_data, initial_guess, 'voigt', bounds)

        elif 2 <= n_detected <= max_peaks_fit:
            # Multi-peak case - use multi-peak deconvolution
            print(f"   Using multi-peak deconvolution for {n_detected} peaks")
            multi_result = self.fit_multi_peak_1d(x_data, y_data, detected_peaks, max_peaks_fit)

            if multi_result and multi_result['success']:
                # Select the peak closest to target
                target_peak = self.select_target_peak(multi_result, target_position)

                if target_peak:
                    # Convert to standard format
                    fit_result = {
                        'success': True,
                        'r_squared': multi_result['r_squared'],
                        'fitted_curve': multi_result['fitted_curve'],
                        'amplitude': target_peak['amplitude'],
                        'center': target_peak['center'],
                        'sigma': target_peak['sigma'],
                        'gamma': target_peak['gamma'],
                        'baseline': multi_result['baseline'],
                        'multi_peak_info': {
                            'n_peaks': multi_result['n_peaks'],
                            'all_peaks': multi_result['peaks']
                        }
                    }

            # Fallback to single peak if multi-peak fails
            if not fit_result:
                print(f"   Multi-peak failed, falling back to single peak")
                initial_guess = self._calculate_initial_guess_1d(x_data, y_data, target_position)
                bounds = self._get_fitting_bounds(initial_guess, x_data, dimension)
                fit_result = self.fit_peak_1d(x_data, y_data, initial_guess, 'voigt', bounds)

        else:
            # Too many peaks - fallback to single peak
            print(f"   Too many peaks ({n_detected}), using single peak fitting")
            initial_guess = self._calculate_initial_guess_1d(x_data, y_data, target_position)
            bounds = self._get_fitting_bounds(initial_guess, x_data, dimension)
            fit_result = self.fit_peak_1d(x_data, y_data, initial_guess, 'voigt', bounds)

        # Quality assessment and iterative optimization fallback
        if fit_result and fit_result.get('success', False):
            r_squared = fit_result.get('r_squared', -float('inf'))

            # If fit quality is poor (R² < 0.5), try iterative optimization
            if r_squared < 0.5:
                print(f"   ⚠️ Poor fit quality (R²={r_squared:.3f}), attempting iterative optimization")
                optimized_result = self.iterative_parameter_optimization(x_data, y_data, target_position, dimension, gui_params)

                if optimized_result and optimized_result.get('r_squared', -float('inf')) > r_squared:
                    return optimized_result
                else:
                    print(f"   📊 Keeping original fit (R²={r_squared:.3f})")
                    return fit_result
            else:
                # Good fit, return as-is
                return fit_result
        else:
            # Fitting failed completely, try iterative optimization as last resort
            print(f"   ❌ Initial fitting failed, attempting iterative optimization")
            optimized_result = self.iterative_parameter_optimization(x_data, y_data, target_position, dimension, gui_params)

            if optimized_result:
                return optimized_result
            else:
                # Complete failure
                return fit_result

    # =================== ENHANCED 2D FITTING ===================

    def fit_peak_voigt_2d(self, peak_x_ppm, peak_y_ppm, assignment="Unknown",
                          use_dynamic_optimization=True, all_peaks_context=None, linewidth_constraints=None): #True was false
        """
        ENHANCED 2D Voigt fitting with dynamic optimization support

        NEW OPTIMIZATION FEATURES:
        - Dynamic baseline optimization to avoid neighboring peak interference
        - Global parameter estimation from well-resolved peaks
        - Comprehensive optimization reporting
        - Backward compatibility with original method

        Parameters:
        - peak_x_ppm: X-axis peak position (ppm)
        - peak_y_ppm: Y-axis peak position (ppm)
        - assignment: Peak assignment/label
        - use_dynamic_optimization: Enable iterative optimization (default: True)
        - all_peaks_context: List of all peak positions [(x1,y1), (x2,y2), ...] for global analysis
        - linewidth_constraints: Dict with 'x' and 'y' constraints for Global Optimization Manager (NEW)

        Returns:
        - Dictionary with comprehensive fitting results or None if failed

        BACKWARD COMPATIBILITY: Set use_dynamic_optimization=False for original behavior
        """
        print(f"Fitting Voigt profiles for peak {assignment} at ({peak_x_ppm:.3f}, {peak_y_ppm:.1f}) ppm")

        # === DYNAMIC OPTIMIZATION OR STANDARD FITTING ===
        if use_dynamic_optimization and self.enhanced_fitter is not None:
            print(f"   🔄 Using dynamic optimization for {assignment}")
            
            # ENHANCEMENT: Pass GUI parameters to enhanced fitter for consistent display
            if hasattr(self.enhanced_fitter, 'set_gui_parameters'):
                self.enhanced_fitter.set_gui_parameters(self.fitting_parameters)
                print(f"   📊 GUI parameters synchronized to enhanced fitter")

            # NEW: Dynamic window optimization
            window_optimization = None
            if hasattr(self, 'optimize_window_dynamically'):
                window_optimization = self.optimize_window_dynamically(
                    peak_x_ppm, peak_y_ppm, assignment,
                    max_iterations=3,  # Limit iterations for performance
                    r2_improvement_threshold=0.03  # Accept modest improvements
                )
                
                if window_optimization.get('success') and window_optimization.get('improvement', 0) > 0.03:
                    # Use optimized windows for this peak
                    optimized_x_window = window_optimization['optimized_x_window'] 
                    optimized_y_window = window_optimization['optimized_y_window']
                    
                    print(f"   🎯 Window optimization successful:")
                    print(f"      GUI → Optimized: X={self.fitting_parameters['fitting_window_x']:.3f}→{optimized_x_window:.3f} ppm")
                    print(f"      GUI → Optimized: Y={self.fitting_parameters['fitting_window_y']:.1f}→{optimized_y_window:.1f} ppm") 
                    print(f"      R² improvement: {window_optimization['improvement']:.3f}")
                    
                    # Extract regions with optimized windows
                    regions = self.extract_peak_region(peak_x_ppm, peak_y_ppm,
                                                     optimized_x_window, optimized_y_window)
                else:
                    print(f"   📊 Using GUI windows (optimization showed no significant improvement)")
                    regions = self.extract_peak_region(peak_x_ppm, peak_y_ppm)
            else:
                regions = self.extract_peak_region(peak_x_ppm, peak_y_ppm)
                
            # ENSURE REGIONS VARIABLE EXISTS for existing code flow:
            if 'regions' not in locals() or regions is None:
                regions = self.extract_peak_region(peak_x_ppm, peak_y_ppm)
                
            # Final check that regions extraction succeeded
            if regions is None:
                print(f"   ❌ Failed to extract peak regions for {assignment}")
                return None

            try:
                # Prepare context for global parameter estimation
                x_context = None
                y_context = None

                if all_peaks_context is not None:
                    x_context = [pos[0] for pos in all_peaks_context]  # Extract x positions
                    y_context = [pos[1] for pos in all_peaks_context]  # Extract y positions
                    print(f"   Using {len(all_peaks_context)} peaks for global parameter estimation")

                # Extract linewidth constraints if provided
                x_linewidth_constraints = linewidth_constraints.get('x') if linewidth_constraints else None
                y_linewidth_constraints = linewidth_constraints.get('y') if linewidth_constraints else None

                # X-dimension fitting with dynamic optimization
                x_fit = self.enhanced_fitter.fit_peak_enhanced(
                    regions['x_ppm_scale'], regions['x_cross_section'],
                    initial_center=peak_x_ppm, nucleus_type='1H',
                    method='iterative_optimization',
                    all_peaks_context=x_context,
                    linewidth_constraints=x_linewidth_constraints
                )

                # Y-dimension fitting with dynamic optimization
                y_fit = self.enhanced_fitter.fit_peak_enhanced(
                    regions['y_ppm_scale'], regions['y_cross_section'],
                    initial_center=peak_y_ppm, nucleus_type='15N',
                    method='iterative_optimization',
                    all_peaks_context=y_context,
                    linewidth_constraints=y_linewidth_constraints
                )

                # Print optimization summary if successful
                if x_fit.get('success', False) and hasattr(self.enhanced_fitter, 'print_optimization_summary'):
                    print(f"   📊 X-dimension optimization summary:")
                    try:
                        self.enhanced_fitter.print_optimization_summary(f"{assignment}_X")
                    except Exception as e:
                        print(f"   Warning: Could not print X optimization report: {e}")

                if y_fit.get('success', False) and hasattr(self.enhanced_fitter, 'print_optimization_summary'):
                    print(f"   📊 Y-dimension optimization summary:")
                    try:
                        self.enhanced_fitter.print_optimization_summary(f"{assignment}_Y")
                    except Exception as e:
                        print(f"   Warning: Could not print Y optimization report: {e}")

            except Exception as e:
                print(f"   ⚠ Dynamic optimization failed ({e}), falling back to standard fitting")
                use_dynamic_optimization = False  # Trigger fallback below

        # === STANDARD FITTING (Original Method) ===
        if not use_dynamic_optimization or self.enhanced_fitter is None:
            if not use_dynamic_optimization:
                print(f"   📈 Using standard fitting for {assignment}")
            else:
                print(f"   📈 Using standard fitting for {assignment} (no enhanced fitter available)")

            # Extract peak regions for standard fitting
            regions = self.extract_peak_region(peak_x_ppm, peak_y_ppm)
            if regions is None:
                return None

            # Enhanced initial guess calculation with improved parameter estimation
            peak_intensity = self.nmr_data[regions['peak_indices'][1], regions['peak_indices'][0]]

            # Better baseline estimation using multiple edge points
            x_baseline = np.mean([np.mean(regions['x_cross_section'][:3]),
                                 np.mean(regions['x_cross_section'][-3:])])
            y_baseline = np.mean([np.mean(regions['y_cross_section'][:5]),
                                 np.mean(regions['y_cross_section'][-5:])])

            # Adaptive width estimation based on data
            x_width_estimate = self._estimate_peak_width(regions['x_ppm_scale'], regions['x_cross_section'], peak_x_ppm)
            y_width_estimate = self._estimate_peak_width(regions['y_ppm_scale'], regions['y_cross_section'], peak_y_ppm)

            # Amplitude estimation with baseline correction
            x_amplitude = abs(peak_intensity - x_baseline) * 1.2  # 20% buffer for fitting
            y_amplitude = abs(peak_intensity - y_baseline) * 1.2

            x_initial_guess = [
                x_amplitude,                       # amplitude
                peak_x_ppm,                        # center
                x_width_estimate * 0.6,            # sigma (Gaussian component)
                x_width_estimate * 0.4,            # gamma (Lorentzian component)
                x_baseline                         # baseline
            ]

            y_initial_guess = [
                y_amplitude,                       # amplitude
                peak_y_ppm,                        # center
                y_width_estimate * 0.6,            # sigma (Gaussian component, typically broader in Y)
                y_width_estimate * 0.4,            # gamma (Lorentzian component)
                y_baseline                         # baseline
            ]

            # Get GUI parameters if available
            gui_params = getattr(self, 'gui_params', None)

            # Fit X-dimension with adaptive strategy (single vs multi-peak)
            print(f"   🔍 Analyzing X-dimension for {assignment}")
            x_fit = self.adaptive_fit_1d(regions['x_ppm_scale'], regions['x_cross_section'],
                                        peak_x_ppm, dimension='x', gui_params=gui_params)
            print(f"   DEBUG: X-fit success={x_fit.get('success', False) if x_fit else False}")

            # Fit Y-dimension with adaptive strategy (single vs multi-peak)
            print(f"   🔍 Analyzing Y-dimension for {assignment}")
            y_fit = self.adaptive_fit_1d(regions['y_ppm_scale'], regions['y_cross_section'],
                                        peak_y_ppm, dimension='y', gui_params=gui_params)
            print(f"   DEBUG: Y-fit success={y_fit.get('success', False) if y_fit else False}")

        # Check fitting quality
        if not (x_fit.get('success', False) and y_fit.get('success', False)):
            print(f"Voigt fitting failed for {assignment}")
            return None

        # Determine overall fitting quality using LOCAL metrics (consistent with enhanced fitter)
        # Use local R-squared if available, otherwise fall back to global
        x_r_squared = x_fit.get('quality_metrics', {}).get('r_squared_local', x_fit['r_squared'])
        y_r_squared = y_fit.get('quality_metrics', {}).get('r_squared_local', y_fit['r_squared'])
        avg_r_squared = (x_r_squared + y_r_squared) / 2

        if avg_r_squared >= 0.9:
            quality = "Excellent"
        elif avg_r_squared >= 0.8:
            quality = "Good"
        elif avg_r_squared >= 0.7:
            quality = "Fair"
        else:
            quality = "Poor"

        # Store detailed results with both local and global quality metrics
        result = {
            'assignment': assignment,
            'peak_position': (peak_x_ppm, peak_y_ppm),
            'x_fit': {
                'ppm_scale': regions['x_ppm_scale'],
                'cross_section': regions['x_cross_section'],
                'fitted_curve': x_fit['fitted_curve'],
                'amplitude': x_fit['amplitude'],
                'center': x_fit['center'],
                'sigma': x_fit['sigma'],
                'gamma': x_fit['gamma'],
                'baseline': x_fit['baseline'],
                'r_squared': x_fit.get('quality_metrics', {}).get('r_squared_local', x_fit['r_squared']),
                # Include quality metrics for local/global assessment
                'quality_metrics': x_fit.get('quality_metrics', {}),
                'r_squared_local': x_fit.get('quality_metrics', {}).get('r_squared_local', x_fit['r_squared']),
                'r_squared_global': x_fit.get('quality_metrics', {}).get('r_squared_global', x_fit['r_squared']),
                # Add window size metadata for display
                'window_size': abs(regions['x_ppm_scale'][-1] - regions['x_ppm_scale'][0]) if len(regions['x_ppm_scale']) > 1 else 'unknown',
                'gui_based': True,
                'gui_window_param': self.fitting_parameters.get('fitting_window_x', 'unknown')
            },
            'y_fit': {
                'ppm_scale': regions['y_ppm_scale'],
                'cross_section': regions['y_cross_section'],
                'fitted_curve': y_fit['fitted_curve'],
                'amplitude': y_fit['amplitude'],
                'center': y_fit['center'],
                'sigma': y_fit['sigma'],
                'gamma': y_fit['gamma'],
                'baseline': y_fit['baseline'],
                'r_squared': y_fit.get('quality_metrics', {}).get('r_squared_local', y_fit['r_squared']),
                # Include quality metrics for local/global assessment
                'quality_metrics': y_fit.get('quality_metrics', {}),
                'r_squared_local': y_fit.get('quality_metrics', {}).get('r_squared_local', y_fit['r_squared']),
                'r_squared_global': y_fit.get('quality_metrics', {}).get('r_squared_global', y_fit['r_squared']),
                # Add window size metadata for display
                'window_size': abs(regions['y_ppm_scale'][-1] - regions['y_ppm_scale'][0]) if len(regions['y_ppm_scale']) > 1 else 'unknown',
                'gui_based': True,
                'gui_window_param': self.fitting_parameters.get('fitting_window_y', 'unknown')
            },
            'fitting_quality': quality,
            'avg_r_squared': avg_r_squared,
            'avg_r_squared_local': (x_r_squared + y_r_squared) / 2,
            'avg_r_squared_global': (x_fit.get('quality_metrics', {}).get('r_squared_global', x_fit['r_squared']) +
                                   y_fit.get('quality_metrics', {}).get('r_squared_global', y_fit['r_squared'])) / 2,
            'gui_window_x': self.fitting_parameters.get('fitting_window_x', 'unknown'),
            'gui_window_y': self.fitting_parameters.get('fitting_window_y', 'unknown'),
            'window_source': window_optimization.get('optimization_type', 'gui_parameters') if window_optimization else 'gui_parameters',
            'window_optimization': window_optimization,
            'timestamp': pd.Timestamp.now()
        }

        # Store in voigt_fits list
        self.voigt_fits.append(result)

        print(f"Voigt fitting completed: {quality} (R² = {avg_r_squared:.3f})")
        return result

    def optimize_window_dynamically(self, peak_x_ppm, peak_y_ppm, assignment, 
                                   initial_x_window=None, initial_y_window=None,
                                   max_iterations=5, r2_improvement_threshold=0.05):
        """
        Intelligent window size optimization starting from GUI parameters.
        
        Uses iterative refinement to find optimal window size balancing:
        - Statistical precision (more data points)
        - Peak isolation (avoiding interference)  
        - Fitting quality (maximizing R²)
        
        Parameters:
        -----------
        peak_x_ppm, peak_y_ppm : float
            Peak position for optimization
        assignment : str
            Peak identifier for logging
        initial_x_window, initial_y_window : float, optional
            Starting window sizes (defaults to GUI parameters)
        max_iterations : int
            Maximum optimization rounds
        r2_improvement_threshold : float
            Minimum R² improvement required to accept new window size
            
        Returns:
        --------
        dict : Optimized window parameters and quality metrics
        """
        import numpy as np
        
        # Phase 1: Initialize with GUI parameters
        base_x_window = initial_x_window or self.fitting_parameters['fitting_window_x']
        base_y_window = initial_y_window or self.fitting_parameters['fitting_window_y']
        
        print(f"🔄 Dynamic window optimization for {assignment}")
        print(f"   Starting windows: X={base_x_window:.3f} ppm, Y={base_y_window:.1f} ppm")
        
        # Phase 2: Baseline quality assessment
        baseline_regions = self.extract_peak_region(peak_x_ppm, peak_y_ppm, 
                                                   base_x_window, base_y_window)
        if not baseline_regions:
            return {'success': False, 'reason': 'baseline_extraction_failed'}
        
        # Perform baseline fit to establish reference quality
        baseline_x_fit = self.adaptive_fit_1d(baseline_regions['x_ppm_scale'], 
                                             baseline_regions['x_cross_section'],
                                             peak_x_ppm, dimension='x')
        baseline_y_fit = self.adaptive_fit_1d(baseline_regions['y_ppm_scale'], 
                                             baseline_regions['y_cross_section'],
                                             peak_y_ppm, dimension='y')
        
        if not (baseline_x_fit.get('success') and baseline_y_fit.get('success')):
            return {'success': False, 'reason': 'baseline_fitting_failed'}
            
        # Extract baseline quality metrics
        baseline_x_r2 = baseline_x_fit.get('quality_metrics', {}).get('r_squared_local', 
                                                                     baseline_x_fit.get('r_squared', 0))
        baseline_y_r2 = baseline_y_fit.get('quality_metrics', {}).get('r_squared_local', 
                                                                     baseline_y_fit.get('r_squared', 0))
        baseline_avg_r2 = (baseline_x_r2 + baseline_y_r2) / 2
        
        print(f"   Baseline quality: R²={baseline_avg_r2:.3f} (X={baseline_x_r2:.3f}, Y={baseline_y_r2:.3f})")
        
        # Phase 3: Analyze spectral crowding
        interference_analysis = self._analyze_peak_interference(peak_x_ppm, peak_y_ppm, 
                                                              base_x_window, base_y_window)
        
        # Phase 4: Determine optimization strategy
        if interference_analysis['isolation_level'] == 'isolated':
            # Strategy: Expand windows for better statistics
            optimization_result = self._optimize_isolated_peak_windows(
                peak_x_ppm, peak_y_ppm, assignment,
                base_x_window, base_y_window,
                baseline_avg_r2, max_iterations, r2_improvement_threshold
            )
        else:
            # Strategy: Contract windows to avoid interference
            optimization_result = self._optimize_crowded_peak_windows(
                peak_x_ppm, peak_y_ppm, assignment,
                base_x_window, base_y_window,
                baseline_avg_r2, max_iterations, r2_improvement_threshold,
                interference_analysis
            )
        
        return optimization_result

    def _analyze_peak_interference(self, peak_x_ppm, peak_y_ppm, x_window, y_window):
        """
        Analyze spectral crowding around target peak.
        
        Returns:
        --------
        dict : Interference analysis with isolation classification
        """
        import numpy as np
        
        # Extract larger region for neighbor analysis (3x current window)
        analysis_x_window = x_window * 3
        analysis_y_window = y_window * 3
        
        analysis_regions = self.extract_peak_region(peak_x_ppm, peak_y_ppm,
                                                   analysis_x_window, analysis_y_window)
        
        if not analysis_regions:
            return {'isolation_level': 'unknown', 'nearby_peaks': []}
        
        # CORRECTED: detect_peaks_1d returns a LIST, not a dict
        try:
            x_peaks_list = self.detect_peaks_1d(analysis_regions['x_ppm_scale'], 
                                              analysis_regions['x_cross_section'],
                                              target_position=peak_x_ppm)
            y_peaks_list = self.detect_peaks_1d(analysis_regions['y_ppm_scale'], 
                                              analysis_regions['y_cross_section'],
                                              target_position=peak_y_ppm)
            
            # Extract positions from peak detection results
            x_peak_positions = []
            y_peak_positions = []
            
            # Handle list of peak dictionaries (actual return format)
            if isinstance(x_peaks_list, list):
                x_peak_positions = [peak.get('position', peak.get('center', 0)) for peak in x_peaks_list 
                                   if isinstance(peak, dict)]
            
            if isinstance(y_peaks_list, list):
                y_peak_positions = [peak.get('position', peak.get('center', 0)) for peak in y_peaks_list 
                                   if isinstance(peak, dict)]
            
        except Exception as e:
            print(f"   Warning: Peak detection failed during interference analysis: {e}")
            x_peak_positions = []
            y_peak_positions = []
        
        # Count peaks within current window boundaries
        x_window_half = x_window / 2
        y_window_half = y_window / 2
        
        x_interferers = [pos for pos in x_peak_positions 
                         if abs(pos - peak_x_ppm) < x_window_half and abs(pos - peak_x_ppm) > 0.005]
        y_interferers = [pos for pos in y_peak_positions 
                         if abs(pos - peak_y_ppm) < y_window_half and abs(pos - peak_y_ppm) > 0.1]
        
        # Classify isolation level
        total_interferers = len(x_interferers) + len(y_interferers)
        
        if total_interferers == 0:
            isolation_level = 'isolated'
        elif total_interferers <= 2:
            isolation_level = 'moderate_interference'  
        else:
            isolation_level = 'heavy_interference'
        
        return {
            'isolation_level': isolation_level,
            'x_interferers': x_interferers,
            'y_interferers': y_interferers,
            'total_interferers': total_interferers,
            'x_peaks_detected': len(x_peak_positions),
            'y_peaks_detected': len(y_peak_positions)
        }

    def _optimize_isolated_peak_windows(self, peak_x_ppm, peak_y_ppm, assignment,
                                       base_x_window, base_y_window, baseline_r2,
                                       max_iterations, r2_threshold):
        """
        Optimize windows for isolated peaks by progressive expansion.
        """
        import numpy as np
        
        print(f"   Strategy: Expanding windows for isolated peak {assignment}")
        
        best_x_window = base_x_window
        best_y_window = base_y_window
        best_r2 = baseline_r2
        iteration_data = []
        
        # Test progressive window expansions
        expansion_factors = [1.0, 1.5, 2.0, 2.5, 3.0]  # Start with GUI size
        
        for iteration, factor in enumerate(expansion_factors):
            if iteration >= max_iterations:
                break
                
            test_x_window = base_x_window * factor
            test_y_window = base_y_window * factor
            
            # Test this window size
            test_regions = self.extract_peak_region(peak_x_ppm, peak_y_ppm,
                                                   test_x_window, test_y_window)
            
            if not test_regions:
                continue
                
            # Fit with test windows
            x_fit = self.adaptive_fit_1d(test_regions['x_ppm_scale'], 
                                       test_regions['x_cross_section'],
                                       peak_x_ppm, dimension='x')
            y_fit = self.adaptive_fit_1d(test_regions['y_ppm_scale'], 
                                       test_regions['y_cross_section'],
                                       peak_y_ppm, dimension='y')
            
            if not (x_fit.get('success') and y_fit.get('success')):
                continue
                
            # Calculate quality metrics
            x_r2 = x_fit.get('quality_metrics', {}).get('r_squared_local', x_fit.get('r_squared', 0))
            y_r2 = y_fit.get('quality_metrics', {}).get('r_squared_local', y_fit.get('r_squared', 0))
            avg_r2 = (x_r2 + y_r2) / 2
            
            iteration_data.append({
                'factor': factor,
                'x_window': test_x_window,
                'y_window': test_y_window,
                'r2': avg_r2,
                'x_points': len(test_regions['x_cross_section']),
                'y_points': len(test_regions['y_cross_section'])
            })
            
            print(f"   Test {iteration+1}: Factor={factor:.1f}x, R²={avg_r2:.3f}, "
                  f"Points=({len(test_regions['x_cross_section'])}, {len(test_regions['y_cross_section'])})")
            
            # Check for improvement
            if avg_r2 > best_r2 + r2_threshold:
                best_x_window = test_x_window
                best_y_window = test_y_window
                best_r2 = avg_r2
                print(f"   ✅ Improvement: R²={avg_r2:.3f} (Δ={avg_r2-baseline_r2:.3f})")
            else:
                # No significant improvement, stop expansion
                print(f"   ⏹ No improvement: R²={avg_r2:.3f}, stopping expansion")
                break
        
        optimization_improvement = best_r2 - baseline_r2
        
        return {
            'success': True,
            'optimization_type': 'isolated_expansion',
            'optimized_x_window': best_x_window,
            'optimized_y_window': best_y_window,
            'optimized_r2': best_r2,
            'baseline_r2': baseline_r2,
            'improvement': optimization_improvement,
            'iterations_tested': len(iteration_data),
            'iteration_data': iteration_data,
            'recommendation': 'expanded' if optimization_improvement > r2_threshold else 'keep_gui'
        }

    def _optimize_crowded_peak_windows(self, peak_x_ppm, peak_y_ppm, assignment,
                                      base_x_window, base_y_window, baseline_r2,
                                      max_iterations, r2_threshold, interference_analysis):
        """
        Optimize windows for crowded peaks by progressive contraction.
        """
        import numpy as np
        
        interferer_count = interference_analysis['total_interferers']
        print(f"   Strategy: Contracting windows for crowded peak {assignment} ({interferer_count} interferers)")
        
        best_x_window = base_x_window
        best_y_window = base_y_window
        best_r2 = baseline_r2
        iteration_data = []
        
        # Test progressive window contractions  
        contraction_factors = [1.0, 0.8, 0.6, 0.5, 0.4]  # Start with GUI size
        
        for iteration, factor in enumerate(contraction_factors):
            if iteration >= max_iterations:
                break
                
            test_x_window = base_x_window * factor
            test_y_window = base_y_window * factor
            
            # Don't contract below minimum reasonable size
            if test_x_window < 0.015 or test_y_window < 0.5:  # Minimum thresholds
                print(f"   ⏹ Minimum window size reached")
                break
                
            test_regions = self.extract_peak_region(peak_x_ppm, peak_y_ppm,
                                                   test_x_window, test_y_window)
            
            if not test_regions:
                continue
                
            # Check if we still have enough data points for reliable fitting
            if len(test_regions['x_cross_section']) < 10 or len(test_regions['y_cross_section']) < 10:
                print(f"   ⏹ Insufficient data points for reliable fitting")
                continue
                
            # Fit with contracted windows
            x_fit = self.adaptive_fit_1d(test_regions['x_ppm_scale'], 
                                       test_regions['x_cross_section'],
                                       peak_x_ppm, dimension='x')
            y_fit = self.adaptive_fit_1d(test_regions['y_ppm_scale'], 
                                       test_regions['y_cross_section'],
                                       peak_y_ppm, dimension='y')
            
            if not (x_fit.get('success') and y_fit.get('success')):
                continue
                
            x_r2 = x_fit.get('quality_metrics', {}).get('r_squared_local', x_fit.get('r_squared', 0))
            y_r2 = y_fit.get('quality_metrics', {}).get('r_squared_local', y_fit.get('r_squared', 0))
            avg_r2 = (x_r2 + y_r2) / 2
            
            iteration_data.append({
                'factor': factor,
                'x_window': test_x_window,
                'y_window': test_y_window,
                'r2': avg_r2,
                'x_points': len(test_regions['x_cross_section']),
                'y_points': len(test_regions['y_cross_section'])
            })
            
            print(f"   Test {iteration+1}: Factor={factor:.1f}x, R²={avg_r2:.3f}, "
                  f"Points=({len(test_regions['x_cross_section'])}, {len(test_regions['y_cross_section'])})")
            
            # Check for improvement (interference removal should increase R²)
            if avg_r2 > best_r2 + r2_threshold:
                best_x_window = test_x_window
                best_y_window = test_y_window
                best_r2 = avg_r2
                print(f"   ✅ Improvement by interference removal: R²={avg_r2:.3f} (Δ={avg_r2-baseline_r2:.3f})")
            
        optimization_improvement = best_r2 - baseline_r2
        
        return {
            'success': True,
            'optimization_type': 'crowded_contraction',
            'optimized_x_window': best_x_window,
            'optimized_y_window': best_y_window,
            'optimized_r2': best_r2,
            'baseline_r2': baseline_r2,
            'improvement': optimization_improvement,
            'interference_removed': interference_analysis,
            'iterations_tested': len(iteration_data),
            'iteration_data': iteration_data,
            'recommendation': 'contracted' if optimization_improvement > r2_threshold else 'keep_gui'
        }

    def extract_peak_region(self, peak_x_ppm, peak_y_ppm, fitting_window_x=None, fitting_window_y=None):
        """
        Extract 1D cross-sections around a peak for fitting

        Returns:
        - x_cross_section: 1D array along x-axis (1H dimension)
        - y_cross_section: 1D array along y-axis (15N/13C dimension)
        - x_ppm_scale: corresponding ppm values
        - y_ppm_scale: corresponding ppm values
        """
        if self.nmr_data is None:
            return None

        # Use default fitting windows if not provided
        if fitting_window_x is None:
            fitting_window_x = self.fitting_parameters['fitting_window_x']
        if fitting_window_y is None:
            fitting_window_y = self.fitting_parameters['fitting_window_y']

        # Find peak position in data points
        x_idx = np.argmin(np.abs(self.ppm_x_axis - peak_x_ppm))
        y_idx = np.argmin(np.abs(self.ppm_y_axis - peak_y_ppm))

        # Calculate fitting window in data points
        # Use absolute values to handle both increasing and decreasing axes
        x_ppm_range = abs(self.ppm_x_axis[0] - self.ppm_x_axis[-1])
        y_ppm_range = abs(self.ppm_y_axis[0] - self.ppm_y_axis[-1])

        x_window_points = max(10, int(fitting_window_x * len(self.ppm_x_axis) / x_ppm_range))
        y_window_points = max(10, int(fitting_window_y * len(self.ppm_y_axis) / y_ppm_range))

        # Define extraction bounds
        x_min = max(0, x_idx - x_window_points//2)
        x_max = min(len(self.ppm_x_axis), x_idx + x_window_points//2)
        y_min = max(0, y_idx - y_window_points//2)
        y_max = min(len(self.ppm_y_axis), y_idx + y_window_points//2)

        # Debug: ensure we have valid bounds
        if x_max <= x_min:
            print(f"Warning: Invalid x bounds ({x_min}, {x_max}), fixing...")
            x_min = max(0, x_idx - 5)
            x_max = min(len(self.ppm_x_axis), x_idx + 5)

        if y_max <= y_min:
            print(f"Warning: Invalid y bounds ({y_min}, {y_max}), fixing...")
            y_min = max(0, y_idx - 5)
            y_max = min(len(self.ppm_y_axis), y_idx + 5)

        # Extract cross-sections
        x_cross_section = self.nmr_data[y_idx, x_min:x_max]  # Horizontal slice
        y_cross_section = self.nmr_data[y_min:y_max, x_idx]  # Vertical slice

        # Corresponding ppm axes
        x_ppm_scale = self.ppm_x_axis[x_min:x_max]
        y_ppm_scale = self.ppm_y_axis[y_min:y_max]

        return {
            'x_cross_section': x_cross_section,
            'y_cross_section': y_cross_section,
            'x_ppm_scale': x_ppm_scale,
            'y_ppm_scale': y_ppm_scale,
            'peak_indices': (x_idx, y_idx),
            'extraction_bounds': (x_min, x_max, y_min, y_max)
        }

    @staticmethod
    def calculate_r_squared(y_true, y_pred):
        """Calculate R-squared value"""
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

    def detect_peaks_full_mode(self, **kwargs):
        """Full peak detection mode"""
        return self.detect_peaks_professional(**kwargs)

    def detect_peaks_inplace_mode(self, **kwargs):
        """In-place fitting mode with enhanced capabilities"""
        # Enhanced detection with better reference handling
        return self.detect_peaks_professional(**kwargs)

    def load_nmr_file(self, nmr_file):
        """Load NMR data file (convenience method for GUI)"""
        if (not hasattr(self, 'peak_list') or
            self.peak_list is None or
            (hasattr(self.peak_list, 'empty') and self.peak_list.empty)):
            # If no peak list is loaded yet, just store the NMR file path
            self.nmr_file_path = nmr_file
            return True
        else:
            # If peak list is already loaded, use the full load_data method
            # Create a temporary peak list file path
            temp_peak_path = getattr(self, 'peak_list_path', None)
            if temp_peak_path:
                return self.load_data(temp_peak_path, nmr_file)
            else:
                # Just load NMR data without peak list
                return self._load_nmr_data_only(nmr_file)

    def _load_nmr_data_only(self, nmr_file):
        """Load only NMR data without peak list"""
        try:
            import nmrglue as ng
            self.nmr_dict, self.nmr_data = ng.pipe.read(nmr_file)
            print(f"Loaded NMR data: {self.nmr_data.shape} from {nmr_file}")

            # Calculate PPM axes
            self._calculate_ppm_axes()

            # Estimate noise level
            self._estimate_noise_level()

            self.nmr_file_path = nmr_file
            return True
        except Exception as e:
            print(f"Error loading NMR data: {e}")
            return False

    def load_peak_list_file(self, peak_list_file):
        """Load peak list file (convenience method for GUI)"""
        try:
            import pandas as pd
            self.peak_list = pd.read_csv(peak_list_file, skipinitialspace=True)
            self.peak_list.columns = self.peak_list.columns.str.strip()
            print(f"Loaded peak list: {len(self.peak_list)} peaks from {peak_list_file}")

            self.peak_list_path = peak_list_file

            # If NMR data was already loaded, we're ready to go
            if hasattr(self, 'nmr_data') and self.nmr_data is not None:
                return True
            # If we have a stored NMR file path, load it now
            elif hasattr(self, 'nmr_file_path'):
                return self.load_data(peak_list_file, self.nmr_file_path)
            else:
                return True  # Peak list loaded, waiting for NMR data

        except Exception as e:
            print(f"Error loading peak list: {e}")
            return False

    def process_peaks(self, **kwargs):
        """Process peaks based on current mode"""
        if self.processing_mode == 'full_detection':
            return self.detect_peaks_full_mode(**kwargs)
        elif self.processing_mode == 'in_place':
            return self.detect_peaks_inplace_mode(**kwargs)
        else:
            raise ValueError(f"Unknown processing mode: {self.processing_mode}")

class EnhancedVoigtIntegrator(VoigtIntegrator):
    """Enhanced integrator with additional advanced features from inplace version"""

    def __init__(self):
        super().__init__()
        self.advanced_parameters = {
            'noise_regions': [],
            'baseline_correction': True,
            'phase_correction': True,
            'apodization': 'exponential',
            'zero_filling': 2,
            'peak_shape_analysis': True
        }
        self.statistics = {
            'processed_spectra': 0,
            'successful_fits': 0,
            'failed_fits': 0,
            'average_quality': 0.0
        }

        # INTEGRATION ENHANCEMENT: Initialize integrated detection-fitting system
        if INTEGRATED_DETECTION_AVAILABLE:
            self.integrated_fitter = create_integrated_fitter(
                self.enhanced_fitter if hasattr(self, 'enhanced_fitter') else None
            )
            print("🚀 Integrated detection-fitting system initialized")
        else:
            self.integrated_fitter = None

        # Integration parameters
        self.integration_mode = 'standard'  # 'standard', 'integrated', 'adaptive'
        self.integration_parameters = {
            'enable_integrated_mode': False,  # Default off for backward compatibility
            'adaptive_thresholds': True,
            'multi_resolution_detection': True,
            'quality_filtering': True,
        }

    def set_advanced_parameters(self, **params):
        """Set advanced processing parameters"""
        self.advanced_parameters.update(params)

    def set_integration_mode(self, mode='standard', **integration_params):
        """
        Set integration mode for detection-fitting workflow (INTEGRATION ENHANCEMENT)

        Args:
            mode: 'standard', 'integrated', or 'adaptive'
            **integration_params: integration-specific parameters
        """
        if mode not in ['standard', 'integrated', 'adaptive']:
            raise ValueError("Mode must be 'standard', 'integrated', or 'adaptive'")

        if mode != 'standard' and not INTEGRATED_DETECTION_AVAILABLE:
            print("⚠️ Integrated detection-fitting not available, falling back to standard mode")
            mode = 'standard'

        self.integration_mode = mode
        self.integration_parameters.update(integration_params)

        if mode == 'integrated':
            self.integration_parameters['enable_integrated_mode'] = True
            print("🚀 Integrated detection-fitting mode enabled")
        elif mode == 'adaptive':
            self.integration_parameters['enable_integrated_mode'] = True
            self.integration_parameters['adaptive_thresholds'] = True
            self.integration_parameters['multi_resolution_detection'] = True
            print("🎯 Adaptive integrated detection-fitting mode enabled")
        else:
            self.integration_parameters['enable_integrated_mode'] = False
            print("📊 Standard detection-fitting mode")

        return self.integration_parameters.copy()

    def get_integration_status(self):
        """Get current integration mode and parameters"""
        return {
            'mode': self.integration_mode,
            'parameters': self.integration_parameters.copy(),
            'integrated_fitter_available': self.integrated_fitter is not None,
            'capabilities': {
                'adaptive_thresholds': INTEGRATED_DETECTION_AVAILABLE,
                'multi_resolution_detection': INTEGRATED_DETECTION_AVAILABLE,
                'quality_scoring': INTEGRATED_DETECTION_AVAILABLE,
                'chemical_shift_context': INTEGRATED_DETECTION_AVAILABLE
            }
        }

    def get_statistics(self):
        """Get processing statistics"""
        return self.statistics.copy()

    def reset_statistics(self):
        """Reset processing statistics"""
        self.statistics = {
            'processed_spectra': 0,
            'successful_fits': 0,
            'failed_fits': 0,
            'average_quality': 0.0
        }

    def update_statistics(self, successful=True, quality=0.0):
        """Update processing statistics"""
        self.statistics['processed_spectra'] += 1
        if successful:
            self.statistics['successful_fits'] += 1
        else:
            self.statistics['failed_fits'] += 1

        # Update average quality
        total = self.statistics['successful_fits']
        if total > 0:
            current_avg = self.statistics['average_quality']
            self.statistics['average_quality'] = ((current_avg * (total - 1)) + quality) / total

    def enhanced_peak_fitting(self, peak_x_ppm, peak_y_ppm, assignment="Unknown"):
        """Enhanced peak fitting with advanced features"""
        result = self.fit_peak_voigt_2d(peak_x_ppm, peak_y_ppm, assignment, use_dynamic_optimization=True) #GM False to True

        if result:
            self.update_statistics(True, result['avg_r_squared'])
        else:
            self.update_statistics(False)

        return result

    def batch_peak_fitting(self, peaks_list, progress_callback=None):
        """Batch fitting with progress tracking"""
        results = []
        total_peaks = len(peaks_list)

        for i, peak in enumerate(peaks_list):
            if progress_callback:
                progress_callback(i, total_peaks, peak.get('assignment', f'Peak_{i+1}'))

            result = self.enhanced_peak_fitting(
                peak['ppm_x'],
                peak['ppm_y'],
                peak.get('assignment', f'Peak_{i+1}')
            )

            if result:
                results.append(result)

        return results

    def optimize_peak_list_globally(self, peak_list: List[Tuple[float, float, str]],
                                   convergence_threshold: float = 0.05,
                                   max_rounds: int = 5) -> Dict[str, Any]:
        """
        GLOBAL OPTIMIZATION: Two-phase optimization of multiple peaks using
        GlobalOptimizationManager

        This method implements the sophisticated two-phase approach:

        Phase 1 (Survey):
        - Fits all peaks with standard methods
        - Collects linewidth statistics from successful peaks
        - Identifies peaks needing improvement

        Phase 2 (Iterative Refinement):
        - Applies linewidth constraints from successful peaks
        - Uses multi-peak deconvolution for complex cases
        - Continues until convergence (<5% improvement per round)

        Args:
            peak_list: List of (x_ppm, y_ppm, peak_id) tuples
            convergence_threshold: Stop when <N% peaks improve per round
            max_rounds: Maximum refinement rounds

        Returns:
            Comprehensive optimization report with all results and statistics

        Example:
            peak_list = [(8.0, 120.0, "Peak1"), (7.9, 118.0, "Peak2")]
            report = integrator.optimize_peak_list_globally(peak_list)
            success_rate = report['optimization_summary']['final_success_rate']
        """
        from lunaNMR.utils.global_optimization_manager import GlobalOptimizationManager

        print("🚀 Starting Global Peak Optimization with Two-Phase Approach")
        print("="*70)

        # Initialize Global Optimization Manager
        optimizer = GlobalOptimizationManager(
            convergence_threshold=convergence_threshold,
            max_rounds=max_rounds
        )

        # Execute two-phase optimization
        optimization_report = optimizer.optimize_peak_list(peak_list, self)

        print("\n✅ Global Optimization Complete!")
        print(f"   Final success rate: {optimization_report['optimization_summary']['final_success_rate']:.1f}%")
        print(f"   Optimization rounds: {optimization_report['optimization_summary']['total_rounds']}")

        return optimization_report
