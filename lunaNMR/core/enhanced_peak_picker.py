# ABOUTME: Prominence-based peak detection with graph clustering and adaptive windows for 2D NMR spectra.
# ABOUTME: Uses scipy.ndimage.maximum_filter, networkx clustering, and S/N analysis for robust candidate detection.
#!/usr/bin/env python3
"""
Enhanced Peak Picker Module

This module provides intelligent peak detection for NMR spectroscopy with improved
SNR filtering, overlapping peak handling, and validation. Works in conjunction
with the enhanced Voigt fitter to provide comprehensive peak analysis.

Key features:
- Intelligent peak detection with adaptive thresholds
- SNR-based filtering with noise estimation
- Overlapping peak detection and separation
- Peak validation based on NMR-specific criteria
- Integration with enhanced Voigt fitting

Author: Guillaume Mas
Date: 2025
"""

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_prominences, peak_widths
from scipy.ndimage import gaussian_filter, maximum_filter, minimum_filter
from sklearn.cluster import DBSCAN
import warnings
warnings.filterwarnings('ignore')

# Import enhanced fitter
try:
    from lunaNMR.core.enhanced_voigt_fitter import EnhancedVoigtFitter
    ENHANCED_FITTING_AVAILABLE = True
except ImportError:
    ENHANCED_FITTING_AVAILABLE = False

class EnhancedPeakPicker:
    """Enhanced peak picker with intelligent detection and validation"""

    def __init__(self):
        self.detection_parameters = {
            'min_snr': 3.0,              # Minimum signal-to-noise ratio
            'min_height_percentile': 5,   # Minimum height as percentile of data
            'min_prominence_factor': 0.1, # Minimum prominence relative to peak height
            'min_width_points': 3,        # Minimum peak width in data points
            'max_width_factor': 0.1,      # Maximum width as fraction of data range
            'noise_estimation_percentile': 10,  # Percentile for noise estimation
            'smoothing_sigma': 1.0,       # Gaussian smoothing parameter
            'overlap_threshold': 0.8,     # Threshold for detecting overlapping peaks
            'validation_fit_threshold': 0.5,  # Minimum R² for peak validation
            'use_centroid': True,         # Use centroid instead of pixel maximum
            'centroid_method': 'top_contour',  # 'top_contour' or 'window'
            'centroid_window_size': 3,    # Half-width of window for 'window' method (in pixels)
            'centroid_intensity_band': 0.05,  # Intensity band for 'top_contour' method (±5%)
            'centroid_window_x_ppm': 0.04,    # Max shift X in ppm (from GUI, default)
            'centroid_window_y_ppm': 0.14     # Max shift Y in ppm (from GUI, default)
        }

        # NMR-specific parameters (corrected ranges)
        self.nmr_ranges = {
            '1H': {'min': 5.5, 'max': 12.0, 'typical_width': 0.02, 'min_snr': 5.0},
            '15N': {'min': 100.0, 'max': 140.0, 'typical_width': 1.5, 'min_snr': 3.0},
            '13C': {'min': 0.0, 'max': 220.0, 'typical_width': 1.0, 'min_snr': 4.0}
        }

        if ENHANCED_FITTING_AVAILABLE:
            self.fitter = EnhancedVoigtFitter()
        else:
            self.fitter = None

        self.last_detection_stats = {}

    def detect_nucleus_type(self, ppm_range):
        """Detect nucleus type based on ppm range"""
        ppm_span = abs(ppm_range[1] - ppm_range[0])
        center_ppm = (ppm_range[0] + ppm_range[1]) / 2

        # Check ¹H range (5.5-12 ppm)
        if 5.0 <= center_ppm <= 13.0 and ppm_span < 20:
            return '1H'

        # Check ¹⁵N range (100-140 ppm)
        elif 90 <= center_ppm <= 150 and ppm_span < 100:
            return '15N'

        # Check ¹³C range (0-220 ppm)
        elif 0 <= center_ppm <= 230 and ppm_span > 50:
            return '13C'

        return '1H'  # Default

    def estimate_noise_level(self, data_2d, method='robust_corners'):
        """
        Estimate noise level from 2D NMR data

        Methods:
        - 'robust_corners': Use corner regions with outlier rejection
        - 'percentile': Use low percentile of entire dataset
        - 'mad': Median Absolute Deviation approach
        """
        try:
            if method == 'robust_corners':
                # Sample from corners of the 2D spectrum
                h, w = data_2d.shape
                corner_size = min(max(10, h//20), 50)  # Adaptive corner size

                # Extract corners
                corners = [
                    data_2d[:corner_size, :corner_size],
                    data_2d[:corner_size, -corner_size:],
                    data_2d[-corner_size:, :corner_size],
                    data_2d[-corner_size:, -corner_size:]
                ]

                corner_data = np.concatenate([c.flatten() for c in corners])

                # Remove outliers using IQR method
                q25, q75 = np.percentile(corner_data, [25, 75])
                iqr = q75 - q25
                lower_bound = q25 - 1.5 * iqr
                upper_bound = q75 + 1.5 * iqr

                clean_data = corner_data[(corner_data >= lower_bound) & (corner_data <= upper_bound)]

                if len(clean_data) > 10:
                    noise_level = np.std(clean_data)
                else:
                    noise_level = np.std(corner_data)

            elif method == 'percentile':
                # Use low percentile as noise estimate
                noise_level = np.percentile(np.abs(data_2d), self.detection_parameters['noise_estimation_percentile'])

            elif method == 'mad':
                # Median Absolute Deviation
                median_val = np.median(data_2d)
                mad = np.median(np.abs(data_2d - median_val))
                noise_level = mad * 1.4826  # Scale factor for normal distribution

            else:
                # Fallback to simple std of corners
                h, w = data_2d.shape
                corners = np.concatenate([
                    data_2d[:h//10, :w//10].flatten(),
                    data_2d[:h//10, -w//10:].flatten(),
                    data_2d[-h//10:, :w//10].flatten(),
                    data_2d[-h//10:, -w//10:].flatten()
                ])
                noise_level = np.std(corners)

            return max(noise_level, 1e-10)  # Avoid zero noise

        except Exception as e:
            return np.std(data_2d) * 0.1  # Fallback

    def preprocess_data_for_detection(self, data_2d, sigma=None, adaptive_smoothing=True):
        """
        Preprocess 2D data for better peak detection with adaptive smoothing
        """
        if sigma is None:
            if adaptive_smoothing:
                sigma = self.get_adaptive_smoothing_parameter(data_2d)
            else:
                sigma = self.detection_parameters['smoothing_sigma']

        # Adaptive Gaussian smoothing based on data characteristics
        if sigma > 0:
            smoothed = gaussian_filter(data_2d, sigma=sigma)
        else:
            smoothed = data_2d

        return smoothed

    def get_adaptive_smoothing_parameter(self, data_2d):
        """
        Calculate adaptive smoothing parameter based on data noise and complexity
        """
        # Estimate noise level
        noise_level = self.estimate_noise_level(data_2d)
        data_max = np.max(data_2d)
        data_range = data_max - np.min(data_2d)

        # Calculate SNR estimate
        snr_estimate = data_max / (noise_level + 1e-10)

        # Adaptive smoothing based on data quality
        base_sigma = self.detection_parameters['smoothing_sigma']

        if snr_estimate > 50:  # High quality data
            adaptive_sigma = base_sigma * 0.5  # Less smoothing for clean data
        elif snr_estimate > 20:  # Good quality data
            adaptive_sigma = base_sigma * 0.8  # Moderate smoothing
        elif snr_estimate > 10:  # Fair quality data
            adaptive_sigma = base_sigma * 1.0  # Standard smoothing
        elif snr_estimate > 5:   # Poor quality data
            adaptive_sigma = base_sigma * 1.5  # More smoothing for noisy data
        else:  # Very noisy data
            adaptive_sigma = base_sigma * 2.0  # Heavy smoothing

        # Additional scaling based on data complexity
        data_complexity = np.std(np.diff(data_2d.flatten())) / (noise_level + 1e-10)
        complexity_factor = min(1.5, max(0.5, 1.0 + 0.1 * np.log10(max(data_complexity, 0.1))))

        adaptive_sigma *= complexity_factor

        # Ensure reasonable bounds
        adaptive_sigma = max(0.3, min(3.0, adaptive_sigma))

        return adaptive_sigma

    def detect_peaks_2d_initial(self, data_2d, ppm_x_axis, ppm_y_axis, nucleus_type=None):
        """
        Initial 2D peak detection using local maxima
        """
        if nucleus_type is None:
            nucleus_type = self.detect_nucleus_type([ppm_x_axis[0], ppm_x_axis[-1]])

        # Get nucleus-specific parameters
        nmr_params = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])

        # Preprocess data
        smoothed_data = self.preprocess_data_for_detection(data_2d)

        # Estimate noise level
        noise_level = self.estimate_noise_level(data_2d)

        # Calculate adaptive dynamic thresholds
        min_snr = nmr_params['min_snr']
        height_threshold = noise_level * min_snr

        # Adaptive percentile threshold based on data quality
        data_snr = np.max(smoothed_data) / (noise_level + 1e-10)
        adaptive_percentile = self.get_adaptive_percentile_threshold(data_snr)
        percentile_threshold = np.percentile(smoothed_data, adaptive_percentile)

        # Use the higher of SNR-based and percentile-based thresholds
        final_threshold = max(height_threshold, percentile_threshold * 0.1)

        # Find local maxima using maximum filter
        neighborhood_size = max(3, min(data_2d.shape) // 50)  # Adaptive neighborhood
        local_max_mask = (smoothed_data == maximum_filter(smoothed_data, size=neighborhood_size))

        # Apply adaptive intensity and prominence thresholds
        intensity_mask = smoothed_data > final_threshold

        # Calculate adaptive prominence threshold
        adaptive_prominence = self.get_adaptive_prominence_threshold(smoothed_data, noise_level)

        # Apply prominence filtering to local maxima
        prominence_mask = self.apply_adaptive_prominence_filter(smoothed_data, local_max_mask, adaptive_prominence)

        # Combine all filters
        peak_candidates = local_max_mask & intensity_mask & prominence_mask

        # Extract peak positions
        peak_indices = np.where(peak_candidates)

        if len(peak_indices[0]) == 0:
            return []

        # Convert to peak list with coordinates and intensities
        initial_peaks = []
        for i, (y_idx, x_idx) in enumerate(zip(peak_indices[0], peak_indices[1])):
            # Calculate centroid position (geometric center) if enabled
            if self.detection_parameters.get('use_centroid', True):
                centroid_method = self.detection_parameters.get('centroid_method', 'top_contour')

                if centroid_method == 'top_contour':
                    # Use top contour method (center of pixels near maximum)
                    centroid_x, centroid_y = self.calculate_top_contour_center(
                        smoothed_data, y_idx, x_idx, ppm_y_axis, ppm_x_axis,
                        intensity_band=self.detection_parameters.get('centroid_intensity_band', 0.05),
                        max_shift_x=self.detection_parameters.get('centroid_window_x_ppm', 0.04),
                        max_shift_y=self.detection_parameters.get('centroid_window_y_ppm', 0.14)
                    )
                else:
                    # Use window method (traditional center of mass in fixed window)
                    centroid_x, centroid_y = self.calculate_peak_centroid(
                        smoothed_data, y_idx, x_idx, ppm_y_axis, ppm_x_axis,
                        window_size=self.detection_parameters.get('centroid_window_size', 3)
                    )

                ppm_x = centroid_x
                ppm_y = centroid_y
            else:
                # Use pixel maximum (original behavior)
                ppm_x = ppm_x_axis[x_idx]
                ppm_y = ppm_y_axis[y_idx]

            peak_info = {
                'id': i,
                'x_idx': x_idx,  # Keep original pixel index for reference
                'y_idx': y_idx,
                'ppm_x': ppm_x,  # Use centroid or pixel position
                'ppm_y': ppm_y,  # Use centroid or pixel position
                'ppm_x_pixel': ppm_x_axis[x_idx],  # Store original pixel position
                'ppm_y_pixel': ppm_y_axis[y_idx],
                'intensity': smoothed_data[y_idx, x_idx],
                'snr': smoothed_data[y_idx, x_idx] / noise_level,
                'nucleus_type': nucleus_type
            }
            initial_peaks.append(peak_info)

        # Sort by intensity (highest first)
        initial_peaks.sort(key=lambda x: x['intensity'], reverse=True)

        return initial_peaks

    def calculate_peak_centroid(self, data_2d, y_idx, x_idx, ppm_y_axis, ppm_x_axis, window_size=3):
        """
        Calculate the centroid (center of mass) of a peak in a local window.

        This provides sub-pixel accuracy by computing the weighted average position
        using intensities as weights. More accurate than using the maximum pixel.

        Parameters:
        -----------
        data_2d : ndarray
            2D NMR data
        y_idx, x_idx : int
            Pixel indices of detected maximum
        ppm_y_axis, ppm_x_axis : ndarray
            PPM axes
        window_size : int, optional
            Half-width of window around peak (in pixels), default=3

        Returns:
        --------
        centroid_x, centroid_y : float
            Centroid position in ppm
        """
        import numpy as np

        # Extract local region around peak
        y_start = max(0, y_idx - window_size)
        y_end = min(data_2d.shape[0], y_idx + window_size + 1)
        x_start = max(0, x_idx - window_size)
        x_end = min(data_2d.shape[1], x_idx + window_size + 1)

        local_data = data_2d[y_start:y_end, x_start:x_end]
        local_y_ppm = ppm_y_axis[y_start:y_end]
        local_x_ppm = ppm_x_axis[x_start:x_end]

        # Use only positive intensities for centroid calculation
        local_data = np.maximum(local_data, 0)

        # Calculate total intensity
        total_intensity = np.sum(local_data)

        if total_intensity == 0:
            # Fallback to pixel position if no signal
            return ppm_x_axis[x_idx], ppm_y_axis[y_idx]

        # Create 2D grids for centroid calculation
        y_grid, x_grid = np.meshgrid(local_y_ppm, local_x_ppm, indexing='ij')

        # Calculate weighted average (center of mass)
        centroid_y = np.sum(y_grid * local_data) / total_intensity
        centroid_x = np.sum(x_grid * local_data) / total_intensity

        return centroid_x, centroid_y

    def calculate_top_contour_center(self, data_2d, y_idx, x_idx, ppm_y_axis, ppm_x_axis,
                                     intensity_band=0.05, max_shift_x=0.04, max_shift_y=0.14):
        """
        Calculate geometric center of pixels within top contour (near maximum intensity).

        This method finds the center of all pixels with intensity close to the maximum
        (within intensity_band), providing a more accurate center than the single maximum
        pixel. Useful for peaks with flat tops or multiple pixels at similar intensity.

        Parameters:
        -----------
        data_2d : ndarray
            2D NMR data
        y_idx, x_idx : int
            Pixel indices of detected maximum
        ppm_y_axis, ppm_x_axis : ndarray
            PPM axes
        intensity_band : float, optional
            Fractional intensity band around maximum (default 0.05 = ±5%)
            Uses pixels with intensity in range [max×(1-band), max×(1+band)]
        max_shift_x, max_shift_y : float, optional
            Maximum allowed shift from original position in ppm (safety constraint)
            Default values should come from GUI centroid_window parameters

        Returns:
        --------
        center_x, center_y : float
            Geometric center position in ppm
        """
        import numpy as np

        # Get original position
        original_x = ppm_x_axis[x_idx]
        original_y = ppm_y_axis[y_idx]
        peak_max = data_2d[y_idx, x_idx]

        # Define search window based on max_shift constraints
        x_ppm_per_pixel = abs(ppm_x_axis[1] - ppm_x_axis[0]) if len(ppm_x_axis) > 1 else 0.004
        y_ppm_per_pixel = abs(ppm_y_axis[1] - ppm_y_axis[0]) if len(ppm_y_axis) > 1 else 0.08

        search_radius_x = int(np.ceil(max_shift_x / x_ppm_per_pixel))
        search_radius_y = int(np.ceil(max_shift_y / y_ppm_per_pixel))

        x_start = max(0, x_idx - search_radius_x)
        x_end = min(data_2d.shape[1], x_idx + search_radius_x + 1)
        y_start = max(0, y_idx - search_radius_y)
        y_end = min(data_2d.shape[0], y_idx + search_radius_y + 1)

        # Extract search region
        search_region = data_2d[y_start:y_end, x_start:x_end]

        # Define intensity range for "top contour"
        # intensity_band = 0.05 means use pixels with 95%-105% of maximum
        intensity_min = peak_max * (1.0 - intensity_band)
        intensity_max = peak_max * (1.0 + intensity_band)  # Allow slight overshoot due to noise

        # Find all pixels within intensity band (the "top contour")
        mask = (search_region >= intensity_min) & (search_region <= intensity_max)

        if not np.any(mask):
            # Fallback to original position if no pixels in band
            return original_x, original_y

        # Get coordinates and intensities of top contour pixels
        local_y_indices, local_x_indices = np.where(mask)
        global_y_indices = local_y_indices + y_start
        global_x_indices = local_x_indices + x_start

        intensities = data_2d[global_y_indices, global_x_indices]

        # Convert to ppm coordinates
        ppm_x_vals = ppm_x_axis[global_x_indices]
        ppm_y_vals = ppm_y_axis[global_y_indices]

        # Calculate intensity-weighted geometric center
        total_intensity = np.sum(intensities)
        center_x = np.sum(ppm_x_vals * intensities) / total_intensity
        center_y = np.sum(ppm_y_vals * intensities) / total_intensity

        # Apply safety constraint: limit maximum shift
        shift_x = center_x - original_x
        shift_y = center_y - original_y

        # Clip shifts to maximum allowed values
        shift_x = np.clip(shift_x, -max_shift_x, max_shift_x)
        shift_y = np.clip(shift_y, -max_shift_y, max_shift_y)

        final_x = original_x + shift_x
        final_y = original_y + shift_y

        return final_x, final_y

    def filter_peaks_by_snr(self, peaks, min_snr=None):
        """Filter peaks based on signal-to-noise ratio"""
        if min_snr is None:
            min_snr = self.detection_parameters['min_snr']

        filtered_peaks = [p for p in peaks if p.get('snr', 0) >= min_snr]
        return filtered_peaks

    def detect_overlapping_peaks(self, peaks, overlap_threshold=None):
        """
        Detect potentially overlapping peaks using distance clustering
        """
        if overlap_threshold is None:
            overlap_threshold = self.detection_parameters['overlap_threshold']

        if len(peaks) < 2:
            return peaks

        # Extract coordinates for clustering
        coordinates = np.array([[p['ppm_x'], p['ppm_y']] for p in peaks])

        # Normalize coordinates for clustering (different scales for x and y)
        x_scale = np.std(coordinates[:, 0]) if len(coordinates) > 1 else 1
        y_scale = np.std(coordinates[:, 1]) if len(coordinates) > 1 else 1

        if x_scale > 0 and y_scale > 0:
            normalized_coords = coordinates.copy()
            normalized_coords[:, 0] /= x_scale
            normalized_coords[:, 1] /= y_scale

            # Use DBSCAN clustering to identify overlapping groups
            try:
                clustering = DBSCAN(eps=overlap_threshold, min_samples=2)
                cluster_labels = clustering.fit_predict(normalized_coords)

                # Mark overlapping peaks
                for i, peak in enumerate(peaks):
                    if cluster_labels[i] >= 0:  # -1 means noise/single peak
                        peak['overlapping'] = True
                        peak['cluster_id'] = cluster_labels[i]
                    else:
                        peak['overlapping'] = False
                        peak['cluster_id'] = -1

            except Exception as e:
                # Mark all as non-overlapping if clustering fails
                for peak in peaks:
                    peak['overlapping'] = False
                    peak['cluster_id'] = -1

        return peaks

    def validate_peak_by_fitting(self, data_2d, ppm_x_axis, ppm_y_axis, peak,
                                window_size_x=0.3, window_size_y=5.0):
        """
        Validate peak by attempting to fit a Voigt profile
        """
        if self.fitter is None:
            # Cannot validate without fitter
            peak['fit_validation'] = 'skipped'
            peak['fit_quality'] = 0.0
            return peak

        try:
            # Extract region around peak for fitting
            x_center_idx = peak['x_idx']
            y_center_idx = peak['y_idx']

            # Calculate window in data points
            x_window_points = max(5, int(window_size_x * len(ppm_x_axis) / abs(ppm_x_axis[-1] - ppm_x_axis[0])))
            y_window_points = max(5, int(window_size_y * len(ppm_y_axis) / abs(ppm_y_axis[-1] - ppm_y_axis[0])))

            # Define extraction bounds
            x_min = max(0, x_center_idx - x_window_points//2)
            x_max = min(len(ppm_x_axis), x_center_idx + x_window_points//2)
            y_min = max(0, y_center_idx - y_window_points//2)
            y_max = min(len(ppm_y_axis), y_center_idx + y_window_points//2)

            # Extract 1D cross-sections for fitting
            x_cross_section = data_2d[y_center_idx, x_min:x_max]
            y_cross_section = data_2d[y_min:y_max, x_center_idx]

            x_ppm_section = ppm_x_axis[x_min:x_max]
            y_ppm_section = ppm_y_axis[y_min:y_max]

            # Attempt fitting in both dimensions
            x_fit = self.fitter.fit_peak_enhanced(x_ppm_section, x_cross_section,
                                                 peak['ppm_x'], peak['nucleus_type'])
            y_fit = self.fitter.fit_peak_enhanced(y_ppm_section, y_cross_section,
                                                 peak['ppm_y'], peak['nucleus_type'])

            # Calculate combined fit quality
            if x_fit['success'] and y_fit['success']:
                avg_r_squared = (x_fit['r_squared'] + y_fit['r_squared']) / 2
                peak['fit_validation'] = 'success'
                peak['fit_quality'] = avg_r_squared
                peak['x_fit_result'] = x_fit
                peak['y_fit_result'] = y_fit

                # Update peak position with fitted centers
                peak['ppm_x_fitted'] = x_fit['center']
                peak['ppm_y_fitted'] = y_fit['center']

            else:
                peak['fit_validation'] = 'failed'
                peak['fit_quality'] = 0.0

        except Exception as e:
            peak['fit_validation'] = 'error'
            peak['fit_quality'] = 0.0

        return peak

    def comprehensive_peak_detection(self, data_2d, ppm_x_axis, ppm_y_axis,
                                   nucleus_type=None, validate_fits=True):
        """
        Hierarchical peak detection with enhanced multi-peak handling

        Enhanced Steps:
        1. Initial peak detection using local maxima
        2. Hierarchical intensity-based detection (dominant → medium → small)
        3. SNR-based filtering with adaptive thresholds
        4. Overlapping peak detection and refinement
        5. Optional fit-based validation
        6. Final quality assessment
        """
        # Step 1: Initial detection
        initial_peaks = self.detect_peaks_2d_initial(data_2d, ppm_x_axis, ppm_y_axis, nucleus_type)

        if len(initial_peaks) == 0:
            return []

        # Step 2: Hierarchical intensity-based detection
        hierarchical_peaks = self.hierarchical_peak_refinement(initial_peaks, data_2d)

        # Step 3: Enhanced SNR filtering with adaptive thresholds
        snr_filtered = self.filter_peaks_by_snr_adaptive(hierarchical_peaks, data_2d)

        # Step 4: Enhanced overlapping peak detection
        overlap_analyzed = self.detect_overlapping_peaks_enhanced(snr_filtered, data_2d)

        # Step 5: Optional fit validation with enhanced criteria
        if validate_fits and len(overlap_analyzed) > 0:
            validated_peaks = []

            for i, peak in enumerate(overlap_analyzed):
                validated_peak = self.validate_peak_by_fitting_enhanced(data_2d, ppm_x_axis, ppm_y_axis, peak)

                # Enhanced validation criteria for small and overlapping peaks
                min_fit_quality = self.get_adaptive_fit_threshold(peak, overlap_analyzed)
                if validated_peak.get('fit_quality', 0) >= min_fit_quality:
                    validated_peaks.append(validated_peak)

            final_peaks = validated_peaks
        else:
            final_peaks = overlap_analyzed

        # Step 6: Final quality assessment and sorting
        for peak in final_peaks:
            # Calculate composite quality score
            snr_score = min(peak.get('snr', 0) / 10, 1.0)  # Normalize to 0-1
            fit_score = peak.get('fit_quality', 0)
            overlap_penalty = 0.8 if peak.get('overlapping', False) else 1.0

            peak['quality_score'] = (snr_score * 0.4 + fit_score * 0.6) * overlap_penalty

        # Sort by quality score
        final_peaks.sort(key=lambda x: x['quality_score'], reverse=True)

        # Store detection statistics
        self.last_detection_stats = {
            'initial_peaks': len(initial_peaks),
            'after_snr_filter': len(snr_filtered),
            'after_overlap_analysis': len(overlap_analyzed),
            'final_peaks': len(final_peaks),
            'validation_enabled': validate_fits,
            'nucleus_type': nucleus_type or 'auto-detected'
        }

        return final_peaks

    def hierarchical_peak_refinement(self, initial_peaks, data_2d):
        """
        Hierarchical peak detection strategy for complex overlapping scenarios

        Strategy stages:
        1. Dominant peaks (>50% max intensity)
        2. Medium peaks (10-50% max intensity)
        3. Small peaks (<10% max intensity)
        4. Global refinement
        """
        if not initial_peaks:
            return initial_peaks

        # Sort peaks by intensity
        sorted_peaks = sorted(initial_peaks, key=lambda p: p.get('intensity', 0), reverse=True)
        max_intensity = sorted_peaks[0].get('intensity', 1.0)

        refined_peaks = []

        # Stage 1: Dominant peaks (>50% max intensity)
        dominant_threshold = max_intensity * 0.5
        dominant_peaks = [p for p in sorted_peaks if p.get('intensity', 0) >= dominant_threshold]
        refined_peaks.extend(dominant_peaks)

        # Stage 2: Medium peaks (10-50% max intensity)
        medium_threshold = max_intensity * 0.1
        medium_peaks = [p for p in sorted_peaks if medium_threshold <= p.get('intensity', 0) < dominant_threshold]

        # Enhanced detection for medium peaks near dominant ones with adaptive distance
        nucleus_type = self.detect_nucleus_type([max_intensity, max_intensity])  # Rough estimate
        data_complexity = len(sorted_peaks) / 10.0  # Simple complexity estimate
        adaptive_shadow_distance = self.get_adaptive_distance_threshold(nucleus_type, data_complexity)

        for peak in medium_peaks:
            # Check if this medium peak is in the shadow of a dominant peak
            is_in_shadow = False
            for dom_peak in dominant_peaks:
                distance = np.sqrt((peak['ppm_x'] - dom_peak['ppm_x'])**2 +
                                 (peak['ppm_y'] - dom_peak['ppm_y'])**2)
                if distance < adaptive_shadow_distance:  # Adaptive shadow detection
                    is_in_shadow = True
                    break

            if not is_in_shadow:
                refined_peaks.append(peak)
            else:
                # Mark as potentially overlapping for special handling
                peak['shadow_peak'] = True
                refined_peaks.append(peak)

        # Stage 3: Small peaks (<10% max intensity) with enhanced sensitivity
        small_peaks = [p for p in sorted_peaks if p.get('intensity', 0) < medium_threshold]

        # Use lower SNR threshold for small peaks
        original_min_snr = self.detection_parameters['min_snr']
        self.detection_parameters['min_snr'] = max(2.0, original_min_snr * 0.6)  # Reduced threshold

        for peak in small_peaks:
            # Enhanced validation for small peaks
            local_noise = self.estimate_local_noise_around_peak(data_2d, peak)
            peak_snr = peak.get('intensity', 0) / (local_noise + 1e-10)

            if peak_snr >= self.detection_parameters['min_snr']:
                peak['small_peak'] = True
                refined_peaks.append(peak)

        # Restore original SNR threshold
        self.detection_parameters['min_snr'] = original_min_snr

        return refined_peaks

    def filter_peaks_by_snr_adaptive(self, peaks, data_2d):
        """Enhanced SNR filtering with adaptive thresholds based on peak characteristics"""
        filtered_peaks = []

        for peak in peaks:
            # Adaptive SNR threshold based on peak type
            if peak.get('small_peak', False):
                min_snr_threshold = max(2.0, self.detection_parameters['min_snr'] * 0.6)
            elif peak.get('shadow_peak', False):
                min_snr_threshold = max(2.5, self.detection_parameters['min_snr'] * 0.8)
            else:
                min_snr_threshold = self.detection_parameters['min_snr']

            # Calculate local SNR
            local_noise = self.estimate_local_noise_around_peak(data_2d, peak)
            peak_snr = peak.get('intensity', 0) / (local_noise + 1e-10)

            if peak_snr >= min_snr_threshold:
                peak['adaptive_snr'] = peak_snr
                filtered_peaks.append(peak)

        return filtered_peaks

    def detect_overlapping_peaks_enhanced(self, peaks, data_2d):
        """Enhanced overlapping peak detection with better handling of complex scenarios"""
        if len(peaks) < 2:
            return peaks

        # Enhanced overlap detection using multiple criteria
        enhanced_peaks = []

        for i, peak in enumerate(peaks):
            peak_enhanced = peak.copy()

            # Check for nearby peaks with multiple distance thresholds
            nearby_peaks = []
            for j, other_peak in enumerate(peaks):
                if i != j:
                    distance = np.sqrt((peak['ppm_x'] - other_peak['ppm_x'])**2 +
                                     (peak['ppm_y'] - other_peak['ppm_y'])**2)

                    # Adaptive overlap threshold based on peak intensities and nucleus type
                    intensity_ratio = min(peak.get('intensity', 1), other_peak.get('intensity', 1)) / \
                                    max(peak.get('intensity', 1), other_peak.get('intensity', 1))

                    # Get nucleus-specific distance threshold
                    nucleus_type = peak.get('nucleus_type', '1H')
                    nucleus_distance_factor = self.get_nucleus_distance_factor(nucleus_type)

                    # Combined adaptive threshold
                    base_threshold = self.detection_parameters['overlap_threshold']
                    intensity_adjustment = (0.5 + 0.5 * intensity_ratio)
                    adaptive_threshold = base_threshold * intensity_adjustment * nucleus_distance_factor

                    if distance < adaptive_threshold:
                        nearby_peaks.append({'peak': other_peak, 'distance': distance, 'intensity_ratio': intensity_ratio})

            # Enhanced overlap classification
            if nearby_peaks:
                peak_enhanced['overlapping'] = True
                peak_enhanced['overlap_complexity'] = len(nearby_peaks)
                peak_enhanced['min_neighbor_distance'] = min(p['distance'] for p in nearby_peaks)
                peak_enhanced['intensity_variation'] = max(p['intensity_ratio'] for p in nearby_peaks)
            else:
                peak_enhanced['overlapping'] = False
                peak_enhanced['overlap_complexity'] = 0

            enhanced_peaks.append(peak_enhanced)

        return enhanced_peaks

    def get_adaptive_percentile_threshold(self, data_snr):
        """
        Calculate adaptive percentile threshold based on data quality
        """
        base_percentile = 95  # Default top 5%

        if data_snr > 50:  # High quality data
            return base_percentile + 2  # Top 3% (97th percentile)
        elif data_snr > 20:  # Good quality data
            return base_percentile + 1  # Top 4% (96th percentile)
        elif data_snr > 10:  # Fair quality data
            return base_percentile      # Top 5% (95th percentile)
        elif data_snr > 5:   # Poor quality data
            return base_percentile - 5  # Top 10% (90th percentile)
        else:  # Very noisy data
            return base_percentile - 15 # Top 20% (80th percentile)

    def get_adaptive_prominence_threshold(self, data, noise_level):
        """
        Calculate adaptive prominence threshold based on data characteristics
        """
        data_max = np.max(data)
        data_snr = data_max / (noise_level + 1e-10)

        # Base prominence factor from detection parameters
        base_prominence = self.detection_parameters['min_prominence_factor']

        # Adaptive prominence based on SNR
        if data_snr > 50:  # High SNR - can detect smaller prominences
            prominence_factor = base_prominence * 0.5
        elif data_snr > 20:  # Good SNR
            prominence_factor = base_prominence * 0.7
        elif data_snr > 10:  # Fair SNR
            prominence_factor = base_prominence * 1.0
        elif data_snr > 5:   # Poor SNR - need higher prominence
            prominence_factor = base_prominence * 1.5
        else:  # Very poor SNR
            prominence_factor = base_prominence * 2.0

        # Additional scaling based on data dynamic range
        data_range = np.max(data) - np.min(data)
        range_factor = max(0.5, min(2.0, data_range / (noise_level * 10)))

        adaptive_prominence = prominence_factor * range_factor

        # Ensure reasonable bounds
        adaptive_prominence = max(0.02, min(0.5, adaptive_prominence))  # 2% to 50%

        return adaptive_prominence

    def apply_adaptive_prominence_filter(self, data, local_max_mask, prominence_threshold):
        """
        Apply adaptive prominence filtering to local maxima
        """
        # Get coordinates of local maxima
        max_coords = np.where(local_max_mask)

        if len(max_coords[0]) == 0:
            return local_max_mask

        # Calculate prominence for each local maximum
        prominence_mask = np.zeros_like(data, dtype=bool)

        for i, (y, x) in enumerate(zip(max_coords[0], max_coords[1])):
            peak_value = data[y, x]

            # Define local neighborhood for prominence calculation
            neighborhood_size = 5
            y_start = max(0, y - neighborhood_size)
            y_end = min(data.shape[0], y + neighborhood_size + 1)
            x_start = max(0, x - neighborhood_size)
            x_end = min(data.shape[1], x + neighborhood_size + 1)

            local_region = data[y_start:y_end, x_start:x_end]

            # Calculate prominence as difference from local minimum
            local_min = np.min(local_region)
            prominence = peak_value - local_min

            # Calculate relative prominence
            relative_prominence = prominence / peak_value if peak_value > 0 else 0

            # Apply adaptive threshold
            if relative_prominence >= prominence_threshold:
                prominence_mask[y, x] = True

        return prominence_mask

    def get_adaptive_distance_threshold(self, nucleus_type, data_complexity=1.0):
        """
        Calculate adaptive distance threshold for shadow peak detection
        """
        # Get nucleus-specific parameters
        nmr_params = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])
        typical_width = nmr_params['typical_width']

        # Base distance threshold
        base_distance = 0.1  # Default 0.1 ppm

        # Nucleus-specific scaling
        if nucleus_type == '1H':
            nucleus_factor = 0.5    # Tighter for 1H due to narrow range
        elif nucleus_type == '15N':
            nucleus_factor = 2.0    # Wider for 15N due to broader range
        elif nucleus_type == '13C':
            nucleus_factor = 1.5    # Medium for 13C
        else:
            nucleus_factor = 1.0

        # Complexity-based adjustment
        complexity_factor = max(0.5, min(2.0, data_complexity))

        # Calculate adaptive distance
        adaptive_distance = base_distance * nucleus_factor * complexity_factor
        adaptive_distance = max(typical_width * 0.5, min(typical_width * 5.0, adaptive_distance))

        return adaptive_distance

    def get_nucleus_distance_factor(self, nucleus_type):
        """
        Get nucleus-specific distance scaling factor for overlap detection
        """
        if nucleus_type == '1H':
            return 0.5  # Tighter overlap detection for 1H (narrow chemical shift range)
        elif nucleus_type == '15N':
            return 2.0  # Wider overlap detection for 15N (broader range)
        elif nucleus_type == '13C':
            return 1.5  # Medium overlap detection for 13C
        else:
            return 1.0  # Default factor

    def validate_peak_by_fitting_enhanced(self, data_2d, ppm_x_axis, ppm_y_axis, peak):
        """Enhanced peak validation with improved criteria for complex scenarios"""
        # Use the existing validation but with enhanced parameters
        validated_peak = self.validate_peak_by_fitting(data_2d, ppm_x_axis, ppm_y_axis, peak)

        # Enhanced validation criteria
        if validated_peak.get('fit_quality', 0) > 0:
            # Adjust quality score based on peak characteristics
            base_quality = validated_peak['fit_quality']

            # Bonus for successfully fitting small peaks
            if peak.get('small_peak', False) and base_quality > 0.3:
                validated_peak['fit_quality'] = min(1.0, base_quality * 1.2)

            # Bonus for resolving overlapping peaks
            if peak.get('overlapping', False) and base_quality > 0.5:
                validated_peak['fit_quality'] = min(1.0, base_quality * 1.1)

            # Penalty for very complex overlaps that may be fitting artifacts
            if peak.get('overlap_complexity', 0) > 3 and base_quality < 0.7:
                validated_peak['fit_quality'] = base_quality * 0.9

        return validated_peak

    def get_adaptive_fit_threshold(self, peak, all_peaks):
        """Get adaptive fit quality threshold based on peak characteristics"""
        base_threshold = self.detection_parameters['validation_fit_threshold']

        # Lower threshold for small peaks (easier to accept)
        if peak.get('small_peak', False):
            return max(0.3, base_threshold * 0.6)

        # Lower threshold for overlapping peaks (more challenging to fit)
        if peak.get('overlapping', False):
            complexity = peak.get('overlap_complexity', 1)
            reduction_factor = max(0.5, 1.0 - 0.1 * complexity)
            return max(0.4, base_threshold * reduction_factor)

        return base_threshold

    def estimate_local_noise_around_peak(self, data_2d, peak):
        """Estimate noise level in the local region around a peak"""
        try:
            # Extract small region around peak for local noise estimation
            x_idx = peak.get('x_idx', 0)
            y_idx = peak.get('y_idx', 0)

            # Define local region (e.g., 10x10 pixels around peak)
            region_size = 10
            x_start = max(0, x_idx - region_size//2)
            x_end = min(data_2d.shape[1], x_idx + region_size//2)
            y_start = max(0, y_idx - region_size//2)
            y_end = min(data_2d.shape[0], y_idx + region_size//2)

            local_region = data_2d[y_start:y_end, x_start:x_end]

            # Estimate noise using robust statistics (excluding the peak itself)
            peak_mask = np.ones_like(local_region, dtype=bool)
            center_x, center_y = local_region.shape[1]//2, local_region.shape[0]//2
            peak_mask[center_y-2:center_y+3, center_x-2:center_x+3] = False  # Exclude peak center

            background_values = local_region[peak_mask]
            if len(background_values) > 5:
                return np.std(background_values)
            else:
                return np.std(local_region) * 0.5  # Conservative estimate
        except:
            # Fallback to global noise estimate
            return np.std(data_2d) * 0.1

    def export_peaks_to_dataframe(self, peaks):
        """Export peak list to pandas DataFrame for analysis"""
        if not peaks:
            return pd.DataFrame()

        # Extract key information for each peak
        peak_data = []
        for i, peak in enumerate(peaks):
            row = {
                'Peak_ID': i + 1,
                'Assignment': f"Peak_{i+1}",
                'X_ppm': peak.get('ppm_x_fitted', peak['ppm_x']),
                'Y_ppm': peak.get('ppm_y_fitted', peak['ppm_y']),
                'Intensity': peak['intensity'],
                'SNR': peak.get('snr', 0),
                'Quality_Score': peak.get('quality_score', 0),
                'Fit_Quality': peak.get('fit_quality', 0),
                'Overlapping': peak.get('overlapping', False),
                'Nucleus_Type': peak.get('nucleus_type', 'Unknown'),
                'Validation': peak.get('fit_validation', 'not_performed')
            }

            # Add fit parameters if available
            if 'x_fit_result' in peak:
                x_fit = peak['x_fit_result']
                row.update({
                    'X_Amplitude': x_fit.get('amplitude', 0),
                    'X_Sigma': x_fit.get('sigma', 0),
                    'X_Gamma': x_fit.get('gamma', 0),
                    'X_R_squared': x_fit.get('r_squared', 0)
                })

            if 'y_fit_result' in peak:
                y_fit = peak['y_fit_result']
                row.update({
                    'Y_Amplitude': y_fit.get('amplitude', 0),
                    'Y_Sigma': y_fit.get('sigma', 0),
                    'Y_Gamma': y_fit.get('gamma', 0),
                    'Y_R_squared': y_fit.get('r_squared', 0)
                })

            peak_data.append(row)

        return pd.DataFrame(peak_data)


if __name__ == "__main__":
    # Test the enhanced peak picker
    print("🧪 Testing Enhanced Peak Picker")
    print("=" * 40)

    # Create synthetic 2D NMR data with multiple peaks
    x_ppm = np.linspace(7.0, 9.0, 100)  # ¹H dimension
    y_ppm = np.linspace(110, 130, 80)   # ¹⁵N dimension

    # Create test data with 3 peaks + noise
    X, Y = np.meshgrid(x_ppm, y_ppm)

    # Synthetic peaks
    peak1 = 1000 * np.exp(-((X - 8.0)**2 / 0.01**2 + (Y - 120)**2 / 2**2))
    peak2 = 800 * np.exp(-((X - 7.5)**2 / 0.015**2 + (Y - 125)**2 / 1.5**2))
    peak3 = 600 * np.exp(-((X - 8.5)**2 / 0.02**2 + (Y - 115)**2 / 3**2))

    # Add noise and baseline
    noise = np.random.normal(0, 100, X.shape)
    baseline = 200

    synthetic_data = peak1 + peak2 + peak3 + noise + baseline

    # Test peak detection
    picker = EnhancedPeakPicker()
    peaks = picker.comprehensive_peak_detection(synthetic_data, x_ppm, y_ppm, '1H', validate_fits=False)

    print(f"✅ Detected {len(peaks)} peaks")
    for i, peak in enumerate(peaks[:3]):  # Show first 3 peaks
        print(f"   Peak {i+1}: ({peak['ppm_x']:.3f}, {peak['ppm_y']:.1f}) ppm, SNR={peak['snr']:.1f}")

    # Test DataFrame export
    df = picker.export_peaks_to_dataframe(peaks)
    print(f"✅ Exported to DataFrame: {len(df)} rows × {len(df.columns)} columns")
