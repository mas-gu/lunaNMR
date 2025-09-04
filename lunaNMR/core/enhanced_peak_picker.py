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
            'validation_fit_threshold': 0.5  # Minimum R² for peak validation
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
            print(f"Noise estimation failed: {e}")
            return np.std(data_2d) * 0.1  # Fallback

    def preprocess_data_for_detection(self, data_2d, sigma=None):
        """
        Preprocess 2D data for better peak detection
        """
        if sigma is None:
            sigma = self.detection_parameters['smoothing_sigma']

        # Light Gaussian smoothing to reduce noise
        if sigma > 0:
            smoothed = gaussian_filter(data_2d, sigma=sigma)
        else:
            smoothed = data_2d

        return smoothed

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

        # Calculate dynamic thresholds
        min_snr = nmr_params['min_snr']
        height_threshold = noise_level * min_snr

        # Additional threshold based on data percentiles
        percentile_threshold = np.percentile(smoothed_data, 95)  # Top 5% of intensities

        # Use the higher of SNR-based and percentile-based thresholds
        final_threshold = max(height_threshold, percentile_threshold * 0.1)

        # Find local maxima using maximum filter
        neighborhood_size = max(3, min(data_2d.shape) // 50)  # Adaptive neighborhood
        local_max_mask = (smoothed_data == maximum_filter(smoothed_data, size=neighborhood_size))

        # Apply intensity threshold
        peak_candidates = local_max_mask & (smoothed_data > final_threshold)

        # Extract peak positions
        peak_indices = np.where(peak_candidates)

        if len(peak_indices[0]) == 0:
            return []

        # Convert to peak list with coordinates and intensities
        initial_peaks = []
        for i, (y_idx, x_idx) in enumerate(zip(peak_indices[0], peak_indices[1])):
            peak_info = {
                'id': i,
                'x_idx': x_idx,
                'y_idx': y_idx,
                'ppm_x': ppm_x_axis[x_idx],
                'ppm_y': ppm_y_axis[y_idx],
                'intensity': smoothed_data[y_idx, x_idx],
                'snr': smoothed_data[y_idx, x_idx] / noise_level,
                'nucleus_type': nucleus_type
            }
            initial_peaks.append(peak_info)

        # Sort by intensity (highest first)
        initial_peaks.sort(key=lambda x: x['intensity'], reverse=True)

        return initial_peaks

    def filter_peaks_by_snr(self, peaks, min_snr=None):
        """Filter peaks based on signal-to-noise ratio"""
        if min_snr is None:
            min_snr = self.detection_parameters['min_snr']

        filtered_peaks = [p for p in peaks if p.get('snr', 0) >= min_snr]

        print(f"SNR filtering: {len(peaks)} → {len(filtered_peaks)} peaks (SNR ≥ {min_snr})")
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
                print(f"Clustering failed: {e}")
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
            print(f"Fit validation failed for peak at ({peak['ppm_x']:.3f}, {peak['ppm_y']:.1f}): {e}")
            peak['fit_validation'] = 'error'
            peak['fit_quality'] = 0.0

        return peak

    def comprehensive_peak_detection(self, data_2d, ppm_x_axis, ppm_y_axis,
                                   nucleus_type=None, validate_fits=True):
        """
        Comprehensive peak detection with all enhancements

        Steps:
        1. Initial peak detection using local maxima
        2. SNR-based filtering
        3. Overlapping peak detection
        4. Optional fit-based validation
        5. Final quality assessment
        """
        print(f"🔍 Starting comprehensive peak detection...")

        # Step 1: Initial detection
        initial_peaks = self.detect_peaks_2d_initial(data_2d, ppm_x_axis, ppm_y_axis, nucleus_type)
        print(f"   Initial detection: {len(initial_peaks)} peaks found")

        if len(initial_peaks) == 0:
            return []

        # Step 2: SNR filtering
        snr_filtered = self.filter_peaks_by_snr(initial_peaks)

        # Step 3: Overlapping peak detection
        overlap_analyzed = self.detect_overlapping_peaks(snr_filtered)

        # Step 4: Optional fit validation
        if validate_fits and len(overlap_analyzed) > 0:
            print(f"   Validating peaks by Voigt fitting...")
            validated_peaks = []

            for i, peak in enumerate(overlap_analyzed):
                if i % max(1, len(overlap_analyzed) // 10) == 0:  # Progress indicator
                    print(f"   Validating peak {i+1}/{len(overlap_analyzed)}")

                validated_peak = self.validate_peak_by_fitting(data_2d, ppm_x_axis, ppm_y_axis, peak)

                # Only keep peaks that pass fit validation
                min_fit_quality = self.detection_parameters['validation_fit_threshold']
                if validated_peak.get('fit_quality', 0) >= min_fit_quality:
                    validated_peaks.append(validated_peak)

            final_peaks = validated_peaks
            print(f"   Fit validation: {len(overlap_analyzed)} → {len(final_peaks)} peaks")
        else:
            final_peaks = overlap_analyzed

        # Step 5: Final quality assessment and sorting
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

        print(f"✅ Peak detection complete: {len(final_peaks)} high-quality peaks")
        return final_peaks

    def get_detection_stats(self):
        """Return statistics from last detection"""
        return self.last_detection_stats.copy()

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


# Convenience function
def detect_peaks_enhanced(data_2d, ppm_x_axis, ppm_y_axis, nucleus_type=None, validate_fits=True):
    """
    Convenience function for enhanced peak detection

    Returns both peak list and DataFrame
    """
    picker = EnhancedPeakPicker()
    peaks = picker.comprehensive_peak_detection(data_2d, ppm_x_axis, ppm_y_axis, nucleus_type, validate_fits)
    df = picker.export_peaks_to_dataframe(peaks)

    return peaks, df


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
