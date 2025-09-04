#!/usr/bin/env python3
"""
Advanced NMR Peak Integration Script - Enhanced Reference-Based Detection

This script implements professional-grade NMR peak integration methods
with enhanced reference-based peak detection capabilities.

Key enhancements over original version:
1. Reference-guided peak detection using known peak positions
2. User-configurable search window parameters (±X ppm, ±Y intensity)
3. User-selectable noise region functionality
4. Adjustable noise threshold multiplier (0.1-10.0 range)
5. Peak propagation for series analysis
6. Improved error handling and validation

Author: Guillaume Mas
Date: 2025
"""

# ========================================================================
# USER CONFIGURATION - MODIFY THESE PARAMETERS AS NEEDED
# ========================================================================

# Input files and directories
PEAK_LIST_FILE = "peak_list/600_T1_KRASB_Q61L_0o0.txt"    # Path to peak list file
NMR_DATA_PATH = "data/600_T1_KRASB_Q61L_0o0.ft"          # Path to NMR data file or directory

# Output settings
OUTPUT_DIR = "inplace_results"                      # Output directory for results
CREATE_VISUALIZATION = True                          # Create visualization plots
SHOW_PLOTS = False                                   # Show plots interactively
SAVE_PPM_PLOTS = True                               # Save plots with PPM axes

# Enhanced Processing parameters
NOISE_THRESHOLD_MULTIPLIER = 3.0                   # Threshold = noise_level * multiplier
SEARCH_RADIUS_PPM = 0.2                            # Search radius in ppm (1H dimension)
SEARCH_RADIUS_PPM_INDIRECT = 3.0                   # Search radius in ppm (15N/13C dimension)
QUALITY_THRESHOLD = 0.7                            # Quality threshold for peak fitting (0-1)

# Reference-based detection parameters
USE_REFERENCE_DETECTION = True                      # Enable reference-based peak detection
SEARCH_WINDOW_X_PPM = 0.2                         # ±X ppm search window (1H dimension)
SEARCH_WINDOW_Y_PPM = 3.0                         # ±Y ppm search window (15N/13C dimension)
RETAIN_REFERENCE_IF_NOT_FOUND = True              # Keep original position if no peak detected

# Noise region configuration
USE_CUSTOM_NOISE_REGIONS = False                   # Use user-defined noise regions
NOISE_REGIONS_PPM = [                             # Noise regions in ppm [x_min, x_max, y_min, y_max]
    [0.5, 1.5, 100, 110],                        # Example region 1
    [11, 12, 130, 140]                           # Example region 2
]

# Threshold control parameters
MIN_THRESHOLD_MULTIPLIER = 0.1                    # Minimum threshold multiplier
MAX_THRESHOLD_MULTIPLIER = 10.0                   # Maximum threshold multiplier
DEFAULT_THRESHOLD_MULTIPLIER = 3.0                # Default threshold multiplier

# Peak fitting parameters
FIT_VOIGT_PROFILES = True                          # Use Voigt profile fitting
MAX_FITTING_ITERATIONS = 100                       # Maximum iterations for fitting
FITTING_TOLERANCE = 1e-6                           # Convergence tolerance

# Visualization options
FIGURE_SIZE = (14, 10)                             # Figure size (width, height)
DPI = 300                                          # Resolution for saved plots
CONTOUR_LEVELS = 20                               # Number of contour levels
COLORMAP = 'viridis'                              # Colormap for visualization
SHOW_PEAK_LABELS = True                           # Show peak assignment labels
LABEL_FONT_SIZE = 8                               # Font size for peak labels

# Quality control
MINIMUM_SNR = 3.0                                 # Minimum signal-to-noise ratio
VALIDATE_PEAK_SHAPES = True                       # Validate peak shape quality
REMOVE_ARTIFACTS = True                           # Remove obvious artifacts

# ========================================================================
# END USER CONFIGURATION
# ========================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import sys
from scipy.optimize import curve_fit, minimize
from scipy.ndimage import maximum_filter, minimum_filter
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

try:
    import nmrglue as ng
except ImportError:
    print("Error: nmrglue not installed. Install with: pip install nmrglue")
    sys.exit(1)

class InPlaceAdvancedNMRIntegrator:
    def __init__(self):
        """Initialize the Enhanced Advanced NMR Integrator"""
        self.peak_list = None
        self.nmr_data = None
        self.nmr_dict = None
        self.noise_level = None
        self.threshold = None

        # PPM axes
        self.ppm_x_axis = None
        self.ppm_y_axis = None

        # Results storage
        self.detected_peaks = []
        self.fitted_peaks = []
        self.integration_results = []

        # Enhanced parameters
        self.search_window_x_ppm = SEARCH_WINDOW_X_PPM
        self.search_window_y_ppm = SEARCH_WINDOW_Y_PPM
        self.threshold_multiplier = DEFAULT_THRESHOLD_MULTIPLIER
        self.noise_regions_ppm = NOISE_REGIONS_PPM
        self.use_custom_noise_regions = USE_CUSTOM_NOISE_REGIONS
        self.use_reference_detection = USE_REFERENCE_DETECTION
        self.consolidation_x_tolerance = 0.05
        self.consolidation_y_tolerance = 2.0

        # Quality tracking
        self.detection_statistics = {
            'total_peaks': 0,
            'detected_peaks': 0,
            'reference_retained': 0,
            'detection_rate': 0.0
        }

    def set_search_window(self, x_ppm, y_ppm):
        """Set search window parameters

        Args:
            x_ppm (float): ±X ppm search window (1H dimension)
            y_ppm (float): ±Y ppm search window (15N/13C dimension)
        """
        if x_ppm <= 0 or y_ppm <= 0:
            raise ValueError("Search window parameters must be positive")
        self.search_window_x_ppm = x_ppm
        self.search_window_y_ppm = y_ppm
        print(f"Search window set to ±{x_ppm} ppm (1H) × ±{y_ppm} ppm (15N/13C)")

    def set_threshold_multiplier(self, multiplier):
        """Set noise threshold multiplier

        Args:
            multiplier (float): Threshold multiplier (0.1-10.0)
        """
        if not (MIN_THRESHOLD_MULTIPLIER <= multiplier <= MAX_THRESHOLD_MULTIPLIER):
            raise ValueError(f"Threshold multiplier must be between {MIN_THRESHOLD_MULTIPLIER} and {MAX_THRESHOLD_MULTIPLIER}")
        self.threshold_multiplier = multiplier
        if self.noise_level is not None:
            self.threshold = self.noise_level * self.threshold_multiplier
        print(f"Threshold multiplier set to {multiplier}")

    def set_noise_regions(self, regions_ppm):
        """Set custom noise regions

        Args:
            regions_ppm (list): List of [x_min, x_max, y_min, y_max] regions in ppm
        """
        if not isinstance(regions_ppm, list) or len(regions_ppm) == 0:
            raise ValueError("Noise regions must be a non-empty list")

        for region in regions_ppm:
            if len(region) != 4:
                raise ValueError("Each noise region must have 4 values: [x_min, x_max, y_min, y_max]")
            if region[0] >= region[1] or region[2] >= region[3]:
                raise ValueError("Invalid region bounds: min values must be less than max values")

        self.noise_regions_ppm = regions_ppm
        self.use_custom_noise_regions = True
        print(f"Custom noise regions set: {len(regions_ppm)} regions")

    def load_data(self, peak_list_file, nmr_file):
        """Load peak list and NMR data"""
        # Load peak list
        try:
            self.peak_list = pd.read_csv(peak_list_file, skipinitialspace=True)
            self.peak_list.columns = self.peak_list.columns.str.strip()
            print(f"Loaded peak list: {len(self.peak_list)} peaks from {peak_list_file}")
        except Exception as e:
            print(f"Error loading peak list: {e}")
            return False

        # Load NMR data
        try:
            self.nmr_dict, self.nmr_data = ng.pipe.read(nmr_file)
            print(f"Loaded NMR data: {self.nmr_data.shape} from {nmr_file}")

            # Calculate PPM axes
            self._calculate_ppm_axes()

            # Estimate noise level
            self._estimate_noise_level()

            return True
        except Exception as e:
            print(f"Error loading NMR data: {e}")
            return False

    def _calculate_ppm_axes(self):
        """Calculate proper PPM axes using nmrglue's built-in conversion"""
        try:
            # Use nmrglue's proper PPM conversion which handles NMRPipe parameters correctly
            uc_f2 = ng.pipe.make_uc(self.nmr_dict, self.nmr_data, dim=1)  # F2 dimension (1H)
            uc_f1 = ng.pipe.make_uc(self.nmr_dict, self.nmr_data, dim=0)  # F1 dimension (15N/13C)

            # Get the full PPM scales
            f2_ppm_scale = uc_f2.ppm_scale()
            f1_ppm_scale = uc_f1.ppm_scale()

            # Create PPM axes (nmrglue handles all the complex NMRPipe parameter conversion)
            self.ppm_x_axis = f2_ppm_scale
            self.ppm_y_axis = f1_ppm_scale

            print(f"F2 (1H) range: {f2_ppm_scale[0]:.2f} to {f2_ppm_scale[-1]:.2f} ppm")
            print(f"F1 (15N/13C) range: {f1_ppm_scale[0]:.1f} to {f1_ppm_scale[-1]:.1f} ppm")

        except Exception as e:
            print(f"Warning: Could not calculate PPM axes properly: {e}")
            # Fallback to simple linear axes
            self.ppm_x_axis = np.linspace(12, 0, self.nmr_data.shape[1])
            self.ppm_y_axis = np.linspace(140, 100, self.nmr_data.shape[0])

    def _estimate_noise_level(self):
        """Estimate noise level in the spectrum using custom or default regions"""
        if self.use_custom_noise_regions and self.noise_regions_ppm:
            # Use custom noise regions
            noise_values = []

            for region in self.noise_regions_ppm:
                x_min_ppm, x_max_ppm, y_min_ppm, y_max_ppm = region

                # Convert ppm to indices
                x_min_idx = np.argmin(np.abs(self.ppm_x_axis - x_max_ppm))  # Note: reversed for NMR
                x_max_idx = np.argmin(np.abs(self.ppm_x_axis - x_min_ppm))
                y_min_idx = np.argmin(np.abs(self.ppm_y_axis - y_max_ppm))  # Note: reversed for NMR
                y_max_idx = np.argmin(np.abs(self.ppm_y_axis - y_min_ppm))

                # Ensure proper ordering
                x_min_idx, x_max_idx = min(x_min_idx, x_max_idx), max(x_min_idx, x_max_idx)
                y_min_idx, y_max_idx = min(y_min_idx, y_max_idx), max(y_min_idx, y_max_idx)

                # Extract region data
                region_data = self.nmr_data[y_min_idx:y_max_idx, x_min_idx:x_max_idx]
                noise_values.extend(region_data.flatten())

                print(f"  Added noise region: {x_min_ppm}-{x_max_ppm} ppm (1H), {y_min_ppm}-{y_max_ppm} ppm (15N/13C)")

            if noise_values:
                self.noise_level = np.std(noise_values)
                print(f"Custom noise estimation from {len(self.noise_regions_ppm)} regions")
            else:
                print("Warning: No valid custom noise regions found, falling back to corner method")
                self._estimate_noise_level_corners()
        else:
            # Use default corner method
            self._estimate_noise_level_corners()

        self.threshold = self.noise_level * self.threshold_multiplier
        print(f"Estimated noise level: {self.noise_level:.2f}")
        print(f"Detection threshold: {self.threshold:.2f} (multiplier: {self.threshold_multiplier})")

    def _estimate_noise_level_corners(self):
        """Estimate noise level using left corner regions (fallback method)"""
        # Use left corner regions to estimate noise
        corner_size = min(50, self.nmr_data.shape[0]//10, self.nmr_data.shape[1]//10)

        corners = [
            self.nmr_data[:corner_size, :corner_size],  # Top-left
            self.nmr_data[-corner_size:, :corner_size], # Bottom-left
        ]

        noise_values = []
        for corner in corners:
            noise_values.extend(corner.flatten())

        self.noise_level = np.std(noise_values)

    def ppm_to_point(self, ppm_x, ppm_y):
        """Convert PPM coordinates to array indices"""
        # Find closest indices
        x_idx = np.argmin(np.abs(self.ppm_x_axis - ppm_x))
        y_idx = np.argmin(np.abs(self.ppm_y_axis - ppm_y))

        # Ensure within bounds
        x_idx = max(0, min(x_idx, len(self.ppm_x_axis) - 1))
        y_idx = max(0, min(y_idx, len(self.ppm_y_axis) - 1))

        return y_idx, x_idx

    def point_to_ppm(self, y_point, x_point):
        """Convert array indices to PPM coordinates"""
        ppm_x = self.ppm_x_axis[x_point]
        ppm_y = self.ppm_y_axis[y_point]
        return ppm_x, ppm_y

    def detect_peaks_reference_based(self):
        """Enhanced reference-based peak detection with Y-peak consolidation."""
        if self.peak_list is None:
            raise ValueError("Peak list must be loaded before detection")
        if self.nmr_data is None:
            raise ValueError("NMR data must be loaded before detection")

        print("Performing reference-based peak detection with Y-consolidation...")

        detected_peaks = []
        reference_retained = 0

        for idx, peak_row in self.peak_list.iterrows():
            ref_x_ppm = float(peak_row['Position_X'])
            ref_y_ppm = float(peak_row['Position_Y'])
            assignment = peak_row.get('Assignment', f'Peak_{idx+1}')

            # Define search window around reference position
            x_min_ppm = ref_x_ppm - self.search_window_x_ppm
            x_max_ppm = ref_x_ppm + self.search_window_x_ppm
            y_min_ppm = ref_y_ppm - self.search_window_y_ppm
            y_max_ppm = ref_y_ppm + self.search_window_y_ppm

            # Convert to array indices
            x_min_idx = np.argmin(np.abs(self.ppm_x_axis - x_max_ppm))
            x_max_idx = np.argmin(np.abs(self.ppm_x_axis - x_min_ppm))
            y_min_idx = np.argmin(np.abs(self.ppm_y_axis - y_max_ppm))
            y_max_idx = np.argmin(np.abs(self.ppm_y_axis - y_min_ppm))

            # Ensure proper bounds
            x_min_idx = max(0, min(x_min_idx, x_max_idx))
            x_max_idx = min(self.nmr_data.shape[1], max(x_min_idx, x_max_idx))
            y_min_idx = max(0, min(y_min_idx, y_max_idx))
            y_max_idx = min(self.nmr_data.shape[0], max(y_min_idx, y_max_idx))

            # Extract search window
            search_window = self.nmr_data[y_min_idx:y_max_idx, x_min_idx:x_max_idx]

            # --- MODIFICATION START ---
            # 1. Find ALL candidate peaks in the window
            candidate_peaks = self._find_all_peaks_in_window(
                search_window, x_min_idx, y_min_idx,
                ref_x_ppm, ref_y_ppm, assignment
            )

            if candidate_peaks:
                # 2. Consolidate the found peaks
                consolidated_peaks = self._consolidate_y_dimension_peaks(
                    candidate_peaks, x_tolerance=self.consolidation_x_tolerance, y_tolerance=self.consolidation_y_tolerance
                )

                # 3. Select the best peak from the consolidated list (closest to reference)
                min_dist = float('inf')
                best_peak = None
                for peak in consolidated_peaks:
                    dist = np.sqrt((peak['ppm_x'] - ref_x_ppm)**2 + (peak['ppm_y'] - ref_y_ppm)**2)
                    if dist < min_dist:
                        min_dist = dist
                        best_peak = peak

                if best_peak:
                    best_peak['distance_from_reference'] = min_dist
                    detected_peaks.append(best_peak)
                else:
                    # This case should not happen if consolidated_peaks is not empty, but as a fallback:
                    if RETAIN_REFERENCE_IF_NOT_FOUND:
                        ref_y_idx, ref_x_idx = self.ppm_to_point(ref_x_ppm, ref_y_ppm)
                        ref_intensity = self.nmr_data[ref_y_idx, ref_x_idx]
                        ref_snr = abs(ref_intensity) / self.noise_level if self.noise_level > 0 else 0
                        retained_peak = {
                            'assignment': assignment, 'ppm_x': ref_x_ppm, 'ppm_y': ref_y_ppm,
                            'x_point': ref_x_idx, 'y_point': ref_y_idx, 'intensity': ref_intensity,
                            'snr': ref_snr, 'detected': False, 'reference_retained': True,
                            'detection_quality': 'Reference'
                        }
                        detected_peaks.append(retained_peak)
                        reference_retained += 1

            else:
                # No peaks found in window, retain reference if desired
                if RETAIN_REFERENCE_IF_NOT_FOUND:
                    ref_y_idx, ref_x_idx = self.ppm_to_point(ref_x_ppm, ref_y_ppm)
                    ref_intensity = self.nmr_data[ref_y_idx, ref_x_idx]
                    ref_snr = abs(ref_intensity) / self.noise_level if self.noise_level > 0 else 0
                    retained_peak = {
                        'assignment': assignment, 'ppm_x': ref_x_ppm, 'ppm_y': ref_y_ppm,
                        'x_point': ref_x_idx, 'y_point': ref_y_idx, 'intensity': ref_intensity,
                        'snr': ref_snr, 'detected': False, 'reference_retained': True,
                        'detection_quality': 'Reference'
                    }
                    detected_peaks.append(retained_peak)
                    reference_retained += 1
            # --- MODIFICATION END ---

        # Update statistics
        self.detection_statistics = {
            'total_peaks': len(self.peak_list),
            'detected_peaks': sum(1 for p in detected_peaks if p['detected']),
            'reference_retained': reference_retained,
            'detection_rate': sum(1 for p in detected_peaks if p['detected']) / len(self.peak_list) * 100
        }

        # Apply Peak Centroid Refinement (always enabled for enhanced accuracy)
        gui_params = getattr(self, 'gui_params', {})
        detected_peaks = self._refine_peaks_with_centroids(detected_peaks, gui_params)

        # Store refined peaks (not original unrefined peaks)
        self.fitted_peaks = detected_peaks

        print(f"Reference-based detection with consolidation completed:")
        print(f"  Total peaks: {self.detection_statistics['total_peaks']}")
        print(f"  Detected peaks: {self.detection_statistics['detected_peaks']}")
        print(f"  Reference retained: {self.detection_statistics['reference_retained']}")
        print(f"  Detection rate: {self.detection_statistics['detection_rate']:.1f}%")

        return detected_peaks

    def _find_best_peak_in_window(self, search_window, x_offset, y_offset, ref_x_ppm, ref_y_ppm, assignment):
        """Find the best peak candidate within a search window"""
        if search_window.size == 0:
            return self._create_failed_detection(ref_x_ppm, ref_y_ppm, assignment, "Empty search window")

        # Find points above threshold in window
        abs_window = np.abs(search_window)
        mask = abs_window > self.threshold

        if not np.any(mask):
            return self._create_failed_detection(ref_x_ppm, ref_y_ppm, assignment, "No points above threshold")

        # Find local maxima
        candidates = []
        coords = np.argwhere(mask)

        for coord in coords:
            local_y, local_x = coord
            global_y = local_y + y_offset
            global_x = local_x + x_offset

            # Check if this is a local maximum
            if self._is_local_maximum(search_window, local_y, local_x):
                intensity = search_window[local_y, local_x]
                ppm_x, ppm_y = self.point_to_ppm(global_y, global_x)
                snr = abs(intensity) / self.noise_level if self.noise_level > 0 else 0

                # Calculate distance from reference
                distance = np.sqrt((ppm_x - ref_x_ppm)**2 + (ppm_y - ref_y_ppm)**2)

                candidates.append({
                    'ppm_x': ppm_x,
                    'ppm_y': ppm_y,
                    'x_point': global_x,
                    'y_point': global_y,
                    'intensity': intensity,
                    'snr': snr,
                    'distance': distance,
                    'quality_score': snr * abs(intensity) / (distance + 0.1)  # Combined score
                })

        if not candidates:
            return self._create_failed_detection(ref_x_ppm, ref_y_ppm, assignment, "No local maxima found")

        # Select best candidate (highest quality score)
        best_candidate = max(candidates, key=lambda x: x['quality_score'])

        # Validate minimum SNR
        if best_candidate['snr'] < MINIMUM_SNR:
            return self._create_failed_detection(ref_x_ppm, ref_y_ppm, assignment, f"SNR too low: {best_candidate['snr']:.1f}")

        # Create successful detection
        return {
            'assignment': assignment,
            'ppm_x': best_candidate['ppm_x'],
            'ppm_y': best_candidate['ppm_y'],
            'x_point': best_candidate['x_point'],
            'y_point': best_candidate['y_point'],
            'intensity': best_candidate['intensity'],
            'snr': best_candidate['snr'],
            'detected': True,
            'reference_retained': False,
            'detection_quality': 'Excellent' if best_candidate['snr'] > 10 else
                               'Good' if best_candidate['snr'] > 5 else 'Fair',
            'distance_from_reference': best_candidate['distance'],
            'quality_score': best_candidate['quality_score']
        }

    def _is_local_maximum(self, window, y, x):
        """Check if a point is a local maximum within a 3x3 neighborhood"""
        h, w = window.shape

        # Define 3x3 neighborhood
        y_min = max(0, y - 1)
        y_max = min(h, y + 2)
        x_min = max(0, x - 1)
        x_max = min(w, x + 2)

        neighborhood = np.abs(window[y_min:y_max, x_min:x_max])
        center_value = abs(window[y, x])

        return center_value >= np.max(neighborhood)

    def _create_failed_detection(self, ref_x_ppm, ref_y_ppm, assignment, reason):
        """Create a failed detection result"""
        return {
            'assignment': assignment,
            'ppm_x': ref_x_ppm,
            'ppm_y': ref_y_ppm,
            'x_point': -1,
            'y_point': -1,
            'intensity': 0.0,
            'snr': 0.0,
            'detected': False,
            'reference_retained': False,
            'detection_quality': 'Failed',
            'failure_reason': reason
        }

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

        # Group peaks by similar X-coordinates (1H dimension)
        x_groups = []

        for peak in peak_info:
            x_pos = peak.get('ppm_x', 0)

            # Find existing group with similar X-coordinate
            assigned_to_group = False
            for group in x_groups:
                group_x_positions = [p.get('ppm_x', 0) for p in group]
                avg_x = sum(group_x_positions) / len(group_x_positions)

                if abs(x_pos - avg_x) <= x_tolerance:
                    group.append(peak)
                    assigned_to_group = True
                    break

            # Create new group if no existing group found
            if not assigned_to_group:
                x_groups.append([peak])

        # For each X-group, consolidate Y-dimension fragments
        consolidated_peaks = []

        for group_idx, group in enumerate(x_groups):
            if len(group) == 1:
                # Single peak in group - no consolidation needed
                consolidated_peaks.extend(group)
                continue

            # Multiple peaks in same X-region - need consolidation
            # Sort peaks in group by Y-coordinate
            group_sorted = sorted(group, key=lambda p: p.get('ppm_y', 0))

            # Sub-group peaks that are close in Y-dimension
            y_subgroups = []
            current_subgroup = [group_sorted[0]]

            for i in range(1, len(group_sorted)):
                current_peak = group_sorted[i]
                previous_peak = group_sorted[i-1]

                y_current = current_peak.get('ppm_y', 0)
                y_previous = previous_peak.get('ppm_y', 0)

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
                strongest_peak = max(subgroup, key=lambda p: p.get('intensity', 0))
                consolidated_peaks.append(strongest_peak)

        return consolidated_peaks

    def _refine_peaks_with_centroids(self, peak_info, gui_params):
        """
        Peak Centroid Detection - Optional post-processing enhancement
        Calculate intensity-weighted centroids for improved coordinate accuracy

        Args:
            peak_info: List of detected peaks with position information
            gui_params: GUI parameters containing centroid settings

        Returns:
            List of peaks with refined centroid coordinates
        """
        if not peak_info:
            return peak_info

        #print(f"         🎯 Peak Centroid Detection: Refining {len(peak_info)} peaks...")

        # Parameters for centroid calculation
        centroid_window_x_ppm = gui_params.get('centroid_window_x_ppm', 0.01)  # window size in ppm (1H dimension)
        centroid_window_y_ppm = gui_params.get('centroid_window_y_ppm', 0.1)   # window size in ppm (15N dimension)
        noise_multiplier = gui_params.get('centroid_noise_multiplier', 5.0)    # noise threshold multiplier

        # Calculate noise level for threshold
        noise_level = np.std(np.abs(self.nmr_data[self.nmr_data != 0]))  # Exclude zeros, use absolute values
        noise_threshold = noise_level * noise_multiplier

        refined_peaks = []
        successful_refinements = 0

        for peak in peak_info:
            try:
                # Skip non-detected peaks
                if not peak.get('detected', False):
                    refined_peaks.append(peak)
                    continue

                # Get original coordinates - handle different field naming conventions
                original_x = peak.get('ppm_x', peak.get('Position_X', 0))
                original_y = peak.get('ppm_y', peak.get('Position_Y', 0))

                # Convert ppm coordinates to matrix indices
                x_idx = np.argmin(np.abs(self.ppm_x_axis - original_x))
                y_idx = np.argmin(np.abs(self.ppm_y_axis - original_y))

                # Convert ppm window sizes to pixel counts
                x_ppm_per_pixel = abs(self.ppm_x_axis[-1] - self.ppm_x_axis[0]) / len(self.ppm_x_axis)
                y_ppm_per_pixel = abs(self.ppm_y_axis[-1] - self.ppm_y_axis[0]) / len(self.ppm_y_axis)
                centroid_window_x_pixels = int(np.ceil(centroid_window_x_ppm / x_ppm_per_pixel))
                centroid_window_y_pixels = int(np.ceil(centroid_window_y_ppm / y_ppm_per_pixel))

                # Define window boundaries (ensure within matrix bounds)
                x_min = max(0, x_idx - centroid_window_x_pixels)
                x_max = min(self.nmr_data.shape[1], x_idx + centroid_window_x_pixels + 1)
                y_min = max(0, y_idx - centroid_window_y_pixels)
                y_max = min(self.nmr_data.shape[0], y_idx + centroid_window_y_pixels + 1)

                # Extract intensity window
                intensity_window = np.abs(self.nmr_data[y_min:y_max, x_min:x_max])  # Use absolute values

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
                if centroid_x_idx < len(self.ppm_x_axis) - 1 and centroid_x_idx >= 0:
                    # Linear interpolation for sub-pixel accuracy
                    x_floor = int(np.floor(centroid_x_idx))
                    x_ceil = int(np.ceil(centroid_x_idx))
                    x_frac = centroid_x_idx - x_floor

                    if x_floor == x_ceil:
                        refined_x_ppm = self.ppm_x_axis[x_floor]
                    else:
                        refined_x_ppm = self.ppm_x_axis[x_floor] * (1 - x_frac) + self.ppm_x_axis[x_ceil] * x_frac
                else:
                    refined_x_ppm = original_x

                if centroid_y_idx < len(self.ppm_y_axis) - 1 and centroid_y_idx >= 0:
                    y_floor = int(np.floor(centroid_y_idx))
                    y_ceil = int(np.ceil(centroid_y_idx))
                    y_frac = centroid_y_idx - y_floor

                    if y_floor == y_ceil:
                        refined_y_ppm = self.ppm_y_axis[y_floor]
                    else:
                        refined_y_ppm = self.ppm_y_axis[y_floor] * (1 - y_frac) + self.ppm_y_axis[y_ceil] * y_frac
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

    def _find_all_peaks_in_window(self, search_window, x_offset, y_offset, ref_x_ppm, ref_y_ppm, assignment):
        """Find all peak candidates within a search window."""
        if search_window.size == 0:
            return []

        # Find points above threshold in window
        abs_window = np.abs(search_window)
        mask = abs_window > self.threshold

        if not np.any(mask):
            return []

        # Find local maxima
        candidates = []
        coords = np.argwhere(mask)

        for coord in coords:
            local_y, local_x = coord
            global_y = local_y + y_offset
            global_x = local_x + x_offset

            # Check if this is a local maximum
            if self._is_local_maximum(search_window, local_y, local_x):
                intensity = search_window[local_y, local_x]
                ppm_x, ppm_y = self.point_to_ppm(global_y, global_x)
                snr = abs(intensity) / self.noise_level if self.noise_level > 0 else 0

                # Validate minimum SNR
                if snr >= MINIMUM_SNR:
                    candidates.append({
                        'ppm_x': ppm_x,
                        'ppm_y': ppm_y,
                        'x_point': global_x,
                        'y_point': global_y,
                        'intensity': intensity,
                        'snr': snr,
                        'assignment': assignment,
                        'detected': True,
                        'reference_retained': False,
                        'detection_quality': 'Excellent' if snr > 10 else 'Good' if snr > 5 else 'Fair'
                    })

        return candidates

    def detect_peaks_professional(self):
        """Professional peak detection - enhanced version with reference option"""
        if self.use_reference_detection and self.peak_list is not None:
            return self.detect_peaks_reference_based()
        else:
            return self.detect_peaks_traditional()

    def detect_peaks_traditional(self):
        """Traditional peak detection using threshold and local maxima (fallback)"""
        print("Detecting peaks using traditional methods...")

        # Find all points above threshold
        mask = np.abs(self.nmr_data) > self.threshold
        coords = np.argwhere(mask)

        detected = []

        # Local maxima detection
        for coord in coords:
            y, x = coord

            # Define local window
            window_size = 3
            y_min = max(0, y - window_size)
            y_max = min(self.nmr_data.shape[0], y + window_size + 1)
            x_min = max(0, x - window_size)
            x_max = min(self.nmr_data.shape[1], x + window_size + 1)

            # Extract local window
            window = np.abs(self.nmr_data[y_min:y_max, x_min:x_max])

            # Check if current point is local maximum
            max_pos = np.unravel_index(np.argmax(window), window.shape)
            if max_pos[0] == y - y_min and max_pos[1] == x - x_min:
                # This is a local maximum
                intensity = self.nmr_data[y, x]
                ppm_x, ppm_y = self.point_to_ppm(y, x)
                snr = np.abs(intensity) / self.noise_level

                if snr >= MINIMUM_SNR:
                    detected.append({
                        'y_point': y,
                        'x_point': x,
                        'ppm_x': ppm_x,
                        'ppm_y': ppm_y,
                        'intensity': intensity,
                        'snr': snr,
                        'detected': True,
                        'detection_quality': 'Traditional'
                    })

        # Sort by intensity
        detected.sort(key=lambda x: abs(x['intensity']), reverse=True)

        # Apply Peak Centroid Refinement (always enabled for enhanced accuracy)
        gui_params = getattr(self, 'gui_params', {})
        detected = self._refine_peaks_with_centroids(detected, gui_params)

        # Store refined peaks (not original unrefined peaks)
        self.detected_peaks = detected

        print(f"Detected {len(detected)} peaks above threshold")
        return detected

    def match_peaks_to_list(self):
        """Match detected peaks to reference list (legacy compatibility)"""
        if self.use_reference_detection:
            print("Using reference-based detection - matching not needed")
            return self.fitted_peaks

        # Traditional matching for non-reference mode
        return self._match_traditional()

    def _match_traditional(self):
        """Traditional peak matching (legacy method)"""
        # Handle DataFrame properly to avoid boolean evaluation errors
        detected_peaks_available = (
            self.detected_peaks is not None and
            (not hasattr(self.detected_peaks, 'empty') or not self.detected_peaks.empty) and
            len(self.detected_peaks) > 0
        )

        if not detected_peaks_available or self.peak_list is None:
            print("No detected peaks or peak list available for matching")
            return []

        matched_peaks = []

        for _, peak_row in self.peak_list.iterrows():
            ref_x = float(peak_row['Position_X'])
            ref_y = float(peak_row['Position_Y'])
            assignment = peak_row.get('Assignment', 'Unknown')

            # Find closest detected peak
            best_match = None
            min_distance = float('inf')

            for detected in self.detected_peaks:
                distance = np.sqrt((detected['ppm_x'] - ref_x)**2 +
                                 (detected['ppm_y'] - ref_y)**2)

                if distance < min_distance and distance < SEARCH_RADIUS_PPM:
                    min_distance = distance
                    best_match = detected

            if best_match:
                matched_peaks.append({
                    'assignment': assignment,
                    'ppm_x': best_match['ppm_x'],
                    'ppm_y': best_match['ppm_y'],
                    'x_point': best_match['x_point'],
                    'y_point': best_match['y_point'],
                    'intensity': best_match['intensity'],
                    'snr': best_match['snr'],
                    'detected': True,
                    'detection_quality': 'Matched'
                })
            else:
                # Not found
                matched_peaks.append({
                    'assignment': assignment,
                    'ppm_x': ref_x,
                    'ppm_y': ref_y,
                    'x_point': -1,
                    'y_point': -1,
                    'intensity': 0,
                    'snr': 0,
                    'detected': False,
                    'detection_quality': 'Not Found'
                })

        self.fitted_peaks = matched_peaks
        return matched_peaks

    def integrate_peaks(self):
        """Integrate detected/matched peaks using enhanced box integration"""
        # Handle DataFrame properly to avoid boolean evaluation errors
        fitted_peaks_available = (
            self.fitted_peaks is not None and
            (not hasattr(self.fitted_peaks, 'empty') or not self.fitted_peaks.empty) and
            len(self.fitted_peaks) > 0
        )

        if not fitted_peaks_available:
            print("No fitted peaks available for integration")
            return []

        print("Integrating peaks using enhanced methods...")

        results = []
        box_size = 3  # Integration box size

        for peak in self.fitted_peaks:
            if not peak.get('detected', False):
                # Skip undetected peaks or use reference position
                result = {
                    'Assignment': peak['assignment'],
                    'Position_X': peak['ppm_x'],
                    'Position_Y': peak['ppm_y'],
                    'Height': 0,
                    'Volume': 0,
                    'SNR': peak.get('snr', 0),
                    'Quality': 'Not Detected',
                    'Integration_Method': 'Reference'
                }
            else:
                # Perform box integration around detected peak
                x_point = peak['x_point']
                y_point = peak['y_point']

                # Define integration box
                x_min = max(0, x_point - box_size)
                x_max = min(self.nmr_data.shape[1], x_point + box_size + 1)
                y_min = max(0, y_point - box_size)
                y_max = min(self.nmr_data.shape[0], y_point + box_size + 1)

                # Extract integration region
                region = self.nmr_data[y_min:y_max, x_min:x_max]

                # Calculate volume (sum) and height (peak value)
                volume = np.sum(region)
                height = peak['intensity']
                snr = peak.get('snr', abs(height) / self.noise_level)

                # Determine quality
                if snr > 10:
                    quality = 'Excellent'
                elif snr > 5:
                    quality = 'Good'
                elif snr > 3:
                    quality = 'Fair'
                else:
                    quality = 'Poor'

                result = {
                    'Assignment': peak['assignment'],
                    'Position_X': peak['ppm_x'],
                    'Position_Y': peak['ppm_y'],
                    'Height': abs(height),
                    'Volume': abs(volume),
                    'SNR': snr,
                    'Quality': quality,
                    'Integration_Method': 'Enhanced_Box',
                    'Detection_Quality': peak.get('detection_quality', 'Unknown'),
                    'Distance_From_Reference': peak.get('distance_from_reference', 0.0)
                }

            results.append(result)

        self.integration_results = results

        # Print summary
        good_quality = sum(1 for r in results if r['Quality'] in ['Excellent', 'Good'])
        detected_count = sum(1 for r in results if r['Integration_Method'] != 'Reference')

        print(f"Integration completed:")
        print(f"  Total peaks: {len(results)}")
        print(f"  Detected and integrated: {detected_count}")
        print(f"  Good quality (≥SNR 5): {good_quality}")

        return results

    def get_detection_statistics(self):
        """Get detailed detection statistics"""
        return self.detection_statistics.copy()

    def export_results(self, output_file):
        """Export integration results to file"""
        # Handle DataFrame properly to avoid boolean evaluation errors
        integration_results_available = (
            self.integration_results is not None and
            (not hasattr(self.integration_results, 'empty') or not self.integration_results.empty) and
            len(self.integration_results) > 0
        )

        if not integration_results_available:
            print("No integration results to export")
            return False

        try:
            df = pd.DataFrame(self.integration_results)
            df.to_csv(output_file, index=False, float_format='%.6f')
            print(f"Results exported to {output_file}")
            return True
        except Exception as e:
            print(f"Export failed: {e}")
            return False

    def create_summary_plot(self, output_file=None):
        """Create summary visualization of detection and integration results"""
        if self.nmr_data is None:
            print("No NMR data available for plotting")
            return

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIGURE_SIZE)

        # Plot 1: Spectrum with detected peaks
        X, Y = np.meshgrid(self.ppm_x_axis, self.ppm_y_axis)
        levels = np.linspace(self.threshold, np.max(np.abs(self.nmr_data)), CONTOUR_LEVELS)

        ax1.contour(X, Y, np.abs(self.nmr_data), levels=levels, cmap=COLORMAP, linewidths=0.5)

        # Plot detected peaks
        for peak in self.fitted_peaks:
            if peak.get('detected', False):
                color = 'green'
                marker = 'o'
            elif peak.get('reference_retained', False):
                color = 'orange'
                marker = 's'
            else:
                color = 'red'
                marker = 'x'

            ax1.scatter(peak['ppm_x'], peak['ppm_y'], c=color, marker=marker, s=50, alpha=0.8)

            if SHOW_PEAK_LABELS:
                ax1.annotate(peak['assignment'], (peak['ppm_x'], peak['ppm_y']),
                           xytext=(3, 3), textcoords='offset points',
                           fontsize=LABEL_FONT_SIZE, alpha=0.8)

        ax1.set_xlabel('¹H Chemical Shift (ppm)')
        ax1.set_ylabel('¹⁵N/¹³C Chemical Shift (ppm)')
        ax1.set_title('Enhanced Peak Detection Results')
        ax1.invert_xaxis()
        ax1.invert_yaxis()
        ax1.grid(True, alpha=0.3)

        # Plot 2: Detection statistics
        stats = self.detection_statistics
        labels = ['Detected', 'Reference\nRetained', 'Failed']
        sizes = [stats['detected_peaks'], stats['reference_retained'],
                stats['total_peaks'] - stats['detected_peaks'] - stats['reference_retained']]
        colors = ['green', 'orange', 'red']

        ax2.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
        ax2.set_title(f'Detection Statistics\n(Total: {stats["total_peaks"]} peaks)')

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
            print(f"Summary plot saved to {output_file}")

        if SHOW_PLOTS:
            plt.show()
        else:
            plt.close()

def main():
    """Main execution function"""
    print("Enhanced NMR Peak Integration - Reference-Based Detection")
    print("=" * 60)

    # Initialize integrator
    integrator = InPlaceAdvancedNMRIntegrator()

    # Load data
    success = integrator.load_data(PEAK_LIST_FILE, NMR_DATA_PATH)
    if not success:
        print("Failed to load data. Exiting.")
        return

    # Set enhanced parameters (examples)
    integrator.set_search_window(SEARCH_WINDOW_X_PPM, SEARCH_WINDOW_Y_PPM)
    integrator.set_threshold_multiplier(DEFAULT_THRESHOLD_MULTIPLIER)

    # Enable custom noise regions if specified
    if USE_CUSTOM_NOISE_REGIONS and NOISE_REGIONS_PPM:
        integrator.set_noise_regions(NOISE_REGIONS_PPM)

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Perform detection
    detected = integrator.detect_peaks_professional()
    print(f"\nDetection completed: {len(detected)} peaks processed")

    # Perform integration
    results = integrator.integrate_peaks()
    print(f"Integration completed: {len(results)} peaks integrated")

    # Export results
    output_file = os.path.join(OUTPUT_DIR, "enhanced_integration_results.csv")
    integrator.export_results(output_file)

    # Create visualization
    if CREATE_VISUALIZATION:
        plot_file = os.path.join(OUTPUT_DIR, "enhanced_detection_summary.png")
        integrator.create_summary_plot(plot_file)

    # Print final statistics
    stats = integrator.get_detection_statistics()
    print(f"\nFinal Statistics:")
    print(f"  Detection rate: {stats['detection_rate']:.1f}%")
    print(f"  Detected peaks: {stats['detected_peaks']}")
    print(f"  Reference retained: {stats['reference_retained']}")

    print(f"\nResults saved to {OUTPUT_DIR}/")
    print("Enhanced peak integration completed successfully!")

if __name__ == "__main__":
    main()
