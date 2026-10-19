#!/usr/bin/env python3
# ABOUTME: Base class of the live integrator: loads a 2D spectrum (NMRPipe or Bruker pdata) and
# ABOUTME: its peak list, builds the ppm axes, detects the nucleus and estimates the noise level.
"""
Spectrum loading for the NMR integrators

The live integrator (core_integrator.EnhancedVoigtIntegrator) inherits its loading,
ppm axes, nucleus detection and noise estimate from this class, and its default
attributes.
"""

# ========================================================================
# DEFAULTS OF THE INTEGRATOR'S ATTRIBUTES
# ========================================================================

# Reference-based detection parameters
USE_REFERENCE_DETECTION = True                      # Enable reference-based peak detection
SEARCH_WINDOW_X_PPM = 0.2                         # ±X ppm search window (1H dimension)
SEARCH_WINDOW_Y_PPM = 3.0                         # ±Y ppm search window (15N/13C dimension)

# Noise region configuration
USE_CUSTOM_NOISE_REGIONS = False                   # Use user-defined noise regions
NOISE_REGIONS_PPM = [                             # Noise regions in ppm [x_min, x_max, y_min, y_max]
    [0.5, 1.5, 100, 110],                        # Example region 1
    [11, 12, 130, 140]                           # Example region 2
]

# Threshold control parameters
DEFAULT_THRESHOLD_MULTIPLIER = 3.0                # Default threshold multiplier

import numpy as np
import pandas as pd

import os
import sys


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

    def _detect_nmr_format(self, file_path):
        """Detect NMR file format from file path.

        Returns:
            'bruker_pdata': Bruker processed data (2rr, 2ri, 2ir, 2ii, 1r, 1i)
            'pipe': NMRPipe format (.ft, .ft2, .pipe, etc.)
        """
        basename = os.path.basename(file_path)

        # Bruker processed data files (no extension, specific names)
        bruker_pdata_names = (
            '2rr', '2ri', '2ir', '2ii', '1r', '1i',
            '3rrr', '3rri', '3rir', '3rii', '3irr', '3iri', '3iir', '3iii'
        )
        if basename in bruker_pdata_names:
            return 'bruker_pdata'

        return 'pipe'

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
            # Detect format and load accordingly
            file_format = self._detect_nmr_format(nmr_file)

            if file_format == 'bruker_pdata':
                # Bruker processed data - use parent directory
                pdata_dir = os.path.dirname(nmr_file)
                self.nmr_dict, self.nmr_data = ng.bruker.read_pdata(dir=pdata_dir)
                self._nmr_format = 'bruker'
                print(f"Loaded Bruker pdata: {self.nmr_data.shape} from {pdata_dir}")
            else:
                # NMRPipe format (default)
                self.nmr_dict, self.nmr_data = ng.pipe.read(nmr_file)
                self._nmr_format = 'pipe'
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
            nmr_format = getattr(self, '_nmr_format', 'pipe')

            if nmr_format == 'bruker':
                # Bruker: use universal dictionary approach
                # strip_fake=True is required for processed data to get correct SW from procs
                udic = ng.bruker.guess_udic(self.nmr_dict, self.nmr_data, strip_fake=True)
                uc_f2 = ng.fileiobase.uc_from_udic(udic, dim=1)  # F2 / 1H
                uc_f1 = ng.fileiobase.uc_from_udic(udic, dim=0)  # F1 / 15N/13C
            else:
                # NMRPipe: use pipe-specific conversion
                uc_f2 = ng.pipe.make_uc(self.nmr_dict, self.nmr_data, dim=1)  # F2 dimension (1H)
                uc_f1 = ng.pipe.make_uc(self.nmr_dict, self.nmr_data, dim=0)  # F1 dimension (15N/13C)

            # Get the full PPM scales
            f2_ppm_scale = uc_f2.ppm_scale()
            f1_ppm_scale = uc_f1.ppm_scale()

            # Create PPM axes (nmrglue handles all the complex parameter conversion)
            self.ppm_x_axis = f2_ppm_scale
            self.ppm_y_axis = f1_ppm_scale
            self._using_fallback_ppm_axes = False

            print(f"F2 (1H) range: {f2_ppm_scale[0]:.2f} to {f2_ppm_scale[-1]:.2f} ppm")
            print(f"F1 (15N/13C) range: {f1_ppm_scale[0]:.1f} to {f1_ppm_scale[-1]:.1f} ppm")

        except Exception as e:
            print(f"Warning: Could not calculate PPM axes properly: {e}")
            # Fallback to simple linear axes
            self.ppm_x_axis = np.linspace(12, 0, self.nmr_data.shape[1])
            self.ppm_y_axis = np.linspace(140, 100, self.nmr_data.shape[0])
            self._using_fallback_ppm_axes = True

    def _detect_nucleus_type(self):
        """
        Automatically detect nucleus type based on spectral dimensions.

        Detection criteria:
          15N-HSQC: X: 5-12 ppm, Y: 100-140 ppm (amide region)
          13C-HSQC: X: -2-4 ppm, Y: 0-80 ppm (aliphatic region, aromatic excluded)

        Returns:
            str: '15N', '13C', or None (if cannot determine)
        """
        if self.ppm_x_axis is None or self.ppm_y_axis is None:
            return None

        # Skip auto-detection if using fallback axes (unreliable)
        if hasattr(self, '_using_fallback_ppm_axes') and self._using_fallback_ppm_axes:
            print("⚠️  Skipping nucleus auto-detection (using fallback PPM axes)")
            return None

        try:
            # Get spectral ranges
            x_min, x_max = float(min(self.ppm_x_axis)), float(max(self.ppm_x_axis))
            y_min, y_max = float(min(self.ppm_y_axis)), float(max(self.ppm_y_axis))

            # Calculate center points for weighted detection
            x_center = (x_min + x_max) / 2.0
            y_center = (y_min + y_max) / 2.0

            # Priority 1: Classic 15N-HSQC (amide region)
            # X: 5-12 ppm (amide protons), Y: 100-140 ppm (backbone/sidechain nitrogens)
            if (5 <= x_center <= 12) and (100 <= y_center <= 140):
                return '15N'

            # Priority 2: Classic 13C-HSQC (aliphatic region)
            # X: -2 to 5 ppm (aliphatic CH), Y: 0-80 ppm (aliphatic carbons)
            if (-2 <= x_center <= 5) and (0 <= y_center <= 80):
                return '13C'

            # Priority 3: Edge cases based on Y-dimension
            if y_center < 90:  # 13C aliphatic likely
                return '13C'
            if y_center > 95:  # 15N backbone likely
                return '15N'

            # Cannot determine with confidence
            return None

        except Exception as e:
            print(f"⚠️  Error during nucleus auto-detection: {e}")
            return None

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
