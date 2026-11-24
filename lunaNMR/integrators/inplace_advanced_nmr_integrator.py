#!/usr/bin/env python3
# ABOUTME: Base class of the live integrator: loads a 2D spectrum (NMRPipe or Bruker pdata) and
# ABOUTME: its peak list, builds the ppm axes, detects the nucleus and estimates the noise level.
"""
Spectrum loading for the NMR integrators

The live integrator (core_integrator.EnhancedVoigtIntegrator) inherits its loading,
ppm axes, nucleus detection and noise estimate from this class, and its default
attributes.
"""

# DEFAULTS OF THE INTEGRATOR'S ATTRIBUTES

# Threshold control parameters
DEFAULT_THRESHOLD_MULTIPLIER = 3.0                # Default threshold multiplier

import numpy as np
import pandas as pd

import os
import sys

from lunaNMR.utils.output_manager import log_info, log_warning, log_error

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
        self.fitted_peaks = []

        # Enhanced parameters
        self.threshold_multiplier = DEFAULT_THRESHOLD_MULTIPLIER

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
        try:
            self.peak_list = pd.read_csv(peak_list_file, skipinitialspace=True)
            self.peak_list.columns = self.peak_list.columns.str.strip()
            log_info(f"Loaded peak list: {len(self.peak_list)} peaks from {peak_list_file}")
        except Exception as e:
            log_error(f"Error loading peak list: {e}")
            return False

        try:
            file_format = self._detect_nmr_format(nmr_file)

            if file_format == 'bruker_pdata':
                # Bruker processed data - use parent directory
                pdata_dir = os.path.dirname(nmr_file)
                self.nmr_dict, self.nmr_data = ng.bruker.read_pdata(dir=pdata_dir)
                self._nmr_format = 'bruker'
                log_info(f"Loaded Bruker pdata: {self.nmr_data.shape} from {pdata_dir}")
            else:
                # NMRPipe format (default)
                self.nmr_dict, self.nmr_data = ng.pipe.read(nmr_file)
                self._nmr_format = 'pipe'
                log_info(f"Loaded NMR data: {self.nmr_data.shape} from {nmr_file}")

            self._calculate_ppm_axes()

            self._estimate_noise_level()

            return True
        except Exception as e:
            log_error(f"Error loading NMR data: {e}")
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

            f2_ppm_scale = uc_f2.ppm_scale()
            f1_ppm_scale = uc_f1.ppm_scale()

            # Create PPM axes (nmrglue handles all the complex parameter conversion)
            self.ppm_x_axis = f2_ppm_scale
            self.ppm_y_axis = f1_ppm_scale
            self._using_fallback_ppm_axes = False

            log_info(f"F2 (1H) range: {f2_ppm_scale[0]:.2f} to {f2_ppm_scale[-1]:.2f} ppm")
            log_info(f"F1 (15N/13C) range: {f1_ppm_scale[0]:.1f} to {f1_ppm_scale[-1]:.1f} ppm")

        except Exception as e:
            log_warning(f"Could not calculate PPM axes properly: {e}")
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
            log_warning("Skipping nucleus auto-detection (using fallback PPM axes)")
            return None

        try:
            x_min, x_max = float(min(self.ppm_x_axis)), float(max(self.ppm_x_axis))
            y_min, y_max = float(min(self.ppm_y_axis)), float(max(self.ppm_y_axis))

            x_center = (x_min + x_max) / 2.0
            y_center = (y_min + y_max) / 2.0

            from lunaNMR.core.ps2d_config import nucleus_of_spectral_window
            return nucleus_of_spectral_window(x_center, y_center)

        except Exception as e:
            log_warning(f"Error during nucleus auto-detection: {e}")
            return None

    def _estimate_noise_level(self):
        """Set the spectrum's noise level and the threshold that follows from it"""
        self.noise_level = self._estimate_noise_level_advanced()

        self.threshold = self.noise_level * self.threshold_multiplier
        log_info(f"Estimated noise level: {self.noise_level:.2f}")
        log_info(f"Detection threshold: {self.threshold:.2f} (multiplier: {self.threshold_multiplier})")
