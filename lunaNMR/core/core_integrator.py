#!/usr/bin/env python3
# ABOUTME: Main integration engine coordinating peak detection, centroid refinement, overlap clustering, and PS2D 2D fitting.
# ABOUTME: Entry point for all spectrum processing workflows (single, batch, series, GUI, CLI).
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
import sys
from typing import List, Dict, Tuple, Any
from scipy.special import wofz  # Faddeeva function for Voigt profiles
from scipy.optimize import curve_fit
import warnings
from lunaNMR.core.ps2d_data_selector import select_data_2d_for_overlap_group
from lunaNMR.core.ps2d_config import get_ps2d_config, set_ps2d_config
warnings.filterwarnings('ignore')

# Import enhanced Voigt fitter
try:
    from lunaNMR.core.enhanced_voigt_fitter import EnhancedVoigtFitter
    ENHANCED_FITTING_AVAILABLE = True
except ImportError:
    ENHANCED_FITTING_AVAILABLE = False
    print("⚠️ Enhanced Voigt fitting not available, using standard fitting")


# Import PS2D-style high-performance fitter (CRITICAL FOR HIGH-QUALITY FITTING)
try:
    from lunaNMR.core.ps2d_style_fitter import (
        MultiStageFitter,
        fit_single_peak_ps2d_style,
        voigt_profile_1d as ps2d_voigt_1d
    )
    PS2D_STYLE_FITTING_AVAILABLE = True
except ImportError as e:
    PS2D_STYLE_FITTING_AVAILABLE = False

# Import PS2D-style multi-peak fitter (EXACT PORT FROM PS2D_SRC)
try:
    from lunaNMR.core.ps2d_multi_peak_fitter import (
        Ps2dMultiPeakFitter,
        fit_overlapping_peaks_ps2d_style
    )
    PS2D_MULTI_PEAK_AVAILABLE = True
except ImportError as e:
    PS2D_MULTI_PEAK_AVAILABLE = False

# Import PS2D data selector (EXACT CLONE OF SPECTRUM.CPP)
try:
    from lunaNMR.core.ps2d_data_selector import Ps2dDataSelector
    PS2D_DATA_SELECTOR_AVAILABLE = True
except ImportError as e:
    PS2D_DATA_SELECTOR_AVAILABLE = False

# Import PS2D exact overlap integration wrapper (EXACT C++ overlap detection)
try:
    from lunaNMR.core.ps2d_exact_overlap_integration import (
        fit_peaks_with_exact_overlap_detection,
        fit_peaks_with_ps2d_exact_overlap,
        OVERLAP_DETECTOR_AVAILABLE,
        MULTI_PEAK_FITTER_AVAILABLE as OVERLAP_FITTER_AVAILABLE
    )
    PS2D_EXACT_OVERLAP_INTEGRATION_AVAILABLE = True
except ImportError as e:
    PS2D_EXACT_OVERLAP_INTEGRATION_AVAILABLE = False

# Import PS2D 2D multi-peak fitter (for closely-spaced overlapping peaks in 2D)
try:
    from lunaNMR.core.ps2d_2d_fitter import Ps2dMultiPeakFitter2D
    from lunaNMR.core.ps2d_data_selector import select_data_2d_for_overlap_group
    PS2D_2D_FITTER_AVAILABLE = True
except ImportError as e:
    PS2D_2D_FITTER_AVAILABLE = False
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
            self.fitted_peaks = []

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
                self._using_fallback_ppm_axes = False

                print(f"PPM axes calculated - X: {self.ppm_x_axis[0]:.2f} to {self.ppm_x_axis[-1]:.2f} ppm")
                print(f"                      Y: {self.ppm_y_axis[0]:.1f} to {self.ppm_y_axis[-1]:.1f} ppm")
            except Exception as e:
                print(f"Warning: Could not calculate PPM axes properly: {e}")
                # Fallback to simple linear axes - ensures plotting/zoom always works
                print("Creating fallback PPM axes...")
                import numpy as np
                self.ppm_x_axis = np.linspace(12, 0, self.nmr_data.shape[1])
                self.ppm_y_axis = np.linspace(140, 100, self.nmr_data.shape[0])
                self._using_fallback_ppm_axes = True
                print(f"Fallback PPM axes created - X: {self.ppm_x_axis[0]:.2f} to {self.ppm_x_axis[-1]:.2f} ppm")
                print(f"                               Y: {self.ppm_y_axis[0]:.1f} to {self.ppm_y_axis[-1]:.1f} ppm")

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
            'min_r_squared': 0.1,  # Lowered to 0.1 for debugging/QC
            'fitting_window_x': 0.2,  # ppm @GM was 0.4 was 0.4 was 0.15
            'fitting_window_y': 2,   # ppm  @GM was 7 was 5 was 3 was 1.5
            'max_expansion_attempts': 3,  # Maximum attempts to expand window on failure
            'failure_expansion_factor': 1.5  # Window expansion factor on fitting failure
        }
        self.processing_mode = 'full_detection'  # 'full_detection' or 'in_place' or 'sn_native'

        # S/N threshold parameters
        self.sn_threshold = 3.0
        self.expected_peak_count = 50
        self.sn_detection_params = {
            'min_snr': 2.0,          # Minimum signal-to-noise ratio
            'max_peaks': 500,        # Maximum peaks to detect
            'noise_estimation': 'corners',  # 'corners' or 'median'
            'peak_separation': 0.01  # Minimum peak separation in ppm
        }

        # Initialize enhanced fitter if available
        if ENHANCED_FITTING_AVAILABLE:
            self.enhanced_fitter = EnhancedVoigtFitter()
            self.enhanced_fitter.parent = self  # Set parent reference for parallel processing
        else:
            self.enhanced_fitter = None

        # NEW: Overlap resolution configuration (optional)
        self.overlap_resolution_enabled = False  # Default OFF for backward compatibility

        # Initialize with nucleus-specific windows (will be overridden by GUI if provided)
        self.detected_nucleus_type = None
        self.gui_window_override = {'x': None, 'y': None}  # Track GUI overrides

        # Nucleus-specific default fitting windows
        # In HSQC: X-axis=1H (direct), Y-axis=15N (indirect)
        self.nucleus_default_windows = {
            'default': {'x': 0.06, 'y': 0.6},
            '15N-HSQC': {'x': 0.06, 'y': 0.6},   # 15N-HSQC: X=1H, Y=15N
            '13C-HSQC': {'x': 0.06, 'y': 0.6},   # 13C-HSQC: X=1H, Y=13C (same window as 15N)
            '1H-15N': {'x': 0.06, 'y': 0.6},     # Alternative naming
            '1H-13C': {'x': 0.06, 'y': 0.6}      # Alternative naming
        }

        # Detection search windows (set by GUI)
        self.search_window_x = 0.01  # Default X search window in ppm
        self.search_window_y = 0.05  # Default Y search window in ppm

        # Integration parameters (for integrated detection-fitting mode)
        self.integration_parameters = {
            'enable_integrated_mode': False  # Disabled by default
        }

    def set_search_window(self, x_ppm, y_ppm):
        """Set search window parameters for peak detection

        Args:
            x_ppm (float): ±X ppm search window (1H dimension)
            y_ppm (float): ±Y ppm search window (15N/13C dimension)
        """
        self.search_window_x = x_ppm
        self.search_window_y = y_ppm
        #print(f"🔍 Search windows set: X=±{x_ppm:.3f} ppm, Y=±{y_ppm:.3f} ppm")

    def configure_overlap_resolution(self, enable: bool = True, config: Dict = None):
        """
        Configure overlap resolution for the integrator

        This enables automatic detection and resolution of overlapping peaks
        during fitting operations.

        Args:
            enable: Enable/disable overlap resolution (default: True)
            config: Optional configuration dictionary with keys:
                - 'jackknife': {'n_resamples': int, 'cv_threshold': float}
                - 'staged_fitting': {'stage2_position_tolerance': float, ...}
                - 'model_selection': {'max_peaks': int, 'criterion': 'AIC'|'BIC'}
                - 'correlation': {'amplitude_correlation_threshold': float}

        Example:
            integrator = VoigtIntegrator()
            integrator.configure_overlap_resolution(
                enable=True,
                config={
                    'jackknife': {'n_resamples': 100},
                    'model_selection': {'max_peaks': 5}
                }
            )

        Note:
            Overlap resolution is disabled by default for backward compatibility.
            Enable it explicitly when you expect overlapping peaks.
        """
        self.overlap_resolution_enabled = enable

        if self.enhanced_fitter is not None:
            # CRITICAL: Call enhanced_fitter's configure_overlap_resolution method
            # This sets the overlap_detection_enabled flag that triggers routing to OverlapResolverEngine
            self.enhanced_fitter.configure_overlap_resolution(enable=enable, config=config)

            status = "enabled" if enable else "disabled"
            #print(f"   ✅ Overlap resolution {status} at integrator level")
        else:
            print(f"   ⚠️  Enhanced fitter not available - cannot configure overlap resolution")

    def enable_overlap_resolution_preset(self, preset: str = 'default'):
        """
        Enable overlap resolution with preset configurations

        Args:
            preset: Configuration preset
                - 'default': Balanced settings (50 resamples)
                - 'fast': Speed-optimized (25 resamples)
                - 'thorough': Accuracy-optimized (100 resamples)

        Example:
            integrator.enable_overlap_resolution_preset('thorough')
        """
        preset_configs = {
            'default': {
                'jackknife': {'n_resamples': 50},
                'model_selection': {'max_peaks': 10}
            },
            'fast': {
                'jackknife': {'n_resamples': 25},
                'staged_fitting': {
                    'stage1_max_iterations': 500,
                    'stage2_max_iterations': 750,
                    'stage3_max_iterations': 1000
                }
            },
            'thorough': {
                'jackknife': {'n_resamples': 100},
                'staged_fitting': {
                    'stage1_max_iterations': 2000,
                    'stage2_max_iterations': 2500,
                    'stage3_max_iterations': 3000
                },
                'model_selection': {'max_peaks': 15}
            }
        }

        if preset not in preset_configs:
            print(f"⚠️  Unknown preset '{preset}'. Using 'default'.")
            preset = 'default'

        self.configure_overlap_resolution(enable=True, config=preset_configs[preset])
        print(f"   📊 Overlap resolution enabled with '{preset}' preset")

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
        #print(f"   X-axis: {ppm_x_axis[0]:.2f} to {ppm_x_axis[-1]:.2f} ppm")
        #print(f"   Y-axis: {ppm_y_axis[0]:.1f} to {ppm_y_axis[-1]:.1f} ppm")

        # Set dummy nmr_dict for compatibility
        self.nmr_dict = {'dummy': 'for_testing'}

        return True

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

    def update_fitting_parameters_from_gui(self, gui_params):
        """
        Update fitting parameters from GUI with proper priority handling

        Priority: GUI > Nucleus-specific > Global defaults
        """
        if not gui_params:
            return

        # Store GUI overrides
        if 'fitting_window_x' in gui_params:
            self.gui_window_override['x'] = gui_params['fitting_window_x']
            self.fitting_parameters['fitting_window_x'] = gui_params['fitting_window_x']
            print(f"   🎛️  GUI window X override: {gui_params['fitting_window_x']:.3f} ppm")

        if 'fitting_window_y' in gui_params:
            self.gui_window_override['y'] = gui_params['fitting_window_y']
            self.fitting_parameters['fitting_window_y'] = gui_params['fitting_window_y']
            print(f"   🎛️  GUI window Y override: {gui_params['fitting_window_y']:.1f} ppm")

        # Update other parameters normally
        for key, value in gui_params.items():
            if key not in ['fitting_window_x', 'fitting_window_y']:
                self.fitting_parameters[key] = value

    def get_adaptive_window_parameters(self, nucleus_type=None):
        """
        Get window parameters with proper priority: GUI > Nucleus-specific > Defaults
        """
        if nucleus_type is None:
            nucleus_type = self.detected_nucleus_type or 'default'

        # Start with nucleus-specific defaults
        nucleus_windows = self.nucleus_default_windows.get(nucleus_type,
                                                          self.nucleus_default_windows['default'])

        # Apply GUI overrides if they exist
        window_x = self.gui_window_override['x'] if self.gui_window_override['x'] is not None else nucleus_windows['x']
        window_y = self.gui_window_override['y'] if self.gui_window_override['y'] is not None else nucleus_windows['y']

        # Update fitting parameters
        self.fitting_parameters['fitting_window_x'] = window_x
        self.fitting_parameters['fitting_window_y'] = window_y

        return window_x, window_y

    def fit_peak_with_adaptive_windows(self, peak_x_ppm, peak_y_ppm, assignment="Unknown",
                                      force_single_peak=False):
        """
        Fit peak with adaptive window sizing and OPTIONAL overlap resolution

        Priority: GUI > Nucleus-specific > Failure-based expansion

        NEW PARAMETER:
            force_single_peak: Force single-peak fitting, skip overlap detection
                              (default: False, allows auto-detection if enabled)

        The overlap resolution behavior is controlled by:
        1. Global enable: self.overlap_resolution_enabled
        2. Per-peak override: force_single_peak parameter

        If overlap resolution is enabled AND not force_single_peak,
        the fitter will automatically detect and resolve overlapping peaks.
        """
        # Detect nucleus type if not already done
        if self.detected_nucleus_type is None:
            self.detected_nucleus_type = self.detect_nucleus_type([self.ppm_y_axis[0], self.ppm_y_axis[-1]])

        # Get initial window parameters
        window_x, window_y = self.get_adaptive_window_parameters(self.detected_nucleus_type)

        # Track window source for diagnostics
        window_source = "GUI" if (self.gui_window_override['x'] is not None or
                                 self.gui_window_override['y'] is not None) else f"{self.detected_nucleus_type}-specific"

        print(f"   🎯 Attempting peak fit: {assignment}")
        print(f"   📊 Initial windows: X={window_x:.3f} ppm, Y={window_y:.1f} ppm ({window_source})")

        # Show overlap resolution status
        if self.overlap_resolution_enabled and not force_single_peak:
            print(f"   🔬 Overlap resolution: ENABLED (will auto-detect overlapping peaks)")
        elif force_single_peak:
            print(f"   🔒 Overlap resolution: DISABLED (forced single-peak mode)")
        else:
            print(f"   📌 Overlap resolution: DISABLED (default mode)")

        # Attempt fitting with progressive window expansion
        for attempt in range(self.fitting_parameters['max_expansion_attempts']):
            try:
                # Extract peak region with current windows
                regions = self.extract_peak_region(peak_x_ppm, peak_y_ppm, window_x, window_y)

                # Attempt fitting using enhanced fitter if available
                if self.enhanced_fitter is not None:
                    # Determine fitting method based on overlap resolution settings
                    use_overlap_detection = self.overlap_resolution_enabled and not force_single_peak

                    if use_overlap_detection:
                        # X-dimension fit with overlap resolution
                        x_result = self.enhanced_fitter.fit_peak_enhanced(
                            regions['x_ppm_scale'], regions['x_cross_section'],
                            peak_x_ppm, nucleus_type=self.detected_nucleus_type
                        )

                        # Y-dimension fit with overlap resolution
                        y_result = self.enhanced_fitter.fit_peak_enhanced(
                            regions['y_ppm_scale'], regions['y_cross_section'],
                            peak_y_ppm, nucleus_type=self.detected_nucleus_type
                        )
                    else:
                        # X-dimension fit (single-peak mode)
                        x_result = self.enhanced_fitter.fit_single_peak_1d(
                            regions['x_ppm_scale'], regions['x_cross_section'],
                            peak_x_ppm, nucleus_type=self.detected_nucleus_type
                        )

                        # Y-dimension fit (single-peak mode)
                        y_result = self.enhanced_fitter.fit_single_peak_1d(
                            regions['y_ppm_scale'], regions['y_cross_section'],
                            peak_y_ppm, nucleus_type=self.detected_nucleus_type
                        )

                    # Evaluate fit quality
                    x_r_squared = x_result.get('r_squared', 0)
                    y_r_squared = y_result.get('r_squared', 0)
                    avg_r_squared = (x_r_squared + y_r_squared) / 2

                    print(f"   📈 Attempt {attempt + 1}: R²_X={x_r_squared:.3f}, R²_Y={y_r_squared:.3f}, Avg={avg_r_squared:.3f}")

                    # Check if fit is acceptable
                    if avg_r_squared >= self.fitting_parameters['min_r_squared']:
                        print(f"   ✅ Successful fit with windows: X={window_x:.3f}, Y={window_y:.1f} ppm")
                        return {
                            'success': True,
                            'x_result': x_result,
                            'y_result': y_result,
                            'r_squared': avg_r_squared,
                            'window_x_used': window_x,
                            'window_y_used': window_y,
                            'expansion_attempts': attempt + 1,
                            'window_source': window_source
                        }

                    # If fit quality is poor and we have expansion attempts left
                    elif attempt < self.fitting_parameters['max_expansion_attempts'] - 1:
                        # Only expand if not GUI-constrained
                        if self.gui_window_override['x'] is None:
                            window_x *= self.fitting_parameters['failure_expansion_factor']
                        if self.gui_window_override['y'] is None:
                            window_y *= self.fitting_parameters['failure_expansion_factor']

                        print(f"   ⚠️  Poor fit (R²={avg_r_squared:.3f}), expanding windows to X={window_x:.3f}, Y={window_y:.1f} ppm")
                        continue

                    else:
                        # Return result even if below threshold for debugging/QC
                        print(f"   ⚠️ Poor fit after {self.fitting_parameters['max_expansion_attempts']} attempts (R²={avg_r_squared:.3f})")
                        return {
                            'success': True,  # Pass through for debugging/QC
                            'x_result': x_result,
                            'y_result': y_result,
                            'r_squared': avg_r_squared,
                            'window_x_used': window_x,
                            'window_y_used': window_y,
                            'expansion_attempts': attempt + 1,
                            'window_source': window_source,
                            'warning': 'poor_fit_quality'
                        }
                else:
                    # Fallback to basic fitting without enhanced fitter
                    print(f"   ⚠️  Enhanced fitter not available, using basic fitting")
                    return {
                        'success': False,
                        'failure_reason': 'no_enhanced_fitter',
                        'window_x_used': window_x,
                        'window_y_used': window_y,
                        'expansion_attempts': 1,
                        'window_source': window_source
                    }

            except Exception as e:
                print(f"   ❌ Fitting error on attempt {attempt + 1}: {str(e)}")
                if attempt < self.fitting_parameters['max_expansion_attempts'] - 1:
                    # Expand windows and try again
                    if self.gui_window_override['x'] is None:
                        window_x *= self.fitting_parameters['failure_expansion_factor']
                    if self.gui_window_override['y'] is None:
                        window_y *= self.fitting_parameters['failure_expansion_factor']
                    continue
                else:
                    return {
                        'success': False,
                        'failure_reason': f'fitting_exception: {str(e)}',
                        'window_x_used': window_x,
                        'window_y_used': window_y,
                        'expansion_attempts': attempt + 1,
                        'window_source': window_source
                    }

    def set_processing_mode(self, mode):
        """Set processing mode: 'full_detection', 'in_place', or 'sn_native'"""
        if mode in ['full_detection', 'in_place', 'sn_native']:
            self.processing_mode = mode
            if mode == 'sn_native':
                print("🎯 S/N native detection mode enabled")
        else:
            raise ValueError("Mode must be 'full_detection', 'in_place', or 'sn_native'")

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
                center - ppm_range * 0.05,  # center can shift ±0.2 ppm for 4.0 ppm window
                0.1,                        # minimum sigma (larger than X)
                0.1,                        # minimum gamma (larger than X)
                baseline - abs(amplitude) * 0.8  # baseline can vary more
            ]
            upper_bounds = [
                amplitude * 3,              # amplitude upper limit
                center + ppm_range * 0.05,  # center can shift ±0.2 ppm for 4.0 ppm window
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
                    lower_bounds, upper_bounds = bounds
                    initial_guess = _clamp_guess_to_bounds(initial_guess, lower_bounds, upper_bounds)
                    popt, pcov = curve_fit(
                        self.voigt_profile_1d, x_data, y_data,
                        p0=initial_guess, bounds=(lower_bounds, upper_bounds),
                        maxfev=self.fitting_parameters['max_iterations']
                    )
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

    def fit_peak_1d_ps2d_style(self, x_data, y_data, target_position, dimension='x', verbose=False):
        """
        Fit 1D Voigt peak using PS2D-style multi-stage fitting

        maximum quality and reliability. This is the RECOMMENDED fitting method.

         expects spectrum to be baseline-corrected during processing (TopSpin/NMRPipe).

        Parameters:
        -----------
        x_data : np.ndarray
            Frequency axis (ppm)
        y_data : np.ndarray
            Intensity data
        target_position : float
            Expected peak position (ppm)
        dimension : str
            'x' for 1H, 'y' for 15N/13C
        verbose : bool
            Print detailed fitting information

        Returns:
        --------
        dict : Fitting results compatible with existing code
        """

        if not PS2D_STYLE_FITTING_AVAILABLE:
            print("⚠️ PS2D-style fitting not available, falling back to standard fitting")
            # Fallback to standard method
            initial_guess = self._calculate_initial_guess_1d(x_data, y_data, target_position)
            bounds = self._get_fitting_bounds(initial_guess, x_data, dimension)
            return self.fit_peak_1d(x_data, y_data, initial_guess, 'voigt', bounds)

        try:
            # Estimate initial linewidth based on dimension
            if dimension.lower() == 'x':
                # 1H dimension - narrower peaks
                initial_linewidth = 0.015  # ppm
            else:
                # 15N/13C dimension - broader peaks
                initial_linewidth = 0.5  # ppm

            # ===================================================================
            # See spectrum.cpp readData2D() line 1100: y.push_back(b);
            # See peakfit.cpp fitGlobal() line 370: yred[j] = y[i];
            # ===================================================================

            # Call PS2D-style fitter with RAW data (NO baseline subtraction)
            result = fit_single_peak_ps2d_style(
                x_data=x_data,
                y_data=y_data,
                peak_position=target_position,
                initial_linewidth=initial_linewidth,
                verbose=verbose
            )

            # Quality threshold for PS2D results (lowered to 0.1 for debugging/QC)
            MIN_R_SQUARED = 0.1

            if result['success'] and result.get('r_squared', 0) >= MIN_R_SQUARED:
                # Convert to standard format expected by rest of code
                return {
                    'success': True,
                    'r_squared': result['r_squared'],
                    'fitted_curve': result['fitted_curve'],  # No baseline adjustment needed
                    'amplitude': result['intensity'],  # PS2D calls it intensity
                    'center': result['pos'],
                    'sigma': result['lw_gauss'] / np.sqrt(8 * np.log(2)),  # Convert FWHM to sigma
                    'gamma': result['lw_lorentz'] / 2.0,  # Convert FWHM to gamma
                    'baseline': 0.0,
                    'profile_type': 'voigt',
                    'method': 'ps2d_style_multistage',
                    'iterations': result['iterations'],
                    'chi2': result['chi2'],
                    'covariance': result['covariance'],
                    'parameters': result['params']
                }
            else:
                # Poor quality or failed - fallback
                if result.get('success'):
                    print(f"⚠️ PS2D fit quality poor (R²={result.get('r_squared', 0):.3f}), falling back to standard")
                else:
                    print("⚠️ PS2D-style fitting failed, falling back to standard")

                initial_guess = self._calculate_initial_guess_1d(x_data, y_data, target_position)
                bounds = self._get_fitting_bounds(initial_guess, x_data, dimension)
                return self.fit_peak_1d(x_data, y_data, initial_guess, 'voigt', bounds)

        except Exception as e:
            print(f"⚠️ PS2D-style fitting error: {e}, falling back to standard")
            initial_guess = self._calculate_initial_guess_1d(x_data, y_data, target_position)
            bounds = self._get_fitting_bounds(initial_guess, x_data, dimension)
            return self.fit_peak_1d(x_data, y_data, initial_guess, 'voigt', bounds)

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
        min_distance = max(1, int(len(y_smooth) / distance_factor))  # Ensure distance ≥ 1 for scipy
        min_prominence = max(signal_max * prominence_threshold, noise_level * 2)

        # Primary peak detection with relaxed criteria
        peaks, properties = find_peaks(y_smooth,
                                     height=min_height,
                                     distance=min_distance,
                                     prominence=min_prominence)

        print(f"   Initial detection: {len(peaks)} peaks (height≥{min_height:.0f}, dist≥{min_distance}, prom≥{min_prominence:.0f})")

        # Secondary detection for overlapping peaks (more sensitive)
        # ENHANCED: Always try sensitive detection for better overlap resolution
        min_height_sensitive = max(signal_max * 0.005, noise_level * 1.5)  # Lower threshold (0.5% of max)
        min_distance_sensitive = max(1, len(y_smooth) // 200)              # Allow very close peaks
        min_prominence_sensitive = max(signal_max * 0.005, noise_level)    # Lower prominence

        peaks_sensitive, _ = find_peaks(y_smooth,
                                      height=min_height_sensitive,
                                      distance=min_distance_sensitive,
                                      prominence=min_prominence_sensitive)

        # Use sensitive detection if it finds more peaks (up to reasonable limit)
        if len(peaks_sensitive) > len(peaks) and len(peaks_sensitive) <= 10:
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
            initial_guess = _clamp_guess_to_bounds(initial_guess, bounds_lower, bounds_upper)
            popt, pcov = curve_fit(
                multi_voigt, x_data, y_data,
                p0=initial_guess,
                bounds=(bounds_lower, bounds_upper),
                maxfev=1000  # Reduced to prevent timeout
            )

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
        # DEBUG: Show gui_params at entry to adaptive_fit_1d
        target_str = f"{target_position:.4f}" if target_position is not None else "None"
        print(f"   [ENTRY] adaptive_fit_1d ({dimension}): target_position={target_str}")
        if gui_params:
            print(f"   ✅ adaptive_fit_1d ({dimension}): gui_params exists with keys={list(gui_params.keys())}")
            print(f"      fix_positions={gui_params.get('fix_positions', 'N/A')}, fix_linewidths={gui_params.get('fix_linewidths', 'N/A')}")
        else:
            print(f"   ⚠️ adaptive_fit_1d ({dimension}): gui_params=None")

        # STANDARD METHOD
        # Get max peaks parameter
        max_peaks_fit = gui_params.get('max_peaks_fit', 4) if gui_params else 4

        # First attempt: use GUI parameters
        detected_peaks = self.detect_peaks_1d(x_data, y_data, target_position, gui_params)
        n_detected = len(detected_peaks)

        print(f"   {dimension.upper()}-dimension: {n_detected} peaks detected")

        # Attempt fitting with current parameters
        fit_result = None

        if n_detected <= 1:
            # Single peak or no peaks - USE PS2D-STYLE FITTING (HIGH QUALITY)
            if PS2D_STYLE_FITTING_AVAILABLE:
                print(f"   🚀 Using PS2D-style high-quality fitting for single peak")
                fit_result = self.fit_peak_1d_ps2d_style(x_data, y_data, target_position, dimension, verbose=False)
            else:
                # Fallback to standard fitting
                initial_guess = self._calculate_initial_guess_1d(x_data, y_data, target_position)
                bounds = self._get_fitting_bounds(initial_guess, x_data, dimension)
                fit_result = self.fit_peak_1d(x_data, y_data, initial_guess, 'voigt', bounds)

        elif 2 <= n_detected <= max_peaks_fit:
            # Multi-peak case - USE PS2D-STYLE MULTI-PEAK FITTING (EXCLUSIVE)
            fit_result = None  # Initialize to avoid variable leakage

            # PS2D multi-peak is now the ONLY method (no fallbacks to old methods)
            # DEBUG: Show decision process
            use_ps2d = gui_params.get('use_ps2d_multi_peak', True) if gui_params else True
            print(f"   🔍 Multi-peak decision: PS2D_available={PS2D_MULTI_PEAK_AVAILABLE}, gui_params={'exists' if gui_params else 'None'}, use_ps2d_multi_peak={use_ps2d}")

            if PS2D_MULTI_PEAK_AVAILABLE:
                print(f"   [PINT] Using PS2D-style multi-peak fitting (5-stage) for {n_detected} peaks in {dimension.upper()}-dimension")

                try:
                    # Get linewidth constraint parameters from GUI (with None check)
                    fix_linewidths = gui_params.get('fix_linewidths', False) if gui_params else False
                    fix_positions = gui_params.get('fix_positions', False) if gui_params else False

                    # Get custom linewidth overrides from GUI (if provided)
                    lw_lorentz_1h = gui_params.get('lw_lorentz_1h', None) if gui_params else None
                    lw_gauss_1h = gui_params.get('lw_gauss_1h', None) if gui_params else None
                    lw_lorentz_15n = gui_params.get('lw_lorentz_15n', None) if gui_params else None
                    lw_gauss_15n = gui_params.get('lw_gauss_15n', None) if gui_params else None

                    # DEBUG: Show constraint values being passed to PS2D
                    print(f"      🔒 Constraints: fix_pos={fix_positions}, fix_lw={fix_linewidths}")
                    print(f"      📏 Custom LW: 1H_L={lw_lorentz_1h}, 1H_G={lw_gauss_1h}, 15N_L={lw_lorentz_15n}, 15N_G={lw_gauss_15n}")

                    # ===================================================================
                    # USE EXACT C++ PS2D OVERLAP DETECTION (ellipsoid + transitive closure)
                    # This implements the complete C++ workflow:
                    # 1. Ellipsoid overlap detection (geometric intersection)
                    # 2. Transitive closure grouping (A overlaps B, B overlaps C → all grouped)
                    # 3. Simultaneous fitting of each overlap group (5-stage deconvolution)
                    # ===================================================================

                    # Check if exact overlap detection is available
                    if PS2D_EXACT_OVERLAP_INTEGRATION_AVAILABLE:
                        print(f"   [EXACT OVERLAP] Using PS2D exact overlap detector (C++ algorithm)")

                        # Call EXACT overlap detection + fitting
                        print(f"   [DEBUG] Passing target_position={target_position:.4f} to PS2D")
                        ps2d_multi_result = fit_peaks_with_exact_overlap_detection(
                            x_data=x_data,
                            y_data=y_data,
                            detected_peaks=detected_peaks,
                            dimension=dimension,
                            target_position=target_position,  # CRITICAL: Pass original target for correct peak selection
                            fix_linewidths=fix_linewidths,
                            fix_positions=fix_positions,
                            lw_lorentz_1h=lw_lorentz_1h,
                            lw_gauss_1h=lw_gauss_1h,
                            lw_lorentz_15n=lw_lorentz_15n,
                            lw_gauss_15n=lw_gauss_15n,
                            use_exact_overlap_detection=True,  # ENABLE exact C++ algorithm
                            verbose=False  # Enable debug output
                        )

                        # Show overlap detection results
                        if ps2d_multi_result.get('n_overlap_groups', 0) > 0:
                            print(f"      ✅ Found {ps2d_multi_result['n_overlap_groups']} overlap groups")
                            print(f"      ✅ Found {ps2d_multi_result['n_isolated_peaks']} isolated peaks")
                    else:
                        # Fallback to standard multi-peak fitting (NO overlap detection)
                        print(f"   [FALLBACK] Exact overlap detector not available - using standard fitting")
                        print(f"   [DEBUG] Passing target_position={target_position:.4f} to PS2D fallback")
                        ps2d_multi_result = fit_overlapping_peaks_ps2d_style(
                            x_data=x_data,
                            y_data=y_data,
                            detected_peaks=detected_peaks,
                            dimension=dimension,
                            target_position=target_position,  # CRITICAL: Pass original target for correct peak selection
                            fix_linewidths=fix_linewidths,
                            fix_positions=fix_positions,
                            lw_lorentz_1h=lw_lorentz_1h,
                            lw_gauss_1h=lw_gauss_1h,
                            lw_lorentz_15n=lw_lorentz_15n,
                            lw_gauss_15n=lw_gauss_15n,
                            verbose=False  # Enable debug output
                        )

                    if ps2d_multi_result['success']:
                        print(f"   ✅ PS2D multi-peak: {len(ps2d_multi_result['peaks'])} peaks, R²={ps2d_multi_result['r_squared']:.3f}")

                        # DEBUG: Show all fitted peaks before selection
                        print(f"      🎯 Target position: {target_position:.4f} ppm")
                        for i, pk in enumerate(ps2d_multi_result['peaks']):
                            dist = abs(pk['position'] - target_position)
                            print(f"         Peak {i+1}: pos={pk['position']:.4f}, intensity={pk['intensity']:.1f}, dist={dist:.4f} ppm")

                        # Find peak closest to target position
                        target_peak = min(ps2d_multi_result['peaks'],
                                        key=lambda p: abs(p['position'] - target_position))
                        print(f"      ✓ Selected peak: pos={target_peak['position']:.4f} (shift={target_peak['position']-target_position:.4f} ppm)")

                        fit_result = {
                            'success': True,
                            'r_squared': ps2d_multi_result['r_squared'],
                            'fitted_curve': ps2d_multi_result['fitted_curve'],  # NO baseline adjustment
                            'amplitude': target_peak['intensity'],
                            'center': target_peak['position'],
                            'sigma': target_peak['lw_gauss'] / np.sqrt(8 * np.log(2)),  # Convert FWHM to sigma
                            'gamma': target_peak['lw_lorentz'] / 2.0,  # Convert FWHM to gamma
                            'baseline': 0.0,
                            'multi_peak_info': {
                                'n_peaks': len(ps2d_multi_result['peaks']),
                                'all_peaks': ps2d_multi_result['peaks'],
                                'method': 'ps2d_multi_peak_5_stage',
                                'convergence_flag': ps2d_multi_result.get('convergence_flag', 0)
                            }
                        }
                    else:
                        print(f"   ⚠️ PS2D multi-peak failed (success=False)")
                        fit_result = None

                except Exception as e:
                    print(f"   ❌ PS2D multi-peak exception: {e}")
                    fit_result = None
            else:
                print(f"   ❌ PS2D multi-peak not available - cannot fit {n_detected} peaks")
                fit_result = None

            # NO FALLBACKS - PS2D is the only method
            # If PS2D failed, the fit_result will be None and fitting fails cleanly

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

    def _safe_extract_params(self, fit_result, default_center=0.0):
        """
        Safely extract fitting parameters from any result format.

        Handles:
        - Flat format: {'amplitude': X, 'center': Y, ...}
        - Nested format: {'x_fit': {'amplitude': X}, 'y_fit': {...}}
        - Missing keys: Returns safe defaults
        - None results: Returns safe defaults

        Args:
            fit_result: Result dict from any fitting method
            default_center: Default center value if missing

        Returns:
            Dict with all required keys: amplitude, center, sigma, gamma, baseline, r_squared, fitted_curve
        """
        if fit_result is None:
            return {
                'amplitude': 1000.0,
                'center': default_center,
                'sigma': 0.01,
                'gamma': 0.01,
                'baseline': 0.0,
                'r_squared': 0.0,
                'fitted_curve': [],
                'success': False,
                'error': 'fit_result_is_none'
            }

        # Try flat format first
        if 'amplitude' in fit_result:
            return {
                'amplitude': fit_result.get('amplitude', 1000.0),
                'center': fit_result.get('center', default_center),
                'sigma': fit_result.get('sigma', 0.01),
                'gamma': fit_result.get('gamma', 0.01),
                'baseline': fit_result.get('baseline', 0.0),
                'r_squared': fit_result.get('r_squared', 0.0),
                'fitted_curve': fit_result.get('fitted_curve', []),
                'success': fit_result.get('success', True),
                'quality_metrics': fit_result.get('quality_metrics', {})
            }

        # Try nested format (x_fit/y_fit structure)
        elif 'x_fit' in fit_result:
            x_fit = fit_result['x_fit']
            return {
                'amplitude': x_fit.get('amplitude', 1000.0),
                'center': x_fit.get('center', default_center),
                'sigma': x_fit.get('sigma', 0.01),
                'gamma': x_fit.get('gamma', 0.01),
                'baseline': x_fit.get('baseline', 0.0),
                'r_squared': fit_result.get('avg_r_squared', x_fit.get('r_squared', 0.0)),
                'fitted_curve': x_fit.get('fitted_curve', []),
                'success': True,
                'quality_metrics': x_fit.get('quality_metrics', {})
            }

        # Try overlap resolution format (has 'peaks' list and 'fitted_params')
        elif 'peaks' in fit_result and isinstance(fit_result['peaks'], list) and len(fit_result['peaks']) > 0:
            # Find peak closest to default_center
            best_peak = min(fit_result['peaks'], key=lambda p: abs(p.get('center', default_center) - default_center))
            return {
                'amplitude': best_peak.get('amplitude', 1000.0),
                'center': best_peak.get('center', default_center),
                'sigma': best_peak.get('sigma', 0.01),
                'gamma': best_peak.get('gamma', 0.01),
                'baseline': fit_result.get('baseline', 0.0),
                'r_squared': fit_result.get('r_squared', 0.0),
                'fitted_curve': fit_result.get('fitted_curve', []),
                'success': fit_result.get('success', True),
                'quality_metrics': {},
                'multi_peak_info': {
                    'n_peaks': fit_result.get('n_peaks', len(fit_result['peaks'])),
                    'all_peaks': fit_result['peaks'],
                    'method': 'overlap_resolution'
                }
            }

        # Unknown format - return defaults with warning
        else:
            print(f"   ⚠️ WARNING: Unknown fit result format, using defaults. Keys: {list(fit_result.keys())}")
            return {
                'amplitude': 1000.0,
                'center': default_center,
                'sigma': 0.01,
                'gamma': 0.01,
                'baseline': 0.0,
                'r_squared': fit_result.get('r_squared', 0.0),
                'fitted_curve': fit_result.get('fitted_curve', []),
                'success': fit_result.get('success', False),
                'error': 'unknown_format'
            }

    def _check_peaks_overlap(self, peak1_pos, peak2_pos,
                             overlap_threshold_x=None,
                             overlap_threshold_y=None) -> bool:
        """
        Check if two peaks' elliptical overlap regions touch or overlap

        Each peak has an elliptical region with semi-axes = overlap_threshold_x/y.
        Two ellipses touch when the normalized elliptical distance <= 2.0.

        Uses proper elliptical distance calculation to handle diagonal overlaps correctly.
        This matches the orange ellipse visualization where touching ellipses should cluster.

        Algorithm:
        1. Compute dx = |x2 - x1|, dy = |y2 - y1|
        2. Normalize: dx_norm = dx/threshold_x, dy_norm = dy/threshold_y
        3. Elliptical distance = sqrt(dx_norm² + dy_norm²)
        4. Overlap if distance <= 2.0 (sum of radii in normalized space)

        This correctly handles:
        - Horizontal overlaps (small dx, large dy OK)
        - Vertical overlaps (small dy, large dx OK)
        - Diagonal overlaps (both dx and dy moderate but distance <= 2.0)

        Uses centralized PS2D configuration (nucleus-adaptive).
        Values automatically scale for 15N vs 13C spectra.

        Parameters:
        -----------
        peak1_pos, peak2_pos : tuple of (x_ppm, y_ppm)
            Peak positions to test
        overlap_threshold_x, overlap_threshold_y : float, optional
            Overlap radii in each dimension (default from PS2D config)

        Returns:
        --------
        bool : True if peaks' elliptical regions overlap
        """
        # Use centralized configuration
        config = get_ps2d_config()

        # Get thresholds (use provided values or read from config)
        if overlap_threshold_x is None:
            overlap_threshold_x = config.overlap_threshold_x
        if overlap_threshold_y is None:
            overlap_threshold_y = config.overlap_threshold_y

        # ALWAYS read multiplier from config (CRITICAL for diagonal overlap detection)
        multiplier = getattr(config, 'overlap_distance_multiplier', 1.0)

        x1, y1 = peak1_pos
        x2, y2 = peak2_pos

        # Two ellipses touch when elliptical distance <= 2.0 * multiplier
        # Uses normalized elliptical distance to handle diagonal overlaps correctly
        # This matches the orange ellipse visualization (proper elliptical regions)
        # Multiplier > 1.0 makes overlap detection more aggressive (catches more diagonal overlaps)
        dx = abs(x2 - x1)
        dy = abs(y2 - y1)

        # Normalize distances by threshold radii (convert to "radius units")
        dx_normalized = dx / overlap_threshold_x
        dy_normalized = dy / overlap_threshold_y

        # Compute elliptical distance
        elliptical_distance = (dx_normalized**2 + dy_normalized**2)**0.5

        # Two ellipses of radius 1.0 each touch when distance <= 2.0
        # Apply multiplier to make detection more/less aggressive
        distance_threshold = 2.0 * multiplier
        overlaps = elliptical_distance <= distance_threshold

        return overlaps

    def identify_overlap_clusters(self, all_peaks_context,
                                   overlap_threshold_x=None,
                                   overlap_threshold_y=None,
                                   max_cluster_diameter_x=0.20,
                                   max_cluster_diameter_y=1.5,
                                   max_cluster_size=None) -> List[List[Tuple[float, float]]]:
        """
        Identify clusters of overlapping peaks using hierarchical clustering
        with spatial and size constraints

        Algorithm:
        1. Start with each peak as singleton cluster
        2. Iteratively merge closest overlapping clusters
        3. Stop merge if diameter exceeds threshold (prevents sprawl)
        4. Split oversized clusters (>max_cluster_size)

        This prevents transitive closure from creating mega-clusters while
        still grouping genuinely overlapping peaks.

        Uses centralized PS2D configuration (nucleus-adaptive).

        Parameters:
        -----------
        all_peaks_context : list of dict
            All peaks with 'x_ppm'/'pos_x' and 'y_ppm'/'pos_y' keys
        overlap_threshold_x, overlap_threshold_y : float, optional
            Overlap detection radii for merging (default from PS2D config)
        max_cluster_diameter_x, max_cluster_diameter_y : float
            Maximum spatial extent of cluster (prevents sprawl)
        max_cluster_size : int, optional
            Maximum peaks per cluster (default from PS2D config)

        Returns:
        --------
        list of clusters : Each cluster is list of (x_ppm, y_ppm) tuples
        """
        # Use centralized configuration
        if overlap_threshold_x is None or overlap_threshold_y is None or max_cluster_size is None:
            config = get_ps2d_config()
            if overlap_threshold_x is None:
                overlap_threshold_x = config.overlap_threshold_x
            if overlap_threshold_y is None:
                overlap_threshold_y = config.overlap_threshold_y
            if max_cluster_size is None:
                max_cluster_size = config.max_cluster_size

        # DIAGNOSTIC: Show overlap threshold values being used for clustering
        config = get_ps2d_config()
        multiplier = getattr(config, 'overlap_distance_multiplier', 1.0)
        print(f"🔍 OVERLAP DETECTION CONFIGURATION:")
        print(f"   overlap_threshold_x (1H): {overlap_threshold_x:.4f} ppm")
        print(f"   overlap_threshold_y (15N/13C): {overlap_threshold_y:.4f} ppm")
        print(f"   overlap_distance_multiplier: {multiplier:.2f} {'(default)' if multiplier == 1.0 else '(INCREASED for diagonal overlaps)'}")
        print(f"   Elliptical distance threshold: {2.0 * multiplier:.2f} (2.0 × {multiplier:.2f})")
        print(f"   max_cluster_size: {max_cluster_size}")

        if not all_peaks_context or len(all_peaks_context) == 0:
            return []

        # Extract peak positions
        peak_positions = []
        for peak in all_peaks_context:
            x = peak.get('x_ppm') or peak.get('pos_x')
            y = peak.get('y_ppm') or peak.get('pos_y')
            if x is not None and y is not None:
                peak_positions.append((x, y))

        if len(peak_positions) == 0:
            return []

        # DIAGNOSTIC: Show first few peaks being clustered
        print(f"📍 Peak positions for clustering (showing first 5 of {len(peak_positions)}):")
        for i, (x, y) in enumerate(peak_positions[:5]):
            print(f"   Peak {i+1}: ({x:.4f}, {y:.3f}) ppm")
        if len(peak_positions) > 5:
            print(f"   ... and {len(peak_positions)-5} more peaks")

        # Start with each peak as its own cluster
        clusters = [[peak] for peak in peak_positions]

        # Hierarchical merging with diameter constraints
        while True:

            # Find closest pair of overlapping clusters
            best_pair = None
            best_distance = float('inf')

            for i in range(len(clusters)):
                for j in range(i + 1, len(clusters)):
                    # Check if any peaks from cluster i overlap with cluster j
                    overlaps = False
                    min_dist = float('inf')

                    for peak_i in clusters[i]:
                        for peak_j in clusters[j]:
                            if self._check_peaks_overlap(peak_i, peak_j,
                                                        overlap_threshold_x,
                                                        overlap_threshold_y):
                                overlaps = True
                                # Compute distance for sorting
                                dist = ((peak_i[0] - peak_j[0])**2 + (peak_i[1] - peak_j[1])**2)**0.5
                                min_dist = min(min_dist, dist)

                    # CRITICAL FIX: Check diameter constraint BEFORE selecting as best pair
                    # This prevents blocking all merges when one pair violates diameter
                    if overlaps:
                        # Test if merge would violate diameter
                        merged = clusters[i] + clusters[j]
                        x_coords = [p[0] for p in merged]
                        y_coords = [p[1] for p in merged]
                        diameter_x = max(x_coords) - min(x_coords)
                        diameter_y = max(y_coords) - min(y_coords)

                        # Only accept as candidate if diameter is OK
                        if (diameter_x <= max_cluster_diameter_x and
                            diameter_y <= max_cluster_diameter_y and
                            min_dist < best_distance):
                            best_distance = min_dist
                            best_pair = (i, j)

            # No more mergeable clusters (that satisfy diameter constraint)
            if best_pair is None:
                break

            # Accept merge (diameter already verified during pair selection)
            i, j = best_pair
            merged = clusters[i] + clusters[j]
            clusters[i] = merged
            clusters.pop(j)

        # Enforce max cluster size by splitting large clusters
        final_clusters = []
        for cluster in clusters:
            if len(cluster) <= max_cluster_size:
                final_clusters.append(cluster)
            else:
                # Split by spatial gaps (thresholds from centralized config)
                # Uses gap detection to split large clusters
                subclusters = self._split_cluster_by_gaps(cluster, max_size=max_cluster_size)
                final_clusters.extend(subclusters)

        # Log clustering results
        n_isolated = sum(1 for c in final_clusters if len(c) == 1)
        n_groups = len(final_clusters) - n_isolated
        max_group_size = max(len(c) for c in final_clusters) if final_clusters else 0

        print(f"📊 Overlap clustering: {len(final_clusters)} clusters from {len(peak_positions)} peaks")
        print(f"   • {n_isolated} isolated peaks")
        print(f"   • {n_groups} overlap groups (max size: {max_group_size})")
        print(f"   🔒 Constraints: diameter ≤ ({max_cluster_diameter_x:.2f}, {max_cluster_diameter_y:.1f}) ppm, size ≤ {max_cluster_size}")

        return final_clusters

    def _split_cluster_by_gaps(self, cluster, gap_threshold_x=None, gap_threshold_y=None, max_size=None):
        """
        Split oversized cluster by finding spatial gaps.
        If no suitable gap exists, force split at median.

        Algorithm:
        1. Sort peaks by F2 (1H) position, find largest gap
        2. If gap ≥ threshold, split there
        3. Otherwise try F1 (15N) position, find largest gap
        4. If gap ≥ threshold, split there
        5. If both fail AND cluster > max_size, force split at median
        6. Recursively split each subcluster if still too large

        Uses centralized PS2D configuration (nucleus-adaptive).

        Parameters:
        -----------
        cluster : list of (x, y) tuples
            Oversized cluster to split
        gap_threshold_x : float, optional
            Minimum gap size in F2 (1H) to trigger split (default from config)
        gap_threshold_y : float, optional
            Minimum gap size in F1 (15N/13C) to trigger split (default from config)
        max_size : int, optional
            Maximum peaks per cluster (default from config)

        Returns:
        --------
        list of subclusters : Each subcluster is list of (x, y) tuples
        """
        # Use centralized configuration
        if gap_threshold_x is None or gap_threshold_y is None or max_size is None:
            config = get_ps2d_config()
            if gap_threshold_x is None:
                gap_threshold_x = config.gap_threshold_x
            if gap_threshold_y is None:
                gap_threshold_y = config.gap_threshold_y
            if max_size is None:
                max_size = config.max_cluster_size

        if len(cluster) <= max_size:
            return [cluster]

        # Sort by F2 (1H) position
        sorted_cluster_f2 = sorted(cluster, key=lambda p: p[0])

        # Find largest gap in F2
        max_gap_f2 = 0
        split_idx_f2 = None

        for i in range(len(sorted_cluster_f2) - 1):
            gap = sorted_cluster_f2[i + 1][0] - sorted_cluster_f2[i][0]
            if gap > max_gap_f2:
                max_gap_f2 = gap
                split_idx_f2 = i + 1

        # Strategy 1: Split at F2 gap if significant
        if max_gap_f2 >= gap_threshold_x and split_idx_f2 is not None:
            subcluster1 = sorted_cluster_f2[:split_idx_f2]
            subcluster2 = sorted_cluster_f2[split_idx_f2:]
            # Recursively split subclusters (pass parameters to maintain consistency)
            result = []
            result.extend(self._split_cluster_by_gaps(subcluster1, gap_threshold_x, gap_threshold_y, max_size))
            result.extend(self._split_cluster_by_gaps(subcluster2, gap_threshold_x, gap_threshold_y, max_size))
            return result

        # Strategy 2: Try splitting by F1 if F2 split didn't work
        sorted_cluster_f1 = sorted(cluster, key=lambda p: p[1])
        max_gap_f1 = 0
        split_idx_f1 = None

        for i in range(len(sorted_cluster_f1) - 1):
            gap = sorted_cluster_f1[i + 1][1] - sorted_cluster_f1[i][1]
            if gap > max_gap_f1:
                max_gap_f1 = gap
                split_idx_f1 = i + 1

        if max_gap_f1 >= gap_threshold_y and split_idx_f1 is not None:
            subcluster1 = sorted_cluster_f1[:split_idx_f1]
            subcluster2 = sorted_cluster_f1[split_idx_f1:]
            # Recursively split subclusters
            result = []
            result.extend(self._split_cluster_by_gaps(subcluster1, gap_threshold_x, gap_threshold_y, max_size))
            result.extend(self._split_cluster_by_gaps(subcluster2, gap_threshold_x, gap_threshold_y, max_size))
            return result

        # Strategy 3: Force split at median if cluster still too large
        # This prevents pathological cases where no gaps exist but cluster > max_size
        if len(cluster) > max_size:
            mid = len(cluster) // 2
            print(f"   ⚠️  No spatial gap found (F2: {max_gap_f2:.3f} < {gap_threshold_x:.3f} ppm, F1: {max_gap_f1:.3f} < {gap_threshold_y:.3f} ppm)")
            print(f"   🔪 Forcing median split: {len(cluster)} peaks → {mid} + {len(cluster)-mid} peaks")
            # Use F2-sorted order for median split (1H dimension is typically more informative)
            subcluster1 = sorted_cluster_f2[:mid]
            subcluster2 = sorted_cluster_f2[mid:]
            # Recursively split subclusters
            result = []
            result.extend(self._split_cluster_by_gaps(subcluster1, gap_threshold_x, gap_threshold_y, max_size))
            result.extend(self._split_cluster_by_gaps(subcluster2, gap_threshold_x, gap_threshold_y, max_size))
            return result

        # Cluster is small enough now
        return [cluster]

    def check_if_peaks_need_2d_fitting(self, peak_x_ppm, peak_y_ppm,
                                        all_peaks_context=None,
                                        overlap_threshold_x=None,
                                        overlap_threshold_y=None) -> Tuple[bool, List]:
        """
        Check if peak requires 2D simultaneous fitting due to overlapping neighbors

        Uses geometric overlap detection: peaks overlap if their elliptical windows intersect.
        Uses centralized PS2D configuration (nucleus-adaptive).

        Parameters:
        -----------
        peak_x_ppm, peak_y_ppm : float
            Target peak position
        all_peaks_context : list of dict, optional
            All detected peaks with positions
        overlap_threshold_x : float, optional
            Overlap radius in X dimension (default from config)
        overlap_threshold_y : float, optional
            Overlap radius in Y dimension (default from config)

        Returns:
        --------
        needs_2d : bool
            True if peak has overlapping neighbors requiring 2D fitting
        overlap_group : list of dict
            Peak dictionaries for all overlapping peaks (including target)
            Each dict contains: 'x_ppm', 'y_ppm', 'intensity' (optional)
        """
        if not PS2D_2D_FITTER_AVAILABLE or all_peaks_context is None:
            return False, []

        # Use centralized configuration
        if overlap_threshold_x is None or overlap_threshold_y is None:
            config = get_ps2d_config()
            if overlap_threshold_x is None:
                overlap_threshold_x = config.overlap_threshold_x
            if overlap_threshold_y is None:
                overlap_threshold_y = config.overlap_threshold_y

        # Find the target peak's full dictionary from context
        target_peak_dict = None
        if all_peaks_context:
            for peak in all_peaks_context:
                px = peak.get('x_ppm') or peak.get('pos_x')
                py = peak.get('y_ppm') or peak.get('pos_y')
                if px and py and abs(px - peak_x_ppm) < 0.001 and abs(py - peak_y_ppm) < 0.01:
                    target_peak_dict = peak
                    break

        # If not found in context, create minimal dict with positions only
        if target_peak_dict is None:
            target_peak_dict = {
                'x_ppm': peak_x_ppm,
                'y_ppm': peak_y_ppm,
                'pos_x': peak_x_ppm,
                'pos_y': peak_y_ppm,
                'intensity': None  # Will trigger fallback re-measurement
            }

        # Start overlap group with target peak dictionary
        overlap_group = [target_peak_dict]

        for other_peak in all_peaks_context:
            other_x = other_peak.get('x_ppm') or other_peak.get('pos_x')
            other_y = other_peak.get('y_ppm') or other_peak.get('pos_y')

            if other_x is None or other_y is None:
                continue

            # Skip self (already added as target)
            if abs(other_x - peak_x_ppm) < 0.001 and abs(other_y - peak_y_ppm) < 0.01:
                continue

            # Elliptical overlap test (normalized distance)
            normalized_dist_sq = ((other_x - peak_x_ppm) / overlap_threshold_x) ** 2 + \
                                 ((other_y - peak_y_ppm) / overlap_threshold_y) ** 2

            if normalized_dist_sq <= 1.0:
                # Overlaps! Append full peak dictionary
                overlap_group.append(other_peak)

        # Need 2D fitting if there are overlapping neighbors
        needs_2d = len(overlap_group) > 1

        return needs_2d, overlap_group

    def _estimate_fwhm_from_slice(self, intensity_slice, ppm_scale, peak_ppm):
        """
        Estimate FWHM (Full Width at Half Maximum) from 1D cross-section

        Parameters:
        -----------
        intensity_slice : np.ndarray
            1D intensity profile
        ppm_scale : np.ndarray
            PPM scale for the slice
        peak_ppm : float
            Expected peak position (ppm)

        Returns:
        --------
        float : Estimated FWHM in ppm
        """
        # Find peak maximum near expected position
        peak_idx = np.argmin(np.abs(ppm_scale - peak_ppm))
        peak_intensity = intensity_slice[peak_idx]

        # Find half-maximum points
        half_max = peak_intensity / 2.0

        # Search left from peak
        left_idx = peak_idx
        while left_idx > 0 and intensity_slice[left_idx] > half_max:
            left_idx -= 1

        # Search right from peak
        right_idx = peak_idx
        while right_idx < len(intensity_slice) - 1 and intensity_slice[right_idx] > half_max:
            right_idx += 1

        # Calculate FWHM
        if left_idx == 0 or right_idx == len(intensity_slice) - 1:
            # Peak extends to edge - use nucleus-adaptive default from config
            # Use 2× min_linewidth as typical FWHM estimate
            ppm_range = abs(ppm_scale[-1] - ppm_scale[0])
            if ppm_range > 5:  # Indirect dimension (15N or 13C)
                config = get_ps2d_config()
                return 2.0 * config.min_linewidth_f1
            else:  # 1H dimension
                return 0.02

        fwhm = abs(ppm_scale[right_idx] - ppm_scale[left_idx])

        # Sanity check: FWHM should be reasonable
        # Use nucleus-adaptive constraints from centralized config
        ppm_range = abs(ppm_scale[-1] - ppm_scale[0])
        if ppm_range > 5:  # Indirect dimension (15N or 13C)
            # Get nucleus-specific minimum from config
            config = get_ps2d_config()
            min_fwhm = config.min_linewidth_f1 * 0.5  # Allow measured FWHM down to 0.5× min
            max_fwhm = 2.0  # Maximum reasonable FWHM
            fwhm = np.clip(fwhm, min_fwhm, max_fwhm)
        else:  # 1H dimension
            fwhm = np.clip(fwhm, 0.002, 0.1)

        return fwhm

    def extract_2d_region_for_overlap_group(self, overlap_group,
                                             radF1=None, radF2=None):
        """
        Extract 2D region covering all peaks in overlap group

        PS2D-compatible implementation: Uses elliptical window radii to define
        bounding box, matching spectrum.cpp:1584-1587 (getRegion function).

        Uses centralized PS2D configuration (nucleus-adaptive).

        Parameters:
        -----------
        overlap_group : list of dict
            Peak dictionaries with 'x_ppm'/'pos_x' and 'y_ppm'/'pos_y' keys
        radF1, radF2 : float, optional
            Elliptical window radii (default from PS2D config)

        Returns:
        --------
        dict : Extracted region with keys:
            'f1_grid', 'f2_grid': 2D coordinate grids (NMRPipe: F1=Y, F2=X)
            'intensity': 2D intensity data
            'f1_ppm', 'f2_ppm': 1D axes
        """
        # Use centralized configuration
        if radF1 is None or radF2 is None:
            config = get_ps2d_config()
            if radF1 is None:
                radF1 = config.radF1
            if radF2 is None:
                radF2 = config.radF2
        # Find bounding box using PS2D approach
        # Each peak contributes region: (pos - rad) to (pos + rad)
        # Final box is union of all peak regions
        x_positions = [p.get('x_ppm') or p.get('pos_x') for p in overlap_group]
        y_positions = [p.get('y_ppm') or p.get('pos_y') for p in overlap_group]

        x_min = min(x_positions) - radF2  # PS2D: peak.f2 - peak.radF2
        x_max = max(x_positions) + radF2  # PS2D: peak.f2 + peak.radF2
        y_min = min(y_positions) - radF1  # PS2D: peak.f1 - peak.radF1
        y_max = max(y_positions) + radF1  # PS2D: peak.f1 + peak.radF1

        # Find indices in full spectrum
        x_mask = (self.ppm_x_axis >= x_min) & (self.ppm_x_axis <= x_max)
        y_mask = (self.ppm_y_axis >= y_min) & (self.ppm_y_axis <= y_max)

        x_indices = np.where(x_mask)[0]
        y_indices = np.where(y_mask)[0]

        if len(x_indices) == 0 or len(y_indices) == 0:
            return None

        # Extract region (NMRPipe convention: data[y, x])
        region_data = self.nmr_data[y_indices.min():y_indices.max()+1,
                                     x_indices.min():x_indices.max()+1]
        region_x = self.ppm_x_axis[x_indices.min():x_indices.max()+1]
        region_y = self.ppm_y_axis[y_indices.min():y_indices.max()+1]

        # Create 2D grids (NMRPipe: F1=Y/15N, F2=X/1H)
        f2_grid, f1_grid = np.meshgrid(region_x, region_y, indexing='xy')

        return {
            'f1_grid': f1_grid,      # Y-axis (15N)
            'f2_grid': f2_grid,      # X-axis (1H)
            'intensity': region_data,
            'f1_ppm': region_y,
            'f2_ppm': region_x
        }

    def _reconstruct_2d_surface(self, region_2d, fitted_peaks, display_multiplier=None):
        """
        Reconstruct 2D Voigt surface from PS2D fitted parameters

        Parameters:
        -----------
        region_2d : dict
            Contains f1_grid, f2_grid, intensity arrays
        fitted_peaks : list of dict
            Fitted peak parameters from PS2D 2D fitter
        display_multiplier : float or None
            Multiplier for fit ellipse radius when clipping individual peaks.
            None = full extent (no clipping)
            2.5 = clip to 2.5× fit ellipse radius (recommended for visualization)

        Returns:
        --------
        tuple : (total_surface, individual_surfaces, baseline)
            total_surface: np.ndarray - Sum of all peaks (without baseline)
            individual_surfaces: list of np.ndarray - Each peak surface separately (without baseline)
            baseline: float - Estimated baseline from region edges
        """
        from scipy.special import wofz

        #print(f"   📊 Reconstructing 2D surface from {len(fitted_peaks)} peaks")

        # Estimate baseline from edges of region (where there are no peaks)
        intensity = region_2d['intensity']
        edge_pixels = np.concatenate([
            intensity[0, :],   # Top edge
            intensity[-1, :],  # Bottom edge
            intensity[:, 0],   # Left edge
            intensity[:, -1]   # Right edge
        ])
        # Use median of edges to avoid outliers
        baseline = np.median(edge_pixels)
        #print(f"      Baseline estimated from edges: {baseline:.2e}")

        # Initialize fitted_surface with zeros (baseline removed per user request 2025-11-19)
        fitted_surface = np.zeros_like(region_2d['intensity'], dtype=np.float64)
        individual_surfaces = []  # Store each peak separately

        # Add each peak's contribution
        for i, peak in enumerate(fitted_peaks):
            # Use PS2D's fitted intensity parameter with EXACT PS2D formula
            # This must match ps2d_style_fitter.py:107 exactly
            # For single peaks from 1D fitting, 'intensity' may not exist - use 'amplitude' as fallback
            fitted_intensity = peak.get('intensity', peak.get('amplitude', 1000))

            # Convert FWHM to sigma/gamma (EXACT PS2D formula)
            SQRT_8LN2 = np.sqrt(8 * np.log(2))
            sigma_f1 = peak['lw_gau_f1'] / SQRT_8LN2
            sigma_f2 = peak['lw_gau_f2'] / SQRT_8LN2
            gamma_f1 = peak['lw_lor_f1'] / 2.0
            gamma_f2 = peak['lw_lor_f2'] / 2.0

            # Prevent division by zero
            sigma_f1 = max(sigma_f1, 1e-10)
            sigma_f2 = max(sigma_f2, 1e-10)

            # Build complex arguments for Faddeeva function (EXACT PS2D formula)
            SQRT_2 = np.sqrt(2)
            z_f1 = ((peak['pos_f1'] - region_2d['f1_grid']) + 1j * gamma_f1) / (sigma_f1 * SQRT_2)
            z_f2 = ((peak['pos_f2'] - region_2d['f2_grid']) + 1j * gamma_f2) / (sigma_f2 * SQRT_2)

            # Compute Faddeeva function
            fade_f1 = np.real(wofz(z_f1))
            fade_f2 = np.real(wofz(z_f2))

            # 2D Voigt using EXACT PS2D formula (ps2d_style_fitter.py:107)
            peak_surface = fitted_intensity * fade_f1 * fade_f2 / (sigma_f1 * sigma_f2 * 2.0 * np.pi)

            # Apply elliptical clipping if requested (for visualization clarity)
            if display_multiplier is not None:
                from .ps2d_config import get_ps2d_config
                config = get_ps2d_config()
                radF1_display = config.radF1 * display_multiplier
                radF2_display = config.radF2 * display_multiplier

                # Create elliptical mask
                ellipse_distance_sq = ((region_2d['f1_grid'] - peak['pos_f1']) / radF1_display) ** 2 + \
                                      ((region_2d['f2_grid'] - peak['pos_f2']) / radF2_display) ** 2
                ellipse_mask = ellipse_distance_sq <= 1.0

                # Zero out peak outside display ellipse
                peak_surface = np.where(ellipse_mask, peak_surface, 0.0)

            peak_max = np.max(peak_surface)
            #print(f"      Peak {i+1}: PS2D intensity={fitted_intensity:.2e}, reconstructed max={peak_max:.2e}")

            # CRITICAL FIX: Override the height field with reconstructed max
            # The calculate_peak_height formula doesn't match actual peak height
            # Use the actual reconstructed maximum as the true height
            fitted_peaks[i]['height'] = peak_max
            fitted_peaks[i]['amplitude'] = peak_max  # Amplitude = height in NMR convention

            # Store individual peak surface (without baseline)
            individual_surfaces.append(peak_surface.copy())

            # Add to total surface
            fitted_surface += peak_surface

        #print(f"   ✅ Fitted surface: min={np.min(fitted_surface):.2e}, max={np.max(fitted_surface):.2e}")
        return fitted_surface, individual_surfaces, baseline

    def _extract_cross_section_for_display(self, peak_x_ppm, peak_y_ppm, dimension, peak_params):
        """
        Extract 1D cross-section and reconstruct fitted Voigt curve for GUI display

        Parameters:
        -----------
        peak_x_ppm : float
            Peak X position in ppm
        peak_y_ppm : float
            Peak Y position in ppm
        dimension : str
            'x' or 'y'
        peak_params : dict
            Fitted peak parameters from PS2D 2D fitter

        Returns:
        --------
        dict with 'ppm_scale', 'cross_section', 'fitted_curve'
        """
        from lunaNMR.core.multi_peak_models import multi_voigt_profile

        # Handle both single peak dict and list of peaks
        all_peaks = peak_params if isinstance(peak_params, list) else [peak_params]

        if dimension == 'x':
            # Extract X cross-section at fitted Y position
            y_idx = np.argmin(np.abs(self.ppm_y_axis - peak_y_ppm))

            # Get fitting window from parameters
            window = self.fitting_parameters.get('fitting_window_x', 0.12)
            x_min, x_max = peak_x_ppm - window, peak_x_ppm + window
            x_mask = (self.ppm_x_axis >= x_min) & (self.ppm_x_axis <= x_max)

            ppm_scale = self.ppm_x_axis[x_mask]
            cross_section = self.nmr_data[y_idx, x_mask]

            # Reconstruct ALL overlapping peaks for accurate visualization
            baseline = np.min(cross_section) if len(cross_section) > 0 else 0
            fitted_curve = np.full_like(ppm_scale, baseline, dtype=np.float64)

            # Add contribution from each fitted peak in overlap group
            for peak in all_peaks:
                # Check if peak is within this cross-section window
                if abs(peak['pos_f2'] - peak_x_ppm) < window * 1.5:
                    # Extract actual amplitude from experimental data at fitted position
                    peak_idx = np.argmin(np.abs(ppm_scale - peak['pos_f2']))
                    amplitude = cross_section[peak_idx] - baseline if peak_idx < len(cross_section) else peak['intensity']

                    params = [amplitude, peak['pos_f2'],
                             peak['lw_gau_f2'], peak['lw_lor_f2'], 0.0]
                    fitted_curve += multi_voigt_profile(ppm_scale, *params) - 0.0

        else:  # dimension == 'y'
            # Extract Y cross-section at fitted X position
            x_idx = np.argmin(np.abs(self.ppm_x_axis - peak_x_ppm))

            # Get fitting window from parameters
            window = self.fitting_parameters.get('fitting_window_y', 1.2)
            y_min, y_max = peak_y_ppm - window, peak_y_ppm + window
            y_mask = (self.ppm_y_axis >= y_min) & (self.ppm_y_axis <= y_max)

            ppm_scale = self.ppm_y_axis[y_mask]
            cross_section = self.nmr_data[y_mask, x_idx]

            # Reconstruct ALL overlapping peaks for accurate visualization
            baseline = np.min(cross_section) if len(cross_section) > 0 else 0
            fitted_curve = np.full_like(ppm_scale, baseline, dtype=np.float64)

            # Add contribution from each fitted peak in overlap group
            for peak in all_peaks:
                # Check if peak is within this cross-section window
                if abs(peak['pos_f1'] - peak_y_ppm) < window * 1.5:
                    # Extract actual amplitude from experimental data at fitted position
                    peak_idx = np.argmin(np.abs(ppm_scale - peak['pos_f1']))
                    amplitude = cross_section[peak_idx] - baseline if peak_idx < len(cross_section) else peak['intensity']

                    params = [amplitude, peak['pos_f1'],
                             peak['lw_gau_f1'], peak['lw_lor_f1'], 0.0]
                    fitted_curve += multi_voigt_profile(ppm_scale, *params) - 0.0

        return {
            'ppm_scale': ppm_scale,
            'cross_section': cross_section,
            'fitted_curve': fitted_curve
        }

    def fit_overlap_group_2d(self, overlap_group, assignment="Unknown", peak_assignments=None,
                             fix_positions=False, fix_linewidths=False):
        """
        Fit overlap group using 2D simultaneous multi-peak fitting

        This is the routing function that calls Ps2dMultiPeakFitter2D when
        peaks are too close to resolve with 1D cross-sections.

        Parameters:
        -----------
        overlap_group : list of dict
            Peak dictionaries with keys: 'x_ppm', 'y_ppm', 'intensity' (optional)
            If 'intensity' is None/missing, will fall back to re-measurement
        assignment : str
            Assignment name for logging (primary peak)
        peak_assignments : list of str, optional
            Assignment names for each peak in overlap_group (for table display)
        fix_positions : bool
            If True, peak positions are held constant during fitting (Stage 2 skipped)
        fix_linewidths : bool
            If True, linewidths are held constant during fitting (Stage 1 skipped)

        Returns:
        --------
        dict : Fitting results for all peaks in group, or None if failed
        """
        if not PS2D_2D_FITTER_AVAILABLE:
            print(f"   ⚠️ 2D fitter not available, cannot fit overlap group")
            return None

        print(f"\n   🎯 2D MULTI-PEAK FITTING triggered for {assignment}")
        print(f"   📍 Overlap group: {len(overlap_group)} peaks")
        for i, peak_dict in enumerate(overlap_group):
            x = peak_dict.get('x_ppm') or peak_dict.get('pos_x')
            y = peak_dict.get('y_ppm') or peak_dict.get('pos_y')
            #print(f"      Peak {i+1}: ({x:.4f}, {y:.4f}) ppm")

        # Extract 2D region
        region = self.extract_2d_region_for_overlap_group(overlap_group)
        if region is None:
            print(f"   ❌ Failed to extract 2D region")
            return None

        # Get nucleus type for logging
        config = get_ps2d_config()
        #print(f"   📦 Extracted region: {region['intensity'].shape}")
        #print(f"      F1 ({config.nucleus_type}): {region['f1_ppm'].min():.3f} - {region['f1_ppm'].max():.3f} ppm")
        #print(f"      F2 (1H): {region['f2_ppm'].min():.3f} - {region['f2_ppm'].max():.3f} ppm")

        # Prepare initial parameters with data-driven linewidth estimates
        initial_peaks = []
        for peak_dict in overlap_group:
            # Extract peak position
            x_ppm = peak_dict.get('x_ppm') or peak_dict.get('pos_x')
            y_ppm = peak_dict.get('y_ppm') or peak_dict.get('pos_y')

            # Try to get detected intensity from peak picker
            detected_intensity = peak_dict.get('intensity')

            # Find nearest grid point (still needed for FWHM estimation)
            x_idx = np.argmin(np.abs(region['f2_ppm'] - x_ppm))
            y_idx = np.argmin(np.abs(region['f1_ppm'] - y_ppm))

            # ====================================================================
            # HEAVY OVERLAP AND TOO-CLOSE DETECTION (DECOUPLED)
            # ====================================================================
            # Two independent concepts:
            # 1. HEAVY OVERLAP: 1D cross-sections contaminated → use typical linewidths
            # 2. TOO CLOSE: Spectral ambiguity → apply L/G and intensity constraints
            # ====================================================================
            config = get_ps2d_config()
            heavily_overlapping = False  # Flag for linewidth estimation strategy
            tooclose = False             # Flag for L/G and intensity constraints

            # Check pairwise distances to all other peaks in cluster
            for other_peak in overlap_group:
                if other_peak is peak_dict:
                    continue  # Skip self-comparison

                other_x = other_peak.get('x_ppm') or other_peak.get('pos_x')
                other_y = other_peak.get('y_ppm') or other_peak.get('pos_y')

                # Calculate elliptical distance (normalized by fitting ellipse radii)
                dx = abs(x_ppm - other_x) / config.radF2
                dy = abs(y_ppm - other_y) / config.radF1
                elliptical_distance = np.sqrt(dx**2 + dy**2)

                # Check against both thresholds (independent criteria)
                if elliptical_distance < config.heavy_overlap_threshold:
                    heavily_overlapping = True

                if elliptical_distance < config.tooclose_threshold:
                    tooclose = True

                # Early exit if both flags are set
                if heavily_overlapping and tooclose:
                    break

            # ====================================================================
            # CONDITIONAL LINEWIDTH ESTIMATION
            # ====================================================================
            if heavily_overlapping:
                # Use spectrum-typical linewidths (not contaminated by overlap)
                #print(f"      ⚠️  Heavy overlap detected for peak at ({x_ppm:.3f}, {y_ppm:.3f})")
                #print(f"         Using typical linewidths instead of measured FWHM")
                fwhm_f1 = config.typical_linewidth_f1
                fwhm_f2 = config.typical_linewidth_f2
            else:
                # Estimate linewidths from 1D cross-sections FWHM (accurate for separated peaks)
                f1_cross = region['intensity'][:, x_idx]  # F1 slice at peak F2 position
                f2_cross = region['intensity'][y_idx, :]  # F2 slice at peak F1 position

                # Estimate FWHM by measuring width at half maximum
                fwhm_f1 = self._estimate_fwhm_from_slice(f1_cross, region['f1_ppm'], y_ppm)
                fwhm_f2 = self._estimate_fwhm_from_slice(f2_cross, region['f2_ppm'], x_ppm)

            # ====================================================================
            # CONVERT FWHM TO VOIGT COMPONENTS
            # ====================================================================
            # FIXED 2025-10-13: Correct FWHM storage convention
            # Bug: Was storing FWHM/2 in lw_gau, but rest of codebase expects full FWHM
            # The variable lw_gau_f1 should contain the Gaussian FWHM, not half of it
            # Use nucleus-adaptive minimum constraints from centralized config

            # NEW (CORRECT): Store full FWHM in lw_gau, small Lorentzian initial guess
            #lw_gau_f1 = max(fwhm_f1, config.min_linewidth_f1)  # lw_gau IS the Gaussian FWHM
            #lw_lor_f1 = max(fwhm_f1 / 10.0, config.min_linewidth_f1 / 5.0)  # Start with small Lorentzian
            #lw_gau_f2 = max(fwhm_f2, config.min_linewidth_f2)  # lw_gau IS the Gaussian FWHM
            #lw_lor_f2 = max(fwhm_f2 / 10.0, config.min_linewidth_f2 / 5.0)  # Start with small Lorentzian

            # OLD (BUG - COMMENTED OUT): Stored FWHM/2, causing confusion throughout codebase
            # Convert FWHM to Gaussian/Lorentzian components (assume 50/50 mix)
            # For Voigt: FWHM ≈ 0.5346*fL + sqrt(0.2166*fL² + fG²)
            # Approximate: set lw_gau ≈ fwhm/2, lw_lor ≈ fwhm/2
            lw_gau_f1 = max(fwhm_f1 / 2.0, config.min_linewidth_f1)  # BUG: Stores FWHM/2
            lw_lor_f1 = max(fwhm_f1 / 25.0, config.min_linewidth_f1) #/2#10
            lw_gau_f2 = max(fwhm_f2 / 2.0, config.min_linewidth_f2)  # BUG: Stores FWHM/2
            lw_lor_f2 = max(fwhm_f2 / 25.0, config.min_linewidth_f2) #10

            # Use detected intensity if available, otherwise fall back to re-measurement
            if detected_intensity is not None and detected_intensity > 0:
                # Use peak picker's measurement directly (no 0.8 scaling needed)
                initial_intensity = detected_intensity
            else:
                # Fallback: re-measure from grid (original behavior)
                local_intensity = region['intensity'][y_idx, x_idx]
                initial_intensity = max(local_intensity * 0.8, 100.0)

            initial_peaks.append({
                'pos_f1': y_ppm,  # NMRPipe: F1=Y=15N
                'pos_f2': x_ppm,  # NMRPipe: F2=X=1H
                'lw_lor_f1': lw_lor_f1, 'lw_gau_f1': lw_gau_f1,  # Data-driven
                'lw_lor_f2': lw_lor_f2, 'lw_gau_f2': lw_gau_f2,  # Data-driven
                'intensity': initial_intensity,
                'heavily_overlapping': heavily_overlapping,  # Flag for linewidth estimation strategy
                'tooclose': tooclose  # Flag for L/G and intensity constraints
            })

        # Log estimated linewidths for diagnostic purposes
        #print(f"   📏 Initial linewidth estimates (FWHM from 1D cross-sections):")
        for i, peak in enumerate(initial_peaks):
            # FIXED 2025-10-13: lw_gau IS the Gaussian FWHM (not half-width)
            #fwhm_f1 = peak['lw_gau_f1']  # NEW: lw_gau IS the FWHM (no compensation needed)
            #fwhm_f2 = peak['lw_gau_f2']  # NEW: lw_gau IS the FWHM (no compensation needed)
            fwhm_f1 = 2.0 * peak['lw_gau_f1']  # OLD: Compensated for FWHM/2 storage bug
            fwhm_f2 = 2.0 * peak['lw_gau_f2']  # OLD: Compensated for FWHM/2 storage bug
            #print(f"      Peak {i+1}: F1={fwhm_f1:.3f} ppm ({config.nucleus_type}), F2={fwhm_f2:.4f} ppm (1H)")

        # Log initial intensity estimates
        #print(f"   📊 Initial intensity estimates:")
        for i, (peak_dict, initial_peak) in enumerate(zip(overlap_group, initial_peaks)):
            detected = peak_dict.get('intensity')
            source = "peak picker" if (detected is not None and detected > 0) else "re-measured"
            #print(f"      Peak {i+1}: {initial_peak['intensity']:.2e} (source: {source})")

        # Estimate fitting time and warn user for large clusters
        n_peaks = len(overlap_group)
        n_points = region['intensity'].size
        # Empirical formula: ~40 seconds per 1000 points per peak² (calibrated from real data)
        estimated_time_sec = 40.0 * (n_points / 1000.0) * (n_peaks**2 / 81.0)

        if estimated_time_sec > 20:
            print(f"   ⏱️  Estimated fitting time: ~{estimated_time_sec:.0f} seconds ({estimated_time_sec/60:.1f} minutes)")
            print(f"      ({n_points:,} data points × {n_peaks} peaks)")

        #print(f"   🔧 Starting 5-stage Levenberg-Marquardt optimization...")

        # CRITICAL: Normalize data to [0,1] range for numerical stability
        # Without this, χ² values explode to 1e15 and damping fails
        max_intensity = np.max(np.abs(region['intensity']))
        if max_intensity == 0:
            max_intensity = 1.0

        normalized_data = region['intensity'] / max_intensity

        # Store original measured heights BEFORE volume estimation
        # These will be used by intensity ratio constraint in PS2D fitter
        original_heights = [peak['intensity'] for peak in initial_peaks]

        # Normalize initial intensity guesses
        # CRITICAL FIX: Convert HEIGHT to VOLUME using linewidth estimates
        # Root cause: Peak picker returns peak HEIGHT (point value at maximum),
        # but Voigt intensity parameter represents VOLUME (integral under 2D surface).
        # For 2D Gaussian: Volume = Height × 2π × σ_F1 × σ_F2
        # We use FWHM estimates from 1D cross-sections to calculate σ = FWHM / (2√(2ln2))
        SQRT_8LN2 = 2.3548200450309493  # 2 * sqrt(2 * ln(2))
        for peak in initial_peaks:
            # Convert FWHM to sigma (Gaussian width parameter)
            # FIXED 2025-10-13: lw_gau IS the Gaussian FWHM (not half-width)
            #fwhm_f1 = peak['lw_gau_f1']  # NEW: lw_gau IS the FWHM (no compensation needed)
            #fwhm_f2 = peak['lw_gau_f2']  # NEW: lw_gau IS the FWHM (no compensation needed)
            fwhm_f1 = 2.0 * peak['lw_gau_f1']  # OLD: Compensated for FWHM/2 storage bug
            fwhm_f2 = 2.0 * peak['lw_gau_f2']  # OLD: Compensated for FWHM/2 storage bug
            sigma_f1 = fwhm_f1 / SQRT_8LN2
            sigma_f2 = fwhm_f2 / SQRT_8LN2

            # Calculate volume from height: Volume = Height × 2π × σ_F1 × σ_F2
            # This is the Gaussian formula, but NMR peaks are Voigt (Gaussian + Lorentzian)
            # Empirically, Voigt volume is ~5× Gaussian volume due to Lorentzian tails
            height = peak['intensity']
            volume_gaussian = height * 2.0 * np.pi * sigma_f1 * sigma_f2
            volume_estimate = volume_gaussian * 5.0  # Voigt correction factor

            # Normalize
            peak['intensity'] = volume_estimate / max_intensity

            # Store original height for intensity ratio constraint
            peak['original_height'] = height

        #print(f"   📊 Data normalization: max = {max_intensity:.2e}")
        sys.stdout.flush()

        # Apply elliptical mask (PS2D approach: union of elliptical windows)
        # Convert overlap_group from dict to (f1, f2) format: (y_ppm, x_ppm)
        peak_positions_f1f2 = [
            (p.get('y_ppm') or p.get('pos_y'), p.get('x_ppm') or p.get('pos_x'))
            for p in overlap_group
        ]

        # Select data inside union of elliptical windows (spectrum.cpp:1010-1020)
        # Use nucleus-adaptive radii from centralized config
        config = get_ps2d_config()
        radF1 = config.radF1  # Nucleus-adaptive F1 ellipse radius
        radF2 = config.radF2  # Nucleus-adaptive F2 ellipse radius
        mask_result = select_data_2d_for_overlap_group(
            region['f1_grid'], region['f2_grid'], region['intensity'],
            peak_positions_f1f2, radF1=radF1, radF2=radF2
        )

        # Reshape mask back to 2D grid
        data_mask = mask_result['mask'].reshape(region['intensity'].shape)

        #print(f"   🎯 Elliptical masking: {mask_result['n_points_selected']}/{region['intensity'].size} points selected")
        #print(f"      (radF1={radF1} ppm, radF2={radF2} ppm)")
        sys.stdout.flush()

        # Fit using 2D multi-peak fitter
        # Constraint parameters controlled by ps2d_2d_fitter.py constructor defaults
        # To enable/disable constraints, change defaults in ps2d_2d_fitter.py:111
        # verbose=True with PROGRESS_PRINT_INTERVAL=200 (prints every 200 iterations)
        fitter = Ps2dMultiPeakFitter2D(verbose=True)
        result = fitter.fit_multi_peak_2d(
            region['f1_grid'], region['f2_grid'], normalized_data,
            initial_peaks,
            fix_positions=fix_positions,    # Pass GUI checkbox state
            fix_linewidths=fix_linewidths,  # Pass GUI checkbox state
            data_mask=data_mask
        )

        # Denormalize fitted intensities AND derived quantities
        # Height/volume/amplitude were calculated with normalized intensity,
        # so they must ALL be scaled back to original units
        if result['success']:
            for peak in result['peaks']:
                peak['intensity'] *= max_intensity
                peak['volume'] *= max_intensity      # Volume = intensity for normalized Voigt
                peak['height'] *= max_intensity      # Height scales linearly with intensity
                peak['amplitude'] *= max_intensity   # Amplitude = height

        if result['success']:
            print(f"   ✅ 2D fitting converged: R² = {result['r_squared']:.4f}")
            print(f"      Iterations: {result['iterations']}, χ² = {result['chi2']:.2e}")

            # Add assignments to each peak for GUI table display
            if peak_assignments and len(peak_assignments) == len(result['peaks']):
                for i, peak_fit in enumerate(result['peaks']):
                    peak_fit['assignment'] = peak_assignments[i]
                    print(f"      Peak {i+1} ({peak_assignments[i]}) fitted: ({peak_fit['pos_f2']:.4f}, {peak_fit['pos_f1']:.4f}) ppm")
            else:
                # Log fitted positions without assignments
                for i, peak_fit in enumerate(result['peaks']):
                    print(f"      Peak {i+1} fitted: ({peak_fit['pos_f2']:.4f}, {peak_fit['pos_f1']:.4f}) ppm")
        else:
            print(f"   ⚠️ 2D fitting did not converge")

        return result

    def fit_peak_voigt_2d(self, peak_x_ppm, peak_y_ppm, assignment="Unknown",
                          use_dynamic_optimization=True, all_peaks_context=None, linewidth_constraints=None): #True was false
        """
        2D simultaneous multi-peak Voigt fitting for all peaks

        FITTING STRATEGY:
        - All peaks are fitted using 2D simultaneous fitting (PS2D algorithm)
        - Isolated peaks are treated as single-peak groups
        - Overlapping peaks are fitted simultaneously for accurate deconvolution
        - No fallback to 1D cross-section fitting

        Parameters:
        - peak_x_ppm: X-axis peak position (ppm)
        - peak_y_ppm: Y-axis peak position (ppm)
        - assignment: Peak assignment/label
        - use_dynamic_optimization: Enable iterative optimization (default: True)
        - all_peaks_context: List of all peak positions [(x1,y1), (x2,y2), ...] for overlap detection
        - linewidth_constraints: Dict with 'x' and 'y' constraints for Global Optimization Manager

        Returns:
        - Dictionary with comprehensive fitting results or None if 2D fitting failed
        """
        print(f"Fitting Voigt profiles for peak {assignment} at ({peak_x_ppm:.3f}, {peak_y_ppm:.1f}) ppm")

        # ====================================================================
        # 2D SIMULTANEOUS FITTING: Always use 2D fitting for all peaks
        # ====================================================================
        # All peaks are fitted using 2D simultaneous multi-peak Voigt fitting.
        # Isolated peaks are treated as single-peak groups.
        if PS2D_2D_FITTER_AVAILABLE and all_peaks_context is not None:
            needs_2d, overlap_group = self.check_if_peaks_need_2d_fitting(
                peak_x_ppm, peak_y_ppm, all_peaks_context
            )

            # Always use 2D fitting (both isolated and overlapping peaks)
            #if len(overlap_group) > 1:
            #    print(f"   🔍 Overlap group: {len(overlap_group)} peaks within elliptical window")
            #else:
            #    print(f"   🔍 Isolated peak: single-peak 2D fitting")
            #print(f"   🎯 Routing to 2D simultaneous fitting...")

            # Extract assignments from all_peaks_context by matching positions
            peak_assignments = []
            for peak_dict in overlap_group:
                x_pos = peak_dict.get('x_ppm') or peak_dict.get('pos_x')
                y_pos = peak_dict.get('y_ppm') or peak_dict.get('pos_y')
                matched_assignment = None
                for peak_ctx in all_peaks_context:
                    ctx_x = peak_ctx.get('x_ppm') or peak_ctx.get('pos_x')
                    ctx_y = peak_ctx.get('y_ppm') or peak_ctx.get('pos_y')
                    # Match by position (within 0.001 ppm tolerance)
                    if abs(ctx_x - x_pos) < 0.001 and abs(ctx_y - y_pos) < 0.01:
                        matched_assignment = peak_ctx.get('assignment', 'Unknown')
                        break
                peak_assignments.append(matched_assignment or 'Unknown')

            # Extract fix_positions and fix_linewidths from GUI parameters
            fix_positions = self.gui_params.get('fix_positions', False)
            fix_linewidths = self.gui_params.get('fix_linewidths', False)

            # Fit entire overlap group using 2D multi-peak fitter
            result_2d = self.fit_overlap_group_2d(overlap_group, assignment,
                                                   peak_assignments=peak_assignments,
                                                   fix_positions=fix_positions,
                                                   fix_linewidths=fix_linewidths)

            if result_2d and result_2d['success']:
                # Find the result for the target peak
                # Match by finding closest position
                best_match_idx = 0
                min_dist = float('inf')
                for i, peak_fit in enumerate(result_2d['peaks']):
                    dist = np.sqrt((peak_fit['pos_f2'] - peak_x_ppm)**2 +
                                  (peak_fit['pos_f1'] - peak_y_ppm)**2)
                    if dist < min_dist:
                        min_dist = dist
                        best_match_idx = i

                target_peak = result_2d['peaks'][best_match_idx]

                # Extract detected intensity from all_peaks_context by matching position
                detected_intensity = None
                if all_peaks_context is not None:
                    for peak_ctx in all_peaks_context:
                        ctx_x = peak_ctx.get('x_ppm') or peak_ctx.get('pos_x')
                        ctx_y = peak_ctx.get('y_ppm') or peak_ctx.get('pos_y')
                        # Match by position (within 0.001 ppm tolerance)
                        if abs(ctx_x - peak_x_ppm) < 0.001 and abs(ctx_y - peak_y_ppm) < 0.01:
                            detected_intensity = peak_ctx.get('intensity')
                            break

                # Extract 2D region data for proper 2D contour visualization
                region_2d = self.extract_2d_region_for_overlap_group(overlap_group)

                # Reconstruct 2D fitted surface from PS2D parameters
                fitted_2d_surface, individual_surfaces, baseline = self._reconstruct_2d_surface(
                    region_2d, result_2d['peaks']
                )

                # Convert to expected format with 2D visualization data
                return {
                    'assignment': assignment,
                    'peak_position': (target_peak['pos_f2'], target_peak['pos_f1']),
                    'x_fit': {
                        'center': target_peak['pos_f2'],
                        'sigma': target_peak['lw_gau_f2'],
                        'gamma': target_peak['lw_lor_f2'],
                        'amplitude': target_peak['amplitude'],
                        'r_squared': result_2d['r_squared'],
                        'success': True,
                        'method': '2d_simultaneous'
                    },
                    'y_fit': {
                        'center': target_peak['pos_f1'],
                        'sigma': target_peak['lw_gau_f1'],
                        'gamma': target_peak['lw_lor_f1'],
                        'amplitude': target_peak['amplitude'],
                        'r_squared': result_2d['r_squared'],
                        'success': True,
                        'method': '2d_simultaneous'
                    },
                    'fitting_quality': self._assign_quality_from_r2(result_2d['r_squared']),
                    'quality': self._assign_quality_from_r2(result_2d['r_squared']),
                    'avg_r_squared': result_2d['r_squared'],
                    'method': '2d_simultaneous_multi_peak',
                    'overlap_group_size': len(overlap_group),
                    # Peak quantitation (from 2D fit)
                    'volume': target_peak['volume'],
                    'height': target_peak['height'],
                    'amplitude': target_peak['amplitude'],
                    'detected_intensity': detected_intensity,  # Intensity from peak detection (before fitting)
                    # 2D visualization data
                    'region_2d': region_2d,
                    'fitted_2d_surface': fitted_2d_surface,
                    'individual_surfaces': individual_surfaces,  # Individual peaks for multi-color visualization
                    'all_peaks': result_2d['peaks'],
                    'baseline': baseline  # Baseline offset for visualization
                }
            else:
                print(f"   ❌ 2D fitting failed - peak not fitted")
                return None

        # === DYNAMIC OPTIMIZATION OR STANDARD FITTING ===
        # CRITICAL: Disable dynamic optimization when PS2D multi-peak is enabled
        # because enhanced_fitter doesn't support Y cross-section re-extraction
        gui_params = getattr(self, 'gui_params', None)
        use_ps2d = gui_params.get('use_ps2d_multi_peak', False) if gui_params else False

        # Initialize window_optimization to avoid reference errors
        window_optimization = None

        if use_ps2d and use_dynamic_optimization:
            print(f"   ℹ️ PS2D multi-peak enabled: using standard fitting (bypassing dynamic optimization)")
            use_dynamic_optimization = False

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
                    max_iterations=50,  # Enhanced: Up to 50 tests for comprehensive exploration
                    r2_improvement_threshold=0.01  # More sensitive to improvements
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

            # DEBUG: Show what gui_params we're passing
            if gui_params:
                print(f"   📋 fit_peak_voigt_2d: gui_params exists, use_ps2d_multi_peak={gui_params.get('use_ps2d_multi_peak', 'N/A')}")
            else:
                print(f"   ⚠️ fit_peak_voigt_2d: gui_params is None!")

            # Fit X-dimension with adaptive strategy (single vs multi-peak)
            print(f"   🔍 Analyzing X-dimension for {assignment}")
            print(f"   [X-FIT] Target X position: {peak_x_ppm:.4f} ppm")
            print(f"   [X-FIT] X cross-section: {len(regions['x_cross_section'])} points from {regions['x_ppm_scale'][0]:.4f} to {regions['x_ppm_scale'][-1]:.4f} ppm")
            x_fit = self.adaptive_fit_1d(regions['x_ppm_scale'], regions['x_cross_section'],
                                        peak_x_ppm, dimension='x', gui_params=gui_params)

            # CRITICAL FIX: Re-extract Y cross-section at the FITTED X position
            # Problem: X fit may shift peak position, so Y cross-section must be updated
            # to ensure it passes through the correct peak in 2D space
            x_params_temp = self._safe_extract_params(x_fit, default_center=peak_x_ppm)
            fitted_x_ppm = x_params_temp.get('center', peak_x_ppm)

            # Check if X position shifted significantly
            x_shift = abs(fitted_x_ppm - peak_x_ppm)
            if x_shift > 0.001:  # Threshold: 0.001 ppm
                print(f"   🔄 X fit shifted peak: {peak_x_ppm:.4f} → {fitted_x_ppm:.4f} ppm (Δ={x_shift:.4f})")
                print(f"   📍 Re-extracting Y cross-section at updated X position")

                # Find new X index for the fitted position
                x_idx_fitted = np.argmin(np.abs(self.ppm_x_axis - fitted_x_ppm))

                # Re-extract Y cross-section at the fitted X position
                if self.ps2d_data_selector is not None:
                    # Use PS2D method for consistency
                    y_cross_full = self.nmr_data[:, x_idx_fitted]
                    y_ppm_full = self.ppm_y_axis

                    print(f"   [Y-REEXTRACT] Extracting Y at fitted X={fitted_x_ppm:.4f} ppm (index {x_idx_fitted})")
                    print(f"   [Y-REEXTRACT] Y window centered at {peak_y_ppm:.4f} ppm with radius 0.6 ppm")

                    # Apply PS2D's elliptical window selection
                    y_selection = self.ps2d_data_selector.select_data_elliptical(
                        y_ppm_full, y_cross_full, fitted_x_ppm, peak_y_ppm, dimension='y'
                    )

                    # Update regions with new Y cross-section
                    regions['y_cross_section'] = y_selection['y_selected']
                    regions['y_ppm_scale'] = y_selection['x_selected']

                    # Show Y range after re-extraction
                    y_min, y_max = y_selection['x_selected'][0], y_selection['x_selected'][-1]
                    print(f"   ✅ Y cross-section updated: {y_selection['n_points_selected']} points")
                    print(f"   [Y-REEXTRACT] Y range: {y_min:.4f} to {y_max:.4f} ppm (expected: {peak_y_ppm-0.6:.4f} to {peak_y_ppm+0.6:.4f})")
                else:
                    # Fallback: simple extraction at fitted X position
                    print(f"   ⚠️ PS2D not available, using simple Y re-extraction")
                    regions['y_cross_section'] = self.nmr_data[:, x_idx_fitted]
            else:
                print(f"   ✓ X position stable (shift={x_shift:.4f} ppm < 0.001 threshold)")

            # Fit Y-dimension with adaptive strategy (single vs multi-peak)
            # Now using the CORRECTED Y cross-section that passes through the fitted X position
            print(f"   🔍 Analyzing Y-dimension for {assignment}")
            print(f"   [Y-FIT] Target Y position: {peak_y_ppm:.4f} ppm")
            print(f"   [Y-FIT] Y cross-section: {len(regions['y_cross_section'])} points from {regions['y_ppm_scale'][0]:.4f} to {regions['y_ppm_scale'][-1]:.4f} ppm")
            y_fit = self.adaptive_fit_1d(regions['y_ppm_scale'], regions['y_cross_section'],
                                        peak_y_ppm, dimension='y', gui_params=gui_params)

            # ITERATIVE REFINEMENT: Re-extract X cross-section at fitted Y position and refine
            y_params_temp = self._safe_extract_params(y_fit, default_center=peak_y_ppm)
            fitted_y_ppm = y_params_temp.get('center', peak_y_ppm)

            print(f"   [Y-FIT] Fitted Y position: {fitted_y_ppm:.4f} ppm (shift from target: {fitted_y_ppm - peak_y_ppm:.4f} ppm)")
            if y_fit and 'multi_peak_info' in y_fit:
                n_peaks = y_fit['multi_peak_info'].get('n_peaks', 0)
                print(f"   [Y-FIT] Multi-peak fit: {n_peaks} peaks detected")
                all_peaks = y_fit['multi_peak_info'].get('all_peaks', [])
                for i, pk in enumerate(all_peaks):
                    dist = abs(pk['position'] - peak_y_ppm)
                    print(f"   [Y-FIT]    Peak {i+1}: pos={pk['position']:.4f} ppm (dist from target: {dist:.4f} ppm)")

            y_shift = abs(fitted_y_ppm - peak_y_ppm)
            if y_shift > 0.01:  # Threshold: 0.01 ppm for Y dimension (10x larger due to 15N scale)
                print(f"   🔄 Y fit shifted peak: {peak_y_ppm:.4f} → {fitted_y_ppm:.4f} ppm (Δ={y_shift:.4f})")
                print(f"   📍 Re-extracting X cross-section at updated Y position for refinement")

                # Find new Y index for the fitted position
                y_idx_fitted = np.argmin(np.abs(self.ppm_y_axis - fitted_y_ppm))

                # Re-extract X cross-section at the fitted Y position
                if self.ps2d_data_selector is not None:
                    x_cross_full = self.nmr_data[y_idx_fitted, :]
                    x_ppm_full = self.ppm_x_axis

                    # Apply PS2D's elliptical window selection
                    x_selection = self.ps2d_data_selector.select_data_elliptical(
                        x_ppm_full, x_cross_full, fitted_x_ppm, fitted_y_ppm, dimension='x'
                    )

                    # Update regions with new X cross-section
                    regions['x_cross_section'] = x_selection['y_selected']
                    regions['x_ppm_scale'] = x_selection['x_selected']
                    print(f"   ✅ X cross-section updated: {x_selection['n_points_selected']} points")

                    # Refine X fit with updated cross-section
                    print(f"   🔄 Refining X-dimension fit with updated cross-section")
                    x_fit = self.adaptive_fit_1d(regions['x_ppm_scale'], regions['x_cross_section'],
                                                fitted_x_ppm, dimension='x', gui_params=gui_params)

                    # Check if refinement changed X position again
                    x_params_refined = self._safe_extract_params(x_fit, default_center=fitted_x_ppm)
                    refined_x_ppm = x_params_refined.get('center', fitted_x_ppm)
                    x_refinement_shift = abs(refined_x_ppm - fitted_x_ppm)

                    if x_refinement_shift > 0.001:
                        print(f"   🔄 X refinement shifted: {fitted_x_ppm:.4f} → {refined_x_ppm:.4f} ppm (Δ={x_refinement_shift:.4f})")
                    else:
                        print(f"   ✓ X position converged (shift={x_refinement_shift:.4f} ppm)")
                else:
                    print(f"   ⚠️ PS2D not available, skipping iterative refinement")
            else:
                print(f"   ✓ Y position stable (shift={y_shift:.4f} ppm < 0.01 threshold)")

        # Safely extract parameters from both fits (handles any format)
        x_params = self._safe_extract_params(x_fit, default_center=peak_x_ppm)
        y_params = self._safe_extract_params(y_fit, default_center=peak_y_ppm)

        # Log final fitted positions vs initial (diagnostic for PS2D integration fix)
        final_x = x_params.get('center', peak_x_ppm)
        final_y = y_params.get('center', peak_y_ppm)
        total_x_shift = abs(final_x - peak_x_ppm)
        total_y_shift = abs(final_y - peak_y_ppm)

        print(f"   📊 Final fitted position: ({final_x:.4f}, {final_y:.4f}) ppm")
        print(f"   📏 Total shift from initial: ΔX={total_x_shift:.4f} ppm, ΔY={total_y_shift:.4f} ppm")

        # Check fitting quality
        if not (x_params.get('success', False) and y_params.get('success', False)):
            print(f"Voigt fitting failed for {assignment}")
            return None

        # Determine overall fitting quality using LOCAL metrics (consistent with enhanced fitter)
        # Use local R-squared if available, otherwise fall back to global
        x_r_squared = x_params.get('quality_metrics', {}).get('r_squared_local', x_params['r_squared'])
        y_r_squared = y_params.get('quality_metrics', {}).get('r_squared_local', y_params['r_squared'])
        avg_r_squared = (x_r_squared + y_r_squared) / 2

        # Assign quality using centralized method (updated thresholds 2025-10-14)
        quality = self._assign_quality_from_r2(avg_r_squared)

        # Extract detected intensity from all_peaks_context by matching position
        detected_intensity = None
        if all_peaks_context is not None:
            for peak_ctx in all_peaks_context:
                ctx_x = peak_ctx.get('x_ppm') or peak_ctx.get('pos_x')
                ctx_y = peak_ctx.get('y_ppm') or peak_ctx.get('pos_y')
                # Match by position (within 0.001 ppm tolerance)
                if abs(ctx_x - peak_x_ppm) < 0.001 and abs(ctx_y - peak_y_ppm) < 0.01:
                    detected_intensity = peak_ctx.get('intensity')
                    break

        # Store detailed results with both local and global quality metrics
        result = {
            'assignment': assignment,
            'peak_position': (peak_x_ppm, peak_y_ppm),
            'x_fit': {
                'ppm_scale': regions['x_ppm_scale'],
                'cross_section': regions['x_cross_section'],
                'fitted_curve': x_params['fitted_curve'],
                'amplitude': x_params['amplitude'],
                'center': x_params['center'],
                'sigma': x_params['sigma'],
                'gamma': x_params['gamma'],
                'baseline': x_params['baseline'],
                'r_squared': x_params.get('quality_metrics', {}).get('r_squared_local', x_params['r_squared']),
                # Include quality metrics for local/global assessment
                'quality_metrics': x_params.get('quality_metrics', {}),
                'r_squared_local': x_params.get('quality_metrics', {}).get('r_squared_local', x_params['r_squared']),
                'r_squared_global': x_params.get('quality_metrics', {}).get('r_squared_global', x_params['r_squared']),
                # Add window size metadata for display
                'window_size': abs(regions['x_ppm_scale'][-1] - regions['x_ppm_scale'][0]) if len(regions['x_ppm_scale']) > 1 else 'unknown',
                'gui_based': True,
                'gui_window_param': self.fitting_parameters.get('fitting_window_x', 'unknown')
            },
            'y_fit': {
                'ppm_scale': regions['y_ppm_scale'],
                'cross_section': regions['y_cross_section'],
                'fitted_curve': y_params['fitted_curve'],
                'amplitude': y_params['amplitude'],
                'center': y_params['center'],
                'sigma': y_params['sigma'],
                'gamma': y_params['gamma'],
                'baseline': y_params['baseline'],
                'r_squared': y_params.get('quality_metrics', {}).get('r_squared_local', y_params['r_squared']),
                # Include quality metrics for local/global assessment
                'quality_metrics': y_params.get('quality_metrics', {}),
                'r_squared_local': y_params.get('quality_metrics', {}).get('r_squared_local', y_params['r_squared']),
                'r_squared_global': y_params.get('quality_metrics', {}).get('r_squared_global', y_params['r_squared']),
                # Add window size metadata for display
                'window_size': abs(regions['y_ppm_scale'][-1] - regions['y_ppm_scale'][0]) if len(regions['y_ppm_scale']) > 1 else 'unknown',
                'gui_based': True,
                'gui_window_param': self.fitting_parameters.get('fitting_window_y', 'unknown')
            },
            'fitting_quality': quality,
            'avg_r_squared': avg_r_squared,
            'avg_r_squared_local': (x_r_squared + y_r_squared) / 2,
            'avg_r_squared_global': (x_params.get('quality_metrics', {}).get('r_squared_global', x_params['r_squared']) +
                                   y_params.get('quality_metrics', {}).get('r_squared_global', y_params['r_squared'])) / 2,
            'detected_intensity': detected_intensity,  # Intensity from peak detection (before fitting)
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

        # Perform baseline fit to establish reference quality (with PS2D fix)
        # Get gui_params from self if available
        gui_params = getattr(self, 'gui_params', None)

        # Fit X-dimension
        baseline_x_fit = self.adaptive_fit_1d(baseline_regions['x_ppm_scale'],
                                             baseline_regions['x_cross_section'],
                                             peak_x_ppm, dimension='x', gui_params=gui_params)

        if not baseline_x_fit.get('success'):
            return {'success': False, 'reason': 'baseline_x_fitting_failed'}

        # Re-extract Y cross-section at fitted X position if needed
        x_params_temp = self._safe_extract_params(baseline_x_fit, default_center=peak_x_ppm)
        fitted_x_ppm = x_params_temp.get('center', peak_x_ppm)
        x_shift = abs(fitted_x_ppm - peak_x_ppm)

        if x_shift > 0.001 and self.ps2d_data_selector is not None:
            print(f"   [BASELINE] X shifted {peak_x_ppm:.4f}→{fitted_x_ppm:.4f}, re-extracting Y")
            x_idx_fitted = np.argmin(np.abs(self.ppm_x_axis - fitted_x_ppm))
            y_cross_full = self.nmr_data[:, x_idx_fitted]
            y_ppm_full = self.ppm_y_axis

            y_selection = self.ps2d_data_selector.select_data_elliptical(
                y_ppm_full, y_cross_full, fitted_x_ppm, peak_y_ppm, dimension='y'
            )

            baseline_regions['y_cross_section'] = y_selection['y_selected']
            baseline_regions['y_ppm_scale'] = y_selection['x_selected']
            print(f"   [BASELINE] Y re-extracted: {y_selection['n_points_selected']} pts")

        # Fit Y-dimension with corrected cross-section
        baseline_y_fit = self.adaptive_fit_1d(baseline_regions['y_ppm_scale'],
                                             baseline_regions['y_cross_section'],
                                             peak_y_ppm, dimension='y', gui_params=gui_params)

        if not baseline_y_fit.get('success'):
            return {'success': False, 'reason': 'baseline_y_fitting_failed'}

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

        # Phase 4: Enhanced comprehensive optimization strategy
        print(f"   Peak interference: {interference_analysis['isolation_level']}")

        # Use enhanced bidirectional optimization for all cases
        optimization_result = self._enhanced_window_optimization(
            peak_x_ppm, peak_y_ppm, assignment,
            base_x_window, base_y_window,
            baseline_avg_r2, max_iterations, r2_improvement_threshold
        )

        # Add interference analysis to results
        optimization_result['interference_analysis'] = interference_analysis

        return optimization_result

    def _enhanced_window_optimization(self, peak_x_ppm, peak_y_ppm, assignment,
                                     base_x_window, base_y_window, baseline_r2,
                                     max_iterations=50, r2_threshold=0.01):
        """
        ENHANCED WINDOW OPTIMIZATION: Bidirectional, Independent X/Y, R²-driven

        Implements comprehensive window optimization with:
        1. Bidirectional exploration (expansion + contraction)
        2. Independent X and Y window optimization
        3. R²-driven adaptive search algorithm
        4. Increment-based exploration (0.1 step size)
        5. Direction reversal when R² degrades
        6. Best-value tracking and restart mechanism

        Args:
            peak_x_ppm, peak_y_ppm: Peak position
            assignment: Peak identifier
            base_x_window, base_y_window: Starting window sizes (user values)
            baseline_r2: Reference R² for comparison
            max_iterations: Maximum optimization tests (default 50)
            r2_threshold: Minimum R² improvement threshold

        Returns:
            dict: Comprehensive optimization results
        """
        import numpy as np

        print(f"   🚀 Enhanced window optimization for {assignment}")
        print(f"   Starting: X={base_x_window:.3f}, Y={base_y_window:.3f}, R²={baseline_r2:.3f}")

        # Search parameters
        step_size = 0.1
        x_bounds = (0.1, 1.0)  # Reasonable X window bounds
        y_bounds = (0.5, 5.0)  # Reasonable Y window bounds
        degradation_threshold = 0.05  # Stop if R² drops by this much

        # Optimization state
        current_x = base_x_window
        current_y = base_y_window
        best_x = base_x_window
        best_y = base_y_window
        best_r2 = baseline_r2

        iteration_history = []
        tests_performed = 0

        # Phase 1: Initial exploration around starting point
        print(f"   📊 Phase 1: Initial exploration (±{step_size})")

        initial_tests = [
            # Test X variations (Y constant)
            (current_x + step_size, current_y, "X+0.1"),
            (current_x - step_size, current_y, "X-0.1"),
            # Test Y variations (X constant)
            (current_x, current_y + step_size, "Y+0.1"),
            (current_x, current_y - step_size, "Y-0.1"),
            # Test diagonal combinations
            (current_x + step_size, current_y + step_size, "X+0.1,Y+0.1"),
            (current_x - step_size, current_y - step_size, "X-0.1,Y-0.1")
        ]

        for test_x, test_y, direction in initial_tests:
            if tests_performed >= max_iterations:
                break

            # Bounds checking
            if not (x_bounds[0] <= test_x <= x_bounds[1] and y_bounds[0] <= test_y <= y_bounds[1]):
                continue

            test_r2 = self._evaluate_window_quality(peak_x_ppm, peak_y_ppm, test_x, test_y)
            tests_performed += 1

            iteration_history.append({
                'test': tests_performed,
                'x_window': test_x,
                'y_window': test_y,
                'r2': test_r2,
                'direction': direction,
                'phase': 'initial_exploration'
            })

            print(f"     Test {tests_performed}: {direction} → X={test_x:.3f}, Y={test_y:.3f}, R²={test_r2:.3f}")

            if test_r2 > best_r2 + r2_threshold:
                best_x, best_y, best_r2 = test_x, test_y, test_r2
                print(f"       ✓ New best: R²={best_r2:.3f} (Δ={test_r2-baseline_r2:.3f})")

        # Phase 2: Adaptive directional search from best point
        print(f"   🎯 Phase 2: Adaptive search from best point (X={best_x:.3f}, Y={best_y:.3f})")

        current_x, current_y = best_x, best_y

        # Determine promising directions from initial exploration
        x_improvements = [h for h in iteration_history if 'X+' in h['direction'] and h['r2'] > baseline_r2]
        y_improvements = [h for h in iteration_history if 'Y+' in h['direction'] and h['r2'] > baseline_r2]

        x_direction = 1 if x_improvements else -1  # Prefer expansion if it helped
        y_direction = 1 if y_improvements else -1

        # Separate X and Y optimization
        for dimension, direction, current_val, bounds_tuple in [
            ('X', x_direction, current_x, x_bounds),
            ('Y', y_direction, current_y, y_bounds)
        ]:
            if tests_performed >= max_iterations:
                break

            print(f"     Optimizing {dimension} dimension (direction: {'+' if direction > 0 else '-'})")

            consecutive_degradations = 0
            search_val = current_val

            while consecutive_degradations < 3 and tests_performed < max_iterations:
                search_val += direction * step_size

                # Bounds check
                if not (bounds_tuple[0] <= search_val <= bounds_tuple[1]):
                    break

                # Test with other dimension at current best
                test_x = search_val if dimension == 'X' else best_x
                test_y = search_val if dimension == 'Y' else best_y

                test_r2 = self._evaluate_window_quality(peak_x_ppm, peak_y_ppm, test_x, test_y)
                tests_performed += 1

                iteration_history.append({
                    'test': tests_performed,
                    'x_window': test_x,
                    'y_window': test_y,
                    'r2': test_r2,
                    'direction': f"{dimension}{'+' if direction > 0 else '-'}{step_size}",
                    'phase': f'{dimension}_optimization'
                })

                trend_emoji = '✅' if test_r2 > best_r2 + r2_threshold else '⚠️' if test_r2 < best_r2 - degradation_threshold else '➡️'
                print(f"       🧪 Test {tests_performed}: {dimension}={search_val:.3f} → R²={test_r2:.3f} {trend_emoji}")

                if test_r2 > best_r2 + r2_threshold:
                    best_x, best_y, best_r2 = test_x, test_y, test_r2
                    consecutive_degradations = 0
                    print(f"         ✅ New optimum: R²={best_r2:.3f} (Δ={test_r2-baseline_r2:.3f})")
                elif test_r2 < best_r2 - degradation_threshold:
                    consecutive_degradations += 1
                    print(f"         ⚠️ Degradation {consecutive_degradations}/3: R²={test_r2:.3f}")
                else:
                    consecutive_degradations += 1
                    print(f"         ➡️ Minimal change: R²={test_r2:.3f}")

            # If search failed in this direction, try opposite direction
            if consecutive_degradations >= 3 and tests_performed < max_iterations:
                print(f"     🔄 Reversing {dimension} direction")
                direction *= -1
                search_val = current_val if dimension == 'X' else current_y
                consecutive_degradations = 0

                for _ in range(5):  # Limited reverse exploration
                    if tests_performed >= max_iterations:
                        break

                    search_val += direction * step_size

                    if not (bounds_tuple[0] <= search_val <= bounds_tuple[1]):
                        break

                    test_x = search_val if dimension == 'X' else best_x
                    test_y = search_val if dimension == 'Y' else best_y

                    test_r2 = self._evaluate_window_quality(peak_x_ppm, peak_y_ppm, test_x, test_y)
                    tests_performed += 1

                    iteration_history.append({
                        'test': tests_performed,
                        'x_window': test_x,
                        'y_window': test_y,
                        'r2': test_r2,
                        'direction': f"{dimension}_reverse{'+' if direction > 0 else '-'}{step_size}",
                        'phase': f'{dimension}_reverse'
                    })

                    print(f"         🔄 Reverse test: {dimension}={search_val:.3f} → R²={test_r2:.3f}")

                    if test_r2 > best_r2 + r2_threshold:
                        best_x, best_y, best_r2 = test_x, test_y, test_r2
                        print(f"         ✅ Reverse breakthrough: R²={best_r2:.3f}")

        # Phase 3: Final convergence around best point
        if tests_performed < max_iterations and best_r2 > baseline_r2 + r2_threshold:
            print(f"   🎯 Phase 3: Fine-tuning around optimal point")

            fine_tune_tests = [
                (best_x + step_size/2, best_y, "fine_X+"),
                (best_x - step_size/2, best_y, "fine_X-"),
                (best_x, best_y + step_size/2, "fine_Y+"),
                (best_x, best_y - step_size/2, "fine_Y-")
            ]

            for test_x, test_y, direction in fine_tune_tests:
                if tests_performed >= max_iterations:
                    break
                if not (x_bounds[0] <= test_x <= x_bounds[1] and y_bounds[0] <= test_y <= y_bounds[1]):
                    continue

                test_r2 = self._evaluate_window_quality(peak_x_ppm, peak_y_ppm, test_x, test_y)
                tests_performed += 1

                print(f"       🎯 Fine-tune: {direction}, R²={test_r2:.3f}")

                if test_r2 > best_r2:
                    best_x, best_y, best_r2 = test_x, test_y, test_r2
                    print(f"         ✅ Fine-tune improved: R²={best_r2:.3f}")

        # Calculate final results
        total_improvement = best_r2 - baseline_r2

        print(f"   🎯 Enhanced optimization complete:")
        print(f"     Final: X={best_x:.3f}, Y={best_y:.3f}, R²={best_r2:.3f}")
        print(f"     Improvement: Δ={total_improvement:.3f} ({tests_performed} tests)")

        return {
            'success': True,
            'optimization_type': 'enhanced_bidirectional',
            'optimized_x_window': best_x,
            'optimized_y_window': best_y,
            'optimized_r2': best_r2,
            'baseline_r2': baseline_r2,
            'improvement': total_improvement,
            'tests_performed': tests_performed,
            'iteration_history': iteration_history,
            'recommendation': 'enhanced_optimal' if total_improvement > r2_threshold else 'keep_gui',
            'search_statistics': {
                'x_range_explored': (min(h['x_window'] for h in iteration_history),
                                   max(h['x_window'] for h in iteration_history)),
                'y_range_explored': (min(h['y_window'] for h in iteration_history),
                                   max(h['y_window'] for h in iteration_history)),
                'best_improvement_test': max(iteration_history, key=lambda x: x['r2'])['test']
            }
        }

    def _evaluate_window_quality(self, peak_x_ppm, peak_y_ppm, x_window, y_window):
        """
        Evaluate R² quality for given window sizes.

        Args:
            peak_x_ppm, peak_y_ppm: Peak position
            x_window, y_window: Window sizes to test

        Returns:
            float: Average R² value for this window configuration
        """
        try:
            # Extract regions with test window sizes
            test_regions = self.extract_peak_region(peak_x_ppm, peak_y_ppm, x_window, y_window)
            if not test_regions:
                return 0.0

            # Perform fits with cross-section re-extraction (PS2D integration fix)
            # Get gui_params from self if available
            gui_params = getattr(self, 'gui_params', None)

            # Fit X-dimension
            x_fit = self.adaptive_fit_1d(test_regions['x_ppm_scale'],
                                       test_regions['x_cross_section'],
                                       peak_x_ppm, dimension='x', gui_params=gui_params)

            if not x_fit.get('success'):
                return 0.0

            # Re-extract Y cross-section at fitted X position if needed
            x_params_temp = self._safe_extract_params(x_fit, default_center=peak_x_ppm)
            fitted_x_ppm = x_params_temp.get('center', peak_x_ppm)
            x_shift = abs(fitted_x_ppm - peak_x_ppm)

            if x_shift > 0.001 and self.ps2d_data_selector is not None:
                # Re-extract Y cross-section at fitted X position
                print(f"     [WINDOW_EVAL] X shifted {peak_x_ppm:.4f}→{fitted_x_ppm:.4f}, re-extracting Y")
                x_idx_fitted = np.argmin(np.abs(self.ppm_x_axis - fitted_x_ppm))
                y_cross_full = self.nmr_data[:, x_idx_fitted]
                y_ppm_full = self.ppm_y_axis

                y_selection = self.ps2d_data_selector.select_data_elliptical(
                    y_ppm_full, y_cross_full, fitted_x_ppm, peak_y_ppm, dimension='y'
                )

                test_regions['y_cross_section'] = y_selection['y_selected']
                test_regions['y_ppm_scale'] = y_selection['x_selected']
                print(f"     [WINDOW_EVAL] Y re-extracted: {y_selection['n_points_selected']} pts")

            # Fit Y-dimension with corrected cross-section
            y_fit = self.adaptive_fit_1d(test_regions['y_ppm_scale'],
                                       test_regions['y_cross_section'],
                                       peak_y_ppm, dimension='y', gui_params=gui_params)

            if not y_fit.get('success'):
                return 0.0

            # Extract R² values
            x_r2 = x_fit.get('quality_metrics', {}).get('r_squared_local', x_fit.get('r_squared', 0))
            y_r2 = y_fit.get('quality_metrics', {}).get('r_squared_local', y_fit.get('r_squared', 0))

            return (x_r2 + y_r2) / 2

        except Exception as e:
            print(f"     ⚠ Window evaluation failed: {e}")
            return 0.0

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

        print(f"   🔍 Strategy: Expanding windows for isolated peak {assignment}")
        print(f"   📊 Baseline: X={base_x_window:.2f}, Y={base_y_window:.2f}, R²={baseline_r2:.3f}")

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
            # Get gui_params from self if available
            gui_params = getattr(self, 'gui_params', None)
            x_fit = self.adaptive_fit_1d(test_regions['x_ppm_scale'],
                                       test_regions['x_cross_section'],
                                       peak_x_ppm, dimension='x', gui_params=gui_params)
            y_fit = self.adaptive_fit_1d(test_regions['y_ppm_scale'],
                                       test_regions['y_cross_section'],
                                       peak_y_ppm, dimension='y', gui_params=gui_params)

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

            print(f"   🧪 Test {iteration+1}: Factor={factor:.1f}x, Windows=({test_x_window:.2f}, {test_y_window:.2f}), "
                  f"R²={avg_r2:.3f}, Points=({len(test_regions['x_cross_section'])}, {len(test_regions['y_cross_section'])})")

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
        print(f"   🔍 Strategy: Contracting windows for crowded peak {assignment} ({interferer_count} interferers)")
        print(f"   📊 Baseline: X={base_x_window:.2f}, Y={base_y_window:.2f}, R²={baseline_r2:.3f}")

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
            # Get gui_params from self if available
            gui_params = getattr(self, 'gui_params', None)
            x_fit = self.adaptive_fit_1d(test_regions['x_ppm_scale'],
                                       test_regions['x_cross_section'],
                                       peak_x_ppm, dimension='x', gui_params=gui_params)
            y_fit = self.adaptive_fit_1d(test_regions['y_ppm_scale'],
                                       test_regions['y_cross_section'],
                                       peak_y_ppm, dimension='y', gui_params=gui_params)

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

            print(f"   🧪 Test {iteration+1}: Factor={factor:.1f}x, Windows=({test_x_window:.2f}, {test_y_window:.2f}), "
                  f"R²={avg_r2:.3f}, Points=({len(test_regions['x_cross_section'])}, {len(test_regions['y_cross_section'])})")

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
        Extract 1D cross-sections around a peak for fitting using PS2D's FIXED elliptical windows.

        - FIXED radius: radF1=0.6 ppm (15N/13C), radF2=0.06 ppm (1H)
        - Elliptical selection: radius² = (F1-pos)²/radF1² + (F2-pos)²/radF2² <= 1.0
        - NO adaptive windows, NO multipliers, NO linewidth-based scaling

        Returns:
        - x_cross_section: 1D array along x-axis (1H dimension)
        - y_cross_section: 1D array along y-axis (15N/13C dimension)
        - x_ppm_scale: corresponding ppm values
        - y_ppm_scale: corresponding ppm values
        """
        if self.nmr_data is None:
            return None

        # PS2D_SRC METHOD: Use FIXED elliptical windows (EXACT clone)
        if self.ps2d_data_selector is not None:
            # Find peak position in data points
            x_idx = np.argmin(np.abs(self.ppm_x_axis - peak_x_ppm))
            y_idx = np.argmin(np.abs(self.ppm_y_axis - peak_y_ppm))

            # Extract FULL horizontal slice (along X/1H dimension at peak Y position)
            x_cross_full = self.nmr_data[y_idx, :]
            x_ppm_full = self.ppm_x_axis

            # Extract FULL vertical slice (along Y/15N dimension at peak X position)
            y_cross_full = self.nmr_data[:, x_idx]
            y_ppm_full = self.ppm_y_axis

            # Apply PS2D's FIXED elliptical window selection (spectrum.cpp lines 1010-1020)
            x_selection = self.ps2d_data_selector.select_data_elliptical(
                x_ppm_full, x_cross_full, peak_x_ppm, peak_y_ppm, dimension='x'
            )
            y_selection = self.ps2d_data_selector.select_data_elliptical(
                y_ppm_full, y_cross_full, peak_x_ppm, peak_y_ppm, dimension='y'
            )

            # Log FIXED window usage (for debugging)
            print(f"   📏 PS2D FIXED windows: X={x_selection['radius_used']:.3f} ppm, Y={y_selection['radius_used']:.1f} ppm")
            print(f"   📊 Data points selected: X={x_selection['n_points_selected']}, Y={y_selection['n_points_selected']}")

            return {
                'x_cross_section': x_selection['y_selected'],
                'y_cross_section': y_selection['y_selected'],
                'x_ppm_scale': x_selection['x_selected'],
                'y_ppm_scale': y_selection['x_selected'],
                'peak_indices': (x_idx, y_idx),
                'extraction_bounds': None,  # Not used in PS2D method
                'ps2d_method': True,  # Flag for logging
                'radF1_used': y_selection['radius_used'],
                'radF2_used': x_selection['radius_used']
            }

        # FALLBACK: OLD ADAPTIVE METHOD (for backward compatibility if PS2D not available)
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

        print(f"   ⚠️ Using FALLBACK adaptive windows (PS2D not available)")

        return {
            'x_cross_section': x_cross_section,
            'y_cross_section': y_cross_section,
            'x_ppm_scale': x_ppm_scale,
            'y_ppm_scale': y_ppm_scale,
            'peak_indices': (x_idx, y_idx),
            'extraction_bounds': (x_min, x_max, y_min, y_max),
            'ps2d_method': False
        }

    @staticmethod
    def calculate_r_squared(y_true, y_pred):
        """Calculate R-squared value"""
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

    def detect_peaks_professional(self, **kwargs):
        """Professional peak detection - enhanced version with reference option"""
        # If we have a peak list loaded, use reference-based detection
        if hasattr(self, 'peak_list') and self.peak_list is not None:
            return self._detect_peaks_reference_based()
        else:
            # No peak list - use fallback
            print("❌ No peak list loaded for standard detection mode")
            print("   Please load a peak list or switch to S/N threshold mode")
            return []

    def _detect_peaks_reference_based(self):
        """Detect all peaks then match closest to peak list (exact backup approach)"""
        if self.peak_list is None:
            raise ValueError("Peak list must be loaded before detection")
        if self.nmr_data is None:
            raise ValueError("NMR data must be loaded before detection")

        print(f"📋 Detection: detect all peaks then match closest to {len(self.peak_list)} reference peaks")

        # Step 1: Detect ALL peaks using threshold (like S/N native mode)
        # Ensure noise level is estimated
        if not hasattr(self, 'noise_level') or self.noise_level is None:
            self._estimate_noise_level_advanced()

        # Use GUI threshold parameters (1H/15N ppm values are search windows)
        threshold_multiplier = 3.0  # Default threshold multiplier
        signal_threshold = self.noise_level * threshold_multiplier

        print(f"   Step 1: Detecting all peaks (threshold={signal_threshold:.2e})")

        # Detect all peaks using the existing threshold method
        all_detected_peaks = self._detect_peaks_by_threshold(signal_threshold)

        print(f"   Found {len(all_detected_peaks)} peaks total")

        # Step 2: Match each reference peak to closest detected peak
        search_window_x_ppm = self.search_window_x  # From GUI: 1H/15N (ppm) first value
        search_window_y_ppm = self.search_window_y  # From GUI: 1H/15N (ppm) second value

        #print(f"   Step 2: Matching to reference peaks (search windows: X=±{search_window_x_ppm:.3f}, Y=±{search_window_y_ppm:.3f} ppm)")

        matched_peaks = []
        used_peaks = set()  # Track which detected peaks have been used

        for idx, peak_row in self.peak_list.iterrows():
            try:
                ref_x_ppm = float(peak_row['Position_X'])
                ref_y_ppm = float(peak_row['Position_Y'])
                assignment = peak_row.get('Assignment', f'Peak_{idx+1}')

                # Find closest detected peak within search window
                best_match = None
                best_distance = float('inf')

                for i, detected_peak in enumerate(all_detected_peaks):
                    if i in used_peaks:  # Skip already used peaks
                        continue

                    # Calculate distance from reference
                    x_distance = abs(detected_peak['x_ppm'] - ref_x_ppm)
                    y_distance = abs(detected_peak['y_ppm'] - ref_y_ppm)

                    # Check if within search window
                    if x_distance <= search_window_x_ppm and y_distance <= search_window_y_ppm:
                        total_distance = np.sqrt(x_distance**2 + y_distance**2)
                        if total_distance < best_distance:
                            best_distance = total_distance
                            best_match = (i, detected_peak)

                if best_match is not None:
                    # Use the matched detected peak
                    peak_idx, detected_peak = best_match
                    used_peaks.add(peak_idx)

                    matched_peak = {
                        'assignment': assignment,
                        'ppm_x': detected_peak['x_ppm'],
                        'ppm_y': detected_peak['y_ppm'],
                        'x_point': detected_peak['indices'][1],
                        'y_point': detected_peak['indices'][0],
                        'intensity': detected_peak['intensity'],
                        'snr': detected_peak['snr'],
                        'detected': True,
                        'detection_quality': 'Matched',
                        'distance_from_reference': best_distance
                    }

                    #if best_distance > 0.001:
                    #    print(f"   {assignment}: matched to peak {best_distance:.4f} ppm away")

                else:
                    # No match found, use reference position
                    ref_x_idx = np.argmin(np.abs(self.ppm_x_axis - ref_x_ppm))
                    ref_y_idx = np.argmin(np.abs(self.ppm_y_axis - ref_y_ppm))
                    ref_intensity = self.nmr_data[ref_y_idx, ref_x_idx]
                    ref_snr = np.abs(ref_intensity) / self.noise_level if self.noise_level > 0 else 0

                    matched_peak = {
                        'assignment': assignment,
                        'ppm_x': ref_x_ppm,
                        'ppm_y': ref_y_ppm,
                        'x_point': ref_x_idx,
                        'y_point': ref_y_idx,
                        'intensity': ref_intensity,
                        'snr': ref_snr,
                        'detected': False,
                        'reference_retained': True,
                        'detection_quality': 'Reference'
                    }

                matched_peaks.append(matched_peak)

            except (ValueError, KeyError) as e:
                print(f"⚠️ Skipping peak {idx+1}: {e}")
                continue

        # Store results
        self.fitted_peaks = matched_peaks

        # Update statistics
        detected_count = sum(1 for p in matched_peaks if p.get('detected', False))
        reference_retained = sum(1 for p in matched_peaks if p.get('reference_retained', False))

        self.detection_statistics = {
            'total_peaks': len(self.peak_list),
            'detected_peaks': detected_count,
            'reference_retained': reference_retained,
            'detection_rate': (detected_count / len(self.peak_list) * 100) if len(self.peak_list) > 0 else 0
        }

        print(f"✅ Detection complete: {detected_count}/{len(self.peak_list)} matched, {reference_retained} references retained")

        return matched_peaks

    def detect_peaks_full_mode(self, **kwargs):
        """Full peak detection mode"""
        return self.detect_peaks_professional(**kwargs)

    def detect_peaks_inplace_mode(self, **kwargs):
        """In-place fitting mode with enhanced capabilities"""
        # Enhanced detection with better reference handling
        return self.detect_peaks_professional(**kwargs)

    def load_nmr_file(self, nmr_file, skip_nucleus_detection=False):
        """
        Load NMR data file (convenience method for GUI)

        Parameters:
        -----------
        nmr_file : str
            Path to NMR data file
        skip_nucleus_detection : bool, optional
            If True, skip auto-detection of nucleus type (default: False)
            Used in series integration to avoid redundant detection for spectra 2+
        """
        # Always load NMR data - spectrum loading should work regardless of peak list
        return self._load_nmr_data_only(nmr_file, skip_nucleus_detection=skip_nucleus_detection)

    def _load_nmr_data_only(self, nmr_file, skip_nucleus_detection=False):
        """
        Load only NMR data without peak list

        Parameters:
        -----------
        nmr_file : str
            Path to NMR data file
        skip_nucleus_detection : bool, optional
            If True, skip auto-detection of nucleus type (default: False)
            Used in series integration to avoid redundant detection for spectra 2+
        """
        try:
            import nmrglue as ng
            self.nmr_dict, self.nmr_data = ng.pipe.read(nmr_file)
            print(f"Loaded NMR data: {self.nmr_data.shape} from {nmr_file}")

            # Calculate PPM axes
            self._calculate_ppm_axes()

            # Auto-detect nucleus type based on spectral dimensions (unless skipped)
            if not skip_nucleus_detection:
                detected_nucleus = self._detect_nucleus_type()
                if detected_nucleus:
                    self.auto_detected_nucleus = detected_nucleus
                    print(f"🔬 Auto-detected nucleus type: {detected_nucleus}-HSQC")
                    print(f"   Spectral range: X={self.ppm_x_axis[0]:.2f}-{self.ppm_x_axis[-1]:.2f} ppm, "
                          f"Y={self.ppm_y_axis[0]:.1f}-{self.ppm_y_axis[-1]:.1f} ppm")
                else:
                    self.auto_detected_nucleus = None
                    print(f"⚠️  Could not auto-detect nucleus type from spectral dimensions")
                    print(f"   Spectral range: X={self.ppm_x_axis[0]:.2f}-{self.ppm_x_axis[-1]:.2f} ppm, "
                          f"Y={self.ppm_y_axis[0]:.1f}-{self.ppm_y_axis[-1]:.1f} ppm")
                    print(f"   Please select nucleus type manually in Expert Mode (PS2D Algorithm Configuration)")
            else:
                print(f"⏭️  Skipping nucleus auto-detection (using reference spectrum configuration)")

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

    def load_spectrum_only(self, nmr_file):
        """Load spectrum without peak list requirement"""
        try:
            success = self._load_nmr_data_only(nmr_file)
            if success:
                self.peak_list = None  # No peak list in S/N mode
                print(f"✅ Spectrum-only loading successful: {nmr_file}")
                print(f"   Data shape: {self.nmr_data.shape}")
                print(f"   X-axis range: {self.ppm_x_axis[0]:.2f} to {self.ppm_x_axis[-1]:.2f} ppm")
                print(f"   Y-axis range: {self.ppm_y_axis[0]:.1f} to {self.ppm_y_axis[-1]:.1f} ppm")
            return success
        except Exception as e:
            print(f"Error in spectrum-only loading: {e}")
            return False

    def process_peaks(self, **kwargs):
        """Process peaks based on current mode"""
        if self.processing_mode == 'full_detection':
            return self.detect_peaks_full_mode(**kwargs)
        elif self.processing_mode == 'in_place':
            return self.detect_peaks_inplace_mode(**kwargs)
        elif self.processing_mode == 'sn_native':
            return self.detect_peaks_sn_native(**kwargs)
        else:
            raise ValueError(f"Unknown processing mode: {self.processing_mode}")

    def detect_peaks_sn_native(self, **kwargs):
        """Native S/N threshold-based peak detection without peak list"""
        print(f"🎯 Starting S/N native detection (threshold={self.sn_threshold}, expected={self.expected_peak_count})")

        if self.nmr_data is None:
            print("❌ No NMR data loaded")
            return []

        # Step 1: Estimate noise level
        noise_level = self._estimate_noise_level_advanced()
        signal_threshold = noise_level * self.sn_threshold

        print(f"   Noise level: {noise_level:.2e}")
        print(f"   S/N threshold: {self.sn_threshold}")
        print(f"   Signal threshold: {signal_threshold:.2e}")

        # Step 2: Detect peaks using threshold
        detected_peaks = self._detect_peaks_by_threshold(signal_threshold)

        # Step 3: Apply expected count cutoff
        if len(detected_peaks) > self.expected_peak_count:
            # Sort by intensity and keep top N
            detected_peaks.sort(key=lambda p: p['intensity'], reverse=True)
            detected_peaks = detected_peaks[:self.expected_peak_count]
            print(f"   Applied count cutoff: {self.expected_peak_count} peaks (from {len(detected_peaks)} detected)")

        print(f"✅ S/N native detection complete: {len(detected_peaks)} peaks")

        # Step 4: Convert to standard format for compatibility
        standardized_peaks = []
        for i, peak in enumerate(detected_peaks):
            peak_data = {
                'assignment': f'Peak_{i+1:03d}',
                'ppm_x': peak['x_ppm'],
                'ppm_y': peak['y_ppm'],
                'intensity': peak['intensity'],
                'snr': peak['snr'],
                'detection_method': 'sn_threshold',
                'detected': True,  # Mark as detected for GUI statistics
                'fitted': False    # No fitting performed in S/N mode
            }
            standardized_peaks.append(peak_data)

        # Step 5: Populate fitted_peaks for GUI compatibility
        self.fitted_peaks = standardized_peaks

        # Step 6: Update detection statistics for GUI status display
        detected_count = len(standardized_peaks)
        self.detection_statistics = {
            'total_peaks': detected_count,  # For S/N detection, all found peaks are the total
            'detected_peaks': detected_count,
            'detection_rate': 100.0,  # 100% since we're finding peaks directly
            'reference_retained': 0,  # No reference peaks in S/N mode
            'method': 'sn_native'
        }
        print(f"   Updated detection statistics: {detected_count} peaks (100% detection rate)")

        return standardized_peaks

    def get_detection_statistics(self):
        """Get detailed detection statistics"""
        return self.detection_statistics.copy()

    def _estimate_noise_level_advanced(self):
        """Advanced noise level estimation for S/N detection"""
        if self.sn_detection_params['noise_estimation'] == 'corners':
            # Use corner regions (more reliable for 2D NMR)
            h, w = self.nmr_data.shape
            corner_size = min(20, h//20, w//20)

            corners = [
                self.nmr_data[:corner_size, :corner_size],
                self.nmr_data[:corner_size, -corner_size:],
                self.nmr_data[-corner_size:, :corner_size],
                self.nmr_data[-corner_size:, -corner_size:]
            ]

            noise_data = np.concatenate([corner.flatten() for corner in corners])
            noise_level = np.std(noise_data)

        else:  # median-based estimation
            # Use median absolute deviation (more robust to outliers)
            flattened = self.nmr_data.flatten()
            median = np.median(flattened)
            mad = np.median(np.abs(flattened - median))
            noise_level = mad * 1.4826  # Convert MAD to std equivalent

        return noise_level

    def _calculate_top_contour_center(self, y_idx, x_idx, intensity_band=0.05,
                                       max_shift_x=0.04, max_shift_y=0.14):
        """
        Calculate sub-pixel peak center using parabolic interpolation.

        Fits a parabola to 3 points around the pixel maximum in both X and Y
        dimensions to find the analytical peak center with sub-pixel precision.
        This is the standard method for peak localization in image processing.

        Algorithm:
        1. Extract 3 points in X direction: I(x-1), I(x), I(x+1)
        2. Fit parabola: I(x) = a*x^2 + b*x + c
        3. Find maximum: x_max = -b/(2*a) = (I_left - I_right) / (2*(I_left - 2*I_center + I_right))
        4. Repeat for Y direction
        5. Apply safety constraints to limit maximum shift

        Parameters:
        -----------
        y_idx, x_idx : int
            Pixel indices of detected maximum
        intensity_band : float
            Not used in parabolic interpolation (kept for API compatibility)
        max_shift_x, max_shift_y : float
            Maximum allowed shift from original position in ppm

        Returns:
        --------
        center_x, center_y : float
            Sub-pixel peak center position in ppm
        """
        # Get original position
        original_x = self.ppm_x_axis[x_idx]
        original_y = self.ppm_y_axis[y_idx]
        peak_max = self.nmr_data[y_idx, x_idx]

        # Define search window based on max_shift constraints
        x_ppm_per_pixel = abs(self.ppm_x_axis[1] - self.ppm_x_axis[0]) if len(self.ppm_x_axis) > 1 else 0.004
        y_ppm_per_pixel = abs(self.ppm_y_axis[1] - self.ppm_y_axis[0]) if len(self.ppm_y_axis) > 1 else 0.08

        search_radius_x = int(np.ceil(max_shift_x / x_ppm_per_pixel))
        search_radius_y = int(np.ceil(max_shift_y / y_ppm_per_pixel))

        x_start = max(0, x_idx - search_radius_x)
        x_end = min(self.nmr_data.shape[1], x_idx + search_radius_x + 1)
        y_start = max(0, y_idx - search_radius_y)
        y_end = min(self.nmr_data.shape[0], y_idx + search_radius_y + 1)

        # Extract search region
        search_region = self.nmr_data[y_start:y_end, x_start:x_end]

        # Define intensity range for "top contour"
        intensity_min = peak_max * (1.0 - intensity_band)
        intensity_max = peak_max * (1.0 + intensity_band)

        # Find all pixels within intensity band
        mask = (search_region >= intensity_min) & (search_region <= intensity_max)

        if not np.any(mask):
            # Fallback to original position
            return original_x, original_y

        # Get coordinates and intensities of top contour pixels
        local_y_indices, local_x_indices = np.where(mask)
        global_y_indices = local_y_indices + y_start
        global_x_indices = local_x_indices + x_start

        intensities = self.nmr_data[global_y_indices, global_x_indices]

        # Convert to ppm coordinates
        ppm_x_vals = self.ppm_x_axis[global_x_indices]
        ppm_y_vals = self.ppm_y_axis[global_y_indices]

        # Parabolic interpolation to find sub-pixel peak maximum
        # Uses 3-point parabola fit in both X and Y dimensions
        # This is the standard method for sub-pixel peak localization

        # Check if we have enough pixels for 3-point interpolation
        if (y_idx > 0 and y_idx < self.nmr_data.shape[0] - 1 and
            x_idx > 0 and x_idx < self.nmr_data.shape[1] - 1):

            # Extract 3 points in X direction (horizontal)
            I_x_left = self.nmr_data[y_idx, x_idx - 1]
            I_x_center = self.nmr_data[y_idx, x_idx]
            I_x_right = self.nmr_data[y_idx, x_idx + 1]

            # Parabolic interpolation formula in X
            # Fits parabola I(x) = a*x^2 + b*x + c to 3 points
            # Maximum at x = -b/(2*a) = (I_left - I_right) / (2*(I_left - 2*I_center + I_right))
            denom_x = 2.0 * (I_x_left - 2.0 * I_x_center + I_x_right)
            if abs(denom_x) > 1e-10:  # Avoid division by zero
                delta_x = (I_x_left - I_x_right) / denom_x  # Fractional pixel shift
                # Clip to reasonable range (±1 pixel)
                delta_x = np.clip(delta_x, -1.0, 1.0)
                x_center_idx = x_idx + delta_x
            else:
                # Parabola is flat (all 3 points equal), use pixel center
                x_center_idx = float(x_idx)

            # Extract 3 points in Y direction (vertical)
            I_y_bottom = self.nmr_data[y_idx - 1, x_idx]
            I_y_center = self.nmr_data[y_idx, x_idx]
            I_y_top = self.nmr_data[y_idx + 1, x_idx]

            # Parabolic interpolation formula in Y
            denom_y = 2.0 * (I_y_bottom - 2.0 * I_y_center + I_y_top)
            if abs(denom_y) > 1e-10:  # Avoid division by zero
                delta_y = (I_y_bottom - I_y_top) / denom_y  # Fractional pixel shift
                # Clip to reasonable range (±1 pixel)
                delta_y = np.clip(delta_y, -1.0, 1.0)
                y_center_idx = y_idx + delta_y
            else:
                # Parabola is flat, use pixel center
                y_center_idx = float(y_idx)

            # Convert interpolated pixel indices to ppm coordinates
            center_x = np.interp(x_center_idx, np.arange(len(self.ppm_x_axis)), self.ppm_x_axis)
            center_y = np.interp(y_center_idx, np.arange(len(self.ppm_y_axis)), self.ppm_y_axis)

        else:
            # Edge case: peak too close to spectrum boundary for 3-point interpolation
            # Use pixel maximum (no sub-pixel refinement)
            center_x = original_x
            center_y = original_y

        # Apply safety constraint: limit maximum shift
        shift_x = center_x - original_x
        shift_y = center_y - original_y

        shift_x = np.clip(shift_x, -max_shift_x, max_shift_x)
        shift_y = np.clip(shift_y, -max_shift_y, max_shift_y)

        final_x = original_x + shift_x
        final_y = original_y + shift_y

        return final_x, final_y

    def _detect_peaks_by_threshold(self, signal_threshold):
        """Detect peaks using signal threshold"""
        from scipy.ndimage import maximum_filter, label

        # Apply threshold to identify peak regions
        peak_mask = self.nmr_data > signal_threshold

        # Use local maxima to find individual peaks
        local_maxima = maximum_filter(self.nmr_data, size=3) == self.nmr_data
        peak_candidates = peak_mask & local_maxima

        # Find connected components (individual peaks)
        labeled_peaks, num_peaks = label(peak_candidates)

        detected_peaks = []

        for peak_id in range(1, num_peaks + 1):
            # Get peak region
            peak_region = labeled_peaks == peak_id
            peak_indices = np.where(peak_region)

            if len(peak_indices[0]) == 0:
                continue

            # Find peak maximum within region
            region_intensities = self.nmr_data[peak_indices]
            max_idx = np.argmax(region_intensities)

            y_idx = peak_indices[0][max_idx]
            x_idx = peak_indices[1][max_idx]

            # Apply top contour centroid (always enabled)
            gui_params = getattr(self, 'gui_params', {})

            # Get pixel position first
            pixel_x_ppm = self.ppm_x_axis[x_idx]
            pixel_y_ppm = self.ppm_y_axis[y_idx]

            # Calculate top contour centroid using GUI parameters
            max_shift_x = gui_params.get('centroid_window_x_ppm', 0.01)
            max_shift_y = gui_params.get('centroid_window_y_ppm', 0.06)

            centroid_x_ppm, centroid_y_ppm = self._calculate_top_contour_center(
                y_idx, x_idx,
                intensity_band=0.05,  # ±5% around maximum
                max_shift_x=max_shift_x,
                max_shift_y=max_shift_y
            )

            # Calculate shift distance
            shift_x = centroid_x_ppm - pixel_x_ppm
            shift_y = centroid_y_ppm - pixel_y_ppm
            shift_distance = np.sqrt(shift_x**2 + shift_y**2)

            # Only print if significant shift (>0.001 ppm)
            #if shift_distance > 0.001:
            #    print(f"   Centroid shift: Δ={shift_distance:.4f} ppm from pixel max")

            x_ppm = centroid_x_ppm
            y_ppm = centroid_y_ppm

            intensity = self.nmr_data[y_idx, x_idx]

            # Calculate S/N ratio
            noise_level = self._estimate_noise_level_advanced()
            snr = intensity / noise_level if noise_level > 0 else 0

            # Apply minimum separation filter
            min_separation = self.sn_detection_params['peak_separation']
            too_close = False

            for existing_peak in detected_peaks:
                x_dist = abs(existing_peak['x_ppm'] - x_ppm)
                y_dist = abs(existing_peak['y_ppm'] - y_ppm)

                if x_dist < min_separation and y_dist < min_separation * 10:  # Y has broader tolerance
                    too_close = True
                    break

            if not too_close and snr >= self.sn_detection_params['min_snr']:
                detected_peaks.append({
                    'x_ppm': x_ppm,
                    'y_ppm': y_ppm,
                    'intensity': intensity,
                    'snr': snr,
                    'indices': (x_idx, y_idx)
                })

        return detected_peaks

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

        if PS2D_DATA_SELECTOR_AVAILABLE:
            self.ps2d_data_selector = Ps2dDataSelector(spectrum_type='15N-HSQC')
        else:
            self.ps2d_data_selector = None

        # Detection statistics for GUI status display
        self.detection_statistics = {
            'total_peaks': 0,
            'detected_peaks': 0,
            'detection_rate': 0.0,
            'reference_retained': 0,
            'method': 'unknown'
        }

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
        Set integration mode for detection-fitting workflow

        Args:
            mode: Only 'standard' is supported (integrated/adaptive modes removed)
            **integration_params: integration-specific parameters
        """
        if mode not in ['standard', 'integrated', 'adaptive']:
            raise ValueError("Mode must be 'standard', 'integrated', or 'adaptive'")

        if mode != 'standard':
            print("⚠️ Integrated/adaptive modes are no longer available, using standard mode")
            mode = 'standard'

        self.integration_mode = mode
        self.integration_parameters.update(integration_params)
        print("📊 Standard detection-fitting mode")

        return self.integration_parameters.copy()

    def get_integration_status(self):
        """Get current integration mode and parameters"""
        return {
            'mode': self.integration_mode,
            'parameters': self.integration_parameters.copy(),
            'integrated_fitter_available': False,
            'capabilities': {
                'adaptive_thresholds': False,
                'multi_resolution_detection': False,
                'quality_scoring': False,
                'chemical_shift_context': False
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

    def _assign_quality_from_r2(self, r_squared):
        """
        Assign quality label based on R² value

        Updated thresholds (2025-10-14):
        - Excellent: R² >= 0.9
        - Good: 0.8 <= R² < 0.9
        - Fair: 0.5 <= R² < 0.8
        - Poor: 0.2 <= R² < 0.5
        - Failed: R² < 0.2
        """
        if r_squared >= 0.9:
            return "Excellent"
        elif r_squared >= 0.8:
            return "Good"
        elif r_squared >= 0.5:
            return "Fair"
        elif r_squared >= 0.2:
            return "Poor"
        else:
            return "Failed"

    def enhanced_peak_fitting(self, peak_x_ppm, peak_y_ppm, assignment="Unknown", linewidth_constraints=None, all_peaks_context=None, use_dynamic_optimization=True):
        """
        Enhanced peak fitting with advanced features and consensus fitting integration

        Parameters
        ----------
        peak_x_ppm : float
            Peak position in X dimension (ppm)
        peak_y_ppm : float
            Peak position in Y dimension (ppm)
        assignment : str, optional
            Peak assignment/label
        linewidth_constraints : dict, optional
            PS2D linewidth reuse constraints from reference spectrum
            Format: {'x': {'sigma_bounds': (min, max), 'gamma_bounds': (min, max)},
                     'y': {'sigma_bounds': (min, max), 'gamma_bounds': (min, max)}}
        all_peaks_context : list of dict, optional
            List of all peaks in the spectrum for overlap detection and 2D routing
            Format: [{'assignment': str, 'x_ppm': float, 'y_ppm': float, 'pos_x': float, 'pos_y': float}, ...]
        use_dynamic_optimization : bool, optional
            Enable iterative optimization (default: True)
        """

        # Check if simplified parameters are enabled (indicating automated mode preference)
        use_consensus = False
        if hasattr(self, 'parameter_manager') and hasattr(self.parameter_manager, 'use_simplified_mode'):
            use_consensus = self.parameter_manager.use_simplified_mode
        elif hasattr(self, 'gui_params') and self.gui_params.get('use_simplified_parameters', False):
            use_consensus = True

        if use_consensus and hasattr(self.voigt_fitter, 'fit_with_consensus'):
            # Use advanced consensus fitting for automated mode
            print(f"🚀 Using consensus fitting for peak at ({peak_x_ppm:.4f}, {peak_y_ppm:.4f})")

            # Extract region data for consensus fitting using existing method
            region_data = self.extract_peak_region(peak_x_ppm, peak_y_ppm, self.fitting_parameters.get('fitting_window_x', 0.15), self.fitting_parameters.get('fitting_window_y', 2.0))
            if region_data is not None:
                # Use the x cross-section for 1D consensus fitting
                x_region = region_data['x_ppm_scale']
                y_region = region_data['x_cross_section']  # Intensity data along x-axis
            else:
                x_region, y_region = None, None
            if x_region is not None and y_region is not None:
                nucleus_type = self.detect_nucleus_type(peak_x_ppm, peak_y_ppm)
                consensus_result = self.voigt_fitter.fit_with_consensus(
                    x_region, y_region, nucleus_type, use_simplified=True
                )

                if consensus_result and consensus_result.get('success', False):
                    # Convert consensus result to standard format
                    result = self.convert_consensus_result_to_standard(consensus_result, peak_x_ppm, peak_y_ppm, assignment)
                    print(f"✅ Consensus fitting successful: R² = {result.get('avg_r_squared', 0):.4f}")

                    self.update_statistics(True, result['avg_r_squared'])
                    return result
                else:
                    print(f"⚠️ Consensus fitting failed, falling back to legacy method")

        # Fallback to legacy fitting method with PS2D linewidth constraints and 2D routing support
        result = self.fit_peak_voigt_2d(
            peak_x_ppm, peak_y_ppm, assignment,
            use_dynamic_optimization=use_dynamic_optimization,
            linewidth_constraints=linewidth_constraints,  # Pass through PS2D constraints
            all_peaks_context=all_peaks_context  # Enable automatic 2D routing for overlapping peaks
        )

        if result:
            self.update_statistics(True, result['avg_r_squared'])
            # Convert nested format to flat format for consistency with consensus results
            result = self.convert_voigt2d_result_to_standard(result)
        else:
            self.update_statistics(False)

        return result


    def detect_nucleus_type(self, peak_x_ppm, peak_y_ppm):
        """Detect nucleus type based on peak position"""
        # Simple heuristic based on chemical shift ranges
        if abs(peak_x_ppm) < 20:  # 1H range approximately 0-15 ppm
            return '1H'
        elif 100 <= abs(peak_x_ppm) <= 200:  # 13C range approximately 0-200 ppm
            return '13C'
        elif abs(peak_y_ppm) > 100:  # 15N range approximately 0-600 ppm
            return '15N'
        else:
            return '1H'  # Default fallback

    def convert_consensus_result_to_standard(self, consensus_result, peak_x_ppm, peak_y_ppm, assignment):
        """Convert consensus fitting result to standard peak result format"""
        try:
            # Extract key metrics from consensus result
            r_squared = consensus_result.get('r_squared', 0)
            parameters = consensus_result.get('parameters', [])
            fitted_curve = consensus_result.get('fitted_curve')

            # Extract fitted parameters (assuming Voigt format: [amp, center, sigma, gamma, baseline])
            if len(parameters) >= 5:
                amplitude = parameters[0]
                center = parameters[1]
                sigma = parameters[2]
                gamma = parameters[3]
                baseline = parameters[4]
                width = sigma + gamma  # Total width
            else:
                # Fallback values
                amplitude, center, sigma, gamma, baseline, width = 1000, peak_x_ppm, 0.01, 0.01, 0, 0.02

            # Create standard result dictionary
            result = {
                'assignment': assignment,
                'peak_x': peak_x_ppm,
                'peak_y': peak_y_ppm,
                'amplitude': amplitude,
                'volume': amplitude * width * np.pi,  # Approximate volume
                'avg_r_squared': r_squared,
                'height': amplitude,
                'width': width,
                'sigma': sigma,
                'gamma': gamma,
                'baseline': baseline,
                'center_x': center,
                'center_y': peak_y_ppm,
                'success': True,
                'method': 'consensus_fitting',
                'quality_class': consensus_result.get('quality_class', 'Fair'),
                'fitted_curve': fitted_curve,
                'original_consensus_result': consensus_result
            }

            return result

        except Exception as e:
            print(f"Failed to convert consensus result: {e}")
            # Return minimal fallback result
            return {
                'assignment': assignment,
                'peak_x': peak_x_ppm,
                'peak_y': peak_y_ppm,
                'amplitude': 1000,
                'volume': 1000,
                'avg_r_squared': 0.5,
                'success': False,
                'method': 'consensus_fitting_failed'
            }

    def convert_voigt2d_result_to_standard(self, voigt2d_result):
        """Convert nested 2D Voigt result to flat standard format

        This ensures enhanced_peak_fitting() always returns a consistent flat format,
        regardless of whether it uses consensus fitting or fit_peak_voigt_2d().

        Args:
            voigt2d_result: Result dict from fit_peak_voigt_2d() with nested x_fit/y_fit structure

        Returns:
            Flattened result dict compatible with consensus format and parallel processing
        """
        try:
            # Extract x and y fit parameters (use x_fit for primary parameters)
            x_fit = voigt2d_result.get('x_fit', {})
            y_fit = voigt2d_result.get('y_fit', {})

            # Extract key parameters from x_fit and y_fit
            amplitude = x_fit.get('amplitude', 1000)
            center_x = x_fit.get('center', 0)
            center_y = y_fit.get('center', 0)
            sigma_x = x_fit.get('sigma', 0.01)
            gamma_x = x_fit.get('gamma', 0.01)
            sigma_y = y_fit.get('sigma', 0.01)
            gamma_y = y_fit.get('gamma', 0.01)
            baseline = x_fit.get('baseline', 0)
            width = sigma_x + gamma_x

            # Get peak positions from original result
            peak_position = voigt2d_result.get('peak_position', (center_x, center_y))

            # CRITICAL FIX: Calculate proper height using PS2D formula
            # For single peak 1D fits, calculate height from fitted Voigt parameters
            # Height = maximum value at peak center, not amplitude
            from lunaNMR.core.ps2d_2d_fitter import calculate_peak_height

            # Calculate proper height (maximum value at peak center)
            # Note: calculate_peak_height expects (lw_lor_f1, lw_gau_f1, lw_lor_f2, lw_gau_f2, intensity)
            height = calculate_peak_height(gamma_y, sigma_y, gamma_x, sigma_x, amplitude)

            # Create flat result dictionary
            avg_r2 = voigt2d_result.get('avg_r_squared', 0)
            result = {
                'assignment': voigt2d_result.get('assignment', 'Unknown'),
                'peak_x': peak_position[0],
                'peak_y': peak_position[1],
                'amplitude': amplitude,
                'volume': amplitude * width * np.pi,  # Approximate volume
                'avg_r_squared': avg_r2,
                'r_squared': avg_r2,  # Alias for GUI compatibility
                'height': height,  # FIXED: Use calculated height instead of amplitude
                'width': width,
                'sigma': sigma_x,
                'gamma': gamma_x,
                'baseline': baseline,
                'center_x': center_x,
                'center_y': center_y,
                'success': True,
                'method': voigt2d_result.get('method', 'voigt_2d_fitting'),
                'fitting_quality': voigt2d_result.get('fitting_quality', 'Unknown'),
                # Preserve original nested structure for advanced analysis
                'x_fit': x_fit,
                'y_fit': y_fit,
                # Preserve additional metadata
                'avg_r_squared_local': voigt2d_result.get('avg_r_squared_local'),
                'avg_r_squared_global': voigt2d_result.get('avg_r_squared_global'),
                'window_optimization': voigt2d_result.get('window_optimization'),
                'timestamp': voigt2d_result.get('timestamp')
            }

            # CRITICAL: Preserve 2D visualization data for simultaneous multi-peak fits
            if 'region_2d' in voigt2d_result:
                result['region_2d'] = voigt2d_result['region_2d']
            if 'fitted_2d_surface' in voigt2d_result:
                result['fitted_2d_surface'] = voigt2d_result['fitted_2d_surface']
            if 'individual_surfaces' in voigt2d_result:
                result['individual_surfaces'] = voigt2d_result['individual_surfaces']
            if 'all_peaks' in voigt2d_result:
                result['all_peaks'] = voigt2d_result['all_peaks']
            if 'overlap_group_size' in voigt2d_result:
                result['overlap_group_size'] = voigt2d_result['overlap_group_size']

            return result

        except Exception as e:
            print(f"Failed to convert voigt2d result: {e}")
            # Return minimal fallback - extract what we can
            return {
                'assignment': voigt2d_result.get('assignment', 'Unknown'),
                'peak_x': voigt2d_result.get('peak_position', (0, 0))[0],
                'peak_y': voigt2d_result.get('peak_position', (0, 0))[1],
                'amplitude': 1000,
                'volume': 1000,
                'avg_r_squared': voigt2d_result.get('avg_r_squared', 0.5),
                'success': False,
                'method': 'voigt_2d_fitting_conversion_failed',
                'error': str(e)
            }

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

    def _extract_direct_intensity_for_ml_comparison(self, ppm_x, ppm_y, local_x, local_y):
        """
        Extract direct intensity from spectrum data to compare with Voigt fitting
        for ML training internal standards. Used in core integrator.
        """
        try:
            if not hasattr(self, 'nmr_data') or self.nmr_data is None:
                return {'intensity': 0.0, 'height': 0.0, 'intensity_ratio': 1.0}

            # Find nearest indices in full spectrum
            x_idx = np.argmin(np.abs(self.ppm_x_axis - ppm_x))
            y_idx = np.argmin(np.abs(self.ppm_y_axis - ppm_y))

            # Extract direct intensity at peak position
            direct_intensity = float(self.nmr_data[y_idx, x_idx])
            direct_height = abs(direct_intensity)

            # Calculate intensity from local data (for comparison with detection/fitting region)
            if len(local_y) > 0:
                local_max = float(np.max(local_y))
                intensity_ratio = direct_height / (local_max + 1e-8) if local_max > 0 else 1.0
            else:
                intensity_ratio = 1.0

            return {
                'intensity': direct_intensity,
                'height': direct_height,
                'intensity_ratio': intensity_ratio,
                'x_idx': int(x_idx),
                'y_idx': int(y_idx),
                'extraction_method': 'spectrum_lookup_core'
            }

        except Exception as e:
            # Silent failure - don't break ML data collection
            return {'intensity': 0.0, 'height': 0.0, 'intensity_ratio': 1.0}

def _clamp_guess_to_bounds(vector, lower_bounds, upper_bounds, margin=1e-6):
    """Clamp parameter vector to be safely inside the provided bounds."""
    clamped = []

    for value, lower, upper in zip(vector, lower_bounds, upper_bounds):
        lower_adj = lower + margin
        upper_adj = upper - margin if np.isfinite(upper) else upper

        if np.isfinite(upper_adj) and upper_adj <= lower_adj:
            lower_adj = lower
            upper_adj = upper

        if np.isfinite(upper_adj):
            clamped.append(min(max(value, lower_adj), upper_adj))
        else:
            clamped.append(max(value, lower_adj))

    return clamped

    def _extract_direct_intensity_for_ml_comparison(self, ppm_x, ppm_y, local_x, local_y):
        """
        Extract direct intensity from spectrum data to compare with Voigt fitting
        for ML training internal standards. Used in core integrator.
        """
        try:
            if not hasattr(self, 'nmr_data') or self.nmr_data is None:
                return {'intensity': 0.0, 'height': 0.0, 'intensity_ratio': 1.0}

            # Find nearest indices in full spectrum
            x_idx = np.argmin(np.abs(self.ppm_x_axis - ppm_x))
            y_idx = np.argmin(np.abs(self.ppm_y_axis - ppm_y))

            # Extract direct intensity at peak position
            direct_intensity = float(self.nmr_data[y_idx, x_idx])
            direct_height = abs(direct_intensity)

            # Calculate intensity from local data (for comparison with detection/fitting region)
            if len(local_y) > 0:
                local_max = float(np.max(local_y))
                intensity_ratio = direct_height / (local_max + 1e-8) if local_max > 0 else 1.0
            else:
                intensity_ratio = 1.0

            return {
                'intensity': direct_intensity,
                'height': direct_height,
                'intensity_ratio': intensity_ratio,
                'x_idx': int(x_idx),
                'y_idx': int(y_idx),
                'extraction_method': 'spectrum_lookup_core'
            }

        except Exception as e:
            # Silent failure - don't break ML data collection
            return {'intensity': 0.0, 'height': 0.0, 'intensity_ratio': 1.0}
