#!/usr/bin/env python3
"""
Series Processing Module

This module handles batch processing of multiple NMR spectra in series,
including progress tracking, result consolidation, and comprehensive analysis.

Classes:
- SeriesProcessor: Main series processing workflow
- SeriesAnalyzer: Analysis and statistics of series results
- BatchResults: Container for batch processing results

Author: Guillaume Mas
Date: 2025
"""

import os
import pandas as pd
import numpy as np
from datetime import datetime
import threading
import time
import json
from pathlib import Path

try:
    from .core_integrator import EnhancedVoigtIntegrator
    from .file_manager import NMRFileManager, DataValidator
except ImportError:
    # Fallback for direct script execution
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.utils.file_manager import NMRFileManager, DataValidator

# Import integration status check
try:
    from lunaNMR.core.core_integrator import INTEGRATED_DETECTION_AVAILABLE
    SERIES_INTEGRATION_AVAILABLE = INTEGRATED_DETECTION_AVAILABLE
except ImportError:
    SERIES_INTEGRATION_AVAILABLE = False

class BatchResults:
    """Container for batch processing results with analysis capabilities"""

    def __init__(self):
        self.results = {}
        self.metadata = {
            'start_time': None,
            'end_time': None,
            'total_spectra': 0,
            'successful_spectra': 0,
            'failed_spectra': 0,
            'processing_mode': 'unknown'
        }
        self.statistics = {}
        self.errors = []

    def add_result(self, spectrum_name, result_data):
        """Add a single spectrum result"""
        self.results[spectrum_name] = result_data

        # Update metadata
        if result_data.get('status') == 'success':
            self.metadata['successful_spectra'] += 1
        else:
            self.metadata['failed_spectra'] += 1
            if 'error' in result_data:
                self.errors.append({
                    'spectrum': spectrum_name,
                    'error': result_data['error']
                })

    def get_summary(self):
        """Get processing summary"""
        total = len(self.results)
        success = self.metadata['successful_spectra']
        failed = self.metadata['failed_spectra']

        if self.metadata['start_time'] and self.metadata['end_time']:
            duration = self.metadata['end_time'] - self.metadata['start_time']
            duration_str = str(duration).split('.')[0]  # Remove microseconds
        else:
            duration_str = "Unknown"

        summary = {
            'total_spectra': total,
            'successful': success,
            'failed': failed,
            'success_rate': (success / total * 100) if total > 0 else 0,
            'duration': duration_str,
            'processing_mode': self.metadata['processing_mode'],
            'error_count': len(self.errors)
        }

        return summary

    def export_results(self, output_folder):
        """Export all results to files"""
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Export summary
        summary_file = os.path.join(output_folder, "batch_summary.json")
        with open(summary_file, 'w') as f:
            json.dump({
                'metadata': self.metadata,
                'summary': self.get_summary(),
                'statistics': self.statistics,
                'errors': self.errors
            }, f, indent=2, default=str)

        # Export individual results
        if self.results:
            # Create consolidated integration results
            all_integrations = []
            detection_stats = []

            for spectrum_name, result in self.results.items():
                # Detection statistics
                detection_stats.append({
                    'spectrum': spectrum_name,
                    'status': result.get('status', 'unknown'),
                    'detected_peaks': result.get('detected_peaks', 0),
                    'total_peaks': result.get('total_peaks', 0),
                    'detection_rate': result.get('detection_rate', 0.0),
                    'noise_level': result.get('noise_level', 0.0),
                    'processing_time': result.get('processing_time', 0.0)
                })

                # Individual integrations
                if result.get('integration_results'):
                    for integration in result['integration_results']:
                        integration_with_spectrum = integration.copy()
                        integration_with_spectrum['spectrum'] = spectrum_name
                        all_integrations.append(integration_with_spectrum)

            # Save detection statistics
            if detection_stats:
                stats_file = os.path.join(output_folder, "detection_statistics.csv")
                pd.DataFrame(detection_stats).to_csv(stats_file, index=False)

            # Save consolidated integrations
            if all_integrations:
                integration_file = os.path.join(output_folder, "all_integrations.csv")
                pd.DataFrame(all_integrations).to_csv(integration_file, index=False, float_format='%.6f')

        return output_folder

class SeriesProcessor:
    """Main series processing workflow manager"""

    def __init__(self):
        self.integrator = EnhancedVoigtIntegrator()
        self.file_manager = NMRFileManager()
        self.validator = DataValidator()

        # Processing parameters
        self.processing_mode = 'full_detection'
        self.reference_peaks = None
        self.output_folder = None

        # Reference spectrum tracking
        self.reference_spectrum_path = None
        self.reference_result = None

        # Progress tracking
        self.progress_callback = None

        # Voigt fitting options for series integration
        self.use_voigt_fitting = True
        self.use_parallel_processing = True  # Default enabled
        self.enhanced_parallel_voigt = True  # Enable enhanced parallel Voigt fitting
        self.parallel_voigt_threshold = 3    # Minimum peaks for parallel processing
        self.use_global_optimization = False
        self.num_integrations = 3  # Default number of integrations per peak

        # PS2D linewidth reuse (C++ peakfit.cpp lines 586-607 logic)
        self.reference_linewidths = {}  # Stores {assignment: {x_sigma, x_gamma, y_sigma, y_gamma}}
        self.use_ps2d_linewidth_reuse = False  # Toggle for PS2D mode

        # INTEGRATION ENHANCEMENT: Integrated detection-fitting options
        self.integration_mode = 'standard'  # 'standard', 'integrated', 'adaptive'
        self.integration_parameters = {
            'enable_integrated_mode': False,  # Default off for backward compatibility
            'adaptive_thresholds': True,
            'multi_resolution_detection': True,
            'quality_filtering': True,
            'detection_confidence_threshold': 0.3,
            'max_integration_iterations': 5,
            'global_parameter_propagation': True
        }

        # GUI integrator reference (set by GUI for unified workflow)
        self.main_gui_integrator = None
        self.status_callback = None
        self.cancelled = False
        self.paused = False

        # Results
        self.batch_results = BatchResults()

    def set_processing_mode(self, mode):
        """Set processing mode"""
        self.processing_mode = mode
        self.integrator.set_processing_mode(mode)
        self.batch_results.metadata['processing_mode'] = mode

    def set_reference_peaks(self, peak_list_path_or_dataframe):
        """Set reference peak list and ensure it is a DataFrame."""
        if isinstance(peak_list_path_or_dataframe, str):
            self.reference_peaks = self.file_manager.load_peak_list(peak_list_path_or_dataframe)
        elif isinstance(peak_list_path_or_dataframe, list):
            print(f"   Converting list of {len(peak_list_path_or_dataframe)} detected peaks to DataFrame.")
            # Standardize keys to ensure DataFrame compatibility
            standardized_peaks = []
            for i, peak in enumerate(peak_list_path_or_dataframe):
                p = {}
                p['Position_X'] = peak.get('ppm_x', peak.get('Position_X'))
                p['Position_Y'] = peak.get('ppm_y', peak.get('Position_Y'))
                p['Assignment'] = peak.get('assignment', peak.get('Assignment', f'Peak_{i+1}'))
                if p['Position_X'] is not None and p['Position_Y'] is not None:
                    standardized_peaks.append(p)
            self.reference_peaks = pd.DataFrame(standardized_peaks)

        elif hasattr(peak_list_path_or_dataframe, 'copy'):
            #self.reference_peaks = peak_list_path_or_dataframe.copy()
## addeed GM
             # This block now correctly handles DataFrames from the peak picker.
             df = peak_list_path_or_dataframe.copy()

             # --- START OF FIX ---
             # Standardize column names from the format used by the peak picker ('X_ppm')
             # to the format expected by the integrator ('Position_X').
             rename_map = {
                 'X_ppm': 'Position_X',
                 'Y_ppm': 'Position_Y'
             }
             df.rename(columns=rename_map, inplace=True)
             # --- END OF FIX ---
## addeed GM
        else:
            raise TypeError(f"Unsupported data type for reference peaks: {type(peak_list_path_or_dataframe)}")

        # Validate reference peaks (which should now be a DataFrame)
        if hasattr(self.reference_peaks, 'columns'):
            issues = self.validator.validate_peak_list_integrity(self.reference_peaks)
            if issues:
                fixed_peaks, fixes = self.validator.auto_fix_peak_list(self.reference_peaks)
                self.reference_peaks = fixed_peaks
                if self.status_callback:
                    self.status_callback(f"Fixed reference peaks: {', '.join(fixes)}")
        else:
            print(f"   Warning: Could not validate reference peaks of type: {type(self.reference_peaks)}")

    def set_callbacks(self, progress_callback=None, status_callback=None):
        """Set callback functions for progress and status updates"""
        self.progress_callback = progress_callback
        self.status_callback = status_callback

    def set_voigt_fitting_options(self, use_voigt_fitting=True, use_parallel_processing=True, use_global_optimization=False, num_integrations=3, **voigt_params):
        """Set Voigt fitting options for series integration"""
        self.use_voigt_fitting = use_voigt_fitting
        self.use_parallel_processing = use_parallel_processing
        self.use_global_optimization = use_global_optimization
        self.num_integrations = num_integrations

        # Store detailed Voigt parameters
        self.voigt_params = voigt_params

        print(f"🔧 Series processor Voigt fitting options:")
        print(f"   Use Voigt fitting: {self.use_voigt_fitting}")
        print(f"   Use parallel processing: {self.use_parallel_processing}")
        print(f"   Use global optimization: {self.use_global_optimization}")
        print(f"   Number of integrations per peak: {self.num_integrations}")

        if voigt_params:
            print(f"📊 Advanced Voigt parameters:")
            for key, value in voigt_params.items():
                print(f"   {key}: {value}")

    def set_integration_mode(self, mode='standard', **integration_params):
        """
        Set integration mode for series processing (INTEGRATION ENHANCEMENT)

        Args:
            mode: 'standard', 'integrated', or 'adaptive'
            **integration_params: integration-specific parameters
        """
        if mode not in ['standard', 'integrated', 'adaptive']:
            raise ValueError("Mode must be 'standard', 'integrated', or 'adaptive'")

        if mode != 'standard' and not SERIES_INTEGRATION_AVAILABLE:
            print("⚠️ Integrated detection-fitting not available for series, using standard mode")
            mode = 'standard'

        self.integration_mode = mode
        self.integration_parameters.update(integration_params)

        if mode == 'integrated':
            self.integration_parameters['enable_integrated_mode'] = True
            print("🚀 Series integrated detection-fitting mode enabled")
        elif mode == 'adaptive':
            self.integration_parameters['enable_integrated_mode'] = True
            self.integration_parameters['adaptive_thresholds'] = True
            self.integration_parameters['multi_resolution_detection'] = True
            self.integration_parameters['global_parameter_propagation'] = True
            print("🎯 Series adaptive integrated detection-fitting mode enabled")
        else:
            self.integration_parameters['enable_integrated_mode'] = False
            print("📊 Series standard detection-fitting mode")

        # Apply integration mode to main integrator if available
        if hasattr(self, 'main_gui_integrator') and self.main_gui_integrator:
            try:
                self.main_gui_integrator.set_integration_mode(mode, **integration_params)
                print("   ✅ Applied to main GUI integrator")
            except AttributeError:
                print("   ⚠️ Main integrator does not support integration mode")

        # Apply to local integrator
        if hasattr(self.integrator, 'set_integration_mode'):
            self.integrator.set_integration_mode(mode, **integration_params)
            print("   ✅ Applied to local integrator")

        return self.integration_parameters.copy()

    def get_integration_status(self):
        """Get current integration mode and parameters for series processing"""
        status = {
            'mode': self.integration_mode,
            'parameters': self.integration_parameters.copy(),
            'series_integration_available': SERIES_INTEGRATION_AVAILABLE,
            'main_integrator_status': None,
            'local_integrator_status': None
        }

        # Get status from main integrator
        if hasattr(self, 'main_gui_integrator') and self.main_gui_integrator:
            if hasattr(self.main_gui_integrator, 'get_integration_status'):
                status['main_integrator_status'] = self.main_gui_integrator.get_integration_status()

        # Get status from local integrator
        if hasattr(self.integrator, 'get_integration_status'):
            status['local_integrator_status'] = self.integrator.get_integration_status()

        return status

    def set_main_gui_integrator(self, gui_integrator):
        """Set reference to main GUI integrator for unified workflow (DEPRECATED - use set_main_gui)"""
        self.main_gui_integrator = gui_integrator

    def set_main_gui(self, main_gui):
        """Set main GUI reference for gold standard workflow"""
        self.main_gui = main_gui
        self.main_gui_integrator = main_gui.integrator  # Backward compatibility
        print(f"✅ Main GUI set for series processing: {type(main_gui)}")

        # Continue with existing logic...
        gui_integrator = main_gui.integrator

        # Sync integration mode with main integrator if available
        if (hasattr(gui_integrator, 'get_integration_status') and
            self.integration_parameters.get('enable_integrated_mode', False)):

            try:
                main_status = gui_integrator.get_integration_status()
                if main_status['mode'] != self.integration_mode:
                    print(f"🔄 Syncing series integration mode with main GUI: {main_status['mode']}")
                    self.integration_mode = main_status['mode']
                    self.integration_parameters.update(main_status['parameters'])
            except Exception as e:
                print(f"⚠️ Could not sync integration mode: {e}")
        print(f"🔧 Series processor will use main GUI integrator for processing")

    def process_series(self, nmr_files, output_folder=None, reference_spectrum_path=None):
        """Process a series of NMR files with proper reference handling"""
        if self.reference_peaks is None or (hasattr(self.reference_peaks, 'empty') and self.reference_peaks.empty):
            raise ValueError("Reference peaks must be set before processing")

        self.output_folder = output_folder or f"series_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        self.reference_spectrum_path = reference_spectrum_path
        self.cancelled = False
        self.paused = False

        # Initialize results
        self.batch_results = BatchResults()
        self.batch_results.metadata['start_time'] = datetime.now()
        self.batch_results.metadata['total_spectra'] = len(nmr_files)
        self.batch_results.metadata['processing_mode'] = self.processing_mode
        self.batch_results.metadata['reference_spectrum'] = os.path.basename(reference_spectrum_path) if reference_spectrum_path else None
        self.batch_results.metadata['processing_files'] = nmr_files.copy()  # Store original file paths

        # Create output directory
        if self.output_folder and not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Process each file
        for i, nmr_file in enumerate(nmr_files):
            if self.cancelled:
                break

            # Handle pause
            while self.paused and not self.cancelled:
                time.sleep(0.1)

            spectrum_name = os.path.basename(nmr_file)
            is_reference = (self.reference_spectrum_path and
                          os.path.abspath(nmr_file) == os.path.abspath(self.reference_spectrum_path))

            # Update progress
            ref_indicator = " [REFERENCE]" if is_reference else ""
            if self.progress_callback:
                self.progress_callback(
                    (i / len(nmr_files)) * 100,
                    f"Processing {spectrum_name}{ref_indicator} ({i+1}/{len(nmr_files)})",
                    f"Started processing {spectrum_name}{ref_indicator}"
                )

            # Process single spectrum
            start_time = time.time()
            result = self.process_single_spectrum(nmr_file)
            processing_time = time.time() - start_time

            # Mark result appropriately
            result['processing_time'] = processing_time
            result['is_reference'] = is_reference
            result['spectrum_type'] = 'reference' if is_reference else 'sample'

            # Enhanced logging for reference spectrum
            if is_reference:
                print(f"📍 REFERENCE SPECTRUM: {spectrum_name}")
                print(f"   Status: {result['status']}")
                if result['status'] == 'failed':
                    print(f"   Error: {result.get('error', 'Unknown error')}")
                else:
                    print(f"   Detected peaks: {result.get('detected_peaks', 0)}/{result.get('total_peaks', 0)}")
                    print(f"   Detection rate: {result.get('detection_rate', 0):.1f}%")
                self.reference_result = result.copy()

            self.batch_results.add_result(spectrum_name, result)

            # Save individual results if requested
            if self.output_folder:
                self.save_individual_results(spectrum_name, result)
                # Note: Individual spectrum plotting disabled due to threading issues on macOS
                # self._save_individual_spectrum_plot(nmr_file, spectrum_name, result)

            # Update progress
            status = "✅" if result['status'] == 'success' else "❌"
            if self.progress_callback:
                self.progress_callback(
                    ((i + 1) / len(nmr_files)) * 100,
                    f"Completed {spectrum_name}{ref_indicator}",
                    f"{status} {spectrum_name}: {result.get('detection_rate', 0):.1f}% detection{ref_indicator}"
                )

        # Finalize results
        self.batch_results.metadata['end_time'] = datetime.now()

        # Summary of reference spectrum processing
        if self.reference_spectrum_path:
            ref_name = os.path.basename(self.reference_spectrum_path)
            if self.reference_result:
                print(f"\n📍 REFERENCE SPECTRUM SUMMARY:")
                print(f"   File: {ref_name}")
                print(f"   Status: {self.reference_result['status']}")
                if self.reference_result['status'] == 'success':
                    print(f"   ✅ Successfully integrated with {self.reference_result.get('detected_peaks', 0)} peaks")
                else:
                    print(f"   ❌ Failed: {self.reference_result.get('error', 'Unknown error')}")
            else:
                print(f"\n⚠️  WARNING: Reference spectrum {ref_name} was not processed!")
        else:
            print(f"\n⚠️  WARNING: No reference spectrum was specified!")

        # Calculate final statistics
        self.calculate_series_statistics()

        # Export results
        if self.output_folder:
            self.batch_results.export_results(self.output_folder)
            # Create comprehensive peak tracking table
            self._create_peak_tracking_table()

        return self.batch_results

    def _extract_direct_intensities_for_series(self, gui_integrator, is_reference=False):
        """Extract direct intensities from spectrum data without Voigt fitting"""
        integration_results = []

        if not hasattr(gui_integrator, 'nmr_data') or gui_integrator.nmr_data is None:
            if is_reference:
                print("   ⚠️ No NMR data available in integrator")
            return integration_results

        nmr_data = gui_integrator.nmr_data

        # Process each reference peak
        for idx, peak_row in self.reference_peaks.iterrows():
            try:
                # Get peak position from reference peaks
                ppm_x = peak_row.get('Position_X', 0.0)
                ppm_y = peak_row.get('Position_Y', 0.0)
                assignment = peak_row.get('Assignment', f'Peak_{idx+1}')

                # Convert PPM to array indices (similar to main_gui._extract_direct_intensities)
                if hasattr(gui_integrator, 'ppm_scale_x') and hasattr(gui_integrator, 'ppm_scale_y'):
                    # Find closest indices
                    x_idx = np.argmin(np.abs(gui_integrator.ppm_scale_x - ppm_x))
                    y_idx = np.argmin(np.abs(gui_integrator.ppm_scale_y - ppm_y))

                    # Extract intensity directly from spectrum
                    if 0 <= y_idx < nmr_data.shape[0] and 0 <= x_idx < nmr_data.shape[1]:
                        intensity = float(nmr_data[y_idx, x_idx])

                        # Create integration result in expected format
                        result = {
                            'assignment': assignment,
                            'ppm_x': ppm_x,
                            'ppm_y': ppm_y,
                            'height': intensity,
                            'volume': intensity * 0.1,  # Approximate volume as height * small factor
                            'snr': intensity / (np.std(nmr_data) + 1e-10),  # Simple SNR calculation
                            'position_x': ppm_x,
                            'position_y': ppm_y,
                            'intensity': intensity,
                            'method': 'direct_extraction'
                        }
                        integration_results.append(result)

                        if is_reference and idx < 3:  # Show first few for debugging
                            print(f"   Peak {assignment}: PPM({ppm_x:.3f}, {ppm_y:.3f}) -> Intensity: {intensity:.2e}")
                    else:
                        if is_reference:
                            print(f"   ⚠️ Peak {assignment}: indices out of bounds (x:{x_idx}, y:{y_idx})")
                        # Still add a zero result to maintain peak count
                        result = {
                            'assignment': assignment,
                            'ppm_x': ppm_x,
                            'ppm_y': ppm_y,
                            'height': 0.0,
                            'volume': 0.0,
                            'snr': 0.0,
                            'position_x': ppm_x,
                            'position_y': ppm_y,
                            'intensity': 0.0,
                            'method': 'direct_extraction'
                        }
                        integration_results.append(result)
                else:
                    if is_reference:
                        print(f"   ⚠️ No PPM scales available in integrator")
                    result = {
                        'assignment': assignment,
                        'ppm_x': ppm_x,
                        'ppm_y': ppm_y,
                        'height': 0.0,
                        'volume': 0.0,
                        'snr': 0.0,
                        'position_x': ppm_x,
                        'position_y': ppm_y,
                        'intensity': 0.0,
                        'method': 'direct_extraction'
                    }
                    integration_results.append(result)

            except Exception as e:
                if is_reference:
                    print(f"   ❌ Error extracting peak {assignment}: {e}")
                # Add failed result
                result = {
                    'assignment': peak_row.get('Assignment', f'Peak_{idx+1}'),
                    'ppm_x': peak_row.get('Position_X', 0.0),
                    'ppm_y': peak_row.get('Position_Y', 0.0),
                    'height': 0.0,
                    'volume': 0.0,
                    'snr': 0.0,
                    'position_x': peak_row.get('Position_X', 0.0),
                    'position_y': peak_row.get('Position_Y', 0.0),
                    'intensity': 0.0,
                    'method': 'direct_extraction'
                }
                integration_results.append(result)

        return integration_results

    def _create_dummy_progress_dialog(self):
        """Create a dummy progress dialog to avoid GUI threading issues"""
        class DummyProgressDialog:
            def __init__(self):
                self.cancelled = False
                self.paused = False
                self.completed_tasks = 0
                self.failed_tasks = 0
                self.progress_var = type('ProgressVar', (object,), {
                    'get': lambda self: 0,
                    'set': lambda self, value: None
                })()

            def update_progress(self, *args, **kwargs):
                # Silently ignore all progress updates
                pass

            def complete(self, message):
                # Silently ignore completion messages
                pass

        return DummyProgressDialog()

    def process_single_spectrum(self, nmr_file, peak_list_path=None):
        """Process a single NMR spectrum with optional peak list (BACKWARD COMPATIBLE)"""
        spectrum_name = os.path.basename(nmr_file)
        is_reference = (self.reference_spectrum_path and
                      os.path.abspath(nmr_file) == os.path.abspath(self.reference_spectrum_path))

        # Check workflow mode if GUI is available
        workflow_mode = "peak_list"  # Default
        if hasattr(self, 'main_gui') and self.main_gui is not None:
            workflow_mode = getattr(self.main_gui, 'workflow_mode', tk.StringVar(value="peak_list")).get()

        # CRITICAL: Use main GUI - it has the perfect working workflow including detect_peaks() method
        if not hasattr(self, 'main_gui') or self.main_gui is None:
            return {
                'status': 'failed',
                'error': 'Main GUI not available - cannot use gold standard workflow',
                'detected_peaks': 0,
                'total_peaks': len(self.reference_peaks) if hasattr(self, 'reference_peaks') and self.reference_peaks is not None else 0,
                'detection_rate': 0.0
            }

        main_gui = self.main_gui
        gui_integrator = main_gui.integrator

        try:
            # Load NMR data into GUI integrator (same as single spectrum)
            success = gui_integrator.load_nmr_file(nmr_file)
            if not success:
                return {
                    'status': 'failed',
                    'error': 'Failed to load NMR data',
                    'detected_peaks': 0,
                    'total_peaks': len(self.reference_peaks) if hasattr(self, 'reference_peaks') and self.reference_peaks is not None else 0,
                    'detection_rate': 0.0
                }

            if workflow_mode == "peak_list":
                # EXISTING WORKFLOW - PRESERVED
                if not hasattr(self, 'reference_peaks') or self.reference_peaks is None:
                    return {
                        'status': 'failed',
                        'error': 'No reference peaks available for peak list mode',
                        'detected_peaks': 0,
                        'total_peaks': 0,
                        'detection_rate': 0.0
                    }

                print("📋 Series processing in Peak List mode")
                # Set reference peaks in GUI integrator (same as single spectrum)
                gui_integrator.peak_list = self.reference_peaks.copy()

            else:
                # NEW S/N WORKFLOW
                print("🎯 Series processing in S/N Threshold mode")
                gui_integrator.peak_list = None

                # Get S/N parameters from main GUI
                sn_threshold = getattr(main_gui, 'sn_threshold', tk.DoubleVar(value=3.0)).get()
                expected_count = getattr(main_gui, 'expected_peak_count', tk.IntVar(value=50)).get()

                # Set S/N parameters
                gui_integrator.sn_threshold = sn_threshold
                gui_integrator.expected_peak_count = expected_count
                gui_integrator.set_processing_mode('sn_native')

                print(f"   S/N Threshold: {sn_threshold}")
                print(f"   Expected Peak Count: {expected_count}")

            # CRITICAL: Apply detailed Voigt parameters to GUI integrator
            if hasattr(self, 'voigt_params') and self.voigt_params:
                # Apply fitting parameters
                if 'fitting_window_x' in self.voigt_params:
                    gui_integrator.fitting_parameters['fitting_window_x'] = self.voigt_params['fitting_window_x']
                if 'fitting_window_y' in self.voigt_params:
                    gui_integrator.fitting_parameters['fitting_window_y'] = self.voigt_params['fitting_window_y']
                if 'min_r_squared' in self.voigt_params:
                    gui_integrator.fitting_parameters['min_r_squared'] = self.voigt_params['min_r_squared']
                if 'max_iterations' in self.voigt_params:
                    gui_integrator.fitting_parameters['max_iterations'] = self.voigt_params['max_iterations']

                # Apply peak detection parameters via gui_params (for multi-peak fitting)
                gui_params = getattr(gui_integrator, 'gui_params', {})
                if 'peak_height_threshold' in self.voigt_params:
                    gui_params['height_threshold'] = self.voigt_params['peak_height_threshold']
                if 'peak_distance_factor' in self.voigt_params:
                    gui_params['distance_factor'] = self.voigt_params['peak_distance_factor']
                if 'peak_prominence_threshold' in self.voigt_params:
                    gui_params['prominence_threshold'] = self.voigt_params['peak_prominence_threshold']
                if 'smoothing_sigma' in self.voigt_params:
                    gui_params['smoothing_sigma'] = self.voigt_params['smoothing_sigma']
                if 'max_peaks_fit' in self.voigt_params:
                    gui_params['max_peaks_fit'] = self.voigt_params['max_peaks_fit']
                if 'max_optimization_iterations' in self.voigt_params:
                    gui_params['max_optimization_iterations'] = self.voigt_params['max_optimization_iterations']

                gui_integrator.gui_params = gui_params

                if is_reference:
                    print(f"   ✅ Applied {len(self.voigt_params)} Voigt parameters to GUI integrator")

            # Apply global optimization, parallel processing, and Voigt fitting settings to GUI
            if hasattr(self, 'use_global_optimization'):
                main_gui.use_global_optimization.set(self.use_global_optimization)
            if hasattr(self, 'use_parallel_processing'):
                main_gui.use_parallel_processing.set(self.use_parallel_processing)
            # CRITICAL FIX: Apply Voigt fitting toggle to BOTH GUI toggles
            if hasattr(self, 'use_voigt_fitting'):
                # Sync both single spectrum and series integration Voigt toggles
                main_gui.use_voigt_fitting.set(self.use_voigt_fitting)
                if hasattr(main_gui, 'series_use_voigt_fitting'):
                    main_gui.series_use_voigt_fitting.set(self.use_voigt_fitting)

            if is_reference:
                print(f"   ✅ Applied settings - Voigt: {getattr(self, 'use_voigt_fitting', True)}, Global opt: {getattr(self, 'use_global_optimization', False)}, Parallel: {getattr(self, 'use_parallel_processing', True)}")

            # STEP 1: Use EXACT same detection method from GUI (same as single spectrum)
            # Call the GUI's detect_peaks() method which handles all modes perfectly
            if is_reference:
                print(f"🔍 Calling GUI detect_peaks() method for {spectrum_name}")

            # This calls the SAME detect_peaks() method that works perfectly in single spectrum
            main_gui.detect_peaks()

            # Check if detection was successful by looking at fitted_peaks or peak_list
            if hasattr(gui_integrator, 'fitted_peaks') and gui_integrator.fitted_peaks is not None:
                detected_peaks = gui_integrator.fitted_peaks
                if is_reference:
                    print(f"   Detection successful: found fitted_peaks with {len(detected_peaks)} peaks")
            elif hasattr(gui_integrator, 'peak_list') and gui_integrator.peak_list is not None:
                detected_peaks = gui_integrator.peak_list
                if is_reference:
                    print(f"   Detection successful: using peak_list with {len(detected_peaks)} peaks")
            else:
                return {
                    'status': 'failed',
                    'error': 'Peak detection failed - no peaks found',
                    'detected_peaks': 0,
                    'total_peaks': len(self.reference_peaks),
                    'detection_rate': 0.0
                }

            # STEP 2: Process peaks based on Voigt fitting setting
            if self.use_voigt_fitting:
                if is_reference:
                    print(f"🔬 Using Voigt profile fitting for {spectrum_name}")

                # Perform Voigt fitting using the integrator's enhanced fitter
                try:
                    fitted_results = []

                    # Process each peak in the reference list
                    for idx, peak_row in self.reference_peaks.iterrows():
                        try:
                            # Get peak position
                            ppm_x = peak_row.get('Position_X', 0.0)
                            ppm_y = peak_row.get('Position_Y', 0.0)
                            assignment = peak_row.get('Assignment', f'Peak_{idx+1}')

                            if is_reference and idx == 0:
                                print(f"   Fitting Voigt profiles for peak {assignment} at ({ppm_x:.3f}, {ppm_y:.3f}) ppm")

                            # PS2D LINEWIDTH REUSE: Build constraints for series spectra
                            linewidth_constraints = None
                            if not is_reference and self.use_ps2d_linewidth_reuse:
                                linewidth_constraints = self._build_linewidth_constraints_from_reference(assignment)
                                if linewidth_constraints and idx == 0:
                                    print(f"   🔒 PS2D Mode: Using FIXED linewidths from reference (C++ peakfit.cpp:586-607)")

                            # Call the integrator's Voigt fitting method WITH constraints
                            if hasattr(gui_integrator, 'enhanced_fitter') and gui_integrator.enhanced_fitter:
                                result = gui_integrator.enhanced_fitter.fit_voigt_profile(
                                    gui_integrator.nmr_data,
                                    ppm_x, ppm_y,
                                    gui_integrator.ppm_scale_x,
                                    gui_integrator.ppm_scale_y,
                                    assignment=assignment,
                                    linewidth_constraints=linewidth_constraints  # NEW: Pass constraints
                                )
                                fitted_results.append(result)
                            else:
                                # Fallback to basic fitting
                                result = {
                                    'assignment': assignment,
                                    'ppm_x': ppm_x,
                                    'ppm_y': ppm_y,
                                    'success': False,
                                    'error': 'Enhanced fitter not available'
                                }
                                fitted_results.append(result)

                        except Exception as e:
                            if is_reference:
                                print(f"   ⚠️ Failed to fit peak {assignment}: {e}")
                            result = {
                                'assignment': peak_row.get('Assignment', f'Peak_{idx+1}'),
                                'ppm_x': peak_row.get('Position_X', 0.0),
                                'ppm_y': peak_row.get('Position_Y', 0.0),
                                'success': False,
                                'error': str(e)
                            }
                            fitted_results.append(result)

                    # Process results
                    successful_fits = len([r for r in fitted_results if r and r.get('success', True)])
                    total_peaks = len(self.reference_peaks)
                    detection_rate = (successful_fits / total_peaks * 100) if total_peaks > 0 else 0

                    # Convert to integration format for compatibility
                    integration_results = self._convert_voigt_to_integration_format(fitted_results, is_reference)

                    if is_reference:
                        print(f"   ✅ Voigt fitting complete: {successful_fits}/{total_peaks} peaks ({detection_rate:.1f}%)")

                    # PS2D LINEWIDTH REUSE: Extract linewidths from reference spectrum
                    if is_reference and self.use_ps2d_linewidth_reuse:
                        self.reference_linewidths = self._extract_reference_linewidths(
                            fitted_results,
                            spectrum_name
                        )

                    return {
                        'status': 'success' if successful_fits > 0 else 'failed',
                        'fitted_results': fitted_results,
                        'detected_peaks': successful_fits,
                        'total_peaks': total_peaks,
                        'detection_rate': detection_rate,
                        'fitted_peaks': self._convert_voigt_to_visualization_format(fitted_results),
                        'integration_results': integration_results,
                        'processing_method': 'voigt_profile_fitting'
                    }

                except Exception as e:
                    if is_reference:
                        print(f"   ❌ Error during Voigt fitting: {e}")
                    return {
                        'status': 'failed',
                        'error': f'Voigt fitting failed: {str(e)}',
                        'detected_peaks': 0,
                        'total_peaks': len(self.reference_peaks),
                        'detection_rate': 0.0
                    }
            else:
                # Direct intensity extraction without Voigt fitting
                if is_reference:
                    print(f"🔍 Using direct intensity extraction for {spectrum_name} (Voigt fitting OFF)")

                # Extract direct intensities from the spectrum data
                integration_results = self._extract_direct_intensities_for_series(gui_integrator, is_reference)

                # Calculate statistics
                successful_extractions = len([r for r in integration_results if r.get('height', 0) > 0])
                total_peaks = len(self.reference_peaks)
                detection_rate = (successful_extractions / total_peaks * 100) if total_peaks > 0 else 0

                if is_reference:
                    print(f"🔍 Direct extraction completed: {successful_extractions}/{total_peaks} peaks ({detection_rate:.1f}%)")
                    print(f"   Integration results count: {len(integration_results)}")

                return {
                    'status': 'success',
                    'detected_peaks': successful_extractions,
                    'total_peaks': total_peaks,
                    'detection_rate': detection_rate,
                    'noise_level': 0.0,  # No noise calculation for direct extraction
                    'integration_results': integration_results,
                    'fitted_peaks': [],  # No fitted peaks for direct extraction
                    'voigt_fits': [],
                    'processing_method': 'direct_intensity_extraction'
                }

        except Exception as e:
            return {
                'status': 'failed',
                'error': str(e),
                'detected_peaks': 0,
                'total_peaks': len(self.reference_peaks) if (
                    self.reference_peaks is not None and
                    (not hasattr(self.reference_peaks, 'empty') or not self.reference_peaks.empty)
                ) else 0,
                'detection_rate': 0.0
            }

    def _extract_reference_linewidths(self, fitted_results, spectrum_name):
        """
        Extract linewidths from reference spectrum for PS2D reuse (C++ lines 586-607)

        C++ PS2D logic:
        - Fit first plane/spectrum with full optimization
        - Extract linewidths (LW F1, LW F2, Voigt params)
        - Reuse these linewidths for all subsequent spectra (fix them)
        - Only fit intensity for subsequent spectra (40% speedup)

        Parameters
        ----------
        fitted_results : list
            Results from Voigt fitting of reference spectrum
        spectrum_name : str
            Name of reference spectrum (for logging)

        Returns
        -------
        dict : Reference linewidths keyed by assignment
        """
        if not fitted_results:
            print(f"   ⚠️ No fitted results to extract linewidths from")
            return {}

        linewidths = {}
        successful_extractions = 0

        for result in fitted_results:
            try:
                assignment = result.get('assignment', 'Unknown')

                # Extract X-dimension linewidths (1H)
                x_fit = result.get('x_fit', {})
                x_sigma = x_fit.get('sigma', 0.015)  # Gaussian width
                x_gamma = x_fit.get('gamma', 0.015)  # Lorentzian width

                # Extract Y-dimension linewidths (15N/13C)
                y_fit = result.get('y_fit', {})
                y_sigma = y_fit.get('sigma', 0.4)    # Gaussian width
                y_gamma = y_fit.get('gamma', 0.4)    # Lorentzian width

                # Store linewidths for this peak
                linewidths[assignment] = {
                    'x_sigma': x_sigma,
                    'x_gamma': x_gamma,
                    'y_sigma': y_sigma,
                    'y_gamma': y_gamma,
                    'x_r_squared': x_fit.get('r_squared', 0),
                    'y_r_squared': y_fit.get('r_squared', 0)
                }

                successful_extractions += 1

            except Exception as e:
                print(f"   ⚠️ Failed to extract linewidths for {result.get('assignment', 'Unknown')}: {e}")

        if successful_extractions > 0:
            print(f"   ✅ PS2D Linewidth Reuse: Extracted linewidths from {successful_extractions} reference peaks")
            print(f"      These linewidths will be FIXED for all series spectra (C++ peakfit.cpp:586-607 logic)")

        return linewidths

    def _build_linewidth_constraints_from_reference(self, assignment):
        """
        Build linewidth constraints from reference values (PS2D C++ style)

        C++ PS2D approach (peakfit.cpp:588-606):
        - Linewidths from plane 1 are COPIED to subsequent planes
        - They are FIXED (not optimized) for planes 2, 3, 4, ...
        - Only intensity is fitted for subsequent planes

        Python implementation:
        - Return very tight bounds around reference values (effectively fixes them)
        - Allow ±1% tolerance for numerical stability
        - Format compatible with enhanced_voigt_fitter's linewidth_constraints parameter

        Parameters
        ----------
        assignment : str
            Peak assignment to get constraints for

        Returns
        -------
        dict or None : Linewidth constraints in format expected by enhanced_voigt_fitter
            {
                'sigma_bounds': (min, max),  # Very tight bounds for X-dimension sigma
                'gamma_bounds': (min, max)   # Very tight bounds for X-dimension gamma
            }
        """
        if not self.reference_linewidths or assignment not in self.reference_linewidths:
            return None

        ref_lw = self.reference_linewidths[assignment]

        # PS2D EXACT approach: FIX linewidths to reference values
        # Use very tight bounds (±1%) to effectively fix them while maintaining numerical stability
        tolerance = 0.01  # 1% tolerance for numerical stability

        # For 1D fitting (X-dimension is primary), use X linewidths
        # The enhanced_voigt_fitter expects sigma_bounds and gamma_bounds as tuples
        x_sigma = ref_lw['x_sigma']
        x_gamma = ref_lw['x_gamma']

        constraints = {
            'sigma_bounds': (x_sigma * (1 - tolerance), x_sigma * (1 + tolerance)),
            'gamma_bounds': (x_gamma * (1 - tolerance), x_gamma * (1 + tolerance))
        }

        return constraints

    def save_individual_results(self, spectrum_name, result):
        """Save individual spectrum results with enhanced compatibility"""
        if not self.output_folder:
            return

        base_name = os.path.splitext(spectrum_name)[0]

        # ENHANCED: Save integration results (works for both modes now)
        if result.get('integration_results'):
            integration_file = os.path.join(self.output_folder, f"{base_name}_integration.csv")
            try:
                df = pd.DataFrame(result['integration_results'])
                # Verify data quality before saving
                if len(df) > 0:
                    # Check for essential columns
                    required_cols = ['Assignment', 'Position_X', 'Position_Y', 'Height', 'Volume', 'SNR', 'Quality']
                    missing_cols = [col for col in required_cols if col not in df.columns]
                    if missing_cols:
                        print(f"   ⚠️ Missing columns in integration results: {missing_cols}")

                    # Check for non-zero values
                    non_zero_heights = (df['Height'] != 0).sum() if 'Height' in df.columns else 0
                    non_zero_volumes = (df['Volume'] != 0).sum() if 'Volume' in df.columns else 0

                    df.to_csv(integration_file, index=False, float_format='%.6f')
                    print(f"   💾 Saved integration results: {len(result['integration_results'])} peaks → {base_name}_integration.csv")
                    print(f"      Non-zero heights: {non_zero_heights}, Non-zero volumes: {non_zero_volumes}")
                else:
                    print(f"   ⚠️ Empty integration results dataframe for {base_name}")
            except Exception as e:
                print(f"   ❌ Failed to save integration results for {base_name}: {e}")
        else:
            print(f"   ⚠️ No integration_results found for {base_name}")

        # ENHANCED: Save Voigt fitting results if available (detailed curve parameters)
        if result.get('fitted_results') or result.get('voigt_fits'):
            voigt_results = result.get('fitted_results', result.get('voigt_fits', []))
            if voigt_results:
                voigt_file = os.path.join(self.output_folder, f"{base_name}_voigt_fits.csv")
                voigt_data = []

                for fit in voigt_results:
                    if isinstance(fit, dict):
                        voigt_data.append({
                            'Assignment': fit.get('assignment', 'Unknown'),
                            'Peak_Number': fit.get('peak_number', 0),
                            'Center_PPM': fit.get('center_ppm', fit.get('ppm_x', 0.0)),
                            'Amplitude': fit.get('amplitude', 0.0),
                            'Sigma': fit.get('sigma', 0.0),
                            'Gamma': fit.get('gamma', 0.0),
                            'R_Squared': fit.get('avg_r_squared', fit.get('r_squared', 0.0)),
                            'Success': fit.get('success', True),
                            'Quality': 'Good' if fit.get('avg_r_squared', 0) > 0.8 else 'Fair' if fit.get('avg_r_squared', 0) > 0.5 else 'Poor'
                        })

                if voigt_data:
                    pd.DataFrame(voigt_data).to_csv(voigt_file, index=False, float_format='%.6f')

        # LEGACY: Save detected peaks (backward compatibility)
        if result.get('fitted_peaks') and not result.get('integration_results'):
            # Only create this if we don't have integration results (avoid duplication)
            peaks_file = os.path.join(self.output_folder, f"{base_name}_detected_peaks.csv")
            peak_data = []

            for peak in result['fitted_peaks']:
                peak_data.append({
                    'Assignment': peak.get('assignment', 'Unknown'),
                    'Position_X': peak.get('ppm_x', 0.0),
                    'Position_Y': peak.get('ppm_y', 0.0),
                    'Detected': peak.get('detected', False),
                    'SNR': peak.get('snr', 0.0),
                    'Quality': peak.get('detection_quality', 'Unknown')
                })

            if peak_data:
                pd.DataFrame(peak_data).to_csv(peaks_file, index=False, float_format='%.6f')

    def calculate_series_statistics(self):
        """Calculate comprehensive series statistics"""
        if not self.batch_results.results or len(self.batch_results.results) == 0:
            return

        # Basic statistics
        detection_rates = []
        noise_levels = []
        processing_times = []

        for result in self.batch_results.results.values():
            if result['status'] == 'success':
                detection_rates.append(result.get('detection_rate', 0.0))
                noise_levels.append(result.get('noise_level', 0.0))
                processing_times.append(result.get('processing_time', 0.0))

        if detection_rates:
            self.batch_results.statistics = {
                'detection_rate': {
                    'mean': np.mean(detection_rates),
                    'median': np.median(detection_rates),
                    'std': np.std(detection_rates),
                    'min': np.min(detection_rates),
                    'max': np.max(detection_rates)
                },
                'noise_level': {
                    'mean': np.mean(noise_levels),
                    'median': np.median(noise_levels),
                    'std': np.std(noise_levels),
                    'min': np.min(noise_levels),
                    'max': np.max(noise_levels)
                },
                'processing_time': {
                    'mean': np.mean(processing_times),
                    'median': np.median(processing_times),
                    'total': np.sum(processing_times)
                }
            }

    def pause(self):
        """Pause processing"""
        self.paused = True

    def resume(self):
        """Resume processing"""
        self.paused = False

    def cancel(self):
        """Cancel processing"""
        self.cancelled = True
        self.paused = False

    def _create_peak_tracking_table(self):
        """Create comprehensive peak tracking table across all spectra"""
        if not self.batch_results.results:
            return

        print("Creating comprehensive peak tracking table...")

        # Get all unique peaks from reference
        if (self.reference_peaks is None or
            (hasattr(self.reference_peaks, 'empty') and self.reference_peaks.empty)):
            return

        # Initialize tracking data structure
        peak_tracking_data = []

        # Get reference peak list
        if hasattr(self.reference_peaks, 'iterrows'):
            # DataFrame
            reference_peak_list = [(idx, row) for idx, row in self.reference_peaks.iterrows()]
        else:
            # List
            reference_peak_list = [(idx, peak) for idx, peak in enumerate(self.reference_peaks)]

        # Process each reference peak
        for peak_idx, peak_data in reference_peak_list:
            if hasattr(peak_data, 'get'):
                # DataFrame row
                assignment = peak_data.get('Assignment', f'Peak_{peak_idx+1}')
                ref_x = peak_data.get('Position_X', 0.0)
                ref_y = peak_data.get('Position_Y', 0.0)
            else:
                # Dictionary
                assignment = peak_data.get('Assignment', f'Peak_{peak_idx+1}')
                ref_x = peak_data.get('Position_X', 0.0)
                ref_y = peak_data.get('Position_Y', 0.0)

            # Initialize row for this peak
            peak_row = {
                'Peak_Number': peak_idx + 1,
                'Assignment': assignment,
                'Reference_X': ref_x,
                'Reference_Y': ref_y
            }

            # Add data for each spectrum
            for spectrum_name, result in self.batch_results.results.items():
                spectrum_key = spectrum_name.replace('.ft', '').replace('.', '_')

                # Initialize default values
                peak_row[f'{spectrum_key}_Detected'] = False
                peak_row[f'{spectrum_key}_Height'] = 0.0
                peak_row[f'{spectrum_key}_Volume'] = 0.0
                peak_row[f'{spectrum_key}_SNR'] = 0.0
                peak_row[f'{spectrum_key}_Quality'] = 'Not Detected'
                peak_row[f'{spectrum_key}_Position_X'] = ref_x
                peak_row[f'{spectrum_key}_Position_Y'] = ref_y

                # ENHANCED: Find matching peak in integration results (works for both modes now)
                if result.get('integration_results'):
                    for integration in result['integration_results']:
                        if integration.get('Assignment') == assignment:
                            peak_row[f'{spectrum_key}_Detected'] = integration.get('Detected', True)
                            peak_row[f'{spectrum_key}_Height'] = integration.get('Height', 0.0)
                            peak_row[f'{spectrum_key}_Volume'] = integration.get('Volume', 0.0)
                            peak_row[f'{spectrum_key}_SNR'] = integration.get('SNR', 0.0)
                            peak_row[f'{spectrum_key}_Quality'] = integration.get('Quality', 'Unknown')
                            peak_row[f'{spectrum_key}_Position_X'] = integration.get('Position_X', ref_x)
                            peak_row[f'{spectrum_key}_Position_Y'] = integration.get('Position_Y', ref_y)

                            # Add Voigt-specific information if available
                            if integration.get('R_Squared') is not None:
                                peak_row[f'{spectrum_key}_R_Squared'] = integration.get('R_Squared', 0.0)
                                peak_row[f'{spectrum_key}_Fit_Success'] = integration.get('Fit_Success', False)

                            break

            peak_tracking_data.append(peak_row)

        # Create DataFrame and save
        if peak_tracking_data:
            tracking_df = pd.DataFrame(peak_tracking_data)

            # Save comprehensive tracking table
            tracking_file = os.path.join(self.output_folder, "comprehensive_peak_tracking.csv")
            tracking_df.to_csv(tracking_file, index=False, float_format='%.6f')

            print(f"📊 Created peak tracking table with {len(peak_tracking_data)} peaks across {len(self.batch_results.results)} spectra")

            # Create simplified intensity matrix
            self._create_intensity_matrix(tracking_df)

            # Create peak presence/absence matrix
            self._create_detection_matrix(tracking_df)
        else:
            print("⚠️ No peak tracking data to export")

    def _create_intensity_matrix(self, tracking_df):
        """Create simplified intensity matrix (Peak vs Spectrum)"""
        intensity_data = []

        # Get spectrum names from columns
        spectrum_columns = [col for col in tracking_df.columns if col.endswith('_Height')]
        spectrum_names = [col.replace('_Height', '') for col in spectrum_columns]

        for _, row in tracking_df.iterrows():
            intensity_row = {
                'Peak_Number': row['Peak_Number'],
                'Assignment': row['Assignment']
            }

            for spectrum_name in spectrum_names:
                intensity_row[spectrum_name] = row.get(f'{spectrum_name}_Height', 0.0)

            intensity_data.append(intensity_row)

        if intensity_data:
            intensity_df = pd.DataFrame(intensity_data)
            intensity_file = os.path.join(self.output_folder, "peak_intensity_matrix.csv")
            intensity_df.to_csv(intensity_file, index=False, float_format='%.2f')
            print(f"Intensity matrix saved to: {intensity_file}")

    def _create_detection_matrix(self, tracking_df):
        """Create peak detection presence/absence matrix"""
        detection_data = []

        # Get spectrum names from columns
        detection_columns = [col for col in tracking_df.columns if col.endswith('_Detected')]
        spectrum_names = [col.replace('_Detected', '') for col in detection_columns]

        for _, row in tracking_df.iterrows():
            detection_row = {
                'Peak_Number': row['Peak_Number'],
                'Assignment': row['Assignment']
            }

            detected_count = 0
            for spectrum_name in spectrum_names:
                detected = row.get(f'{spectrum_name}_Detected', False)
                detection_row[spectrum_name] = 1 if detected else 0
                if detected:
                    detected_count += 1

            # Add summary statistics
            detection_row['Total_Detected'] = detected_count
            detection_row['Detection_Rate'] = detected_count / len(spectrum_names) if spectrum_names else 0
            detection_row['Reliability'] = 'High' if detection_row['Detection_Rate'] >= 0.8 else \
                                         'Medium' if detection_row['Detection_Rate'] >= 0.5 else 'Low'

            detection_data.append(detection_row)

        if detection_data:
            detection_df = pd.DataFrame(detection_data)
            detection_file = os.path.join(self.output_folder, "peak_detection_matrix.csv")
            detection_df.to_csv(detection_file, index=False)
            print(f"Detection matrix saved to: {detection_file}")

    def _save_individual_spectrum_plot(self, nmr_file, spectrum_name, result):
        """Save individual spectrum plot with peak overlays"""
        try:
            import matplotlib.pyplot as plt
            try:
                from .visualization import SpectrumPlotter
            except ImportError:
                # Fallback for direct script execution
                from lunaNMR.gui.visualization import SpectrumPlotter

            # Create plots directory
            plots_dir = os.path.join(self.output_folder, "individual_plots")
            os.makedirs(plots_dir, exist_ok=True)

            # Only create plots for successfully processed spectra
            if result.get('status') != 'success':
                return

            # Create figure and axis
            fig, ax = plt.subplots(1, 1, figsize=(12, 10))
            plotter = SpectrumPlotter(fig, ax)

            # Load the integrator temporarily to get the data
            temp_integrator = self.integrator.__class__()
            temp_integrator.peak_list = self.reference_peaks

            # Try to load NMR data
            try:
                temp_integrator.load_nmr_file(nmr_file)
            except:
                # If loading fails, skip this plot
                plt.close(fig)
                return

            # Set fitted peaks from result
            if result.get('fitted_peaks'):
                temp_integrator.fitted_peaks = result['fitted_peaks']

            # Plot spectrum
            plotter.plot_spectrum(temp_integrator)
            plotter.plot_peaks(temp_integrator, show_detected=True, show_assigned=True)

            # Add title and annotations
            is_reference = result.get('is_reference', False)
            ref_text = " [REFERENCE]" if is_reference else ""
            detection_rate = result.get('detection_rate', 0)

            ax.set_title(f'{spectrum_name}{ref_text}\n'
                        f'Detection: {result.get("detected_peaks", 0)}/{result.get("total_peaks", 0)} peaks '
                        f'({detection_rate:.1f}%)',
                        fontsize=14, fontweight='bold')

            # Add processing information
            info_text = (f"Processing Mode: {self.processing_mode}\n"
                        f"Noise Level: {result.get('noise_level', 0):.2e}\n"
                        f"Status: {result.get('status', 'Unknown')}")

            ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
                   verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3",
                   facecolor="white", alpha=0.8), fontsize=9)

            # Save plot
            base_name = spectrum_name.replace('.ft', '').replace('.', '_')
            plot_file = os.path.join(plots_dir, f"{base_name}_spectrum.png")
            fig.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close(fig)

            print(f"Individual plot saved: {plot_file}")

        except Exception as e:
            print(f"Warning: Failed to save individual plot for {spectrum_name}: {e}")

    def _process_single_spectrum_with_voigt_fitting(self, integrator, nmr_file, is_reference):
        """
        Process a single spectrum with complete Voigt fitting workflow
        This replicates the same workflow as the GUI's fit_all_peaks method
        """
        try:
            spectrum_name = os.path.basename(nmr_file)

            if is_reference:
                print(f"🔬 Voigt fitting for reference spectrum: {spectrum_name}")

            # Ensure integrator has the peak list
            if not hasattr(integrator, 'peak_list') or integrator.peak_list is None:
                integrator.peak_list = self.reference_peaks.copy()

            # Import ParallelPeakFitter for complete Voigt fitting
            try:
                from lunaNMR.processors.parallel_fitting import ParallelPeakFitter
            except ImportError as e:
                print(f"❌ Failed to import ParallelPeakFitter: {e}")
                # Fallback to basic integration
                return {
                    'status': 'failed',
                    'error': 'ParallelPeakFitter not available - falling back to basic integration',
                    'detected_peaks': 0,
                    'total_peaks': len(self.reference_peaks),
                    'detection_rate': 0.0
                }

            # Check if global optimization is enabled
            if self.use_global_optimization:
                return self._process_with_global_optimization(integrator, spectrum_name, is_reference)
            else:
                return self._process_with_linear_fitting(integrator, spectrum_name, is_reference)

        except Exception as e:
            return {
                'status': 'failed',
                'error': f"Voigt fitting error: {str(e)}",
                'detected_peaks': 0,
                'total_peaks': len(self.reference_peaks),
                'detection_rate': 0.0
            }

    def _process_with_linear_fitting(self, integrator, spectrum_name, is_reference):
        """Process with linear fitting workflow (parallel or sequential)"""
        try:
            from lunaNMR.processors.parallel_fitting import ParallelPeakFitter

            peak_list = integrator.peak_list
            total_count = len(peak_list)

            # Build all_peaks_context for 2D overlap detection and routing
            all_peaks_context = []
            for _, row in peak_list.iterrows():
                all_peaks_context.append({
                    'assignment': str(row.get('Assignment', 'Unknown')),
                    'x_ppm': float(row['Position_X']),
                    'y_ppm': float(row['Position_Y']),
                    'pos_x': float(row['Position_X']),
                    'pos_y': float(row['Position_Y'])
                })

            if is_reference:
                print(f"   🎯 2D overlap detection enabled with {len(all_peaks_context)} peaks context")

            if self.use_parallel_processing:
                if is_reference:
                    print(f"   Using parallel processing with {ParallelPeakFitter(integrator).max_workers} workers")

                # Create parallel fitter
                parallel_fitter = ParallelPeakFitter(integrator)

                # Define progress callback
                def progress_callback(completed, total, current_assignment):
                    if is_reference:
                        print(f"   Progress: {completed}/{total} peaks fitted")

                # Run parallel fitting with all_peaks_context
                fitted_results = parallel_fitter.fit_peaks_parallel(peak_list, progress_callback,
                                                                    all_peaks_context=all_peaks_context)
            else:
                if is_reference:
                    print(f"   Using sequential processing")

                # Sequential fitting
                fitted_results = []
                for i, (peak_idx, peak_row) in enumerate(peak_list.iterrows()):
                    peak_number = i + 1
                    assignment = peak_row.get('Assignment', f'Peak_{peak_number}')
                    peak_x = float(peak_row['Position_X'])
                    peak_y = float(peak_row['Position_Y'])

                    try:
                        result = integrator.enhanced_peak_fitting(peak_x, peak_y, assignment,
                                                                 all_peaks_context=all_peaks_context)
                        if result:
                            result['peak_number'] = peak_number
                            fitted_results.append(result)
                    except Exception as e:
                        if is_reference:
                            print(f"   ❌ Peak {peak_number} ({assignment}): {str(e)}")

            # Process results
            successful_fits = len([r for r in fitted_results if r.get('success', True)])
            detection_rate = (successful_fits / total_count * 100) if total_count > 0 else 0

            if is_reference:
                print(f"   ✅ Voigt fitting completed: {successful_fits}/{total_count} successful ({detection_rate:.1f}%)")

            # Store fitted results in integrator using our storage fix
            if fitted_results and hasattr(integrator, 'fitted_peaks'):
                integrator.fitted_peaks = fitted_results

            # CRITICAL FIX: Convert Voigt fitting results to integration format for export compatibility
            integration_results = self._convert_voigt_to_integration_format(fitted_results, is_reference)

            return {
                'status': 'success',
                'detected_peaks': successful_fits,
                'total_peaks': total_count,
                'detection_rate': detection_rate,
                'fitted_results': fitted_results,
                'fitted_peaks': self._convert_voigt_to_visualization_format(fitted_results),
                'voigt_fits': fitted_results,  # Backward compatibility
                'integration_results': integration_results  # For export compatibility
            }

        except Exception as e:
            return {
                'status': 'failed',
                'error': f"Linear fitting error: {str(e)}",
                'detected_peaks': 0,
                'total_peaks': len(self.reference_peaks),
                'detection_rate': 0.0
            }

    def _process_with_global_optimization(self, integrator, spectrum_name, is_reference):
        """Process with global optimization workflow (hybrid approach)"""
        try:
            if is_reference:
                print(f"   Using global optimization (2-stage hybrid approach)")

            peak_list = integrator.peak_list
            total_count = len(peak_list)

            # Stage 1: Linear analysis
            if is_reference:
                print(f"   Stage 1: Linear analysis for baseline statistics")
            linear_results = self._process_with_linear_fitting(integrator, spectrum_name, False)

            if linear_results['status'] != 'success' or not linear_results.get('fitted_results'):
                if is_reference:
                    print(f"   ❌ Stage 1 failed - falling back to linear results")
                return linear_results

            # Stage 2: Global optimization (simplified for series processing)
            if is_reference:
                print(f"   Stage 2: Constraint-based refinement")

            # For series processing, we'll use the linear results but apply global constraints
            # This is a simplified version of the full global optimization
            fitted_results = linear_results['fitted_results']

            # Apply global constraints and statistics
            refined_results = self._apply_global_constraints(fitted_results, is_reference)

            successful_fits = len([r for r in refined_results if r.get('success', True)])
            detection_rate = (successful_fits / total_count * 100) if total_count > 0 else 0

            if is_reference:
                print(f"   ✅ Global optimization completed: {successful_fits}/{total_count} successful ({detection_rate:.1f}%)")

            # Store fitted results in integrator
            if refined_results and hasattr(integrator, 'fitted_peaks'):
                integrator.fitted_peaks = refined_results

            # CRITICAL FIX: Convert Voigt fitting results to integration format for export compatibility
            integration_results = self._convert_voigt_to_integration_format(refined_results, is_reference)

            return {
                'status': 'success',
                'detected_peaks': successful_fits,
                'total_peaks': total_count,
                'detection_rate': detection_rate,
                'fitted_results': refined_results,
                'fitted_peaks': refined_results,
                'voigt_fits': refined_results,  # Backward compatibility
                'integration_results': integration_results  # For export compatibility
            }

        except Exception as e:
            return {
                'status': 'failed',
                'error': f"Global optimization error: {str(e)}",
                'detected_peaks': 0,
                'total_peaks': len(self.reference_peaks),
                'detection_rate': 0.0
            }

    def _apply_global_constraints(self, fitted_results, is_reference):
        """Apply global constraints to fitted results (simplified version)"""
        if not fitted_results:
            return fitted_results

        # Extract statistics from successful fits
        successful_results = [r for r in fitted_results if r.get('success', True) and 'avg_r_squared' in r]

        if not successful_results:
            return fitted_results

        # Calculate global statistics
        r_squared_values = [r['avg_r_squared'] for r in successful_results]
        if r_squared_values:
            median_r_squared = sorted(r_squared_values)[len(r_squared_values)//2]
            min_acceptable_r_squared = max(0.5, median_r_squared * 0.7)  # 70% of median, minimum 0.5

            if is_reference:
                print(f"   Global constraint: R² ≥ {min_acceptable_r_squared:.3f}")

            # Apply constraints
            refined_results = []
            for result in fitted_results:
                if result.get('success', True) and result.get('avg_r_squared', 0) >= min_acceptable_r_squared:
                    refined_results.append(result)
                else:
                    # Mark as failed due to global constraints
                    result_copy = result.copy() if isinstance(result, dict) else result
                    if isinstance(result_copy, dict):
                        result_copy['success'] = False
                        result_copy['constraint_failed'] = True
                    refined_results.append(result_copy)

            return refined_results

        return fitted_results

    def _convert_voigt_to_integration_format(self, fitted_results, is_reference=False):
        """
        Convert Voigt fitting results to integration format for export compatibility
        This ensures both Voigt and integration-only modes produce the same export format
        """
        if not fitted_results:
            return []

        integration_results = []

        for fit_result in fitted_results:
            if not isinstance(fit_result, dict):
                continue

            try:
                # Extract key information from Voigt fit result
                assignment = fit_result.get('assignment', 'Unknown')
                peak_number = fit_result.get('peak_number', 0)
                success = fit_result.get('success', True)

                # Add debug output for problematic results
                if is_reference:
                    print(f"   Converting result for {assignment}: success={success}")

                # Convert Voigt parameters to integration values
                if success and fit_result.get('success', True):
                    # Extract data from Voigt fitting result - check multiple possible field names
                    # Position X: try various field names used in Voigt fitting
                    position_x = 0.0
                    if fit_result.get('x_fit', {}).get('center') is not None:
                        position_x = float(fit_result['x_fit']['center'])
                    elif fit_result.get('center_ppm') is not None:
                        position_x = float(fit_result['center_ppm'])
                    elif isinstance(fit_result.get('peak_position'), (list, tuple)) and len(fit_result['peak_position']) > 0:
                        position_x = float(fit_result['peak_position'][0])

                    # Position Y: For NMR, Y is typically the intensity coordinate (amplitude)
                    position_y = 0.0
                    if fit_result.get('y_fit', {}).get('center') is not None:
                        position_y = float(fit_result['y_fit']['center'])
                    elif isinstance(fit_result.get('peak_position'), (list, tuple)) and len(fit_result['peak_position']) > 1:
                        position_y = float(fit_result['peak_position'][1])

                    # Height from Voigt amplitude (prefer X dimension for NMR)
                    height = 0.0
                    if fit_result.get('x_fit', {}).get('amplitude') is not None:
                        height = abs(float(fit_result['x_fit']['amplitude']))
                    elif fit_result.get('amplitude') is not None:
                        height = abs(float(fit_result['amplitude']))
                    elif fit_result.get('y_fit', {}).get('amplitude') is not None:
                        height = abs(float(fit_result['y_fit']['amplitude']))

                    # Calculate volume from Voigt parameters (proper integration)
                    x_sigma = fit_result.get('x_fit', {}).get('sigma', 0.0)
                    x_gamma = fit_result.get('x_fit', {}).get('gamma', 0.0)
                    y_sigma = fit_result.get('y_fit', {}).get('sigma', 0.0)
                    y_gamma = fit_result.get('y_fit', {}).get('gamma', 0.0)

                    if height > 0 and (x_sigma > 0 or x_gamma > 0):
                        # For NMR, primary integration is in X dimension
                        x_fwhm = 2 * np.sqrt(2 * np.log(2)) * x_sigma + 2 * x_gamma  # Voigt FWHM approximation
                        if y_sigma > 0 or y_gamma > 0:
                            y_fwhm = 2 * np.sqrt(2 * np.log(2)) * y_sigma + 2 * y_gamma
                            volume = abs(height) * x_fwhm * y_fwhm * self.num_integrations
                        else:
                            volume = abs(height) * x_fwhm * self.num_integrations  # 1D integration
                    else:
                        volume = abs(height) * 0.02 * self.num_integrations  # Conservative fallback

                    # Use R² as quality indicator for SNR approximation
                    r_squared = fit_result.get('avg_r_squared', fit_result.get('r_squared', 0.0))
                    if isinstance(r_squared, (list, tuple)):
                        r_squared = float(r_squared[0]) if len(r_squared) > 0 else 0.0
                    else:
                        r_squared = float(r_squared) if r_squared is not None else 0.0

                    # Convert R² to realistic SNR values
                    if r_squared > 0.9:
                        snr = 50 + (r_squared - 0.9) * 500  # 50-100 for excellent fits
                    elif r_squared > 0.7:
                        snr = 10 + (r_squared - 0.7) * 200  # 10-50 for good fits
                    else:
                        snr = max(1.0, r_squared * 14)  # 1-10 for poor fits

                    quality = 'Excellent' if r_squared > 0.9 else 'Good' if r_squared > 0.8 else 'Fair' if r_squared > 0.5 else 'Poor'
                    detected = True

                    if is_reference:
                        print(f"     Extracted: pos_x={position_x:.3f}, pos_y={position_y:.1f}, height={height:.2e}, volume={volume:.2e}, R²={r_squared:.3f}")

                else:
                    # Failed fit - use reference values or defaults
                    if isinstance(fit_result.get('peak_position'), (list, tuple)) and len(fit_result['peak_position']) >= 2:
                        position_x = float(fit_result['peak_position'][0])
                        position_y = float(fit_result['peak_position'][1])
                    else:
                        position_x = 0.0
                        position_y = 0.0
                    height = 0.0
                    volume = 0.0
                    snr = 0.0
                    quality = 'Not Fitted'
                    detected = False

                # Create integration-compatible result
                integration_result = {
                    'Assignment': assignment,
                    'Peak_Number': peak_number,
                    'Position_X': float(position_x),
                    'Position_Y': float(position_y),
                    'Height': float(height),
                    'Volume': float(volume),
                    'SNR': float(snr),
                    'Quality': quality,
                    'Detected': detected,
                    # Additional Voigt-specific data for advanced users
                    'R_Squared': fit_result.get('avg_r_squared', 0.0),
                    'Sigma': fit_result.get('sigma', 0.0),
                    'Gamma': fit_result.get('gamma', 0.0),
                    'Fit_Success': success
                }

                integration_results.append(integration_result)

            except Exception as e:
                if is_reference:
                    print(f"   Warning: Failed to convert fit result for {fit_result.get('assignment', 'unknown')}: {e}")
                continue

        if is_reference and integration_results:
            print(f"   📊 Converted {len(integration_results)} Voigt results to integration format")
            if len(integration_results) > 0:
                sample = integration_results[0]
                print(f"   Sample converted result:")
                for key, value in sample.items():
                    print(f"     {key}: {value} ({type(value)})")

        return integration_results

    def _convert_voigt_to_visualization_format(self, fitted_results):
        """
        Convert Voigt fitting results to visualization format for peak plotting
        The visualization code expects peaks with fields: ppm_x, ppm_y, detected, assignment
        """
        if not fitted_results:
            return []

        visualization_peaks = []

        for fit_result in fitted_results:
            if not isinstance(fit_result, dict):
                continue

            try:
                # Extract key information from Voigt fit result
                assignment = fit_result.get('assignment', 'Unknown')
                success = fit_result.get('success', True)

                if success and fit_result.get('success', True):
                    # Extract positions from Voigt fitting result
                    ppm_x = 0.0
                    ppm_y = 0.0

                    # Try to get position from x_fit/y_fit structure
                    if fit_result.get('x_fit', {}).get('center') is not None:
                        ppm_x = float(fit_result['x_fit']['center'])
                    elif fit_result.get('ppm_x') is not None:
                        ppm_x = float(fit_result['ppm_x'])
                    elif fit_result.get('center_x') is not None:
                        ppm_x = float(fit_result['center_x'])

                    if fit_result.get('y_fit', {}).get('center') is not None:
                        ppm_y = float(fit_result['y_fit']['center'])
                    elif fit_result.get('ppm_y') is not None:
                        ppm_y = float(fit_result['ppm_y'])
                    elif fit_result.get('center_y') is not None:
                        ppm_y = float(fit_result['center_y'])

                    # Create visualization peak in expected format
                    viz_peak = {
                        'ppm_x': ppm_x,
                        'ppm_y': ppm_y,
                        'detected': True,  # All Voigt fitted peaks are "detected"
                        'assignment': assignment,
                        'peak_number': fit_result.get('peak_number', len(visualization_peaks) + 1),
                        'success': success,
                        # Keep original Voigt data for reference
                        'voigt_data': fit_result
                    }

                    visualization_peaks.append(viz_peak)

            except (ValueError, TypeError, KeyError) as e:
                print(f"⚠️ Warning: Failed to convert Voigt result for {assignment}: {e}")
                continue

        print(f"✅ Converted {len(fitted_results)} Voigt results → {len(visualization_peaks)} visualization peaks")
        return visualization_peaks

class SeriesAnalyzer:
    """Analysis and visualization of series processing results"""

    def __init__(self):
        self.results = None

    def load_results(self, batch_results):
        """Load batch results for analysis"""
        self.results = batch_results

    def analyze_detection_trends(self):
        """Analyze detection rate trends across series"""
        if self.results is None or not hasattr(self.results, 'results') or not self.results.results:
            return None

        detection_data = []
        for spectrum_name, result in self.results.results.items():
            if result['status'] == 'success':
                detection_data.append({
                    'spectrum': spectrum_name,
                    'detection_rate': result.get('detection_rate', 0.0),
                    'detected_peaks': result.get('detected_peaks', 0),
                    'total_peaks': result.get('total_peaks', 0),
                    'noise_level': result.get('noise_level', 0.0)
                })

        if not detection_data:
            return None

        df = pd.DataFrame(detection_data)

        analysis = {
            'data': df,
            'trends': {
                'best_spectrum': df.loc[df['detection_rate'].idxmax()]['spectrum'],
                'worst_spectrum': df.loc[df['detection_rate'].idxmin()]['spectrum'],
                'average_detection': df['detection_rate'].mean(),
                'detection_consistency': 1.0 - (df['detection_rate'].std() / df['detection_rate'].mean()) if df['detection_rate'].mean() > 0 else 0
            },
            'correlations': {}
        }

        # Analyze correlations
        if len(df) > 2:
            # Correlation between noise level and detection rate
            noise_detection_corr = df['noise_level'].corr(df['detection_rate'])
            analysis['correlations']['noise_vs_detection'] = noise_detection_corr

        return analysis

    def generate_quality_report(self):
        """Generate comprehensive quality assessment report"""
        if not self.results:
            return None

        summary = self.results.get_summary()

        # Quality thresholds
        thresholds = {
            'excellent_detection': 90.0,
            'good_detection': 70.0,
            'fair_detection': 50.0,
            'min_success_rate': 80.0
        }

        report = {
            'overall_grade': 'Unknown',
            'summary': summary,
            'quality_metrics': {},
            'recommendations': [],
            'issues': []
        }

        # Calculate quality metrics
        if summary['success_rate'] >= thresholds['min_success_rate']:
            if 'statistics' in self.results.statistics and 'detection_rate' in self.results.statistics:
                avg_detection = self.results.statistics['detection_rate']['mean']

                if avg_detection >= thresholds['excellent_detection']:
                    report['overall_grade'] = 'Excellent'
                elif avg_detection >= thresholds['good_detection']:
                    report['overall_grade'] = 'Good'
                elif avg_detection >= thresholds['fair_detection']:
                    report['overall_grade'] = 'Fair'
                else:
                    report['overall_grade'] = 'Poor'
                    report['issues'].append(f"Low average detection rate: {avg_detection:.1f}%")

                report['quality_metrics']['average_detection_rate'] = avg_detection
                report['quality_metrics']['detection_consistency'] = 1.0 - (
                    self.results.statistics['detection_rate']['std'] / avg_detection
                ) if avg_detection > 0 else 0
            else:
                report['overall_grade'] = 'Poor'
                report['issues'].append("No detection statistics available")
        else:
            report['overall_grade'] = 'Failed'
            report['issues'].append(f"Low success rate: {summary['success_rate']:.1f}%")

        # Generate recommendations
        if report['overall_grade'] in ['Poor', 'Fair']:
            report['recommendations'].append("Consider adjusting detection parameters")
            report['recommendations'].append("Review noise threshold settings")
            report['recommendations'].append("Check data quality and preprocessing")

        if summary['error_count'] > 0:
            report['recommendations'].append(f"Investigate {summary['error_count']} processing errors")

        return report

    def export_analysis(self, output_file):
        """Export analysis results to file"""
        if not self.results:
            return False

        analysis_data = {
            'timestamp': datetime.now().isoformat(),
            'summary': self.results.get_summary(),
            'statistics': self.results.statistics,
            'trends': self.analyze_detection_trends(),
            'quality_report': self.generate_quality_report()
        }

        try:
            with open(output_file, 'w') as f:
                json.dump(analysis_data, f, indent=2, default=str)
            return True
        except Exception as e:
            print(f"Failed to export analysis: {e}")
            return False
