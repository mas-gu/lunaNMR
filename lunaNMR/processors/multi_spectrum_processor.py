"""
Independent Multi-Spectrum Processor for NMR Series Analysis
Combines complete GUI decoupling with retrocompatibility features

Author: Guillaume Mas
Date: 2025
"""

import os
import pandas as pd
import numpy as np
from datetime import datetime
import threading
import time
from typing import List, Dict, Any, Optional
from datetime import datetime

# Independent imports - no GUI dependencies
from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
from lunaNMR.utils.file_manager import NMRFileManager

class MultiSpectrumProcessor:
    """
    Completely independent multi-spectrum processor.
    Has its own integrator instance and creates all necessary output files.
    """

    def __init__(self, voigt_params):
        """Initialize with independent integrator instance"""
        # Create OWN integrator instance - completely independent from GUI
        self.integrator = EnhancedVoigtIntegrator()
        self.file_manager = NMRFileManager()
        self.voigt_params = voigt_params
        self.reference_peaks = None
        self.processing_active = False
        self.progress_callback = None

        # Configure the internal integrator with provided parameters
        self.integrator.fitting_parameters.update(self.voigt_params)
        self.integrator.gui_params = self.voigt_params.copy()

        # Initialize output folder for retrocompatibility
        self.output_folder = None

    def process_nmr_series(self, nmr_files: List[str], reference_peaks: pd.DataFrame,
                          output_folder: str, peak_source_mode: str = 'reference',
                          progress_callback: Optional[callable] = None) -> Dict[str, Any]:
        """
        Main entry point for processing NMR series.
        Returns comprehensive results in both new and legacy-compatible formats.
        """
        self.reference_peaks = reference_peaks.copy()
        self.output_folder = output_folder
        self.progress_callback = progress_callback
        self.processing_active = True

        try:
            # Initialize comprehensive batch results structure
            batch_results = self._initialize_comprehensive_batch_results(
                nmr_files, output_folder, peak_source_mode
            )

            total_spectra = len(nmr_files)
            successful_spectra = 0
            total_peaks_processed = 0
            total_successful_fits = 0

            # Process each spectrum independently
            for i, nmr_file in enumerate(nmr_files):
                if not self.processing_active:
                    break

                spectrum_name = os.path.basename(nmr_file)
                spectrum_key = os.path.splitext(spectrum_name)[0]

                # Update progress
                progress = (i / total_spectra) * 100
                if self.progress_callback:
                    self.progress_callback(
                        progress,
                        f"Processing spectrum {i+1}/{total_spectra}",
                        f"Loading {spectrum_name}"
                    )

                try:
                    # Process single spectrum using independent integrator
                    spectrum_result = self._process_single_spectrum_independent(
                        nmr_file, spectrum_name, i+1, total_spectra
                    )

                    if spectrum_result['success']:
                        successful_spectra += 1
                        total_successful_fits += spectrum_result['successful_fits']

                    total_peaks_processed += spectrum_result['total_peaks']
                    batch_results['results'][spectrum_name] = spectrum_result

                except Exception as e:
                    print(f" Failed to process {spectrum_name}: {e}")
                    batch_results['results'][spectrum_name] = {
                        'success': False,
                        'error': str(e),
                        'fitted_results': [],
                        'total_peaks': len(self.reference_peaks),
                        'successful_fits': 0,
                        'spectrum_file': spectrum_name
                    }

            # Finalize results with comprehensive statistics
            batch_results['summary'] = {
                'total_spectra': total_spectra,
                'successful': successful_spectra,
                'failed': total_spectra - successful_spectra,
                'success_rate': (successful_spectra / total_spectra * 100) if total_spectra > 0 else 0,
                'total_peaks_processed': total_peaks_processed,
                'total_successful_fits': total_successful_fits,
                'overall_detection_rate': (total_successful_fits / total_peaks_processed * 100) if total_peaks_processed > 0 else 0
            }

            # Create all output files for downstream compatibility
            self._create_comprehensive_output_files(batch_results)

            return batch_results

        except Exception as e:
            print(f" Multi-spectrum processing failed: {e}")
            return {'error': str(e), 'results': {}, 'summary': {}}
        finally:
            self.processing_active = False

    def _process_single_spectrum_independent(self, nmr_file: str, spectrum_name: str,
                                           spectrum_number: int, total_spectra: int) -> Dict[str, Any]:
        """Process a single spectrum using the independent integrator"""

        # REFINED STATE MANAGEMENT: Preserve essential integrator functionality
        # Save critical components before clearing
        essential_components = {}
        preserve_attrs = ['enhanced_fitter', 'fitting_parameters', 'gui_params', 'threshold_multiplier', 'processing_mode']

        for attr in preserve_attrs:
            if hasattr(self.integrator, attr):
                essential_components[attr] = getattr(self.integrator, attr)

        # Clear only problematic state attributes
        problematic_attrs = ['peak_list_path', 'nmr_file_path']
        for attr in problematic_attrs:
            if hasattr(self.integrator, attr):
                delattr(self.integrator, attr)

        # Clear NMR data only (let load_nmr_file handle proper loading)
        self.integrator.nmr_data = None
        self.integrator.nmr_dict = None
        self.integrator.ppm_x_axis = None
        self.integrator.ppm_y_axis = None

        # Always reset peak_list to None to ensure clean loading
        self.integrator.peak_list = None

        # Restore essential components
        for attr, value in essential_components.items():
            setattr(self.integrator, attr, value)
            print(f"   ♻️ Restored {attr}")

        # USE PROPER LOADING: Let integrator handle loading correctly
        print(f"   🔄 Loading NMR file using integrator's load method: {nmr_file}")

        # First attempt: Use integrator's load method
        load_success = False
        try:
            load_success = self.integrator.load_nmr_file(nmr_file)
            if load_success:
                print(f"   ✅ NMR file loaded successfully via load_nmr_file()")
            else:
                print(f"   ⚠️ load_nmr_file() returned False, trying fallback")
        except Exception as e:
            print(f"   ⚠️ load_nmr_file() failed: {e}, trying fallback")

        # Fallback: Direct loading only if regular loading fails
        if not load_success or self.integrator.nmr_data is None:
            print(f"   🔄 Attempting direct loading fallback")
            try:
                import nmrglue as ng
                self.integrator.nmr_dict, self.integrator.nmr_data = ng.pipe.read(nmr_file)

                # Only calculate axes if not already done
                if not hasattr(self.integrator, 'ppm_x_axis') or self.integrator.ppm_x_axis is None:
                    self.integrator._calculate_ppm_axes()

                # Only estimate noise if not already done
                if not hasattr(self.integrator, 'noise_level'):
                    self.integrator._estimate_noise_level()

                self.integrator.nmr_file_path = nmr_file
                print(f"   ✅ Direct loading successful: {self.integrator.nmr_data.shape}")

            except Exception as e:
                print(f"   ❌ Direct loading also failed: {e}")
                raise IOError(f"All loading methods failed for: {nmr_file}")
##


        # VERIFICATION: Ensure NMR data is properly loaded
        if (self.integrator.nmr_data is None or
            self.integrator.ppm_x_axis is None or
            self.integrator.ppm_y_axis is None):
            print(f"❌ NMR data verification failed for {spectrum_name}")
            print(f"   nmr_data: {type(self.integrator.nmr_data)}")
            print(f"   ppm_x_axis: {type(self.integrator.ppm_x_axis)}")
            print(f"   ppm_y_axis: {type(self.integrator.ppm_y_axis)}")
            raise IOError(f"NMR data not properly loaded after direct loading: {nmr_file}")

        print(f"✅ NMR data verified for {spectrum_name}: {self.integrator.nmr_data.shape}")

        # Set peak list in independent integrator
        self.integrator.peak_list = self.reference_peaks.copy()

        fitted_results = []
        successful_fits = 0
        total_peaks = len(self.reference_peaks)

        # Process each reference peak
        for peak_idx, peak_row in self.reference_peaks.iterrows():
            if not self.processing_active:
                break

            # Update progress for individual peaks
            peak_progress = ((spectrum_number - 1) / total_spectra +
                           (peak_idx + 1) / (len(self.reference_peaks) * total_spectra)) * 100

            assignment = peak_row.get('Assignment', f'Peak_{peak_idx + 1}')
            peak_x = float(peak_row['Position_X'])
            peak_y = float(peak_row['Position_Y'])

            if self.progress_callback:
                self.progress_callback(
                    peak_progress,
                    f"Spectrum {spectrum_number}/{total_spectra}",
                    f"Fitting {assignment} at ({peak_x:.3f}, {peak_y:.1f})"
                )

            # Perform fitting using independent integrator
            result = self.integrator.enhanced_peak_fitting(peak_x, peak_y, assignment)

            if result:
                # Add comprehensive metadata
                result['peak_number'] = peak_idx + 1
                result['spectrum_file'] = spectrum_name
                result['processing_mode'] = 'independent_series'
                fitted_results.append(result)
                successful_fits += 1

        # Calculate spectrum-level statistics
        success_rate = (successful_fits / total_peaks * 100) if total_peaks > 0 else 0

        return {
            'success': successful_fits > 0,
            'fitted_results': fitted_results,
            'total_peaks': total_peaks,
            'successful_fits': successful_fits,
            'success_rate': success_rate,
            'spectrum_file': spectrum_name,
            'integration_results': self._convert_to_integration_format(fitted_results)
        }

    def _create_tidy_results_file(self, batch_results: Dict[str, Any]):
        """
        Creates a single, tidy CSV file from all series results.
        This format is ideal for external analysis and plotting.
        """
        if not batch_results.get('results'):
            return

        tidy_data = []
        for spectrum_name, result_data in batch_results['results'].items():
            if not result_data.get('success', False):
                continue

            # Use the standardized integration results
            integration_results = result_data.get('integration_results', [])
            for peak in integration_results:
                row = {
                    'spectrum_name': os.path.splitext(spectrum_name)[0],
                    'assignment': peak.get('assignment', 'Unknown'),
                    'peak_number': peak.get('peak_number', 0),
                    'ppm_x': peak.get('ppm_x', 0.0),
                    'ppm_y': peak.get('ppm_y', 0.0),
                    'height': peak.get('height', 0.0),
                    'volume': peak.get('volume', 0.0),
                    'snr': peak.get('snr', 0.0),
                    'quality': peak.get('quality', 'Unknown'),
                    'r_squared': peak.get('r_squared', 0.0)
                }
                tidy_data.append(row)

        if not tidy_data:
            return

        tidy_df = pd.DataFrame(tidy_data)
        tidy_file = os.path.join(self.output_folder, "series_analysis_tidy.csv")
        tidy_df.to_csv(tidy_file, index=False, float_format='%.6f')
        print(f"✨ Created tidy results file for easy analysis: {tidy_file}")

    def _create_comprehensive_output_files(self, batch_results: Dict[str, Any]):
        """
        Creates all output files, including the new tidy format.
        (This is an update to the existing method).
        """
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # --- Existing code ---
        tracking_df = self._create_peak_tracking_table(batch_results)
        if not tracking_df.empty:
            tracking_file = os.path.join(self.output_folder, "comprehensive_peak_tracking.csv")
            tracking_df.to_csv(tracking_file, index=False, float_format='%.6f')
            self._create_intensity_matrix(tracking_df)
            self._create_detection_matrix(tracking_df)
            self._create_summary_statistics_file(batch_results)

        # --- ADD THIS CALL ---
        self._create_tidy_results_file(batch_results)

## gm added

    def _convert_to_integration_format(self, fitted_results: List[Dict]) -> List[Dict]:
        """
        Convert detailed Voigt fitting results to the new standardized peak data format.
        This format is used by all downstream components, including the spectrum viewer.
        FIXED: Uses reference coordinates when Voigt coordinates are missing/zero.
        """
        standardized_results = []
        for i, fit_result in enumerate(fitted_results):
            if not isinstance(fit_result, dict):
                continue

            # Get assignment for reference lookup
            assignment = fit_result.get('assignment', f'Peak_{i+1}')

            # Find reference coordinates from the original peak list
            ref_x, ref_y = 0.0, 0.0
            if self.reference_peaks is not None:
                try:
                    matching_ref = self.reference_peaks[
                        self.reference_peaks['Assignment'] == assignment
                    ]
                    if not matching_ref.empty:
                        ref_x = float(matching_ref['Position_X'].iloc[0])
                        ref_y = float(matching_ref['Position_Y'].iloc[0])
                        print(f"   🎯 Found reference coordinates for {assignment}: ({ref_x:.3f}, {ref_y:.1f})")
                except Exception as e:
                    print(f"   ⚠️ Error getting reference coordinates for {assignment}: {e}")

            # Initialize with reference coordinates as defaults
            peak_data = {
                'assignment': assignment,
                'peak_number': fit_result.get('peak_number', i + 1),
                'ppm_x': ref_x,  # Use reference coordinates as default
                'ppm_y': ref_y,  # Use reference coordinates as default
                'height': 0.0,
                'volume': 0.0,
                'snr': 0.0,
                'quality': 'Failed',
                'r_squared': 0.0,
                'detected': True,   # Always mark as detected since it was processed
                'fitted': False,
                'voigt_fit_data': None,

                # Add legacy format fields for backward compatibility
                'Assignment': assignment,
                'Position_X': ref_x,
                'Position_Y': ref_y,
                'Height': 0.0,
                'Volume': 0.0,
                'SNR': 0.0,
                'Quality': 'Failed',
                'R_Squared': 0.0,
                'Peak_Number': fit_result.get('peak_number', i + 1)
            }

            #success = fit_result.get('success', False)
            # Check multiple success indicators - be more permissive
            success = (fit_result.get('success', False) or
                      fit_result.get('fitted', False) or
                      bool(fit_result.get('x_fit')) or
                      bool(fit_result.get('y_fit')) or
                      fit_result.get('snr', 0) > 0 or
                      fit_result.get('avg_r_squared', 0) > 0)
            if success:
                # Extract data from the complex fit_result dictionary
                x_fit = fit_result.get('x_fit', {})
                y_fit = fit_result.get('y_fit', {})

                # Get fitted coordinates with robust fallback
                fitted_x = None
                fitted_y = None

                # Try multiple ways to get fitted coordinates
                if x_fit and 'center' in x_fit and x_fit['center'] is not None:
                    fitted_x = float(x_fit['center'])
                elif 'ppm_x' in fit_result and fit_result['ppm_x'] is not None:
                    fitted_x = float(fit_result['ppm_x'])
                elif 'peak_position' in fit_result and fit_result['peak_position']:
                    try:
                        fitted_x = float(fit_result['peak_position'][0])
                    except:
                        pass

                if y_fit and 'center' in y_fit and y_fit['center'] is not None:
                    fitted_y = float(y_fit['center'])
                elif 'ppm_y' in fit_result and fit_result['ppm_y'] is not None:
                    fitted_y = float(fit_result['ppm_y'])
                elif 'peak_position' in fit_result and fit_result['peak_position']:
                    try:
                        fitted_y = float(fit_result['peak_position'][1])
                    except:
                        pass

                # Use fitted coordinates if available and valid, otherwise keep reference
                final_x = fitted_x if fitted_x is not None and fitted_x != 0.0 else ref_x
                final_y = fitted_y if fitted_y is not None and fitted_y != 0.0 else ref_y

                # Height extraction - try multiple sources
                height = 0.0
                if x_fit and 'amplitude' in x_fit:
                    height = float(x_fit['amplitude'])
                elif y_fit and 'amplitude' in y_fit:
                    height = float(y_fit['amplitude'])
                elif 'height' in fit_result:
                    height = float(fit_result['height'])
                elif 'peak_intensity' in fit_result:
                    height = float(fit_result['peak_intensity'])

                # SNR extraction - try multiple sources
                snr = 0.0
                if 'snr' in fit_result:
                    snr = float(fit_result['snr'])
                elif 'peak_snr' in fit_result:
                    snr = float(fit_result['peak_snr'])
                elif height > 0 and 'noise_level' in fit_result:
                    snr = height / fit_result['noise_level']

                # R-squared extraction - try multiple sources
                r_squared = 0.0
                if 'avg_r_squared' in fit_result:
                    r_squared = float(fit_result['avg_r_squared'])
                elif 'r_squared' in fit_result:
                    r_squared = float(fit_result['r_squared'])
                elif x_fit and 'r_squared' in x_fit:
                    r_squared = float(x_fit['r_squared'])
                elif y_fit and 'r_squared' in y_fit:
                    r_squared = float(y_fit['r_squared'])

                # Quality determination with fallbacks
                quality = 'Failed'
                if r_squared > 0.90:
                    quality = 'Excellent'
                elif r_squared > 0.80:
                    quality = 'Good'
                elif r_squared > 0.60:
                    quality = 'Fair'
                elif snr > 5.0:
                    quality = 'Poor'
                elif 'fitting_quality' in fit_result:
                    quality = str(fit_result['fitting_quality'])

                print(f"   📊 Peak {assignment}: height={height:.1e}, snr={snr:.1f}, r²={r_squared:.3f}, quality={quality}")

                # Update peak data with fitted values
                peak_data.update({
                    'ppm_x': final_x,
                    'ppm_y': final_y,
                    'height': height,
                    'snr': snr,
                    'quality': quality,
                    'r_squared': r_squared,
                    'fitted': True,
                    'voigt_fit_data': {
                        'x_fit': x_fit,
                        'y_fit': y_fit
                    },

                    # Update legacy fields as well
                    'Position_X': final_x,
                    'Position_Y': final_y,
                    'Height': height,
                    'SNR': snr,
                    'Quality': quality,
                    'R_Squared': r_squared
                })

                # Calculate volume from Voigt parameters
                volume = 0.0
                if x_fit and y_fit and height > 0:
                    x_sigma = x_fit.get('sigma', 0)
                    x_gamma = x_fit.get('gamma', 0)
                    y_sigma = y_fit.get('sigma', 0)
                    y_gamma = y_fit.get('gamma', 0)

                    if (x_sigma > 0 or x_gamma > 0) and (y_sigma > 0 or y_gamma > 0):
                        x_fwhm = 2 * np.sqrt(2 * np.log(2)) * x_sigma + 2 * x_gamma
                        y_fwhm = 2 * np.sqrt(2 * np.log(2)) * y_sigma + 2 * y_gamma
                        volume = abs(height) * x_fwhm * y_fwhm

                peak_data['volume'] = volume
                peak_data['Volume'] = volume  # Legacy field

                print(f"   ✅ Peak {assignment}: fitted=({final_x:.3f}, {final_y:.1f}), ref=({ref_x:.3f}, {ref_y:.1f})")

            else:
                print(f"   ❌ Peak {assignment}: fit failed, using reference coordinates ({ref_x:.3f}, {ref_y:.1f})")

            standardized_results.append(peak_data)

        return standardized_results

##

##

    def _initialize_comprehensive_batch_results(self, nmr_files: List[str],
                                              output_folder: str, peak_source_mode: str) -> Dict[str, Any]:
        """Initialize comprehensive batch results structure"""
        return {
            'metadata': {
                'processing_started': datetime.now().isoformat(),
                'peak_source_mode': peak_source_mode,
                'total_spectra': len(nmr_files),
                'output_folder': output_folder,
                'reference_peaks_count': len(self.reference_peaks)
            },
            'results': {},
            'summary': {}
        }

    def _create_comprehensive_output_files(self, batch_results: Dict[str, Any]):
        """Create all output files needed for downstream GUI components"""
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Create tracking DataFrame
        tracking_df = self._create_peak_tracking_table(batch_results)

        if not tracking_df.empty:
            # Main comprehensive tracking file
            tracking_file = os.path.join(self.output_folder, "comprehensive_peak_tracking.csv")
            tracking_df.to_csv(tracking_file, index=False, float_format='%.6f')

            # Create intensity matrix
            self._create_intensity_matrix(tracking_df)

            # Create detection matrix
            self._create_detection_matrix(tracking_df)

            # Create summary statistics file
            self._create_summary_statistics_file(batch_results)

    def _create_peak_tracking_table(self, batch_results: Dict[str, Any]) -> pd.DataFrame:
        """Create comprehensive peak tracking table"""
        all_peaks_data = []

        for peak_idx, ref_peak in self.reference_peaks.iterrows():
            peak_row = {
                'Peak_Number': peak_idx + 1,
                'Assignment': ref_peak.get('Assignment', f'Peak_{peak_idx + 1}'),
                'Reference_X': ref_peak['Position_X'],
                'Reference_Y': ref_peak['Position_Y']
            }

            # Add data for each spectrum
            for spectrum_name, result_data in batch_results['results'].items():
                spectrum_key = os.path.splitext(spectrum_name)[0]
                found_peak = False

                if result_data.get('success', False):
                    integration_results = result_data.get('integration_results', [])
                    for peak in integration_results:
                        if peak.get('Assignment') == peak_row['Assignment']:
                            peak_row[f'{spectrum_key}_Detected'] = True
                            peak_row[f'{spectrum_key}_Height'] = peak.get('Height', 0.0)
                            peak_row[f'{spectrum_key}_Volume'] = peak.get('Volume', 0.0)
                            peak_row[f'{spectrum_key}_SNR'] = peak.get('SNR', 0.0)
                            peak_row[f'{spectrum_key}_Quality'] = peak.get('Quality', 'Poor')
                            peak_row[f'{spectrum_key}_Position_X'] = peak.get('Position_X', 0.0)
                            peak_row[f'{spectrum_key}_Position_Y'] = peak.get('Position_Y', 0.0)
                            peak_row[f'{spectrum_key}_R_Squared'] = peak.get('R_Squared', 0.0)
                            found_peak = True
                            break

                if not found_peak:
                    peak_row[f'{spectrum_key}_Detected'] = False
                    peak_row[f'{spectrum_key}_Height'] = 0.0
                    peak_row[f'{spectrum_key}_Volume'] = 0.0
                    peak_row[f'{spectrum_key}_SNR'] = 0.0
                    peak_row[f'{spectrum_key}_Quality'] = 'Not Detected'
                    peak_row[f'{spectrum_key}_Position_X'] = 0.0
                    peak_row[f'{spectrum_key}_Position_Y'] = 0.0
                    peak_row[f'{spectrum_key}_R_Squared'] = 0.0

            all_peaks_data.append(peak_row)

        return pd.DataFrame(all_peaks_data)

    def _create_intensity_matrix(self, tracking_df: pd.DataFrame):
        """Create peak intensity matrix file"""
        intensity_data = tracking_df[['Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y']].copy()

        for col in tracking_df.columns:
            if col.endswith('_Height'):
                spectrum_name = col.replace('_Height', '')
                intensity_data[spectrum_name] = tracking_df[col]

        intensity_file = os.path.join(self.output_folder, "peak_intensity_matrix.csv")
        intensity_data.to_csv(intensity_file, index=False, float_format='%.6f')

    def _create_detection_matrix(self, tracking_df: pd.DataFrame):
        """Create peak detection matrix file"""
        detection_data = tracking_df[['Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y']].copy()

        for col in tracking_df.columns:
            if col.endswith('_Detected'):
                spectrum_name = col.replace('_Detected', '')
                detection_data[spectrum_name] = tracking_df[col].astype(int)

        detection_file = os.path.join(self.output_folder, "peak_detection_matrix.csv")
        detection_data.to_csv(detection_file, index=False)

    def _create_summary_statistics_file(self, batch_results: Dict[str, Any]):
        """Create summary statistics file"""
        summary_data = []

        for spectrum_name, result_data in batch_results['results'].items():
            summary_row = {
                'Spectrum': spectrum_name,
                'Success': result_data.get('success', False),
                'Total_Peaks': result_data.get('total_peaks', 0),
                'Successful_Fits': result_data.get('successful_fits', 0),
                'Success_Rate': result_data.get('success_rate', 0.0)
            }
            summary_data.append(summary_row)

        summary_df = pd.DataFrame(summary_data)
        summary_file = os.path.join(self.output_folder, "processing_summary.csv")
        summary_df.to_csv(summary_file, index=False, float_format='%.2f')

    def cancel_processing(self):
        """Cancel current processing operation"""
        self.processing_active = False
        print("⏹️ Multi-spectrum processing cancelled")
