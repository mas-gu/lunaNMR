#!/usr/bin/env python3
"""
Enhanced Series NMR Integration Script - Reference-Based Peak Propagation

This script implements advanced batch processing of multiple NMR spectra using
enhanced reference-based peak detection with intelligent peak propagation.

Key enhancements over original version:
1. Reference-guided peak detection with configurable search windows
2. Intelligent peak propagation across series
3. User-configurable noise regions and threshold controls
4. Enhanced quality tracking and statistics
5. Improved error handling and validation
6. Advanced series analysis capabilities

Author: Guillaume Mas
Date: 2025
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import glob
from typing import List, Dict, Tuple, Optional
import nmrglue as ng
import time
from datetime import datetime

# Import the enhanced integrator
from lunaNMR.integrators.inplace_advanced_nmr_integrator import InPlaceAdvancedNMRIntegrator

# ========================================================================
# USER CONFIGURATION - MODIFY THESE PARAMETERS AS NEEDED
# ========================================================================

# Input settings
PEAK_LIST_FILE = "peak_list/KRASB/WT/600_KRASB_GTP.txt"    # Initial peak list for reference
NMR_DATA_DIRECTORY = "data/KRASB/WT/WT_600/600_T1_WT/"                # Directory containing .ft files
REFERENCE_SPECTRUM = "600_T1_KRASB_WT_0o0.ft"                          # Reference spectrum filename

# File pattern matching
FILE_PATTERN = "*.ft"                                 # Pattern to match NMR files
EXCLUDE_FILES = []                                    # List of files to exclude

# Output settings
OUTPUT_DIR = "series_results/KRASB/WT/WT_600/600_T1_WT"     # Output directory for series results
CREATE_INDIVIDUAL_PLOTS = True                       # Create plot for each spectrum
CREATE_PDF_PLOTS = True                              # Save plots as PDF format
CREATE_SERIES_SUMMARY = True                         # Create series summary plots
SAVE_REFERENCE_PEAKS = True                          # Save refined peak positions from reference

# Enhanced Processing parameters
USE_REFERENCE_DETECTION = True                       # Enable reference-based detection
SEARCH_WINDOW_X_PPM = 0.2                           # ±X ppm search window (1H dimension)
SEARCH_WINDOW_Y_PPM = 3.0                           # ±Y ppm search window (15N/13C dimension)
PROPAGATE_SUCCESSFUL_DETECTIONS = True              # Propagate successful detections forward
ADAPTIVE_SEARCH_WINDOW = True                       # Adaptively adjust search windows

# Reference spectrum parameters
MIN_REFERENCE_SNR = 3.0                             # Minimum SNR for reference peaks
MIN_REFERENCE_QUALITY = 0.5                         # Minimum quality for reference peaks
REFERENCE_THRESHOLD_MULTIPLIER = 3.0                # Threshold multiplier for reference

# Series propagation parameters
PEAK_POSITION_TOLERANCE = 0.05                      # Tolerance for position updates (ppm)
SEARCH_RADIUS_SERIES = 0.3                          # Search radius for series spectra (ppm)
MIN_SERIES_DETECTION = 3                            # Minimum spectra where peak must be found
PROPAGATION_CONFIDENCE_THRESHOLD = 0.7              # Confidence threshold for propagation

# Enhanced Quality control
FAIL_THRESHOLD = 0.3                                # Fraction of peaks that must be detected
OUTLIER_Z_SCORE = 3.0                               # Z-score threshold for outlier detection
QUALITY_ASSESSMENT_WINDOW = 5                       # Window for quality assessment
MIN_DETECTION_RATE = 0.4                           # Minimum detection rate to continue

# Noise handling for series
USE_ADAPTIVE_NOISE_ESTIMATION = True                # Adapt noise estimation per spectrum
NOISE_CONSISTENCY_CHECK = True                      # Check noise level consistency
MAX_NOISE_VARIATION = 2.0                          # Maximum allowed noise variation factor

# Visualization
SERIES_COLORMAP = 'tab10'                           # Colormap for series visualization
SHOW_SERIES_TRENDS = True                           # Show integration trends across series
FIGURE_SIZE_SERIES = (16, 10)                      # Figure size for series plots
COLORMAP = 'viridis'                                # Colormap for contour plots
CONTOUR_LEVELS = 20                                 # Number of contour levels
DPI = 300                                           # Resolution for saved plots
SHOW_PEAK_LABELS = True                             # Show peak assignment labels
LABEL_FONT_SIZE = 8                                 # Font size for peak labels

# ========================================================================

class InPlaceSeriesNMRIntegrator(InPlaceAdvancedNMRIntegrator):
    """Enhanced Series NMR integrator with reference-based peak propagation"""

    def __init__(self):
        """Initialize the Enhanced Series NMR Integrator"""
        super().__init__()

        # Series-specific attributes
        self.reference_peaks = None
        self.reference_spectrum_file = None
        self.series_results = []
        self.spectrum_files = []
        self.failed_spectra = []
        self.series_summary = None

        # Enhanced tracking
        self.peak_detection_rates = {}
        self.peak_quality_trends = {}
        self.peak_propagation_history = {}
        self.noise_level_history = {}
        self.detection_confidence_scores = {}

        # Series processing state
        self.current_spectrum_index = 0
        self.total_spectra = 0
        self.processing_start_time = None

        # Enhanced configuration
        self.series_config = {
            'use_reference_detection': USE_REFERENCE_DETECTION,
            'search_window_x_ppm': SEARCH_WINDOW_X_PPM,
            'search_window_y_ppm': SEARCH_WINDOW_Y_PPM,
            'propagate_successful_detections': PROPAGATE_SUCCESSFUL_DETECTIONS,
            'adaptive_search_window': ADAPTIVE_SEARCH_WINDOW,
            'min_reference_snr': MIN_REFERENCE_SNR,
            'min_reference_quality': MIN_REFERENCE_QUALITY,
            'reference_threshold_multiplier': REFERENCE_THRESHOLD_MULTIPLIER
        }

    def find_nmr_files(self, directory: str, pattern: str = "*.ft") -> List[str]:
        """Find all NMR files in directory matching pattern"""
        search_path = os.path.join(directory, pattern)
        files = glob.glob(search_path)

        # Filter out excluded files
        if EXCLUDE_FILES:
            files = [f for f in files if os.path.basename(f) not in EXCLUDE_FILES]

        # Sort files naturally
        files.sort()

        print(f"Found {len(files)} NMR files in {directory}")
        return files

    def setup_reference_spectrum(self, spectrum_file: str, peak_list_file: str) -> bool:
        """Setup reference spectrum with enhanced detection"""
        print(f"Setting up reference spectrum: {spectrum_file}")

        # Load reference data
        success = self.load_data(peak_list_file, spectrum_file)
        if not success:
            print(f"Failed to load reference spectrum: {spectrum_file}")
            return False

        # Configure enhanced parameters for reference
        self.set_search_window(
            self.series_config['search_window_x_ppm'],
            self.series_config['search_window_y_ppm']
        )
        self.set_threshold_multiplier(self.series_config['reference_threshold_multiplier'])
        self.use_reference_detection = self.series_config['use_reference_detection']

        # Perform enhanced detection
        print("Performing enhanced reference detection...")
        detected_peaks = self.detect_peaks_professional()

        if not detected_peaks:
            print("No peaks detected in reference spectrum")
            return False

        # Filter reference peaks by quality
        quality_filtered_peaks = []
        for peak in detected_peaks:
            if (peak.get('snr', 0) >= self.series_config['min_reference_snr'] and
                peak.get('detected', False)):
                quality_filtered_peaks.append(peak)

        print(f"Reference spectrum analysis:")
        print(f"  Total detected: {len(detected_peaks)}")
        print(f"  Quality filtered: {len(quality_filtered_peaks)}")

        if len(quality_filtered_peaks) < len(self.peak_list) * 0.3:
            print("Warning: Low detection rate in reference spectrum")

        # Store reference peaks for series propagation
        self.reference_peaks = quality_filtered_peaks
        self.reference_spectrum_file = spectrum_file

        # Store reference noise level for consistency checking
        self.noise_level_history[spectrum_file] = self.noise_level

        # Perform integration on reference
        integration_results = self.integrate_peaks()

        # Save reference results
        if SAVE_REFERENCE_PEAKS:
            self._save_reference_results()

        print(f"Reference spectrum setup completed: {len(self.reference_peaks)} reference peaks")
        return True

    def _save_reference_results(self):
        """Save reference spectrum results"""
        os.makedirs(OUTPUT_DIR, exist_ok=True)

        # Save reference peak positions
        if self.reference_peaks:
            ref_df = pd.DataFrame([
                {
                    'Assignment': peak['assignment'],
                    'Reference_X': peak['ppm_x'],
                    'Reference_Y': peak['ppm_y'],
                    'Reference_SNR': peak.get('snr', 0),
                    'Detection_Quality': peak.get('detection_quality', 'Unknown'),
                    'Quality_Score': peak.get('quality_score', 0)
                }
                for peak in self.reference_peaks
            ])

            ref_file = os.path.join(OUTPUT_DIR, "reference_peaks.csv")
            ref_df.to_csv(ref_file, index=False, float_format='%.6f')
            print(f"Reference peaks saved to {ref_file}")

        # Save reference integration results
        if hasattr(self, 'integration_results') and self.integration_results:
            ref_integration_file = os.path.join(OUTPUT_DIR, "reference_integration.csv")
            ref_df = pd.DataFrame(self.integration_results)
            ref_df.to_csv(ref_integration_file, index=False, float_format='%.6f')
            print(f"Reference integration saved to {ref_integration_file}")

    def process_spectrum_in_series(self, spectrum_file: str, spectrum_name: str) -> Dict:
        """Process individual spectrum using reference-based detection with propagation"""
        print(f"Processing spectrum: {spectrum_name}")

        # Store current peak list temporarily
        temp_peak_list = self.peak_list

        try:
            # Load only NMR data
            dic, data = ng.pipe.read(spectrum_file)
            self.nmr_data = data
            self.nmr_dict = dic
            self._calculate_ppm_axes()
            self._estimate_noise_level()

            # Check noise consistency
            if NOISE_CONSISTENCY_CHECK and self.reference_spectrum_file in self.noise_level_history:
                ref_noise = self.noise_level_history[self.reference_spectrum_file]
                noise_ratio = self.noise_level / ref_noise
                if noise_ratio > MAX_NOISE_VARIATION or noise_ratio < (1/MAX_NOISE_VARIATION):
                    print(f"  Warning: Noise level variation detected (ratio: {noise_ratio:.2f})")

            # Store noise level
            self.noise_level_history[spectrum_name] = self.noise_level

        except Exception as e:
            print(f"  Failed to load spectrum: {spectrum_name} - {e}")
            return {'file': spectrum_file, 'status': 'failed', 'error': 'load_failed'}

        # Create propagated peak list from reference detections
        propagated_peak_list = self._create_propagated_peak_list()
        if propagated_peak_list is None:
            return {'file': spectrum_file, 'status': 'failed', 'error': 'no_reference_peaks'}

        # Use propagated peak list for this spectrum
        self.peak_list = propagated_peak_list

        # Perform enhanced detection with propagated positions
        try:
            detected_peaks = self.detect_peaks_professional()

            if not detected_peaks:
                print(f"  No peaks detected in {spectrum_name}")
                return {'file': spectrum_file, 'status': 'failed', 'error': 'no_detection'}

            # Perform integration
            integration_results = self.integrate_peaks()

            # Calculate detection statistics
            detection_stats = self.get_detection_statistics()
            detected_count = detection_stats['detected_peaks']
            total_count = detection_stats['total_peaks']
            detection_rate = detection_stats['detection_rate'] / 100.0

            print(f"  Detection: {detected_count}/{total_count} peaks ({detection_rate:.1%})")

            # Check if spectrum meets quality threshold
            if detection_rate < FAIL_THRESHOLD:
                print(f"  Warning: Low detection rate in {spectrum_name} ({detection_rate:.1%})")
                status = 'poor_quality'
            else:
                status = 'success'

            # Update propagation history
            self._update_propagation_history(spectrum_name, detected_peaks)

            # Store results
            result = {
                'file': spectrum_file,
                'name': spectrum_name,
                'status': status,
                'detected_peaks': detected_count,
                'total_peaks': total_count,
                'detection_rate': detection_rate,
                'noise_level': self.noise_level,
                'integration_results': integration_results,
                'fitted_peaks': detected_peaks.copy(),
                'detection_statistics': detection_stats
            }

            return result

        except Exception as e:
            print(f"  Processing failed for {spectrum_name}: {e}")
            return {'file': spectrum_file, 'status': 'failed', 'error': str(e)}

        finally:
            # Restore original peak list
            self.peak_list = temp_peak_list

    def _create_propagated_peak_list(self) -> Optional[pd.DataFrame]:
        """Create peak list using propagated positions from reference detections"""
        if not self.reference_peaks:
            return None

        propagated_data = []
        for ref_peak in self.reference_peaks:
            # Use detected position from reference as new reference
            propagated_data.append({
                'Assignment': ref_peak['assignment'],
                'Position_X': ref_peak['ppm_x'],  # Use detected position
                'Position_Y': ref_peak['ppm_y']   # Use detected position
            })

        return pd.DataFrame(propagated_data)

    def _update_propagation_history(self, spectrum_name: str, detected_peaks: List[Dict]):
        """Update peak propagation history for analysis"""
        if spectrum_name not in self.peak_propagation_history:
            self.peak_propagation_history[spectrum_name] = {}

        for peak in detected_peaks:
            assignment = peak['assignment']
            if assignment not in self.peak_propagation_history[spectrum_name]:
                self.peak_propagation_history[spectrum_name][assignment] = {
                    'detected': peak.get('detected', False),
                    'position_x': peak['ppm_x'],
                    'position_y': peak['ppm_y'],
                    'snr': peak.get('snr', 0),
                    'quality': peak.get('detection_quality', 'Unknown')
                }

    def process_series(self, data_directory: str, peak_list_file: str,
                      reference_spectrum: Optional[str] = None) -> bool:
        """Process entire series with enhanced reference-based detection"""
        print("Enhanced Series NMR Integration - Reference-Based Detection")
        print("=" * 70)

        import time
        self.processing_start_time = time.time()

        # Find all NMR files
        self.spectrum_files = self.find_nmr_files(data_directory, FILE_PATTERN)
        if not self.spectrum_files:
            print(f"No NMR files found in {data_directory}")
            return False

        self.total_spectra = len(self.spectrum_files)

        # Determine reference spectrum
        if reference_spectrum:
            ref_file = os.path.join(data_directory, reference_spectrum)
            if not os.path.exists(ref_file):
                print(f"Reference spectrum not found: {ref_file}")
                return False
        else:
            ref_file = self.spectrum_files[0]
            print(f"Using first spectrum as reference: {os.path.basename(ref_file)}")

        # Setup reference spectrum
        success = self.setup_reference_spectrum(ref_file, peak_list_file)
        if not success:
            print("Failed to setup reference spectrum")
            return False

        # Process each spectrum in series
        print(f"\nProcessing {self.total_spectra} spectra...")

        for i, spectrum_file in enumerate(self.spectrum_files):
            self.current_spectrum_index = i
            spectrum_name = os.path.basename(spectrum_file)

            print(f"\n[{i+1}/{self.total_spectra}] Processing: {spectrum_name}")

            # Skip reference spectrum (already processed)
            if spectrum_file == ref_file:
                print("  Skipping reference spectrum")
                continue

            # Process spectrum
            result = self.process_spectrum_in_series(spectrum_file, spectrum_name)

            if result['status'] == 'failed':
                self.failed_spectra.append(result)
                print(f"  ❌ Failed: {result.get('error', 'Unknown error')}")
            elif result['status'] == 'poor_quality':
                self.series_results.append(result)
                print(f"  ⚠️ Poor quality: {result['detection_rate']:.1%} detection rate")
            else:
                self.series_results.append(result)
                print(f"  ✅ Success: {result['detection_rate']:.1%} detection rate")

            # Check if we should continue based on recent performance
            if i > 5 and self._should_stop_processing():
                print("  Stopping processing due to poor detection rates")
                break

        # Generate series analysis
        self._analyze_series_results()

        # Create output directory
        os.makedirs(OUTPUT_DIR, exist_ok=True)

        # Export results
        self._export_series_results()

        # Create visualizations
        if CREATE_SERIES_SUMMARY:
            self._create_series_summary_plots()

        # Print final summary
        self._print_final_summary()

        return True

    def _should_stop_processing(self) -> bool:
        """Determine if processing should stop based on recent performance"""
        if len(self.series_results) < 3:
            return False

        # Check recent detection rates
        recent_results = self.series_results[-3:]
        recent_detection_rates = [r['detection_rate'] for r in recent_results]
        avg_recent_rate = np.mean(recent_detection_rates)

        if avg_recent_rate < MIN_DETECTION_RATE:
            print(f"  Recent average detection rate too low: {avg_recent_rate:.1%}")
            return True

        return False

    def _analyze_series_results(self):
        """Analyze results across the entire series"""
        print("\nAnalyzing series results...")

        if not self.series_results:
            print("No successful results to analyze")
            return

        # Combine all integration results
        all_integration_results = []
        for result in self.series_results:
            if 'integration_results' in result:
                for integration in result['integration_results']:
                    integration['Spectrum'] = result['name']
                    integration['Detection_Rate'] = result['detection_rate']
                    all_integration_results.append(integration)

        if not all_integration_results:
            print("No integration results to analyze")
            return

        # Convert to DataFrame for analysis
        df = pd.DataFrame(all_integration_results)

        # Analyze peak detection across series
        peak_detection_count = {}
        for assignment in df['Assignment'].unique():
            peak_data = df[df['Assignment'] == assignment]
            detected_count = len(peak_data[peak_data['Integration_Method'] != 'Reference'])
            total_spectra = len(self.series_results)
            detection_rate = detected_count / total_spectra if total_spectra > 0 else 0

            peak_detection_count[assignment] = {
                'detected_count': detected_count,
                'total_spectra': total_spectra,
                'detection_rate': detection_rate,
                'mean_snr': peak_data['SNR'].mean(),
                'std_snr': peak_data['SNR'].std()
            }

        self.peak_detection_rates = peak_detection_count

        # Create series summary
        self.series_summary = {
            'total_spectra': len(self.spectrum_files),
            'successful_spectra': len(self.series_results),
            'failed_spectra': len(self.failed_spectra),
            'mean_detection_rate': np.mean([r['detection_rate'] for r in self.series_results]),
            'std_detection_rate': np.std([r['detection_rate'] for r in self.series_results]),
            'peak_statistics': peak_detection_count
        }

        print(f"Series analysis completed:")
        print(f"  Successful spectra: {self.series_summary['successful_spectra']}")
        print(f"  Failed spectra: {self.series_summary['failed_spectra']}")
        print(f"  Mean detection rate: {self.series_summary['mean_detection_rate']:.1%}")


    def _export_series_results(self):
        """Export comprehensive series results"""
        print("Exporting enhanced series results...")

        # Export individual spectrum results
        if self.series_results:
            # Create summary DataFrame
            summary_data = []
            for result in self.series_results:
                summary_data.append({
                    'Spectrum': result['name'],
                    'Status': result['status'],
                    'Detected_Peaks': result['detected_peaks'],
                    'Total_Peaks': result['total_peaks'],
                    'Detection_Rate': result['detection_rate'],
                    'Noise_Level': result['noise_level']
                })

            summary_df = pd.DataFrame(summary_data)
            summary_file = os.path.join(OUTPUT_DIR, "enhanced_series_summary.csv")
            summary_df.to_csv(summary_file, index=False, float_format='%.6f')
            print(f"Series summary saved to {summary_file}")

            # Export consolidated integration results
            all_integrations = []
            for result in self.series_results:
                if 'integration_results' in result:
                    for integration in result['integration_results']:
                        integration['Spectrum'] = result['name']
                        integration['Detection_Rate'] = result['detection_rate']
                        all_integrations.append(integration)

            if all_integrations:
                integration_df = pd.DataFrame(all_integrations)
                integration_file = os.path.join(OUTPUT_DIR, "enhanced_series_integrations.csv")
                integration_df.to_csv(integration_file, index=False, float_format='%.6f')
                print(f"Consolidated integrations saved to {integration_file}")

        # Export peak detection statistics
        if self.peak_detection_rates:
            peak_stats_data = []
            for assignment, stats in self.peak_detection_rates.items():
                peak_stats_data.append({
                    'Assignment': assignment,
                    'Detected_Count': stats['detected_count'],
                    'Total_Spectra': stats['total_spectra'],
                    'Detection_Rate': stats['detection_rate'],
                    'Mean_SNR': stats['mean_snr'],
                    'Std_SNR': stats['std_snr']
                })

            peak_stats_df = pd.DataFrame(peak_stats_data)
            peak_stats_file = os.path.join(OUTPUT_DIR, "enhanced_peak_detection_statistics.csv")
            peak_stats_df.to_csv(peak_stats_file, index=False, float_format='%.6f')
            print(f"Peak detection statistics saved to {peak_stats_file}")

        # Export failed spectra information
        if self.failed_spectra:
            failed_df = pd.DataFrame(self.failed_spectra)
            failed_file = os.path.join(OUTPUT_DIR, "failed_spectra.csv")
            failed_df.to_csv(failed_file, index=False)
            print(f"Failed spectra log saved to {failed_file}")


    def _create_series_summary_plots(self):
        """Create enhanced series summary visualizations"""
        print("Creating enhanced series summary plots...")

        if not self.series_results:
            print("No results available for plotting")
            return

        # Create multi-panel summary figure
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=FIGURE_SIZE_SERIES)

        # Plot 1: Detection rates across series
        spectrum_names = [r['name'] for r in self.series_results]
        detection_rates = [r['detection_rate'] for r in self.series_results]

        ax1.plot(range(len(detection_rates)), detection_rates, 'bo-', markersize=4)
        ax1.axhline(y=FAIL_THRESHOLD, color='r', linestyle='--', alpha=0.7, label=f'Fail Threshold ({FAIL_THRESHOLD:.1%})')
        ax1.set_xlabel('Spectrum Index')
        ax1.set_ylabel('Detection Rate')
        ax1.set_title('Enhanced Detection Rates Across Series')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # Plot 2: Peak detection statistics
        if self.peak_detection_rates:
            assignments = list(self.peak_detection_rates.keys())[:20]  # Top 20 peaks
            rates = [self.peak_detection_rates[a]['detection_rate'] for a in assignments]

            ax2.barh(range(len(assignments)), rates, color='skyblue', alpha=0.7)
            ax2.set_yticks(range(len(assignments)))
            ax2.set_yticklabels(assignments, fontsize=8)
            ax2.set_xlabel('Detection Rate')
            ax2.set_title('Peak Detection Rates (Top 20)')
            ax2.grid(True, alpha=0.3)

        # Plot 3: Noise level consistency
        noise_levels = [r['noise_level'] for r in self.series_results]
        ax3.plot(range(len(noise_levels)), noise_levels, 'go-', markersize=4)
        ax3.set_xlabel('Spectrum Index')
        ax3.set_ylabel('Noise Level')
        ax3.set_title('Noise Level Consistency')
        ax3.grid(True, alpha=0.3)

        # Plot 4: Quality distribution
        statuses = [r['status'] for r in self.series_results]
        status_counts = pd.Series(statuses).value_counts()

        ax4.pie(status_counts.values, labels=status_counts.index, autopct='%1.1f%%', startangle=90)
        ax4.set_title('Processing Quality Distribution')

        plt.tight_layout()

        # Save plot
        plot_file = os.path.join(OUTPUT_DIR, "enhanced_series_summary.png")
        plt.savefig(plot_file, dpi=DPI, bbox_inches='tight')

        if CREATE_PDF_PLOTS:
            pdf_file = os.path.join(OUTPUT_DIR, "enhanced_series_summary.pdf")
            plt.savefig(pdf_file, bbox_inches='tight')

        plt.close()
        print(f"Enhanced series summary plot saved to {plot_file}")

    def _print_final_summary(self):
        """Print comprehensive final summary"""
        elapsed_time = time.time() - self.processing_start_time if self.processing_start_time else 0

        print("\n" + "=" * 70)
        print("ENHANCED SERIES PROCESSING SUMMARY")
        print("=" * 70)

        if self.series_summary:
            print(f"📊 Total Spectra: {self.series_summary['total_spectra']}")
            print(f"✅ Successfully Processed: {self.series_summary['successful_spectra']}")
            print(f"❌ Failed: {self.series_summary['failed_spectra']}")
            print(f"📈 Mean Detection Rate: {self.series_summary['mean_detection_rate']:.1%} ± {self.series_summary['std_detection_rate']:.1%}")

        if self.reference_peaks:
            print(f"🎯 Reference Peaks: {len(self.reference_peaks)}")

        if self.peak_detection_rates:
            high_quality_peaks = sum(1 for stats in self.peak_detection_rates.values()
                                   if stats['detection_rate'] >= 0.8)
            reliable_peaks = sum(1 for stats in self.peak_detection_rates.values()
                               if stats['detection_rate'] >= 0.6)

            print(f"🏆 High-quality peaks (≥80% detection): {high_quality_peaks}")
            print(f"✔️ Reliable peaks (≥60% detection): {reliable_peaks}")
            print(f"⚠️ Poor peaks (<60% detection): {len(self.peak_detection_rates) - reliable_peaks}")

        print(f"⏱️ Processing time: {elapsed_time:.1f} seconds")
        print(f"📁 Results saved to: {OUTPUT_DIR}")

        print("\n🔬 Enhanced Reference-Based Detection completed successfully!")
        print("=" * 70)

def main():
    """Main execution function for enhanced series processing"""
    print("Enhanced NMR Series Integration - Reference-Based Peak Propagation")
    print("=" * 70)

    # Initialize enhanced series integrator
    integrator = InPlaceSeriesNMRIntegrator()

    # Process the series
    success = integrator.process_series(
        data_directory=NMR_DATA_DIRECTORY,
        peak_list_file=PEAK_LIST_FILE,
        reference_spectrum=REFERENCE_SPECTRUM
    )

    if success:
        print("\n🎉 Enhanced series processing completed successfully!")
    else:
        print("\n❌ Enhanced series processing failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()
