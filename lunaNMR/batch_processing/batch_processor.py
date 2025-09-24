#!/usr/bin/env python3
"""
Core Batch Processing Engine for LunaNMR

This module provides the core functionality for batch processing multiple NMR spectra
for automated analysis and ML training data generation. It operates independently
from the GUI interface and integrates seamlessly with existing ML data collection.

Features:
- Automated spectrum loading and processing
- Integration with existing ML training data collection
- Robust error handling with individual spectrum failure recovery
- Progress tracking and detailed logging
- Support for multiple nucleus types (1H, 15N, 13C)
"""

import os
import sys
import time
import logging
from pathlib import Path
from typing import List, Dict, Optional, Union, Any
import numpy as np
from datetime import datetime
import traceback

# Import lunaNMR core components (READ-ONLY usage)
try:
    from ..core.core_integrator import EnhancedVoigtIntegrator
    from ..core.enhanced_voigt_fitter import EnhancedVoigtFitter
    CORE_AVAILABLE = True
except ImportError:
    CORE_AVAILABLE = False

# Note: Processor imports will be done dynamically in methods to avoid import issues

class BatchProcessor:
    """
    Core batch processing engine for automated NMR spectrum analysis.

    This class handles the automated processing of multiple NMR spectra,
    including S/N detection, peak fitting, and ML training data collection.
    It operates completely independently from the GUI interface.
    """

    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize the batch processor.

        Args:
            config: Optional configuration dictionary with processing parameters
        """
        # Set up logging
        self.logger = logging.getLogger(__name__)
        self._setup_logging()

        # Default configuration
        self.config = {
            'sn_thresholds': {
                '15N1H': 2.5,
                '13C1H': 2.0,
                '1H': 2.0,
                'default': 2.2
            },
            'expected_peaks': {
                '15N1H': 100,
                '13C1H': 80,
                '1H': 60,
                'default': 150
            },
            'quality_threshold': 0.6,
            'max_fitting_attempts': 3,
            'skip_on_error': True,
            'detailed_logging': True,
            'progress_interval': 5,
            # Add missing iteration parameters
            'max_iterations': 1000,                    # Curve fitting max iterations
            'max_optimization_iterations': 50,         # Parameter optimization iterations
            'baseline_max_iter': 50,                   # Baseline correction iterations
            'convergence_tolerance': 1e-8,             # Fitting convergence tolerance
            'processing_parameters': {
                'noise_threshold': 3.0,
                'search_window_x': 0.2,
                'search_window_y': 3.0,
                'fitting_window_x': 0.3,
                'fitting_window_y': 8.0,
                'use_reference_detection': True
            }
        }

        # Default file handling configuration (matching ConfigManager structure)
        self.default_file_handling = {
            'extensions': ['.ft2', '.ft', '.pipe', '.ucsf', '.nmrpipe'],
            'recursive_search': True,
            'pattern_filters': [],
            'exclude_patterns': ['*backup*', '*temp*']
        }

        # Update with provided config
        if config:
            # Handle nested config structure from aggressive_config.json
            if 'processing' in config:
                processing_config = config['processing']
                # Extract processing parameters to flat structure
                if 'sn_thresholds' in processing_config:
                    self.config['sn_thresholds'] = processing_config['sn_thresholds']
                if 'expected_peaks' in processing_config:
                    self.config['expected_peaks'] = processing_config['expected_peaks']
                if 'quality_threshold' in processing_config:
                    self.config['quality_threshold'] = processing_config['quality_threshold']
                if 'max_fitting_attempts' in processing_config:
                    self.config['max_fitting_attempts'] = processing_config['max_fitting_attempts']
                if 'skip_on_error' in processing_config:
                    self.config['skip_on_error'] = processing_config['skip_on_error']
                if 'max_iterations' in processing_config:
                    self.config['max_iterations'] = processing_config['max_iterations']
                if 'max_iter' in processing_config:
                    self.config['max_iterations'] = processing_config['max_iter']
                if 'max_optimization_iterations' in processing_config:
                    self.config['max_optimization_iterations'] = processing_config['max_optimization_iterations']
                if 'baseline_max_iter' in processing_config:
                    self.config['baseline_max_iter'] = processing_config['baseline_max_iter']
                if 'convergence_tolerance' in processing_config:
                    self.config['convergence_tolerance'] = processing_config['convergence_tolerance']

            if 'processing_parameters' in config:
                self.config['processing_parameters'].update(config['processing_parameters'])

            # Handle other top-level config sections
            for key in ['file_handling', 'optimization', 'logging', 'output']:
                if key in config:
                    self.config[key] = config[key]

        # Ensure file_handling section exists with defaults
        if 'file_handling' not in self.config:
            self.config['file_handling'] = self.default_file_handling
        else:
            # Merge file_handling defaults with provided config
            for key, default_value in self.default_file_handling.items():
                if key not in self.config['file_handling']:
                    self.config['file_handling'][key] = default_value

        # Processing statistics
        self.stats = {
            'total_files': 0,
            'processed_files': 0,
            'failed_files': 0,
            'total_peaks_detected': 0,
            'total_peaks_fitted': 0,
            'total_ml_samples': 0,
            'processing_start_time': None,
            'processing_end_time': None,
            'failed_files_list': []
        }

        # Initialize core integrator
        if CORE_AVAILABLE:
            self.integrator = EnhancedVoigtIntegrator()

            # Configure integrator with batch config parameters
            max_iter = self.config.get('max_iterations', 1000)
            self.integrator.fitting_parameters['max_iterations'] = max_iter
            if hasattr(self.integrator, 'enhanced_fitter') and self.integrator.enhanced_fitter:
                self.integrator.enhanced_fitter.fitting_parameters['max_iterations'] = max_iter
                self.logger.info(f"Configured max_iterations: {max_iter}")

            proc_param_config = self.config.get('processing_parameters', {})

            search_window_x = proc_param_config.get('search_window_x')
            search_window_y = proc_param_config.get('search_window_y')
            if search_window_x and search_window_y and hasattr(self.integrator, 'set_search_window'):
                try:
                    self.integrator.set_search_window(float(search_window_x), float(search_window_y))
                    self.logger.info(
                        f"Configured search window: ±{search_window_x} ppm (1H) × ±{search_window_y} ppm (15N/13C)"
                    )
                except Exception as exc:
                    self.logger.warning(f"Unable to apply search window overrides: {exc}")

            fitting_window_x = proc_param_config.get('fitting_window_x')
            if fitting_window_x is not None:
                try:
                    window_x_value = float(fitting_window_x)
                    self.integrator.fitting_parameters['fitting_window_x'] = window_x_value
                    if hasattr(self.integrator, 'enhanced_fitter') and self.integrator.enhanced_fitter:
                        self.integrator.enhanced_fitter.fitting_parameters['fitting_window_x'] = window_x_value
                    self.logger.info(f"Configured fitting window X: ±{window_x_value} ppm")
                except (TypeError, ValueError) as exc:
                    self.logger.warning(f"Invalid fitting_window_x override ({fitting_window_x}): {exc}")

            fitting_window_y = proc_param_config.get('fitting_window_y')
            if fitting_window_y is not None:
                try:
                    window_y_value = float(fitting_window_y)
                    self.integrator.fitting_parameters['fitting_window_y'] = window_y_value
                    if hasattr(self.integrator, 'enhanced_fitter') and self.integrator.enhanced_fitter:
                        self.integrator.enhanced_fitter.fitting_parameters['fitting_window_y'] = window_y_value
                    self.logger.info(f"Configured fitting window Y: ±{window_y_value} ppm")
                except (TypeError, ValueError) as exc:
                    self.logger.warning(f"Invalid fitting_window_y override ({fitting_window_y}): {exc}")

            if hasattr(self.integrator, 'ml_data_collector') and self.integrator.ml_data_collector:
                try:
                    self.integrator.ml_data_collector.update_window_settings(proc_param_config)
                except AttributeError:
                    self.logger.debug("ML data collector does not support window settings update (legacy version)")

            self.logger.info("Batch processor initialized with ML data collection")
        else:
            self.integrator = None
            self.logger.error("Core components not available - batch processing disabled")

    def _setup_logging(self):
        """Set up logging system for batch processing."""
        # Create logger
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)

        # Create console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)

        # Create file handler
        log_file = f"batch_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)

        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        console_handler.setFormatter(formatter)
        file_handler.setFormatter(formatter)

        # Add handlers
        logger.addHandler(console_handler)
        logger.addHandler(file_handler)

    def detect_nucleus_type(self, filename: str, spectrum_data: Optional[np.ndarray] = None,
                          f1_range: Optional[tuple] = None, f2_range: Optional[tuple] = None,
                          return_details: bool = False) -> Union[str, Dict[str, Any]]:
        """
        Enhanced nucleus type detection using chemical shifts and filename patterns.

        Args:
            filename: Name of the spectrum file
            spectrum_data: Optional spectrum data (maintained for backward compatibility)
            f1_range: Optional tuple of (min, max) chemical shifts for F1 dimension
            f2_range: Optional tuple of (min, max) chemical shifts for F2 dimension
            return_details: If True, return detailed detection results, else just nucleus type string

        Returns:
            If return_details=False: Detected nucleus type string ('15N1H', '13C1H', '1H', or 'unknown')
            If return_details=True: Dictionary with detailed detection results and confidence scores
        """
        detection_results = []

        # Method 1: Chemical shift analysis (highest priority if available)
        if f1_range is not None and f2_range is not None:
            try:
                chem_result = self._detect_from_chemical_shifts(f1_range, f2_range)
                detection_results.append(chem_result)

                # If high confidence chemical shift detection, use it immediately
                if chem_result['confidence'] >= 0.8:
                    if return_details:
                        return {
                            'nucleus_type': chem_result['nucleus_type'],
                            'confidence': chem_result['confidence'],
                            'primary_method': 'chemical_shift',
                            'all_methods': [chem_result],
                            'filename': filename,
                            'chemical_shifts': {'f1_range': f1_range, 'f2_range': f2_range}
                        }
                    else:
                        return chem_result['nucleus_type']
            except Exception as e:
                self.logger.debug(f"Chemical shift detection failed for {filename}: {e}")

        # Method 2: Filename pattern matching (always performed for compatibility)
        try:
            filename_result = self._detect_from_filename(filename)
            detection_results.append(filename_result)
        except Exception as e:
            self.logger.debug(f"Filename detection failed for {filename}: {e}")

        # Method 3: Combine results and resolve conflicts
        final_result = self._resolve_detection_conflicts(detection_results, filename)

        # Log detection details for debugging
        if len(detection_results) > 1:
            methods_used = [r['method'] for r in detection_results]
            confidences = [r['confidence'] for r in detection_results]
            self.logger.debug(f"Detection for {filename}: methods={methods_used}, confidences={confidences}, final={final_result['nucleus_type']}")

        # Return result in requested format
        if return_details:
            return {
                'nucleus_type': final_result['nucleus_type'],
                'confidence': final_result['confidence'],
                'primary_method': final_result['method'],
                'all_methods': detection_results,
                'filename': filename,
                'chemical_shifts': {'f1_range': f1_range, 'f2_range': f2_range} if f1_range and f2_range else None
            }
        else:
            # Backward compatibility: return just the nucleus type string
            return final_result['nucleus_type']

    def _resolve_detection_conflicts(self, detection_results: List[Dict[str, Any]], filename: str) -> Dict[str, Any]:
        """
        Resolve conflicts between different detection methods.

        Args:
            detection_results: List of detection results from different methods
            filename: Filename for logging purposes

        Returns:
            Final detection result with highest confidence/priority
        """
        if not detection_results:
            self.logger.warning(f"Could not determine nucleus type for {filename}, using default")
            return {
                'nucleus_type': 'unknown',
                'confidence': 0.0,
                'method': 'default'
            }

        # Filter out unknown results if we have any confident detections
        confident_results = [r for r in detection_results if r['nucleus_type'] != 'unknown' and r['confidence'] > 0.5]

        if confident_results:
            # Check for conflicts between high-confidence methods
            chemical_shift_results = [r for r in confident_results if r['method'] == 'chemical_shift']
            filename_results = [r for r in confident_results if r['method'] == 'filename']

            # If we have both chemical shift and filename results, check for conflicts
            if chemical_shift_results and filename_results:
                chem_nucleus = chemical_shift_results[0]['nucleus_type']
                filename_nucleus = filename_results[0]['nucleus_type']

                if chem_nucleus != filename_nucleus:
                    self.logger.warning(f"Detection conflict for {filename}: "
                                      f"Chemical shifts suggest {chem_nucleus} "
                                      f"(confidence: {chemical_shift_results[0]['confidence']:.2f}), "
                                      f"filename suggests {filename_nucleus} "
                                      f"(confidence: {filename_results[0]['confidence']:.2f})")

                    # Trust chemical shifts over filename (they are more reliable)
                    if chemical_shift_results[0]['confidence'] >= 0.7:
                        self.logger.info(f"Using chemical shift detection: {chem_nucleus}")
                        return chemical_shift_results[0]

            # Return the result with highest confidence
            best_result = max(confident_results, key=lambda x: x['confidence'])
            return best_result
        else:
            # No confident results, return the best of what we have
            best_result = max(detection_results, key=lambda x: x['confidence'])
            if best_result['nucleus_type'] == 'unknown':
                self.logger.warning(f"Could not determine nucleus type for {filename}, using default")
            return best_result

    def _detect_from_chemical_shifts(self, f1_range: tuple, f2_range: tuple) -> Dict[str, Any]:
        """
        Detect nucleus type from F1/F2 chemical shift ranges.

        Args:
            f1_range: Tuple of (min, max) chemical shifts for F1 dimension
            f2_range: Tuple of (min, max) chemical shifts for F2 dimension

        Returns:
            Dictionary with detection results
        """
        f1_min, f1_max = f1_range
        f2_min, f2_max = f2_range

        # Chemical shift signatures for different nucleus types
        signatures = {
            '15N1H': {
                'f1_range': (90, 140),      # 15N backbone/side chains
                'f2_range': (5, 12),        # Amide protons
                'confidence_base': 0.9
            },
            '13C1H': {
                'f1_ranges': [(0, 60), (100, 180)],  # Aliphatic OR aromatic
                'f2_range': (-1, 10),       # Broad 1H range
                'confidence_base': 0.8
            },
            '1H': {
                'f1_range': (-1, 12),       # 1H range
                'f2_range': (-1, 12),       # Same range (homonuclear)
                'confidence_base': 0.7
            }
        }

        detection_results = []

        for nucleus_type, sig in signatures.items():
            confidence = 0.0

            if nucleus_type == '15N1H':
                # Check if ranges match 15N1H signature
                f1_match = sig['f1_range'][0] <= f1_min <= f1_max <= sig['f1_range'][1]
                f2_match = sig['f2_range'][0] <= f2_min <= f2_max <= sig['f2_range'][1]

                if f1_match and f2_match:
                    confidence = sig['confidence_base']
                elif f1_match or f2_match:
                    confidence = sig['confidence_base'] * 0.6

            elif nucleus_type == '13C1H':
                # Check aliphatic OR aromatic carbon ranges
                f1_aliphatic = sig['f1_ranges'][0][0] <= f1_min <= f1_max <= sig['f1_ranges'][0][1]
                f1_aromatic = sig['f1_ranges'][1][0] <= f1_min <= f1_max <= sig['f1_ranges'][1][1]
                f2_match = sig['f2_range'][0] <= f2_min <= f2_max <= sig['f2_range'][1]

                if (f1_aliphatic or f1_aromatic) and f2_match:
                    confidence = sig['confidence_base']
                elif f1_aliphatic or f1_aromatic or f2_match:
                    confidence = sig['confidence_base'] * 0.5

            elif nucleus_type == '1H':
                # Check if F1 and F2 ranges are similar (homonuclear)
                f1_match = sig['f1_range'][0] <= f1_min <= f1_max <= sig['f1_range'][1]
                f2_match = sig['f2_range'][0] <= f2_min <= f2_max <= sig['f2_range'][1]

                # Calculate range similarity
                f1_center = (f1_min + f1_max) / 2
                f2_center = (f2_min + f2_max) / 2
                range_similarity = 1.0 - min(abs(f1_center - f2_center) / 12.0, 1.0)

                if f1_match and f2_match and range_similarity > 0.8:
                    confidence = sig['confidence_base'] * range_similarity

            if confidence > 0:
                detection_results.append({
                    'nucleus_type': nucleus_type,
                    'confidence': confidence,
                    'method': 'chemical_shift',
                    'details': {
                        'f1_range': f1_range,
                        'f2_range': f2_range,
                        'signature_match': sig
                    }
                })

        # Return the highest confidence result, or unknown if none found
        if detection_results:
            best_result = max(detection_results, key=lambda x: x['confidence'])
            return best_result
        else:
            return {
                'nucleus_type': 'unknown',
                'confidence': 0.0,
                'method': 'chemical_shift',
                'details': {'f1_range': f1_range, 'f2_range': f2_range}
            }

    def _detect_from_filename(self, filename: str) -> Dict[str, Any]:
        """
        Detect nucleus type from filename patterns (original method refactored).

        Args:
            filename: Name of the spectrum file

        Returns:
            Dictionary with detection results
        """
        filename_lower = filename.lower()

        # Primary filename patterns
        if '15n' in filename_lower or 'hsqc' in filename_lower:
            return {
                'nucleus_type': '15N1H',
                'confidence': 0.9,
                'method': 'filename',
                'details': {'pattern_found': '15n or hsqc in filename'}
            }
        elif '13c' in filename_lower or 'hmqc' in filename_lower:
            return {
                'nucleus_type': '13C1H',
                'confidence': 0.9,
                'method': 'filename',
                'details': {'pattern_found': '13c or hmqc in filename'}
            }
        elif '1h' in filename_lower and ('15n' not in filename_lower and '13c' not in filename_lower):
            return {
                'nucleus_type': '1H',
                'confidence': 0.9,
                'method': 'filename',
                'details': {'pattern_found': '1h without 15n/13c in filename'}
            }

        # Extended pattern search
        nucleus_patterns = {
            '15N1H': ['n15', 'nitrogen', '15n1h', 'hn', 'hsqc'],
            '13C1H': ['c13', 'carbon', '13c1h', 'hc', 'hmqc'],
            '1H': ['proton', '1h', 'h1']
        }

        for nucleus_type, patterns in nucleus_patterns.items():
            for pattern in patterns:
                if pattern in filename_lower:
                    return {
                        'nucleus_type': nucleus_type,
                        'confidence': 0.8,
                        'method': 'filename',
                        'details': {'pattern_found': f'{pattern} in filename'}
                    }

        # No patterns found
        return {
            'nucleus_type': 'unknown',
            'confidence': 0.0,
            'method': 'filename',
            'details': {'patterns_searched': list(nucleus_patterns.keys())}
        }

    def _extract_chemical_shift_ranges_from_integrator(self) -> Optional[Dict[str, tuple]]:
        """
        Extract chemical shift ranges from loaded spectrum in integrator.

        Returns:
            Dictionary with f1_range and f2_range tuples, or None if not available
        """
        try:
            if not hasattr(self.integrator, 'spectrum_data') or self.integrator.spectrum_data is None:
                return None

            ranges = {}

            # Method 1: Check for direct ppm scale attributes
            if hasattr(self.integrator, 'f1_ppm') and self.integrator.f1_ppm is not None:
                f1_ppm = self.integrator.f1_ppm
                ranges['f1_range'] = (float(min(f1_ppm)), float(max(f1_ppm)))

            if hasattr(self.integrator, 'f2_ppm') and self.integrator.f2_ppm is not None:
                f2_ppm = self.integrator.f2_ppm
                ranges['f2_range'] = (float(min(f2_ppm)), float(max(f2_ppm)))

            # Method 2: Check for ppm_scale attribute (for F2/X dimension)
            if 'f2_range' not in ranges and hasattr(self.integrator, 'ppm_scale'):
                ppm_scale = self.integrator.ppm_scale
                if ppm_scale is not None and len(ppm_scale) > 0:
                    ranges['f2_range'] = (float(min(ppm_scale)), float(max(ppm_scale)))

            # Method 3: Try to extract from spectrum metadata attributes
            for attr_name in ['ppm_x', 'ppm_f2', 'x_ppm']:
                if 'f2_range' not in ranges and hasattr(self.integrator, attr_name):
                    ppm_data = getattr(self.integrator, attr_name)
                    if ppm_data is not None and len(ppm_data) > 0:
                        ranges['f2_range'] = (float(min(ppm_data)), float(max(ppm_data)))
                        break

            for attr_name in ['ppm_y', 'ppm_f1', 'y_ppm']:
                if 'f1_range' not in ranges and hasattr(self.integrator, attr_name):
                    ppm_data = getattr(self.integrator, attr_name)
                    if ppm_data is not None and len(ppm_data) > 0:
                        ranges['f1_range'] = (float(min(ppm_data)), float(max(ppm_data)))
                        break

            # Method 4: Parse from recent log output (fallback)
            if len(ranges) < 2:
                ranges_from_log = self._parse_chemical_shifts_from_log()
                if ranges_from_log:
                    ranges.update(ranges_from_log)

            # Ensure we have both dimensions
            if 'f1_range' in ranges and 'f2_range' in ranges:
                self.logger.debug(f"Extracted chemical shift ranges: F1={ranges['f1_range']}, F2={ranges['f2_range']}")
                return ranges
            else:
                self.logger.debug(f"Incomplete chemical shift ranges: {ranges}")
                return None

        except Exception as e:
            self.logger.debug(f"Could not extract chemical shift ranges: {e}")
            return None

    def _parse_chemical_shifts_from_log(self) -> Optional[Dict[str, tuple]]:
        """
        Parse chemical shift ranges from log output as fallback.

        Returns:
            Dictionary with chemical shift ranges or None
        """
        try:
            # This is a fallback method that could parse from captured log output
            # For now, return None and rely on the direct attribute methods
            return None
        except:
            return None

    def get_optimal_parameters(self, nucleus_type: str, spectrum_file: Path) -> Dict[str, Any]:
        """
        Get optimal parameters for a given nucleus type and spectrum.

        Args:
            nucleus_type: Type of nucleus ('15N1H', '13C1H', etc.)
            spectrum_file: Path to the spectrum file

        Returns:
            Dictionary with optimal processing parameters
        """
        # Get base parameters from config
        sn_threshold = self.config['sn_thresholds'].get(nucleus_type,
                                                       self.config['sn_thresholds']['default'])
        expected_peaks = self.config['expected_peaks'].get(nucleus_type,
                                                          self.config['expected_peaks']['default'])

        parameters = {
            'nucleus_type': nucleus_type,
            'sn_threshold': sn_threshold,
            'expected_peak_count': expected_peaks,
            'quality_threshold': self.config['quality_threshold'],
            'max_attempts': self.config['max_fitting_attempts']
        }

        self.logger.debug(f"Parameters for {nucleus_type}: {parameters}")
        return parameters

    def process_single_spectrum(self, spectrum_file: Path, parameters: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process a single spectrum file.

        Args:
            spectrum_file: Path to the spectrum file
            parameters: Processing parameters

        Returns:
            Dictionary with processing results
        """
        result = {
            'filename': spectrum_file.name,
            'success': False,
            'peaks_detected': 0,
            'peaks_fitted': 0,
            'ml_samples_collected': 0,
            'processing_time': 0,
            'error': None
        }

        start_time = time.time()

        try:
            self.logger.info(f"Processing: {spectrum_file.name}")

            # Record ML sample count at start of this spectrum
            spectrum_start_ml_samples = len(self.integrator.ml_data_collector.session_data) if \
                                        (self.integrator.ml_data_collector) else 0

            # Load spectrum
            if self._load_spectrum(spectrum_file):
                # Set up integrator parameters
                self.integrator.sn_threshold = parameters['sn_threshold']
                self.integrator.expected_peak_count = parameters['expected_peak_count']

                # Step 1: Perform S/N detection to find peaks
                detected_peaks = self.integrator.detect_peaks_sn_native()

                if detected_peaks and len(detected_peaks) > 0:
                    result['peaks_detected'] = len(detected_peaks)
                    self.logger.info(f"  Detected {len(detected_peaks)} peaks")

                    # Step 2: Perform proper Voigt fitting on detected peaks
                    fitted_count = self._perform_voigt_fitting(detected_peaks, parameters)
                    result['peaks_fitted'] = fitted_count

                    # Check ML data collection for this spectrum
                    if self.integrator.ml_data_collector:
                        spectrum_end_ml_samples = len(self.integrator.ml_data_collector.session_data)
                        result['ml_samples_collected'] = spectrum_end_ml_samples - spectrum_start_ml_samples

                    result['success'] = True
                    self.logger.info(f"  Successfully processed: {fitted_count} peaks fitted, "
                                   f"{result['ml_samples_collected']} ML samples collected")

                else:
                    self.logger.warning(f"  No peaks detected in {spectrum_file.name}")
                    result['error'] = "No peaks detected"

            else:
                result['error'] = "Failed to load spectrum"

        except Exception as e:
            result['error'] = str(e)
            self.logger.error(f"  Error processing {spectrum_file.name}: {e}")

        result['processing_time'] = time.time() - start_time
        return result

    def _load_spectrum(self, spectrum_file: Path) -> bool:
        """
        Load a spectrum file using the core integrator.

        Args:
            spectrum_file: Path to the spectrum file

        Returns:
            True if loading successful, False otherwise
        """
        try:
            # Use the core integrator's load_spectrum_only method for batch processing
            # This method loads spectra without requiring peak lists
            if hasattr(self.integrator, 'load_spectrum_only'):
                success = self.integrator.load_spectrum_only(str(spectrum_file))
                return success
            elif hasattr(self.integrator, 'load_nmr_file'):
                # Fallback to load_nmr_file method
                success = self.integrator.load_nmr_file(str(spectrum_file))
                return success
            else:
                self.logger.error("No suitable spectrum loading method available")
                return False

        except Exception as e:
            self.logger.error(f"Failed to load {spectrum_file.name}: {e}")
            return False

    def _perform_voigt_fitting(self, detected_peaks: List[Dict], parameters: Dict[str, Any]) -> int:
        """
        Perform comprehensive Voigt fitting on detected peaks with parallel processing and enhanced ML data collection.

        This method implements the complete Voigt fitting pipeline:
        1. Parallel or sequential Voigt profile fitting
        2. Progress tracking and detailed logging
        3. Comprehensive error handling and recovery
        4. Enhanced ML data collection integration
        5. Quality assessment and filtering

        Args:
            detected_peaks: List of detected peak dictionaries from S/N detection
            parameters: Processing parameters including quality thresholds

        Returns:
            Number of successfully fitted peaks
        """
        if not detected_peaks:
            self.logger.warning("  No detected peaks provided for Voigt fitting")
            return 0

        total_peaks = len(detected_peaks)
        self.logger.info(f"  🚀 Starting comprehensive Voigt fitting on {total_peaks} detected peaks")

        try:
            # Enhanced parameter setup
            use_parallel = total_peaks >= 10  # Use parallel processing for 10+ peaks
            quality_threshold = parameters.get('quality_threshold', 0.8)
            max_attempts = parameters.get('max_attempts', 3)
            nucleus_type = parameters.get('nucleus_type', 'unknown')

            self.logger.info(f"     Parallel processing: {use_parallel}")
            self.logger.info(f"     Quality threshold: {quality_threshold}")
            self.logger.info(f"     Max fitting attempts: {max_attempts}")
            self.logger.info(f"     Nucleus type: {nucleus_type}")

            # Prepare peak data for fitting
            import pandas as pd
            fitted_results = []
            failed_peaks = []
            processing_start_time = time.time()

            # Use direct enhanced_peak_fitting calls for guaranteed ML data collection
            # with fallback to SingleSpectrumProcessor for compatibility
            try:
                self.logger.info(f"     Using direct enhanced_peak_fitting (guaranteed ML collection)")
                fitted_results = self._fit_peaks_with_enhanced_peak_fitting(
                    detected_peaks, parameters, use_parallel
                )
            except Exception as e:
                self.logger.warning(f"     Direct fitting failed ({e}), falling back to SingleSpectrumProcessor")
                fitted_results = self._fit_peaks_with_single_spectrum_processor(
                    detected_peaks, parameters, use_parallel
                )

            # Debug: Check first result structure
            if fitted_results and len(fitted_results) > 0:
                first_result = fitted_results[0]
                if first_result:
                    self.logger.info(f"     DEBUG: First result keys: {list(first_result.keys())}")
                    self.logger.info(f"     DEBUG: Success field: {first_result.get('success', 'MISSING')}")
                    avg_r2 = first_result.get('avg_r_squared', 'MISSING')
                    r2 = first_result.get('r_squared', 'MISSING')
                    R2 = first_result.get('R2', 'MISSING')
                    self.logger.info(f"     DEBUG: R² fields - avg_r_squared: {avg_r2}, r_squared: {r2}, R2: {R2}")

            # Process results and collect statistics
            processing_time = time.time() - processing_start_time

            if fitted_results:
                # More robust quality filtering with multiple R² field names
                successful_fits = []
                failed_fits = []

                for r in fitted_results:
                    if r is None:
                        failed_fits.append(r)
                        continue

                    # Check for success flag (more permissive)
                    success = r.get('success', True)  # Default to True if missing
                    if not success:
                        failed_fits.append(r)
                        continue

                    # Check R² with multiple possible field names
                    r_squared = r.get('avg_r_squared', 0) or r.get('r_squared', 0) or r.get('R2', 0)

                    if r_squared >= quality_threshold:
                        successful_fits.append(r)
                    else:
                        failed_fits.append(r)

                # Store results in integrator for ML data collection
                self.integrator.fitted_peaks = successful_fits

                success_count = len(successful_fits)
                failure_count = len(failed_fits)
                success_rate = (success_count / total_peaks) * 100 if total_peaks > 0 else 0

                self.logger.info(f"  ✅ Voigt fitting completed in {processing_time:.2f}s")
                self.logger.info(f"     Successful fits: {success_count}/{total_peaks} ({success_rate:.1f}%)")
                self.logger.info(f"     Failed/low-quality fits: {failure_count}")

                if failure_count > 0:
                    self.logger.debug(f"     Failed peaks: {[f.get('assignment', 'Unknown') for f in failed_fits if f]}")

                # ML data collection is now handled automatically by SingleSpectrumProcessor
                # No manual collection needed - it happens during the fitting process

                return success_count

            else:
                self.logger.warning(f"  ❌ Voigt fitting failed - no results returned")
                return 0

        except Exception as e:
            self.logger.error(f"  ❌ Voigt fitting failed with error: {e}")
            import traceback
            self.logger.debug(f"     Full traceback: {traceback.format_exc()}")
            return 0

    def _fit_peaks_with_single_spectrum_processor(self, detected_peaks: List[Dict], parameters: Dict[str, Any], use_parallel: bool) -> List[Dict]:
        """
        Use the EXACT same SingleSpectrumProcessor as GUI "Fit All Peaks" and Series Integration.

        This is the identical implementation used by:
        - GUI main_gui.py._fit_all_peaks_with_new_processor()
        - Series Integration processors
        - All parallel processing capabilities
        """
        try:
            import pandas as pd

            # Dynamic import of processor components for maximum compatibility
            try:
                from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor
                from lunaNMR.utils.parameter_manager import NMRParameterManager
                processor_available = True
            except ImportError as e:
                self.logger.warning(f"     ⚠️ SingleSpectrumProcessor not available ({e}), using fallback")
                return self._fit_peaks_fallback_method(detected_peaks, parameters)

            self.logger.info(f"     🚀 Initializing SingleSpectrumProcessor (parallel={use_parallel})")

            # Convert detected peaks to the exact DataFrame format expected by SingleSpectrumProcessor
            peak_list_data = []
            for peak in detected_peaks:
                peak_list_data.append({
                    'Assignment': peak['assignment'],
                    'Position_X': peak['ppm_x'],
                    'Position_Y': peak['ppm_y'],
                    'X_HZ': 0,  # Placeholder - same as GUI
                    'Y_HZ': 0   # Placeholder - same as GUI
                })

            peak_list = pd.DataFrame(peak_list_data)

            # CRITICAL: Store peak list in integrator exactly like GUI does
            # This is essential for ML data collection to work properly
            self.integrator.peak_list = peak_list

            # Create parameter manager (same as GUI setup)
            param_manager = NMRParameterManager()

            # CRITICAL: Sync parameter manager exactly like GUI's update_from_gui_variables
            # This ensures the SingleSpectrumProcessor has the correct configuration
            current_params = {
                'use_parallel_processing': use_parallel,
                'use_global_optimization': False,  # Standard batch processing
                'quality_threshold': parameters.get('quality_threshold', 0.6),  # Use batch threshold
                'max_attempts': parameters.get('max_attempts', 3),
                'sn_threshold': parameters.get('sn_threshold', 2.2)
            }

            processing_overrides = self.config.get('processing_parameters', {})

            numeric_overrides = {
                key: processing_overrides[key]
                for key in ['search_window_x', 'search_window_y', 'fitting_window_x', 'fitting_window_y', 'noise_threshold']
                if key in processing_overrides
            }
            for key, value in numeric_overrides.items():
                try:
                    current_params[key] = float(value)
                except (TypeError, ValueError):
                    self.logger.warning(f"Invalid numeric override for {key}: {value}")

            if 'use_reference_detection' in processing_overrides:
                current_params['use_reference_detection'] = bool(processing_overrides['use_reference_detection'])

            param_manager.current_params.update(current_params)

            # Initialize SingleSpectrumProcessor (EXACT same as GUI)
            single_spectrum_processor = SingleSpectrumProcessor(self.integrator, param_manager)

            # Set up progress callback (same pattern as GUI)
            def progress_callback(progress, task, log_msg=None, failed=False):
                if progress % 10 == 0 or failed:  # Log every 10% or on failure
                    status = "❌ FAILED" if failed else "✅ Processing"
                    self.logger.info(f"     {status}: {progress:.1f}% - {task}")
                    if log_msg:
                        self.logger.debug(f"       {log_msg}")

            # Set processing options (EXACT same as GUI)
            processing_options = {
                'use_parallel': use_parallel,
                'use_global_optimization': False  # Standard fitting for batch
            }

            self.logger.info(f"     📋 Processing {len(peak_list)} peaks with options: {processing_options}")

            # Record ML samples before processing (exactly like GUI)
            ml_samples_before = len(self.integrator.ml_data_collector.session_data) if hasattr(self.integrator, 'ml_data_collector') and self.integrator.ml_data_collector else 0

            # Process all peaks using the EXACT same method as GUI
            fitted_results = single_spectrum_processor.process_peak_list(
                peak_list,
                processing_options,
                progress_callback
            )

            # Check ML data collection results (exactly like GUI monitoring)
            ml_samples_after = len(self.integrator.ml_data_collector.session_data) if hasattr(self.integrator, 'ml_data_collector') and self.integrator.ml_data_collector else 0
            ml_samples_collected = ml_samples_after - ml_samples_before

            if ml_samples_collected > 0:
                self.logger.info(f"     🤖 ML data automatically collected: {ml_samples_collected} new samples")
            else:
                self.logger.warning(f"     ⚠️ No ML data collected automatically - this indicates an issue")

            # Store fitted results in integrator exactly like GUI does
            if fitted_results:
                # Get comprehensive summary like GUI
                summary = single_spectrum_processor.get_processing_summary(fitted_results, len(peak_list))

                # Standardize results exactly like GUI's _complete_new_single_spectrum_fitting
                standardized_results = []
                for result in fitted_results:
                    standardized_result = {
                        'assignment': result.get('Assignment', result.get('assignment', '')),
                        'ppm_x': float(result.get('Position_X', result.get('ppm_x', 0))),
                        'ppm_y': float(result.get('Position_Y', result.get('ppm_y', 0))),
                        'detected': True,
                        'fitted': True
                    }
                    # Preserve additional fitting data
                    for key, value in result.items():
                        if key not in ['Assignment', 'Position_X', 'Position_Y']:
                            standardized_result[key] = value
                    standardized_results.append(standardized_result)

                # Store in integrator exactly like GUI
                self.integrator.fitted_peaks = standardized_results
                self.logger.info(f"     🔧 Standardized {len(standardized_results)} fitting results for compatibility")

            self.logger.info(f"     ✅ SingleSpectrumProcessor completed: {len(fitted_results)} results")
            return fitted_results

        except Exception as e:
            self.logger.error(f"     ❌ SingleSpectrumProcessor failed: {e}")
            import traceback
            self.logger.debug(f"     Full traceback: {traceback.format_exc()}")
            return []

    def _fit_peaks_with_enhanced_peak_fitting(self, detected_peaks: List[Dict], parameters: Dict[str, Any], use_parallel: bool) -> List[Dict]:
        """
        Direct enhanced_peak_fitting implementation with guaranteed ML data collection.

        This method uses the exact same enhanced_peak_fitting calls as GUI individual peak fitting,
        ensuring 100% ML data collection compatibility with optional multiprocessing support.
        """
        import time
        import os

        total_peaks = len(detected_peaks)
        fitted_results = []

        # Record ML samples before processing
        ml_samples_before = len(self.integrator.ml_data_collector.session_data) if hasattr(self.integrator, 'ml_data_collector') and self.integrator.ml_data_collector else 0

        self.logger.info(f"     🎯 Direct enhanced_peak_fitting: {total_peaks} peaks")
        self.logger.info(f"     🔄 Parallel processing: {use_parallel}")

        if use_parallel and total_peaks > 10:
            # Use multiprocessing for large peak sets
            fitted_results = self._fit_peaks_parallel_direct(detected_peaks, parameters)
        else:
            # Sequential processing for small peak sets or when parallel is disabled
            fitted_results = self._fit_peaks_sequential_direct(detected_peaks, parameters)

        # Check ML data collection results
        ml_samples_after = len(self.integrator.ml_data_collector.session_data) if hasattr(self.integrator, 'ml_data_collector') and self.integrator.ml_data_collector else 0
        ml_samples_collected = ml_samples_after - ml_samples_before

        success_count = len([r for r in fitted_results if r and r.get('success', True)])

        self.logger.info(f"     ✅ Direct fitting completed: {success_count}/{total_peaks} successful")
        self.logger.info(f"     🤖 ML samples collected: {ml_samples_collected}")

        return fitted_results

    def _fit_peaks_sequential_direct(self, detected_peaks: List[Dict], parameters: Dict[str, Any]) -> List[Dict]:
        """Sequential direct enhanced_peak_fitting calls"""
        fitted_results = []
        total_peaks = len(detected_peaks)

        for i, peak in enumerate(detected_peaks):
            try:
                # Progress logging every 20%
                if i % max(1, total_peaks // 5) == 0:
                    progress = (i / total_peaks) * 100
                    self.logger.info(f"     📈 Sequential progress: {progress:.0f}% ({i}/{total_peaks})")

                # Direct enhanced_peak_fitting call (same as GUI)
                result = self.integrator.enhanced_peak_fitting(
                    peak['ppm_x'],
                    peak['ppm_y'],
                    peak['assignment']
                )

                if result:
                    # Add batch-specific metadata
                    result['peak_number'] = i + 1
                    result['processing_mode'] = 'batch_direct_sequential'
                    result['success'] = True
                    fitted_results.append(result)

            except Exception as e:
                self.logger.debug(f"     ❌ Failed to fit {peak['assignment']}: {e}")
                continue

        return fitted_results

    def _fit_peaks_parallel_direct(self, detected_peaks: List[Dict], parameters: Dict[str, Any]) -> List[Dict]:
        """Parallel direct enhanced_peak_fitting using threading for ML data collection compatibility"""
        import threading
        import queue
        from concurrent.futures import ThreadPoolExecutor, as_completed

        fitted_results = []
        total_peaks = len(detected_peaks)

        # Use threading instead of multiprocessing to maintain ML data collector access
        max_workers = min(os.cpu_count() or 4, 8)  # Reasonable limit

        self.logger.info(f"     🚀 Threading parallel: {max_workers} workers")

        def fit_single_peak(peak_data):
            """Fit a single peak using direct enhanced_peak_fitting"""
            peak, index = peak_data
            try:
                result = self.integrator.enhanced_peak_fitting(
                    peak['ppm_x'],
                    peak['ppm_y'],
                    peak['assignment']
                )

                if result:
                    result['peak_number'] = index + 1
                    result['processing_mode'] = 'batch_direct_parallel'
                    result['success'] = True
                    return result
                return None

            except Exception as e:
                return None

        # Process peaks in parallel using ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all peaks for processing
            future_to_peak = {
                executor.submit(fit_single_peak, (peak, i)): i
                for i, peak in enumerate(detected_peaks)
            }

            completed = 0
            for future in as_completed(future_to_peak):
                completed += 1

                # Progress logging every 20%
                if completed % max(1, total_peaks // 5) == 0:
                    progress = (completed / total_peaks) * 100
                    self.logger.info(f"     📈 Parallel progress: {progress:.0f}% ({completed}/{total_peaks})")

                try:
                    result = future.result()
                    if result:
                        fitted_results.append(result)
                except Exception as e:
                    continue

        # Sort results by peak number to maintain order
        fitted_results.sort(key=lambda x: x.get('peak_number', 0))

        return fitted_results

    def _fit_peaks_fallback_method(self, detected_peaks: List[Dict], parameters: Dict[str, Any]) -> List[Dict]:
        """
        Fallback method for when SingleSpectrumProcessor is not available.
        Maintains backward compatibility and essential functionality.
        """
        try:
            import pandas as pd

            self.logger.info(f"     🔄 Using fallback sequential fitting method")

            # Convert detected peaks to peak list format
            peak_list_data = []
            for peak in detected_peaks:
                peak_list_data.append({
                    'Assignment': peak['assignment'],
                    'Position_X': peak['ppm_x'],
                    'Position_Y': peak['ppm_y'],
                    'X_HZ': 0,
                    'Y_HZ': 0
                })

            # Create temporary peak list for compatibility
            self.integrator.peak_list = pd.DataFrame(peak_list_data)

            # Switch to fitting mode
            self.integrator.set_processing_mode('full_detection')

            # Fit peaks sequentially (reliable fallback)
            fitted_results = []
            total_peaks = len(detected_peaks)

            for i, peak in enumerate(detected_peaks):
                try:
                    # Progress logging
                    if i % max(1, total_peaks // 10) == 0:
                        progress = (i / total_peaks) * 100
                        self.logger.info(f"     Fallback progress: {progress:.1f}% ({i+1}/{total_peaks}) - {peak['assignment']}")

                    # Use enhanced_peak_fitting method
                    result = self.integrator.enhanced_peak_fitting(
                        peak['ppm_x'],
                        peak['ppm_y'],
                        peak['assignment']
                    )

                    if result:
                        fitted_results.append(result)
                    else:
                        self.logger.debug(f"     Failed to fit {peak['assignment']}")

                except Exception as e:
                    self.logger.debug(f"     Error fitting {peak['assignment']}: {e}")
                    continue

            self.logger.info(f"     ✅ Fallback method completed: {len(fitted_results)} results")
            return fitted_results

        except Exception as e:
            self.logger.error(f"     ❌ Fallback method failed: {e}")
            return []

    def _collect_ml_data_from_fit(self, fit_result: Dict, nucleus_type: str):
        """
        Collect enhanced ML training data from successful fit results.
        """
        try:
            # Check if result is valid (if it has essential fields, assume success)
            if not fit_result or not fit_result.get('avg_r_squared', 0):
                return

            # Extract fitting data for ML collection
            if hasattr(self.integrator, 'ml_data_collector') and self.integrator.ml_data_collector:
                # FORCE ML data collection for this successful fit
                assignment = fit_result.get('assignment', 'Unknown')
                r_squared = fit_result.get('avg_r_squared', 0)

                self.logger.debug(f"     Forcing ML data collection for {assignment}: R²={r_squared:.3f}")

                # Create synthetic training data from the fit result
                try:
                    # Get peak position from fit result (handle both dict and tuple formats)
                    peak_pos = fit_result.get('peak_position', None)
                    ppm_x, ppm_y = None, None

                    if isinstance(peak_pos, dict):
                        ppm_x = peak_pos.get('ppm_x')
                        ppm_y = peak_pos.get('ppm_y')
                    elif isinstance(peak_pos, (tuple, list)) and len(peak_pos) >= 2:
                        ppm_x, ppm_y = peak_pos[0], peak_pos[1]

                    if ppm_x is not None and ppm_y is not None:
                        # Get fit data
                        x_fit = fit_result.get('x_fit', {})
                        y_fit = fit_result.get('y_fit', {})

                        # Create minimal training sample with available data
                        sample_data = {
                            'assignment': assignment,
                            'ppm_x': ppm_x,
                            'ppm_y': ppm_y,
                            'r_squared': r_squared,
                            'nucleus_type': nucleus_type,
                            'fitting_quality': fit_result.get('fitting_quality', 'batch_processed'),
                            'source': 'batch_processor_forced'
                        }

                        # Try to collect with any available data structure
                        if x_fit and isinstance(x_fit, dict):
                            x_data = x_fit.get('x_data', [])
                            y_data = x_fit.get('y_data', [])
                            fit_params = x_fit.get('fit_params', {})

                            if x_data and y_data:
                                self.integrator.ml_data_collector.collect_training_sample(
                                    x_data=x_data,
                                    y_data=y_data,
                                    fit_params=fit_params,
                                    fit_result={'r_squared': r_squared},
                                    peak_info=sample_data,
                                    optimization_info={'method': 'batch_forced_collection'}
                                )
                                self.logger.debug(f"     ✅ ML data collected for {assignment} (full data)")
                            else:
                                # Fallback: create minimal sample
                                self.integrator.ml_data_collector.collect_training_sample(
                                    x_data=[ppm_x],  # Minimal data
                                    y_data=[ppm_y],
                                    fit_params={'r_squared': r_squared},
                                    fit_result={'r_squared': r_squared},
                                    peak_info=sample_data,
                                    optimization_info={'method': 'batch_minimal_collection'}
                                )
                                self.logger.debug(f"     ✅ ML data collected for {assignment} (minimal)")
                        else:
                            # Create minimal sample
                            self.integrator.ml_data_collector.collect_training_sample(
                                x_data=[ppm_x],
                                y_data=[ppm_y],
                                fit_params={'r_squared': r_squared},
                                fit_result={'r_squared': r_squared},
                                peak_info=sample_data,
                                optimization_info={'method': 'batch_minimal_collection'}
                            )
                            self.logger.debug(f"     ✅ ML data collected for {assignment} (minimal fallback)")
                    else:
                        self.logger.debug(f"     ⚠️ No position data for {assignment}")

                except Exception as e:
                    self.logger.debug(f"     ❌ Failed to force ML collection for {assignment}: {e}")

                self.logger.debug(f"     ML data for {assignment}: R²={r_squared:.3f}, Success={fit_result.get('success', False)}")

        except Exception as e:
            self.logger.debug(f"     ML data collection enhancement failed: {e}")
            pass

    def find_spectrum_files(self, folder_path: Union[str, Path]) -> List[Path]:
        """
        Find all NMR spectrum files in a folder.

        Args:
            folder_path: Path to the folder containing spectrum files

        Returns:
            List of Path objects for found spectrum files
        """
        folder = Path(folder_path)
        if not folder.exists():
            raise FileNotFoundError(f"Folder not found: {folder}")

        spectrum_files = []

        # Get file handling configuration with fallbacks
        file_config = self.config.get('file_handling', {})
        extensions = file_config.get('extensions', self.config.get('file_extensions', ['.ft2', '.ft', '.pipe', '.ucsf', '.nmrpipe']))
        recursive_search = file_config.get('recursive_search', False)
        pattern_filters = file_config.get('pattern_filters', [])
        exclude_patterns = file_config.get('exclude_patterns', [])

        # Determine search pattern (recursive or not)
        if recursive_search:
            self.logger.debug(f"Using recursive search in {folder}")
        else:
            self.logger.debug(f"Using non-recursive search in {folder}")

        # Find files by extension
        for ext in extensions:
            if recursive_search:
                pattern = f"**/*{ext}"
                found_files = list(folder.glob(pattern))
            else:
                pattern = f"*{ext}"
                found_files = list(folder.glob(pattern))
            spectrum_files.extend(found_files)
        self.logger.info(f"Found {len(spectrum_files)} spectrum files: {[f.name for f in spectrum_files]}")
        self.logger.debug(f"Found {len(spectrum_files)} files with extensions {extensions}")

        # Apply pattern filters (include only files matching these patterns)
        if pattern_filters:
            filtered_files = []
            for file_path in spectrum_files:
                for pattern in pattern_filters:
                    if self._matches_pattern(file_path.name, pattern):
                        filtered_files.append(file_path)
                        break  # File matches at least one pattern, include it
            spectrum_files = filtered_files
            self.logger.debug(f"After pattern filtering ({pattern_filters}): {len(spectrum_files)} files")

        # Apply exclude patterns (remove files matching these patterns)
        if exclude_patterns:
            excluded_files = []
            for file_path in spectrum_files:
                should_exclude = False
                for pattern in exclude_patterns:
                    if self._matches_pattern(file_path.name, pattern):
                        should_exclude = True
                        break
                if not should_exclude:
                    excluded_files.append(file_path)
            spectrum_files = excluded_files
            self.logger.debug(f"After exclusion filtering ({exclude_patterns}): {len(spectrum_files)} files")

        # Remove duplicates (can happen with recursive search and multiple extensions)
        spectrum_files = list(set(spectrum_files))

        # Sort by name for consistent processing order
        spectrum_files.sort(key=lambda x: x.name.lower())

        self.logger.info(f"Found {len(spectrum_files)} spectrum files in {folder}")

        # Log file discovery details in debug mode
        if self.logger.isEnabledFor(logging.DEBUG) and spectrum_files:
            self.logger.debug("Found spectrum files:")
            for i, file_path in enumerate(spectrum_files[:10]):  # Show first 10
                self.logger.debug(f"  {i+1:3d}. {file_path}")
            if len(spectrum_files) > 10:
                self.logger.debug(f"  ... and {len(spectrum_files) - 10} more files")

        return spectrum_files

    def _matches_pattern(self, filename: str, pattern: str) -> bool:
        """
        Check if a filename matches a glob-style pattern.

        Args:
            filename: Name of the file to check
            pattern: Glob-style pattern (supports * and ?)

        Returns:
            True if filename matches pattern, False otherwise
        """
        import fnmatch
        return fnmatch.fnmatch(filename.lower(), pattern.lower())

    def _apply_auto_optimization(self, spectrum_file: Path, nucleus_type: str, base_parameters: Dict[str, Any]) -> Dict[str, Any]:
        """
        Apply automatic parameter optimization for a specific spectrum.

        Args:
            spectrum_file: Path to the spectrum file
            nucleus_type: Detected nucleus type
            base_parameters: Base parameters to optimize

        Returns:
            Optimized parameters dictionary
        """
        try:
            # Import ParameterOptimizer
            from .parameter_optimizer import ParameterOptimizer

            # Create optimizer with config
            opt_config = self.config.get('optimization', {})
            optimizer = ParameterOptimizer(opt_config)

            self.logger.info(f"  🎯 Applying auto-optimization for {spectrum_file.name} ({nucleus_type})")

            # Check if we have the necessary spectrum data (should already be loaded)
            if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
                # Try to load spectrum data if not available
                if not self._load_spectrum(spectrum_file):
                    self.logger.warning(f"  ⚠️ Could not load spectrum for optimization, using base parameters")
                    return base_parameters

                # Check again after loading
                if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
                    self.logger.warning(f"  ⚠️ No spectrum data available for optimization, using base parameters")
                    return base_parameters

            # Extract spectrum characteristics
            spectrum_data = self.integrator.nmr_data
            f1_ppm = getattr(self.integrator, 'ppm_y_axis', None)  # Y-axis is F1 (15N/13C)
            f2_ppm = getattr(self.integrator, 'ppm_x_axis', None)  # X-axis is F2 (1H)

            if f1_ppm is None or f2_ppm is None:
                self.logger.warning(f"  ⚠️ Missing PPM scales for optimization, using base parameters")
                return base_parameters

            # Perform optimization
            optimized_params = optimizer.optimize_parameters(
                self.integrator, spectrum_data, f2_ppm, f1_ppm, nucleus_type
            )

            # Merge optimized parameters with base parameters
            updated_parameters = base_parameters.copy()

            # Update S/N threshold if optimized
            if 'sn_threshold' in optimized_params:
                updated_parameters['sn_threshold'] = optimized_params['sn_threshold']
                self.logger.info(f"  📈 Optimized S/N threshold: {optimized_params['sn_threshold']:.2f}")

            # Update expected peak count if optimized
            if 'expected_peak_count' in optimized_params:
                updated_parameters['expected_peak_count'] = optimized_params['expected_peak_count']
                self.logger.info(f"  📊 Optimized expected peaks: {optimized_params['expected_peak_count']}")

            # Log optimization strategy used
            strategy = opt_config.get('optimization_strategy', 'balanced')
            self.logger.info(f"  🎛️ Used {strategy} optimization strategy")

            return updated_parameters

        except ImportError:
            self.logger.warning("  ⚠️ ParameterOptimizer not available, using base parameters")
            return base_parameters
        except Exception as e:
            self.logger.warning(f"  ⚠️ Auto-optimization failed ({e}), using base parameters")
            return base_parameters

    def process_folder(self, folder_path: Union[str, Path],
                      nucleus_type: Optional[str] = None,
                      auto_optimize: bool = False) -> Dict[str, Any]:
        """
        Process all spectrum files in a folder.

        Args:
            folder_path: Path to folder containing spectrum files
            nucleus_type: Optional fixed nucleus type, if None will auto-detect
            auto_optimize: Whether to use automatic parameter optimization

        Returns:
            Dictionary with complete processing results
        """
        if not CORE_AVAILABLE:
            raise RuntimeError("Core components not available - cannot process spectra")

        self.logger.info(f"Starting batch processing of folder: {folder_path}")
        self.stats['processing_start_time'] = datetime.now()

        try:
            # Find spectrum files
            spectrum_files = self.find_spectrum_files(folder_path)
            if not spectrum_files:
                raise ValueError(f"No spectrum files found in {folder_path}")

            self.stats['total_files'] = len(spectrum_files)

            # Process each spectrum
            for i, spectrum_file in enumerate(spectrum_files):
                try:
                    # Progress reporting
                    if i % self.config['progress_interval'] == 0:
                        progress = (i / len(spectrum_files)) * 100
                        self.logger.info(f"Progress: {progress:.1f}% ({i}/{len(spectrum_files)})")

                    # Enhanced nucleus type detection with chemical shift analysis
                    if nucleus_type:
                        # User provided nucleus type - use it directly (backward compatibility)
                        detected_nucleus = nucleus_type
                        detection_details = None
                    else:
                        # Attempt to extract chemical shift ranges for enhanced detection
                        chemical_shift_ranges = None
                        try:
                            # Load spectrum temporarily to extract chemical shift information
                            if self._load_spectrum(spectrum_file):
                                chemical_shift_ranges = self._extract_chemical_shift_ranges_from_integrator()
                                if chemical_shift_ranges:
                                    self.logger.debug(f"Extracted chemical shifts for {spectrum_file.name}: {chemical_shift_ranges}")
                        except Exception as e:
                            self.logger.debug(f"Could not extract chemical shifts for {spectrum_file.name}: {e}")

                        # Perform enhanced nucleus detection
                        if chemical_shift_ranges and 'f1_range' in chemical_shift_ranges and 'f2_range' in chemical_shift_ranges:
                            # Use enhanced detection with chemical shifts
                            detection_details = self.detect_nucleus_type(
                                spectrum_file.name,
                                f1_range=chemical_shift_ranges['f1_range'],
                                f2_range=chemical_shift_ranges['f2_range'],
                                return_details=True
                            )
                            detected_nucleus = detection_details['nucleus_type']

                            # Log enhanced detection results
                            confidence = detection_details['confidence']
                            method = detection_details['primary_method']
                            self.logger.info(f"Enhanced detection for {spectrum_file.name}: {detected_nucleus} "
                                           f"(confidence: {confidence:.2f}, method: {method})")

                            if confidence >= 0.8:
                                self.logger.debug(f"High confidence detection using {method}")
                            elif len(detection_details['all_methods']) > 1:
                                self.logger.debug(f"Multiple detection methods used: {[m['method'] for m in detection_details['all_methods']]}")
                        else:
                            # Fallback to filename-only detection (backward compatibility)
                            detected_nucleus = self.detect_nucleus_type(spectrum_file.name)
                            detection_details = None
                            self.logger.debug(f"Using filename-only detection for {spectrum_file.name}: {detected_nucleus}")

                    # Get optimal parameters
                    parameters = self.get_optimal_parameters(detected_nucleus, spectrum_file)

                    # Apply auto-optimization if requested
                    if auto_optimize or self.config.get('optimization', {}).get('enable_auto_optimization', False):
                        parameters = self._apply_auto_optimization(spectrum_file, detected_nucleus, parameters)

                    # Process the spectrum
                    result = self.process_single_spectrum(spectrum_file, parameters)

                    # Add enhanced detection metadata to result
                    if detection_details:
                        result['nucleus_detection'] = {
                            'detected_nucleus': detected_nucleus,
                            'confidence': detection_details['confidence'],
                            'method': detection_details['primary_method'],
                            'chemical_shifts': detection_details.get('chemical_shifts'),
                            'all_methods': [m['method'] for m in detection_details['all_methods']]
                        }
                    else:
                        result['nucleus_detection'] = {
                            'detected_nucleus': detected_nucleus,
                            'confidence': 0.8 if nucleus_type else 0.0,  # User-provided vs unknown
                            'method': 'user_provided' if nucleus_type else 'filename_fallback',
                            'chemical_shifts': None,
                            'all_methods': []
                        }

                    # Update statistics
                    if result['success']:
                        self.stats['processed_files'] += 1
                        self.stats['total_peaks_detected'] += result['peaks_detected']
                        self.stats['total_peaks_fitted'] += result['peaks_fitted']
                        self.stats['total_ml_samples'] += result['ml_samples_collected']
                    else:
                        self.stats['failed_files'] += 1
                        self.stats['failed_files_list'].append({
                            'filename': result['filename'],
                            'error': result['error']
                        })

                except Exception as e:
                    self.logger.error(f"Unexpected error processing {spectrum_file.name}: {e}")
                    self.stats['failed_files'] += 1
                    self.stats['failed_files_list'].append({
                        'filename': spectrum_file.name,
                        'error': str(e)
                    })

                    if not self.config['skip_on_error']:
                        raise

            self.stats['processing_end_time'] = datetime.now()
            self._generate_summary_report()

            return self.stats

        except Exception as e:
            self.logger.error(f"Batch processing failed: {e}")
            self.stats['processing_end_time'] = datetime.now()
            raise

    def _generate_summary_report(self):
        """Generate and log a summary report of the batch processing."""
        if self.stats['processing_start_time'] and self.stats['processing_end_time']:
            processing_time = self.stats['processing_end_time'] - self.stats['processing_start_time']
            processing_seconds = processing_time.total_seconds()
        else:
            processing_seconds = 0

        self.logger.info("\n" + "="*60)
        self.logger.info("BATCH PROCESSING SUMMARY")
        self.logger.info("="*60)
        self.logger.info(f"Total files found: {self.stats['total_files']}")
        self.logger.info(f"Successfully processed: {self.stats['processed_files']}")
        self.logger.info(f"Failed processing: {self.stats['failed_files']}")
        self.logger.info(f"Total peaks detected: {self.stats['total_peaks_detected']}")
        self.logger.info(f"Total peaks fitted: {self.stats['total_peaks_fitted']}")

        # Get actual ML samples from collector (session-wide count)
        actual_ml_samples = len(self.integrator.ml_data_collector.session_data) if \
                           (hasattr(self, 'integrator') and hasattr(self.integrator, 'ml_data_collector') and
                            self.integrator.ml_data_collector) else 0

        self.logger.info(f"ML training samples collected: {actual_ml_samples}")
        self.logger.info(f"Processing time: {processing_seconds:.1f} seconds")

        if self.stats['failed_files'] > 0:
            self.logger.info(f"\nFailed files ({self.stats['failed_files']}):")
            for failed_file in self.stats['failed_files_list']:
                self.logger.info(f"  - {failed_file['filename']}: {failed_file['error']}")

        success_rate = (self.stats['processed_files'] / self.stats['total_files'] * 100) if self.stats['total_files'] > 0 else 0
        self.logger.info(f"\nSuccess rate: {success_rate:.1f}%")

        if actual_ml_samples > 0:
            self.logger.info(f"ML training data generation: SUCCESSFUL")
            self.logger.info(f"Ready for Phase 2 ML model development")
        else:
            self.logger.warning(f"No ML training samples collected - check spectrum quality")

        self.logger.info("="*60)

    def get_processing_statistics(self) -> Dict[str, Any]:
        """
        Get current processing statistics.

        Returns:
            Dictionary with current processing statistics
        """
        return self.stats.copy()

    def _detect_nucleus_type(self, file_path: Path) -> str:
        """
        Detect nucleus type from file path using enhanced detection methods.

        Args:
            file_path: Path to the spectrum file

        Returns:
            Detected nucleus type string
        """
        try:
            # Load spectrum temporarily for chemical shift analysis
            if self._load_spectrum(file_path):
                chemical_shift_ranges = self._extract_chemical_shift_ranges_from_integrator()
                if chemical_shift_ranges and 'f1_range' in chemical_shift_ranges and 'f2_range' in chemical_shift_ranges:
                    # Use enhanced detection with chemical shifts
                    detection_result = self.detect_nucleus_type(
                        file_path.name,
                        f1_range=chemical_shift_ranges['f1_range'],
                        f2_range=chemical_shift_ranges['f2_range']
                    )
                    return detection_result

            # Fallback to filename detection
            return self.detect_nucleus_type(file_path.name)
        except Exception as e:
            self.logger.debug(f"Nucleus detection failed for {file_path}: {e}")
            return 'unknown'

    def _get_nucleus_parameters(self, nucleus_type: str) -> Dict[str, Any]:
        """
        Get processing parameters for a specific nucleus type.

        Args:
            nucleus_type: Type of nucleus ('15N1H', '13C1H', etc.)

        Returns:
            Dictionary with nucleus-specific parameters
        """
        # Get S/N threshold for this nucleus type
        sn_thresholds = self.config.get('sn_thresholds', {})
        sn_threshold = sn_thresholds.get(nucleus_type, sn_thresholds.get('default', 2.2))

        # Get expected peak count for this nucleus type
        expected_peaks = self.config.get('expected_peaks', {})
        expected_peak_count = expected_peaks.get(nucleus_type, expected_peaks.get('default', 50))

        return {
            'nucleus_type': nucleus_type,
            'sn_threshold': float(sn_threshold),
            'expected_peak_count': int(expected_peak_count),
            'quality_threshold': float(self.config.get('quality_threshold', 0.8)),
            'max_attempts': int(self.config.get('max_fitting_attempts', 3)),
            'max_iterations': int(self.config.get('max_iterations', 1000)),
            'convergence_tolerance': float(self.config.get('convergence_tolerance', 1e-8))
        }

    def _prepare_processing_parameters(self, nucleus_type: str) -> Dict[str, Any]:
        """
        Prepare complete processing parameters for a nucleus type.

        Args:
            nucleus_type: Type of nucleus

        Returns:
            Complete parameter dictionary for processing
        """
        base_params = self._get_nucleus_parameters(nucleus_type)

        # Add additional processing parameters
        processing_params = {
            **base_params,
            'skip_on_error': self.config.get('skip_on_error', True),
            'max_optimization_iterations': int(self.config.get('max_optimization_iterations', 50)),
            'baseline_max_iter': int(self.config.get('baseline_max_iter', 50)),
            'detailed_logging': self.config.get('detailed_logging', True)
        }

        return processing_params

if __name__ == "__main__":
    # Basic command-line usage
    import argparse

    parser = argparse.ArgumentParser(description="Batch process NMR spectra for ML training")
    parser.add_argument("folder", help="Path to folder containing NMR spectra")
    parser.add_argument("--nucleus", choices=['15N1H', '13C1H', '1H'],
                       help="Nucleus type (if not specified, will auto-detect)")
    parser.add_argument("--auto-optimize", action="store_true",
                       help="Enable automatic parameter optimization")

    args = parser.parse_args()

    # Create and run batch processor
    processor = BatchProcessor()
    try:
        results = processor.process_folder(args.folder, args.nucleus, args.auto_optimize)
        print(f"\nBatch processing completed successfully!")
        print(f"Processed {results['processed_files']}/{results['total_files']} files")
        print(f"ML training samples collected: {results['total_ml_samples']}")
    except Exception as e:
        print(f"Batch processing failed: {e}")
        sys.exit(1)
