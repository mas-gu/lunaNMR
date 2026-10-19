#!/usr/bin/env python3
"""
Core Batch Processing Engine for LunaNMR

This module provides the core functionality for batch processing multiple NMR spectra
for automated analysis. It operates independently from the GUI interface.

Features:
- Automated spectrum loading and processing
- Robust error handling with individual spectrum failure recovery
- Progress tracking and detailed logging
- Support for multiple nucleus types (1H, 15N, 13C)
"""

import sys
import time
import logging
from pathlib import Path
from typing import List, Dict, Optional, Union, Any
import numpy as np
from datetime import datetime

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
    including S/N detection and peak fitting.
    It operates completely independently from the GUI interface.
    """

    def __init__(self, config: Optional[Dict] = None, output_dir: Union[str, Path, None] = None):
        """
        Initialize the batch processor.

        Args:
            config: Optional configuration dictionary with processing parameters
            output_dir: Folder for the batch log and each spectrum's fitting results CSV
                (default: the current directory)
        """
        self.output_dir = Path(output_dir) if output_dir else Path.cwd()
        # Set up logging
        self.logger = logging.getLogger(__name__)
        self._setup_logging(self.output_dir)

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
            # Search/fitting window overrides; without them a spectrum is fitted with
            # the windows the GUI's Fit Spectrum uses
            'processing_parameters': {}
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
            'processing_start_time': None,
            'processing_end_time': None,
            'failed_files_list': [],
            # One entry per spectrum: what was detected, fitted and written
            'files': []
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

            self.logger.info("Batch processor initialized")
        else:
            self.integrator = None
            self.logger.error("Core components not available - batch processing disabled")

    def _setup_logging(self, log_dir: Path):
        """Set up logging system for batch processing, logging to a file in `log_dir`.

        The logger is shared by every processor in the process, so the handlers a previous
        processor added are closed first: otherwise each line is printed once per processor
        and every earlier log file stays open.
        """
        # Create logger
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        for handler in list(logger.handlers):
            logger.removeHandler(handler)
            handler.close()

        # Create console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)

        # Create file handler
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / f"batch_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
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
        except Exception:
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
            'nucleus_type': parameters.get('nucleus_type'),
            'sn_threshold': parameters['sn_threshold'],
            'expected_peak_count': parameters['expected_peak_count'],
            'peaks_detected': 0,
            'peaks_fitted': 0,
            'peaks_above_quality_threshold': 0,
            'results_file': None,
            'fitted_results': [],
            'processing_time': 0,
            'error': None
        }

        start_time = time.time()

        try:
            self.logger.info(f"Processing: {spectrum_file.name}")

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

                    # Step 2: Fit the detected peaks, and write their results as the GUI does
                    fitted_results = self._fit_detected_peaks(detected_peaks)
                    result['fitted_results'] = fitted_results
                    if fitted_results:
                        result['results_file'] = str(self._write_fitting_results(spectrum_file, fitted_results))

                    fitted = [r for r in fitted_results if r.get('success')]
                    quality_threshold = parameters.get('quality_threshold', 0.0)
                    result['peaks_fitted'] = len(fitted)
                    result['peaks_above_quality_threshold'] = sum(
                        1 for r in fitted if (r.get('avg_r_squared') or r.get('r_squared') or 0) >= quality_threshold)

                    if fitted:
                        result['success'] = True
                        self.logger.info(f"  Fitted {len(fitted)}/{len(detected_peaks)} peaks "
                                         f"({result['peaks_above_quality_threshold']} with R² >= {quality_threshold})")
                    else:
                        result['error'] = f"No peak fitted ({len(detected_peaks)} detected)"
                        self.logger.warning(f"  No peak fitted in {spectrum_file.name}")

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

    def _write_fitting_results(self, spectrum_file: Path, fitted_results: List[Dict]) -> Path:
        """Write the spectrum's fit results into the output folder as the GUI's
        `<spectrum>_fitting_results.csv`."""
        from lunaNMR.processors.single_spectrum_processor import write_fitting_results_csv

        self.output_dir.mkdir(parents=True, exist_ok=True)
        output_file = self.output_dir / f"{spectrum_file.stem}_fitting_results.csv"
        write_fitting_results_csv(fitted_results, output_file)
        self.logger.info(f"  Saved fitting results: {output_file}")
        return output_file

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
            elif hasattr(self.integrator, 'load_nmr_file'):
                # Fallback to load_nmr_file method
                success = self.integrator.load_nmr_file(str(spectrum_file))
            else:
                self.logger.error("No suitable spectrum loading method available")
                return False
            if success:
                self._apply_fit_all_parameters()
            return success

        except Exception as e:
            self.logger.error(f"Failed to load {spectrum_file.name}: {e}")
            return False

    def _fit_all_parameter_manager(self, nucleus: str):
        """
        The parameters the GUI's Fit Spectrum fits with, for a spectrum of `nucleus`.

        These are the main window's defaults as its sync_parameters_to_integrator writes
        them into its parameter manager, and the PS2D configuration of the nucleus it
        auto-detects on load. A config's `processing_parameters` replace the search and
        fitting windows, and its `max_iterations` the fit's iteration limit.
        """
        from lunaNMR.core.ps2d_config import set_ps2d_config
        from lunaNMR.utils.parameter_manager import NMRParameterManager

        ps2d_config = set_ps2d_config(nucleus)
        param_manager = NMRParameterManager()
        param_manager.current_params.update({
            'search_window_x': ps2d_config.search_window_f2,
            'search_window_y': ps2d_config.search_window_f1,
            'centroid_window_x_ppm': 0.01,
            'centroid_window_y_ppm': 0.1,
            'centroid_noise_multiplier': 2.0,
            'fix_linewidths': False,
            'fix_positions': False,
            'use_parallel_processing': True,
            'use_adaptive_optimization': False,
            'lw_lorentz_1h': None,
            'lw_gauss_1h': None,
            'lw_lorentz_15n': None,
            'lw_gauss_15n': None,
            'max_iterations': int(self.config.get('max_iterations', 1000)),
        })

        processing_overrides = self.config.get('processing_parameters', {})
        for key in ['search_window_x', 'search_window_y', 'fitting_window_x', 'fitting_window_y']:
            if key in processing_overrides:
                try:
                    param_manager.current_params[key] = float(processing_overrides[key])
                except (TypeError, ValueError):
                    self.logger.warning(f"Invalid numeric override for {key}: {processing_overrides[key]}")
        return param_manager

    def _apply_fit_all_parameters(self):
        """
        Put the GUI's Fit Spectrum parameters on the integrator for the loaded spectrum.

        The GUI syncs them before it detects as well as before it fits: the detection's
        centroid windows are among them.
        """
        from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor

        nucleus = getattr(self.integrator, 'auto_detected_nucleus', None) or '15N'
        self.fit_processor = SingleSpectrumProcessor(
            self.integrator, self._fit_all_parameter_manager(nucleus))
        self.fit_processor.sync_parameters_to_integrator()

    def _fit_detected_peaks(self, detected_peaks: List[Dict]) -> List[Dict]:
        """
        Fit the detected peaks as the GUI's Fit Spectrum does after a detection: the
        detected peaks are the peak list, clustered afresh and fitted by
        SingleSpectrumProcessor.process_peak_list.

        Returns:
            One fit result per peak, as the GUI's fit returns them
        """
        from lunaNMR.processors.single_spectrum_processor import detected_peak_list

        peak_list = detected_peak_list(detected_peaks)
        self.logger.info(f"  Fitting {len(peak_list)} detected peaks")

        # Fit Spectrum clusters the peaks it is given, never the clusters of a fit before
        if self.integrator.enhanced_fitter is not None:
            self.integrator.enhanced_fitter.series_params = None

        processing_options = {
            'use_parallel': self.fit_processor.parameter_manager.current_params['use_parallel_processing'],
            'use_global_optimization': False,
            'use_voigt_fitting': True,
            'skip_series_params': True
        }
        fitted_results, _ = self.fit_processor.process_peak_list(peak_list, processing_options)
        return fitted_results

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
                    self.logger.warning("  ⚠️ Could not load spectrum for optimization, using base parameters")
                    return base_parameters

                # Check again after loading
                if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
                    self.logger.warning("  ⚠️ No spectrum data available for optimization, using base parameters")
                    return base_parameters

            # Extract spectrum characteristics
            spectrum_data = self.integrator.nmr_data
            f1_ppm = getattr(self.integrator, 'ppm_y_axis', None)  # Y-axis is F1 (15N/13C)
            f2_ppm = getattr(self.integrator, 'ppm_x_axis', None)  # X-axis is F2 (1H)

            if f1_ppm is None or f2_ppm is None:
                self.logger.warning("  ⚠️ Missing PPM scales for optimization, using base parameters")
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
                    self.stats['files'].append({key: value for key, value in result.items()
                                                if key not in ('fitted_results', 'nucleus_detection')})
                    if result['success']:
                        self.stats['processed_files'] += 1
                        self.stats['total_peaks_detected'] += result['peaks_detected']
                        self.stats['total_peaks_fitted'] += result['peaks_fitted']
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
        self.logger.info(f"Processing time: {processing_seconds:.1f} seconds")

        if self.stats['failed_files'] > 0:
            self.logger.info(f"\nFailed files ({self.stats['failed_files']}):")
            for failed_file in self.stats['failed_files_list']:
                self.logger.info(f"  - {failed_file['filename']}: {failed_file['error']}")

        success_rate = (self.stats['processed_files'] / self.stats['total_files'] * 100) if self.stats['total_files'] > 0 else 0
        self.logger.info(f"\nSuccess rate: {success_rate:.1f}%")
        self.logger.info("="*60)

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

if __name__ == "__main__":
    # Basic command-line usage
    import argparse

    parser = argparse.ArgumentParser(description="Batch process NMR spectra")
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
        print("\nBatch processing completed successfully!")
        print(f"Processed {results['processed_files']}/{results['total_files']} files")
    except Exception as e:
        print(f"Batch processing failed: {e}")
        sys.exit(1)
