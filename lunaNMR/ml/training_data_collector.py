#!/usr/bin/env python3
"""
Enhanced ML Training Data Collection Module for LunaNMR - Version 2.0

This module transparently collects comprehensive training data from successful peak fits
for advanced ML model development. Operates completely in the background
without affecting user workflows or performance.

Enhancements in Version 2.0:
- Initial parameter estimates and optimization trajectory data
- Multi-method fitting results and consensus analysis
- Physics-based parameter validation and constraints
- Advanced quality metrics with residual analysis
- Optimization difficulty and convergence information
- Enhanced spectral analysis (noise, baseline, peak shape)
- Chemical environment and coupling predictions
- Cross-validation and uncertainty quantification

Features:
- Quality-filtered data collection (R² > 0.8)
- Comprehensive feature extraction (spectral, chemical, context, optimization)
- Seamless integration with all fitting methods
- Graceful error handling with silent operation
- Automatic batch saving and storage management
- Rich optimization context for advanced ML models
"""

import numpy as np
import pickle
import json
import os
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any
import logging
from pathlib import Path
from scipy import stats

class MLTrainingDataCollector:
    """
    Enhanced ML training data collector with comprehensive feature extraction.

    Captures not just the final fit results, but the entire optimization context:
    - Initial estimates vs final parameters
    - Convergence behavior and optimization difficulty
    - Multi-method comparison results
    - Physics-based parameter validation
    - Advanced spectral and chemical analysis
    """

    def __init__(self, storage_path: Optional[str] = None):
        """
        Initialize enhanced ML training data collector.

        Args:
            storage_path: Directory for storing collected data
                         Defaults to 'ml_training_data/' in project root
        """
        # Set up storage path
        if storage_path is None:
            # Default to project root / ml_training_data
            project_root = Path(__file__).parent.parent.parent.parent
            storage_path = project_root / "ml_training_data"

        self.storage_path = Path(storage_path)
        self.storage_path.mkdir(exist_ok=True)

        # Enhanced collection configuration v2.0
        self.collection_version = '2.0'
        self.enhanced_features_enabled = True
        self.collect_rejected_samples = True  # NEW: Collect hard cases for ML learning

        # Session data buffer
        self.session_data = []
        self.batch_size = 100  # Save after collecting 100 samples

        # Initialize logging
        self.logger = logging.getLogger(__name__)

        # Track window settings used for data collection (optional metadata)
        self.window_settings = {}

        # Quality thresholds
        self.min_r_squared = 0.8  # Threshold for "good" vs "reject" labels
        self.min_r_squared_absolute = 0.3  # Absolute minimum to collect (even rejects)
        self.max_samples_per_session = 1000  # Prevent memory issues

        # Physics constraints for validation
        self.physics_constraints = {
            '1H': {'width_range': (0.001, 0.1), 'shift_range': (0, 15)},
            '15N': {'width_range': (0.1, 10), 'shift_range': (80, 180)},
            '13C': {'width_range': (0.1, 5), 'shift_range': (0, 220)}
        }

        # Enhanced feature schema v2.0
        self.feature_schema = {
            'spectral_features': [
                'peak_height', 'baseline_level', 'dynamic_range', 'snr', 'spectral_range',
                'data_points', 'peak_width_fwhm', 'peak_width_moment', 'peak_width_10percent',
                'peak_width_50percent', 'peak_width_90percent', 'peak_width_ratio_90_10',
                'peak_asymmetry', 'peak_fronting', 'peak_tailing', 'baseline_drift',
                'baseline_slope', 'baseline_curvature', 'noise_variance', 'noise_correlation',
                'noise_level_robust', 'intensity_mean', 'intensity_std', 'intensity_skewness',
                'intensity_kurtosis', 'spectral_resolution', 'data_density'
            ],
            'chemical_features': [
                'nucleus_type', 'chemical_shift_center', 'chemical_region', 'shift_start',
                'shift_end', 'shift_range', 'expected_multiplicity', 'coupling_environment',
                'relaxation_estimate', 'chemical_shift_anisotropy', 'j_coupling_estimate'
            ],
            'context_features': [
                'detection_confidence', 'estimated_amplitude', 'estimated_width',
                'nearby_peaks_count', 'peak_complexity', 'noise_level',
                'overlap_severity', 'interference_sources', 'peak_isolation_score',
                'baseline_stability', 'data_quality_score'
            ],
            'optimization_features': [
                'initial_parameters', 'parameter_bounds_used', 'convergence_iterations',
                'optimization_method', 'converged', 'parameter_uncertainties',
                'gradient_norm', 'hessian_condition', 'optimization_difficulty',
                'parameter_correlation_matrix', 'cost_function_landscape',
                'initial_cost', 'final_cost', 'parameter_changes'
            ],
            'multi_method_results': [
                'methods_tested', 'best_method_r_squared', 'worst_method_r_squared',
                'consensus_parameters', 'parameter_stability', 'method_agreement_score',
                'robust_parameter_estimates', 'r_squared_range'
            ],
            'physics_validation': [
                'parameter_validity', 'spectroscopic_constraints',
                'shape_realism_score', 'linewidth_consistency', 'chemical_shift_validity'
            ],
            'target_parameters': [
                'amplitude', 'center', 'sigma', 'gamma', 'baseline',
                'parameter_uncertainties', 'confidence_intervals'
            ],
            'quality_metrics': [
                'r_squared', 'rmse_normalized', 'parameters_physical', 'fit_success',
                'residual_analysis', 'outlier_detection', 'cross_validation_score',
                'goodness_of_fit_tests', 'model_selection_criteria'
            ]
        }

        # Enhanced statistics tracking
        self.collection_stats = {
            'total_attempts': 0,
            'quality_filtered': 0,
            'successfully_collected': 0,
            'batch_saves': 0,
            'errors': 0,
            'convergence_failures': 0,
            'physics_violations': 0,
            'multi_method_comparisons': 0,
            'enhanced_features_extracted': 0,
            'optimization_trajectories_saved': 0,
            'hard_cases_collected': 0  # NEW: Track rejected samples for ML learning
        }

        self.logger.info(f"Enhanced ML Training Data Collector v2.0 initialized - Storage: {self.storage_path}")

    def collect_training_sample(self,
                               x_data: np.ndarray,
                               y_data: np.ndarray,
                               fit_result: Dict[str, Any],
                               context: Dict[str, Any],
                               nucleus_type: str,
                               initial_params: Optional[Dict[str, Any]] = None,
                               optimization_info: Optional[Dict[str, Any]] = None,
                               alternative_results: Optional[List[Dict]] = None) -> bool:
        """
        Collect enhanced training sample with comprehensive optimization context.

        Args:
            x_data: X-axis spectral data (ppm)
            y_data: Y-axis intensity data
            fit_result: Final fitting result dictionary
            context: Detection and spectral context
            nucleus_type: '1H', '15N', '13C', etc.
            initial_params: Initial parameter estimates (Enhanced v2.0)
            optimization_info: Convergence and optimization details (Enhanced v2.0)
            alternative_results: Results from alternative fitting methods (Enhanced v2.0)

        Returns:
            True if sample was collected successfully
        """
        self.collection_stats['total_attempts'] += 1

        try:
            # Enhanced quality handling - collect both good and bad fits with labels
            r_squared = fit_result.get('r_squared', 0.0)

            # Determine quality label
            if r_squared >= self.min_r_squared:
                quality_label = 'good'
            elif r_squared >= self.min_r_squared_absolute and self.collect_rejected_samples:
                quality_label = 'reject'  # Hard cases for ML learning
            else:
                # Either truly terrible fits OR rejected samples disabled - skip
                self.collection_stats['quality_filtered'] += 1
                return False

            # Convergence validation (still collect difficult cases for learning)
            converged = optimization_info.get('converged', True) if optimization_info else True
            if not converged:
                self.collection_stats['convergence_failures'] += 1

            # Prevent memory issues
            if len(self.session_data) >= self.max_samples_per_session:
                self.logger.warning("Max samples per session reached, skipping collection")
                return False

            # Extract enhanced features with error handling
            try:
                spectral_features = self._extract_enhanced_spectral_features(x_data, y_data)
                chemical_features = self._extract_enhanced_chemical_features(x_data, nucleus_type)
                context_features = self._extract_enhanced_context_features(context, x_data, y_data)

                # Enhanced features (v2.0)
                optimization_features = {}
                if initial_params or optimization_info:
                    optimization_features = self._extract_optimization_features(
                        initial_params or {}, fit_result, optimization_info or {}
                    )
                    self.collection_stats['optimization_trajectories_saved'] += 1

                multi_method_features = {}
                if alternative_results:
                    multi_method_features = self._analyze_multi_method_results(fit_result, alternative_results)
                    self.collection_stats['multi_method_comparisons'] += 1

                physics_features = self._extract_physics_validation(fit_result, nucleus_type)
                quality_metrics = self._compute_enhanced_quality_metrics(
                    x_data, y_data, fit_result, optimization_info or {}
                )

                self.collection_stats['enhanced_features_extracted'] += 1

            except Exception as e:
                self.logger.debug(f"Enhanced feature extraction failed: {e}")
                self.collection_stats['errors'] += 1
                return False

            # Validate target parameters
            required_params = ['amplitude', 'center', 'sigma', 'gamma', 'baseline']
            for param in required_params:
                if param not in fit_result or not np.isfinite(fit_result[param]):
                    self.logger.debug(f"Invalid target parameter {param}: {fit_result.get(param, 'missing')}")
                    self.collection_stats['errors'] += 1
                    return False

            # Create enhanced training sample
            training_sample = {
                'timestamp': datetime.now().isoformat(),
                'collection_version': '2.0',
                'sample_id': f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{len(self.session_data)}",

                # Enhanced features
                'spectral_features': spectral_features,
                'chemical_features': chemical_features,
                'context_features': context_features,
                'optimization_features': optimization_features,
                'multi_method_results': multi_method_features,
                'physics_validation': physics_features,

                # Target parameters with uncertainty
                'target_parameters': {
                    'amplitude': float(fit_result['amplitude']),
                    'center': float(fit_result['center']),
                    'sigma': float(fit_result['sigma']),
                    'gamma': float(fit_result['gamma']),
                    'baseline': float(fit_result['baseline']),
                    'parameter_uncertainties': optimization_info.get('parameter_errors', {}) if optimization_info else {}
                },

                # Enhanced quality metrics with ML learning label
                'quality_metrics': {**quality_metrics, 'label': quality_label},

                # Optimization trajectory (if available)
                'optimization_trajectory': {
                    'parameter_evolution': optimization_info.get('parameter_history', []) if optimization_info else [],
                    'cost_evolution': optimization_info.get('cost_history', []) if optimization_info else [],
                    'gradient_evolution': optimization_info.get('gradient_history', []) if optimization_info else []
                },

                # Raw data with preprocessing info
                'raw_data': {
                    'x_data': self._downsample_data(x_data).tolist(),
                    'y_data': self._downsample_data(y_data).tolist(),
                    'original_length': len(x_data),
                    'preprocessing_applied': context.get('preprocessing_steps', []),
                    'data_quality_flags': self._assess_data_quality_flags(x_data, y_data)
                }
            }

            # Add to session buffer with label tracking
            self.session_data.append(training_sample)
            self.collection_stats['successfully_collected'] += 1
            if quality_label == 'reject':
                self.collection_stats['hard_cases_collected'] = self.collection_stats.get('hard_cases_collected', 0) + 1

            # Trigger batch save if threshold reached
            if len(self.session_data) >= self.batch_size:
                self._save_batch()

            return True

        except Exception as e:
            # Silent error handling - never break workflows
            self.logger.debug(f"Enhanced ML data collection failed: {e}")
            self.collection_stats['errors'] += 1
            return False

    def update_window_settings(self, settings: Dict[str, Any]) -> None:
        """Record the detection/fitting window settings for metadata export."""

        if not isinstance(settings, dict):
            return

        sanitized: Dict[str, Any] = {}
        for key in ['search_window_x', 'search_window_y', 'fitting_window_x', 'fitting_window_y']:
            if key in settings and settings[key] is not None:
                try:
                    sanitized[key] = float(settings[key])
                except (TypeError, ValueError):
                    self.logger.debug(f"Invalid window setting for {key}: {settings[key]}")

        if sanitized:
            self.window_settings.update(sanitized)

    def _extract_enhanced_spectral_features(self, x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, float]:
        """Extract comprehensive spectral characteristics as ML features."""

        # Basic spectral statistics
        peak_height = float(np.max(y_data))
        baseline_level = float(np.percentile(y_data, 10))
        dynamic_range = peak_height - baseline_level

        # Advanced noise analysis
        noise_analysis = self._analyze_noise_characteristics(y_data)

        # Signal-to-noise ratio (robust estimation)
        signal = peak_height - baseline_level
        noise_region = y_data[y_data < np.percentile(y_data, 25)]
        noise = float(np.std(noise_region)) if len(noise_region) > 1 else signal * 0.01
        snr = signal / (noise + 1e-8)

        # Enhanced peak width estimation (multiple methods)
        peak_width_estimates = self._estimate_peak_width_multiple_methods(x_data, y_data)

        # Peak asymmetry and tailing analysis
        asymmetry_metrics = self._analyze_peak_asymmetry(x_data, y_data)

        # Baseline characteristics
        baseline_metrics = self._analyze_baseline_characteristics(x_data, y_data)

        # Spectral resolution and density
        if len(x_data) > 1:
            spectral_resolution = float(abs(x_data[1] - x_data[0]))
            data_density = float(len(x_data) / abs(x_data[-1] - x_data[0]))
        else:
            spectral_resolution = 0.001
            data_density = 1000.0

        # Intensity distribution analysis
        intensity_metrics = self._analyze_intensity_distribution(y_data)

        # Combine all spectral features
        spectral_features = {
            'peak_height': peak_height,
            'baseline_level': baseline_level,
            'dynamic_range': dynamic_range,
            'snr': float(snr),
            'spectral_range': float(x_data[-1] - x_data[0]) if len(x_data) > 1 else 0.0,
            'data_points': int(len(x_data)),
            'spectral_resolution': spectral_resolution,
            'data_density': data_density,
            **noise_analysis,
            **peak_width_estimates,
            **asymmetry_metrics,
            **baseline_metrics,
            **intensity_metrics
        }

        return spectral_features

    def _extract_enhanced_chemical_features(self, x_data: np.ndarray, nucleus_type: str) -> Dict[str, Any]:
        """Extract comprehensive chemical context features."""

        chemical_shift_center = float(np.mean(x_data))
        shift_range = float(x_data[-1] - x_data[0]) if len(x_data) > 1 else 0.0

        # Chemical region classification (enhanced)
        region = self._classify_chemical_region_enhanced(chemical_shift_center, nucleus_type)

        # Expected multiplicity prediction
        multiplicity = self._predict_multiplicity(chemical_shift_center, nucleus_type)

        # Coupling environment assessment
        coupling_env = self._assess_coupling_environment(chemical_shift_center, nucleus_type)

        # Relaxation time estimates
        relaxation_estimates = self._estimate_relaxation_properties(chemical_shift_center, nucleus_type)

        # Chemical shift anisotropy estimation
        csa_estimate = self._estimate_chemical_shift_anisotropy(chemical_shift_center, nucleus_type)

        # J-coupling estimates
        j_coupling_estimate = self._estimate_j_coupling(chemical_shift_center, nucleus_type)

        return {
            'nucleus_type': nucleus_type,
            'chemical_shift_center': chemical_shift_center,
            'shift_range': shift_range,
            'chemical_region': region,
            'shift_start': float(x_data[0]) if len(x_data) > 0 else 0.0,
            'shift_end': float(x_data[-1]) if len(x_data) > 0 else 0.0,
            'expected_multiplicity': multiplicity,
            'coupling_environment': coupling_env,
            'relaxation_estimate': relaxation_estimates,
            'chemical_shift_anisotropy': csa_estimate,
            'j_coupling_estimate': j_coupling_estimate
        }

    def _extract_enhanced_context_features(self, context: Dict[str, Any], x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, Any]:
        """Extract comprehensive contextual features from fitting context."""

        # Basic context features
        detection_confidence = float(context.get('detection_confidence', 0.5))
        estimated_amplitude = float(context.get('estimated_amplitude', np.max(y_data) - np.min(y_data)))
        estimated_width = float(context.get('estimated_width', 0.02))

        # Enhanced overlap analysis
        overlap_analysis = self._analyze_peak_overlap_enhanced(context, x_data, y_data)

        # Interference source identification
        interference_sources = self._identify_interference_sources(context)

        # Peak isolation scoring
        isolation_score = self._calculate_peak_isolation_score(context, x_data)

        # Baseline stability assessment
        baseline_stability = self._assess_baseline_stability_enhanced(y_data)

        # Data quality scoring
        data_quality = self._assess_data_quality_comprehensive(x_data, y_data, context)

        # Overlapping peaks information
        overlapping_peaks = context.get('overlapping_peaks', [])
        nearby_peaks_count = len(overlapping_peaks) if overlapping_peaks else 0

        # Peak complexity classification (enhanced)
        complexity = self._classify_peak_complexity_enhanced(context, nearby_peaks_count)

        # Noise level estimation (advanced)
        noise_level = self._estimate_noise_level_advanced(y_data, context)

        return {
            'detection_confidence': detection_confidence,
            'estimated_amplitude': estimated_amplitude,
            'estimated_width': estimated_width,
            'nearby_peaks_count': int(nearby_peaks_count),
            'peak_complexity': complexity,
            'noise_level': noise_level,
            'overlap_severity': overlap_analysis.get('severity', 0.0),
            'interference_sources': interference_sources,
            'peak_isolation_score': isolation_score,
            'baseline_stability': baseline_stability,
            'data_quality_score': data_quality
        }

    def _extract_optimization_features(self, initial_params: Dict[str, Any],
                                     final_result: Dict[str, Any],
                                     optimization_info: Dict[str, Any]) -> Dict[str, Any]:
        """Extract optimization and convergence features."""

        features = {}

        # Initial vs final parameter comparison
        param_changes = {}
        for param in ['amplitude', 'center', 'sigma', 'gamma', 'baseline']:
            if param in initial_params and param in final_result:
                initial_val = float(initial_params[param])
                final_val = float(final_result[param])
                if abs(initial_val) > 1e-10:
                    param_changes[f'{param}_relative_change'] = abs(final_val - initial_val) / abs(initial_val)
                    param_changes[f'{param}_improvement_ratio'] = final_val / initial_val
                else:
                    param_changes[f'{param}_relative_change'] = 0.0
                    param_changes[f'{param}_improvement_ratio'] = 1.0

                param_changes[f'{param}_absolute_change'] = abs(final_val - initial_val)

        features['parameter_changes'] = param_changes
        features['initial_parameters'] = initial_params

        # Parameter bounds information
        features['parameter_bounds_used'] = optimization_info.get('bounds_used', {})

        # Convergence metrics
        features['convergence_iterations'] = int(optimization_info.get('iterations', 0))
        features['converged'] = bool(optimization_info.get('converged', False))
        features['optimization_method'] = optimization_info.get('method', 'unknown')
        features['final_cost'] = float(optimization_info.get('final_cost', 0.0))
        features['initial_cost'] = float(optimization_info.get('initial_cost', features['final_cost']))

        # Parameter uncertainty estimates
        features['parameter_uncertainties'] = optimization_info.get('parameter_errors', {})

        # Optimization difficulty assessment
        features['optimization_difficulty'] = self._assess_optimization_difficulty(optimization_info)

        # Gradient and Hessian information
        features['gradient_norm'] = float(optimization_info.get('gradient_norm', 0.0))
        features['hessian_condition'] = float(optimization_info.get('hessian_condition', 1.0))

        # Parameter correlation matrix (if available)
        features['parameter_correlation_matrix'] = optimization_info.get('correlation_matrix', {})

        # Cost function landscape information
        features['cost_function_landscape'] = {
            'local_minima_detected': optimization_info.get('local_minima', 0),
            'cost_function_curvature': optimization_info.get('curvature', 0.0),
            'optimization_path_length': len(optimization_info.get('parameter_history', []))
        }

        return features

    def _analyze_multi_method_results(self, primary_result: Dict[str, Any],
                                    alternative_results: List[Dict]) -> Dict[str, Any]:
        """Analyze results from multiple fitting methods."""

        analysis = {}

        # Method comparison
        methods_tested = [result.get('method', 'unknown') for result in alternative_results]
        analysis['methods_tested'] = methods_tested

        # Parameter consensus analysis
        consensus_params = {}
        stability_metrics = {}

        for param in ['amplitude', 'center', 'sigma', 'gamma', 'baseline']:
            if param in primary_result:
                values = [primary_result[param]]
                values.extend([result.get(param, 0) for result in alternative_results if param in result])

                if len(values) > 1:
                    consensus_params[f'{param}_mean'] = float(np.mean(values))
                    consensus_params[f'{param}_std'] = float(np.std(values))
                    consensus_params[f'{param}_range'] = float(np.max(values) - np.min(values))
                    consensus_params[f'{param}_median'] = float(np.median(values))

                    # Parameter stability assessment
                    if consensus_params[f'{param}_mean'] != 0:
                        stability_metrics[f'{param}_cv'] = consensus_params[f'{param}_std'] / abs(consensus_params[f'{param}_mean'])
                    else:
                        stability_metrics[f'{param}_cv'] = 0.0

        analysis['consensus_parameters'] = consensus_params
        analysis['parameter_stability'] = stability_metrics

        # Method agreement scoring
        r_squares = [primary_result.get('r_squared', 0)]
        r_squares.extend([result.get('r_squared', 0) for result in alternative_results])

        analysis['method_agreement_score'] = float(np.std(r_squares)) if len(r_squares) > 1 else 0.0
        analysis['best_method_r_squared'] = float(max(r_squares))
        analysis['worst_method_r_squared'] = float(min(r_squares))
        analysis['r_squared_range'] = analysis['best_method_r_squared'] - analysis['worst_method_r_squared']

        # Robust parameter estimates (using median)
        robust_estimates = {}
        for param in ['amplitude', 'center', 'sigma', 'gamma', 'baseline']:
            if param in consensus_params:
                robust_estimates[f'{param}_robust'] = consensus_params[f'{param}_median']

        analysis['robust_parameter_estimates'] = robust_estimates

        return analysis

    def _extract_physics_validation(self, fit_result: Dict[str, Any], nucleus_type: str) -> Dict[str, Any]:
        """Validate parameters against spectroscopic physics."""

        validation = {}

        # Parameter range validation
        param_valid = {}
        nucleus_clean = nucleus_type.replace('1H', '1H').replace('15N', '15N').replace('13C', '13C')

        for param, value in fit_result.items():
            if param in ['amplitude', 'center', 'sigma', 'gamma', 'baseline']:
                param_valid[f'{param}_valid'] = self._validate_parameter_physics(param, value, nucleus_clean)

        validation['parameter_validity'] = param_valid

        # Peak shape realism
        sigma = fit_result.get('sigma', 0.01)
        gamma = fit_result.get('gamma', 0.01)
        validation['shape_realism_score'] = self._assess_peak_shape_realism(sigma, gamma, nucleus_clean)

        # Linewidth consistency
        validation['linewidth_consistency'] = self._check_linewidth_consistency(sigma, gamma, nucleus_clean)

        # Chemical shift validity
        center = fit_result.get('center', 0.0)
        validation['chemical_shift_validity'] = self._validate_chemical_shift(center, nucleus_clean)

        # Spectroscopic constraints
        validation['spectroscopic_constraints'] = self._check_spectroscopic_constraints(fit_result, nucleus_clean)

        return validation

    def _compute_enhanced_quality_metrics(self, x_data: np.ndarray, y_data: np.ndarray,
                                        fit_result: Dict[str, Any],
                                        optimization_info: Dict[str, Any]) -> Dict[str, Any]:
        """Compute comprehensive quality metrics."""

        metrics = {}

        # Basic quality metrics
        metrics['r_squared'] = float(fit_result.get('r_squared', 0.0))
        metrics['rmse_normalized'] = float(fit_result.get('rmse_normalized', 0.0))
        metrics['fit_success'] = bool(fit_result.get('success', False))

        # Enhanced residual analysis
        residuals = optimization_info.get('residuals', [])
        if len(residuals) > 0:
            residual_metrics = self._analyze_residuals_comprehensive(residuals)
            metrics.update(residual_metrics)

        # Outlier detection
        metrics['outlier_fraction'] = self._detect_outlier_fraction(y_data, fit_result)

        # Cross-validation estimate (if available)
        metrics['cross_validation_score'] = float(optimization_info.get('cv_score', 0.0))

        # Goodness of fit tests
        gof_tests = self._perform_goodness_of_fit_tests(y_data, fit_result, residuals)
        metrics['goodness_of_fit_tests'] = gof_tests

        # Model selection criteria
        n_params = 5  # Voigt model has 5 parameters
        n_data = len(y_data)
        if metrics['r_squared'] > 0:
            sse = np.sum(np.array(residuals)**2) if residuals else 0
            metrics['aic'] = n_data * np.log(sse/n_data) + 2 * n_params
            metrics['bic'] = n_data * np.log(sse/n_data) + n_params * np.log(n_data)

        return metrics

    # Helper methods for enhanced feature extraction

    def _analyze_noise_characteristics(self, y_data: np.ndarray) -> Dict[str, float]:
        """Comprehensive noise analysis."""

        # Robust noise estimation using MAD
        noise_robust = float(np.median(np.abs(y_data - np.median(y_data))) * 1.4826)

        # Noise autocorrelation
        if len(y_data) > 10:
            autocorr = np.corrcoef(y_data[:-1], y_data[1:])[0, 1]
            noise_correlation = float(autocorr) if np.isfinite(autocorr) else 0.0
        else:
            noise_correlation = 0.0

        return {
            'noise_level_robust': noise_robust,
            'noise_correlation': noise_correlation,
            'noise_variance': float(np.var(y_data))
        }

    def _estimate_peak_width_multiple_methods(self, x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, float]:
        """Multiple peak width estimation methods."""

        # FWHM method
        fwhm = self._calculate_fwhm(x_data, y_data)

        # Moment-based width
        moment_width = self._calculate_moment_width(x_data, y_data)

        # Peak width at different height fractions
        width_10 = self._calculate_width_at_fraction(x_data, y_data, 0.1)
        width_50 = self._calculate_width_at_fraction(x_data, y_data, 0.5)  # FWHM
        width_90 = self._calculate_width_at_fraction(x_data, y_data, 0.9)

        return {
            'peak_width_fwhm': fwhm,
            'peak_width_moment': moment_width,
            'peak_width_10percent': width_10,
            'peak_width_50percent': width_50,
            'peak_width_90percent': width_90,
            'peak_width_ratio_90_10': width_90 / (width_10 + 1e-8)
        }

    def _analyze_peak_asymmetry(self, x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, float]:
        """Comprehensive peak asymmetry analysis."""

        peak_idx = np.argmax(y_data)

        if 0 < peak_idx < len(y_data) - 1:
            # Traditional asymmetry
            left_area = float(np.sum(y_data[:peak_idx]))
            right_area = float(np.sum(y_data[peak_idx:]))
            total_area = left_area + right_area
            asymmetry = (right_area - left_area) / (total_area + 1e-8)

            # Tailing factor
            tailing = self._calculate_peak_tailing(x_data, y_data)

            # Fronting factor
            fronting = self._calculate_peak_fronting(x_data, y_data)
        else:
            asymmetry = 0.0
            tailing = 1.0
            fronting = 1.0

        return {
            'peak_asymmetry': float(np.clip(asymmetry, -1.0, 1.0)),
            'peak_tailing': tailing,
            'peak_fronting': fronting
        }

    def _analyze_baseline_characteristics(self, x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, float]:
        """Comprehensive baseline analysis."""

        # Baseline drift estimation
        baseline_points = y_data[:len(y_data)//10]  # First 10% of points
        if len(baseline_points) > 1:
            baseline_drift = float(np.std(baseline_points))
            baseline_slope = float(np.polyfit(x_data[:len(baseline_points)], baseline_points, 1)[0])
        else:
            baseline_drift = 0.0
            baseline_slope = 0.0

        # Baseline curvature (using more points)
        n_baseline = min(len(y_data)//5, 50)  # Use up to 20% or 50 points
        if n_baseline > 3:
            try:
                poly_coeffs = np.polyfit(x_data[:n_baseline], y_data[:n_baseline], 2)
                baseline_curvature = float(poly_coeffs[0])  # Quadratic coefficient
            except:
                baseline_curvature = 0.0
        else:
            baseline_curvature = 0.0

        return {
            'baseline_drift': baseline_drift,
            'baseline_slope': baseline_slope,
            'baseline_curvature': baseline_curvature
        }

    def _analyze_intensity_distribution(self, y_data: np.ndarray) -> Dict[str, float]:
        """Analyze intensity distribution characteristics."""

        # Statistical moments
        intensity_mean = float(np.mean(y_data))
        intensity_std = float(np.std(y_data))
        intensity_skew = self._calculate_skewness(y_data)
        intensity_kurtosis = float(stats.kurtosis(y_data)) if len(y_data) > 3 else 0.0

        return {
            'intensity_mean': intensity_mean,
            'intensity_std': intensity_std,
            'intensity_skewness': intensity_skew,
            'intensity_kurtosis': intensity_kurtosis
        }

    # Additional helper methods (implementations)

    def _calculate_fwhm(self, x_data: np.ndarray, y_data: np.ndarray) -> float:
        """Calculate Full Width at Half Maximum."""
        peak_height = np.max(y_data)
        baseline = np.min(y_data)
        half_max = baseline + (peak_height - baseline) / 2

        above_half = y_data >= half_max
        if np.sum(above_half) >= 2:
            indices = np.where(above_half)[0]
            width_points = indices[-1] - indices[0] + 1
            if len(x_data) > 1:
                return width_points * abs(x_data[1] - x_data[0])

        return 0.02  # Default fallback

    def _calculate_moment_width(self, x_data: np.ndarray, y_data: np.ndarray) -> float:
        """Calculate width using second moment."""
        if len(x_data) != len(y_data) or len(x_data) < 3:
            return 0.02

        # Center of mass
        total_intensity = np.sum(y_data)
        if total_intensity == 0:
            return 0.02

        center_of_mass = np.sum(x_data * y_data) / total_intensity

        # Second moment (variance)
        variance = np.sum(y_data * (x_data - center_of_mass)**2) / total_intensity

        return float(2 * np.sqrt(variance))  # Width as 2*sigma

    def _calculate_width_at_fraction(self, x_data: np.ndarray, y_data: np.ndarray, fraction: float) -> float:
        """Calculate peak width at given height fraction."""
        peak_height = np.max(y_data)
        baseline = np.min(y_data)
        threshold = baseline + (peak_height - baseline) * fraction

        above_threshold = y_data >= threshold
        if np.sum(above_threshold) >= 2:
            indices = np.where(above_threshold)[0]
            width_points = indices[-1] - indices[0] + 1
            if len(x_data) > 1:
                return width_points * abs(x_data[1] - x_data[0])

        return 0.02  # Default fallback

    def _calculate_peak_tailing(self, x_data: np.ndarray, y_data: np.ndarray) -> float:
        """Calculate peak tailing factor."""
        peak_idx = np.argmax(y_data)
        if peak_idx == 0 or peak_idx == len(y_data) - 1:
            return 1.0

        # Find 10% height points
        peak_height = y_data[peak_idx]
        baseline = np.min(y_data)
        ten_percent = baseline + (peak_height - baseline) * 0.1

        # Find leading and tailing edges
        leading_idx = np.where(y_data[:peak_idx] >= ten_percent)[0]
        tailing_idx = np.where(y_data[peak_idx:] >= ten_percent)[0]

        if len(leading_idx) > 0 and len(tailing_idx) > 0:
            leading_width = peak_idx - leading_idx[0]
            tailing_width = tailing_idx[-1]
            if leading_width > 0:
                tailing_ratio = tailing_width / leading_width
                if tailing_ratio <= 0:
                    return 1.0  # Treat degenerate shapes as symmetric
                return float(tailing_ratio)

        return 1.0

    def _calculate_peak_fronting(self, x_data: np.ndarray, y_data: np.ndarray) -> float:
        """Calculate peak fronting factor."""
        tailing = self._calculate_peak_tailing(x_data, y_data)

        if tailing <= 0:
            return 1.0  # Fall back to neutral fronting when tailing is invalid

        return 1.0 / tailing

    def _classify_chemical_region_enhanced(self, shift: float, nucleus_type: str) -> str:
        """Enhanced chemical shift region classification."""

        if '1H' in nucleus_type:
            if shift < 1.0:
                return 'aliphatic_methyl'
            elif shift < 2.0:
                return 'aliphatic_upfield'
            elif shift < 4.0:
                return 'aliphatic_downfield'
            elif shift < 5.5:
                return 'alpha_protons'
            elif shift < 6.0:
                return 'oxygen_bearing'
            elif shift < 7.0:
                return 'olefinic'
            elif shift < 8.0:
                return 'aromatic'
            elif shift < 10.0:
                return 'aromatic_downfield'
            else:
                return 'highly_deshielded'

        elif '15N' in nucleus_type:
            if shift < 100:
                return 'aliphatic_amines'
            elif shift < 110:
                return 'amines_upfield'
            elif shift < 125:
                return 'amides_backbone'
            elif shift < 140:
                return 'amides_sidechain'
            elif shift < 160:
                return 'aromatics_heterocycles'
            elif shift < 180:
                return 'imines_imidazoles'
            else:
                return 'highly_deshielded'

        elif '13C' in nucleus_type:
            if shift < 20:
                return 'methyl_carbons'
            elif shift < 50:
                return 'aliphatic'
            elif shift < 80:
                return 'oxygenated_carbons'
            elif shift < 100:
                return 'anomeric_carbons'
            elif shift < 140:
                return 'aromatic_CH'
            elif shift < 160:
                return 'aromatic_quaternary'
            elif shift < 180:
                return 'carbonyl_esters'
            elif shift < 210:
                return 'carbonyl_ketones'
            else:
                return 'highly_deshielded'

        return 'unknown_region'

    def _predict_multiplicity(self, shift: float, nucleus_type: str) -> str:
        """Predict expected multiplicity based on chemical shift."""

        # Simplified multiplicity prediction
        if '1H' in nucleus_type:
            if 0.5 <= shift <= 2.0:
                return 'methyl_or_methylene'
            elif 2.0 < shift <= 4.5:
                return 'methine_or_heteroatom'
            elif 6.5 <= shift <= 8.5:
                return 'aromatic'
            else:
                return 'variable'
        elif '15N' in nucleus_type:
            return 'singlet'  # Most 15N are singlets in HSQC
        elif '13C' in nucleus_type:
            return 'singlet'  # Most 13C are singlets in HSQC

        return 'unknown'

    def _assess_coupling_environment(self, shift: float, nucleus_type: str) -> str:
        """Assess coupling environment."""

        # Simplified coupling assessment
        if '1H' in nucleus_type:
            if 0.5 <= shift <= 4.5:
                return 'aliphatic_coupled'
            elif 6.5 <= shift <= 8.5:
                return 'aromatic_coupled'
            else:
                return 'isolated_or_exchangeable'

        return 'decoupled'

    def _estimate_relaxation_properties(self, shift: float, nucleus_type: str) -> Dict[str, float]:
        """Estimate relaxation properties."""

        # Simplified relaxation estimates based on nucleus and environment
        if '1H' in nucleus_type:
            return {"T1_estimate": 2.0, "T2_estimate": 0.5}  # Typical for 1H
        elif '15N' in nucleus_type:
            return {"T1_estimate": 5.0, "T2_estimate": 0.1}  # Typical for 15N
        elif '13C' in nucleus_type:
            return {"T1_estimate": 10.0, "T2_estimate": 0.2}  # Typical for 13C

        return {"T1_estimate": 1.0, "T2_estimate": 0.1}

    def _estimate_chemical_shift_anisotropy(self, shift: float, nucleus_type: str) -> float:
        """Estimate chemical shift anisotropy."""

        # Simplified CSA estimates
        if '15N' in nucleus_type:
            return 160.0  # Typical CSA for 15N
        elif '13C' in nucleus_type:
            if shift > 100:  # Aromatic
                return 180.0
            else:  # Aliphatic
                return 25.0

        return 0.0

    def _estimate_j_coupling(self, shift: float, nucleus_type: str) -> Dict[str, float]:
        """Estimate J-coupling constants."""

        # Simplified J-coupling estimates
        if '1H' in nucleus_type:
            return {"J_HH_typical": 7.0, "J_HH_range": [6.0, 8.0]}
        elif '15N' in nucleus_type:
            return {"J_NH_typical": 94.0, "J_NH_range": [90.0, 100.0]}
        elif '13C' in nucleus_type:
            return {"J_CH_typical": 140.0, "J_CH_range": [125.0, 160.0]}

        return {}

    # Additional helper methods continue...
    # (Many more methods would be implemented for complete functionality)

    def _analyze_peak_overlap_enhanced(self, context: Dict, x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, float]:
        """Enhanced peak overlap analysis."""

        overlapping_peaks = context.get('overlapping_peaks', [])
        if not overlapping_peaks:
            return {'severity': 0.0, 'type': 'isolated'}

        # Analyze overlap characteristics
        peak_separations = []
        for other_peak in overlapping_peaks:
            if 'center' in other_peak:
                separation = abs(other_peak['center'] - np.mean(x_data))
                peak_separations.append(separation)

        if peak_separations:
            min_separation = min(peak_separations)
            avg_separation = np.mean(peak_separations)

            # Estimate typical peak width
            peak_width = self._calculate_fwhm(x_data, y_data)

            # Overlap severity based on separation vs width
            if min_separation < peak_width / 2:
                severity = 0.9  # Severe overlap
                overlap_type = 'severe'
            elif min_separation < peak_width:
                severity = 0.6  # Moderate overlap
                overlap_type = 'moderate'
            elif min_separation < 2 * peak_width:
                severity = 0.3  # Mild overlap
                overlap_type = 'mild'
            else:
                severity = 0.1  # Minimal overlap
                overlap_type = 'minimal'
        else:
            severity = 0.0
            overlap_type = 'isolated'

        return {
            'severity': severity,
            'type': overlap_type,
            'nearby_peak_count': len(overlapping_peaks),
            'min_separation': min(peak_separations) if peak_separations else float('inf'),
            'avg_separation': np.mean(peak_separations) if peak_separations else float('inf')
        }

    def _identify_interference_sources(self, context: Dict) -> List[str]:
        """Identify sources of interference."""

        sources = []

        # Check for various interference types
        if context.get('baseline_issues', False):
            sources.append('baseline_drift')

        if context.get('noise_level', 0) > 0.1:
            sources.append('high_noise')

        overlapping_peaks = context.get('overlapping_peaks', [])
        if len(overlapping_peaks) > 2:
            sources.append('multiple_overlaps')
        elif len(overlapping_peaks) > 0:
            sources.append('peak_overlap')

        if context.get('peak_complexity', 'simple') in ['complex', 'very_complex']:
            sources.append('peak_complexity')

        if context.get('solvent_peaks', False):
            sources.append('solvent_interference')

        return sources

    def _calculate_peak_isolation_score(self, context: Dict, x_data: np.ndarray) -> float:
        """Calculate how well-isolated the peak is."""

        overlapping_peaks = context.get('overlapping_peaks', [])

        if not overlapping_peaks:
            return 1.0  # Perfectly isolated

        # Calculate isolation based on peak separations
        peak_center = np.mean(x_data)
        peak_width = abs(x_data[-1] - x_data[0])

        min_distance = float('inf')
        for other_peak in overlapping_peaks:
            if 'center' in other_peak:
                distance = abs(other_peak['center'] - peak_center)
                min_distance = min(min_distance, distance)

        if min_distance == float('inf'):
            return 1.0

        # Isolation score based on distance relative to peak width
        isolation_ratio = min_distance / (peak_width + 1e-8)
        isolation_score = 1.0 / (1.0 + np.exp(-2 * (isolation_ratio - 1)))  # Sigmoid function

        return float(np.clip(isolation_score, 0.0, 1.0))

    def _assess_baseline_stability_enhanced(self, y_data: np.ndarray) -> float:
        """Enhanced baseline stability assessment."""

        # Use first and last 10% of data points for baseline assessment
        n_baseline = max(len(y_data) // 10, 3)

        baseline_start = y_data[:n_baseline]
        baseline_end = y_data[-n_baseline:]

        # Calculate stability metrics
        start_std = np.std(baseline_start)
        end_std = np.std(baseline_end)
        start_mean = np.mean(baseline_start)
        end_mean = np.mean(baseline_end)

        # Stability based on consistency between start and end
        if start_mean != 0:
            baseline_drift = abs(end_mean - start_mean) / abs(start_mean)
        else:
            baseline_drift = 0.0

        # Stability based on local variations
        avg_std = (start_std + end_std) / 2
        if start_mean != 0:
            relative_noise = avg_std / abs(start_mean)
        else:
            relative_noise = 0.0

        # Combine metrics into stability score (0-1, higher is more stable)
        stability = 1.0 / (1.0 + baseline_drift + relative_noise)

        return float(np.clip(stability, 0.0, 1.0))

    def _assess_data_quality_comprehensive(self, x_data: np.ndarray, y_data: np.ndarray, context: Dict) -> float:
        """Comprehensive data quality assessment."""

        quality_factors = []

        # Resolution factor
        if len(x_data) > 1:
            resolution = abs(x_data[1] - x_data[0])
            if resolution < 0.01:
                quality_factors.append(1.0)  # High resolution
            elif resolution < 0.05:
                quality_factors.append(0.8)  # Good resolution
            else:
                quality_factors.append(0.6)  # Lower resolution
        else:
            quality_factors.append(0.5)

        # Signal-to-noise factor
        snr = self._calculate_snr_robust(y_data)
        if snr > 10:
            quality_factors.append(1.0)
        elif snr > 5:
            quality_factors.append(0.8)
        elif snr > 2:
            quality_factors.append(0.6)
        else:
            quality_factors.append(0.3)

        # Dynamic range factor
        dynamic_range = np.max(y_data) - np.min(y_data)
        if dynamic_range > 1000:
            quality_factors.append(1.0)
        elif dynamic_range > 100:
            quality_factors.append(0.8)
        else:
            quality_factors.append(0.6)

        # Data completeness (no missing values, reasonable range)
        if np.all(np.isfinite(y_data)) and np.all(np.isfinite(x_data)):
            quality_factors.append(1.0)
        else:
            quality_factors.append(0.3)

        # Context-based quality adjustments
        if context.get('preprocessing_applied', False):
            quality_factors.append(0.9)  # Slight penalty for heavy preprocessing

        return float(np.mean(quality_factors))

    def _calculate_snr_robust(self, y_data: np.ndarray) -> float:
        """Calculate robust signal-to-noise ratio."""

        signal = np.max(y_data) - np.median(y_data)

        # Use MAD for robust noise estimation
        noise = np.median(np.abs(y_data - np.median(y_data))) * 1.4826

        if noise > 0:
            return signal / noise
        else:
            return float('inf')

    def _classify_peak_complexity_enhanced(self, context: Dict, nearby_peaks_count: int) -> str:
        """Enhanced peak complexity classification."""

        # Base complexity from nearby peaks
        if nearby_peaks_count == 0:
            base_complexity = 'simple'
        elif nearby_peaks_count <= 2:
            base_complexity = 'moderate'
        elif nearby_peaks_count <= 5:
            base_complexity = 'complex'
        else:
            base_complexity = 'very_complex'

        # Adjust based on other factors
        if context.get('baseline_issues', False):
            if base_complexity == 'simple':
                base_complexity = 'moderate'
            elif base_complexity == 'moderate':
                base_complexity = 'complex'

        if context.get('peak_asymmetry', 0) > 0.5:
            if base_complexity == 'simple':
                base_complexity = 'moderate'

        return base_complexity

    def _estimate_noise_level_advanced(self, y_data: np.ndarray, context: Dict) -> float:
        """Advanced noise level estimation."""

        # Multiple noise estimation methods
        methods = []

        # Method 1: MAD-based
        mad_noise = np.median(np.abs(y_data - np.median(y_data))) * 1.4826
        methods.append(mad_noise)

        # Method 2: Percentile-based
        low_percentile_std = np.std(y_data[y_data <= np.percentile(y_data, 25)])
        methods.append(low_percentile_std)

        # Method 3: Corner-based (if data is 2D-like)
        if len(y_data) > 20:
            corner_size = len(y_data) // 10
            corner_noise = np.std(y_data[:corner_size])
            methods.append(corner_noise)

        # Method 4: Context-provided
        context_noise = context.get('noise_level', 0)
        if context_noise > 0:
            methods.append(context_noise)

        # Return robust average
        valid_methods = [m for m in methods if m > 0 and np.isfinite(m)]
        if valid_methods:
            return float(np.median(valid_methods))
        else:
            return 0.01  # Default fallback

    def _assess_optimization_difficulty(self, opt_info: Dict) -> float:
        """Assess optimization difficulty score."""

        difficulty_factors = []

        # Iteration count factor
        iterations = opt_info.get('iterations', 0)
        if iterations > 100:
            difficulty_factors.append(0.9)  # High difficulty
        elif iterations > 50:
            difficulty_factors.append(0.6)  # Moderate difficulty
        elif iterations > 10:
            difficulty_factors.append(0.3)  # Low difficulty
        else:
            difficulty_factors.append(0.1)  # Very easy

        # Convergence factor
        if not opt_info.get('converged', True):
            difficulty_factors.append(0.8)  # Failed convergence indicates difficulty

        # Gradient norm factor
        gradient_norm = opt_info.get('gradient_norm', 0)
        if gradient_norm > 1.0:
            difficulty_factors.append(0.7)
        elif gradient_norm > 0.1:
            difficulty_factors.append(0.4)
        else:
            difficulty_factors.append(0.2)

        # Hessian condition number
        condition = opt_info.get('hessian_condition', 1.0)
        if condition > 1000:
            difficulty_factors.append(0.8)  # Ill-conditioned
        elif condition > 100:
            difficulty_factors.append(0.5)
        else:
            difficulty_factors.append(0.2)

        return float(np.mean(difficulty_factors)) if difficulty_factors else 0.5

    def _validate_parameter_physics(self, param: str, value: float, nucleus_type: str) -> bool:
        """Validate parameter against physical constraints."""

        if not np.isfinite(value):
            return False

        # Get constraints for nucleus type
        constraints = self.physics_constraints.get(nucleus_type, {})

        if param == 'amplitude' and value <= 0:
            return False
        elif param in ['sigma', 'gamma']:
            width_range = constraints.get('width_range', (0, float('inf')))
            return width_range[0] <= value <= width_range[1]
        elif param == 'center':
            shift_range = constraints.get('shift_range', (-float('inf'), float('inf')))
            return shift_range[0] <= value <= shift_range[1]
        elif param == 'baseline':
            return True  # Baseline can be any finite value

        return True

    def _assess_peak_shape_realism(self, sigma: float, gamma: float, nucleus_type: str) -> float:
        """Assess how realistic the peak shape parameters are."""

        if sigma <= 0 or gamma < 0:
            return 0.0

        # Voigt profile mixing parameter
        mixing_ratio = gamma / (sigma + gamma + 1e-10)

        # Realistic mixing ratios for different nuclei
        if '1H' in nucleus_type:
            # 1H typically more Gaussian
            if 0.1 <= mixing_ratio <= 0.7:
                return 1.0
            elif mixing_ratio < 0.1:
                return 0.8  # Very Gaussian
            else:
                return 0.6  # Very Lorentzian
        elif '15N' in nucleus_type or '13C' in nucleus_type:
            # Heteronuclei can be more Lorentzian
            if 0.2 <= mixing_ratio <= 0.8:
                return 1.0
            else:
                return 0.7

        return 0.8  # Default reasonable score

    def _check_linewidth_consistency(self, sigma: float, gamma: float, nucleus_type: str) -> bool:
        """Check if linewidths are consistent with nucleus type."""

        total_width = sigma + gamma

        if '1H' in nucleus_type:
            return 0.001 <= total_width <= 0.1
        elif '15N' in nucleus_type:
            return 0.1 <= total_width <= 10
        elif '13C' in nucleus_type:
            return 0.1 <= total_width <= 5

        return True  # Unknown nucleus, assume valid

    def _validate_chemical_shift(self, center: float, nucleus_type: str) -> bool:
        """Validate chemical shift against typical ranges."""

        if '1H' in nucleus_type:
            return -2 <= center <= 15
        elif '15N' in nucleus_type:
            return 0 <= center <= 250
        elif '13C' in nucleus_type:
            return -10 <= center <= 250

        return True  # Unknown nucleus, assume valid

    def _check_spectroscopic_constraints(self, fit_result: Dict, nucleus_type: str) -> Dict[str, bool]:
        """Check various spectroscopic constraints."""

        constraints = {}

        # Amplitude positivity
        constraints['amplitude_positive'] = fit_result.get('amplitude', 0) > 0

        # Width positivity
        constraints['sigma_positive'] = fit_result.get('sigma', 0) > 0
        constraints['gamma_nonnegative'] = fit_result.get('gamma', 0) >= 0

        # Chemical shift reasonableness
        center = fit_result.get('center', 0)
        constraints['shift_reasonable'] = self._validate_chemical_shift(center, nucleus_type)

        # Peak shape reasonableness
        sigma = fit_result.get('sigma', 0.01)
        gamma = fit_result.get('gamma', 0.01)
        constraints['shape_reasonable'] = self._assess_peak_shape_realism(sigma, gamma, nucleus_type) > 0.5

        return constraints

    def _analyze_residuals_comprehensive(self, residuals: List) -> Dict[str, float]:
        """Comprehensive residual analysis."""

        if not residuals or len(residuals) < 3:
            return {}

        residuals = np.array(residuals)

        analysis = {}

        # Basic statistics
        analysis['residual_mean'] = float(np.mean(residuals))
        analysis['residual_std'] = float(np.std(residuals))
        analysis['residual_max'] = float(np.max(np.abs(residuals)))

        # Distribution analysis
        analysis['residual_skewness'] = self._calculate_skewness(residuals)
        analysis['residual_kurtosis'] = float(stats.kurtosis(residuals)) if len(residuals) > 3 else 0.0

        # Autocorrelation analysis
        analysis['residual_autocorrelation'] = self._calculate_residual_autocorrelation(residuals)

        # Runs test for randomness
        analysis['runs_test_pvalue'] = self._runs_test(residuals)

        return analysis

    def _calculate_residual_autocorrelation(self, residuals: List) -> float:
        """Calculate autocorrelation of residuals."""

        if len(residuals) < 2:
            return 0.0

        residuals = np.array(residuals)

        # Lag-1 autocorrelation
        if len(residuals) > 1:
            autocorr = np.corrcoef(residuals[:-1], residuals[1:])[0, 1]
            return float(autocorr) if np.isfinite(autocorr) else 0.0

        return 0.0

    def _runs_test(self, residuals: List) -> float:
        """Perform runs test for residual randomness."""

        if len(residuals) < 10:
            return 1.0  # Can't perform test

        residuals = np.array(residuals)
        median = np.median(residuals)

        # Convert to binary sequence
        binary = (residuals > median).astype(int)

        # Count runs
        runs = 1
        for i in range(1, len(binary)):
            if binary[i] != binary[i-1]:
                runs += 1

        # Approximate p-value (simplified)
        n = len(binary)
        expected_runs = (2 * n - 1) / 3

        if runs > expected_runs:
            p_value = 2 * (1 - stats.norm.cdf(abs(runs - expected_runs) / np.sqrt(expected_runs)))
        else:
            p_value = 2 * stats.norm.cdf(abs(runs - expected_runs) / np.sqrt(expected_runs))

        return float(p_value)

    def _detect_outlier_fraction(self, y_data: np.ndarray, fit_result: Dict) -> float:
        """Detect fraction of outlier points."""

        # Use residuals if available, otherwise estimate
        residuals = fit_result.get('residuals', [])

        if not residuals:
            return 0.0  # Can't assess without residuals

        residuals = np.array(residuals)

        # Use MAD for robust outlier detection
        median_residual = np.median(residuals)
        mad = np.median(np.abs(residuals - median_residual)) * 1.4826

        # Points beyond 3 MADs are outliers
        outliers = np.abs(residuals - median_residual) > 3 * mad

        return float(np.sum(outliers) / len(residuals))

    def _perform_goodness_of_fit_tests(self, y_data: np.ndarray, fit_result: Dict, residuals: List) -> Dict[str, float]:
        """Perform various goodness-of-fit tests."""

        tests = {}

        if not residuals or len(residuals) < 10:
            return tests

        residuals = np.array(residuals)

        # Shapiro-Wilk test for normality of residuals
        try:
            if len(residuals) >= 3:
                stat, p_value = stats.shapiro(residuals)
                tests['shapiro_wilk_pvalue'] = float(p_value)
        except:
            tests['shapiro_wilk_pvalue'] = 1.0

        # Kolmogorov-Smirnov test
        try:
            stat, p_value = stats.kstest(residuals, 'norm')
            tests['ks_test_pvalue'] = float(p_value)
        except:
            tests['ks_test_pvalue'] = 1.0

        # Jarque-Bera test
        try:
            if len(residuals) >= 8:
                stat, p_value = stats.jarque_bera(residuals)
                tests['jarque_bera_pvalue'] = float(p_value)
        except:
            tests['jarque_bera_pvalue'] = 1.0

        return tests

    def _assess_data_quality_flags(self, x_data: np.ndarray, y_data: np.ndarray) -> List[str]:
        """Assess data quality and return flags for issues."""

        flags = []

        # Check for missing or invalid data
        if np.any(~np.isfinite(x_data)):
            flags.append('invalid_x_data')
        if np.any(~np.isfinite(y_data)):
            flags.append('invalid_y_data')

        # Check for insufficient data
        if len(x_data) < 10:
            flags.append('insufficient_data_points')

        # Check for zero or negative intensities
        if np.any(y_data <= 0):
            flags.append('zero_or_negative_intensities')

        # Check for very low dynamic range
        dynamic_range = np.max(y_data) - np.min(y_data)
        if dynamic_range < 10:
            flags.append('low_dynamic_range')

        # Check for potential baseline issues
        baseline_variation = np.std(y_data[:len(y_data)//10])  # First 10%
        signal_level = np.max(y_data)
        if baseline_variation > signal_level * 0.1:
            flags.append('baseline_instability')

        # Check for potential saturation
        if np.sum(y_data == np.max(y_data)) > len(y_data) * 0.05:  # > 5% at max
            flags.append('potential_saturation')

        return flags

    def _calculate_skewness(self, data: np.ndarray) -> float:
        """Calculate skewness of distribution."""
        if len(data) < 3:
            return 0.0

        mean_val = np.mean(data)
        std_val = np.std(data)

        if std_val == 0:
            return 0.0

        skewness = np.mean(((data - mean_val) / std_val) ** 3)
        return float(np.clip(skewness, -3.0, 3.0))  # Clip extreme values

    def _downsample_data(self, data: np.ndarray, max_points: int = 200) -> np.ndarray:
        """Downsample data for efficient storage while preserving shape."""
        if len(data) <= max_points:
            return data

        # Use uniform sampling to preserve overall shape
        indices = np.linspace(0, len(data) - 1, max_points, dtype=int)
        return data[indices]

    def _save_batch(self) -> None:
        """Save collected training data to storage."""
        if not self.session_data:
            return

        try:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = self.storage_path / f"enhanced_training_batch_{timestamp}.pkl"

            # Ensure directory exists
            filename.parent.mkdir(exist_ok=True)

            # Save with compression
            with open(filename, 'wb') as f:
                pickle.dump(self.session_data, f, protocol=pickle.HIGHEST_PROTOCOL)

            # Enhanced metadata
            metadata = {
                'timestamp': datetime.now().isoformat(),
                'collection_version': '2.0',
                'samples_count': len(self.session_data),
                'collection_stats': self.collection_stats.copy(),
                'feature_schema': self.feature_schema,
                'window_settings': self.window_settings,
                'enhancements': [
                    'optimization_trajectory', 'multi_method_comparison',
                    'physics_validation', 'enhanced_quality_metrics',
                    'comprehensive_spectral_analysis', 'chemical_environment_prediction'
                ],
                'quality_threshold': self.min_r_squared
            }

            metadata_file = filename.with_suffix('.json')
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)

            self.logger.info(f"Saved enhanced ML training batch: {len(self.session_data)} samples to {filename}")

            # Update statistics
            self.collection_stats['batch_saves'] += 1

            # Clear session data
            self.session_data.clear()

        except Exception as e:
            self.logger.error(f"Failed to save enhanced training batch: {e}")
            # Don't clear data on save failure - try again later

    # Backward compatibility methods for existing tests
    def _extract_spectral_features(self, x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, float]:
        """
        Backward compatibility method for v1.0 tests.

        Maps new enhanced features to old expected format.
        """
        enhanced_features = self._extract_enhanced_spectral_features(x_data, y_data)

        # Map new features to old expected names
        return {
            'peak_height': enhanced_features.get('peak_height', 0.0),
            'baseline_level': enhanced_features.get('baseline_level', 0.0),
            'snr': enhanced_features.get('snr', 0.0),
            'peak_width_estimate': enhanced_features.get('peak_width_fwhm', 0.0)  # Use FWHM as primary width
        }

    def _extract_chemical_features(self, x_data: np.ndarray, nucleus_type: str) -> Dict[str, Any]:
        """
        Backward compatibility method for v1.0 tests.

        Maps new enhanced chemical features to old expected format.
        """
        enhanced_features = self._extract_enhanced_chemical_features(x_data, nucleus_type)

        # Map new features to old expected names
        return {
            'nucleus_type': enhanced_features.get('nucleus_type', nucleus_type),
            'chemical_shift_center': enhanced_features.get('chemical_shift', 0.0),
            'chemical_region': enhanced_features.get('chemical_region', 'unknown')
        }

    def finalize_session(self) -> Dict[str, Any]:
        """Finalize collection session and save remaining data."""

        # Save any remaining data
        if self.session_data:
            self._save_batch()

        # Return session statistics
        return {
            'session_stats': self.collection_stats.copy(),
            'storage_path': str(self.storage_path),
            'total_files': len(list(self.storage_path.glob("*training_batch_*.pkl"))),
            'collection_version': self.collection_version,
            'enhanced_features': self.enhanced_features_enabled
        }

    def get_collection_stats(self) -> Dict[str, Any]:
        """Get current collection statistics."""
        stats = self.collection_stats.copy()
        stats['pending_samples'] = len(self.session_data)
        stats['storage_path'] = str(self.storage_path)
        stats['collection_version'] = self.collection_version

        if stats['total_attempts'] > 0:
            stats['collection_rate'] = stats['successfully_collected'] / stats['total_attempts']
            stats['quality_filter_rate'] = stats['quality_filtered'] / stats['total_attempts']
            stats['error_rate'] = stats['errors'] / stats['total_attempts']
            stats['convergence_failure_rate'] = stats['convergence_failures'] / stats['total_attempts']

        return stats

    def __del__(self):
        """Cleanup on object destruction."""
        try:
            if hasattr(self, 'session_data') and self.session_data:
                self._save_batch()
        except:
            pass  # Silent cleanup
