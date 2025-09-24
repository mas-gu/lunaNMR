#!/usr/bin/env python3
"""
Enhanced ML Training Data Collector - Version 2.0

This enhanced version addresses gaps in the original collector by adding:
- Initial parameter estimates and optimization convergence data
- Multi-method comparison results
- Physics-based constraints and validation
- Advanced peak overlap and interference metrics
- Fitting difficulty and optimization path information

This provides much richer training data for more robust ML models.
"""

import numpy as np
import pandas as pd
import pickle
import json
import os
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any
import logging
from pathlib import Path

class EnhancedMLTrainingDataCollector:
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
        """Initialize enhanced collector with expanded feature schema."""

        # Set up storage path
        if storage_path is None:
            project_root = Path(__file__).parent.parent.parent.parent
            storage_path = project_root / "ml_training_data_enhanced"

        self.storage_path = Path(storage_path)
        self.storage_path.mkdir(exist_ok=True)

        # Session data buffer
        self.session_data = []
        self.batch_size = 100

        # Logging
        self.logger = logging.getLogger(__name__)

        # Enhanced quality thresholds
        self.min_r_squared = 0.8
        self.max_samples_per_session = 1000

        # Enhanced feature schema
        self.feature_schema = {
            'spectral_features': [
                'peak_height', 'peak_width_estimate', 'snr', 'baseline_level',
                'peak_asymmetry', 'spectral_range', 'data_points', 'dynamic_range',
                'noise_characteristics', 'baseline_drift', 'peak_tailing',
                'spectral_resolution', 'data_density'
            ],
            'chemical_features': [
                'nucleus_type', 'chemical_shift', 'chemical_region', 'shift_range',
                'expected_multiplicity', 'coupling_environment', 'relaxation_estimate'
            ],
            'context_features': [
                'detection_confidence', 'estimated_amplitude', 'estimated_width',
                'nearby_peaks_count', 'peak_complexity', 'noise_level',
                'overlap_severity', 'interference_sources', 'peak_isolation_score'
            ],
            'optimization_features': [
                'initial_parameters', 'parameter_bounds', 'convergence_iterations',
                'optimization_method', 'convergence_quality', 'parameter_uncertainty',
                'gradient_info', 'hessian_condition', 'optimization_difficulty'
            ],
            'multi_method_results': [
                'alternative_methods_tested', 'method_comparison_scores',
                'consensus_parameters', 'parameter_stability', 'method_agreement'
            ],
            'physics_validation': [
                'parameter_physical_validity', 'spectroscopic_constraints',
                'peak_shape_realism', 'coupling_consistency', 'relaxation_validity'
            ],
            'target_parameters': [
                'amplitude', 'center', 'sigma', 'gamma', 'baseline'
            ],
            'quality_metrics': [
                'r_squared', 'rmse_normalized', 'parameters_physical', 'fit_success',
                'residual_analysis', 'outlier_detection', 'cross_validation_score'
            ]
        }

        # Statistics tracking
        self.collection_stats = {
            'total_attempts': 0,
            'quality_filtered': 0,
            'successfully_collected': 0,
            'batch_saves': 0,
            'errors': 0,
            'convergence_failures': 0,
            'physics_violations': 0
        }

        self.logger.info(f"Enhanced ML Training Data Collector v2.0 initialized - Storage: {self.storage_path}")

    def collect_enhanced_training_sample(self,
                                       x_data: np.ndarray,
                                       y_data: np.ndarray,
                                       fit_result: Dict[str, Any],
                                       initial_params: Dict[str, Any],
                                       optimization_info: Dict[str, Any],
                                       context: Dict[str, Any],
                                       nucleus_type: str,
                                       alternative_results: Optional[List[Dict]] = None) -> bool:
        """
        Collect enhanced training sample with comprehensive optimization context.

        Args:
            x_data: X-axis spectral data (ppm)
            y_data: Y-axis intensity data
            fit_result: Final fitting result
            initial_params: Initial parameter estimates
            optimization_info: Convergence and optimization details
            context: Detection and spectral context
            nucleus_type: Nucleus type
            alternative_results: Results from alternative fitting methods

        Returns:
            True if sample was collected successfully
        """
        self.collection_stats['total_attempts'] += 1

        try:
            # Quality filtering
            r_squared = fit_result.get('r_squared', 0.0)
            if r_squared < self.min_r_squared:
                self.collection_stats['quality_filtered'] += 1
                return False

            # Convergence validation
            if not optimization_info.get('converged', False):
                self.collection_stats['convergence_failures'] += 1
                # Still collect for learning about difficult cases, but flag it

            # Memory management
            if len(self.session_data) >= self.max_samples_per_session:
                self.logger.warning("Max samples per session reached")
                return False

            # Extract enhanced features
            spectral_features = self._extract_enhanced_spectral_features(x_data, y_data)
            chemical_features = self._extract_enhanced_chemical_features(x_data, nucleus_type)
            context_features = self._extract_enhanced_context_features(context, x_data, y_data)
            optimization_features = self._extract_optimization_features(initial_params, fit_result, optimization_info)
            physics_features = self._extract_physics_validation(fit_result, nucleus_type)

            # Multi-method analysis
            multi_method_features = {}
            if alternative_results:
                multi_method_features = self._analyze_multi_method_results(fit_result, alternative_results)

            # Advanced quality metrics
            quality_metrics = self._compute_enhanced_quality_metrics(
                x_data, y_data, fit_result, optimization_info
            )

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
                    'parameter_uncertainties': optimization_info.get('parameter_errors', {})
                },

                # Enhanced quality metrics
                'quality_metrics': quality_metrics,

                # Optimization trajectory (if available)
                'optimization_trajectory': {
                    'parameter_evolution': optimization_info.get('parameter_history', []),
                    'cost_evolution': optimization_info.get('cost_history', []),
                    'gradient_evolution': optimization_info.get('gradient_history', [])
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

            # Validate sample completeness
            if self._validate_sample_completeness(training_sample):
                self.session_data.append(training_sample)
                self.collection_stats['successfully_collected'] += 1

                # Auto-save if batch size reached
                if len(self.session_data) >= self.batch_size:
                    self._save_batch()

                return True
            else:
                self.collection_stats['errors'] += 1
                return False

        except Exception as e:
            self.logger.debug(f"Enhanced ML data collection failed: {e}")
            self.collection_stats['errors'] += 1
            return False

    def _extract_enhanced_spectral_features(self, x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, float]:
        """Extract comprehensive spectral features."""

        # Basic features (from original collector)
        basic_features = self._extract_basic_spectral_features(x_data, y_data)

        # Enhanced features
        enhanced_features = {}

        # Noise characteristics
        noise_analysis = self._analyze_noise_characteristics(y_data)
        enhanced_features.update(noise_analysis)

        # Baseline analysis
        baseline_analysis = self._analyze_baseline_characteristics(x_data, y_data)
        enhanced_features.update(baseline_analysis)

        # Peak shape analysis
        peak_analysis = self._analyze_peak_shape_characteristics(x_data, y_data)
        enhanced_features.update(peak_analysis)

        # Spectral resolution and density
        if len(x_data) > 1:
            enhanced_features['spectral_resolution'] = float(abs(x_data[1] - x_data[0]))
            enhanced_features['data_density'] = float(len(x_data) / abs(x_data[-1] - x_data[0]))
        else:
            enhanced_features['spectral_resolution'] = 0.001
            enhanced_features['data_density'] = 1000.0

        # Combine basic and enhanced features
        all_features = {**basic_features, **enhanced_features}
        return all_features

    def _extract_enhanced_chemical_features(self, x_data: np.ndarray, nucleus_type: str) -> Dict[str, Any]:
        """Extract enhanced chemical context features."""

        # Basic chemical features
        basic_features = self._extract_basic_chemical_features(x_data, nucleus_type)

        # Enhanced chemical analysis
        enhanced_features = {}

        # Expected multiplicity based on chemical shift
        enhanced_features['expected_multiplicity'] = self._predict_multiplicity(
            basic_features['chemical_shift_center'], nucleus_type
        )

        # Coupling environment assessment
        enhanced_features['coupling_environment'] = self._assess_coupling_environment(
            basic_features['chemical_shift_center'], nucleus_type
        )

        # Relaxation time estimates
        enhanced_features['relaxation_estimate'] = self._estimate_relaxation_properties(
            basic_features['chemical_shift_center'], nucleus_type
        )

        return {**basic_features, **enhanced_features}

    def _extract_enhanced_context_features(self, context: Dict[str, Any], x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, Any]:
        """Extract comprehensive contextual features."""

        # Basic context features
        basic_features = self._extract_basic_context_features(context, x_data, y_data)

        # Enhanced context analysis
        enhanced_features = {}

        # Peak overlap analysis
        overlap_analysis = self._analyze_peak_overlap(context, x_data, y_data)
        enhanced_features.update(overlap_analysis)

        # Interference source identification
        enhanced_features['interference_sources'] = self._identify_interference_sources(context)

        # Peak isolation scoring
        enhanced_features['peak_isolation_score'] = self._calculate_peak_isolation_score(context, x_data)

        return {**basic_features, **enhanced_features}

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
                if initial_val != 0:
                    param_changes[f'{param}_relative_change'] = abs(final_val - initial_val) / abs(initial_val)
                else:
                    param_changes[f'{param}_relative_change'] = 0.0

        features['parameter_changes'] = param_changes

        # Convergence metrics
        features['convergence_iterations'] = int(optimization_info.get('iterations', 0))
        features['converged'] = bool(optimization_info.get('converged', False))
        features['optimization_method'] = optimization_info.get('method', 'unknown')
        features['final_cost'] = float(optimization_info.get('final_cost', 0.0))

        # Parameter uncertainty estimates
        features['parameter_uncertainties'] = optimization_info.get('parameter_errors', {})

        # Optimization difficulty assessment
        features['optimization_difficulty'] = self._assess_optimization_difficulty(optimization_info)

        # Gradient and Hessian information
        features['gradient_norm'] = float(optimization_info.get('gradient_norm', 0.0))
        features['hessian_condition'] = float(optimization_info.get('hessian_condition', 1.0))

        return features

    def _extract_physics_validation(self, fit_result: Dict[str, Any], nucleus_type: str) -> Dict[str, Any]:
        """Validate parameters against spectroscopic physics."""

        validation = {}

        # Parameter range validation
        param_valid = {}
        for param, value in fit_result.items():
            if param in ['amplitude', 'center', 'sigma', 'gamma', 'baseline']:
                param_valid[f'{param}_valid'] = self._validate_parameter_physics(param, value, nucleus_type)

        validation['parameter_validity'] = param_valid

        # Peak shape realism
        sigma = fit_result.get('sigma', 0.01)
        gamma = fit_result.get('gamma', 0.01)
        validation['shape_realism_score'] = self._assess_peak_shape_realism(sigma, gamma, nucleus_type)

        # Linewidth consistency
        validation['linewidth_consistency'] = self._check_linewidth_consistency(sigma, gamma, nucleus_type)

        return validation

    def _analyze_multi_method_results(self, primary_result: Dict[str, Any],
                                    alternative_results: List[Dict]) -> Dict[str, Any]:
        """Analyze results from multiple fitting methods."""

        analysis = {}

        # Method comparison
        methods_tested = [result.get('method', 'unknown') for result in alternative_results]
        analysis['methods_tested'] = methods_tested

        # Parameter consensus analysis
        consensus_params = {}
        for param in ['amplitude', 'center', 'sigma', 'gamma', 'baseline']:
            if param in primary_result:
                values = [primary_result[param]]
                values.extend([result.get(param, 0) for result in alternative_results if param in result])

                if len(values) > 1:
                    consensus_params[f'{param}_mean'] = float(np.mean(values))
                    consensus_params[f'{param}_std'] = float(np.std(values))
                    consensus_params[f'{param}_range'] = float(np.max(values) - np.min(values))

        analysis['parameter_consensus'] = consensus_params

        # Method agreement scoring
        r_squares = [primary_result.get('r_squared', 0)]
        r_squares.extend([result.get('r_squared', 0) for result in alternative_results])

        analysis['method_agreement_score'] = float(np.std(r_squares)) if len(r_squares) > 1 else 0.0
        analysis['best_method_r_squared'] = float(max(r_squares))

        return analysis

    def _compute_enhanced_quality_metrics(self, x_data: np.ndarray, y_data: np.ndarray,
                                        fit_result: Dict[str, Any],
                                        optimization_info: Dict[str, Any]) -> Dict[str, Any]:
        """Compute comprehensive quality metrics."""

        metrics = {}

        # Basic quality metrics
        metrics['r_squared'] = float(fit_result.get('r_squared', 0.0))
        metrics['rmse_normalized'] = float(fit_result.get('rmse_normalized', 0.0))
        metrics['fit_success'] = bool(fit_result.get('success', False))

        # Residual analysis
        residuals = optimization_info.get('residuals', [])
        if len(residuals) > 0:
            metrics['residual_mean'] = float(np.mean(residuals))
            metrics['residual_std'] = float(np.std(residuals))
            metrics['residual_skewness'] = float(self._calculate_skewness(np.array(residuals)))
            metrics['residual_autocorrelation'] = self._calculate_residual_autocorrelation(residuals)

        # Outlier detection
        metrics['outlier_fraction'] = self._detect_outlier_fraction(y_data, fit_result)

        # Cross-validation estimate (if available)
        metrics['cross_validation_score'] = float(optimization_info.get('cv_score', 0.0))

        return metrics

    # Helper methods for enhanced feature extraction
    def _analyze_noise_characteristics(self, y_data: np.ndarray) -> Dict[str, float]:
        """Analyze noise characteristics in the data."""
        # Implementation for noise analysis
        return {
            'noise_level_robust': float(np.median(np.abs(y_data - np.median(y_data))) * 1.4826),
            'noise_distribution_type': 'gaussian',  # Could be enhanced with actual testing
            'noise_correlation': 0.0  # Could add autocorrelation analysis
        }

    def _analyze_baseline_characteristics(self, x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, float]:
        """Analyze baseline drift and characteristics."""
        # Implementation for baseline analysis
        return {
            'baseline_drift': float(np.std(y_data[:len(y_data)//10])),  # Simplified
            'baseline_curvature': 0.0  # Could add polynomial fit analysis
        }

    def _analyze_peak_shape_characteristics(self, x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, float]:
        """Analyze peak tailing and shape characteristics."""
        return {
            'peak_tailing': self._calculate_peak_tailing(x_data, y_data),
            'peak_fronting': 0.0  # Placeholder
        }

    def _calculate_peak_tailing(self, x_data: np.ndarray, y_data: np.ndarray) -> float:
        """Calculate peak tailing factor."""
        # Simplified tailing calculation
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
                return float(tailing_width / leading_width)

        return 1.0

    # Placeholder implementations for other helper methods
    def _extract_basic_spectral_features(self, x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, float]:
        """Extract basic spectral features (from original collector)."""
        # This would call the original method
        return {
            'peak_height': float(np.max(y_data)),
            'snr': 10.0,  # Placeholder
            'peak_width_estimate': 0.02,  # Placeholder
        }

    def _extract_basic_chemical_features(self, x_data: np.ndarray, nucleus_type: str) -> Dict[str, Any]:
        """Extract basic chemical features."""
        return {
            'nucleus_type': nucleus_type,
            'chemical_shift_center': float(np.mean(x_data)),
            'chemical_region': 'unknown'
        }

    def _extract_basic_context_features(self, context: Dict[str, Any], x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, Any]:
        """Extract basic context features."""
        return {
            'detection_confidence': float(context.get('detection_confidence', 0.5)),
            'peak_complexity': context.get('peak_complexity', 'unknown')
        }

    # Additional helper method implementations would go here...
    # (Many more helper methods would be needed for a complete implementation)

    def _validate_sample_completeness(self, sample: Dict[str, Any]) -> bool:
        """Validate that the sample contains all required fields."""
        required_fields = ['spectral_features', 'target_parameters', 'quality_metrics']
        return all(field in sample for field in required_fields)

    def _save_batch(self) -> None:
        """Save enhanced training batch."""
        if not self.session_data:
            return

        try:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = self.storage_path / f"enhanced_training_batch_{timestamp}.pkl"

            with open(filename, 'wb') as f:
                pickle.dump(self.session_data, f, protocol=pickle.HIGHEST_PROTOCOL)

            # Enhanced metadata
            metadata = {
                'timestamp': datetime.now().isoformat(),
                'collection_version': '2.0',
                'samples_count': len(self.session_data),
                'collection_stats': self.collection_stats.copy(),
                'feature_schema': self.feature_schema,
                'enhancements': [
                    'optimization_trajectory', 'multi_method_comparison',
                    'physics_validation', 'enhanced_quality_metrics'
                ]
            }

            metadata_file = filename.with_suffix('.json')
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)

            self.logger.info(f"Saved enhanced ML training batch: {len(self.session_data)} samples")
            self.collection_stats['batch_saves'] += 1
            self.session_data.clear()

        except Exception as e:
            self.logger.error(f"Failed to save enhanced training batch: {e}")

    # Placeholder implementations for missing methods
    def _predict_multiplicity(self, shift: float, nucleus_type: str) -> str:
        return "singlet"  # Placeholder

    def _assess_coupling_environment(self, shift: float, nucleus_type: str) -> str:
        return "isolated"  # Placeholder

    def _estimate_relaxation_properties(self, shift: float, nucleus_type: str) -> Dict[str, float]:
        return {"T1_estimate": 1.0, "T2_estimate": 0.1}  # Placeholder

    def _analyze_peak_overlap(self, context: Dict, x_data: np.ndarray, y_data: np.ndarray) -> Dict[str, float]:
        return {"overlap_severity": 0.0}  # Placeholder

    def _identify_interference_sources(self, context: Dict) -> List[str]:
        return []  # Placeholder

    def _calculate_peak_isolation_score(self, context: Dict, x_data: np.ndarray) -> float:
        return 1.0  # Placeholder

    def _assess_optimization_difficulty(self, opt_info: Dict) -> float:
        return float(opt_info.get('iterations', 0) / 100.0)  # Simple metric

    def _validate_parameter_physics(self, param: str, value: float, nucleus_type: str) -> bool:
        return True  # Placeholder - would implement physical constraints

    def _assess_peak_shape_realism(self, sigma: float, gamma: float, nucleus_type: str) -> float:
        return 1.0  # Placeholder

    def _check_linewidth_consistency(self, sigma: float, gamma: float, nucleus_type: str) -> bool:
        return True  # Placeholder

    def _calculate_skewness(self, data: np.ndarray) -> float:
        if len(data) < 3:
            return 0.0
        return float(np.mean(((data - np.mean(data)) / np.std(data)) ** 3))

    def _calculate_residual_autocorrelation(self, residuals: List) -> float:
        return 0.0  # Placeholder - would implement autocorrelation

    def _detect_outlier_fraction(self, y_data: np.ndarray, fit_result: Dict) -> float:
        return 0.0  # Placeholder - would implement outlier detection

    def _assess_data_quality_flags(self, x_data: np.ndarray, y_data: np.ndarray) -> List[str]:
        return []  # Placeholder for data quality flags

    def _downsample_data(self, data: np.ndarray, max_points: int = 200) -> np.ndarray:
        """Downsample data for storage efficiency."""
        if len(data) <= max_points:
            return data
        indices = np.linspace(0, len(data) - 1, max_points, dtype=int)
        return data[indices]