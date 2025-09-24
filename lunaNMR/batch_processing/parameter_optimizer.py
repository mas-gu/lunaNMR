#!/usr/bin/env python3
"""
Parameter Optimization Engine for Batch Processing

This module provides intelligent parameter optimization for different spectrum types
and characteristics. It automatically tunes S/N thresholds and other parameters
to maximize the quality and quantity of ML training data collected.

Features:
- Spectrum quality assessment
- Adaptive S/N threshold optimization
- Peak count targeting for optimal training data
- Parameter feedback loops based on fitting success
"""

import numpy as np
import logging
from typing import Dict, List, Tuple, Any, Optional
from pathlib import Path

class ParameterOptimizer:
    """
    Intelligent parameter optimization for batch NMR processing.

    This class analyzes spectrum characteristics and optimizes processing
    parameters to maximize the quality and quantity of ML training data.
    """

    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize the parameter optimizer.

        Args:
            config: Optional configuration dictionary
        """
        self.logger = logging.getLogger(__name__)

        # Default optimization configuration
        self.config = {
            'sn_test_range': {
                '15N1H': [1, 1.5, 2.0, 3.0, 5.0],
                '13C1H': [1.0, 2.0, 3.0, 4.0, 5.0],
                '1H': [1.5, 1.8, 2.0, 2.5],
                'default': [1.2, 1.5, 1.8, 2.0, 2.5]
            },
            'optimal_peak_ranges': {
                '15N1H': (100, 300),    # Target range for 15N1H spectra
                '13C1H': (80, 200),    # Target range for 13C1H spectra
                '1H': (15, 40),       # Target range for 1H spectra
                'default': (100, 400)
            },
            'quality_thresholds': {
                'excellent': 0.95,
                'good': 0.85,
                'acceptable': 0.75,
                'poor': 0.60
            },
            'optimization_strategy': 'balanced',  # 'conservative', 'balanced', 'aggressive'
            'max_optimization_attempts': 5,
            'noise_estimation_method': 'corners'
        }

        # Update with provided config using deep merge
        if config:
            self.config = self._merge_configs(self.config, config)

        # Optimization history for learning
        self.optimization_history = []

    def _merge_configs(self, base_config: Dict[str, Any], override_config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Recursively merge two configuration dictionaries.

        Args:
            base_config: Base configuration dictionary
            override_config: Configuration to merge/override with

        Returns:
            Merged configuration dictionary
        """
        merged = base_config.copy()

        for key, value in override_config.items():
            if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
                # Recursively merge nested dictionaries
                merged[key] = self._merge_configs(merged[key], value)
            else:
                # Override value
                merged[key] = value

        return merged

    def analyze_spectrum_characteristics(self, spectrum_data: np.ndarray,
                                       ppm_x: np.ndarray,
                                       ppm_y: np.ndarray,
                                       nucleus_type: str) -> Dict[str, Any]:
        """
        Analyze spectrum characteristics to guide parameter optimization.

        Args:
            spectrum_data: 2D spectrum data
            ppm_x: X-axis chemical shift values
            ppm_y: Y-axis chemical shift values
            nucleus_type: Type of nucleus

        Returns:
            Dictionary with spectrum characteristics
        """
        characteristics = {
            'nucleus_type': nucleus_type,
            'shape': spectrum_data.shape,
            'dynamic_range': 0,
            'noise_level': 0,
            'signal_level': 0,
            'snr_estimate': 0,
            'peak_density_estimate': 0,
            'baseline_stability': 0,
            'quality_assessment': 'unknown'
        }

        try:
            # Basic spectrum statistics
            characteristics['dynamic_range'] = np.max(spectrum_data) - np.min(spectrum_data)

            # Noise level estimation using corner regions
            characteristics['noise_level'] = self._estimate_noise_level(spectrum_data)

            # Signal level estimation
            characteristics['signal_level'] = np.percentile(spectrum_data, 95)

            # SNR estimate
            if characteristics['noise_level'] > 0:
                characteristics['snr_estimate'] = characteristics['signal_level'] / characteristics['noise_level']
            else:
                characteristics['snr_estimate'] = float('inf')

            # Peak density estimation (rough)
            characteristics['peak_density_estimate'] = self._estimate_peak_density(spectrum_data)

            # Baseline stability assessment
            characteristics['baseline_stability'] = self._assess_baseline_stability(spectrum_data)

            # Overall quality assessment
            characteristics['quality_assessment'] = self._assess_overall_quality(characteristics)

            self.logger.debug(f"Spectrum characteristics: {characteristics}")

        except Exception as e:
            self.logger.warning(f"Error analyzing spectrum characteristics: {e}")

        return characteristics

    def _estimate_noise_level(self, spectrum_data: np.ndarray) -> float:
        """
        Estimate noise level using corner regions of the spectrum.

        Args:
            spectrum_data: 2D spectrum data

        Returns:
            Estimated noise level
        """
        try:
            h, w = spectrum_data.shape
            corner_size = min(h, w) // 8  # Use corner regions

            if corner_size < 3:
                # Fallback for very small spectra
                return np.std(spectrum_data) * 0.1

            # Extract corner regions
            corners = [
                spectrum_data[:corner_size, :corner_size],           # Top-left
                spectrum_data[:corner_size, -corner_size:],          # Top-right
                spectrum_data[-corner_size:, :corner_size],          # Bottom-left
                spectrum_data[-corner_size:, -corner_size:]          # Bottom-right
            ]

            # Calculate noise as standard deviation of corner regions
            corner_data = np.concatenate([corner.flatten() for corner in corners])
            noise_level = np.std(corner_data)

            return float(noise_level)

        except Exception as e:
            self.logger.warning(f"Error estimating noise level: {e}")
            return float(np.std(spectrum_data) * 0.1)

    def _estimate_peak_density(self, spectrum_data: np.ndarray) -> float:
        """
        Estimate peak density in the spectrum.

        Args:
            spectrum_data: 2D spectrum data

        Returns:
            Estimated peak density (peaks per unit area)
        """
        try:
            # Simple peak density estimation based on local maxima
            # This is a rough estimate for parameter optimization

            # Use median + 3*std as a threshold for potential peaks
            threshold = np.median(spectrum_data) + 3 * np.std(spectrum_data)

            # Count pixels above threshold
            above_threshold = np.sum(spectrum_data > threshold)

            # Calculate density relative to spectrum size
            total_pixels = spectrum_data.size
            density = above_threshold / total_pixels

            return float(density)

        except Exception as e:
            self.logger.warning(f"Error estimating peak density: {e}")
            return 0.01  # Default low density

    def _assess_baseline_stability(self, spectrum_data: np.ndarray) -> float:
        """
        Assess baseline stability across the spectrum.

        Args:
            spectrum_data: 2D spectrum data

        Returns:
            Baseline stability score (0-1, higher is more stable)
        """
        try:
            # Sample baseline regions (bottom 10th percentile)
            baseline_threshold = np.percentile(spectrum_data, 10)
            baseline_points = spectrum_data[spectrum_data <= baseline_threshold]

            if len(baseline_points) < 10:
                return 0.5  # Default moderate stability

            # Calculate coefficient of variation for baseline
            baseline_std = np.std(baseline_points)
            baseline_mean = np.mean(baseline_points)

            if baseline_mean > 0:
                cv = baseline_std / baseline_mean
                # Convert to stability score (lower CV = higher stability)
                stability = np.exp(-cv)  # Exponential decay from CV
                return float(np.clip(stability, 0, 1))
            else:
                return 0.5

        except Exception as e:
            self.logger.warning(f"Error assessing baseline stability: {e}")
            return 0.5

    def _assess_overall_quality(self, characteristics: Dict[str, Any]) -> str:
        """
        Assess overall spectrum quality based on characteristics.

        Args:
            characteristics: Dictionary with spectrum characteristics

        Returns:
            Quality assessment string
        """
        try:
            score = 0
            factors = 0

            # SNR contribution
            if characteristics['snr_estimate'] > 10:
                score += 3
            elif characteristics['snr_estimate'] > 5:
                score += 2
            elif characteristics['snr_estimate'] > 2:
                score += 1
            factors += 3

            # Baseline stability contribution
            if characteristics['baseline_stability'] > 0.8:
                score += 2
            elif characteristics['baseline_stability'] > 0.6:
                score += 1
            factors += 2

            # Dynamic range contribution
            if characteristics['dynamic_range'] > 50000:
                score += 2
            elif characteristics['dynamic_range'] > 10000:
                score += 1
            factors += 2

            # Calculate overall score
            if factors > 0:
                overall_score = score / factors
            else:
                overall_score = 0.5

            # Map to quality categories
            if overall_score >= 0.9:
                return 'excellent'
            elif overall_score >= 0.7:
                return 'good'
            elif overall_score >= 0.5:
                return 'acceptable'
            else:
                return 'poor'

        except Exception as e:
            self.logger.warning(f"Error assessing overall quality: {e}")
            return 'unknown'

    def optimize_sn_threshold(self, integrator, characteristics: Dict[str, Any]) -> float:
        """
        Optimize S/N threshold based on spectrum characteristics.

        Args:
            integrator: The EnhancedVoigtIntegrator instance
            characteristics: Dictionary with spectrum characteristics

        Returns:
            Optimized S/N threshold
        """
        nucleus_type = characteristics['nucleus_type']
        test_thresholds = self.config['sn_test_range'].get(nucleus_type,
                                                           self.config['sn_test_range']['default'])
        optimal_range = self.config['optimal_peak_ranges'].get(nucleus_type,
                                                              self.config['optimal_peak_ranges']['default'])

        best_threshold = test_thresholds[len(test_thresholds)//2]  # Default to middle value
        best_score = 0

        self.logger.debug(f"Optimizing S/N threshold for {nucleus_type}, target range: {optimal_range}")

        for threshold in test_thresholds:
            try:
                # Set threshold and test detection
                integrator.sn_threshold = threshold
                detected_peaks = integrator.detect_peaks_sn_native()

                if detected_peaks is None:
                    continue

                peak_count = len(detected_peaks)

                # Score based on how close we are to optimal range
                score = self._score_peak_count(peak_count, optimal_range, characteristics)

                self.logger.debug(f"  Threshold {threshold}: {peak_count} peaks, score: {score:.3f}")

                if score > best_score:
                    best_score = score
                    best_threshold = threshold

            except Exception as e:
                self.logger.debug(f"  Threshold {threshold} failed: {e}")
                continue

        self.logger.info(f"Optimized S/N threshold: {best_threshold} (score: {best_score:.3f})")
        return best_threshold

    def _score_peak_count(self, peak_count: int, optimal_range: Tuple[int, int],
                         characteristics: Dict[str, Any]) -> float:
        """
        Score a peak count based on how well it fits the optimal range.

        Args:
            peak_count: Number of detected peaks
            optimal_range: Tuple of (min, max) optimal peak counts
            characteristics: Spectrum characteristics

        Returns:
            Score (0-1, higher is better)
        """
        min_peaks, max_peaks = optimal_range

        if min_peaks <= peak_count <= max_peaks:
            # Peak count is in optimal range
            # Score based on position within range (favor middle)
            range_size = max_peaks - min_peaks
            if range_size > 0:
                position = (peak_count - min_peaks) / range_size
                # Quadratic score favoring middle of range
                score = 1.0 - 4 * (position - 0.5) ** 2
            else:
                score = 1.0
        elif peak_count < min_peaks:
            # Too few peaks - score based on how close we are
            if min_peaks > 0:
                score = peak_count / min_peaks * 0.7  # Max 70% score for too few
            else:
                score = 0.5
        else:  # peak_count > max_peaks
            # Too many peaks - score based on excess
            excess_factor = (peak_count - max_peaks) / max_peaks
            score = max(0.1, 0.8 * np.exp(-excess_factor))  # Exponential decay

        # Adjust score based on spectrum quality
        quality_multiplier = {
            'excellent': 1.0,
            'good': 0.95,
            'acceptable': 0.85,
            'poor': 0.7,
            'unknown': 0.8
        }.get(characteristics.get('quality_assessment', 'unknown'), 0.8)

        return float(score * quality_multiplier)

    def optimize_expected_peak_count(self, detected_peak_count: int,
                                   nucleus_type: str,
                                   characteristics: Dict[str, Any]) -> int:
        """
        Optimize expected peak count based on detected peaks and spectrum characteristics.

        Args:
            detected_peak_count: Number of peaks detected with current parameters
            nucleus_type: Type of nucleus
            characteristics: Spectrum characteristics

        Returns:
            Optimized expected peak count
        """
        optimal_range = self.config['optimal_peak_ranges'].get(nucleus_type,
                                                              self.config['optimal_peak_ranges']['default'])

        # Strategy based on optimization approach
        strategy = self.config['optimization_strategy']

        if strategy == 'conservative':
            # Use lower end of optimal range
            target = min(detected_peak_count, optimal_range[1])
        elif strategy == 'aggressive':
            # Use higher end, potentially above detected
            target = max(detected_peak_count, optimal_range[0])
        else:  # balanced
            # Target middle of optimal range, adjusted for detected count
            target_middle = (optimal_range[0] + optimal_range[1]) // 2
            if detected_peak_count < optimal_range[0]:
                target = optimal_range[0]
            elif detected_peak_count > optimal_range[1]:
                target = min(detected_peak_count, optimal_range[1] + 10)
            else:
                target = min(detected_peak_count, target_middle)

        self.logger.debug(f"Expected peak count optimization: detected={detected_peak_count}, "
                         f"target={target}, strategy={strategy}")

        return target

    def optimize_parameters(self, integrator, spectrum_data: np.ndarray,
                          ppm_x: np.ndarray, ppm_y: np.ndarray,
                          nucleus_type: str) -> Dict[str, Any]:
        """
        Perform comprehensive parameter optimization for a spectrum.

        Args:
            integrator: The EnhancedVoigtIntegrator instance
            spectrum_data: 2D spectrum data
            ppm_x: X-axis chemical shift values
            ppm_y: Y-axis chemical shift values
            nucleus_type: Type of nucleus

        Returns:
            Dictionary with optimized parameters
        """
        self.logger.info(f"Starting parameter optimization for {nucleus_type} spectrum")

        try:
            # Analyze spectrum characteristics
            characteristics = self.analyze_spectrum_characteristics(
                spectrum_data, ppm_x, ppm_y, nucleus_type
            )

            # Optimize S/N threshold
            optimal_sn_threshold = self.optimize_sn_threshold(integrator, characteristics)

            # Test detection with optimized threshold
            integrator.sn_threshold = optimal_sn_threshold
            detected_peaks = integrator.detect_peaks_sn_native()
            detected_count = len(detected_peaks) if detected_peaks else 0

            # Optimize expected peak count
            optimal_expected_peaks = self.optimize_expected_peak_count(
                detected_count, nucleus_type, characteristics
            )

            # Compile optimized parameters
            optimized_params = {
                'nucleus_type': nucleus_type,
                'sn_threshold': optimal_sn_threshold,
                'expected_peak_count': optimal_expected_peaks,
                'characteristics': characteristics,
                'optimization_method': 'adaptive',
                'detected_peaks_with_params': detected_count
            }

            # Record optimization history for learning
            self.optimization_history.append({
                'nucleus_type': nucleus_type,
                'characteristics': characteristics,
                'parameters': optimized_params,
                'success': detected_count > 0
            })

            self.logger.info(f"Parameter optimization complete: S/N={optimal_sn_threshold}, "
                           f"expected={optimal_expected_peaks}, detected={detected_count}")

            return optimized_params

        except Exception as e:
            self.logger.error(f"Parameter optimization failed: {e}")
            # Return safe default parameters
            return {
                'nucleus_type': nucleus_type,
                'sn_threshold': 2.2,
                'expected_peak_count': 50,
                'characteristics': {},
                'optimization_method': 'default_fallback',
                'detected_peaks_with_params': 0
            }

    def get_optimization_statistics(self) -> Dict[str, Any]:
        """
        Get statistics about parameter optimization performance.

        Returns:
            Dictionary with optimization statistics
        """
        if not self.optimization_history:
            return {
                'total_optimizations': 0,
                'success_rate': 0,
                'average_detected_peaks': 0,
                'nucleus_type_breakdown': {}
            }

        total = len(self.optimization_history)
        successful = sum(1 for opt in self.optimization_history if opt['success'])

        # Calculate average detected peaks
        detected_counts = [opt['parameters']['detected_peaks_with_params']
                         for opt in self.optimization_history if opt['success']]
        avg_detected = np.mean(detected_counts) if detected_counts else 0

        # Nucleus type breakdown
        nucleus_breakdown = {}
        for opt in self.optimization_history:
            nucleus = opt['nucleus_type']
            if nucleus not in nucleus_breakdown:
                nucleus_breakdown[nucleus] = {'total': 0, 'successful': 0}
            nucleus_breakdown[nucleus]['total'] += 1
            if opt['success']:
                nucleus_breakdown[nucleus]['successful'] += 1

        return {
            'total_optimizations': total,
            'success_rate': successful / total if total > 0 else 0,
            'average_detected_peaks': float(avg_detected),
            'nucleus_type_breakdown': nucleus_breakdown
        }
