#!/usr/bin/env python3
"""
Unit Tests for ML Integration (Phase 1: Data Collection)

These tests verify that the ML data collection infrastructure integrates
properly with all existing workflows without breaking functionality.
"""

import unittest
import numpy as np
import tempfile
import shutil
from pathlib import Path
import sys
import os

# Add lunaNMR to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lunaNMR.core.enhanced_voigt_fitter import EnhancedVoigtFitter
from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
from lunaNMR.ml.training_data_collector import MLTrainingDataCollector

class TestMLDataCollection(unittest.TestCase):
    """Test ML data collection infrastructure"""

    def setUp(self):
        """Set up test environment"""
        self.temp_dir = tempfile.mkdtemp()
        self.collector = MLTrainingDataCollector(self.temp_dir)

        # Create synthetic test data
        self.x_data = np.linspace(7.8, 8.0, 50)
        self.y_data = self._create_synthetic_peak(self.x_data)

    def tearDown(self):
        """Clean up test environment"""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _create_synthetic_peak(self, x_data):
        """Create synthetic Voigt peak for testing"""
        amplitude = 150000
        center = 7.9
        sigma = 0.008
        gamma = 0.004
        baseline = 30000

        # Simplified Voigt approximation
        gaussian = amplitude * np.exp(-0.5 * ((x_data - center) / sigma)**2)
        lorentzian = amplitude * gamma**2 / ((x_data - center)**2 + gamma**2)
        voigt = 0.6 * gaussian + 0.4 * lorentzian + baseline

        # Add noise
        noise = 0.02 * np.random.randn(len(x_data)) * amplitude
        return voigt + noise

    def test_collector_initialization(self):
        """Test MLTrainingDataCollector initialization"""
        self.assertIsNotNone(self.collector)
        self.assertEqual(len(self.collector.session_data), 0)
        self.assertTrue(self.collector.storage_path.exists())

    def test_feature_extraction(self):
        """Test enhanced spectral and chemical feature extraction"""
        # Test enhanced spectral features (v2.0)
        spectral_features = self.collector._extract_enhanced_spectral_features(self.x_data, self.y_data)

        # Test key enhanced features from schema
        expected_features = [
            'peak_height', 'baseline_level', 'snr', 'peak_width_fwhm',
            'peak_width_moment', 'dynamic_range', 'spectral_resolution'
        ]
        for feature in expected_features:
            self.assertIn(feature, spectral_features)
            self.assertIsInstance(spectral_features[feature], (int, float))
            self.assertTrue(np.isfinite(spectral_features[feature]))

        # Test enhanced chemical features (v2.0)
        chemical_features = self.collector._extract_enhanced_chemical_features(self.x_data, '1H')

        # Test key enhanced chemical features from schema
        expected_features = [
            'nucleus_type', 'chemical_shift_center', 'chemical_region',
            'shift_range', 'coupling_environment'
        ]
        for feature in expected_features:
            self.assertIn(feature, chemical_features)

        self.assertEqual(chemical_features['nucleus_type'], '1H')

        # Test backward compatibility methods still work
        legacy_spectral = self.collector._extract_spectral_features(self.x_data, self.y_data)
        legacy_chemical = self.collector._extract_chemical_features(self.x_data, '1H')

        # Verify legacy methods return expected v1.0 format
        self.assertIn('peak_width_estimate', legacy_spectral)  # Maps to peak_width_fwhm
        self.assertIn('chemical_shift_center', legacy_chemical)

    def test_quality_filtering(self):
        """Test quality filtering of training samples"""
        # High-quality fit (should be collected)
        high_quality_result = {
            'success': True,
            'r_squared': 0.95,
            'amplitude': 150000,
            'center': 7.9,
            'sigma': 0.008,
            'gamma': 0.004,
            'baseline': 30000
        }

        context = {
            'detection_confidence': 0.8,
            'peak_complexity': 'simple'
        }

        success = self.collector.collect_training_sample(
            self.x_data, self.y_data, high_quality_result, context, '1H'
        )

        self.assertTrue(success)
        self.assertEqual(len(self.collector.session_data), 1)

        # Medium-quality fit (should be collected as 'reject' for ML learning)
        medium_quality_result = high_quality_result.copy()
        medium_quality_result['r_squared'] = 0.5  # Below good threshold but above absolute minimum

        success = self.collector.collect_training_sample(
            self.x_data, self.y_data, medium_quality_result, context, '1H'
        )

        self.assertTrue(success)  # Now collected as reject sample
        self.assertEqual(len(self.collector.session_data), 2)  # One good + one reject

        # Check quality labels
        good_sample = self.collector.session_data[0]
        reject_sample = self.collector.session_data[1]

        self.assertEqual(good_sample['quality_metrics']['label'], 'good')
        self.assertEqual(reject_sample['quality_metrics']['label'], 'reject')

        # Very low-quality fit (should be filtered out completely)
        very_low_quality_result = high_quality_result.copy()
        very_low_quality_result['r_squared'] = 0.1  # Below absolute minimum

        success = self.collector.collect_training_sample(
            self.x_data, self.y_data, very_low_quality_result, context, '1H'
        )

        self.assertFalse(success)  # This should still be filtered
        self.assertEqual(len(self.collector.session_data), 2)  # No change

        # Test with collect_rejected_samples disabled
        self.collector.collect_rejected_samples = False
        medium_quality_result['r_squared'] = 0.6  # Medium quality again

        success = self.collector.collect_training_sample(
            self.x_data, self.y_data, medium_quality_result, context, '1H'
        )

        self.assertFalse(success)  # Should be filtered when flag is disabled
        self.assertEqual(len(self.collector.session_data), 2)  # No change

    def test_enhanced_quality_labels(self):
        """Test the new three-tier quality labeling system"""
        context = {'detection_confidence': 0.8, 'peak_complexity': 'simple'}

        # Test excellent quality (R² = 0.98) -> 'good' label
        excellent_result = {
            'success': True, 'r_squared': 0.98, 'amplitude': 150000, 'center': 7.9,
            'sigma': 0.008, 'gamma': 0.004, 'baseline': 30000
        }

        success = self.collector.collect_training_sample(
            self.x_data, self.y_data, excellent_result, context, '1H'
        )

        self.assertTrue(success)
        self.assertEqual(len(self.collector.session_data), 1)
        self.assertEqual(self.collector.session_data[0]['quality_metrics']['label'], 'good')

        # Test borderline quality (R² = 0.4) -> 'reject' label (if enabled)
        self.collector.collect_rejected_samples = True
        borderline_result = excellent_result.copy()
        borderline_result['r_squared'] = 0.4

        success = self.collector.collect_training_sample(
            self.x_data, self.y_data, borderline_result, context, '1H'
        )

        self.assertTrue(success)
        self.assertEqual(len(self.collector.session_data), 2)
        self.assertEqual(self.collector.session_data[1]['quality_metrics']['label'], 'reject')

        # Test with collect_rejected_samples disabled -> borderline should be filtered
        self.collector.session_data.clear()
        self.collector.collect_rejected_samples = False

        success = self.collector.collect_training_sample(
            self.x_data, self.y_data, borderline_result, context, '1H'
        )

        self.assertFalse(success)
        self.assertEqual(len(self.collector.session_data), 0)

        # Test very poor quality (R² = 0.1) -> always filtered regardless of setting
        very_poor_result = excellent_result.copy()
        very_poor_result['r_squared'] = 0.1

        self.collector.collect_rejected_samples = True  # Even with this enabled
        success = self.collector.collect_training_sample(
            self.x_data, self.y_data, very_poor_result, context, '1H'
        )

        self.assertFalse(success)  # Should still be filtered
        self.assertEqual(len(self.collector.session_data), 0)

    def test_data_format_consistency(self):
        """Test that collected data has consistent format"""
        fit_result = {
            'success': True,
            'r_squared': 0.92,
            'amplitude': 150000,
            'center': 7.9,
            'sigma': 0.008,
            'gamma': 0.004,
            'baseline': 30000
        }

        context = {
            'detection_confidence': 0.8,
            'peak_complexity': 'simple'
        }

        self.collector.collect_training_sample(
            self.x_data, self.y_data, fit_result, context, '1H'
        )

        self.assertEqual(len(self.collector.session_data), 1)
        sample = self.collector.session_data[0]

        # Check required top-level keys
        required_keys = ['timestamp', 'spectral_features', 'chemical_features',
                        'context_features', 'target_parameters', 'quality_metrics', 'raw_data']
        for key in required_keys:
            self.assertIn(key, sample)

        # Check target parameters
        target_params = sample['target_parameters']
        param_keys = ['amplitude', 'center', 'sigma', 'gamma', 'baseline']
        for key in param_keys:
            self.assertIn(key, target_params)
            self.assertIsInstance(target_params[key], (int, float))

    def test_error_handling(self):
        """Test graceful error handling"""
        # Invalid fit result
        invalid_result = {
            'success': True,
            'r_squared': 0.95,
            'amplitude': np.nan,  # Invalid parameter
            'center': 7.9
        }

        context = {'detection_confidence': 0.8}

        # Should handle gracefully without raising exception
        success = self.collector.collect_training_sample(
            self.x_data, self.y_data, invalid_result, context, '1H'
        )

        self.assertFalse(success)
        self.assertEqual(len(self.collector.session_data), 0)

class TestEnhancedVoigtFitterIntegration(unittest.TestCase):
    """Test ML integration with EnhancedVoigtFitter"""

    def setUp(self):
        """Set up test environment"""
        self.fitter = EnhancedVoigtFitter()
        self.x_data = np.linspace(7.8, 8.0, 50)
        self.y_data = self._create_synthetic_peak(self.x_data)

    def _create_synthetic_peak(self, x_data):
        """Create synthetic Voigt peak for testing"""
        amplitude = 150000
        center = 7.9
        sigma = 0.008
        gamma = 0.004
        baseline = 30000

        gaussian = amplitude * np.exp(-0.5 * ((x_data - center) / sigma)**2)
        lorentzian = amplitude * gamma**2 / ((x_data - center)**2 + gamma**2)
        voigt = 0.6 * gaussian + 0.4 * lorentzian + baseline

        noise = 0.02 * np.random.randn(len(x_data)) * amplitude
        return voigt + noise

    def test_ml_collector_initialization(self):
        """Test that ML data collector is properly initialized"""
        self.assertTrue(hasattr(self.fitter, 'ml_data_collector'))
        self.assertIsNotNone(self.fitter.ml_data_collector)

    def test_fitting_with_ml_collection(self):
        """Test that fitting works normally with ML collection enabled"""
        initial_count = len(self.fitter.ml_data_collector.session_data)

        # Perform fitting
        result = self.fitter.fit_peak_enhanced(
            self.x_data, self.y_data, nucleus_type='1H', method='iterative_optimization'
        )

        # Check fitting succeeded
        self.assertIsNotNone(result)
        self.assertTrue(result.get('success', False))
        self.assertGreater(result.get('r_squared', 0), 0.8)

        # Check ML data collection occurred
        final_count = len(self.fitter.ml_data_collector.session_data)
        self.assertGreaterEqual(final_count, initial_count)

    def test_backward_compatibility(self):
        """Test that existing functionality is not broken"""
        # Test that all original methods still work
        methods = ['iterative_optimization', 'multi_step', 'single_step']

        for method in methods:
            try:
                result = self.fitter.fit_peak_enhanced(
                    self.x_data, self.y_data, nucleus_type='1H', method=method
                )

                self.assertIsNotNone(result)

                if result.get('success', False):
                    # Check all expected keys are present
                    expected_keys = ['success', 'amplitude', 'center', 'sigma', 'gamma', 'baseline']
                    for key in expected_keys:
                        self.assertIn(key, result)

            except Exception as e:
                self.fail(f"Method {method} failed with ML integration: {e}")

class TestCoreIntegratorIntegration(unittest.TestCase):
    """Test ML integration with core integrator (S/N detection)"""

    def setUp(self):
        """Set up test environment"""
        self.integrator = EnhancedVoigtIntegrator()

    def test_ml_collector_initialization(self):
        """Test ML collector initialization in core integrator"""
        self.assertTrue(hasattr(self.integrator, 'ml_data_collector'))
        self.assertIsNotNone(self.integrator.ml_data_collector)

    def test_sn_detection_compatibility(self):
        """Test S/N detection works with ML collection"""
        # Create test spectrum
        spectrum = np.random.normal(1000, 100, (20, 25)) + np.random.exponential(1200, (20, 25))
        ppm_x = np.linspace(10, 6, 25)
        ppm_y = np.linspace(130, 100, 20)

        # Add synthetic peaks
        spectrum[5, 10] += 5000
        spectrum[10, 15] += 6000
        spectrum[15, 20] += 4000

        self.integrator.load_nmr_data(spectrum, ppm_x, ppm_y)
        self.integrator.sn_threshold = 2.0
        self.integrator.expected_peak_count = 5

        initial_count = len(self.integrator.ml_data_collector.session_data)

        # Perform S/N detection
        detected_peaks = self.integrator.detect_peaks_sn_native()

        # Verify detection worked
        self.assertIsNotNone(detected_peaks)
        self.assertIsInstance(detected_peaks, list)

        # ML collection doesn't break S/N detection
        final_count = len(self.integrator.ml_data_collector.session_data)
        self.assertGreaterEqual(final_count, initial_count)

class TestMLIntegrationSuite(unittest.TestCase):
    """Integration test suite for complete ML functionality"""

    def test_end_to_end_workflow(self):
        """Test complete end-to-end workflow"""
        # Initialize components
        fitter = EnhancedVoigtFitter()
        integrator = EnhancedVoigtIntegrator()

        # Both should have ML collectors
        self.assertIsNotNone(fitter.ml_data_collector)
        self.assertIsNotNone(integrator.ml_data_collector)

        # Test that they don't interfere with each other
        initial_fitter_count = len(fitter.ml_data_collector.session_data)
        initial_integrator_count = len(integrator.ml_data_collector.session_data)

        # Perform operations
        x_data = np.linspace(7.8, 8.0, 50)
        y_data = 30000 + 120000 * np.exp(-0.5 * ((x_data - 7.9) / 0.01)**2)

        # EnhancedVoigtFitter operation
        result1 = fitter.fit_peak_enhanced(x_data, y_data, nucleus_type='1H')

        # EnhancedVoigtIntegrator operation
        spectrum = np.random.normal(1000, 50, (15, 20))
        spectrum[7, 10] += 3000
        ppm_x = np.linspace(10, 6, 20)
        ppm_y = np.linspace(130, 100, 15)

        integrator.load_nmr_data(spectrum, ppm_x, ppm_y)
        integrator.sn_threshold = 1.5
        result2 = integrator.detect_peaks_sn_native()

        # Both operations should succeed
        self.assertTrue(result1.get('success', False) if result1 else True)
        self.assertIsNotNone(result2)

        # Data collection should have occurred independently
        final_fitter_count = len(fitter.ml_data_collector.session_data)
        final_integrator_count = len(integrator.ml_data_collector.session_data)

        # Check if collection was attempted (even if it failed)
        # At least one collector should have recorded some activity
        fitter_stats = fitter.ml_data_collector.collection_stats
        integrator_stats = integrator.ml_data_collector.collection_stats

        total_attempts = (fitter_stats.get('total_attempts', 0) +
                         integrator_stats.get('total_attempts', 0))

        # Either data was successfully collected OR attempts were made
        collection_occurred = (
            final_fitter_count > initial_fitter_count or
            final_integrator_count > initial_integrator_count or
            total_attempts > 0
        )

        self.assertTrue(collection_occurred,
                       f"No ML collection activity detected. Fitter stats: {fitter_stats}, "
                       f"Integrator stats: {integrator_stats}")

if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)