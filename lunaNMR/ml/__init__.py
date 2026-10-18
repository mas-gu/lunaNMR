"""
Machine Learning Module for LunaNMR

This module provides ML enhancement capabilities for peak fitting and analysis,
including parameter prediction, statistical fallback, and adaptive learning.
"""

# Core data structures
from .feature_extractor import SpectrumFeatures, OptimalParameters, FeatureExtractor

# Predictors
from .ml_predictor import MLPredictor
from .stats_predictor import StatisticalPredictor

# Decision and management
from .decision_engine import DecisionEngine
from .model_manager import ModelManager

# Local adaptation
from .training_collector import TrainingCollector
from .adaptation_engine import AdaptationEngine

# Peak classifier for ML-based peak detection
from .peak_feature_extractor import PeakFeatures, PeakFeatureExtractor
from .peak_classifier import PeakClassifier, MLPeakDetector

# Version and metadata
__version__ = "1.0.0"

__all__ = [
    # Data structures
    'SpectrumFeatures',
    'OptimalParameters',
    'FeatureExtractor',
    # Predictors
    'MLPredictor',
    'StatisticalPredictor',
    # Decision and management
    'DecisionEngine',
    'ModelManager',
    # Training and adaptation
    'TrainingCollector',
    'AdaptationEngine',
    # Peak classifier
    'PeakFeatures',
    'PeakFeatureExtractor',
    'PeakClassifier',
    'MLPeakDetector',
]