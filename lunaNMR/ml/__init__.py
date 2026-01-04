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

# Existing training collector (PS2D)
from .ps2d_training_collector import PS2DTrainingDataCollector

# Local adaptation
from .training_collector import TrainingCollector
from .adaptation_engine import AdaptationEngine

# Comprehensive training data collection
from .storage import get_ml_data_dir, get_training_data_path
from .comprehensive_collector import (
    ComprehensiveTrainingCollector,
    ComprehensiveTrainingSample,
    PeakTrainingData,
)

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
    'PS2DTrainingDataCollector',
    'TrainingCollector',
    'AdaptationEngine',
    # Comprehensive collection
    'get_ml_data_dir',
    'get_training_data_path',
    'ComprehensiveTrainingCollector',
    'ComprehensiveTrainingSample',
    'PeakTrainingData',
    # Peak classifier
    'PeakFeatures',
    'PeakFeatureExtractor',
    'PeakClassifier',
    'MLPeakDetector',
]