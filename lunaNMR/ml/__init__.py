"""
Machine Learning Module for LunaNMR

This module provides ML enhancement capabilities for peak fitting and analysis.
Currently implements training data collection for future ML model development.

Phase 1: Data Collection Infrastructure
- Transparent collection of high-quality fit results
- Feature extraction for spectral, chemical, and context data
- Seamless integration with all existing workflows
- Zero impact on user experience
"""

# Version and metadata
__version__ = "0.1.0"
__phase__ = "Data Collection"

# Import core components when available
try:
    from .training_data_collector import MLTrainingDataCollector
    ML_COLLECTION_AVAILABLE = True
except ImportError:
    ML_COLLECTION_AVAILABLE = False

__all__ = ['MLTrainingDataCollector', 'ML_COLLECTION_AVAILABLE']