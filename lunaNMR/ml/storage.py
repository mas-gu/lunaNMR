# ABOUTME: Local storage paths for ML training data within the project
# ABOUTME: Provides unified storage location for GUI and CLI data collection

from pathlib import Path
from typing import Optional


def _get_project_root() -> Path:
    """Get the lunaNMR_v1o0 project root directory."""
    # This file is at: lunaNMR_v1o0/lunaNMR/ml/storage.py
    # Project root is 3 levels up
    return Path(__file__).parent.parent.parent


def get_ml_data_dir(custom_dir: Optional[Path] = None) -> Path:
    """
    Get the directory for ML training data.

    Parameters
    ----------
    custom_dir : Path, optional
        Custom directory to use instead of default

    Returns
    -------
    Path
        Directory path for ML data storage

    Notes
    -----
    Default location: lunaNMR_v1o0/ml_training_data/
    """
    if custom_dir is not None:
        ml_dir = Path(custom_dir)
        ml_dir.mkdir(parents=True, exist_ok=True)
        return ml_dir

    # Store locally within the project
    project_root = _get_project_root()
    ml_dir = project_root / "ml_training_data"
    ml_dir.mkdir(parents=True, exist_ok=True)

    return ml_dir


def get_training_data_path(custom_dir: Optional[Path] = None) -> Path:
    """Get path to the main training data JSON file."""
    return get_ml_data_dir(custom_dir) / "training_data.json"


def get_session_log_dir(custom_dir: Optional[Path] = None) -> Path:
    """Get path to session logs directory."""
    log_dir = get_ml_data_dir(custom_dir) / "session_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir


def get_stats_history_path(custom_dir: Optional[Path] = None) -> Path:
    """Get path to statistical predictor history file."""
    return get_ml_data_dir(custom_dir) / "stats_history.json"


def get_metadata_path(custom_dir: Optional[Path] = None) -> Path:
    """Get path to training metadata file."""
    return get_ml_data_dir(custom_dir) / "training_metadata.json"
