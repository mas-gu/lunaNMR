# ABOUTME: Training data collector for PS2D 2D multi-peak fitting results
# ABOUTME: Stores successful 2D Voigt fits as JSON for future ML model development

import json
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional
import logging


class PS2DTrainingDataCollector:
    """
    Collects training data from PS2D 2D multi-peak fitting results.

    Stores successful 2D Voigt fits with their parameters, quality metrics,
    and spectrum metadata for future machine learning model development.
    """

    def __init__(self, output_dir: Optional[str] = None):
        """
        Initialize the PS2D training data collector.

        Parameters:
        -----------
        output_dir : str, optional
            Directory to save training data. If None, uses './ml_training_data'
        """
        self.logger = logging.getLogger(__name__)

        # Set output directory
        if output_dir is None:
            output_dir = Path.cwd() / 'ml_training_data'
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Session data storage
        self.session_data = []
        self.session_start_time = datetime.now()

        self.logger.info(f"PS2D training data collector initialized: {self.output_dir}")

    def collect_ps2d_fit(self,
                         fit_result: Dict[str, Any],
                         spectrum_metadata: Optional[Dict[str, Any]] = None) -> bool:
        """
        Collect a PS2D fit result for training data.

        Parameters:
        -----------
        fit_result : dict
            PS2D fit result dictionary containing:
            - 'success': bool
            - 'pragmatic_acceptance': bool
            - 'n_peaks': int
            - 'peaks': list of peak parameter dicts
            - 'r_squared': float
            - 'chi2': float
            - 'iterations': int
        spectrum_metadata : dict, optional
            Additional spectrum information:
            - 'spectrum_name': str
            - 'nucleus_f1': str (e.g., '15N')
            - 'nucleus_f2': str (e.g., '1H')
            - 'temperature': float
            - 'field_strength': float
            - etc.

        Returns:
        --------
        bool : True if data was collected, False if rejected
        """
        # Only collect successful fits
        if not fit_result.get('success', False):
            self.logger.debug("Rejecting fit: not successful")
            return False

        # Ensure minimum quality threshold (R² > 0.3 as per PS2D algorithm)
        r_squared = fit_result.get('r_squared', 0.0)
        if r_squared < 0.3:
            self.logger.debug(f"Rejecting fit: R² = {r_squared:.3f} < 0.3")
            return False

        # Extract peak parameters (convert numpy arrays to lists for JSON)
        peaks_data = []
        for peak in fit_result.get('peaks', []):
            peak_dict = {
                'pos_f1': float(peak['pos_f1']),
                'lw_lor_f1': float(peak['lw_lor_f1']),
                'lw_gau_f1': float(peak['lw_gau_f1']),
                'pos_f2': float(peak['pos_f2']),
                'lw_lor_f2': float(peak['lw_lor_f2']),
                'lw_gau_f2': float(peak['lw_gau_f2']),
                'intensity': float(peak['intensity']),
                'volume': float(peak.get('volume', 0.0)),
                'height': float(peak.get('height', 0.0)),
                'amplitude': float(peak.get('amplitude', 0.0))
            }
            peaks_data.append(peak_dict)

        # Create training sample
        training_sample = {
            'timestamp': datetime.now().isoformat(),
            'fit_quality': {
                'r_squared': float(r_squared),
                'chi2': float(fit_result.get('chi2', 0.0)),
                'iterations': int(fit_result.get('iterations', 0)),
                'formal_convergence': bool(fit_result.get('formal_convergence', False)),
                'pragmatic_acceptance': True
            },
            'peaks': peaks_data,
            'n_peaks': len(peaks_data),
            'metadata': spectrum_metadata if spectrum_metadata else {}
        }

        # Add to session data
        self.session_data.append(training_sample)

        self.logger.info(f"Collected PS2D fit: {len(peaks_data)} peaks, R² = {r_squared:.3f}")
        return True

    def save_session(self, filename: Optional[str] = None) -> Path:
        """
        Save collected training data to JSON file.

        Parameters:
        -----------
        filename : str, optional
            Output filename. If None, generates timestamped filename.

        Returns:
        --------
        Path : Path to saved file
        """
        if not self.session_data:
            self.logger.warning("No training data to save")
            return None

        # Generate filename if not provided
        if filename is None:
            timestamp = self.session_start_time.strftime('%Y%m%d_%H%M%S')
            filename = f"ps2d_training_data_{timestamp}.json"

        output_path = self.output_dir / filename

        # Prepare output data
        output_data = {
            'metadata': {
                'session_start': self.session_start_time.isoformat(),
                'session_end': datetime.now().isoformat(),
                'n_samples': len(self.session_data),
                'collector_version': '1.0.0'
            },
            'samples': self.session_data
        }

        # Save to JSON
        with open(output_path, 'w') as f:
            json.dump(output_data, f, indent=2)

        self.logger.info(f"Saved {len(self.session_data)} training samples to {output_path}")
        return output_path

    def get_session_summary(self) -> Dict[str, Any]:
        """
        Get summary statistics of collected training data.

        Returns:
        --------
        dict : Summary containing sample counts, quality statistics, etc.
        """
        if not self.session_data:
            return {
                'n_samples': 0,
                'message': 'No training data collected'
            }

        # Calculate statistics
        r_squared_values = [s['fit_quality']['r_squared'] for s in self.session_data]
        peak_counts = [s['n_peaks'] for s in self.session_data]

        summary = {
            'n_samples': len(self.session_data),
            'r_squared': {
                'mean': float(np.mean(r_squared_values)),
                'median': float(np.median(r_squared_values)),
                'min': float(np.min(r_squared_values)),
                'max': float(np.max(r_squared_values))
            },
            'peaks_per_fit': {
                'mean': float(np.mean(peak_counts)),
                'median': float(np.median(peak_counts)),
                'min': int(np.min(peak_counts)),
                'max': int(np.max(peak_counts)),
                'total_peaks': int(np.sum(peak_counts))
            },
            'session_duration': str(datetime.now() - self.session_start_time)
        }

        return summary

    def clear_session(self):
        """Clear all collected training data from current session."""
        n_samples = len(self.session_data)
        self.session_data = []
        self.session_start_time = datetime.now()
        self.logger.info(f"Cleared {n_samples} training samples from session")

    def __len__(self) -> int:
        """Return number of collected training samples."""
        return len(self.session_data)

    def __repr__(self) -> str:
        """String representation of collector."""
        return f"PS2DTrainingDataCollector(n_samples={len(self.session_data)}, output_dir={self.output_dir})"
