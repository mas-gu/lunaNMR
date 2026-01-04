#!/usr/bin/env python3
# ABOUTME: CLI for automated ML training data collection from NMR spectra
# ABOUTME: Processes spectra through full pipeline (PASS1, Adaptive, PASS1-bis, PASS2) and collects training data

"""
Build Training Data CLI
========================

Automated ML training data collection by processing NMR spectra through the
complete lunaNMR pipeline, exactly as in the GUI.

Modes:
------
1. Independent Mode: Process individual spectra with their own peak lists
2. Series Mode: Process reference spectrum first, then apply parameters to series

Usage:
------
    # Independent mode
    python build_training_data.py --mode independent --input-dir ./spectra

    # Series mode
    python build_training_data.py --mode series --input-dir ./series_data

    # Custom output directory
    python build_training_data.py --mode independent --input-dir ./spectra --output-dir ./my_training_data

    # Resume from checkpoint
    python build_training_data.py --mode independent --input-dir ./spectra --resume

Expected Input Structure:
-------------------------
Independent mode:
    input_dir/
    ├── spectrum1.ft
    ├── spectrum1_peaks.txt
    ├── spectrum2.ft
    ├── spectrum2_peaks.txt
    └── ...

Series mode:
    input_dir/
    ├── reference.ft          # Reference spectrum (first alphabetically or named 'reference')
    ├── peaks.txt             # Peak list for all spectra
    ├── series_001.ft         # Time series spectra
    ├── series_002.ft
    └── ...
"""

import argparse
import json
import logging
import os
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor
from lunaNMR.utils.parameter_manager import NMRParameterManager
from lunaNMR.ml.comprehensive_collector import ComprehensiveTrainingCollector
from lunaNMR.ml.storage import get_ml_data_dir


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)


class TrainingDataBuilder:
    """
    Automated training data collection from NMR spectra.

    Processes spectra through the full lunaNMR pipeline and collects
    comprehensive training data for ML model development.
    """

    SUPPORTED_EXTENSIONS = {'.ft', '.ft2', '.ft1', '.2rr', '.2ii', '.fid', '.ucsf'}
    PEAK_LIST_PATTERNS = ['*_peaks.txt', '*_peaklist.txt', '*.txt', '*.csv']

    def __init__(
        self,
        input_dir: Path,
        output_dir: Optional[Path] = None,
        mode: str = 'independent',
        min_r2: float = 0.70,
        use_parallel: bool = True,
        resume: bool = False,
        skip_adaptive: bool = False,
    ):
        """
        Initialize the training data builder.

        Parameters
        ----------
        input_dir : Path
            Directory containing spectra and peak lists
        output_dir : Path, optional
            Custom output directory. If None, uses platform default.
        mode : str
            'independent' or 'series'
        min_r2 : float
            Minimum R² threshold for accepting samples
        use_parallel : bool
            Use parallel processing
        resume : bool
            Resume from checkpoint if available
        skip_adaptive : bool
            Skip adaptive parameter optimization (use 1.0x multiplier)
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir) if output_dir else None
        self.mode = mode
        self.min_r2 = min_r2
        self.use_parallel = use_parallel
        self.resume = resume
        self.skip_adaptive = skip_adaptive

        # Initialize collector
        self.collector = ComprehensiveTrainingCollector(
            storage_dir=self.output_dir,
            min_r2=min_r2,
            auto_save=True,
        )

        # Processing state
        self.checkpoint_file = self._get_checkpoint_path()
        self.processed_spectra: List[str] = []
        self.failed_spectra: List[Dict] = []
        self.reference_params: Optional[Dict] = None
        self.reference_statistics: Optional[Dict] = None

        # Load checkpoint if resuming
        if resume and self.checkpoint_file.exists():
            self._load_checkpoint()

    def _get_checkpoint_path(self) -> Path:
        """Get checkpoint file path."""
        base_dir = self.output_dir if self.output_dir else get_ml_data_dir()
        return base_dir / 'build_checkpoint.json'

    def _load_checkpoint(self):
        """Load processing checkpoint."""
        try:
            with open(self.checkpoint_file, 'r') as f:
                checkpoint = json.load(f)

            self.processed_spectra = checkpoint.get('processed', [])
            self.failed_spectra = checkpoint.get('failed', [])
            self.reference_params = checkpoint.get('reference_params')
            self.reference_statistics = checkpoint.get('reference_statistics')

            logger.info(f"Resumed from checkpoint: {len(self.processed_spectra)} already processed")
        except Exception as e:
            logger.warning(f"Could not load checkpoint: {e}")

    def _save_checkpoint(self):
        """Save processing checkpoint."""
        checkpoint = {
            'timestamp': datetime.now().isoformat(),
            'mode': self.mode,
            'input_dir': str(self.input_dir),
            'processed': self.processed_spectra,
            'failed': self.failed_spectra,
            'reference_params': self.reference_params,
            'reference_statistics': self.reference_statistics,
        }

        self.checkpoint_file.parent.mkdir(parents=True, exist_ok=True)
        with open(self.checkpoint_file, 'w') as f:
            json.dump(checkpoint, f, indent=2)

    def discover_spectra(self) -> List[Tuple[Path, Path]]:
        """
        Discover spectrum and peak list pairs in input directory.

        Returns
        -------
        list of tuples
            List of (spectrum_path, peak_list_path) pairs
        """
        pairs = []

        # Find all spectrum files
        spectrum_files = []
        for ext in self.SUPPORTED_EXTENSIONS:
            spectrum_files.extend(self.input_dir.glob(f'*{ext}'))
            spectrum_files.extend(self.input_dir.glob(f'**/*{ext}'))

        spectrum_files = sorted(set(spectrum_files))
        logger.info(f"Found {len(spectrum_files)} spectrum files")

        if self.mode == 'series':
            # Series mode: find reference and shared peak list
            pairs = self._discover_series_spectra(spectrum_files)
        else:
            # Independent mode: match each spectrum with its peak list
            pairs = self._discover_independent_spectra(spectrum_files)

        return pairs

    def _discover_independent_spectra(self, spectrum_files: List[Path]) -> List[Tuple[Path, Path]]:
        """Discover spectrum-peak list pairs for independent mode."""
        pairs = []

        for spectrum_path in spectrum_files:
            peak_list = self._find_peak_list_for_spectrum(spectrum_path)
            if peak_list:
                pairs.append((spectrum_path, peak_list))
            else:
                logger.warning(f"No peak list found for {spectrum_path.name}")

        return pairs

    def _discover_series_spectra(self, spectrum_files: List[Path]) -> List[Tuple[Path, Path]]:
        """Discover spectra for series mode."""
        pairs = []

        # Find shared peak list
        peak_list = self._find_series_peak_list()
        if not peak_list:
            logger.error("No shared peak list found for series mode")
            return []

        # Sort spectra - reference first
        reference = None
        series_spectra = []

        for spectrum_path in spectrum_files:
            name_lower = spectrum_path.stem.lower()
            if 'reference' in name_lower or 'ref' in name_lower:
                reference = spectrum_path
            else:
                series_spectra.append(spectrum_path)

        # If no explicit reference, use first alphabetically
        if reference is None and series_spectra:
            reference = series_spectra.pop(0)

        # Add reference first
        if reference:
            pairs.append((reference, peak_list))

        # Add series spectra
        for spectrum_path in sorted(series_spectra):
            pairs.append((spectrum_path, peak_list))

        return pairs

    def _find_peak_list_for_spectrum(self, spectrum_path: Path) -> Optional[Path]:
        """Find peak list file for a spectrum."""
        stem = spectrum_path.stem
        parent = spectrum_path.parent

        # Try various naming patterns
        candidates = [
            parent / f"{stem}_peaks.txt",
            parent / f"{stem}_peaklist.txt",
            parent / f"{stem}.txt",
            parent / f"{stem}.csv",
            parent / f"{stem}_peaks.csv",
        ]

        for candidate in candidates:
            if candidate.exists():
                return candidate

        # Try any .txt or .csv in same directory
        for pattern in ['*.txt', '*.csv']:
            matches = list(parent.glob(pattern))
            if len(matches) == 1:
                return matches[0]

        return None

    def _find_series_peak_list(self) -> Optional[Path]:
        """Find shared peak list for series mode."""
        # Look for common peak list names
        candidates = [
            self.input_dir / 'peaks.txt',
            self.input_dir / 'peaklist.txt',
            self.input_dir / 'peak_list.txt',
            self.input_dir / 'peaks.csv',
        ]

        for candidate in candidates:
            if candidate.exists():
                return candidate

        # Find any single peak list file
        txt_files = list(self.input_dir.glob('*.txt'))
        csv_files = list(self.input_dir.glob('*.csv'))

        # Filter out spectrum-related files
        peak_files = [
            f for f in txt_files + csv_files
            if 'peak' in f.stem.lower() or not any(
                f.stem.lower().endswith(ext.strip('.'))
                for ext in self.SUPPORTED_EXTENSIONS
            )
        ]

        if len(peak_files) == 1:
            return peak_files[0]

        return None

    def process_spectrum(
        self,
        spectrum_path: Path,
        peak_list_path: Path,
        is_reference: bool = False,
        series_index: int = 0,
    ) -> bool:
        """
        Process a single spectrum through the full pipeline.

        Parameters
        ----------
        spectrum_path : Path
            Path to spectrum file
        peak_list_path : Path
            Path to peak list file
        is_reference : bool
            Whether this is a series reference spectrum
        series_index : int
            Index in series (0 for reference)

        Returns
        -------
        bool
            True if processing succeeded
        """
        spectrum_name = spectrum_path.name
        logger.info(f"Processing: {spectrum_name}")

        start_time = time.time()
        timing_info = {}

        try:
            # Initialize integrator
            integrator = EnhancedVoigtIntegrator()
            param_manager = NMRParameterManager()

            # Load spectrum
            if not integrator.load_nmr_file(str(spectrum_path)):
                raise RuntimeError(f"Failed to load spectrum: {spectrum_path}")

            # Get spectrum metadata
            field_strength_mhz = self._extract_field_strength(integrator)
            nucleus_type = self._detect_nucleus_type(integrator)
            ppm_ranges = self._get_ppm_ranges(integrator)

            # Load peak list
            peak_list = self._load_peak_list(peak_list_path)
            if peak_list is None or len(peak_list) == 0:
                raise RuntimeError(f"Failed to load peak list: {peak_list_path}")

            integrator.peak_list = peak_list

            # Configure detection (match GUI)
            integrator.search_window_x = 0.05
            integrator.search_window_y = 0.5
            integrator.processing_mode = 'full_detection'
            integrator.noise_level = integrator._estimate_noise_level_advanced()

            # Run peak detection
            timing_info['detection_start'] = time.time()
            detected_peaks = integrator.process_peaks()
            timing_info['detection'] = time.time() - timing_info['detection_start']

            if detected_peaks is None:
                raise RuntimeError("Peak detection failed")

            # Convert detected peaks to DataFrame with intensity (like multi_spectrum_processor)
            # This ensures Height/Intensity columns are populated from detection
            detected_peak_data = []
            for peak in detected_peaks:
                detected_peak_data.append({
                    'Assignment': peak.get('assignment', 'Unknown'),
                    'Position_X': peak.get('ppm_x', 0),
                    'Position_Y': peak.get('ppm_y', 0),
                    'Height': peak.get('intensity', 0),
                    'Intensity': peak.get('intensity', 0)
                })
            detected_df = pd.DataFrame(detected_peak_data)

            # Update integrator.peak_list with detected DataFrame (includes intensity)
            integrator.peak_list = detected_df
            peak_list = detected_df  # Use detected DataFrame for fitting

            # Create processor
            processor = SingleSpectrumProcessor(integrator, param_manager)

            # Processing options
            processing_options = {
                'use_parallel': self.use_parallel,
                'use_global_optimization': False,
                'skip_ml_prediction': True,  # We're collecting data, not using predictions
                'skip_adaptive': self.skip_adaptive,  # Pass to processor
            }

            # For series mode, use reference statistics for non-reference spectra
            pre_learned_stats = None
            if self.mode == 'series' and not is_reference and self.reference_statistics:
                pre_learned_stats = self.reference_statistics

            # Run full processing pipeline
            timing_info['fitting_start'] = time.time()
            fitted_results, learned_statistics = processor.process_peak_list(
                peak_list,
                processing_options,
                progress_callback=self._progress_callback,
                pre_learned_statistics=pre_learned_stats,
            )
            timing_info['fitting'] = time.time() - timing_info['fitting_start']

            # For series reference, store learned statistics
            if is_reference and learned_statistics:
                self.reference_statistics = learned_statistics

                # Also get adaptive optimization results from processor
                if hasattr(processor, 'integrator') and hasattr(processor.integrator, 'enhanced_fitter'):
                    fitter = processor.integrator.enhanced_fitter
                    if hasattr(fitter, 'optimal_params'):
                        self.reference_params = fitter.optimal_params

            # Get adaptive results and pass results from parallel processor
            adaptive_results = None
            pass1_results = None
            pass1bis_results = None
            pass2_results = None
            pp_timing_info = None
            if hasattr(processor, 'integrator') and hasattr(processor.integrator, 'enhanced_fitter'):
                fitter = processor.integrator.enhanced_fitter
                # Get results from parallel processor
                if hasattr(fitter, 'parallel_processor'):
                    pp = fitter.parallel_processor
                    adaptive_results = getattr(pp, 'optimal_params', None)
                    pass1_results = getattr(pp, 'pass1_results', None)
                    pass1bis_results = getattr(pp, 'pass1bis_results', None)
                    pass2_results = getattr(pp, 'pass2_results', None)
                    pp_timing_info = getattr(pp, 'timing_info', None)
                # Fallback: check fitter directly
                elif hasattr(fitter, 'optimal_params'):
                    adaptive_results = fitter.optimal_params

            # Collect cluster info
            cluster_info = self._extract_cluster_info(processor, integrator)

            # Calculate timing - merge parallel processor timing with local timing
            total_time = time.time() - start_time
            timing_info['total'] = total_time
            # Add pass-level timing from parallel processor
            if pp_timing_info:
                timing_info['pass1'] = pp_timing_info.get('pass1', 0)
                timing_info['adaptive'] = pp_timing_info.get('adaptive', 0)
                timing_info['pass1bis'] = pp_timing_info.get('pass1bis', 0)
                timing_info['pass2'] = pp_timing_info.get('pass2', 0)

            # Get detection parameters and statistics
            detection_params = {
                'search_window_x': getattr(integrator, 'search_window_x', 0.01),
                'search_window_y': getattr(integrator, 'search_window_y', 0.05),
                'noise_threshold': getattr(integrator, 'threshold_multiplier', 3.0),
                'noise_level': getattr(integrator, 'noise_level', 0),
            }
            detection_statistics = None
            if hasattr(integrator, 'get_detection_statistics'):
                detection_statistics = integrator.get_detection_statistics()

            # Collect training data
            success = self.collector.collect_spectrum_processing(
                spectrum_data=integrator.nmr_data,
                peak_list=peak_list.to_dict('records'),
                fit_results=fitted_results,
                pass1_statistics=learned_statistics,
                adaptive_results=adaptive_results,
                pass1_results=pass1_results,
                pass1bis_results=pass1bis_results,
                pass2_results=pass2_results,
                cluster_info=cluster_info,
                timing_info=timing_info,
                spectrum_name=spectrum_name,
                nucleus_type=nucleus_type,
                field_strength_mhz=field_strength_mhz,
                ppm_ranges=ppm_ranges,
                processing_mode='parallel' if self.use_parallel else 'sequential',
                is_reference=is_reference,
                series_index=series_index,
                detection_params=detection_params,
                detection_statistics=detection_statistics,
            )

            if success:
                self.processed_spectra.append(str(spectrum_path))
                logger.info(f"  Collected: {len(fitted_results)} peaks, R²={np.mean([r.get('r_squared', 0) for r in fitted_results if r.get('success')]):.3f}")
            else:
                logger.warning(f"  Rejected (below R² threshold)")

            return success

        except Exception as e:
            logger.error(f"  Failed: {e}")
            self.failed_spectra.append({
                'spectrum': str(spectrum_path),
                'error': str(e),
                'timestamp': datetime.now().isoformat(),
            })
            return False

        finally:
            self._save_checkpoint()

    def _extract_field_strength(self, integrator) -> float:
        """Extract field strength from spectrum metadata."""
        # Try NMRPipe FDF2OBS
        if hasattr(integrator, 'dic') and integrator.dic:
            field = integrator.dic.get('FDF2OBS', None)
            if field:
                return float(field)

        # Default to 600 MHz
        return 600.0

    def _detect_nucleus_type(self, integrator) -> str:
        """Detect nucleus type from spectrum.

        Uses multiple strategies in order of reliability:
        1. Filename (user's naming convention)
        2. Spectrum metadata (axis labels)
        3. CoreIntegrator's built-in detection (spectral ranges)
        4. Default to 15N
        """
        # 1. Check filename first (most reliable for user-named files)
        if hasattr(integrator, 'spectrum_path') and integrator.spectrum_path:
            filename = str(integrator.spectrum_path).upper()
            if '13C' in filename:
                return '13C'
            if '15N' in filename:
                return '15N'

        # 2. Try to get from spectrum metadata (axis labels)
        if hasattr(integrator, 'dic') and integrator.dic:
            label = str(integrator.dic.get('FDF1LABEL', '')).upper()
            if '13C' in label:
                return '13C'
            if '15N' in label or label == 'N':
                return '15N'

        # 3. Use CoreIntegrator's built-in detection method (spectral ranges)
        if hasattr(integrator, '_detect_nucleus_type'):
            detected = integrator._detect_nucleus_type()
            if detected:
                return detected

        return '15N'  # Default assumption

    def _get_ppm_ranges(self, integrator) -> Dict:
        """Get PPM ranges from integrator."""
        ranges = {'f1': (0, 0), 'f2': (0, 0)}

        # ppm_y_axis is F1 (15N), ppm_x_axis is F2 (1H)
        if hasattr(integrator, 'ppm_y_axis') and integrator.ppm_y_axis is not None:
            ranges['f1'] = (float(min(integrator.ppm_y_axis)), float(max(integrator.ppm_y_axis)))
        if hasattr(integrator, 'ppm_x_axis') and integrator.ppm_x_axis is not None:
            ranges['f2'] = (float(min(integrator.ppm_x_axis)), float(max(integrator.ppm_x_axis)))

        return ranges

    def _load_peak_list(self, path: Path) -> Optional[pd.DataFrame]:
        """Load peak list from file."""
        try:
            # Try different separators
            for sep in [',', '\t', ' ', ';']:
                try:
                    df = pd.read_csv(path, sep=sep)
                    if len(df.columns) >= 2:
                        break
                except:
                    continue

            # Clean column names
            df.columns = [c.strip() for c in df.columns]

            # Standardize column names
            column_mapping = {
                'Position F1 (ppm)': 'Position_Y',
                'Position F2 (ppm)': 'Position_X',
                'F1': 'Position_Y',
                'F2': 'Position_X',
                '15N': 'Position_Y',
                '1H': 'Position_X',
                'ppm_f1': 'Position_Y',
                'ppm_f2': 'Position_X',
            }

            for old, new in column_mapping.items():
                if old in df.columns and new not in df.columns:
                    df = df.rename(columns={old: new})

            # Ensure required columns exist
            if 'Position_X' not in df.columns or 'Position_Y' not in df.columns:
                # Try to use first two numeric columns
                numeric_cols = df.select_dtypes(include=[np.number]).columns
                if len(numeric_cols) >= 2:
                    df = df.rename(columns={
                        numeric_cols[0]: 'Position_X',
                        numeric_cols[1]: 'Position_Y',
                    })

            return df

        except Exception as e:
            logger.error(f"Failed to load peak list: {e}")
            return None

    def _extract_cluster_info(self, processor, integrator) -> Dict:
        """Extract clustering information from processor."""
        cluster_info = {
            'isolated_clusters': [],
            'multi_peak_clusters': [],
        }

        # Try to get cluster info from processor
        if hasattr(processor, '_last_cluster_info'):
            cluster_info = processor._last_cluster_info

        return cluster_info

    def _progress_callback(self, percent: int, task: str, detail: str = ''):
        """Progress callback for processing."""
        # Minimal logging during batch processing
        pass

    def run(self) -> Dict:
        """
        Run the training data collection.

        Returns
        -------
        dict
            Summary of processing results
        """
        logger.info("=" * 60)
        logger.info("Build Training Data CLI")
        logger.info("=" * 60)
        logger.info(f"Mode: {self.mode}")
        logger.info(f"Input: {self.input_dir}")
        logger.info(f"Output: {self.output_dir or 'default'}")
        logger.info(f"Min R²: {self.min_r2}")
        logger.info(f"Skip adaptive: {self.skip_adaptive}")
        logger.info("=" * 60)

        # Discover spectra
        pairs = self.discover_spectra()
        if not pairs:
            logger.error("No spectrum-peak list pairs found!")
            return {'success': False, 'error': 'No input data found'}

        logger.info(f"Found {len(pairs)} spectrum-peak list pairs")

        # Filter already processed
        if self.resume:
            pairs = [
                (s, p) for s, p in pairs
                if str(s) not in self.processed_spectra
            ]
            logger.info(f"After resume filter: {len(pairs)} remaining")

        # Process spectra
        start_time = time.time()

        for i, (spectrum_path, peak_list_path) in enumerate(pairs):
            is_reference = (self.mode == 'series' and i == 0)
            series_index = i

            logger.info(f"\n[{i+1}/{len(pairs)}] {'REFERENCE: ' if is_reference else ''}{spectrum_path.name}")

            self.process_spectrum(
                spectrum_path,
                peak_list_path,
                is_reference=is_reference,
                series_index=series_index,
            )

        elapsed = time.time() - start_time

        # Summary
        logger.info("\n" + "=" * 60)
        logger.info("SUMMARY")
        logger.info("=" * 60)
        logger.info(f"Total time: {elapsed:.1f}s")
        logger.info(f"Processed: {len(self.processed_spectra)}")
        logger.info(f"Failed: {len(self.failed_spectra)}")
        logger.info(f"Collected samples: {self.collector.n_collected}")
        logger.info(f"Rejected samples: {self.collector.n_rejected}")

        if self.failed_spectra:
            logger.info("\nFailed spectra:")
            for f in self.failed_spectra:
                logger.info(f"  - {f['spectrum']}: {f['error']}")

        # Clean up checkpoint on success
        if not self.failed_spectra and self.checkpoint_file.exists():
            self.checkpoint_file.unlink()
            logger.info("\nCheckpoint cleaned up")

        return {
            'success': True,
            'processed': len(self.processed_spectra),
            'failed': len(self.failed_spectra),
            'collected': self.collector.n_collected,
            'rejected': self.collector.n_rejected,
            'elapsed': elapsed,
        }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Build ML training data from NMR spectra',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        '--input-dir', '-i',
        type=Path,
        required=True,
        help='Directory containing spectra and peak lists'
    )

    parser.add_argument(
        '--output-dir', '-o',
        type=Path,
        default=None,
        help='Output directory (default: platform-specific location)'
    )

    parser.add_argument(
        '--mode', '-m',
        choices=['independent', 'series'],
        default='independent',
        help='Processing mode (default: independent)'
    )

    parser.add_argument(
        '--min-r2',
        type=float,
        default=0.70,
        help='Minimum R² threshold for accepting samples (default: 0.70)'
    )

    parser.add_argument(
        '--no-parallel',
        action='store_true',
        help='Disable parallel processing'
    )

    parser.add_argument(
        '--resume',
        action='store_true',
        help='Resume from checkpoint if available'
    )

    parser.add_argument(
        '--skip-adaptive',
        action='store_true',
        help='Skip adaptive parameter optimization (use 1.0x multiplier)'
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Validate input directory
    if not args.input_dir.exists():
        logger.error(f"Input directory does not exist: {args.input_dir}")
        sys.exit(1)

    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Suppress lunaNMR output for cleaner logs
    import lunaNMR.utils.output_manager as om
    om.log_info = lambda *a, **k: None
    om.log_progress = lambda *a, **k: None
    om.log_warning = lambda *a, **k: None

    # Create builder and run
    builder = TrainingDataBuilder(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        mode=args.mode,
        min_r2=args.min_r2,
        use_parallel=not args.no_parallel,
        resume=args.resume,
        skip_adaptive=args.skip_adaptive,
    )

    result = builder.run()

    if result['success']:
        logger.info("\nTraining data collection complete!")
        sys.exit(0)
    else:
        logger.error(f"\nFailed: {result.get('error', 'Unknown error')}")
        sys.exit(1)


if __name__ == '__main__':
    main()
