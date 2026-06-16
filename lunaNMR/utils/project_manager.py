# ABOUTME: Manages save/load of lunaNMR project bundles (.lunaNMR directories)
# ABOUTME: Serializes GUI state, peak lists, fit results, and file references

import json
import logging
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any, List, Tuple
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _safe_bundle_name(name: str) -> str:
    """A filesystem-safe folder name for a saved analysis: keep word chars,
    dot and dash; collapse anything else to '_'. Never empty."""
    import re
    safe = re.sub(r'[^\w.\-]+', '_', str(name)).strip('._')
    return safe or 'analysis'


class ProjectManager:
    """Manages save/load of lunaNMR project bundles.

    A project bundle is a directory ending with .lunaNMR containing:
    - project.json: Master manifest with metadata and file references
    - gui_state.json: All GUI parameters
    - peak_list.csv: Peak list DataFrame
    - fit_results/: Fit results with numpy arrays
    - series_results/: Series integration results (if applicable)
    - dynamixs/: DynamiXs analysis state (if applicable)
    """

    SCHEMA_VERSION = "1.1"
    PROJECT_EXTENSION = ".lunaNMR"

    # List of GUI parameters to serialize
    # From main_window.py lines 290-462
    GUI_PARAMS = [
        # File paths (stored as references)
        'current_nmr_file', 'current_peak_file', 'current_nmr_folder',
        'current_peak_folder', 'series_output_folder',

        # Workflow configuration
        'workflow_mode', 'sn_threshold', 'expected_peak_count',

        # Display options
        'show_detected', 'show_assigned', 'show_fitted_curves', 'show_ellipses',

        # Contour display
        'contour_levels', 'contour_min', 'contour_increment',

        # Zoom parameters
        'zoom_x_range', 'zoom_y_range',
        # Note: saved_xlim, saved_ylim are tuples - handle separately

        # Processing parameters
        'use_voigt_fitting', 'use_parallel_processing', 'use_adaptive_optimization',
        'use_ps2d_multi_peak', 'fix_linewidths', 'fix_positions',

        # Peak detection parameters
        'noise_threshold', 'search_window_x', 'search_window_y',
        'detection_square_size', 'detection_rectangle_y',
        'centroid_window_x_ppm', 'centroid_window_y_ppm',
        'use_centroid_refinement', 'centroid_noise_multiplier',
        'auto_add_dummy_peaks',

        # Peak detection for series
        'height_threshold', 'distance_factor', 'prominence_threshold',
        'smoothing_sigma', 'max_peaks_fit', 'max_optimization_iterations',

        # Fitting quality parameters
        'min_r_squared', 'max_iterations', 'fitting_window_x', 'fitting_window_y',
        'use_global_optimization',

        # PS2D algorithm configuration
        'nucleus_type', 'ps2d_radF1', 'ps2d_radF2', 'ps2d_max_iterations',
        'ps2d_max_cluster_size', 'ps2d_overlap_x', 'ps2d_overlap_y',

        # Custom linewidth parameters
        'use_custom_linewidths', 'lw_lorentz_1h', 'lw_gauss_1h',
        'lw_lorentz_15n', 'lw_gauss_15n',

        # Series integration parameters
        'use_ps2d_linewidth_reuse', 'enable_cascade_drift_limit',

        # Shift peak list
        'shift_1h_value', 'shift_15n_value',

        # 3D Voigt analysis
        'show_exp_3d', 'show_fit_3d', 'show_individual_3d', 'show_peak_labels_3d',
        'show_resid_3d', 'limit_peak_display_3d', 'residual_mode_3d',
        'color_scheme_3d', 'intensity_scale_3d',
    ]

    def __init__(self, main_window):
        """Initialize ProjectManager.

        Args:
            main_window: Reference to LunaNMRMainWindow for state access
        """
        self.main_window = main_window

    # =========================================================================
    # SAVE OPERATIONS
    # =========================================================================

    def save_project(self, project_path: Path,
                     embed_region_data: bool = False) -> Tuple[bool, Dict[str, Any]]:
        """Save complete project state to bundle directory.

        Args:
            project_path: Path ending with .lunaNMR
            embed_region_data: When True, store each fit region's raw intensity
                grid in the bundle so fit-detail plots open without the original
                spectrum. When False (default), only the region bounds are stored
                and intensity is re-sliced from the spectrum on load.

        Returns:
            Tuple of (success: bool, summary: Dict[str, Any])
        """
        project_path = Path(project_path)
        summary = {
            'project_name': project_path.stem,
            'saved_items': [],
            'dynamixs': {}
        }

        try:
            # Import QApplication for processEvents
            try:
                from PySide6.QtWidgets import QApplication
                process_events = lambda: QApplication.processEvents()
            except ImportError:
                process_events = lambda: None  # No-op if not in GUI context

            # Create bundle directory
            project_path.mkdir(parents=True, exist_ok=True)

            # Save components and track what was saved
            self._save_manifest(project_path)
            summary['saved_items'].append('Project manifest')
            process_events()

            self._save_gui_state(project_path)
            summary['saved_items'].append('GUI state')
            process_events()

            if self._save_peak_list(project_path):
                peak_count = len(getattr(self.main_window, 'peaks', []))
                if peak_count > 0:
                    summary['saved_items'].append(f'Peak list ({peak_count} peaks)')
            process_events()

            if self._save_fit_results(project_path, embed_region_data):
                summary['saved_items'].append('Fit results')
            process_events()

            if self._save_series_results(project_path):
                summary['saved_items'].append('Series results')
            process_events()

            # Save DynamiXs state and get detailed summary
            dynamixs_summary = self._save_dynamixs_state_with_summary(project_path)
            if dynamixs_summary:
                summary['dynamixs'] = dynamixs_summary
            process_events()

            if self._save_kd_state(project_path):
                summary['saved_items'].append('Kd titration')
            process_events()

            logger.info(f"Project saved to {project_path}")
            return True, summary

        except Exception as e:
            logger.error(f"Failed to save project: {e}")
            return False, summary

    def _save_manifest(self, bundle_path: Path) -> bool:
        """Save project.json with metadata and file references."""
        manifest = {
            "schema_version": self.SCHEMA_VERSION,
            "project_name": bundle_path.stem,
            "created_at": datetime.now().isoformat(),
            "modified_at": datetime.now().isoformat(),
            "lunaNMR_version": "1.0",

            "file_references": {
                "nmr_spectrum": self._get_safe_attr('current_nmr_file'),
                "peak_list_source": self._get_safe_attr('current_peak_file'),
                "output_folder": self._get_safe_attr('series_output_folder'),
            },

            "workflow_state": {
                "workflow_mode": getattr(self.main_window, 'workflow_mode', 'peak_list'),
                "has_peak_list": self._has_peak_list(),
                "has_fit_results": self._has_fit_results(),
                "has_series_results": self._has_series_results(),
            }
        }

        manifest_path = bundle_path / "project.json"
        with open(manifest_path, 'w') as f:
            json.dump(manifest, f, indent=2, default=str)

        return True

    def _save_gui_state(self, bundle_path: Path) -> bool:
        """Save all GUI parameters to gui_state.json."""
        params = self._collect_gui_params()

        gui_state_path = bundle_path / "gui_state.json"
        with open(gui_state_path, 'w') as f:
            json.dump(params, f, indent=2, default=self._json_serializer)

        return True

    def _save_peak_list(self, bundle_path: Path) -> bool:
        """Save peak list DataFrame to CSV with fitted positions and R² values.

        When fit results exist, uses fitted positions (center_x, center_y) instead
        of original reference positions for better accuracy.
        """
        fit_results = getattr(self.main_window, 'last_fitting_results', None)

        # If we have fit results, build peak list from fitted data (more accurate)
        if fit_results:
            rows = []
            for fit in fit_results:
                assignment = fit.get('assignment')
                if assignment is None:
                    continue

                row = {
                    'Assignment': assignment,
                    'Position_X': fit.get('center_x', fit.get('peak_x', 0)),
                    'Position_Y': fit.get('center_y', fit.get('peak_y', 0)),
                    'Height': fit.get('height', fit.get('amplitude', '')),
                    'Intensity': fit.get('volume', fit.get('intensity', '')),
                    'R_Squared': fit.get('r_squared'),
                    'Fitting_Quality': fit.get('fitting_quality', ''),
                    'Method': fit.get('method', ''),
                }
                rows.append(row)

            if rows:
                peak_list = pd.DataFrame(rows)
                peak_list_path = bundle_path / "peak_list.csv"
                peak_list.to_csv(peak_list_path, index=False)
                logger.info(f"Saved {len(rows)} fitted peaks to peak_list.csv")
                return True

        # Fallback: save reference peak list if no fit results
        if not self._has_peak_list():
            return True  # Nothing to save

        peak_list = self.main_window.integrator.peak_list.copy()
        peak_list_path = bundle_path / "peak_list.csv"
        peak_list.to_csv(peak_list_path, index=False)
        logger.info(f"Saved {len(peak_list)} reference peaks to peak_list.csv")

        return True

    def _save_fit_results(self, bundle_path: Path,
                          embed_region_data: bool = False) -> bool:
        """Save fit results, separating numpy arrays from scalars.

        Fit results from "Fit Spectrum" (PS2D) contain:
        - Scalar params: r_squared, amplitude, center_x/y, sigma_x/y, gamma_x/y, etc.
        - Per-peak params: all_peaks (list of dicts)
        - Numpy arrays: region_2d (dict with arrays), fitted_2d_surface, individual_surfaces

        The surfaces and the region coordinate grids are reconstructable from the
        per-peak params and the spectrum, so they are never stored. Each region's
        bounds (indices into the full spectrum) are stored so the raw intensity can
        be re-sliced on load. When embed_region_data is True, the raw intensity is
        also stored so the project opens without the original spectrum.
        """
        if not self._has_fit_results():
            return True  # Nothing to save

        # Create fit_results directory structure
        fit_dir = bundle_path / "fit_results"
        fit_dir.mkdir(exist_ok=True)
        arrays_dir = fit_dir / "arrays"
        arrays_dir.mkdir(exist_ok=True)

        # Get fit results from main_window
        fit_results = getattr(self.main_window, 'last_fitting_results', None)
        if fit_results is None or len(fit_results) == 0:
            return True

        integrator = getattr(self.main_window, 'integrator', None)

        # Serialize each fit result
        serialized_results = []
        for idx, fit_result in enumerate(fit_results):
            if not isinstance(fit_result, dict):
                continue

            peak_id = f"peak_{idx:03d}"
            serialized = self._serialize_fit_result(
                fit_result, peak_id, arrays_dir, embed_region_data, integrator)
            serialized_results.append(serialized)

        # Save scalars to JSON (compact format for faster saving)
        fits_path = fit_dir / "fits.json"
        with open(fits_path, 'w') as f:
            json.dump(serialized_results, f, default=self._json_serializer)

        logger.info(f"Saved {len(serialized_results)} fit results to {fit_dir}")
        return True

    def _serialize_fit_result(self, fit_result: dict, peak_id: str, arrays_dir: Path,
                              embed_region_data: bool, integrator) -> dict:
        """Serialize single fit result, dropping reconstructable arrays.

        fitted_2d_surface and individual_surfaces are always dropped (rebuilt from
        all_peaks on load). region_2d is replaced by region_bounds, plus an
        intensity .npz when embed_region_data is True.
        """
        serialized = {}

        DROP_KEYS = ('fitted_2d_surface', 'individual_surfaces')

        for key, value in fit_result.items():
            if key in DROP_KEYS:
                # Reconstructed on load from all_peaks + region grid
                continue

            elif key == 'region_2d' and isinstance(value, dict):
                bounds = self._compute_region_bounds(value, integrator)
                if bounds is not None:
                    serialized['region_bounds'] = bounds
                if embed_region_data:
                    region_file = arrays_dir / f"{peak_id}_region_2d.npz"
                    self._save_region_intensity(value, region_file)
                    serialized['region_2d'] = f"arrays/{peak_id}_region_2d.npz"

            elif key == 'all_peaks' and isinstance(value, list):
                # List of peak parameter dicts - serialize as-is (no arrays inside)
                serialized[key] = value

            elif isinstance(value, np.ndarray):
                # Other numpy arrays - convert to list for JSON
                serialized[key] = value.tolist()

            elif isinstance(value, (np.integer, np.floating)):
                # Convert numpy scalars
                serialized[key] = float(value) if isinstance(value, np.floating) else int(value)

            elif hasattr(value, '_mock_name'):
                # Skip mock objects (testing)
                continue

            else:
                # Regular values (str, int, float, bool, dict, list)
                serialized[key] = value

        return serialized

    def _compute_region_bounds(self, region_2d: dict, integrator):
        """Compute [y0, y1, x0, x1] index bounds of a region in the full spectrum.

        Returns None when the spectrum axes or region axes are unavailable.
        """
        if integrator is None:
            return None
        ppm_y = getattr(integrator, 'ppm_y_axis', None)
        ppm_x = getattr(integrator, 'ppm_x_axis', None)
        f1_ppm = region_2d.get('f1_ppm')
        f2_ppm = region_2d.get('f2_ppm')
        if ppm_y is None or ppm_x is None or f1_ppm is None or f2_ppm is None:
            return None

        f1_ppm = np.asarray(f1_ppm)
        f2_ppm = np.asarray(f2_ppm)
        if f1_ppm.size == 0 or f2_ppm.size == 0:
            return None

        y0 = int(np.argmin(np.abs(np.asarray(ppm_y) - f1_ppm[0])))
        x0 = int(np.argmin(np.abs(np.asarray(ppm_x) - f2_ppm[0])))
        return [y0, y0 + int(f1_ppm.size), x0, x0 + int(f2_ppm.size)]

    def _save_region_intensity(self, region_2d: dict, filepath: Path):
        """Save only the raw intensity grid and 1D axes of a region."""
        arrays_to_save = {}
        for key in ('intensity', 'f1_ppm', 'f2_ppm'):
            value = region_2d.get(key)
            if value is not None:
                arrays_to_save[key] = np.asarray(value)
        np.savez_compressed(filepath, **arrays_to_save)

    # =========================================================================
    # LOAD OPERATIONS
    # =========================================================================

    def load_project(self, project_path: Path) -> Tuple[bool, List[str], Dict[str, Any]]:
        """Load project state from bundle directory.

        Args:
            project_path: Path to .lunaNMR bundle

        Returns:
            Tuple of (success: bool, missing_files: List[str], summary: Dict[str, Any])
        """
        project_path = Path(project_path)
        missing_files = []
        summary = {
            'project_name': project_path.stem,
            'loaded_items': [],
            'dynamixs': {}
        }

        try:
            # Load and validate manifest
            manifest = self._load_manifest(project_path)
            if manifest is None:
                return False, ["project.json not found or invalid"], summary

            summary['loaded_items'].append('Project manifest')

            # Check file references
            missing_files = self._validate_file_references(manifest)

            # Load components and track what was loaded
            self._load_gui_state(project_path)
            summary['loaded_items'].append('GUI state')

            if self._load_peak_list(project_path):
                peak_count = len(getattr(self.main_window, 'peaks', []))
                if peak_count > 0:
                    summary['loaded_items'].append(f'Peak list ({peak_count} peaks)')

            if self._load_fit_results(project_path):
                summary['loaded_items'].append('Fit results')

            if self._load_series_results(project_path):
                summary['loaded_items'].append('Series results')

            # Load DynamiXs state and get summary
            dynamixs_summary = self._load_dynamixs_state_with_summary(project_path)
            if dynamixs_summary:
                summary['dynamixs'] = dynamixs_summary

            if self._load_kd_state(project_path):
                summary['loaded_items'].append('Kd titration')

            logger.info(f"Project loaded from {project_path}")
            return True, missing_files, summary

        except Exception as e:
            logger.error(f"Failed to load project: {e}")
            return False, [str(e)], summary

    def _load_manifest(self, bundle_path: Path) -> Optional[Dict]:
        """Load and validate project.json."""
        manifest_path = bundle_path / "project.json"

        if not manifest_path.exists():
            return None

        with open(manifest_path, 'r') as f:
            manifest = json.load(f)

        # Check schema version
        schema_version = manifest.get('schema_version', '1.0')
        if schema_version != self.SCHEMA_VERSION:
            logger.warning(f"Project uses schema version {schema_version}, "
                          f"current is {self.SCHEMA_VERSION}")

        return manifest

    def _validate_file_references(self, manifest: Dict) -> List[str]:
        """Check which referenced files exist, return list of missing paths."""
        missing = []
        file_refs = manifest.get('file_references', {})

        for ref_name, ref_path in file_refs.items():
            if ref_path and not Path(ref_path).exists():
                missing.append(ref_path)

        return missing

    def get_missing_files_structured(self, project_path: Path) -> List[Dict]:
        """Get structured list of missing files with type info.

        Args:
            project_path: Path to .lunaNMR bundle

        Returns:
            List of dicts with 'path' and 'type' keys
        """
        project_path = Path(project_path)
        missing_files = []

        # Load manifest
        manifest = self._load_manifest(project_path)
        if manifest is None:
            return missing_files

        file_refs = manifest.get('file_references', {})

        # Map reference names to file types for the dialog
        type_mapping = {
            'nmr_spectrum': 'nmr_spectrum',
            'peak_list_source': 'peak_list',
            'output_folder': 'output_folder',
        }

        for ref_name, ref_path in file_refs.items():
            if ref_path and not Path(ref_path).exists():
                file_type = type_mapping.get(ref_name, ref_name)
                missing_files.append({
                    'path': ref_path,
                    'type': file_type,
                    'ref_name': ref_name,
                })

        return missing_files

    def apply_path_remapping(self, remapped_paths: Dict[str, str]):
        """Apply path remapping to GUI state after loading.

        Args:
            remapped_paths: Dict mapping old paths to new paths
        """
        # Remap file paths in main_window
        path_attrs = ['current_nmr_file', 'current_peak_file', 'series_output_folder']

        for attr in path_attrs:
            current_value = getattr(self.main_window, attr, None)
            if current_value and current_value in remapped_paths:
                setattr(self.main_window, attr, remapped_paths[current_value])
                logger.debug(f"Remapped {attr}: {current_value} -> {remapped_paths[current_value]}")

    def _load_gui_state(self, bundle_path: Path) -> bool:
        """Load GUI parameters from gui_state.json."""
        gui_state_path = bundle_path / "gui_state.json"

        if not gui_state_path.exists():
            logger.warning("gui_state.json not found, using defaults")
            return False

        with open(gui_state_path, 'r') as f:
            params = json.load(f)

        self._apply_gui_params(params)
        return True

    def _load_peak_list(self, bundle_path: Path) -> bool:
        """Load peak list DataFrame from CSV."""
        peak_list_path = bundle_path / "peak_list.csv"

        if not peak_list_path.exists():
            logger.debug("peak_list.csv not found")
            return False

        df = pd.read_csv(peak_list_path)
        self.main_window.integrator.peak_list = df

        return True

    def _load_fit_results(self, bundle_path: Path) -> bool:
        """Load fit results including numpy arrays.

        Reconstructs the fit results list with numpy arrays loaded from .npz files.
        """
        fit_dir = bundle_path / "fit_results"
        fits_path = fit_dir / "fits.json"

        if not fits_path.exists():
            logger.debug("fits.json not found")
            return False

        with open(fits_path, 'r') as f:
            serialized_results = json.load(f)

        # Deserialize each fit result
        fit_results = []
        for serialized in serialized_results:
            fit_result = self._deserialize_fit_result(serialized, fit_dir)
            fit_results.append(fit_result)

        # Store in main_window
        self.main_window.last_fitting_results = fit_results

        logger.info(f"Loaded {len(fit_results)} fit results from {fit_dir}")
        return True

    def _deserialize_fit_result(self, serialized: dict, fit_dir: Path) -> dict:
        """Reconstruct fit result, loading embedded region intensity if present."""
        result = {}

        for key, value in serialized.items():
            if key == 'region_2d' and isinstance(value, str) and value.endswith('.npz'):
                # Load embedded region intensity from npz
                npz_path = fit_dir / value
                if npz_path.exists():
                    result[key] = self._load_region_2d(npz_path)
                else:
                    logger.warning(f"Missing array file: {npz_path}")
                    result[key] = None

            else:
                # Regular value - use as-is
                result[key] = value

        return result

    def _load_region_2d(self, npz_path: Path) -> dict:
        """Load region_2d dict from npz file."""
        region_2d = {}
        with np.load(npz_path) as data:
            for key in data.files:
                arr = data[key]
                # Convert 0-d arrays back to scalars if needed
                if arr.ndim == 0:
                    region_2d[key] = arr.item()
                else:
                    region_2d[key] = arr
        return region_2d

    def reconstruct_fit_arrays(self) -> None:
        """Rebuild region grids and fit surfaces for loaded fit results.

        Must run after the spectrum is loaded so off-mode regions can be
        re-sliced. Results whose region cannot be resolved (no embedded data and
        no spectrum) are left without surfaces; downstream plotting degrades
        gracefully.
        """
        results = getattr(self.main_window, 'last_fitting_results', None)
        integrator = getattr(self.main_window, 'integrator', None)
        if not results or integrator is None:
            return

        for result in results:
            if not isinstance(result, dict) or 'all_peaks' not in result:
                continue
            try:
                region_2d = self._resolve_region_2d(result, integrator)
                if region_2d is None:
                    continue
                surface, individual, baseline = integrator._reconstruct_2d_surface(
                    region_2d, result['all_peaks'])
            except Exception as e:
                logger.warning(f"Could not reconstruct fit surface for "
                               f"{result.get('assignment', '?')}: {e}")
                continue
            result['region_2d'] = region_2d
            result['fitted_2d_surface'] = surface
            result['individual_surfaces'] = individual
            result.setdefault('baseline', baseline)

    def _resolve_region_2d(self, result: dict, integrator):
        """Return a region_2d dict with intensity and grids, or None."""
        region = result.get('region_2d')
        if isinstance(region, dict) and region.get('intensity') is not None \
                and region.get('f1_ppm') is not None and region.get('f2_ppm') is not None:
            return self._build_region_2d(region['f1_ppm'], region['f2_ppm'],
                                         region['intensity'])

        bounds = result.get('region_bounds')
        if bounds is None or len(bounds) != 4:
            return None
        nmr_data = getattr(integrator, 'nmr_data', None)
        ppm_y = getattr(integrator, 'ppm_y_axis', None)
        ppm_x = getattr(integrator, 'ppm_x_axis', None)
        if nmr_data is None or ppm_y is None or ppm_x is None:
            return None

        y0, y1, x0, x1 = bounds
        intensity = nmr_data[y0:y1, x0:x1]
        if intensity.size == 0:
            return None
        return self._build_region_2d(np.asarray(ppm_y)[y0:y1],
                                     np.asarray(ppm_x)[x0:x1], intensity)

    def _build_region_2d(self, f1_ppm, f2_ppm, intensity) -> dict:
        """Assemble a region_2d dict, rebuilding the coordinate grids."""
        f1_ppm = np.asarray(f1_ppm)
        f2_ppm = np.asarray(f2_ppm)
        f2_grid, f1_grid = np.meshgrid(f2_ppm, f1_ppm, indexing='xy')
        return {
            'f1_ppm': f1_ppm,
            'f2_ppm': f2_ppm,
            'f1_grid': f1_grid,
            'f2_grid': f2_grid,
            'intensity': np.asarray(intensity),
        }

    # =========================================================================
    # UTILITY METHODS
    # =========================================================================

    def _collect_gui_params(self) -> Dict[str, Any]:
        """Extract all GUI parameters from main_window as dict."""
        params = {}

        for param_name in self.GUI_PARAMS:
            if hasattr(self.main_window, param_name):
                value = getattr(self.main_window, param_name)
                # Skip mock objects (testing) - only serialize real values
                if hasattr(value, '_mock_name'):
                    continue
                params[param_name] = value

        # Handle tuple params separately (saved_xlim, saved_ylim)
        for tuple_param in ['saved_xlim', 'saved_ylim']:
            if hasattr(self.main_window, tuple_param):
                value = getattr(self.main_window, tuple_param)
                # Skip mock objects
                if hasattr(value, '_mock_name'):
                    continue
                if value is not None:
                    params[tuple_param] = list(value)  # Convert tuple to list for JSON

        return params

    def _apply_gui_params(self, params: Dict[str, Any]) -> None:
        """Apply dict of parameters to main_window GUI."""
        for param_name, value in params.items():
            if param_name in self.GUI_PARAMS or param_name in ['saved_xlim', 'saved_ylim']:
                # Convert lists back to tuples for xlim/ylim
                if param_name in ['saved_xlim', 'saved_ylim'] and value is not None:
                    value = tuple(value)

                setattr(self.main_window, param_name, value)

    def _json_serializer(self, obj):
        """Custom JSON serializer for non-standard types."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        if isinstance(obj, Path):
            return str(obj)
        # Handle mock objects (for testing) - serialize as None
        if hasattr(obj, '_mock_name'):
            return None
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

    def _get_safe_attr(self, attr_name: str, default=None):
        """Get attribute from main_window, returning default for mock objects."""
        value = getattr(self.main_window, attr_name, default)
        # Return None for mock objects (testing)
        if hasattr(value, '_mock_name'):
            return default
        return value

    def _has_peak_list(self) -> bool:
        """Check if main_window has a peak list."""
        integrator = getattr(self.main_window, 'integrator', None)
        if integrator is None:
            return False
        peak_list = getattr(integrator, 'peak_list', None)
        return peak_list is not None and len(peak_list) > 0

    def _has_fit_results(self) -> bool:
        """Check if main_window has fit results."""
        return (getattr(self.main_window, 'current_voigt_result', None) is not None or
                getattr(self.main_window, 'last_fitting_results', None) is not None)

    def _has_series_results(self) -> bool:
        """Check if main_window has series integration results."""
        # Check for multiple saved series
        saved_series = getattr(self.main_window, 'saved_series', None)
        if saved_series and len(saved_series) > 0:
            return True
        # Fallback to single batch_results for backward compatibility
        return getattr(self.main_window, 'batch_results', None) is not None

    # =========================================================================
    # SERIES RESULTS SAVE/LOAD
    # =========================================================================

    def _save_series_results(self, bundle_path: Path) -> bool:
        """Save series integration results (multiple series) to JSON.

        Supports multiple named series in saved_series dict.
        Falls back to single batch_results for backward compatibility.

        Each series is saved in its own subdirectory:
        series_results/{series_name}/batch_results.json
        """
        # Import QApplication for processEvents
        try:
            from PySide6.QtWidgets import QApplication
            process_events = lambda: QApplication.processEvents()
        except ImportError:
            process_events = lambda: None  # No-op if not in GUI context

        # Create series_results directory
        series_dir = bundle_path / "series_results"
        series_dir.mkdir(exist_ok=True)

        # Get all saved series
        saved_series = getattr(self.main_window, 'saved_series', {}) or {}

        # For backward compatibility, if no saved_series but batch_results exists,
        # save it as "default" series
        if not saved_series:
            batch_results = getattr(self.main_window, 'batch_results', None)
            if batch_results is not None:
                saved_series = {'default': batch_results}

        if not saved_series:
            return True  # Nothing to save

        # Save manifest listing all series
        manifest = {
            'series_count': len(saved_series),
            'series_list': []
        }

        # Save each series
        for series_name, batch_results in saved_series.items():
            process_events()  # Keep GUI responsive
            # Create subdirectory for this series
            series_subdir = series_dir / series_name
            series_subdir.mkdir(exist_ok=True)

            # Build serializable data structure
            data = {
                'results': {},
                'metadata': self._serialize_metadata(batch_results.metadata if hasattr(batch_results, 'metadata') else {}),
                'statistics': batch_results.statistics if hasattr(batch_results, 'statistics') else {},
                'errors': batch_results.errors if hasattr(batch_results, 'errors') else [],
            }

            # Serialize each spectrum result
            if hasattr(batch_results, 'results'):
                for name, result in batch_results.results.items():
                    data['results'][name] = self._serialize_spectrum_result(result)

            # Write to JSON (compact format for faster saving)
            batch_file = series_subdir / "batch_results.json"
            with open(batch_file, 'w') as f:
                json.dump(data, f, default=self._json_serializer)

            # Add to manifest
            spectrum_count = len(batch_results.results) if hasattr(batch_results, 'results') else 0
            manifest['series_list'].append({
                'name': series_name,
                'spectrum_count': spectrum_count,
            })

            logger.info(f"Saved series '{series_name}': {spectrum_count} spectra")

        # Write manifest (compact format for faster saving)
        manifest_file = series_dir / "series_manifest.json"
        with open(manifest_file, 'w') as f:
            json.dump(manifest, f)

        logger.info(f"Saved {len(saved_series)} series integrations")
        return True

    def _serialize_metadata(self, metadata: dict) -> dict:
        """Serialize metadata, converting datetime objects to ISO strings."""
        serialized = {}
        for key, value in metadata.items():
            if hasattr(value, 'isoformat'):  # datetime object
                serialized[key] = value.isoformat()
            else:
                serialized[key] = value
        return serialized

    def _serialize_spectrum_result(self, result: dict) -> dict:
        """Serialize a single spectrum result.

        Handles numpy arrays and scalars in integration_results.
        """
        serialized = {}
        for key, value in result.items():
            if isinstance(value, np.ndarray):
                serialized[key] = value.tolist()
            elif isinstance(value, list):
                # Handle integration_results list
                serialized[key] = [
                    self._serialize_integration_result(r) if isinstance(r, dict) else r
                    for r in value
                ]
            elif isinstance(value, (np.integer, np.floating)):
                serialized[key] = float(value) if isinstance(value, np.floating) else int(value)
            else:
                serialized[key] = value
        return serialized

    def _serialize_integration_result(self, result: dict) -> dict:
        """Serialize a single integration result dict."""
        serialized = {}
        for key, value in result.items():
            if isinstance(value, np.ndarray):
                serialized[key] = value.tolist()
            elif isinstance(value, (np.integer, np.floating)):
                serialized[key] = float(value) if isinstance(value, np.floating) else int(value)
            else:
                serialized[key] = value
        return serialized

    def _load_series_results(self, bundle_path: Path) -> bool:
        """Load series integration results (multiple series) and reconstruct BatchResults.

        Supports multiple named series from series_manifest.json.
        Falls back to single batch_results.json for backward compatibility.
        """
        series_dir = bundle_path / "series_results"

        if not series_dir.exists():
            return True  # No series results to load

        from lunaNMR.processors.multi_spectrum_processor import BatchResults

        # Initialize saved_series dict
        self.main_window.saved_series = {}

        # Check for manifest (new format with multiple series)
        manifest_file = series_dir / "series_manifest.json"
        if manifest_file.exists():
            with open(manifest_file, 'r') as f:
                manifest = json.load(f)

            # Load each series from its subdirectory
            for series_info in manifest.get('series_list', []):
                series_name = series_info.get('name', 'unnamed')
                series_subdir = series_dir / series_name
                batch_file = series_subdir / "batch_results.json"

                if not batch_file.exists():
                    logger.warning(f"Series '{series_name}' data not found, skipping")
                    continue

                with open(batch_file, 'r') as f:
                    data = json.load(f)

                # Reconstruct BatchResults object
                batch_results = BatchResults()
                batch_results.metadata = data.get('metadata', {})
                batch_results.metadata['series_name'] = series_name
                batch_results.statistics = data.get('statistics', {})
                batch_results.errors = data.get('errors', [])

                # Restore results dict
                for name, result in data.get('results', {}).items():
                    batch_results.results[name] = self._deserialize_spectrum_result(result)

                # Add to saved_series
                self.main_window.saved_series[series_name] = batch_results
                logger.info(f"Loaded series '{series_name}': {len(batch_results.results)} spectra")

            # Set first series as current batch_results for backward compatibility
            if self.main_window.saved_series:
                first_name = list(self.main_window.saved_series.keys())[0]
                self.main_window.batch_results = self.main_window.saved_series[first_name]

            logger.info(f"Loaded {len(self.main_window.saved_series)} series integrations")

        else:
            # Fallback: try old format (single batch_results.json)
            batch_file = series_dir / "batch_results.json"
            if batch_file.exists():
                with open(batch_file, 'r') as f:
                    data = json.load(f)

                batch_results = BatchResults()
                batch_results.metadata = data.get('metadata', {})
                batch_results.statistics = data.get('statistics', {})
                batch_results.errors = data.get('errors', [])

                for name, result in data.get('results', {}).items():
                    batch_results.results[name] = self._deserialize_spectrum_result(result)

                # Store in both places
                self.main_window.batch_results = batch_results
                self.main_window.saved_series['default'] = batch_results

                logger.info(f"Loaded legacy series results: {len(batch_results.results)} spectra")

        return True

    def _deserialize_spectrum_result(self, result: dict) -> dict:
        """Deserialize spectrum result, converting lists back to numpy arrays.

        integration_results may contain region_2d with arrays that were
        serialized to lists. Convert them back to numpy arrays for plotting.
        """
        deserialized = {}
        for key, value in result.items():
            if key == 'integration_results' and isinstance(value, list):
                deserialized[key] = [
                    self._deserialize_integration_result(r) if isinstance(r, dict) else r
                    for r in value
                ]
            else:
                deserialized[key] = value
        return deserialized

    def _deserialize_integration_result(self, result: dict) -> dict:
        """Deserialize integration result, converting region_2d lists to arrays."""
        deserialized = {}
        for key, value in result.items():
            if key == 'region_2d' and isinstance(value, dict):
                deserialized[key] = self._deserialize_region_2d(value)
            elif key == 'fitted_2d_surface' and isinstance(value, list):
                deserialized[key] = np.array(value)
            elif key == 'individual_surfaces' and isinstance(value, list):
                deserialized[key] = [
                    np.array(s) if isinstance(s, list) else s for s in value
                ]
            else:
                deserialized[key] = value
        return deserialized

    def _deserialize_region_2d(self, region_2d: dict) -> dict:
        """Convert region_2d dict lists back to numpy arrays."""
        deserialized = {}
        # Keys expected to be numpy arrays
        array_keys = ['f1_ppm', 'f2_ppm', 'f1_grid', 'f2_grid', 'intensity']
        for key, value in region_2d.items():
            if key in array_keys and isinstance(value, list):
                deserialized[key] = np.array(value)
            else:
                deserialized[key] = value
        return deserialized

    # =========================================================================
    # DYNAMIXS STATE SAVE/LOAD
    # =========================================================================

    def _save_dynamixs_state(self, bundle_path: Path) -> bool:
        """Save DynamiXs dialog state, file references, and result files.

        DynamiXs state is stored in:
        - dynamixs/state.json: Parameter values for all pages
        - dynamixs/file_refs.json: File paths for input/output
        - dynamixs/results/: Copied result files (CSV, TXT)
        - dynamixs/results_manifest.json: Mapping of original paths to bundled files
        """

        # IMPORTANT: If DynamiXsDialog is currently open, get fresh state from it
        # (state is only transferred to main_window when dialog closes)
        dialog = getattr(self.main_window, 'dynamixs_dialog', None)
        if dialog is not None and hasattr(dialog, 'get_state'):
            try:
                state = dialog.get_state()
                file_refs = dialog.get_file_refs()
                # Also update main_window so it's in sync
                self.main_window.dynamixs_state = state
                self.main_window.dynamixs_file_refs = file_refs
                logger.info("Captured DynamiXs state from open dialog")
            except Exception as e:
                logger.warning(f"Could not get state from open DynamiXsDialog: {e}")
                state = getattr(self.main_window, 'dynamixs_state', None)
                file_refs = getattr(self.main_window, 'dynamixs_file_refs', None)
        else:
            state = getattr(self.main_window, 'dynamixs_state', None)
            file_refs = getattr(self.main_window, 'dynamixs_file_refs', None)

        # Handle mock objects from testing
        if state is not None and not isinstance(state, dict):
            state = None
        if file_refs is not None and not isinstance(file_refs, dict):
            file_refs = None

        if not state and not file_refs:
            return True  # Nothing to save

        # Create dynamixs directory
        dynamixs_dir = bundle_path / "dynamixs"
        dynamixs_dir.mkdir(exist_ok=True)

        # Save state
        if state:
            state_file = dynamixs_dir / "state.json"
            with open(state_file, 'w') as f:
                json.dump(state, f, indent=2, default=self._json_serializer)

            # Save analysis metadata separately for easy access
            if 'analysis_metadata' in state and state['analysis_metadata']:
                meta_file = dynamixs_dir / "analysis_meta.json"
                with open(meta_file, 'w') as f:
                    json.dump(state['analysis_metadata'], f, indent=2, default=self._json_serializer)

        # Save file references
        if file_refs:
            refs_file = dynamixs_dir / "file_refs.json"
            with open(refs_file, 'w') as f:
                json.dump(file_refs, f, indent=2, default=self._json_serializer)

        # Copy result files from output directories
        self._save_dynamixs_results(dynamixs_dir, state)

        logger.info(f"Saved DynamiXs state to {dynamixs_dir}")
        return True

    def _save_dynamixs_state_with_summary(self, bundle_path: Path) -> Dict[str, Any]:
        """Save DynamiXs state and return a summary of what was saved.

        Returns:
            Dict with summary information about what was saved.
        """
        summary = {
            't1t2_fitting': {},
            'spectral_density': {},
            'model_free': {},
            'methyl_t2': {},
        }

        # Get current state (from dialog if open, otherwise from main_window)
        dialog = getattr(self.main_window, 'dynamixs_dialog', None)
        if dialog is not None and hasattr(dialog, 'get_state'):
            try:
                state = dialog.get_state()
            except Exception:
                state = getattr(self.main_window, 'dynamixs_state', None)
        else:
            state = getattr(self.main_window, 'dynamixs_state', None)

        if not state:
            return {}

        # Call the actual save method
        self._save_dynamixs_state(bundle_path)

        # Extract summary from state
        # T1/T2 Fitting
        t1t2_state = state.get('t1t2', {})
        if t1t2_state:
            session_results = t1t2_state.get('session_results', {})
            fitted_experiments = t1t2_state.get('fitted_experiments', [])
            if fitted_experiments:
                summary['t1t2_fitting'] = {
                    'fitted_experiments': fitted_experiments,
                    'output_dir': t1t2_state.get('output_dir'),
                }

        # Spectral Density
        spectral_state = state.get('spectral', {})
        if spectral_state:
            session_results = spectral_state.get('session_results', {})
            if session_results and session_results.get('analysis_complete'):
                summary['spectral_density'] = {
                    'analysis_complete': True,
                    'output_dir': session_results.get('output_dir'),
                    'results_file': session_results.get('results_file'),
                }

        # Model Free (Integrated Analysis)
        integrated_state = state.get('integrated', {})
        if integrated_state:
            session_results = integrated_state.get('session_results', {})
            if session_results and session_results.get('analysis_complete'):
                summary['model_free'] = {
                    'analysis_complete': True,
                    'output_dir': session_results.get('output_dir'),
                    'json_folder': session_results.get('json_folder'),
                }

        # Methyl T2 Fitting
        methyl_state = state.get('methyl_t2', {})
        if methyl_state and methyl_state.get('last_results_file'):
            summary['methyl_t2'] = {
                'analysis_complete': True,
                'output_dir': methyl_state.get('output_dir'),
                'results_file': methyl_state.get('last_results_file'),
            }

        # Remove empty entries
        summary = {k: v for k, v in summary.items() if v}

        return summary

    def _save_dynamixs_results(self, dynamixs_dir: Path, state: dict) -> bool:
        """Copy DynamiXs result files into the project bundle.

        Args:
            dynamixs_dir: Path to dynamixs/ folder in bundle
            state: DynamiXs state dict containing output_dir paths
        """
        import shutil
        import glob

        if not state:
            return True

        results_dir = dynamixs_dir / "results"
        results_manifest = {}

        # Patterns that identify DynamiXs result files
        result_patterns = [
            '_fit_results', '_basic', '_detailed', '_results',
            '_T2_from_T1rho', '_spectral_density', '_series_tidy',
            '_data', '_intensities', 'hetNOE_ratios'
        ]

        # Collect output directories from state
        output_dirs = {}
        if 'integrated' in state and 'output_dir' in state['integrated']:
            output_dirs['integrated'] = state['integrated']['output_dir']
        # T1T2 page uses input file directory or explicit output_dir
        if 't1t2' in state and 'output_dir' in state.get('t1t2', {}):
            output_dirs['t1t2'] = state['t1t2']['output_dir']
        # Spectral Density page output_dir
        if 'spectral' in state and 'output_dir' in state.get('spectral', {}):
            output_dirs['spectral'] = state['spectral']['output_dir']
        # Methyl T2 page output_dir
        if 'methyl_t2' in state and 'output_dir' in state.get('methyl_t2', {}):
            output_dirs['methyl_t2'] = state['methyl_t2']['output_dir']

        # Also check if DynamiXsDialog is open and get output dirs from it
        dialog = getattr(self.main_window, 'dynamixs_dialog', None)
        if dialog and hasattr(dialog, 'get_output_directories'):
            try:
                dialog_dirs = dialog.get_output_directories()
                output_dirs.update(dialog_dirs)
            except Exception:
                pass  # Dialog may not be fully initialized

        # Copy result files from each output directory
        for page_name, output_dir in output_dirs.items():
            if not output_dir or not Path(output_dir).is_dir():
                continue

            # Find result files
            all_files = []
            all_files.extend(glob.glob(str(Path(output_dir) / "*.csv")))
            all_files.extend(glob.glob(str(Path(output_dir) / "*.txt")))

            # Filter to only result files
            result_files = [f for f in all_files
                           if any(p in Path(f).name for p in result_patterns)]

            if result_files:
                # Create results directory
                results_dir.mkdir(exist_ok=True)
                page_results_dir = results_dir / page_name
                page_results_dir.mkdir(exist_ok=True)

                # Copy each file
                page_manifest = {}
                for src_path in result_files:
                    src = Path(src_path)
                    dst = page_results_dir / src.name
                    shutil.copy2(src, dst)
                    page_manifest[src.name] = str(src_path)

                results_manifest[page_name] = {
                    'original_dir': output_dir,
                    'files': page_manifest
                }

            # Copy JSON folders (contain T1/T2 fit data for FitViewer)
            json_folders = glob.glob(str(Path(output_dir) / "json*"))
            for json_src in json_folders:
                json_src_path = Path(json_src)
                if json_src_path.is_dir():
                    results_dir.mkdir(exist_ok=True)
                    page_results_dir = results_dir / page_name
                    page_results_dir.mkdir(exist_ok=True)
                    json_dst = page_results_dir / json_src_path.name
                    if json_dst.exists():
                        shutil.rmtree(json_dst)
                    shutil.copytree(json_src_path, json_dst)
                    # Track in manifest
                    if page_name not in results_manifest:
                        results_manifest[page_name] = {
                            'original_dir': output_dir,
                            'files': {}
                        }
                    results_manifest[page_name]['json_folder'] = {
                        'original': str(json_src_path),
                        'bundled': str(json_dst)
                    }

        # Save manifest
        if results_manifest:
            manifest_file = dynamixs_dir / "results_manifest.json"
            with open(manifest_file, 'w') as f:
                json.dump(results_manifest, f, indent=2)

        return True

    def _load_dynamixs_state(self, bundle_path: Path) -> bool:
        """Load DynamiXs state, file references, and restore result files.

        Restores state to main_window.dynamixs_state and main_window.dynamixs_file_refs.
        Result files are copied from bundle to a restored location.
        When DynamiXsDialog is opened, it reads from these attributes.
        """
        dynamixs_dir = bundle_path / "dynamixs"

        if not dynamixs_dir.exists():
            return True  # No DynamiXs state to load

        # Load state
        state_file = dynamixs_dir / "state.json"
        if state_file.exists():
            with open(state_file, 'r') as f:
                self.main_window.dynamixs_state = json.load(f)

        # Load analysis metadata (may have been edited separately)
        meta_file = dynamixs_dir / "analysis_meta.json"
        if meta_file.exists():
            with open(meta_file, 'r') as f:
                metadata = json.load(f)
                # Merge into state (metadata file takes precedence if edited)
                if hasattr(self.main_window, 'dynamixs_state') and self.main_window.dynamixs_state:
                    self.main_window.dynamixs_state['analysis_metadata'] = metadata
                    if 'analysis_name' in metadata:
                        self.main_window.dynamixs_state['analysis_name'] = metadata['analysis_name']

        # Load file references
        refs_file = dynamixs_dir / "file_refs.json"
        if refs_file.exists():
            with open(refs_file, 'r') as f:
                self.main_window.dynamixs_file_refs = json.load(f)

        # Restore result files
        self._load_dynamixs_results(bundle_path)

        logger.info(f"Loaded DynamiXs state from {dynamixs_dir}")
        return True

    def _load_dynamixs_results(self, bundle_path: Path) -> bool:
        """Restore DynamiXs result files from bundle.

        Result files are copied to a 'dynamixs_results' folder next to the project.
        The state is updated to point to the new location.
        """
        import shutil

        dynamixs_dir = bundle_path / "dynamixs"
        results_dir = dynamixs_dir / "results"
        manifest_file = dynamixs_dir / "results_manifest.json"

        if not results_dir.exists() or not manifest_file.exists():
            return True  # No results to restore

        # Load manifest
        with open(manifest_file, 'r') as f:
            manifest = json.load(f)

        # Create restore directory next to the project bundle
        restore_base = bundle_path.parent / f"{bundle_path.stem}_dynamixs_results"

        # Copy result files and update state
        state = getattr(self.main_window, 'dynamixs_state', {}) or {}

        for page_name, page_data in manifest.items():
            page_results_dir = results_dir / page_name
            if not page_results_dir.exists():
                continue

            # Create restore directory for this page
            restore_dir = restore_base / page_name
            restore_dir.mkdir(parents=True, exist_ok=True)

            # Copy CSV/TXT files
            for filename in page_data.get('files', {}).keys():
                src = page_results_dir / filename
                dst = restore_dir / filename
                if src.exists():
                    shutil.copy2(src, dst)

            # Copy JSON folders (T1/T2 fit data)
            json_info = page_data.get('json_folder')
            restored_json_folder = None
            if json_info and 'bundled' in json_info:
                bundled_json = Path(json_info['bundled'])
                # The bundled path is absolute from save time, recalculate relative to current bundle
                json_folder_name = bundled_json.name
                src_json = page_results_dir / json_folder_name
                if src_json.exists():
                    dst_json = restore_dir / json_folder_name
                    if dst_json.exists():
                        shutil.rmtree(dst_json)
                    shutil.copytree(src_json, dst_json)
                    restored_json_folder = str(dst_json)

            # Update state to point to restored location
            if page_name in state:
                state[page_name]['output_dir'] = str(restore_dir)
                # Update session_results json_folder if present
                if 'session_results' in state[page_name] and restored_json_folder:
                    sr = state[page_name]['session_results']
                    if isinstance(sr, dict) and 'json_folder' in sr:
                        sr['json_folder'] = restored_json_folder
                # Methyl T2 tracks its viewer/results paths at the top level
                if page_name == 'methyl_t2':
                    if restored_json_folder:
                        state[page_name]['last_json_folder'] = restored_json_folder
                    rf = state[page_name].get('last_results_file')
                    if rf and (restore_dir / Path(rf).name).exists():
                        state[page_name]['last_results_file'] = str(restore_dir / Path(rf).name)

        # Update main_window state
        if state:
            self.main_window.dynamixs_state = state

        # Store the restore base path for reference
        self.main_window.dynamixs_results_dir = str(restore_base)

        return True

    def _load_dynamixs_state_with_summary(self, bundle_path: Path) -> Dict[str, Any]:
        """Load DynamiXs state and return a summary of what was loaded.

        Returns:
            Dict with summary information about what was loaded.
        """
        summary = {
            't1t2_fitting': {},
            'spectral_density': {},
            'model_free': {},
            'methyl_t2': {},
        }

        # Call the actual load method
        self._load_dynamixs_state(bundle_path)

        # Get loaded state
        state = getattr(self.main_window, 'dynamixs_state', None)
        if not state:
            return {}

        # Extract summary from state
        # T1/T2 Fitting
        t1t2_state = state.get('t1t2', {})
        if t1t2_state:
            fitted_experiments = t1t2_state.get('fitted_experiments', [])
            if fitted_experiments:
                summary['t1t2_fitting'] = {
                    'fitted_experiments': fitted_experiments,
                    'output_dir': t1t2_state.get('output_dir'),
                }

        # Spectral Density
        spectral_state = state.get('spectral', {})
        if spectral_state:
            session_results = spectral_state.get('session_results', {})
            if session_results and session_results.get('analysis_complete'):
                summary['spectral_density'] = {
                    'analysis_complete': True,
                    'output_dir': spectral_state.get('output_dir'),
                }

        # Model Free (Integrated Analysis)
        integrated_state = state.get('integrated', {})
        if integrated_state:
            session_results = integrated_state.get('session_results', {})
            if session_results and session_results.get('analysis_complete'):
                summary['model_free'] = {
                    'analysis_complete': True,
                    'output_dir': integrated_state.get('output_dir'),
                }

        # Methyl T2 Fitting
        methyl_state = state.get('methyl_t2', {})
        if methyl_state and methyl_state.get('last_results_file'):
            summary['methyl_t2'] = {
                'analysis_complete': True,
                'output_dir': methyl_state.get('output_dir'),
                'results_file': methyl_state.get('last_results_file'),
            }

        # Remove empty entries
        summary = {k: v for k, v in summary.items() if v}

        return summary

    # =========================================================================
    # KD / TITRATION STATE SAVE/LOAD
    # =========================================================================

    def _save_kd_state(self, bundle_path: Path) -> bool:
        """Save the named Kd analyses into the bundle.

        Each analysis is self-contained under kd/analyses/<name>/:
        - fit_data.json: the Kd fit JSON (raw per-point series + params + fits)
        - meta.json: page parameters and comparison-view selections
        The user adds analyses from the Kd page ("Save analysis to project");
        they are held on main_window.kd_analyses until the project is saved.
        """
        # Auto-capture the fit currently on screen (upsert by series name) so saving
        # the project always keeps the visible fit, not only ones explicitly saved.
        dialog = getattr(self.main_window, 'kd_titration_dialog', None)
        if dialog is not None and hasattr(dialog, 'ensure_current_saved'):
            try:
                dialog.ensure_current_saved()
            except Exception as e:
                logger.warning(f"Could not auto-capture current Kd fit: {e}")

        analyses = getattr(self.main_window, 'kd_analyses', None)
        if not isinstance(analyses, dict) or not analyses:
            return False

        analyses_dir = bundle_path / "kd" / "analyses"
        analyses_dir.mkdir(parents=True, exist_ok=True)
        kept = set()
        for name, entry in analyses.items():
            if not isinstance(entry, dict):
                continue
            safe = _safe_bundle_name(name)
            kept.add(safe)
            adir = analyses_dir / safe
            adir.mkdir(exist_ok=True)
            with open(adir / "fit_data.json", 'w') as f:
                json.dump(entry.get('fit_data', {}), f, indent=2,
                          default=self._json_serializer)
            with open(adir / "meta.json", 'w') as f:
                json.dump(entry.get('meta', {}), f, indent=2,
                          default=self._json_serializer)

        # Prune folders for analyses the user deleted from the session.
        import shutil
        for sub in analyses_dir.iterdir():
            if sub.is_dir() and sub.name not in kept:
                shutil.rmtree(sub, ignore_errors=True)

        logger.info(f"Saved {len(kept)} Kd analysis(es) to {analyses_dir}")
        return True

    def _load_kd_state(self, bundle_path: Path) -> bool:
        """Load the named Kd analyses from the bundle into main_window.kd_analyses.

        Each entry's fit JSON is referenced by its bundled path (fit_data_path)
        so the viewer opens the copy inside the project.
        """
        analyses_dir = bundle_path / "kd" / "analyses"
        if not analyses_dir.is_dir():
            return False

        loaded = {}
        for adir in sorted(analyses_dir.iterdir()):
            fit_path = adir / "fit_data.json"
            if not adir.is_dir() or not fit_path.exists():
                continue
            meta_path = adir / "meta.json"
            try:
                # Skip a corrupt analysis rather than aborting the whole project load.
                fit_data = json.loads(fit_path.read_text())
                meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
            except (ValueError, OSError) as e:
                logger.warning(f"Skipping unreadable Kd analysis '{adir.name}': {e}")
                continue
            loaded[adir.name] = {
                'fit_data': fit_data, 'meta': meta, 'fit_data_path': str(fit_path),
            }
        if not loaded:
            return False

        self.main_window.kd_analyses = loaded
        logger.info(f"Loaded {len(loaded)} Kd analysis(es) from {analyses_dir}")
        return True

    # =========================================================================
    # BUNDLE INVENTORY (project browser)
    # =========================================================================

    def _path_size(self, path: Path) -> int:
        """Total size in bytes of a file or directory (recursive)."""
        path = Path(path)
        if path.is_file():
            return path.stat().st_size
        if path.is_dir():
            return sum(p.stat().st_size for p in path.rglob('*') if p.is_file())
        return 0

    def inventory(self, project_path: Path) -> List[Dict[str, Any]]:
        """Describe a bundle's contents grouped by category.

        Returns a list of category dicts, each:
            {id, label, size, items: [{id, label, paths, size, removable}]}
        Only categories/items that exist on disk are included. `paths` are
        bundle-relative strings; `removable` marks items the browser may delete.
        """
        project_path = Path(project_path)
        categories: List[Dict[str, Any]] = []

        def add(cat_id, label, specs):
            items = []
            for it_id, it_label, rel_paths, removable in specs:
                present = [p for p in rel_paths if (project_path / p).exists()]
                if not present:
                    continue
                size = sum(self._path_size(project_path / p) for p in present)
                items.append({'id': it_id, 'label': it_label, 'paths': present,
                              'size': size, 'removable': removable})
            if items:
                categories.append({'id': cat_id, 'label': label,
                                   'size': sum(i['size'] for i in items),
                                   'items': items})

        add('manifest', 'Project manifest',
            [('manifest', 'project.json', ['project.json'], False)])
        add('gui_state', 'GUI state',
            [('gui_state', 'gui_state.json', ['gui_state.json'], True)])
        add('peak_list', 'Peak list',
            [('peak_list', 'peak_list.csv', ['peak_list.csv'], True)])
        add('fit_results', 'Fit results', [
            ('fit_results/fits', 'Fit parameters (fits.json)',
             ['fit_results/fits.json'], True),
            ('fit_results/arrays', 'Embedded region data',
             ['fit_results/arrays'], True),
        ])
        add('series_results', 'Series results',
            [('series_results', 'All series results', ['series_results'], True)])

        dynamixs_specs = [
            ('dynamixs/state', 'Parameters (state)',
             ['dynamixs/state.json', 'dynamixs/file_refs.json',
              'dynamixs/analysis_meta.json'], True),
        ]
        results_root = project_path / 'dynamixs' / 'results'
        if results_root.is_dir():
            for sub in sorted(results_root.iterdir()):
                if sub.is_dir():
                    dynamixs_specs.append(
                        (f'dynamixs/results/{sub.name}', f'{sub.name} results',
                         [f'dynamixs/results/{sub.name}'], True))
        add('dynamixs', 'DynamiXs', dynamixs_specs)

        kd_specs = []
        kd_analyses_root = project_path / 'kd' / 'analyses'
        if kd_analyses_root.is_dir():
            for sub in sorted(kd_analyses_root.iterdir()):
                if sub.is_dir():
                    kd_specs.append((f'kd/analyses/{sub.name}', sub.name,
                                     [f'kd/analyses/{sub.name}'], True))
        add('kd', 'Kd analyses', kd_specs)

        return categories

    def remove_bundle_paths(self, project_path: Path, rel_paths: List[str]) -> int:
        """Delete bundle-relative paths; return bytes freed.

        Raises ValueError if a path resolves outside the bundle or is the bundle
        root itself (guards against deleting unintended files).
        """
        import shutil
        project_path = Path(project_path).resolve()
        freed = 0
        for rel in rel_paths:
            target = (project_path / rel).resolve()
            if target == project_path or project_path not in target.parents:
                raise ValueError(f"Refusing to delete path outside bundle: {rel}")
            if not target.exists():
                continue
            freed += self._path_size(target)
            if target.is_dir():
                shutil.rmtree(target)
            else:
                target.unlink()
        return freed
