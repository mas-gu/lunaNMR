#!/usr/bin/env python3
"""
Background Worker Threads for DynamiXs v2.0

This module provides QThread-based workers for running analysis operations
in the background without blocking the GUI. All workers emit signals for
progress updates, completion, and error handling.

Usage:
    worker = T1T2FittingWorker(params)
    worker.progress.connect(self.update_progress)
    worker.finished.connect(self.on_complete)
    worker.error.connect(self.on_error)
    worker.start()
"""

import os
import tempfile
import traceback
from typing import Any, Callable, Dict, Optional
from dataclasses import dataclass, field

import pandas as pd
from PySide6.QtCore import QThread, Signal, QObject


# =============================================================================
# BASE WORKER CLASS
# =============================================================================

class AnalysisWorker(QThread):
    """
    Base class for background analysis workers.

    Provides common signal infrastructure for progress reporting,
    completion notification, and error handling.

    Signals:
        progress: Emits progress messages (str)
        progress_value: Emits progress percentage (int, 0-100)
        finished: Emits results on successful completion (object)
        error: Emits error message on failure (str)
        status: Emits status updates (str)
    """

    # Signals for communication with GUI
    progress = Signal(str)           # Progress message
    progress_value = Signal(int)     # Progress percentage (0-100)
    finished = Signal(object)        # Result object on completion
    error = Signal(str)              # Error message
    status = Signal(str)             # Status update

    def __init__(self, parent: Optional[QObject] = None):
        """Initialize the worker thread."""
        super().__init__(parent)
        self._cancelled = False

    def cancel(self):
        """Request cancellation of the analysis."""
        self._cancelled = True

    def is_cancelled(self) -> bool:
        """Check if cancellation was requested."""
        return self._cancelled

    def run(self):
        """
        Execute the analysis. Override in subclasses.

        This method runs in a separate thread. Use signals to communicate
        with the GUI thread.
        """
        raise NotImplementedError("Subclasses must implement run()")

    def emit_progress(self, message: str, percentage: Optional[int] = None):
        """
        Emit a progress update.

        Args:
            message: Progress message text
            percentage: Optional progress percentage (0-100)
        """
        self.progress.emit(message)
        if percentage is not None:
            self.progress_value.emit(percentage)

    def emit_error(self, error: Exception):
        """
        Emit an error with full traceback.

        Args:
            error: The exception that occurred
        """
        error_msg = f"{str(error)}\n\n{traceback.format_exc()}"
        self.error.emit(error_msg)


# =============================================================================
# GENERIC FUNCTION WORKER
# =============================================================================

class FunctionWorker(AnalysisWorker):
    """
    Worker that executes an arbitrary function in a background thread.

    This is useful for wrapping existing analysis functions without
    creating custom worker classes.

    Example:
        def my_analysis(params, progress_callback):
            progress_callback("Starting...")
            result = do_analysis(params)
            progress_callback("Complete!")
            return result

        worker = FunctionWorker(my_analysis, params)
        worker.finished.connect(handle_result)
        worker.start()
    """

    def __init__(
        self,
        func: Callable,
        *args,
        progress_callback: bool = True,
        **kwargs
    ):
        """
        Initialize the function worker.

        Args:
            func: The function to execute
            *args: Positional arguments for the function
            progress_callback: If True, pass progress callback as first arg
            **kwargs: Keyword arguments for the function
        """
        super().__init__()
        self._func = func
        self._args = args
        self._kwargs = kwargs
        self._use_callback = progress_callback

    def run(self):
        """Execute the wrapped function."""
        try:
            if self._use_callback:
                # Pass progress callback as first argument
                result = self._func(
                    self.emit_progress,
                    *self._args,
                    **self._kwargs
                )
            else:
                result = self._func(*self._args, **self._kwargs)

            if not self._cancelled:
                self.finished.emit(result)

        except Exception as e:
            self.emit_error(e)


# =============================================================================
# T1/T2 FITTING WORKER
# =============================================================================

@dataclass
class T1T2FittingParams:
    """Parameters for T1/T2 fitting analysis.

    Data can be provided in three ways (checked in order):
    1. input_df: DataFrame already in fitting format (pivot with Assignment + delay columns)
    2. series_name: Name of series to look up from lunaNMR saved_series
    3. input_file: Path to CSV file

    When input_df or series_name is used, a temporary CSV file is created
    for the fitting module.
    """
    input_file: str = ""  # CSV file path (legacy mode)
    output_dir: str = ""
    experiment_type: str = "T1"  # "T1" or "T2"
    results_prefix: str = ""
    multicore: bool = True
    initial_amplitude: float = 1.0
    initial_decay: float = 1000.0
    n_bootstrap: int = 1000
    error_method: str = "analytical"  # "analytical" or "bootstrap"
    json_folder: str = ""
    field_name: str = ""
    field_freq: float = 0.0
    # New fields for memory-based data flow
    input_df: Optional[pd.DataFrame] = field(default=None, repr=False)
    series_name: str = ""


class T1T2FittingWorker(AnalysisWorker):
    """
    Worker for T1/T2 exponential decay fitting analysis.

    This worker runs the T1 or T2 fitting scripts in a background thread,
    reporting progress and returning the fit results.
    """

    def __init__(self, params: T1T2FittingParams, parent: Optional[QObject] = None):
        """
        Initialize the T1/T2 fitting worker.

        Args:
            params: Fitting parameters
            parent: Parent QObject
        """
        super().__init__(parent)
        self.params = params

    def run(self):
        """Execute the T1/T2 fitting analysis."""
        temp_csv_path = None
        try:
            self.status.emit("Initializing fitting...")
            self.emit_progress("Loading data...", 5)

            # Import the fitting module
            if self.params.multicore:
                from dynamiXs_T1_T2 import fitmulti__Tx_NMRRE as fitting_module
            else:
                from dynamiXs_T1_T2 import fit_Tx_NMRRE as fitting_module

            # Determine input CSV file path
            # Priority: input_df > input_file
            if self.params.input_df is not None:
                # Write DataFrame to temp file for fitting module
                fd, temp_csv_path = tempfile.mkstemp(suffix='.csv', prefix='dynamixs_')
                os.close(fd)
                self.params.input_df.to_csv(temp_csv_path, index=False)
                input_csv_file = temp_csv_path
                self.emit_progress("Loaded data from memory...", 8)
            elif self.params.input_file:
                input_csv_file = self.params.input_file
            else:
                raise ValueError("No input data provided: need input_df or input_file")

            self.emit_progress("Starting fit...", 10)

            # Build output file paths
            output_prefix = os.path.join(self.params.output_dir, self.params.results_prefix)
            results_txt_file = f"{output_prefix}_fit_results.txt"

            # Prepare parameters dictionary matching run_analysis_with_params() expected keys
            fit_params = {
                'input_csv_file': input_csv_file,
                'output_prefix': output_prefix,
                'results_txt_file': results_txt_file,
                'experiment_type': self.params.experiment_type,
                'time_units': 'ms',
                'signal_units': 'Intensity',
                'initial_A': self.params.initial_amplitude,
                'initial_t2': self.params.initial_decay,
                'n_bootstrap': self.params.n_bootstrap,
                'error_method': self.params.error_method,
                'n_plots_per_figure': 20,
            }

            # Add JSON export parameters if provided
            if self.params.json_folder:
                fit_params['json_folder'] = self.params.json_folder
                fit_params['field_name'] = self.params.field_name
                fit_params['field_freq'] = self.params.field_freq

            # Create progress callback that emits signals
            def progress_callback(completed, total, residue_name, message):
                if self._cancelled:
                    return
                # Calculate progress: 10% for setup, 85% for fitting, 5% for finalization
                fitting_progress = int((completed / total) * 85) + 10
                self.emit_progress(f"Fitting {residue_name} ({completed}/{total})", fitting_progress)

            # Run the analysis with progress callback
            result = fitting_module.run_analysis_with_params(fit_params, progress_callback=progress_callback)

            if not self._cancelled:
                self.emit_progress("Fitting complete!", 100)
                # Augment result with session tracking fields
                result['field_name'] = self.params.field_name or 'field1'
                result['experiment_type'] = self.params.experiment_type
                result['output_dir'] = self.params.output_dir
                result['json_folder'] = self.params.json_folder
                result['field_freq'] = self.params.field_freq
                self.finished.emit(result)

        except Exception as e:
            self.emit_error(e)
        finally:
            # Clean up temp file if created
            if temp_csv_path and os.path.exists(temp_csv_path):
                try:
                    os.remove(temp_csv_path)
                except OSError:
                    pass  # Best effort cleanup


# =============================================================================
# METHYL T2 BI-EXP WORKER
# =============================================================================

@dataclass
class MethylT2FittingParams:
    """Parameters for the shared-amplitude methyl bi-exp T2 fitter.

    Model: I(t) = 0.5 * A * (exp(-t/T2a) + exp(-t/T2b)).
    Same data-source rules as T1T2FittingParams: input_df > input_file.
    """
    input_file: str = ""
    output_dir: str = ""
    results_prefix: str = ""
    json_folder: str = ""
    field_name: str = "field1"
    field_freq: float = 0.0
    initial_A: Optional[float] = None  # None -> data-driven default = span(y)
    initial_t2_a: float = 100.0  # ms (slow component)
    initial_t2_b: float = 20.0   # ms (fast component)
    n_bootstrap: int = 1000
    error_method: str = "analytical"
    input_df: Optional[pd.DataFrame] = field(default=None, repr=False)
    series_name: str = ""


class MethylT2FittingWorker(AnalysisWorker):
    """Worker for the methyl bi-exp T2 fitter.

    Calls dynamiXs_T1_T2.fit_methyl_T2.run_methyl_t2_analysis_with_params on a
    background thread, reporting per-residue progress.
    """

    def __init__(self, params: MethylT2FittingParams, parent: Optional[QObject] = None):
        super().__init__(parent)
        self.params = params

    def run(self):
        temp_csv_path = None
        try:
            self.status.emit("Initializing methyl T2 fit...")
            self.emit_progress("Loading data...", 5)

            from dynamiXs_T1_T2 import fit_methyl_T2 as fitting_module

            if self.params.input_df is not None:
                fd, temp_csv_path = tempfile.mkstemp(suffix=".csv", prefix="dynamixs_methyl_")
                os.close(fd)
                self.params.input_df.to_csv(temp_csv_path, index=False)
                input_csv = temp_csv_path
            elif self.params.input_file:
                input_csv = self.params.input_file
            else:
                raise ValueError("No input data provided: need input_df or input_file")

            self.emit_progress("Starting fit...", 10)

            output_prefix = os.path.join(self.params.output_dir, self.params.results_prefix)
            results_txt_file = f"{output_prefix}_fit_results.txt"

            fit_params = {
                "input_csv_file": input_csv,
                "output_prefix": output_prefix,
                "results_txt_file": results_txt_file,
                "json_folder": self.params.json_folder,
                "field_name": self.params.field_name,
                "field_freq": self.params.field_freq,
                "time_units": "ms",
                "signal_units": "Intensity",
                "initial_A": self.params.initial_A,
                "initial_t2_a": self.params.initial_t2_a,
                "initial_t2_b": self.params.initial_t2_b,
                "n_bootstrap": self.params.n_bootstrap,
                "error_method": self.params.error_method,
            }

            def progress_callback(completed, total, residue_name, message):
                if self._cancelled:
                    return
                pct = int((completed / total) * 85) + 10
                self.emit_progress(f"Fitting {residue_name} ({completed}/{total})", pct)

            result = fitting_module.run_methyl_t2_analysis_with_params(
                fit_params, progress_callback=progress_callback
            )

            if not self._cancelled:
                self.emit_progress("Fitting complete!", 100)
                result["field_name"] = self.params.field_name
                result["experiment_type"] = "methylT2"
                result["json_folder"] = self.params.json_folder
                result["field_freq"] = self.params.field_freq
                self.finished.emit(result)

        except Exception as e:
            self.emit_error(e)
        finally:
            if temp_csv_path and os.path.exists(temp_csv_path):
                try:
                    os.remove(temp_csv_path)
                except OSError:
                    pass


# =============================================================================
# SPECTRAL DENSITY WORKER
# =============================================================================

@dataclass
class SpectralDensityParams:
    """Parameters for spectral density analysis."""
    input_file: str  # Field 1 input file
    input_file2: str = ""  # Field 2 input file (for dual-field)
    output_dir: str = ""
    results_prefix: str = "spectral_density"
    field1_freq: float = 600.0
    field2_freq: float = 700.0
    dual_field: bool = True
    use_087: bool = True
    r_nh: float = 1.015
    csa: float = -172.0
    monte_carlo: bool = False
    n_samples: int = 1000


class SpectralDensityWorker(AnalysisWorker):
    """
    Worker for spectral density mapping analysis.

    Runs the Farrow 1995 spectral density calculations in a background
    thread, supporting both single-field and dual-field analysis.
    """

    def __init__(self, params: SpectralDensityParams, parent: Optional[QObject] = None):
        """
        Initialize the spectral density worker.

        Args:
            params: Analysis parameters
            parent: Parent QObject
        """
        super().__init__(parent)
        self.params = params

    def run(self):
        """Execute the spectral density analysis."""
        import os
        try:
            self.status.emit("Initializing spectral density analysis...")
            self.emit_progress("Loading data...", 10)

            # Convert physical parameters to proper units
            rNH_meters = self.params.r_nh * 1e-10  # Angstroms to meters
            csaN_units = self.params.csa * 1e-6   # ppm to dimensionless

            # Build output paths
            output_prefix = os.path.join(self.params.output_dir, self.params.results_prefix)

            self.emit_progress("Computing spectral densities...", 30)

            if self.params.dual_field:
                # Dual-field analysis requires two input files
                if not self.params.input_file2:
                    raise ValueError("Dual-field analysis requires two input files")

                if self.params.use_087:
                    from dynamiXs_density_functions.ZZ_multi_2fields_density_087 import DualFieldSpectralDensityAnalysis
                else:
                    from dynamiXs_density_functions.ZZ_multi_2fields_density import DualFieldSpectralDensityAnalysis

                # Create analyzer
                analyzer = DualFieldSpectralDensityAnalysis(
                    field1_freq=self.params.field1_freq,
                    field2_freq=self.params.field2_freq,
                    rNH=rNH_meters,
                    csaN=csaN_units
                )

                self.emit_progress("Running dual-field analysis...", 50)

                # Run analysis
                results_df = analyzer.analyze_dual_field_csv(
                    csv_file1=self.params.input_file,
                    csv_file2=self.params.input_file2,
                    use_monte_carlo_errors=self.params.monte_carlo,
                    n_monte_carlo=self.params.n_samples,
                    use_multiprocessing=True
                )

                # Save results
                self.emit_progress("Saving results...", 80)
                basic_csv = f"{output_prefix}_basic.csv"
                detailed_csv = f"{output_prefix}_detailed.csv"
                plots_pdf = f"{output_prefix}_plots.pdf"

                results_df.to_csv(basic_csv, index=False)
                analyzer.save_dual_field_results(results_df, detailed_csv)
                analyzer.plot_dual_field_results(results_df, save_plots=True, plot_filename=plots_pdf)

                result = {
                    'n_residues': len(results_df),
                    'output_file': detailed_csv,
                    'plots_file': plots_pdf
                }

            else:
                # Single-field analysis
                if self.params.use_087:
                    from dynamiXs_density_functions.ZZ_multi_density_087 import ReducedSpectralDensityAnalysis
                else:
                    from dynamiXs_density_functions.ZZ_multi_density import ReducedSpectralDensityAnalysis

                # Create analyzer
                analyzer = ReducedSpectralDensityAnalysis(
                    spectrometer_frequency=self.params.field1_freq,
                    rNH=rNH_meters,
                    csaN=csaN_units
                )

                self.emit_progress("Running single-field analysis...", 50)

                # Run analysis
                results_df = analyzer.analyze_csv(
                    csv_file=self.params.input_file,
                    use_monte_carlo_errors=self.params.monte_carlo,
                    n_monte_carlo=self.params.n_samples,
                    use_multiprocessing=True
                )

                # Save results
                self.emit_progress("Saving results...", 80)
                output_csv = f"{output_prefix}_results.csv"
                plots_pdf = f"{output_prefix}_plots.pdf"

                results_df.to_csv(output_csv, index=False)
                analyzer.plot_results(results_df, save_plots=True, plot_filename=plots_pdf)

                result = {
                    'n_residues': len(results_df),
                    'output_file': output_csv,
                    'plots_file': plots_pdf
                }

            if not self._cancelled:
                self.emit_progress("Analysis complete!", 100)
                self.finished.emit(result)

        except Exception as e:
            self.emit_error(e)


# =============================================================================
# INTEGRATED ANALYSIS WORKER
# =============================================================================

@dataclass
class IntegratedAnalysisParams:
    """Parameters for integrated model-free analysis."""
    # Field 1 data
    field1_t1_file: str = ""
    field1_t2_file: str = ""
    field1_noe_sat_file: str = ""
    field1_noe_unsat_file: str = ""
    field1_freq: float = 600.0

    # Field 2 data (optional for dual-field)
    field2_t1_file: str = ""
    field2_t2_file: str = ""
    field2_noe_sat_file: str = ""
    field2_noe_unsat_file: str = ""
    field2_freq: float = 700.0
    dual_field: bool = True

    # Analysis settings
    analysis_method: str = "dual_087"  # single, single_087, dual, dual_087
    r_nh: float = 1.015
    csa: float = -172.0

    # Fitting parameters
    initial_amplitude: float = 5.0
    initial_t1: float = 800.0
    initial_t2: float = 100.0
    n_bootstrap: int = 1000
    error_method: str = "analytical"  # "analytical" or "bootstrap"
    n_monte_carlo: int = 50

    # Output settings
    output_dir: str = ""
    results_prefix: str = "analysis"
    json_folder: str = ""


class IntegratedAnalysisWorker(AnalysisWorker):
    """
    Worker for complete integrated model-free analysis.

    This worker orchestrates the full analysis pipeline:
    1. T1/T2 fitting for both fields
    2. hetNOE calculation
    3. Spectral density mapping
    4. Model-free parameter extraction
    """

    def __init__(self, params: IntegratedAnalysisParams, parent: Optional[QObject] = None):
        """
        Initialize the integrated analysis worker.

        Args:
            params: Analysis parameters
            parent: Parent QObject
        """
        super().__init__(parent)
        self.params = params

    def run(self):
        """Execute the complete integrated analysis pipeline."""
        try:
            self.status.emit("Starting integrated analysis...")

            # Import the integrated analysis pipeline
            from dynamiXs_integrated.integrated_analysis import (
                IntegratedAnalysisPipeline,
                IntegratedAnalysisParameters
            )

            # Create pipeline parameters and set attributes
            pipeline_params = IntegratedAnalysisParameters()

            # Field 1 files and frequency
            pipeline_params.field1_t1_file = self.params.field1_t1_file
            pipeline_params.field1_t2_file = self.params.field1_t2_file
            pipeline_params.field1_noe_sat_file = self.params.field1_noe_sat_file
            pipeline_params.field1_noe_unsat_file = self.params.field1_noe_unsat_file
            pipeline_params.field1_freq_mhz = self.params.field1_freq

            # Field 2 files and frequency
            pipeline_params.field2_t1_file = self.params.field2_t1_file
            pipeline_params.field2_t2_file = self.params.field2_t2_file
            pipeline_params.field2_noe_sat_file = self.params.field2_noe_sat_file
            pipeline_params.field2_noe_unsat_file = self.params.field2_noe_unsat_file
            pipeline_params.field2_freq_mhz = self.params.field2_freq
            pipeline_params.enable_dual_field = self.params.dual_field

            # Analysis method
            pipeline_params.analysis_method = self.params.analysis_method

            # Physical parameters
            pipeline_params.rNH_angstrom = self.params.r_nh
            pipeline_params.csaN_ppm = self.params.csa

            # Fitting parameters
            pipeline_params.t1_initial_amplitude = self.params.initial_amplitude
            pipeline_params.t1_initial_time = self.params.initial_t1
            pipeline_params.t1_bootstrap_iterations = self.params.n_bootstrap
            pipeline_params.t2_initial_amplitude = self.params.initial_amplitude
            pipeline_params.t2_initial_time = self.params.initial_t2
            pipeline_params.t2_bootstrap_iterations = self.params.n_bootstrap
            pipeline_params.error_method = self.params.error_method
            pipeline_params.monte_carlo_iterations = self.params.n_monte_carlo

            # Output settings
            pipeline_params.output_dir = self.params.output_dir
            pipeline_params.output_prefix = self.params.results_prefix
            pipeline_params.json_folder = self.params.json_folder

            # Define progress callback
            def progress_cb(msg, pct=None):
                if not self._cancelled:
                    self.emit_progress(msg, pct)

            # Create pipeline with progress callback and run
            pipeline = IntegratedAnalysisPipeline(pipeline_params, progress_callback=progress_cb)
            result = pipeline.run_complete_analysis()

            if not self._cancelled:
                self.emit_progress("Integrated analysis complete!", 100)
                self.finished.emit(result)

        except Exception as e:
            self.emit_error(e)


# =============================================================================
# CPMG ANALYSIS WORKER
# =============================================================================

@dataclass
class CPMGAnalysisParams:
    """Parameters for CPMG relaxation dispersion analysis."""
    input_file: str
    output_dir: str
    results_prefix: str
    field_freq: float = 700.0
    single_field: bool = True
    n_bootstrap: int = 1000


class CPMGAnalysisWorker(AnalysisWorker):
    """
    Worker for CPMG relaxation dispersion analysis.
    """

    def __init__(self, params: CPMGAnalysisParams, parent: Optional[QObject] = None):
        """
        Initialize the CPMG analysis worker.

        Args:
            params: Analysis parameters
            parent: Parent QObject
        """
        super().__init__(parent)
        self.params = params

    def run(self):
        """Execute the CPMG analysis."""
        try:
            self.status.emit("Initializing CPMG analysis...")
            self.emit_progress("Loading data...", 10)

            # Import the appropriate CPMG module
            if self.params.single_field:
                from dynamiXs_cpmg import cpmg_RD as cpmg_module
            else:
                from dynamiXs_cpmg import cpmg_RD_bootstrap as cpmg_module

            self.emit_progress("Fitting CPMG data...", 30)

            def progress_cb(msg, pct=None):
                if not self._cancelled:
                    self.emit_progress(msg, pct)

            result = cpmg_module.run_analysis(
                input_file=self.params.input_file,
                output_dir=self.params.output_dir,
                prefix=self.params.results_prefix,
                field_freq=self.params.field_freq,
                n_bootstrap=self.params.n_bootstrap,
                progress_callback=progress_cb
            )

            if not self._cancelled:
                self.emit_progress("CPMG analysis complete!", 100)
                self.finished.emit(result)

        except Exception as e:
            self.emit_error(e)


# =============================================================================
# FILE LOADING WORKER
# =============================================================================

class FileLoadWorker(AnalysisWorker):
    """
    Worker for loading and validating data files.

    Use this for loading large CSV files or performing validation
    without blocking the GUI.
    """

    def __init__(
        self,
        file_path: str,
        validator: Optional[Callable] = None,
        parent: Optional[QObject] = None
    ):
        """
        Initialize the file load worker.

        Args:
            file_path: Path to file to load
            validator: Optional validation function
            parent: Parent QObject
        """
        super().__init__(parent)
        self.file_path = file_path
        self.validator = validator

    def run(self):
        """Load and validate the file."""
        try:
            import pandas as pd

            self.emit_progress(f"Loading {self.file_path}...", 20)

            # Load the file
            df = pd.read_csv(self.file_path)

            self.emit_progress("Validating data...", 60)

            # Run validation if provided
            if self.validator:
                validation_result = self.validator(df)
                if not validation_result.get('valid', True):
                    raise ValueError(validation_result.get('error', 'Validation failed'))

            self.emit_progress("File loaded successfully!", 100)
            self.finished.emit({
                'file_path': self.file_path,
                'data': df,
                'n_rows': len(df),
                'columns': list(df.columns)
            })

        except Exception as e:
            self.emit_error(e)


# =============================================================================
# BATCH PROCESSING WORKER
# =============================================================================

class BatchProcessingWorker(AnalysisWorker):
    """
    Worker for processing multiple files or residues in batch.

    Provides progress updates for each item in the batch.
    """

    def __init__(
        self,
        items: list,
        process_func: Callable,
        parent: Optional[QObject] = None
    ):
        """
        Initialize the batch processing worker.

        Args:
            items: List of items to process
            process_func: Function to call for each item
            parent: Parent QObject
        """
        super().__init__(parent)
        self.items = items
        self.process_func = process_func

    def run(self):
        """Process all items in the batch."""
        try:
            results = []
            total = len(self.items)

            for i, item in enumerate(self.items):
                if self._cancelled:
                    break

                progress = int((i / total) * 100)
                self.emit_progress(f"Processing item {i+1}/{total}...", progress)

                result = self.process_func(item)
                results.append(result)

            if not self._cancelled:
                self.emit_progress("Batch processing complete!", 100)
                self.finished.emit(results)

        except Exception as e:
            self.emit_error(e)


# =============================================================================
# SPECTRA PROCESSING WORKER (LunaNMR Integration)
# =============================================================================

@dataclass
class SpectraProcessingParams:
    """Parameters for NMR spectra processing via LunaNMR."""
    experiment_folders: Dict[str, str]  # {experiment_type: folder_path}
    peak_list_path: str
    field_mhz: float
    output_dir: str
    intensity_column: str = 'height'  # 'height' or 'volume'
    use_voigt: bool = True
    use_parallel: bool = True
    fix_positions: bool = False
    fix_linewidths: bool = False


class SpectraProcessingWorker(AnalysisWorker):
    """
    Worker for processing NMR spectra folders using LunaNMR.

    This worker:
    1. Runs LunaNMR series integration on each experiment folder
    2. Converts the tidy CSV output to DynamiXs format
    3. Returns paths to generated DynamiXs-compatible CSV files
    """

    def __init__(self, params: SpectraProcessingParams, parent: Optional[QObject] = None):
        """
        Initialize the spectra processing worker.

        Args:
            params: Processing parameters
            parent: Parent QObject
        """
        super().__init__(parent)
        self.params = params

    def run(self):
        """Execute the spectra processing pipeline."""
        import os
        import sys
        from pathlib import Path

        try:
            self.status.emit("Initializing LunaNMR integration...")
            self.emit_progress("Loading LunaNMR modules...", 5)

            # Add lunaNMR to path
            lunaNMR_root = Path(__file__).parent.parent.parent
            if str(lunaNMR_root) not in sys.path:
                sys.path.insert(0, str(lunaNMR_root))

            from lunaNMR.utils.relaxation_data_bridge import RelaxationDataBridge
            from lunaNMR.utils.delay_extractor import DelayExtractor

            # Initialize bridge
            bridge = RelaxationDataBridge(self.params.peak_list_path)
            extractor = DelayExtractor()

            # Create output directory
            os.makedirs(self.params.output_dir, exist_ok=True)

            results = {
                'processed_experiments': [],
                'output_files': {},
                'errors': []
            }

            # Process each experiment type
            experiments_to_process = list(self.params.experiment_folders.items())
            total_experiments = len(experiments_to_process)

            for idx, (exp_type, folder_path) in enumerate(experiments_to_process):
                if self._cancelled:
                    break

                base_progress = int((idx / total_experiments) * 80) + 10
                self.emit_progress(f"Processing {exp_type}...", base_progress)

                try:
                    if exp_type in ['T1', 'T2']:
                        # T1/T2 processing: run series integration then convert
                        output_file = self._process_t1_t2_folder(
                            exp_type, folder_path, bridge, extractor, base_progress
                        )
                        results['output_files'][exp_type] = output_file
                        results['processed_experiments'].append(exp_type)

                    elif exp_type in ['hetNOE_sat', 'hetNOE_unsat']:
                        # hetNOE: just extract intensities (no series fitting needed)
                        output_file = self._process_hetnoe_folder(
                            exp_type, folder_path, bridge, base_progress
                        )
                        results['output_files'][exp_type] = output_file
                        results['processed_experiments'].append(exp_type)

                except Exception as e:
                    error_msg = f"{exp_type}: {str(e)}"
                    results['errors'].append(error_msg)
                    self.emit_progress(f"Error processing {exp_type}: {str(e)}", base_progress)

            # Generate hetNOE ratio file if both sat and unsat were processed
            if ('hetNOE_sat' in results['output_files'] and
                'hetNOE_unsat' in results['output_files']):
                try:
                    self.emit_progress("Calculating hetNOE ratios...", 95)
                    noe_output = self._calculate_hetnoe_ratio(
                        results['output_files']['hetNOE_sat'],
                        results['output_files']['hetNOE_unsat']
                    )
                    results['output_files']['hetNOE'] = noe_output
                except Exception as e:
                    results['errors'].append(f"hetNOE ratio: {str(e)}")

            if not self._cancelled:
                self.emit_progress("Processing complete!", 100)
                self.finished.emit(results)

        except Exception as e:
            self.emit_error(e)

    def _process_t1_t2_folder(
        self,
        exp_type: str,
        folder_path: str,
        bridge: 'RelaxationDataBridge',
        extractor: 'DelayExtractor',
        base_progress: int
    ) -> str:
        """
        Process a T1 or T2 folder.

        For now, this creates a mock tidy CSV from the folder contents.
        Full LunaNMR integration would call MultiSpectrumProcessor here.
        """
        import os
        import pandas as pd

        # Scan folder for files with delays
        files_with_delays = extractor.scan_folder(folder_path)

        if not files_with_delays:
            raise ValueError(f"No NMR files with delay patterns found in {folder_path}")

        self.emit_progress(f"Found {len(files_with_delays)} {exp_type} spectra", base_progress + 2)

        # For now, create a placeholder indicating integration is needed
        # Full implementation would call:
        #   from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
        #   processor = MultiSpectrumProcessor(...)
        #   processor.process_nmr_series(...)

        # Create a stub tidy CSV for testing (in production, this comes from LunaNMR)
        # This allows the GUI integration to proceed while full series integration is wired up

        output_tidy = os.path.join(self.params.output_dir, f"{exp_type}_series_tidy.csv")
        output_dynamixs = os.path.join(self.params.output_dir, f"{exp_type}_data.csv")

        # Check if tidy file already exists (from previous LunaNMR run)
        if os.path.exists(output_tidy):
            self.emit_progress(f"Converting existing {exp_type} tidy CSV...", base_progress + 5)
            result = bridge.convert_tidy_to_dynamixs_format(
                tidy_csv_path=output_tidy,
                output_csv_path=output_dynamixs,
                intensity_column=self.params.intensity_column
            )
            self.emit_progress(
                f"{exp_type}: {result['residue_count']} residues, {result['delay_count']} delays",
                base_progress + 8
            )
            return output_dynamixs
        else:
            # No existing tidy file - would need to run series integration
            raise ValueError(
                f"Series integration not yet run for {exp_type}. "
                f"Please run 'Fit Series' from LunaNMR first, or use the standalone "
                f"series integration dialog."
            )

    def _process_hetnoe_folder(
        self,
        exp_type: str,
        folder_path: str,
        bridge: 'RelaxationDataBridge',
        base_progress: int
    ) -> str:
        """
        Process a hetNOE folder (sat or unsat).

        For hetNOE, we just need single-spectrum integration, not series.
        """
        import os

        output_file = os.path.join(self.params.output_dir, f"{exp_type}_intensities.csv")

        # Check if integration results already exist
        # In production, would check for existing peak integration results
        tidy_file = os.path.join(folder_path, "integration_results.csv")

        if os.path.exists(tidy_file):
            # Convert from tidy format
            import pandas as pd
            df = pd.read_csv(tidy_file)

            # Extract assignment and intensity columns
            with open(output_file, 'w') as f:
                for _, row in df.iterrows():
                    assignment = str(row.get('assignment', row.get('Assignment', '')))
                    intensity = row.get(self.params.intensity_column,
                                        row.get('Height', row.get('height', 0)))
                    residue = bridge.format_residue_id(assignment)
                    f.write(f"{residue},{intensity}\n")

            return output_file
        else:
            raise ValueError(
                f"Integration not yet run for {exp_type}. "
                f"Please run peak integration from LunaNMR first."
            )

    def _calculate_hetnoe_ratio(self, sat_file: str, unsat_file: str) -> str:
        """Calculate hetNOE ratios from sat and unsat intensity files."""
        import os
        import pandas as pd

        # Load intensities
        sat_data = {}
        unsat_data = {}

        with open(sat_file, 'r') as f:
            for line in f:
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    sat_data[parts[0]] = float(parts[1])

        with open(unsat_file, 'r') as f:
            for line in f:
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    unsat_data[parts[0]] = float(parts[1])

        # Calculate ratios
        output_file = os.path.join(self.params.output_dir, "hetNOE_ratios.csv")

        with open(output_file, 'w') as f:
            f.write("Residue,hetNOE,hetNOE_err\n")
            for residue in sorted(sat_data.keys(), key=lambda x: int(x) if x.isdigit() else 0):
                if residue in unsat_data and unsat_data[residue] > 0:
                    ratio = sat_data[residue] / unsat_data[residue]
                    # Estimate 2% error
                    error = ratio * 0.02
                    f.write(f"{residue},{ratio:.5f},{error:.5f}\n")

        return output_file
