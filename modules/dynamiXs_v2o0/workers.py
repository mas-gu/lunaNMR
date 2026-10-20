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
from typing import Callable, Dict, Optional
from dataclasses import dataclass, field

import pandas as pd
from PySide6.QtCore import QThread, Signal, QObject


# BASE WORKER CLASS

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


# GENERIC FUNCTION WORKER

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


# T1/T2 FITTING WORKER

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
    input_file: str = ""  # CSV file path (single-file mode)
    output_dir: str = ""
    experiment_type: str = "T1"  # "T1" or "T2"
    results_prefix: str = ""
    multicore: bool = True
    initial_amplitude: float = 1.0
    initial_decay: float = 1000.0  # in time_units
    time_units: str = "ms"  # units of the series' delays
    n_bootstrap: int = 1000
    error_method: str = "analytical"  # "analytical" or "bootstrap"
    json_folder: str = ""
    field_name: str = ""
    # Starts the fit JSON's name, before the field (<run_prefix>_<field>_<exp>), as the
    # CLI's --prefix does: runs sharing a JSON folder keep their own.
    run_prefix: str = ""
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

            if self.params.multicore:
                from dynamiXs_T1_T2 import fitmulti__Tx_NMRRE as fitting_module
            else:
                from dynamiXs_T1_T2 import fit_Tx_NMRRE as fitting_module

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

            output_prefix = os.path.join(self.params.output_dir, self.params.results_prefix)
            results_txt_file = f"{output_prefix}_fit_results.txt"

            fit_params = {
                'input_csv_file': input_csv_file,
                'output_prefix': output_prefix,
                'results_txt_file': results_txt_file,
                'experiment_type': self.params.experiment_type,
                'time_units': self.params.time_units,
                'signal_units': 'Intensity',
                'initial_A': self.params.initial_amplitude,
                'initial_t2': self.params.initial_decay,
                'n_bootstrap': self.params.n_bootstrap,
                'error_method': self.params.error_method,
                'n_plots_per_figure': 20,
            }

            if self.params.json_folder:
                fit_params['json_folder'] = self.params.json_folder
                fit_params['field_name'] = "_".join(
                    p for p in (self.params.run_prefix, self.params.field_name) if p)
                fit_params['field_freq'] = self.params.field_freq

            def progress_callback(completed, total, residue_name, message):
                if self._cancelled:
                    return
                # Calculate progress: 10% for setup, 85% for fitting, 5% for finalization
                fitting_progress = int((completed / total) * 85) + 10
                self.emit_progress(f"Fitting {residue_name} ({completed}/{total})", fitting_progress)

            result = fitting_module.run_analysis_with_params(fit_params, progress_callback=progress_callback)

            if not self._cancelled:
                self.emit_progress("Fitting complete!", 100)
                # Augment result with session tracking fields
                result['field_name'] = self.params.field_name or 'field1'
                result['experiment_type'] = self.params.experiment_type
                result['output_dir'] = self.params.output_dir
                result['json_folder'] = self.params.json_folder
                result['field_freq'] = self.params.field_freq
                # What every file this run wrote starts with.
                result['file_prefix'] = self.params.results_prefix
                self.finished.emit(result)

        except Exception as e:
            self.emit_error(e)
        finally:
            if temp_csv_path and os.path.exists(temp_csv_path):
                try:
                    os.remove(temp_csv_path)
                except OSError:
                    pass  # Best effort cleanup


# METHYL T2 BI-EXP WORKER

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
    initial_t2_a: float = 100.0  # slow component, in time_units
    initial_t2_b: float = 20.0   # fast component, in time_units
    time_units: str = "ms"  # units of the series' delays
    n_bootstrap: int = 1000
    error_method: str = "analytical"
    input_df: Optional[pd.DataFrame] = field(default=None, repr=False)
    series_name: str = ""


@dataclass
class KdTitrationFittingParams:
    """Parameters for the Kd / titration binding fit (CSP and intensity ratio)."""
    input_file: str = ""
    output_dir: str = ""
    output_prefix: str = "kd"
    concentrations: Optional[list] = None   # one [L] per titration point (sorted order)
    conc_units: str = "absolute"            # 'absolute' or 'equivalents' of [P]0
    intensity_scales: Optional[list] = None  # per-point height/volume scale (scan-count correction)
    protein_conc: float = 50.0              # [P]0
    alpha: float = 0.14                     # 15N -> 1H CSP scaling
    observables: Optional[list] = None      # subset of ['csp', 'intensity']
    intensity_value: str = "height"         # 'height' or 'volume'
    n_bootstrap: int = 0
    thresholds: Optional[dict] = None       # kd_params.THRESHOLD_KEYS -> value
    input_df: Optional[pd.DataFrame] = field(default=None, repr=False)
    series_name: str = ""


class KdTitrationFittingWorker(AnalysisWorker):
    """Worker for the Kd titration fitter (dynamiXs_Kd.kd_fit)."""

    def __init__(self, params: KdTitrationFittingParams, parent: Optional[QObject] = None):
        super().__init__(parent)
        self.params = params

    def run(self):
        import sys
        temp_csv_path = None
        try:
            self.status.emit("Initializing Kd titration fit...")
            self.emit_progress("Loading data...", 5)

            kd_dir = os.path.join(os.path.dirname(__file__), "dynamiXs_Kd")
            if kd_dir not in sys.path:
                sys.path.insert(0, kd_dir)
            import kd_fit

            if self.params.input_df is not None:
                fd, temp_csv_path = tempfile.mkstemp(suffix=".csv", prefix="dynamixs_kd_")
                os.close(fd)
                self.params.input_df.to_csv(temp_csv_path, index=False)
                input_csv = temp_csv_path
            elif self.params.input_file:
                input_csv = self.params.input_file
            else:
                raise ValueError("No input data provided: need input_df or input_file")

            self.emit_progress("Starting fit...", 10)
            fit_params = {
                "input_csv_file": input_csv,
                "output_dir": self.params.output_dir,
                "output_prefix": self.params.output_prefix,
                "concentrations": self.params.concentrations,
                "conc_units": self.params.conc_units,
                "intensity_scales": self.params.intensity_scales,
                "protein_conc": self.params.protein_conc,
                "alpha": self.params.alpha,
                "observables": self.params.observables or ["csp", "intensity"],
                "intensity_value": self.params.intensity_value,
                "n_bootstrap": self.params.n_bootstrap,
                **(self.params.thresholds or {}),
            }

            def progress_callback(completed, total, residue_name, message):
                if self._cancelled:
                    return
                pct = int((completed / total) * 85) + 10
                self.emit_progress(f"Fitting {residue_name} ({completed}/{total})", pct)

            result = kd_fit.run_kd_analysis_with_params(fit_params, progress_callback=progress_callback)

            if not self._cancelled:
                self.emit_progress("Fitting complete!", 100)
                result["analysis"] = "Kd_titration"
                self.finished.emit(result)
        except Exception as e:
            self.emit_error(e)
        finally:
            if temp_csv_path and os.path.exists(temp_csv_path):
                try:
                    os.remove(temp_csv_path)
                except OSError:
                    pass


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
                "time_units": self.params.time_units,
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


# SPECTRAL DENSITY WORKER

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

            rNH_meters = self.params.r_nh * 1e-10  # Angstroms to meters
            csaN_units = self.params.csa * 1e-6   # ppm to dimensionless

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

                analyzer = DualFieldSpectralDensityAnalysis(
                    field1_freq=self.params.field1_freq,
                    field2_freq=self.params.field2_freq,
                    rNH=rNH_meters,
                    csaN=csaN_units
                )

                self.emit_progress("Running dual-field analysis...", 50)

                results_df = analyzer.analyze_dual_field_csv(
                    csv_file1=self.params.input_file,
                    csv_file2=self.params.input_file2,
                    use_monte_carlo_errors=self.params.monte_carlo,
                    n_monte_carlo=self.params.n_samples,
                    use_multiprocessing=True
                )

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

                analyzer = ReducedSpectralDensityAnalysis(
                    spectrometer_frequency=self.params.field1_freq,
                    rNH=rNH_meters,
                    csaN=csaN_units
                )

                self.emit_progress("Running single-field analysis...", 50)

                results_df = analyzer.analyze_csv(
                    csv_file=self.params.input_file,
                    use_monte_carlo_errors=self.params.monte_carlo,
                    n_monte_carlo=self.params.n_samples,
                    use_multiprocessing=True
                )

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


# INTEGRATED ANALYSIS WORKER

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

    # Units of every relaxation series' delays (the page states one for all four)
    time_units: str = "ms"


def pipeline_units(params: IntegratedAnalysisParams) -> Dict[str, str]:
    """The pipeline's per-series delay units, each the page's one choice."""
    return {f"field{f}_t{t}_units": params.time_units for f in (1, 2) for t in (1, 2)}


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

            from dynamiXs_integrated.integrated_analysis import (
                IntegratedAnalysisPipeline,
                IntegratedAnalysisParameters
            )

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
            for name, units in pipeline_units(self.params).items():
                setattr(pipeline_params, name, units)

            # Output settings
            pipeline_params.output_dir = self.params.output_dir
            pipeline_params.output_prefix = self.params.results_prefix
            pipeline_params.json_folder = self.params.json_folder

            def progress_cb(msg, pct=None):
                if not self._cancelled:
                    self.emit_progress(msg, pct)

            pipeline = IntegratedAnalysisPipeline(pipeline_params, progress_callback=progress_cb)
            result = pipeline.run_complete_analysis()

            if not self._cancelled:
                self.emit_progress("Integrated analysis complete!", 100)
                self.finished.emit(result)

        except Exception as e:
            self.emit_error(e)
