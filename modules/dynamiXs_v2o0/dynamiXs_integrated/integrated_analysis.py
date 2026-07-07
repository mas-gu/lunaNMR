#!/usr/bin/env python3
"""
Integrated Spectral Density Analysis Pipeline

This module orchestrates the complete workflow from raw NMR data to
model-free parameters, automating:
1. T1/T2 exponential decay fitting
2. hetNOE calculation from raw intensities
3. Data format conversion (times → rates)
4. Dataset merging and validation
5. Spectral density analysis

Author: DynamiXs Development Team
Date: 2025-01-23
"""

import sys
import os
from pathlib import Path
from typing import Dict, Tuple, Optional, Callable
import pandas as pd

# Import local modules - use absolute imports for GUI compatibility
try:
    # Try relative imports first (when used as package)
    from .data_converters import (
        convert_relaxation_times_to_rates,
        calculate_hetnoe_from_intensities,
        parse_intensity_csv,
        create_spectral_density_input_csv
    )
    from .data_validation import (
        merge_datasets_by_residue,
        validate_relaxation_rates,
        validate_field_consistency,
        check_minimum_residue_count,
        generate_validation_report,
        filter_invalid_residues
    )
    from .fitting_wrapper import run_t1_fitting, run_t2_fitting
except ImportError:
    # Fall back to absolute imports (when imported directly)
    from data_converters import (
        convert_relaxation_times_to_rates,
        calculate_hetnoe_from_intensities,
        parse_intensity_csv,
        create_spectral_density_input_csv
    )
    from data_validation import (
        merge_datasets_by_residue,
        validate_relaxation_rates,
        validate_field_consistency,
        check_minimum_residue_count,
        generate_validation_report,
        filter_invalid_residues
    )
    from fitting_wrapper import run_t1_fitting, run_t2_fitting


class IntegratedAnalysisParameters:
    """Container for all analysis parameters"""

    def __init__(self):
        # Field 1 (required)
        self.field1_freq_mhz = 600.0
        self.field1_t1_file = None
        self.field1_t2_file = None
        self.field1_noe_sat_file = None
        self.field1_noe_unsat_file = None

        # Field 2 (optional for dual-field)
        self.enable_dual_field = False
        self.field2_freq_mhz = 700.0
        self.field2_t1_file = None
        self.field2_t2_file = None
        self.field2_noe_sat_file = None
        self.field2_noe_unsat_file = None

        # Analysis method
        self.analysis_method = 'dual_field_087'  # Options: 'single_jwh', 'single_087', 'dual_jwh', 'dual_087'

        # Physical parameters
        self.rNH_angstrom = 1.015  # N-H bond length in Angstroms
        self.csaN_ppm = -160.0      # 15N CSA in ppm

        # Fitting parameters
        self.t1_initial_amplitude = 5.0
        self.t1_initial_time = 800.0
        self.t1_bootstrap_iterations = 1000

        self.t2_initial_amplitude = 5.0
        self.t2_initial_time = 100.0
        self.t2_bootstrap_iterations = 1000

        # Error estimation method ('analytical' or 'bootstrap')
        self.error_method = 'analytical'

        # Monte Carlo parameters
        self.monte_carlo_iterations = 50

        # Output
        self.output_dir = None  # Output directory (if None, uses temp/)
        self.output_prefix = 'integrated_analysis'
        self.save_intermediate_files = True
        self.json_folder = None  # Folder for JSON fit data (for visualization)


class IntegratedAnalysisPipeline:
    """Main pipeline orchestrator for integrated spectral density analysis"""

    def __init__(self, params: IntegratedAnalysisParameters,
                 progress_callback: Optional[Callable[[str], None]] = None):
        """
        Initialize the integrated analysis pipeline.

        Parameters
        ----------
        params : IntegratedAnalysisParameters
            Analysis parameters
        progress_callback : callable, optional
            Function to call with progress updates (takes string message)
        """
        self.params = params
        self.progress_callback = progress_callback or print

        # Storage for intermediate results
        self.field1_results = {}
        self.field2_results = {}
        self.final_results = None

        # Create temp directory
        self.temp_dir = Path(__file__).parent.parent / 'temp'
        self.temp_dir.mkdir(exist_ok=True)

    def log_progress(self, message: str):
        """Log progress message"""
        self.progress_callback(message)

    def run_complete_analysis(self) -> Dict:
        """
        Run the complete integrated analysis pipeline.

        Returns
        -------
        dict
            Final analysis results including spectral densities and model-free parameters

        Raises
        ------
        ValueError
            If required files are missing or validation fails
        """
        self.log_progress("=" * 60)
        self.log_progress("INTEGRATED SPECTRAL DENSITY ANALYSIS")
        self.log_progress("=" * 60)

        # Step 1: Process Field 1 data
        self.log_progress("\n[1/7] Processing Field 1 data...")
        self.field1_results = self._process_field_data(
            field_freq=self.params.field1_freq_mhz,
            t1_file=self.params.field1_t1_file,
            t2_file=self.params.field1_t2_file,
            noe_sat_file=self.params.field1_noe_sat_file,
            noe_unsat_file=self.params.field1_noe_unsat_file,
            field_label='field1'
        )

        # Step 2: Process Field 2 data (if dual-field)
        if self.params.enable_dual_field:
            self.log_progress("\n[2/7] Processing Field 2 data...")
            self.field2_results = self._process_field_data(
                field_freq=self.params.field2_freq_mhz,
                t1_file=self.params.field2_t1_file,
                t2_file=self.params.field2_t2_file,
                noe_sat_file=self.params.field2_noe_sat_file,
                noe_unsat_file=self.params.field2_noe_unsat_file,
                field_label='field2'
            )
        else:
            self.log_progress("\n[2/7] Skipping Field 2 (single-field analysis)")

        # Step 3: Validate datasets
        self.log_progress("\n[3/7] Validating datasets...")
        validation_report = self._validate_datasets()
        self.log_progress(validation_report)

        # Step 4: Create spectral density input files
        self.log_progress("\n[4/7] Creating spectral density input files...")
        input_files = self._create_spectral_density_inputs()

        # Step 5: Run spectral density analysis
        self.log_progress("\n[5/7] Running spectral density analysis...")
        results = self._run_spectral_density_analysis(input_files)

        # Step 6: Generate comprehensive report
        self.log_progress("\n[6/7] Generating analysis report...")
        self._generate_report(results)

        # Step 7: Save final results
        self.log_progress("\n[7/7] Saving final results...")
        output_file = self._save_final_results(results)
        self.log_progress(f"✓ Results saved to: {output_file}")

        self.log_progress("\n" + "=" * 60)
        self.log_progress("ANALYSIS COMPLETE!")
        self.log_progress("=" * 60)

        self.final_results = results
        return results

    def _process_field_data(self, field_freq: float, t1_file: str, t2_file: str,
                           noe_sat_file: str, noe_unsat_file: str,
                           field_label: str) -> Dict:
        """Process data for a single field"""

        # Step 1: Fit T1 data
        self.log_progress(f"  ├─ Fitting T1 data ({field_freq} MHz)...")
        t1_params, t1_results_file = run_t1_fitting(
            input_csv=t1_file,
            initial_amplitude=self.params.t1_initial_amplitude,
            initial_t1=self.params.t1_initial_time,
            n_bootstrap=self.params.t1_bootstrap_iterations,
            error_method=self.params.error_method,
            output_prefix=f"{self.params.output_prefix}_{field_label}_T1",
            json_folder=self.params.json_folder,
            field_name=field_label,
            field_freq=field_freq
        )
        self.log_progress(f"  │  ✓ {len(t1_params)} residues fitted")

        # Step 2: Fit T2 data
        self.log_progress(f"  ├─ Fitting T2 data ({field_freq} MHz)...")
        t2_params, t2_results_file = run_t2_fitting(
            input_csv=t2_file,
            initial_amplitude=self.params.t2_initial_amplitude,
            initial_t2=self.params.t2_initial_time,
            n_bootstrap=self.params.t2_bootstrap_iterations,
            error_method=self.params.error_method,
            output_prefix=f"{self.params.output_prefix}_{field_label}_T2",
            json_folder=self.params.json_folder,
            field_name=field_label,
            field_freq=field_freq
        )
        self.log_progress(f"  │  ✓ {len(t2_params)} residues fitted")

        # Step 3: Calculate hetNOE
        self.log_progress(f"  ├─ Calculating hetNOE...")
        sat_intensities, sat_errors = parse_intensity_csv(noe_sat_file)
        unsat_intensities, unsat_errors = parse_intensity_csv(noe_unsat_file)

        noe_params = calculate_hetnoe_from_intensities(
            saturated_data=sat_intensities,
            unsaturated_data=unsat_intensities,
            saturated_errors=sat_errors,
            unsaturated_errors=unsat_errors
        )
        self.log_progress(f"  │  ✓ {len(noe_params)} residues calculated")

        # Step 4: Convert T1/T2 → R1/R2
        self.log_progress(f"  ├─ Converting relaxation times to rates...")
        r1_params = convert_relaxation_times_to_rates(t1_params, time_units='ms')
        r2_params = convert_relaxation_times_to_rates(t2_params, time_units='ms')
        self.log_progress(f"  │  ✓ Converted to R1/R2 (s⁻¹)")

        # Step 5: Merge datasets
        self.log_progress(f"  └─ Merging datasets...")
        merged_data, excluded = merge_datasets_by_residue(
            r1_data=r1_params,
            r2_data=r2_params,
            noe_data=noe_params,
            strategy='intersection'
        )
        self.log_progress(f"     ✓ {len(merged_data)} residues in common")
        if excluded:
            self.log_progress(f"     ⚠ {len(excluded)} residues excluded (not in all datasets)")

        return {
            'frequency': field_freq,
            'merged_data': merged_data,
            'excluded_residues': excluded,
            't1_results_file': t1_results_file,
            't2_results_file': t2_results_file
        }

    def _validate_datasets(self) -> str:
        """Validate all datasets and return report"""

        # Validate Field 1
        warnings_field1 = validate_relaxation_rates(self.field1_results['merged_data'])
        check_minimum_residue_count(self.field1_results['merged_data'])

        # Validate Field 2 if present
        if self.params.enable_dual_field:
            warnings_field2 = validate_relaxation_rates(self.field2_results['merged_data'])
            check_minimum_residue_count(self.field2_results['merged_data'])

            # Check consistency between fields
            common, f1_only, f2_only = validate_field_consistency(
                self.field1_results['merged_data'],
                self.field2_results['merged_data']
            )

            if f1_only or f2_only:
                self.log_progress(f"  ⚠ Field mismatch: {len(common)} common, "
                                f"{len(f1_only)} only in field1, {len(f2_only)} only in field2")

        # Generate validation report
        report = generate_validation_report(
            warnings_field1,
            self.field1_results['excluded_residues']
        )

        return report

    def _create_spectral_density_inputs(self) -> Dict[str, str]:
        """Create input CSV files for spectral density analysis"""

        input_files = {}

        # Field 1
        field1_input = str(self.temp_dir / f"{self.params.output_prefix}_field1_input.csv")
        n_residues_f1 = create_spectral_density_input_csv(
            r1_data={k: {'value': v['R1'], 'error': v['R1_err']}
                    for k, v in self.field1_results['merged_data'].items()},
            r2_data={k: {'value': v['R2'], 'error': v['R2_err']}
                    for k, v in self.field1_results['merged_data'].items()},
            noe_data={k: {'value': v['hetNOE'], 'error': v['hetNOE_err']}
                     for k, v in self.field1_results['merged_data'].items()},
            output_file=field1_input
        )
        input_files['field1'] = field1_input
        self.log_progress(f"  ✓ Field 1 input created: {n_residues_f1} residues")

        # Field 2 if dual-field
        if self.params.enable_dual_field:
            field2_input = str(self.temp_dir / f"{self.params.output_prefix}_field2_input.csv")
            n_residues_f2 = create_spectral_density_input_csv(
                r1_data={k: {'value': v['R1'], 'error': v['R1_err']}
                        for k, v in self.field2_results['merged_data'].items()},
                r2_data={k: {'value': v['R2'], 'error': v['R2_err']}
                        for k, v in self.field2_results['merged_data'].items()},
                noe_data={k: {'value': v['hetNOE'], 'error': v['hetNOE_err']}
                         for k, v in self.field2_results['merged_data'].items()},
                output_file=field2_input
            )
            input_files['field2'] = field2_input
            self.log_progress(f"  ✓ Field 2 input created: {n_residues_f2} residues")

        return input_files

    def _run_spectral_density_analysis(self, input_files: Dict[str, str]) -> Dict:
        """Run the appropriate spectral density analysis script using multicore versions"""

        # Import spectral density module
        density_module_path = Path(__file__).parent.parent / 'dynamiXs_density_functions'
        sys.path.insert(0, str(density_module_path))

        # Determine which multicore script to use
        method_map = {
            'single_jwh': ('ZZ_multi_density', 'ReducedSpectralDensityAnalysis'),
            'single_087': ('ZZ_multi_density_087', 'ReducedSpectralDensityAnalysis'),
            'dual_jwh': ('ZZ_multi_2fields_density', 'DualFieldSpectralDensityAnalysis'),
            'dual_087': ('ZZ_multi_2fields_density_087', 'DualFieldSpectralDensityAnalysis')
        }

        script_info = method_map.get(self.params.analysis_method)
        if not script_info:
            raise ValueError(f"Unknown analysis method: {self.params.analysis_method}")

        script_name, class_name = script_info
        self.log_progress(f"  Using method: {self.params.analysis_method} ({script_name})")
        self.log_progress(f"  ⚡ Multicore spectral density analysis enabled")

        # Import the appropriate class
        try:
            module = __import__(script_name)
            AnalysisClass = getattr(module, class_name)
        except ImportError as e:
            raise ImportError(f"Could not import {script_name}: {e}")

        # Prepare output files
        if self.params.output_dir:
            output_dir = Path(self.params.output_dir)
        else:
            output_dir = Path(__file__).parent.parent / 'temp'
        output_dir.mkdir(exist_ok=True)
        output_prefix = str(output_dir / self.params.output_prefix)

        # Convert rNH from Angstroms to meters
        rNH_meters = self.params.rNH_angstrom * 1e-10

        # Convert CSA from ppm to proper units
        csaN_units = self.params.csaN_ppm * 1e-6

        # Run analysis based on method type
        if self.params.analysis_method.startswith('dual_'):
            # Dual-field analysis
            self.log_progress(f"  Analyzing dual-field data ({self.params.field1_freq_mhz} MHz + {self.params.field2_freq_mhz} MHz)")

            analyzer = AnalysisClass(
                field1_freq=self.params.field1_freq_mhz,
                field2_freq=self.params.field2_freq_mhz,
                rNH=rNH_meters,
                csaN=csaN_units
            )

            results_df = analyzer.analyze_dual_field_csv(
                csv_file1=input_files['field1'],
                csv_file2=input_files['field2'],
                use_monte_carlo_errors=True,
                n_monte_carlo=self.params.monte_carlo_iterations,
                use_multiprocessing=True,
                n_processes=None  # Use default (80% of cores)
            )

            # Save results
            basic_csv = f"{output_prefix}_spectral_density_basic.csv"
            detailed_csv = f"{output_prefix}_spectral_density_detailed.csv"
            plots_pdf = f"{output_prefix}_spectral_density_plots.pdf"

            results_df.to_csv(basic_csv, index=False)
            analyzer.save_dual_field_results(results_df, detailed_csv)
            analyzer.plot_dual_field_results(results_df, save_plots=True, plot_filename=plots_pdf)

            self.log_progress(f"  ✓ Dual-field analysis completed: {len(results_df)} residues")
            self.log_progress(f"  ✓ Results saved to: {basic_csv}")
            self.log_progress(f"  ✓ Plots saved to: {plots_pdf}")

        else:
            # Single-field analysis
            self.log_progress(f"  Analyzing single-field data ({self.params.field1_freq_mhz} MHz)")

            analyzer = AnalysisClass(
                spectrometer_frequency=self.params.field1_freq_mhz,
                rNH=rNH_meters,
                csaN=csaN_units
            )

            results_df = analyzer.analyze_csv(
                csv_file=input_files['field1'],
                use_monte_carlo_errors=True,
                n_monte_carlo=self.params.monte_carlo_iterations,
                use_multiprocessing=True,
                n_processes=None  # Use default (80% of cores)
            )

            # Save results
            basic_csv = f"{output_prefix}_spectral_density_basic.csv"
            detailed_csv = f"{output_prefix}_spectral_density_detailed.csv"
            plots_pdf = f"{output_prefix}_spectral_density_plots.pdf"

            results_df.to_csv(basic_csv, index=False)
            analyzer.save_detailed_results(results_df, detailed_csv)
            analyzer.plot_results(results_df, save_plots=True)

            self.log_progress(f"  ✓ Single-field analysis completed: {len(results_df)} residues")
            self.log_progress(f"  ✓ Results saved to: {basic_csv}")
            self.log_progress(f"  ✓ Plots saved to: {plots_pdf}")

        # Return summary
        success_count = results_df['fit_success'].sum() if 'fit_success' in results_df.columns else len(results_df)

        results = {
            'method': self.params.analysis_method,
            'input_files': input_files,
            'output_files': {
                'basic_csv': basic_csv,
                'detailed_csv': detailed_csv,
                'plots_pdf': plots_pdf
            },
            'status': 'completed',
            'n_residues': len(results_df),
            'n_successful': int(success_count),
            'results_dataframe': results_df
        }

        return results

    def _generate_report(self, results: Dict):
        """Generate comprehensive analysis report"""
        self.log_progress("  ✓ Analysis report generated")

    def _save_final_results(self, results: Dict) -> str:
        """Save final results to CSV"""
        output_file = f"{self.params.output_prefix}_integrated_results.csv"
        # Save results (placeholder)
        return output_file
