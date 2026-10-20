#!/usr/bin/env python3
"""
Reduced Spectral Density Mapping Analysis Script with Multiprocessing Support

This script imports R1, R1err, R2, R2err, hetNOE, and hetNOEerr from a CSV file
and performs reduced spectral density calculations including J(0), J(wN), J(wH),
S2, Rex, and te parameters with proper error propagation.

omega_h_factor 0.87 evaluates J at 0.87 wH instead of wH - J(0.87wH), for a more
accurate treatment of cross-correlation effects and dipolar-CSA interference - and
names its column JwH_087.

MULTICORE VERSION: Modified to utilize multiple CPU cores for parallel processing.
"""

#####


#####

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

# Multiprocessing imports
import multiprocessing
from multiprocessing import cpu_count
from dynamiXs_density_functions.model_free_search import (
    SINGLE_FIELD_STARTS, S2_BOUNDS, TE_BOUNDS, REX_BOUNDS, TE_SCALE, best_least_squares)

from dynamiXs_density_functions.spectral_density_physics import (
    spectral_density, omega_h, omega_n, dipolar_constant, csa_constant,
    reduced_spectral_densities, reduced_spectral_density_errors, relaxation_deviations,
    overall_correlation_time)

# The J(wH) column for each proton-frequency factor J(wH) can be evaluated at
JWH_COLUMN = {1.0: 'JwH', 0.87: 'JwH_087'}


class ReducedSpectralDensityAnalysis:
    """
    Class for performing reduced spectral density mapping analysis
    """
    
    def __init__(self, spectrometer_frequency=600.0, rNH=1.015e-10, csaN=-160.0e-6,
                 omega_h_factor=1.0):
        """
        Initialize the analysis with experimental parameters
        
        Parameters:
        -----------
        spectrometer_frequency : float
            Spectrometer frequency in MHz (default: 600.0)
        rNH : float
            N-H bond length in meters (default: 1.015e-10)
        csaN : float
            15N CSA in ppm converted to frequency units (default: -160.0e-6)
        omega_h_factor : float
            J(wH) is evaluated at omega_h_factor * wH: 1.0, or 0.87 for J(0.87wH)
        """
        if omega_h_factor not in JWH_COLUMN:
            raise ValueError(f"omega_h_factor must be one of {sorted(JWH_COLUMN)}, "
                             f"not {omega_h_factor}")
        self.spectrometer_frequency = spectrometer_frequency
        self.rNH = rNH
        self.csaN = csaN
        
        self.omegaH = self._calculate_omegaH()
        self.omegaN = self._calculate_omegaN()
        
        self.d_factor = self._calculate_d_factor()
        self.c_factor = self._calculate_c_factor()

        self.omega_h_factor = omega_h_factor
        self.omegaH_eff = self.omegaH * omega_h_factor
        self.jwh_column = JWH_COLUMN[omega_h_factor]
        self.jwh_label = 'J(ωH)' if omega_h_factor == 1.0 else 'J(0.87ωH)'
        # Appended to plot titles and summaries of the 0.87 analysis
        self.variant_note = '' if omega_h_factor == 1.0 else ' (J(0.87ωH))'
        self.file_suffix = '' if omega_h_factor == 1.0 else '_087'
        
    def _calculate_omegaH(self):
        """Calculate 1H Larmor frequency in rad/s"""
        return omega_h(self.spectrometer_frequency)
    
    def _calculate_omegaN(self):
        """Calculate 15N Larmor frequency in rad/s"""
        return omega_n(self.spectrometer_frequency)
    
    def _calculate_d_factor(self):
        """Dipolar coupling constant d² (Farrow et al. 1995)"""
        return dipolar_constant(self.rNH)
    
    def _calculate_c_factor(self):
        """Calculate CSA constant"""
        return csa_constant(self.omegaN, self.csaN)
    
    def calculate_isotropic_spectral_density(self, s2, tc, te, omega):
        """J(w) of the Lipari-Szabo model-free form"""
        return spectral_density(s2, tc, te, omega)
    
    def estimate_overall_correlation_time(self, r1_array, r2_array):
        """Overall correlation time from R1/R2 data by the Fushman method"""
        return overall_correlation_time(r1_array, r2_array, abs(self.omegaN))

    def _single_field_residuals(self, s2, te, rex, t1, t2, noe, tc_fixed):
        """Relative deviations of the predicted T1, T2 and (when measured) NOE from the
        measured ones, at fixed tc - the fit's objective, summed in squares."""
        return relaxation_deviations(s2, te, rex, t1, t2, noe, tc_fixed,
                                     self.d_factor, self.c_factor,
                                     self.omegaH_eff, self.omegaH, abs(self.omegaN))

    def _fit_single_dataset(self, r1, r2, noe, tc_fixed, niter):
        """Fit S2, te and Rex for one residue at fixed tc.

        Bounded least squares from a grid of starts (model_free_search): a parameter on a
        bound can leave it, and the answer does not depend on a random seed. `niter` caps
        the residual evaluations of each start. The fit carries no uncertainty of its own -
        the residuals are relative deviations, not weighted by measurement errors - so the
        errors are NaN here and come from fit_model_free_individual's Monte-Carlo refits.
        """
        t1 = 1.0 / r1
        t2 = 1.0 / r2
        fit = best_least_squares(
            lambda p: self._single_field_residuals(p[0], p[1], p[2], t1, t2, noe, tc_fixed),
            SINGLE_FIELD_STARTS,
            lower=(S2_BOUNDS[0], TE_BOUNDS[0], REX_BOUNDS[0]),
            upper=(S2_BOUNDS[1], TE_BOUNDS[1], REX_BOUNDS[1]),
            scale=(1.0, TE_SCALE, 1.0), max_nfev=niter)
        s2_best, te_best, rex_best = fit['params']
        return {
            'S2': s2_best,
            'tc': tc_fixed * 1e9,  # Convert to ns
            'te': te_best * 1e12,  # Convert to ps
            'Rex': rex_best,
            'S2_err': np.nan,
            'te_err': np.nan,
            'Rex_err': np.nan,
            'fit_success': True,
            'chi2': fit['chi2'],
            'iterations': fit['nfev'],
            'at_iteration_cap': not fit['converged']
        }

    def fit_model_free_individual(self, r1, r2, noe, tc_fixed,
                                r1_err=None, r2_err=None, noe_err=None,
                                n_monte_carlo=100, niter=10_000):
        """
        Fit individual residue model-free parameters using fixed tc
        Now includes proper error propagation via Monte Carlo sampling
        
        Parameters:
        -----------
        r1, r2, noe : float
            Experimental relaxation parameters
        tc_fixed : float
            Fixed overall correlation time (seconds)
        r1_err, r2_err, noe_err : float, optional
            Experimental uncertainties. If provided, enables Monte Carlo error propagation
        n_monte_carlo : int
            Number of Monte Carlo samples for error propagation (default: 100)
        niter : int
            Number of iterations per fit (default: 10_000)
            
        Returns:
        --------
        dict : Fitted parameters (S2, te, Rex) with errors
        """
        # First, fit the actual experimental data
        best_fit = self._fit_single_dataset(r1, r2, noe, tc_fixed, niter)
        
        # If no errors provided or fit failed, return original behavior
        if (r1_err is None or r2_err is None or noe_err is None or 
            not best_fit['fit_success']):
            return best_fit
        
        # Monte Carlo error propagation
        print(f"  Running Monte Carlo error analysis ({n_monte_carlo} samples)...", end='', flush=True)
        mc_results = {'S2': [], 'te': [], 'Rex': [], 'chi2': []}
        
        # Use reduced iterations for MC samples since we start from a good guess
        mc_niter = max(1000, niter // 10)
        
        for i in range(n_monte_carlo):
            # Sample within experimental uncertainties
            r1_sample = np.random.normal(r1, r1_err)
            r2_sample = np.random.normal(r2, r2_err) 
            noe_sample = np.random.normal(noe, noe_err)
            
            r1_sample = max(0.1, r1_sample)  # R1 must be positive
            r2_sample = max(r1_sample * 1.1, r2_sample)  # R2 > R1
            # NOE can be negative but typically between -0.5 and 1.0
            noe_sample = np.clip(noe_sample, -0.5, 1.0)
            
            mc_fit = self._fit_single_dataset(r1_sample, r2_sample, noe_sample, 
                                            tc_fixed, mc_niter)
            
            if mc_fit['fit_success']:
                mc_results['S2'].append(mc_fit['S2'])
                mc_results['te'].append(mc_fit['te']) 
                mc_results['Rex'].append(mc_fit['Rex'])
                mc_results['chi2'].append(mc_fit['chi2'])
        
        print(" done")
        
        n_successful = len(mc_results['S2'])
        if n_successful < 10:  # Need minimum samples for reliable errors
            print(f"    Warning: Only {n_successful}/{n_monte_carlo} MC fits succeeded")
            # Fall back to ensemble-based errors
            return best_fit
        
        best_fit['S2_err'] = np.std(mc_results['S2'])
        best_fit['te_err'] = np.std(mc_results['te'])
        best_fit['Rex_err'] = np.std(mc_results['Rex'])
        
        best_fit['S2_95CI'] = np.percentile(mc_results['S2'], [2.5, 97.5])
        best_fit['te_95CI'] = np.percentile(mc_results['te'], [2.5, 97.5])
        best_fit['Rex_95CI'] = np.percentile(mc_results['Rex'], [2.5, 97.5])
        
        best_fit['mc_success_rate'] = n_successful / n_monte_carlo
        
        return best_fit

    
    def propagate_errors(self, r1, r1_err, r2, r2_err, noe, noe_err):
        """
        Spectral densities and their first-order errors

        Parameters:
        -----------
        r1, r1_err : float
            R1 value and error
        r2, r2_err : float
            R2 value and error
        noe, noe_err : float
            NOE value and error

        Returns:
        --------
        dict : Spectral density values and errors
        """
        j0, jwn, jwh = reduced_spectral_densities(r1, r2, noe, self.d_factor, self.c_factor)
        j0_err, jwn_err, jwh_err = reduced_spectral_density_errors(
            r1, r1_err, r2, r2_err, noe, noe_err, self.d_factor, self.c_factor)
        return {
            'J0': j0, 'J0_err': j0_err,
            'JwN': jwn, 'JwN_err': jwn_err,
            self.jwh_column: jwh, f'{self.jwh_column}_err': jwh_err
        }
    
    def analyze_csv(self, csv_file, residue_col='Residue', use_monte_carlo_errors=True,
                    n_monte_carlo=100, use_multiprocessing=True, n_processes=None):
        """
        Analyze relaxation data from CSV file with multiprocessing support
        
        Parameters:
        -----------
        csv_file : str
            Path to CSV file containing columns:
            R1, R1err, R2, R2err, hetNOE, hetNOEerr
        residue_col : str
            Name of residue identifier column
        use_monte_carlo_errors : bool
            Whether to use Monte Carlo error propagation (default: True)
        n_monte_carlo : int
            Number of Monte Carlo samples for error analysis (default: 100)
        use_multiprocessing : bool
            Whether to use multiprocessing for residue analysis (default: True)
        n_processes : int, optional
            Number of processes for parallel processing (default: 80% of CPU cores)
            
        Returns:
        --------
        pandas.DataFrame : Analysis results
        """
        data = pd.read_csv(csv_file)
        
        required_cols = ['R1', 'R1err', 'R2', 'R2err', 'hetNOE', 'hetNOEerr']
        missing_cols = [col for col in required_cols if col not in data.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        if n_processes is None:
            n_processes = max(1, int(cpu_count() * 0.8))
        
        print(f"\nProcessing {len(data)} residues from {csv_file}...")
        if self.omega_h_factor != 1.0:
            print(f"Using J(0.87ωH) approach with factor = {self.omega_h_factor}")
        if use_monte_carlo_errors:
            print(f"Using Monte Carlo error propagation with {n_monte_carlo} samples per residue")
        if use_multiprocessing:
            print(f"Using multiprocessing with {n_processes} processes")
        
        # First pass: collect valid R1, R2 data for overall tc estimation
        valid_r1 = []
        valid_r2 = []
        
        for idx, row in data.iterrows():
            r1, r2 = row['R1'], row['R2']
            if not (pd.isna(r1) or pd.isna(r2) or r1 <= 0 or r2 <= 0):
                valid_r1.append(r1)
                valid_r2.append(r2)
        
        if len(valid_r1) < 3:
            print("ERROR: Not enough valid R1/R2 data to estimate overall correlation time!")
            return pd.DataFrame()
        
        tc_overall = self.estimate_overall_correlation_time(np.array(valid_r1), np.array(valid_r2))
        print(f"Estimated overall correlation time: {tc_overall*1e9:.1f} ns")
        
        # Process residues: choose between parallel and sequential processing
        if use_multiprocessing and len(data) > 1:
            return self._process_residues_parallel(data, tc_overall, residue_col, 
                                                 use_monte_carlo_errors, n_monte_carlo, n_processes)
        else:
            return self._process_residues_sequential(data, tc_overall, residue_col, 
                                                   use_monte_carlo_errors, n_monte_carlo)
    
    def _process_residues_sequential(self, data, tc_overall, residue_col, 
                                   use_monte_carlo_errors, n_monte_carlo):
        """
        Process residues sequentially (original behavior)
        """
        results = []
        
        # Counters for data quality reporting
        total_residues = len(data)
        skipped_nan = 0
        skipped_zero_negative = 0
        processed = 0
        
        for idx, row in data.iterrows():
            r1 = row['R1']
            r1_err = row['R1err'] 
            r2 = row['R2']
            r2_err = row['R2err']
            noe = row['hetNOE']
            noe_err = row['hetNOEerr']
            
            residue_id = row[residue_col] if residue_col in data.columns else f"Index_{idx}"
            
            if any(pd.isna([r1, r1_err, r2, r2_err, noe, noe_err])):
                skipped_nan += 1
                print(f"  Skipping {residue_id}: Contains NaN values")
                continue
            
            # Skip if any core measurements are zero or negative (indicating missing data)
            if (r1 <= 0 or r2 <= 0 or r1_err <= 0 or r2_err <= 0 or 
                noe == 0 or noe_err <= 0):
                skipped_zero_negative += 1
                reason = []
                if r1 <= 0: reason.append(f"R1={r1}")
                if r2 <= 0: reason.append(f"R2={r2}")
                if r1_err <= 0: reason.append(f"R1err={r1_err}")
                if r2_err <= 0: reason.append(f"R2err={r2_err}")
                if noe == 0: reason.append(f"hetNOE={noe}")
                if noe_err <= 0: reason.append(f"hetNOEerr={noe_err}")
                print(f"  Skipping {residue_id}: Invalid values - {', '.join(reason)}")
                continue
            
            # Additional sanity checks for physical reasonableness
            warnings = []
            if r1 > 10: warnings.append(f"unusually high R1={r1:.2f}")
            if r2 > 100: warnings.append(f"unusually high R2={r2:.2f}")
            if noe < 0.2 or noe > 1.2: warnings.append(f"unusual NOE={noe:.3f}")
            if r1_err/r1 > 0.5: warnings.append(f"large R1 error ({r1_err/r1*100:.1f}%)")
            if r2_err/r2 > 0.5: warnings.append(f"large R2 error ({r2_err/r2*100:.1f}%)")
            
            if warnings:
                print(f"  Warning for {residue_id}: {', '.join(warnings)}")
                
            j_results = self.propagate_errors(r1, r1_err, r2, r2_err, noe, noe_err)
            
            # Fit model-free parameters using fixed overall tc
            if use_monte_carlo_errors:
                # New behavior: include experimental errors
                mf_results = self.fit_model_free_individual(
                    r1, r2, noe, tc_overall, 
                    r1_err=r1_err, r2_err=r2_err, noe_err=noe_err,
                    n_monte_carlo=n_monte_carlo
                )
            else:
                # Original behavior: no experimental error propagation
                mf_results = self.fit_model_free_individual(r1, r2, noe, tc_overall)
            
            # Compile results
            result_row = {
                'Index': idx,
                'R1': r1, 'R1_err': r1_err,
                'R2': r2, 'R2_err': r2_err, 
                'hetNOE': noe, 'hetNOE_err': noe_err,
                **j_results,
                **mf_results
            }
            
            if residue_col in data.columns:
                result_row['Residue'] = row[residue_col]
                
            results.append(result_row)
            processed += 1
        
        # Report data quality summary
        print(f"\nData Quality Summary{self.variant_note}:")
        print(f"  Total residues in file: {total_residues}")
        print(f"  Processed successfully: {processed}")
        print(f"  Skipped (NaN values): {skipped_nan}")
        print(f"  Skipped (zero/negative values): {skipped_zero_negative}")
        print(f"  Success rate: {processed/total_residues*100:.1f}%")
        
        if processed == 0:
            print("  WARNING: No valid data found to process!")
            return pd.DataFrame()
        
        return pd.DataFrame(results)
    
    def _process_residues_parallel(self, data, tc_overall, residue_col, 
                                 use_monte_carlo_errors, n_monte_carlo, n_processes):
        """
        Process residues in parallel using multiprocessing
        """
        # Prepare analyzer state for serialization
        analyzer_state = {
            'spectrometer_frequency': self.spectrometer_frequency,
            'rNH': self.rNH,
            'csaN': self.csaN,
            'omega_h_factor': self.omega_h_factor
        }
        
        worker_args = []
        for idx, row in data.iterrows():
            worker_args.append((
                idx, row.to_dict(), analyzer_state, tc_overall, 
                use_monte_carlo_errors, n_monte_carlo, residue_col
            ))
        
        print("Starting parallel processing...")
        with multiprocessing.get_context('spawn').Pool(processes=n_processes) as pool:
            worker_results = pool.map(_worker_process_residue_single, worker_args)
        
        results = []
        skipped_nan = 0
        skipped_invalid = 0
        errors = 0
        
        for worker_result in worker_results:
            if worker_result['status'] == 'success':
                results.append(worker_result['result'])
                if worker_result['warnings']:
                    print(f"  Warning for {worker_result['residue_id']}: {', '.join(worker_result['warnings'])}")
            elif worker_result['status'] == 'skipped_nan':
                skipped_nan += 1
                print(f"  Skipping {worker_result['residue_id']}: Contains NaN values")
            elif worker_result['status'] == 'skipped_invalid':
                skipped_invalid += 1
                print(f"  Skipping {worker_result['residue_id']}: Invalid values")
            elif worker_result['status'] == 'error':
                errors += 1
                print(f"  Error processing {worker_result['residue_id']}: {worker_result['error']}")
        
        # Report data quality summary
        total_residues = len(data)
        processed = len(results)
        
        print(f"\nData Quality Summary{self.variant_note}:")
        print(f"  Total residues in file: {total_residues}")
        print(f"  Processed successfully: {processed}")
        print(f"  Skipped (NaN values): {skipped_nan}")
        print(f"  Skipped (invalid values): {skipped_invalid}")
        print(f"  Processing errors: {errors}")
        print(f"  Success rate: {processed/total_residues*100:.1f}%")
        
        if processed == 0:
            print("  WARNING: No valid data found to process!")
            return pd.DataFrame()
        
        return pd.DataFrame(results)
    
    def plot_results(self, results_df, save_plots=True, plot_filename=None):
        """
        Generate plots of the analysis results
        
        Parameters:
        -----------
        results_df : pandas.DataFrame
            Results from analyze_csv()
        save_plots : bool
            Whether to save plots to files
        """
        fig = Figure(figsize=(15, 10))
        axes = fig.subplots(2, 3)
        fig.suptitle(f'Reduced Spectral Density Analysis Results{self.variant_note}', fontsize=16)
        
        success_mask = results_df['fit_success'] == True
        good_data = results_df[success_mask]
        
        if len(good_data) == 0:
            print("No successful fits to plot!")
            return
        
        residues = good_data.index if 'Residue' not in good_data.columns else good_data['Residue']
        
        # J(0)
        axes[0,0].errorbar(residues, good_data['J0'], yerr=good_data['J0_err'], 
                          fmt='o', capsize=3)
        axes[0,0].set_title('J(0)')
        axes[0,0].set_ylabel('J(0) (ns/rad²)')
        
        # J(ωN)  
        axes[0,1].errorbar(residues, good_data['JwN'], yerr=good_data['JwN_err'],
                          fmt='o', capsize=3, color='orange')
        axes[0,1].set_title('J(ωN)')
        axes[0,1].set_ylabel('J(ωN) (ns/rad²)')
        
        # J(ωH)
        axes[0,2].errorbar(residues, good_data[self.jwh_column],
                          yerr=good_data[f'{self.jwh_column}_err'],
                          fmt='o', capsize=3, color='green') 
        axes[0,2].set_title(self.jwh_label)
        axes[0,2].set_ylabel(f'{self.jwh_label} (ns/rad²)')
        
        # S²
        axes[1,0].errorbar(residues, good_data['S2'], yerr=good_data['S2_err'],
                          fmt='o', capsize=3, color='blue')         
        axes[1,0].set_title(f'Order Parameter S²{self.variant_note}')
        axes[1,0].set_ylabel('S²')
        axes[1,0].set_ylim(0, 1)
        
        # τe  
        axes[1,1].errorbar(residues, good_data['te'], yerr=good_data['te_err'],
                          fmt='o', capsize=3, color='purple')
        axes[1,1].set_title(f'Internal Correlation Time τe{self.variant_note}')
        axes[1,1].set_ylim(0, 100) 
        axes[1,1].set_ylabel('τe (ps)')
        
        # Rex
        axes[1,2].errorbar(residues, good_data['Rex'], yerr=good_data['Rex_err'],
                          fmt='o', capsize=3, color='brown')
        axes[1,2].set_title(f'Chemical Exchange Rex{self.variant_note}')
        axes[1,2].set_ylabel('Rex (s⁻¹)')
        
        for ax in axes.flat:
            ax.set_xlabel('Residue')
            ax.grid(True, alpha=0.3)
            
        fig.tight_layout()

        if save_plots:
            filename = plot_filename or f'ZZ_rsdm_multicore_analysis_results{self.file_suffix}.pdf'
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Plots saved as '{filename}'")

    def save_detailed_results(self, results_df, filename=None):
        """
        Save results with confidence intervals if Monte Carlo was used
        
        Parameters:
        -----------
        results_df : pandas.DataFrame
            Results from analyze_csv()
        filename : str
            Output filename (default: ZZ_multicore_detailed_results[_087].csv)
        """
        filename = filename or f'ZZ_multicore_detailed_results{self.file_suffix}.csv'
        if 'S2_95CI' in results_df.columns:
            results_df['S2_CI'] = results_df['S2_95CI'].apply(
                lambda x: f"[{x[0]:.3f}, {x[1]:.3f}]" if not pd.isna(x).any() else ""
            )
            results_df['te_CI'] = results_df['te_95CI'].apply(
                lambda x: f"[{x[0]:.1f}, {x[1]:.1f}]" if not pd.isna(x).any() else ""
            )
            results_df['Rex_CI'] = results_df['Rex_95CI'].apply(
                lambda x: f"[{x[0]:.2f}, {x[1]:.2f}]" if not pd.isna(x).any() else ""
            )
            
            # Drop the raw CI arrays
            results_df = results_df.drop(columns=['S2_95CI', 'te_95CI', 'Rex_95CI'])
        
        results_df.to_csv(filename, index=False)
        print(f"Detailed results{self.variant_note} saved to '{filename}'")


# Worker functions for multiprocessing (must be defined at module level)
def _worker_process_residue_single(args):
    """
    Worker function for multiprocessing - processes a single residue
    This function is defined at module level to be pickle-able by multiprocessing
    """
    (idx, residue_data, analyzer_state, tc_overall, use_monte_carlo_errors, 
     n_monte_carlo, residue_col) = args
    
    try:
        # Reconstruct analyzer object from state
        analyzer = ReducedSpectralDensityAnalysis(
            spectrometer_frequency=analyzer_state['spectrometer_frequency'],
            rNH=analyzer_state['rNH'],
            csaN=analyzer_state['csaN'],
            omega_h_factor=analyzer_state['omega_h_factor']
        )
        
        r1 = residue_data['R1']
        r1_err = residue_data['R1err']
        r2 = residue_data['R2']
        r2_err = residue_data['R2err']
        noe = residue_data['hetNOE']
        noe_err = residue_data['hetNOEerr']
        
        residue_id = residue_data[residue_col] if residue_col in residue_data else f"Index_{idx}"
        
        if any(pd.isna([r1, r1_err, r2, r2_err, noe, noe_err])):
            return {
                'status': 'skipped_nan',
                'residue_id': residue_id,
                'result': None,
                'warnings': [],
                'error': None
            }
        
        if (r1 <= 0 or r2 <= 0 or r1_err <= 0 or r2_err <= 0 or 
            noe == 0 or noe_err <= 0):
            return {
                'status': 'skipped_invalid',
                'residue_id': residue_id,
                'result': None,
                'warnings': [],
                'error': None
            }
        
        # Physical reasonableness checks
        warnings = []
        if r1 > 10: warnings.append(f"unusually high R1={r1:.2f}")
        if r2 > 100: warnings.append(f"unusually high R2={r2:.2f}")
        if noe < 0.2 or noe > 1.2: warnings.append(f"unusual NOE={noe:.3f}")
        if r1_err/r1 > 0.5: warnings.append(f"large R1 error ({r1_err/r1*100:.1f}%)")
        if r2_err/r2 > 0.5: warnings.append(f"large R2 error ({r2_err/r2*100:.1f}%)")
        
        j_results = analyzer.propagate_errors(r1, r1_err, r2, r2_err, noe, noe_err)
        
        if use_monte_carlo_errors:
            # Use sequential Monte Carlo in worker processes to avoid nested multiprocessing
            mf_results = analyzer.fit_model_free_individual(
                r1, r2, noe, tc_overall, 
                r1_err=r1_err, r2_err=r2_err, noe_err=noe_err,
                n_monte_carlo=n_monte_carlo
            )
        else:
            mf_results = analyzer.fit_model_free_individual(r1, r2, noe, tc_overall)
        
        # Compile results
        result_row = {
            'Index': idx,
            'R1': r1, 'R1_err': r1_err,
            'R2': r2, 'R2_err': r2_err, 
            'hetNOE': noe, 'hetNOE_err': noe_err,
            **j_results,
            **mf_results
        }
        
        if residue_col in residue_data:
            result_row['Residue'] = residue_data[residue_col]
        
        return {
            'status': 'success',
            'residue_id': residue_id,
            'result': result_row,
            'warnings': warnings,
            'error': None
        }
        
    except Exception as e:
        return {
            'status': 'error',
            'residue_id': residue_id if 'residue_id' in locals() else f"Index_{idx}",
            'result': None,
            'warnings': [],
            'error': str(e)
        }




def main():
    """
    Example usage of the ReducedSpectralDensityAnalysis class with multiprocessing
    """
    analyzer = ReducedSpectralDensityAnalysis(
        spectrometer_frequency=700.0,  # MHz
        rNH=1.023e-10,                   # meters
        csaN=-160.0e-6                   # ppm in frequency units
    )
    
    # Analyze CSV file with multiprocessing support
    try:
        print("Running analysis with multiprocessing and Monte Carlo error propagation...")
        results = analyzer.analyze_csv('data_in_700.csv',
                                     use_monte_carlo_errors=True,
                                     n_monte_carlo=100,
                                     use_multiprocessing=True,
                                     n_processes=None)  # Use 80% of available cores
        
        if len(results) == 0:
            print("No data to analyze after filtering!")
            return
        
        print("\nReduced Spectral Density Analysis Results (Multiprocessing):")
        print("=" * 70)
        
        if 'Residue' in results.columns:
            display_cols = ['Residue', 'J0', 'JwN', 'JwH', 'S2', 'tc', 'te', 'Rex']
        else:
            display_cols = ['Index', 'J0', 'JwN', 'JwH', 'S2', 'tc', 'te', 'Rex']
        
        print(results[display_cols].round(4))
        
        results.to_csv('ZZ_rsdm_multicore.csv', index=False)
        
        analyzer.save_detailed_results(results)
        
        analyzer.plot_results(results)
        
        success_rate = results['fit_success'].mean() * 100
        print("\nModel-Free Fitting Results:")
        print(f"  Fit success rate: {success_rate:.1f}%")
        print(f"  Multiprocessing enabled: {success_rate > 0}")
        
        if success_rate > 0:
            successful = results[results['fit_success']]
            print(f"  Successfully fitted residues: {len(successful)}")
            print(f"  Average S²: {successful['S2'].mean():.3f} ± {successful['S2'].std():.3f}")
            print(f"  Average τc: {successful['tc'].mean():.1f} ± {successful['tc'].std():.1f} ns")
            print(f"  Average τe: {successful['te'].mean():.1f} ± {successful['te'].std():.1f} ps")
            print(f"  Average Rex: {successful['Rex'].mean():.1f} ± {successful['Rex'].std():.1f} s⁻¹")
            
            # Report Monte Carlo success if used
            if 'mc_success_rate' in successful.columns:
                avg_mc_success = successful['mc_success_rate'].mean()
                print(f"  Average Monte Carlo success rate: {avg_mc_success*100:.1f}%")
        else:
            print("  No successful model-free fits obtained!")
            
    except FileNotFoundError:
        print("Please provide a CSV file named 'data_in.csv' with columns:")
        print("R1, R1err, R2, R2err, hetNOE, hetNOEerr, and optionally Residue")
        
        example_data = pd.DataFrame({
            'Residue': [f'A{i}' for i in range(1, 11)],
            'R1': np.random.normal(1.5, 0.2, 10),
            'R1err': np.random.normal(0.05, 0.01, 10),
            'R2': np.random.normal(15.0, 2.0, 10), 
            'R2err': np.random.normal(0.5, 0.1, 10),
            'hetNOE': np.random.normal(0.8, 0.1, 10),
            'hetNOEerr': np.random.normal(0.03, 0.005, 10)
        })
        example_data.to_csv('example_relaxation_data.csv', index=False)
        print("\nExample data file 'example_relaxation_data.csv' created.")
        print("You can test the script with this example data.")


if __name__ == "__main__":
    main()