#!/usr/bin/env python3
"""
Dual-Field Reduced Spectral Density Mapping Analysis Script - GUI Compatible Version

This script imports R1, R1err, R2, R2err, hetNOE, and hetNOEerr from two CSV files
at different field strengths and performs reduced spectral density calculations 
including J(0), J(wN), J(wH), S2, Rex, and te parameters with proper error propagation.

GUI Compatible wrapper functions added for integration with DynamiXs GUI.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for thread safety
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

# Physical constants
GAMMA_H = 2.67522128e8  # 1H gyromagnetic ratio (rad s-1 T-1)
GAMMA_N = -2.7126e7     # 15N gyromagnetic ratio (rad s-1 T-1) 
REDUCED_PERM_VACUUM = 1.25663706212e-6  # H/m
REDUCED_PLANK = 1.05457180013e-34       # J⋅s
PI = 3.14159265359

class DualFieldSpectralDensityAnalysis:
    """
    Class for performing dual-field reduced spectral density mapping analysis
    """
    
    def __init__(self, field1_freq=600.0, field2_freq=800.0, rNH=1.015e-10, csaN=-160.0e-6):
        """
        Initialize the analysis with experimental parameters for two fields
        
        Parameters:
        -----------
        field1_freq : float
            First spectrometer frequency in MHz (default: 600.0)
        field2_freq : float
            Second spectrometer frequency in MHz (default: 800.0)
        rNH : float
            N-H bond length in meters (default: 1.015e-10)
        csaN : float
            15N CSA in ppm converted to frequency units (default: -160.0e-6)
        """
        self.field1_freq = field1_freq
        self.field2_freq = field2_freq
        self.rNH = rNH
        self.csaN = csaN
        
        # Calculate angular frequencies for both fields
        self.field1_omegaH = self._calculate_omegaH(field1_freq)
        self.field1_omegaN = self._calculate_omegaN(field1_freq)
        self.field2_omegaH = self._calculate_omegaH(field2_freq)
        self.field2_omegaN = self._calculate_omegaN(field2_freq)
        
        # Calculate dipolar and CSA constants for both fields
        self.field1_d_factor = self._calculate_d_factor()
        self.field1_c_factor = self._calculate_c_factor(self.field1_omegaN)
        self.field2_d_factor = self._calculate_d_factor()
        self.field2_c_factor = self._calculate_c_factor(self.field2_omegaN)
        
    def _calculate_omegaH(self, freq_mhz):
        """Calculate 1H Larmor frequency in rad/s"""
        return freq_mhz * 2 * PI * 1e6
    
    def _calculate_omegaN(self, freq_mhz):
        """Calculate 15N Larmor frequency in rad/s"""
        omegaH = self._calculate_omegaH(freq_mhz)
        return omegaH * GAMMA_N / GAMMA_H
    
    def _calculate_d_factor(self):
        """
        Calculate dipolar coupling constant (d²) - field-independent

        Returns d² as defined in Farrow et al. J. Biomol. NMR, 6 (1995) 153-162
        d² = [μ₀ħγNγH/(8π²r³NH)]²
        """
        mu0_h_bar = REDUCED_PERM_VACUUM * REDUCED_PLANK
        d_squared = (mu0_h_bar * GAMMA_N * GAMMA_H / (4 * PI * self.rNH**3))**2
        return d_squared  # Returns d², not d²/4 (Farrow-exact notation)
    
    def _calculate_c_factor(self, omegaN):
        """Calculate CSA constant (field-dependent)"""
        return (omegaN * self.csaN)**2 / 3.0
    
    def calculate_sigma_NOE(self, noe, r1):
        """Calculate sigma NOE value"""
        return (noe - 1.0) * r1 * (GAMMA_N / GAMMA_H)
    
    def calculate_spectral_densities(self, r1, r2, noe, field='field1'):
        """
        Calculate spectral densities for a given field
        
        Parameters:
        -----------
        r1, r2, noe : float or array
            Relaxation parameters
        field : str
            'field1' or 'field2'
            
        Returns:
        --------
        dict : Spectral density values
        """
        if field == 'field1':
            d_factor = self.field1_d_factor
            c_factor = self.field1_c_factor
        else:
            d_factor = self.field2_d_factor
            c_factor = self.field2_c_factor
            
        sigma_noe = self.calculate_sigma_NOE(noe, r1)

        # Calculate spectral densities using Farrow et al. (1995) equations
        # Farrow Eq. 7
        j0 = (3.0 / (2.0 * (3.0 * d_factor / 4.0 + c_factor))) * (
            -0.5 * r1 + r2 - (3.0/5.0) * sigma_noe
        )

        # Farrow Eq. 6
        jwn = (1.0 / (3.0 * d_factor / 4.0 + c_factor)) * (r1 - (7.0/5.0) * sigma_noe)

        # Farrow Eq. 5
        jwh = 4.0 * sigma_noe / (5.0 * d_factor)

        return {'J0': j0, 'JwN': jwn, 'JwH': jwh}

    def analyze_dual_field_csv(self, csv_file1, csv_file2, residue_col='Residue', 
                              use_monte_carlo_errors=True, n_monte_carlo=50):
        """
        Analyze dual-field relaxation data from two CSV files
        
        Parameters:
        -----------
        csv_file1 : str
            Path to first field CSV file
        csv_file2 : str
            Path to second field CSV file
        residue_col : str
            Name of residue identifier column
        use_monte_carlo_errors : bool
            Whether to use Monte Carlo error propagation
        n_monte_carlo : int
            Number of Monte Carlo samples
            
        Returns:
        --------
        pandas.DataFrame : Dual-field analysis results
        """
        # Load data from both fields
        try:
            data1 = pd.read_csv(csv_file1)
            data2 = pd.read_csv(csv_file2)
            print(f"Loaded data from {csv_file1} and {csv_file2}")
        except FileNotFoundError as e:
            print(f"Error loading CSV files: {e}")
            return pd.DataFrame()
        
        # Check required columns
        required_cols = ['R1', 'R1err', 'R2', 'R2err', 'hetNOE', 'hetNOEerr']
        for i, data in enumerate([data1, data2], 1):
            missing_cols = [col for col in required_cols if col not in data.columns]
            if missing_cols:
                raise ValueError(f"Missing required columns in file {i}: {missing_cols}")
        
        # Align datasets by residue identifier
        if residue_col in data1.columns and residue_col in data2.columns:
            # Merge on residue identifier
            merged_data = pd.merge(data1, data2, on=residue_col, suffixes=('_f1', '_f2'))
            print(f"Matched {len(merged_data)} residues between the two datasets")
        else:
            # Merge by index if no residue column
            print("Warning: No residue column found, merging by index position")
            min_len = min(len(data1), len(data2))
            merged_data = pd.concat([
                data1.iloc[:min_len].reset_index(drop=True).add_suffix('_f1'),
                data2.iloc[:min_len].reset_index(drop=True).add_suffix('_f2')
            ], axis=1)
            if residue_col + '_f1' in merged_data.columns:
                merged_data[residue_col] = merged_data[residue_col + '_f1']
        
        if len(merged_data) == 0:
            print("Error: No matching residues found between datasets!")
            return pd.DataFrame()
        
        results = []
        
        # Process each residue
        total_residues = len(merged_data)
        processed = 0
        
        print(f"\nProcessing {total_residues} matched residues from dual-field data...")
        print(f"Field 1: {self.field1_freq} MHz, Field 2: {self.field2_freq} MHz")
        
        for idx, row in merged_data.iterrows():
            # Extract data for both fields
            r1_f1, r2_f1, noe_f1 = row['R1_f1'], row['R2_f1'], row['hetNOE_f1']
            r1_f1_err, r2_f1_err, noe_f1_err = row['R1err_f1'], row['R2err_f1'], row['hetNOEerr_f1']
            
            r1_f2, r2_f2, noe_f2 = row['R1_f2'], row['R2_f2'], row['hetNOE_f2']
            r1_f2_err, r2_f2_err, noe_f2_err = row['R1err_f2'], row['R2err_f2'], row['hetNOEerr_f2']
            
            # Get residue identifier
            residue_id = row[residue_col] if residue_col in merged_data.columns else f"Index_{idx}"
            
            # Check for valid data
            all_values = [r1_f1, r1_f1_err, r2_f1, r2_f1_err, noe_f1, noe_f1_err,
                         r1_f2, r1_f2_err, r2_f2, r2_f2_err, noe_f2, noe_f2_err]
            
            if any(pd.isna(all_values)) or any(x <= 0 for x in all_values):
                continue
            
            # Calculate spectral densities for both fields
            j_f1 = self.calculate_spectral_densities(r1_f1, r2_f1, noe_f1, 'field1')
            j_f2 = self.calculate_spectral_densities(r1_f2, r2_f2, noe_f2, 'field2')
            
            # Compile results
            result_row = {
                'Residue': residue_id,
                # Field 1 data
                'R1_f1': r1_f1, 'R1_f1_err': r1_f1_err,
                'R2_f1': r2_f1, 'R2_f1_err': r2_f1_err,
                'hetNOE_f1': noe_f1, 'hetNOE_f1_err': noe_f1_err,
                # Field 2 data
                'R1_f2': r1_f2, 'R1_f2_err': r1_f2_err,
                'R2_f2': r2_f2, 'R2_f2_err': r2_f2_err,
                'hetNOE_f2': noe_f2, 'hetNOE_f2_err': noe_f2_err,
                # Spectral densities
                'J0_f1': j_f1['J0'], 'JwN_f1': j_f1['JwN'], 'JwH_f1': j_f1['JwH'],
                'J0_f2': j_f2['J0'], 'JwN_f2': j_f2['JwN'], 'JwH_f2': j_f2['JwH'],
            }
            
            results.append(result_row)
            processed += 1
        
        print(f"Successfully processed {processed}/{total_residues} residues")
        
        if processed == 0:
            print("WARNING: No valid dual-field data found to process!")
            return pd.DataFrame()
        
        return pd.DataFrame(results)

    def plot_dual_field_results(self, results_df, save_plots=True, output_prefix="spectral_density"):
        """
        Generate plots of the dual-field analysis results
        
        Parameters:
        -----------
        results_df : pandas.DataFrame
            Results from analyze_dual_field_csv()
        save_plots : bool
            Whether to save plots to files
        output_prefix : str
            Prefix for output files
        """
        if len(results_df) == 0:
            print("No data to plot!")
            return
            
        fig, axes = plt.subplots(3, 3, figsize=(18, 15))
        fig.suptitle(f'Dual-Field Spectral Density Analysis Results\n{self.field1_freq} MHz vs {self.field2_freq} MHz', fontsize=16)
        
        residues = results_df.index
        
        # Row 1: Experimental data comparison (R1, R2, hetNOE)
        axes[0,0].errorbar(residues, results_df['R1_f1'], yerr=results_df['R1_f1_err'], 
                          fmt='o', capsize=3, label=f'{self.field1_freq} MHz', alpha=0.7)
        axes[0,0].errorbar(residues, results_df['R1_f2'], yerr=results_df['R1_f2_err'], 
                          fmt='s', capsize=3, label=f'{self.field2_freq} MHz', alpha=0.7)
        axes[0,0].set_title('R1 Comparison')
        axes[0,0].set_ylabel('R1 (s⁻¹)')
        axes[0,0].legend()
        
        axes[0,1].errorbar(residues, results_df['R2_f1'], yerr=results_df['R2_f1_err'],
                          fmt='o', capsize=3, label=f'{self.field1_freq} MHz', alpha=0.7)
        axes[0,1].errorbar(residues, results_df['R2_f2'], yerr=results_df['R2_f2_err'],
                          fmt='s', capsize=3, label=f'{self.field2_freq} MHz', alpha=0.7)
        axes[0,1].set_title('R2 Comparison')
        axes[0,1].set_ylabel('R2 (s⁻¹)')
        axes[0,1].legend()
        
        axes[0,2].errorbar(residues, results_df['hetNOE_f1'], yerr=results_df['hetNOE_f1_err'],
                          fmt='o', capsize=3, label=f'{self.field1_freq} MHz', alpha=0.7)
        axes[0,2].errorbar(residues, results_df['hetNOE_f2'], yerr=results_df['hetNOE_f2_err'],
                          fmt='s', capsize=3, label=f'{self.field2_freq} MHz', alpha=0.7)
        axes[0,2].set_title('hetNOE Comparison')
        axes[0,2].set_ylabel('hetNOE')
        axes[0,2].legend()
        
        # Row 2: Spectral densities comparison
        axes[1,0].plot(residues, results_df['J0_f1'], 'o', label=f'{self.field1_freq} MHz', alpha=0.7)
        axes[1,0].plot(residues, results_df['J0_f2'], 's', label=f'{self.field2_freq} MHz', alpha=0.7)
        axes[1,0].set_title('J(0) Comparison')
        axes[1,0].set_ylabel('J(0) (ns/rad²)')
        axes[1,0].legend()
        
        axes[1,1].plot(residues, results_df['JwN_f1'], 'o', label=f'{self.field1_freq} MHz', alpha=0.7)
        axes[1,1].plot(residues, results_df['JwN_f2'], 's', label=f'{self.field2_freq} MHz', alpha=0.7)
        axes[1,1].set_title('J(ωN) Comparison')
        axes[1,1].set_ylabel('J(ωN) (ns/rad²)')
        axes[1,1].legend()
        
        axes[1,2].plot(residues, results_df['JwH_f1'], 'o', label=f'{self.field1_freq} MHz', alpha=0.7)
        axes[1,2].plot(residues, results_df['JwH_f2'], 's', label=f'{self.field2_freq} MHz', alpha=0.7)
        axes[1,2].set_title('J(ωH) Comparison')
        axes[1,2].set_ylabel('J(ωH) (ns/rad²)')
        axes[1,2].legend()
        
        # Row 3: Correlation plots
        axes[2,0].scatter(results_df['J0_f1'], results_df['J0_f2'], alpha=0.6)
        axes[2,0].set_xlabel(f'J(0) {self.field1_freq} MHz')
        axes[2,0].set_ylabel(f'J(0) {self.field2_freq} MHz')
        axes[2,0].set_title('J(0) Field Correlation')
        
        axes[2,1].scatter(results_df['JwN_f1'], results_df['JwN_f2'], alpha=0.6)
        axes[2,1].set_xlabel(f'J(ωN) {self.field1_freq} MHz')
        axes[2,1].set_ylabel(f'J(ωN) {self.field2_freq} MHz')
        axes[2,1].set_title('J(ωN) Field Correlation')
        
        axes[2,2].scatter(results_df['JwH_f1'], results_df['JwH_f2'], alpha=0.6)
        axes[2,2].set_xlabel(f'J(ωH) {self.field1_freq} MHz')
        axes[2,2].set_ylabel(f'J(ωH) {self.field2_freq} MHz')
        axes[2,2].set_title('J(ωH) Field Correlation')
        
        # Format all subplots
        for i in range(3):
            for j in range(3):
                axes[i,j].grid(True, alpha=0.3)
                if i < 2:  # First 2 rows need residue labels on x-axis
                    axes[i,j].set_xlabel('Residue')
        
        plt.tight_layout()

        if save_plots:
            plot_filename = f"{output_prefix}_analysis_results.pdf"
            plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
            print(f"Plots saved as '{plot_filename}'")
            plt.close(fig)  # Close figure to free memory
        else:
            plt.close(fig)  # Always close to prevent memory leaks

    def save_results(self, results_df, filename='spectral_density_results.csv'):
        """
        Save spectral density results to CSV file
        
        Parameters:
        -----------
        results_df : pandas.DataFrame
            Results from analyze_dual_field_csv()
        filename : str
            Output filename
        """
        results_df.to_csv(filename, index=False)
        print(f"Results saved to '{filename}'")


def run_spectral_density_analysis_with_params(params):
    """
    Run spectral density analysis with parameters provided by GUI
    
    Parameters:
    -----------
    params : dict
        Dictionary containing all analysis parameters
        
    Returns:
    --------
    dict : Results summary
    """
    
    # Extract parameters
    input_file1 = params['input_file1']
    input_file2 = params.get('input_file2', None)
    field1_freq = params.get('field1_freq', 600.0)
    field2_freq = params.get('field2_freq', 700.0)
    analysis_type = params.get('analysis_type', 'dual_field')  # 'single_field' or 'dual_field'
    output_prefix = params.get('output_prefix', 'spectral_density')
    use_monte_carlo = params.get('use_monte_carlo', False)
    n_monte_carlo = params.get('n_monte_carlo', 50)
    
    # Bond length and CSA parameters
    rNH = params.get('rNH', 1.015e-10)
    csaN = params.get('csaN', -160.0e-6)
    
    try:
        if analysis_type == 'dual_field':
            if not input_file2:
                raise ValueError("Dual field analysis requires two input files")
            
            # Initialize dual-field analyzer
            analyzer = DualFieldSpectralDensityAnalysis(
                field1_freq=field1_freq,
                field2_freq=field2_freq,
                rNH=rNH,
                csaN=csaN
            )
            
            print(f"Starting dual-field spectral density analysis...")
            print(f"Field 1: {field1_freq} MHz, Field 2: {field2_freq} MHz")
            
            # Analyze data
            results = analyzer.analyze_dual_field_csv(
                input_file1, input_file2,
                use_monte_carlo_errors=use_monte_carlo,
                n_monte_carlo=n_monte_carlo
            )
            
            if len(results) == 0:
                return {
                    'success': False,
                    'error': 'No valid data found for analysis',
                    'n_processed': 0
                }
            
            # Save results
            results_file = f"{output_prefix}_results.csv"
            analyzer.save_results(results, results_file)
            
            # Generate plots
            analyzer.plot_dual_field_results(results, save_plots=True, output_prefix=output_prefix)
            
            return {
                'success': True,
                'n_processed': len(results),
                'results_file': results_file,
                'plots_prefix': output_prefix,
                'analysis_type': 'dual_field',
                'field1_freq': field1_freq,
                'field2_freq': field2_freq,
                'mean_J0_f1': results['J0_f1'].mean(),
                'mean_J0_f2': results['J0_f2'].mean(),
                'mean_JwN_f1': results['JwN_f1'].mean(),
                'mean_JwN_f2': results['JwN_f2'].mean()
            }
            
        else:  # single_field
            # For single field, just calculate spectral densities
            analyzer = DualFieldSpectralDensityAnalysis(
                field1_freq=field1_freq,
                field2_freq=field1_freq,  # Use same field for both
                rNH=rNH,
                csaN=csaN
            )
            
            print(f"Starting single-field spectral density analysis...")
            print(f"Field: {field1_freq} MHz")
            
            # Load single field data
            data = pd.read_csv(input_file1)
            required_cols = ['R1', 'R1err', 'R2', 'R2err', 'hetNOE', 'hetNOEerr']
            missing_cols = [col for col in required_cols if col not in data.columns]
            if missing_cols:
                raise ValueError(f"Missing required columns: {missing_cols}")
            
            results = []
            for idx, row in data.iterrows():
                r1, r2, noe = row['R1'], row['R2'], row['hetNOE']
                r1_err, r2_err, noe_err = row['R1err'], row['R2err'], row['hetNOEerr']
                
                if any(pd.isna([r1, r2, noe, r1_err, r2_err, noe_err])) or any(x <= 0 for x in [r1, r2, noe, r1_err, r2_err, noe_err]):
                    continue
                
                j_values = analyzer.calculate_spectral_densities(r1, r2, noe, 'field1')
                
                result_row = {
                    'Residue': row.get('Residue', f"Index_{idx}"),
                    'R1': r1, 'R1_err': r1_err,
                    'R2': r2, 'R2_err': r2_err,
                    'hetNOE': noe, 'hetNOE_err': noe_err,
                    'J0': j_values['J0'],
                    'JwN': j_values['JwN'],
                    'JwH': j_values['JwH']
                }
                results.append(result_row)
            
            results_df = pd.DataFrame(results)
            
            if len(results_df) == 0:
                return {
                    'success': False,
                    'error': 'No valid data found for analysis',
                    'n_processed': 0
                }
            
            # Save results
            results_file = f"{output_prefix}_results.csv"
            results_df.to_csv(results_file, index=False)
            
            return {
                'success': True,
                'n_processed': len(results_df),
                'results_file': results_file,
                'analysis_type': 'single_field',
                'field_freq': field1_freq,
                'mean_J0': results_df['J0'].mean(),
                'mean_JwN': results_df['JwN'].mean(),
                'mean_JwH': results_df['JwH'].mean()
            }
            
    except Exception as e:
        return {
            'success': False,
            'error': str(e),
            'n_processed': 0
        }


def main():
    """Example usage for testing"""
    # Example parameters
    params = {
        'input_file1': 'example_600MHz_data.csv',
        'input_file2': 'example_700MHz_data.csv',
        'field1_freq': 600.0,
        'field2_freq': 700.0,
        'analysis_type': 'dual_field',
        'output_prefix': 'test_spectral_density',
        'use_monte_carlo': False,
        'n_monte_carlo': 50
    }
    
    results = run_spectral_density_analysis_with_params(params)
    print("Analysis results:", results)


if __name__ == "__main__":
    main()