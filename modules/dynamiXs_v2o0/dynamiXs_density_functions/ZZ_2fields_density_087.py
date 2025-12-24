#!/usr/bin/env python3
"""
Dual-Field Reduced Spectral Density Mapping Analysis Script with J(0.87ωH)

This script imports R1, R1err, R2, R2err, hetNOE, and hetNOEerr from two CSV files
at different field strengths and performs reduced spectral density calculations 
including J(0), J(wN), J(0.87wH), S2, Rex, and te parameters with proper error propagation.

The dual-field approach provides better separation of chemical exchange (Rex) from 
internal dynamics (te) and more reliable model-free parameter determination.

This version uses J(0.87ωH) instead of J(ωH) for more accurate treatment of
cross-correlation effects and dipolar-CSA interference in the relaxation mechanisms.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend (must be before pyplot import)
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

# Cross-correlation factor for J(0.87ωH)
OMEGA_H_FACTOR = 0.87  # Factor for more accurate spectral density calculation

class DualFieldSpectralDensityAnalysis:
    """
    Class for performing dual-field reduced spectral density mapping analysis
    using J(0.87ωH) for improved accuracy
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
        
        # Calculate effective proton frequencies with 0.87 factor
        self.field1_omegaH_eff = self.field1_omegaH * OMEGA_H_FACTOR
        self.field2_omegaH_eff = self.field2_omegaH * OMEGA_H_FACTOR
        
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
        Calculate spectral densities for a given field using J(0.87ωH)
        
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

        # Farrow Eq. 5: J(0.87ωH) calculation
        jwh_087 = 4.0 * sigma_noe / (5.0 * d_factor)

        return {'J0': j0, 'JwN': jwn, 'JwH_087': jwh_087}
    
    def calculate_isotropic_spectral_density(self, s2, tc, te, omega):
        """
        Calculate J(w) using Lipari-Szabo model-free approach
        """
        if te == 0:
            j = s2 * tc / (1.0 + (omega * tc)**2)
        else:
            teff = 1.0 / (1.0/tc + 1.0/te)
            j = (s2 * tc / (1.0 + (omega * tc)**2) + 
                 (1.0 - s2) * teff / (1.0 + (omega * teff)**2))
        
        return j * 0.4  # 2/5 prefactor
    
    def estimate_overall_correlation_time(self, r1_data1, r2_data1, r1_data2, r2_data2):
        """
        Estimate overall correlation time from dual-field R1/R2 data
        Uses the higher field data for better sensitivity
        """
        # Use higher field data for tc estimation (better sensitivity)
        if self.field2_freq > self.field1_freq:
            r1_array, r2_array = r1_data2, r2_data2
            omegaN = abs(self.field2_omegaN)
        else:
            r1_array, r2_array = r1_data1, r2_data1  
            omegaN = abs(self.field1_omegaN)
            
        t1_array = 1.0 / r1_array
        t2_array = 1.0 / r2_array
        
        # Fushman method: tc = [1/(2*omegaN)] * sqrt[(6*T1/T2) - 7]
        ratio = 6.0 * t1_array / t2_array
        valid_mask = ratio > 7.0
        
        if np.sum(valid_mask) < 3:
            tc_est = 1.0 / (2.0 * omegaN) * np.sqrt(np.mean(ratio) - 7.0) if np.mean(ratio) > 7 else 10e-9
        else:
            tc_values = 1.0 / (2.0 * omegaN) * np.sqrt(ratio[valid_mask] - 7.0)
            tc_est = np.median(tc_values)
            
        return max(1e-12, min(100e-9, tc_est))

    def _fit_dual_field_dataset(self, r1_f1, r2_f1, noe_f1, r1_f2, r2_f2, noe_f2, tc_fixed, niter):
        """
        Fit dual-field dataset for a single residue using J(0.87ωH)
        
        Parameters:
        -----------
        r1_f1, r2_f1, noe_f1 : float
            Field 1 relaxation parameters
        r1_f2, r2_f2, noe_f2 : float  
            Field 2 relaxation parameters
        tc_fixed : float
            Fixed overall correlation time
        niter : int
            Number of iterations
            
        Returns:
        --------
        dict : Fitted parameters
        """
        from random import random, randint
        from math import exp
        
        # Convert to T1, T2 for both fields
        t1_f1, t2_f1 = 1.0 / r1_f1, 1.0 / r2_f1
        t1_f2, t2_f2 = 1.0 / r1_f2, 1.0 / r2_f2
        
        # Constants for both fields
        A1, A2 = self.field1_d_factor, self.field2_d_factor  # Same for both fields
        C1, C2 = self.field1_c_factor, self.field2_c_factor  # Different for each field
        gammaHN = GAMMA_H / GAMMA_N
        
        # Initialize parameter ensemble
        ensemble = []
        
        # Parameter starting values
        s2_start = [x * 0.05 for x in range(1, 21)]  # 0.05 to 1.0
        te_start = [x * 1e-12 for x in [0.1, 0.5, 1, 2, 5, 10, 20, 50]]  # 1ps to 1ns
        rex1_start = [float(x) for x in range(16)]  # 0 to 15 s-1 for field 1
        rex2_start = [float(x) for x in range(16)]  # 0 to 15 s-1 for field 2
        
        # Create initial ensemble - now includes Rex for both fields
        for s2 in s2_start:
            for te in te_start:
                for rex1 in rex1_start:
                    for rex2 in rex2_start:
                        ensemble.append((None, s2, te, rex1, rex2))
        
        # Calculate initial scores
        for k, (score, s2, te, rex1, rex2) in enumerate(ensemble):
            # Calculate spectral densities for both fields using 0.87ωH
            # Field 1
            jH1 = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, self.field1_omegaH_eff)
            jN1 = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, abs(self.field1_omegaN))
            jHpN1 = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, self.field1_omegaH + abs(self.field1_omegaN))
            jHmN1 = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, self.field1_omegaH - abs(self.field1_omegaN))
            j01 = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, 0.0)
            
            # Field 2
            jH2 = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, self.field2_omegaH_eff)
            jN2 = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, abs(self.field2_omegaN))
            jHpN2 = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, self.field2_omegaH + abs(self.field2_omegaN))
            jHmN2 = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, self.field2_omegaH - abs(self.field2_omegaN))
            j02 = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, 0.0)
            
            # Calculate predicted relaxation parameters for both fields
            # Field 1
            r1_calc_f1 = (A1 / 4.0) * (3*jN1 + 6*jHpN1 + jHmN1) + C1 * jN1
            r2_calc_f1 = 0.5 * (A1 / 4.0) * (4*j01 + 3*jN1 + 6*jHpN1 + 6*jH1 + jHmN1) + C1 * (2*j01/3.0 + 0.5*jN1) + rex1

            # Field 2
            r1_calc_f2 = (A2 / 4.0) * (3*jN2 + 6*jHpN2 + jHmN2) + C2 * jN2
            r2_calc_f2 = 0.5 * (A2 / 4.0) * (4*j02 + 3*jN2 + 6*jHpN2 + 6*jH2 + jHmN2) + C2 * (2*j02/3.0 + 0.5*jN2) + rex2
            
            t1_calc_f1, t2_calc_f1 = 1.0 / r1_calc_f1, 1.0 / r2_calc_f1
            t1_calc_f2, t2_calc_f2 = 1.0 / r1_calc_f2, 1.0 / r2_calc_f2
            
            # Calculate combined score (chi-squared for both fields)
            d1_f1 = (t1_calc_f1 - t1_f1) / t1_f1
            d2_f1 = (t2_calc_f1 - t2_f1) / t2_f1
            d1_f2 = (t1_calc_f2 - t1_f2) / t1_f2
            d2_f2 = (t2_calc_f2 - t2_f2) / t2_f2
            
            score = d1_f1*d1_f1 + d2_f1*d2_f1 + d1_f2*d1_f2 + d2_f2*d2_f2
            
            # Add NOE contributions for both fields
            if noe_f1 is not None and not pd.isna(noe_f1):
                noe_calc_f1 = 1.0 + ((A1 / 4.0) * gammaHN * t1_calc_f1 * (6*jHpN1 - jHmN1))
                dn_f1 = (noe_calc_f1 - noe_f1) / (0.1 if noe_f1 == 0 else abs(noe_f1))
                score += dn_f1*dn_f1

            if noe_f2 is not None and not pd.isna(noe_f2):
                noe_calc_f2 = 1.0 + ((A2 / 4.0) * gammaHN * t1_calc_f2 * (6*jHpN2 - jHmN2))
                dn_f2 = (noe_calc_f2 - noe_f2) / (0.1 if noe_f2 == 0 else abs(noe_f2))
                score += dn_f2*dn_f2
                
            ensemble[k] = (score, s2, te, rex1, rex2, t1_calc_f1, t2_calc_f1, t1_calc_f2, t2_calc_f2)
        
        # Sort and keep best 10
        ensemble.sort()
        ensemble = ensemble[:10]
        ensemble_size = len(ensemble)
        
        # Monte Carlo optimization
        for i in range(niter):
            f = exp(-10.0 * i / float(niter))  # Cooling factor
            
            ensemble.sort()
            if ensemble[0][0] < 1e-10:  # Converged
                break
                
            prevScore, s2, te, rex1, rex2, t1_calc_f1, t2_calc_f1, t1_calc_f2, t2_calc_f2 = ensemble[-1]
            
            # Mutate parameters
            d = ((random() - 0.382) * f) + 1.0
            s2_new = max(0.01, min(1.0, s2 * d))
            
            d = ((random() - 0.382) * f) + 1.0  
            te_new = max(1e-15, min(10e-10, te * d))
            
            d = ((random() - 0.382) * f) + 1.0
            rex1_new = max(0.0, min(50.0, rex1 * d))
            
            d = ((random() - 0.382) * f) + 1.0
            rex2_new = max(0.0, min(50.0, rex2 * d))
            
            # Calculate new score with dual-field data using 0.87ωH
            # Field 1
            jH1 = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, self.field1_omegaH_eff)
            jN1 = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, abs(self.field1_omegaN))
            jHpN1 = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, self.field1_omegaH + abs(self.field1_omegaN))
            jHmN1 = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, self.field1_omegaH - abs(self.field1_omegaN))
            j01 = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, 0.0)
            
            # Field 2
            jH2 = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, self.field2_omegaH_eff)
            jN2 = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, abs(self.field2_omegaN))
            jHpN2 = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, self.field2_omegaH + abs(self.field2_omegaN))
            jHmN2 = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, self.field2_omegaH - abs(self.field2_omegaN))
            j02 = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, 0.0)
            
            # Calculate relaxation parameters
            r1_calc_f1 = (A1 / 4.0) * (3*jN1 + 6*jHpN1 + jHmN1) + C1 * jN1
            r2_calc_f1 = 0.5 * (A1 / 4.0) * (4*j01 + 3*jN1 + 6*jHpN1 + 6*jH1 + jHmN1) + C1 * (2*j01/3.0 + 0.5*jN1) + rex1_new
            r1_calc_f2 = (A2 / 4.0) * (3*jN2 + 6*jHpN2 + jHmN2) + C2 * jN2
            r2_calc_f2 = 0.5 * (A2 / 4.0) * (4*j02 + 3*jN2 + 6*jHpN2 + 6*jH2 + jHmN2) + C2 * (2*j02/3.0 + 0.5*jN2) + rex2_new
            
            t1_calc_f1, t2_calc_f1 = 1.0 / r1_calc_f1, 1.0 / r2_calc_f1
            t1_calc_f2, t2_calc_f2 = 1.0 / r1_calc_f2, 1.0 / r2_calc_f2
            
            d1_f1 = (t1_calc_f1 - t1_f1) / t1_f1
            d2_f1 = (t2_calc_f1 - t2_f1) / t2_f1
            d1_f2 = (t1_calc_f2 - t1_f2) / t1_f2
            d2_f2 = (t2_calc_f2 - t2_f2) / t2_f2
            
            score = d1_f1*d1_f1 + d2_f1*d2_f1 + d1_f2*d1_f2 + d2_f2*d2_f2
            
            # Add NOE contributions
            if noe_f1 is not None and not pd.isna(noe_f1):
                noe_calc_f1 = 1.0 + ((A1 / 4.0) * gammaHN * t1_calc_f1 * (6*jHpN1 - jHmN1))
                dn_f1 = (noe_calc_f1 - noe_f1) / (0.1 if noe_f1 == 0 else abs(noe_f1))
                score += dn_f1*dn_f1

            if noe_f2 is not None and not pd.isna(noe_f2):
                noe_calc_f2 = 1.0 + ((A2 / 4.0) * gammaHN * t1_calc_f2 * (6*jHpN2 - jHmN2))
                dn_f2 = (noe_calc_f2 - noe_f2) / (0.1 if noe_f2 == 0 else abs(noe_f2))
                score += dn_f2*dn_f2
            
            # Accept or reject
            ratio = exp(prevScore - score)
            if ratio > 1.0:
                ensemble[-1] = (score, s2_new, te_new, rex1_new, rex2_new, t1_calc_f1, t2_calc_f1, t1_calc_f2, t2_calc_f2)
            else:
                k = randint(0, ensemble_size - 1)
                ensemble[-1] = ensemble[k]
        
        # Return best result
        ensemble.sort()
        if len(ensemble) > 0:
            score, s2_best, te_best, rex1_best, rex2_best, t1_calc_f1, t2_calc_f1, t1_calc_f2, t2_calc_f2 = ensemble[0]
            
            # Calculate errors from ensemble spread
            s2_values = [e[1] for e in ensemble[:10]]
            te_values = [e[2] for e in ensemble[:10]]
            rex1_values = [e[3] for e in ensemble[:10]]
            rex2_values = [e[4] for e in ensemble[:10]]
            
            return {
                'S2': s2_best,
                'tc': tc_fixed * 1e9,  # Convert to ns
                'te': te_best * 1e12,  # Convert to ps
                'Rex_field1': rex1_best,
                'Rex_field2': rex2_best,
                'S2_err': np.std(s2_values),
                'te_err': np.std(te_values) * 1e12,
                'Rex_field1_err': np.std(rex1_values),
                'Rex_field2_err': np.std(rex2_values),
                'fit_success': True,
                'chi2': score,
                'iterations': i
            }
        else:
            return {
                'S2': np.nan, 'tc': np.nan, 'te': np.nan, 
                'Rex_field1': np.nan, 'Rex_field2': np.nan,
                'S2_err': np.nan, 'te_err': np.nan, 
                'Rex_field1_err': np.nan, 'Rex_field2_err': np.nan,
                'fit_success': False, 'chi2': np.inf, 'iterations': niter
            }

    def fit_dual_field_model_free(self, r1_f1, r2_f1, noe_f1, r1_f2, r2_f2, noe_f2, tc_fixed,
                                r1_f1_err=None, r2_f1_err=None, noe_f1_err=None,
                                r1_f2_err=None, r2_f2_err=None, noe_f2_err=None,
                                n_monte_carlo=100, niter=10_000):
        """
        Fit dual-field model-free parameters with error propagation using J(0.87ωH)
        
        Parameters:
        -----------
        r1_f1, r2_f1, noe_f1 : float
            Field 1 relaxation parameters
        r1_f2, r2_f2, noe_f2 : float
            Field 2 relaxation parameters  
        tc_fixed : float
            Fixed overall correlation time
        *_err : float, optional
            Experimental uncertainties for Monte Carlo error propagation
        n_monte_carlo : int
            Number of Monte Carlo samples
        niter : int
            Iterations per fit
            
        Returns:
        --------
        dict : Fitted parameters with errors
        """
        # First, fit the actual experimental data
        best_fit = self._fit_dual_field_dataset(r1_f1, r2_f1, noe_f1, r1_f2, r2_f2, noe_f2, tc_fixed, niter)
        
        # Check if errors are provided for Monte Carlo
        errors_provided = all(err is not None for err in [r1_f1_err, r2_f1_err, noe_f1_err,
                                                         r1_f2_err, r2_f2_err, noe_f2_err])

        if not errors_provided or not best_fit['fit_success']:
            return best_fit

        # Monte Carlo error propagation
        print(f"  Running dual-field Monte Carlo error analysis ({n_monte_carlo} samples)...", end='', flush=True)
        mc_results = {'S2': [], 'te': [], 'Rex_field1': [], 'Rex_field2': [], 'chi2': []}
        
        mc_niter = max(1000, niter // 10)
        
        for i in range(n_monte_carlo):
            # Sample within experimental uncertainties for both fields
            r1_f1_sample = max(0.1, np.random.normal(r1_f1, r1_f1_err))
            r2_f1_sample = max(r1_f1_sample * 1.1, np.random.normal(r2_f1, r2_f1_err))
            noe_f1_sample = np.clip(np.random.normal(noe_f1, noe_f1_err), -0.5, 1.0)
            
            r1_f2_sample = max(0.1, np.random.normal(r1_f2, r1_f2_err))
            r2_f2_sample = max(r1_f2_sample * 1.1, np.random.normal(r2_f2, r2_f2_err))
            noe_f2_sample = np.clip(np.random.normal(noe_f2, noe_f2_err), -0.5, 1.0)
            
            # Fit the synthetic dual-field dataset
            mc_fit = self._fit_dual_field_dataset(r1_f1_sample, r2_f1_sample, noe_f1_sample,
                                                r1_f2_sample, r2_f2_sample, noe_f2_sample,
                                                tc_fixed, mc_niter)
            
            if mc_fit['fit_success']:
                mc_results['S2'].append(mc_fit['S2'])
                mc_results['te'].append(mc_fit['te'])
                mc_results['Rex_field1'].append(mc_fit['Rex_field1'])
                mc_results['Rex_field2'].append(mc_fit['Rex_field2'])
                mc_results['chi2'].append(mc_fit['chi2'])
        
        print(" done")
        
        # Calculate errors from Monte Carlo distribution
        n_successful = len(mc_results['S2'])
        if n_successful < 10:
            print(f"    Warning: Only {n_successful}/{n_monte_carlo} dual-field MC fits succeeded")
            return best_fit
        
        # Update errors with Monte Carlo results
        best_fit['S2_err'] = np.std(mc_results['S2'])
        best_fit['te_err'] = np.std(mc_results['te'])
        best_fit['Rex_field1_err'] = np.std(mc_results['Rex_field1'])
        best_fit['Rex_field2_err'] = np.std(mc_results['Rex_field2'])
        
        # Add confidence intervals
        best_fit['S2_95CI'] = np.percentile(mc_results['S2'], [2.5, 97.5])
        best_fit['te_95CI'] = np.percentile(mc_results['te'], [2.5, 97.5])
        best_fit['Rex_field1_95CI'] = np.percentile(mc_results['Rex_field1'], [2.5, 97.5])
        best_fit['Rex_field2_95CI'] = np.percentile(mc_results['Rex_field2'], [2.5, 97.5])
        
        best_fit['mc_success_rate'] = n_successful / n_monte_carlo
        
        return best_fit
    
    def propagate_dual_field_errors(self, r1_f1, r1_f1_err, r2_f1, r2_f1_err, noe_f1, noe_f1_err,
                                  r1_f2, r1_f2_err, r2_f2, r2_f2_err, noe_f2, noe_f2_err):
        """
        Propagate errors through dual-field spectral density calculations using J(0.87ωH)
        
        Parameters:
        -----------
        r1_f1, r1_f1_err, r2_f1, r2_f1_err, noe_f1, noe_f1_err : float
            Field 1 relaxation parameters and errors
        r1_f2, r1_f2_err, r2_f2, r2_f2_err, noe_f2, noe_f2_err : float
            Field 2 relaxation parameters and errors
            
        Returns:
        --------
        dict : Spectral density values and errors for both fields
        """
        # Calculate spectral densities for both fields
        j_f1 = self.calculate_spectral_densities(r1_f1, r2_f1, noe_f1, 'field1')
        j_f2 = self.calculate_spectral_densities(r1_f2, r2_f2, noe_f2, 'field2')
        
        # For error propagation, use analytical derivatives (simplified version)
        # In practice, Monte Carlo is more robust for complex error propagation
        
        # Field 1 errors (using same approach as single field)
        d1 = self.field1_d_factor
        c1 = self.field1_c_factor
        gamma_ratio = GAMMA_N / GAMMA_H
        
        sigma_noe_f1 = (noe_f1 - 1.0) * r1_f1 * gamma_ratio
        dsigma_dr1_f1 = (noe_f1 - 1.0) * gamma_ratio
        dsigma_dnoe_f1 = r1_f1 * gamma_ratio
        sigma_noe_f1_err = np.sqrt((dsigma_dr1_f1 * r1_f1_err)**2 + (dsigma_dnoe_f1 * noe_f1_err)**2)
        
        # J(0) field 1
        factor_j0_f1 = 3.0 / (2.0 * (3.0 * d1 / 4.0 + c1))
        dj0_dr1_f1 = factor_j0_f1 * (-0.5 - (3.0/5.0) * dsigma_dr1_f1)
        dj0_dr2_f1 = factor_j0_f1
        dj0_dnoe_f1 = factor_j0_f1 * (-(3.0/5.0) * dsigma_dnoe_f1)
        j0_f1_err = np.sqrt((dj0_dr1_f1 * r1_f1_err)**2 + (dj0_dr2_f1 * r2_f1_err)**2 + (dj0_dnoe_f1 * noe_f1_err)**2)
        
        # Field 2 errors (similar calculation)
        d2 = self.field2_d_factor
        c2 = self.field2_c_factor
        
        sigma_noe_f2 = (noe_f2 - 1.0) * r1_f2 * gamma_ratio
        dsigma_dr1_f2 = (noe_f2 - 1.0) * gamma_ratio
        dsigma_dnoe_f2 = r1_f2 * gamma_ratio
        sigma_noe_f2_err = np.sqrt((dsigma_dr1_f2 * r1_f2_err)**2 + (dsigma_dnoe_f2 * noe_f2_err)**2)
        
        factor_j0_f2 = 3.0 / (2.0 * (3.0 * d2 / 4.0 + c2))
        dj0_dr1_f2 = factor_j0_f2 * (-0.5 - (3.0/5.0) * dsigma_dr1_f2)
        dj0_dr2_f2 = factor_j0_f2
        dj0_dnoe_f2 = factor_j0_f2 * (-(3.0/5.0) * dsigma_dnoe_f2)
        j0_f2_err = np.sqrt((dj0_dr1_f2 * r1_f2_err)**2 + (dj0_dr2_f2 * r2_f2_err)**2 + (dj0_dnoe_f2 * noe_f2_err)**2)
    
########################    ########################    ########################  
########################    ########################    ########################  
########################    ########################    ########################  
  
    
        # J(ωN) error for field 1
        factor_jwn_f1 = 1.0 / (3.0 * d1 / 4.0 + c1)
        djwn_dr1_f1 = factor_jwn_f1 * (1.0 - (7.0/5.0) * (noe_f1 - 1.0) * gamma_ratio)
        djwn_dnoe_f1 = factor_jwn_f1 * (-(7.0/5.0) * r1_f1 * gamma_ratio)
        jwn_f1_err = np.sqrt((djwn_dr1_f1 * r1_f1_err)**2 + (djwn_dnoe_f1 * noe_f1_err)**2)

        # J(0.87ωH) error for field 1  
        factor_jwh_f1 = 4.0 / (5.0 * d1)
        djwh_dr1_f1 = factor_jwh_f1 * (noe_f1 - 1.0) * gamma_ratio
        djwh_dnoe_f1 = factor_jwh_f1 * r1_f1 * gamma_ratio
        jwh_087_f1_err = np.sqrt((djwh_dr1_f1 * r1_f1_err)**2 + (djwh_dnoe_f1 * noe_f1_err)**2) 
        
        # J(ωN) error for field 2
        factor_jwn_f2 = 1.0 / (3.0 * d2 / 4.0 + c2)
        djwn_dr1_f2 = factor_jwn_f2 * (1.0 - (7.0/5.0) * (noe_f2 - 1.0) * gamma_ratio)
        djwn_dnoe_f2 = factor_jwn_f2 * (-(7.0/5.0) * r1_f2 * gamma_ratio)
        jwn_f2_err = np.sqrt((djwn_dr1_f2 * r1_f2_err)**2 + (djwn_dnoe_f2 * noe_f2_err)**2)

        # J(0.87ωH) error for field 2  
        factor_jwh_f2 = 4.0 / (5.0 * d2)
        djwh_dr1_f2 = factor_jwh_f2 * (noe_f2 - 1.0) * gamma_ratio
        djwh_dnoe_f2 = factor_jwh_f2 * r1_f2 * gamma_ratio
        jwh_087_f2_err = np.sqrt((djwh_dr1_f2 * r1_f2_err)**2 + (djwh_dnoe_f2 * noe_f2_err)**2)
               
        
########################    ########################    ########################          
########################    ########################    ########################         
########################    ########################    ########################    
        
        return {
            'Field1': {
                'J0': j_f1['J0'], 'J0_err': j0_f1_err,
                'JwN': j_f1['JwN'], 'JwN_err': jwn_f1_err,
                'JwH_087': j_f1['JwH_087'], 'JwH_087_err': jwh_087_f1_err                         
            },
            'Field2': {
                'J0': j_f2['J0'], 'J0_err': j0_f2_err,             
                'JwN': j_f2['JwN'], 'JwN_err': jwn_f2_err,
                'JwH_087': j_f2['JwH_087'], 'JwH_087_err': jwh_087_f2_err
            }
        }
    
    def analyze_dual_field_csv(self, csv_file1, csv_file2, residue_col='Residue', 
                              use_monte_carlo_errors=True, n_monte_carlo=100):
        """
        Analyze dual-field relaxation data from two CSV files using J(0.87ωH)
        
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

            # Sort numerically by residue ID (extract numeric part for proper ordering)
            import re
            def extract_numeric(res_id):
                """Extract numeric part from residue ID for sorting (e.g., '3.LYS' -> 3, '10' -> 10)"""
                match = re.search(r'\d+', str(res_id))
                return int(match.group()) if match else 0

            merged_data['_sort_key'] = merged_data[residue_col].apply(extract_numeric)
            merged_data = merged_data.sort_values('_sort_key').drop('_sort_key', axis=1).reset_index(drop=True)
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
        
        # Counters for data quality reporting
        total_residues = len(merged_data)
        skipped_nan = 0
        skipped_zero_negative = 0
        processed = 0
        
        print(f"\nProcessing {total_residues} matched residues from dual-field data...")
        print(f"Field 1: {self.field1_freq} MHz, Field 2: {self.field2_freq} MHz")
        print(f"Using J(0.87ωH) approach with factor = {OMEGA_H_FACTOR}")
        if use_monte_carlo_errors:
            print(f"Using Monte Carlo error propagation with {n_monte_carlo} samples per residue")
        
        # First pass: collect valid data for overall tc estimation
        valid_r1_f1, valid_r2_f1 = [], []
        valid_r1_f2, valid_r2_f2 = [], []
        
        for idx, row in merged_data.iterrows():
            # Extract data for both fields
            r1_f1, r2_f1 = row['R1_f1'], row['R2_f1']
            r1_f2, r2_f2 = row['R1_f2'], row['R2_f2']
            
            if not any(pd.isna([r1_f1, r2_f1, r1_f2, r2_f2])) and all(x > 0 for x in [r1_f1, r2_f1, r1_f2, r2_f2]):
                valid_r1_f1.append(r1_f1)
                valid_r2_f1.append(r2_f1)
                valid_r1_f2.append(r1_f2)
                valid_r2_f2.append(r2_f2)
        
        if len(valid_r1_f1) < 3:
            print("ERROR: Not enough valid dual-field R1/R2 data to estimate overall correlation time!")
            return pd.DataFrame()
        
        # Estimate overall correlation time using dual-field data
        tc_overall = self.estimate_overall_correlation_time(
            np.array(valid_r1_f1), np.array(valid_r2_f1),
            np.array(valid_r1_f2), np.array(valid_r2_f2)
        )
        print(f"Estimated overall correlation time: {tc_overall*1e9:.1f} ns")
        
        # Second pass: process individual residues
        for idx, row in merged_data.iterrows():
            # Extract data for both fields
            r1_f1 = row['R1_f1']
            r1_f1_err = row['R1err_f1']
            r2_f1 = row['R2_f1']
            r2_f1_err = row['R2err_f1']
            noe_f1 = row['hetNOE_f1']
            noe_f1_err = row['hetNOEerr_f1']
            
            r1_f2 = row['R1_f2']
            r1_f2_err = row['R1err_f2']
            r2_f2 = row['R2_f2']
            r2_f2_err = row['R2err_f2']
            noe_f2 = row['hetNOE_f2']
            noe_f2_err = row['hetNOEerr_f2']
            
            # Get residue identifier
            residue_id = row[residue_col] if residue_col in merged_data.columns else f"Index_{idx}"
            
            # Check for NaN values in either field
            all_values = [r1_f1, r1_f1_err, r2_f1, r2_f1_err, noe_f1, noe_f1_err,
                         r1_f2, r1_f2_err, r2_f2, r2_f2_err, noe_f2, noe_f2_err]
            
            if any(pd.isna(all_values)):
                skipped_nan += 1
                print(f"  Skipping {residue_id}: Contains NaN values in one or both fields")
                continue
            
            # Check for invalid values
            invalid_conditions = [
                r1_f1 <= 0, r2_f1 <= 0, r1_f1_err <= 0, r2_f1_err <= 0, noe_f1 == 0, noe_f1_err <= 0,
                r1_f2 <= 0, r2_f2 <= 0, r1_f2_err <= 0, r2_f2_err <= 0, noe_f2 == 0, noe_f2_err <= 0
            ]
            
            if any(invalid_conditions):
                skipped_zero_negative += 1
                print(f"  Skipping {residue_id}: Invalid values in dual-field data")
                continue
            
            # Physical reasonableness checks
            warnings = []
            if r1_f1 > 10 or r1_f2 > 10: warnings.append("unusually high R1")
            if r2_f1 > 100 or r2_f2 > 100: warnings.append("unusually high R2")
            if noe_f1 < 0.2 or noe_f1 > 1.2 or noe_f2 < 0.2 or noe_f2 > 1.2: warnings.append("unusual NOE")
            
            if warnings:
                print(f"  Warning for {residue_id}: {', '.join(warnings)}")
            
            # Calculate spectral densities for both fields with error propagation
            j_results = self.propagate_dual_field_errors(
                r1_f1, r1_f1_err, r2_f1, r2_f1_err, noe_f1, noe_f1_err,
                r1_f2, r1_f2_err, r2_f2, r2_f2_err, noe_f2, noe_f2_err
            )
            
            # Fit dual-field model-free parameters with Monte Carlo error propagation
            mf_results = self.fit_dual_field_model_free(
                r1_f1, r2_f1, noe_f1, r1_f2, r2_f2, noe_f2, tc_overall,
                r1_f1_err=r1_f1_err, r2_f1_err=r2_f1_err, noe_f1_err=noe_f1_err,
                r1_f2_err=r1_f2_err, r2_f2_err=r2_f2_err, noe_f2_err=noe_f2_err,
                n_monte_carlo=n_monte_carlo
            )
            
            # Compile results
            result_row = {
                'Index': idx,
                # Field 1 data
                'R1_f1': r1_f1, 'R1_f1_err': r1_f1_err,
                'R2_f1': r2_f1, 'R2_f1_err': r2_f1_err,
                'hetNOE_f1': noe_f1, 'hetNOE_f1_err': noe_f1_err,
                # Field 2 data
                'R1_f2': r1_f2, 'R1_f2_err': r1_f2_err,
                'R2_f2': r2_f2, 'R2_f2_err': r2_f2_err,
                'hetNOE_f2': noe_f2, 'hetNOE_f2_err': noe_f2_err,
                # Spectral densities for both fields (now with J(0.87ωH))
                'J0_f1': j_results['Field1']['J0'], 'J0_f1_err': j_results['Field1']['J0_err'],
                'JwN_f1': j_results['Field1']['JwN'], 'JwN_f1_err': j_results['Field1']['JwN_err'],
                'JwH_087_f1': j_results['Field1']['JwH_087'], 'JwH_087_f1_err': j_results['Field1']['JwH_087_err'],
                'J0_f2': j_results['Field2']['J0'], 'J0_f2_err': j_results['Field2']['J0_err'],
                'JwN_f2': j_results['Field2']['JwN'], 'JwN_f2_err': j_results['Field2']['JwN_err'],
                'JwH_087_f2': j_results['Field2']['JwH_087'], 'JwH_087_f2_err': j_results['Field2']['JwH_087_err'],
                # Model-free parameters
                **mf_results
            }
            
            # Add residue info if available
            if residue_col in merged_data.columns:
                result_row['Residue'] = row[residue_col]
            
            results.append(result_row)
            processed += 1
        
        # Report data quality summary
        print(f"\nDual-Field Data Quality Summary (J(0.87ωH) Analysis):")
        print(f"  Total matched residues: {total_residues}")
        print(f"  Processed successfully: {processed}")
        print(f"  Skipped (NaN values): {skipped_nan}")
        print(f"  Skipped (zero/negative values): {skipped_zero_negative}")
        print(f"  Success rate: {processed/total_residues*100:.1f}%")
        
        if processed == 0:
            print("  WARNING: No valid dual-field data found to process!")
            return pd.DataFrame()
        
        return pd.DataFrame(results)
    
    def plot_dual_field_results(self, results_df, save_plots=True):
        """
        Generate plots of the dual-field analysis results using J(0.87ωH)
        
        Parameters:
        -----------
        results_df : pandas.DataFrame
                Results from analyze_dual_field_csv()
        save_plots : bool
                Whether to save plots to files
        """
        fig, axes = plt.subplots(4, 3, figsize=(18, 20))
        fig.suptitle(f'Dual-Field Spectral Density Analysis Results (J(0.87ωH))\n{self.field1_freq} MHz vs {self.field2_freq} MHz', fontsize=16)
        
        # Filter successful fits
        success_mask = results_df['fit_success'] == True
        good_data = results_df[success_mask]

        if len(good_data) == 0:
                print("No successful fits to plot!")
                return

        residues = good_data.index if 'Residue' not in good_data.columns else good_data['Residue']

        # Detect missing residues for grey bar plotting
        import re
        def extract_numeric(res_id):
                match = re.search(r'\d+', str(res_id))
                return int(match.group()) if match else 0

        residue_numbers = [extract_numeric(r) for r in residues]
        if len(residue_numbers) > 0:
                max_res = max(residue_numbers)
                all_residues_range = set(range(1, max_res + 1))  # Always start from residue 1
                present_residues = set(residue_numbers)
                missing_residues = sorted(all_residues_range - present_residues)
        else:
                missing_residues = []

        # Helper function to add grey bars for missing residues
        def add_missing_residue_bars(ax, missing_residues):
                """Add grey vertical bars for missing residues spanning full y-axis range"""
                if len(missing_residues) == 0:
                        return
                y_min, y_max = ax.get_ylim()
                for res_num in missing_residues:
                        ax.axvspan(res_num - 0.5, res_num + 0.5,
                                  facecolor='lightgrey', alpha=0.3, zorder=0)

        # Row 1: Experimental data comparison (R1, R2, hetNOE)
        # R1 comparison
        axes[0,0].set_ylim(0, 2)
        add_missing_residue_bars(axes[0,0], missing_residues)
        axes[0,0].errorbar(residue_numbers, good_data['R1_f1'], yerr=good_data['R1_f1_err'],
                                          fmt='o', capsize=3, label=f'{self.field1_freq} MHz', alpha=0.7)
        axes[0,0].errorbar(residue_numbers, good_data['R1_f2'], yerr=good_data['R1_f2_err'],
                                          fmt='s', capsize=3, label=f'{self.field2_freq} MHz', alpha=0.7)
        axes[0,0].set_title('R1 Comparison')
        axes[0,0].set_ylabel('R1 (s⁻¹)')
        axes[0,0].legend()

        # R2 comparison
        axes[0,1].set_ylim(0, 30)
        add_missing_residue_bars(axes[0,1], missing_residues)
        axes[0,1].errorbar(residue_numbers, good_data['R2_f1'], yerr=good_data['R2_f1_err'],
                                          fmt='o', capsize=3, label=f'{self.field1_freq} MHz', alpha=0.7)
        axes[0,1].errorbar(residue_numbers, good_data['R2_f2'], yerr=good_data['R2_f2_err'],
                                          fmt='s', capsize=3, label=f'{self.field2_freq} MHz', alpha=0.7)
        axes[0,1].set_title('R2 Comparison')
        axes[0,1].set_ylabel('R2 (s⁻¹)')
        axes[0,1].legend()

        # hetNOE comparison
        axes[0,2].set_ylim(0, 1)
        add_missing_residue_bars(axes[0,2], missing_residues)
        axes[0,2].errorbar(residue_numbers, good_data['hetNOE_f1'], yerr=good_data['hetNOE_f1_err'],
                                          fmt='o', capsize=3, label=f'{self.field1_freq} MHz', alpha=0.7)
        axes[0,2].errorbar(residue_numbers, good_data['hetNOE_f2'], yerr=good_data['hetNOE_f2_err'],
                                          fmt='s', capsize=3, label=f'{self.field2_freq} MHz', alpha=0.7)
        axes[0,2].set_title('hetNOE Comparison')
        axes[0,2].set_ylabel('hetNOE')
        axes[0,2].legend()
        
        # Row 2: Spectral densities comparison (J0, JwN, J(0.87ωH))
        # J(0) comparison
        axes[1,0].set_ylim(0, 0.000000008)
        add_missing_residue_bars(axes[1,0], missing_residues)
        axes[1,0].errorbar(residue_numbers, good_data['J0_f1'], yerr=good_data['J0_f1_err'],
                                          fmt='o', capsize=3, label=f'{self.field1_freq} MHz', alpha=0.7)
        axes[1,0].errorbar(residue_numbers, good_data['J0_f2'], yerr=good_data['J0_f2_err'],
                                          fmt='s', capsize=3, label=f'{self.field2_freq} MHz', alpha=0.7)
        axes[1,0].set_title('J(0) Comparison')
        axes[1,0].set_ylabel('J(0) (ns/rad²)')
        axes[1,0].legend()

        # J(ωN) comparison
        axes[1,1].set_ylim(0, 0.0000000004)
        add_missing_residue_bars(axes[1,1], missing_residues)
        axes[1,1].errorbar(residue_numbers, good_data['JwN_f1'], yerr=good_data['JwN_f1_err'],
                                          fmt='o', capsize=3, label=f'{self.field1_freq} MHz', alpha=0.7)
        axes[1,1].errorbar(residue_numbers, good_data['JwN_f2'], yerr=good_data['JwN_f2_err'],
                                          fmt='s', capsize=3, label=f'{self.field2_freq} MHz', alpha=0.7)
        axes[1,1].set_title('J(ωN) Comparison')
        axes[1,1].set_ylabel('J(ωN) (ns/rad²)')
        axes[1,1].legend()

        # J(0.87ωH) comparison
        axes[1,2].set_ylim(0, 0.00000000002)
        add_missing_residue_bars(axes[1,2], missing_residues)
        axes[1,2].errorbar(residue_numbers, good_data['JwH_087_f1'], yerr=good_data['JwH_087_f1_err'],
                                          fmt='o', capsize=3, label=f'{self.field1_freq} MHz', alpha=0.7)
        axes[1,2].errorbar(residue_numbers, good_data['JwH_087_f2'], yerr=good_data['JwH_087_f2_err'],
                                          fmt='s', capsize=3, label=f'{self.field2_freq} MHz', alpha=0.7)
        axes[1,2].set_title('J(0.87ωH) Comparison')
        axes[1,2].set_ylabel('J(0.87ωH) (ns/rad²)')
        axes[1,2].legend()
        
        # Row 3: Model-free parameters (S², τe, Rex comparison)
        # S²
        axes[2,0].set_ylim(0, 1)
        add_missing_residue_bars(axes[2,0], missing_residues)
        axes[2,0].errorbar(residue_numbers, good_data['S2'], yerr=good_data['S2_err'],
                                          fmt='o', capsize=3, color='blue')
        axes[2,0].set_title('Order Parameter S² (Dual-Field, J(0.87ωH))')
        axes[2,0].set_ylabel('S²')

        # τe
        axes[2,1].set_ylim(0, 100)
        add_missing_residue_bars(axes[2,1], missing_residues)
        axes[2,1].errorbar(residue_numbers, good_data['te'], yerr=good_data['te_err'],
                                          fmt='o', capsize=3, color='purple')
        axes[2,1].set_title('Internal Correlation Time τe (Dual-Field, J(0.87ωH))')
        axes[2,1].set_ylabel('τe (ps)')

        # Rex comparison
        axes[2,2].set_ylim(0, 20)
        add_missing_residue_bars(axes[2,2], missing_residues)
        axes[2,2].errorbar(residue_numbers, good_data['Rex_field1'], yerr=good_data['Rex_field1_err'],
                                          fmt='o', capsize=3, label=f'Rex {self.field1_freq} MHz', alpha=0.7)
        axes[2,2].errorbar(residue_numbers, good_data['Rex_field2'], yerr=good_data['Rex_field2_err'],
                                          fmt='s', capsize=3, label=f'Rex {self.field2_freq} MHz', alpha=0.7)
        axes[2,2].set_title('Chemical Exchange Rex Comparison (J(0.87ωH))')
        axes[2,2].set_ylabel('Rex (s⁻¹)')
        axes[2,2].legend()
        
        # Row 4: Field-dependence analysis
        # Rex field dependence
        axes[3,0].scatter(good_data['Rex_field1'], good_data['Rex_field2'], alpha=0.6)
        axes[3,0].plot([0, good_data[['Rex_field1', 'Rex_field2']].max().max()],
                                  [0, good_data[['Rex_field1', 'Rex_field2']].max().max()], 'k--', alpha=0.5)
        axes[3,0].set_xlabel(f'Rex {self.field1_freq} MHz (s⁻¹)')
        axes[3,0].set_ylabel(f'Rex {self.field2_freq} MHz (s⁻¹)')
        axes[3,0].set_title('Rex Field Dependence (J(0.87ωH))')

        # Rex ratio vs field ratio (should be quadratic for chemical exchange)
        # Option 3: Filter by minimum threshold + outlier removal
        field_ratio = (self.field2_freq / self.field1_freq)**2
        min_rex = 0.5  # s⁻¹, minimum Rex threshold for meaningful scaling analysis

        # Filter by minimum threshold at both fields
        mask = (good_data['Rex_field1'] > min_rex) & (good_data['Rex_field2'] > min_rex)
        filtered_rex_data = good_data[mask].copy()

        if len(filtered_rex_data) > 0:
            # Calculate ratios without offset (not needed after filtering)
            rex_ratio = filtered_rex_data['Rex_field2'] / filtered_rex_data['Rex_field1']

            # Remove outliers (ratio > 3× expected or < 0.3× expected)
            ratio_mask = (rex_ratio < 3 * field_ratio) & (rex_ratio > 0.3 * field_ratio)
            final_filtered_data = filtered_rex_data[ratio_mask]
            final_rex_ratio = rex_ratio[ratio_mask]

            # Plot filtered data
            axes[3,1].scatter(np.full(len(final_filtered_data), field_ratio), final_rex_ratio, alpha=0.6)
            axes[3,1].axhline(y=field_ratio, color='red', linestyle='--',
                             label=f'Expected ratio = {field_ratio:.2f}')

            # Add info about filtering
            n_total = len(good_data)
            n_filtered = len(final_filtered_data)
            axes[3,1].text(0.05, 0.95, f'n = {n_filtered}/{n_total}\n(Rex > {min_rex} s⁻¹)',
                          transform=axes[3,1].transAxes, verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=9)
        else:
            # No data passed filtering
            axes[3,1].text(0.5, 0.5, f'No residues with Rex > {min_rex} s⁻¹\nat both fields',
                          transform=axes[3,1].transAxes, ha='center', va='center',
                          bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

        axes[3,1].set_xlabel('Field Ratio² (B₂²/B₁²)')
        axes[3,1].set_ylabel('Rex Ratio (Rex₂/Rex₁)')
        axes[3,1].set_title('Rex Field Scaling (J(0.87ωH))')
        axes[3,1].legend()
        
        # Chi-squared distribution
        axes[3,2].hist(good_data['chi2'], bins=20, alpha=0.7, edgecolor='black')
        axes[3,2].set_xlabel('χ² (Fit Quality)')
        axes[3,2].set_ylabel('Count')
        axes[3,2].set_title('Fit Quality Distribution (J(0.87ωH))')
        
        # Format all subplots
        for i in range(3):  # First 3 rows need residue labels
                for j in range(3):
                        axes[i,j].set_xlabel('Residue')
                        axes[i,j].grid(True, alpha=0.3)

        # Fourth row has different x-labels
        axes[3,0].grid(True, alpha=0.3)
        axes[3,1].grid(True, alpha=0.3)
        axes[3,2].grid(True, alpha=0.3)

        # Set x-axis ticks to show every 10 residues (1, 10, 20, 30, ...)
        if len(residue_numbers) > 0:
                max_res = max(residue_numbers)

                # Create ticks: start at 1, then every 10 (10, 20, 30, ...)
                tick_positions = [1] + list(range(10, max_res + 1, 10))

                # Apply to first 3 rows (plots with residue x-axis)
                for i in range(3):
                        for j in range(3):
                                axes[i,j].set_xticks(tick_positions)
                                axes[i,j].tick_params(axis='x', rotation=0)

        plt.tight_layout()

        if save_plots:
                plt.savefig('dual_field_rsdm_analysis_results_087.pdf', dpi=300, bbox_inches='tight')
                print("Plots saved as 'dual_field_rsdm_analysis_results_087.pdf'")
                plt.close(fig)  # Close figure to free memory
        else:
                plt.close(fig)  # Always close to prevent memory leaks

    def save_dual_field_results(self, results_df, filename='dual_field_detailed_results_087.csv'):
        """
        Save dual-field results with confidence intervals (J(0.87ωH) version)
        
        Parameters:
        -----------
        results_df : pandas.DataFrame
            Results from analyze_dual_field_csv()
        filename : str
            Output filename
        """
        # Check if Monte Carlo results are available
        if 'S2_95CI' in results_df.columns:
            # Format confidence intervals for readability
            results_df['S2_CI'] = results_df['S2_95CI'].apply(
                lambda x: f"[{x[0]:.3f}, {x[1]:.3f}]" if not pd.isna(x).any() else ""
            )
            results_df['te_CI'] = results_df['te_95CI'].apply(
                lambda x: f"[{x[0]:.1f}, {x[1]:.1f}]" if not pd.isna(x).any() else ""
            )
            results_df['Rex_field1_CI'] = results_df['Rex_field1_95CI'].apply(
                lambda x: f"[{x[0]:.2f}, {x[1]:.2f}]" if not pd.isna(x).any() else ""
            )
            results_df['Rex_field2_CI'] = results_df['Rex_field2_95CI'].apply(
                lambda x: f"[{x[0]:.2f}, {x[1]:.2f}]" if not pd.isna(x).any() else ""
            )
            
            # Drop the raw CI arrays
            results_df = results_df.drop(columns=['S2_95CI', 'te_95CI', 'Rex_field1_95CI', 'Rex_field2_95CI'])
        
        results_df.to_csv(filename, index=False)
        print(f"Detailed dual-field results (J(0.87ωH)) saved to '{filename}'")


def main():
    """
    Example usage of the DualFieldSpectralDensityAnalysis class with J(0.87ωH)
    """
    # Initialize analysis with dual-field experimental parameters
    analyzer = DualFieldSpectralDensityAnalysis(
        field1_freq=700.093,    # First field in MHz
        field2_freq=600.133,    # Second field in MHz  
        rNH=1.023e-10,        # N-H bond length in meters
        csaN=-160.0e-6        # 15N CSA in ppm
    )
    
    # Analyze dual-field CSV files
    try:
        print("Running dual-field analysis with J(0.87ωH) and Monte Carlo error propagation...")
        print(f"Using proton frequency correction factor: {OMEGA_H_FACTOR}")
        
        # Specify your two CSV files here
        csv_file1 = 'data_in_WT_700.csv'  # Field 1 data file
        csv_file2 = 'data_in_WT_600.csv'  # Field 2 data file
        
        results = analyzer.analyze_dual_field_csv(
            csv_file1, csv_file2,
            use_monte_carlo_errors=True,
            n_monte_carlo=50
        )
        
        # Check if any data was processed
        if len(results) == 0:
            print("No dual-field data to analyze after filtering!")
            return
        
        # Display results
        print("\nDual-Field Reduced Spectral Density Analysis Results (J(0.87ωH)):")
        print("=" * 80)
        
        # Show key results
        if 'Residue' in results.columns:
            display_cols = ['Residue', 'S2', 'tc', 'te', 'Rex_field1', 'Rex_field2', 'chi2']
        else:
            display_cols = ['Index', 'S2', 'tc', 'te', 'Rex_field1', 'Rex_field2', 'chi2']
        
        print(results[display_cols].round(4))
        
        # Save basic results
        results.to_csv('dual_field_rsdm_results_087.csv', index=False)
        
        # Save detailed results with confidence intervals
        analyzer.save_dual_field_results(results)
        
        # Generate plots
        analyzer.plot_dual_field_results(results)
        
        # Print summary statistics
        success_rate = results['fit_success'].mean() * 100
        print(f"\nDual-Field Model-Free Fitting Results (J(0.87ωH)):")
        print(f"  Fit success rate: {success_rate:.1f}%")
        print(f"  Fields analyzed: {analyzer.field1_freq} MHz and {analyzer.field2_freq} MHz")
        print(f"  Proton frequency factor: {OMEGA_H_FACTOR} (J(0.87ωH))")
        
        if success_rate > 0:
            successful = results[results['fit_success']]
            print(f"  Successfully fitted residues: {len(successful)}")
            print(f"  Average S²: {successful['S2'].mean():.3f} ± {successful['S2'].std():.3f}")
            print(f"  Average τc: {successful['tc'].mean():.1f} ± {successful['tc'].std():.1f} ns")
            print(f"  Average τe: {successful['te'].mean():.1f} ± {successful['te'].std():.1f} ps")
            print(f"  Average Rex (field 1): {successful['Rex_field1'].mean():.1f} ± {successful['Rex_field1'].std():.1f} s⁻¹")
            print(f"  Average Rex (field 2): {successful['Rex_field2'].mean():.1f} ± {successful['Rex_field2'].std():.1f} s⁻¹")
            
            # Analyze Rex field dependence
            rex_f1_mean = successful['Rex_field1'].mean()
            rex_f2_mean = successful['Rex_field2'].mean()
            expected_ratio = (analyzer.field2_freq / analyzer.field1_freq)**2
            observed_ratio = rex_f2_mean / rex_f1_mean if rex_f1_mean > 0 else np.nan
            
            print(f"\nRex Field Dependence Analysis (J(0.87ωH)):")
            print(f"  Expected Rex ratio (B₂²/B₁²): {expected_ratio:.2f}")
            print(f"  Observed Rex ratio: {observed_ratio:.2f}")
            print(f"  Correlation between fields: {successful['Rex_field1'].corr(successful['Rex_field2']):.3f}")
            
            # Report Monte Carlo success if used
            if 'mc_success_rate' in successful.columns:
                avg_mc_success = successful['mc_success_rate'].mean()
                print(f"  Average Monte Carlo success rate: {avg_mc_success*100:.1f}%")
                
            # Identify residues with significant exchange
            high_rex_threshold = 2.0  # s⁻¹
            high_rex_residues = successful[
                (successful['Rex_field1'] > high_rex_threshold) | 
                (successful['Rex_field2'] > high_rex_threshold)
            ]
            
            if len(high_rex_residues) > 0:
                print(f"\nResidues with significant chemical exchange (Rex > {high_rex_threshold} s⁻¹):")
                if 'Residue' in high_rex_residues.columns:
                    exchange_list = high_rex_residues['Residue'].tolist()
                else:
                    exchange_list = high_rex_residues.index.tolist()
                print(f"  {len(exchange_list)} residues: {exchange_list}")
                
        else:
            print("  No successful dual-field model-free fits obtained!")
            
    except FileNotFoundError:
        print("Please provide two CSV files with dual-field relaxation data:")
        print("File 1 (e.g., 'data_600MHz.csv') with columns: R1, R1err, R2, R2err, hetNOE, hetNOEerr")
        print("File 2 (e.g., 'data_800MHz.csv') with columns: R1, R1err, R2, R2err, hetNOE, hetNOEerr")
        print("Both files should also include a 'Residue' column for matching")
        
        # Create example dual-field data
        np.random.seed(42)  # For reproducible example
        n_residues = 15
        
        # Generate realistic dual-field relaxation data
        # Field 1 (600 MHz) - typically lower Rex values
        example_data1 = pd.DataFrame({
            'Residue': [f'A{i}' for i in range(1, n_residues + 1)],
            'R1': np.random.normal(1.3, 0.3, n_residues),
            'R1err': np.random.normal(0.05, 0.01, n_residues),
            'R2': np.random.normal(12.0, 3.0, n_residues), 
            'R2err': np.random.normal(0.4, 0.1, n_residues),
            'hetNOE': np.random.normal(0.75, 0.15, n_residues),
            'hetNOEerr': np.random.normal(0.03, 0.005, n_residues)
        })
        
        # Field 2 (800 MHz) - typically higher Rex values due to field dependence
        field_ratio_squared = (700.0/600.0)**2
        example_data2 = pd.DataFrame({
            'Residue': [f'A{i}' for i in range(1, n_residues + 1)],
            'R1': np.random.normal(1.1, 0.25, n_residues),  # Slightly different R1
            'R1err': np.random.normal(0.04, 0.008, n_residues),
            'R2': np.random.normal(14.0, 4.0, n_residues),  # Higher R2 at higher field
            'R2err': np.random.normal(0.5, 0.12, n_residues),
            'hetNOE': np.random.normal(0.78, 0.12, n_residues),  # Slightly different NOE
            'hetNOEerr': np.random.normal(0.025, 0.004, n_residues)
        })
        
        # Ensure all values are positive and physically reasonable
        for df in [example_data1, example_data2]:
            df['R1'] = np.clip(df['R1'], 0.5, 3.0)
            df['R1err'] = np.clip(df['R1err'], 0.01, 0.1)
            df['R2'] = np.clip(df['R2'], df['R1'] * 1.5, 50.0)  # R2 > R1
            df['R2err'] = np.clip(df['R2err'], 0.1, 1.0)
            df['hetNOE'] = np.clip(df['hetNOE'], 0.3, 1.0)
            df['hetNOEerr'] = np.clip(df['hetNOEerr'], 0.01, 0.05)
        
        example_data1.to_csv('example_600MHz_data_087.csv', index=False)
        example_data2.to_csv('example_800MHz_data_087.csv', index=False)
        
        print("\nExample dual-field data files created for J(0.87ωH) analysis:")
        print("  'example_600MHz_data_087.csv' - Field 1 data")
        print("  'example_800MHz_data_087.csv' - Field 2 data")
        print("You can test the J(0.87ωH) script with these example files.")
        print("\nTo use your own data, ensure both CSV files have the same residue identifiers")
        print("and modify the field frequencies in the main() function as needed.")


if __name__ == "__main__":
    main()