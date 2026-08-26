#!/usr/bin/env python3
"""
Reduced Spectral Density Mapping Analysis Script with Multiprocessing Support and J(0.87ωH)

This script imports R1, R1err, R2, R2err, hetNOE, and hetNOEerr from a CSV file
and performs reduced spectral density calculations including J(0), J(wN), J(0.87wH),
S2, Rex, and te parameters with proper error propagation.

MULTICORE VERSION: Modified to utilize multiple CPU cores for parallel processing.
J(0.87ωH) VERSION: Uses J(0.87ωH) instead of J(ωH) for more accurate treatment of
cross-correlation effects and dipolar-CSA interference in the relaxation mechanisms.
"""

#####

import os
os.environ['OMP_NUM_THREADS'] = str(int(os.cpu_count() * 0.8))
os.environ['MKL_NUM_THREADS'] = str(int(os.cpu_count() * 0.8))

#####

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend (must be before pyplot import)
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

# Multiprocessing imports
from multiprocessing import Pool, cpu_count
from functools import partial
import pickle

# Numba JIT compilation for performance (optional but recommended)
try:
    from numba import jit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    # Fallback: decorator that does nothing
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator

# Physical constants
GAMMA_H = 2.67522128e8  # 1H gyromagnetic ratio (rad s-1 T-1)
GAMMA_N = -2.7126e7     # 15N gyromagnetic ratio (rad s-1 T-1) 
REDUCED_PERM_VACUUM = 1.25663706212e-6  # H/m
REDUCED_PLANK = 1.05457180013e-34       # J⋅s
PI = 3.14159265359

# Cross-correlation factor for J(0.87ωH)
OMEGA_H_FACTOR = 0.87  # Factor for more accurate spectral density calculation

# JIT-compiled spectral density calculation for performance
@jit(nopython=False, cache=True, forceobj=True)
def _calculate_spectral_density_jit(s2, tc, te, omega):
    """
    Calculate J(w) using Lipari-Szabo model-free approach (JIT-compiled version)
    """
    if te == 0:
        j = s2 * tc / (1.0 + (omega * tc)**2)
    else:
        teff = 1.0 / (1.0/tc + 1.0/te)
        j = (s2 * tc / (1.0 + (omega * tc)**2) +
             (1.0 - s2) * teff / (1.0 + (omega * teff)**2))
    return j * 0.4  # 2/5 prefactor


class ReducedSpectralDensityAnalysis:
    """
    Class for performing reduced spectral density mapping analysis using J(0.87ωH)
    """
    
    def __init__(self, spectrometer_frequency=600.0, rNH=1.015e-10, csaN=-160.0e-6):
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
        """
        self.spectrometer_frequency = spectrometer_frequency
        self.rNH = rNH
        self.csaN = csaN
        
        # Calculate angular frequencies
        self.omegaH = self._calculate_omegaH()
        self.omegaN = self._calculate_omegaN()
        
        # Calculate effective proton frequency with 0.87 factor
        self.omegaH_eff = self.omegaH * OMEGA_H_FACTOR
        
        # Calculate dipolar and CSA constants
        self.d_factor = self._calculate_d_factor()
        self.c_factor = self._calculate_c_factor()
        
    def _calculate_omegaH(self):
        """Calculate 1H Larmor frequency in rad/s"""
        return self.spectrometer_frequency * 2 * PI * 1e6
    
    def _calculate_omegaN(self):
        """Calculate 15N Larmor frequency in rad/s"""
        return self.omegaH * GAMMA_N / GAMMA_H
    
    def _calculate_d_factor(self):
        """
        Calculate dipolar coupling constant (d²)

        Returns d² as defined in Farrow et al. J. Biomol. NMR, 6 (1995) 153-162
        d² = [μ₀ħγNγH/(8π²r³NH)]²
        """
        mu0_h_bar = REDUCED_PERM_VACUUM * REDUCED_PLANK
        d_squared = (mu0_h_bar * GAMMA_N * GAMMA_H / (4 * PI * self.rNH**3))**2
        return d_squared  # Returns d², not d²/4 (Farrow-exact notation)
    
    def _calculate_c_factor(self):
        """Calculate CSA constant"""
        return (self.omegaN * self.csaN)**2 / 3.0
    
    def calculate_sigma_NOE(self, noe, r1):
        """
        Calculate sigma NOE value
        
        Parameters:
        -----------
        noe : array-like
            Heteronuclear NOE values
        r1 : array-like
            R1 relaxation rates (s-1)
            
        Returns:
        --------
        array : sigma NOE values
        """
        return (noe - 1.0) * r1 * (GAMMA_N / GAMMA_H)
    
    def calculate_J0(self, noe, r1, r2):
        """
        Calculate J(0) spectral density

        Implements Farrow et al. (1995) Equation 7:
        J(0) = [R₂ - (3d²/8 + c²/2)J(ωN) - (13d²/8)J(ω̃H)] / (d²/2 + 2c²/3)

        Parameters:
        -----------
        noe : array-like
            Heteronuclear NOE values
        r1 : array-like
            R1 relaxation rates (s-1)
        r2 : array-like
            R2 relaxation rates (s-1)

        Returns:
        --------
        array : J(0) values
        """
        sigma_noe = self.calculate_sigma_NOE(noe, r1)
        d = self.d_factor  # d = d² (Farrow notation)
        c = self.c_factor  # c = c²

        # Farrow Eq. 7 with algebraically simplified form
        j0 = (3.0 / (2.0 * (3.0 * d / 4.0 + c))) * (
            -0.5 * r1 + r2 - (3.0/5.0) * sigma_noe
        )
        return j0
    
    def calculate_JwN(self, noe, r1, r2):
        """
        Calculate J(wN) spectral density

        Implements Farrow et al. (1995) Equation 6:
        J(ωN) = [R₁ - (7d²/4)J(ω̃H)] / [(3d²/4) + c²]

        Parameters:
        -----------
        noe : array-like
            Heteronuclear NOE values
        r1 : array-like
            R1 relaxation rates (s-1)
        r2 : array-like
            R2 relaxation rates (s-1)

        Returns:
        --------
        array : J(wN) values
        """
        sigma_noe = self.calculate_sigma_NOE(noe, r1)
        d = self.d_factor  # d = d² (Farrow notation)
        c = self.c_factor  # c = c²

        # Farrow Eq. 6: denominator has (3d²/4) + c²
        jwn = (1.0 / (3.0 * d / 4.0 + c)) * (r1 - (7.0/5.0) * sigma_noe)
        return jwn
    
    def calculate_JwH_087(self, noe, r1):
        """
        Calculate J(0.87ωH) spectral density

        Implements Farrow et al. (1995) Equation 5:
        J(0.87ωH) = [4/(5d²)] × (γN/γH) × (NOE - 1) × R₁

        Parameters:
        -----------
        noe : array-like
            Heteronuclear NOE values
        r1 : array-like
            R1 relaxation rates (s-1)

        Returns:
        --------
        array : J(0.87ωH) values
        """
        sigma_noe = self.calculate_sigma_NOE(noe, r1)
        d = self.d_factor  # d = d² (Farrow notation)

        # Farrow Eq. 5: J(0.87ωH) = [4/(5d²)] × σNOE
        jwh_087 = 4.0 * sigma_noe / (5.0 * d)
        return jwh_087
    
    def calculate_isotropic_spectral_density(self, s2, tc, te, omega):
        """
        Calculate J(w) using Lipari-Szabo model-free approach
        Delegates to JIT-compiled function for performance
        """
        return _calculate_spectral_density_jit(s2, tc, te, omega)
    
    def estimate_overall_correlation_time(self, r1_array, r2_array):
        """
        Estimate overall correlation time from R1/R2 data using Fushman method
        
        Parameters:
        -----------
        r1_array : array
            Array of R1 values
        r2_array : array  
            Array of R2 values
            
        Returns:
        --------
        float : Overall correlation time in seconds
        """
        t1_array = 1.0 / r1_array
        t2_array = 1.0 / r2_array
        
        # Fushman method: tc = [1/(2*omegaN)] * sqrt[(6*T1/T2) - 7]
        ratio = 6.0 * t1_array / t2_array
        valid_mask = ratio > 7.0  # Only use residues where ratio > 7
        
        if np.sum(valid_mask) < 3:
            # Fallback: use simple average if not enough valid residues
            tc_est = 1.0 / (2.0 * abs(self.omegaN)) * np.sqrt(np.mean(ratio) - 7.0) if np.mean(ratio) > 7 else 10e-9
        else:
            tc_values = 1.0 / (2.0 * abs(self.omegaN)) * np.sqrt(ratio[valid_mask] - 7.0)
            tc_est = np.median(tc_values)  # Use median for robustness
            
        return max(1e-12, min(100e-9, tc_est))  # Constrain to reasonable range

    def _fit_single_dataset(self, r1, r2, noe, tc_fixed, niter):
        """
        Original fitting function using J(0.87ωH), now separated for clarity
        This is the existing fit_model_free_individual code modified for J(0.87ωH)
        """
        from random import random, randint
        from math import exp
        
        # Convert to T1, T2
        t1 = 1.0 / r1
        t2 = 1.0 / r2
        
        # Constants
        A = self.d_factor
        C = self.c_factor  
        gammaHN = GAMMA_H / GAMMA_N
        
        # Initialize parameter ensemble
        ensemble = []
        
        # Parameter starting values - RESTORED to original for accuracy
        s2_start = [x * 0.05 for x in range(1, 21)]  # 0.05 to 1.0 (20 values)
        te_start = [x * 1e-12 for x in [0.1, 0.5, 1, 2, 5, 10, 20, 50]]  # 8 values
        rex_start = [float(x) for x in range(16)]  # 0 to 15 s-1 (16 values)
        
        # Create initial ensemble
        for s2 in s2_start:
            for te in te_start:
                for rex in rex_start:
                    ensemble.append((None, s2, te, rex))
        
        # Calculate initial scores
        for k, (score, s2, te, rex) in enumerate(ensemble):
            # Calculate spectral densities using J(0.87ωH)
            jH = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, self.omegaH_eff)
            jN = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, abs(self.omegaN))
            jHpN = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, self.omegaH + abs(self.omegaN))
            jHmN = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, self.omegaH - abs(self.omegaN))
            j0 = self.calculate_isotropic_spectral_density(s2, tc_fixed, te, 0.0)
            
            # Calculate predicted relaxation parameters
            r1_calc = (A / 4.0) * (3*jN + 6*jHpN + jHmN) + C * jN
            r2_calc = 0.5 * (A / 4.0) * (4*j0 + 3*jN + 6*jHpN + 6*jH + jHmN) + C * (2*j0/3.0 + 0.5*jN) + rex
            
            t1_calc = 1.0 / r1_calc
            t2_calc = 1.0 / r2_calc
            
            # Calculate score (chi-squared)
            d1 = (t1_calc - t1) / t1
            d2 = (t2_calc - t2) / t2
            score = d1*d1 + d2*d2
            
            # Add NOE contribution if available
            if noe is not None and not pd.isna(noe):
                noe_calc = 1.0 + ((A / 4.0) * gammaHN * t1 * (6*jHpN - jHmN))
                dn = (noe_calc - noe) / (0.1 if noe == 0 else abs(noe))  # Avoid division by zero
                score += dn*dn
                
            ensemble[k] = (score, s2, te, rex, t1_calc, t2_calc)
        
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

            prevScore, s2, te, rex, t1_calc, t2_calc = ensemble[-1]
            
            # Mutate parameters
            d = ((random() - 0.382) * f) + 1.0
            s2_new = max(0.01, min(1.0, s2 * d))
            
            d = ((random() - 0.382) * f) + 1.0  
            te_new = max(1e-15, min(10e-10, te * d)) #GM 10e-9 to 10e-10
            
            d = ((random() - 0.382) * f) + 1.0
            rex_new = max(0.0, min(50.0, rex * d))
            
            # Calculate new score using J(0.87ωH)
            jH = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, self.omegaH_eff)
            jN = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, abs(self.omegaN))
            jHpN = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, self.omegaH + abs(self.omegaN))
            jHmN = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, self.omegaH - abs(self.omegaN))
            j0 = self.calculate_isotropic_spectral_density(s2_new, tc_fixed, te_new, 0.0)
            
            r1_calc = (A / 4.0) * (3*jN + 6*jHpN + jHmN) + C * jN
            r2_calc = 0.5 * (A / 4.0) * (4*j0 + 3*jN + 6*jHpN + 6*jH + jHmN) + C * (2*j0/3.0 + 0.5*jN) + rex_new
            
            t1_calc = 1.0 / r1_calc
            t2_calc = 1.0 / r2_calc
            
            d1 = (t1_calc - t1) / t1
            d2 = (t2_calc - t2) / t2
            score = d1*d1 + d2*d2
            
            if noe is not None and not pd.isna(noe):
                noe_calc = 1.0 + ((A / 4.0) * gammaHN * t1 * (6*jHpN - jHmN))
                dn = (noe_calc - noe) / (0.1 if noe == 0 else abs(noe))
                score += dn*dn
            
            # Accept or reject
            ratio = exp(prevScore - score)
            if ratio > 1.0:
                ensemble[-1] = (score, s2_new, te_new, rex_new, t1_calc, t2_calc)
            else:
                k = randint(0, ensemble_size - 1)
                ensemble[-1] = ensemble[k]
        
        # Return best result
        ensemble.sort()
        if len(ensemble) > 0:
            score, s2_best, te_best, rex_best, t1_calc, t2_calc = ensemble[0]
            
            # Calculate errors from ensemble spread
            s2_values = [e[1] for e in ensemble[:100]]  # Top 100 results
            te_values = [e[2] for e in ensemble[:10]]
            rex_values = [e[3] for e in ensemble[:100]]
            
            return {
                'S2': s2_best,
                'tc': tc_fixed * 1e9,  # Convert to ns
                'te': te_best * 1e12,  # Convert to ps
                'Rex': rex_best,
                'S2_err': np.std(s2_values),
                'te_err': np.std(te_values) * 1e12,
                'Rex_err': np.std(rex_values), 
                'fit_success': True,
                'chi2': score,
                'iterations': i
            }
        else:
            return {
                'S2': np.nan, 'tc': np.nan, 'te': np.nan, 'Rex': np.nan,
                'S2_err': np.nan, 'te_err': np.nan, 'Rex_err': np.nan,
                'fit_success': False, 'chi2': np.inf, 'iterations': niter
            }

    def fit_model_free_individual(self, r1, r2, noe, tc_fixed,
                                r1_err=None, r2_err=None, noe_err=None,
                                n_monte_carlo=100, niter=10_000):
        """
        Fit individual residue model-free parameters using fixed tc and J(0.87ωH)
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
            
            # Ensure physical validity
            r1_sample = max(0.1, r1_sample)  # R1 must be positive
            r2_sample = max(r1_sample * 1.1, r2_sample)  # R2 > R1
            # NOE can be negative but typically between -0.5 and 1.0
            noe_sample = np.clip(noe_sample, -0.5, 1.0)
            
            # Fit the synthetic dataset with reduced iterations
            mc_fit = self._fit_single_dataset(r1_sample, r2_sample, noe_sample, 
                                            tc_fixed, mc_niter)
            
            if mc_fit['fit_success']:
                mc_results['S2'].append(mc_fit['S2'])
                mc_results['te'].append(mc_fit['te']) 
                mc_results['Rex'].append(mc_fit['Rex'])
                mc_results['chi2'].append(mc_fit['chi2'])
        
        print(" done")
        
        # Calculate errors from Monte Carlo distribution
        n_successful = len(mc_results['S2'])
        if n_successful < 10:  # Need minimum samples for reliable errors
            print(f"    Warning: Only {n_successful}/{n_monte_carlo} MC fits succeeded")
            # Fall back to ensemble-based errors
            return best_fit
        
        # Update errors with Monte Carlo results
        best_fit['S2_err'] = np.std(mc_results['S2'])
        best_fit['te_err'] = np.std(mc_results['te'])
        best_fit['Rex_err'] = np.std(mc_results['Rex'])
        
        # Add confidence intervals
        best_fit['S2_95CI'] = np.percentile(mc_results['S2'], [2.5, 97.5])
        best_fit['te_95CI'] = np.percentile(mc_results['te'], [2.5, 97.5])
        best_fit['Rex_95CI'] = np.percentile(mc_results['Rex'], [2.5, 97.5])
        
        # Add MC success rate
        best_fit['mc_success_rate'] = n_successful / n_monte_carlo
        
        return best_fit

    
    def propagate_errors(self, r1, r1_err, r2, r2_err, noe, noe_err):
        """
        Propagate errors through spectral density calculations using J(0.87ωH)
        
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
        # Calculate partial derivatives for error propagation
        d = self.d_factor
        c = self.c_factor
        gamma_ratio = GAMMA_N / GAMMA_H
        
        # Sigma NOE and its error
        sigma_noe = (noe - 1.0) * r1 * gamma_ratio
        dsigma_dr1 = (noe - 1.0) * gamma_ratio
        dsigma_dnoe = r1 * gamma_ratio
        sigma_noe_err = np.sqrt((dsigma_dr1 * r1_err)**2 + (dsigma_dnoe * noe_err)**2)
        
        # J(0) and its error (Farrow-exact notation with d²/4)
        j0 = (3.0 / (2.0 * (3.0 * d / 4.0 + c))) * (-0.5 * r1 + r2 - (3.0/5.0) * sigma_noe)
        factor_j0 = 3.0 / (2.0 * (3.0 * d / 4.0 + c))
        dj0_dr1 = factor_j0 * (-0.5 - (3.0/5.0) * dsigma_dr1)
        dj0_dr2 = factor_j0
        dj0_dnoe = factor_j0 * (-(3.0/5.0) * dsigma_dnoe)
        j0_err = np.sqrt((dj0_dr1 * r1_err)**2 + (dj0_dr2 * r2_err)**2 + (dj0_dnoe * noe_err)**2)
        
        # J(wN) and its error (Farrow-exact notation with d²/4)
        jwn = (1.0 / (3.0 * d / 4.0 + c)) * (r1 - (7.0/5.0) * sigma_noe)
        factor_jwn = 1.0 / (3.0 * d / 4.0 + c)
        djwn_dr1 = factor_jwn * (1.0 - (7.0/5.0) * dsigma_dr1)
        djwn_dnoe = factor_jwn * (-(7.0/5.0) * dsigma_dnoe)
        jwn_err = (np.sqrt((djwn_dr1 * r1_err)**2 + (djwn_dnoe * noe_err)**2))/2 ## /2 GM 
        
        # J(0.87ωH) and its error (Farrow-exact notation with factor 4)
        jwh_087 = 4.0 * sigma_noe / (5.0 * d)
        factor_jwh = 4.0 / (5.0 * d)
        djwh_dr1 = factor_jwh * dsigma_dr1
        djwh_dnoe = factor_jwh * dsigma_dnoe
        jwh_087_err = np.sqrt((djwh_dr1 * r1_err)**2 + (djwh_dnoe * noe_err)**2)
        
        return {
            'J0': j0, 'J0_err': j0_err,
            'JwN': jwn, 'JwN_err': jwn_err,
            'JwH_087': jwh_087, 'JwH_087_err': jwh_087_err
        }
    
    def analyze_csv(self, csv_file, residue_col='Residue', use_monte_carlo_errors=True,
                    n_monte_carlo=100, use_multiprocessing=True, n_processes=None):
        """
        Analyze relaxation data from CSV file with multiprocessing support and J(0.87ωH)
        
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
        # Load data
        data = pd.read_csv(csv_file)
        
        # Check required columns
        required_cols = ['R1', 'R1err', 'R2', 'R2err', 'hetNOE', 'hetNOEerr']
        missing_cols = [col for col in required_cols if col not in data.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # Determine number of processes
        if n_processes is None:
            n_processes = max(1, int(cpu_count() * 0.8))
        
        print(f"\nProcessing {len(data)} residues from {csv_file}...")
        print(f"Using J(0.87ωH) approach with factor = {OMEGA_H_FACTOR}")
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
        
        # Estimate overall correlation time
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
        Process residues sequentially (original behavior) using J(0.87ωH)
        """
        results = []
        
        # Counters for data quality reporting
        total_residues = len(data)
        skipped_nan = 0
        skipped_zero_negative = 0
        processed = 0
        
        # Process individual residues
        for idx, row in data.iterrows():
            # Extract data
            r1 = row['R1']
            r1_err = row['R1err'] 
            r2 = row['R2']
            r2_err = row['R2err']
            noe = row['hetNOE']
            noe_err = row['hetNOEerr']
            
            # Get residue identifier for reporting
            residue_id = row[residue_col] if residue_col in data.columns else f"Index_{idx}"
            
            # Skip if any values are NaN
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
                
            # Calculate spectral densities with error propagation using J(0.87ωH)
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
            
            # Add residue info if available
            if residue_col in data.columns:
                result_row['Residue'] = row[residue_col]
                
            results.append(result_row)
            processed += 1
        
        # Report data quality summary
        print(f"\nData Quality Summary (J(0.87ωH) Analysis):")
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
        Process residues in parallel using multiprocessing with J(0.87ωH)
        """
        # Prepare analyzer state for serialization
        analyzer_state = {
            'spectrometer_frequency': self.spectrometer_frequency,
            'rNH': self.rNH,
            'csaN': self.csaN
        }
        
        # Prepare arguments for worker processes
        worker_args = []
        for idx, row in data.iterrows():
            worker_args.append((
                idx, row.to_dict(), analyzer_state, tc_overall, 
                use_monte_carlo_errors, n_monte_carlo, residue_col
            ))
        
        # Process residues in parallel
        print("Starting parallel processing...")
        with Pool(processes=n_processes) as pool:
            worker_results = pool.map(_worker_process_residue_single_087, worker_args)
        
        # Collect results and count statistics
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
        
        print(f"\nData Quality Summary (J(0.87ωH) Analysis):")
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
        Generate plots of the analysis results using J(0.87ωH)
        
        Parameters:
        -----------
        results_df : pandas.DataFrame
            Results from analyze_csv()
        save_plots : bool
            Whether to save plots to files
        """
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle('Reduced Spectral Density Analysis Results (J(0.87ωH))', fontsize=16)
        
        # Filter successful fits
        success_mask = results_df['fit_success'] == True
        good_data = results_df[success_mask]
        
        if len(good_data) == 0:
            print("No successful fits to plot!")
            return
        
        # Plot spectral densities
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
        
        # J(0.87ωH)
        axes[0,2].errorbar(residues, good_data['JwH_087'], yerr=good_data['JwH_087_err'],
                          fmt='o', capsize=3, color='green') 
        axes[0,2].set_title('J(0.87ωH)')
        axes[0,2].set_ylabel('J(0.87ωH) (ns/rad²)')
        
        # S²
        axes[1,0].errorbar(residues, good_data['S2'], yerr=good_data['S2_err'],
                          fmt='o', capsize=3, color='blue')         
        axes[1,0].set_title('Order Parameter S² (J(0.87ωH))')
        axes[1,0].set_ylabel('S²')
        axes[1,0].set_ylim(0, 1)
        
        # τe  
        axes[1,1].errorbar(residues, good_data['te'], yerr=good_data['te_err'],
                          fmt='o', capsize=3, color='purple')
        axes[1,1].set_title('Internal Correlation Time τe (J(0.87ωH))')
        axes[1,1].set_ylim(0, 100) 
        axes[1,1].set_ylabel('τe (ps)')
        
        # Rex
        axes[1,2].errorbar(residues, good_data['Rex'], yerr=good_data['Rex_err'],
                          fmt='o', capsize=3, color='brown')
        axes[1,2].set_title('Chemical Exchange Rex (J(0.87ωH))')
        axes[1,2].set_ylabel('Rex (s⁻¹)')
        
        # Format all subplots
        for ax in axes.flat:
            ax.set_xlabel('Residue')
            ax.grid(True, alpha=0.3)
            
        plt.tight_layout()

        if save_plots:
            filename = plot_filename or 'ZZ_rsdm_multicore_analysis_results_087.pdf'
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Plots saved as '{filename}'")
            plt.close(fig)  # Close figure to free memory
        else:
            plt.close(fig)  # Always close to prevent memory leaks

    def save_detailed_results(self, results_df, filename='ZZ_multicore_detailed_results_087.csv'):
        """
        Save results with confidence intervals if Monte Carlo was used (J(0.87ωH) version)
        
        Parameters:
        -----------
        results_df : pandas.DataFrame
            Results from analyze_csv()
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
            results_df['Rex_CI'] = results_df['Rex_95CI'].apply(
                lambda x: f"[{x[0]:.2f}, {x[1]:.2f}]" if not pd.isna(x).any() else ""
            )
            
            # Drop the raw CI arrays
            results_df = results_df.drop(columns=['S2_95CI', 'te_95CI', 'Rex_95CI'])
        
        results_df.to_csv(filename, index=False)
        print(f"Detailed results (J(0.87ωH)) saved to '{filename}'")


# Worker functions for multiprocessing (must be defined at module level)
def _worker_process_residue_single_087(args):
    """
    Worker function for multiprocessing - processes a single residue using J(0.87ωH)
    This function is defined at module level to be pickle-able by multiprocessing
    """
    (idx, residue_data, analyzer_state, tc_overall, use_monte_carlo_errors, 
     n_monte_carlo, residue_col) = args
    
    try:
        # Reconstruct analyzer object from state
        analyzer = ReducedSpectralDensityAnalysis(
            spectrometer_frequency=analyzer_state['spectrometer_frequency'],
            rNH=analyzer_state['rNH'],
            csaN=analyzer_state['csaN']
        )
        
        # Extract data
        r1 = residue_data['R1']
        r1_err = residue_data['R1err']
        r2 = residue_data['R2']
        r2_err = residue_data['R2err']
        noe = residue_data['hetNOE']
        noe_err = residue_data['hetNOEerr']
        
        # Get residue identifier
        residue_id = residue_data[residue_col] if residue_col in residue_data else f"Index_{idx}"
        
        # Check for NaN values
        if any(pd.isna([r1, r1_err, r2, r2_err, noe, noe_err])):
            return {
                'status': 'skipped_nan',
                'residue_id': residue_id,
                'result': None,
                'warnings': [],
                'error': None
            }
        
        # Check for invalid values
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
        
        # Calculate spectral densities with error propagation using J(0.87ωH)
        j_results = analyzer.propagate_errors(r1, r1_err, r2, r2_err, noe, noe_err)
        
        # Fit model-free parameters
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
        
        # Add residue info if available
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
    Example usage of the ReducedSpectralDensityAnalysis class with multiprocessing and J(0.87ωH)
    """
    # Initialize analysis with experimental parameters
    analyzer = ReducedSpectralDensityAnalysis(
        spectrometer_frequency=700.0,  # MHz
        rNH=1.023e-10,                   # meters
        csaN=-160.0e-6                   # ppm in frequency units
    )
    
    # Analyze CSV file with multiprocessing support
    try:
        print("Running analysis with multiprocessing, J(0.87ωH), and Monte Carlo error propagation...")
        print(f"Using proton frequency correction factor: {OMEGA_H_FACTOR}")
        
        results = analyzer.analyze_csv('data_in_700.csv',
                                     use_monte_carlo_errors=True,
                                     n_monte_carlo=100,
                                     use_multiprocessing=True,
                                     n_processes=None)  # Use 80% of available cores
        
        # Check if any data was processed
        if len(results) == 0:
            print("No data to analyze after filtering!")
            return
        
        # Display results
        print("\nReduced Spectral Density Analysis Results (Multiprocessing + J(0.87ωH)):")
        print("=" * 80)
        
        # Show results with residue info if available
        if 'Residue' in results.columns:
            display_cols = ['Residue', 'J0', 'JwN', 'JwH_087', 'S2', 'tc', 'te', 'Rex']
        else:
            display_cols = ['Index', 'J0', 'JwN', 'JwH_087', 'S2', 'tc', 'te', 'Rex']
        
        print(results[display_cols].round(4))
        
        # Save basic results
        results.to_csv('ZZ_rsdm_multicore_087.csv', index=False)
        
        # Save detailed results with confidence intervals
        analyzer.save_detailed_results(results)
        
        # Generate plots
        analyzer.plot_results(results)
        
        # Print summary statistics
        success_rate = results['fit_success'].mean() * 100
        print(f"\nModel-Free Fitting Results (J(0.87ωH)):")
        print(f"  Fit success rate: {success_rate:.1f}%")
        print(f"  Multiprocessing enabled: {success_rate > 0}")
        print(f"  Proton frequency factor: {OMEGA_H_FACTOR} (J(0.87ωH))")
        
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
        
        # Create example data
        example_data = pd.DataFrame({
            'Residue': [f'A{i}' for i in range(1, 11)],
            'R1': np.random.normal(1.5, 0.2, 10),
            'R1err': np.random.normal(0.05, 0.01, 10),
            'R2': np.random.normal(15.0, 2.0, 10), 
            'R2err': np.random.normal(0.5, 0.1, 10),
            'hetNOE': np.random.normal(0.8, 0.1, 10),
            'hetNOEerr': np.random.normal(0.03, 0.005, 10)
        })
        example_data.to_csv('example_relaxation_data_087.csv', index=False)
        print("\nExample data file 'example_relaxation_data_087.csv' created for J(0.87ωH) analysis.")
        print("You can test the J(0.87ωH) script with this example data.")


if __name__ == "__main__":
    main()