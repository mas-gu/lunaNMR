# ABOUTME: Physics shared by the single- and dual-field spectral density analyses: constants, Larmor
# ABOUTME: frequencies, dipolar/CSA constants, reduced spectral densities (Farrow 1995) and predicted relaxation.

"""
One copy of every formula the spectral density and model-free analyses use, so the
single- and dual-field fits cannot disagree on one. d and c are the Farrow (1995)
constants d^2 and c^2: `dipolar_constant` returns d^2, not d^2/4.
"""

import numpy as np
import pandas as pd

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


@jit(cache=True, forceobj=True)
def spectral_density(s2, tc, te, omega):
    """
    J(w) of the Lipari-Szabo model-free form, with the 2/5 prefactor.
    """
    if te == 0:
        j = s2 * tc / (1.0 + (omega * tc)**2)
    else:
        teff = 1.0 / (1.0/tc + 1.0/te)
        j = (s2 * tc / (1.0 + (omega * tc)**2) +
             (1.0 - s2) * teff / (1.0 + (omega * teff)**2))
    return j * 0.4  # 2/5 prefactor


def omega_h(freq_mhz):
    """1H Larmor frequency in rad/s"""
    return freq_mhz * 2 * PI * 1e6


def omega_n(freq_mhz):
    """15N Larmor frequency in rad/s (negative: GAMMA_N is)"""
    return omega_h(freq_mhz) * GAMMA_N / GAMMA_H


def dipolar_constant(r_nh):
    """
    Dipolar coupling constant d², field-independent

    d² as defined in Farrow et al. J. Biomol. NMR, 6 (1995) 153-162
    d² = [μ₀ħγNγH/(4πr³NH)]² = [μ₀hγNγH/(8π²r³NH)]²
    """
    mu0_h_bar = REDUCED_PERM_VACUUM * REDUCED_PLANK
    return (mu0_h_bar * GAMMA_N * GAMMA_H / (4 * PI * r_nh**3))**2


def csa_constant(omega_n_rad, csa_n):
    """CSA constant c² (field-dependent)"""
    return (omega_n_rad * csa_n)**2 / 3.0


def sigma_noe(noe, r1):
    """Cross-relaxation rate from the hetNOE and R1"""
    return (noe - 1.0) * r1 * (GAMMA_N / GAMMA_H)


def reduced_spectral_densities(r1, r2, noe, d, c):
    """
    J(0), J(wN) and J(wH) from R1, R2 and the hetNOE - Farrow et al. (1995) Eqs. 5-7:

    J(0)  = [R₂ - (3d²/8 + c²/2)J(ωN) - (13d²/8)J(ω̃H)] / (d²/2 + 2c²/3)
    J(ωN) = [R₁ - (7d²/4)J(ω̃H)] / [(3d²/4) + c²]
    J(ω̃H) = [4/(5d²)] × (γN/γH) × (NOE - 1) × R₁

    J(0) and J(ωN) are evaluated in their algebraically simplified forms.
    """
    sigma = sigma_noe(noe, r1)
    j0 = (3.0 / (2.0 * (3.0 * d / 4.0 + c))) * (-0.5 * r1 + r2 - (3.0/5.0) * sigma)
    jwn = (1.0 / (3.0 * d / 4.0 + c)) * (r1 - (7.0/5.0) * sigma)
    jwh = 4.0 * sigma / (5.0 * d)
    return j0, jwn, jwh


def reduced_spectral_density_errors(r1, r1_err, r2, r2_err, noe, noe_err, d, c):
    """
    First-order errors of J(0), J(wN) and J(wH): sqrt(sum((dJ/dx * sigma_x)^2)) over
    R1, R2 and the hetNOE.
    """
    gamma_ratio = GAMMA_N / GAMMA_H
    dsigma_dr1 = (noe - 1.0) * gamma_ratio
    dsigma_dnoe = r1 * gamma_ratio

    factor_j0 = 3.0 / (2.0 * (3.0 * d / 4.0 + c))
    dj0_dr1 = factor_j0 * (-0.5 - (3.0/5.0) * dsigma_dr1)
    dj0_dr2 = factor_j0
    dj0_dnoe = factor_j0 * (-(3.0/5.0) * dsigma_dnoe)
    j0_err = np.sqrt((dj0_dr1 * r1_err)**2 + (dj0_dr2 * r2_err)**2 + (dj0_dnoe * noe_err)**2)

    factor_jwn = 1.0 / (3.0 * d / 4.0 + c)
    djwn_dr1 = factor_jwn * (1.0 - (7.0/5.0) * dsigma_dr1)
    djwn_dnoe = factor_jwn * (-(7.0/5.0) * dsigma_dnoe)
    jwn_err = np.sqrt((djwn_dr1 * r1_err)**2 + (djwn_dnoe * noe_err)**2)

    factor_jwh = 4.0 / (5.0 * d)
    djwh_dr1 = factor_jwh * dsigma_dr1
    djwh_dnoe = factor_jwh * dsigma_dnoe
    jwh_err = np.sqrt((djwh_dr1 * r1_err)**2 + (djwh_dnoe * noe_err)**2)

    return j0_err, jwn_err, jwh_err


def relaxation_deviations(s2, te, rex, t1, t2, noe, tc, d, c, omega_h_j, omega_h_rad, omega_n_abs):
    """Relative deviations of the predicted T1, T2 and (when measured) NOE from the
    measured ones, at fixed tc - one field's share of a model-free fit's objective.

    J(wH) is evaluated at `omega_h_j` (0.87 wH in the 0.87 variants); J(wH +/- wN) at the
    true frequencies. The NOE is predicted from the T1 the model predicts, so the
    prediction depends on the parameters alone.
    """
    gammaHN = GAMMA_H / GAMMA_N
    jH = spectral_density(s2, tc, te, omega_h_j)
    jN = spectral_density(s2, tc, te, omega_n_abs)
    jHpN = spectral_density(s2, tc, te, omega_h_rad + omega_n_abs)
    jHmN = spectral_density(s2, tc, te, omega_h_rad - omega_n_abs)
    j0 = spectral_density(s2, tc, te, 0.0)

    r1_calc = (d / 4.0) * (3*jN + 6*jHpN + jHmN) + c * jN
    r2_calc = 0.5 * (d / 4.0) * (4*j0 + 3*jN + 6*jHpN + 6*jH + jHmN) + c * (2*j0/3.0 + 0.5*jN) + rex
    t1_calc = 1.0 / r1_calc

    deviations = [(t1_calc - t1) / t1, (1.0 / r2_calc - t2) / t2]
    if noe is not None and not pd.isna(noe):
        noe_calc = 1.0 + ((d / 4.0) * gammaHN * t1_calc * (6*jHpN - jHmN))
        deviations.append((noe_calc - noe) / (0.1 if noe == 0 else abs(noe)))
    return deviations


def overall_correlation_time(r1_array, r2_array, omega_n_abs):
    """
    Overall correlation time from R1/R2 by the Fushman method,
    tc = [1/(2*omegaN)] * sqrt[(6*T1/T2) - 7], median over the residues where the ratio
    exceeds 7; the mean ratio when fewer than three do, 10 ns when even that is <= 7.
    Constrained to 1 ps - 100 ns.
    """
    t1_array = 1.0 / r1_array
    t2_array = 1.0 / r2_array

    ratio = 6.0 * t1_array / t2_array
    valid_mask = ratio > 7.0

    if np.sum(valid_mask) < 3:
        tc_est = 1.0 / (2.0 * omega_n_abs) * np.sqrt(np.mean(ratio) - 7.0) if np.mean(ratio) > 7 else 10e-9
    else:
        tc_values = 1.0 / (2.0 * omega_n_abs) * np.sqrt(ratio[valid_mask] - 7.0)
        tc_est = np.median(tc_values)  # Use median for robustness

    return max(1e-12, min(100e-9, tc_est))
