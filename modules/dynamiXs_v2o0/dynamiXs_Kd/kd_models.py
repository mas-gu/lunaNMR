# ABOUTME: 1:1 binding isotherm models for Kd titration fitting (CSP and intensity).
# ABOUTME: Full quadratic fraction-bound (accounts for protein concentration P0).

import numpy as np


def fraction_bound(L, P, Kd):
    """Fraction of protein bound for a 1:1 interaction (full quadratic isotherm).

    Accounts for ligand depletion (free ligand != total) via the exact solution
    of P0*L0 = (Kd + free)·bound. Valid for any [P], [L]; reduces to the
    hyperbola L/(L+Kd) only when L >> P.

    Parameters
    ----------
    L : float or array
        Total ligand concentration.
    P : float
        Total protein concentration (P0).
    Kd : float
        Dissociation constant (same units as L and P).
    """
    L = np.asarray(L, dtype=float)
    s = P + L + Kd
    disc = np.maximum(s * s - 4.0 * P * L, 0.0)
    fb = (s - np.sqrt(disc)) / (2.0 * P)
    return fb if fb.ndim else float(fb)


def csp_model(L, dd_max, Kd, P):
    """Observed CSP for a 1:1 titration: Δδ_obs = Δδ_max · fraction_bound."""
    return dd_max * fraction_bound(L, P, Kd)


def intensity_decay(L, I0, Kd):
    """Intensity decay for a titration: I = I0 · exp(-L / Kd).

    Matches the lab's KD scripts (fitmulti__KD_NMRRE.py: A·exp(-x/t2)). I0 is the
    fitted reference intensity at L=0; Kd is the decay constant (the [L] where I
    falls to I0/e). Fitting I0 anchors the curve to the data's first point instead
    of forcing a binding isotherm onto an intensity-loss profile.
    """
    L = np.asarray(L, dtype=float)
    out = I0 * np.exp(-L / Kd)
    return out if out.ndim else float(out)


def compute_csp(dH, dN, alpha=0.14):
    """Combined chemical-shift perturbation magnitude.

    CSP = sqrt(Δδ_H^2 + (alpha · Δδ_N)^2), where alpha scales the 15N shift to
    1H-equivalent ppm (default 0.14).
    """
    dH = np.asarray(dH, dtype=float)
    dN = np.asarray(dN, dtype=float)
    csp = np.sqrt(dH * dH + (alpha * dN) ** 2)
    return csp if csp.ndim else float(csp)
