# ABOUTME: Bounded least-squares search for model-free parameters, from a grid of starts.
# ABOUTME: Shared by every model-free variant; each one supplies its own residuals.

import itertools

import numpy as np
from scipy.optimize import least_squares

# S2, te (s) and Rex (s-1) starts. Rex starts at 0 AND above it, and least squares can move
# a parameter off a bound, so a fit is not held at Rex = 0 by where it began.
S2_STARTS = (0.3, 0.6, 0.85, 0.95)
TE_STARTS = (5e-12, 5e-11, 3e-10)
REX_STARTS = (0.0, 1.0, 5.0)
SINGLE_FIELD_STARTS = tuple(itertools.product(S2_STARTS, TE_STARTS, REX_STARTS))
DUAL_FIELD_STARTS = tuple((s2, te, rex, rex * ratio)
                          for s2, te, rex in SINGLE_FIELD_STARTS for ratio in (1.0, 1.36))

S2_BOUNDS = (0.01, 1.0)
TE_BOUNDS = (1e-15, 10e-10)
REX_BOUNDS = (0.0, 50.0)

# te is ~1e-11 s; the search works on parameters of order one.
TE_SCALE = 1e-11


def best_least_squares(residuals, starts, lower, upper, scale, max_nfev):
    """Minimise sum(residuals(p)**2) within [lower, upper] from every start, keeping the
    lowest. Deterministic: the same inputs give the same answer.

    Returns {'params', 'chi2', 'converged', 'nfev'} for the best start. `converged` is
    False when that start stopped on max_nfev rather than on a convergence test.
    """
    lo, hi, scale = (np.asarray(v, dtype=float) for v in (lower, upper, scale))
    best = None
    for start in starts:
        x0 = np.clip(np.asarray(start, dtype=float), lo, hi) / scale
        sol = least_squares(lambda q: residuals(q * scale), x0,
                            bounds=(lo / scale, hi / scale), max_nfev=max_nfev)
        if best is None or sol.cost < best.cost:
            best = sol
    return {'params': best.x * scale, 'chi2': 2.0 * best.cost,
            'converged': bool(best.status > 0), 'nfev': int(best.nfev)}
