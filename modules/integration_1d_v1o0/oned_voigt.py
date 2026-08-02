# ABOUTME: Pure 1D Voigt maths for the 1D integration module - profile, height, FWHM, R2.
# ABOUTME: Depends only on numpy and scipy; no lunaNMR imports, so it is testable in isolation.

import numpy as np
from scipy.special import wofz

# Smallest Gaussian width the profile will evaluate. sigma reaches zero whenever an
# optimiser probes the pure-Lorentzian edge of the parameter space, and the profile
# divides by sigma twice.
SIGMA_FLOOR = 1e-9

# FWHM of a Gaussian of unit sigma: 2*sqrt(2*ln2).
GAUSS_FWHM_FACTOR = 2.0 * np.sqrt(2.0 * np.log(2.0))

# Olivero & Longbothum pseudo-Voigt width coefficients (accurate to ~0.02%).
_OLIVERO_LINEAR = 0.5346
_OLIVERO_QUADRATIC = 0.2166


def voigt_1d(x, area, center, sigma, gamma, baseline=0.0):
    """Area-normalised Voigt profile evaluated on `x`.

    `area` is the integrated intensity of the profile, i.e. the quantity reported
    as the peak volume - not the peak maximum. Use voigt_height() for the maximum.

    sigma is the Gaussian standard deviation, gamma the Lorentzian half-width at
    half-maximum, both in the units of `x` (ppm).
    """
    sigma = max(float(sigma), SIGMA_FLOOR)
    x = np.asarray(x, dtype=float)

    z = ((x - center) + 1j * gamma) / (sigma * np.sqrt(2.0))
    profile = area * np.real(wofz(z)) / (sigma * np.sqrt(2.0 * np.pi))

    return profile + baseline


def voigt_height(area, sigma, gamma):
    """Peak maximum of an area-normalised Voigt, evaluated analytically at the centre.

    The Gaussian relation height = area/(sigma*sqrt(2pi)) holds only at gamma=0;
    any Lorentzian content lowers the maximum by Re(wofz(i*gamma/(sigma*sqrt2))).
    Seeding an `area` from a measured height must invert this function, not the
    Gaussian shortcut, or the seed is low by a factor of ~2 at typical widths.
    """
    sigma = max(float(sigma), SIGMA_FLOOR)

    z = 1j * gamma / (sigma * np.sqrt(2.0))

    return area * np.real(wofz(z)) / (sigma * np.sqrt(2.0 * np.pi))


def voigt_fwhm(sigma, gamma):
    """Full width at half maximum of a Voigt profile, in the units of the axis.

    Olivero-Longbothum approximation combining the Gaussian and Lorentzian widths;
    exact in both limits and within ~0.02% between them.
    """
    gaussian_fwhm = GAUSS_FWHM_FACTOR * float(sigma)
    lorentzian_fwhm = 2.0 * float(gamma)

    return (_OLIVERO_LINEAR * lorentzian_fwhm
            + np.sqrt(_OLIVERO_QUADRATIC * lorentzian_fwhm ** 2 + gaussian_fwhm ** 2))


def r_squared(y, y_fit):
    """Coefficient of determination of `y_fit` against observed `y`.

    Flat data has no variance to explain, so the usual ratio is 0/0; an exact
    reproduction of it scores 1.0 and anything else 0.0.
    """
    y = np.asarray(y, dtype=float)
    y_fit = np.asarray(y_fit, dtype=float)

    ss_res = float(np.sum((y - y_fit) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))

    if ss_tot == 0.0:
        return 1.0 if ss_res == 0.0 else 0.0

    return 1.0 - ss_res / ss_tot
