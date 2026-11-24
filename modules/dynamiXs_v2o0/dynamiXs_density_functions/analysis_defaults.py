# ABOUTME: What a spectral density or model-free analysis starts with when the user states nothing:
# ABOUTME: the one declaration the DynamiXs pages, their workers, the pipeline and the CLI all read.

"""
Nothing is imported here, so any front end can read these without loading an analysis.
"""

# 15N chemical shift anisotropy, in ppm
CSA_N_PPM = -172.0

# J at the proton frequency is evaluated at this multiple of wH: 0.87 for J(0.87 wH)
OMEGA_H_FACTOR = 0.87
# That choice as the front ends spell it: a yes/no, and the ending of a method name
USE_087 = OMEGA_H_FACTOR == 0.87
METHOD_VARIANT = '087' if USE_087 else 'jwh'

# Monte-Carlo samples per residue for the errors of a spectral density analysis, and
# of the model-free pipeline, which runs one for every residue it then fits
DENSITY_MONTE_CARLO_SAMPLES = 1000
MODELFREE_MONTE_CARLO_SAMPLES = 50
