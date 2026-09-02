# ABOUTME: The seed a peak's intensity is built from, and what the fit says about its bounds.
# ABOUTME: Real synthetic planes and the real integrator - no mocks anywhere.

"""
Two things the 2D path could not say, and one it decided badly.

`initial_intensity` falls back to `max(sample * 0.8, 100.0)` when the picked point
carries no positive value. 100.0 is an absolute constant in the spectrum's own
arbitrary units, and the whole CLUSTER's intensity interval is derived from that seed
(`[max_initial * 0.001, max_initial * 5.0]`), so a floored seed puts the answer outside
the box the optimiser may search. It then converges immediately against a bound, which
`formal_convergence` reports as success.

Measured on a real titration before this was fixed: 4 of 69 picked points hold a
non-positive sample at the reference and 13 of 69 at the 2:1 point; one cluster had
both its members floored, giving normalised bounds nine orders of magnitude from the
truth at R-squared 0.0208.
"""

import numpy as np
import pytest

from lunaNMR.core.ps2d_config import set_ps2d_config
from lunaNMR.core.core_integrator import VoigtIntegrator
from lunaNMR.core.ps2d_style_fitter import voigt_profile_2d


LW_LOR_F1, LW_GAU_F1 = 0.18, 0.06
LW_LOR_F2, LW_GAU_F2 = 0.018, 0.006

# Amplitudes at the scale a real spectrum carries. At amplitude 1.0 the engine's
# absolute 100.0 seed floor is LARGER than the data and the defect inverts, so a test
# written at unit amplitude passes whatever the code does.
AMPLITUDE = 5.0e8


@pytest.fixture(autouse=True)
def nucleus():
    set_ps2d_config('15N')


def _integrator(peaks, dip_at=None):
    """A 15N-HSQC-like patch, optionally with one sample driven negative."""
    ppm_f2 = np.linspace(7.80, 7.60, 128)
    ppm_f1 = np.linspace(122.4, 120.6, 96)
    f2_grid, f1_grid = np.meshgrid(ppm_f2, ppm_f1)
    data = np.zeros_like(f1_grid)
    for f1, f2, amplitude in peaks:
        data += voigt_profile_2d(f1_grid, f2_grid, f1, LW_GAU_F1, f2,
                                 LW_GAU_F2, LW_LOR_F1, LW_LOR_F2, amplitude)
    if dip_at is not None:
        f1, f2 = dip_at
        yi = int(np.argmin(np.abs(ppm_f1 - f1)))
        xi = int(np.argmin(np.abs(ppm_f2 - f2)))
        data[yi, xi] = -0.5 * abs(data[yi, xi])

    integrator = VoigtIntegrator()
    integrator.nmr_data = data
    integrator.ppm_x_axis = ppm_f2
    integrator.ppm_y_axis = ppm_f1
    return integrator


def test_a_negative_sample_at_the_picked_point_does_not_decide_the_answer():
    """A noise dip, a dispersive wiggle or a coordinate a few ppb out is enough to make
    that one sample non-positive. The peak either side of it is untouched."""
    truth = (121.50, 7.700, AMPLITUDE)
    clean = _integrator([truth])
    dipped = _integrator([truth], dip_at=(truth[0], truth[1]))

    group = [{'x_ppm': truth[1], 'y_ppm': truth[0]}]
    good = clean.fit_overlap_group_2d(group, assignment='clean')
    poor = dipped.fit_overlap_group_2d(group, assignment='dipped')

    assert good and good['success'] and poor and poor['success']
    assert poor['peaks'][0]['height'] == pytest.approx(
        good['peaks'][0]['height'], rel=0.25), (
        'one bad sample decided the answer instead of the data around it')


def test_a_fit_reports_the_boxes_it_was_confined_to():
    """A fit that stopped on a wall is otherwise indistinguishable from one that found
    a minimum. The intensity is reported as a FRACTION of its interval, not as bounds:
    the interval is in the fitter's normalised units while the intensity beside it has
    been scaled back, so a bound would be in different units from the value."""
    truth = (121.50, 7.700, AMPLITUDE)
    result = _integrator([truth]).fit_overlap_group_2d(
        [{'x_ppm': truth[1], 'y_ppm': truth[0]}], assignment='onTarget')

    assert result and result['success']
    peak = result['peaks'][0]
    for key in ('pos_bounds_f1', 'pos_bounds_f2'):
        lower, upper = peak[key]
        assert lower < peak[key.replace('_bounds', '')] < upper

    fraction = peak['intensity_fraction']
    assert 0.0 < fraction < 1.0, 'a peak on its own maximum is not on either bound'


def test_a_weak_cluster_member_pins_on_the_floor_its_neighbours_set():
    """The intensity floor is CLUSTER-RELATIVE - 0.1% of the brightest member's seed -
    so a weak member among bright ones is driven onto it and stays. Its volume is then
    a product of its neighbours' seeds, and the cluster R-squared its neighbours carry
    cannot show it."""
    truth = [(121.50, 7.700, AMPLITUDE)]
    integrator = _integrator(truth)
    group = [{'x_ppm': 7.700, 'y_ppm': 121.50},
             {'x_ppm': 7.720, 'y_ppm': 121.50}]      # nothing here at all
    result = integrator.fit_overlap_group_2d(
        group, assignment='pair', peak_assignments=['bright', 'empty'])

    assert result and result['success']
    fractions = [pk['intensity_fraction'] for pk in result['peaks']]
    assert all(f is not None for f in fractions)
    assert min(fractions) < max(fractions), (
        'the empty member must be distinguishable from the real one, which the '
        'cluster R-squared cannot do')


def test_a_pinned_fit_is_reported_as_pinned():
    """`peak_reached_a_bound` turns the boxes the fit reports into the two questions a
    caller actually has. Reported, never gated on - what to do with a peak held at a
    bound is the caller's decision, and a fading peak is still a measurement."""
    from lunaNMR.core.core_integrator import peak_reached_a_bound

    settled = {'pos_f1': 121.5, 'pos_bounds_f1': (121.4, 121.6),
               'pos_f2': 7.70, 'pos_bounds_f2': (7.69, 7.71),
               'intensity_fraction': 0.5}
    assert peak_reached_a_bound(settled) == (False, False)

    onWall = dict(settled, pos_f1=121.4)
    assert peak_reached_a_bound(onWall)[0] is True

    onFloor = dict(settled, intensity_fraction=0.0)
    assert peak_reached_a_bound(onFloor)[1] is True

    # Silence is not "no": a fit that did not report its box says nothing either way.
    assert peak_reached_a_bound({'intensity_fraction': None}) == (None, None)
