# ABOUTME: The integrator hands the 2D fitter one position margin per peak of a cluster.
# ABOUTME: Protection is all-or-nothing: a partial set would leave territories overlapping.

import numpy as np

from lunaNMR.core.core_integrator import VoigtIntegrator


def _integrator(territories=None):
    # _territory_margins_for reads one attribute and nothing else, so skip __init__
    # rather than build a whole integrator around it.
    integrator = VoigtIntegrator.__new__(VoigtIntegrator)
    integrator.titration_pos_margins = territories
    return integrator


TERRITORIES = {'3': (0.30, 0.030), '4': (0.15, 0.015), '5': (0.02, 0.002)}


def test_no_territories_leaves_the_shared_margin_in_charge():
    """Relaxation mode sets none and must keep the behaviour it always had."""
    assert _integrator()._territory_margins_for(['3', '4']) is None


def test_margins_come_back_in_peak_order():
    """The fitter indexes them by peak, so the order is the contract."""
    assert _integrator(TERRITORIES)._territory_margins_for(['5', '3']) == \
           [(0.02, 0.002), (0.30, 0.030)]


def test_numeric_assignments_still_match():
    """A peak list stores assignments as whatever the DataFrame held - float, numpy
    int, str - while the territories are keyed from that same list elsewhere. Both
    sides go through str(), so the two representations meet."""
    territories = {'3.0': (0.30, 0.030), '4': (0.15, 0.015)}
    assert _integrator(territories)._territory_margins_for([3.0, np.int64(4)]) == \
           [(0.30, 0.030), (0.15, 0.015)]


def test_one_missing_peak_disables_the_whole_cluster():
    """A peak left on the shared margin can still reach a neighbour that gave up
    room for it, so a partial set is not a weaker guarantee - it is none."""
    assert _integrator(TERRITORIES)._territory_margins_for(['3', 'unseen']) is None
