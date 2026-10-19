# ABOUTME: The fitted peaks of a cluster go back to its members by distance in linewidths.
# ABOUTME: Raw ppm counts a 15N ppm as a 1H ppm and gives a member its neighbour's fit.

"""
A 2D cluster fit returns its peaks in no particular identity, and each is handed back
to the member it lies closest to (one fit per member). In raw ppm, a fit 0.1 ppm off in
15N - a quarter of a 15N linewidth - counted as further than one 0.03 ppm off in 1H,
three quarters of a 1H linewidth. Two members a linewidth apart in 1H, each with its
fit a quarter-linewidth off in 15N, then swapped fits.
"""

import pytest


@pytest.fixture
def nucleus():
    from lunaNMR.core.ps2d_config import set_ps2d_config
    yield set_ps2d_config
    set_ps2d_config('15N')


def _fit(x, y):
    return {'pos_f2': x, 'pos_f1': y}


def test_each_member_takes_the_fit_closest_in_linewidths(nucleus):
    from lunaNMR.core.parallel_voigt_processor import match_fits_to_members
    nucleus('15N')
    members = [(8.00, 120.0), (8.03, 120.1)]
    fits = [_fit(8.03, 120.0), _fit(8.00, 120.1)]

    assert match_fits_to_members(members, fits) == {0: 1, 1: 0}


def test_a_13c_cluster_is_measured_in_13c_linewidths(nucleus):
    from lunaNMR.core.parallel_voigt_processor import match_fits_to_members
    from lunaNMR.core.ps2d_config import distance_in_linewidths
    config = nucleus('13C')
    members = [(1.00, 20.0), (1.02, 20.0)]
    fits = [_fit(1.00, 20.0 + 0.4 * config.radF1), _fit(1.02, 20.0 + 0.4 * config.radF1)]

    assert match_fits_to_members(members, fits) == {0: 0, 1: 1}
    assert distance_in_linewidths(0.02, 0.4 * config.radF1) == pytest.approx(
        ((0.02 / config.radF2) ** 2 + 0.4 ** 2) ** 0.5)


def test_no_fits_maps_nothing():
    from lunaNMR.core.parallel_voigt_processor import match_fits_to_members
    assert match_fits_to_members([(8.0, 120.0)], []) == {}
