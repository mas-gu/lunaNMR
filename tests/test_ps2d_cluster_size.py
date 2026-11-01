# ABOUTME: A PS2D overlap cluster holds at most 20 peaks, whichever nucleus is selected.
# ABOUTME: Larger groups are split, so the limit decides which peaks are fitted together.

import pytest

from lunaNMR.core.core_integrator import VoigtIntegrator
from lunaNMR.core.ps2d_config import set_ps2d_config

MAX_CLUSTER_SIZE = 20


@pytest.mark.parametrize("nucleus", ["15N", "13C"])
def test_default_max_cluster_size(nucleus):
    assert set_ps2d_config(nucleus).max_cluster_size == MAX_CLUSTER_SIZE


@pytest.fixture(autouse=True)
def fresh_config():
    set_ps2d_config('15N')
    yield
    set_ps2d_config('15N')


def _chain(count, spacing=0.02):
    """`count` peaks in a row along 1H, each touching the next; 0.58 ppm across at 30."""
    return [{'x_ppm': 7.5 + index * spacing, 'y_ppm': 121.0,
             'pos_x': 7.5 + index * spacing, 'pos_y': 121.0} for index in range(count)]


def test_a_chain_is_split_at_the_cap():
    clusters = VoigtIntegrator().identify_overlap_clusters(_chain(30))
    assert max(len(cluster) for cluster in clusters) <= MAX_CLUSTER_SIZE
