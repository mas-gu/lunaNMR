# ABOUTME: A PS2D overlap cluster holds at most 20 peaks, whichever nucleus is selected.
# ABOUTME: Larger groups are split, so the limit decides which peaks are fitted together.

import pytest

from lunaNMR.core.ps2d_config import set_ps2d_config

MAX_CLUSTER_SIZE = 20


@pytest.mark.parametrize("nucleus", ["15N", "13C"])
def test_default_max_cluster_size(nucleus):
    assert set_ps2d_config(nucleus).max_cluster_size == MAX_CLUSTER_SIZE
