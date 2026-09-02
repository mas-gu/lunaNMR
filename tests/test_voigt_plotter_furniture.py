# ABOUTME: The 3D Voigt view's furniture and its per-assignment colour identity.
# ABOUTME: Ported from the CcpNmr PS2D plugin, where each decision was measured.

"""
What the 3D fit view draws around the data, and what colour a residue gets.

These are the two things a render showed were wrong and code review did not: the
default Axes3D cage is a heavier structure than the data inside it, and the same
residue was given different colours in different panels of the same figure.
"""

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pytest

from lunaNMR.gui.components.voigt_analysis_plotter import (
    COMPONENT_COLOURS, furniture_colours, style_3d_furniture, indirect_axis_label)


SCHEMES = {
    'onWhite': {'background': 'white', 'experimental': {'color': 'black'}},
    'onBlack': {'background': 'black', 'experimental': {'color': 'silver'}},
}


@pytest.fixture
def axes3d():
    figure = plt.figure()
    yield figure.add_subplot(projection='3d')
    plt.close(figure)


@pytest.mark.parametrize('name', sorted(SCHEMES))
def test_the_cage_is_stripped_back(axes3d, name):
    """Three filled panes, a grid on all three and near-black axis lines together
    outweigh the surface they surround. The grid goes ENTIRELY, not softened: each
    axis draws its grid on the floor AND the two back walls, and with the panes
    unfilled the wall lines render as grey strokes hanging in the air."""
    style_3d_furniture(axes3d, SCHEMES[name])

    # `Axes3D.grid(False)` records itself on the AXES, not the axis: the per-axis
    # `gridOn` tick kwarg is already False by default on a 3D projection and so says
    # nothing about whether the call was made.
    assert axes3d._draw_grid is False
    for axis in (axes3d.xaxis, axes3d.yaxis, axes3d.zaxis):
        assert axis.pane.fill is False
        # The Line2D, not _axinfo['axisline']: Axis3D builds its line from that dict
        # in __init__ and never reads it again, so setting the dict changes nothing
        # on screen and nothing says so.
        assert axis.line.get_color() != 'k'
        assert axis.line.get_linewidth() == pytest.approx(0.8)


def test_the_neutrals_follow_the_scheme_they_sit_on(axes3d):
    """A neutral has to be legible against what is behind it, and a preset on black
    needs the opposite treatment from one on white."""
    light = furniture_colours(SCHEMES['onWhite'])
    dark = furniture_colours(SCHEMES['onBlack'])
    assert light['muted'] != dark['muted']
    for colours in (light, dark):
        # Ordered: the data mesh is the substrate, the furniture is quieter still.
        assert colours['data'] != colours['muted'] != colours['faint']


def test_one_palette_so_a_residue_keeps_its_colour():
    """`_get_peak_color` indexes with `idx % len(colors)`, so palettes of different
    lengths gave the same residue a different colour in the cross-sections from the
    one it had in the 3D view - four palettes of lengths 5, 5, 8 and 8."""
    import inspect
    from lunaNMR.gui.components import voigt_analysis_plotter as module

    source = inspect.getsource(module)
    assert "colors = ['red', 'orange'" not in source
    assert "colors = ['orange', 'purple'" not in source
    assert source.count('colors = COMPONENT_COLOURS') == 4
    assert len(COMPONENT_COLOURS) == len(set(COMPONENT_COLOURS))


def test_the_indirect_axis_names_the_nucleus_it_was_given():
    """The label was the literal '15N/13C Chemical Shift (ppm)' - a guess printed on
    the axis whose whole job is to say what is being measured."""
    from lunaNMR.core.ps2d_config import set_ps2d_config

    set_ps2d_config('15N')
    nitrogen = indirect_axis_label()
    set_ps2d_config('13C')
    carbon = indirect_axis_label()
    set_ps2d_config('15N')

    assert '15' in nitrogen and '13' not in nitrogen
    assert '13' in carbon and '15' not in carbon
    assert nitrogen != carbon
