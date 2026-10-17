# ABOUTME: Draws a spectrum's positive AND negative peaks as contours, negatives dashed.
# ABOUTME: Shared by every spectrum view so a negative peak is never hidden or shown as positive.

"""
Negative peaks are real data: a flexible residue in a saturated hetNOE plane, side-chain
NH2 groups there, CH2 in an edited 13C-HSQC, a folded peak. Contouring |data| draws them
as positive peaks; contouring the signed data at positive levels only does not draw them
at all. Both happened here, in different views.

The negative levels mirror the positive ones, so a peak of either sign at the same height
gets the same number of contours. They are dashed - the convention, and independent of
colour, so it survives a per-spectrum colour in an overlay.
"""

import numpy as np


def draw_signed_contours(ax, x, y, data, levels, negative_colors=None, **style):
    """Contour `data` at `levels` and at their mirror images below zero.

    :param levels: positive, increasing contour levels
    :param negative_colors: colour(s) for the negative contours; by default the same
        `colors` as the positive ones. Needed when the positives use a `cmap`.
    :param style: passed to `ax.contour` for both signs (colors/cmap, linewidths, alpha)
    :return: (positive ContourSet, negative ContourSet or None when no point of the
        data reaches the lowest negative level)
    """
    levels = np.asarray(levels, dtype=float)
    positive = ax.contour(x, y, data, levels=levels, **style)

    if not np.nanmin(data) <= -levels[0]:
        return positive, None

    negative_style = dict(style)
    if negative_colors is not None:
        negative_style.pop('cmap', None)
        negative_style['colors'] = negative_colors
    negative = ax.contour(x, y, data, levels=-levels[::-1], linestyles='dashed',
                          **negative_style)
    return positive, negative
