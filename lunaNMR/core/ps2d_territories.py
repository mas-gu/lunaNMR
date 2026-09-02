# ABOUTME: Per-peak position bounds narrow enough that no two peaks can reach one position.
# ABOUTME: Keeps a tracked peak from taking a neighbour's position, and so its identity.

"""
Position bounds for cascade tracking.

A titration seeds each spectrum at the previous spectrum's fitted positions and lets
every peak move by one per-step window. Clustered peaks are within the overlap
threshold of each other by definition - 0.04 ppm in 1H - while the window reaches
0.03 ppm, so a clustered peak is licensed to travel more than halfway to its
neighbour in a single step. Peaks in one cluster are also parameters of a single
optimisation, which nothing stops from putting two of them on the same maximum.

The remedy is to give each peak its own bound instead of one window for all: no
further than halfway to its nearest neighbour. Then the reachable regions are
disjoint, and a peak can neither meet another nor cross one - so the arrangement of
the peak list is preserved for the whole series, not merely checked afterwards.

Half is the same fraction `_scoreTracking` already calls untrusted: a peak that moved
more than halfway to a rival could equally have arrived AT the rival. This makes that
region unreachable rather than reportable.
"""

# Just under half, so two peaks driven at each other stop short of meeting. Bounds
# are hard clips: at exactly half, both peaks pinned at their bound would land on the
# same face and be coincident, which is the failure this exists to prevent.
TERRITORY_FRACTION = 0.49


def territory_margins(positions, margin_f1, margin_f2):
    """Per-peak position margins whose reachable regions do not overlap.

    Distance is measured in units of the requested window, which is what makes the
    two dimensions comparable: 0.02 ppm is most of a 1H window and a rounding error
    in 15N. Each peak keeps the SHAPE of the requested window, scaled down until it
    clears its nearest neighbour.

    :param positions: [(f1, f2)] ppm per peak, in peak order
    :param margin_f1: the per-step window in F1 (indirect), ppm
    :param margin_f2: the per-step window in F2 (1H), ppm
    :return: [(margin_f1, margin_f2)] per peak, parallel to `positions`
    """
    if margin_f1 <= 0 or margin_f2 <= 0:
        raise ValueError(f'a window must be positive, not {margin_f1}/{margin_f2}')

    scaled = [(f1 / margin_f1, f2 / margin_f2) for f1, f2 in positions]

    margins = []
    for i, (u, v) in enumerate(scaled):
        # The window itself is the unit box in these coordinates, so a peak never
        # gets more than the caller asked for however far away its neighbours are.
        room = 1.0
        for j, (other_u, other_v) in enumerate(scaled):
            if j == i:
                continue
            gap = max(abs(u - other_u), abs(v - other_v))
            room = min(room, TERRITORY_FRACTION * gap)
        margins.append((room * margin_f1, room * margin_f2))
    return margins
