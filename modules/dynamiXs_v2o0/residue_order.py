# ABOUTME: The order of residues in DynamiXs and Kd: by residue number, names without one last.
# ABOUTME: The sort key for lists and tables, and the x positions of a plot against residue number.

import re


def residue_sort_key(name):
    """Order residues by sequence number, not amino-acid letter: 'K14' before
    'A17'. Names without a number sort last, alphabetically among themselves; a dummy
    peak is one of them, since its index is not a residue number."""
    text = str(name)
    m = None if text.lower().startswith('dummy') else re.search(r'\d+', text)
    return (int(m.group()) if m else float('inf'), text)


def residue_axis_positions(names):
    """The x position of each name on a residue-number axis: its own number, and for
    names without one the positions after the highest number, in sorted order."""
    keys = [residue_sort_key(name) for name in names]
    highest = max((number for number, _ in keys if number != float('inf')), default=0)
    unnumbered = sorted({text for number, text in keys if number == float('inf')})
    return [highest + 1 + unnumbered.index(text) if number == float('inf') else number
            for number, text in keys]
