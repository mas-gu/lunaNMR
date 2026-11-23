# ABOUTME: The order of residues in lunaNMR: by residue number, a name without a number last.
# ABOUTME: The DynamiXs and Kd modules import nothing from lunaNMR and carry the same rule.

import re


def residue_sort_key(name):
    """Order residues by sequence number, not amino-acid letter: 'K14' before
    'A17'. Names without a number sort last, alphabetically among themselves; a dummy
    peak is one of them, since its index is not a residue number."""
    text = str(name)
    m = None if text.lower().startswith('dummy') else re.search(r'\d+', text)
    return (int(m.group()) if m else float('inf'), text)
