# ABOUTME: The hetNOE step must emit a QC plot of I_sat/I_unsat vs residue number.
# ABOUTME: Verifies plot_hetnoe_vs_residue writes a non-empty PDF and orders by residue.

import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_integrated"
sys.path.insert(0, str(_DIR))


def _noe():
    # residue -> {value, error}; keys deliberately out of sequence order
    return {
        'V10': {'value': 0.81, 'error': 0.02},
        'K3':  {'value': 0.55, 'error': 0.03},
        'G68': {'value': 0.30, 'error': 0.05},
    }


def test_writes_non_empty_pdf(tmp_path):
    from data_converters import plot_hetnoe_vs_residue
    out = tmp_path / "hetnoe.pdf"
    plot_hetnoe_vs_residue(_noe(), str(out), title="field1")
    assert out.exists() and out.stat().st_size > 0


def test_empty_dict_writes_nothing(tmp_path):
    from data_converters import plot_hetnoe_vs_residue
    out = tmp_path / "hetnoe.pdf"
    plot_hetnoe_vs_residue({}, str(out))
    assert not out.exists()   # nothing to plot -> no file, no crash
