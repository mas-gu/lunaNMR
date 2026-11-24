# ABOUTME: The T1/T2 fit PDFs draw every residue in order, by either fitter: a failed one with its
# ABOUTME: data and the words that it was not fitted. Real 700 MHz T1 matrix, one residue left with two delays.

import importlib
import re
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_DIR))
REPO_ROOT = Path(__file__).resolve().parents[3]

REAL_T1 = (REPO_ROOT.parent / "data_test" / "data_JA1_relaxation" / "run_20260825_1405"
           / "series" / "700_T1" / "peak_intensity_matrix.csv")

pytestmark = pytest.mark.skipif(not REAL_T1.exists(), reason="data_test not present")


def _pdf_texts(path):
    """The text strings drawn in a matplotlib PDF, page order kept."""
    raw = Path(path).read_bytes()
    texts = []
    for stream in re.finditer(rb'stream\r?\n(.*?)\r?\nendstream', raw, re.S):
        try:
            content = zlib.decompress(stream.group(1))
        except zlib.error:
            continue
        for array, single in re.findall(rb'\[(.*?)\]\s*TJ|\((.*?)\)\s*Tj', content, re.S):
            text = b''.join(re.findall(rb'\(((?:[^()\\]|\\.)*)\)', array)) if array else single
            texts.append(text.decode('latin-1'))
    return texts


@pytest.fixture(scope="module", params=["fit_Tx_NMRRE", "fitmulti__Tx_NMRRE"])
def run(request, tmp_path_factory):
    """The real T1 matrix with its first residue measured at its two shortest delays only:
    two points cannot carry a fit of two parameters with one to spare, so it fails before
    fitting. Two residues of this matrix fit to a time constant beyond the window and fail
    after it, their data still attached to the result."""
    import matplotlib
    matplotlib.use("Agg")
    from delay_parser import parse_delay_columns
    fitter = importlib.import_module(request.param)
    folder = tmp_path_factory.mktemp("failed_residues")
    matrix = pd.read_csv(REAL_T1)
    delays = list(matrix.columns[4:])
    x, keep, _ = parse_delay_columns(delays)
    kept = [delays[i] for i in np.flatnonzero(keep)]
    shortest_two = {kept[i] for i in np.argsort(x)[:2]}
    for column in delays:
        if column not in shortest_two:
            matrix.loc[0, column] = np.nan
    truncated = matrix.loc[0, 'Assignment']
    matrix.to_csv(folder / 'peak_intensity_matrix.csv', index=False)
    summary = fitter.run_analysis_with_params({
        'input_csv_file': str(folder / 'peak_intensity_matrix.csv'),
        'output_prefix': str(folder / 'T1'),
        'results_txt_file': str(folder / 'T1_results.txt'),
        'experiment_type': 'T1', 'time_units': 'ms', 'error_method': 'analytical'})
    return folder, summary, truncated, list(matrix['Assignment'])


def test_a_residue_fails_before_fitting_and_two_after(run):
    _, summary, truncated, _ = run
    assert summary['failed_residues'] == [truncated, '10ValH', '43GluH']
    assert summary['n_fitted'] == summary['n_total'] - len(summary['failed_residues'])


def test_the_pdfs_draw_every_residue_in_order(run):
    folder, _, _, residues = run
    pdfs = sorted(folder.glob('T1_fit_results_fig*.pdf'), key=lambda p: int(p.stem.rsplit('fig', 1)[1]))
    assert len(pdfs) == int(np.ceil(len(residues) / 20))
    texts = [text for pdf in pdfs for text in _pdf_texts(pdf)]
    assert [text[len('Residue: '):] for text in texts if text.startswith('Residue: ')] == residues


def test_a_failed_residue_is_drawn_as_not_fitted(run):
    folder, summary, _, residues = run
    texts = [text for pdf in sorted(folder.glob('T1_fit_results_fig*.pdf'))
             for text in _pdf_texts(pdf)]
    assert texts.count('T1: not fitted') == len(summary['failed_residues'])
    assert texts.count('Fit') == len(residues) - len(summary['failed_residues'])
    assert not [text for text in texts if 'nan' in text.lower()]
