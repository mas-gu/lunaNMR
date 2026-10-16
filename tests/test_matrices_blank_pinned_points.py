# ABOUTME: A point whose fit stopped on its intensity floor is written empty in the intensity and volume matrices.
# ABOUTME: A position-pinned or decayed point keeps its value - neither makes the intensity wrong.

"""
The matrices carry no fit flags, and the relaxation analyses read only the matrices. A
height taken from the intensity floor - cluster-relative, so a product of the
NEIGHBOURS' seeds - would otherwise be fitted as a decay point. Blanking it is
point-wise: the residue keeps every other point.

A position that stopped on its margin is not blanked here. Measured on KRASB 600 hetNOE,
most such peaks stop on the 0.020 ppm 15N margin, a quarter of one 15N point; widened
to 0.05 ppm the freed peaks' heights moved by a median 0.0% (q10-q90 -3% to +5%).
Blanking them halved the hetNOE residue set (82 -> 41).
"""

import pandas as pd

from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor


def _peak(assignment, height, position_pinned=False, intensity_pinned=False):
    return {'Assignment': assignment, 'Height': height, 'Volume': 2 * height,
            'detected_intensity': height, 'SNR': 10.0, 'Quality': 'Good',
            'Position_X': 8.0, 'Position_Y': 120.0, 'R_Squared': 0.98,
            'position_pinned': position_pinned, 'intensity_pinned': intensity_pinned}


def _write_matrices(tmp_path):
    p = MultiSpectrumProcessor.__new__(MultiSpectrumProcessor)
    p.output_folder = str(tmp_path)
    p.reference_peaks = pd.DataFrame({
        'Assignment': ['Healthy', 'Decayed', 'PosPinned', 'IntPinned', 'Legacy'],
        'Position_X': [8.0] * 5, 'Position_Y': [120.0] * 5})
    batch = {'results': {
        'T2_10ms.ft': {'success': True, 'integration_results': [
            _peak('Healthy', 1.0e9), _peak('Decayed', 1.0e9),
            _peak('PosPinned', 1.0e9), _peak('IntPinned', 1.0e9),
            _peak('Legacy', 1.0e9, position_pinned=None, intensity_pinned=None)]},
        'T2_90ms.ft': {'success': True, 'integration_results': [
            _peak('Healthy', 4.0e8), _peak('Decayed', 2.0e6),
            _peak('PosPinned', 3.0e8, position_pinned=True),
            _peak('IntPinned', 1.0e6, intensity_pinned=True),
            _peak('Legacy', 5.0e8, position_pinned=None, intensity_pinned=None)]},
    }}
    tracking = p._create_peak_tracking_table(batch)
    p._create_intensity_matrix(tracking)
    p._create_volume_matrix(tracking)
    p._create_detected_intensity_matrix(tracking)
    read = lambda name: pd.read_csv(tmp_path / name).set_index('Assignment')
    return (read('peak_intensity_matrix.csv'), read('peak_volume_matrix.csv'),
            read('peak_intensity_detected_matrix.csv'))


def test_an_intensity_pinned_point_is_empty_in_the_intensity_and_volume_matrices(tmp_path):
    intensity, volume, _ = _write_matrices(tmp_path)
    for matrix in (intensity, volume):
        assert pd.isna(matrix.loc['IntPinned', 'T2_90ms'])


def test_a_position_pinned_point_keeps_its_intensity(tmp_path):
    intensity, volume, _ = _write_matrices(tmp_path)
    assert intensity.loc['PosPinned', 'T2_90ms'] == 3.0e8
    assert volume.loc['PosPinned', 'T2_90ms'] == 6.0e8


def test_blanking_is_point_wise(tmp_path):
    intensity, volume, _ = _write_matrices(tmp_path)
    for matrix in (intensity, volume):
        assert matrix.loc['PosPinned', 'T2_10ms'] > 0
        assert matrix.loc['IntPinned', 'T2_10ms'] > 0


def test_a_decayed_point_keeps_its_value(tmp_path):
    """Dropping a low point from a decay invents a plateau."""
    intensity, _, _ = _write_matrices(tmp_path)
    assert intensity.loc['Decayed', 'T2_90ms'] == 2.0e6


def test_an_unknown_flag_is_not_read_as_pinned(tmp_path):
    """A fit that did not report its bounds is unknown, and unknown is kept."""
    intensity, _, _ = _write_matrices(tmp_path)
    assert intensity.loc['Legacy', 'T2_90ms'] == 5.0e8


def test_the_detected_intensity_matrix_is_left_alone(tmp_path):
    """It is read from the spectrum at the peak, not from the fit, so a fit's bound says
    nothing about it."""
    _, _, detected = _write_matrices(tmp_path)
    assert detected.loc['IntPinned', 'T2_90ms'] == 1.0e6


def test_every_output_the_series_summary_names_is_written(tmp_path):
    """`series --format json` lists its outputs by name and keeps only those that exist,
    so a misspelt name is not an error - the file just never appears in the summary."""
    from lunaNMR.cli import SERIES_OUTPUTS
    _write_matrices(tmp_path)
    p = MultiSpectrumProcessor.__new__(MultiSpectrumProcessor)
    p.output_folder = str(tmp_path)
    p.reference_peaks = pd.DataFrame({'Assignment': ['A'], 'Position_X': [8.0],
                                      'Position_Y': [120.0]})
    p._create_detection_matrix(p._create_peak_tracking_table(
        {'results': {'T2_10ms.ft': {'success': True,
                                    'integration_results': [_peak('A', 1.0e9)]}}}))
    matrices = [n for n in SERIES_OUTPUTS if n.endswith('_matrix.csv')]
    assert matrices
    assert [n for n in matrices if not (tmp_path / n).exists()] == []


def _tracking_csv(tmp_path, pinned_peak):
    p = MultiSpectrumProcessor.__new__(MultiSpectrumProcessor)
    p.output_folder = str(tmp_path)
    p.reference_peaks = pd.DataFrame({'Assignment': ['Healthy', 'Pinned'],
                                      'Position_X': [8.0, 8.5], 'Position_Y': [120.0, 121.0]})
    batch = {'results': {
        'tit_0o0.ft': {'success': True, 'integration_results': [
            _peak('Healthy', 1.0e9), _peak('Pinned', 1.0e9)]},
        'tit_1o0.ft': {'success': True, 'integration_results': [
            _peak('Healthy', 6.0e8), pinned_peak]},
    }}
    path = tmp_path / 'comprehensive_peak_tracking.csv'
    p._format_dataframe_for_csv(p._create_peak_tracking_table(batch)).to_csv(path, index=False)
    return path


def _load_tracking(path):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]
                           / 'modules' / 'dynamiXs_v2o0' / 'dynamiXs_Kd'))
    from kd_input import load_titration_tracking
    return load_titration_tracking(str(path))


def test_kd_reading_the_tracking_csv_keeps_the_intensity_of_a_position_pinned_point(tmp_path):
    """A position on its margin is not a shift anyone measured, so CSP must not see it;
    the intensity does not depend on it, so I/I0 keeps it."""
    path = _tracking_csv(tmp_path, _peak('Pinned', 7.0e8, position_pinned=True))
    _, residues = _load_tracking(path)
    pinned = residues['Pinned']
    assert pinned['height'][1] == 7.0e8 and pinned['volume'][1] == 1.4e9
    assert pd.isna(pinned['ppm_x'][1]) and pd.isna(pinned['ppm_y'][1])
    assert pd.notna(pinned['ppm_x'][0]), 'blanking is point-wise'


def test_kd_reading_the_tracking_csv_blanks_an_intensity_pinned_point(tmp_path):
    """The tracking CSV keeps every value and marks the bounds in `{pt}_Position_Pinned`
    and `{pt}_Intensity_Pinned`. Written by the real formatter, read by the real Kd
    loader: the flag survives the CSV as text ('1.0'), which is where a loader can
    silently miss it."""
    path = _tracking_csv(tmp_path, _peak('Pinned', 1.0e6, intensity_pinned=True))
    points, residues = _load_tracking(path)
    pinned, healthy = residues['Pinned'], residues['Healthy']
    for key in ('height', 'volume', 'ppm_x', 'ppm_y'):
        assert pd.isna(pinned[key][1]), f'{key} taken from a bound was read as data'
        assert pd.notna(pinned[key][0]), 'blanking is point-wise'
        assert pd.notna(healthy[key][1])
