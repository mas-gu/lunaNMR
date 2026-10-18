# ABOUTME: The series tracking table matches each reference peak to its own result.
# ABOUTME: Numeric assignments that differ after the decimal point stay distinct peaks.

"""
The table compared str(int(float(a))) on both sides, so 45.1 and 45.2 both read '45' and
the second peak got the first one's height, volume and position - in the tracking table
and in every matrix built from it.
"""

import pandas as pd


def _tracking(reference, results):
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    processor = MultiSpectrumProcessor({})
    processor.reference_peaks = pd.DataFrame(
        {'Assignment': reference, 'Position_X': 8.0, 'Position_Y': 120.0})
    batch = {'results': {'s1.ft': {'success': True, 'integration_results': [
        {'Assignment': a, 'Height': h, 'Volume': h / 10} for a, h in results]}}}
    return processor._create_peak_tracking_table(batch)


def test_assignments_differing_after_the_decimal_point_are_distinct():
    table = _tracking(['45.1', '45.2', '12NH'], [('45.1', 100.0), ('45.2', 999.0), ('12NH', 5.0)])
    assert table['s1_Height'].tolist() == [100.0, 999.0, 5.0]


def test_an_integer_assignment_matches_its_float_spelling():
    table = _tracking([142], [('142.0', 7.0)])
    assert table['s1_Height'].tolist() == [7.0]
