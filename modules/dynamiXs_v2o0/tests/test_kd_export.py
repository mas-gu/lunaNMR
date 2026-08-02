# ABOUTME: Tests the wide CSV + JSON writers for the Kd vs-sequence (ref->point) data.
# ABOUTME: Excel-ready wide CSV (blank=absent) + a structured JSON mirror (null=absent).

import sys
import csv
import json
from pathlib import Path

_KD_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_Kd"
sys.path.insert(0, str(_KD_DIR))


def _fits():
    return [
        {'residue': 'A17', 'series': {
            'ppm_x': [8.00, 8.02, 8.05], 'ppm_y': [120.0, 120.1, 120.3],
            'height': [100.0, 80.0, 50.0]}},
        {'residue': 'K14', 'series': {
            'ppm_x': [7.50, 7.51, 0.0], 'ppm_y': [110.0, 110.05, 0.0],
            'height': [200.0, 150.0, 0.0]}},
    ]


def test_export_writes_csv_and_json(tmp_path):
    from kd_export import export_ref_vs_point
    base = tmp_path / "X_intensity_ref_vs_point"
    paths = export_ref_vs_point(str(base), _fits(), [0.0, 1.0, 2.0], 'intensity')
    assert sorted(Path(p).suffix for p in paths) == ['.csv', '.json']

    with open(str(base) + ".csv") as fh:
        rows = list(csv.reader(fh))
    assert rows[0] == ['Residue', '0.0', '1.0', '2.0']
    assert [r[0] for r in rows[1:]] == ['K14', 'A17']       # sequence order
    a17 = next(r for r in rows[1:] if r[0] == 'A17')
    assert a17[1] == '1.0'                                    # I(0)/I(0)
    k14 = next(r for r in rows[1:] if r[0] == 'K14')
    assert k14[3] == ''                                       # absent -> blank


def test_json_mirror_structure(tmp_path):
    from kd_export import export_ref_vs_point
    base = tmp_path / "X_csp_ref_vs_point"
    export_ref_vs_point(str(base), _fits(), [0.0, 1.0, 2.0], 'csp')
    data = json.loads((Path(str(base) + ".json")).read_text())
    assert data['observable'] == 'csp'
    assert data['unit'] == 'ppm'
    assert data['labels'] == [0.0, 1.0, 2.0]
    k14 = next(r for r in data['residues'] if r['residue'] == 'K14')
    assert k14['values'][0] == 0.0        # CSP reference vs itself
    assert k14['values'][2] is None       # absent -> null


def _global_intensity_fits():
    return [
        {'residue': 'A17', 'intensity': {'success': True,
            'L': [0.0, 1.0, 2.0], 'obs': [1.0, 0.6, 0.4]}},
        {'residue': 'K14', 'intensity': {'success': True,
            'L': [0.0, 1.0, 2.0], 'obs': [1.0, 0.5, 0.3]}},
    ]


def _global_intensity():
    return {'intensity': {'success': True, 'Kd': 1.5, 'Kd_err': 0.2,
                          'I0': {'A17': 1.0, 'K14': 1.0},
                          'I_inf': {'A17': 0.3, 'K14': 0.2}}}


def test_export_global_fit_intensity(tmp_path):
    from kd_export import export_global_fit
    base = tmp_path / "X_intensity_global_fit"
    paths = export_global_fit(str(base), _global_intensity_fits(),
                              _global_intensity(), 'intensity', P0=50.0)
    assert sorted(Path(p).suffix for p in paths) == ['.csv', '.json']

    with open(str(base) + ".csv") as fh:
        rows = list(csv.reader(fh))
    assert rows[0] == ['Residue', 'I0', 'I_inf', 'global_Kd', 'R2']
    assert [r[0] for r in rows[1:]] == ['K14', 'A17']      # sequence order
    a17 = next(r for r in rows[1:] if r[0] == 'A17')
    assert float(a17[3]) == 1.5                             # shared Kd column

    data = json.loads((Path(str(base) + ".json")).read_text())
    assert data['observable'] == 'intensity'
    assert data['global_Kd'] == 1.5
    assert data['global_Kd_err'] == 0.2
    k14 = next(r for r in data['residues'] if r['residue'] == 'K14')
    assert k14['I0'] == 1.0 and k14['I_inf'] == 0.2


def test_export_global_fit_csp_header(tmp_path):
    from kd_export import export_global_fit
    base = tmp_path / "X_csp_global_fit"
    fits = [{'residue': 'A17', 'csp': {'success': True,
             'L': [0.0, 1.0, 2.0], 'obs': [0.0, 0.05, 0.08]}}]
    g = {'csp': {'success': True, 'Kd': 2.0, 'dd_max': {'A17': 0.1}}}
    export_global_fit(str(base), fits, g, 'csp', P0=50.0)
    with open(str(base) + ".csv") as fh:
        rows = list(csv.reader(fh))
    assert rows[0] == ['Residue', 'dd_max', 'global_Kd', 'R2']
    assert rows[1][0] == 'A17'


def test_export_global_fit_absent_returns_empty(tmp_path):
    from kd_export import export_global_fit
    base = tmp_path / "X_intensity_global_fit"
    paths = export_global_fit(str(base), _global_intensity_fits(),
                              {}, 'intensity', P0=50.0)
    assert paths == []
    assert not (Path(str(base) + ".csv")).exists()
