# ABOUTME: Tests the importable binding-parameters JSON (schema, load/dump, discovery).
# ABOUTME: One schema shared with fit-JSON metadata, so old fit JSONs are importable too.

import sys
import json
from pathlib import Path

_KD_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_Kd"
sys.path.insert(0, str(_KD_DIR))


def test_normalize_fills_defaults():
    from kd_params import normalize_params
    p = normalize_params({})
    assert p['protein_conc'] == 50.0
    assert p['alpha'] == 0.14
    assert p['observables'] == ['csp', 'intensity']
    assert p['intensity_value'] == 'height'
    assert p['n_bootstrap'] == 0


def test_normalize_observable_string_alias():
    from kd_params import normalize_params
    assert normalize_params({'observable': 'csp+intensity'})['observables'] == ['csp', 'intensity']
    assert normalize_params({'observable': 'intensity'})['observables'] == ['intensity']
    assert normalize_params({'observable': 'csp'})['observables'] == ['csp']


def test_bootstrap_alias_and_coercion():
    from kd_params import normalize_params
    p = normalize_params({'bootstrap': '500', 'protein_conc': '25', 'alpha': '0.2'})
    assert p['n_bootstrap'] == 500 and p['protein_conc'] == 25.0 and p['alpha'] == 0.2


def test_roundtrip_dump_load(tmp_path):
    from kd_params import dump_params_json, load_params, normalize_params
    params = normalize_params({
        'points': [0.0, 1.0, 2.0], 'concentrations': [0.0, 6.25, 12.5],
        'intensity_scales': [1.2, 1.0, 0.5], 'protein_conc': 50.0, 'alpha': 0.14,
        'observables': ['csp'], 'intensity_value': 'volume', 'n_bootstrap': 200})
    path = tmp_path / "kd_params.json"
    dump_params_json(path, params)
    assert load_params(path) == params


def test_load_from_fit_json_metadata(tmp_path):
    """A saved fit JSON (params under 'metadata') must load as params — this is how
    old datasets become importable without re-fitting."""
    from kd_params import load_params
    fit = {'metadata': {'analysis': 'Kd_titration', 'protein_conc': 30.0, 'alpha': 0.1,
                        'concentrations': [0.0, 5.0], 'points': [0.0, 1.0],
                        'intensity_scales': [1.0, 0.5], 'intensity_value': 'height',
                        'observables': ['csp', 'intensity'], 'n_bootstrap': 100},
           'fits': [], 'global': {}}
    path = tmp_path / "X_kd_fit_data.json"
    path.write_text(json.dumps(fit))
    p = load_params(path)
    assert p['protein_conc'] == 30.0 and p['concentrations'] == [0.0, 5.0]
    assert p['n_bootstrap'] == 100


def test_combo_index_roundtrip():
    from kd_params import observables_to_combo_index, combo_index_to_observables
    for idx, obs in [(0, ['csp', 'intensity']), (1, ['csp']), (2, ['intensity'])]:
        assert observables_to_combo_index(obs) == idx
        assert combo_index_to_observables(idx) == obs


def test_find_params_source_prefers_sibling(tmp_path):
    from kd_params import find_params_source
    csv = tmp_path / "series_analysis_tidy.csv"
    csv.write_text("x\n")
    sibling = tmp_path / "kd_params.json"
    sibling.write_text("{}")
    assert find_params_source(str(csv)) == str(sibling)


def test_find_params_source_none_when_absent(tmp_path):
    from kd_params import find_params_source
    csv = tmp_path / "series_analysis_tidy.csv"
    csv.write_text("x\n")
    assert find_params_source(str(csv)) is None
