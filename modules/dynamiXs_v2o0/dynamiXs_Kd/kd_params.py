# ABOUTME: Importable binding-parameters JSON for Kd titration (concentrations, scales, P0, ...).
# ABOUTME: One schema shared with fit-JSON metadata, so a saved fit JSON is also importable.

import glob
import json
import os

PARAMS_SUFFIX = "_kd_params.json"       # written per analysis, next to results
SIBLING_NAME = "kd_params.json"         # hand-placed next to an input CSV

# Canonical keys (identical to the fit JSON's metadata) with their defaults.
_DEFAULTS = {
    'points': None,
    'concentrations': None,
    'intensity_scales': None,
    'protein_conc': 50.0,
    'alpha': 0.14,
    'observables': ['csp', 'intensity'],
    'intensity_value': 'height',
    'n_bootstrap': 0,
}

_OBS_BY_INDEX = {0: ['csp', 'intensity'], 1: ['csp'], 2: ['intensity']}


def _parse_observables(value):
    """Accept a canonical list or a hand-written string ('csp+intensity', 'csp',
    'intensity', 'both') and return the canonical ordered list."""
    if isinstance(value, str):
        tokens = value.replace('+', ' ').replace(',', ' ').split()
        if 'both' in tokens:
            tokens = ['csp', 'intensity']
    else:
        tokens = list(value or [])
    tokens = {str(t).strip().lower() for t in tokens}
    obs = [o for o in ('csp', 'intensity') if o in tokens]
    return obs or ['csp', 'intensity']


def _float_or_none_list(value):
    if value is None:
        return None
    return [float(v) for v in value]


def normalize_params(d):
    """Coerce a raw params dict (from JSON or the GUI) to the canonical schema,
    filling defaults and accepting the 'observable'/'bootstrap' hand-edit aliases."""
    d = dict(d or {})
    obs_source = d['observables'] if 'observables' in d else d.get('observable')
    boot = d.get('n_bootstrap', d.get('bootstrap', _DEFAULTS['n_bootstrap']))
    return {
        'points': _float_or_none_list(d.get('points')),
        'concentrations': _float_or_none_list(d.get('concentrations')),
        'intensity_scales': _float_or_none_list(d.get('intensity_scales')),
        'protein_conc': float(d.get('protein_conc', _DEFAULTS['protein_conc'])),
        'alpha': float(d.get('alpha', _DEFAULTS['alpha'])),
        'observables': _parse_observables(obs_source),
        'intensity_value': str(d.get('intensity_value', _DEFAULTS['intensity_value'])),
        'n_bootstrap': int(boot),
    }


def load_params(path):
    """Load params from a standalone *_kd_params.json OR a saved *_kd_fit_data.json
    (params live under its 'metadata'). Returns the normalized params dict."""
    with open(path) as fh:
        raw = json.load(fh)
    if isinstance(raw, dict) and isinstance(raw.get('metadata'), dict):
        raw = raw['metadata']           # a fit JSON — params are in metadata
    return normalize_params(raw)


def dump_params_json(path, params):
    """Write normalized params to a JSON file (indent=2, hand-editable)."""
    with open(path, 'w') as fh:
        json.dump(normalize_params(params), fh, indent=2)


def observables_to_combo_index(observables):
    """Map an observables list to the GUI 'Observable' combo index (0/1/2)."""
    obs = set(_parse_observables(observables))
    for idx, ref in _OBS_BY_INDEX.items():
        if obs == set(ref):
            return idx
    return 0


def combo_index_to_observables(index):
    """Map the GUI 'Observable' combo index to an observables list."""
    return list(_OBS_BY_INDEX.get(int(index), _OBS_BY_INDEX[0]))


def _newest(pattern):
    matches = glob.glob(pattern)
    return max(matches, key=os.path.getmtime) if matches else None


def find_params_source(csv_path):
    """Locate a params source near an input titration CSV, or None. Search order:
    a sibling kd_params.json, then the newest *_kd_params.json beside the CSV, then
    inside a sibling 'kd_analysis' output folder (newest *_kd_params.json, else the
    newest *_kd_fit_data.json)."""
    d = os.path.dirname(os.path.abspath(csv_path))
    candidates = [
        os.path.join(d, SIBLING_NAME),
        _newest(os.path.join(d, "*" + PARAMS_SUFFIX)),
        _newest(os.path.join(d, "kd_analysis", "*" + PARAMS_SUFFIX)),
        _newest(os.path.join(d, "kd_analysis", "*_kd_fit_data.json")),
    ]
    for c in candidates:
        if c and os.path.isfile(c):
            return c
    return None
