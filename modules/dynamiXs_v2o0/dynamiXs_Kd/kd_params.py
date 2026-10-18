# ABOUTME: Importable binding-parameters JSON for Kd titration (concentrations, scales, P0, ...).
# ABOUTME: One schema shared with fit-JSON metadata, so a saved fit JSON is also importable.

import glob
import json
import os

PARAMS_SUFFIX = "_kd_params.json"       # written per analysis, next to results
SIBLING_NAME = "kd_params.json"         # hand-placed next to an input CSV

# Canonical keys (identical to the fit JSON's metadata) with their defaults.
_DEFAULTS = {
    # The input CSV these settings were written for (its file name). A stored file is
    # only ever applied to that input: a folder can hold several titrations.
    'input': None,
    'points': None,
    'concentrations': None,
    'intensity_scales': None,
    'protein_conc': 50.0,
    'alpha': 0.14,
    'observables': ['csp', 'intensity'],
    'intensity_value': 'height',
    'n_bootstrap': 0,
    # 'equivalents' multiplies every concentration by protein_conc. Titration filenames
    # commonly encode equivalents, and read as absolute they give a confident Kd wrong by
    # the factor P0.
    'conc_units': 'absolute',
    # Survey thresholds. Persisted so a dataset remembers how it was judged; CLI flags
    # override. Defaults are measured, not assumed — see kd_survey.py for the evidence.
    # Multiples of the trimmed CSP spread a residue must exceed at the last
    # titration point to enter the shared CSP fit. 2.0 is a common stricter variant.
    'csp_sigma_multiple': 1.0,
    # Robust z (median/MAD on log10 Kd) beyond which a residue leaves the shared
    # CSP fit. 0 or non-finite disables the gate.
    'kd_outlier_z': 3.0,
    'noise_quantile': 0.25,
    'dd_runaway_ratio': 10.0,
    'ref_max_ratio': 10.0,
}

# The survey thresholds above: how a dataset was judged, restored with it by every front end.
THRESHOLD_KEYS = ('csp_sigma_multiple', 'kd_outlier_z', 'noise_quantile',
                  'dd_runaway_ratio', 'ref_max_ratio')

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
    source = d.get('input') or (os.path.basename(str(d['input_csv_file']))
                                if d.get('input_csv_file') else None)
    return {
        'input': source,
        'points': _float_or_none_list(d.get('points')),
        'concentrations': _float_or_none_list(d.get('concentrations')),
        'intensity_scales': _float_or_none_list(d.get('intensity_scales')),
        'protein_conc': float(d.get('protein_conc', _DEFAULTS['protein_conc'])),
        'alpha': float(d.get('alpha', _DEFAULTS['alpha'])),
        'observables': _parse_observables(obs_source),
        'intensity_value': str(d.get('intensity_value', _DEFAULTS['intensity_value'])),
        'n_bootstrap': int(boot),
        'conc_units': str(d.get('conc_units') or _DEFAULTS['conc_units']),
        'csp_sigma_multiple': float(d.get('csp_sigma_multiple',
                                          _DEFAULTS['csp_sigma_multiple'])),
        'kd_outlier_z': float(d.get('kd_outlier_z', _DEFAULTS['kd_outlier_z'])),
        'noise_quantile': float(d.get('noise_quantile', _DEFAULTS['noise_quantile'])),
        'dd_runaway_ratio': float(d.get('dd_runaway_ratio',
                                        _DEFAULTS['dd_runaway_ratio'])),
        'ref_max_ratio': float(d.get('ref_max_ratio', _DEFAULTS['ref_max_ratio'])),
    }


def thresholds(params=None):
    """The survey thresholds a params dict records, schema defaults filling the rest."""
    p = normalize_params(params)
    return {k: p[k] for k in THRESHOLD_KEYS}


def load_params(path):
    """Load params from a standalone *_kd_params.json OR a saved *_kd_fit_data.json
    (params live under its 'metadata'). Returns the normalized params dict."""
    with open(path) as fh:
        raw = json.load(fh)
    if isinstance(raw, dict) and isinstance(raw.get('metadata'), dict):
        raw = raw['metadata']           # a fit JSON — params are in metadata
    return normalize_params(raw)


def recorded_keys(path):
    """The keys a params source actually states, before load_params fills the defaults
    in: a hand-placed file that says nothing about conc_units has not declared 'absolute'."""
    with open(path) as fh:
        raw = json.load(fh)
    if isinstance(raw, dict) and isinstance(raw.get('metadata'), dict):
        raw = raw['metadata']
    return set(raw) if isinstance(raw, dict) else set()


def dump_params_json(path, params):
    """Write normalized params to a JSON file (indent=2, hand-editable). A conc_units of
    None means nobody stated the units: the file then records none, rather than the
    schema default every reader would take for a statement."""
    out = normalize_params(params)
    if 'conc_units' in params and params['conc_units'] is None:
        del out['conc_units']
    with open(path, 'w') as fh:
        json.dump(out, fh, indent=2)


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


def points_mismatch(stored, points):
    """Why the titration points a params source records are not `points`, or None."""
    recorded = stored.get('points')
    if recorded is not None and [float(p) for p in recorded] != [float(p) for p in points]:
        return f"its titration points {recorded} are not this input's {list(points)}"
    return None


def _mismatch(path, stored, csv_path, points):
    """Why a stored params file does not describe this input, or None if it does."""
    reason = points_mismatch(stored, points)
    if reason:
        return reason
    name = os.path.basename(csv_path)
    if stored.get('input') is not None:
        return None if stored['input'] == name else f"it was written for {stored['input']}"
    if os.path.basename(path) != SIBLING_NAME:
        return "it does not record which input it was written for"
    return None


def find_params_source(csv_path, points=None, rejected=None):
    """Locate a params source near an input titration CSV, or None. Search order:
    a sibling kd_params.json, then the newest *_kd_params.json beside the CSV, then
    inside a sibling 'kd_analysis' output folder (newest *_kd_params.json, else the
    newest *_kd_fit_data.json).

    Given the input's `points`, only a file that describes THIS input is returned: one
    that loads, whose recorded points match, and that names this input - a hand-placed
    kd_params.json is the user's own statement for the folder and need not. Each file
    passed over is appended to `rejected` as {'path', 'reason'}, so the caller can say so.

    A run writes its machine-readable output under <out>/data/, so the kd_analysis
    candidates are checked at both depths: the bare paths only ever match output
    predating that layout, and dropping them would make an old analysis unreadable.
    """
    d = os.path.dirname(os.path.abspath(csv_path))
    analysis = os.path.join(d, "kd_analysis")
    candidates = [
        os.path.join(d, SIBLING_NAME),
        _newest(os.path.join(d, "*" + PARAMS_SUFFIX)),
        _newest(os.path.join(analysis, "data", "*" + PARAMS_SUFFIX)),
        _newest(os.path.join(analysis, "*" + PARAMS_SUFFIX)),
        _newest(os.path.join(analysis, "data", "*_kd_fit_data.json")),
        _newest(os.path.join(analysis, "*_kd_fit_data.json")),
    ]
    for c in candidates:
        if not (c and os.path.isfile(c)):
            continue
        if points is None:
            return c
        try:
            reason = _mismatch(c, load_params(c), csv_path, points)
        except Exception as exc:
            reason = f"it could not be read: {exc}"
        if reason is None:
            return c
        if rejected is not None:
            rejected.append({'path': c, 'reason': reason})
    return None
