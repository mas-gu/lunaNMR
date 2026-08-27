# ABOUTME: Per-residue evidence and a suggested residue selection for a Kd titration.
# ABOUTME: Surveys without fitting a global Kd, so a human picks which residues to fit.

import numpy as np

from kd_fit import _is_dummy, fit_residue_csp, fit_residue_intensity
from kd_input import csp_series, intensity_ratio_series
from kd_models import residue_sort_key

# A residue whose fitted plateau sits this far above its own largest observed CSP has
# not been constrained by the data — the titration never approached saturation for it.
# Cited as evidence, never used to drop a residue: the fraction flagged varies from 17%
# to 43% between two datasets, so it describes the experiment, not the residue.
DD_RUNAWAY_RATIO = 10.0

# Residues below this quantile of the dataset's own max-CSP distribution are treated as
# non-movers. Derived per dataset rather than fixed, because the noise floor is a
# property of the spectra.
NO_MOTION_QUANTILE = 0.25


def survey_residues(residues, L, P0, alpha=0.14, value='height',
                    noise_quantile=NO_MOTION_QUANTILE,
                    dd_runaway_ratio=DD_RUNAWAY_RATIO, ref_max_ratio=None):
    """Per-residue evidence for whether a residue is worth fitting.

    Runs provisional per-residue fits (needed for the plateau evidence) but never a
    global fit — the quotable Kd only ever comes from residues a human selected.

    Returns a list of dicts ordered by residue sequence number, each carrying the
    observables, the evidence flags, and a `verdict` of 'unusable' | 'check' | 'ok'.
    """
    L = np.asarray(L, dtype=float)
    rows = []
    for name, res in residues.items():
        csp = np.asarray(csp_series(res, alpha=alpha), dtype=float)
        ratio = np.asarray(intensity_ratio_series(
            res, value=value, ref_max_ratio=ref_max_ratio), dtype=float)
        finite_csp = csp[np.isfinite(csp)]
        max_csp = float(np.max(finite_csp)) if finite_csp.size else float('nan')

        csp_fit = fit_residue_csp(L, csp, P0)
        dd_ratio = float('nan')
        if csp_fit.get('success') and np.isfinite(max_csp) and max_csp > 0:
            dd_ratio = float(csp_fit['dd_max']) / max_csp

        finite_ratio = ratio[np.isfinite(ratio)]
        rows.append({
            'residue': name,
            'max_csp': max_csp,
            # I/I0 at the last point that survived, i.e. the largest perturbation
            # actually measured for this residue.
            'intensity_final': float(finite_ratio[-1]) if finite_ratio.size else float('nan'),
            'dd_ratio': dd_ratio,
            'n_failed_points': int(res.get('n_failed_points', 0)),
            'n_points_csp': int(np.count_nonzero(np.isfinite(csp))),
            'n_points_intensity': int(np.count_nonzero(np.isfinite(ratio))),
            'csp_r_squared': csp_fit.get('r_squared', float('nan')),
            'intensity_r_squared': fit_residue_intensity(
                L, ratio, P0).get('r_squared', float('nan')),
        })

    _add_verdicts(rows, noise_quantile=noise_quantile,
                  dd_runaway_ratio=dd_runaway_ratio)
    rows.sort(key=lambda r: residue_sort_key(r['residue']))
    return rows


def _add_verdicts(rows, noise_quantile=NO_MOTION_QUANTILE,
                  dd_runaway_ratio=DD_RUNAWAY_RATIO):
    """Attach a verdict and its reasons, using a noise floor derived from these rows."""
    movers = [r['max_csp'] for r in rows if np.isfinite(r['max_csp'])]
    floor = float(np.quantile(movers, noise_quantile)) if movers else 0.0

    for r in rows:
        r['noise_floor'] = floor
        unusable, check = [], []
        if _is_dummy(r['residue']):
            unusable.append('dummy')
        # intensity_ratio_series blanks every point of a residue whose reference is
        # unusable, so no surviving intensity point means the reference was rejected.
        if r['n_points_intensity'] == 0:
            unusable.append('reference_unusable')
        if r['n_points_csp'] < 2:
            unusable.append('no_position_data')

        # Ordered by how strongly each argues against fitting: a residue that never
        # moved is unfittable however clean it is, while one failed point of seven
        # leaves six good ones.
        if np.isfinite(r['max_csp']) and r['max_csp'] <= floor:
            check.append('no_motion')
        if r['n_failed_points']:
            check.append('failed_points')
        if np.isfinite(r['dd_ratio']) and r['dd_ratio'] > dd_runaway_ratio:
            check.append('plateau_unconstrained')

        r['reason_keys'] = unusable + check
        r['verdict'] = 'unusable' if unusable else ('check' if check else 'ok')
    return rows


_REASON_TEXT = {
    'dummy': lambda r: 'placeholder, not a residue',
    'reference_unusable': lambda r: 'reference point unusable',
    'no_position_data': lambda r: 'no usable position data',
    'no_motion': lambda r: (f"CSP {r['max_csp']:.4f} ppm at or below the noise floor "
                            f"({r['noise_floor']:.4f})"),
    'failed_points': lambda r: f"{r['n_failed_points']} point(s) lost to a failed fit",
    'plateau_unconstrained': lambda r: (f"plateau {r['dd_ratio']:.0f}x its own largest "
                                        f"CSP — not constrained by this titration"),
}


def reason_text(row):
    """Render a row's reason keys. Kept apart from the predicates so wording and
    firing can change independently."""
    return [_REASON_TEXT[k](row) for k in row['reason_keys']]


def prior_states(text, known):
    """Which residues an existing selection file selects, as name -> bool.

    A line counts only if its first token names a residue we actually surveyed, so
    the prose header cannot be mistaken for a deselected residue.
    """
    states = {}
    for line in text.splitlines():
        stripped = line.lstrip()
        body = stripped.lstrip('#').strip().split('#', 1)[0].split()
        if body and body[0] in known:
            states[body[0]] = not stripped.startswith('#')
    return states


def format_residues_file(rows, prior=None):
    """Render a survey as an editable selection file, preserving human decisions.

    Unusable residues are commented out. Everything else stays selected, with any
    evidence against it as a trailing comment, so the debatable calls are put in front
    of the reader rather than made for them.

    `prior` carries the selections read from an existing file. A residue the human has
    already ruled on keeps that ruling and is marked as an override where it disagrees
    with the current suggestion — re-running a survey must not quietly undo an edit.
    """
    prior = prior or {}
    lines_out = []
    for r in rows:
        suggested = r['verdict'] != 'unusable'
        selected = prior.get(r['residue'], suggested)
        why = reason_text(r)
        if r['residue'] in prior and selected != suggested:
            why.insert(0, 'survey suggests ' + ('excluding' if suggested is False
                                                else 'including'))
            mark = '[kept] ' if selected else '[dropped] '
        else:
            mark = ''
        if selected and r['verdict'] == 'unusable':
            # Selection expresses interest; a mechanical exclusion expresses data
            # validity. Wanting the residue cannot rescue it, so say so here rather
            # than letting the fit drop it silently.
            why.insert(0, 'data still unusable, will not fit')
        note = mark + '; '.join(why)
        if selected:
            lines_out.append(f"{r['residue']:<10} # {note}" if note else r['residue'])
        else:
            # The reason sits behind its own '#', so deleting the leading one leaves a
            # valid entry with a trailing comment rather than a corrupt residue name.
            lines_out.append(f"# {r['residue']:<10} # excluded: {note}")

    n_selected = sum(1 for r in rows
                     if prior.get(r['residue'], r['verdict'] != 'unusable'))
    header = [
        '# lunaNMR Kd residue selection',
        '#',
        '# Residues listed here are fitted. Comment a line out with # to exclude it;',
        '# delete the # to put one back. Trailing comments are evidence, not exclusions.',
        f'# {n_selected} of {len(rows)} residues selected.',
        '#',
    ]
    return '\n'.join(header + lines_out) + '\n'


def parse_residues_file(text):
    """Read a selection file or a comma-separated list into residue names.

    Everything from a '#' onwards is a comment, so a fully commented line selects
    nothing and a trailing comment leaves its residue selected.
    """
    if '\n' not in text and ',' in text:
        return [n.strip() for n in text.split(',') if n.strip()]
    names = []
    for line in text.splitlines():
        name = line.split('#', 1)[0].strip()
        if name:
            names.append(name)
    return names


_TABLE_COLUMNS = ['residue', 'max_csp', 'intensity_final', 'dd_ratio', 'n_failed_points',
                  'n_points_csp', 'n_points_intensity', 'csp_r_squared',
                  'intensity_r_squared', 'verdict']


def survey_table(rows):
    """The survey as a header plus rows, for the evidence CSV."""
    header = _TABLE_COLUMNS + ['reasons']
    out = [[r[c] for c in _TABLE_COLUMNS] + ['; '.join(reason_text(r))] for r in rows]
    return header, out


def run_kd_survey_with_params(params):
    """Survey a titration: write the editable selection file and the evidence table.

    Deliberately computes no global Kd. The number people quote comes only from a run
    over residues someone chose.
    """
    import csv
    import os

    from kd_fit import load_inputs

    data_in = load_inputs(params)
    rows = survey_residues(
        data_in['residues'], data_in['L'], data_in['P0'],
        alpha=data_in['alpha'], value=data_in['value'],
        noise_quantile=float(params.get('noise_quantile', NO_MOTION_QUANTILE)),
        dd_runaway_ratio=float(params.get('dd_runaway_ratio', DD_RUNAWAY_RATIO)),
        ref_max_ratio=params.get('ref_max_ratio'))

    os.makedirs(params['output_dir'], exist_ok=True)
    prefix = params.get('output_prefix', 'kd')
    # The selection is the durable artefact of the workflow, so it defaults beside the
    # INPUT rather than inside --out: a scratch output path changes between runs and the
    # merge below can only protect a file we actually find again.
    residues_file = params.get('selection_file') or os.path.join(
        os.path.dirname(os.path.abspath(params['input_csv_file'])),
        f"{prefix}_residues.txt")
    os.makedirs(os.path.dirname(residues_file) or '.', exist_ok=True)
    # Merge rather than overwrite: a re-survey refreshes the evidence but must never
    # discard a selection someone made from the last one.
    prior = {}
    if os.path.exists(residues_file):
        with open(residues_file) as fh:
            prior = prior_states(fh.read(), {r['residue'] for r in rows})
    with open(residues_file, 'w') as fh:
        fh.write(format_residues_file(rows, prior=prior))

    data_dir = os.path.join(params['output_dir'], 'data')
    os.makedirs(data_dir, exist_ok=True)
    survey_file = os.path.join(data_dir, f"{prefix}_survey.csv")
    header, table = survey_table(rows)
    with open(survey_file, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(table)

    # Persist the thresholds beside the INPUT, where find_params_source looks. Named
    # <prefix>_kd_params.json rather than the bare sibling kd_params.json so a
    # hand-placed file still wins the lookup and is never clobbered.
    from kd_params import PARAMS_SUFFIX, dump_params_json
    params_file = os.path.join(
        os.path.dirname(os.path.abspath(params['input_csv_file'])),
        f"{prefix}{PARAMS_SUFFIX}")
    # Persist the run's own params, not the loaded ones: data_in['concs'] has already
    # had equivalents multiplied by P0, so storing it beside a conc_units of
    # 'equivalents' would convert twice on read-back, and storing it as 'absolute'
    # would silently contradict a --conc-units the next run restates. Recording what
    # the caller gave keeps the numbers and their units describing the same thing.
    # normalize_params keeps the schema keys and drops the rest, so every setting a run
    # was judged with round-trips without being enumerated here.
    dump_params_json(params_file, {**params, 'points': data_in['points']})

    counts = {v: sum(1 for r in rows if r['verdict'] == v)
              for v in ('ok', 'check', 'unusable')}
    return {'rows': rows, 'residues_file': residues_file, 'survey_file': survey_file,
            'params_file': params_file,
            'n_total': len(rows), 'n_selected': len(rows) - counts['unusable'],
            'counts': counts, 'points': data_in['points'], 'L': data_in['L']}
