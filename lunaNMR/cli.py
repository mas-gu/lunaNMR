# ABOUTME: Unified command-line entry point for lunaNMR (`python -m lunaNMR <subcommand>`).
# ABOUTME: Dispatches to headless analysis engines; heavy deps are imported lazily per subcommand.

"""Command-line interface for lunaNMR.

Subcommands:
  kd        Fit binding affinity (Kd) from a titration series (CSP / intensity).
  series    Process a multi-spectrum series / titration.
  dynamixs  T1/T2 and methyl-T2 relaxation fitting.
  export    Render figures / reports from a fit JSON (headless).
  project   Inspect / prune a .lunaNMR project bundle.
  batch     Batch peak detection + Voigt/PS2D fitting over a folder of spectra.

Each subcommand imports its own dependencies lazily so that unrelated subcommands
(and ``import lunaNMR``) stay headless — no Qt, no display required for ``kd``.
"""

import argparse
import glob
import os
import re
import sys


def _float_list(text):
    """Parse a comma-separated list of floats, e.g. '0,10,25.5'."""
    return [float(x) for x in text.split(',') if x.strip() != '']


def _str_list(text):
    """Parse a comma-separated list of stripped strings, e.g. 'csp,intensity'."""
    return [x.strip() for x in text.split(',') if x.strip() != '']


def _choice_list(*allowed):
    """A comma-separated list whose members are checked, since argparse's own `choices`
    validates the whole string rather than the items inside it. A misspelling otherwise
    reached the renderer and produced no figures without saying why."""
    allowed_set = set(allowed)

    def parse(text):
        values = _str_list(text)
        unknown = [v for v in values if v not in allowed_set]
        if unknown:
            raise argparse.ArgumentTypeError(
                f"unknown value(s) {', '.join(unknown)}; choose from "
                f"{', '.join(sorted(allowed_set))}")
        return values
    return parse


# Bumped only when a key changes meaning or disappears. Additions do not bump it: a
# consumer that reads the keys it knows keeps working across them.
SCHEMA_VERSION = 1


def _json_safe(obj):
    """Replace non-finite floats with None.

    json.dumps writes bare NaN and Infinity, which are not JSON — jq, Go and Rust all
    reject them. `mean_t2` is NaN when nothing fitted, so the output became unparseable
    on exactly the degenerate runs worth inspecting.
    """
    import math
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def _emit(args, summary, *human_lines):
    """Print the run summary: a single JSON object under --format json, else human lines.

    `ok`, `schema_version` and the full `command` are injected here rather than written
    into every call site, so a summary cannot be added without them.
    """
    import json
    if getattr(args, 'format', 'text') == 'json':
        envelope = {'ok': True, 'schema_version': SCHEMA_VERSION,
                    'command': _full_command(args)}
        envelope.update(summary)
        print(json.dumps(_json_safe(envelope), allow_nan=False))
    else:
        for line in human_lines:
            print(line)


def _emit_error(args, exc):
    """Report a failure in the caller's own format.

    Under --format json the documented usage is `json.loads(proc.stdout)`; every failure
    path wrote to stderr and left stdout empty, so that parse raised on every error and
    the exception said nothing about the cause. Text mode keeps the bare stderr line.
    """
    import json
    if getattr(args, 'format', 'text') == 'json':
        print(json.dumps({'ok': False, 'schema_version': SCHEMA_VERSION,
                          'command': _full_command(args),
                          'error': {'type': type(exc).__name__, 'message': str(exc)}}))
    else:
        print(f"error: {exc}", file=sys.stderr)


import contextlib


@contextlib.contextmanager
def _engine_stdout(args):
    """Under --format json, send all engine chatter to stderr so the CLI's own stdout
    carries only the JSON summary.

    Redirects at both levels: the Python `sys.stdout` (so parent prints and pytest's
    capsys are covered) and the underlying fd 1 via os.dup2 (so worker processes spawned
    by the parallel processor, which inherit fd 1, also write to stderr rather than
    leaking onto stdout).
    """
    if getattr(args, 'format', 'text') != 'json':
        yield
        return
    sys.stdout.flush()
    saved_fd = os.dup(1)
    os.dup2(2, 1)
    try:
        with contextlib.redirect_stdout(sys.stderr):
            yield
    finally:
        sys.stdout.flush()
        os.dup2(saved_fd, 1)
        os.close(saved_fd)


def _full_command(args):
    """The full subcommand path, e.g. 'dynamixs density'.

    argparse records the nested choice under its own dest, so a nested handler reported
    only 'dynamixs' in a dry-run while its real run emitted 'dynamixs t1t2'. The one key
    every summary carries has to discriminate, or it cannot be used to route.
    """
    parts = [getattr(args, 'command', None),
             getattr(args, 'dynamixs_command', None),
             getattr(args, 'export_command', None),
             getattr(args, 'project_command', None)]
    return ' '.join(p for p in parts if p)


def _dry_run(args, inputs, planned):
    """Validate that required input paths exist and report the plan without running.

    `inputs` is a list of (label, path); `planned` a dict of planned outputs. Returns 0
    if every input exists, else 1.
    """
    missing = [p for _, p in inputs if not p or not os.path.exists(p)]
    summary = {
        'command': _full_command(args),
        'dry_run': True,
        'inputs': {label: path for label, path in inputs},
        'planned_outputs': planned,
        'missing_inputs': missing,
    }
    human = [f"[dry-run] {summary['command']}"]
    for label, path in inputs:
        ok = bool(path) and os.path.exists(path)
        human.append(f"  input  {label}: {path} [{'OK' if ok else 'MISSING'}]")
    for key, val in planned.items():
        human.append(f"  output {key}: {val}")
    _emit(args, summary, *human)
    return 1 if missing else 0


def _add_modules_path(*parts):
    """Put a dynamiXs module directory on sys.path so its top-level modules import.

    The dynamiXs fitters live in modules/ (sibling of the lunaNMR package) and use
    bare sibling imports (e.g. `from kd_models import ...`), so their directory must
    be importable directly. Mirrors modules/dynamiXs_v2o0/workers.py.
    """
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'modules', *parts))
    if path not in sys.path:
        sys.path.insert(0, path)
    # Also expose on PYTHONPATH so multiprocessing-spawn children (macOS default) and
    # subprocesses can import these sibling modules — the density analyzers fan out to a
    # Pool whose workers re-import the module by name.
    existing = os.environ.get('PYTHONPATH', '')
    if path not in existing.split(os.pathsep):
        os.environ['PYTHONPATH'] = path + (os.pathsep + existing if existing else '')
    return path


# Documented per-subcommand defaults for the delay units. They differ on purpose and are
# NOT unified: on t1t2/methyl-t2 the value only labels the output, while on t1rho and
# modelfree it rescales, so changing either default silently moves published numbers.
# Kept as constants because the flags default to None -- the handler has to tell "the
# user asked for seconds" from "the user said nothing", to know whether a sidecar may
# supply the answer.
_TIME_UNITS_DEFAULT = 's'          # dynamixs t1t2 / methyl-t2 (labels only)
_T1RHO_TIME_UNITS_DEFAULT = 'ms'   # dynamixs t1rho (rescales)
_MODELFREE_UNITS_DEFAULT = 'ms'    # dynamixs modelfree --f{1,2}-t{1,2}-units (rescales)

_SIDECAR_NAME = 'series_metadata.json'


def _series_sidecar(input_path):
    """The `series` run-metadata written beside a matrix, or None.

    `series` records the delay units it normalised to and which spectra had no
    parseable delay. Both were previously re-asserted by hand on the command line
    while the answer sat next to the input file.
    """
    import json
    path = os.path.join(os.path.dirname(os.path.abspath(input_path)), _SIDECAR_NAME)
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as fh:
            meta = json.load(fh)
    except (OSError, ValueError):
        return None            # a damaged sidecar must not block a runnable input
    if not isinstance(meta, dict):
        return None
    meta['path'] = path
    return meta


def _unparsed_spectra(sidecar):
    """The spectra whose filename carried no delay, named rather than counted."""
    return [c.get('spectrum') for c in (sidecar.get('columns') or [])
            if c.get('value') is None]


def _check_unparsed_delays(sidecar, allow):
    """Refuse a relaxation fit whose series lost spectra to unparseable filenames.

    Those columns carry no time, so they are dropped from the fit: the run measures
    fewer points than were acquired and says so nowhere the caller is looking. Two
    hetNOE planes are the legitimate case, hence the escape hatch.
    """
    n = int(sidecar.get('n_value_unparsed') or 0)
    if not n or allow:
        return
    names = [s for s in _unparsed_spectra(sidecar) if s]
    raise ValueError(
        f"{n} spectrum(s) in this series have no parseable delay, so they carry no "
        f"time and cannot be fitted: {', '.join(names) or 'see the sidecar'}. "
        f"Recorded in {sidecar['path']}. Rename them with an explicit unit "
        f"(_300ms, _2.4s, _2400us) and re-run `series`, or pass "
        f"--allow-unparsed-delays to fit the remaining points without them."
    )


def _resolve_time_units(flag, sidecar, default, flag_name='--time-units'):
    """Delay units: an explicit flag, else the sidecar's, else the documented default.

    A flag that contradicts the sidecar is refused rather than preferred. One of the
    two is wrong, and picking either way is how a T1 in seconds meets a T2 in
    milliseconds and puts R1 out by 1000x.
    """
    recorded = (sidecar or {}).get('value_units')
    if flag is None:
        return recorded or default
    if recorded and recorded != flag:
        raise ValueError(
            f"{flag_name} says {flag!r} but {sidecar['path']} records the delays as "
            f"{recorded!r}. One of them is wrong; guessing puts R1 out by 1000x. "
            f"Drop the flag to use the recorded units, or re-run `series` if the "
            f"sidecar is stale."
        )
    return flag


def _sidecar_summary(sidecar):
    """The sidecar fields worth reporting: which file was read, the units it supplied,
    and how many spectra it says had no delay. None when there was no sidecar, so a
    consumer can tell "checked and clean" from "nothing to check against"."""
    if not sidecar:
        return None
    return {'path': sidecar['path'], 'value_units': sidecar.get('value_units'),
            'n_value_unparsed': int(sidecar.get('n_value_unparsed') or 0)}


_PERSISTED_SETTINGS = ('conc_units', 'csp_sigma_multiple', 'kd_outlier_z',
                       'noise_quantile',
                       'dd_runaway_ratio', 'ref_max_ratio')

# The list-valued inputs, mapped schema key -> argparse dest. They round-trip too, so a
# survey's concentrations do not have to be retyped on the fit, but they differ from the
# settings above: an unrecorded one is stored as null rather than being absent, so the
# restore tests the value and not just the key.
_PERSISTED_INPUTS = {'concentrations': 'conc', 'intensity_scales': 'intensity_scale'}


def _apply_thresholds(args, params):
    """Resolve persisted binding params: an explicit flag beats a value persisted beside
    the input, which beats the schema default.

    Concentrations are stored as the caller typed them, next to the conc_units that
    describe them, so restoring the pair converts equivalents exactly once however the
    two are mixed between flags and the stored file.
    """
    from kd_params import find_params_source, load_params
    stored = {}
    src = find_params_source(args.input)
    if src:
        try:
            stored = load_params(src)
        except Exception:
            src, stored = None, {}
    for key in _PERSISTED_SETTINGS:
        flag = getattr(args, key, None)
        if flag is not None:
            params[key] = flag
        elif key in stored:
            params[key] = stored[key]
    for key, dest in _PERSISTED_INPUTS.items():
        flag = getattr(args, dest, None)
        if flag is not None:
            params[key] = flag
        elif stored.get(key) is not None:
            params[key] = stored[key]
    return src


_DATA_SUBDIR = 'data'   # CSV/JSON companions sit one level under the figures


def _data_dir(out):
    """Machine-readable companions live under <out>/data; the top level stays figures."""
    d = os.path.join(out, _DATA_SUBDIR)
    os.makedirs(d, exist_ok=True)
    return d


def _read_residue_selection(value):
    """Resolve --residues: a selection file if it names one, else a comma-separated list."""
    from kd_survey import parse_residues_file
    text = value
    if os.path.isfile(value):
        with open(value) as fh:
            text = fh.read()
    names = parse_residues_file(text)
    if not names:
        raise ValueError(f"no residues selected in {value}")
    return names


def _draw_vs_sequence(ax, rows, key, ylabel, floor=None):
    """Per-residue bars against sequence, greying out residues the survey rejected."""
    import numpy as np
    names = [r['residue'] for r in rows]
    vals = [r[key] if np.isfinite(r[key]) else 0.0 for r in rows]
    colors = ['#d0d0d0' if r['verdict'] == 'unusable'
              else ('#e8a33d' if r['verdict'] == 'check' else 'steelblue') for r in rows]
    ax.bar(range(len(rows)), vals, color=colors)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(names, rotation=90, fontsize=5)
    ax.set_ylabel(ylabel)
    if floor is not None and np.isfinite(floor):
        ax.axhline(floor, color='red', ls=':', lw=0.8,
                   label=f"noise floor {floor:.4f}")
        ax.legend(fontsize=7)


def _run_kd_survey(args):
    """Survey a titration and render the vs-sequence figures. Writes no fit JSON."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from kd_survey import run_kd_survey_with_params

    params = {
        'input_csv_file': args.input,
        'output_dir': args.out,
        'output_prefix': args.prefix,
        'protein_conc': args.p0,
        'alpha': args.alpha,
        'intensity_value': args.intensity_from,
    }
    if args.selection is not None:
        params['selection_file'] = args.selection
    _apply_thresholds(args, params)

    with _engine_stdout(args):
        result = run_kd_survey_with_params(params)
    rows = result['rows']

    matplotlib.rcParams['pdf.fonttype'] = 42
    movers = [r['max_csp'] for r in rows if np.isfinite(r['max_csp'])]
    floor = float(np.quantile(movers, 0.25)) if movers else None
    figures = []
    for key, ylabel, tag, ref in (
            ('max_csp', 'max CSP (ppm)', 'csp', floor),
            ('intensity_final', 'I/I$_0$ at last point', 'intensity', None)):
        fig, ax = plt.subplots(figsize=(max(8.0, len(rows) * 0.16), 4.0))
        _draw_vs_sequence(ax, rows, key, ylabel, floor=ref)
        ax.set_title(f"{tag} vs sequence — {result['n_selected']}/{result['n_total']} selected")
        fig.tight_layout()
        path = os.path.join(args.out, f"{args.prefix}_{tag}_vs_sequence.pdf")
        fig.savefig(path)
        plt.close(fig)
        figures.append(path)

    counts = result['counts']
    _emit(args,
          {'command': 'kd survey', 'n_total': result['n_total'],
           'n_selected': result['n_selected'], 'counts': counts,
           'residues_file': result['residues_file'],
           'survey_file': result['survey_file'], 'figures': figures},
          f"Survey complete: {result['n_selected']}/{result['n_total']} residues selected "
          f"({counts['ok']} ok, {counts['check']} to check, {counts['unusable']} unusable)",
          f"  Selection: {result['residues_file']}",
          f"  Evidence:  {result['survey_file']}",
          *(f"  Figure:    {p}" for p in figures),
          "  Edit the selection file, then re-run with --residues to fit.")
    return 0


def _run_kd(args):
    """Wrap dynamiXs_Kd.kd_fit.run_kd_analysis_with_params for the CLI."""
    if args.dry_run:
        # The survey writes its selection file beside the INPUT, not into --out, so a
        # plan naming only --out sends the caller to the wrong place for the one file
        # it is then told to edit.
        planned = {'output_dir': args.out}
        if args.survey:
            planned['residues_file'] = os.path.join(
                os.path.dirname(os.path.abspath(args.input)), f"{args.prefix}_residues.txt")
            planned['survey_file'] = os.path.join(args.out, 'data',
                                                  f"{args.prefix}_survey.csv")
        return _dry_run(args, [('input', args.input)], planned)
    _add_modules_path('dynamiXs_v2o0', 'dynamiXs_Kd')
    if args.survey:
        return _run_kd_survey(args)
    import kd_fit

    params = {
        'input_csv_file': args.input,
        'output_dir': args.out,
        'output_prefix': args.prefix,
        'protein_conc': args.p0,
        'alpha': args.alpha,
        'observables': args.observable,
        'intensity_value': args.intensity_from,
        'n_bootstrap': args.bootstrap,
    }
    if args.residues is not None:
        params['residues'] = _read_residue_selection(args.residues)
    params_source = _apply_thresholds(args, params)

    with _engine_stdout(args):
        result = kd_fit.run_kd_analysis_with_params(params)
    for w in result.get('quality_warnings') or []:
        print('\n' + '!' * 72, file=sys.stderr)
        print(w, file=sys.stderr)
        print('!' * 72 + '\n', file=sys.stderr)
    _emit(args,
          {'command': 'kd', 'n_fitted': result['n_fitted'], 'n_total': result['n_total'],
           'quality_warnings': result.get('quality_warnings') or [],
           'conc_units': params.get('conc_units', 'absolute'),
           'params_source': params_source,
           'params_file': result.get('params_file'),
           'json_file': result['json_file'], 'results_file': result['results_file']},
          f"Kd analysis complete: {result['n_fitted']}/{result['n_total']} residues fitted",
          *([f"  Settings restored from: {params_source}"] if params_source else []),
          f"  JSON:    {result['json_file']}",
          f"  Params:  {result.get('params_file')}",
          f"  Results: {result['results_file']}")
    return 0


def _run_dynamixs_t1t2(args):
    """Wrap dynamiXs_T1_T2.fit_Tx_NMRRE.run_analysis_with_params (T1/T2 relaxation)."""
    if args.dry_run:
        return _dry_run(args, [('input', args.input)],
                        {'output_dir': args.out, 'experiment': args.exp})
    _add_modules_path('dynamiXs_v2o0', 'dynamiXs_T1_T2')
    from fit_Tx_NMRRE import run_analysis_with_params
    sidecar = _series_sidecar(args.input)
    if sidecar:
        _check_unparsed_delays(sidecar, args.allow_unparsed)
    time_units = _resolve_time_units(args.time_units, sidecar, _TIME_UNITS_DEFAULT)
    os.makedirs(args.out, exist_ok=True)
    params = {
        'input_csv_file': args.input,
        'output_prefix': os.path.join(args.out, args.prefix),
        'results_txt_file': os.path.join(args.out, f"{args.prefix}_fit_results.txt"),
        'experiment_type': args.exp,
        'time_units': time_units,
        'error_method': args.error_method,
        'n_bootstrap': args.bootstrap,
        'field_name': args.field_name,
        'field_freq': args.field_freq,
        'json_folder': None if args.no_json else args.out,
    }
    with _engine_stdout(args):
        result = run_analysis_with_params(params)
    n_total = result['n_fitted'] + result.get('n_excluded', 0)
    human = [f"{args.exp} analysis complete: {result['n_fitted']}/{n_total} residues fitted, "
             f"mean {args.exp} = {result['mean_t2']:.2f} {time_units}",
             f"  Results: {result['results_file']}"]
    if result.get('json_file'):
        human.append(f"  JSON:    {result['json_file']}")
    marginal = result.get('n_window_marginal')
    if marginal:
        human.append(f"  NOTE: {marginal} residues have t_max/T < 2 — the delay window is "
                     f"too short to pin the time constant down for those")
    _emit(args,
          {'command': 'dynamixs t1t2', 'experiment': args.exp, 'n_fitted': result['n_fitted'],
           'time_units': time_units, 'series_metadata': _sidecar_summary(sidecar),
           'dropped_columns': result.get('dropped_columns') or [],
           # n_fitted alone is not a rate. The engine already returns the excluded
           # count; dropping it here left 40 survivors of 200 reading as 40.
           'n_excluded': result.get('n_excluded', 0),
           'n_total': result['n_fitted'] + result.get('n_excluded', 0),
           'mean_t2': result['mean_t2'], 'results_file': result['results_file'],
           'json_file': result.get('json_file'),
           'baseline_fixed': result.get('baseline_fixed', True),
           'n_window_marginal': marginal},
          *human)
    return 0 if result['n_fitted'] else 1


def _run_dynamixs_methyl(args):
    """Wrap dynamiXs_T1_T2.fit_methyl_T2.run_methyl_t2_analysis_with_params (bi-exp methyl T2)."""
    if args.dry_run:
        return _dry_run(args, [('input', args.input)], {'output_dir': args.out})
    _add_modules_path('dynamiXs_v2o0', 'dynamiXs_T1_T2')
    from fit_methyl_T2 import run_methyl_t2_analysis_with_params
    sidecar = _series_sidecar(args.input)
    if sidecar:
        _check_unparsed_delays(sidecar, args.allow_unparsed)
    time_units = _resolve_time_units(args.time_units, sidecar, _TIME_UNITS_DEFAULT)
    os.makedirs(args.out, exist_ok=True)
    params = {
        'input_csv_file': args.input,
        'output_prefix': os.path.join(args.out, args.prefix),
        'results_txt_file': os.path.join(args.out, f"{args.prefix}_fit_results.txt"),
        'json_folder': args.out,
        'field_name': args.field_name,
        'field_freq': args.field_freq,
        'time_units': time_units,
        'error_method': args.error_method,
        'n_bootstrap': args.bootstrap,
    }
    with _engine_stdout(args):
        result = run_methyl_t2_analysis_with_params(params)
    _emit(args,
          {'command': 'dynamixs methyl-t2', 'n_fitted': result['n_fitted'],
           'time_units': time_units, 'series_metadata': _sidecar_summary(sidecar),
           'dropped_columns': result.get('dropped_columns') or [],
           'n_total': result['n_total'], 'results_file': result['results_file'],
           'json_file': result['json_file']},
          f"Methyl T2 analysis complete: {result['n_fitted']}/{result['n_total']} residues fitted",
          f"  Results: {result['results_file']}",
          f"  JSON:    {result['json_file']}")
    return 0 if result['n_fitted'] else 1


def _natural_key(text):
    """Sort key that orders embedded numbers numerically (s_2 before s_10)."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', text)]


def _discover_spectra(spectra, extensions=None):
    """Resolve --spectra (a folder or a glob) to a naturally-sorted list of spectrum files.

    `extensions` defaults to NMRFileManager.supported_nmr_formats, resolved at call time
    rather than restated here: a literal default is a second list, and the one that used
    to sit here had fallen behind by an extension, so NMRPipe 1D spectra were invisible
    to every caller that took the default.
    """
    if extensions is None:
        from lunaNMR.utils.file_manager import NMRFileManager
        extensions = NMRFileManager().supported_nmr_formats
    if os.path.isdir(spectra):
        files = []
        for ext in extensions:
            files.extend(glob.glob(os.path.join(spectra, f'*.{ext}')))
    else:
        files = glob.glob(spectra)
    return sorted(files, key=lambda f: _natural_key(os.path.basename(f)))


def _default_series_params(parallel=False):
    """Default series/voigt parameters, mirroring the GUI's getattr fallbacks.

    The GUI assembles this nested dict from the main window with these same defaults
    (gui/dialogs/series_integration_dialog.py::_get_voigt_parameters); a headless run
    just uses the defaults directly. `parallel` enables the two-pass parallel processor.
    """
    return {
        'detection_params': {
            'search_window_x': 0.070,
            'search_window_y': 0.070,
            'noise_threshold': 3.0,
        },
        'gui_params': {
            'fix_positions': False,
            'fix_linewidths': False,
            'use_parallel_processing': parallel,
            'use_centroid_refinement': True,
            'centroid_window_x_ppm': 0.02,
            'centroid_window_y_ppm': 1.0,
            'centroid_noise_multiplier': 2.0,
            'use_ps2d_multi_peak': True,
            'use_ps2d_linewidth_reuse': False,
            'collect_training_data': False,
            'height_threshold': 0.1,
            'distance_factor': 2.0,
            'prominence_threshold': 0.05,
            'smoothing_sigma': 1.0,
            'max_peaks_fit': 50,
            'max_optimization_iterations': 50,
        },
        'fitting_params': {
            'min_r_squared': 0.5,
            'max_iterations': 1000,
            'fitting_window_x': 0.2,
            'fitting_window_y': 2.0,
        },
        'processing_options': {
            'use_parallel_processing': parallel,
            'use_global_optimization': False,
            'enable_cascade_drift_limit': True,
            'rerun_adaptive_per_spectrum': False,
            'lock_cluster_assignments': False,
            'use_original_reference_for_detection': False,
        },
    }


def _serialise_checks(result):
    """Findings are NamedTuples and the peak list carries its coordinates; neither belongs
    in a JSON summary as-is."""
    out = {}
    for key, value in result.items():
        if key == 'findings':
            out[key] = [f._asdict() for f in value]
        elif key == 'peak_list' and isinstance(value, dict):
            out[key] = {'path': value.get('path'), 'n_peaks': len(value.get('peaks', [])),
                        'assignments': value.get('assignments', [])}
        elif key == 'experiments':
            out[key] = [_serialise_checks(e) for e in value]
        elif isinstance(value, dict) and 'findings' in value:
            out[key] = {k: v for k, v in value.items() if k != 'findings'}
        else:
            out[key] = value
    return out


def _format_experiment(check):
    """The per-experiment block of the human-readable report."""
    lines = [f"\n### {check['experiment']}  —  {check['n_spectra']} spectra"]
    delays = check.get('delays')
    if delays and delays['values']:
        lines.append(f"    delays     : {delays['values']}")
        if delays['repeats']:
            lines.append(f"    repeats    : {delays['repeats']}")
    if delays and delays['unparsed']:
        lines.append(f"    unparsed   : {len(delays['unparsed'])} spectra -> stem-named columns")
    peak_list = check.get('peak_list')
    if peak_list:
        lines.append(f"    peak list  : {os.path.basename(peak_list['path'])} "
                     f"({len(peak_list['peaks'])} peaks)")
    for spec in check['spectra']:
        flag = ''
        offset = max(abs(spec['dx']), abs(spec['dy'])) / 0.070
        if offset > 0.5:
            flag = '   <-- OFFSET'
        elif offset > 0.15:
            flag = '   (minor)'
        if spec['decayed']:
            flag += '   (decayed)'
        lines.append(f"    {spec['spectrum']:34s} shift=({spec['dx']:+.4f},{spec['dy']:+.3f}) "
                     f"capture={spec['capture']}/{spec['n_peaks']} "
                     f"medS/N={spec['median_snr']:.0f} weak={spec['weak']}{flag}")
    if check.get('hetnoe'):
        h = check['hetnoe']
        lines.append(f"    hetNOE sat/unsat = {h['ratio']:.3f}   SATURATED = {h['saturated']}"
                     f"   UNSATURATED = {h['unsaturated']}")
    if check.get('subseries'):
        sub_ = check['subseries']
        lines.append(f"    repeat acquisitions, shared delays {sub_['shared']}:")
        for delay, rec in sorted(sub_['per_delay'].items()):
            lines.append(f"        {delay:g}: {rec['second']}/{rec['first']} = {rec['ratio']:.3f}"
                         f"  (n={rec['n_strong']})")
    return lines


def _format_findings(findings):
    if not findings:
        return ['', 'NO ISSUES — safe to run the pipeline.']
    lines = ['']
    for severity in ('FAIL', 'WARN'):
        for f in findings:
            if f.severity == severity:
                lines.append(f'  [{severity}] {f.message}')
    n_fail = sum(1 for f in findings if f.severity == 'FAIL')
    n_warn = len(findings) - n_fail
    lines.append(f'\n{n_fail} FAIL, {n_warn} WARN — review before running the pipeline.')
    return lines


def _findings_exit_code(findings, strict=False):
    """0 when it is safe to run the pipeline, 1 when it is not.

    FAIL always gates: `diagnose && series` exists so the shell stops there, and
    returning 0 regardless made it sail past the exact failure the command is for.
    WARN does not, by default -- it includes routine conditions like a folder with no
    peak list, and a gate that fires on those gets deleted rather than heeded.
    """
    severities = {f.severity for f in findings}
    if 'FAIL' in severities:
        return 1
    return 1 if (strict and 'WARN' in severities) else 0


def _run_diagnose(args):
    """Read-only pre-flight over a whole dataset: registration, capture, delays, peak lists."""
    from lunaNMR.validation.spectra_check import check_dataset
    root = os.path.abspath(args.root)
    if not os.path.isdir(root):
        raise FileNotFoundError(f"dataset root not found: {root}")
    result = check_dataset(root, quick=args.quick, sample=args.sample, series_mode=args.mode)
    human = [f'DIAGNOSTIC — {root}', '=' * 72]
    for check in result['experiments']:
        human.extend(_format_experiment(check))
    if result.get('assignments'):
        human.append('\n### assignment sets')
        for name, n in result['assignments']['per_experiment'].items():
            human.append(f'    {name:28s} {n} residues')
        human.append(f"    common to all: {result['assignments']['common']}   "
                     f"union: {result['assignments']['union']}")
    human.append('\n' + '=' * 72)
    human.extend(_format_findings(result['findings']))
    _emit(args, {'command': 'diagnose', **_serialise_checks(result)}, *human)
    return _findings_exit_code(result['findings'], strict=args.strict)


def _validate_series(parser, args):
    """`--deep` inspects the spectra instead of running them, so it needs `--dry-run`."""
    if getattr(args, 'deep', False) and not args.dry_run:
        parser.error('--deep only applies to --dry-run (it inspects the spectra rather than '
                     'processing them)')


def _run_series(args):
    """Wrap MultiSpectrumProcessor.process_nmr_series for a headless series/titration run."""
    import matplotlib
    matplotlib.use('Agg')  # headless guard before any plotting-capable import
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    from lunaNMR.utils.file_manager import NMRFileManager

    file_manager = NMRFileManager()
    nmr_files = _discover_spectra(args.spectra, file_manager.supported_nmr_formats)

    if args.dry_run:
        peaks_ok = os.path.exists(args.peaks)
        missing = ([] if nmr_files else ['spectra']) + ([] if peaks_ok else [args.peaks])
        summary = {'command': 'series', 'dry_run': True, 'spectra_found': len(nmr_files),
                   'peaks': args.peaks, 'peaks_exists': peaks_ok, 'mode': args.mode,
                   'peak_source': args.peak_source, 'parallel': args.parallel,
                   'output_dir': args.out, 'missing_inputs': missing}
        deep_lines, deep_code = [], 0
        if getattr(args, 'deep', False) and nmr_files and peaks_ok:
            from lunaNMR.validation.spectra_check import check_experiment
            checks = check_experiment(os.path.dirname(nmr_files[0]), args.peaks,
                                      quick=args.quick, sample=args.sample,
                                      series_mode=args.mode)
            summary['checks'] = _serialise_checks(checks)
            deep_lines = _format_experiment(checks) + _format_findings(checks['findings'])
            # --deep runs the same checks as `diagnose`, so a FAIL means the same
            # thing: do not run this series. A dry-run that reports one and exits 0
            # is the gap `diagnose` had.
            deep_code = _findings_exit_code(checks['findings'])
        _emit(args, summary,
              "[dry-run] series",
              f"  spectra found: {len(nmr_files)}",
              f"  peaks: {args.peaks} [{'OK' if peaks_ok else 'MISSING'}]",
              f"  mode={args.mode} peak-source={args.peak_source} parallel={args.parallel}",
              f"  output: {args.out}",
              *deep_lines)
        return 1 if missing else deep_code

    if not nmr_files:
        _emit_error(args, FileNotFoundError(f"No spectrum files found in {args.spectra}"))
        return 1
    reference_peaks = file_manager.load_peak_list(args.peaks)
    if reference_peaks.empty:
        _emit_error(args, ValueError(f"Peak list is empty: {args.peaks}"))
        return 1

    os.makedirs(args.out, exist_ok=True)
    print(f"Processing {len(nmr_files)} spectra ({args.mode} mode, {args.peak_source} peaks)...",
          file=sys.stderr)
    series_params = _default_series_params(parallel=args.parallel)
    if args.peak_source == 'independent':
        # Independent mode means a full detect+fit per spectrum, which is two options
        # rather than one. The GUI sets both (series_integration_dialog, the
        # independent_radio branch); the CLI accepted the same flag name and left them
        # at their cascade defaults, so the two ran different algorithms in silence.
        series_params['processing_options'].update(
            rerun_adaptive_per_spectrum=True,
            use_original_reference_for_detection=True)
    processor = MultiSpectrumProcessor(series_params)
    with _engine_stdout(args):
        result = processor.process_nmr_series(
            nmr_files, reference_peaks, args.out,
            peak_source_mode=args.peak_source, series_mode=args.mode, extract_delays=True,
        )
    if getattr(result, 'errors', None):
        for err in result.errors:
            print(f"  error: {err}", file=sys.stderr)
    # process_nmr_series reports failure by returning a result with no successful
    # spectra (or an empty result on a top-level exception), never by raising.
    results = getattr(result, 'results', None) or {}
    n_success = result.metadata.get('successful_spectra', 0)
    if not n_success:
        _emit_error(args, RuntimeError(
            f"Series produced no successful fits: 0 of {len(results)} spectra "
            f"(all failed, or none loaded)"))
        return 1
    output_folder = result.metadata.get('output_folder', args.out)
    # Naming them beats globbing for them, and series_analysis_tidy.csv is the only
    # place per-peak fit quality reaches the series output — a caller that misses it
    # has no way to gate residues.
    outputs = {name: os.path.join(output_folder, name)
               for name in ('peak_intensity_matrix.csv', 'peak_volume_matrix.csv',
                            'peak_detected_matrix.csv', 'comprehensive_peak_tracking.csv',
                            'series_analysis_tidy.csv', 'series_metadata.json')
               if os.path.exists(os.path.join(output_folder, name))}
    _emit(args,
          {'command': 'series', 'spectra_fitted': n_success, 'spectra_total': len(results),
           'output_folder': output_folder, 'outputs': outputs, 'parallel': args.parallel},
          f"Series analysis complete: {n_success}/{len(results)} spectra fitted",
          f"  Output: {output_folder}",
          *(f"  {name}" for name in sorted(outputs)))
    return 0


def _human_size(n):
    """Format a byte count as a short human-readable string."""
    size = float(n)
    for unit in ('B', 'KB', 'MB', 'GB'):
        if size < 1024 or unit == 'GB':
            return f"{size:.0f} {unit}" if unit == 'B' else f"{size:.1f} {unit}"
        size /= 1024


def _project_manager():
    """A headless ProjectManager: inventory/remove only touch the filesystem, so a
    duck-typed session object suffices (no Qt main window needed)."""
    import types
    from lunaNMR.utils.project_manager import ProjectManager
    return ProjectManager(types.SimpleNamespace())


def _run_project_inventory(args):
    if not os.path.isdir(args.bundle):
        _emit_error(args, FileNotFoundError(f"Not a project bundle: {args.bundle}"))
        return 1
    categories = _project_manager().inventory(args.bundle)

    def _lines():
        """The listing, with the bundle-relative paths `project remove` actually takes.

        They were in the inventory data all along and only the English label was
        printed, so there was no way to learn a valid argument from the tool itself.
        """
        if not categories:
            return ["(empty bundle — no recognized categories)"]
        out = []
        for cat in categories:
            out.append(f"{cat['label']}  [{_human_size(cat['size'])}]")
            for item in cat['items']:
                if item['removable']:
                    out.append(f"  - {item['label']}  [{_human_size(item['size'])}]")
                    out.append(f"      {' '.join(item['paths'])}")
                else:
                    out.append(f"  - {item['label']}  [{_human_size(item['size'])}]"
                               f"  (protected)")
        return out

    _emit(args, {'command': 'project inventory', 'bundle': args.bundle,
                 'categories': categories}, *_lines())
    return 0


def _run_project_remove(args):
    if not os.path.isdir(args.bundle):
        _emit_error(args, FileNotFoundError(f"Not a project bundle: {args.bundle}"))
        return 1
    manager = _project_manager()
    if args.dry_run:
        # Deletion here is immediate and has no undo, so the preview must apply the
        # same containment check the real run does -- a preview that accepts a path
        # the run would refuse is worse than none.
        try:
            targets = manager.resolve_bundle_paths(args.bundle, args.paths)
        except ValueError as exc:
            _emit_error(args, exc)
            return 1
        _emit(args, {'command': 'project remove', 'dry_run': True,
                     'bundle': args.bundle, 'paths': args.paths,
                     'existing': [p for p, exists, _ in targets if exists],
                     'missing': [p for p, exists, _ in targets if not exists],
                     'would_free': sum(size for _, _, size in targets)},
              f"[dry-run] project remove — nothing deleted",
              *(f"  {'would remove' if exists else 'not present'}: {rel}"
                f"{'  [' + _human_size(size) + ']' if exists else ''}"
                for rel, exists, size in targets),
              f"  would free {_human_size(sum(s for _, _, s in targets))}")
        return 0
    try:
        freed = manager.remove_bundle_paths(args.bundle, args.paths)
    except ValueError as exc:
        _emit_error(args, exc)
        return 1
    _emit(args, {'command': 'project remove', 'bundle': args.bundle,
                 'paths': args.paths, 'freed': freed},
          f"Removed {len(args.paths)} path(s), freed {_human_size(freed)}")
    return 0


def _safe_name(text):
    """Filesystem-safe version of a residue label."""
    return re.sub(r'[^\w.+-]+', '_', str(text))


# Shared y-axis (with equal top/bottom margin) for every intensity curve panel, so
# raw-scale residues (I0 in the hundreds) and already-normalized ones plot on the same
# 0-1 I/I(0) scale instead of each auto-scaling to its own data range.
_INTENSITY_YLIM = (-0.05, 1.05)


def _build_kd_panel(residue, fit, L, y, Ld, obs, P0):
    """Build one residue's curve-fit panel dict. Intensity panels are normalized to
    I/I(0) (dividing by the fitted I0 amplitude) so every intensity panel shares a
    meaningful 0-1 y-axis regardless of the input data's raw scale."""
    from kd_models import csp_model, intensity_decay
    if obs == 'csp':
        yc = csp_model(Ld, fit['dd_max'], fit['Kd'], P0)
        ylabel = 'CSP (ppm)'
    else:
        yc = intensity_decay(Ld, fit['I0'], fit['I_inf'], fit['Kd'])
        y = y / fit['I0']
        yc = yc / fit['I0']
        ylabel = 'I / I(0)'
    return {'residue': residue, 'L': L, 'y': y, 'Ld': Ld, 'yc': yc, 'ylabel': ylabel,
            'kd': fit['Kd'], 'r2': fit.get('r_squared', float('nan'))}


def _draw_kd_panel(ax, panel, ylim=None):
    """Draw one residue's observed points + fitted binding curve into an axis."""
    ax.plot(panel['L'], panel['y'], 'o', color='#1f77b4', ms=4)
    ax.plot(panel['Ld'], panel['yc'], '-', color='#d62728', lw=1.2)
    ax.set_title(f"{panel['residue']}  Kd={panel['kd']:.2g}  R²={panel['r2']:.2f}", fontsize=8)
    ax.set_xlabel('[ligand]', fontsize=7)
    ax.set_ylabel(panel['ylabel'], fontsize=7)
    ax.tick_params(labelsize=6)
    if ylim is not None:
        ax.set_ylim(*ylim)


def _draw_ref_bars(ax, names, vals, obs, ref_label, cmp_label, threshold=None,
                   threshold_label=None):
    """Bar chart of a reference→point observable per residue. Residues absent at either
    point (NaN) draw as a full-height grey bar so the gap stays visible. I/I₀ is a ratio
    (axis fixed 0–1); CSP is in ppm (auto-scaled)."""
    import numpy as np
    finite = [v for v in vals if np.isfinite(v)]
    if obs == 'intensity':
        data_max = max(finite + [1.0]) if finite else 1.0
        top = 1.0 if data_max <= 1.0 else data_max * 1.05
        color, ylabel, sym = 'indianred', 'Intensity ratio I/I₀', 'I/I₀'
    else:
        top = (max(finite) * 1.1) if finite else 1.0
        color, ylabel, sym = 'seagreen', 'CSP (ppm)', 'CSP'
    if threshold is not None and np.isfinite(threshold):
        top = max(top, float(threshold) * 1.15)
    top = top if top > 0 else 1.0          # all-zero observable -> avoid a singular axis
    x = np.arange(len(names))
    heights = [v if np.isfinite(v) else top for v in vals]
    colors = [color if np.isfinite(v) else '#d0d0d0' for v in vals]
    ax.bar(x, heights, color=colors)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_ylim(0, top)
    ax.set_ylabel(ylabel)
    if threshold is not None and np.isfinite(threshold):
        # The gate's own threshold, read from the fit rather than recomputed, so the line
        # and the pool can never disagree. It is defined at the last titration point, so
        # on earlier pages it shows how far the shifts still are from clearing it.
        ax.axhline(float(threshold), color='red', ls='--', lw=1.2,
                   label=threshold_label or f"threshold {threshold:.4f}")
        ax.legend(fontsize=7, loc='upper right')

    def _lab(v):
        return f"{v:g}" if isinstance(v, (int, float)) else str(v)
    ax.set_title(f"{sym}:  {_lab(ref_label)} → {_lab(cmp_label)}  (ref → point)")


def _draw_kd_bars(ax, names, kds, errs, obs, global_kd):
    """Per-residue Kd bar chart. Residues with no successful fit (NaN) draw as a
    full-height grey bar; the shared global Kd (if any) is a red dashed line."""
    import numpy as np
    finite = [v for v in kds if np.isfinite(v)]
    top = (max(finite) * 1.15) if finite else 1.0
    top = top if top > 0 else 1.0          # all-zero Kd -> avoid a singular axis
    x = np.arange(len(names))
    heights = [v if np.isfinite(v) else top for v in kds]
    colors = ['steelblue' if np.isfinite(v) else '#d0d0d0' for v in kds]
    yerr = [e if np.isfinite(v) else 0.0 for v, e in zip(kds, errs)]
    ax.bar(x, heights, yerr=yerr, color=colors, ecolor='#555555', capsize=2)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_ylim(0, top)
    ax.set_ylabel('Kd')
    lbl = 'CSP' if obs == 'csp' else 'intensity (apparent)'
    if global_kd is not None and np.isfinite(global_kd):
        gl = 'global Kd' if obs == 'csp' else 'global apparent Kd'
        ax.axhline(global_kd, color='red', ls='--', label=f"{gl}={global_kd:.3g}")
        ax.legend()
    ax.set_title(f"Kd vs residue  ({lbl})")


def _draw_global_panel(ax, panel, ylim=None):
    """One residue's observed points + the shared-Kd global-model curve. Title carries
    the per-residue R² of the data against the global curve (goodness under one Kd)."""
    ax.plot(panel['L'], panel['y'], 'o', color='black', ms=4)
    ax.plot(panel['Ld'], panel['yc'], '-', color='#1f77b4', lw=1.2)   # blue = global
    ax.set_title(f"{panel['residue']}  R²(global)={panel['r2']:.2f}", fontsize=8)
    ax.set_xlabel('[ligand]', fontsize=7)
    ax.set_ylabel(panel['ylabel'], fontsize=7)
    ax.tick_params(labelsize=6)
    if ylim is not None:
        ax.set_ylim(*ylim)


def _run_export_kd(args):
    """Render CSP / intensity fit figures + a summary from a self-contained kd fit JSON."""
    import json
    import csv
    if args.dry_run:
        return _dry_run(args, [('json', args.json)], {'output_dir': args.out})
    if not os.path.isfile(args.json):
        _emit_error(args, FileNotFoundError(f"Fit JSON not found: {args.json}"))
        return 1
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType, so text stays editable in Illustrator
    import matplotlib.pyplot as plt
    import numpy as np
    _add_modules_path('dynamiXs_v2o0', 'dynamiXs_Kd')
    from kd_models import csp_model, intensity_decay, ref_point_values

    unknown_kind = [k for k in args.kind if k not in ('curves', 'ref-bars', 'kd-bars',
                                                      'global-fit')]
    if unknown_kind:
        print(f"Unknown --kind value(s): {','.join(unknown_kind)} "
              "(use curves, ref-bars, kd-bars and/or global-fit)", file=sys.stderr)
        return 1

    with open(args.json) as fh:
        data = json.load(fh)
    fits = data.get('fits', [])
    meta = data.get('metadata', {})
    P0 = meta.get('protein_conc')
    observables = args.observable or [o for o in ('csp', 'intensity')
                                      if any(f.get(o) for f in fits)]
    os.makedirs(args.out, exist_ok=True)
    tag = f'{args.prefix}_' if args.prefix else ''

    # Collect each residue's fit into per-observable panel lists (+ the summary rows).
    summary_rows = []
    panels = {obs: [] for obs in observables}
    for f in fits:
        residue = f.get('residue', 'peak')
        for obs in observables:
            fit = f.get(obs)
            if not fit or not fit.get('success'):
                continue
            summary_rows.append({
                'residue': residue, 'observable': obs,
                'Kd': fit.get('Kd'), 'Kd_err': fit.get('Kd_err'),
                'r_squared': fit.get('r_squared'),
            })
            L = np.asarray(fit['L'], dtype=float)
            y = np.asarray(fit['obs'], dtype=float)
            good = np.isfinite(L) & np.isfinite(y)
            if good.sum() < 2:
                continue
            Ld = np.linspace(float(L[good].min()), float(L[good].max()), 200)
            panels[obs].append(_build_kd_panel(residue, fit, L[good], y[good], Ld, obs, P0))

    summary_path = os.path.join(_data_dir(args.out), f'{tag}summary.csv')
    with open(summary_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['residue', 'observable', 'Kd', 'Kd_err', 'r_squared'])
        writer.writeheader()
        writer.writerows(summary_rows)

    unknown = [fmt for fmt in args.fig_format if fmt not in ('pdf', 'png')]
    if unknown:
        print(f"Unknown --fig-format value(s): {','.join(unknown)} (use pdf and/or png)",
              file=sys.stderr)
        return 1

    outputs = []
    if not args.summary_only and 'curves' in args.kind:
        cols = min(4, max(1, args.per_page))
        rows = max(1, -(-args.per_page // cols))  # ceil
        for obs in observables:
            plist = panels[obs]
            if not plist:
                continue
            ylim = _INTENSITY_YLIM if obs == 'intensity' else None
            if 'png' in args.fig_format:
                obs_dir = os.path.join(args.out, f'{tag}{obs}')
                os.makedirs(obs_dir, exist_ok=True)
                for p in plist:
                    fig, ax = plt.subplots(figsize=(4, 3))
                    _draw_kd_panel(ax, p, ylim=ylim)
                    fig.tight_layout()
                    fig.savefig(os.path.join(obs_dir, f"{_safe_name(p['residue'])}.png"), dpi=120)
                    plt.close(fig)
                outputs.append(obs_dir)
            if 'pdf' in args.fig_format:
                from matplotlib.backends.backend_pdf import PdfPages
                pdf_path = os.path.join(args.out, f"{tag}{obs}_fits.pdf")
                with PdfPages(pdf_path) as pdf:
                    for start in range(0, len(plist), args.per_page):
                        page = plist[start:start + args.per_page]
                        fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 3))
                        axes = np.atleast_1d(axes).flatten()
                        for ax, p in zip(axes, page):
                            _draw_kd_panel(ax, p, ylim=ylim)
                        for ax in axes[len(page):]:
                            ax.axis('off')
                        fig.suptitle(f"{obs.upper()} binding fits", fontsize=12)
                        fig.tight_layout(rect=[0, 0, 1, 0.98])
                        pdf.savefig(fig)
                        plt.close(fig)
                outputs.append(pdf_path)

    # Reference→point: the observable per residue between point 0 and each later
    # titration point. The data (wide CSV + JSON) is always written so the figures are
    # reproducible in Excel; the bars PDF is written only when requested.
    if not args.summary_only and 'ref-bars' in args.kind:
        labels = meta.get('concentrations') or meta.get('points') or []
        alpha = float(meta.get('alpha', 0.14))
        value = meta.get('intensity_value', 'height')
        # Ref→point is model-free (computed from the raw series), so default to BOTH
        # observables even when only one was fitted; --observable still restricts.
        ref_observables = args.observable or ['csp', 'intensity']
        if len(labels) >= 2 and any(f.get('series') for f in fits):
            from kd_export import export_ref_vs_point
            for obs in ref_observables:
                outputs.extend(export_ref_vs_point(
                    os.path.join(_data_dir(args.out), f"{tag}{obs}_ref_vs_point"),
                    fits, labels, obs, alpha=alpha, value=value))
            if 'pdf' in args.fig_format:
                from matplotlib.backends.backend_pdf import PdfPages
                for obs in ref_observables:
                    pdf_path = os.path.join(args.out, f"{tag}{obs}_ref_vs_point.pdf")
                    with PdfPages(pdf_path) as pdf:
                        for j in range(1, len(labels)):
                            names, vals = ref_point_values(fits, 0, j, obs,
                                                           alpha=alpha, value=value)
                            fig, ax = plt.subplots(figsize=(11, 6))
                            sig = (meta.get('csp_significance') or {}) if obs == 'csp' else {}
                            thr = sig.get('threshold')
                            thr_lab = None
                            if thr is not None:
                                last = labels[-1]
                                last = f"{last:g}" if isinstance(last, (int, float)) else last
                                thr_lab = (f"{sig.get('multiple', 1.0):g}σ significance "
                                           f"threshold = {thr:.4f}  (defined at {last})")
                            _draw_ref_bars(ax, names, vals, obs, labels[0], labels[j],
                                           threshold=thr, threshold_label=thr_lab)
                            fig.tight_layout()
                            pdf.savefig(fig)
                            plt.close(fig)
                    outputs.append(pdf_path)

    # Kd-per-residue bars (+ the shared global-Kd line) for each observable. PDF only.
    if not args.summary_only and 'kd-bars' in args.kind and 'pdf' in args.fig_format:
        from matplotlib.backends.backend_pdf import PdfPages
        global_fit = data.get('global', {}) or {}
        for obs in observables:
            # Only residues whose Kd is actually reported — i.e. those in the shared
            # fit's pool. Drawing every residue put insignificant and unmeasured ones
            # beside the ones that count, with nothing saying which was which.
            pool = (global_fit.get(obs, {}) or {}).get('residues')
            names, kds, errs = [], [], []
            for f in fits:
                fit = f.get(obs) or {}
                ok = fit.get('success') and isinstance(fit.get('Kd'), (int, float))
                if pool is not None and f.get('residue') not in pool:
                    continue
                names.append(f.get('residue', 'peak'))
                kds.append(float(fit['Kd']) if ok else float('nan'))
                e = fit.get('Kd_err')
                errs.append(float(e) if ok and isinstance(e, (int, float)) else 0.0)
            if not any(np.isfinite(v) for v in kds):
                continue
            gkd = (global_fit.get(obs, {}) or {}).get('Kd')
            pdf_path = os.path.join(args.out, f"{tag}{obs}_kd_vs_residue.pdf")
            with PdfPages(pdf_path) as pdf:
                fig, ax = plt.subplots(figsize=(11, 6))
                _draw_kd_bars(ax, names, kds, errs, obs, gkd)
                excl = data.get('metadata', {}).get('csp_pool_excluded') or {}
                if obs == 'csp' and excl:
                    ax.set_title(f"{ax.get_title()}  —  {len(names)} of "
                                 f"{len(names) + len(excl)} residues "
                                 f"(see csp_pool_excluded for the rest)", fontsize=9)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
            outputs.append(pdf_path)

    # Global shared-Kd per-residue params (amplitudes + shared Kd + per-residue R^2) as
    # Excel-ready CSV + JSON, so the global-fit figure is reproducible as data.
    if not args.summary_only and 'global-fit' in args.kind:
        from kd_export import export_global_fit
        global_fit = data.get('global', {}) or {}
        for obs in observables:
            outputs.extend(export_global_fit(
                os.path.join(_data_dir(args.out), f"{tag}{obs}_global_fit"),
                fits, global_fit, obs, P0))

    # Global shared-Kd fit over the data: per-residue observed points + the single-Kd
    # model curve, so one can judge how well one Kd fits every residue. PDF only.
    if not args.summary_only and 'global-fit' in args.kind and 'pdf' in args.fig_format:
        from matplotlib.backends.backend_pdf import PdfPages
        global_fit = data.get('global', {}) or {}
        cols = min(4, max(1, args.per_page))
        rows = max(1, -(-args.per_page // cols))
        for obs in observables:
            g = global_fit.get(obs) or {}
            amp = g.get('dd_max') if obs == 'csp' else g.get('I0')
            if not g.get('success') or not amp:
                continue
            kd = g.get('Kd')
            gpanels = []
            for f in fits:
                res = f.get('residue', 'peak')
                fit = f.get(obs) or {}
                if res not in amp or not fit.get('L'):
                    continue
                L = np.asarray(fit['L'], dtype=float)
                y = np.asarray(fit['obs'], dtype=float)
                good = np.isfinite(L) & np.isfinite(y)
                if good.sum() < 2:
                    continue
                Lg = np.linspace(0.0, float(L[good].max()) * 1.05, 200)
                y_plot = y[good]
                if obs == 'csp':
                    yc = csp_model(Lg, g['dd_max'][res], kd, P0)
                    yhat = csp_model(L[good], g['dd_max'][res], kd, P0)
                    ylabel = 'CSP (ppm)'
                else:
                    yc = intensity_decay(Lg, g['I0'][res], g['I_inf'][res], kd)
                    yhat = intensity_decay(L[good], g['I0'][res], g['I_inf'][res], kd)
                    ylabel = 'I / I(0)'
                ss_res = float(np.sum((y[good] - yhat) ** 2))
                ss_tot = float(np.sum((y[good] - np.mean(y[good])) ** 2))
                r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')
                if obs == 'intensity':
                    # Normalize to I/I(0) so this panel shares the same 0-1 axis as
                    # every other intensity panel, regardless of raw data scale.
                    y_plot = y_plot / g['I0'][res]
                    yc = yc / g['I0'][res]
                gpanels.append({'residue': res, 'L': L[good], 'y': y_plot, 'Ld': Lg,
                                'yc': yc, 'ylabel': ylabel, 'r2': r2})
            if not gpanels:
                continue
            ylim = _INTENSITY_YLIM if obs == 'intensity' else None
            pdf_path = os.path.join(args.out, f"{tag}{obs}_global_fit.pdf")
            with PdfPages(pdf_path) as pdf:
                for start in range(0, len(gpanels), args.per_page):
                    page = gpanels[start:start + args.per_page]
                    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 3))
                    axes = np.atleast_1d(axes).flatten()
                    for ax, p in zip(axes, page):
                        _draw_global_panel(ax, p, ylim=ylim)
                    for ax in axes[len(page):]:
                        ax.axis('off')
                    fig.suptitle(f"Global shared-Kd fit ({obs}): one Kd={kd:.4g} for all residues",
                                 fontsize=12)
                    fig.tight_layout(rect=[0, 0, 1, 0.98])
                    pdf.savefig(fig)
                    plt.close(fig)
            outputs.append(pdf_path)

    _emit(args,
          {'command': 'export kd', 'n_fits': len(summary_rows),
           'summary_csv': summary_path, 'figures': outputs},
          f"Exported {len(summary_rows)} fit(s)",
          f"  Summary: {summary_path}",
          *[f"  Figures: {o}" for o in outputs])
    return 0


def _run_dynamixs_hetnoe(args):
    """Compute heteronuclear NOE (I_sat / I_unsat) per residue, as the GUI pipeline does."""
    import csv
    if args.dry_run:
        return _dry_run(args, [('sat', args.sat), ('unsat', args.unsat)], {'output_dir': args.out})
    _add_modules_path('dynamiXs_v2o0')
    from dynamiXs_integrated.data_converters import (parse_intensity_csv,
                                                     calculate_hetnoe_from_intensities,
                                                     plot_hetnoe_vs_residue)
    for label, path in (('--sat', args.sat), ('--unsat', args.unsat)):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"{label} file not found: {path}")
    with _engine_stdout(args):
        sat_i, sat_e = parse_intensity_csv(args.sat)
        unsat_i, unsat_e = parse_intensity_csv(args.unsat)
        noe = calculate_hetnoe_from_intensities(saturated_data=sat_i, unsaturated_data=unsat_i,
                                                saturated_errors=sat_e, unsaturated_errors=unsat_e)
    os.makedirs(args.out, exist_ok=True)
    out_csv = os.path.join(args.out, f"{args.prefix}_hetnoe.csv")
    with open(out_csv, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['Residue', 'hetNOE', 'hetNOEerr'])
        for res, rec in noe.items():
            writer.writerow([res, rec['value'], rec['error']])
    out_pdf = os.path.join(args.out, f"{args.prefix}_hetnoe.pdf")
    with _engine_stdout(args):
        plot_hetnoe_vs_residue(noe, out_pdf, title=args.prefix)
    # Which error columns were actually supplied decides whether the reported hetNOE
    # error is measured or estimated, and the two differ by ~3x. A caller weighting the
    # model-free fit by these needs to know which it got.
    supplied = {(True, True): 'both', (True, False): 'sat',
                (False, True): 'unsat', (False, False): 'none'}[
                    (sat_e is not None, unsat_e is not None)]
    _emit(args,
          {'command': 'dynamixs hetnoe', 'n_residues': len(noe),
           'errors_supplied': supplied,
           'n_sat': len(sat_i), 'n_unsat': len(unsat_i),
           'n_common': len(set(sat_i) & set(unsat_i)),
           'n_dropped_nonpositive_ref': len(set(sat_i) & set(unsat_i)) - len(noe),
           'output': out_csv, 'plot': out_pdf},
          f"hetNOE complete: {len(noe)} residues "
          f"({len(sat_i)} sat, {len(unsat_i)} unsat, errors: {supplied})",
          *(["  NOTE: no error columns supplied — hetNOE errors are a flat 2% estimate, "
             "~3x tighter than a realistic floor, which over-weights hetNOE downstream"]
            if supplied == 'none' else []),
          f"  Output: {out_csv}",
          f"  Plot:   {out_pdf}")
    return 0


def _run_dynamixs_t1rho(args):
    """Fit a T1rho series and convert it to R2, so T1/T1rho/hetNOE reaches model-free.

    R1rho = R1 cos^2(theta) + R2 sin^2(theta), so R2 needs R1 as well as the spin-lock
    geometry. The tilt angle is computed per residue: residues sit at different offsets
    from the spin-lock carrier, so one nominal angle is not enough.
    """
    inputs = [('input', args.input), ('t1', args.t1), ('peaks', args.peaks)]
    if args.dry_run:
        return _dry_run(args, inputs, {'output_dir': args.out})
    for label, path in inputs:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"--{label} file not found: {path}")
    import matplotlib
    matplotlib.use('Agg')
    import json
    _add_modules_path('dynamiXs_v2o0', 'dynamiXs_T1_T2')
    from fit_Tx_NMRRE import run_analysis_with_params
    from lunaNMR.utils.file_manager import NMRFileManager
    from lunaNMR.utils.t1rho_calculator import r2_table_from_fits
    import pandas as pd

    # Both series are fitted here, so both sidecars matter; they must agree with each
    # other as well as with the flag, or the tilt correction mixes two time bases.
    sidecars = {label: _series_sidecar(path)
                for label, path in (('--input', args.input), ('--t1', args.t1))}
    for label, sidecar in sidecars.items():
        if sidecar:
            _check_unparsed_delays(sidecar, args.allow_unparsed)
    time_units = args.time_units
    for label, sidecar in sidecars.items():
        time_units = _resolve_time_units(time_units, sidecar, _T1RHO_TIME_UNITS_DEFAULT,
                                         flag_name=f'--time-units (against {label})')
    os.makedirs(args.out, exist_ok=True)

    def _fit(input_csv, label):
        """Both series are mono-exponential decays, so the shared T1/T2 fitter is used."""
        params = {
            'input_csv_file': input_csv,
            'output_prefix': os.path.join(args.out, f"{args.prefix}_{label}"),
            'results_txt_file': os.path.join(args.out, f"{args.prefix}_{label}_fit_results.txt"),
            'experiment_type': 'T2',
            'time_units': time_units,
            'error_method': args.error_method,
            'n_bootstrap': args.bootstrap,
            'field_name': label,
            'field_freq': args.field_freq,
            'json_folder': args.out,
        }
        with _engine_stdout(args):
            result = run_analysis_with_params(params)
        data = json.load(open(result['json_file']))
        # Keep only the residues the fitter itself judged reliable. It prints
        # "unreliable (no decay in window), excluded" for the rest and leaves them out
        # of n_fitted, but they were still in the JSON -- and the only downstream filter
        # is finite-and-positive, which a T2 of 2.7e9 ms passes.
        fits = {str(f['residue']): {'value': f['t2'], 'error': f.get('t2_err', 0.0)}
                for f in data['fits'] if f.get('success', True)}
        return fits, result

    t1_fits, t1_res = _fit(os.path.abspath(args.t1), 'T1')
    rho_fits, rho_res = _fit(os.path.abspath(args.input), 'T1rho')

    peaks_df = NMRFileManager().load_peak_list(os.path.abspath(args.peaks))
    peak_list = pd.DataFrame({'residue': peaks_df['Assignment'].astype(str),
                              'N_ppm': peaks_df['Position_Y'].astype(float)})

    table = r2_table_from_fits(t1_fits, rho_fits, peak_list,
                               omega1_hz=args.omega1, carrier_ppm=args.carrier,
                               theta_deg=args.theta, spec_freq_mhz=args.field_freq,
                               time_units=time_units)
    out_csv = os.path.join(args.out, f"{args.prefix}_r2_from_t1rho.csv")
    table.to_csv(out_csv, index=False)
    n_unmatched = int(table.attrs.get('n_shift_unmatched', 0))
    if n_unmatched:
        print(f"  WARNING: {n_unmatched}/{len(table)} residues had no 15N shift in "
              f"--peaks, so their tilt angle fell back to the nominal "
              f"{args.theta} deg (no per-residue correction): "
              f"{', '.join(table.attrs.get('shift_unmatched', [])[:10])}",
              file=sys.stderr)
    _emit(args,
          {'command': 'dynamixs t1rho', 'n_residues': len(table),
           'n_shift_unmatched': n_unmatched,
           'shift_unmatched': table.attrs.get('shift_unmatched', []),
           'time_units': time_units,
           'series_metadata': {k: _sidecar_summary(v) for k, v in sidecars.items()},
           'n_t1_fitted': len(t1_fits), 'n_t1rho_fitted': len(rho_fits),
           'theta_nominal_deg': args.theta, 'omega1_hz': args.omega1,
           'r2_table': out_csv,
           't1_results': t1_res['results_file'], 't1rho_results': rho_res['results_file']},
          f"T1rho complete: R2 derived for {len(table)} residues",
          f"  theta (nominal): {args.theta} deg, omega1: {args.omega1} Hz "
          f"(per-residue tilt applied)",
          f"  R2 table: {out_csv}")
    return 0 if len(table) else 1


def _run_dynamixs_density(args):
    """Reduced spectral density mapping (Farrow 1995) from an R1/R2/hetNOE table."""
    if args.dual and not args.input2:
        raise ValueError("dual-field density needs a second table (--input2)")
    if args.dual and args.field2_freq is None:
        raise ValueError("dual-field density needs --field2-freq (the second "
                         "spectrometer 1H frequency in MHz); there is no safe default")
    inputs = [('input', args.input)] + ([('input2', args.input2)] if args.dual else [])
    if args.dry_run:
        return _dry_run(args, inputs, {'output_dir': args.out})
    import matplotlib
    matplotlib.use('Agg')
    from importlib import import_module
    _add_modules_path('dynamiXs_v2o0')
    if not os.path.isfile(args.input):
        raise FileNotFoundError(f"Input not found: {args.input}")
    os.makedirs(args.out, exist_ok=True)
    r_nh_m = args.rnh * 1e-10       # Angstrom -> metre
    csa_n = args.csa * 1e-6         # ppm -> dimensionless
    prefix = os.path.join(args.out, args.prefix)
    suffix = '_087' if args.use_087 else ''

    with _engine_stdout(args):
        if args.dual:
            cls = getattr(import_module(f'dynamiXs_density_functions.ZZ_multi_2fields_density{suffix}'),
                          'DualFieldSpectralDensityAnalysis')
            analyzer = cls(field1_freq=args.field1_freq, field2_freq=args.field2_freq,
                           rNH=r_nh_m, csaN=csa_n)
            df = analyzer.analyze_dual_field_csv(
                csv_file1=args.input, csv_file2=args.input2,
                use_monte_carlo_errors=args.monte_carlo, n_monte_carlo=args.n_samples,
                use_multiprocessing=not args.no_parallel)
            basic, detailed = f"{prefix}_basic.csv", f"{prefix}_detailed.csv"
            df.to_csv(basic, index=False)
            analyzer.save_dual_field_results(df, detailed)
            outputs = [basic, detailed]
            if not args.no_plot:
                plots = f"{prefix}_plots.pdf"
                analyzer.plot_dual_field_results(df, save_plots=True, plot_filename=plots)
                outputs.append(plots)
        else:
            cls = getattr(import_module(f'dynamiXs_density_functions.ZZ_multi_density{suffix}'),
                          'ReducedSpectralDensityAnalysis')
            analyzer = cls(spectrometer_frequency=args.field1_freq, rNH=r_nh_m, csaN=csa_n)
            df = analyzer.analyze_csv(
                csv_file=args.input, use_monte_carlo_errors=args.monte_carlo,
                n_monte_carlo=args.n_samples, use_multiprocessing=not args.no_parallel)
            results = f"{prefix}_results.csv"
            df.to_csv(results, index=False)
            outputs = [results]
            if not args.no_plot:
                # The single-field plotter writes a hard-coded filename to the CWD, so run
                # it from the output directory and record whatever PDF it produced.
                cwd = os.getcwd()
                try:
                    os.chdir(args.out)
                    before = set(f for f in os.listdir('.') if f.endswith('.pdf'))
                    analyzer.plot_results(df, save_plots=True)
                    new = [f for f in os.listdir('.') if f.endswith('.pdf') and f not in before]
                finally:
                    os.chdir(cwd)
                outputs += [os.path.join(args.out, f) for f in new]

    _emit(args,
          {'command': 'dynamixs density', 'n_residues': len(df),
           'dual': args.dual, 'outputs': outputs},
          f"Spectral density complete: {len(df)} residues",
          *[f"  {o}" for o in outputs])
    return 0


def _validate_modelfree(parser, args):
    """R2 comes from a T2 series or from a precomputed table, never both."""
    for field in ('f1', 'f2'):
        if getattr(args, f'{field}_t2') and getattr(args, f'{field}_r2_table'):
            parser.error(f'--{field}-t2 and --{field}-r2-table are alternatives; give one')
    if not (args.f1_t2 or args.f1_r2_table):
        parser.error('one of --f1-t2 or --f1-r2-table is required')
    # There is no conventional second field, so a default is a guess -- and analysing
    # an 800 MHz dataset at 700 surfaces as a tau_c mismatch that the QC table blames
    # on the data rather than on the flag.
    if args.dual and args.field2_freq is None:
        parser.error('--dual needs --field2-freq (the second spectrometer 1H frequency '
                     'in MHz); there is no safe default for it')


def _run_dynamixs_modelfree(args):
    """Integrated model-free pipeline: T1/T2 fit -> hetNOE -> density -> Lipari-Szabo."""
    f1 = [('f1-t1', args.f1_t1), ('f1-t2', args.f1_t2 or args.f1_r2_table),
          ('f1-noe-sat', args.f1_noe_sat), ('f1-noe-unsat', args.f1_noe_unsat)]
    f2 = [('f2-t1', args.f2_t1), ('f2-t2', args.f2_t2 or args.f2_r2_table),
          ('f2-noe-sat', args.f2_noe_sat), ('f2-noe-unsat', args.f2_noe_unsat)]
    if args.dual and not all(v for _, v in f2):
        missing = ', '.join(f'--{label}' for label, value in f2 if not value)
        raise ValueError(f"dual-field model-free needs all four --f2-* files; "
                         f"missing: {missing}")
    if args.dry_run:
        return _dry_run(args, f1 + (f2 if args.dual else []), {'output_dir': args.out})
    import matplotlib
    matplotlib.use('Agg')
    _add_modules_path('dynamiXs_v2o0')
    _add_modules_path('dynamiXs_v2o0', 'dynamiXs_T1_T2')
    from dynamiXs_integrated.integrated_analysis import (IntegratedAnalysisPipeline,
                                                         IntegratedAnalysisParameters)
    # Each relaxation series carries its own units. `series` already recorded them, so
    # the flag is only needed for hand-built tables -- and a flag that contradicts the
    # sidecar is refused rather than believed: a T1 read as seconds against a T2 in
    # milliseconds puts R1 out by 1000x, and nothing downstream notices.
    units, sidecars = {}, {}
    for field, exp, path in (('f1', 't1', args.f1_t1), ('f1', 't2', args.f1_t2),
                             ('f2', 't1', args.f2_t1), ('f2', 't2', args.f2_t2)):
        key = f'{field}_{exp}_units'
        if not path:
            units[key] = getattr(args, key) or _MODELFREE_UNITS_DEFAULT
            continue
        sidecar = _series_sidecar(path)
        if sidecar:
            _check_unparsed_delays(sidecar, args.allow_unparsed)
            sidecars[f'{field}-{exp}'] = _sidecar_summary(sidecar)
        units[key] = _resolve_time_units(getattr(args, key), sidecar,
                                         _MODELFREE_UNITS_DEFAULT,
                                         flag_name=f'--{field}-{exp}-units')
    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)
    p = IntegratedAnalysisParameters()
    p.field1_t1_file = os.path.abspath(args.f1_t1)
    p.field1_t2_file = os.path.abspath(args.f1_t2) if args.f1_t2 else None
    p.field1_r2_table_file = os.path.abspath(args.f1_r2_table) if args.f1_r2_table else None
    p.field1_noe_sat_file = os.path.abspath(args.f1_noe_sat)
    p.field1_noe_unsat_file = os.path.abspath(args.f1_noe_unsat)
    p.field1_freq_mhz = args.field1_freq
    p.field1_t1_units = units['f1_t1_units']
    p.field1_t2_units = units['f1_t2_units']
    p.enable_dual_field = args.dual
    if args.dual:
        p.field2_t1_file = os.path.abspath(args.f2_t1)
        p.field2_t2_file = os.path.abspath(args.f2_t2) if args.f2_t2 else None
        p.field2_r2_table_file = os.path.abspath(args.f2_r2_table) if args.f2_r2_table else None
        p.field2_noe_sat_file = os.path.abspath(args.f2_noe_sat)
        p.field2_noe_unsat_file = os.path.abspath(args.f2_noe_unsat)
        p.field2_freq_mhz = args.field2_freq
        p.field2_t1_units = units['f2_t1_units']
        p.field2_t2_units = units['f2_t2_units']
    # Field count must match --dual (the density step branches on the method prefix);
    # keep only the user's 087/jwh variant preference.
    variant = 'jwh' if (args.method or '').endswith('jwh') else '087'
    p.analysis_method = ('dual_' if args.dual else 'single_') + variant
    p.rNH_angstrom = args.rnh
    p.csaN_ppm = args.csa
    p.t1_initial_amplitude = args.init_amp
    p.t2_initial_amplitude = args.init_amp
    p.t1_initial_time = args.init_t1
    p.t2_initial_time = args.init_t2
    p.t1_bootstrap_iterations = args.n_bootstrap
    p.t2_bootstrap_iterations = args.n_bootstrap
    p.error_method = args.error_method
    p.monte_carlo_iterations = args.n_monte_carlo
    p.output_dir = out_dir
    p.output_prefix = args.prefix

    # Under --format json the engine's own prints already go to stderr via
    # _engine_stdout, but this callback bypasses it. Nulling it discarded
    # validate_relaxation_rates -- the only check that catches a units error -- in
    # exactly the mode the runbooks tell an agent to use.
    cb = ((lambda msg, *a, **k: print(msg, file=sys.stderr))
          if getattr(args, 'format', 'text') == 'json' else None)
    # Run from the output dir: the single-field density plotter and the intermediate
    # fit-result txt files are written to the CWD with relative/hard-coded names.
    cwd = os.getcwd()
    try:
        os.chdir(out_dir)
        with _engine_stdout(args):
            result = IntegratedAnalysisPipeline(p, progress_callback=cb).run_complete_analysis()
    finally:
        os.chdir(cwd)
    out_files = result.get('output_files', {})
    _emit(args,
          {'command': 'dynamixs modelfree', 'method': result.get('method'),
           'time_units': units, 'series_metadata': sidecars or None,
           # Counts per category, so a caller can threshold on them. The residue lists
           # are in the human report; a summary that carries neither is the black hole.
           'validation': result.get('validation_warnings') or {},
           'residue_counts': result.get('residue_counts') or {},
           'n_residues': result.get('n_residues'), 'n_successful': result.get('n_successful'),
           'output_files': out_files},
          f"Model-free complete: {result.get('n_successful')}/{result.get('n_residues')} residues",
          *[f"  {k}: {v}" for k, v in out_files.items()])
    return 0


def _run_batch(batch_argv):
    """Delegate to the existing batch CLI, passing through all of its flags."""
    from lunaNMR.batch_processing.cli_interface import CLIInterface
    return CLIInterface().main(batch_argv)


def build_parser():
    parser = argparse.ArgumentParser(
        prog='lunaNMR',
        description='LunaNMR command-line interface for headless NMR analysis.',
        # `--help` is the one entry point every caller already runs, so it is where the
        # documentation an agent needs before running anything has to be named.
        epilog=(
            'Documentation:\n'
            '  docs/CLI.md                human reference for every subcommand and flag\n'
            '  docs/CLI_AGENT.md          machine contract, output shapes, and the\n'
            '                             silent-corruption gotchas -- read this before\n'
            '                             driving the CLI programmatically\n'
            '  docs/CLI_AGENTS_DEEP/      long-form runbooks: phase structure, physical\n'
            '                             QC bands, and worked relaxation/affinity flows\n'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--version', action='version', version=f'lunaNMR {_version()}')
    sub = parser.add_subparsers(dest='command', metavar='<subcommand>')

    # Shared flags for the analysis subcommands: output format + input validation.
    # Split because `diagnose` reads spectra, fits nothing and writes nothing: a
    # --dry-run of a read-only command has nothing to mean, and it was accepted and
    # silently ignored.
    fmt = argparse.ArgumentParser(add_help=False)
    fmt.add_argument('--format', choices=['text', 'json'], default='text',
                     help='Run-summary output format (default: text)')
    common = argparse.ArgumentParser(add_help=False, parents=[fmt])
    common.add_argument('--dry-run', action='store_true', dest='dry_run',
                        help='Validate inputs and print the plan without running')

    kd = sub.add_parser('kd', parents=[common],
                        help='Fit Kd from a titration series (CSP / intensity)')
    kd.add_argument('--input', required=True,
                    help='Titration CSV (tidy / comprehensive_peak_tracking / intensity matrix)')
    kd.add_argument('--out', required=True, help='Output directory for the fit JSON + results')
    kd.add_argument('--prefix', default='kd', help='Output filename prefix (default: kd)')
    kd.add_argument('--p0', type=float, required=True, help='Total protein concentration [P]0')
    kd.add_argument('--conc', type=_float_list, default=None,
                    help='Comma-separated ligand concentrations (default: CSV point labels)')
    kd.add_argument('--alpha', type=float, default=0.14,
                    help='CSP N/H scaling factor (default: 0.14)')
    kd.add_argument('--observable', type=_choice_list('csp', 'intensity'),
                    default=['csp', 'intensity'],
                    help='Comma-separated observables to fit: csp,intensity (default: both)')
    kd.add_argument('--intensity-from', choices=['height', 'volume'], default='height',
                    dest='intensity_from', help='Intensity source for the ratio (default: height)')
    kd.add_argument('--bootstrap', type=int, default=0,
                    help='Bootstrap iterations for error estimates (default: 0)')
    kd.add_argument('--intensity-scale', type=_float_list, default=None,
                    help='Comma-separated per-point height/volume scale factors')
    kd.add_argument('--survey', action='store_true',
                    help='Survey residues and write an editable selection file + '
                         'vs-sequence figures. Writes no fit JSON and no global Kd')
    kd.add_argument('--conc-units', choices=['absolute', 'equivalents'], default=None,
                    dest='conc_units',
                    help='Whether concentrations (parsed from spectrum names or given via '
                         '--conc) are absolute or equivalents of --p0. Equivalents are '
                         'multiplied by --p0 (default: absolute)')
    kd.add_argument('--csp-sigma-multiple', type=float, default=None,
                    dest='csp_sigma_multiple',
                    help='Multiples of the trimmed CSP spread at the last titration '
                         'point a residue must exceed to enter the shared CSP fit '
                         '(default 1.0)')
    kd.add_argument('--kd-outlier-z', type=float, default=None, dest='kd_outlier_z',
                    help='Robust z (median/MAD on log10 Kd) beyond which a residue '
                         'leaves the shared CSP fit (default 3.0; 0 disables)')
    kd.add_argument('--selection', default=None,
                    help='Path for the survey selection file (default: beside --input)')
    kd.add_argument('--noise-quantile', type=float, default=None, dest='noise_quantile',
                    help='Quantile of the max-CSP distribution below which a residue '
                         'counts as a non-mover (default 0.25; persisted in the params JSON)')
    kd.add_argument('--dd-runaway-ratio', type=float, default=None, dest='dd_runaway_ratio',
                    help='Flag a residue whose fitted plateau exceeds its own largest CSP '
                         'by more than this (default 10.0; 17%% of one dataset, 43%% of another)')
    kd.add_argument('--ref-max-ratio', type=float, default=None, dest='ref_max_ratio',
                    help='Reject a reference intensity this many times below the residue '
                         "own series max (default 10.0; legitimate values top out at 1.30)")
    kd.add_argument('--residues', default=None,
                    help='Restrict the fit to these residues: a selection file from '
                         '--survey (one name per line, # comments out) or a '
                         'comma-separated list')
    kd.set_defaults(func=_run_kd)

    series = sub.add_parser('series', parents=[common],
                            help='Process a multi-spectrum series/titration')
    series.add_argument('--spectra', required=True, help='Folder or glob of spectrum files')
    series.add_argument('--peaks', required=True, help='Reference peak-list file (Assignment, Position_X, Position_Y)')
    series.add_argument('--out', required=True, help='Output folder for the series_results CSVs')
    series.add_argument('--mode', choices=['time', 'titration'], default='time',
                        help='Series type: time (relaxation) or titration (default: time)')
    series.add_argument('--peak-source', choices=['reference', 'cascade', 'detected', 'independent'],
                        default='reference', dest='peak_source',
                        help='Peak position source across spectra (default: reference)')
    series.add_argument('--deep', action='store_true',
                        help='With --dry-run: also check registration, capture, delays and '
                             'the peak list against the spectra (reads them; seconds, not instant)')
    series.add_argument('--quick', action='store_true',
                        help='Coarser registration grid for --deep')
    series.add_argument('--sample', action='store_true',
                        help='With --deep: assess only the first and last spectrum')
    series.add_argument('--parallel', action='store_true',
                        help='Use the two-pass parallel processor (~2.7x faster)')
    series.set_defaults(func=_run_series, _validate=_validate_series)

    diagnose = sub.add_parser('diagnose', parents=[fmt],
                              help='Read-only pre-flight over a dataset: registration, capture, '
                                   'delays, peak lists, and the cross-experiment residue set')
    diagnose.add_argument('root', help='Dataset root containing the experiment folders')
    diagnose.add_argument('--quick', action='store_true', help='Coarser registration grid')
    diagnose.add_argument('--sample', action='store_true',
                          help='Assess only the first and last spectrum of each experiment')
    diagnose.add_argument('--mode', choices=['time', 'titration'], default='time',
                          help='How filenames encode the series value (default: time)')
    diagnose.add_argument('--strict', action='store_true',
                          help='Also exit 1 on WARN findings (default: FAIL only)')
    diagnose.set_defaults(func=_run_diagnose)

    dx = sub.add_parser('dynamixs', help='Relaxation fitting: T1/T2 and methyl-T2')
    dx_sub = dx.add_subparsers(dest='dynamixs_command', metavar='<kind>')

    t1t2 = dx_sub.add_parser('t1t2', parents=[common],
                             help='Mono-exponential T1 or T2 relaxation fit')
    _add_relaxation_flags(t1t2)
    t1t2.add_argument('--exp', choices=['T1', 'T2'], required=True, help='Experiment type')
    t1t2.add_argument('--no-json', action='store_true', help='Skip writing the JSON fit data')
    t1t2.set_defaults(func=_run_dynamixs_t1t2)

    methyl = dx_sub.add_parser('methyl-t2', parents=[common],
                               help='Bi-exponential (Tugarinov-Kay) methyl T2 fit')
    _add_relaxation_flags(methyl)
    methyl.set_defaults(func=_run_dynamixs_methyl)

    hetnoe = dx_sub.add_parser('hetnoe', parents=[common],
                               help='Heteronuclear NOE = I_sat / I_unsat per residue')
    hetnoe.add_argument('--sat', required=True, help='Saturated intensities (residue,intensity[,err])')
    hetnoe.add_argument('--unsat', required=True, help='Unsaturated intensities (residue,intensity[,err])')
    hetnoe.add_argument('--out', required=True, help='Output directory')
    hetnoe.add_argument('--prefix', default='field1', help='Output filename prefix (default: field1)')
    hetnoe.set_defaults(func=_run_dynamixs_hetnoe)

    t1rho = dx_sub.add_parser('t1rho', parents=[common],
                              help='Fit a T1rho series and convert it to R2 (needs T1)')
    t1rho.add_argument('--input', required=True, help='T1rho relaxation matrix')
    t1rho.add_argument('--t1', required=True, help='T1 relaxation matrix (R1 is needed for the tilt correction)')
    t1rho.add_argument('--peaks', required=True, help='Peak list, for each residue\'s 15N shift')
    t1rho.add_argument('--omega1', type=float, required=True, help='Spin-lock field strength in Hz (cnst27)')
    t1rho.add_argument('--carrier', type=float, required=True, help='15N carrier position in ppm')
    t1rho.add_argument('--theta', type=float, required=True, help='Nominal tilt angle in degrees (cnst28)')
    t1rho.add_argument('--field-freq', type=float, default=600.0, dest='field_freq',
                       help='1H spectrometer frequency in MHz (default: 600)')
    t1rho.add_argument('--out', required=True, help='Output directory')
    t1rho.add_argument('--prefix', default='field1', help='Output filename prefix')
    t1rho.add_argument('--time-units', choices=['ms', 's', 'us'], default=None, dest='time_units',
                       help='Units of the delay values (RESCALES the fitted rates). '
                            'Default: taken from series_metadata.json beside the input, else ms')
    t1rho.add_argument('--allow-unparsed-delays', action='store_true', dest='allow_unparsed',
                       help='Fit anyway when the series sidecar records spectra whose '
                            'filename carried no delay (those points are dropped)')
    t1rho.add_argument('--error-method', choices=['analytical', 'bootstrap'], default='analytical',
                       dest='error_method')
    t1rho.add_argument('--bootstrap', type=int, default=1000)
    t1rho.set_defaults(func=_run_dynamixs_t1rho)

    density = dx_sub.add_parser('density', parents=[common],
                                help='Reduced spectral density mapping from an R1/R2/hetNOE table')
    density.add_argument('--input', required=True,
                         help='Field-1 table: Residue,R1,R1err,R2,R2err,hetNOE,hetNOEerr')
    density.add_argument('--input2', help='Field-2 table (dual-field)')
    density.add_argument('--out', required=True, help='Output directory')
    density.add_argument('--prefix', default='spectral_density', help='Output prefix')
    density.add_argument('--dual', action='store_true', help='Dual-field analysis (needs --input2)')
    density.add_argument('--no-087', action='store_false', dest='use_087',
                         help='Evaluate J at full omega_H instead of 0.87*omega_H (default: 0.87)')
    density.add_argument('--monte-carlo', action='store_true', dest='monte_carlo',
                         help='Monte-Carlo error propagation (slower)')
    density.add_argument('--n-samples', type=int, default=1000, dest='n_samples',
                         help='Monte-Carlo samples (default: 1000)')
    density.add_argument('--no-plot', action='store_true', dest='no_plot', help='Skip the PDF plots')
    density.add_argument('--no-parallel', action='store_true', dest='no_parallel',
                         help='Disable multiprocessing')
    _add_density_flags(density)
    density.set_defaults(func=_run_dynamixs_density)

    mf = dx_sub.add_parser('modelfree', parents=[common],
                           help='Integrated model-free: T1/T2 fit -> hetNOE -> density -> Lipari-Szabo')
    mf.add_argument('--f1-t1', required=True, dest='f1_t1', help='Field-1 T1 relaxation matrix')
    mf.add_argument('--f1-t2', dest='f1_t2', help='Field-1 T2 relaxation matrix')
    mf.add_argument('--f1-r2-table', dest='f1_r2_table',
                    help='Field-1 per-residue R2 table (from `dynamixs t1rho`), instead of --f1-t2')
    mf.add_argument('--f1-noe-sat', required=True, dest='f1_noe_sat', help='Field-1 NOE saturated intensities')
    mf.add_argument('--f1-noe-unsat', required=True, dest='f1_noe_unsat', help='Field-1 NOE unsaturated intensities')
    mf.add_argument('--f2-t1', dest='f2_t1', help='Field-2 T1 matrix (dual-field)')
    mf.add_argument('--f2-t2', dest='f2_t2', help='Field-2 T2 matrix (dual-field)')
    mf.add_argument('--f2-r2-table', dest='f2_r2_table',
                    help='Field-2 per-residue R2 table, instead of --f2-t2')
    mf.add_argument('--f2-noe-sat', dest='f2_noe_sat', help='Field-2 NOE saturated (dual-field)')
    mf.add_argument('--f2-noe-unsat', dest='f2_noe_unsat', help='Field-2 NOE unsaturated (dual-field)')
    for f in ('f1', 'f2'):
        for exp in ('t1', 't2'):
            mf.add_argument(f'--{f}-{exp}-units', choices=['ms', 's', 'us'], default=None,
                            dest=f'{f}_{exp}_units',
                            help=f'Delay units of the {f.upper()} {exp.upper()} series '
                                 f'(RESCALES the rates). Default: from that series\''
                                 f' series_metadata.json, else ms')
    mf.add_argument('--allow-unparsed-delays', action='store_true', dest='allow_unparsed',
                    help='Fit anyway when a series sidecar records spectra whose '
                         'filename carried no delay (those points are dropped)')
    mf.add_argument('--out', required=True, help='Output directory')
    mf.add_argument('--prefix', default='modelfree', help='Output prefix')
    mf.add_argument('--dual', action='store_true', help='Dual-field (needs the f2 files)')
    mf.add_argument('--method', choices=['single_jwh', 'single_087', 'dual_jwh', 'dual_087'],
                    default=None,
                    help='Density/model-free variant; field count follows --dual (default: 087)')
    mf.add_argument('--init-amp', type=float, default=None, dest='init_amp',
                    help='Initial fit amplitude (default: derived from the data)')
    mf.add_argument('--init-t1', type=float, default=800.0, dest='init_t1', help='Initial T1 (ms)')
    mf.add_argument('--init-t2', type=float, default=100.0, dest='init_t2', help='Initial T2 (ms)')
    mf.add_argument('--n-bootstrap', type=int, default=1000, dest='n_bootstrap')
    mf.add_argument('--error-method', choices=['analytical', 'bootstrap'], default='analytical',
                    dest='error_method')
    mf.add_argument('--n-monte-carlo', type=int, default=50, dest='n_monte_carlo')
    _add_density_flags(mf)
    mf.set_defaults(func=_run_dynamixs_modelfree, _validate=_validate_modelfree)

    dx.set_defaults(func=lambda a: (dx.print_help(sys.stderr) or 2))

    export = sub.add_parser('export', help='Render figures / reports from a fit JSON (headless)')
    export_sub = export.add_subparsers(dest='export_command', metavar='<kind>')
    ex_kd = export_sub.add_parser('kd', parents=[common],
                                  help='CSP / intensity fit figures + summary from a kd fit JSON')
    ex_kd.add_argument('--json', required=True, help='kd fit JSON (…_kd_fit_data.json)')
    ex_kd.add_argument('--out', required=True, help='Output directory for figures + summary.csv')
    ex_kd.add_argument('--observable', type=_choice_list('csp', 'intensity'), default=None,
                       help='Comma-separated observables to render (default: those present)')
    ex_kd.add_argument('--fig-format', type=_choice_list('pdf', 'png'),
                       default=['pdf'], dest='fig_format',
                       help='pdf (multi-page grid per observable) and/or png (one file per '
                            'residue); comma-separated for both, e.g. pdf,png (default: pdf)')
    ex_kd.add_argument('--per-page', type=int, default=20, dest='per_page',
                       help='Panels per PDF page (default: 20 = 5x4, like T1/T2)')
    ex_kd.add_argument('--kind', type=_str_list,
                       default=['curves', 'ref-bars', 'kd-bars', 'global-fit'],
                       help="What to render: 'curves' (per-residue binding fits), "
                            "'ref-bars' (reference→point observable bars, one page per "
                            "point), 'kd-bars' (per-residue Kd + global-Kd line), and/or "
                            "'global-fit' (per-residue data + the shared-Kd global curve); "
                            "kd-bars is PDF only, ref-bars/global-fit also write CSV+JSON; "
                            "comma-separated (default: curves,ref-bars,kd-bars,global-fit)")
    ex_kd.add_argument('--summary-only', action='store_true', dest='summary_only',
                       help='Write only summary.csv, no figures')
    ex_kd.add_argument('--prefix', default='',
                       help='Output filename prefix, e.g. the sample name (default: none, '
                            'i.e. summary.csv/<obs>_fits.pdf/... unprefixed)')
    ex_kd.set_defaults(func=_run_export_kd)
    export.set_defaults(func=lambda a: (export.print_help(sys.stderr) or 2))

    proj = sub.add_parser('project', help='Inspect / prune a .lunaNMR project bundle')
    proj_sub = proj.add_subparsers(dest='project_command', metavar='<action>')
    inv = proj_sub.add_parser('inventory', parents=[fmt],
                              help='List a bundle\'s contents by category')
    inv.add_argument('bundle', help='Path to a .lunaNMR bundle directory')
    inv.set_defaults(func=_run_project_inventory)
    rm = proj_sub.add_parser('remove', parents=[common],
                             help='Delete bundle-relative paths from a bundle')
    rm.add_argument('bundle', help='Path to a .lunaNMR bundle directory')
    rm.add_argument('paths', nargs='+', help='Bundle-relative path(s) to delete')
    rm.set_defaults(func=_run_project_remove)
    proj.set_defaults(func=lambda a: (proj.print_help(sys.stderr) or 2))

    # `batch` is intercepted in main() before argparse so the batch CLI owns all of its
    # own flags (including -h). This entry exists only so it appears in the top-level help.
    sub.add_parser('batch', add_help=False,
                   help='Batch detect + Voigt/PS2D fit over a folder (pass -h for its flags)')

    return parser


def _add_relaxation_flags(p):
    """Shared flags for the dynamixs t1t2 / methyl-t2 subcommands."""
    p.add_argument('--input', required=True, help='Relaxation series CSV (LunaNMR or DynamiXs format)')
    p.add_argument('--out', required=True, help='Output directory for results + JSON')
    p.add_argument('--prefix', default='field1', help='Output filename prefix (default: field1)')
    p.add_argument('--field-name', default='field1', dest='field_name',
                   help='Field label used in the JSON filename (default: field1)')
    p.add_argument('--field-freq', type=float, default=600.0, dest='field_freq',
                   help='Spectrometer field frequency in MHz (default: 600)')
    p.add_argument('--time-units', choices=['ms', 's', 'us'], default=None, dest='time_units',
                   help='Units of the delay values (labels output; does not rescale). '
                        'Default: taken from series_metadata.json beside the input, else s')
    p.add_argument('--allow-unparsed-delays', action='store_true', dest='allow_unparsed',
                   help='Fit anyway when the series sidecar records spectra whose '
                        'filename carried no delay (those points are dropped)')
    p.add_argument('--error-method', choices=['analytical', 'bootstrap'], default='analytical',
                   dest='error_method', help='Error estimation method (default: analytical)')
    p.add_argument('--bootstrap', type=int, default=1000,
                   help='Bootstrap iterations when --error-method bootstrap (default: 1000)')


def _add_density_flags(p):
    """Physical/spectrometer flags shared by the density and modelfree subcommands."""
    p.add_argument('--field1-freq', type=float, default=600.0, dest='field1_freq',
                   help='Field-1 1H frequency in MHz (default: 600)')
    p.add_argument('--field2-freq', type=float, default=None, dest='field2_freq',
                   help='Field-2 1H frequency in MHz. Required with --dual on modelfree; '
                        'no default, because guessing it reads as a tau_c mismatch')
    p.add_argument('--rnh', type=float, default=1.015,
                   help='N-H bond length in Angstrom (default: 1.015)')
    p.add_argument('--csa', type=float, default=-172.0,
                   help='15N CSA in ppm (default: -172)')


def _version():
    from lunaNMR import __version__
    return __version__


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    # Delegate `batch` before argparse so the batch CLI parses its own flags cleanly.
    if argv and argv[0] == 'batch':
        return _run_batch(argv[1:])

    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, 'command', None):
        parser.print_help(sys.stderr)
        return 2
    validate = getattr(args, '_validate', None)
    if validate:
        validate(parser, args)
    try:
        return args.func(args)
    except (FileNotFoundError, ValueError, RuntimeError, KeyError, TypeError) as exc:
        # Expected bad-input failures from the wrapped engines (missing file, bad
        # concentrations, malformed CSV/JSON missing an expected column/key, an
        # unparseable delay label): report cleanly instead of dumping a traceback.
        _emit_error(args, exc)
        return 1


if __name__ == '__main__':
    sys.exit(main())
