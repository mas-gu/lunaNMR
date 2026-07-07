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


def _emit(args, summary, *human_lines):
    """Print the run summary: a single JSON object under --format json, else human lines."""
    import json
    if getattr(args, 'format', 'text') == 'json':
        print(json.dumps(summary))
    else:
        for line in human_lines:
            print(line)


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


def _dry_run(args, inputs, planned):
    """Validate that required input paths exist and report the plan without running.

    `inputs` is a list of (label, path); `planned` a dict of planned outputs. Returns 0
    if every input exists, else 1.
    """
    missing = [p for _, p in inputs if not os.path.exists(p)]
    summary = {
        'command': getattr(args, 'command', None),
        'dry_run': True,
        'inputs': {label: path for label, path in inputs},
        'planned_outputs': planned,
        'missing_inputs': missing,
    }
    human = [f"[dry-run] {summary['command']}"]
    for label, path in inputs:
        human.append(f"  input  {label}: {path} [{'OK' if os.path.exists(path) else 'MISSING'}]")
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
    return path


def _run_kd(args):
    """Wrap dynamiXs_Kd.kd_fit.run_kd_analysis_with_params for the CLI."""
    if args.dry_run:
        return _dry_run(args, [('input', args.input)], {'output_dir': args.out})
    _add_modules_path('dynamiXs_v2o0', 'dynamiXs_Kd')
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
    if args.conc is not None:
        params['concentrations'] = args.conc
    if args.intensity_scale is not None:
        params['intensity_scales'] = args.intensity_scale

    with _engine_stdout(args):
        result = kd_fit.run_kd_analysis_with_params(params)
    _emit(args,
          {'command': 'kd', 'n_fitted': result['n_fitted'], 'n_total': result['n_total'],
           'json_file': result['json_file'], 'results_file': result['results_file']},
          f"Kd analysis complete: {result['n_fitted']}/{result['n_total']} residues fitted",
          f"  JSON:    {result['json_file']}",
          f"  Results: {result['results_file']}")
    return 0


def _run_dynamixs_t1t2(args):
    """Wrap dynamiXs_T1_T2.fit_Tx_NMRRE.run_analysis_with_params (T1/T2 relaxation)."""
    if args.dry_run:
        return _dry_run(args, [('input', args.input)],
                        {'output_dir': args.out, 'experiment': args.exp})
    _add_modules_path('dynamiXs_v2o0', 'dynamiXs_T1_T2')
    from fit_Tx_NMRRE import run_analysis_with_params
    os.makedirs(args.out, exist_ok=True)
    params = {
        'input_csv_file': args.input,
        'output_prefix': os.path.join(args.out, args.prefix),
        'results_txt_file': os.path.join(args.out, f"{args.prefix}_fit_results.txt"),
        'experiment_type': args.exp,
        'error_method': args.error_method,
        'n_bootstrap': args.bootstrap,
        'field_name': args.field_name,
        'field_freq': args.field_freq,
        'json_folder': None if args.no_json else args.out,
    }
    with _engine_stdout(args):
        result = run_analysis_with_params(params)
    human = [f"{args.exp} analysis complete: {result['n_fitted']} residues fitted, "
             f"mean {args.exp} = {result['mean_t2']:.2f} ms",
             f"  Results: {result['results_file']}"]
    if result.get('json_file'):
        human.append(f"  JSON:    {result['json_file']}")
    _emit(args,
          {'command': 'dynamixs t1t2', 'experiment': args.exp, 'n_fitted': result['n_fitted'],
           'mean_t2': result['mean_t2'], 'results_file': result['results_file'],
           'json_file': result.get('json_file')},
          *human)
    return 0


def _run_dynamixs_methyl(args):
    """Wrap dynamiXs_T1_T2.fit_methyl_T2.run_methyl_t2_analysis_with_params (bi-exp methyl T2)."""
    if args.dry_run:
        return _dry_run(args, [('input', args.input)], {'output_dir': args.out})
    _add_modules_path('dynamiXs_v2o0', 'dynamiXs_T1_T2')
    from fit_methyl_T2 import run_methyl_t2_analysis_with_params
    os.makedirs(args.out, exist_ok=True)
    params = {
        'input_csv_file': args.input,
        'output_prefix': os.path.join(args.out, args.prefix),
        'results_txt_file': os.path.join(args.out, f"{args.prefix}_fit_results.txt"),
        'json_folder': args.out,
        'field_name': args.field_name,
        'field_freq': args.field_freq,
        'error_method': args.error_method,
        'n_bootstrap': args.bootstrap,
    }
    with _engine_stdout(args):
        result = run_methyl_t2_analysis_with_params(params)
    _emit(args,
          {'command': 'dynamixs methyl-t2', 'n_fitted': result['n_fitted'],
           'n_total': result['n_total'], 'results_file': result['results_file'],
           'json_file': result['json_file']},
          f"Methyl T2 analysis complete: {result['n_fitted']}/{result['n_total']} residues fitted",
          f"  Results: {result['results_file']}",
          f"  JSON:    {result['json_file']}")
    return 0


def _natural_key(text):
    """Sort key that orders embedded numbers numerically (s_2 before s_10)."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', text)]


def _discover_spectra(spectra, extensions=('ft', 'ser', 'ft2', 'ft3', 'pipe', 'ucsf')):
    """Resolve --spectra (a folder or a glob) to a naturally-sorted list of spectrum files.

    `extensions` defaults to NMRFileManager.supported_nmr_formats (passed explicitly by
    the series handler) so discovery never diverges from what the loader accepts.
    """
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
            'search_window_x': 0.08,
            'search_window_y': 0.8,
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
        _emit(args,
              {'command': 'series', 'dry_run': True, 'spectra_found': len(nmr_files),
               'peaks': args.peaks, 'peaks_exists': peaks_ok, 'mode': args.mode,
               'peak_source': args.peak_source, 'parallel': args.parallel,
               'output_dir': args.out, 'missing_inputs': missing},
              "[dry-run] series",
              f"  spectra found: {len(nmr_files)}",
              f"  peaks: {args.peaks} [{'OK' if peaks_ok else 'MISSING'}]",
              f"  mode={args.mode} peak-source={args.peak_source} parallel={args.parallel}",
              f"  output: {args.out}")
        return 0 if not missing else 1

    if not nmr_files:
        print(f"No spectrum files found in {args.spectra}", file=sys.stderr)
        return 1
    reference_peaks = file_manager.load_peak_list(args.peaks)
    if reference_peaks.empty:
        print(f"Peak list is empty: {args.peaks}", file=sys.stderr)
        return 1

    os.makedirs(args.out, exist_ok=True)
    print(f"Processing {len(nmr_files)} spectra ({args.mode} mode, {args.peak_source} peaks)...",
          file=sys.stderr)
    processor = MultiSpectrumProcessor(_default_series_params(parallel=args.parallel))
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
        print("Series produced no successful fits (all spectra failed or none loaded)",
              file=sys.stderr)
        return 1
    output_folder = result.metadata.get('output_folder', args.out)
    _emit(args,
          {'command': 'series', 'spectra_fitted': n_success, 'spectra_total': len(results),
           'output_folder': output_folder, 'parallel': args.parallel},
          f"Series analysis complete: {n_success}/{len(results)} spectra fitted",
          f"  Output: {output_folder}")
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
        print(f"Not a project bundle: {args.bundle}", file=sys.stderr)
        return 1
    categories = _project_manager().inventory(args.bundle)
    if not categories:
        print("(empty bundle — no recognized categories)")
        return 0
    for cat in categories:
        print(f"{cat['label']}  [{_human_size(cat['size'])}]")
        for item in cat['items']:
            lock = '' if item['removable'] else '  (protected)'
            print(f"  - {item['label']}  [{_human_size(item['size'])}]{lock}")
    return 0


def _run_project_remove(args):
    if not os.path.isdir(args.bundle):
        print(f"Not a project bundle: {args.bundle}", file=sys.stderr)
        return 1
    try:
        freed = _project_manager().remove_bundle_paths(args.bundle, args.paths)
    except ValueError as exc:
        print(f"Refused: {exc}", file=sys.stderr)
        return 1
    print(f"Removed {len(args.paths)} path(s), freed {_human_size(freed)}")
    return 0


def _safe_name(text):
    """Filesystem-safe version of a residue label."""
    return re.sub(r'[^\w.+-]+', '_', str(text))


def _run_export_kd(args):
    """Render CSP / intensity fit figures + a summary from a self-contained kd fit JSON."""
    import json
    import csv
    if args.dry_run:
        return _dry_run(args, [('json', args.json)], {'output_dir': args.out})
    if not os.path.isfile(args.json):
        print(f"Fit JSON not found: {args.json}", file=sys.stderr)
        return 1
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    _add_modules_path('dynamiXs_v2o0', 'dynamiXs_Kd')
    from kd_models import csp_model, intensity_decay

    with open(args.json) as fh:
        data = json.load(fh)
    fits = data.get('fits', [])
    P0 = data.get('metadata', {}).get('protein_conc')
    observables = args.observable or [o for o in ('csp', 'intensity')
                                      if any(f.get(o) for f in fits)]
    os.makedirs(args.out, exist_ok=True)

    summary_rows = []
    n_figs = 0
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
            if args.summary_only:
                continue
            L = np.asarray(fit['L'], dtype=float)
            y = np.asarray(fit['obs'], dtype=float)
            good = np.isfinite(L) & np.isfinite(y)
            if good.sum() < 2:
                continue
            Ld = np.linspace(float(L[good].min()), float(L[good].max()), 200)
            if obs == 'csp':
                yc = csp_model(Ld, fit['dd_max'], fit['Kd'], P0)
                ylabel = 'CSP (ppm)'
            else:
                yc = intensity_decay(Ld, fit['I0'], fit['I_inf'], fit['Kd'])
                ylabel = 'I / I(0)' if max(y[good]) <= 1.5 else 'Intensity'
            fig, ax = plt.subplots(figsize=(4, 3))
            ax.plot(L[good], y[good], 'o', color='#1f77b4', label='observed')
            ax.plot(Ld, yc, '-', color='#d62728',
                    label=f"Kd={fit['Kd']:.2g}, R²={fit.get('r_squared', float('nan')):.3f}")
            ax.set_xlabel('[ligand]')
            ax.set_ylabel(ylabel)
            ax.set_title(f"{residue} ({obs})")
            ax.legend(fontsize=7)
            fig.tight_layout()
            obs_dir = os.path.join(args.out, obs)
            os.makedirs(obs_dir, exist_ok=True)
            fig.savefig(os.path.join(obs_dir, f"{_safe_name(residue)}.png"), dpi=120)
            plt.close(fig)
            n_figs += 1

    summary_path = os.path.join(args.out, 'summary.csv')
    with open(summary_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['residue', 'observable', 'Kd', 'Kd_err', 'r_squared'])
        writer.writeheader()
        writer.writerows(summary_rows)

    _emit(args,
          {'command': 'export kd', 'n_fits': len(summary_rows), 'n_figures': n_figs,
           'summary_csv': summary_path},
          f"Exported {len(summary_rows)} fit(s), {n_figs} figure(s)",
          f"  Summary: {summary_path}")
    return 0


def _run_batch(batch_argv):
    """Delegate to the existing batch CLI, passing through all of its flags."""
    from lunaNMR.batch_processing.cli_interface import CLIInterface
    return CLIInterface().main(batch_argv)


def build_parser():
    parser = argparse.ArgumentParser(
        prog='lunaNMR',
        description='LunaNMR command-line interface for headless NMR analysis.',
    )
    parser.add_argument('--version', action='version', version=f'lunaNMR {_version()}')
    sub = parser.add_subparsers(dest='command', metavar='<subcommand>')

    # Shared flags for the analysis subcommands: output format + input validation.
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('--format', choices=['text', 'json'], default='text',
                        help='Run-summary output format (default: text)')
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
    kd.add_argument('--observable', type=_str_list, default=['csp', 'intensity'],
                    help='Comma-separated observables to fit: csp,intensity (default: both)')
    kd.add_argument('--intensity-from', choices=['height', 'volume'], default='height',
                    dest='intensity_from', help='Intensity source for the ratio (default: height)')
    kd.add_argument('--bootstrap', type=int, default=0,
                    help='Bootstrap iterations for error estimates (default: 0)')
    kd.add_argument('--intensity-scale', type=_float_list, default=None,
                    help='Comma-separated per-point height/volume scale factors')
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
    series.add_argument('--parallel', action='store_true',
                        help='Use the two-pass parallel processor (~2.7x faster)')
    series.set_defaults(func=_run_series)

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

    dx.set_defaults(func=lambda a: (dx.print_help(sys.stderr) or 2))

    export = sub.add_parser('export', help='Render figures / reports from a fit JSON (headless)')
    export_sub = export.add_subparsers(dest='export_command', metavar='<kind>')
    ex_kd = export_sub.add_parser('kd', parents=[common],
                                  help='CSP / intensity fit figures + summary from a kd fit JSON')
    ex_kd.add_argument('--json', required=True, help='kd fit JSON (…_kd_fit_data.json)')
    ex_kd.add_argument('--out', required=True, help='Output directory for figures + summary.csv')
    ex_kd.add_argument('--observable', type=_str_list, default=None,
                       help='Comma-separated observables to render (default: those present)')
    ex_kd.add_argument('--summary-only', action='store_true', dest='summary_only',
                       help='Write only summary.csv, no figures')
    ex_kd.set_defaults(func=_run_export_kd)
    export.set_defaults(func=lambda a: (export.print_help(sys.stderr) or 2))

    proj = sub.add_parser('project', help='Inspect / prune a .lunaNMR project bundle')
    proj_sub = proj.add_subparsers(dest='project_command', metavar='<action>')
    inv = proj_sub.add_parser('inventory', help='List a bundle\'s contents by category')
    inv.add_argument('bundle', help='Path to a .lunaNMR bundle directory')
    inv.set_defaults(func=_run_project_inventory)
    rm = proj_sub.add_parser('remove', help='Delete bundle-relative paths from a bundle')
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
    p.add_argument('--error-method', choices=['analytical', 'bootstrap'], default='analytical',
                   dest='error_method', help='Error estimation method (default: analytical)')
    p.add_argument('--bootstrap', type=int, default=1000,
                   help='Bootstrap iterations when --error-method bootstrap (default: 1000)')


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
    try:
        return args.func(args)
    except (FileNotFoundError, ValueError, RuntimeError, KeyError) as exc:
        # Expected bad-input failures from the wrapped engines (missing file, bad
        # concentrations, malformed CSV/JSON missing an expected column/key): report
        # cleanly instead of dumping a traceback.
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
