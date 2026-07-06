# ABOUTME: Unified command-line entry point for lunaNMR (`python -m lunaNMR <subcommand>`).
# ABOUTME: Dispatches to headless analysis engines; heavy deps are imported lazily per subcommand.

"""Command-line interface for lunaNMR.

Subcommands:
  kd      Fit binding affinity (Kd) from a titration series (CSP / intensity).
  batch   Batch peak detection + Voigt/PS2D fitting over a folder of spectra.

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

    result = kd_fit.run_kd_analysis_with_params(params)
    print(f"Kd analysis complete: {result['n_fitted']}/{result['n_total']} residues fitted")
    print(f"  JSON:    {result['json_file']}")
    print(f"  Results: {result['results_file']}")
    return 0


def _run_dynamixs_t1t2(args):
    """Wrap dynamiXs_T1_T2.fit_Tx_NMRRE.run_analysis_with_params (T1/T2 relaxation)."""
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
    result = run_analysis_with_params(params)
    print(f"{args.exp} analysis complete: {result['n_fitted']} residues fitted, "
          f"mean {args.exp} = {result['mean_t2']:.2f} ms")
    print(f"  Results: {result['results_file']}")
    if result.get('json_file'):
        print(f"  JSON:    {result['json_file']}")
    return 0


def _run_dynamixs_methyl(args):
    """Wrap dynamiXs_T1_T2.fit_methyl_T2.run_methyl_t2_analysis_with_params (bi-exp methyl T2)."""
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
    result = run_methyl_t2_analysis_with_params(params)
    print(f"Methyl T2 analysis complete: {result['n_fitted']}/{result['n_total']} residues fitted")
    print(f"  Results: {result['results_file']}")
    print(f"  JSON:    {result['json_file']}")
    return 0


def _natural_key(text):
    """Sort key that orders embedded numbers numerically (s_2 before s_10)."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', text)]


def _discover_spectra(spectra):
    """Resolve --spectra (a folder or a glob) to a naturally-sorted list of spectrum files."""
    if os.path.isdir(spectra):
        files = []
        for ext in ('ft', 'ft2', 'fid'):
            files.extend(glob.glob(os.path.join(spectra, f'*.{ext}')))
    else:
        files = glob.glob(spectra)
    return sorted(files, key=lambda f: _natural_key(os.path.basename(f)))


def _default_series_params():
    """Default series/voigt parameters, mirroring the GUI's getattr fallbacks.

    The GUI assembles this nested dict from the main window with these same defaults
    (gui/dialogs/series_integration_dialog.py::_get_voigt_parameters); a headless run
    just uses the defaults directly.
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
            'use_parallel_processing': False,
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
            'use_parallel_processing': False,
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

    nmr_files = _discover_spectra(args.spectra)
    if not nmr_files:
        print(f"No spectrum files (.ft/.ft2/.fid) found in {args.spectra}", file=sys.stderr)
        return 1
    reference_peaks = NMRFileManager().load_peak_list(args.peaks)
    if reference_peaks is None or reference_peaks.empty:
        print(f"Failed to load peak list: {args.peaks}", file=sys.stderr)
        return 1

    os.makedirs(args.out, exist_ok=True)
    print(f"Processing {len(nmr_files)} spectra ({args.mode} mode, {args.peak_source} peaks)...")
    processor = MultiSpectrumProcessor(_default_series_params())
    result = processor.process_nmr_series(
        nmr_files, reference_peaks, args.out,
        peak_source_mode=args.peak_source, series_mode=args.mode, extract_delays=True,
    )
    if getattr(result, 'errors', None):
        for err in result.errors:
            print(f"  error: {err}", file=sys.stderr)
    output_folder = result.metadata.get('output_folder', args.out)
    print(f"Series analysis complete: {len(getattr(result, 'results', {}) or {})} spectra processed")
    print(f"  Output: {output_folder}")
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

    kd = sub.add_parser('kd', help='Fit Kd from a titration series (CSP / intensity)')
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

    series = sub.add_parser('series', help='Process a multi-spectrum series/titration')
    series.add_argument('--spectra', required=True, help='Folder or glob of spectrum files (.ft/.ft2/.fid)')
    series.add_argument('--peaks', required=True, help='Reference peak-list file (Assignment, Position_X, Position_Y)')
    series.add_argument('--out', required=True, help='Output folder for the series_results CSVs')
    series.add_argument('--mode', choices=['time', 'titration'], default='time',
                        help='Series type: time (relaxation) or titration (default: time)')
    series.add_argument('--peak-source', choices=['reference', 'cascade', 'detected', 'independent'],
                        default='reference', dest='peak_source',
                        help='Peak position source across spectra (default: reference)')
    series.set_defaults(func=_run_series)

    dx = sub.add_parser('dynamixs', help='Relaxation fitting: T1/T2 and methyl-T2')
    dx_sub = dx.add_subparsers(dest='dynamixs_command', metavar='<kind>')

    t1t2 = dx_sub.add_parser('t1t2', help='Mono-exponential T1 or T2 relaxation fit')
    _add_relaxation_flags(t1t2)
    t1t2.add_argument('--exp', choices=['T1', 'T2'], required=True, help='Experiment type')
    t1t2.add_argument('--no-json', action='store_true', help='Skip writing the JSON fit data')
    t1t2.set_defaults(func=_run_dynamixs_t1t2)

    methyl = dx_sub.add_parser('methyl-t2', help='Bi-exponential (Tugarinov-Kay) methyl T2 fit')
    _add_relaxation_flags(methyl)
    methyl.set_defaults(func=_run_dynamixs_methyl)

    dx.set_defaults(func=lambda a: (dx.print_help(sys.stderr) or 2))

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
    p.add_argument('--bootstrap', type=int, default=0,
                   help='Bootstrap iterations when --error-method bootstrap (default: 0)')


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
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
