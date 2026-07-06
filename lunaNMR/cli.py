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
import os
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
