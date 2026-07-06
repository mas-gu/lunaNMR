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


def _run_kd(args):
    """Wrap dynamiXs_Kd.kd_fit.run_kd_analysis_with_params for the CLI."""
    # kd_fit uses top-level sibling imports (from kd_models import ...), so its
    # directory must be on sys.path (mirrors modules/dynamiXs_v2o0/workers.py).
    kd_dir = os.path.join(os.path.dirname(__file__), '..', 'modules',
                          'dynamiXs_v2o0', 'dynamiXs_Kd')
    kd_dir = os.path.abspath(kd_dir)
    if kd_dir not in sys.path:
        sys.path.insert(0, kd_dir)
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

    # `batch` is intercepted in main() before argparse so the batch CLI owns all of its
    # own flags (including -h). This entry exists only so it appears in the top-level help.
    sub.add_parser('batch', add_help=False,
                   help='Batch detect + Voigt/PS2D fit over a folder (pass -h for its flags)')

    return parser


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
