# ABOUTME: Builds `python -m lunaNMR` argument vectors and runs them in a subprocess for the MCP server.
# ABOUTME: The CLI's JSON envelope is the contract; a run that prints none becomes an error envelope.

import asyncio
import json
import re
import sys

STDERR_TAIL_LINES = 40
MESSAGE_CHARS = 1500
RESIDUE_NAME = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.+\-]{0,63}')


class InvalidResidueName(ValueError):
    """A residue name that could be read as something else (a path, a comment, a flag)."""


def selection_text(residues):
    """The text of a selection file naming these residues, one per line.

    The CLI reads a `--residues` value as a file when one exists, so names never travel on the
    command line: they are written to a file this server owns.
    """
    if not residues:
        raise InvalidResidueName('no residue names given')
    for name in residues:
        if not RESIDUE_NAME.fullmatch(name):
            raise InvalidResidueName(f'not a residue name: {name[:40]!r}')
    return '\n'.join(residues) + '\n'


def _flag(args, name, enabled):
    if enabled:
        args.append(name)


def _value(args, name, value):
    if value is not None:
        args += [name, str(value)]


def _csv(values):
    return ','.join(str(value) for value in values)


def diagnose_args(path, quick=False, mode='time', sample=False, strict=False):
    args = ['diagnose', str(path), '--mode', mode, '--format', 'json']
    _flag(args, '--quick', quick)
    _flag(args, '--sample', sample)
    _flag(args, '--strict', strict)
    return args


def plan_series_args(spectra, peaks, out, mode='time', peak_source='reference', deep=True,
                     quick=False, sample=False, params_file=None):
    args = ['series', '--dry-run', '--spectra', str(spectra), '--peaks', str(peaks),
            '--out', str(out), '--mode', mode, '--peak-source', peak_source,
            '--format', 'json']
    _flag(args, '--deep', deep)
    _flag(args, '--quick', quick)
    _flag(args, '--sample', sample)
    _value(args, '--params', params_file)
    return args


def start_series_args(spectra, peaks, out, mode='time', peak_source='reference', parallel=False,
                      params_file=None):
    args = ['series', '--spectra', str(spectra), '--peaks', str(peaks), '--out', str(out),
            '--mode', mode, '--peak-source', peak_source, '--format', 'json']
    _flag(args, '--parallel', parallel)
    _value(args, '--params', params_file)
    return args


def _kd_args(input_file, p0, conc_units, out, conc, observables, intensity_scales, prefix,
             alpha, intensity_from, bootstrap, csp_sigma_multiple, kd_outlier_z,
             noise_quantile, dd_runaway_ratio, ref_max_ratio):
    args = ['kd', '--input', str(input_file), '--out', str(out), '--p0', str(p0),
            '--conc-units', conc_units, '--format', 'json']
    if conc is not None:
        args += ['--conc', _csv(conc)]
    if observables is not None:
        args += ['--observable', _csv(observables)]
    if intensity_scales is not None:
        args += ['--intensity-scale', _csv(intensity_scales)]
    _value(args, '--prefix', prefix)
    _value(args, '--alpha', alpha)
    _value(args, '--intensity-from', intensity_from)
    _value(args, '--bootstrap', bootstrap)
    _value(args, '--csp-sigma-multiple', csp_sigma_multiple)
    _value(args, '--kd-outlier-z', kd_outlier_z)
    _value(args, '--noise-quantile', noise_quantile)
    _value(args, '--dd-runaway-ratio', dd_runaway_ratio)
    _value(args, '--ref-max-ratio', ref_max_ratio)
    return args


def kd_survey_args(input_file, p0, conc_units, out, conc=None, intensity_scales=None,
                   selection_file=None, observables=None, prefix=None, alpha=None,
                   intensity_from=None, bootstrap=None, csp_sigma_multiple=None,
                   kd_outlier_z=None, noise_quantile=None, dd_runaway_ratio=None,
                   ref_max_ratio=None):
    args = _kd_args(input_file, p0, conc_units, out, conc, observables, intensity_scales, prefix,
                    alpha, intensity_from, bootstrap, csp_sigma_multiple, kd_outlier_z,
                    noise_quantile, dd_runaway_ratio, ref_max_ratio)
    args.insert(1, '--survey')
    _value(args, '--selection', selection_file)
    return args


def kd_fit_args(input_file, p0, conc_units, out, conc=None, residues_file=None,
                observables=None, intensity_scales=None, prefix=None, alpha=None,
                intensity_from=None, bootstrap=None, csp_sigma_multiple=None,
                kd_outlier_z=None, noise_quantile=None, dd_runaway_ratio=None,
                ref_max_ratio=None):
    args = _kd_args(input_file, p0, conc_units, out, conc, observables, intensity_scales, prefix,
                    alpha, intensity_from, bootstrap, csp_sigma_multiple, kd_outlier_z,
                    noise_quantile, dd_runaway_ratio, ref_max_ratio)
    _value(args, '--residues', residues_file)
    return args


def export_kd_args(fit_json, out, kinds=None, observables=None, fig_format=None, per_page=None,
                   prefix=None, summary_only=False):
    args = ['export', 'kd', '--json', str(fit_json), '--out', str(out), '--format', 'json']
    if kinds is not None:
        args += ['--kind', _csv(kinds)]
    if observables is not None:
        args += ['--observable', _csv(observables)]
    if fig_format is not None:
        args += ['--fig-format', _csv(fig_format)]
    _value(args, '--per-page', per_page)
    _value(args, '--prefix', prefix)
    _flag(args, '--summary-only', summary_only)
    return args


def peaks_shift_args(peaks, out, dx=None, dy=None, auto=False, spectrum=None, quick=False):
    args = ['peaks', 'shift', '--peaks', str(peaks), '--out', str(out), '--format', 'json']
    _value(args, '--dx', dx)
    _value(args, '--dy', dy)
    _flag(args, '--auto', auto)
    _value(args, '--spectrum', spectrum)
    _flag(args, '--quick', quick)
    return args


def integrate_1d_args(spectra, out, peaks=None, detect=False, prefix=None, min_snr=None,
                      value=None):
    args = ['integrate-1d', '--spectra', str(spectra), '--out', str(out), '--format', 'json']
    if peaks is not None:
        args += ['--peaks', _csv(peaks)]
    _flag(args, '--detect', detect)
    _value(args, '--prefix', prefix)
    _value(args, '--min-snr', min_snr)
    _value(args, '--value', value)
    return args


def trimmed(envelope):
    """The envelope without the parts that grow with the dataset and that a client does not
    act on: the echoed parameter block and each peak list's assignment names (the count stays)."""
    envelope = {key: value for key, value in envelope.items() if key != 'params'}
    peak_lists = [envelope.get('checks', {}).get('peak_list')]
    peak_lists += [experiment.get('peak_list') for experiment in envelope.get('experiments', [])]
    for peak_list in peak_lists:
        if isinstance(peak_list, dict):
            peak_list.pop('assignments', None)
    return envelope


async def run_cli(args):
    """Run one CLI command and return its JSON envelope.

    A subprocess keeps the engines' prints, and any spawn worker's, off the server's
    stdout, which carries the MCP protocol.
    """
    process = await asyncio.create_subprocess_exec(
        sys.executable, '-m', 'lunaNMR', *args,
        stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE)
    stdout, stderr = await process.communicate()
    try:
        return trimmed(json.loads(stdout.decode()))
    except ValueError:
        tail = '\n'.join(stderr.decode(errors='replace').splitlines()[-STDERR_TAIL_LINES:])
        message = f'exit {process.returncode}, no JSON on stdout. {tail}'
        return {'ok': False,
                'error': {'type': 'NoJsonOutput', 'message': message[-MESSAGE_CHARS:]}}
