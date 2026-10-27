# ABOUTME: Builds `python -m lunaNMR` argument vectors and runs them in a subprocess for the MCP server.
# ABOUTME: The CLI's JSON envelope is the contract; a run that prints none becomes an error envelope.

import asyncio
import json
import sys

STDERR_TAIL_LINES = 40


def diagnose_args(path, quick=False):
    args = ['diagnose', str(path), '--format', 'json']
    if quick:
        args.append('--quick')
    return args


def plan_series_args(spectra, peaks, out, mode='time', peak_source='reference', deep=True):
    args = ['series', '--dry-run', '--spectra', str(spectra), '--peaks', str(peaks),
            '--out', str(out), '--mode', mode, '--peak-source', peak_source,
            '--format', 'json']
    if deep:
        args.append('--deep')
    return args


def kd_survey_args(input_file, p0, conc_units, out, conc=None):
    args = ['kd', '--survey', '--input', str(input_file), '--out', str(out),
            '--p0', str(p0), '--conc-units', conc_units, '--format', 'json']
    if conc is not None:
        args += ['--conc', ','.join(str(value) for value in conc)]
    return args


def export_kd_args(fit_json, out, kinds=None):
    args = ['export', 'kd', '--json', str(fit_json), '--out', str(out), '--format', 'json']
    if kinds is not None:
        args += ['--kind', ','.join(kinds)]
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
        tail = '\n'.join(stderr.decode().splitlines()[-STDERR_TAIL_LINES:])
        return {'ok': False,
                'error': {'type': 'NoJsonOutput',
                          'message': f'exit {process.returncode}, no JSON on stdout. {tail}'}}
