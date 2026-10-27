# ABOUTME: The lunaNMR MCP server: typed tools over the headless CLI, confined to one data folder.
# ABOUTME: Run it as `lunanmr-mcp --root <dataset folder>`; it speaks MCP over stdio.

import json
import sys
from typing import Literal

from mcp.server.mcpserver import MCPServer
from mcp.types import CallToolResult, TextContent, ToolAnnotations

from lunaNMR.mcp import check, cli_call
from lunaNMR.mcp.command_line import build_parser
from lunaNMR.mcp.data_root import DataRoot, PathRefused
from lunaNMR.mcp.jobs import JOB_FOLDER, JobManager, NotAJob, NotCancellable, TooManyJobs


MAX_RESULT_CHARS = 200_000


def _result(envelope):
    """A tool result whose error flag follows the envelope's `ok`, so a client treats a
    failed run as a failure and still reads the reason. Anything larger than a client should
    be handed is replaced by an error: the text comes from files the dataset controls."""
    text = json.dumps(envelope)
    if len(text) > MAX_RESULT_CHARS:
        envelope = {'ok': False, 'error': {
            'type': 'ResultTooLarge',
            'message': f'the result was {len(text)} characters; read the files in the output '
                       f'folder instead'}}
        text = json.dumps(envelope)
    return CallToolResult(
        content=[TextContent(type='text', text=text)],
        structured_content=envelope,
        is_error=not envelope.get('ok', False))


def _failure(error):
    return _result({'ok': False, 'error': {'type': type(error).__name__, 'message': str(error)}})


READS_ONLY = ToolAnnotations(read_only_hint=True, open_world_hint=False)
# Rewrites the residue list and settings file beside its input (merging, so decisions survive).
REWRITES_BESIDE_INPUT = ToolAnnotations(
    read_only_hint=False, destructive_hint=True, idempotent_hint=True, open_world_hint=False)
CANCELS_A_JOB = ToolAnnotations(
    read_only_hint=False, destructive_hint=True, idempotent_hint=True, open_world_hint=False)
WRITES_NEW_FOLDER = ToolAnnotations(
    read_only_hint=False, destructive_hint=False, idempotent_hint=False, open_world_hint=False)


async def _run(args):
    return _result(await cli_call.run_cli(args))


def build_server(root, max_running_jobs=2):
    server = MCPServer(
        'lunaNMR',
        instructions='NMR peak fitting and analysis on the data folder this server was '
                     'started with. Every path is relative to that folder.')
    jobs = JobManager(max_running_jobs)

    def _job_envelope(status):
        status['job'] = str(root.relative(status['job']))
        return _result({'ok': True, **status})

    @server.tool(annotations=READS_ONLY)
    async def diagnose(path: str = '.', quick: bool = False,
                       mode: Literal['time', 'titration'] = 'time') -> CallToolResult:
        """Read-only pre-flight over a dataset folder: per-spectrum registration offset and
        capture rate, delays, peak-list hygiene, repeat-acquisition scale. Run it before
        any fit. Changes nothing. `path` is the folder that holds the spectra. Set `mode` to
        'titration' when the filenames carry ligand ratios (0o5, 1o0), not time delays."""
        try:
            folder = root.existing(path)
            root.refuse_escaping_links(folder)
        except PathRefused as error:
            return _failure(error)
        return await _run(cli_call.diagnose_args(folder, quick=quick, mode=mode))

    @server.tool(annotations=READS_ONLY)
    async def plan_series(spectra: str, peaks: str, out: str,
                          mode: Literal['time', 'titration'] = 'time',
                          peak_source: Literal['reference', 'cascade', 'detected',
                                               'independent'] = 'reference',
                          deep: bool = True) -> CallToolResult:
        """Dry run of a series fit: checks the spectra folder, the peak list and the
        parameters, runs the per-spectrum registration and capture checks, and reports
        the plan. Writes nothing. `out` is the folder the real run would create and must
        not exist yet."""
        try:
            spectra_folder = root.existing(spectra)
            root.refuse_escaping_links(spectra_folder)
            peak_list = root.existing(peaks)
            out_folder = root.new_output(out)
        except PathRefused as error:
            return _failure(error)
        return await _run(cli_call.plan_series_args(
            spectra_folder, peak_list, out_folder, mode=mode, peak_source=peak_source,
            deep=deep))

    @server.tool(annotations=REWRITES_BESIDE_INPUT)
    async def kd_survey(input: str, p0: float, conc_units: Literal['absolute', 'equivalents'],
                        out: str, conc: list[float] | None = None,
                        intensity_scales: list[float] | None = None) -> CallToolResult:
        """Survey a titration before fitting a Kd: rates every residue ok / check /
        unusable and proposes a residue selection. `input` is series_analysis_tidy.csv.
        `p0` is the protein concentration and `conc_units` says whether the ligand
        concentrations are absolute or equivalents; `conc` gives them per titration point
        when the point labels are not numbers. `intensity_scales` is one factor per titration
        point applied to every peak's height and volume (a scan-count correction, e.g. 0.5 for
        a point acquired with twice the scans); it never changes positions. Figures and the
        survey table go to `out`
        (new folder); the residue list and a settings file are written beside `input`."""
        try:
            input_file = root.existing(input)
            root.refuse_escaping_links(input_file.parent)
            out_folder = root.new_output(out)
        except PathRefused as error:
            return _failure(error)
        return await _run(cli_call.kd_survey_args(
            input_file, p0, conc_units, out_folder, conc=conc,
            intensity_scales=intensity_scales))

    @server.tool(annotations=WRITES_NEW_FOLDER)
    async def export_kd(fit_json: str, out: str, kinds: list[Literal[
            'curves', 'ref-bars', 'kd-bars', 'global-fit']] | None = None) -> CallToolResult:
        """Figures and tables from a Kd fit JSON (`…_kd_fit_data.json`). `kinds` picks
        what to render (default: all four). `out` is a new folder."""
        try:
            fit_file = root.existing(fit_json)
            out_folder = root.new_output(out)
        except PathRefused as error:
            return _failure(error)
        return await _run(cli_call.export_kd_args(
            fit_file, out_folder, kinds=kinds))

    @server.tool(annotations=WRITES_NEW_FOLDER)
    async def start_series(spectra: str, peaks: str, out: str,
                           mode: Literal['time', 'titration'] = 'time',
                           peak_source: Literal['reference', 'cascade', 'detected',
                                                'independent'] = 'reference',
                           parallel: bool = False) -> CallToolResult:
        """Start a series fit in the background and return its job id at once; poll
        `job_status` until it finishes. Run `plan_series` first. `out` is a new folder.
        Many peaks per spectrum may justify `parallel`; on a short list it is slower."""
        try:
            spectra_folder = root.existing(spectra)
            root.refuse_escaping_links(spectra_folder)
            peak_list = root.existing(peaks)
            out_folder = root.new_output(out)
            started = jobs.start('series', cli_call.start_series_args(
                spectra_folder, peak_list, out_folder, mode=mode, peak_source=peak_source,
                parallel=parallel), out_folder)
        except (PathRefused, TooManyJobs) as error:
            return _failure(error)
        return _job_envelope(jobs.status(started))

    @server.tool(annotations=WRITES_NEW_FOLDER)
    async def kd_fit(input: str, p0: float, conc_units: Literal['absolute', 'equivalents'],
                     out: str, conc: list[float] | None = None,
                     residues: list[str] | None = None, residues_file: str | None = None,
                     observables: list[Literal['csp', 'intensity']] | None = None,
                     intensity_scales: list[float] | None = None) -> CallToolResult:
        """Start a Kd fit in the background; poll `job_status`. Fit all residues, or only
        those named in `residues` or listed in `residues_file` (the selection file
        `kd_survey` writes beside its input) -- not both. `observables` defaults to both.
        The finished job reports each shared Kd with its `reliable` flag and how many
        residues are in the shared fit. `intensity_scales` is the same scan-count correction
        as in `kd_survey` and must match the survey's. `out` is a new folder."""
        if residues is not None and residues_file is not None:
            return _failure(ValueError('give residues or residues_file, not both'))
        try:
            input_file = root.existing(input)
            root.refuse_escaping_links(input_file.parent)
            selection = root.existing(residues_file) if residues_file is not None else None
            named = cli_call.selection_text(residues) if residues is not None else None
            out_folder = root.new_output(out)
            if named is not None:
                selection = out_folder / JOB_FOLDER / 'residues.txt'
            started = jobs.start('kd', cli_call.kd_fit_args(
                input_file, p0, conc_units, out_folder, conc=conc, residues_file=selection,
                observables=observables, intensity_scales=intensity_scales), out_folder,
                files=None if named is None else {'residues.txt': named})
        except (PathRefused, TooManyJobs, cli_call.InvalidResidueName) as error:
            return _failure(error)
        return _job_envelope(jobs.status(started))

    @server.tool(annotations=READS_ONLY)
    async def job_status(job: str) -> CallToolResult:
        """State of a job started by `start_series` or `kd_fit`: running (with the last
        log lines), done, failed or cancelled. A finished job carries the command's summary
        and what its exit code hides: pinned and undetected points for a series, the shared
        Kd and its pool for a Kd fit. Works for jobs from an earlier server run."""
        try:
            return _job_envelope(jobs.status(root.existing(job)))
        except (PathRefused, NotAJob) as error:
            return _failure(error)

    @server.tool(annotations=CANCELS_A_JOB)
    async def cancel_job(job: str) -> CallToolResult:
        """Stop a running job started by this server. The output folder keeps whatever was
        written so far; start a new run in a new folder."""
        try:
            return _job_envelope(await jobs.cancel(root.existing(job)))
        except (PathRefused, NotAJob, NotCancellable) as error:
            return _failure(error)

    return server


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.check:
        sys.exit(check.report(check.run_checks()))
    if args.root is None:
        parser.error('--root is required unless --check is given')
    try:
        root = DataRoot(args.root)
    except PathRefused as error:
        parser.error(str(error))
    build_server(root).run()


if __name__ == '__main__':
    main()
