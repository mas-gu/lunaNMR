# ABOUTME: MCP tools for a dataset's spectra: the pre-flight check and the series fit.
# ABOUTME: Every folder they scan is checked for links that leave the data root.

from typing import Literal

from mcp.types import CallToolResult

from lunaNMR.mcp import cli_call
from lunaNMR.mcp.context import REFUSALS
from lunaNMR.mcp.results import READS_ONLY, WRITES_NEW_FOLDER, failure


def register(server, ctx):
    root = ctx.root

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
        except REFUSALS as error:
            return failure(error)
        return await ctx.run(cli_call.diagnose_args(folder, quick=quick, mode=mode))

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
        except REFUSALS as error:
            return failure(error)
        return await ctx.run(cli_call.plan_series_args(
            spectra_folder, peak_list, out_folder, mode=mode, peak_source=peak_source,
            deep=deep))

    @server.tool(annotations=WRITES_NEW_FOLDER)
    async def start_series(spectra: str, peaks: str, out: str,
                           mode: Literal['time', 'titration'] = 'time',
                           peak_source: Literal['reference', 'cascade', 'detected',
                                                'independent'] = 'reference',
                           parallel: bool = False) -> CallToolResult:
        """Start a series fit in the background; a quick one returns finished, a slow one
        returns a job to poll with `job_status`. Run `plan_series` first. `out` is a new
        folder. Many peaks per spectrum may justify `parallel`; on a short list it is slower."""
        try:
            spectra_folder = root.existing(spectra)
            root.refuse_escaping_links(spectra_folder)
            peak_list = root.existing(peaks)
            out_folder = root.new_output(out)
            return await ctx.start_job('series', cli_call.start_series_args(
                spectra_folder, peak_list, out_folder, mode=mode, peak_source=peak_source,
                parallel=parallel), out_folder)
        except REFUSALS as error:
            return failure(error)
