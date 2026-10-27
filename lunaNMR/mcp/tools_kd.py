# ABOUTME: MCP tools for a binding titration: survey the residues, fit the Kd, export figures.
# ABOUTME: Residue names are written to a selection file the server owns, never passed as paths.

from typing import Literal

from mcp.types import CallToolResult

from lunaNMR.mcp import cli_call
from lunaNMR.mcp.context import REFUSALS
from lunaNMR.mcp.jobs import JOB_FOLDER
from lunaNMR.mcp.results import (READS_ONLY, REWRITES_BESIDE_INPUT, WRITES_NEW_FOLDER,
                                 failure)


def register(server, ctx):
    root = ctx.root

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
        survey table go to `out` (new folder); the residue list and a settings file are
        written beside `input`."""
        try:
            input_file = root.existing(input)
            root.refuse_escaping_links(input_file.parent)
            out_folder = root.new_output(out)
        except REFUSALS as error:
            return failure(error)
        return await ctx.run(cli_call.kd_survey_args(
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
        except REFUSALS as error:
            return failure(error)
        return await ctx.run(cli_call.export_kd_args(fit_file, out_folder, kinds=kinds))

    @server.tool(annotations=WRITES_NEW_FOLDER)
    async def kd_fit(input: str, p0: float, conc_units: Literal['absolute', 'equivalents'],
                     out: str, conc: list[float] | None = None,
                     residues: list[str] | None = None, residues_file: str | None = None,
                     observables: list[Literal['csp', 'intensity']] | None = None,
                     intensity_scales: list[float] | None = None) -> CallToolResult:
        """Start a Kd fit in the background; a quick one returns finished, a slow one
        returns a job to poll with `job_status`. Fit all residues, or only
        those named in `residues` or listed in `residues_file` (the selection file
        `kd_survey` writes beside its input) -- not both. `observables` defaults to both.
        The finished job reports each shared Kd with its `reliable` flag and how many
        residues are in the shared fit. `intensity_scales` is the same scan-count correction
        as in `kd_survey` and must match the survey's. `out` is a new folder."""
        if residues is not None and residues_file is not None:
            return failure(ValueError('give residues or residues_file, not both'))
        try:
            input_file = root.existing(input)
            root.refuse_escaping_links(input_file.parent)
            selection = root.existing(residues_file) if residues_file is not None else None
            named = cli_call.selection_text(residues) if residues is not None else None
            out_folder = root.new_output(out)
            if named is not None:
                selection = out_folder / JOB_FOLDER / 'residues.txt'
            return await ctx.start_job('kd', cli_call.kd_fit_args(
                input_file, p0, conc_units, out_folder, conc=conc, residues_file=selection,
                observables=observables, intensity_scales=intensity_scales), out_folder,
                files=None if named is None else {'residues.txt': named})
        except REFUSALS as error:
            return failure(error)
