# ABOUTME: MCP tools for a binding titration: survey the residues, fit the Kd, export figures.
# ABOUTME: Residue names are written to a selection file the server owns, never passed as paths.

import inspect
from typing import Literal

from mcp.types import CallToolResult

from lunaNMR.mcp import cli_call
from lunaNMR.mcp.context import REFUSALS
from lunaNMR.mcp.jobs import JOB_FOLDER
from lunaNMR.mcp.names import safe_name
from lunaNMR.mcp.results import DRY_RUN, REWRITES_BESIDE_INPUT, WRITES_NEW_FOLDER, failure

THRESHOLDS = """`alpha` scales the N shift against H in the CSP (default 0.14). `intensity_from` picks
        height or volume for the intensity ratio. `bootstrap` sets the bootstrap iterations for
        errors. `csp_sigma_multiple`, `kd_outlier_z`, `noise_quantile`, `dd_runaway_ratio` and
        `ref_max_ratio` tune which residues enter the shared fit and which are flagged; leave
        them unset for the defaults. `prefix` names the output files."""


KD_SURVEY_DESCRIPTION = inspect.cleandoc("""Survey a titration before fitting a Kd: rates every residue ok / check /
        unusable and proposes a residue selection. `input` is series_analysis_tidy.csv.
        `p0` is the protein concentration and `conc_units` says whether the ligand
        concentrations are absolute or equivalents; `conc` gives them per titration point
        when the point labels are not numbers. `intensity_scales` is one factor per titration
        point applied to every peak's height and volume (a scan-count correction, e.g. 0.5 for
        a point acquired with twice the scans); it never changes positions. Figures and the
        survey table go to `out` (new folder); the residue list and a settings file are
        written beside `input`, or the residue list to `selection_file` (a file inside the data
        folder). """ + THRESHOLDS + " " + DRY_RUN)

KD_FIT_DESCRIPTION = inspect.cleandoc("""Start a Kd fit in the background; a quick one returns finished, a slow one
        returns a job to poll with `job_status`. Fit all residues, or only
        those named in `residues` or listed in `residues_file` (the selection file
        `kd_survey` writes beside its input) -- not both. `observables` defaults to both.
        The finished job reports each shared Kd with its `reliable` flag and how many
        residues are in the shared fit. `intensity_scales` is the same scan-count correction
        as in `kd_survey` and must match the survey's. """ + THRESHOLDS + """ `out` is a new folder. """ + DRY_RUN)


EXPORT_KD_DESCRIPTION = inspect.cleandoc("""
    Figures and tables from a Kd fit JSON (`…_kd_fit_data.json`). `kinds` picks what to render
    (default: all four) and `observables` which observable(s). `fig_format` is pdf, png or both,
    `per_page` the panels per page, `prefix` names the files, and `summary_only` writes the
    summary table without figures. `out` is a new folder.""") + " " + DRY_RUN


def register(server, ctx):
    root = ctx.root

    @server.tool(annotations=REWRITES_BESIDE_INPUT, description=KD_SURVEY_DESCRIPTION)
    async def kd_survey(input: str, p0: float, conc_units: Literal['absolute', 'equivalents'],
                        out: str, conc: list[float] | None = None,
                        intensity_scales: list[float] | None = None,
                        selection_file: str | None = None,
                        observables: list[Literal['csp', 'intensity']] | None = None,
                        prefix: str | None = None, alpha: float | None = None,
                        intensity_from: Literal['height', 'volume'] | None = None,
                        bootstrap: int | None = None, csp_sigma_multiple: float | None = None,
                        kd_outlier_z: float | None = None, noise_quantile: float | None = None,
                        dd_runaway_ratio: float | None = None,
                        ref_max_ratio: float | None = None,
                        dry_run: bool = False) -> CallToolResult:
        try:
            input_file = root.existing(input)
            root.refuse_escaping_links(input_file.parent)
            out_folder = root.new_output(out)
            selection = root.writable_file(selection_file) if selection_file is not None else None
            name = safe_name(prefix, 'prefix') if prefix is not None else None
        except REFUSALS as error:
            return failure(error)
        return await ctx.run(cli_call.kd_survey_args(
            input_file, p0, conc_units, out_folder, conc=conc, intensity_scales=intensity_scales,
            selection_file=selection, observables=observables, prefix=name, alpha=alpha,
            intensity_from=intensity_from, bootstrap=bootstrap,
            csp_sigma_multiple=csp_sigma_multiple, kd_outlier_z=kd_outlier_z,
            noise_quantile=noise_quantile, dd_runaway_ratio=dd_runaway_ratio,
            ref_max_ratio=ref_max_ratio, dry_run=dry_run))

    @server.tool(annotations=WRITES_NEW_FOLDER, description=EXPORT_KD_DESCRIPTION)
    async def export_kd(fit_json: str, out: str,
                        kinds: list[Literal['curves', 'ref-bars', 'kd-bars', 'global-fit']] | None = None,
                        observables: list[Literal['csp', 'intensity']] | None = None,
                        fig_format: list[Literal['pdf', 'png']] | None = None,
                        per_page: int | None = None, prefix: str | None = None,
                        summary_only: bool = False, dry_run: bool = False) -> CallToolResult:
        try:
            fit_file = root.existing(fit_json)
            out_folder = root.new_output(out)
            name = safe_name(prefix, 'prefix') if prefix is not None else None
            if per_page is not None and per_page < 1:
                raise ValueError('per_page must be at least 1')
        except (*REFUSALS, ValueError) as error:
            return failure(error)
        return await ctx.run(cli_call.export_kd_args(
            fit_file, out_folder, kinds=kinds, observables=observables, fig_format=fig_format,
            per_page=per_page, prefix=name, summary_only=summary_only, dry_run=dry_run))

    @server.tool(annotations=WRITES_NEW_FOLDER, description=KD_FIT_DESCRIPTION)
    async def kd_fit(input: str, p0: float, conc_units: Literal['absolute', 'equivalents'],
                     out: str, conc: list[float] | None = None,
                     residues: list[str] | None = None, residues_file: str | None = None,
                     observables: list[Literal['csp', 'intensity']] | None = None,
                     intensity_scales: list[float] | None = None,
                     prefix: str | None = None, alpha: float | None = None,
                     intensity_from: Literal['height', 'volume'] | None = None,
                     bootstrap: int | None = None, csp_sigma_multiple: float | None = None,
                     kd_outlier_z: float | None = None, noise_quantile: float | None = None,
                     dd_runaway_ratio: float | None = None,
                     ref_max_ratio: float | None = None,
                     dry_run: bool = False) -> CallToolResult:
        if residues is not None and residues_file is not None:
            return failure(ValueError('give residues or residues_file, not both'))
        try:
            input_file = root.existing(input)
            root.refuse_escaping_links(input_file.parent)
            selection = root.existing(residues_file) if residues_file is not None else None
            named = cli_call.selection_text(residues) if residues is not None else None
            name = safe_name(prefix, 'prefix') if prefix is not None else None
            out_folder = root.new_output(out)
            if named is not None:
                selection = out_folder / JOB_FOLDER / 'residues.txt'
            if dry_run and named is not None:
                selection = None
            return await ctx.launch('kd', cli_call.kd_fit_args(
                input_file, p0, conc_units, out_folder, conc=conc, residues_file=selection,
                observables=observables, intensity_scales=intensity_scales, prefix=name,
                alpha=alpha, intensity_from=intensity_from, bootstrap=bootstrap,
                csp_sigma_multiple=csp_sigma_multiple, kd_outlier_z=kd_outlier_z,
                noise_quantile=noise_quantile, dd_runaway_ratio=dd_runaway_ratio,
                ref_max_ratio=ref_max_ratio, dry_run=dry_run), out_folder, dry_run,
                files=None if named is None else {'residues.txt': named})
        except REFUSALS as error:
            return failure(error)
