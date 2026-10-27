# ABOUTME: MCP tools for .lunaNMR project bundles: build one from a run, list its contents, export it.
# ABOUTME: Removing from a bundle and overwriting an export are not offered; both destroy data.

from mcp.types import CallToolResult

from lunaNMR.mcp import cli_call
from lunaNMR.mcp.context import REFUSALS
from lunaNMR.mcp.results import READS_ONLY, WRITES_NEW_FOLDER, describe, failure


def register(server, ctx):
    root = ctx.root

    @server.tool(annotations=WRITES_NEW_FOLDER)
    @describe
    async def project_build(run_dir: str, out: str, spectrum: str, peaks: str,
                            dry_run: bool = False) -> CallToolResult:
        """Assemble a `.lunaNMR` project bundle `out` (a new folder, name it `*.lunaNMR`) from
        a run folder holding series, Kd and relaxation outputs, found anywhere beneath
        `run_dir`. `spectrum` and `peaks` name the spectrum and peak list the bundle records for
        viewing; they are required, and must be inside the data folder, so the bundle never
        takes them from paths written in a run's own metadata. @DRY_RUN@"""
        try:
            run_folder = root.existing(run_dir)
            root.refuse_escaping_links(run_folder)
            spectrum_file = root.existing(spectrum)
            peak_list = root.existing(peaks)
            out_folder = root.new_output(out)
        except REFUSALS as error:
            return failure(error)
        return await ctx.run(cli_call.project_build_args(
            run_folder, out_folder, spectrum_file, peak_list, dry_run=dry_run))

    @server.tool(annotations=READS_ONLY)
    async def project_inventory(bundle: str) -> CallToolResult:
        """List what a `.lunaNMR` project bundle holds, by category (series, Kd analyses, DynamiXs
        runs, fit results) with the size of each item. Read-only: nothing in the bundle changes."""
        try:
            bundle_folder = root.existing(bundle)
            root.refuse_escaping_links(bundle_folder)
        except REFUSALS as error:
            return failure(error)
        return await ctx.run(cli_call.project_inventory_args(bundle_folder))

    @server.tool(annotations=WRITES_NEW_FOLDER)
    @describe
    async def project_export(bundle: str, out: str, dry_run: bool = False) -> CallToolResult:
        """Copy a bundle's contents to a plain folder `out` (new) with a manifest. The bundle
        itself is not changed. @DRY_RUN@"""
        try:
            bundle_folder = root.existing(bundle)
            root.refuse_escaping_links(bundle_folder)
            out_folder = root.new_output(out)
        except REFUSALS as error:
            return failure(error)
        return await ctx.run(cli_call.project_export_args(
            bundle_folder, out_folder, dry_run=dry_run))
