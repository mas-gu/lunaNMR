# ABOUTME: MCP tools for peak lists and 1D spectra: shift a list into register, integrate a 1D series.
# ABOUTME: A shifted list is written to a new file only; a prefix is a name, never a path.

from typing import Literal

from mcp.types import CallToolResult

from lunaNMR.mcp import cli_call
from lunaNMR.mcp.context import REFUSALS
from lunaNMR.mcp.names import safe_name
from lunaNMR.mcp.results import WRITES_NEW_FOLDER, failure


def register(server, ctx):
    root = ctx.root

    @server.tool(annotations=WRITES_NEW_FOLDER)
    async def peaks_shift(peaks: str, out: str, dx: float | None = None, dy: float | None = None,
                          auto: bool = False, spectrum: str | None = None,
                          quick: bool = False) -> CallToolResult:
        """Apply a rigid ppm offset to a peak list and write it to a new file `out` (the
        original is never changed). Give `dx` (1H) and `dy` (15N) in ppm, or set `auto` and
        name a `spectrum` to measure the offset against it. An automatic shift that does not
        register the list is refused and writes nothing. `quick` uses a coarser grid, which
        quantises the 15N offset to 0.03 ppm and misreports smaller real offsets."""
        try:
            peak_list = root.existing(peaks)
            out_file = root.new_file(out)
            spectrum_path = root.existing(spectrum) if spectrum is not None else None
        except REFUSALS as error:
            return failure(error)
        return await ctx.run(cli_call.peaks_shift_args(
            peak_list, out_file, dx=dx, dy=dy, auto=auto, spectrum=spectrum_path, quick=quick))

    @server.tool(annotations=WRITES_NEW_FOLDER)
    async def integrate_1d(spectra: str, out: str, peaks: list[float] | None = None,
                           detect: bool = False, prefix: str | None = None,
                           min_snr: float | None = None,
                           value: Literal['height', 'area', 'fit_area'] | None = None
                           ) -> CallToolResult:
        """Measure peaks across a series of 1D spectra in the folder `spectra`: give the peak
        positions in ppm as `peaks`, or set `detect` to pick them in the first spectrum above
        `min_snr` times the noise. `value` chooses the observable of the wide table: `area`
        (the region sum, default) stays proportional when linewidths change; `fit_area` is
        diagnostic only, because it extrapolates tails the data never sampled. `prefix` names
        the output files. `out` is a new folder. The result lists the files written."""
        try:
            spectra_folder = root.existing(spectra)
            root.refuse_escaping_links(spectra_folder)
            out_folder = root.new_output(out)
            name = safe_name(prefix, 'prefix') if prefix is not None else None
            return await ctx.start_job('integrate-1d', cli_call.integrate_1d_args(
                spectra_folder, out_folder, peaks=peaks, detect=detect, prefix=name,
                min_snr=min_snr, value=value), out_folder)
        except REFUSALS as error:
            return failure(error)
