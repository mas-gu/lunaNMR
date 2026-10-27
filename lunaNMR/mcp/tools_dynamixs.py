# ABOUTME: MCP tools for relaxation and dynamics: T1/T2, methyl T2, hetNOE, T1rho, spectral density, model-free.
# ABOUTME: Each is a job that waits briefly; input files are guarded paths and names never become paths.

import inspect
from typing import Literal

from mcp.types import CallToolResult

from lunaNMR.mcp import cli_call
from lunaNMR.mcp.context import REFUSALS
from lunaNMR.mcp.names import safe_name
from lunaNMR.mcp.results import WRITES_NEW_FOLDER, describe, failure

TimeUnits = Literal['ms', 's', 'us']
ErrorMethod = Literal['analytical', 'bootstrap']

DELAYS = """Delay values are read from the matrix column headers or its series_metadata.json;
`time_units` states them when neither does, and RESCALES the fitted rates, so a wrong unit gives
a confidently wrong number. `allow_unparsed_delays` fits anyway when some spectra carried no
delay (those points are dropped)."""


def register(server, ctx):
    root = ctx.root

    def analysis_tool(function):
        """Register `function` as a tool whose docstring has the shared delay text filled in."""
        function.__doc__ = function.__doc__.replace('@DELAYS@', DELAYS)
        return server.tool(annotations=WRITES_NEW_FOLDER)(describe(function))

    def _optional(path):
        return root.existing(path) if path is not None else None

    def _name(value, what):
        return safe_name(value, what) if value is not None else None

    @analysis_tool
    async def dynamixs_t1t2(input: str, out: str, exp: Literal['T1', 'T2'],
                            prefix: str | None = None, field_name: str | None = None,
                            field_freq: float | None = None, time_units: TimeUnits | None = None,
                            allow_unparsed_delays: bool = False,
                            error_method: ErrorMethod | None = None,
                            bootstrap: int | None = None, dry_run: bool = False) -> CallToolResult:
        """Fit T1 or T2 (`exp`) to a series matrix (`peak_intensity_matrix.csv` from a series
        run; a tidy CSV is refused) with a mono-exponential. The finished job reports how many
        residues fail the usual gate (fitted, finite T, error below the value) and names them.
        In the result, `mean_t2` is the mean time constant of whichever experiment was fitted.
        `field_freq` is the 1H frequency in MHz, `error_method` analytical or bootstrap with
        `bootstrap` iterations. @DELAYS@ `out` is a new folder. @DRY_RUN@"""
        try:
            matrix = root.existing(input)
            out_folder = root.new_output(out)
            names = (_name(prefix, 'prefix'), _name(field_name, 'field_name'))
            return await ctx.launch('t1t2', cli_call.dynamixs_t1t2_args(
                matrix, out_folder, exp, prefix=names[0], field_name=names[1],
                field_freq=field_freq, time_units=time_units,
                allow_unparsed_delays=allow_unparsed_delays, error_method=error_method,
                bootstrap=bootstrap, dry_run=dry_run), out_folder, dry_run)
        except REFUSALS as error:
            return failure(error)

    @analysis_tool
    async def dynamixs_methyl_t2(input: str, out: str, prefix: str | None = None,
                                 field_name: str | None = None, field_freq: float | None = None,
                                 time_units: TimeUnits | None = None,
                                 allow_unparsed_delays: bool = False,
                                 error_method: ErrorMethod | None = None,
                                 bootstrap: int | None = None, dry_run: bool = False) -> CallToolResult:
        """Fit the bi-exponential methyl T2 (Tugarinov-Kay shared-amplitude form) to a series
        matrix of 13CH3 intensities. T2a, the slow component, is the methyl T2 to report. The
        result counts the fits whose bi-exponential could not be identified (`n_bi_exp_unidentifiable`):
        those residues collapsed to a single exponential, so read it before trusting the
        slow/fast split. @DELAYS@ `out` is a new folder. @DRY_RUN@"""
        try:
            matrix = root.existing(input)
            out_folder = root.new_output(out)
            names = (_name(prefix, 'prefix'), _name(field_name, 'field_name'))
            return await ctx.launch('methyl-t2', cli_call.dynamixs_methyl_t2_args(
                matrix, out_folder, prefix=names[0], field_name=names[1], field_freq=field_freq,
                time_units=time_units, allow_unparsed_delays=allow_unparsed_delays,
                error_method=error_method, bootstrap=bootstrap, dry_run=dry_run), out_folder, dry_run)
        except REFUSALS as error:
            return failure(error)

    @analysis_tool
    async def dynamixs_hetnoe(sat: str, unsat: str, out: str,
                              prefix: str | None = None, dry_run: bool = False) -> CallToolResult:
        """Compute the heteronuclear NOE, I_sat / I_unsat per residue, from two headerless files
        `residue,intensity[,error]`. Decide which plane is saturated by the measured intensity
        (the saturated one is lower, a ratio of 0.75-0.85 for a rigid fold), never by filename.
        Without an error column a flat error is assumed, which is too tight. The result gives
        the median hetNOE and how many residues fall outside the possible range. `out` is a new
        folder. @DRY_RUN@"""
        try:
            sat_file = root.existing(sat)
            unsat_file = root.existing(unsat)
            out_folder = root.new_output(out)
            name = _name(prefix, 'prefix')
            return await ctx.launch('hetnoe', cli_call.dynamixs_hetnoe_args(
                sat_file, unsat_file, out_folder, prefix=name, dry_run=dry_run), out_folder, dry_run)
        except REFUSALS as error:
            return failure(error)

    @analysis_tool
    async def dynamixs_t1rho(input: str, t1: str, peaks: str, omega1: float, carrier: float,
                             theta: float, out: str, prefix: str | None = None,
                             field_freq: float | None = None, time_units: TimeUnits | None = None,
                             allow_unparsed_delays: bool = False,
                             error_method: ErrorMethod | None = None,
                             bootstrap: int | None = None, dry_run: bool = False) -> CallToolResult:
        """Fit a T1rho series and convert it to R2 per residue, correcting the spin-lock tilt
        angle from each residue's offset (R1rho = R1 cos^2 theta + R2 sin^2 theta). `input` is
        the T1rho matrix, `t1` the T1 matrix (R1 is needed for the tilt), `peaks` the peak list
        (each residue's 15N shift). `omega1` is the spin-lock field in Hz, `carrier` the 15N
        carrier in ppm and `theta` the nominal tilt angle in degrees: these are acquisition
        parameters recoverable from neither the spectra nor the filenames, so ask for them
        rather than guessing. The R2 table it writes feeds model-free as `f1_r2_table`.
        @DELAYS@ `out` is a new folder. @DRY_RUN@"""
        try:
            rho = root.existing(input)
            t1_matrix = root.existing(t1)
            peak_list = root.existing(peaks)
            out_folder = root.new_output(out)
            name = _name(prefix, 'prefix')
            return await ctx.launch('t1rho', cli_call.dynamixs_t1rho_args(
                rho, t1_matrix, peak_list, omega1, carrier, theta, out_folder, prefix=name,
                field_freq=field_freq, time_units=time_units,
                allow_unparsed_delays=allow_unparsed_delays, error_method=error_method,
                bootstrap=bootstrap, dry_run=dry_run), out_folder, dry_run)
        except REFUSALS as error:
            return failure(error)

    @analysis_tool
    async def dynamixs_density(input: str, out: str, input2: str | None = None,
                               prefix: str | None = None, dual: bool = False,
                               use_087: bool = True, monte_carlo: bool = False,
                               n_samples: int | None = None, plot: bool = True,
                               parallel: bool = True, field1_freq: float | None = None,
                               field2_freq: float | None = None, rnh: float | None = None,
                               csa: float | None = None, dry_run: bool = False) -> CallToolResult:
        """Reduced spectral density mapping, J(0), J(wN), J(0.87 wH), from a table with columns
        Residue,R1,R1err,R2,R2err,hetNOE,hetNOEerr (a model-free run writes one as
        `<prefix>_field1_input.csv`). For two fields give `input2` and set `dual`.
        `use_087` uses the 0.87 wH approximation (default), `monte_carlo` adds error samples
        (`n_samples`), `plot` and `parallel` can be switched off to save time. `field1_freq` and
        `field2_freq` are 1H frequencies in MHz, `rnh` the N-H bond length in Angstrom and `csa`
        the 15N CSA in ppm. `out` is a new folder. @DRY_RUN@"""
        try:
            table = root.existing(input)
            table2 = _optional(input2)
            out_folder = root.new_output(out)
            name = _name(prefix, 'prefix')
            return await ctx.launch('density', cli_call.dynamixs_density_args(
                table, out_folder, input2=table2, prefix=name, dual=dual, use_087=use_087,
                monte_carlo=monte_carlo, n_samples=n_samples, plot=plot, parallel=parallel,
                field1_freq=field1_freq, field2_freq=field2_freq, rnh=rnh, csa=csa, dry_run=dry_run), out_folder, dry_run)
        except REFUSALS as error:
            return failure(error)

    @analysis_tool
    async def dynamixs_modelfree(
            out: str, f1_t1: str, f1_noe_sat: str, f1_noe_unsat: str,
            f1_t2: str | None = None, f1_r2_table: str | None = None,
            f2_t1: str | None = None, f2_t2: str | None = None, f2_r2_table: str | None = None,
            f2_noe_sat: str | None = None, f2_noe_unsat: str | None = None,
            f1_t1_units: TimeUnits | None = None, f1_t2_units: TimeUnits | None = None,
            f2_t1_units: TimeUnits | None = None, f2_t2_units: TimeUnits | None = None,
            allow_unparsed_delays: bool = False, prefix: str | None = None, dual: bool = False,
            method: Literal['single_jwh', 'single_087', 'dual_jwh', 'dual_087'] | None = None,
            init_amp: float | None = None, init_t1: float | None = None,
            init_t2: float | None = None, n_bootstrap: int | None = None,
            error_method: ErrorMethod | None = None, n_monte_carlo: int | None = None,
            field1_freq: float | None = None, field2_freq: float | None = None,
            rnh: float | None = None, csa: float | None = None, dry_run: bool = False) -> CallToolResult:
        """Model-free (Lipari-Szabo) analysis from T1, T2 and hetNOE: fits T1/T2, builds the
        spectral density and fits S2, tau_e and Rex per residue by bounded least squares.
        Field 1 needs `f1_t1`, `f1_noe_sat`, `f1_noe_unsat` and one of `f1_t2` or `f1_r2_table`
        (an R2 table from `dynamixs_t1rho`; not both). For two fields give the `f2_` inputs, set
        `dual` and both `field1_freq` and `field2_freq` (MHz). Run the single-field fits too:
        they are the control for whether to believe the dual one. Each series' delay units
        default to ms; give `f1_t1_units` etc. when a series is in seconds, because a wrong unit
        silently rescales every rate. The result reports the residues kept at each stage, any
        fits that hit the iteration cap, and the validation checks; a dual-field run takes about
        40 seconds, so expect a job. `out` is a new folder. @DRY_RUN@"""
        try:
            paths = {name: _optional(value) for name, value in (
                ('f1_t1', f1_t1), ('f1_t2', f1_t2), ('f1_r2_table', f1_r2_table),
                ('f1_noe_sat', f1_noe_sat), ('f1_noe_unsat', f1_noe_unsat), ('f2_t1', f2_t1),
                ('f2_t2', f2_t2), ('f2_r2_table', f2_r2_table), ('f2_noe_sat', f2_noe_sat),
                ('f2_noe_unsat', f2_noe_unsat))}
            out_folder = root.new_output(out)
            name = _name(prefix, 'prefix')
            return await ctx.launch('modelfree', cli_call.dynamixs_modelfree_args(
                out_folder, **paths, f1_t1_units=f1_t1_units, f1_t2_units=f1_t2_units,
                f2_t1_units=f2_t1_units, f2_t2_units=f2_t2_units,
                allow_unparsed_delays=allow_unparsed_delays, prefix=name, dual=dual,
                method=method, init_amp=init_amp, init_t1=init_t1, init_t2=init_t2,
                n_bootstrap=n_bootstrap, error_method=error_method,
                n_monte_carlo=n_monte_carlo, field1_freq=field1_freq,
                field2_freq=field2_freq, rnh=rnh, csa=csa, dry_run=dry_run), out_folder, dry_run)
        except REFUSALS as error:
            return failure(error)
