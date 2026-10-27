# ABOUTME: Keeps the MCP server and the CLI in step: every subcommand and every flag is exposed or excluded on purpose.
# ABOUTME: A flag added to the CLI later fails here until a tool carries it or this file says why it does not.

import argparse
import inspect

import pytest

pytest.importorskip("mcp")

from lunaNMR.cli import build_parser
from lunaNMR.mcp import cli_call

# The subcommands no MCP tool runs, and why.
NOT_EXPOSED_SUBCOMMANDS = {
    "project remove": "deletes bundle contents irreversibly; the Project Browser does it with confirmation",
    "batch": "an unmaintained pass-through to another parser with no flags of its own",
}

# The flags no MCP tool sets, and why.
NOT_EXPOSED_FLAGS = {
    ("dynamixs t1t2", "--no-json"): "the JSON fit file is what the server reads results from",
    ("project export", "--force"): "overwrites an existing folder; exports go to a new one",
}

# The argument builders that run each subcommand (a tool may use two: plan and start).
BUILDERS = {
    "diagnose": [cli_call.diagnose_args],
    "series": [cli_call.plan_series_args, cli_call.start_series_args],
    "kd": [cli_call.kd_survey_args, cli_call.kd_fit_args],
    "export kd": [cli_call.export_kd_args],
    "integrate-1d": [cli_call.integrate_1d_args],
    "peaks shift": [cli_call.peaks_shift_args],
    "dynamixs t1t2": [cli_call.dynamixs_t1t2_args],
    "dynamixs methyl-t2": [cli_call.dynamixs_methyl_t2_args],
    "dynamixs hetnoe": [cli_call.dynamixs_hetnoe_args],
    "dynamixs t1rho": [cli_call.dynamixs_t1rho_args],
    "dynamixs density": [cli_call.dynamixs_density_args],
    "dynamixs modelfree": [cli_call.dynamixs_modelfree_args],
    "project build": [cli_call.project_build_args],
    "project inventory": [cli_call.project_inventory_args],
    "project export": [cli_call.project_export_args],
}
# Parameters that switch a default-on CLI behaviour off: the flag appears when they are False.
SWITCH_OFF = {cli_call.dynamixs_density_args: {"use_087", "plot", "parallel"}}


def _subparsers(parser):
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action.choices
    return {}


def _all_subcommands():
    found = {}
    for name, parser in _subparsers(build_parser()).items():
        inner = _subparsers(parser)
        if inner:
            found.update({f"{name} {sub}": subparser for sub, subparser in inner.items()})
        else:
            found[name] = parser
    return found


def _flags_of(parser):
    return {option for action in parser._actions for option in action.option_strings
            if option.startswith("--") and option != "--help"}


def _flags_a_builder_can_emit(builder):
    """Call the builder with every optional parameter set and collect the flags it emits."""
    switch_off = SWITCH_OFF.get(builder, set())
    lists = {"conc", "intensity_scales", "peaks", "kinds", "observables", "fig_format"}
    arguments = {}
    for name, parameter in inspect.signature(builder).parameters.items():
        if name in switch_off:
            arguments[name] = False
        elif parameter.default is False:
            arguments[name] = True
        else:
            arguments[name] = ["x"] if name in lists else "x"
    argv = builder(**arguments)
    return {token for token in argv if isinstance(token, str) and token.startswith("--")}


def test_every_subcommand_is_exposed_or_excluded_on_purpose():
    subcommands = set(_all_subcommands())
    accounted = set(BUILDERS) | set(NOT_EXPOSED_SUBCOMMANDS)
    assert subcommands - accounted == set(), "a CLI subcommand has no MCP tool and no stated reason"
    assert accounted - subcommands == set(), "this file names a subcommand the CLI no longer has"


@pytest.mark.parametrize("subcommand", sorted(BUILDERS))
def test_every_flag_of_an_exposed_subcommand_is_carried_or_excluded_on_purpose(subcommand):
    parser_flags = _flags_of(_all_subcommands()[subcommand])
    carried = set().union(*(_flags_a_builder_can_emit(b) for b in BUILDERS[subcommand]))
    excluded = {flag for (name, flag) in NOT_EXPOSED_FLAGS if name == subcommand}
    missing = parser_flags - carried - excluded
    assert missing == set(), f"{subcommand}: no MCP tool can set {sorted(missing)}"
    assert not (carried & excluded), f"{subcommand}: {sorted(carried & excluded)} is both carried and excluded"
    # --dry-run and --format are the server's to set; every other emitted flag must be real
    assert carried - parser_flags - {"--format", "--dry-run"} == set()
