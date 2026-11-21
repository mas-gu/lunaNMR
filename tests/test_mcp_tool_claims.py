# ABOUTME: What an MCP tool tells the model is what its CLI command does: no parameter the command
# ABOUTME: ignores, no file it does not write, and delay units described as the CLI's own help does.

import argparse

import pytest

pytest.importorskip("mcp")

from lunaNMR.cli import build_parser
from mcp_helpers import list_tools


@pytest.fixture(scope="module")
def tools(tmp_path_factory):
    return {tool.name: tool for tool in list_tools(tmp_path_factory.mktemp("root"))}


def _cli_help(*subcommand, flag):
    """The help text the CLI gives `flag` of a subcommand."""
    parser = build_parser()
    for name in subcommand:
        choices = next(action.choices for action in parser._actions
                       if isinstance(action, argparse._SubParsersAction))
        parser = choices[name]
    return next(action.help for action in parser._actions if flag in action.option_strings)


def test_the_survey_offers_no_parameter_a_survey_does_not_read(tools):
    """`kd --survey` fits nothing: the observables to fit and the bootstrap count of the
    fit's errors reach no part of it."""
    survey = tools["kd_survey"]
    assert not {"observables", "bootstrap"} & set(survey.input_schema["properties"])
    assert "bootstrap" not in survey.description


def test_the_fit_keeps_the_parameters_a_fit_reads(tools):
    fit = tools["kd_fit"]
    assert {"observables", "bootstrap"} <= set(fit.input_schema["properties"])
    assert "`bootstrap`" in fit.description


@pytest.mark.parametrize("tool, subcommand", [
    ("dynamixs_t1t2", "t1t2"),
    ("dynamixs_methyl_t2", "methyl-t2"),
    ("dynamixs_t1rho", "t1rho"),
])
def test_time_units_are_described_as_the_cli_describes_them(tools, tool, subcommand):
    """t1t2 and methyl-t2 report time constants in the units given: the units label the
    output. t1rho converts to rates, which the units rescale."""
    rescales = "RESCALES" in _cli_help("dynamixs", subcommand, flag="--time-units")
    description = tools[tool].description
    assert ("RESCALES" in description) is rescales
    assert ("does not rescale" in description) is not rescales


def test_a_project_export_promises_no_file_it_does_not_write(tools):
    assert "manifest" not in tools["project_export"].description
