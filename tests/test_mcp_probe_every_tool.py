# ABOUTME: Feeds a path outside the data root into every string parameter of every MCP tool.
# ABOUTME: Each must be refused at once, naming that value, with nothing outside read, written or leaked.

import json

import pytest

pytest.importorskip("mcp")

from mcp_helpers import (BINDING, call_result, list_tools, relaxation_data_present,
                         stage_binding_dataset, stage_relaxation_dataset)

pytestmark = pytest.mark.skipif(
    not (BINDING.is_dir() and relaxation_data_present()),
    reason="the binding and relaxation test data are both needed")

REFUSALS = {"PathRefused", "InvalidName", "InvalidResidueName"}
SECRET = "SECRET_MARKER_DO_NOT_LEAK"

# A valid value for each required string parameter, so a refusal can only come from the hostile one.
VALID = {
    "path": "series", "out": "new_out", "spectra": "series", "peaks": "series/DNAJA1.txt",
    "input": "tidy.csv", "fit_json": "tidy.csv", "t1": "tidy.csv", "sat": "tidy.csv",
    "unsat": "tidy.csv", "f1_t1": "tidy.csv", "f1_noe_sat": "tidy.csv", "f1_noe_unsat": "tidy.csv",
    "run_dir": "series", "spectrum": "series/2D_DNAJA1_0o0.ft", "bundle": "series",
    "job": "series"}


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    base = tmp_path_factory.mktemp("probe")
    root = base / "root"
    root.mkdir()
    stage_binding_dataset(root)
    stage_relaxation_dataset(root)
    outside = base / "outside"
    outside.mkdir()
    (outside / "secret.txt").write_text(SECRET)
    return root, outside


def _tree(folder):
    return sorted((str(path.relative_to(folder)), path.stat().st_size)
                  for path in folder.rglob("*"))


def _string_params(tool):
    """(name, is_list, required) for every free-text parameter: not an enum, not a number."""
    schema = tool.input_schema
    found = []
    for name, spec in schema.get("properties", {}).items():
        branches = spec.get("anyOf", [spec])
        for branch in branches:
            if branch.get("type") == "string" and "enum" not in branch:
                found.append((name, False, name in schema.get("required", [])))
            elif branch.get("type") == "array":
                items = branch.get("items", {})
                if items.get("type") == "string" and "enum" not in items:
                    found.append((name, True, name in schema.get("required", [])))
    return found


def _valid_value(name, spec):
    """Some valid value for a required parameter of any type, so only the hostile one can fail."""
    branch = next(b for b in spec.get("anyOf", [spec]) if b.get("type") != "null")
    if "enum" in branch:
        return branch["enum"][0]
    kind = branch.get("type")
    if kind == "string":
        assert name in VALID, f"give {name} a valid value in VALID"
        return VALID[name]
    return {"number": 1, "integer": 1, "boolean": False, "array": []}[kind]


def _baseline(tool):
    schema = tool.input_schema
    return {name: _valid_value(name, schema["properties"][name])
            for name in schema.get("required", [])}


def _cases():
    import tempfile
    from pathlib import Path
    tools = list_tools(Path(tempfile.gettempdir()))
    return [(tool.name, param, is_list) for tool in tools
            for param, is_list, _ in _string_params(tool)]


def test_every_string_parameter_is_covered():
    cases = _cases()
    assert len(cases) > 60, "the probe found suspiciously few string parameters"
    assert {name for name, _, _ in cases} >= {"diagnose", "kd_fit", "dynamixs_modelfree",
                                              "project_build", "job_status"}


@pytest.mark.parametrize("tool_name,param,is_list", _cases())
def test_a_hostile_value_is_refused_at_once_and_touches_nothing(world, tool_name, param, is_list):
    root, outside = world
    hostile = str(outside / "secret.txt")
    tool = next(t for t in list_tools(root) if t.name == tool_name)
    arguments = _baseline(tool)
    arguments[param] = [hostile] if is_list else hostile

    before_root, before_outside = _tree(root), _tree(outside)
    result = call_result(root, tool_name, arguments)
    text = result.content[0].text
    envelope = json.loads(text)

    assert result.is_error is True, text[:300]
    assert envelope["error"]["type"] in REFUSALS, text[:300]
    assert hostile[:30] in envelope["error"]["message"], text[:300]
    assert SECRET not in text
    assert _tree(root) == before_root and _tree(outside) == before_outside
