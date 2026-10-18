# ABOUTME: `project build` finds a run's outputs through symlinked folders under run_dir,
# ABOUTME: and a symlink that loops back into the run does not make it walk forever.

"""
A run folder assembled by linking result folders in, rather than copying them, gave a
bundle with nothing in it: the search did not descend into a symlinked folder. Real
DNAJA1 titration run, linked in, built, then read back with the real loader.
"""

import types

import pytest


@pytest.fixture
def linked_run(tmp_path, dnaja1_titration_series):
    """A run folder holding only links: the series folder, and one back to the run."""
    run = tmp_path / "run"
    run.mkdir()
    (run / "series").symlink_to(dnaja1_titration_series.parent, target_is_directory=True)
    (run / "loop").symlink_to(run, target_is_directory=True)
    return run


def test_a_symlinked_series_folder_is_bundled(linked_run, tmp_path):
    from lunaNMR.cli import main
    from lunaNMR.utils.project_manager import ProjectManager
    bundle = tmp_path / "linked.lunaNMR"
    assert main(["project", "build", str(linked_run), "--out", str(bundle)]) == 0
    window = types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=None))
    assert ProjectManager(window).load_project(bundle)
    assert len(window.saved_series) == 1
    (series,) = window.saved_series.values()
    assert len(series.results) == 3


def test_a_folder_reached_twice_is_found_once(linked_run):
    from lunaNMR.utils.project_builder import find_series
    assert len(find_series(linked_run)) == 1
