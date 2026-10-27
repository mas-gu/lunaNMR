# ABOUTME: The MCP server's data-root guard: every path an AI client names must stay inside one folder.
# ABOUTME: A new run may not land on an existing folder, so a client cannot overwrite earlier results.

import pytest

from lunaNMR.mcp.data_root import DataRoot, PathRefused


@pytest.fixture
def root(tmp_path):
    folder = tmp_path / "dataset"
    (folder / "series").mkdir(parents=True)
    return DataRoot(folder)


def test_existing_path_inside_the_root_resolves(root):
    assert root.existing("series") == (root.path / "series").resolve()


def test_absolute_path_inside_the_root_resolves(root):
    inside = root.path / "series"
    assert root.existing(str(inside)) == inside.resolve()


def test_dotdot_escape_is_refused(root):
    with pytest.raises(PathRefused, match="outside the data root"):
        root.existing("../elsewhere")


def test_absolute_path_outside_the_root_is_refused(root, tmp_path):
    with pytest.raises(PathRefused, match="outside the data root"):
        root.existing(str(tmp_path))


def test_symlink_pointing_outside_is_refused(root, tmp_path):
    outside = tmp_path / "outside"
    outside.mkdir()
    (root.path / "link").symlink_to(outside)
    with pytest.raises(PathRefused, match="outside the data root"):
        root.existing("link")


def test_missing_path_is_refused_with_its_name(root):
    with pytest.raises(PathRefused, match="does not exist: nothing_here"):
        root.existing("nothing_here")


def test_new_output_inside_the_root_resolves(root):
    assert root.new_output("run_1") == (root.path / "run_1").resolve()


def test_new_output_may_not_be_an_existing_folder(root):
    with pytest.raises(PathRefused, match="already exists"):
        root.new_output("series")


def test_new_output_outside_the_root_is_refused(root, tmp_path):
    with pytest.raises(PathRefused, match="outside the data root"):
        root.new_output(str(tmp_path / "run_1"))


def test_the_root_must_be_an_existing_folder(tmp_path):
    with pytest.raises(PathRefused, match="data root"):
        DataRoot(tmp_path / "missing")


def test_a_new_file_inside_the_root_resolves(root):
    assert root.new_file("shifted.txt") == (root.path / "shifted.txt").resolve()


def test_a_new_file_may_not_overwrite_an_existing_one(root):
    (root.path / "taken.txt").write_text("x")
    with pytest.raises(PathRefused, match="already exists"):
        root.new_file("taken.txt")


def test_a_new_file_needs_an_existing_folder_inside_the_root(root, tmp_path):
    with pytest.raises(PathRefused, match="folder does not exist"):
        root.new_file("no_such_folder/out.txt")
    with pytest.raises(PathRefused, match="outside the data root"):
        root.new_file(str(tmp_path / "out.txt"))
