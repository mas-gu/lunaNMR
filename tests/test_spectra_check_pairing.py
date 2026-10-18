# ABOUTME: `diagnose` checks each spectrum against the peak list of the same name.
# ABOUTME: A folder holding several states, each with its own list, must not borrow the first one.

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

BINDING = Path(__file__).resolve().parents[2] / "data_test" / "data_binding"
FIRST = BINDING / "ERDJ6_HSPA8"
SECOND = BINDING / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not (FIRST.exists() and SECOND.exists()),
                                reason="data_test not present")


def _link(folder, name, target):
    (folder / name).symlink_to(target)


@pytest.fixture(scope="module")
def mixed(tmp_path_factory):
    """Two unrelated proteins in one folder, each spectrum beside its own list, plus a
    spectrum with no list of its name. Sorted, `first.txt` comes before `second.txt`, so
    a folder-wide list would judge `second.ft` against the other protein's peaks."""
    from lunaNMR.validation.spectra_check import check_dataset
    root = tmp_path_factory.mktemp("root")
    d = root / "states"
    d.mkdir()
    _link(d, "first.ft", FIRST / "2D_ERDJ6_0o0.ft")
    _link(d, "first.txt", FIRST / "ERDJ6.txt")
    _link(d, "second.ft", SECOND / "2D_DNAJA1_0o0.ft")
    _link(d, "second.txt", SECOND / "DNAJA1.txt")
    _link(d, "unpaired.ft", SECOND / "2D_DNAJA1_HSPA8_3o0.ft")
    return check_dataset(root, quick=True)


def _spectra_by_list(result):
    return {Path(e["peak_list"]["path"]).name: [s["spectrum"] for s in e["spectra"]]
            for e in result["experiments"]}


class TestSeveralListsInOneFolder:
    def test_each_spectrum_is_checked_against_its_own_list(self, mixed):
        assert _spectra_by_list(mixed) == {"first.txt": ["first.ft"],
                                           "second.txt": ["second.ft"]}

    def test_paired_spectra_register(self, mixed):
        for experiment in mixed["experiments"]:
            for spec in experiment["spectra"]:
                assert spec["capture_rate"] > 0.85, spec
        assert not [f for f in mixed["findings"] if f.check == "registration"]

    def test_a_spectrum_without_its_own_list_is_a_warning_naming_it_and_the_lists(self, mixed):
        warns = [f for f in mixed["findings"] if f.check == "peak_list"]
        assert len(warns) == 1
        assert warns[0].severity == "WARN"
        for name in ("unpaired.ft", "first.txt", "second.txt"):
            assert name in warns[0].message

    def test_an_unpaired_spectrum_is_not_checked_against_some_list(self, mixed):
        checked = [s for spectra in _spectra_by_list(mixed).values() for s in spectra]
        assert "unpaired.ft" not in checked


class TestASingleListServesTheFolder:
    def test_every_spectrum_uses_it_whatever_its_name(self, tmp_path):
        from lunaNMR.validation.spectra_check import check_dataset
        d = tmp_path / "titration"
        d.mkdir()
        _link(d, "point_a.ft", SECOND / "2D_DNAJA1_0o0.ft")
        _link(d, "point_b.ft", SECOND / "2D_DNAJA1_HSPA8_3o0.ft")
        _link(d, "assignment.txt", SECOND / "DNAJA1.txt")
        result = check_dataset(tmp_path, quick=True)
        assert _spectra_by_list(result) == {"assignment.txt": ["point_a.ft", "point_b.ft"]}
        assert not [f for f in result["findings"] if f.check == "peak_list"]
