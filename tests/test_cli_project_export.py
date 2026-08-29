# ABOUTME: `project export` extracts a .lunaNMR bundle into a plain folder, headlessly.
# ABOUTME: It must copy everything, including files inventory does not classify.

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))


def _bundle(tmp_path):
    """Shaped like a real bundle, including a loose file under series_results/ —
    inventory enumerates only directories there, so that file is unclassified."""
    b = tmp_path / "proj.lunaNMR"
    (b / "fit_results" / "arrays").mkdir(parents=True)
    (b / "series_results" / "run_A").mkdir(parents=True)
    (b / "project.json").write_text('{"schema_version": "1.1"}')
    (b / "peak_list.csv").write_text("Assignment,Position_X,Position_Y\nK3,8.0,120.0\n")
    (b / "fit_results" / "fits.json").write_text("[]")
    (b / "fit_results" / "arrays" / "peak_000.npz").write_bytes(b"x" * 512)
    (b / "series_results" / "run_A" / "peak_intensity_matrix.csv").write_text("a,b\n1,2\n")
    (b / "series_results" / "batch_results.json").write_text("{}")   # loose: unclassified
    return b


class TestEverythingComesOut:
    """A partial export is worse than none: the caller cannot tell which files were
    left behind, and the bundle is the only copy."""

    def test_every_file_is_copied(self, tmp_path):
        from lunaNMR.cli import main
        b, out = _bundle(tmp_path), tmp_path / "exported"
        assert main(["project", "export", str(b), "--out", str(out)]) == 0
        want = {str(p.relative_to(b)) for p in b.rglob("*") if p.is_file()}
        got = {str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()}
        assert got == want

    def test_a_file_inventory_does_not_classify_is_still_copied(self, tmp_path):
        """inventory only walks directories under series_results/, so a loose file
        there is invisible to it. Driving the copy off inventory would drop it."""
        from lunaNMR.cli import main
        b, out = _bundle(tmp_path), tmp_path / "exported"
        main(["project", "export", str(b), "--out", str(out)])
        assert (out / "series_results" / "batch_results.json").exists()

    def test_a_file_no_category_covers_is_reported(self, tmp_path, capsys):
        """Rather than hidden: a file inventory cannot account for is copied AND named,
        so the gap stays visible instead of being silently absorbed. (The loose
        series_results file that first exposed this is now classified; a stray file at
        the bundle root still is not, and the guarantee is what matters.)"""
        from lunaNMR.cli import main
        b, out = _bundle(tmp_path), tmp_path / "exported"
        (b / "scratch_notes.txt").write_text("nothing claims this")
        main(["project", "export", str(b), "--out", str(out), "--format", "json"])
        payload = json.loads(capsys.readouterr().out)
        assert "scratch_notes.txt" in payload["unclassified"]
        assert (out / "scratch_notes.txt").exists()

    def test_content_survives(self, tmp_path):
        from lunaNMR.cli import main
        b, out = _bundle(tmp_path), tmp_path / "exported"
        main(["project", "export", str(b), "--out", str(out)])
        assert (out / "peak_list.csv").read_text().startswith("Assignment")
        assert (out / "fit_results" / "arrays" / "peak_000.npz").stat().st_size == 512


class TestItRefusesToDestroy:
    def test_an_existing_destination_is_refused(self, tmp_path):
        from lunaNMR.cli import main
        b, out = _bundle(tmp_path), tmp_path / "exported"
        out.mkdir()
        (out / "precious.txt").write_text("do not clobber")
        assert main(["project", "export", str(b), "--out", str(out)]) == 1
        assert (out / "precious.txt").read_text() == "do not clobber"

    def test_force_allows_it(self, tmp_path):
        from lunaNMR.cli import main
        b, out = _bundle(tmp_path), tmp_path / "exported"
        out.mkdir()
        assert main(["project", "export", str(b), "--out", str(out), "--force"]) == 0
        assert (out / "peak_list.csv").exists()

    def test_a_missing_bundle_is_an_error(self, tmp_path):
        from lunaNMR.cli import main
        assert main(["project", "export", str(tmp_path / "nope.lunaNMR"),
                     "--out", str(tmp_path / "o")]) == 1

    def test_dry_run_writes_nothing(self, tmp_path):
        from lunaNMR.cli import main
        b, out = _bundle(tmp_path), tmp_path / "exported"
        assert main(["project", "export", str(b), "--out", str(out), "--dry-run"]) == 0
        assert not out.exists()


class TestInventoryAccountsForEveryFile:
    """inventory's docstring says it describes a bundle's contents. It enumerated only
    DIRECTORIES under series_results/, so a loose file there appeared in no category —
    invisible in the Project Browser, and untargetable by `project remove`, which takes
    the paths inventory publishes.
    """

    def test_a_loose_series_file_is_listed(self, tmp_path):
        import types
        from lunaNMR.utils.project_manager import ProjectManager
        b = _bundle(tmp_path)
        cats = ProjectManager(types.SimpleNamespace()).inventory(b)
        paths = {p for c in cats for i in c['items'] for p in i['paths']}
        assert "series_results/batch_results.json" in paths

    def test_it_is_removable(self, tmp_path):
        """Listing it without making it removable would leave the same dead end."""
        import types
        from lunaNMR.utils.project_manager import ProjectManager
        b = _bundle(tmp_path)
        cats = ProjectManager(types.SimpleNamespace()).inventory(b)
        entry = [i for c in cats for i in c['items']
                 if "series_results/batch_results.json" in i['paths']]
        assert entry and entry[0]['removable'] is True

    def test_export_then_reports_nothing_unclassified(self, tmp_path, capsys):
        from lunaNMR.cli import main
        b = _bundle(tmp_path)
        main(["project", "export", str(b), "--out", str(tmp_path / "e"),
              "--format", "json"])
        assert json.loads(capsys.readouterr().out)["unclassified"] == []
