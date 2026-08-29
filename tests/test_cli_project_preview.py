# ABOUTME: `project remove` is the only destructive subcommand and had no preview.
# ABOUTME: inventory printed English labels; remove takes bundle-relative paths.

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))


def _bundle(tmp_path):
    bundle = tmp_path / "proj.lunaNMR"
    (bundle / "fit_results" / "arrays").mkdir(parents=True)
    (bundle / "project.json").write_text('{"schema_version": "1.1"}')
    (bundle / "gui_state.json").write_text('{}')
    (bundle / "fit_results" / "fits.json").write_text('[]')
    (bundle / "fit_results" / "arrays" / "peak_000.npz").write_bytes(b"x" * 5000)
    return bundle


class TestInventoryNamesWhatRemoveTakes:
    """remove takes bundle-relative paths; inventory printed only English labels, so
    there was no way to learn a valid argument from the tool itself. The paths were in
    the inventory data all along and were dropped at the print.
    """

    def test_json_carries_the_removable_paths(self, tmp_path, capsys):
        from lunaNMR.cli import main
        assert main(["project", "inventory", str(_bundle(tmp_path)), "--format", "json"]) == 0
        payload = json.loads(capsys.readouterr().out)
        paths = {p for cat in payload["categories"]
                 for item in cat["items"] if item["removable"] for p in item["paths"]}
        assert "fit_results/fits.json" in paths

    def test_the_text_listing_still_reads_as_before(self, tmp_path, capsys):
        from lunaNMR.cli import main
        assert main(["project", "inventory", str(_bundle(tmp_path))]) == 0
        out = capsys.readouterr().out
        assert "Fit results" in out


class TestRemoveHasAPreview:
    """One command deletes, immediately, with no undo. It should be possible to see what
    that means before it happens."""

    def test_dry_run_reports_without_deleting(self, tmp_path, capsys):
        from lunaNMR.cli import main
        bundle = _bundle(tmp_path)
        assert main(["project", "remove", str(bundle),
                     "fit_results/fits.json", "--dry-run"]) == 0
        assert (bundle / "fit_results" / "fits.json").exists()
        assert "fit_results/fits.json" in capsys.readouterr().out

    def test_dry_run_refuses_the_same_paths_the_real_run_would(self, tmp_path):
        """An escape must fail identically in preview, or the preview lies."""
        from lunaNMR.cli import main
        bundle = _bundle(tmp_path)
        assert main(["project", "remove", str(bundle), "../secret.txt", "--dry-run"]) == 1

    def test_without_dry_run_it_still_deletes(self, tmp_path):
        from lunaNMR.cli import main
        bundle = _bundle(tmp_path)
        assert main(["project", "remove", str(bundle), "fit_results/fits.json"]) == 0
        assert not (bundle / "fit_results" / "fits.json").exists()
