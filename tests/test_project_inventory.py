# ABOUTME: Tests the project-bundle inventory + safe removal used by the project browser.
# ABOUTME: Inventory/removal run headless; one test drives the browser dialog under Qt.

import os
import types
from unittest.mock import patch

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from lunaNMR.utils.project_manager import ProjectManager


def _pm():
    return ProjectManager(types.SimpleNamespace())


def _make_bundle(tmp_path):
    bundle = tmp_path / "proj.lunaNMR"
    (bundle / "fit_results" / "arrays").mkdir(parents=True)
    (bundle / "dynamixs" / "results" / "methyl_t2").mkdir(parents=True)
    (bundle / "kd" / "analyses" / "A1_assi").mkdir(parents=True)

    (bundle / "project.json").write_text('{"schema_version": "1.1"}')
    (bundle / "gui_state.json").write_text('{"a": 1}')
    (bundle / "fit_results" / "fits.json").write_text('[]')
    (bundle / "fit_results" / "arrays" / "peak_000_region_2d.npz").write_bytes(b"x" * 5000)
    (bundle / "dynamixs" / "state.json").write_text('{}')
    (bundle / "dynamixs" / "results" / "methyl_t2" / "r.csv").write_text("a,b\n")
    (bundle / "kd" / "analyses" / "A1_assi" / "fit_data.json").write_text("{}")
    (bundle / "kd" / "analyses" / "A1_assi" / "meta.json").write_text("{}")
    return bundle


def _cat(inv, cat_id):
    return next(c for c in inv if c['id'] == cat_id)


def _item(inv, cat_id, item_id):
    return next(i for i in _cat(inv, cat_id)['items'] if i['id'] == item_id)


def test_inventory_lists_present_categories(tmp_path):
    bundle = _make_bundle(tmp_path)
    inv = _pm().inventory(bundle)
    ids = {c['id'] for c in inv}
    assert {'manifest', 'gui_state', 'fit_results',
            'dynamixs', 'kd'} <= ids
    # peak_list.csv was never created → category absent
    assert 'peak_list' not in ids


def test_inventory_sizes_and_embedded_region_flag(tmp_path):
    bundle = _make_bundle(tmp_path)
    inv = _pm().inventory(bundle)
    arrays = _item(inv, 'fit_results', 'fit_results/arrays')
    assert arrays['size'] >= 5000
    assert arrays['removable'] is True
    # Manifest is core and must not be removable
    assert _item(inv, 'manifest', 'manifest')['removable'] is False


def test_inventory_lists_dynamixs_result_subdirs(tmp_path):
    bundle = _make_bundle(tmp_path)
    inv = _pm().inventory(bundle)
    item_ids = {i['id'] for i in _cat(inv, 'dynamixs')['items']}
    assert 'dynamixs/results/methyl_t2' in item_ids


def test_remove_bundle_paths_deletes_and_reports_size(tmp_path):
    bundle = _make_bundle(tmp_path)
    pm = _pm()
    freed = pm.remove_bundle_paths(bundle, ['fit_results/arrays'])
    assert freed >= 5000
    assert not (bundle / "fit_results" / "arrays").exists()
    # Other content untouched
    assert (bundle / "fit_results" / "fits.json").exists()
    # Category disappears from inventory once emptied
    inv = pm.inventory(bundle)
    arrays_items = [i for i in _cat(inv, 'fit_results')['items']
                    if i['id'] == 'fit_results/arrays']
    assert arrays_items == []


def test_remove_bundle_paths_rejects_escape(tmp_path):
    bundle = _make_bundle(tmp_path)
    sibling = tmp_path / "secret.txt"
    sibling.write_text("keep me")
    pm = _pm()
    with pytest.raises(ValueError):
        pm.remove_bundle_paths(bundle, ['../secret.txt'])
    assert sibling.exists()


def test_remove_bundle_paths_rejects_bundle_root(tmp_path):
    bundle = _make_bundle(tmp_path)
    pm = _pm()
    with pytest.raises(ValueError):
        pm.remove_bundle_paths(bundle, ['.'])
    assert bundle.exists()


def test_project_browser_dialog_populates_and_removes(tmp_path):
    from PySide6.QtWidgets import QApplication, QMessageBox
    from lunaNMR.gui.dialogs.project_browser_dialog import ProjectBrowserDialog

    app = QApplication.instance() or QApplication([])
    bundle = _make_bundle(tmp_path)
    dlg = ProjectBrowserDialog(parent=None, main_window=None, project_path=bundle)
    assert dlg.tree.topLevelItemCount() > 0

    # Locate the "Embedded region data" leaf and select it.
    target = None
    for i in range(dlg.tree.topLevelItemCount()):
        cat = dlg.tree.topLevelItem(i)
        for j in range(cat.childCount()):
            if cat.child(j).text(0) == "Embedded region data":
                target = cat.child(j)
    assert target is not None
    dlg.tree.setCurrentItem(target)
    assert dlg.remove_btn.isEnabled()

    with patch.object(QMessageBox, "question", return_value=QMessageBox.Yes):
        dlg._on_remove()

    assert not (bundle / "fit_results" / "arrays").exists()
    dlg.deleteLater()


def test_series_results_listed_individually(tmp_path):
    # Series results should list each series by name (like Kd analyses), so they can
    # be removed individually — not one all-or-nothing "All series results" item.
    bundle = tmp_path / "proj.lunaNMR"
    for name in ("A1_HSPA1A", "B4_HSPA8"):
        d = bundle / "series_results" / name
        d.mkdir(parents=True)
        (d / "series_analysis_tidy.csv").write_text("spectrum_name,assignment\n0,R1\n")
    inv = _pm().inventory(bundle)
    series = next(c for c in inv if c['id'] == 'series_results')
    labels = {it['label'] for it in series['items']}
    assert labels == {"A1_HSPA1A", "B4_HSPA8"}
    assert all(it['removable'] for it in series['items'])
    assert all(it['id'].startswith('series_results/') for it in series['items'])
