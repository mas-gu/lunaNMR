# ABOUTME: Project bundle browser — shows .lunaNMR contents by category with sizes.
# ABOUTME: Lets the user delete categories/items; deletions are immediate and irreversible.

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QTreeWidget, QTreeWidgetItem,
    QFileDialog, QMessageBox, QHeaderView,
)

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.utils.project_manager import ProjectManager


def _human_size(n: int) -> str:
    """Format a byte count as a short human-readable string."""
    size = float(n)
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if size < 1024 or unit == 'TB':
            return f"{int(size)} {unit}" if unit == 'B' else f"{size:.1f} {unit}"
        size /= 1024


class ProjectBrowserDialog(BaseDialog):
    """Browse a .lunaNMR bundle by category and remove items to clean it up."""

    def __init__(self, parent=None, main_window=None, project_path=None):
        super().__init__(
            parent=parent,
            title="Project Contents",
            default_size=(640, 520),
            min_size=(520, 400),
            modal=False,
        )
        self.main_window = main_window
        self.pm = (main_window._get_project_manager()
                   if main_window is not None and hasattr(main_window, '_get_project_manager')
                   else ProjectManager(main_window))
        self.project_path = Path(project_path) if project_path else None

        self._build_ui()
        if self.project_path:
            self._set_project(self.project_path)

    def _build_ui(self):
        layout = QVBoxLayout()

        path_row = QHBoxLayout()
        self.path_label = QLabel("No project selected")
        self.path_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        path_row.addWidget(self.path_label, stretch=1)
        browse_btn = QPushButton("Browse…")
        browse_btn.clicked.connect(self._on_browse)
        path_row.addWidget(browse_btn)
        layout.addLayout(path_row)

        self.tree = QTreeWidget()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(["Item", "Size"])
        self.tree.header().setSectionResizeMode(0, QHeaderView.Stretch)
        self.tree.header().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.tree.itemSelectionChanged.connect(self._update_remove_enabled)
        self.tree.itemSelectionChanged.connect(self._update_open_enabled)
        layout.addWidget(self.tree, stretch=1)

        btn_row = QHBoxLayout()
        self.total_label = QLabel("Total: —")
        btn_row.addWidget(self.total_label, stretch=1)
        self.open_analysis_btn = QPushButton("Open analysis")
        self.open_analysis_btn.setEnabled(False)
        self.open_analysis_btn.clicked.connect(self._on_open_analysis)
        btn_row.addWidget(self.open_analysis_btn)
        self.remove_btn = QPushButton("Remove selected")
        self.remove_btn.setEnabled(False)
        self.remove_btn.clicked.connect(self._on_remove)
        btn_row.addWidget(self.remove_btn)
        refresh_btn = QPushButton("Refresh")
        refresh_btn.clicked.connect(self._populate)
        btn_row.addWidget(refresh_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

        self.setLayout(layout)

    def _on_browse(self):
        start = str(self.project_path.parent) if self.project_path else ""
        folder = QFileDialog.getExistingDirectory(
            self, "Select a .lunaNMR project bundle", start)
        if folder:
            self._set_project(folder)

    def _set_project(self, folder):
        self.project_path = Path(folder)
        self.path_label.setText(str(self.project_path))
        self._populate()

    def _populate(self):
        self.tree.clear()
        if not self.project_path or not self.project_path.exists():
            self.total_label.setText("Total: —")
            self._update_remove_enabled()
            return

        total = 0
        for cat in self.pm.inventory(self.project_path):
            cat_paths = [p for it in cat['items'] if it['removable'] for p in it['paths']]
            cat_node = QTreeWidgetItem([cat['label'], _human_size(cat['size'])])
            cat_node.setData(0, Qt.UserRole, {
                'label': cat['label'], 'paths': cat_paths,
                'removable': bool(cat_paths), 'size': cat['size'],
            })
            for it in cat['items']:
                child = QTreeWidgetItem([it['label'], _human_size(it['size'])])
                child.setData(0, Qt.UserRole, {
                    'id': it.get('id'),
                    'label': it['label'], 'paths': it['paths'],
                    'removable': it['removable'], 'size': it['size'],
                })
                cat_node.addChild(child)
            self.tree.addTopLevelItem(cat_node)
            cat_node.setExpanded(True)
            total += cat['size']

        self.total_label.setText(f"Total: {_human_size(total)}")
        self._update_remove_enabled()

    def _selected_data(self):
        node = self.tree.currentItem()
        return node.data(0, Qt.UserRole) if node is not None else None

    def _update_remove_enabled(self):
        data = self._selected_data()
        self.remove_btn.setEnabled(bool(data and data.get('removable') and data.get('paths')))

    def _openable_kind(self, data):
        """'kd' / 'dynamixs' if the selected item is a reopenable saved analysis, else None."""
        item_id = str((data or {}).get('id') or '')
        if item_id.startswith('kd/analyses/'):
            return 'kd'
        if item_id.startswith('dynamixs/analyses/'):
            return 'dynamixs'
        return None

    def _update_open_enabled(self):
        can_open = self._openable_kind(self._selected_data()) is not None and self.main_window is not None
        self.open_analysis_btn.setEnabled(can_open)

    def _on_open_analysis(self):
        data = self._selected_data()
        kind = self._openable_kind(data)
        if kind is None or self.main_window is None:
            return
        opener = (self.main_window.open_kd_analysis if kind == 'kd'
                  else self.main_window.open_dynamixs_analysis)
        try:
            opened = opener(data['label'])
        except Exception as e:
            QMessageBox.critical(self, "Open failed", str(e))
            return
        if not opened:
            # The browser can point at any bundle, but open reads the loaded
            # project's analyses — so a non-loaded bundle has nothing to open.
            QMessageBox.information(
                self, "Not loaded",
                f"'{data['label']}' isn't in the loaded project. Open its project "
                "first (File ▸ Open), then reopen the analysis.")

    def _on_remove(self):
        data = self._selected_data()
        if not data or not data.get('removable') or not data.get('paths'):
            return
        resp = QMessageBox.question(
            self, "Remove from project",
            f"Remove '{data['label']}' ({_human_size(data['size'])})?\n\n"
            "This permanently deletes files from the bundle and cannot be undone.",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if resp != QMessageBox.Yes:
            return
        try:
            self.pm.remove_bundle_paths(self.project_path, data['paths'])
        except Exception as e:
            QMessageBox.critical(self, "Removal failed", str(e))
            return
        self._populate()
