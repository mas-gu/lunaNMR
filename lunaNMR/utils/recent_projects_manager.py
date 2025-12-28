# ABOUTME: Manages recently opened projects using QSettings for persistence
# ABOUTME: Provides add/get/clear operations for a bounded list of recent project paths

"""
RecentProjectsManager - Persistent Recent Projects Storage

Uses Qt's QSettings for cross-platform storage of recently opened/saved projects.
Maintains a bounded list of the most recent project paths.
"""

from pathlib import Path
from typing import List
from PySide6.QtCore import QSettings


class RecentProjectsManager:
    """Manages a list of recently opened projects using QSettings.

    Stores up to MAX_RECENT project paths, with the most recent first.
    Automatically removes duplicates and non-existent paths.
    """

    MAX_RECENT = 5
    SETTINGS_KEY = "recentProjects/paths"

    def __init__(self):
        """Initialize the manager with QSettings."""
        self._settings = QSettings("LunaNMR", "LunaNMR_v1")

    def get_recent_projects(self) -> List[Path]:
        """Get list of recent project paths.

        Returns:
            List of Path objects for recent projects, most recent first.
            Only includes paths that still exist on disk.
        """
        paths_str = self._settings.value(self.SETTINGS_KEY, [])

        # Handle case where QSettings returns a string (single item)
        if isinstance(paths_str, str):
            paths_str = [paths_str] if paths_str else []

        # Convert to Path objects and filter non-existent
        paths = []
        for p in paths_str:
            if p:
                path = Path(p)
                if path.exists():
                    paths.append(path)

        return paths[:self.MAX_RECENT]

    def add_recent_project(self, project_path: Path) -> None:
        """Add a project to the recent list.

        Args:
            project_path: Path to the project file/folder

        The project is added at the front of the list. If it already
        exists in the list, it's moved to the front. The list is
        trimmed to MAX_RECENT items.
        """
        project_path = Path(project_path).resolve()

        # Get current list
        current = self.get_recent_projects()

        # Remove if already present (will re-add at front)
        current = [p for p in current if p.resolve() != project_path]

        # Add at front
        current.insert(0, project_path)

        # Trim to max size
        current = current[:self.MAX_RECENT]

        # Save as strings
        paths_str = [str(p) for p in current]
        self._settings.setValue(self.SETTINGS_KEY, paths_str)
        self._settings.sync()

    def clear_recent_projects(self) -> None:
        """Clear all recent projects."""
        self._settings.setValue(self.SETTINGS_KEY, [])
        self._settings.sync()

    def remove_project(self, project_path: Path) -> None:
        """Remove a specific project from the recent list.

        Args:
            project_path: Path to remove
        """
        project_path = Path(project_path).resolve()
        current = self.get_recent_projects()
        current = [p for p in current if p.resolve() != project_path]

        paths_str = [str(p) for p in current]
        self._settings.setValue(self.SETTINGS_KEY, paths_str)
        self._settings.sync()
