# ABOUTME: Undo/redo system for peak manipulation using Qt's QUndoStack
# ABOUTME: Supports add, delete, move operations with configurable history depth

"""
UndoManager for Peak Manipulation

Provides undo/redo functionality for:
- Add peak (Shift+click)
- Delete peak (Ctrl+click, Delete key, or middle-click + Delete)
- Move peak (middle-click then click new position)

Uses Qt's QUndoStack for:
- Automatic menu integration (Edit -> Undo/Redo)
- Keyboard shortcuts (Cmd+Z, Cmd+Shift+Z)
- Clean/dirty state tracking
- Command merging (consecutive moves of same peak)
"""

from PySide6.QtGui import QUndoCommand, QUndoStack
from PySide6.QtCore import QObject, Signal
import copy
import logging
import pandas as pd

logger = logging.getLogger(__name__)


class DeletePeakCommand(QUndoCommand):
    """Command to delete a single peak."""

    def __init__(self, main_window, peak_id: int, peak_type: str, description: str = None):
        """
        Args:
            main_window: Reference to LunaNMRMainWindow
            peak_id: Stable peak identifier (not row index)
            peak_type: 'reference' or 'detected'
            description: Optional custom description for undo menu
        """
        desc = description or f"Delete {peak_type} peak"
        super().__init__(desc)
        self.main_window = main_window
        self.peak_id = peak_id
        self.peak_type = peak_type
        self.snapshot = None  # Stored on first redo()

    def redo(self):
        """Execute deletion (called on first push and on redo)."""
        integrator = self.main_window.integrator

        if self.peak_type == 'reference':
            # Find row by peak_id, store snapshot, delete
            if integrator.peak_list is not None and 'peak_id' in integrator.peak_list.columns:
                mask = integrator.peak_list['peak_id'] == self.peak_id
                if mask.any():
                    self.snapshot = integrator.peak_list[mask].iloc[0].to_dict()
                    integrator.peak_list = integrator.peak_list[~mask]
                    logger.info(f"Deleted reference peak (peak_id={self.peak_id})")

        elif self.peak_type == 'detected':
            # Find in list by peak_id, store snapshot, delete
            for i, peak in enumerate(integrator.fitted_peaks):
                if peak.get('peak_id') == self.peak_id:
                    self.snapshot = copy.deepcopy(peak)
                    del integrator.fitted_peaks[i]
                    logger.info(f"Deleted detected peak (peak_id={self.peak_id})")
                    break

        self._refresh_ui()

    def undo(self):
        """Restore the deleted peak."""
        if self.snapshot is None:
            logger.warning(f"Cannot undo delete: no snapshot for peak_id={self.peak_id}")
            return

        integrator = self.main_window.integrator

        if self.peak_type == 'reference':
            # Re-insert row into DataFrame
            new_row = pd.DataFrame([self.snapshot])
            integrator.peak_list = pd.concat([integrator.peak_list, new_row], ignore_index=True)
            logger.info(f"Restored reference peak (peak_id={self.peak_id})")

        elif self.peak_type == 'detected':
            # Re-insert into list
            integrator.fitted_peaks.append(copy.deepcopy(self.snapshot))
            logger.info(f"Restored detected peak (peak_id={self.peak_id})")

        self._refresh_ui()

    def _refresh_ui(self):
        """Update UI after change."""
        self.main_window.update_spectrum_plot(preserve_zoom=True)
        if hasattr(self.main_window, 'peak_navigator'):
            self.main_window.peak_navigator.refresh_peak_list()
        if hasattr(self.main_window, 'update_statistics_panel'):
            self.main_window.update_statistics_panel()


class DeleteMultiplePeaksCommand(QUndoCommand):
    """Command to delete multiple peaks at once (e.g., from rectangle selection)."""

    def __init__(self, main_window, peaks_info: list):
        """
        Args:
            main_window: Reference to LunaNMRMainWindow
            peaks_info: List of dicts with 'peak_id' and 'type' keys
        """
        count = len(peaks_info)
        super().__init__(f"Delete {count} peaks")
        self.main_window = main_window
        self.peaks_info = peaks_info
        self.snapshots = []  # Populated on first redo

    def redo(self):
        """Delete all peaks, storing snapshots for undo."""
        self.snapshots = []
        integrator = self.main_window.integrator

        for info in self.peaks_info:
            peak_id = info['peak_id']
            peak_type = info['type']

            if peak_type == 'reference':
                if integrator.peak_list is not None and 'peak_id' in integrator.peak_list.columns:
                    mask = integrator.peak_list['peak_id'] == peak_id
                    if mask.any():
                        snapshot = integrator.peak_list[mask].iloc[0].to_dict()
                        self.snapshots.append({'type': 'reference', 'data': snapshot})
                        integrator.peak_list = integrator.peak_list[~mask]

            elif peak_type == 'detected':
                for i, peak in enumerate(integrator.fitted_peaks):
                    if peak.get('peak_id') == peak_id:
                        self.snapshots.append({'type': 'detected', 'data': copy.deepcopy(peak)})
                        del integrator.fitted_peaks[i]
                        break

        logger.info(f"Deleted {len(self.snapshots)} peaks")
        self._refresh_ui()

    def undo(self):
        """Restore all deleted peaks."""
        integrator = self.main_window.integrator

        for snapshot in self.snapshots:
            if snapshot['type'] == 'reference':
                new_row = pd.DataFrame([snapshot['data']])
                integrator.peak_list = pd.concat([integrator.peak_list, new_row], ignore_index=True)
            elif snapshot['type'] == 'detected':
                integrator.fitted_peaks.append(copy.deepcopy(snapshot['data']))

        logger.info(f"Restored {len(self.snapshots)} peaks")
        self._refresh_ui()

    def _refresh_ui(self):
        """Update UI after change."""
        self.main_window.update_spectrum_plot(preserve_zoom=True)
        if hasattr(self.main_window, 'peak_navigator'):
            self.main_window.peak_navigator.refresh_peak_list()
        if hasattr(self.main_window, 'update_statistics_panel'):
            self.main_window.update_statistics_panel()


class MovePeakCommand(QUndoCommand):
    """Command to move a peak to a new position."""

    # Command ID for merging consecutive moves of the same peak
    MERGE_ID = 1001

    def __init__(self, main_window, peak_id: int, peak_type: str,
                 old_x: float, old_y: float, new_x: float, new_y: float):
        """
        Args:
            main_window: Reference to LunaNMRMainWindow
            peak_id: Stable peak identifier
            peak_type: 'reference' or 'detected'
            old_x, old_y: Original position (ppm)
            new_x, new_y: New position (ppm)
        """
        super().__init__(f"Move {peak_type} peak")
        self.main_window = main_window
        self.peak_id = peak_id
        self.peak_type = peak_type
        self.old_x, self.old_y = old_x, old_y
        self.new_x, self.new_y = new_x, new_y

    def redo(self):
        """Move peak to new position."""
        self._set_position(self.new_x, self.new_y)
        logger.info(f"Moved {self.peak_type} peak (peak_id={self.peak_id}) to ({self.new_x:.3f}, {self.new_y:.1f})")

    def undo(self):
        """Move peak back to original position."""
        self._set_position(self.old_x, self.old_y)
        logger.info(f"Restored {self.peak_type} peak (peak_id={self.peak_id}) to ({self.old_x:.3f}, {self.old_y:.1f})")

    def _set_position(self, x: float, y: float):
        """Set peak position in data structure."""
        integrator = self.main_window.integrator

        if self.peak_type == 'reference':
            if integrator.peak_list is not None and 'peak_id' in integrator.peak_list.columns:
                mask = integrator.peak_list['peak_id'] == self.peak_id
                integrator.peak_list.loc[mask, 'Position_X'] = x
                integrator.peak_list.loc[mask, 'Position_Y'] = y

        elif self.peak_type == 'detected':
            for peak in integrator.fitted_peaks:
                if peak.get('peak_id') == self.peak_id:
                    peak['ppm_x'] = x
                    peak['ppm_y'] = y
                    break

        self.main_window.update_spectrum_plot(preserve_zoom=True)

    def id(self):
        """Return command ID for merging."""
        return self.MERGE_ID

    def mergeWith(self, other):
        """Merge consecutive moves of the same peak.

        This combines multiple small moves into a single undo operation.
        """
        if not isinstance(other, MovePeakCommand):
            return False
        if other.peak_id != self.peak_id or other.peak_type != self.peak_type:
            return False
        # Keep our old position, take their new position
        self.new_x, self.new_y = other.new_x, other.new_y
        logger.debug(f"Merged move commands for peak_id={self.peak_id}")
        return True


class AddPeakCommand(QUndoCommand):
    """Command to add a new peak."""

    def __init__(self, main_window, peak_data: dict, peak_type: str):
        """
        Args:
            main_window: Reference to LunaNMRMainWindow
            peak_data: Dict with peak data including 'peak_id'
            peak_type: 'reference' or 'detected'
        """
        assignment = peak_data.get('Assignment', peak_data.get('assignment', 'Unknown'))
        super().__init__(f"Add peak '{assignment}'")
        self.main_window = main_window
        self.peak_data = copy.deepcopy(peak_data)
        self.peak_type = peak_type
        self.peak_id = peak_data.get('peak_id')

    def redo(self):
        """Add the peak."""
        integrator = self.main_window.integrator

        if self.peak_type == 'reference':
            new_row = pd.DataFrame([self.peak_data])
            integrator.peak_list = pd.concat([integrator.peak_list, new_row], ignore_index=True)
            logger.info(f"Added reference peak (peak_id={self.peak_id})")

        elif self.peak_type == 'detected':
            integrator.fitted_peaks.append(copy.deepcopy(self.peak_data))
            logger.info(f"Added detected peak (peak_id={self.peak_id})")

        self._refresh_ui()

    def undo(self):
        """Remove the added peak."""
        integrator = self.main_window.integrator

        if self.peak_type == 'reference':
            if integrator.peak_list is not None and 'peak_id' in integrator.peak_list.columns:
                mask = integrator.peak_list['peak_id'] == self.peak_id
                integrator.peak_list = integrator.peak_list[~mask]
                logger.info(f"Removed reference peak (peak_id={self.peak_id})")

        elif self.peak_type == 'detected':
            integrator.fitted_peaks = [
                p for p in integrator.fitted_peaks
                if p.get('peak_id') != self.peak_id
            ]
            logger.info(f"Removed detected peak (peak_id={self.peak_id})")

        self._refresh_ui()

    def _refresh_ui(self):
        """Update UI after change."""
        self.main_window.update_spectrum_plot(preserve_zoom=True)
        if hasattr(self.main_window, 'peak_navigator'):
            self.main_window.peak_navigator.refresh_peak_list()
        if hasattr(self.main_window, 'update_statistics_panel'):
            self.main_window.update_statistics_panel()


class PeakUndoManager(QObject):
    """Manager for peak manipulation undo/redo operations.

    Wraps QUndoStack with additional convenience methods and signals.
    """

    # Emitted when undo/redo state changes (for updating menu state)
    can_undo_changed = Signal(bool)
    can_redo_changed = Signal(bool)

    # Stack limits
    MAX_UNDO_DEPTH = 30
    POST_SAVE_KEEP = 20

    def __init__(self, main_window, parent=None):
        """
        Args:
            main_window: Reference to LunaNMRMainWindow
            parent: Optional QObject parent
        """
        super().__init__(parent)
        self.main_window = main_window
        self.undo_stack = QUndoStack(self)
        self.undo_stack.setUndoLimit(self.MAX_UNDO_DEPTH)

        # Connect signals
        self.undo_stack.canUndoChanged.connect(self.can_undo_changed.emit)
        self.undo_stack.canRedoChanged.connect(self.can_redo_changed.emit)

    def push(self, command: QUndoCommand):
        """Push a command onto the undo stack.

        The command's redo() method is called automatically.
        """
        self.undo_stack.push(command)

    def undo(self):
        """Undo the last command."""
        if self.undo_stack.canUndo():
            self.undo_stack.undo()
            return True
        return False

    def redo(self):
        """Redo the last undone command."""
        if self.undo_stack.canRedo():
            self.undo_stack.redo()
            return True
        return False

    def can_undo(self) -> bool:
        """Check if undo is available."""
        return self.undo_stack.canUndo()

    def can_redo(self) -> bool:
        """Check if redo is available."""
        return self.undo_stack.canRedo()

    def clear(self):
        """Clear all undo/redo history."""
        self.undo_stack.clear()
        logger.info("Undo stack cleared")

    def on_file_load(self):
        """Called when a new file is loaded. Clears undo history."""
        self.clear()

    def on_major_operation(self):
        """Called before major operations like 'Fit All Peaks'.

        Clears undo history since fitted_peaks will be completely replaced.
        """
        self.clear()

    def on_save(self):
        """Called when project is saved.

        Marks the current state as clean. History is preserved.
        """
        self.undo_stack.setClean()
        logger.info("Undo stack marked as clean (saved)")

    def is_clean(self) -> bool:
        """Check if current state matches last saved state."""
        return self.undo_stack.isClean()

    def create_undo_action(self, parent):
        """Create a QAction for Undo menu item.

        The action is automatically enabled/disabled and updates its text.
        """
        return self.undo_stack.createUndoAction(parent, "&Undo")

    def create_redo_action(self, parent):
        """Create a QAction for Redo menu item.

        The action is automatically enabled/disabled and updates its text.
        """
        return self.undo_stack.createRedoAction(parent, "&Redo")

    # Convenience methods for creating commands

    def delete_peak(self, peak_id: int, peak_type: str):
        """Create and push a delete peak command."""
        cmd = DeletePeakCommand(self.main_window, peak_id, peak_type)
        self.push(cmd)

    def delete_multiple_peaks(self, peaks_info: list):
        """Create and push a delete multiple peaks command."""
        if not peaks_info:
            return
        cmd = DeleteMultiplePeaksCommand(self.main_window, peaks_info)
        self.push(cmd)

    def move_peak(self, peak_id: int, peak_type: str,
                  old_x: float, old_y: float, new_x: float, new_y: float):
        """Create and push a move peak command."""
        cmd = MovePeakCommand(self.main_window, peak_id, peak_type,
                              old_x, old_y, new_x, new_y)
        self.push(cmd)

    def add_peak(self, peak_data: dict, peak_type: str):
        """Create and push an add peak command."""
        cmd = AddPeakCommand(self.main_window, peak_data, peak_type)
        self.push(cmd)
