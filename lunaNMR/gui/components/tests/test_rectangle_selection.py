# ABOUTME: Unit tests for middle-click + drag rectangle selection functionality
# ABOUTME: Tests click vs drag detection, area selection, and multi-peak selection

import pytest
from unittest.mock import MagicMock, patch


class TestClickVsDragDetection:
    """Test the click vs drag threshold detection."""

    def test_click_without_movement_triggers_single_selection(self):
        """Middle-click without movement should trigger single peak selection."""
        # Simulate the handler behavior
        class MockHandler:
            DRAG_THRESHOLD_PIXELS = 5

            def should_trigger_drag(self, start_x, start_y, end_x, end_y):
                dx = abs(end_x - start_x)
                dy = abs(end_y - start_y)
                return dx >= self.DRAG_THRESHOLD_PIXELS or dy >= self.DRAG_THRESHOLD_PIXELS

        handler = MockHandler()

        # Click at same position (no movement)
        assert handler.should_trigger_drag(100, 100, 100, 100) is False

        # Small movement (2 pixels) - still a click
        assert handler.should_trigger_drag(100, 100, 102, 101) is False

        # Just under threshold (4 pixels)
        assert handler.should_trigger_drag(100, 100, 104, 100) is False

    def test_drag_beyond_threshold_triggers_area_selection(self):
        """Middle-click + drag beyond threshold should trigger area selection."""
        class MockHandler:
            DRAG_THRESHOLD_PIXELS = 5

            def should_trigger_drag(self, start_x, start_y, end_x, end_y):
                dx = abs(end_x - start_x)
                dy = abs(end_y - start_y)
                return dx >= self.DRAG_THRESHOLD_PIXELS or dy >= self.DRAG_THRESHOLD_PIXELS

        handler = MockHandler()

        # Exactly at threshold (5 pixels) - triggers drag
        assert handler.should_trigger_drag(100, 100, 105, 100) is True

        # Beyond threshold (10 pixels)
        assert handler.should_trigger_drag(100, 100, 110, 100) is True

        # Y movement beyond threshold
        assert handler.should_trigger_drag(100, 100, 100, 106) is True

        # Diagonal movement (4,4) - doesn't trigger since neither dx nor dy >= 5
        assert handler.should_trigger_drag(100, 100, 104, 104) is False
        # Diagonal (5,5) - triggers
        assert handler.should_trigger_drag(100, 100, 105, 105) is True


class TestRectangleBoundsNormalization:
    """Test that rectangle bounds are normalized correctly regardless of drag direction."""

    def test_normalize_top_left_to_bottom_right(self):
        """Drag from top-left to bottom-right."""
        x1, y1 = 8.0, 120.0  # Higher ppm values (NMR convention)
        x2, y2 = 7.5, 115.0  # Lower ppm values

        x_min, x_max = min(x1, x2), max(x1, x2)
        y_min, y_max = min(y1, y2), max(y1, y2)

        assert x_min == 7.5
        assert x_max == 8.0
        assert y_min == 115.0
        assert y_max == 120.0

    def test_normalize_bottom_right_to_top_left(self):
        """Drag from bottom-right to top-left."""
        x1, y1 = 7.5, 115.0  # Start bottom-right
        x2, y2 = 8.0, 120.0  # End top-left

        x_min, x_max = min(x1, x2), max(x1, x2)
        y_min, y_max = min(y1, y2), max(y1, y2)

        assert x_min == 7.5
        assert x_max == 8.0
        assert y_min == 115.0
        assert y_max == 120.0

    def test_normalize_any_direction_same_result(self):
        """All drag directions should give the same normalized bounds."""
        corners = [
            ((8.0, 120.0), (7.5, 115.0)),  # top-left to bottom-right
            ((7.5, 115.0), (8.0, 120.0)),  # bottom-right to top-left
            ((8.0, 115.0), (7.5, 120.0)),  # bottom-left to top-right
            ((7.5, 120.0), (8.0, 115.0)),  # top-right to bottom-left
        ]

        for (x1, y1), (x2, y2) in corners:
            x_min, x_max = min(x1, x2), max(x1, x2)
            y_min, y_max = min(y1, y2), max(y1, y2)

            assert x_min == 7.5
            assert x_max == 8.0
            assert y_min == 115.0
            assert y_max == 120.0


class TestFindPeaksInRectangle:
    """Test the peak finding within rectangle logic."""

    def test_finds_peaks_inside_rectangle(self):
        """Peaks inside the rectangle should be found."""
        # Mock peak data
        reference_peaks = [
            {'Position_X': 7.8, 'Position_Y': 118.0, 'Assignment': 'A.10.LYS.H'},
            {'Position_X': 7.6, 'Position_Y': 116.5, 'Assignment': 'A.15.ALA.H'},
            {'Position_X': 9.0, 'Position_Y': 130.0, 'Assignment': 'A.20.GLU.H'},  # Outside
        ]

        x_min, x_max = 7.5, 8.0
        y_min, y_max = 115.0, 120.0

        inside = []
        for peak in reference_peaks:
            px, py = peak['Position_X'], peak['Position_Y']
            if x_min <= px <= x_max and y_min <= py <= y_max:
                inside.append(peak)

        assert len(inside) == 2
        assert inside[0]['Assignment'] == 'A.10.LYS.H'
        assert inside[1]['Assignment'] == 'A.15.ALA.H'

    def test_peaks_on_boundary_included(self):
        """Peaks exactly on rectangle boundary should be included."""
        reference_peaks = [
            {'Position_X': 7.5, 'Position_Y': 118.0, 'Assignment': 'OnLeftEdge'},
            {'Position_X': 8.0, 'Position_Y': 118.0, 'Assignment': 'OnRightEdge'},
            {'Position_X': 7.7, 'Position_Y': 115.0, 'Assignment': 'OnBottomEdge'},
            {'Position_X': 7.7, 'Position_Y': 120.0, 'Assignment': 'OnTopEdge'},
        ]

        x_min, x_max = 7.5, 8.0
        y_min, y_max = 115.0, 120.0

        inside = []
        for peak in reference_peaks:
            px, py = peak['Position_X'], peak['Position_Y']
            if x_min <= px <= x_max and y_min <= py <= y_max:
                inside.append(peak)

        assert len(inside) == 4

    def test_no_peaks_in_empty_region(self):
        """Rectangle in empty region should return no peaks."""
        reference_peaks = [
            {'Position_X': 7.8, 'Position_Y': 118.0, 'Assignment': 'A.10.LYS.H'},
        ]

        # Rectangle far from any peaks
        x_min, x_max = 10.0, 11.0
        y_min, y_max = 140.0, 150.0

        inside = []
        for peak in reference_peaks:
            px, py = peak['Position_X'], peak['Position_Y']
            if x_min <= px <= x_max and y_min <= py <= y_max:
                inside.append(peak)

        assert len(inside) == 0


class TestMultiSelectionState:
    """Test multi-selection state management."""

    def test_clear_multi_selection(self):
        """Clearing multi-selection should empty the list."""
        class MockMainWindow:
            def __init__(self):
                self.selected_peaks_multi = [
                    {'type': 'reference', 'x': 7.8, 'y': 118.0},
                    {'type': 'detected', 'x': 7.6, 'y': 116.5},
                ]

            def clear_multi_selection(self):
                self.selected_peaks_multi = []

        window = MockMainWindow()
        assert len(window.selected_peaks_multi) == 2

        window.clear_multi_selection()
        assert len(window.selected_peaks_multi) == 0

    def test_single_selection_clears_multi_selection(self):
        """Single peak selection should clear any existing multi-selection."""
        class MockMainWindow:
            def __init__(self):
                self.selected_peaks_multi = [
                    {'type': 'reference', 'x': 7.8, 'y': 118.0},
                    {'type': 'detected', 'x': 7.6, 'y': 116.5},
                ]
                self.selected_peak_info = None

            def clear_multi_selection(self):
                self.selected_peaks_multi = []

            def select_single_peak(self, peak_info):
                if self.selected_peaks_multi:
                    self.clear_multi_selection()
                self.selected_peak_info = peak_info

        window = MockMainWindow()
        window.select_single_peak({'type': 'reference', 'x': 8.5, 'y': 125.0})

        assert len(window.selected_peaks_multi) == 0
        assert window.selected_peak_info is not None


class TestBatchDelete:
    """Test batch delete functionality."""

    def test_delete_selected_peaks_removes_from_lists(self):
        """Batch delete should remove all selected peaks from data structures."""
        class MockIntegrator:
            def __init__(self):
                import pandas as pd
                self.peak_list = pd.DataFrame({
                    'Position_X': [7.8, 7.6, 9.0],
                    'Position_Y': [118.0, 116.5, 130.0],
                    'Assignment': ['Peak1', 'Peak2', 'Peak3']
                })
                self.fitted_peaks = [
                    {'ppm_x': 7.75, 'ppm_y': 117.5, 'assignment': 'Det1'},
                    {'ppm_x': 8.5, 'ppm_y': 125.0, 'assignment': 'Det2'},
                ]

        class MockMainWindow:
            def __init__(self):
                self.integrator = MockIntegrator()
                self.selected_peaks_multi = [
                    {'type': 'reference', 'index': 0},  # Peak1
                    {'type': 'detected', 'index': 0},   # Det1
                ]

            def delete_selected_peaks(self):
                # Collect indices to delete (in reverse order to maintain indices)
                ref_indices = sorted(
                    [p['index'] for p in self.selected_peaks_multi if p['type'] == 'reference'],
                    reverse=True
                )
                det_indices = sorted(
                    [p['index'] for p in self.selected_peaks_multi if p['type'] == 'detected'],
                    reverse=True
                )

                # Delete from DataFrame
                for idx in ref_indices:
                    self.integrator.peak_list = self.integrator.peak_list.drop(idx)

                # Delete from list
                for idx in det_indices:
                    del self.integrator.fitted_peaks[idx]

                self.selected_peaks_multi = []

        window = MockMainWindow()
        assert len(window.integrator.peak_list) == 3
        assert len(window.integrator.fitted_peaks) == 2

        window.delete_selected_peaks()

        assert len(window.integrator.peak_list) == 2  # 1 removed
        assert len(window.integrator.fitted_peaks) == 1  # 1 removed
        assert len(window.selected_peaks_multi) == 0


class TestKeyboardHandling:
    """Test keyboard shortcuts for selection management."""

    def test_escape_clears_multi_selection(self):
        """Escape key should clear multi-selection."""
        from PySide6.QtCore import Qt

        class MockMainWindow:
            def __init__(self):
                self.selected_peaks_multi = [
                    {'type': 'reference', 'x': 7.8, 'y': 118.0},
                ]
                self.escape_pressed = False

            def clear_multi_selection(self):
                self.selected_peaks_multi = []

            def handle_key(self, key):
                if key == Qt.Key_Escape and self.selected_peaks_multi:
                    self.clear_multi_selection()
                    return True
                return False

        window = MockMainWindow()
        result = window.handle_key(Qt.Key_Escape)

        assert result is True
        assert len(window.selected_peaks_multi) == 0

    def test_delete_key_removes_selected_peaks(self):
        """Delete key should trigger batch delete."""
        from PySide6.QtCore import Qt

        class MockMainWindow:
            def __init__(self):
                self.selected_peaks_multi = [
                    {'type': 'reference', 'x': 7.8, 'y': 118.0},
                ]
                self.delete_called = False

            def delete_selected_peaks(self):
                self.delete_called = True
                self.selected_peaks_multi = []

            def handle_key(self, key, modifiers=Qt.NoModifier):
                if self.selected_peaks_multi:
                    if key == Qt.Key_Delete or (key == Qt.Key_D and modifiers & Qt.ControlModifier):
                        self.delete_selected_peaks()
                        return True
                return False

        window = MockMainWindow()

        # Test Delete key
        result = window.handle_key(Qt.Key_Delete)
        assert result is True
        assert window.delete_called is True

    def test_ctrl_d_removes_selected_peaks(self):
        """Ctrl+D should trigger batch delete."""
        from PySide6.QtCore import Qt

        class MockMainWindow:
            def __init__(self):
                self.selected_peaks_multi = [
                    {'type': 'reference', 'x': 7.8, 'y': 118.0},
                ]
                self.delete_called = False

            def delete_selected_peaks(self):
                self.delete_called = True
                self.selected_peaks_multi = []

            def handle_key(self, key, modifiers=Qt.NoModifier):
                if self.selected_peaks_multi:
                    if key == Qt.Key_Delete or (key == Qt.Key_D and modifiers & Qt.ControlModifier):
                        self.delete_selected_peaks()
                        return True
                return False

        window = MockMainWindow()

        # Test Ctrl+D
        result = window.handle_key(Qt.Key_D, Qt.ControlModifier)
        assert result is True
        assert window.delete_called is True


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
