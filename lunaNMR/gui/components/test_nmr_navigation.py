#!/usr/bin/env python3
# ABOUTME: Unit tests for NMR navigation handler zoom/pan algorithms
# ABOUTME: Tests the mathematical correctness of zoom-at-cursor and pan with inverted axes

"""
Tests for NMRNavigationHandler algorithms.

These tests verify the mathematical correctness of:
- Zoom-at-cursor with inverted NMR axes
- Pan direction handling with inverted axes
"""

import pytest
import numpy as np
import sys
import os

# Add parent directory to path for standalone testing
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

from lunaNMR.gui.components.nmr_navigation_handler import calculate_zoom_limits, calculate_pan_limits


class TestZoomAtCursorMath:
    """Test the zoom-at-cursor algorithm with inverted NMR axes."""

    def test_zoom_in_preserves_cursor_position(self):
        """Cursor position should remain fixed after zoom in."""

        # Initial limits (inverted NMR axes: x_min > x_max)
        x_lim = (10.0, 6.0)  # High ppm on left
        y_lim = (130.0, 110.0)  # High ppm on bottom

        # Cursor at center of view
        x_cursor = 8.0
        y_cursor = 120.0

        # Zoom in by factor of 2
        new_x_lim, new_y_lim = calculate_zoom_limits(
            x_lim, y_lim, x_cursor, y_cursor, factor=2.0
        )

        # Cursor should still be at same relative position in new view
        # For inverted axis: rel = (cursor - max) / (min - max)
        old_rel_x = (x_cursor - x_lim[1]) / (x_lim[0] - x_lim[1])
        new_rel_x = (x_cursor - new_x_lim[1]) / (new_x_lim[0] - new_x_lim[1])

        old_rel_y = (y_cursor - y_lim[1]) / (y_lim[0] - y_lim[1])
        new_rel_y = (y_cursor - new_y_lim[1]) / (new_y_lim[0] - new_y_lim[1])

        assert abs(old_rel_x - new_rel_x) < 0.001, f"X relative position changed: {old_rel_x} -> {new_rel_x}"
        assert abs(old_rel_y - new_rel_y) < 0.001, f"Y relative position changed: {old_rel_y} -> {new_rel_y}"

        # Range should be halved (zoomed in 2x)
        old_x_range = abs(x_lim[0] - x_lim[1])
        new_x_range = abs(new_x_lim[0] - new_x_lim[1])
        assert abs(new_x_range - old_x_range / 2.0) < 0.001

    def test_zoom_out_preserves_cursor_position(self):
        """Cursor position should remain fixed after zoom out."""
        # Initial limits (inverted)
        x_lim = (10.0, 6.0)
        y_lim = (130.0, 110.0)

        # Cursor off-center
        x_cursor = 9.0
        y_cursor = 125.0

        # Zoom out by factor of 0.5
        new_x_lim, new_y_lim = calculate_zoom_limits(
            x_lim, y_lim, x_cursor, y_cursor, factor=0.5
        )

        # Cursor relative position should be preserved
        old_rel_x = (x_cursor - x_lim[1]) / (x_lim[0] - x_lim[1])
        new_rel_x = (x_cursor - new_x_lim[1]) / (new_x_lim[0] - new_x_lim[1])

        assert abs(old_rel_x - new_rel_x) < 0.001

        # Range should be doubled (zoomed out 2x)
        old_x_range = abs(x_lim[0] - x_lim[1])
        new_x_range = abs(new_x_lim[0] - new_x_lim[1])
        assert abs(new_x_range - old_x_range * 2.0) < 0.001


class TestPanMath:
    """Test the pan algorithm with inverted NMR axes."""

    def test_pan_right_decreases_ppm(self):
        """Dragging right should shift view to lower ppm (right side of NMR spectrum)."""
        # Initial limits (inverted)
        x_lim = (10.0, 6.0)  # View from 10 to 6 ppm
        y_lim = (130.0, 110.0)

        # Drag from left to right in data coordinates
        # Start at 9 ppm, end at 7 ppm (moved 2 ppm to the right visually)
        start_x, start_y = 9.0, 120.0
        end_x, end_y = 7.0, 120.0

        new_x_lim, new_y_lim = calculate_pan_limits(
            x_lim, y_lim, start_x, start_y, end_x, end_y
        )

        # View should shift to lower ppm values
        # Original: 10-6, after panning right: should be ~8-4
        assert new_x_lim[0] < x_lim[0], "X min should decrease (shift right)"
        assert new_x_lim[1] < x_lim[1], "X max should decrease (shift right)"

    def test_pan_preserves_range(self):
        """Pan should not change the zoom level (range should be constant)."""
        x_lim = (10.0, 6.0)
        y_lim = (130.0, 110.0)

        start_x, start_y = 8.0, 120.0
        end_x, end_y = 7.5, 118.0

        new_x_lim, new_y_lim = calculate_pan_limits(
            x_lim, y_lim, start_x, start_y, end_x, end_y
        )

        old_x_range = abs(x_lim[0] - x_lim[1])
        new_x_range = abs(new_x_lim[0] - new_x_lim[1])
        old_y_range = abs(y_lim[0] - y_lim[1])
        new_y_range = abs(new_y_lim[0] - new_y_lim[1])

        assert abs(old_x_range - new_x_range) < 0.001
        assert abs(old_y_range - new_y_range) < 0.001


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
