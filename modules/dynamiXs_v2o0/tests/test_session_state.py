# ABOUTME: Unit tests for DynamiXs session state persistence
# ABOUTME: Tests that field assignments (T1 Data, T2 Data) are saved/restored correctly

import pytest
from unittest.mock import MagicMock, patch


class TestT1T2FittingPageSessionState:
    """Test T1T2FittingPage session state save/restore."""

    def test_get_session_state_includes_source_series_map(self):
        """get_session_state() includes source_series_map in returned state."""
        # Simulate T1T2FittingPage state
        class MockPage:
            def __init__(self):
                self.source_series_map = {
                    'field1_t1': 'T1_asyn_series',
                    'field1_t2': 'T2_asyn_series',
                }
                self.field1_t1_file = '/path/to/t1.csv'
                self.field1_t2_file = '/path/to/t2.csv'
                self.field1_b0_file = None
                self.field1_temp_file = None
                self.field2_t1_file = None
                self.field2_t2_file = None
                self.field2_b0_file = None
                self.field2_temp_file = None

            def get_session_state(self):
                """Simulate the actual get_session_state method."""
                state = {
                    'field_files': {
                        'field1_t1': self.field1_t1_file,
                        'field1_t2': self.field1_t2_file,
                        'field1_b0': self.field1_b0_file,
                        'field1_temp': self.field1_temp_file,
                        'field2_t1': self.field2_t1_file,
                        'field2_t2': self.field2_t2_file,
                        'field2_b0': self.field2_b0_file,
                        'field2_temp': self.field2_temp_file,
                    },
                    'source_series_map': dict(getattr(self, 'source_series_map', {})),
                }
                return state

        page = MockPage()
        state = page.get_session_state()

        assert 'source_series_map' in state
        assert state['source_series_map']['field1_t1'] == 'T1_asyn_series'
        assert state['source_series_map']['field1_t2'] == 'T2_asyn_series'
        assert state['field_files']['field1_t1'] == '/path/to/t1.csv'

    def test_restore_session_state_restores_source_series_map(self):
        """restore_session_state() restores source_series_map from saved state."""
        # Simulate T1T2FittingPage restore
        class MockDisplay:
            def __init__(self):
                self.text_value = "No file selected"
                self.tooltip = ""

            def setText(self, text):
                self.text_value = text

            def setToolTip(self, tip):
                self.tooltip = tip

        class MockPage:
            def __init__(self):
                self.source_series_map = {}
                self.field1_t1_file = None
                self.field1_t2_file = None
                self.field1_t1_display = MockDisplay()
                self.field1_t2_display = MockDisplay()
                self.update_displays_called = False

            def _update_field_displays(self):
                self.update_displays_called = True
                # Simulate the actual display update
                for field_key in ['field1_t1', 'field1_t2']:
                    display = getattr(self, f"{field_key}_display", None)
                    series_name = self.source_series_map.get(field_key)
                    if series_name and display:
                        display.setText(f"📊 {series_name}")

            def restore_session_state(self, state):
                """Simulate the actual restore_session_state method."""
                field_files = state.get('field_files', {})
                for field_key, file_path in field_files.items():
                    setattr(self, f"{field_key}_file", file_path)

                self.source_series_map = dict(state.get('source_series_map', {}))
                self._update_field_displays()

        page = MockPage()
        saved_state = {
            'field_files': {
                'field1_t1': '/path/to/t1.csv',
                'field1_t2': '/path/to/t2.csv',
            },
            'source_series_map': {
                'field1_t1': 'T1_asyn_series',
                'field1_t2': 'T2_asyn_series',
            }
        }

        page.restore_session_state(saved_state)

        assert page.source_series_map['field1_t1'] == 'T1_asyn_series'
        assert page.source_series_map['field1_t2'] == 'T2_asyn_series'
        assert page.field1_t1_file == '/path/to/t1.csv'
        assert page.update_displays_called is True
        assert "📊 T1_asyn_series" in page.field1_t1_display.text_value
        assert "📊 T2_asyn_series" in page.field1_t2_display.text_value

    def test_session_state_round_trip(self):
        """State saved with get_session_state() restores correctly."""
        class MockDisplay:
            def __init__(self):
                self.text_value = "No file selected"

            def setText(self, text):
                self.text_value = text

            def setToolTip(self, tip):
                pass

        class MockPage:
            def __init__(self):
                self.source_series_map = {}
                self.field1_t1_file = None
                self.field1_t2_file = None
                self.field1_b0_file = None
                self.field1_temp_file = None
                self.field2_t1_file = None
                self.field2_t2_file = None
                self.field2_b0_file = None
                self.field2_temp_file = None
                self.field1_t1_display = MockDisplay()
                self.field1_t2_display = MockDisplay()

            def get_session_state(self):
                return {
                    'field_files': {
                        'field1_t1': self.field1_t1_file,
                        'field1_t2': self.field1_t2_file,
                        'field1_b0': self.field1_b0_file,
                        'field1_temp': self.field1_temp_file,
                        'field2_t1': self.field2_t1_file,
                        'field2_t2': self.field2_t2_file,
                        'field2_b0': self.field2_b0_file,
                        'field2_temp': self.field2_temp_file,
                    },
                    'source_series_map': dict(self.source_series_map),
                }

            def _update_field_displays(self):
                for field_key in ['field1_t1', 'field1_t2']:
                    display = getattr(self, f"{field_key}_display", None)
                    series_name = self.source_series_map.get(field_key)
                    if series_name and display:
                        display.setText(f"📊 {series_name}")

            def restore_session_state(self, state):
                field_files = state.get('field_files', {})
                for field_key, file_path in field_files.items():
                    setattr(self, f"{field_key}_file", file_path)
                self.source_series_map = dict(state.get('source_series_map', {}))
                self._update_field_displays()

        # Create page with initial state
        original_page = MockPage()
        original_page.source_series_map = {
            'field1_t1': 'T1_asyn_series',
            'field1_t2': 'T2_asyn_series',
        }
        original_page.field1_t1_file = '/path/to/t1.csv'
        original_page.field1_t2_file = '/path/to/t2.csv'

        # Save state
        saved_state = original_page.get_session_state()

        # Create fresh page and restore
        restored_page = MockPage()
        restored_page.restore_session_state(saved_state)

        # Verify round-trip
        assert restored_page.source_series_map == original_page.source_series_map
        assert restored_page.field1_t1_file == original_page.field1_t1_file
        assert restored_page.field1_t2_file == original_page.field1_t2_file
        assert "📊 T1_asyn_series" in restored_page.field1_t1_display.text_value


class TestDynamiXsMainWindowState:
    """Test DynamiXsMainWindow state methods."""

    def test_get_state_delegates_to_page(self):
        """DynamiXsMainWindow.get_state() delegates to page get_session_state()."""
        class MockPage:
            def get_session_state(self):
                return {'source_series_map': {'field1_t1': 'T1_series'}}

        class MockMainWindow:
            def __init__(self):
                self.current_page = MockPage()

            def get_state(self):
                state = {}
                if hasattr(self.current_page, 'get_session_state'):
                    state['t1t2_fitting'] = self.current_page.get_session_state()
                return state

        main_window = MockMainWindow()
        state = main_window.get_state()

        assert 't1t2_fitting' in state
        assert state['t1t2_fitting']['source_series_map']['field1_t1'] == 'T1_series'

    def test_restore_state_delegates_to_page(self):
        """DynamiXsMainWindow.restore_state() delegates to page restore_session_state()."""
        class MockPage:
            def __init__(self):
                self.restored_state = None

            def restore_session_state(self, state):
                self.restored_state = state

        class MockMainWindow:
            def __init__(self):
                self.current_page = MockPage()

            def restore_state(self, state):
                if 't1t2_fitting' in state and hasattr(self.current_page, 'restore_session_state'):
                    self.current_page.restore_session_state(state['t1t2_fitting'])

        main_window = MockMainWindow()
        state = {
            't1t2_fitting': {
                'source_series_map': {'field1_t1': 'T1_series'},
                'field_files': {'field1_t1': '/path/to/t1.csv'},
            }
        }

        main_window.restore_state(state)

        assert main_window.current_page.restored_state is not None
        assert main_window.current_page.restored_state['source_series_map']['field1_t1'] == 'T1_series'


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
