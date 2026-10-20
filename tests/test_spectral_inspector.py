# ABOUTME: Tests for the SpectralInspector data model (InspectorPeak, InspectorSpectrum, InspectorGroup)
# ABOUTME: Mostly display-free; a few canvas-rendering tests run under an offscreen Qt app.

import pytest


class TestInspectorPeak:
    """Tests for InspectorPeak dataclass and dict conversion."""

    def test_from_dict_standard_keys(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        d = {'ppm_x': 8.5, 'ppm_y': 120.3, 'assignment': 'A23-HN', 'quality': 'Excellent'}
        peak = InspectorPeak.from_dict(d)
        assert peak.ppm_x == pytest.approx(8.5)
        assert peak.ppm_y == pytest.approx(120.3)
        assert peak.assignment == 'A23-HN'
        assert peak.quality == 'Excellent'

    def test_from_dict_alias_keys(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        d = {'center_x': 7.2, 'center_y': 115.0}
        peak = InspectorPeak.from_dict(d)
        assert peak.ppm_x == pytest.approx(7.2)
        assert peak.ppm_y == pytest.approx(115.0)

    def test_from_dict_peak_x_alias(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        d = {'peak_x': 6.1, 'peak_y': 110.5, 'Quality': 'Good'}
        peak = InspectorPeak.from_dict(d)
        assert peak.ppm_x == pytest.approx(6.1)
        assert peak.quality == 'Good'

    def test_from_dict_position_aliases(self):
        """Main window peak-list columns (Position_X/Position_Y) are recognized."""
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        d = {'Assignment': 1, 'Position_X': 8.0407, 'Position_Y': 123.1998}
        peak = InspectorPeak.from_dict(d)
        assert peak.ppm_x == pytest.approx(8.0407)
        assert peak.ppm_y == pytest.approx(123.1998)
        assert peak.assignment == '1'

    def test_from_dict_nan_assignment_is_empty(self):
        """A blank Assignment cell (pandas NaN) must not become the string 'nan'."""
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        peak = InspectorPeak.from_dict(
            {'Assignment': float('nan'), 'Position_X': 8.0, 'Position_Y': 120.0}
        )
        assert peak.assignment == ''

    def test_from_dict_nan_quality_falls_back(self):
        """A blank Quality cell (NaN) must fall back to 'Unknown', not 'nan'."""
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        peak = InspectorPeak.from_dict(
            {'ppm_x': 8.0, 'ppm_y': 120.0, 'quality': float('nan')}
        )
        assert peak.quality == 'Unknown'

    def test_from_dict_nan_coord_defaults_zero(self):
        """A blank coordinate cell (NaN) falls back to 0.0, not a NaN peak."""
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        peak = InspectorPeak.from_dict(
            {'Position_X': 8.0, 'Position_Y': float('nan')}
        )
        assert peak.ppm_x == pytest.approx(8.0)
        assert peak.ppm_y == pytest.approx(0.0)

    def test_to_dict_standard_and_alias_keys(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        peak = InspectorPeak(ppm_x=8.0, ppm_y=118.0, assignment='B12-HN', quality='Good')
        d = peak.to_dict()
        assert d['ppm_x'] == pytest.approx(8.0)
        assert d['ppm_y'] == pytest.approx(118.0)
        assert d['center_x'] == pytest.approx(8.0)
        assert d['center_y'] == pytest.approx(118.0)
        assert d['assignment'] == 'B12-HN'
        assert d['quality'] == 'Good'
        assert 'peak_id' in d

    def test_unique_ids(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        assert InspectorPeak().peak_id != InspectorPeak().peak_id

    def test_from_dict_missing_coords_defaults_zero(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        peak = InspectorPeak.from_dict({})
        assert peak.ppm_x == pytest.approx(0.0)
        assert peak.ppm_y == pytest.approx(0.0)
        assert peak.quality == 'Unknown'

    def test_from_dict_zero_coord_preserved(self):
        """ppm_x=0.0 must not be silently discarded by falsy or-chain."""
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        peak = InspectorPeak.from_dict({'ppm_x': 0.0, 'ppm_y': 0.0})
        assert peak.ppm_x == pytest.approx(0.0)
        assert peak.ppm_y == pytest.approx(0.0)

    def test_from_dict_zero_coord_with_alias_fallback(self):
        """Explicit zero on first key must not fall through to alias key."""
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        # ppm_x is explicitly 0.0 — center_x should NOT override it
        peak = InspectorPeak.from_dict({'ppm_x': 0.0, 'center_x': 7.5})
        assert peak.ppm_x == pytest.approx(0.0)

    def test_from_dict_empty_assignment_not_overridden_by_alias(self):
        """Explicit empty assignment must not fall through to Assignment alias."""
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        peak = InspectorPeak.from_dict({'assignment': '', 'Assignment': 'B3'})
        assert peak.assignment == ''


class TestInspectorSpectrum:
    """Tests for InspectorSpectrum dataclass."""

    def test_defaults(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorSpectrum
        spec = InspectorSpectrum(display_name='test', file_path='/some/path')
        assert spec.visible is True
        assert spec.peaks == []
        assert spec.loaded is False
        assert spec.integrator is None

    def test_unique_ids(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorSpectrum
        assert InspectorSpectrum().spectrum_id != InspectorSpectrum().spectrum_id

    def test_default_color_is_valid_hex(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorSpectrum
        spec = InspectorSpectrum()
        assert spec.color.startswith('#')
        assert len(spec.color) == 7


class TestColorCycling:
    """Tests for color assignment cycling."""

    def test_cycles_through_palette(self):
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector, DEFAULT_COLORS
        inspector = SpectralInspector.__new__(SpectralInspector)
        inspector._color_index = 0
        colors = [inspector._next_color() for _ in range(len(DEFAULT_COLORS))]
        assert colors == DEFAULT_COLORS

    def test_wraps_around(self):
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector, DEFAULT_COLORS
        inspector = SpectralInspector.__new__(SpectralInspector)
        inspector._color_index = len(DEFAULT_COLORS) - 1
        last = inspector._next_color()
        first = inspector._next_color()
        assert last == DEFAULT_COLORS[-1]
        assert first == DEFAULT_COLORS[0]


class TestInspectorDataModelOperations:
    """Tests for group/spectrum management operations."""

    def _make_inspector(self):
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._color_index = 0
        return insp

    def test_all_spectra_empty(self):
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        insp = self._make_inspector()
        assert insp._all_spectra() == []

    def test_all_spectra_flattens_groups(self):
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = self._make_inspector()
        g1 = InspectorGroup(name='A')
        g2 = InspectorGroup(name='B')
        s1 = InspectorSpectrum(display_name='s1')
        s2 = InspectorSpectrum(display_name='s2')
        s3 = InspectorSpectrum(display_name='s3')
        g1.spectra.extend([s1, s2])
        g2.spectra.append(s3)
        insp.groups = [g1, g2]
        assert len(insp._all_spectra()) == 3

    def test_find_spectrum_across_groups(self):
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = self._make_inspector()
        g = InspectorGroup(name='A')
        spec = InspectorSpectrum(display_name='target')
        g.spectra.append(spec)
        insp.groups = [g]
        found = insp._find_spectrum(spec.spectrum_id)
        assert found is not None
        assert found.display_name == 'target'

    def test_find_spectrum_missing_returns_none(self):
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        insp = self._make_inspector()
        assert insp._find_spectrum('nonexistent-id') is None

    def test_peaks_from_dicts(self):
        """Loading a peak list from dicts populates InspectorPeak list."""
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorSpectrum, InspectorPeak
        spec = InspectorSpectrum(display_name='s')
        raw = [
            {'ppm_x': 8.1, 'ppm_y': 120.0, 'assignment': 'A1', 'quality': 'Excellent'},
            {'center_x': 7.5, 'center_y': 115.5, 'Assignment': 'B2'},
        ]
        spec.peaks = [InspectorPeak.from_dict(d) for d in raw]
        assert len(spec.peaks) == 2
        assert spec.peaks[0].assignment == 'A1'
        assert spec.peaks[1].ppm_x == pytest.approx(7.5)


class TestSpectralInspectorMutations:
    """Tests for group/spectrum add/remove operations on SpectralInspector."""

    def _make_inspector(self):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._color_index = 0
        insp._contour_levels = 10
        insp._contour_min_level = 0.05
        insp._contour_increment = 1.3
        insp._canvas = MagicMock()
        insp._toolbar_buttons = {}
        insp._status_bar_enabled = False
        return insp

    def test_add_group_appends_to_groups(self):
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector, InspectorGroup
        insp = self._make_inspector()
        result = insp._add_group('WT')
        assert len(insp.groups) == 1
        assert isinstance(result, InspectorGroup)
        assert result.name == 'WT'

    def test_add_group_returns_inspector_group(self):
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        insp = self._make_inspector()
        g = insp._add_group('Mutant')
        assert g.name == 'Mutant'
        assert g in insp.groups

    def test_add_multiple_groups(self):
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        insp = self._make_inspector()
        insp._add_group('A')
        insp._add_group('B')
        assert len(insp.groups) == 2
        assert insp.groups[0].name == 'A'
        assert insp.groups[1].name == 'B'

    def test_remove_group_removes_from_groups(self):
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        insp = self._make_inspector()
        g = insp._add_group('X')
        insp._remove_group(g.group_id)
        assert len(insp.groups) == 0

    def test_remove_group_missing_id_is_noop(self):
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        insp = self._make_inspector()
        insp._add_group('X')
        insp._remove_group('nonexistent-id')
        assert len(insp.groups) == 1

    def test_remove_spectrum_removes_from_group(self):
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorSpectrum
        )
        insp = self._make_inspector()
        g = insp._add_group('G')
        spec = InspectorSpectrum(display_name='s1')
        g.spectra.append(spec)
        insp._remove_spectrum(spec.spectrum_id)
        assert len(g.spectra) == 0

    def test_remove_spectrum_missing_id_is_noop(self):
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorSpectrum
        )
        insp = self._make_inspector()
        g = insp._add_group('G')
        spec = InspectorSpectrum(display_name='s1')
        g.spectra.append(spec)
        insp._remove_spectrum('nonexistent-id')
        assert len(g.spectra) == 1

    def test_load_spectrum_nonexistent_path_returns_none(self):
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        insp = self._make_inspector()
        g = insp._add_group('G')
        result = insp._load_spectrum_file('/nonexistent/path/to/spectrum.ft2', g)
        assert result is None

    def test_load_spectrum_populates_group(self):
        """Successful load adds InspectorSpectrum to group with loaded=True."""
        from unittest.mock import patch, MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector, InspectorSpectrum
        insp = self._make_inspector()
        g = insp._add_group('G')

        fake = MagicMock()
        fake.nmr_data = object()
        fake.ppm_x_axis = [1.0, 2.0]
        fake.ppm_y_axis = [100.0, 110.0]
        fake.load_nmr_file.return_value = True

        with patch('lunaNMR.core.core_integrator.EnhancedVoigtIntegrator', return_value=fake):
            result = insp._load_spectrum_file(__file__, g)

        assert result is not None
        assert isinstance(result, InspectorSpectrum)
        assert result.loaded is True
        assert result in g.spectra

    def test_load_spectrum_assigns_color_from_palette(self):
        """Each loaded spectrum gets the next palette color."""
        from unittest.mock import patch, MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector, DEFAULT_COLORS
        insp = self._make_inspector()
        g = insp._add_group('G')

        def make_fake():
            fake = MagicMock()
            fake.nmr_data = object()
            fake.ppm_x_axis = [1.0]
            fake.ppm_y_axis = [100.0]
            fake.load_nmr_file.return_value = True
            return fake

        with patch('lunaNMR.core.core_integrator.EnhancedVoigtIntegrator', side_effect=make_fake):
            s1 = insp._load_spectrum_file(__file__, g)
            s2 = insp._load_spectrum_file(__file__, g)

        assert s1 is not None and s2 is not None
        assert s1.color == DEFAULT_COLORS[0]
        assert s2.color == DEFAULT_COLORS[1]


class TestContourLevelCalculation:
    """Tests for InspectorCanvas._calculate_contour_levels (pure function, no display)."""

    def _calc(self, max_intensity, num_levels, min_level_fraction, increment):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas
        return InspectorCanvas._calculate_contour_levels(
            max_intensity, num_levels, min_level_fraction, increment
        )

    def test_first_level_is_fraction_of_max(self):
        levels = self._calc(max_intensity=100.0, num_levels=5, min_level_fraction=0.05, increment=1.3)
        assert levels[0] == pytest.approx(5.0)

    def test_levels_are_multiplicative(self):
        levels = self._calc(max_intensity=1000.0, num_levels=5, min_level_fraction=0.01, increment=2.0)
        assert levels[0] == pytest.approx(10.0)
        assert levels[1] == pytest.approx(20.0)
        assert levels[2] == pytest.approx(40.0)

    def test_levels_stop_before_exceeding_max(self):
        levels = self._calc(max_intensity=100.0, num_levels=20, min_level_fraction=0.5, increment=2.0)
        assert all(lv <= 100.0 for lv in levels)

    def test_num_levels_respected_when_range_allows(self):
        levels = self._calc(max_intensity=1e6, num_levels=8, min_level_fraction=0.001, increment=1.5)
        assert len(levels) == 8

    def test_empty_when_min_exceeds_max(self):
        levels = self._calc(max_intensity=10.0, num_levels=5, min_level_fraction=2.0, increment=1.3)
        assert levels == []

    def test_single_level(self):
        levels = self._calc(max_intensity=100.0, num_levels=1, min_level_fraction=0.05, increment=1.3)
        assert len(levels) == 1
        assert levels[0] == pytest.approx(5.0)


class TestSpectralInspectorContourParams:
    """Tests for contour parameter storage on SpectralInspector."""

    def _make_inspector(self):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._color_index = 0
        insp._contour_levels = 10
        insp._contour_min_level = 0.05
        insp._contour_increment = 1.3
        insp._canvas = MagicMock()
        insp._toolbar_buttons = {}
        insp._status_bar_enabled = False
        return insp

    def test_contour_params_stored_on_update(self):
        insp = self._make_inspector()
        insp._on_contour_update_requested(15, 0.03, 1.5)
        assert insp._contour_levels == 15
        assert insp._contour_min_level == pytest.approx(0.03)
        assert insp._contour_increment == pytest.approx(1.5)

    def test_visibility_change_updates_spectrum_model(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorGroup, InspectorSpectrum
        insp = self._make_inspector()
        g = insp._add_group('G')
        spec = InspectorSpectrum(display_name='s', visible=True)
        g.spectra.append(spec)
        insp._on_spectrum_visibility(spec.spectrum_id, False)
        assert spec.visible is False

    def test_color_change_updates_spectrum_model(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorGroup, InspectorSpectrum
        insp = self._make_inspector()
        g = insp._add_group('G')
        spec = InspectorSpectrum(display_name='s', color='#aabbcc')
        g.spectra.append(spec)
        insp._on_spectrum_color(spec.spectrum_id, '#ff0000')
        assert spec.color == '#ff0000'


class TestPeakManagement:
    """Tests for peak add/delete operations on SpectralInspector."""

    def _make_inspector(self):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector, InspectorGroup, InspectorSpectrum
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._color_index = 0
        insp._canvas = MagicMock()
        insp._status_bar_enabled = False
        g = InspectorGroup(name='G')
        spec = InspectorSpectrum(display_name='s1')
        g.spectra.append(spec)
        insp.groups.append(g)
        return insp, g, spec

    def test_add_peak_returns_inspector_peak(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        insp, g, spec = self._make_inspector()
        peak = insp._add_peak(spec.spectrum_id, 8.5, 120.3)
        assert isinstance(peak, InspectorPeak)
        assert peak.ppm_x == pytest.approx(8.5)
        assert peak.ppm_y == pytest.approx(120.3)

    def test_add_peak_appends_to_spectrum(self):
        insp, g, spec = self._make_inspector()
        insp._add_peak(spec.spectrum_id, 7.1, 115.0)
        assert len(spec.peaks) == 1

    def test_add_peak_nonexistent_spectrum_returns_none(self):
        insp, g, spec = self._make_inspector()
        result = insp._add_peak('nonexistent-id', 8.0, 120.0)
        assert result is None

    def test_delete_peak_removes_from_spectrum(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        insp, g, spec = self._make_inspector()
        peak = InspectorPeak(ppm_x=8.0, ppm_y=120.0)
        spec.peaks.append(peak)
        insp._delete_peak(spec.spectrum_id, peak.peak_id)
        assert len(spec.peaks) == 0

    def test_delete_peak_nonexistent_id_is_noop(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        insp, g, spec = self._make_inspector()
        peak = InspectorPeak(ppm_x=8.0, ppm_y=120.0)
        spec.peaks.append(peak)
        insp._delete_peak(spec.spectrum_id, 'bad-peak-id')
        assert len(spec.peaks) == 1


class TestPeakListParsing:
    """Tests for peak list file parsing and loading."""

    def test_parse_json_peak_list(self, tmp_path):
        import json
        from lunaNMR.gui.dialogs.spectral_inspector import _parse_peak_list_file
        data = [{'ppm_x': 8.1, 'ppm_y': 120.0, 'assignment': 'A1', 'quality': 'Good'}]
        f = tmp_path / 'peaks.json'
        f.write_text(json.dumps(data))
        result = _parse_peak_list_file(str(f))
        assert len(result) == 1
        assert result[0]['ppm_x'] == pytest.approx(8.1)

    def test_parse_csv_peak_list(self, tmp_path):
        from lunaNMR.gui.dialogs.spectral_inspector import _parse_peak_list_file
        f = tmp_path / 'peaks.csv'
        f.write_text('ppm_x,ppm_y,assignment,quality\n8.5,120.3,A1,Excellent\n7.2,115.0,B2,Good\n')
        result = _parse_peak_list_file(str(f))
        assert len(result) == 2
        assert float(result[0]['ppm_x']) == pytest.approx(8.5)

    def test_parse_txt_main_format(self, tmp_path):
        """A .txt in the main window's format (Assignment, Position_X, Position_Y)."""
        from lunaNMR.gui.dialogs.spectral_inspector import _parse_peak_list_file
        f = tmp_path / 'peaks.txt'
        # Header has spaces after commas, exactly like the real peak lists
        f.write_text('Assignment, Position_X, Position_Y\n'
                     '1,8.040740,123.199796\n'
                     '2,8.537622,123.723142\n')
        result = _parse_peak_list_file(str(f))
        assert len(result) == 2
        # Column names are stripped, matching the main window's loader
        assert float(result[0]['Position_X']) == pytest.approx(8.040740)
        assert float(result[1]['Position_Y']) == pytest.approx(123.723142)

    def test_load_txt_peak_list_populates_peaks(self, tmp_path):
        """End-to-end: a main-format .txt loads into a spectrum's peaks."""
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        g = InspectorGroup(name='G')
        spec = InspectorSpectrum(display_name='s')
        g.spectra.append(spec)
        insp.groups.append(g)

        f = tmp_path / 'peaks.txt'
        f.write_text('Assignment, Position_X, Position_Y\n'
                     '1,8.040740,123.199796\n'
                     '2,8.537622,123.723142\n')
        count = insp._load_peak_list_from_file(str(f), spec.spectrum_id)
        assert count == 2
        assert spec.peaks[0].ppm_x == pytest.approx(8.040740)
        assert spec.peaks[0].ppm_y == pytest.approx(123.199796)
        assert spec.peaks[1].assignment == '2'

    def test_load_peak_list_from_file_populates_peaks(self, tmp_path):
        import json
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._color_index = 0
        g = InspectorGroup(name='G')
        spec = InspectorSpectrum(display_name='s')
        g.spectra.append(spec)
        insp.groups.append(g)

        data = [
            {'ppm_x': 8.1, 'ppm_y': 120.0, 'assignment': 'A1'},
            {'center_x': 7.5, 'center_y': 115.5, 'Assignment': 'B2'},
        ]
        f = tmp_path / 'peaks.json'
        f.write_text(json.dumps(data))

        count = insp._load_peak_list_from_file(str(f), spec.spectrum_id)
        assert count == 2
        assert len(spec.peaks) == 2

    def test_load_peak_list_nonexistent_file_returns_zero(self):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        g = InspectorGroup(name='G')
        spec = InspectorSpectrum(display_name='s')
        g.spectra.append(spec)
        insp.groups.append(g)
        count = insp._load_peak_list_from_file('/nonexistent/peaks.json', spec.spectrum_id)
        assert count == 0


class TestFindNearestPeak:
    """Tests for InspectorCanvas._find_nearest_peak_in_spectra (pure function)."""

    def _make_spectra_with_peaks(self):
        from lunaNMR.gui.dialogs.spectral_inspector import (
            InspectorSpectrum, InspectorPeak, InspectorGroup
        )
        spec = InspectorSpectrum(display_name='s', visible=True)
        spec.peaks = [
            InspectorPeak(ppm_x=8.0, ppm_y=120.0),
            InspectorPeak(ppm_x=7.5, ppm_y=115.0),
        ]
        return [spec]

    def test_finds_nearest_peak(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas
        spectra = self._make_spectra_with_peaks()
        result = InspectorCanvas._find_nearest_peak_in_spectra(
            8.01, 120.5, spectra, xlim_range=3.0, ylim_range=30.0
        )
        assert result is not None
        peak, sid = result
        assert peak.ppm_x == pytest.approx(8.0)

    def test_returns_none_when_no_peaks(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas, InspectorSpectrum
        spec = InspectorSpectrum(display_name='empty', visible=True)
        result = InspectorCanvas._find_nearest_peak_in_spectra(
            8.0, 120.0, [spec], xlim_range=3.0, ylim_range=30.0
        )
        assert result is None

    def test_returns_none_when_beyond_threshold(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas
        spectra = self._make_spectra_with_peaks()
        # Click far from any peak (50% of axis range away)
        result = InspectorCanvas._find_nearest_peak_in_spectra(
            5.0, 90.0, spectra, xlim_range=3.0, ylim_range=30.0
        )
        assert result is None

    def test_skips_invisible_spectra(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas, InspectorSpectrum, InspectorPeak
        spec = InspectorSpectrum(display_name='s', visible=False)
        spec.peaks = [InspectorPeak(ppm_x=8.0, ppm_y=120.0)]
        result = InspectorCanvas._find_nearest_peak_in_spectra(
            8.0, 120.0, [spec], xlim_range=3.0, ylim_range=30.0
        )
        assert result is None


class TestPerSpectrumContour:
    """Tests for per-spectrum contour parameters."""

    def _make_inspector(self):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._color_index = 0
        insp._canvas = MagicMock()
        insp._toolbar_buttons = {}
        insp._status_bar_enabled = False
        insp._contour_levels = 10
        insp._contour_min_level = 0.05
        insp._contour_increment = 1.3
        g = InspectorGroup(name='G')
        s1 = InspectorSpectrum(display_name='s1')
        s2 = InspectorSpectrum(display_name='s2')
        g.spectra.extend([s1, s2])
        insp.groups.append(g)
        return insp, s1, s2

    def test_spectrum_has_contour_defaults(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorSpectrum
        spec = InspectorSpectrum()
        assert spec.contour_levels == 10
        assert spec.contour_min_level == pytest.approx(0.05)
        assert spec.contour_increment == pytest.approx(1.3)

    def test_update_spectrum_contour(self):
        insp, s1, s2 = self._make_inspector()
        insp._update_spectrum_contour(s1.spectrum_id, 15, 0.03, 1.5)
        assert s1.contour_levels == 15
        assert s1.contour_min_level == pytest.approx(0.03)
        assert s1.contour_increment == pytest.approx(1.5)
        # s2 unchanged
        assert s2.contour_levels == 10

    def test_update_spectrum_contour_nonexistent_is_noop(self):
        insp, s1, s2 = self._make_inspector()
        insp._update_spectrum_contour('bad-id', 5, 1.0, 2.0)  # no error raised
        assert s1.contour_levels == 10

    def test_apply_all_contour_updates_all_spectra(self):
        insp, s1, s2 = self._make_inspector()
        insp._on_contour_update_requested(20, 0.08, 1.6)
        assert s1.contour_levels == 20
        assert s2.contour_levels == 20
        assert s1.contour_min_level == pytest.approx(0.08)
        assert s2.contour_increment == pytest.approx(1.6)


class TestPropagateSettings:
    """Tests for settings propagation between spectra."""

    def _make_inspector(self, n=3):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._color_index = 0
        insp._canvas = MagicMock()
        insp._toolbar_buttons = {}
        insp._status_bar_enabled = False
        g = InspectorGroup(name='G')
        spectra = [InspectorSpectrum(display_name=f's{i}') for i in range(n)]
        g.spectra.extend(spectra)
        insp.groups.append(g)
        return insp, spectra

    def test_propagate_contour_levels(self):
        insp, (s1, s2, s3) = self._make_inspector()
        s1.contour_levels = 20
        insp._propagate_settings(s1.spectrum_id, [s2.spectrum_id, s3.spectrum_id],
                                  ['contour_levels'])
        assert s2.contour_levels == 20
        assert s3.contour_levels == 20

    def test_propagate_multiple_fields(self):
        insp, (s1, s2, s3) = self._make_inspector()
        s1.contour_levels = 15
        s1.contour_min_level = 0.03
        s1.color = '#ff0000'
        insp._propagate_settings(s1.spectrum_id, [s2.spectrum_id],
                                  ['contour_levels', 'contour_min_level', 'color'])
        assert s2.contour_levels == 15
        assert s2.contour_min_level == pytest.approx(0.03)
        assert s2.color == '#ff0000'

    def test_propagate_source_unchanged(self):
        insp, (s1, s2, s3) = self._make_inspector()
        s1.contour_levels = 20
        insp._propagate_settings(s1.spectrum_id, [s2.spectrum_id], ['contour_levels'])
        assert s1.contour_levels == 20  # source not modified

    def test_propagate_skips_missing_target_id(self):
        insp, (s1, s2, s3) = self._make_inspector()
        s1.contour_levels = 20
        insp._propagate_settings(s1.spectrum_id, ['nonexistent', s2.spectrum_id],
                                  ['contour_levels'])  # no error, s2 still updated
        assert s2.contour_levels == 20

    def test_propagate_source_not_found_is_noop(self):
        insp, (s1, s2, s3) = self._make_inspector()
        insp._propagate_settings('bad-id', [s2.spectrum_id], ['contour_levels'])
        assert s2.contour_levels == 10  # unchanged

    def test_propagate_ignores_unknown_fields(self):
        insp, (s1, s2, s3) = self._make_inspector()
        insp._propagate_settings(s1.spectrum_id, [s2.spectrum_id],
                                  ['nonexistent_field'])  # no error
        assert s2.contour_levels == 10  # unchanged


class TestInspectorGroup:
    """Tests for InspectorGroup dataclass."""

    def test_defaults(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorGroup
        group = InspectorGroup(name='WT')
        assert group.name == 'WT'
        assert group.spectra == []

    def test_unique_ids(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorGroup
        assert InspectorGroup().group_id != InspectorGroup().group_id

    def test_append_spectrum(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorGroup, InspectorSpectrum
        group = InspectorGroup(name='Mutant')
        s1 = InspectorSpectrum(display_name='spec1')
        s2 = InspectorSpectrum(display_name='spec2')
        group.spectra.append(s1)
        group.spectra.append(s2)
        assert len(group.spectra) == 2
        assert group.spectra[0].display_name == 'spec1'

    def test_find_spectrum_by_id(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorGroup, InspectorSpectrum
        group = InspectorGroup(name='Test')
        spec = InspectorSpectrum(display_name='myspec')
        group.spectra.append(spec)
        found = next((s for s in group.spectra if s.spectrum_id == spec.spectrum_id), None)
        assert found is not None
        assert found.display_name == 'myspec'


class TestSpectrumDragDrop:
    """Tests for drag-drop spectrum movement between groups (data model only)."""

    def _make_inspector(self):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._color_index = 0
        insp._contour_levels = 10
        insp._contour_min_level = 0.05
        insp._contour_increment = 1.3
        insp._canvas = MagicMock()
        insp._toolbar_buttons = {}
        insp._status_bar_enabled = False
        insp._library_panel = MagicMock()
        g1 = InspectorGroup(name='G1')
        g2 = InspectorGroup(name='G2')
        s1 = InspectorSpectrum(display_name='s1')
        g1.spectra.append(s1)
        insp.groups.extend([g1, g2])
        return insp, g1, g2, s1

    def test_move_spectrum_removes_from_source_group(self):
        insp, g1, g2, s1 = self._make_inspector()
        insp._on_spectrum_moved(s1.spectrum_id, g2.group_id)
        assert s1 not in g1.spectra

    def test_move_spectrum_adds_to_target_group(self):
        insp, g1, g2, s1 = self._make_inspector()
        insp._on_spectrum_moved(s1.spectrum_id, g2.group_id)
        assert s1 in g2.spectra

    def test_move_to_same_group_is_noop(self):
        insp, g1, g2, s1 = self._make_inspector()
        insp._on_spectrum_moved(s1.spectrum_id, g1.group_id)
        assert s1 in g1.spectra
        assert s1 not in g2.spectra

    def test_move_nonexistent_spectrum_is_noop(self):
        insp, g1, g2, s1 = self._make_inspector()
        insp._on_spectrum_moved('bad-id', g2.group_id)
        assert len(g1.spectra) == 1
        assert len(g2.spectra) == 0

    def test_move_to_nonexistent_group_is_noop(self):
        insp, g1, g2, s1 = self._make_inspector()
        insp._on_spectrum_moved(s1.spectrum_id, 'bad-group-id')
        assert s1 in g1.spectra


class TestActiveSpectrumTracking:
    """Tests that the controller tracks the active spectrum id."""

    def _make_inspector(self):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._active_spectrum_id = None
        insp._canvas = MagicMock()
        insp._status_bar_enabled = False
        g = InspectorGroup(name='G')
        spec = InspectorSpectrum(display_name='s1')
        g.spectra.append(spec)
        insp.groups.append(g)
        return insp, spec

    def test_active_id_tracked_on_change(self):
        insp, spec = self._make_inspector()
        insp._on_active_spectrum_changed(spec.spectrum_id)
        assert insp._active_spectrum_id == spec.spectrum_id

    def test_active_change_forwards_to_canvas(self):
        insp, spec = self._make_inspector()
        insp._on_active_spectrum_changed(spec.spectrum_id)
        insp._canvas.set_active_spectrum.assert_called_once_with(spec.spectrum_id)


class TestShiftPeakList:
    """Tests for the Shift Peak List model methods and handlers."""

    def _make_inspector(self, with_active=True):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum, InspectorPeak
        )
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._canvas = MagicMock()
        insp._status_bar_enabled = False
        g = InspectorGroup(name='G')
        spec = InspectorSpectrum(display_name='s1')
        spec.peaks = [
            InspectorPeak(ppm_x=8.0, ppm_y=120.0),
            InspectorPeak(ppm_x=7.5, ppm_y=115.0),
        ]
        g.spectra.append(spec)
        insp.groups.append(g)
        insp._active_spectrum_id = spec.spectrum_id if with_active else None
        return insp, spec

    def test_shift_moves_peaks(self):
        insp, spec = self._make_inspector()
        insp._shift_spectrum_peaks(spec.spectrum_id, 0.1, -1.0)
        assert spec.peaks[0].ppm_x == pytest.approx(8.1)
        assert spec.peaks[0].ppm_y == pytest.approx(119.0)
        assert spec.peaks[1].ppm_x == pytest.approx(7.6)

    def test_shift_tracks_applied_offset_per_peak(self):
        insp, spec = self._make_inspector()
        insp._shift_spectrum_peaks(spec.spectrum_id, 0.1, -1.0)
        insp._shift_spectrum_peaks(spec.spectrum_id, 0.05, 0.5)
        for peak in spec.peaks:
            assert peak.applied_shift_h == pytest.approx(0.15)
            assert peak.applied_shift_n == pytest.approx(-0.5)

    def test_shift_returns_peak_count(self):
        insp, spec = self._make_inspector()
        assert insp._shift_spectrum_peaks(spec.spectrum_id, 0.1, 0.0) == 2

    def test_shift_missing_spectrum_returns_zero(self):
        insp, spec = self._make_inspector()
        assert insp._shift_spectrum_peaks('bad-id', 0.1, 0.0) == 0

    def test_reset_restores_original_positions(self):
        insp, spec = self._make_inspector()
        insp._shift_spectrum_peaks(spec.spectrum_id, 0.1, -1.0)
        insp._reset_spectrum_shift(spec.spectrum_id)
        assert spec.peaks[0].ppm_x == pytest.approx(8.0)
        assert spec.peaks[0].ppm_y == pytest.approx(120.0)
        assert spec.peaks[1].ppm_x == pytest.approx(7.5)

    def test_reset_zeros_applied_offset(self):
        insp, spec = self._make_inspector()
        insp._shift_spectrum_peaks(spec.spectrum_id, 0.1, -1.0)
        insp._reset_spectrum_shift(spec.spectrum_id)
        for peak in spec.peaks:
            assert peak.applied_shift_h == pytest.approx(0.0)
            assert peak.applied_shift_n == pytest.approx(0.0)

    def test_reset_leaves_peaks_added_after_shift_untouched(self):
        """Regression: a peak added after a shift carries no applied shift, so
        Reset must leave it where the user placed it (not subtract the net shift)."""
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        insp, spec = self._make_inspector()
        insp._shift_spectrum_peaks(spec.spectrum_id, 0.1, -1.0)   # shift existing peaks
        # user adds a new peak at its true position AFTER the shift
        new_peak = insp._add_peak(spec.spectrum_id, 7.0, 110.0)
        insp._reset_spectrum_shift(spec.spectrum_id)
        # existing peaks restored...
        assert spec.peaks[0].ppm_x == pytest.approx(8.0)
        # ...but the later-added peak is untouched
        assert new_peak.ppm_x == pytest.approx(7.0)
        assert new_peak.ppm_y == pytest.approx(110.0)

    def test_apply_handler_shifts_active(self):
        insp, spec = self._make_inspector()
        insp._on_shift_apply_requested(0.1, -1.0)
        assert spec.peaks[0].ppm_x == pytest.approx(8.1)
        insp._canvas.update_plot.assert_called()

    def test_apply_handler_no_active_is_noop(self):
        insp, spec = self._make_inspector(with_active=False)
        insp._on_shift_apply_requested(0.1, -1.0)
        assert spec.peaks[0].ppm_x == pytest.approx(8.0)  # unchanged

    def test_reset_handler_restores_active(self):
        insp, spec = self._make_inspector()
        insp._on_shift_apply_requested(0.1, -1.0)
        insp._on_shift_reset_requested()
        assert spec.peaks[0].ppm_x == pytest.approx(8.0)
        assert spec.peaks[0].applied_shift_h == pytest.approx(0.0)


class TestPeakMarkersAndZoomForwarding:
    """Tests that peak-marker and reset-zoom handlers forward to the canvas."""

    def _make_inspector(self):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        insp = SpectralInspector.__new__(SpectralInspector)
        insp._canvas = MagicMock()
        insp._status_bar_enabled = False
        return insp

    def test_peak_markers_toggle_forwards(self):
        insp = self._make_inspector()
        insp._on_peak_markers_toggled(False)
        insp._canvas.set_peaks_visible.assert_called_once_with(False)

    def test_reset_zoom_forwards(self):
        insp = self._make_inspector()
        insp._on_reset_zoom_requested()
        insp._canvas.reset_zoom.assert_called_once()


class TestCanvasZoomPreservation:
    """Canvas rendering tests (offscreen Qt) for zoom-vs-empty-state behavior."""

    def _canvas(self):
        from PySide6.QtWidgets import QApplication
        QApplication.instance() or QApplication([])
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas
        return InspectorCanvas()

    def _spec(self):
        import numpy as np
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorSpectrum
        integ = MagicMock()
        integ.nmr_data = np.random.rand(50, 60)
        integ.ppm_x_axis = np.linspace(4.0, 12.0, 60)
        integ.ppm_y_axis = np.linspace(100.0, 132.0, 50)
        return InspectorSpectrum(display_name='s', integrator=integ, loaded=True)

    def test_full_range_on_first_draw(self):
        canvas = self._canvas()
        spec = self._spec()
        canvas.update_plot([spec])
        xlim = canvas._plot.axes.get_xlim()
        # Inverted NMR axis over the data range, NOT matplotlib's default (0,1)
        assert xlim[0] == pytest.approx(12.0, abs=0.5)
        assert xlim[1] == pytest.approx(4.0, abs=0.5)

    def test_view_refits_after_hide_then_show(self):
        """Regression: hiding all spectra then showing must refit to full range,
        not lock the view to the empty state's (0,1) range."""
        canvas = self._canvas()
        spec = self._spec()
        canvas.update_plot([spec])            # full range
        spec.visible = False
        canvas.update_plot([spec])            # empty state → axes at (0,1)
        spec.visible = True
        canvas.update_plot([spec])            # must refit, not preserve (0,1)
        xlim = canvas._plot.axes.get_xlim()
        assert not (xlim[0] == pytest.approx(0.0) and xlim[1] == pytest.approx(1.0))
        assert xlim[0] == pytest.approx(12.0, abs=0.5)
        assert xlim[1] == pytest.approx(4.0, abs=0.5)

    def test_user_zoom_preserved_across_visibility_toggle(self):
        """A real zoom must survive a redraw that keeps at least one spectrum visible."""
        canvas = self._canvas()
        s1 = self._spec()
        s2 = self._spec()
        canvas.update_plot([s1, s2])
        canvas._plot.axes.set_xlim(9.0, 7.0)   # user zooms in
        canvas._plot.axes.set_ylim(125.0, 115.0)
        s2.visible = False
        canvas.update_plot([s1, s2])           # s1 still visible → keep zoom
        xlim = canvas._plot.axes.get_xlim()
        assert xlim[0] == pytest.approx(9.0, abs=0.1)
        assert xlim[1] == pytest.approx(7.0, abs=0.1)

    def test_reset_zoom_refits_to_full(self):
        canvas = self._canvas()
        spec = self._spec()
        canvas.update_plot([spec])
        canvas._plot.axes.set_xlim(9.0, 7.0)   # zoomed in
        canvas.reset_zoom()
        xlim = canvas._plot.axes.get_xlim()
        assert xlim[0] == pytest.approx(12.0, abs=0.5)
        assert xlim[1] == pytest.approx(4.0, abs=0.5)

    def test_failed_contour_does_not_lock_blank_view(self):
        """A degenerate contour (increment<1 → non-increasing levels → ax.contour
        raises) must not lock the blank (0,1) view; a later good draw must refit."""
        canvas = self._canvas()
        spec = self._spec()
        canvas.update_plot([spec])                 # good draw
        spec.contour_increment = 0.5               # degenerate → every contour fails
        canvas.update_plot([spec])
        assert canvas._has_plotted is False
        spec.contour_increment = 1.3               # fixed
        canvas.update_plot([spec])
        xlim = canvas._plot.axes.get_xlim()
        assert not (xlim[0] == pytest.approx(0.0) and xlim[1] == pytest.approx(1.0))
        assert xlim[0] == pytest.approx(12.0, abs=0.5)


class TestZoomBox:
    """Canvas rendering tests (offscreen Qt) for the middle-drag zoom box."""

    def _canvas(self):
        from PySide6.QtWidgets import QApplication
        QApplication.instance() or QApplication([])
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas
        return InspectorCanvas()

    def _spec(self):
        import numpy as np
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorSpectrum
        integ = MagicMock()
        integ.nmr_data = np.random.rand(50, 60)
        integ.ppm_x_axis = np.linspace(4.0, 12.0, 60)
        integ.ppm_y_axis = np.linspace(100.0, 132.0, 50)
        return InspectorSpectrum(display_name='s', integrator=integ, loaded=True)

    def test_zoom_box_is_wired_to_handler(self):
        canvas = self._canvas()
        assert canvas._nav_handler.on_area_select == canvas._on_zoom_box
        assert canvas._nav_handler.on_rect_drag == canvas._on_zoom_rubberband

    def test_zoom_box_sets_independent_inverted_limits(self):
        canvas = self._canvas()
        canvas.update_plot([self._spec()])
        canvas._on_zoom_box(9.0, 118.0, 7.0, 122.0)
        xlim = canvas._plot.axes.get_xlim()
        ylim = canvas._plot.axes.get_ylim()
        # Inverted NMR orientation: high ppm first
        assert xlim[0] == pytest.approx(9.0) and xlim[1] == pytest.approx(7.0)
        assert ylim[0] == pytest.approx(122.0) and ylim[1] == pytest.approx(118.0)
        assert canvas._has_plotted is True

    def test_zoom_box_reshapes_aspect(self):
        """The box's x-span:y-span ratio need not match the data's — that's the point."""
        canvas = self._canvas()
        canvas.update_plot([self._spec()])
        canvas._on_zoom_box(9.0, 118.0, 8.5, 130.0)   # narrow x, wide y
        xspan = abs(canvas._plot.axes.get_xlim()[1] - canvas._plot.axes.get_xlim()[0])
        yspan = abs(canvas._plot.axes.get_ylim()[1] - canvas._plot.axes.get_ylim()[0])
        assert xspan == pytest.approx(0.5)
        assert yspan == pytest.approx(12.0)

    def test_zoom_box_ignores_degenerate_span(self):
        canvas = self._canvas()
        canvas.update_plot([self._spec()])
        before = canvas._plot.axes.get_xlim()
        canvas._on_zoom_box(8.0, 120.0, 8.0, 121.0)   # zero x-span
        assert canvas._plot.axes.get_xlim() == pytest.approx(before)

    def test_zoom_box_ignores_none_coords(self):
        canvas = self._canvas()
        canvas.update_plot([self._spec()])
        before = canvas._plot.axes.get_xlim()
        canvas._on_zoom_box(8.0, 120.0, None, 121.0)  # off-axes release
        assert canvas._plot.axes.get_xlim() == pytest.approx(before)

    def test_rubberband_draws_then_clears_on_finalize(self):
        canvas = self._canvas()
        canvas.update_plot([self._spec()])
        canvas._on_zoom_rubberband(9.0, 118.0, 7.0, 122.0)
        assert canvas._rubberband is not None
        assert canvas._rubberband in canvas._plot.axes.patches
        canvas._on_zoom_box(9.0, 118.0, 7.0, 122.0)
        assert canvas._rubberband is None


class TestExportPdf:
    """Canvas rendering tests (offscreen Qt) for vector-PDF export."""

    def _canvas(self):
        from PySide6.QtWidgets import QApplication
        QApplication.instance() or QApplication([])
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas
        return InspectorCanvas()

    def _spec(self):
        import numpy as np
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorSpectrum, InspectorPeak
        integ = MagicMock()
        integ.nmr_data = np.random.rand(50, 60)
        integ.ppm_x_axis = np.linspace(4.0, 12.0, 60)
        integ.ppm_y_axis = np.linspace(100.0, 132.0, 50)
        spec = InspectorSpectrum(display_name='S1', integrator=integ, loaded=True)
        spec.peaks = [InspectorPeak(ppm_x=8.0, ppm_y=120.0, assignment='23HN')]
        return spec

    def test_export_writes_valid_pdf(self, tmp_path):
        canvas = self._canvas()
        canvas.update_plot([self._spec()])
        out = tmp_path / 'export.pdf'
        canvas.export_pdf(str(out))
        assert out.exists()
        assert out.read_bytes()[:5] == b'%PDF-'   # a real PDF, not a raster image

    def test_export_uses_truetype_and_restores_rcparams(self, tmp_path):
        """During save pdf.fonttype must be 42 (editable text in Illustrator),
        and the global rcParam must be restored afterwards."""
        import matplotlib as mpl
        canvas = self._canvas()
        canvas.update_plot([self._spec()])

        mpl.rcParams['pdf.fonttype'] = 3          # simulate the un-editable default
        captured = {}
        orig_savefig = canvas._plot.figure.savefig

        def spy(*args, **kwargs):
            captured['fonttype'] = mpl.rcParams['pdf.fonttype']
            return orig_savefig(*args, **kwargs)

        canvas._plot.figure.savefig = spy
        try:
            canvas.export_pdf(str(tmp_path / 'x.pdf'))
        finally:
            canvas._plot.figure.savefig = orig_savefig

        assert captured['fonttype'] == 42          # TrueType while saving
        assert mpl.rcParams['pdf.fonttype'] == 3   # restored, no global side-effect

    def test_export_pdf_contains_editable_label_text(self, tmp_path):
        """The peak assignment should appear as recoverable text in the PDF."""
        canvas = self._canvas()
        canvas.update_plot([self._spec()])
        out = tmp_path / 'labels.pdf'
        canvas.export_pdf(str(out))
        raw = out.read_bytes()
        # Vector output: embedded fonts (text stays text) and NO rasterized image.
        # (`/Image` alone appears in the standard /ProcSet list, so check the real
        # rasterization marker `/Subtype /Image` instead.)
        assert b'/Type /Font' in raw or b'/Font' in raw
        assert b'/Subtype /Image' not in raw

    def test_export_restores_state_when_savefig_raises(self):
        """A failed save must still restore rcParams and the selection ring."""
        import matplotlib as mpl
        canvas = self._canvas()
        spec = self._spec()
        canvas.update_plot([spec])
        # select a peak so there is a highlight to restore
        canvas._selected_peak = spec.peaks[0]
        canvas._selected_spectrum_id = spec.spectrum_id

        mpl.rcParams['pdf.fonttype'] = 3            # a non-42 baseline to detect restore

        def boom(*args, **kwargs):
            raise IOError('disk full')

        canvas._plot.figure.savefig = boom
        with pytest.raises(IOError):
            canvas.export_pdf('/nonexistent/dir/out.pdf')

        assert mpl.rcParams['pdf.fonttype'] == 3    # restored despite the error
        assert canvas._highlight is not None        # selection ring redrawn


class TestMinLevelScaleButtons:
    """Offscreen-Qt tests for the ×2 / ÷2 buttons: emit a scale factor, don't set-all."""

    def _panel(self):
        from PySide6.QtWidgets import QApplication
        QApplication.instance() or QApplication([])
        from lunaNMR.gui.dialogs.spectral_inspector import SpectrumLibraryPanel
        return SpectrumLibraryPanel()

    def test_times2_emits_factor_2(self):
        panel = self._panel()
        captured = []
        panel.min_scale_requested.connect(lambda f: captured.append(f))
        panel._on_min_times2()
        assert captured == [pytest.approx(2.0)]

    def test_div2_emits_factor_half(self):
        panel = self._panel()
        captured = []
        panel.min_scale_requested.connect(lambda f: captured.append(f))
        panel._on_min_div2()
        assert captured == [pytest.approx(0.5)]

    def test_scale_does_not_touch_the_spinbox(self):
        """×2/÷2 is a per-spectrum relative scale — it must NOT set the global spinbox
        (which would flatten every spectrum to one value)."""
        panel = self._panel()
        panel._min_spin.setValue(0.05)
        panel._on_min_times2()
        assert panel._min_spin.value() == pytest.approx(0.05)   # unchanged


class TestMinScaleControl:
    """Controller-level ×2/÷2: scales each spectrum's OWN min level, preserving divergence."""

    def _make(self):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._canvas = MagicMock()
        insp._status_bar_enabled = False
        g = InspectorGroup(name='G')
        s1 = InspectorSpectrum(display_name='s1')
        s2 = InspectorSpectrum(display_name='s2')
        s1.contour_min_level = 0.05
        s2.contour_min_level = 0.10   # deliberately different
        g.spectra.extend([s1, s2])
        insp.groups.append(g)
        return insp, s1, s2

    def test_times2_scales_each_individually(self):
        insp, s1, s2 = self._make()
        insp._on_min_scale_requested(2.0)
        assert s1.contour_min_level == pytest.approx(0.10)
        assert s2.contour_min_level == pytest.approx(0.20)   # divergence preserved

    def test_div2_scales_each_down(self):
        insp, s1, s2 = self._make()
        insp._on_min_scale_requested(0.5)
        assert s1.contour_min_level == pytest.approx(0.025)
        assert s2.contour_min_level == pytest.approx(0.05)

    def test_clamps_at_ceiling(self):
        insp, s1, s2 = self._make()
        s1.contour_min_level = 0.8
        insp._on_min_scale_requested(2.0)
        assert s1.contour_min_level == pytest.approx(1.0)     # 1.6 clamped to 1.0

    def test_clamps_at_floor(self):
        insp, s1, s2 = self._make()
        s1.contour_min_level = 0.0001
        insp._on_min_scale_requested(0.5)
        assert s1.contour_min_level == pytest.approx(0.0001)

    def test_scale_refreshes_canvas(self):
        insp, s1, s2 = self._make()
        insp._on_min_scale_requested(2.0)
        insp._canvas.update_plot.assert_called()

    def test_scale_no_spectra_is_noop(self):
        insp, s1, s2 = self._make()
        insp.groups = []
        insp._on_min_scale_requested(2.0)   # must not raise
        insp._canvas.update_plot.assert_not_called()


class TestContourPropagateDialog:
    """Offscreen-Qt tests for the propagate-settings picker dialog."""

    def _specs(self, n=3):
        from PySide6.QtWidgets import QApplication
        QApplication.instance() or QApplication([])
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorSpectrum
        return [InspectorSpectrum(display_name=f's{i}') for i in range(n)]

    def _dlg(self, specs, active_id=None, on_apply=None):
        from lunaNMR.gui.dialogs.spectral_inspector import ContourPropagateDialog
        return ContourPropagateDialog(specs, active_id=active_id, on_apply=on_apply)

    def test_source_defaults_to_active(self):
        specs = self._specs()
        dlg = self._dlg(specs, active_id=specs[1].spectrum_id)
        assert dlg._source_combo.currentData() == specs[1].spectrum_id

    def test_targets_exclude_source(self):
        specs = self._specs()
        dlg = self._dlg(specs, active_id=specs[0].spectrum_id)
        targets = dlg.selected_targets()
        assert specs[0].spectrum_id not in targets
        assert set(targets) == {specs[1].spectrum_id, specs[2].spectrum_id}

    def test_default_fields_are_all_three(self):
        specs = self._specs()
        dlg = self._dlg(specs)
        assert dlg.selected_fields() == [
            'contour_levels', 'contour_min_level', 'contour_increment'
        ]

    def test_apply_invokes_callback_with_selections(self):
        specs = self._specs()
        captured = {}
        dlg = self._dlg(specs, active_id=specs[1].spectrum_id,
                        on_apply=lambda s, t, f: captured.update(source=s, targets=t, fields=f))
        dlg._apply()
        assert captured['source'] == specs[1].spectrum_id
        assert specs[1].spectrum_id not in captured['targets']
        assert 'contour_min_level' in captured['fields']

    def test_switching_source_rechecks_former_source(self):
        """Changing the source re-enables and re-checks the previous source as a target."""
        from PySide6.QtCore import Qt
        specs = self._specs()
        dlg = self._dlg(specs, active_id=specs[0].spectrum_id)
        # switch source from s0 to s1
        dlg._source_combo.setCurrentIndex(1)
        targets = dlg.selected_targets()
        assert specs[1].spectrum_id not in targets   # new source excluded
        assert specs[0].spectrum_id in targets       # former source re-checked


class TestApplyPropagation:
    """Controller-level propagation (model + refresh)."""

    def _make(self):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._canvas = MagicMock()
        insp._status_bar_enabled = False
        g = InspectorGroup(name='G')
        s1 = InspectorSpectrum(display_name='s1')
        s2 = InspectorSpectrum(display_name='s2')
        s1.contour_levels = 22
        s1.contour_min_level = 0.09
        g.spectra.extend([s1, s2])
        insp.groups.append(g)
        return insp, s1, s2

    def test_apply_propagation_copies_and_refreshes(self):
        insp, s1, s2 = self._make()
        insp._apply_propagation(s1.spectrum_id, [s2.spectrum_id],
                                ['contour_levels', 'contour_min_level'])
        assert s2.contour_levels == 22
        assert s2.contour_min_level == pytest.approx(0.09)
        insp._canvas.update_plot.assert_called()

    def test_apply_propagation_no_targets_is_noop(self):
        insp, s1, s2 = self._make()
        insp._apply_propagation(s1.spectrum_id, [], ['contour_levels'])
        assert s2.contour_levels == 10   # unchanged default

    def test_propagate_requested_needs_two_spectra(self):
        from unittest.mock import patch
        insp, s1, s2 = self._make()
        insp.groups[0].spectra = [s1]          # only one spectrum
        insp._active_spectrum_id = s1.spectrum_id
        with patch('lunaNMR.gui.dialogs.spectral_inspector.ContourPropagateDialog') as Dlg:
            insp._on_propagate_requested()
            Dlg.assert_not_called()

    def test_propagate_requested_opens_dialog_with_two(self):
        from unittest.mock import patch
        insp, s1, s2 = self._make()
        insp._active_spectrum_id = s1.spectrum_id
        with patch('lunaNMR.gui.dialogs.spectral_inspector.ContourPropagateDialog') as Dlg:
            insp._on_propagate_requested()
            Dlg.assert_called_once()


class TestPeakStateRoundtrip:
    """InspectorPeak.to_state / from_state (pure, no display)."""

    def test_roundtrip_preserves_all_fields(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        p = InspectorPeak(ppm_x=8.04, ppm_y=120.3, assignment='A23', quality='Good')
        p.applied_shift_h = 0.1
        p.applied_shift_n = -0.5
        p2 = InspectorPeak.from_state(p.to_state())
        assert p2.peak_id == p.peak_id
        assert p2.ppm_x == pytest.approx(8.04)
        assert p2.ppm_y == pytest.approx(120.3)
        assert p2.assignment == 'A23'
        assert p2.quality == 'Good'
        assert p2.applied_shift_h == pytest.approx(0.1)
        assert p2.applied_shift_n == pytest.approx(-0.5)

    def test_from_state_generates_id_when_missing(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorPeak
        p = InspectorPeak.from_state({'ppm_x': 1.0, 'ppm_y': 2.0})
        assert p.peak_id       # non-empty uuid


class TestInspectorStateRoundtrip:
    """Full get_state / load_state round-trip (offscreen Qt)."""

    def _inspector(self):
        from PySide6.QtWidgets import QApplication
        QApplication.instance() or QApplication([])
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        return SpectralInspector(parent=None)

    def _populate(self, w):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorSpectrum, InspectorPeak
        g = w._add_group('WT')
        w._library_panel.add_group(g)
        s = InspectorSpectrum(display_name='A', file_path='/tmp/nope.ft', color='#abcdef',
                              visible=False, contour_levels=22, contour_min_level=0.09,
                              contour_increment=1.5)
        s.peaks = [InspectorPeak(ppm_x=8.0, ppm_y=120.0, assignment='23HN')]
        g.spectra.append(s)
        w._library_panel.add_spectrum(s, g.group_id)
        w._add_toolbar_spectrum_button(s)
        w._active_spectrum_id = s.spectrum_id
        w._canvas._h_trace_on = True
        w._canvas._trace_scale = 0.4
        return s

    def test_state_is_json_serializable(self):
        import json
        w = self._inspector()
        self._populate(w)
        json.dumps(w.get_state())   # must not raise

    def test_roundtrip_preserves_model(self):
        w = self._inspector()
        s = self._populate(w)
        state = w.get_state()
        w2 = self._inspector()
        w2.load_state(state)
        s2 = w2._all_spectra()[0]
        assert [g.name for g in w2.groups] == ['WT']
        assert (s2.display_name, s2.color, s2.visible) == ('A', '#abcdef', False)
        assert s2.contour_levels == 22
        assert s2.contour_min_level == pytest.approx(0.09)
        assert s2.contour_increment == pytest.approx(1.5)
        assert s2.spectrum_id == s.spectrum_id            # stable id
        assert w2._active_spectrum_id == s.spectrum_id
        assert s2.peaks[0].assignment == '23HN'
        assert w2._canvas._h_trace_on is True
        assert w2._canvas._trace_scale == pytest.approx(0.4)
        assert w2._library_panel._h_trace_cb.isChecked() is True
        assert w2._library_panel._trace_scale_spin.value() == pytest.approx(0.4)
        assert w2._library_panel._trace_scale_slider.value() == 40

    def test_load_state_replaces_existing_content(self):
        w = self._inspector()
        old = w._add_group('Old')
        w._library_panel.add_group(old)
        w.load_state({'groups': [{'group_id': 'g1', 'name': 'New', 'spectra': []}]})
        assert [g.name for g in w.groups] == ['New']

    def test_missing_file_degrades_gracefully(self):
        w = self._inspector()
        w.load_state({'groups': [{'group_id': 'g', 'name': 'G', 'spectra': [
            {'spectrum_id': 's1', 'display_name': 'x', 'file_path': '/nonexistent.ft',
             'peaks': []}]}]})
        s = w._all_spectra()[0]
        assert s.loaded is False
        assert s.integrator is None
        assert s.display_name == 'x'      # still shown in the tree


class TestSpectrumReorder:
    """Reorder spectra within a group (bottom-to-top draw stacking)."""

    def _make(self):
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = SpectralInspector.__new__(SpectralInspector)
        g = InspectorGroup(name='G')
        s = [InspectorSpectrum(display_name=n) for n in ('A', 'B', 'C')]
        g.spectra.extend(s)
        insp.groups = [g]
        return insp, s

    def _order(self, insp):
        return [sp.display_name for sp in insp.groups[0].spectra]

    def test_move_down_swaps_with_next(self):
        insp, s = self._make()
        assert insp._reorder_spectrum(s[0].spectrum_id, +1) is True
        assert self._order(insp) == ['B', 'A', 'C']

    def test_move_up_swaps_with_previous(self):
        insp, s = self._make()
        assert insp._reorder_spectrum(s[2].spectrum_id, -1) is True
        assert self._order(insp) == ['A', 'C', 'B']

    def test_move_up_at_top_is_noop(self):
        insp, s = self._make()
        assert insp._reorder_spectrum(s[0].spectrum_id, -1) is False
        assert self._order(insp) == ['A', 'B', 'C']

    def test_move_down_at_bottom_is_noop(self):
        insp, s = self._make()
        assert insp._reorder_spectrum(s[2].spectrum_id, +1) is False
        assert self._order(insp) == ['A', 'B', 'C']

    def test_reorder_missing_spectrum_is_false(self):
        insp, s = self._make()
        assert insp._reorder_spectrum('bad-id', +1) is False

    def test_reorder_handler_updates_toolbar_order(self):
        from PySide6.QtWidgets import QApplication
        QApplication.instance() or QApplication([])
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorSpectrum
        )
        insp = SpectralInspector(parent=None)
        g = insp._add_group('G')
        insp._library_panel.add_group(g)
        a = InspectorSpectrum(display_name='A')
        b = InspectorSpectrum(display_name='B')
        for sp in (a, b):
            g.spectra.append(sp)
            insp._library_panel.add_spectrum(sp, g.group_id)
            insp._add_toolbar_spectrum_button(sp)
        insp._on_spectrum_reorder_requested(a.spectrum_id, +1)   # A toward front
        assert [sp.display_name for sp in insp._all_spectra()] == ['B', 'A']
        assert list(insp._toolbar_actions.keys()) == [b.spectrum_id, a.spectrum_id]


class TestProjectPersistence:
    """ProjectManager save/load/inventory for the Spectral Inspector bundle folder."""

    def _pm(self, main_window):
        from lunaNMR.utils.project_manager import ProjectManager
        return ProjectManager(main_window)

    def test_save_writes_state_json(self, tmp_path):
        from types import SimpleNamespace
        mw = SimpleNamespace(_spectral_inspector=None,
                             spectral_inspector_state={'schema': 1,
                                                       'groups': [{'group_id': 'g', 'name': 'G',
                                                                   'spectra': []}]})
        assert self._pm(mw)._save_spectral_inspector_state(tmp_path) is True
        assert (tmp_path / 'spectral_inspector' / 'state.json').exists()

    def test_save_no_content_prunes_stale_folder(self, tmp_path):
        from types import SimpleNamespace
        folder = tmp_path / 'spectral_inspector'
        folder.mkdir()
        (folder / 'state.json').write_text('{}')
        mw = SimpleNamespace(_spectral_inspector=None, spectral_inspector_state=None)
        assert self._pm(mw)._save_spectral_inspector_state(tmp_path) is False
        assert not folder.exists()

    def test_save_empty_groups_returns_false(self, tmp_path):
        from types import SimpleNamespace
        mw = SimpleNamespace(_spectral_inspector=None,
                             spectral_inspector_state={'groups': []})
        assert self._pm(mw)._save_spectral_inspector_state(tmp_path) is False

    def test_load_reads_state(self, tmp_path):
        import json
        from types import SimpleNamespace
        folder = tmp_path / 'spectral_inspector'
        folder.mkdir()
        (folder / 'state.json').write_text(json.dumps(
            {'schema': 1, 'groups': [{'group_id': 'g', 'name': 'G', 'spectra': []}]}))
        mw = SimpleNamespace(spectral_inspector_state=None)
        assert self._pm(mw)._load_spectral_inspector_state(tmp_path) is True
        assert mw.spectral_inspector_state['groups'][0]['name'] == 'G'

    def test_load_missing_returns_false_and_clears(self, tmp_path):
        from types import SimpleNamespace
        mw = SimpleNamespace(spectral_inspector_state='stale')
        assert self._pm(mw)._load_spectral_inspector_state(tmp_path) is False
        assert mw.spectral_inspector_state is None

    def test_save_captures_live_inspector(self, tmp_path):
        import json
        from types import SimpleNamespace
        from PySide6.QtWidgets import QApplication
        QApplication.instance() or QApplication([])
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector, InspectorSpectrum
        insp = SpectralInspector(parent=None)
        g = insp._add_group('G')
        insp._library_panel.add_group(g)
        g.spectra.append(InspectorSpectrum(display_name='A'))
        mw = SimpleNamespace(_spectral_inspector=insp, spectral_inspector_state=None)
        assert self._pm(mw)._save_spectral_inspector_state(tmp_path) is True
        data = json.loads((tmp_path / 'spectral_inspector' / 'state.json').read_text())
        assert data['groups'][0]['spectra'][0]['display_name'] == 'A'
        assert mw.spectral_inspector_state is not None      # captured onto main_window too

    def test_inventory_lists_spectral_inspector(self, tmp_path):
        from types import SimpleNamespace
        folder = tmp_path / 'spectral_inspector'
        folder.mkdir()
        (folder / 'state.json').write_text('{}')
        inv = self._pm(SimpleNamespace()).inventory(tmp_path)
        assert 'spectral_inspector' in [c['id'] for c in inv]


class TestTraceExtraction:
    """Pure cross-section extraction + scaling (no display needed)."""

    def _grid(self):
        import numpy as np
        data = np.array([[0, 1, 2, 3],
                         [4, 5, 6, 7],
                         [8, 9, 10, 11]], dtype=float)   # (ny=3, nx=4)
        ppm_y = [130.0, 120.0, 110.0]
        ppm_x = [10.0, 9.0, 8.0, 7.0]
        return data, ppm_x, ppm_y

    def test_horizontal_trace_picks_nearest_row(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas
        data, px, py = self._grid()
        xarr, row = InspectorCanvas._horizontal_trace(data, px, py, 119.0)  # nearest 120 → row 1
        assert list(row) == [4, 5, 6, 7]
        assert list(xarr) == [10, 9, 8, 7]

    def test_vertical_trace_picks_nearest_column(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas
        data, px, py = self._grid()
        yarr, col = InspectorCanvas._vertical_trace(data, px, py, 8.9)   # nearest 9 → col 1
        assert list(col) == [1, 5, 9]
        assert list(yarr) == [130, 120, 110]

    def test_offset_intensity_maps_amplitude(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas
        off = InspectorCanvas._offset_intensity([0.0, 2.0], max_int=4.0, scale=0.1,
                                                axis_range=30.0, baseline=100.0)
        assert off[0] == pytest.approx(100.0)        # zero intensity → baseline
        assert off[1] == pytest.approx(101.5)        # 100 + 0.1*30*(2/4)

    def test_offset_intensity_zero_max_is_safe(self):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas
        off = InspectorCanvas._offset_intensity([2.0], max_int=0.0, scale=1.0,
                                                axis_range=10.0, baseline=0.0)
        assert off[0] == pytest.approx(20.0)         # denom falls back to 1.0, no div-by-zero


class TestTraceRendering:
    """Offscreen-Qt tests for trace toggles + blitted rendering."""

    def _canvas(self):
        from PySide6.QtWidgets import QApplication
        QApplication.instance() or QApplication([])
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas
        return InspectorCanvas()

    def _spec(self, color='#ff0000'):
        import numpy as np
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorSpectrum
        integ = MagicMock()
        integ.nmr_data = np.random.rand(30, 40)
        integ.ppm_x_axis = np.linspace(4.0, 12.0, 40)
        integ.ppm_y_axis = np.linspace(100.0, 132.0, 30)
        return InspectorSpectrum(display_name='s', color=color, integrator=integ, loaded=True)

    def test_set_traces_updates_state(self):
        canvas = self._canvas()
        canvas.set_traces(horizontal=True)
        assert canvas._h_trace_on and not canvas._v_trace_on
        canvas.set_traces(vertical=True)
        assert canvas._h_trace_on and canvas._v_trace_on

    def test_set_trace_scale_updates_state(self):
        canvas = self._canvas()
        canvas.set_trace_scale(0.3)
        assert canvas._trace_scale == pytest.approx(0.3)

    def test_horizontal_only_draws_one_line_per_spectrum(self):
        canvas = self._canvas()
        canvas.update_plot([self._spec()])
        canvas.set_traces(horizontal=True, vertical=False)
        canvas._draw_traces_at(8.0, 116.0)
        visible = [ln for ln in canvas._trace_lines if ln.get_visible()]
        assert len(visible) == 1
        assert visible[0].get_color() == '#ff0000'      # trace color == spectrum color

    def test_both_traces_two_lines_per_spectrum(self):
        canvas = self._canvas()
        canvas.update_plot([self._spec()])
        canvas.set_traces(horizontal=True, vertical=True)
        canvas._draw_traces_at(8.0, 116.0)
        visible = [ln for ln in canvas._trace_lines if ln.get_visible()]
        assert len(visible) == 2

    def test_two_spectra_use_their_own_colors(self):
        canvas = self._canvas()
        canvas.update_plot([self._spec('#ff0000'), self._spec('#00ff00')])
        canvas.set_traces(horizontal=True, vertical=False)
        canvas._draw_traces_at(8.0, 116.0)
        colors = {ln.get_color() for ln in canvas._trace_lines if ln.get_visible()}
        assert colors == {'#ff0000', '#00ff00'}

    def test_toggling_off_hides_all_trace_lines(self):
        canvas = self._canvas()
        canvas.update_plot([self._spec()])
        canvas.set_traces(horizontal=True, vertical=True)
        canvas._draw_traces_at(8.0, 116.0)
        canvas.set_traces(horizontal=False, vertical=False)
        assert all(not ln.get_visible() for ln in canvas._trace_lines)

    def test_horizontal_trace_bumps_up_for_positive_peak(self):
        """On the inverted y-axis, a positive peak must plot ABOVE the baseline
        (toward the top = lower ppm_y), not below."""
        import numpy as np
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorSpectrum
        canvas = self._canvas()
        integ = MagicMock()
        data = np.zeros((30, 40))
        data[10, 20] = 1.0                       # single positive peak
        integ.nmr_data = data
        integ.ppm_x_axis = np.linspace(4.0, 12.0, 40)
        integ.ppm_y_axis = np.linspace(100.0, 132.0, 30)
        spec = InspectorSpectrum(display_name='s', integrator=integ, loaded=True)
        canvas.update_plot([spec])
        assert canvas._plot.axes.get_ylim()[0] > canvas._plot.axes.get_ylim()[1]  # inverted

        cy = float(integ.ppm_y_axis[10])
        canvas.set_traces(horizontal=True, vertical=False)
        canvas._draw_traces_at(float(integ.ppm_x_axis[20]), cy)
        line = [ln for ln in canvas._trace_lines if ln.get_visible()][0]
        ydata = np.asarray(line.get_ydata())
        assert ydata.min() < cy                  # peak bumps toward the top
        assert ydata.max() == pytest.approx(cy)  # zero-intensity points on the baseline

    def test_update_plot_resets_trace_pool(self):
        canvas = self._canvas()
        canvas.update_plot([self._spec()])
        canvas.set_traces(horizontal=True)
        canvas._draw_traces_at(8.0, 116.0)
        assert canvas._trace_lines            # created
        canvas.update_plot([self._spec()])    # ax.clear() removed the stale artists
        assert canvas._trace_lines == []      # pool dropped (no refs to removed artists)


class TestTracePanelSignals:
    """Offscreen-Qt tests for the 1D Traces panel controls."""

    def _panel(self):
        from PySide6.QtWidgets import QApplication
        QApplication.instance() or QApplication([])
        from lunaNMR.gui.dialogs.spectral_inspector import SpectrumLibraryPanel
        return SpectrumLibraryPanel()

    def test_horizontal_checkbox_emits(self):
        panel = self._panel()
        captured = []
        panel.horizontal_trace_toggled.connect(lambda b: captured.append(b))
        panel._h_trace_cb.setChecked(True)
        assert captured == [True]

    def test_vertical_checkbox_emits(self):
        panel = self._panel()
        captured = []
        panel.vertical_trace_toggled.connect(lambda b: captured.append(b))
        panel._v_trace_cb.setChecked(True)
        assert captured == [True]

    def test_scale_slider_maps_to_fraction(self):
        panel = self._panel()
        captured = []
        panel.trace_scale_changed.connect(lambda f: captured.append(f))
        panel._trace_scale_slider.setValue(50)
        assert captured[-1] == pytest.approx(0.5)

    def test_scale_slider_reaches_amplification(self):
        """The slider must go past 1× so weak traces can be blown up."""
        panel = self._panel()
        captured = []
        panel.trace_scale_changed.connect(lambda f: captured.append(f))
        panel._trace_scale_slider.setValue(panel._trace_scale_slider.maximum())
        assert captured[-1] >= 10.0

    def test_scale_spinbox_emits_and_moves_slider(self):
        panel = self._panel()
        captured = []
        panel.trace_scale_changed.connect(lambda f: captured.append(f))
        panel._trace_scale_spin.setValue(2.5)
        assert captured == [pytest.approx(2.5)]
        assert panel._trace_scale_slider.value() == 250

    def test_slider_updates_spinbox_once(self):
        panel = self._panel()
        captured = []
        panel.trace_scale_changed.connect(lambda f: captured.append(f))
        panel._trace_scale_slider.setValue(300)
        assert panel._trace_scale_spin.value() == pytest.approx(3.0)
        assert captured == [pytest.approx(3.0)]   # no echo from the spin box

    def test_set_scale_silently_syncs_both_without_emitting(self):
        panel = self._panel()
        captured = []
        panel.trace_scale_changed.connect(lambda f: captured.append(f))
        panel.set_trace_scale_silently(4.0)
        assert panel._trace_scale_slider.value() == 400
        assert panel._trace_scale_spin.value() == pytest.approx(4.0)
        assert captured == []


class TestTraceControllerForwarding:
    """Controller trace handlers forward to the canvas."""

    def _make(self):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
        insp = SpectralInspector.__new__(SpectralInspector)
        insp._canvas = MagicMock()
        insp._status_bar_enabled = False
        return insp

    def test_horizontal_forwards(self):
        insp = self._make()
        insp._on_horizontal_trace_toggled(True)
        insp._canvas.set_traces.assert_called_once_with(horizontal=True)

    def test_vertical_forwards(self):
        insp = self._make()
        insp._on_vertical_trace_toggled(False)
        insp._canvas.set_traces.assert_called_once_with(vertical=False)

    def test_scale_forwards(self):
        insp = self._make()
        insp._on_trace_scale_changed(0.3)
        insp._canvas.set_trace_scale.assert_called_once_with(0.3)


class TestExportPdfController:
    """Controller-level PDF export (guard, extension, error handling)."""

    def _make(self, with_spectrum=True):
        from unittest.mock import MagicMock
        from lunaNMR.gui.dialogs.spectral_inspector import (
            SpectralInspector, InspectorGroup, InspectorSpectrum
        )
        insp = SpectralInspector.__new__(SpectralInspector)
        insp.groups = []
        insp._canvas = MagicMock()
        insp._status_bar_enabled = False
        if with_spectrum:
            g = InspectorGroup(name='G')
            s = InspectorSpectrum(display_name='s', loaded=True, integrator=object())
            g.spectra.append(s)
            insp.groups.append(g)
        return insp

    def test_no_content_skips_dialog(self):
        from unittest.mock import patch
        insp = self._make(with_spectrum=False)
        with patch('lunaNMR.gui.dialogs.spectral_inspector.QFileDialog') as fd:
            insp._on_export_pdf()
            fd.getSaveFileName.assert_not_called()
        insp._canvas.export_pdf.assert_not_called()

    def test_appends_pdf_extension(self):
        from unittest.mock import patch
        insp = self._make()
        with patch('lunaNMR.gui.dialogs.spectral_inspector.QFileDialog.getSaveFileName',
                   return_value=('/tmp/out', '')):
            insp._on_export_pdf()
        insp._canvas.export_pdf.assert_called_once()
        assert insp._canvas.export_pdf.call_args[0][0] == '/tmp/out.pdf'

    def test_cancel_is_noop(self):
        from unittest.mock import patch
        insp = self._make()
        with patch('lunaNMR.gui.dialogs.spectral_inspector.QFileDialog.getSaveFileName',
                   return_value=('', '')):
            insp._on_export_pdf()
        insp._canvas.export_pdf.assert_not_called()

    def test_export_error_is_caught(self):
        from unittest.mock import patch
        insp = self._make()
        insp._canvas.export_pdf.side_effect = RuntimeError('boom')
        with patch('lunaNMR.gui.dialogs.spectral_inspector.QFileDialog.getSaveFileName',
                   return_value=('/tmp/out.pdf', '')):
            insp._on_export_pdf()   # must not raise
        insp._canvas.export_pdf.assert_called_once()
