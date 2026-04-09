# ABOUTME: Tests for the SpectralInspector data model (InspectorPeak, InspectorSpectrum, InspectorGroup)
# ABOUTME: GUI tests are excluded — data model tests run without a display.

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
        insp._contour_min_pct = 5.0
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

    def _calc(self, max_intensity, num_levels, min_pct, increment):
        from lunaNMR.gui.dialogs.spectral_inspector import InspectorCanvas
        return InspectorCanvas._calculate_contour_levels(
            max_intensity, num_levels, min_pct, increment
        )

    def test_first_level_is_min_pct_of_max(self):
        levels = self._calc(max_intensity=100.0, num_levels=5, min_pct=5.0, increment=1.3)
        assert levels[0] == pytest.approx(5.0)

    def test_levels_are_multiplicative(self):
        levels = self._calc(max_intensity=1000.0, num_levels=5, min_pct=1.0, increment=2.0)
        assert levels[0] == pytest.approx(10.0)
        assert levels[1] == pytest.approx(20.0)
        assert levels[2] == pytest.approx(40.0)

    def test_levels_stop_before_exceeding_max(self):
        levels = self._calc(max_intensity=100.0, num_levels=20, min_pct=50.0, increment=2.0)
        assert all(lv <= 100.0 for lv in levels)

    def test_num_levels_respected_when_range_allows(self):
        levels = self._calc(max_intensity=1e6, num_levels=8, min_pct=0.1, increment=1.5)
        assert len(levels) == 8

    def test_empty_when_min_exceeds_max(self):
        levels = self._calc(max_intensity=10.0, num_levels=5, min_pct=200.0, increment=1.3)
        assert levels == []

    def test_single_level(self):
        levels = self._calc(max_intensity=100.0, num_levels=1, min_pct=5.0, increment=1.3)
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
        insp._contour_min_pct = 5.0
        insp._contour_increment = 1.3
        insp._canvas = MagicMock()
        insp._toolbar_buttons = {}
        insp._status_bar_enabled = False
        return insp

    def test_contour_params_stored_on_update(self):
        insp = self._make_inspector()
        insp._on_contour_update_requested(15, 3.0, 1.5)
        assert insp._contour_levels == 15
        assert insp._contour_min_pct == pytest.approx(3.0)
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
        insp._contour_min_pct = 5.0
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
        assert spec.contour_min_pct == pytest.approx(5.0)
        assert spec.contour_increment == pytest.approx(1.3)

    def test_update_spectrum_contour(self):
        insp, s1, s2 = self._make_inspector()
        insp._update_spectrum_contour(s1.spectrum_id, 15, 3.0, 1.5)
        assert s1.contour_levels == 15
        assert s1.contour_min_pct == pytest.approx(3.0)
        assert s1.contour_increment == pytest.approx(1.5)
        # s2 unchanged
        assert s2.contour_levels == 10

    def test_update_spectrum_contour_nonexistent_is_noop(self):
        insp, s1, s2 = self._make_inspector()
        insp._update_spectrum_contour('bad-id', 5, 1.0, 2.0)  # no error raised
        assert s1.contour_levels == 10

    def test_apply_all_contour_updates_all_spectra(self):
        insp, s1, s2 = self._make_inspector()
        insp._on_contour_update_requested(20, 8.0, 1.6)
        assert s1.contour_levels == 20
        assert s2.contour_levels == 20
        assert s1.contour_min_pct == pytest.approx(8.0)
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
        s1.contour_min_pct = 3.0
        s1.color = '#ff0000'
        insp._propagate_settings(s1.spectrum_id, [s2.spectrum_id],
                                  ['contour_levels', 'contour_min_pct', 'color'])
        assert s2.contour_levels == 15
        assert s2.contour_min_pct == pytest.approx(3.0)
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
        insp._contour_min_pct = 5.0
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
