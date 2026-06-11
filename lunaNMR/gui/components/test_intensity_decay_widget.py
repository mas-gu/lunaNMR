# ABOUTME: Unit tests for IntensityDecayWidget component
# ABOUTME: Tests delay extraction from spectrum names and plot data generation

import pytest

# Test delay extraction logic (standalone function, no Qt)


class TestDelayExtraction:
    """Test extracting delay values from spectrum names."""

    def test_numeric_spectrum_name(self):
        """Spectrum name is already a delay value."""
        from lunaNMR.gui.components.intensity_decay_widget import extract_delay_from_spectrum_name
        assert extract_delay_from_spectrum_name("50") == 50.0
        assert extract_delay_from_spectrum_name("100") == 100.0
        assert extract_delay_from_spectrum_name("1000") == 1000.0

    def test_numeric_with_suffix(self):
        """Spectrum name with duplicate suffix (e.g., 50_2)."""
        from lunaNMR.gui.components.intensity_decay_widget import extract_delay_from_spectrum_name
        assert extract_delay_from_spectrum_name("50_2") == 50.0
        assert extract_delay_from_spectrum_name("100_3") == 100.0

    def test_filename_with_ms_suffix(self):
        """Filenames with _XXms pattern."""
        from lunaNMR.gui.components.intensity_decay_widget import extract_delay_from_spectrum_name
        assert extract_delay_from_spectrum_name("T1asyn8_50ms") == 50.0
        assert extract_delay_from_spectrum_name("T1asyn8_100ms.ft") == 100.0
        assert extract_delay_from_spectrum_name("relax_250ms") == 250.0

    def test_filename_with_s_suffix(self):
        """Filenames with _Xs pattern (seconds)."""
        from lunaNMR.gui.components.intensity_decay_widget import extract_delay_from_spectrum_name
        assert extract_delay_from_spectrum_name("T1_1s") == 1000.0
        assert extract_delay_from_spectrum_name("T1_2s.ft") == 2000.0
        assert extract_delay_from_spectrum_name("T1_0.5s") == 500.0

    def test_unparseable_returns_none(self):
        """Unparseable names return None."""
        from lunaNMR.gui.components.intensity_decay_widget import extract_delay_from_spectrum_name
        assert extract_delay_from_spectrum_name("spectrum_001") is None
        assert extract_delay_from_spectrum_name("test") is None


class TestTitrationExtraction:
    """Titration points are read from the _<value> suffix via the o-for-dot convention."""

    def test_o_convention_suffix(self):
        from lunaNMR.gui.components.intensity_decay_widget import extract_titration_from_spectrum_name
        assert extract_titration_from_spectrum_name("sample_1o0.ft") == 1.0
        assert extract_titration_from_spectrum_name("sample_0o5.ft") == 0.5
        assert extract_titration_from_spectrum_name("titr_0o25.ft") == 0.25
        assert extract_titration_from_spectrum_name("sample_10o0") == 10.0

    def test_non_titration_returns_none(self):
        from lunaNMR.gui.components.intensity_decay_widget import extract_titration_from_spectrum_name
        # No o-decimal token -> not a titration name (avoids hijacking index/time files)
        assert extract_titration_from_spectrum_name("spectrum_001") is None
        assert extract_titration_from_spectrum_name("T1_50ms.ft") is None
        assert extract_titration_from_spectrum_name("reference") is None

    def test_index_style_not_treated_as_titration(self):
        """An integer index suffix must NOT be parsed as a titration point."""
        from lunaNMR.gui.components.intensity_decay_widget import extract_titration_from_spectrum_name
        assert extract_titration_from_spectrum_name("experiment_002.ft") is None

    def test_collect_decay_data_titration_mode(self):
        from lunaNMR.gui.components.intensity_decay_widget import collect_decay_data
        spectra = [
            {'name': 'titr_1o0.ft', 'fitted_peaks': [
                {'assignment': 'A.2.ASP.H', 'volume': 1000000}
            ]},
            {'name': 'titr_0o5.ft', 'fitted_peaks': [
                {'assignment': 'A.2.ASP.H', 'volume': 800000}
            ]},
        ]
        x_values, volumes, indices, mode = collect_decay_data(spectra, 'A.2.ASP.H')
        assert mode == 'titration'
        assert x_values == [1.0, 0.5]
        assert volumes == [1000000, 800000]


class TestDecayDataGeneration:
    """Test generating decay plot data from spectra."""

    def test_collect_decay_data_basic(self):
        """Collect volume vs delay from spectra."""
        from lunaNMR.gui.components.intensity_decay_widget import collect_decay_data

        # Mock spectra data
        spectra = [
            {'name': '50', 'fitted_peaks': [
                {'assignment': 'A.2.ASP.H', 'volume': 1000000}
            ]},
            {'name': '100', 'fitted_peaks': [
                {'assignment': 'A.2.ASP.H', 'volume': 800000}
            ]},
            {'name': '200', 'fitted_peaks': [
                {'assignment': 'A.2.ASP.H', 'volume': 500000}
            ]},
        ]

        delays, volumes, indices, mode = collect_decay_data(spectra, 'A.2.ASP.H')

        assert delays == [50.0, 100.0, 200.0]
        assert volumes == [1000000, 800000, 500000]
        assert indices == [0, 1, 2]
        assert mode == 'delay'

    def test_collect_decay_data_missing_peak(self):
        """Handle spectra where peak is missing."""
        from lunaNMR.gui.components.intensity_decay_widget import collect_decay_data

        spectra = [
            {'name': '50', 'fitted_peaks': [
                {'assignment': 'A.2.ASP.H', 'volume': 1000000}
            ]},
            {'name': '100', 'fitted_peaks': [
                {'assignment': 'A.3.VAL.H', 'volume': 900000}  # Different peak
            ]},
            {'name': '200', 'fitted_peaks': [
                {'assignment': 'A.2.ASP.H', 'volume': 500000}
            ]},
        ]

        delays, volumes, indices, mode = collect_decay_data(spectra, 'A.2.ASP.H')

        # Should only include spectra 0 and 2
        assert delays == [50.0, 200.0]
        assert volumes == [1000000, 500000]
        assert indices == [0, 2]
        assert mode == 'delay'

    def test_collect_decay_data_unparseable_delay(self):
        """Skip spectra with unparseable delay."""
        from lunaNMR.gui.components.intensity_decay_widget import collect_decay_data

        spectra = [
            {'name': '50', 'fitted_peaks': [
                {'assignment': 'A.2.ASP.H', 'volume': 1000000}
            ]},
            {'name': 'reference', 'fitted_peaks': [  # No delay in name
                {'assignment': 'A.2.ASP.H', 'volume': 900000}
            ]},
        ]

        delays, volumes, indices, mode = collect_decay_data(spectra, 'A.2.ASP.H')

        # Should only include spectrum 0
        assert delays == [50.0]
        assert volumes == [1000000]
        assert indices == [0]
        assert mode == 'delay'


class TestIndexExtraction:
    """Tests for spectrum index extraction from filenames."""

    def test_extract_index_basic(self):
        """Extract 3-digit index from filename."""
        from lunaNMR.gui.components.intensity_decay_widget import extract_index_from_spectrum_name

        assert extract_index_from_spectrum_name('03_2D_NR_ATP_ref_noCa_001.ft') == 1
        assert extract_index_from_spectrum_name('experiment_042.ucsf') == 42
        assert extract_index_from_spectrum_name('scan_0123.ft') == 123

    def test_extract_index_no_match(self):
        """Return None for filenames without 3+ digit suffix."""
        from lunaNMR.gui.components.intensity_decay_widget import extract_index_from_spectrum_name

        assert extract_index_from_spectrum_name('spectrum.ft') is None
        assert extract_index_from_spectrum_name('scan_01.ft') is None  # Only 2 digits
        assert extract_index_from_spectrum_name('') is None

    def test_collect_decay_data_index_fallback(self):
        """Fall back to index mode when no delays parseable."""
        from lunaNMR.gui.components.intensity_decay_widget import collect_decay_data

        spectra = [
            {'name': 'experiment_001.ft', 'fitted_peaks': [
                {'assignment': 'A.2.ASP.H', 'volume': 1000000}
            ]},
            {'name': 'experiment_002.ft', 'fitted_peaks': [
                {'assignment': 'A.2.ASP.H', 'volume': 800000}
            ]},
            {'name': 'experiment_003.ft', 'fitted_peaks': [
                {'assignment': 'A.2.ASP.H', 'volume': 500000}
            ]},
        ]

        x_values, volumes, indices, mode = collect_decay_data(spectra, 'A.2.ASP.H')

        assert x_values == [1.0, 2.0, 3.0]
        assert volumes == [1000000, 800000, 500000]
        assert indices == [0, 1, 2]
        assert mode == 'index'

    def test_delay_mode_takes_priority(self):
        """Delay mode should be used when delays are parseable."""
        from lunaNMR.gui.components.intensity_decay_widget import collect_decay_data

        # Delay patterns take priority over index extraction
        # Use proper delay format (pure numeric or _XXms pattern)
        spectra = [
            {'name': '50', 'fitted_peaks': [
                {'assignment': 'A.2.ASP.H', 'volume': 1000000}
            ]},
        ]

        x_values, volumes, indices, mode = collect_decay_data(spectra, 'A.2.ASP.H')

        assert mode == 'delay'
        assert x_values == [50.0]


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
