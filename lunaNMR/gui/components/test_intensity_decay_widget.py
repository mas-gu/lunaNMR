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

        delays, volumes, indices = collect_decay_data(spectra, 'A.2.ASP.H')

        assert delays == [50.0, 100.0, 200.0]
        assert volumes == [1000000, 800000, 500000]
        assert indices == [0, 1, 2]

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

        delays, volumes, indices = collect_decay_data(spectra, 'A.2.ASP.H')

        # Should only include spectra 0 and 2
        assert delays == [50.0, 200.0]
        assert volumes == [1000000, 500000]
        assert indices == [0, 2]

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

        delays, volumes, indices = collect_decay_data(spectra, 'A.2.ASP.H')

        # Should only include spectrum 0
        assert delays == [50.0]
        assert volumes == [1000000]
        assert indices == [0]


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
