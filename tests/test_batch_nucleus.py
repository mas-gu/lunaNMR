# ABOUTME: `lunaNMR batch` takes each spectrum's nucleus from the integrator's own detection on
# ABOUTME: load, the one the PS2D config uses. Real 15N-HSQC and 13C methyl spectra.

"""
batch looked for the spectral ranges under attribute names the integrator does not have,
so it fell back to the filename: `2D_DNAJA1_0o0.ft` names no nucleus, and a 15N-HSQC ran
as 'unknown' with the `default` threshold and peak count instead of the 15N ones.
"""

from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parents[2] / "data_test"
HSQC_15N = DATA / "data_binding" / "DNAJA1_HSPA8" / "2D_DNAJA1_0o0.ft"
METHYL_13C = DATA / "methyl" / "sampleNBD" / "001_T2_NBD_3ms.ft"


@pytest.fixture
def processor(tmp_path):
    from lunaNMR.batch_processing.batch_processor import BatchProcessor
    return BatchProcessor(output_dir=tmp_path)


@pytest.mark.skipif(not HSQC_15N.exists(), reason="data_test not present")
def test_a_15n_hsqc_is_15n1h(processor):
    assert processor.spectrum_nucleus(HSQC_15N) == '15N1H'
    assert processor.integrator.auto_detected_nucleus == '15N'


@pytest.mark.skipif(not METHYL_13C.exists(), reason="data_test not present")
def test_a_13c_methyl_spectrum_is_13c1h(processor):
    assert processor.spectrum_nucleus(METHYL_13C) == '13C1H'
    assert processor.integrator.auto_detected_nucleus == '13C'


@pytest.mark.skipif(not HSQC_15N.exists(), reason="data_test not present")
def test_a_spectrum_that_does_not_load_is_unknown(processor, tmp_path):
    broken = tmp_path / "broken.ft"
    broken.write_bytes(b"not a spectrum")
    assert processor.spectrum_nucleus(broken) == 'unknown'


@pytest.mark.skipif(not HSQC_15N.exists(), reason="data_test not present")
def test_the_15n_entries_of_the_config_are_used(processor, tmp_path):
    """The 15N1H threshold and peak count, not the default ones."""
    from lunaNMR.batch_processing.batch_processor import BatchProcessor
    config = {"processing": {"sn_thresholds": {"15N1H": 1000.0, "default": 2.0},
                             "expected_peaks": {"15N1H": 7, "default": 150}}}
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    (spectra / HSQC_15N.name).symlink_to(HSQC_15N)
    stats = BatchProcessor(config, output_dir=tmp_path / "out").process_folder(spectra)
    run = stats['files'][0]
    assert (run['nucleus_type'], run['sn_threshold'], run['expected_peak_count']) == \
        ('15N1H', 1000.0, 7)
