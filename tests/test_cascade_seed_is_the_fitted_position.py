# ABOUTME: A titration seeds each spectrum from where the previous one FITTED every peak, alone or
# ABOUTME: in a cluster; in a time series a cluster member keeps the position it was seeded from.

import contextlib
import io
from pathlib import Path

import pandas as pd
import pytest

from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
TITRATION_POINTS = ("0o0", "1o0", "3o0")


def _processor(titration=True):
    proc = MultiSpectrumProcessor.__new__(MultiSpectrumProcessor)
    proc.titration_tracking = titration
    proc.reference_peaks = pd.DataFrame({'Assignment': ['A'], 'Position_X': [8.000],
                                         'Position_Y': [120.00]})
    return proc


class TestATitrationSeedIsTheFittedPosition:

    def test_a_cluster_member_seeds_from_its_fit_not_from_its_own_seed(self):
        """A cluster member's result keeps the position it was seeded from in peak_x/peak_y
        and its fitted position in center_x/center_y."""
        member = {'assignment': 'A', 'success': True, 'peak_x': 8.000, 'peak_y': 120.00,
                  'center_x': 8.015, 'center_y': 120.20, 'height': 5.0, 'volume': 9.0}
        seed = _processor()._convert_fitted_results_to_reference_df([member])
        assert seed['Position_X'].iloc[0] == 8.015
        assert seed['Position_Y'].iloc[0] == 120.20

    def test_a_peak_fitted_alone_seeds_from_its_fit(self):
        alone = {'assignment': 'A', 'success': True, 'peak_x': 8.015, 'peak_y': 120.20,
                 'center_x': 8.015, 'center_y': 120.20, 'height': 5.0, 'volume': 9.0}
        seed = _processor()._convert_fitted_results_to_reference_df([alone])
        assert seed['Position_X'].iloc[0] == 8.015
        assert seed['Position_Y'].iloc[0] == 120.20

    def test_a_result_without_a_fitted_centre_seeds_from_the_position_it_carries(self):
        bare = {'assignment': 'A', 'success': True, 'peak_x': 8.010, 'peak_y': 120.10}
        seed = _processor()._convert_fitted_results_to_reference_df([bare])
        assert seed['Position_X'].iloc[0] == 8.010
        assert seed['Position_Y'].iloc[0] == 120.10


class TestATimeSeriesHoldsAClusterMemberWhereItWasSeeded:
    """Relaxation peaks do not move. Measured in cascade mode, letting cluster members
    follow their fits moved a 600 MHz T1 series further from the reference-mode values
    (median 0.95% to 1.95%) at unchanged fit errors."""

    def test_a_cluster_member_keeps_the_position_it_was_seeded_from(self):
        member = {'assignment': 'A', 'success': True, 'peak_x': 8.000, 'peak_y': 120.00,
                  'center_x': 8.015, 'center_y': 120.20, 'height': 5.0, 'volume': 9.0}
        seed = _processor(titration=False)._convert_fitted_results_to_reference_df([member])
        assert seed['Position_X'].iloc[0] == 8.000
        assert seed['Position_Y'].iloc[0] == 120.00


@pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")
def test_every_peak_of_a_real_titration_is_seeded_where_it_was_fitted(tmp_path, monkeypatch):
    """Three points of the DNAJA1 titration: each seed list handed to the next spectrum holds,
    for every peak that seeds from its own fit, the position that fit reported."""
    handed_on = []
    convert = MultiSpectrumProcessor._convert_fitted_results_to_reference_df

    def recording(self, fitted_results):
        seeds = convert(self, fitted_results)
        handed_on.append((fitted_results, seeds))
        return seeds

    monkeypatch.setattr(MultiSpectrumProcessor, '_convert_fitted_results_to_reference_df', recording)
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    for point in TITRATION_POINTS:
        source = next(DNAJA1.glob(f"2D_*_{point}.ft"))
        (spectra / source.name).symlink_to(source)
    from lunaNMR.cli import main
    with contextlib.redirect_stdout(io.StringIO()):
        assert main(["series", "--spectra", str(spectra), "--peaks", str(DNAJA1 / "DNAJA1.txt"),
                     "--mode", "titration", "--out", str(tmp_path / "out")]) == 0

    assert len(handed_on) == len(TITRATION_POINTS) - 1
    clustered = 0
    for fitted_results, seeds in handed_on:
        assert len(seeds) == len(fitted_results)
        for result, (_, seed) in zip(fitted_results, seeds.iterrows()):
            if not (result and result.get('success') and 'center_x' in result):
                continue
            clustered += result.get('method') == '2d_simultaneous_multi_peak'
            assert seed['Position_X'] == pytest.approx(result['center_x'], abs=1e-9), result['assignment']
            assert seed['Position_Y'] == pytest.approx(result['center_y'], abs=1e-9), result['assignment']
    assert clustered > 0
