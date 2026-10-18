# ABOUTME: A series fits each spectrum with the PS2D config of the nucleus it detected.
# ABOUTME: Real NBD methyl 13C-HSQC and DNAJA1 15N-HSQC spectra, run back to back in one process.

"""
Headless, the detected nucleus was logged and never applied: only the GUI called
set_ps2d_config. A 13C-HSQC was fitted with the 15N windows (radF1 0.4 against 0.15),
and 70% of its heights moved by more than 5%. The parallel workers also copied only eight
of the config's keys, so the 13C linewidth limits and the expert overrides never reached them.
"""

from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parents[2] / "data_test"
METHYL = DATA / "methyl" / "sampleNBD"
DNAJA1 = DATA / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not (METHYL.exists() and DNAJA1.exists()),
                                reason="data_test not present")


@pytest.fixture(autouse=True)
def fresh_config():
    from lunaNMR.core.ps2d_config import set_ps2d_config
    set_ps2d_config('15N')
    yield
    set_ps2d_config('15N')


def _series(tmp_path, folder, spectrum, peaks, name):
    from lunaNMR.cli import main
    spectra = tmp_path / name
    spectra.mkdir()
    (spectra / spectrum).symlink_to(folder / spectrum)
    assert main(["series", "--spectra", str(spectra), "--peaks", str(folder / peaks),
                 "--mode", "time", "--peak-source", "reference",
                 "--out", str(tmp_path / f"{name}_out")]) == 0


def test_a_13c_series_is_fitted_with_the_13c_config(tmp_path):
    from lunaNMR.core.ps2d_config import PS2DConfig, get_ps2d_config
    _series(tmp_path, METHYL, "001_T2_NBD_3ms.ft", "001_NBDADP_NBD.txt", "methyl")
    config = get_ps2d_config()
    assert config.nucleus_type == '13C'
    assert config.radF1 == PS2DConfig('13C').radF1


def test_a_15n_series_after_it_returns_to_15n(tmp_path):
    from lunaNMR.core.ps2d_config import PS2DConfig, get_ps2d_config
    _series(tmp_path, METHYL, "001_T2_NBD_3ms.ft", "001_NBDADP_NBD.txt", "methyl")
    _series(tmp_path, DNAJA1, "2D_DNAJA1_0o0.ft", "DNAJA1.txt", "amide")
    assert get_ps2d_config().nucleus_type == '15N'
    assert get_ps2d_config().radF1 == PS2DConfig('15N').radF1


def test_a_worker_receives_the_whole_config():
    from lunaNMR.core.parallel_voigt_processor import ParallelVoigtProcessor
    from lunaNMR.core.ps2d_config import get_ps2d_config, set_ps2d_config
    set_ps2d_config('13C')
    get_ps2d_config().max_drift_f1 = 0.123
    shipped = ParallelVoigtProcessor(object())._serialize_ps2d_config()
    assert shipped == dict(vars(get_ps2d_config()))
