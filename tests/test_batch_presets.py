# ABOUTME: A nucleus preset of `lunaNMR batch` detects with its own threshold and peak count.
# ABOUTME: The real CLI on the DNAJA1 15N-HSQC and a 13C methyl spectrum, read from the json envelope.

"""
batch reads each spectrum's threshold and peak count from the entry of its nucleus
(`15N1H`, `13C1H`), falling back to `default` only for an unknown one. The nucleus
presets set only `default`, so a 15N-HSQC run with --preset 15N1H detected with the
base config's 15N1H peak count (100), never the preset's (60).
"""

import contextlib
import io
import json
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parents[2] / "data_test"
HSQC_15N = DATA / "data_binding" / "DNAJA1_HSPA8" / "2D_DNAJA1_0o0.ft"
METHYL_13C = DATA / "methyl" / "sampleNBD" / "001_T2_NBD_3ms.ft"

pytestmark = pytest.mark.skipif(not HSQC_15N.exists(), reason="data_test not present")


def _batch_run(tmp_path, spectrum, preset):
    """The envelope's entry for `spectrum` after one `lunaNMR batch --preset` run."""
    from lunaNMR.cli import main
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    (spectra / spectrum.name).symlink_to(spectrum)
    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(io.StringIO()):
        code = main(["batch", "--preset", preset, str(spectra), "--output", str(tmp_path / "out"),
                     "--format", "json"])
    assert code == 0
    return json.loads(stdout.getvalue())["files"][0]


@pytest.mark.parametrize("spectrum, preset, expected", [
    (HSQC_15N, "15N1H", ("15N1H", 2.5, 60)),
    (METHYL_13C, "13C1H", ("13C1H", 2.0, 40)),
])
def test_a_nucleus_preset_detects_with_its_own_values(tmp_path, spectrum, preset, expected):
    run = _batch_run(tmp_path, spectrum, preset)
    assert (run["nucleus_type"], run["sn_threshold"], run["expected_peak_count"]) == expected


@pytest.mark.parametrize("preset, nucleus, expected", [
    ("15N1H", "15N1H", (2.5, 60)),
    ("13C1H", "13C1H", (2.0, 40)),
    ("1H", "1H", (2.0, 30)),
    ("conservative", "15N1H", (3.0, 30)),
    ("aggressive", "13C1H", (1.5, 150)),
])
def test_each_preset_sets_the_entry_of_the_nucleus_it_is_read_for(preset, nucleus, expected):
    from lunaNMR.batch_processing.config_manager import ConfigManager
    processing = ConfigManager().get_preset_config(preset)['processing']
    assert (processing['sn_thresholds'][nucleus], processing['expected_peaks'][nucleus]) == expected
