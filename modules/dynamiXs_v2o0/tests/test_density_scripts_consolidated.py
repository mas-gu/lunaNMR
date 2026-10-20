# ABOUTME: The spectral density functions are the analyses something runs, sharing one physics module.
# ABOUTME: Standalone copies that nothing imported are gone, so a formula cannot drift between them.

from pathlib import Path

import pytest

DENSITY_DIR = Path(__file__).resolve().parents[1] / "dynamiXs_density_functions"


@pytest.mark.parametrize("name", ["ZZ_density", "ZZ_2fields_density", "ZZ_2fields_density_087",
                                  "ZZ_2fields_density_claude_rex_mcmc_error"])
def test_the_uncalled_standalone_script_is_gone(name):
    assert not (DENSITY_DIR / f"{name}.py").exists()
