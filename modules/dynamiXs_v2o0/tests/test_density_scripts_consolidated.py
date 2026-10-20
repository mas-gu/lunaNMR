# ABOUTME: The spectral density functions are the analyses something runs, sharing one physics module.
# ABOUTME: Standalone copies that nothing imported are gone, so a formula cannot drift between them.

from pathlib import Path

import pytest

DENSITY_DIR = Path(__file__).resolve().parents[1] / "dynamiXs_density_functions"


@pytest.mark.parametrize("name", ["ZZ_density", "ZZ_2fields_density", "ZZ_2fields_density_087",
                                  "ZZ_2fields_density_claude_rex_mcmc_error"])
def test_the_uncalled_standalone_script_is_gone(name):
    assert not (DENSITY_DIR / f"{name}.py").exists()


ANALYSES = ["ZZ_multi_density", "ZZ_multi_density_087", "ZZ_multi_2fields_density",
            "ZZ_multi_2fields_density_087"]


@pytest.mark.parametrize("name", ANALYSES)
def test_the_analysis_takes_its_physics_from_the_shared_module(name):
    import ast
    tree = ast.parse((DENSITY_DIR / f"{name}.py").read_text())
    defined = {t.id for node in tree.body if isinstance(node, ast.Assign)
               for t in node.targets if isinstance(t, ast.Name)}
    functions = {node.name for node in tree.body if isinstance(node, ast.FunctionDef)}
    assert not defined & {'GAMMA_H', 'GAMMA_N', 'REDUCED_PERM_VACUUM', 'REDUCED_PLANK', 'PI'}
    assert '_calculate_spectral_density_jit' not in functions
    imports = {alias.name for node in tree.body if isinstance(node, ast.ImportFrom)
               and node.module == 'dynamiXs_density_functions.spectral_density_physics'
               for alias in node.names}
    assert {'spectral_density', 'reduced_spectral_density_errors',
            'relaxation_deviations', 'overall_correlation_time'} <= imports
