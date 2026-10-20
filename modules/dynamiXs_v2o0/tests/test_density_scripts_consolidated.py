# ABOUTME: One single-field and one dual-field spectral density analysis, sharing one physics module;
# ABOUTME: J(0.87wH) is a factor of each, not a copy, so no formula can drift between variants.

from pathlib import Path

import pytest

DENSITY_DIR = Path(__file__).resolve().parents[1] / "dynamiXs_density_functions"


@pytest.mark.parametrize("name", ["ZZ_density", "ZZ_2fields_density", "ZZ_2fields_density_087",
                                  "ZZ_2fields_density_claude_rex_mcmc_error",
                                  "ZZ_multi_density_087", "ZZ_multi_2fields_density_087"])
def test_the_copy_is_gone(name):
    assert not (DENSITY_DIR / f"{name}.py").exists()


ANALYSES = ["ZZ_multi_density", "ZZ_multi_2fields_density"]


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


@pytest.fixture
def analyses():
    import sys
    sys.path.insert(0, str(DENSITY_DIR.parent))
    from dynamiXs_density_functions.ZZ_multi_density import ReducedSpectralDensityAnalysis
    from dynamiXs_density_functions.ZZ_multi_2fields_density import DualFieldSpectralDensityAnalysis
    return {'single': lambda f: ReducedSpectralDensityAnalysis(spectrometer_frequency=600.0,
                                                               omega_h_factor=f),
            'dual': lambda f: DualFieldSpectralDensityAnalysis(field1_freq=600.0, field2_freq=700.0,
                                                               omega_h_factor=f)}


@pytest.mark.parametrize("kind", ["single", "dual"])
@pytest.mark.parametrize("factor, column", [(1.0, 'JwH'), (0.87, 'JwH_087')])
def test_the_factor_names_the_jwh_column(analyses, kind, factor, column):
    assert analyses[kind](factor).jwh_column == column


@pytest.mark.parametrize("kind", ["single", "dual"])
def test_a_factor_with_no_column_name_is_refused(analyses, kind):
    with pytest.raises(ValueError, match="omega_h_factor"):
        analyses[kind](0.9)
