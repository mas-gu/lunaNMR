# ABOUTME: Tests slim fit-result saving — surfaces dropped, region intensity optional.
# ABOUTME: Round-trips a real 2D fit result through ProjectManager save/load/reconstruct.

import json
import types
from pathlib import Path

import numpy as np
import pytest

from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
from lunaNMR.utils.project_manager import ProjectManager


def _make_integrator():
    """Synthetic 2D spectrum with descending ppm axes (NMR convention)."""
    integ = EnhancedVoigtIntegrator()
    ny, nx = 100, 200
    integ.ppm_y_axis = np.linspace(130.0, 110.0, ny)   # F1 / 15N
    integ.ppm_x_axis = np.linspace(10.0, 6.0, nx)       # F2 / 1H
    rng = np.random.default_rng(0)
    integ.nmr_data = rng.normal(0.0, 1.0, size=(ny, nx)).astype(np.float64)
    return integ


def _make_fit_result(integ):
    """Build a real 2D multi-peak result using the integrator's own primitives."""
    peaks = [
        {'pos_f1': 120.0, 'pos_f2': 8.00, 'x_ppm': 8.00, 'y_ppm': 120.0,
         'lw_gau_f1': 0.5, 'lw_gau_f2': 0.02, 'lw_lor_f1': 0.3, 'lw_lor_f2': 0.015,
         'intensity': 1.0e6, 'amplitude': 1.0e6, 'volume': 1.0e6, 'height': 1.0e6},
        {'pos_f1': 120.4, 'pos_f2': 8.03, 'x_ppm': 8.03, 'y_ppm': 120.4,
         'lw_gau_f1': 0.5, 'lw_gau_f2': 0.02, 'lw_lor_f1': 0.3, 'lw_lor_f2': 0.015,
         'intensity': 8.0e5, 'amplitude': 8.0e5, 'volume': 8.0e5, 'height': 8.0e5},
    ]
    region_2d = integ.extract_2d_region_for_overlap_group(peaks, radF1=1.0, radF2=0.1)
    assert region_2d is not None
    surface, individual, baseline = integ._reconstruct_2d_surface(region_2d, peaks)
    return {
        'assignment': 'A1',
        'method': '2d_simultaneous_multi_peak',
        'avg_r_squared': 0.95,
        'all_peaks': peaks,
        'region_2d': region_2d,
        'fitted_2d_surface': surface,
        'individual_surfaces': individual,
        'baseline': baseline,
        'intensity': 1.0e6, 'volume': 1.0e6, 'height': 1.0e6,
    }


@pytest.fixture
def saved_off(tmp_path):
    """Save one fit result with region data NOT embedded (default)."""
    integ = _make_integrator()
    result = _make_fit_result(integ)
    original = {
        'intensity': result['region_2d']['intensity'].copy(),
        'surface': result['fitted_2d_surface'].copy(),
        'individual': [s.copy() for s in result['individual_surfaces']],
    }
    mw = types.SimpleNamespace(last_fitting_results=[result], integrator=integ)
    pm = ProjectManager(mw)
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    pm._save_fit_results(bundle, embed_region_data=False)
    return bundle, integ, original


def test_off_save_writes_no_surface_arrays(saved_off):
    bundle, _, _ = saved_off
    arrays = bundle / "fit_results" / "arrays"
    files = list(arrays.glob("*.npz")) if arrays.exists() else []
    names = [f.name for f in files]
    assert not any("surface" in n for n in names), names
    assert not any("region_2d" in n for n in names), names


def test_off_fits_json_drops_surfaces_keeps_bounds(saved_off):
    bundle, _, _ = saved_off
    data = json.loads((bundle / "fit_results" / "fits.json").read_text())
    rec = data[0]
    assert 'fitted_2d_surface' not in rec
    assert 'individual_surfaces' not in rec
    assert 'all_peaks' in rec
    assert 'region_bounds' in rec and len(rec['region_bounds']) == 4


def test_off_roundtrip_reconstructs_from_spectrum(saved_off):
    bundle, integ, original = saved_off
    mw = types.SimpleNamespace(last_fitting_results=None, integrator=integ)
    pm = ProjectManager(mw)
    pm._load_fit_results(bundle)
    pm.reconstruct_fit_arrays()

    res = mw.last_fitting_results[0]
    np.testing.assert_array_equal(res['region_2d']['intensity'], original['intensity'])
    np.testing.assert_allclose(res['fitted_2d_surface'], original['surface'], rtol=1e-9, atol=1e-6)
    for got, exp in zip(res['individual_surfaces'], original['individual']):
        np.testing.assert_allclose(got, exp, rtol=1e-9, atol=1e-6)


def test_on_embed_then_load_without_spectrum(tmp_path):
    integ = _make_integrator()
    result = _make_fit_result(integ)
    original_surface = result['fitted_2d_surface'].copy()
    original_intensity = result['region_2d']['intensity'].copy()

    mw = types.SimpleNamespace(last_fitting_results=[result], integrator=integ)
    pm = ProjectManager(mw)
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    pm._save_fit_results(bundle, embed_region_data=True)

    arrays = bundle / "fit_results" / "arrays"
    assert any("region_2d" in f.name for f in arrays.glob("*.npz"))

    # Reload with an integrator that has NO spectrum loaded.
    empty = EnhancedVoigtIntegrator()
    mw2 = types.SimpleNamespace(last_fitting_results=None, integrator=empty)
    pm2 = ProjectManager(mw2)
    pm2._load_fit_results(bundle)
    pm2.reconstruct_fit_arrays()

    res = mw2.last_fitting_results[0]
    np.testing.assert_array_equal(res['region_2d']['intensity'], original_intensity)
    np.testing.assert_allclose(res['fitted_2d_surface'], original_surface, rtol=1e-9, atol=1e-6)


def test_off_without_spectrum_degrades_gracefully(saved_off):
    bundle, _, _ = saved_off
    empty = EnhancedVoigtIntegrator()  # nmr_data is None
    mw = types.SimpleNamespace(last_fitting_results=None, integrator=empty)
    pm = ProjectManager(mw)
    pm._load_fit_results(bundle)
    pm.reconstruct_fit_arrays()  # must not raise

    res = mw.last_fitting_results[0]
    assert 'fitted_2d_surface' not in res or res.get('fitted_2d_surface') is None
