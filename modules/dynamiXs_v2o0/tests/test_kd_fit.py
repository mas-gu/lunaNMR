# ABOUTME: Tests the Kd titration fitters: per-residue and global shared-Kd.
# ABOUTME: Recovers an injected Kd from synthetic CSP/intensity isotherms.

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_KD_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_Kd"
sys.path.insert(0, str(_KD_DIR))


L = np.array([0, 5, 10, 20, 40, 80, 160, 320], dtype=float)
P0 = 50.0


class TestPerResidueCsp:
    def test_recovers_kd_and_ddmax(self):
        from kd_models import csp_model
        from kd_fit import fit_residue_csp
        y = csp_model(L, dd_max=0.3, Kd=15.0, P=P0)
        r = fit_residue_csp(L, y, P0)
        assert r['success']
        assert r['Kd'] == pytest.approx(15.0, rel=1e-3)
        assert r['dd_max'] == pytest.approx(0.3, rel=1e-3)
        assert r['r_squared'] > 0.999

    def test_too_few_points_fails_gracefully(self):
        from kd_fit import fit_residue_csp
        r = fit_residue_csp([0.0, 10.0], [0.0, 0.1], P0)
        assert r['success'] is False


class TestPerResidueIntensity:
    def test_recovers_kd_with_intensity_loss(self):
        from kd_models import intensity_model
        from kd_fit import fit_residue_intensity
        y = intensity_model(L, baseline=1.0, amp=-0.7, Kd=25.0, P=P0)
        r = fit_residue_intensity(L, y, P0)
        assert r['success']
        assert r['Kd'] == pytest.approx(25.0, rel=1e-3)
        assert r['amp'] == pytest.approx(-0.7, rel=1e-3)


class TestGlobalKd:
    def test_shared_kd_across_residues(self):
        from kd_models import csp_model
        from kd_fit import fit_global_kd_csp
        kd_true = 18.0
        residues = {
            'A': csp_model(L, 0.30, kd_true, P0),
            'B': csp_model(L, 0.12, kd_true, P0),
            'C': csp_model(L, 0.45, kd_true, P0),
        }
        g = fit_global_kd_csp(residues, L, P0)
        assert g['success']
        assert g['Kd'] == pytest.approx(kd_true, rel=1e-3)
        assert g['dd_max']['A'] == pytest.approx(0.30, rel=1e-3)
        assert g['dd_max']['C'] == pytest.approx(0.45, rel=1e-3)


class TestRunAnalysis:
    def test_end_to_end_writes_json(self, tmp_path):
        from kd_models import csp_model
        from kd_fit import run_kd_analysis_with_params
        # synthetic tidy CSV: one mover with known Kd, 5 points
        pts = [0.0, 10.0, 25.0, 60.0, 150.0]
        ddH = csp_model(np.array(pts), 0.2, 20.0, P0)  # use CSP as 1H shift directly
        rows = []
        for p, d in zip(pts, ddH):
            rows.append((str(p), 'R1', 8.0 + d, 120.0, 1000.0, 2000.0))
        df = pd.DataFrame(rows, columns=['spectrum_name', 'assignment',
                                         'ppm_x', 'ppm_y', 'height', 'volume'])
        csv = tmp_path / 'series_analysis_tidy.csv'
        df.to_csv(csv, index=False)

        params = {
            'input_csv_file': str(csv),
            'output_dir': str(tmp_path),
            'output_prefix': 'kd_test',
            'concentrations': pts,
            'protein_conc': P0,
            'alpha': 0.14,
            'observables': ['csp'],
            'n_bootstrap': 0,
        }
        result = run_kd_analysis_with_params(params)
        assert result['n_fitted'] >= 1
        import json
        data = json.loads(Path(result['json_file']).read_text())
        r1 = next(f for f in data['fits'] if f['residue'] == 'R1')
        assert r1['csp']['Kd'] == pytest.approx(20.0, rel=1e-2)
        assert 'global' in data
