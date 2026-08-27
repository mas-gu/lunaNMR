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
    def test_recovers_decay_constant_I0_and_plateau(self):
        from kd_models import intensity_decay
        from kd_fit import fit_residue_intensity
        y = intensity_decay(L, I0=1.0, I_inf=0.05, Kd=25.0)   # plateau above zero
        r = fit_residue_intensity(L, y, P0)
        assert r['success']
        assert r['Kd'] == pytest.approx(25.0, rel=1e-2)
        assert r['I0'] == pytest.approx(1.0, rel=1e-2)
        assert r['I_inf'] == pytest.approx(0.05, rel=1e-2)    # tracks the floor, not 0

    def test_bootstrap_gives_finite_errors(self):
        from kd_models import intensity_decay
        from kd_fit import fit_residue_intensity
        y = intensity_decay(L, I0=1.0, I_inf=0.05, Kd=25.0) + np.array(
            [0.02, -0.01, 0.0, 0.01, -0.02, 0.0, 0.01, -0.01])  # mild noise
        r = fit_residue_intensity(L, y, n_bootstrap=200)
        assert r['success']
        assert np.isfinite(r['Kd_err']) and r['Kd_err'] > 0
        assert np.isfinite(r['I0_err']) and np.isfinite(r['I_inf_err'])


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

    def test_shared_intensity_decay_across_residues(self):
        from kd_models import intensity_decay
        from kd_fit import fit_global_kd_intensity
        kd_true = 22.0
        residues = {
            'A': intensity_decay(L, 1.00, 0.05, kd_true),
            'B': intensity_decay(L, 0.80, 0.20, kd_true),
            'C': intensity_decay(L, 1.20, 0.00, kd_true),
        }
        g = fit_global_kd_intensity(residues, L)
        assert g['success']
        assert g['Kd'] == pytest.approx(kd_true, rel=1e-3)
        assert g['I0']['A'] == pytest.approx(1.00, rel=1e-3)
        assert g['I_inf']['B'] == pytest.approx(0.20, rel=1e-3)
        assert g['n_residues'] == 3
        # a shared-Kd standard error is reported and finite for clean data
        assert g['Kd_err'] is not None and np.isfinite(g['Kd_err']) and g['Kd_err'] >= 0

    def test_global_intensity_bootstrap_gives_finite_error(self):
        # Residual-resampling bootstrap (same convention as fit_residue_intensity's
        # own n_bootstrap), scaled to the joint multi-residue fit: pool every residue's
        # residuals, resample, refit the whole joint problem, collect the shared Kd.
        from kd_models import intensity_decay
        from kd_fit import fit_global_kd_intensity
        kd_true = 22.0
        noise = np.array([0.02, -0.01, 0.0, 0.01, -0.02, 0.0, 0.01, -0.01])
        residues = {
            'A': intensity_decay(L, 1.00, 0.05, kd_true) + noise,
            'B': intensity_decay(L, 0.80, 0.20, kd_true) - noise,
            'C': intensity_decay(L, 1.20, 0.00, kd_true) + noise[::-1],
        }
        g = fit_global_kd_intensity(residues, L, n_bootstrap=30)
        assert g['success']
        assert np.isfinite(g['Kd_err']) and g['Kd_err'] > 0

    def test_bootstrap_global_intensity_runs_without_crashing(self):
        # Regression test for a parameter-shadowing bug: the bootstrap-count
        # argument `n` was being overwritten by the `for i, n in enumerate(names)`
        # loop, so `for _ in range(n)` later raised TypeError on a residue-name
        # string and the bare except swallowed it, always returning nan. Calling
        # the helper directly (not through the try/except in fit_global_kd_intensity)
        # would have surfaced that as a nan here, before the fix.
        from kd_models import intensity_decay
        from kd_fit import fit_global_kd_intensity, _bootstrap_global_intensity
        kd_true = 22.0
        residues = {
            'A': intensity_decay(L, 1.00, 0.05, kd_true),
            'B': intensity_decay(L, 0.80, 0.20, kd_true),
            'C': intensity_decay(L, 1.20, 0.00, kd_true),
        }
        g = fit_global_kd_intensity(residues, L)
        names = list(residues.keys())
        series = {n: np.asarray(residues[n], dtype=float) for n in names}
        masks = {n: ~np.isnan(series[n]) for n in names}
        popt = np.array([g['Kd']] + [v for n in names for v in (g['I0'][n], g['I_inf'][n])])
        lo = [1e-9] + [0.0] * (2 * len(names))
        hi = [np.inf] * (1 + 2 * len(names))
        result = _bootstrap_global_intensity(L, series, masks, names, popt, lo, hi, 20)
        assert np.isfinite(result)

    def test_global_intensity_bootstrap_falls_back_to_analytic_on_failure(self, monkeypatch):
        # If every bootstrap refit fails, Kd_err must fall back to the analytic
        # covariance error rather than surface a non-finite value.
        from kd_models import intensity_decay
        import kd_fit
        kd_true = 22.0
        residues = {
            'A': intensity_decay(L, 1.00, 0.05, kd_true),
            'B': intensity_decay(L, 0.80, 0.20, kd_true),
            'C': intensity_decay(L, 1.20, 0.00, kd_true),
        }
        monkeypatch.setattr(kd_fit, '_bootstrap_global_intensity',
                            lambda *a, **k: float('nan'))
        g = kd_fit.fit_global_kd_intensity(residues, L, n_bootstrap=10)
        assert g['success']
        assert np.isfinite(g['Kd_err'])

    def test_global_csp_reports_kd_err(self):
        from kd_models import csp_model
        from kd_fit import fit_global_kd_csp
        residues = {n: csp_model(L, dd, 18.0, P0)
                    for n, dd in (('A', 0.30), ('B', 0.12), ('C', 0.45))}
        g = fit_global_kd_csp(residues, L, P0)
        assert g['success']
        assert 'Kd_err' in g

    def test_global_intensity_fewer_than_two_residues_fails(self):
        from kd_models import intensity_decay
        from kd_fit import fit_global_kd_intensity
        g = fit_global_kd_intensity({'A': intensity_decay(L, 1.0, 0.05, 20.0)}, L)
        assert not g['success']

    def test_global_csp_kd_bounded_to_titration_range(self):
        # A titration from 5-320 has essentially no power to resolve a Kd of 1e-8: the
        # model looks saturated at every single point, indistinguishable from Kd=1e-6 or
        # Kd=1e-3. Without a bound the optimizer can run the shared Kd to that nonsensical
        # extreme (with a misleadingly small formal error); it must stay within one decade
        # of the tested concentration range instead.
        from kd_models import csp_model
        from kd_fit import fit_global_kd_csp
        residues = {n: csp_model(L, dd, 1e-8, P0) for n, dd in (('A', 0.30), ('B', 0.12))}
        g = fit_global_kd_csp(residues, L, P0)
        assert g['success']
        nonzero = L[L > 0]
        lo = nonzero.min() / 10.0
        assert g['Kd'] >= lo
        assert g['reliable'] is False

    def test_global_csp_kd_bounded_above_titration_range(self):
        # The mirror case: essentially flat/no-signal data (true Kd far above the highest
        # tested concentration) must not send the shared Kd to an arbitrarily huge value.
        from kd_models import csp_model
        from kd_fit import fit_global_kd_csp
        residues = {n: csp_model(L, dd, 1e8, P0) for n, dd in (('A', 0.30), ('B', 0.12))}
        g = fit_global_kd_csp(residues, L, P0)
        assert g['success']
        hi = L.max() * 10.0
        assert g['Kd'] <= hi
        assert g['reliable'] is False

    def test_global_intensity_reliable_true_for_well_determined_fit(self):
        """AFFINITY_PLAYBOOK and AFFINITY_PROMPT tell an agent never to report a Kd
        without its `reliable` flag. Only the CSP global had one, so the instruction was
        unfollowable for obs=intensity."""
        from kd_models import intensity_decay
        from kd_fit import fit_global_kd_intensity
        residues = {n: intensity_decay(L, i0, inf, 22.0)
                    for n, i0, inf in (('A', 1.0, 0.05), ('B', 0.8, 0.20))}
        g = fit_global_kd_intensity(residues, L)
        assert g['success']
        assert g['reliable'] is True

    def test_global_intensity_flat_data_is_not_reliable(self):
        """Data that never decays does not push Kd out of the window — the model fits it
        by driving I_inf to I0, which cancels the exponential and leaves Kd free near its
        seed. So this case is caught by the relative error, not the window: the two gates
        are not redundant, and dropping either would pass a Kd nothing measured."""
        from kd_models import intensity_decay
        from kd_fit import fit_global_kd_intensity, _resolvable_window
        residues = {n: intensity_decay(L, i0, 0.0, 1e8)
                    for n, i0 in (('A', 1.0), ('B', 0.8))}
        g = fit_global_kd_intensity(residues, L)
        assert g['success']
        lo, hi = _resolvable_window(L)
        assert lo <= g['Kd'] <= hi                     # inside the window ...
        assert g['Kd_err'] > 0.3 * g['Kd']             # ... but unconstrained
        assert g['reliable'] is False

    def test_global_intensity_below_the_titration_range_is_not_reliable(self):
        """The mirror case: decayed at every measured point, so any smaller Kd fits
        equally well."""
        from kd_models import intensity_decay
        from kd_fit import fit_global_kd_intensity
        residues = {n: intensity_decay(L, i0, 0.0, 1e-6)
                    for n, i0 in (('A', 1.0), ('B', 0.8))}
        g = fit_global_kd_intensity(residues, L)
        assert g['success']
        assert g['reliable'] is False

    def test_global_csp_reliable_true_for_well_determined_fit(self):
        from kd_models import csp_model
        from kd_fit import fit_global_kd_csp
        residues = {n: csp_model(L, dd, 18.0, P0)
                    for n, dd in (('A', 0.30), ('B', 0.12), ('C', 0.45))}
        g = fit_global_kd_csp(residues, L, P0)
        assert g['success']
        assert g['reliable'] is True


class TestJsonSafe:
    def test_nan_inf_become_none(self):
        from kd_fit import json_safe
        out = json_safe({'a': float('nan'), 'b': [1.0, float('inf')], 'c': {'d': 2.0}})
        assert out['a'] is None
        assert out['b'] == [1.0, None]
        assert out['c']['d'] == 2.0
        import json
        json.loads(json.dumps(out))   # round-trips as strict JSON (no bare NaN)


class TestIntensityScaling:
    def test_per_point_scale_affects_intensity_not_csp(self, tmp_path):
        from kd_models import csp_model
        from kd_fit import run_kd_analysis_with_params
        import json
        pts = [0.0, 10.0, 25.0, 60.0, 150.0]
        ddH = csp_model(np.array(pts), 0.2, 20.0, P0)
        rows = []
        for p, d in zip(pts, ddH):
            rows.append((str(p), 'R1', 8.0 + d, 120.0, 1000.0 - 5 * p, 2000.0))
        df = pd.DataFrame(rows, columns=['spectrum_name', 'assignment',
                                         'ppm_x', 'ppm_y', 'height', 'volume'])
        csv = tmp_path / 'series_analysis_tidy.csv'
        df.to_csv(csv, index=False)

        def run(scales, out):
            r = run_kd_analysis_with_params({
                'input_csv_file': str(csv), 'output_dir': str(out), 'output_prefix': 'k',
                'concentrations': pts, 'protein_conc': P0, 'alpha': 0.14,
                'observables': ['csp', 'intensity'], 'intensity_scales': scales})
            return json.loads(Path(r['json_file']).read_text())['fits'][0]

        base = run(None, tmp_path / 'a')
        scaled = run([2.0, 1.0, 1.0, 1.0, 1.0], tmp_path / 'b')   # double the reference point
        # CSP identical (positions untouched)
        assert scaled['csp']['Kd'] == pytest.approx(base['csp']['Kd'], rel=1e-6)
        # intensity ratio changed: reference doubled -> ratios halved
        assert scaled['intensity']['obs'][1] == pytest.approx(base['intensity']['obs'][1] / 2.0, rel=1e-6)

    def test_scale_count_mismatch_raises(self, tmp_path):
        from kd_models import csp_model
        from kd_fit import run_kd_analysis_with_params
        df = pd.DataFrame([(str(p), 'R1', 8.0, 120.0, 1000.0, 2000.0) for p in (0.0, 1.0, 2.0)],
                          columns=['spectrum_name', 'assignment', 'ppm_x', 'ppm_y', 'height', 'volume'])
        csv = tmp_path / 'series_analysis_tidy.csv'
        df.to_csv(csv, index=False)
        with pytest.raises(ValueError):
            run_kd_analysis_with_params({
                'input_csv_file': str(csv), 'output_dir': str(tmp_path), 'output_prefix': 'k',
                'concentrations': [0.0, 1.0, 2.0], 'protein_conc': P0,
                'observables': ['intensity'], 'intensity_scales': [1.0, 2.0]})  # only 2 scales


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

    def test_global_fit_excludes_low_r2_residues(self, tmp_path):
        # A residue whose own decay is poorly fit (low R²) must NOT pool into the global
        # shared-Kd fit — otherwise one garbage residue hijacks the shared Kd.
        from kd_models import intensity_decay
        from kd_fit import run_kd_analysis_with_params
        import json
        pts = [0.0, 10.0, 25.0, 60.0, 150.0, 300.0]
        rows = []
        for name, (i0, iinf) in {'A1': (1000.0, 50.0), 'K2': (800.0, 200.0),
                                 'G3': (1200.0, 0.0)}.items():
            hs = intensity_decay(np.array(pts), i0, iinf, 40.0)   # clean, shared Kd=40
            rows += [(str(p), name, 8.0, 120.0, h, 2 * h) for p, h in zip(pts, hs)]
        # garbage residue: near-flat noisy intensity -> poor exp-decay fit (low R²)
        bad = [500.0, 505.0, 495.0, 500.0, 510.0, 490.0]
        rows += [(str(p), 'BAD', 8.0, 120.0, h, 2 * h) for p, h in zip(pts, bad)]
        df = pd.DataFrame(rows, columns=['spectrum_name', 'assignment',
                                         'ppm_x', 'ppm_y', 'height', 'volume'])
        csv = tmp_path / 'series_analysis_tidy.csv'
        df.to_csv(csv, index=False)
        result = run_kd_analysis_with_params({
            'input_csv_file': str(csv), 'output_dir': str(tmp_path), 'output_prefix': 'g',
            'concentrations': pts, 'protein_conc': P0,
            'observables': ['intensity'], 'n_bootstrap': 0})
        data = json.loads(Path(result['json_file']).read_text())
        gi = data['global']['intensity']
        assert gi['Kd'] == pytest.approx(40.0, rel=0.1)   # not dragged by BAD
        assert 'BAD' not in gi['I0']                       # low-R² residue excluded
        assert set(gi['I0']) == {'A1', 'K2', 'G3'}

    def test_end_to_end_writes_global_intensity(self, tmp_path):
        from kd_models import intensity_decay
        from kd_fit import run_kd_analysis_with_params
        import json
        pts = [0.0, 10.0, 25.0, 60.0, 150.0, 300.0]
        # two residues sharing a decay constant, differing amplitudes
        specs = {'R1': (1000.0, 50.0, 40.0), 'R2': (800.0, 200.0, 40.0)}
        rows = []
        for name, (i0, iinf, kd) in specs.items():
            heights = intensity_decay(np.array(pts), i0, iinf, kd)
            for p, h in zip(pts, heights):
                rows.append((str(p), name, 8.0, 120.0, h, 2 * h))
        df = pd.DataFrame(rows, columns=['spectrum_name', 'assignment',
                                         'ppm_x', 'ppm_y', 'height', 'volume'])
        csv = tmp_path / 'series_analysis_tidy.csv'
        df.to_csv(csv, index=False)
        result = run_kd_analysis_with_params({
            'input_csv_file': str(csv), 'output_dir': str(tmp_path), 'output_prefix': 'ki',
            'concentrations': pts, 'protein_conc': P0, 'alpha': 0.14,
            'observables': ['intensity'], 'n_bootstrap': 0})
        data = json.loads(Path(result['json_file']).read_text())
        gi = data['global']['intensity']
        assert gi['success']
        assert gi['Kd'] == pytest.approx(40.0, rel=5e-2)
        # results.txt records the apparent-Kd line
        txt = Path(result['results_file']).read_text()
        assert 'apparent Kd (intensity decay)' in txt
        # ... and the reliability verdict beside it. The runbooks say never to quote a
        # Kd without stating whether `reliable` was true, which a reader of this file
        # could not do while the flag existed only in the JSON.
        assert 'reliable' in txt.lower()

    def _tidy_at_precision(self, tmp_path, decimals):
        """A CSP titration whose fitted positions carry `decimals` decimal places."""
        from kd_models import csp_model
        pts = [0.0, 10.0, 25.0, 60.0, 150.0, 300.0]
        rows = []
        for name, dd in (('R1', 0.25), ('R2', 0.18)):
            dh = csp_model(np.array(pts), dd, 30.0, P0)
            for p, d in zip(pts, dh):
                rows.append((str(p), name, round(8.0 + d, decimals),
                             round(120.0 + 5 * d, decimals), 1000.0, 2000.0))
        df = pd.DataFrame(rows, columns=['spectrum_name', 'assignment',
                                         'ppm_x', 'ppm_y', 'height', 'volume'])
        csv = tmp_path / 'series_analysis_tidy.csv'
        df.to_csv(csv, index=False)
        return csv

    def _run(self, csv, tmp_path, observables):
        from kd_fit import run_kd_analysis_with_params
        import json
        result = run_kd_analysis_with_params({
            'input_csv_file': str(csv), 'output_dir': str(tmp_path / f"o_{observables[0]}"),
            'output_prefix': 'p', 'protein_conc': P0, 'alpha': 0.14,
            'observables': observables, 'n_bootstrap': 0})
        return json.loads(Path(result['json_file']).read_text()), result

    def test_coarse_positions_warn_because_csp_cannot_survive_them(self, tmp_path):
        """Series runs before positions were written to 4 decimals rounded them to
        0.1 ppm. CSP is the difference between two positions, so a 0.1 ppm grid
        quantises it at five times the significance threshold — and the fit still
        returns a finite Kd with success=True and exit 0."""
        data, _ = self._run(self._tidy_at_precision(tmp_path, 1), tmp_path, ['csp'])
        assert data['metadata']['position_precision_warning']
        assert any('precision' in w.lower() for w in data['metadata']['quality_warnings'])

    def test_four_decimal_positions_do_not_warn(self, tmp_path):
        data, _ = self._run(self._tidy_at_precision(tmp_path, 4), tmp_path, ['csp'])
        assert not data['metadata']['position_precision_warning']

    def test_an_intensity_only_fit_on_coarse_positions_is_fine(self, tmp_path):
        """Intensity never reads a position, so the same CSV is perfectly valid for it.
        Warning there would train the reader to ignore the warning."""
        data, _ = self._run(self._tidy_at_precision(tmp_path, 1), tmp_path, ['intensity'])
        assert not data['metadata']['position_precision_warning']
        assert not [w for w in data['metadata']['quality_warnings']
                    if 'precision' in w.lower()]

    def test_embeds_raw_series_per_residue(self, tmp_path):
        # The comparison view (arbitrary reference point) and exact-reopen need the
        # raw per-point positions/intensities, not just the point-0-relative obs.
        from kd_models import csp_model
        from kd_fit import run_kd_analysis_with_params
        pts = [0.0, 10.0, 25.0, 60.0, 150.0]
        ddH = csp_model(np.array(pts), 0.2, 20.0, P0)
        rows = [(str(p), 'R1', 8.0 + d, 120.0 + i, 1000.0 - 3 * i, 2000.0 - 5 * i)
                for i, (p, d) in enumerate(zip(pts, ddH))]
        df = pd.DataFrame(rows, columns=['spectrum_name', 'assignment',
                                         'ppm_x', 'ppm_y', 'height', 'volume'])
        csv = tmp_path / 'series_analysis_tidy.csv'
        df.to_csv(csv, index=False)
        result = run_kd_analysis_with_params({
            'input_csv_file': str(csv), 'output_dir': str(tmp_path),
            'output_prefix': 'kd_test', 'concentrations': pts, 'protein_conc': P0,
            'alpha': 0.14, 'observables': ['csp', 'intensity'], 'n_bootstrap': 0})
        import json
        data = json.loads(Path(result['json_file']).read_text())
        r1 = next(f for f in data['fits'] if f['residue'] == 'R1')
        s = r1['series']
        assert len(s['ppm_x']) == len(pts) == len(s['ppm_y'])
        assert len(s['height']) == len(pts) == len(s['volume'])
        assert s['ppm_y'][0] == pytest.approx(120.0)
        assert s['ppm_y'][-1] == pytest.approx(124.0)   # raw N positions preserved
        assert s['height'][-1] == pytest.approx(1000.0 - 3 * 4)
