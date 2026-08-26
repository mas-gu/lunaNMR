# ABOUTME: The relaxation fit must not carry a free baseline: it is unidentifiable at realistic sampling.
# ABOUTME: C is fixed at zero by default; a shared-baseline diagnostic reports when that assumption is wrong.

import sys
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_DIR))

from fit_Tx_NMRRE import (fit_single_residue, shared_baseline_test,                          save_results, DEGENERATE_T2_OVER_TMAX)

T_TRUE = 900.0
DELAYS = np.array([0, 100, 150, 300, 600, 900, 1200, 1500], float)   # a real truncated T1 grid: t_max/T = 1.7


def _decay(t=T_TRUE, amp=1e7, offset=0.0, noise=0.0, seed=0):
    y = amp * np.exp(-DELAYS / t) + offset
    if noise:
        y = y + np.random.default_rng(seed).normal(0, noise * amp, DELAYS.size)
    return y


class TestDefaultIsNoBaseline:
    def test_c_is_fixed_at_zero(self):
        res = fit_single_residue(DELAYS, _decay(), 'R1', n_bootstrap=0)
        assert res['C'] == 0.0
        assert res['baseline_fixed'] is True

    def test_c_stays_zero_even_when_the_data_has_an_offset(self):
        """The default must not chase an offset. Whether one is real is decided by the
        shared-baseline test over all residues, not per residue where it is unidentifiable."""
        res = fit_single_residue(DELAYS, _decay(offset=0.15e7), 'R1', n_bootstrap=0)
        assert res['C'] == 0.0

    def test_the_offset_model_is_still_reachable_and_still_works(self):
        """Capability preserved for the shared-baseline diagnostic; only the default changed."""
        res = fit_single_residue(DELAYS, _decay(offset=0.15e7), 'R1',
                                 n_bootstrap=0, fit_baseline=True)
        assert res['baseline_fixed'] is False
        assert res['C'] == pytest.approx(0.15e7, rel=0.05)

    def test_t_is_recovered_when_the_baseline_really_is_zero(self):
        res = fit_single_residue(DELAYS, _decay(noise=0.03, seed=1), 'R1', n_bootstrap=0)
        assert res['t2'] == pytest.approx(T_TRUE, rel=0.10)


class TestPrecision:
    """The reason for the change: at t_max/T ~ 1.7 the three-parameter fit is unidentifiable."""

    def test_fixing_the_baseline_sharpens_t_on_a_truncated_series(self):
        errs_free, errs_fixed = [], []
        for seed in range(25):
            y = _decay(noise=0.03, seed=seed)
            errs_free.append(fit_single_residue(DELAYS, y, 'R', n_bootstrap=0,
                                                fit_baseline=True)['t2_err'])
            errs_fixed.append(fit_single_residue(DELAYS, y, 'R', n_bootstrap=0)['t2_err'])
        free, fixed = np.nanmedian(errs_free), np.nanmedian(errs_fixed)
        assert fixed < free / 2.0, f'expected >2x sharper, got free={free:.1f} fixed={fixed:.1f}'

    def test_and_does_not_bias_t_when_the_baseline_is_zero(self):
        """Precision bought by biasing the centre would be worse than the status quo."""
        got = [fit_single_residue(DELAYS, _decay(noise=0.03, seed=s), 'R',
                                  n_bootstrap=0)['t2'] for s in range(25)]
        assert np.median(got) == pytest.approx(T_TRUE, rel=0.05)


class TestBootstrapRespectsIt:
    def test_bootstrap_does_not_silently_re_free_the_baseline(self):
        """bootstrap_errors rebuilds params via make_params; without care it re-frees C."""
        res = fit_single_residue(DELAYS, _decay(noise=0.03, seed=2), 'R1',
                                 n_bootstrap=40, error_method='bootstrap')
        assert res['C'] == 0.0
        assert res['C_err'] == 0.0


class TestOutputContract:
    def test_a_fixed_baseline_is_distinguishable_from_an_unmeasured_one(self):
        """C_err = NaN means the covariance was singular. A fixed C must not look like that."""
        res = fit_single_residue(DELAYS, _decay(), 'R1', n_bootstrap=0)
        assert res['C_err'] == 0.0 and not np.isnan(res['C_err'])
        assert res['baseline_fixed'] is True

    def test_results_table_carries_the_success_column(self, tmp_path):
        """fitting_wrapper reads row['Success']; the single-core writer omitted it."""
        out = tmp_path / 'r.txt'
        save_results([{'residue': 'R1', 'A': 1e7, 't2': 900.0, 'C': 0.0,
                           'A_err': 1e5, 't2_err': 40.0, 'C_err': 0.0,
                           'reliable': True}], str(out), 'T1')
        lines = out.read_text().splitlines()
        header = lines[0].split('\t')
        assert 'Success' in header
        assert lines[1].split('\t')[header.index('Success')] == 'Yes'


class TestSharedBaselineDiagnostic:
    """A blanket C=0 is only honest if the case where it is wrong announces itself."""

    def _series(self, f, n=90, noise=0.03):
        rng = np.random.default_rng(7)
        out = {}
        for i in range(n):
            amp, t = 1e7 * (0.5 + rng.random()), T_TRUE * (0.8 + 0.4 * rng.random())
            y = amp * ((1 - f) * np.exp(-DELAYS / t) + f)
            out[f'R{i}'] = (DELAYS, y + rng.normal(0, noise * amp, DELAYS.size))
        return out

    def test_silent_when_there_is_no_shared_baseline(self):
        r = shared_baseline_test(self._series(0.0))
        assert r['significant'] is False
        assert abs(r['f']) < 0.03

    def test_detects_a_real_shared_baseline(self):
        """A real truncated T1 series sits at f = 0.15, with 70% of residues
        improving (p = 5e-07)."""
        r = shared_baseline_test(self._series(0.15))
        assert r['significant'] is True
        assert r['f'] == pytest.approx(0.15, abs=0.04)
        assert r['frac_improved'] > 0.6

    def test_reports_the_cost_of_assuming_zero(self):
        """Aggregated as a median of per-residue scaled SSR, so a handful of rows that do
        not fit at all cannot set the answer."""
        r = shared_baseline_test(self._series(0.15))
        assert r['profile_at_zero'] > r['profile_min']
        assert 0.0 <= r['p_value'] <= 1.0


class TestMeasurabilityGuard:
    def test_flags_a_window_too_short_to_measure_t(self):
        """t_max/T < 2 is where the fit stops being able to separate decay from offset."""
        res = fit_single_residue(DELAYS, _decay(t=1400.0), 'R1', n_bootstrap=0)
        assert res['window_ratio'] == pytest.approx(1500.0 / res['t2'], rel=0.01)
        assert res['window_marginal'] is True

    def test_quiet_on_a_well_sampled_series(self):
        res = fit_single_residue(DELAYS, _decay(t=400.0), 'R1', n_bootstrap=0)
        assert res['window_marginal'] is False

    def test_the_degeneracy_constant_is_defined_once(self):
        import fitmulti__Tx_NMRRE as fm
        assert fm.DEGENERATE_T2_OVER_TMAX is DEGENERATE_T2_OVER_TMAX


class TestProvenanceSurvivesToDisk:
    """The fitter computed these and both writers dropped them, so nothing outside the
    process could see them. Assert the round trip, not just the in-memory dict."""

    def _run(self, tmp_path, tau, delays):
        import pandas as pd
        from fit_Tx_NMRRE import run_analysis_with_params
        rng = np.random.default_rng(3)
        cols = [str(d) for d in delays]
        rows = [[i + 1, f"R{i+1}", 8.0, 120.0]
                + [1e6 * np.exp(-d / tau) * (1 + rng.normal(0, 0.02)) for d in delays]
                for i in range(12)]
        csv = tmp_path / 'm.csv'
        pd.DataFrame(rows, columns=["Peak_Number", "Assignment", "Reference_X",
                                    "Reference_Y"] + cols).to_csv(csv, index=False)
        return run_analysis_with_params({
            'input_csv_file': str(csv), 'output_prefix': str(tmp_path / 'f'),
            'results_txt_file': str(tmp_path / 'f.txt'), 'experiment_type': 'T1',
            'json_folder': str(tmp_path)})

    def test_json_records_carry_baseline_and_window(self, tmp_path):
        import glob, json
        self._run(tmp_path, 900.0, [0, 100, 150, 300, 600, 900, 1200, 1500])
        data = json.load(open(glob.glob(str(tmp_path / '*_fit_data.json'))[0]))
        rec = data['fits'][0]
        for key in ('baseline_fixed', 'window_ratio', 'window_marginal'):
            assert key in rec, f'{key} never reaches the JSON'
        assert rec['baseline_fixed'] is True

    def test_results_table_carries_the_window_ratio(self, tmp_path):
        self._run(tmp_path, 900.0, [0, 100, 150, 300, 600, 900, 1200, 1500])
        header = (tmp_path / 'f.txt').read_text().splitlines()[0].split('\t')
        assert 'WindowRatio' in header

    def test_a_marginal_window_is_reported_in_the_run_summary(self, tmp_path):
        """t_max 1500 against T ~ 1400 is the marginal-window regime real data shows."""
        r = self._run(tmp_path, 1400.0, [0, 100, 150, 300, 600, 900, 1200, 1500])
        assert r['n_window_marginal'] >= 1
        assert r['baseline_fixed'] is True


class TestSharedBaselineRobustness:
    """Three defects found by a peer running it on real data: it ignored the package's
    dummy_* convention, silently collapsed duplicate residue keys, and reported high
    confidence (F = 9.7, p = 0.002) from a profile whose minimum was not determined."""

    def _flat_series(self, n=40, noise=0.03):
        """No baseline, so the SSR profile over f has no real minimum."""
        rng = np.random.default_rng(11)
        out = {}
        for i in range(n):
            amp, t = 1e7 * (0.5 + rng.random()), T_TRUE * (0.8 + 0.4 * rng.random())
            out[f'R{i}'] = (DELAYS, amp * np.exp(-DELAYS / t)
                            + rng.normal(0, noise * amp, DELAYS.size))
        return out

    def test_dummy_rows_are_excluded_and_counted(self):
        s = self._flat_series(20)
        rng = np.random.default_rng(2)
        for k in ('dummy_001', 'dummy_002'):
            s[k] = (DELAYS, rng.normal(0, 1e5, DELAYS.size))
        r = shared_baseline_test(s)
        assert r['n_excluded_dummy'] == 2
        assert r['n_residues'] == 20

    def test_duplicate_residue_names_are_not_silently_collapsed(self):
        """A dict cannot hold them, so a sequence of triples must be accepted; the raw
        file had dummy_001 twice, which is why two routes disagreed on n."""
        s = self._flat_series(10)
        items = [(k, x, y) for k, (x, y) in s.items()]
        items.append(items[0])                      # same residue twice
        r = shared_baseline_test(items)
        assert r['n_residues'] == 11

    def test_a_flat_profile_is_not_reported_as_significant(self):
        """The failure mode: an argmin that moved when 3 of 153 residues were removed,
        reported at p = 0.002. Robust aggregation plus a per-residue vote must veto it."""
        r = shared_baseline_test(self._flat_series())
        assert r['significant'] is False

    def test_a_real_shared_baseline_is_still_well_determined(self):
        rng = np.random.default_rng(7)
        s = {}
        for i in range(40):
            amp, t = 1e7 * (0.5 + rng.random()), T_TRUE * (0.8 + 0.4 * rng.random())
            s[f'R{i}'] = (DELAYS, amp * (0.85 * np.exp(-DELAYS / t) + 0.15)
                          + rng.normal(0, 0.02 * amp, DELAYS.size))
        r = shared_baseline_test(s)
        assert r['significant'] is True and r['well_determined'] is True
        assert r['f'] == pytest.approx(0.15, abs=0.04)

    def test_the_per_residue_vote_is_reported_so_a_caller_can_judge(self):
        """A pooled F-test turns a 1% SSR gain into p ~ 0 across hundreds of residues.
        The fraction of residues that individually improve is the honest statement."""
        r = shared_baseline_test(self._flat_series())
        assert 0 <= r['n_improved'] <= r['n_residues']
        assert r['frac_improved'] == pytest.approx(r['n_improved'] / r['n_residues'], abs=1e-9)


class TestGridQuantisation:
    """The returned f was pinned to the search grid, giving a +/-0.0125 floor at default
    settings — 8% of the value at f ~ 0.15, and enough to make two independent series
    report a spuriously identical estimate."""

    def _series(self, f, n=90, noise=0.03):
        rng = np.random.default_rng(21)
        out = {}
        for i in range(n):
            amp, t = 1e7 * (0.5 + rng.random()), T_TRUE * (0.8 + 0.4 * rng.random())
            out[f'R{i}'] = (DELAYS, amp * ((1 - f) * np.exp(-DELAYS / t) + f)
                            + rng.normal(0, noise * amp, DELAYS.size))
        return out

    def test_the_quantisation_floor_is_small_enough_to_quote_f_to_two_digits(self):
        """f IS grid-quantised; the contract is that the step is small and declared.
        At the old default of 25 points the floor was 0.0125, 8% of a typical f."""
        r = shared_baseline_test(self._series(0.17))
        assert r['f_grid_step'] == pytest.approx(0.01, abs=1e-9)

    def test_the_quantisation_floor_is_reported(self):
        r = shared_baseline_test(self._series(0.15))
        assert r['f_grid_step'] > 0

    def test_a_coarse_and_a_fine_grid_agree_within_the_step(self):
        """The refinement must make the answer insensitive to how the grid was chosen."""
        coarse = shared_baseline_test(self._series(0.17), n_grid=25)
        fine = shared_baseline_test(self._series(0.17), n_grid=121)
        assert abs(coarse['f'] - fine['f']) < coarse['f_grid_step']


class TestVerdictIsSelfExplaining:
    """A caller should not have to read fields in the right order to avoid a wrong
    conclusion. On real 800 T1rho every field looks green — well_determined True,
    frac_improved 60%, p 0.0088 — and only the SIGN of f stops the row."""

    def _series(self, f, n=90, noise=0.03, seed=31):
        rng = np.random.default_rng(seed)
        out = {}
        for i in range(n):
            amp, t = 1e7 * (0.5 + rng.random()), T_TRUE * (0.8 + 0.4 * rng.random())
            out[f'R{i}'] = (DELAYS, amp * ((1 - f) * np.exp(-DELAYS / t) + f)
                            + rng.normal(0, noise * amp, DELAYS.size))
        return out

    def test_a_negative_baseline_says_why_it_was_rejected(self):
        r = shared_baseline_test(self._series(-0.12))
        assert r['significant'] is False
        assert r['f'] < 0
        assert 'negative' in r['reason'].lower(), r['reason']

    def test_well_determined_describes_precision_not_existence(self):
        """The two genuinely disagree on that row, so the docs must not conflate them."""
        r = shared_baseline_test(self._series(-0.12))
        assert r['well_determined'] is True and r['significant'] is False

    def test_a_real_baseline_gives_an_empty_reason(self):
        r = shared_baseline_test(self._series(0.17))
        assert r['significant'] is True and r['reason'] == ''

    def test_an_undetectable_baseline_is_rejected_with_a_measured_reason(self):
        """A baseline far below the noise is caught by the per-residue vote rather than by
        the size gate: at 3% noise a true f = 0.004 leaves only 44% of residues improving.
        Either gate is a correct rejection; what matters is that it names a number."""
        r = shared_baseline_test(self._series(0.004))
        assert r['significant'] is False
        assert r['reason'] and any(ch.isdigit() for ch in r['reason']), r['reason']

    def test_the_size_gate_is_reachable(self):
        """Real 800 T2 lands at f = 0.0100 against a 0.01 step and is stopped here."""
        r = shared_baseline_test(self._series(0.0), n_grid=61)
        if 0 <= r['f'] <= r['f_grid_step']:
            assert 'grid step' in r['reason'], r['reason']
        assert r['significant'] is False

    def test_grid_step_is_not_a_floating_point_smear(self):
        r = shared_baseline_test(self._series(0.15))
        assert repr(r['f_grid_step']) == '0.01'


class TestDepthCarriesItsConfound:
    """profile_depth is not scale-free in t_max/T — truncating a real series moved it
    6.25% -> 15.71% and another 25.81% -> 15.08%. It orders within a series, not across,
    so window_ratio is reported alongside it."""

    def _series(self, f, delays, n=60, noise=0.03, seed=41):
        rng = np.random.default_rng(seed)
        return {f'R{i}': (delays, (a := 1e7 * (0.5 + rng.random()))
                          * ((1 - f) * np.exp(-delays / (T_TRUE * (0.8 + 0.4 * rng.random()))) + f)
                          + rng.normal(0, noise * a, delays.size)) for i in range(n)}

    def test_window_ratio_is_reported(self):
        r = shared_baseline_test(self._series(0.0, DELAYS))
        assert np.isfinite(r['window_ratio_at_f']) and r['window_ratio_at_f'] > 0

    def test_it_tracks_the_window_actually_used(self):
        """A series truncated to half its delays must report a smaller t_max/T."""
        full = shared_baseline_test(self._series(0.0, DELAYS))
        short = shared_baseline_test(self._series(0.0, DELAYS[:5]))
        assert short['window_ratio_at_f'] < full['window_ratio_at_f']

    def test_depth_is_not_stable_under_truncation(self):
        """Documenting the limitation rather than pretending it is absent: a caller
        comparing two depths must check window_ratio first."""
        full = shared_baseline_test(self._series(0.15, DELAYS))
        short = shared_baseline_test(self._series(0.15, DELAYS[:5]))
        assert abs(full['profile_depth'] - short['profile_depth']) > 0.01


class TestWindowRatioSaysWhichModelItIsUnder:
    """Two fields with the same meaning under different models diverged 37-43% on exactly
    the rows where a baseline is claimed, and straddled the t_max/T < 2 line — opposite
    readings of whether the window was adequate. The names now carry the conditioning."""

    def _series(self, f, n=60, noise=0.03, seed=53):
        rng = np.random.default_rng(seed)
        return {f'R{i}': (DELAYS, (a := 1e7 * (0.5 + rng.random()))
                          * ((1 - f) * np.exp(-DELAYS / (T_TRUE * (0.8 + 0.4 * rng.random()))) + f)
                          + rng.normal(0, noise * a, DELAYS.size)) for i in range(n)}

    def test_both_conditionings_are_reported(self):
        r = shared_baseline_test(self._series(0.15))
        assert np.isfinite(r['window_ratio_at_f'])
        assert np.isfinite(r['window_ratio_at_zero'])

    def test_they_diverge_when_a_baseline_is_claimed(self):
        """A plateau absorbs part of the tail, so T shortens and the ratio rises."""
        r = shared_baseline_test(self._series(0.15))
        assert r['window_ratio_at_f'] > r['window_ratio_at_zero'] * 1.1

    def test_divergence_tracks_the_size_of_f_not_the_verdict(self):
        """Keyed on |f| deliberately. The earlier version asserted agreement whenever the
        series was REJECTED, which held only because no fixture was both rejected and far
        from zero — a real T2 series is rejected on sign at f = -0.06 and diverges 11%.
        The relationship is monotone in |f| and signed by f: 0% at f = 0, -11% at -0.06,
        +43% at +0.15."""
        near = shared_baseline_test(self._series(0.0))
        assert near['window_ratio_at_f'] == pytest.approx(near['window_ratio_at_zero'], rel=0.05)

        # rejected, but a long way from zero: it must still diverge, and downward
        far = shared_baseline_test(self._series(-0.12))
        assert far['significant'] is False
        rel = (far['window_ratio_at_f'] - far['window_ratio_at_zero']) / far['window_ratio_at_zero']
        assert rel < -0.05, f'expected a negative f to shorten the ratio, got {rel:+.1%}' 
