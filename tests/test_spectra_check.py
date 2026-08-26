# ABOUTME: Pre-flight checks on a relaxation dataset - registration, capture, delays, peak lists.
# ABOUTME: Anything about labels or identity must use the pipeline's own parsers, so it cannot drift.

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

ng = pytest.importorskip("nmrglue")

from lunaNMR.utils.delay_extractor import DelayExtractor
from lunaNMR.validation.spectra_check import (
    Finding, assess_registration, check_dataset, check_experiment,
)

# 15N: 19.7 ppm over 512 pts = 0.0385 ppm/pt.  1H: 10 ppm over 1024 = 0.0098 ppm/pt.
_AXES = {0: dict(sw=1200.0, obs=60.8, car=120.0, size=512, label='15N', encoding='states'),
         1: dict(sw=6000.0, obs=600.1, car=8.0, size=1024, label='1H', encoding='direct')}


def _spectrum(path, peaks, *, shift=(0.0, 0.0), amps=None, noise=0.0, seed=0):
    """Write a real NMRPipe file with a Gaussian at each peak, optionally displaced."""
    udic = ng.fileiobase.create_blank_udic(2)
    for d, ax in _AXES.items():
        udic[d].update(dict(sw=ax['sw'], obs=ax['obs'], car=ax['car'] * ax['obs'],
                            size=ax['size'], label=ax['label'], encoding=ax['encoding'],
                            complex=False, time=False, freq=True))
    dic = ng.pipe.create_dic(udic)
    data = np.zeros((_AXES[0]['size'], _AXES[1]['size']), 'float32')
    if noise:
        data += np.random.default_rng(seed).normal(0, noise, data.shape).astype('float32')
    uy = ng.pipe.make_uc(dic, data, dim=0)
    ux = ng.pipe.make_uc(dic, data, dim=1)
    amps = amps if amps is not None else [1e6] * len(peaks)
    yy, xx = np.mgrid[0:data.shape[0], 0:data.shape[1]]
    for (_, x, y), amp in zip(peaks, amps):
        iy, ix = uy(y + shift[1], 'ppm'), ux(x + shift[0], 'ppm')
        data += (amp * np.exp(-((yy - iy) ** 2) / 8.0 - ((xx - ix) ** 2) / 8.0)).astype('float32')
    ng.pipe.write(str(path), dic, data, overwrite=True)
    return path


def _peaks(n=6):
    return [(f'{10 + 3 * i}ValH', 8.0 + 0.35 * i, 118.0 + 1.7 * i) for i in range(n)]


def _write_list(path, peaks, *, lead_ws=False, dup=False, dummy=False):
    lines = ['Assignment, Position_X, Position_Y']
    for name, x, y in peaks:
        lines.append(f"{' ' if lead_ws else ''}{name}, {x:.4f}, {y:.4f}")
    if dup:
        n, x, y = peaks[0]
        lines.append(f"{n}, {x:.4f}, {y:.4f}")
    if dummy:
        lines.append("dummy_1, 9.5, 130.0")
    path.write_text('\n'.join(lines) + '\n')
    return path


def _folder(tmp_path, name, spectra, peaks, **kw):
    """One experiment folder: spectra named as given, plus its own peak list."""
    d = tmp_path / name
    d.mkdir()
    for fname, opts in spectra.items():
        _spectrum(d / fname, peaks, **opts)
    _write_list(d / f'{name}_assignment.txt', peaks, **kw)
    return d


def _sev(result, check):
    return [f.severity for f in result['findings'] if f.check == check]


class TestDelayReporting:
    """The delays reported must be the delays `series` will actually parse."""

    @pytest.mark.parametrize('stem', ['T1_2400ms', 'T1_2s', 'T1_2400', 'T1_300msb', 'vd_2400us',
                                      'T2_8ms_2', 'T1_500us'])
    def test_matches_the_pipeline_parser(self, tmp_path, stem):
        peaks = _peaks()
        d = _folder(tmp_path, 'T1', {f'{stem}.ft': {}}, peaks)
        got = check_experiment(d, d / 'T1_assignment.txt', quick=True)
        expected = DelayExtractor(mode='time').extract_value(f'{stem}.ft')
        assert got['delays']['parsed'].get(f'{stem}.ft') == expected

    def test_unparsed_spectra_are_named_not_silently_dropped(self, tmp_path):
        peaks = _peaks()
        d = _folder(tmp_path, 'T1', {'T1_50ms.ft': {}, 'T1_2400.ft': {}, 'T1_100ms.ft': {}}, peaks)
        r = check_experiment(d, d / 'T1_assignment.txt', quick=True)
        assert r['delays']['n_unparsed'] == 1
        assert 'T1_2400.ft' in r['delays']['unparsed']
        assert 'WARN' in _sev(r, 'delay_parsing')


class TestPeakListHygiene:
    """Hygiene is judged on what the loader returns, not on the raw bytes."""

    def test_leading_whitespace_is_not_a_failure(self, tmp_path):
        """The loader strips it, so flagging it would be a false positive."""
        peaks = _peaks()
        d = _folder(tmp_path, 'T1', {'T1_50ms.ft': {}}, peaks, lead_ws=True)
        r = check_experiment(d, d / 'T1_assignment.txt', quick=True)
        assert 'FAIL' not in _sev(r, 'peak_list')
        assert all(a == a.strip() for a in r['peak_list']['assignments'])

    def test_duplicate_assignments_fail(self, tmp_path):
        peaks = _peaks()
        d = _folder(tmp_path, 'T1', {'T1_50ms.ft': {}}, peaks, dup=True)
        r = check_experiment(d, d / 'T1_assignment.txt', quick=True)
        assert 'FAIL' in _sev(r, 'peak_list')

    def test_dummy_entries_warn(self, tmp_path):
        peaks = _peaks()
        d = _folder(tmp_path, 'T1', {'T1_50ms.ft': {}}, peaks, dummy=True)
        r = check_experiment(d, d / 'T1_assignment.txt', quick=True)
        assert 'WARN' in _sev(r, 'peak_list')


class TestSubSeriesScale:
    """Repeat acquisitions are found by delay collision, so naming cannot hide them."""

    def test_infix_b_repeats_are_compared(self, tmp_path):
        """The `_b_` convention: the old prefix split collapsed these into one group."""
        peaks = _peaks()
        d = _folder(tmp_path, 'T2', {
            'T2_A1_51ms.ft': {}, 'T2_A1_102ms.ft': {}, 'T2_A1_170ms.ft': {},
            'T2_A1_b_51ms.ft': {'amps': [5e5] * 6}, 'T2_A1_b_102ms.ft': {'amps': [5e5] * 6},
        }, peaks)
        r = check_experiment(d, d / 'T2_assignment.txt', quick=True)
        assert r['subseries'] is not None
        assert sorted(r['subseries']['shared']) == [51.0, 102.0]

    def test_prefix_repeats_are_compared(self, tmp_path):
        """The `bT2_` convention, which the old prefix split did handle."""
        peaks = _peaks()
        d = _folder(tmp_path, 'T2', {
            'T2_8ms.ft': {}, 'T2_102ms.ft': {},
            'bT2_8ms.ft': {'amps': [5e5] * 6}, 'bT2_102ms.ft': {'amps': [5e5] * 6},
        }, peaks)
        r = check_experiment(d, d / 'T2_assignment.txt', quick=True)
        assert sorted(r['subseries']['shared']) == [8.0, 102.0]

    def test_matched_repeats_raise_nothing(self, tmp_path):
        peaks = _peaks()
        d = _folder(tmp_path, 'T2', {
            'T2_51ms.ft': {}, 'T2_102ms.ft': {},
            'T2_b_51ms.ft': {}, 'T2_b_102ms.ft': {},
        }, peaks)
        r = check_experiment(d, d / 'T2_assignment.txt', quick=True)
        assert r['subseries']['ratio'] == pytest.approx(1.0, abs=0.02)
        assert _sev(r, 'subseries') == []

    def test_mismatched_repeats_fail_and_report_the_scale(self, tmp_path):
        peaks = _peaks()
        d = _folder(tmp_path, 'T2', {
            'T2_51ms.ft': {}, 'T2_102ms.ft': {},
            'T2_b_51ms.ft': {'amps': [8e5] * 6}, 'T2_b_102ms.ft': {'amps': [8e5] * 6},
        }, peaks)
        r = check_experiment(d, d / 'T2_assignment.txt', quick=True)
        assert r['subseries']['ratio'] == pytest.approx(0.8, abs=0.02)
        assert r['subseries']['scale'] == pytest.approx(1.25, abs=0.04)
        assert 'FAIL' in _sev(r, 'subseries')


class TestRegistration:
    """The one check the pipeline has no equivalent for."""

    def test_recovers_a_known_offset(self, tmp_path):
        peaks = _peaks()
        f = _spectrum(tmp_path / 's.ft', peaks, shift=(0.05, 0.40))
        dx, dy = assess_registration(f, peaks)
        assert dx == pytest.approx(0.05, abs=0.006)
        assert dy == pytest.approx(0.40, abs=0.05)

    def test_aligned_list_reports_no_shift(self, tmp_path):
        peaks = _peaks()
        f = _spectrum(tmp_path / 's.ft', peaks)
        dx, dy = assess_registration(f, peaks)
        assert abs(dx) < 0.006 and abs(dy) < 0.05

    def test_large_offset_fails(self, tmp_path):
        peaks = _peaks()
        d = tmp_path / 'T1'
        d.mkdir()
        _spectrum(d / 'T1_50ms.ft', peaks, shift=(0.06, 0.0))
        _write_list(d / 'T1_assignment.txt', peaks)
        r = check_experiment(d, d / 'T1_assignment.txt')
        assert 'FAIL' in _sev(r, 'registration')


class TestCaptureRate:
    """Only the reference point can distinguish missing peaks from a decayed tail."""

    def test_missing_peaks_at_the_reference_point_warn(self, tmp_path):
        peaks = _peaks(6)
        d = tmp_path / 'T1'
        d.mkdir()
        # half the list is absent from the very first point, before anything can decay
        _spectrum(d / 'T1_50ms.ft', peaks, amps=[1e6, 1e6, 1e6, 0, 0, 0], noise=1e3)
        _spectrum(d / 'T1_100ms.ft', peaks, amps=[1e6] * 6, noise=1e3)
        _write_list(d / 'T1_assignment.txt', peaks)
        r = check_experiment(d, d / 'T1_assignment.txt', quick=True)
        assert 'WARN' in _sev(r, 'capture')

    def test_a_decayed_tail_raises_nothing(self, tmp_path):
        """Every relaxation series ends in the noise; warning on it buries the real cases."""
        peaks = _peaks(6)
        d = tmp_path / 'T1'
        d.mkdir()
        _spectrum(d / 'T1_50ms.ft', peaks, amps=[1e6] * 6, noise=1e3)
        _spectrum(d / 'T1_2400ms.ft', peaks, amps=[4e3] * 6, noise=1e3)
        _write_list(d / 'T1_assignment.txt', peaks)
        r = check_experiment(d, d / 'T1_assignment.txt', quick=True)
        assert _sev(r, 'capture') == []
        tail = [s for s in r['spectra'] if s['spectrum'] == 'T1_2400ms.ft'][0]
        assert tail['decayed'] is True and tail['is_reference'] is False

    def test_a_tie_at_the_shortest_delay_resolves_the_way_series_labels_columns(self, tmp_path):
        """Both files sit at 8 ms; the reference must be the one owning the bare label."""
        peaks = _peaks(6)
        d = tmp_path / 'T2'
        d.mkdir()
        _spectrum(d / 'T2_8ms.ft', peaks, amps=[1e6] * 6, noise=1e3)
        _spectrum(d / 'bT2_8ms.ft', peaks, amps=[1e6] * 6, noise=1e3)
        _spectrum(d / 'T2_102ms.ft', peaks, amps=[5e5] * 6, noise=1e3)
        _write_list(d / 'T2_assignment.txt', peaks)
        r = check_experiment(d, d / 'T2_assignment.txt', quick=True)
        assert r['reference_spectrum'] == 'T2_8ms.ft'

    def test_the_reference_point_is_the_shortest_delay_not_the_first_filename(self, tmp_path):
        peaks = _peaks(6)
        d = tmp_path / 'T1'
        d.mkdir()
        _spectrum(d / 'T1_100ms.ft', peaks, amps=[1e6] * 6, noise=1e3)
        _spectrum(d / 'T1_50ms.ft', peaks, amps=[1e6] * 6, noise=1e3)
        _write_list(d / 'T1_assignment.txt', peaks)
        r = check_experiment(d, d / 'T1_assignment.txt', quick=True)
        assert r['reference_spectrum'] == 'T1_50ms.ft'


class TestHetnoePlanes:
    def test_saturated_plane_identified_by_intensity_ratio(self, tmp_path):
        peaks = _peaks()
        d = _folder(tmp_path, 'hetNOE', {
            'plane_001.ft': {},
            'plane_002.ft': {'amps': [8e5] * 6},
        }, peaks)
        r = check_experiment(d, d / 'hetNOE_assignment.txt', quick=True)
        assert r['hetnoe']['saturated'] == 'plane_002.ft'
        assert r['hetnoe']['unsaturated'] == 'plane_001.ft'
        assert r['hetnoe']['ratio'] == pytest.approx(0.8, abs=0.02)


class TestDataset:
    def test_cross_field_assignment_intersection(self, tmp_path):
        a, b = _peaks(6), _peaks(6)[:5]
        _folder(tmp_path, 'T1', {'T1_50ms.ft': {}}, a)
        _folder(tmp_path, 'T2', {'T2_8ms.ft': {}}, b)
        r = check_dataset(tmp_path, quick=True)
        assert r['assignments']['common'] == 5
        assert r['assignments']['union'] == 6
        assert 'WARN' in [f.severity for f in r['findings'] if f.check == 'assignments']

    def test_reports_each_experiment(self, tmp_path):
        peaks = _peaks()
        _folder(tmp_path, 'T1', {'T1_50ms.ft': {}}, peaks)
        _folder(tmp_path, 'T2', {'T2_8ms.ft': {}}, peaks)
        r = check_dataset(tmp_path, quick=True)
        assert sorted(e['experiment'] for e in r['experiments']) == ['T1', 'T2']
