# ABOUTME: T1rho -> R2 conversion and its CLI surface, so a T1/T1rho/hetNOE dataset
# ABOUTME: reaches model-free without anyone re-deriving the tilt-angle algebra by hand.

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from lunaNMR.utils.t1rho_calculator import (
    calculate_R2_from_R1_R1rho, calculate_per_residue_theta, r2_table_from_fits)

TRUE_R1, TRUE_R2 = 1.5, 12.0
DELAYS_T1 = [0, 100, 300, 600, 900, 1200]
DELAYS_RHO = [4, 8, 16, 24, 40, 64, 84, 104]


def _peaks(n=6):
    return [(f'{10 + 3 * i}ValH', 8.0 + 0.3 * i, 118.0 + 2.0 * i) for i in range(n)]


def _matrix(path, peaks, delays, rate, amplitude=1e7):
    """A clean mono-exponential decay at a known rate, in peak_intensity_matrix shape."""
    rows = []
    for i, (name, x, y) in enumerate(peaks, 1):
        row = {'Peak_Number': i, 'Assignment': name, 'Reference_X': x, 'Reference_Y': y}
        for d in delays:
            row[str(d)] = amplitude * np.exp(-rate * d / 1000.0)
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def _peak_file(path, peaks):
    lines = ['Assignment, Position_X, Position_Y']
    lines += [f'{n}, {x:.4f}, {y:.4f}' for n, x, y in peaks]
    path.write_text('\n'.join(lines) + '\n')
    return path


class TestTiltAngleAlgebra:
    def test_r2_inverts_the_r1rho_relation(self):
        """R1rho = R1 cos^2(theta) + R2 sin^2(theta); the conversion must invert it."""
        theta = 55.0
        r1rho = (TRUE_R1 * np.cos(np.radians(theta)) ** 2
                 + TRUE_R2 * np.sin(np.radians(theta)) ** 2)
        assert calculate_R2_from_R1_R1rho(TRUE_R1, r1rho, theta) == pytest.approx(TRUE_R2)

    def test_on_resonance_r2_equals_r1rho(self):
        """At theta = 90 the spin-lock is fully transverse, so no R1 correction applies."""
        assert calculate_R2_from_R1_R1rho(TRUE_R1, 9.9, 90.0) == pytest.approx(9.9)

    def test_a_vanishing_tilt_is_refused_rather_than_dividing_by_zero(self):
        with pytest.raises(ValueError):
            calculate_R2_from_R1_R1rho(TRUE_R1, 9.9, 0.0)

    def test_theta_varies_with_chemical_shift(self):
        """The whole point of the per-residue correction: residues off the spin-lock
        carrier are tilted differently, so one nominal angle is not enough."""
        kw = dict(carrier_ppm=118.0, omega1_hz=1500.0, theta_nominal_deg=45.0,
                  spec_freq_mhz=600.0)
        near = calculate_per_residue_theta(residue_ppm=118.0, **kw)
        far = calculate_per_residue_theta(residue_ppm=132.0, **kw)
        assert abs(near - far) > 1.0


class TestR2Table:
    def test_recovers_r2_from_fitted_t1_and_t1rho(self):
        peaks = _peaks()
        theta = 90.0   # on-resonance: R2 == R1rho, so the expected answer is exact
        t1 = {n: {'value': 1000.0 / TRUE_R1, 'error': 1.0} for n, _, _ in peaks}
        rho = {n: {'value': 1000.0 / TRUE_R2, 'error': 0.5} for n, _, _ in peaks}
        pl = pd.DataFrame([{'residue': n, 'N_ppm': y} for n, _, y in peaks])
        table = r2_table_from_fits(t1, rho, pl, omega1_hz=2000.0, carrier_ppm=118.0,
                                   theta_deg=theta, spec_freq_mhz=600.0)
        assert len(table) == len(peaks)
        assert table['R1'].iloc[0] == pytest.approx(TRUE_R1, rel=1e-6)
        assert table['R2'].iloc[0] == pytest.approx(TRUE_R2, rel=0.05)

    def test_carries_the_columns_the_density_step_needs(self):
        peaks = _peaks(3)
        t1 = {n: {'value': 660.0, 'error': 10.0} for n, _, _ in peaks}
        rho = {n: {'value': 80.0, 'error': 2.0} for n, _, _ in peaks}
        pl = pd.DataFrame([{'residue': n, 'N_ppm': y} for n, _, y in peaks])
        table = r2_table_from_fits(t1, rho, pl, omega1_hz=2000.0, carrier_ppm=118.0,
                                   theta_deg=90.0, spec_freq_mhz=600.0)
        for col in ('residue', 'R1', 'R1_err', 'R1rho', 'R1rho_err', 'theta', 'R2', 'R2_err'):
            assert col in table.columns


def _run(*argv):
    return subprocess.run([sys.executable, '-m', 'lunaNMR', *argv],
                          cwd=str(ROOT), capture_output=True, text=True)


class TestCli:
    def _inputs(self, tmp_path, n=12):
        peaks = _peaks(n)
        return (_matrix(tmp_path / 't1.csv', peaks, DELAYS_T1, TRUE_R1),
                _matrix(tmp_path / 'rho.csv', peaks, DELAYS_RHO, TRUE_R2),
                _peak_file(tmp_path / 'peaks.txt', peaks))

    def test_writes_an_r2_table(self, tmp_path):
        """With a spin-lock far stronger than every residue's offset, the tilt is ~90 deg
        for all of them and R2 must come back as R1rho — an answer known independently of
        the code under test. Off-resonance behaviour is covered by TestTiltAngleAlgebra."""
        t1, rho, pl = self._inputs(tmp_path)
        p = _run('dynamixs', 't1rho', '--input', str(rho), '--t1', str(t1),
                 '--peaks', str(pl), '--omega1', '50000', '--carrier', '118',
                 '--theta', '90', '--field-freq', '600',
                 '--out', str(tmp_path / 'out'), '--format', 'json')
        assert p.returncode == 0, p.stderr
        out = json.loads(p.stdout)
        table = pd.read_csv(out['r2_table'])
        assert len(table) == 12
        assert table['R2'].median() == pytest.approx(TRUE_R2, rel=0.05)
        assert table['R1'].median() == pytest.approx(TRUE_R1, rel=0.05)

    def test_dry_run_validates_inputs(self, tmp_path):
        t1, rho, pl = self._inputs(tmp_path)
        p = _run('dynamixs', 't1rho', '--input', str(rho), '--t1', str(tmp_path / 'nope.csv'),
                 '--peaks', str(pl), '--omega1', '2000', '--carrier', '118',
                 '--theta', '90', '--field-freq', '600',
                 '--out', str(tmp_path / 'out'), '--format', 'json', '--dry-run')
        assert p.returncode == 1

    def test_modelfree_accepts_the_r2_table_in_place_of_a_t2_matrix(self, tmp_path):
        t1, rho, pl = self._inputs(tmp_path)
        r2 = _run('dynamixs', 't1rho', '--input', str(rho), '--t1', str(t1), '--peaks', str(pl),
                  '--omega1', '2000', '--carrier', '118', '--theta', '90', '--field-freq', '600',
                  '--out', str(tmp_path / 'out'), '--format', 'json')
        table = json.loads(r2.stdout)['r2_table']
        peaks = _peaks(12)
        sat = tmp_path / 'sat.csv'; unsat = tmp_path / 'unsat.csv'
        sat.write_text('\n'.join(f'{n},{0.8e6},{0.035e6}' for n, _, _ in peaks) + '\n')
        unsat.write_text('\n'.join(f'{n},{1.0e6},{0.044e6}' for n, _, _ in peaks) + '\n')
        p = _run('dynamixs', 'modelfree', '--f1-t1', str(t1), '--f1-r2-table', table,
                 '--f1-noe-sat', str(sat), '--f1-noe-unsat', str(unsat),
                 '--field1-freq', '600', '--out', str(tmp_path / 'mf'),
                 '--prefix', 'rho', '--format', 'json')
        assert p.returncode == 0, p.stderr
        assert json.loads(p.stdout)['n_successful'] >= 1

    def test_t2_matrix_and_r2_table_are_mutually_exclusive(self, tmp_path):
        t1, rho, pl = self._inputs(tmp_path)
        p = _run('dynamixs', 'modelfree', '--f1-t1', str(t1), '--f1-t2', str(rho),
                 '--f1-r2-table', str(rho), '--f1-noe-sat', str(rho),
                 '--f1-noe-unsat', str(rho), '--field1-freq', '600',
                 '--out', str(tmp_path / 'mf'))
        assert p.returncode != 0


class TestOffResonanceMatters:
    """If the per-residue tilt were ignored, an off-resonance residue's R2 would be wrong
    by tens of percent — which is the reason this correction exists."""

    def test_a_residue_far_off_carrier_gets_a_different_r2(self):
        peaks = [('10ValH', 8.0, 118.0), ('40ValH', 8.5, 140.0)]
        t1 = {n: {'value': 1000.0 / TRUE_R1, 'error': 1.0} for n, _, _ in peaks}
        rho = {n: {'value': 1000.0 / TRUE_R2, 'error': 0.5} for n, _, _ in peaks}
        pl = pd.DataFrame([{'residue': n, 'N_ppm': y} for n, _, y in peaks])
        table = r2_table_from_fits(t1, rho, pl, omega1_hz=2000.0, carrier_ppm=118.0,
                                   theta_deg=90.0, spec_freq_mhz=600.0).set_index('residue')
        # identical R1 and R1rho, different chemical shift -> different tilt -> different R2
        assert table.loc['10ValH', 'theta'] == pytest.approx(90.0, abs=0.5)
        assert table.loc['40ValH', 'theta'] < 70.0
        assert table.loc['40ValH', 'R2'] > table.loc['10ValH', 'R2'] * 1.1
