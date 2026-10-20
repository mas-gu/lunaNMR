# ABOUTME: The 2D fitter and the integrator's loading report through output_manager, not print().
# ABOUTME: A print() reaches the terminal only; the GUI's progress log and a worker's callback never see it.

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from lunaNMR.core import ps2d_2d_fitter
from lunaNMR.core.ps2d_config import set_ps2d_config
from lunaNMR.core.ps2d_2d_fitter import Ps2dMultiPeakFitter2D
from lunaNMR.core.ps2d_style_fitter import voigt_profile_2d
from lunaNMR.integrators.inplace_advanced_nmr_integrator import InPlaceAdvancedNMRIntegrator
from lunaNMR.utils.output_manager import OutputManager


@pytest.fixture
def messages():
    """Every (message, failed) the current thread's output_manager callback receives."""
    received = []
    OutputManager.set_callback(lambda progress, task, msg, failed: received.append((msg, failed)))
    try:
        yield received
    finally:
        OutputManager.set_callback(None)


def _plane():
    ppm_f2 = np.linspace(7.80, 7.60, 64)
    ppm_f1 = np.linspace(122.4, 120.6, 48)
    f2_grid, f1_grid = np.meshgrid(ppm_f2, ppm_f1)
    intensity = voigt_profile_2d(f1_grid, f2_grid, 121.5, 0.06, 7.70, 0.006, 0.18, 0.018, 1.0)
    return f1_grid, f2_grid, intensity


def test_a_failed_region_sum_is_reported_as_a_warning(messages, monkeypatch, capsys):
    """The region sum is measured beside the analytic volume; when it cannot be, the fit
    says so on the output_manager channel and reports no region volume."""
    set_ps2d_config('15N')

    def unmeasurable(*args, **kwargs):
        raise ValueError("mask does not match the data")

    monkeypatch.setattr(ps2d_2d_fitter, 'region_sum_volumes', unmeasurable)
    seed = {'pos_f1': 121.5, 'lw_lor_f1': 0.18, 'lw_gau_f1': 0.06,
            'pos_f2': 7.70, 'lw_lor_f2': 0.018, 'lw_gau_f2': 0.006, 'intensity': 1.0}
    result = Ps2dMultiPeakFitter2D(verbose=False).fit_multi_peak_2d(*_plane(), [seed])

    assert result['success']
    assert ("Warning: region-sum volume unavailable: mask does not match the data", False) \
        in messages
    assert "region-sum volume unavailable" not in capsys.readouterr().out


def test_an_unreadable_peak_list_is_reported_as_an_error(messages, tmp_path, capsys):
    missing = tmp_path / "no_such_peaks.csv"
    assert InPlaceAdvancedNMRIntegrator().load_data(str(missing), str(tmp_path / "x.ft2")) is False

    [(message, failed)] = messages
    assert message.startswith("Error loading peak list: ") and str(missing) in message
    assert failed is True
    assert capsys.readouterr().out == ""


def test_importing_the_base_integrator_first_leaves_the_live_one_built_on_it():
    """lunaNMR.utils imports lunaNMR.core, whose integrator imports this base class, so
    the base reaches output_manager only from inside its methods. Imported at module
    level, the live integrator would find the base half-initialised and fall back to a
    stand-in class whenever the base was imported first."""
    code = ("import lunaNMR.integrators.inplace_advanced_nmr_integrator as base\n"
            "from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator\n"
            "print(issubclass(EnhancedVoigtIntegrator, base.InPlaceAdvancedNMRIntegrator))\n")
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                       cwd=Path(__file__).resolve().parents[1])
    assert p.returncode == 0, p.stderr[-600:]
    assert p.stdout.strip().splitlines() == ["True"], p.stdout
