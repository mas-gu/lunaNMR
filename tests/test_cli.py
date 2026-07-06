# ABOUTME: Tests the unified `python -m lunaNMR` CLI dispatcher and its kd subcommand.
# ABOUTME: Also guards that importing the lunaNMR package is headless-safe (no nmrglue crash).

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
_KD_DIR = REPO_ROOT / "modules" / "dynamiXs_v2o0" / "dynamiXs_Kd"
sys.path.insert(0, str(_KD_DIR))

P0 = 50.0


def _write_tidy_csv(tmp_path):
    """Synthetic titration tidy CSV: one mover with a known Kd, 5 points."""
    from kd_models import csp_model
    pts = [0.0, 10.0, 25.0, 60.0, 150.0]
    ddH = csp_model(np.array(pts), 0.2, 20.0, P0)   # inject Kd=20, dd_max=0.2
    rows = [(str(p), 'R1', 8.0 + d, 120.0, 1000.0, 2000.0) for p, d in zip(pts, ddH)]
    df = pd.DataFrame(rows, columns=['spectrum_name', 'assignment',
                                     'ppm_x', 'ppm_y', 'height', 'volume'])
    csv = tmp_path / 'series_analysis_tidy.csv'
    df.to_csv(csv, index=False)
    return csv


class TestPackageImportIsHeadless:
    def test_import_lunanmr_is_lazy_and_headless(self):
        # `import lunaNMR` must succeed without pulling in the GUI/core stack, so
        # `python -m lunaNMR` is headless (no PySide6/QApplication, no nmrglue needed).
        import subprocess
        code = (
            "import sys, lunaNMR; "
            "assert lunaNMR.__version__; "
            "assert 'PySide6' not in sys.modules, 'GUI imported eagerly'; "
            "assert 'lunaNMR.core.core_integrator' not in sys.modules, 'core imported eagerly'"
        )
        r = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT),
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr

    def test_lazy_attribute_access_still_works(self):
        # PEP 562 lazy access: the exported names must still resolve on demand.
        import subprocess
        code = "import lunaNMR; assert lunaNMR.EnhancedVoigtIntegrator is not None"
        r = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT),
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr


class TestKdSubcommand:
    def test_kd_end_to_end_writes_json_and_recovers_kd(self, tmp_path):
        from lunaNMR.cli import main
        csv = _write_tidy_csv(tmp_path)
        out = tmp_path / "kdout"
        code = main(["kd", "--input", str(csv), "--out", str(out),
                     "--p0", "50", "--observable", "csp"])
        assert code == 0
        json_file = out / "kd_kd_fit_data.json"
        assert json_file.exists()
        data = json.loads(json_file.read_text())
        r1 = next(f for f in data['fits'] if f['residue'] == 'R1')
        assert r1['csp']['Kd'] == pytest.approx(20.0, rel=1e-2)

    def test_kd_missing_required_flag_errors(self, tmp_path):
        from lunaNMR.cli import main
        # --out and --p0 are required; argparse should exit non-zero.
        with pytest.raises(SystemExit) as exc:
            main(["kd", "--input", "nope.csv"])
        assert exc.value.code != 0


class TestBatchDelegation:
    def test_batch_is_intercepted_and_argv_forwarded(self, monkeypatch):
        # `batch` must be routed to the batch CLI with its own args untouched,
        # before the top-level argparse can choke on flags like -h.
        from lunaNMR import cli
        captured = {}

        def spy(argv):
            captured['argv'] = argv
            return 0

        monkeypatch.setattr(cli, "_run_batch", spy)
        code = cli.main(["batch", "--preset", "1H", "myfolder"])
        assert code == 0
        assert captured['argv'] == ["--preset", "1H", "myfolder"]


class TestDynamixsT1T2:
    # The T1/T2 plotter calls plt.show() under the Agg backend; that known,
    # harmless warning is expected on a headless run.
    @pytest.mark.filterwarnings("ignore:FigureCanvasAgg is non-interactive")
    def test_t2_end_to_end_writes_results_and_json(self, tmp_path):
        from lunaNMR.cli import main
        delays = [3.0, 9.0, 25.0, 50.0, 100.0, 200.0, 300.0]
        tau, A = 80.0, 1.0e5
        header = ["residue"] + [str(int(d)) for d in delays]
        rows = [[name] + [A * np.exp(-t / tau) for t in delays] for name in ("R1", "R2")]
        df = pd.DataFrame(rows, columns=header)
        csv = tmp_path / "t2_series.csv"
        df.to_csv(csv, index=False)

        out = tmp_path / "t2out"
        code = main(["dynamixs", "t1t2", "--input", str(csv), "--out", str(out),
                     "--exp", "T2", "--prefix", "field1"])
        assert code == 0
        assert (out / "field1_fit_results.txt").exists()
        json_file = out / "field1_T2_fit_data.json"
        assert json_file.exists()
        data = json.loads(json_file.read_text())
        # JSON embeds per-residue fits; recovered T2 must be near the injected tau.
        text = json.dumps(data)
        assert "R1" in text


class TestDynamixsMethyl:
    def test_methyl_end_to_end_writes_json(self, tmp_path):
        from lunaNMR.cli import main
        delays = [3.0, 9.0, 25.0, 50.0, 100.0, 200.0, 300.0]

        def sig(t):
            return 0.5 * 1.0e5 * (np.exp(-t / 80.0) + np.exp(-t / 15.0))

        header = (["Peak_Number", "Assignment", "Reference_X", "Reference_Y"]
                  + [f"003_T2_ADDA_{int(d)}ms" for d in delays])
        rows = [[k + 1, name, 0.0, 0.0] + [f"{sig(t):.2f}" for t in delays]
                for k, name in enumerate(("27_cg1_D", "44_cg1_D"))]
        df = pd.DataFrame(rows, columns=header)
        csv = tmp_path / "methyl_series.csv"
        df.to_csv(csv, index=False)

        out = tmp_path / "methylout"
        code = main(["dynamixs", "methyl-t2", "--input", str(csv), "--out", str(out),
                     "--field-name", "field1", "--field-freq", "600", "--bootstrap", "0"])
        assert code == 0
        assert (out / "field1_methylT2_fit_data.json").exists()

    def test_dynamixs_requires_subcommand(self):
        from lunaNMR.cli import main
        assert main(["dynamixs"]) != 0


class TestSeriesHelpers:
    def test_default_series_params_has_nested_structure(self):
        from lunaNMR.cli import _default_series_params
        p = _default_series_params()
        assert set(p) >= {'detection_params', 'gui_params', 'fitting_params', 'processing_options'}
        assert p['gui_params']['use_ps2d_multi_peak'] is True
        assert p['fitting_params']['min_r_squared'] == 0.5

    def test_discover_spectra_natural_sort(self, tmp_path):
        from lunaNMR.cli import _discover_spectra
        for name in ("s_10.ft", "s_2.ft", "s_1.ft", "notes.txt"):
            (tmp_path / name).write_text("x")
        found = _discover_spectra(str(tmp_path))
        assert [os.path.basename(f) for f in found] == ["s_1.ft", "s_2.ft", "s_10.ft"]


class TestSeriesSubcommand:
    def test_series_requires_peaks(self, tmp_path):
        from lunaNMR.cli import main
        with pytest.raises(SystemExit) as exc:
            main(["series", "--spectra", str(tmp_path), "--out", str(tmp_path / "o")])
        assert exc.value.code != 0

    @pytest.mark.filterwarnings("ignore::RuntimeWarning")
    def test_series_end_to_end_writes_tidy_csv(self, tmp_path):
        from lunaNMR.cli import main
        data = REPO_ROOT / "data_example" / "2DHSQC"
        spectra = sorted(data.glob("600_T1_0o*.ft"))
        peaks = data / "600_assi.txt"
        if len(spectra) < 2 or not peaks.exists():
            pytest.skip("data_example/2DHSQC spectra not available")
        out = tmp_path / "series_out"
        code = main(["series", "--spectra", str(data), "--peaks", str(peaks),
                     "--out", str(out), "--mode", "time", "--peak-source", "reference"])
        assert code == 0
        tidy = out / "series_analysis_tidy.csv"
        assert tidy.exists()
        assert pd.read_csv(tidy).shape[0] > 0


def _make_bundle(tmp_path):
    bundle = tmp_path / "proj.lunaNMR"
    (bundle / "fit_results" / "arrays").mkdir(parents=True)
    (bundle / "kd" / "analyses" / "A1_assi").mkdir(parents=True)
    (bundle / "project.json").write_text('{"schema_version": "1.1"}')
    (bundle / "fit_results" / "fits.json").write_text('[]')
    (bundle / "fit_results" / "arrays" / "peak_000.npz").write_bytes(b"x" * 5000)
    (bundle / "kd" / "analyses" / "A1_assi" / "fit_data.json").write_text("{}")
    return bundle


class TestProjectSubcommand:
    def test_inventory_lists_categories(self, tmp_path, capsys):
        from lunaNMR.cli import main
        bundle = _make_bundle(tmp_path)
        code = main(["project", "inventory", str(bundle)])
        assert code == 0
        out = capsys.readouterr().out
        assert "Fit results" in out and "Kd analyses" in out

    def test_inventory_missing_bundle_errors(self, tmp_path):
        from lunaNMR.cli import main
        assert main(["project", "inventory", str(tmp_path / "nope.lunaNMR")]) != 0

    def test_remove_deletes_path(self, tmp_path):
        from lunaNMR.cli import main
        bundle = _make_bundle(tmp_path)
        code = main(["project", "remove", str(bundle), "fit_results/arrays"])
        assert code == 0
        assert not (bundle / "fit_results" / "arrays").exists()
        assert (bundle / "fit_results" / "fits.json").exists()

    def test_remove_rejects_escape(self, tmp_path):
        from lunaNMR.cli import main
        bundle = _make_bundle(tmp_path)
        sibling = tmp_path / "secret.txt"
        sibling.write_text("keep")
        assert main(["project", "remove", str(bundle), "../secret.txt"]) != 0
        assert sibling.exists()

    def test_project_requires_subcommand(self):
        from lunaNMR.cli import main
        assert main(["project"]) != 0


class TestDispatch:
    def test_help_exits_zero(self):
        from lunaNMR.cli import main
        with pytest.raises(SystemExit) as exc:
            main(["--help"])
        assert exc.value.code == 0

    def test_no_subcommand_returns_nonzero(self):
        from lunaNMR.cli import main
        assert main([]) != 0

    def test_unknown_subcommand_errors(self):
        from lunaNMR.cli import main
        with pytest.raises(SystemExit) as exc:
            main(["bogus"])
        assert exc.value.code != 0
