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

    def test_submodule_access_after_bare_import(self):
        # `import lunaNMR` then `lunaNMR.core` / `lunaNMR.utils` must still resolve,
        # as it did when __init__ imported them eagerly.
        import subprocess
        code = "import lunaNMR; assert lunaNMR.utils is not None; assert lunaNMR.core is not None"
        r = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT),
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr

    def test_unknown_attribute_raises_attributeerror(self):
        # A missing/unavailable name must surface as AttributeError (so hasattr works),
        # not ModuleNotFoundError from a lazy import.
        import lunaNMR
        with pytest.raises(AttributeError):
            lunaNMR.DefinitelyNotARealAttribute
        assert hasattr(lunaNMR, "AlsoNotReal") is False


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
        # Recovered T2 must be near the injected tau (80 ms).
        r1 = next(f for f in data['fits'] if f['residue'] == 'R1')
        assert r1['t2'] == pytest.approx(80.0, rel=0.05)


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
        json_file = out / "field1_methylT2_fit_data.json"
        assert json_file.exists()
        # Recovered slow/fast components must be near the injected bi-exp (80/15 ms).
        r = next(f for f in json.loads(json_file.read_text())['fits'] if f['residue'] == '27_cg1_D')
        assert r['t2_a'] == pytest.approx(80.0, rel=0.05)
        assert r['t2_b'] == pytest.approx(15.0, rel=0.1)

    def test_dynamixs_requires_subcommand(self):
        from lunaNMR.cli import main
        assert main(["dynamixs"]) != 0

    @pytest.mark.filterwarnings("ignore:FigureCanvasAgg is non-interactive")
    def test_t1t2_reads_lunanmr_matrix_with_o_decimal_headers(self, tmp_path):
        # peak_intensity_matrix.csv from a series run has metadata columns plus
        # delay columns labelled like '600_T1_0o0' and dummy_XXX rows. The T1/T2
        # loader must parse the 'o'-decimal headers and drop the dummy rows.
        from lunaNMR.cli import main
        delays = [0.0, 0.1, 0.2, 0.4, 0.8, 1.0, 1.5, 2.0]
        labels = [f"600_T1_{str(d).replace('.', 'o')}" for d in delays]
        A, tau = 1.0e6, 0.5
        header = ["Peak_Number", "Assignment", "Reference_X", "Reference_Y"] + labels
        rows = [
            [1, "3.0", 8.3, 128.3] + [A * np.exp(-t / tau) for t in delays],
            [2, "4.0", 8.6, 122.8] + [A * np.exp(-t / tau) for t in delays],
            [3, "dummy_001", 0.0, 0.0] + [0.0] * len(delays),
        ]
        csv = tmp_path / "peak_intensity_matrix.csv"
        pd.DataFrame(rows, columns=header).to_csv(csv, index=False)

        out = tmp_path / "t1out"
        code = main(["dynamixs", "t1t2", "--input", str(csv), "--out", str(out), "--exp", "T1"])
        assert code == 0
        data = json.loads((out / "field1_T1_fit_data.json").read_text())
        residues = {f["residue"] for f in data["fits"]}
        assert residues == {"3.0", "4.0"}          # dummy row excluded
        # Recovered T1 must be near the injected tau (0.5), not a degenerate huge value.
        assert all(f["t2"] == pytest.approx(0.5, rel=0.1) for f in data["fits"])


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

    def test_discover_spectra_matches_loader_formats(self, tmp_path):
        # Discovery must track NMRFileManager.supported_nmr_formats: include ucsf/pipe/ft3,
        # exclude fid (which the loader rejects).
        from lunaNMR.cli import _discover_spectra
        for name in ("a.ucsf", "b.pipe", "c.ft2", "d.ft3", "junk.fid"):
            (tmp_path / name).write_text("x")
        found = {os.path.basename(f) for f in _discover_spectra(str(tmp_path))}
        assert {"a.ucsf", "b.pipe", "c.ft2", "d.ft3"} <= found
        assert "junk.fid" not in found


class TestRelaxationFlags:
    def test_bootstrap_default_is_usable(self):
        # --error-method bootstrap without --bootstrap must not silently run 0 iterations.
        from lunaNMR.cli import build_parser
        args = build_parser().parse_args(
            ["dynamixs", "t1t2", "--input", "x", "--out", "y", "--exp", "T2"])
        assert args.bootstrap >= 1000


class TestSeriesSubcommand:
    def test_series_requires_peaks(self, tmp_path):
        from lunaNMR.cli import main
        with pytest.raises(SystemExit) as exc:
            main(["series", "--spectra", str(tmp_path), "--out", str(tmp_path / "o")])
        assert exc.value.code != 0

    @pytest.mark.filterwarnings("ignore::RuntimeWarning")
    def test_series_returns_nonzero_when_no_spectra_processed(self, tmp_path):
        from lunaNMR.cli import main
        peaks = REPO_ROOT / "data_example" / "2DHSQC" / "600_assi.txt"
        if not peaks.exists():
            pytest.skip("data_example peak list not available")
        (tmp_path / "bad.ft").write_bytes(b"not a real spectrum")
        code = main(["series", "--spectra", str(tmp_path), "--peaks", str(peaks),
                     "--out", str(tmp_path / "o")])
        assert code == 1

    # Running the real PS2D fitter headlessly emits benign core-engine warnings
    # (numba object-mode fallback / cache) unrelated to the CLI under test.
    @pytest.mark.filterwarnings("ignore::RuntimeWarning")
    @pytest.mark.filterwarnings("ignore::numba.core.errors.NumbaWarning")
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


class TestExportKd:
    def _make_kd_json(self, tmp_path, observable="csp"):
        from lunaNMR.cli import main
        csv = _write_tidy_csv(tmp_path)
        kdout = tmp_path / "kd"
        assert main(["kd", "--input", str(csv), "--out", str(kdout),
                     "--p0", "50", "--observable", observable]) == 0
        return kdout / "kd_kd_fit_data.json"

    def test_export_writes_pdf_and_summary_by_default(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs)]) == 0
        assert (figs / "summary.csv").exists()
        pdf = figs / "csp_fits.pdf"
        assert pdf.exists() and pdf.stat().st_size > 0
        assert not (figs / "csp").exists()          # no per-residue PNG dir in pdf mode

    def test_export_png_format_writes_per_residue(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs_png"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs),
                     "--fig-format", "png"]) == 0
        assert len(list((figs / "csp").glob("*.png"))) >= 1
        assert not (figs / "csp_fits.pdf").exists()

    def test_export_both_formats(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs_both"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs),
                     "--fig-format", "pdf,png"]) == 0
        assert (figs / "csp_fits.pdf").exists()
        assert len(list((figs / "csp").glob("*.png"))) >= 1

    def test_export_summary_only_skips_figures(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs2"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs),
                     "--summary-only"]) == 0
        assert (figs / "summary.csv").exists()
        assert not (figs / "csp").exists()

    def test_export_missing_json_errors(self, tmp_path):
        from lunaNMR.cli import main
        assert main(["export", "kd", "--json", str(tmp_path / "no.json"),
                     "--out", str(tmp_path / "f")]) != 0

    def test_export_requires_subcommand(self):
        from lunaNMR.cli import main
        assert main(["export"]) != 0


class TestErrorHandling:
    def test_kd_missing_input_file_returns_clean_error(self, tmp_path, capsys):
        from lunaNMR.cli import main
        code = main(["kd", "--input", str(tmp_path / "nope.csv"),
                     "--out", str(tmp_path / "o"), "--p0", "50"])
        assert code == 1
        assert capsys.readouterr().err.strip() != ""

    def test_series_missing_peaks_file_returns_clean_error(self, tmp_path):
        from lunaNMR.cli import main
        (tmp_path / "spec.ft").write_bytes(b"\x00")  # a file so discovery passes
        code = main(["series", "--spectra", str(tmp_path),
                     "--peaks", str(tmp_path / "nope.txt"), "--out", str(tmp_path / "o")])
        assert code == 1

    def test_dynamixs_missing_input_returns_clean_error(self, tmp_path):
        from lunaNMR.cli import main
        code = main(["dynamixs", "methyl-t2", "--input", str(tmp_path / "nope.csv"),
                     "--out", str(tmp_path / "o")])
        assert code == 1

    def test_kd_malformed_csv_returns_clean_error(self, tmp_path):
        from lunaNMR.cli import main
        # tidy layout (has spectrum_name) but the assignment column is misnamed →
        # the loader raises KeyError on groupby('assignment').
        df = pd.DataFrame({"spectrum_name": ["0", "1"], "wrong": ["R1", "R1"],
                           "ppm_x": [8.0, 8.1]})
        csv = tmp_path / "bad.csv"
        df.to_csv(csv, index=False)
        code = main(["kd", "--input", str(csv), "--out", str(tmp_path / "o"), "--p0", "50"])
        assert code == 1

    def test_export_kd_malformed_fit_returns_clean_error(self, tmp_path):
        from lunaNMR.cli import main
        # A fit marked success but missing dd_max must not crash with a raw KeyError.
        j = tmp_path / "f.json"
        j.write_text(json.dumps({
            "metadata": {"protein_conc": 50.0},
            "fits": [{"residue": "R1",
                      "csp": {"success": True, "Kd": 20.0, "L": [0, 1, 2], "obs": [0, 0.1, 0.2]}}],
        }))
        code = main(["export", "kd", "--json", str(j), "--out", str(tmp_path / "o")])
        assert code == 1


class TestPackaging:
    def test_pyproject_declares_console_script(self):
        import tomllib
        cfg = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
        assert cfg["project"]["name"] == "lunaNMR"
        assert cfg["project"]["scripts"]["lunanmr"] == "lunaNMR.cli:main"


class TestSeriesParallel:
    def test_default_params_parallel_toggle(self):
        from lunaNMR.cli import _default_series_params
        on = _default_series_params(parallel=True)
        assert on['gui_params']['use_parallel_processing'] is True
        assert on['processing_options']['use_parallel_processing'] is True
        assert _default_series_params()['gui_params']['use_parallel_processing'] is False

    def test_series_accepts_parallel_flag(self):
        from lunaNMR.cli import build_parser
        args = build_parser().parse_args(
            ["series", "--spectra", "x", "--peaks", "p", "--out", "o", "--parallel"])
        assert args.parallel is True


class TestOutputFormat:
    def test_kd_format_json_is_clean_stdout(self, tmp_path, capsys):
        from lunaNMR.cli import main
        csv = _write_tidy_csv(tmp_path)
        code = main(["kd", "--input", str(csv), "--out", str(tmp_path / "o"),
                     "--p0", "50", "--observable", "csp", "--format", "json"])
        assert code == 0
        # Entire stdout must be valid JSON — engine chatter goes to stderr.
        data = json.loads(capsys.readouterr().out)
        assert data["command"] == "kd"
        assert data["n_fitted"] >= 1
        assert data["json_file"].endswith(".json")


class TestDryRun:
    def test_kd_dry_run_validates_without_running(self, tmp_path):
        from lunaNMR.cli import main
        csv = _write_tidy_csv(tmp_path)
        out = tmp_path / "o"
        code = main(["kd", "--input", str(csv), "--out", str(out), "--p0", "50", "--dry-run"])
        assert code == 0
        assert not (out / "kd_kd_fit_data.json").exists()  # nothing was executed

    def test_kd_dry_run_missing_input_returns_1(self, tmp_path):
        from lunaNMR.cli import main
        code = main(["kd", "--input", str(tmp_path / "nope.csv"),
                     "--out", str(tmp_path / "o"), "--p0", "50", "--dry-run"])
        assert code == 1

    def test_dry_run_json_plan(self, tmp_path):
        from lunaNMR.cli import main
        csv = _write_tidy_csv(tmp_path)
        import io
        import contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = main(["kd", "--input", str(csv), "--out", str(tmp_path / "o"),
                         "--p0", "50", "--dry-run", "--format", "json"])
        assert code == 0
        data = json.loads(buf.getvalue())
        assert data["dry_run"] is True


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
