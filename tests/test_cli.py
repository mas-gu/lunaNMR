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
        json_file = out / "data" / "kd_kd_fit_data.json"
        assert json_file.exists()
        data = json.loads(json_file.read_text())
        r1 = next(f for f in data['fits'] if f['residue'] == 'R1')
        assert r1['csp']['Kd'] == pytest.approx(20.0, rel=1e-2)

    def test_kd_missing_required_flag_errors(self, tmp_path):
        from lunaNMR.cli import main
        # --out and --p0 are required; argparse should exit non-zero.
        with pytest.raises(SystemExit) as exc:
            main(["kd", "--input", "nope.csv"])
        assert exc.value.code == 2


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

    def test_t1t2_zero_fit_returns_nonzero(self, tmp_path):
        # No decay within the delay window -> every residue is unreliable and
        # excluded, so nothing is fitted and the command must report failure.
        from lunaNMR.cli import main
        delays = [3.0, 9.0, 25.0, 50.0, 100.0]
        header = ["residue"] + [str(int(d)) for d in delays]
        rows = [[name] + [1.0e5 * np.exp(-t / 1.0e8) for t in delays] for name in ("R1", "R2")]
        pd.DataFrame(rows, columns=header).to_csv(tmp_path / "flat.csv", index=False)
        code = main(["dynamixs", "t1t2", "--input", str(tmp_path / "flat.csv"),
                     "--out", str(tmp_path / "o"), "--exp", "T2", "--prefix", "field1"])
        assert code == 1


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

    def test_methyl_zero_fit_returns_nonzero(self, tmp_path):
        # Fewer time points than the bi-exp has parameters -> every fit fails ->
        # nothing is fitted, so the command must report failure.
        from lunaNMR.cli import main
        delays = [3.0, 100.0]   # 2 points < 3 bi-exp params
        header = (["Peak_Number", "Assignment", "Reference_X", "Reference_Y"]
                  + [f"003_T2_ADDA_{int(d)}ms" for d in delays])
        rows = [[k + 1, name, 0.0, 0.0] + [1.0e5, 5.0e4]
                for k, name in enumerate(("27_cg1_D", "44_cg1_D"))]
        pd.DataFrame(rows, columns=header).to_csv(tmp_path / "m.csv", index=False)
        code = main(["dynamixs", "methyl-t2", "--input", str(tmp_path / "m.csv"),
                     "--out", str(tmp_path / "o"), "--field-name", "field1",
                     "--field-freq", "600", "--bootstrap", "0"])
        assert code == 1

    def test_dynamixs_requires_subcommand(self):
        from lunaNMR.cli import main
        assert main(["dynamixs"]) == 2

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
        # 1H/15N search window starting value is 0.070 ppm in both dimensions.
        assert p['detection_params']['search_window_x'] == 0.070
        assert p['detection_params']['search_window_y'] == 0.070

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
        for name in ("a.ucsf", "b.pipe", "c.ft2", "d.ft3", "e.ft1", "junk.fid"):
            (tmp_path / name).write_text("x")
        found = {os.path.basename(f) for f in _discover_spectra(str(tmp_path))}
        assert {"a.ucsf", "b.pipe", "c.ft2", "d.ft3", "e.ft1"} <= found
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
        assert exc.value.code == 2

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


class TestIndependentModeMatchesTheGui:
    """--peak-source independent names the GUI's Independent Mode, but the GUI sets two
    processing options for it that the CLI left at their defaults. Same flag, different
    algorithm, and neither side said so.
    """

    def _options_for(self, monkeypatch, tmp_path, peak_source):
        """Capture the processing_options the CLI hands the processor."""
        import lunaNMR.processors.multi_spectrum_processor as msp
        captured = {}

        class _Spy:
            def __init__(self, params):
                captured.update(params['processing_options'])

            def process_nmr_series(self, *a, **kw):
                raise SystemExit(0)      # stop before any fitting happens

        monkeypatch.setattr(msp, 'MultiSpectrumProcessor', _Spy)
        (tmp_path / "s_8ms.ft").write_bytes(b"")
        peaks = tmp_path / "p.txt"
        peaks.write_text("Assignment\tPosition_X\tPosition_Y\nK3\t8.0\t120.0\n")
        from lunaNMR.cli import main
        with pytest.raises(SystemExit):
            main(["series", "--spectra", str(tmp_path), "--peaks", str(peaks),
                  "--out", str(tmp_path / "o"), "--peak-source", peak_source])
        return captured

    def test_independent_reruns_adaptive_and_uses_the_original_reference(
            self, monkeypatch, tmp_path):
        opts = self._options_for(monkeypatch, tmp_path, "independent")
        assert opts['rerun_adaptive_per_spectrum'] is True
        assert opts['use_original_reference_for_detection'] is True

    @pytest.mark.parametrize("source", ["reference", "cascade", "detected"])
    def test_the_other_sources_are_unchanged(self, monkeypatch, tmp_path, source):
        opts = self._options_for(monkeypatch, tmp_path, source)
        assert opts['rerun_adaptive_per_spectrum'] is False
        assert opts['use_original_reference_for_detection'] is False


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
        assert main(["project", "inventory", str(tmp_path / "nope.lunaNMR")]) == 1

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
        assert main(["project", "remove", str(bundle), "../secret.txt"]) == 1
        assert sibling.exists()

    def test_project_requires_subcommand(self):
        from lunaNMR.cli import main
        assert main(["project"]) == 2


class TestExportKd:
    def _make_kd_json(self, tmp_path, observable="csp"):
        from lunaNMR.cli import main
        csv = _write_tidy_csv(tmp_path)
        kdout = tmp_path / "kd"
        assert main(["kd", "--input", str(csv), "--out", str(kdout),
                     "--p0", "50", "--observable", observable]) == 0
        return kdout / "data" / "kd_kd_fit_data.json"

    def test_export_writes_pdf_and_summary_by_default(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs)]) == 0
        assert (figs / "data" / "summary.csv").exists()
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
        assert (figs / "data" / "summary.csv").exists()
        assert not (figs / "csp").exists()

    def test_export_ref_vs_point_pdf_by_default(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs_ref"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs)]) == 0
        # both the per-residue curves and the new ref->point bars are emitted
        assert (figs / "csp_fits.pdf").exists()
        ref = figs / "csp_ref_vs_point.pdf"
        assert ref.exists() and ref.stat().st_size > 0

    def test_export_kind_curves_only_skips_ref_bars(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs_curves"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs),
                     "--kind", "curves"]) == 0
        assert (figs / "csp_fits.pdf").exists()
        assert not (figs / "csp_ref_vs_point.pdf").exists()

    def test_export_ref_bars_cover_both_observables_by_default(self, tmp_path):
        # Ref->point bars are model-free (from the raw series), so a CSP-only fit still
        # gets BOTH csp and intensity ref-bars by default (parity with the GUI).
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs_bothref"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs)]) == 0
        assert (figs / "csp_ref_vs_point.pdf").exists()
        assert (figs / "intensity_ref_vs_point.pdf").exists()

    def test_export_observable_flag_limits_ref_bars(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs_reffilter"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs),
                     "--observable", "csp"]) == 0
        assert (figs / "csp_ref_vs_point.pdf").exists()
        assert not (figs / "intensity_ref_vs_point.pdf").exists()

    def test_export_kind_ref_bars_only_skips_curves(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs_refonly"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs),
                     "--kind", "ref-bars"]) == 0
        assert (figs / "csp_ref_vs_point.pdf").exists()
        assert not (figs / "csp_fits.pdf").exists()

    def test_export_kd_vs_residue_pdf_by_default(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs_kdbars"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs)]) == 0
        pdf = figs / "csp_kd_vs_residue.pdf"
        assert pdf.exists() and pdf.stat().st_size > 0

    def test_export_kind_kd_bars_only_skips_curves_and_ref(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs_kdonly"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs),
                     "--kind", "kd-bars"]) == 0
        assert (figs / "csp_kd_vs_residue.pdf").exists()
        assert not (figs / "csp_fits.pdf").exists()
        assert not (figs / "csp_ref_vs_point.pdf").exists()

    def test_export_prefix_is_prepended_to_all_outputs(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs_prefixed"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs),
                     "--prefix", "DNAJA1_HSPA8"]) == 0
        assert (figs / "data" / "DNAJA1_HSPA8_summary.csv").exists()
        assert (figs / "DNAJA1_HSPA8_csp_fits.pdf").exists()
        assert (figs / "DNAJA1_HSPA8_csp_ref_vs_point.pdf").exists()
        assert (figs / "data" / "DNAJA1_HSPA8_csp_ref_vs_point.csv").exists()
        assert (figs / "DNAJA1_HSPA8_csp_kd_vs_residue.pdf").exists()
        # no unprefixed duplicates left behind
        assert not (figs / "data" / "summary.csv").exists()
        assert not (figs / "csp_fits.pdf").exists()
        assert not (figs / "csp_ref_vs_point.pdf").exists()
        assert not (figs / "csp_kd_vs_residue.pdf").exists()

    def test_export_no_prefix_keeps_unprefixed_names(self, tmp_path):
        # Omitting --prefix must stay fully backward compatible.
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs_noprefix"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs)]) == 0
        assert (figs / "data" / "summary.csv").exists()
        assert (figs / "csp_fits.pdf").exists()

    def _make_multi_kd_json(self, tmp_path):
        """A 3-residue titration (intensity decay, shared Kd) so a global fit fires."""
        from kd_models import intensity_decay
        from lunaNMR.cli import main
        pts = [0.0, 10.0, 25.0, 60.0, 150.0, 300.0]
        specs = {'A1': (1000.0, 50.0), 'K2': (800.0, 200.0), 'G3': (1200.0, 0.0)}
        rows = []
        for name, (i0, iinf) in specs.items():
            hs = intensity_decay(np.array(pts), i0, iinf, 40.0)
            rows += [(str(p), name, 8.0, 120.0, h, 2 * h) for p, h in zip(pts, hs)]
        df = pd.DataFrame(rows, columns=['spectrum_name', 'assignment',
                                         'ppm_x', 'ppm_y', 'height', 'volume'])
        csv = tmp_path / "series_analysis_tidy.csv"
        df.to_csv(csv, index=False)
        out = tmp_path / "kdm"
        assert main(["kd", "--input", str(csv), "--out", str(out), "--p0", "50",
                     "--conc", "0,10,25,60,150,300", "--observable", "intensity"]) == 0
        return out / "data" / "kd_kd_fit_data.json"

    def test_export_global_fit_pdf_by_default(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_multi_kd_json(tmp_path)
        figs = tmp_path / "figs_global"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs)]) == 0
        pdf = figs / "intensity_global_fit.pdf"
        assert pdf.exists() and pdf.stat().st_size > 0

    def test_export_kind_global_fit_only_skips_others(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_multi_kd_json(tmp_path)
        figs = tmp_path / "figs_globalonly"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs),
                     "--kind", "global-fit"]) == 0
        assert (figs / "intensity_global_fit.pdf").exists()
        assert not (figs / "intensity_fits.pdf").exists()
        assert not (figs / "intensity_kd_vs_residue.pdf").exists()

    def test_export_prefix_applies_to_global_fit(self, tmp_path):
        from lunaNMR.cli import main
        jf = self._make_multi_kd_json(tmp_path)
        figs = tmp_path / "figs_globalprefixed"
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs),
                     "--prefix", "DNAJB12_HSPA6"]) == 0
        assert (figs / "DNAJB12_HSPA6_intensity_global_fit.pdf").exists()
        assert not (figs / "intensity_global_fit.pdf").exists()

    def test_export_kd_sets_pdf_fonttype_42_for_illustrator(self, tmp_path):
        # pdf.fonttype=3 (matplotlib's default) embeds Type 3/bitmap fonts that Illustrator
        # doesn't recognize as real text; 42 (TrueType) keeps labels editable on import.
        import matplotlib
        from lunaNMR.cli import main
        jf = self._make_kd_json(tmp_path, "csp")
        figs = tmp_path / "figs_fonttype"
        matplotlib.rcParams['pdf.fonttype'] = 3
        assert main(["export", "kd", "--json", str(jf), "--out", str(figs)]) == 0
        assert matplotlib.rcParams['pdf.fonttype'] == 42

    def test_export_missing_json_errors(self, tmp_path):
        from lunaNMR.cli import main
        assert main(["export", "kd", "--json", str(tmp_path / "no.json"),
                     "--out", str(tmp_path / "f")]) != 0

    def test_export_requires_subcommand(self):
        from lunaNMR.cli import main
        assert main(["export"]) == 2


class TestIntensityPanelAxisSharing:
    """Every intensity curve panel (raw-scale or already-normalized) must plot as a true
    I/I(0) ratio and share one 0-1-ish y-axis, so raw-scale residues (I0 in the hundreds)
    don't get clipped to invisible flat lines by a blind 0-1 clamp."""

    def test_build_kd_panel_normalizes_raw_scale_intensity_to_i0(self):
        from lunaNMR.cli import _build_kd_panel
        fit = {'I0': 500.0, 'I_inf': 50.0, 'Kd': 40.0, 'r_squared': 0.9}
        L = np.array([0.0, 10.0, 50.0])
        y = np.array([500.0, 300.0, 100.0])
        Ld = np.linspace(0.0, 50.0, 5)
        panel = _build_kd_panel('R1', fit, L, y, Ld, 'intensity', P0=None)
        assert panel['ylabel'] == 'I / I(0)'
        np.testing.assert_allclose(panel['y'], y / 500.0)
        assert panel['yc'].max() <= 1.0 + 1e-9

    def test_build_kd_panel_leaves_csp_unnormalized(self):
        from lunaNMR.cli import _build_kd_panel
        fit = {'dd_max': 0.5, 'Kd': 40.0, 'r_squared': 0.9}
        L = np.array([0.0, 10.0, 50.0])
        y = np.array([0.0, 0.2, 0.4])
        Ld = np.linspace(0.0, 50.0, 5)
        panel = _build_kd_panel('R1', fit, L, y, Ld, 'csp', P0=50.0)
        assert panel['ylabel'] == 'CSP (ppm)'
        np.testing.assert_allclose(panel['y'], y)

    def test_draw_kd_panel_applies_shared_ylim(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from lunaNMR.cli import _draw_kd_panel, _INTENSITY_YLIM
        fig, ax = plt.subplots()
        panel = {'residue': 'R1', 'L': np.array([0.0, 1.0]), 'y': np.array([1.0, 0.5]),
                'Ld': np.array([0.0, 1.0]), 'yc': np.array([1.0, 0.5]),
                'ylabel': 'I / I(0)', 'kd': 40.0, 'r2': 0.9}
        _draw_kd_panel(ax, panel, ylim=_INTENSITY_YLIM)
        assert ax.get_ylim() == _INTENSITY_YLIM
        plt.close(fig)

    def test_draw_kd_panel_no_ylim_autoscales(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from lunaNMR.cli import _draw_kd_panel
        fig, ax = plt.subplots()
        panel = {'residue': 'R1', 'L': np.array([0.0, 1.0]), 'y': np.array([10.0, 5.0]),
                'Ld': np.array([0.0, 1.0]), 'yc': np.array([10.0, 5.0]),
                'ylabel': 'CSP (ppm)', 'kd': 40.0, 'r2': 0.9}
        _draw_kd_panel(ax, panel)
        assert ax.get_ylim() != (0.0, 1.0)
        plt.close(fig)

    def test_export_intensity_curves_pdf_uses_shared_ylim_for_raw_scale_data(self, tmp_path):
        # End-to-end: a raw-scale (I0=500-ish) intensity fit must not crash and must
        # still write a real PDF (the actual per-axis ylim isn't inspectable once saved,
        # but the unit tests above cover the normalization + ylim application directly).
        from kd_models import intensity_decay
        from lunaNMR.cli import main
        pts = [0.0, 10.0, 25.0, 60.0, 150.0, 300.0]
        hs = intensity_decay(np.array(pts), 500.0, 50.0, 40.0)
        rows = [(str(p), 'A1', 8.0, 120.0, h, 2 * h) for p, h in zip(pts, hs)]
        df = pd.DataFrame(rows, columns=['spectrum_name', 'assignment',
                                         'ppm_x', 'ppm_y', 'height', 'volume'])
        csv = tmp_path / "series_analysis_tidy.csv"
        df.to_csv(csv, index=False)
        kdout = tmp_path / "kd"
        assert main(["kd", "--input", str(csv), "--out", str(kdout), "--p0", "50",
                     "--conc", ",".join(str(p) for p in pts),
                     "--observable", "intensity"]) == 0
        figs = tmp_path / "figs"
        assert main(["export", "kd", "--json", str(kdout / "data" / "kd_kd_fit_data.json"),
                     "--out", str(figs), "--observable", "intensity"]) == 0
        pdf = figs / "intensity_fits.pdf"
        assert pdf.exists() and pdf.stat().st_size > 0


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
        assert not (out / "data" / "kd_kd_fit_data.json").exists()  # nothing was executed

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


class TestADryRunAgreesWithTheRealRun:
    """A dry-run exists so a caller can find out whether the real run will work. One
    that exits 0 where the real run exits 1 is worse than none: it certifies a command
    that cannot run.
    """

    def _f1(self, tmp_path):
        for name in ("t1", "t2", "sat", "unsat"):
            (tmp_path / f"{name}.csv").write_text("8,1000\n51,500\n")
        return ["--f1-t1", str(tmp_path / "t1.csv"), "--f1-t2", str(tmp_path / "t2.csv"),
                "--f1-noe-sat", str(tmp_path / "sat.csv"),
                "--f1-noe-unsat", str(tmp_path / "unsat.csv")]

    def test_density_dual_without_input2_fails_both_ways(self, tmp_path):
        """--dual needs --input2, but that was checked after the dry-run returned."""
        from lunaNMR.cli import main
        table = tmp_path / "f1.csv"
        table.write_text("Residue,R1,R1err,R2,R2err,hetNOE,hetNOEerr\n1,1.0,0.1,10.0,0.5,0.8,0.02\n")
        argv = ["dynamixs", "density", "--dual", "--input", str(table),
                "--out", str(tmp_path / "o")]
        assert main(argv + ["--dry-run"]) == main(argv) == 1

    def test_density_dual_dry_run_names_the_missing_flag(self, tmp_path, capsys):
        from lunaNMR.cli import main
        table = tmp_path / "f1.csv"
        table.write_text("Residue,R1,R1err,R2,R2err,hetNOE,hetNOEerr\n1,1.0,0.1,10.0,0.5,0.8,0.02\n")
        main(["dynamixs", "density", "--dual", "--input", str(table),
              "--out", str(tmp_path / "o"), "--dry-run"])
        assert "--input2" in capsys.readouterr().err

    def test_modelfree_dual_dry_run_names_the_missing_flags(self, tmp_path, capsys):
        """It already exited 1, but by way of os.path.exists(None) raising TypeError:
        'stat: path should be string, bytes, os.PathLike or integer, not NoneType' --
        which names neither --dual nor the files it wanted."""
        from lunaNMR.cli import main
        assert main(["dynamixs", "modelfree", "--dual", *self._f1(tmp_path),
                     "--field2-freq", "800",   # required with --dual; not what this tests
                     "--out", str(tmp_path / "o"), "--dry-run"]) == 1
        assert "--f2" in capsys.readouterr().err

    def test_a_dry_run_plan_names_the_full_subcommand(self, tmp_path, capsys):
        """The one universal key has to discriminate: nested handlers reported
        'dynamixs' in dry-run and 'dynamixs density' in the real run."""
        from lunaNMR.cli import main
        table = tmp_path / "f1.csv"
        table.write_text("Residue,R1,R1err,R2,R2err,hetNOE,hetNOEerr\n1,1.0,0.1,10.0,0.5,0.8,0.02\n")
        main(["dynamixs", "density", "--input", str(table), "--out", str(tmp_path / "o"),
              "--dry-run", "--format", "json"])
        assert json.loads(capsys.readouterr().out)["command"] == "dynamixs density"


class TestDynamixsHetnoe:
    def test_hetnoe_ratio(self, tmp_path):
        from lunaNMR.cli import main
        (tmp_path / "sat.csv").write_text("R1,800,16\nR2,1000,20\n")
        (tmp_path / "unsat.csv").write_text("R1,1000,20\nR2,1000,20\n")
        out = tmp_path / "h"
        assert main(["dynamixs", "hetnoe", "--sat", str(tmp_path / "sat.csv"),
                     "--unsat", str(tmp_path / "unsat.csv"), "--out", str(out)]) == 0
        df = pd.read_csv(out / "field1_hetnoe.csv")
        assert float(df.loc[df.Residue == "R1", "hetNOE"].iloc[0]) == pytest.approx(0.8)
        # QC plot written alongside the CSV
        assert (out / "field1_hetnoe.pdf").exists()

    def test_hetnoe_missing_file_errors(self, tmp_path):
        from lunaNMR.cli import main
        (tmp_path / "unsat.csv").write_text("R1,1000\n")
        assert main(["dynamixs", "hetnoe", "--sat", str(tmp_path / "nope.csv"),
                     "--unsat", str(tmp_path / "unsat.csv"), "--out", str(tmp_path / "h")]) == 1

    def test_hetnoe_dry_run(self, tmp_path):
        from lunaNMR.cli import main
        (tmp_path / "sat.csv").write_text("R1,800\n")
        (tmp_path / "unsat.csv").write_text("R1,1000\n")
        out = tmp_path / "h"
        assert main(["dynamixs", "hetnoe", "--sat", str(tmp_path / "sat.csv"),
                     "--unsat", str(tmp_path / "unsat.csv"),
                     "--out", str(out), "--dry-run"]) == 0
        assert not out.exists()   # dry-run must not touch the filesystem


class TestDynamixsDensity:
    def _table(self, path, n=4):
        rows = [(f"{10 + i}.VAL", 1.4 + 0.01 * i, 0.06, 8.0 + 0.1 * i, 0.3,
                 0.70 + 0.01 * i, 0.014) for i in range(n)]
        pd.DataFrame(rows, columns=["Residue", "R1", "R1err", "R2", "R2err",
                                    "hetNOE", "hetNOEerr"]).to_csv(path, index=False)

    def test_single_field_density(self, tmp_path):
        from lunaNMR.cli import main
        t = tmp_path / "f1.csv"
        self._table(t)
        out = tmp_path / "dens"
        code = main(["dynamixs", "density", "--input", str(t), "--out", str(out),
                     "--field1-freq", "600", "--no-parallel", "--no-plot"])
        assert code == 0
        df = pd.read_csv(out / "spectral_density_results.csv")
        assert len(df) == 4
        assert {"J0", "JwN", "S2"} <= set(df.columns)
        assert (df["S2"].between(0, 1.05)).all()   # order parameters are physical

    def test_density_dry_run(self, tmp_path):
        from lunaNMR.cli import main
        t = tmp_path / "f1.csv"
        self._table(t)
        out = tmp_path / "dens"
        assert main(["dynamixs", "density", "--input", str(t), "--out", str(out),
                     "--field1-freq", "600", "--dry-run"]) == 0
        assert not out.exists()


class TestDynamixsModelfree:
    N = 12   # pipeline requires >= 10 residues

    def _matrix(self, path, delays, tau):
        # Noise is not decoration here. A perfectly noiseless decay fits exactly, so the
        # fit error is zero, and the spectral-density step rejects a zero error as
        # unweightable (`r2_err <= 0` -> "invalid values"). Real data always has noise;
        # a noiseless matrix is not a case the pipeline can or should handle.
        rng = np.random.default_rng(0)
        header = ["Peak_Number", "Assignment", "Reference_X", "Reference_Y"] + [str(d) for d in delays]
        rows = [[i + 1, f"R{i + 1}", 8.0, 120.0]
                + [1e5 * np.exp(-t / tau) * (1.0 + rng.normal(0, 0.01)) for t in delays]
                for i in range(self.N)]
        pd.DataFrame(rows, columns=header).to_csv(path, index=False)

    def _noe(self, path, value):
        path.write_text("".join(f"R{i + 1},{value}\n" for i in range(self.N)))

    def test_single_field_modelfree_end_to_end(self, tmp_path, monkeypatch):
        from lunaNMR.cli import main
        monkeypatch.chdir(tmp_path)               # pipeline writes intermediate txt to CWD
        delays = [10, 30, 60, 100, 150]           # ms
        self._matrix(tmp_path / "t1.csv", delays, 800.0)   # T1 ~800 ms
        self._matrix(tmp_path / "t2.csv", delays, 100.0)   # T2 ~100 ms
        self._noe(tmp_path / "sat.csv", 800)
        self._noe(tmp_path / "unsat.csv", 1000)
        out = tmp_path / "mf"
        code = main(["dynamixs", "modelfree",
                     "--f1-t1", str(tmp_path / "t1.csv"), "--f1-t2", str(tmp_path / "t2.csv"),
                     "--f1-noe-sat", str(tmp_path / "sat.csv"),
                     "--f1-noe-unsat", str(tmp_path / "unsat.csv"),
                     "--field1-freq", "600", "--method", "single_087", "--out", str(out)])
        assert code == 0
        import glob
        basic = glob.glob(str(out / "*basic*.csv"))
        assert basic, "no spectral-density output written"
        df = pd.read_csv(basic[0])
        assert "S2" in df.columns and len(df) >= 1

    def test_validation_warnings_are_not_discarded_under_json(self, tmp_path, monkeypatch):
        """validate_relaxation_rates is the only check that catches a units error, and
        it reports through the progress callback. Under --format json that callback was
        a null lambda -- the mode both prompts prescribe, so the check ran into nothing.

        R2 < R1 is unphysical for isotropic tumbling, so a T2 longer than T1 is a case
        it must flag."""
        from lunaNMR.cli import main
        monkeypatch.chdir(tmp_path)
        delays = [10, 30, 60, 100, 150]
        self._matrix(tmp_path / "t1.csv", delays, 100.0)    # T1 shorter ...
        self._matrix(tmp_path / "t2.csv", delays, 800.0)    # ... than T2: R2 < R1
        self._noe(tmp_path / "sat.csv", 800)
        self._noe(tmp_path / "unsat.csv", 1000)
        argv = ["dynamixs", "modelfree",
                "--f1-t1", str(tmp_path / "t1.csv"), "--f1-t2", str(tmp_path / "t2.csv"),
                "--f1-noe-sat", str(tmp_path / "sat.csv"),
                "--f1-noe-unsat", str(tmp_path / "unsat.csv"),
                "--field1-freq", "600", "--out", str(tmp_path / "mf"),
                "--format", "json"]
        import contextlib, io
        out_buf, err_buf = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out_buf), contextlib.redirect_stderr(err_buf):
            code = main(argv)
        assert code == 0, err_buf.getvalue()[-800:]
        summary = json.loads(out_buf.getvalue())
        assert "validation" in summary
        assert summary["validation"]["r2_less_than_r1"] > 0

    def test_the_validation_report_reaches_stderr_not_a_black_hole(self, tmp_path, monkeypatch):
        from lunaNMR.cli import main
        monkeypatch.chdir(tmp_path)
        delays = [10, 30, 60, 100, 150]
        self._matrix(tmp_path / "t1.csv", delays, 100.0)
        self._matrix(tmp_path / "t2.csv", delays, 800.0)
        self._noe(tmp_path / "sat.csv", 800)
        self._noe(tmp_path / "unsat.csv", 1000)
        import contextlib, io
        err_buf = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err_buf):
            main(["dynamixs", "modelfree",
                  "--f1-t1", str(tmp_path / "t1.csv"), "--f1-t2", str(tmp_path / "t2.csv"),
                  "--f1-noe-sat", str(tmp_path / "sat.csv"),
                  "--f1-noe-unsat", str(tmp_path / "unsat.csv"),
                  "--field1-freq", "600", "--out", str(tmp_path / "mf"),
                  "--format", "json"])
        assert "Validat" in err_buf.getvalue()

    def test_dual_requires_the_second_field_frequency(self, tmp_path):
        """--field2-freq defaults to 700, so an 800 MHz second field was analysed at
        700 and surfaced as a tau_c mismatch that the QC table blames on the data. A
        second field frequency is never a safe guess; there is no conventional value."""
        from lunaNMR.cli import main
        for name in ("t1", "t2", "sat", "unsat"):
            (tmp_path / f"f1_{name}.csv").write_text("8,1000\n")
            (tmp_path / f"f2_{name}.csv").write_text("8,1000\n")
        argv = ["dynamixs", "modelfree", "--dual",
                *sum([[f"--f1-{f}", str(tmp_path / f"f1_{n}.csv")]
                      for f, n in (("t1", "t1"), ("t2", "t2"),
                                   ("noe-sat", "sat"), ("noe-unsat", "unsat"))], []),
                *sum([[f"--f2-{f}", str(tmp_path / f"f2_{n}.csv")]
                      for f, n in (("t1", "t1"), ("t2", "t2"),
                                   ("noe-sat", "sat"), ("noe-unsat", "unsat"))], []),
                "--field1-freq", "600", "--out", str(tmp_path / "o"), "--dry-run"]
        with pytest.raises(SystemExit) as exc:
            main(argv)
        assert exc.value.code == 2

    def test_single_field_does_not_require_it(self, tmp_path):
        """Only --dual reads a second frequency, so it must not become mandatory."""
        from lunaNMR.cli import main
        for name in ("t1", "t2", "sat", "unsat"):
            (tmp_path / f"{name}.csv").write_text("8,1000\n")
        assert main(["dynamixs", "modelfree",
                     "--f1-t1", str(tmp_path / "t1.csv"), "--f1-t2", str(tmp_path / "t2.csv"),
                     "--f1-noe-sat", str(tmp_path / "sat.csv"),
                     "--f1-noe-unsat", str(tmp_path / "unsat.csv"),
                     "--field1-freq", "600", "--out", str(tmp_path / "o"),
                     "--dry-run"]) == 0

    def test_modelfree_t1_units_seconds(self, tmp_path, monkeypatch):
        # A T1 series labelled in seconds must be honoured via --f1-t1-units s.
        # Without it the hard-coded ms conversion makes R1 ~1000x too large.
        from lunaNMR.cli import main
        monkeypatch.chdir(tmp_path)
        self._matrix(tmp_path / "t1.csv", [0.02, 0.05, 0.1, 0.2, 0.4, 0.8], 0.6)  # T1 ~0.6 s
        self._matrix(tmp_path / "t2.csv", [10, 30, 60, 100, 150], 100.0)          # T2 ~100 ms
        self._noe(tmp_path / "sat.csv", 800)
        self._noe(tmp_path / "unsat.csv", 1000)
        out = tmp_path / "mf"
        code = main(["dynamixs", "modelfree",
                     "--f1-t1", str(tmp_path / "t1.csv"), "--f1-t2", str(tmp_path / "t2.csv"),
                     "--f1-noe-sat", str(tmp_path / "sat.csv"),
                     "--f1-noe-unsat", str(tmp_path / "unsat.csv"),
                     "--field1-freq", "600", "--f1-t1-units", "s", "--out", str(out)])
        assert code == 0
        import glob
        df = pd.read_csv(glob.glob(str(out / "*basic*.csv"))[0])
        assert df["R1"].median() < 10   # physical (~1.7 s^-1), not ~1700 under wrong units

    def test_modelfree_default_method_writes_only_to_out(self, tmp_path, monkeypatch):
        # With no --method and single field (--dual off), the method must derive to
        # single_087 (not the dual default that crashes), and every file the pipeline
        # emits must land under --out even though the CWD is elsewhere.
        from lunaNMR.cli import main
        cwd = tmp_path / "cwd"
        cwd.mkdir()
        monkeypatch.chdir(cwd)
        delays = [10, 30, 60, 100, 150]
        self._matrix(tmp_path / "t1.csv", delays, 800.0)
        self._matrix(tmp_path / "t2.csv", delays, 100.0)
        self._noe(tmp_path / "sat.csv", 800)
        self._noe(tmp_path / "unsat.csv", 1000)
        out = tmp_path / "mf"
        code = main(["dynamixs", "modelfree",
                     "--f1-t1", str(tmp_path / "t1.csv"), "--f1-t2", str(tmp_path / "t2.csv"),
                     "--f1-noe-sat", str(tmp_path / "sat.csv"),
                     "--f1-noe-unsat", str(tmp_path / "unsat.csv"),
                     "--field1-freq", "600", "--out", str(out)])
        assert code == 0
        import glob
        assert glob.glob(str(out / "*basic*.csv")), "density output not under --out"
        # The CWD must stay clean: no stray intermediate txt / plot files leaked there.
        assert not glob.glob(str(cwd / "*.txt"))
        assert not glob.glob(str(cwd / "*.pdf"))

    def test_modelfree_dry_run(self, tmp_path):
        from lunaNMR.cli import main
        for name in ("t1.csv", "t2.csv", "sat.csv", "unsat.csv"):
            (tmp_path / name).write_text("x\n")
        out = tmp_path / "mf"
        code = main(["dynamixs", "modelfree",
                     "--f1-t1", str(tmp_path / "t1.csv"), "--f1-t2", str(tmp_path / "t2.csv"),
                     "--f1-noe-sat", str(tmp_path / "sat.csv"),
                     "--f1-noe-unsat", str(tmp_path / "unsat.csv"),
                     "--field1-freq", "600", "--out", str(out), "--dry-run"])
        assert code == 0
        assert not out.exists()


class TestDispatch:
    def test_help_exits_zero(self):
        from lunaNMR.cli import main
        with pytest.raises(SystemExit) as exc:
            main(["--help"])
        assert exc.value.code == 0

    def test_no_subcommand_returns_nonzero(self):
        from lunaNMR.cli import main
        assert main([]) == 2

    def test_unknown_subcommand_errors(self):
        from lunaNMR.cli import main
        with pytest.raises(SystemExit) as exc:
            main(["bogus"])
        assert exc.value.code == 2
