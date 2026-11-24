# ABOUTME: A copy of the lunaNMR package alone, as `pip install .` leaves it: what needs the
# ABOUTME: analysis modules beside it says to install in editable mode, the rest still runs.

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "data_example" / "2DHSQC"
ERROR_TYPE = "AnalysisModulesNotFound"
REMEDY = "pip install -e ."


def _package_copy(site):
    shutil.copytree(REPO / "lunaNMR", site / "lunaNMR",
                    ignore=shutil.ignore_patterns("__pycache__", ".DS_Store"))
    return site


@pytest.fixture(scope="module")
def package_only(tmp_path_factory):
    """A folder holding the lunaNMR package and no `modules/` beside it."""
    return _package_copy(tmp_path_factory.mktemp("site_packages"))


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    """Files that exist, so a command gets as far as needing its analysis module."""
    folder = tmp_path_factory.mktemp("inputs")
    (folder / "matrix.csv").write_text("Residue,10,20,40\n1AlaH,100,80,50\n")
    (folder / "peaks.txt").write_text("Assignment,Position_X,Position_Y\n1AlaH,8.0,120.0\n")
    (folder / "fit.json").write_text("{}")
    run = folder / "run"
    (run / "data").mkdir(parents=True)
    (run / "data" / "kd_kd_fit_data.json").write_text('{"metadata": {}, "residues": {}}')
    (folder / "spectra_1d").mkdir()
    return folder


@pytest.fixture(scope="module")
def example_series(tmp_path_factory):
    """Two spectra of the shipped T1 series and a peak list, in a folder of their own."""
    names = ["600_T1_0o0.ft", "600_T1_0o1.ft", "r_600_peaks.txt"]
    if not all((EXAMPLE / name).exists() for name in names):
        pytest.skip("data_example/2DHSQC not available")
    folder = tmp_path_factory.mktemp("example_series")
    for name in names:
        (folder / name).symlink_to(EXAMPLE / name)
    return folder


def _python(site, *argv, cwd):
    env = dict(os.environ, PYTHONPATH=str(site), QT_QPA_PLATFORM="offscreen")
    return subprocess.run([sys.executable, *argv], cwd=str(cwd), env=env,
                          capture_output=True, text=True, timeout=600)


def _cli(site, *argv, cwd):
    return _python(site, "-m", "lunaNMR", *argv, cwd=cwd)


def _commands_needing_modules(inputs):
    matrix, out = str(inputs / "matrix.csv"), str(inputs / "out")
    return {
        "dynamixs t1t2": ["dynamixs", "t1t2", "--input", matrix, "--out", out, "--exp", "T1"],
        "dynamixs methyl-t2": ["dynamixs", "methyl-t2", "--input", matrix, "--out", out],
        "dynamixs hetnoe": ["dynamixs", "hetnoe", "--sat", matrix, "--unsat", matrix,
                            "--out", out],
        "dynamixs t1rho": ["dynamixs", "t1rho", "--input", matrix, "--t1", matrix,
                           "--peaks", str(inputs / "peaks.txt"), "--omega1", "2000",
                           "--carrier", "118", "--theta", "90", "--out", out],
        "dynamixs density": ["dynamixs", "density", "--input", matrix, "--out", out,
                             "--field1-freq", "600"],
        "dynamixs modelfree": ["dynamixs", "modelfree", "--f1-t1", matrix, "--f1-t2", matrix,
                               "--f1-noe-sat", matrix, "--f1-noe-unsat", matrix,
                               "--field1-freq", "600", "--out", out],
        "kd": ["kd", "--input", matrix, "--out", out, "--p0", "50",
               "--conc-units", "absolute"],
        "export kd": ["export", "kd", "--json", str(inputs / "fit.json"), "--out", out],
        "integrate-1d": ["integrate-1d", "--spectra", str(inputs / "spectra_1d"),
                         "--out", out, "--detect"],
        "project build": ["project", "build", str(inputs / "run"),
                          "--out", str(inputs / "built.lunaNMR")],
    }


COMMANDS_NEEDING_MODULES = ["dynamixs t1t2", "dynamixs methyl-t2", "dynamixs hetnoe",
                            "dynamixs t1rho", "dynamixs density", "dynamixs modelfree",
                            "kd", "export kd", "integrate-1d", "project build"]


def test_the_copy_is_the_package_that_runs(package_only, tmp_path):
    done = _python(package_only, "-c", "import lunaNMR; print(lunaNMR.__file__)", cwd=tmp_path)
    assert done.returncode == 0, done.stderr[-2000:]
    assert Path(done.stdout.strip()) == package_only / "lunaNMR" / "__init__.py"


class TestACommandThatNeedsTheModules:

    @pytest.mark.parametrize("name", COMMANDS_NEEDING_MODULES)
    def test_it_stops_on_one_line_naming_the_editable_install(self, package_only, inputs,
                                                              tmp_path, name):
        done = _cli(package_only, *_commands_needing_modules(inputs)[name], cwd=tmp_path)
        assert done.returncode == 1, done.stdout[-1000:] + done.stderr[-2000:]
        assert "Traceback" not in done.stderr
        errors = [line for line in done.stderr.splitlines() if line.startswith("error: ")]
        assert len(errors) == 1, done.stderr[-2000:]
        assert REMEDY in errors[0]
        assert str(package_only / "modules") in errors[0]

    @pytest.mark.parametrize("name", COMMANDS_NEEDING_MODULES)
    def test_under_json_it_prints_the_error_envelope(self, package_only, inputs, tmp_path,
                                                     name):
        done = _cli(package_only, *_commands_needing_modules(inputs)[name],
                    "--format", "json", cwd=tmp_path)
        assert done.returncode == 1, done.stdout[-1000:] + done.stderr[-2000:]
        assert "Traceback" not in done.stderr
        envelope = json.loads(done.stdout)
        assert envelope["ok"] is False
        assert envelope["command"] == name
        assert envelope["error"]["type"] == ERROR_TYPE
        assert REMEDY in envelope["error"]["message"]
        assert "internal" not in envelope["error"]

    @pytest.mark.parametrize("name", COMMANDS_NEEDING_MODULES)
    def test_its_dry_run_stops_the_same_way(self, package_only, inputs, tmp_path, name):
        done = _cli(package_only, *_commands_needing_modules(inputs)[name], "--dry-run",
                    "--format", "json", cwd=tmp_path)
        assert done.returncode == 1, done.stdout[-1000:] + done.stderr[-2000:]
        assert json.loads(done.stdout)["error"]["type"] == ERROR_TYPE

    def test_a_failed_build_leaves_no_bundle(self, package_only, inputs, tmp_path):
        _cli(package_only, *_commands_needing_modules(inputs)["project build"], cwd=tmp_path)
        assert not (inputs / "built.lunaNMR").exists()


class TestWhatNeedsNoModulesStillRuns:

    @pytest.mark.parametrize("argv", [["--help"], ["--version"], ["series", "--help"],
                                      ["dynamixs", "t1t2", "--help"], ["kd", "--help"],
                                      ["integrate-1d", "--help"]],
                             ids=lambda argv: " ".join(argv))
    def test_help(self, package_only, tmp_path, argv):
        done = _cli(package_only, *argv, cwd=tmp_path)
        assert done.returncode == 0, done.stderr[-2000:]
        assert done.stdout.strip()

    def test_diagnose_reads_a_dataset(self, package_only, example_series, tmp_path):
        done = _cli(package_only, "diagnose", str(example_series), "--quick",
                    "--format", "json", cwd=tmp_path)
        assert "Traceback" not in done.stderr
        envelope = json.loads(done.stdout)
        assert envelope["command"] == "diagnose"
        assert envelope.get("error", {}).get("type") in (None, "PreflightFailed")
        assert envelope["experiments"], envelope

    def test_peaks_shift_writes_a_list(self, package_only, example_series, tmp_path):
        done = _cli(package_only, "peaks", "shift", "--peaks",
                    str(example_series / "r_600_peaks.txt"), "--out",
                    str(tmp_path / "shifted.txt"), "--dx", "0.01", "--format", "json",
                    cwd=tmp_path)
        assert done.returncode == 0, done.stdout[-1000:] + done.stderr[-2000:]
        assert json.loads(done.stdout)["n_peaks"] == 6
        assert (tmp_path / "shifted.txt").is_file()

    def test_batch_fits_a_spectrum(self, package_only, example_series, tmp_path):
        spectra = tmp_path / "spectra"
        spectra.mkdir()
        (spectra / "600_T1_0o0.ft").symlink_to(example_series / "600_T1_0o0.ft")
        done = _cli(package_only, "batch", str(spectra), "--output", str(tmp_path / "out"),
                    "--format", "json", cwd=tmp_path)
        assert done.returncode == 0, done.stdout[-1000:] + done.stderr[-2000:]
        envelope = json.loads(done.stdout)
        assert envelope["processed_files"] == 1
        assert envelope["total_peaks_fitted"] > 0

    def test_a_series_runs_and_is_built_into_a_bundle(self, package_only, example_series,
                                                      tmp_path):
        run = tmp_path / "run"
        series = _cli(package_only, "series", "--spectra", str(example_series),
                      "--peaks", str(example_series / "r_600_peaks.txt"),
                      "--out", str(run / "series"), "--format", "json", cwd=tmp_path)
        assert series.returncode == 0, series.stdout[-1000:] + series.stderr[-2000:]
        assert json.loads(series.stdout)["spectra_fitted"] == 2
        assert (run / "series" / "peak_intensity_matrix.csv").is_file()

        bundle = tmp_path / "built.lunaNMR"
        build = _cli(package_only, "project", "build", str(run), "--out", str(bundle),
                     "--format", "json", cwd=tmp_path)
        assert build.returncode == 0, build.stdout[-1000:] + build.stderr[-2000:]
        assert json.loads(build.stdout)["series"] == ["series"]

        inventory = _cli(package_only, "project", "inventory", str(bundle),
                         "--format", "json", cwd=tmp_path)
        assert inventory.returncode == 0, inventory.stderr[-2000:]
        assert json.loads(inventory.stdout)["ok"] is True


class TestTheInstallationCheck:

    def test_it_fails_on_the_missing_modules_and_names_the_editable_install(
            self, package_only, tmp_path):
        script = package_only / "lunaNMR" / "validation" / "verify_installation.py"
        done = _python(package_only, str(script), cwd=tmp_path)
        assert done.returncode == 1, done.stdout[-3000:]
        assert "INSTALLATION INCOMPLETE" in done.stdout
        assert "[FAIL] modules/dynamiXs_v2o0" in done.stdout
        assert "[FAIL] modules/integration_1d_v1o0" in done.stdout
        report = done.stdout.split("INSTALLATION REPORT")[1]
        assert REMEDY in report
        assert str(package_only / "modules") in report

    def test_it_passes_the_modules_of_this_checkout(self, tmp_path):
        script = REPO / "lunaNMR" / "validation" / "verify_installation.py"
        env = dict(os.environ, QT_QPA_PLATFORM="offscreen")
        done = subprocess.run([sys.executable, str(script)], cwd=str(tmp_path), env=env,
                              capture_output=True, text=True, timeout=600)
        assert done.returncode == 0, done.stdout[-3000:]
        assert "[ok] modules/dynamiXs_v2o0" in done.stdout
        assert "[ok] modules/integration_1d_v1o0" in done.stdout

    def test_the_mcp_check_fails_on_the_missing_modules(self, package_only, tmp_path):
        pytest.importorskip("mcp")
        done = _python(package_only, "-m", "lunaNMR.mcp.server", "--check", cwd=tmp_path)
        assert done.returncode == 1, done.stdout + done.stderr
        failed = [line for line in done.stdout.splitlines() if line.startswith("FAIL")]
        assert len(failed) == 1, done.stdout
        assert failed[0].startswith("FAIL modules: ")
        assert REMEDY in failed[0]
        assert "Traceback" not in done.stderr


class TestTheGui:
    """`lunaNMR.gui` imports the three module dialogs, so the window cannot open without
    the modules: the import is where the message is given."""

    def _import_error(self, site, cwd):
        done = _python(site, "-c", "import lunaNMR.gui", cwd=cwd)
        assert done.returncode == 1, done.stdout[-1000:]
        return done.stderr.strip().splitlines()[-1]

    def test_importing_it_names_the_editable_install(self, package_only, tmp_path):
        last = self._import_error(package_only, tmp_path)
        assert ERROR_TYPE in last
        assert REMEDY in last
        assert str(package_only / "modules" / "dynamiXs_v2o0") in last

    def test_a_modules_folder_without_1d_integration_is_named(self, tmp_path):
        site = _package_copy(tmp_path / "site")
        (site / "modules").mkdir()
        (site / "modules" / "dynamiXs_v2o0").symlink_to(REPO / "modules" / "dynamiXs_v2o0")
        last = self._import_error(site, tmp_path)
        assert ERROR_TYPE in last
        assert REMEDY in last
        assert str(site / "modules" / "integration_1d_v1o0") in last

    def test_the_lazy_main_export_carries_the_message(self, package_only, tmp_path):
        done = _python(package_only, "-c", "import lunaNMR; lunaNMR.main()", cwd=tmp_path)
        assert done.returncode == 1
        assert REMEDY in done.stderr.strip().splitlines()[-1]


class TestTheOneFunction:

    def test_it_finds_the_modules_of_this_checkout(self):
        from lunaNMR.analysis_modules import module_folder
        assert Path(module_folder("dynamiXs_v2o0", "dynamiXs_Kd")) == (
            REPO / "modules" / "dynamiXs_v2o0" / "dynamiXs_Kd")

    def test_a_folder_that_is_not_there_is_an_import_error_with_the_remedy(self):
        from lunaNMR.analysis_modules import AnalysisModulesNotFound, module_folder
        with pytest.raises(AnalysisModulesNotFound) as raised:
            module_folder("no_such_module")
        assert isinstance(raised.value, ImportError)
        assert REMEDY in str(raised.value)
        assert str(REPO / "modules" / "no_such_module") in str(raised.value)
        assert "\n" not in str(raised.value)

    def test_the_cli_puts_the_same_folder_on_the_import_path(self, monkeypatch):
        from lunaNMR.analysis_modules import module_folder
        from lunaNMR.cli import _add_modules_path
        monkeypatch.setattr(sys, "path", list(sys.path))
        monkeypatch.setenv("PYTHONPATH", "")
        folder = module_folder("integration_1d_v1o0")
        assert _add_modules_path("integration_1d_v1o0") == folder
        assert folder in sys.path
        assert folder in os.environ["PYTHONPATH"].split(os.pathsep)
