# ABOUTME: Importing a DynamiXs fitter or plot script leaves the process's matplotlib backend alone,
# ABOUTME: and its plots draw from a worker thread under QtAgg without pyplot, as the GUI workers call them.

import os
import subprocess
import sys
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parent.parent
RUN = (Path(__file__).resolve().parents[4] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")

_PATHS = [str(_DIR), str(_DIR / "dynamiXs_T1_T2"), str(_DIR / "dynamiXs_density_functions"),
          str(_DIR / "dynamiXs_plot")]

_QT_PROCESS = f"""
import sys
sys.path[:0] = {_PATHS!r}
import matplotlib
matplotlib.use('QtAgg')
from PySide6.QtWidgets import QApplication
app = QApplication(sys.argv)
"""

FITTERS = ["fit_Tx_NMRRE", "fitmulti__Tx_NMRRE"]
SINGLE_FIELD = ["ZZ_multi_density", "ZZ_multi_density_087"]
DUAL_FIELD = ["ZZ_multi_2fields_density", "ZZ_multi_2fields_density_087"]
PLOT_SCRIPTS = ["compare_datasets", "compare_datasets_087", "plot_dual_field_datasets",
                "plot_dual_field_datasets_087"]


def _run(code, cwd):
    return subprocess.run([sys.executable, "-c", _QT_PROCESS + code], capture_output=True,
                          text=True, cwd=str(cwd), timeout=300,
                          env={**os.environ, "QT_QPA_PLATFORM": "offscreen"})


@pytest.mark.parametrize("module", FITTERS + SINGLE_FIELD + DUAL_FIELD + PLOT_SCRIPTS)
def test_importing_the_module_keeps_the_qt_backend(module, tmp_path):
    run = _run(f"import {module}\nprint(matplotlib.get_backend())", tmp_path)

    assert run.returncode == 0, run.stderr
    assert run.stdout.strip().lower() == "qtagg"


def _assert_drawn_off_the_main_thread(run, written):
    assert run.returncode == 0, run.stderr
    assert "outside of the main thread" not in run.stderr
    assert "pyplot figures: 0" in run.stdout
    assert written and all(p.stat().st_size > 0 for p in written)


_IN_A_WORKER = """
import threading
import matplotlib.pyplot as plt
worker = threading.Thread(target=draw)
worker.start()
worker.join()
print("pyplot figures:", len(plt.get_fignums()))
"""


@pytest.mark.skipif(not RUN.exists(), reason="data_test not present")
@pytest.mark.parametrize("module", FITTERS)
def test_fit_figures_draw_from_a_worker_thread(module, tmp_path):
    code = f"""
import numpy as np, pandas as pd
import {module} as fitter
m = pd.read_csv({str(RUN / 'series/600_T1/peak_intensity_matrix.csv')!r})
fit = pd.read_csv({str(RUN / 'fits/600_T1/600_T1_fit_results.txt')!r}, sep='\\t').set_index('Residue')
cols = [c for c in m.columns[4:] if '_' not in c]
x = np.array([float(c) for c in cols])
results = [{{'x': x, 'y': row[cols].to_numpy(float), 'residue': row['Assignment'],
             'A': fit.loc[row['Assignment'], 'A'], 't2': fit.loc[row['Assignment'], 'T1'],
             't2_err': fit.loc[row['Assignment'], 'T1_err'], 'C': fit.loc[row['Assignment'], 'C'],
             'success': True}} for _, row in m.head(25).iterrows()]
def draw():
    fitter.create_plots(results, 'T1', 20, 'T1', 'ms', 'Intensity')
""" + _IN_A_WORKER
    run = _run(code, tmp_path)

    _assert_drawn_off_the_main_thread(run, sorted(tmp_path.glob("T1_fit_results_fig*.pdf")))
    assert len(list(tmp_path.glob("T1_fit_results_fig*.pdf"))) == 2


def _density_frame(detailed_csv, jwh):
    """A model-free detailed results table, with the column names each script's plot reads."""
    return f"""
d = pd.read_csv({str(detailed_csv)!r})
if {jwh!r}:
    d = d.rename(columns=lambda c: c.replace('JwH_087', 'JwH'))
for f in ('1', '2'):
    for s in ('', '_err'):
        if f'Rex_f{{f}}{{s}}' in d:
            d[f'Rex_field{{f}}{{s}}'] = d[f'Rex_f{{f}}{{s}}']
"""


@pytest.mark.skipif(not RUN.exists(), reason="data_test not present")
@pytest.mark.parametrize("module", SINGLE_FIELD + DUAL_FIELD)
def test_density_figures_draw_from_a_worker_thread(module, tmp_path):
    dual = module in DUAL_FIELD
    source = RUN / ("modelfree/dual_raw/dual_raw_spectral_density_detailed.csv" if dual
                    else "modelfree/f600/f600_spectral_density_detailed.csv")
    cls = "DualFieldSpectralDensityAnalysis" if dual else "ReducedSpectralDensityAnalysis"
    method = "plot_dual_field_results" if dual else "plot_results"
    code = (f"import pandas as pd\nimport {module} as script\n"
            + _density_frame(source, not module.endswith("_087"))
            + f"def draw():\n    script.{cls}().{method}(d, save_plots=True)\n"
            + _IN_A_WORKER)
    run = _run(code, tmp_path)

    _assert_drawn_off_the_main_thread(run, list(tmp_path.glob("*.pdf")))
    assert len(list(tmp_path.glob("*.pdf"))) == 1
