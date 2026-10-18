# ABOUTME: Shared pytest fixtures for the lunaNMR test suite: a short real DNAJA1 titration run
# ABOUTME: once per session through `lunaNMR series`, and a series rebuilt in memory from a tidy CSV.

import json
from pathlib import Path

import pandas as pd
import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
TITRATION_POINTS = ("0o0", "0o5", "3o0")


@pytest.fixture(scope="session")
def dnaja1_titration_series(tmp_path_factory):
    """Output folder of a 3-point DNAJA1 titration run with the CLI defaults; the run's
    JSON summary is written beside it as summary.json."""
    if not DNAJA1.exists():
        pytest.skip("data_test not present")
    import contextlib
    import io
    from lunaNMR.cli import main
    root = tmp_path_factory.mktemp("dnaja1_titration")
    (root / "spectra").mkdir()
    for point in TITRATION_POINTS:
        source = next(DNAJA1.glob(f"2D_*_{point}.ft"))
        (root / "spectra" / source.name).symlink_to(source)
    summary = io.StringIO()
    with contextlib.redirect_stdout(summary):
        code = main(["series", "--spectra", str(root / "spectra"), "--peaks", str(DNAJA1 / "DNAJA1.txt"),
                     "--mode", "titration", "--out", str(root / "out"), "--format", "json"])
    assert code == 0
    (root / "summary.json").write_text(summary.getvalue())
    return root / "out"


@pytest.fixture
def series_from_tidy():
    """Build a series as the GUI holds it (BatchResults) from a series folder's tidy CSV:
    one result per spectrum, carrying the per-peak fields the tidy CSV was written from.
    Returns (BatchResults, the tidy frame)."""
    from lunaNMR.processors.multi_spectrum_processor import BatchResults

    def build(tidy_csv):
        meta = json.loads((Path(tidy_csv).parent / "series_metadata.json").read_text())
        spectrum_of = {c['column']: c['spectrum'] for c in meta['columns']}
        tidy = pd.read_csv(tidy_csv, dtype={'spectrum_name': str})
        br = BatchResults()
        for column, rows in tidy.groupby('spectrum_name'):
            peaks = [{k: (None if pd.isna(v) else v) for k, v in row.items()
                      if k != 'spectrum_name'} for row in rows.to_dict('records')]
            br.add_result(spectrum_of[column], {'status': 'success', 'success': True,
                                                'integration_results': peaks})
        return br, tidy
    return build


@pytest.fixture
def open_overlay_viewer():
    """Open the overlay viewer on a series output folder, built as the series dialog
    builds it: each stored result with its spectrum's name and file."""
    import os
    from lunaNMR.gui.dialogs.multi_spectrum_viewer_dialog import MultiSpectrumViewerDialog

    def open_viewer(folder):
        stored = json.loads((Path(folder) / "batch_results.json").read_text())
        data_folder = stored['metadata']['data_folder']
        all_results = []
        for name, result in stored['results'].items():
            result = dict(result)
            result['spectrum_name'] = name
            result['spectrum_file'] = os.path.join(data_folder, name)
            all_results.append(result)
        return MultiSpectrumViewerDialog(all_results=all_results, output_folder=str(folder))
    return open_viewer
