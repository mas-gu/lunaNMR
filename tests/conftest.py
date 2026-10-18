# ABOUTME: Shared pytest fixtures for the lunaNMR test suite.
# ABOUTME: A short real DNAJA1 titration series, run once per session through `lunaNMR series`.

from pathlib import Path

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
