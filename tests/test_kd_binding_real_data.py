# ABOUTME: Kd titration contract tests against the real DNAJA1/ERDJ6 binding series.
# ABOUTME: Covers dummy/failed-point exclusion, dd_max runaway flagging, determinism and unit invariants.

import json
import math
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_BINDING = Path(__file__).resolve().parents[2] / "data_test" / "data_binding"

DNAJA1 = DATA_BINDING / "DNAJA1_HSPA8"
ERDJ6 = DATA_BINDING / "ERDJ6_HSPA8"

# 50 uM protein titrated to 0, 0.125, 0.25, 0.5, 1, 2, 3 equivalents.
P0 = "50"
CONC = "0,6.25,12.5,25,50,100,150"
CONC_X10 = "0,62.5,125,250,500,1000,1500"

# DNAJA1.txt carries 69 rows, one of which is the placeholder dummy_001.
DNAJA1_ROWS = 69
DNAJA1_REAL = 68
ERDJ6_REAL = 60

# D36 loses 5 of its 7 points, leaving too few to fit a three-parameter decay.
UNFITTABLE_RESIDUE = "D36"
# The intensity model has three free parameters (I0, I_inf, Kd).
MIN_POINTS = 4


def _partially_failed(data):
    """A residue that lost some points but kept enough to fit — the case that
    distinguishes point-wise exclusion from dropping the whole residue."""
    lost = data["metadata"].get("n_failed_points") or {}
    for res, n in sorted(lost.items()):
        if 0 < n <= 7 - MIN_POINTS:
            fit = _by_residue(data)[res].get("intensity") or {}
            if fit.get("L"):
                return res, n
    pytest.fail(f"no partially-failed residue available; n_failed_points={lost}")


def _require(dataset):
    if not dataset.is_dir() or not list(dataset.glob("*.ft")):
        pytest.skip(f"{dataset.name} binding data not available")


def _series(tmp_dir, dataset, peaks_name):
    """Run the real series pipeline once and return its tidy CSV."""
    from lunaNMR.cli import main
    out = tmp_dir / f"{dataset.name}_series"
    code = main(["series", "--spectra", str(dataset), "--peaks", str(dataset / peaks_name),
                 "--out", str(out), "--mode", "titration", "--peak-source", "cascade"])
    assert code == 0, f"series failed on {dataset.name}"
    tidy = out / "series_analysis_tidy.csv"
    # An exit code can be swallowed by a wrapper; a file on disk cannot.
    assert tidy.exists(), f"no tidy CSV produced for {dataset.name}"
    assert not pd.read_csv(tidy).empty
    return tidy


@pytest.fixture(scope="session")
def dnaja1_tidy(tmp_path_factory):
    _require(DNAJA1)
    return _series(tmp_path_factory.mktemp("dnaja1"), DNAJA1, "DNAJA1.txt")


@pytest.fixture
def own_tidy(dnaja1_tidy, tmp_path):
    """A private copy of the tidy CSV. The selection and params files default to
    siblings of the input, so tests that survey must not share an input directory —
    otherwise one test's hand-edit leaks into the next and results depend on order."""
    import shutil
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    local = data_dir / dnaja1_tidy.name
    shutil.copy(dnaja1_tidy, local)
    return local


@pytest.fixture(scope="session")
def erdj6_tidy(tmp_path_factory):
    _require(ERDJ6)
    return _series(tmp_path_factory.mktemp("erdj6"), ERDJ6, "ERDJ6.txt")


def _fit(tidy, out_dir, conc=CONC, p0=P0, extra=()):
    """Run `lunaNMR kd` and return (fit JSON dict, results.txt text)."""
    from lunaNMR.cli import main
    args = ["kd", "--input", str(tidy), "--out", str(out_dir), "--prefix", "t",
            "--p0", p0, "--conc", conc, "--observable", "csp,intensity", *extra]
    assert main(args) == 0
    # Machine artefacts live one level down; the human-readable report stays on top.
    data = json.loads((out_dir / "data" / "t_kd_fit_data.json").read_text())
    return data, (out_dir / "t_kd_results.txt").read_text()


def _by_residue(data):
    return {f["residue"]: f for f in data["fits"]}


class TestDummyResiduesAreExcluded:
    """dummy_* placeholders are excluded by every other fitter in the package
    (fit_Tx_NMRRE.py:225-232). The Kd path must match, and must report the count."""

    def test_only_real_residues_are_fitted(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        assert len(data["fits"]) == DNAJA1_REAL

    def test_no_dummy_residue_survives_into_the_fits(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        assert [r for r in _by_residue(data) if r.lower().startswith("dummy")] == []

    def test_no_dummy_row_in_the_results_table(self, dnaja1_tidy, tmp_path):
        _, results = _fit(dnaja1_tidy, tmp_path / "kd")
        assert "dummy" not in results.lower()

    def test_the_exclusion_count_is_reported_not_silent(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        assert data["metadata"].get("n_excluded_dummy") == 1

    def test_a_dataset_with_no_dummies_is_untouched(self, erdj6_tidy, tmp_path):
        """A filter that changes a dataset with nothing to filter matches too much."""
        data, _ = _fit(erdj6_tidy, tmp_path / "kd")
        assert len(data["fits"]) == ERDJ6_REAL
        assert data["metadata"].get("n_excluded_dummy", 0) == 0


class TestFailedPointsAreExcluded:
    """The series pipeline already labels unusable fits quality='Failed' and writes a
    0.1 height sentinel (once, a negative one). The Kd path consumes them unfiltered."""

    def test_the_series_really_does_label_them(self, dnaja1_tidy):
        """Guards the premise: if this stops holding, the tests below are meaningless."""
        df = pd.read_csv(dnaja1_tidy)
        failed = df[df["quality"].astype(str).str.lower() == "failed"]
        assert len(failed) > 0
        assert (failed["snr"] == 0).all()
        assert (df[df["height"] < 0]["quality"].astype(str).str.lower() == "failed").all()

    def test_failed_points_do_not_reach_the_fit(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        res, lost = _partially_failed(data)
        fit = _by_residue(data)[res]["intensity"]
        assert len(fit["L"]) == 7 - lost

    def test_no_negative_intensity_reaches_a_fit(self, dnaja1_tidy, tmp_path):
        """Says only that no negative survives into fitted data. It does NOT prove the
        exclusion is point-wise — a residue-wise drop would satisfy it too. Point-wise
        is proven by TestExclusionIsPointWiseNotResidueWise."""
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        for res, f in _by_residue(data).items():
            obs = (f.get("intensity") or {}).get("obs") or []
            assert all(v >= 0 for v in obs if v is not None), f"{res} fitted a negative intensity"

    def test_the_exclusion_count_is_reported_not_silent(self, dnaja1_tidy, tmp_path):
        """Reported per residue, so a reader can see which residues lost what."""
        lost = _fit(dnaja1_tidy, tmp_path / "kd")[0]["metadata"].get("n_failed_points")
        assert lost and sum(lost.values()) > 0

    def test_a_residue_left_with_too_few_points_says_so(self, dnaja1_tidy, tmp_path):
        """D36 keeps 2 of 7 points — fewer than the model's free parameters. Failing
        with a stated reason is correct; silently returning a fit would not be."""
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        fit = _by_residue(data)[UNFITTABLE_RESIDUE]["intensity"]
        assert fit["success"] is False
        assert fit.get("reason")


class TestTheGlobalPoolIsAuditable:
    """The shared Kd is pooled from a filtered subset, but only its size is recorded.
    Which residues produced the headline number cannot be recovered from the output."""

    @pytest.mark.parametrize("observable", ["csp", "intensity"])
    def test_the_pooled_residues_are_named(self, dnaja1_tidy, tmp_path, observable):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        block = data["global"][observable]
        assert "residues" in block
        assert len(block["residues"]) == block["n_residues"]
        assert set(block["residues"]) <= set(_by_residue(data))


class TestReferencePointIntegrity:
    """I/I0 divides by the zero-equivalent point. intensity_ratio_series rejects a
    reference of <= 0 but nothing else, so a failed reference silently rescales the
    whole residue: I48 has a reference height of 430.4 against a series median of
    3.4e+11 and its ratios reach 3.8e+07."""

    def test_the_series_really_does_contain_a_broken_reference(self, dnaja1_tidy):
        """Guards the premise for the assertions below."""
        df = pd.read_csv(dnaja1_tidy)
        ref = df[df["spectrum_name"] == 0.0]
        assert (ref["height"] < 1e3).any()

    def test_an_intensity_ratio_never_exceeds_one_by_orders_of_magnitude(
            self, dnaja1_tidy, tmp_path):
        """Binding attenuates intensity; a ratio of 1e+07 is a broken denominator."""
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        for res, f in _by_residue(data).items():
            obs = [v for v in ((f.get("intensity") or {}).get("obs") or []) if v is not None]
            if obs:
                assert max(obs) < 10, f"{res}: max I/I0 = {max(obs):.4g}"

    def test_a_failed_reference_point_is_not_used_as_a_denominator(
            self, dnaja1_tidy, tmp_path):
        """I48's reference is quality='Failed'; the residue cannot be normalised."""
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        fit = _by_residue(data).get("I48", {}).get("intensity") or {}
        assert not fit.get("success")

    def test_quality_labels_alone_do_not_catch_this(self, dnaja1_tidy):
        """T5 (height 3.0) and E43 (height 409.0) are labelled 'Good' at the reference
        point, so filtering on quality='Failed' is necessary but not sufficient."""
        df = pd.read_csv(dnaja1_tidy)
        ref = df[df["spectrum_name"] == 0.0]
        tiny_but_good = ref[(ref["height"] < 1e3)
                            & (ref["quality"].astype(str).str.lower() != "failed")]
        assert not tiny_but_good.empty


class TestNonPositiveIntensitiesAreRejected:
    """kd_models.peak_present already defines the convention (finite and > 0) so that
    CSP and intensity treat the same points as missing. The fit path must honour it at
    every point, not only the reference — excluded, never clamped to zero."""

    def test_a_negative_height_is_excluded_rather_than_clamped(self, dnaja1_tidy, tmp_path):
        """D36 at 0.125 eq has height -3.163e+09. Excluding says 'not measured';
        clamping to zero would claim the intensity was measured as zero."""
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        for res, f in _by_residue(data).items():
            fit = f.get("intensity") or {}
            for conc, val in zip(fit.get("L") or [], fit.get("obs") or []):
                assert val > 0, f"{res} fitted I/I0={val} at [L]={conc}"


class TestExclusionIsPointWiseNotResidueWise:
    """A residue with one unusable point keeps its remaining good points. Dropping the
    whole residue would silently discard measured data and nothing would report it."""

    def test_a_residue_with_bad_points_is_still_fitted(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        res, _lost = _partially_failed(data)
        assert res in _by_residue(data)

    def test_it_keeps_every_point_that_was_good(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        res, lost = _partially_failed(data)
        assert len(_by_residue(data)[res]["intensity"]["L"]) == 7 - lost

    def test_residues_with_one_bad_point_are_not_wholesale_dropped(self, dnaja1_tidy, tmp_path):
        """17 residues carry >=1 Failed point; only a residue failing at every point
        has nothing left to fit."""
        df = pd.read_csv(dnaja1_tidy)
        q = df["quality"].astype(str).str.lower()
        partial = {r for r in df[q == "failed"]["assignment"].unique()
                   if not (q[df["assignment"] == r] == "failed").all()
                   and not str(r).lower().startswith("dummy")}
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        assert partial - set(_by_residue(data)) == set()


def _selection_file(tmp_path, lines):
    path = tmp_path / "residues.txt"
    path.write_text("\n".join(lines) + "\n")
    return path


class TestResidueSelection:
    """--residues is how the human's step-4 choice reaches the fitter."""

    def test_only_the_listed_residues_are_fitted(self, dnaja1_tidy, tmp_path):
        sel = _selection_file(tmp_path, ["K3", "E4", "T5"])
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd", extra=["--residues", str(sel)])
        assert set(_by_residue(data)) == {"K3", "E4", "T5"}

    def test_an_inline_list_works_too(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd", extra=["--residues", "K3,E4"])
        assert set(_by_residue(data)) == {"K3", "E4"}

    def test_the_global_pool_is_drawn_only_from_the_selection(self, dnaja1_tidy, tmp_path):
        sel = _selection_file(tmp_path, ["K24", "K44", "K28", "E4", "E74", "G42"])
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd", extra=["--residues", str(sel)])
        pooled = set(data["global"]["csp"].get("residues") or [])
        assert pooled <= {"K24", "K44", "K28", "E4", "E74", "G42"}

    def test_a_commented_residue_is_left_out(self, dnaja1_tidy, tmp_path):
        sel = _selection_file(tmp_path, ["K3", "# E41 excluded: dd_max runaway", "E4"])
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd", extra=["--residues", str(sel)])
        assert set(_by_residue(data)) == {"K3", "E4"}

    def test_uncommenting_puts_it_back(self, dnaja1_tidy, tmp_path):
        """The round trip is the whole point: the file is the human's working document."""
        sel = tmp_path / "residues.txt"

        sel.write_text("K3\n# E4\nT5\n")
        without, _ = _fit(dnaja1_tidy, tmp_path / "kd_without", extra=["--residues", str(sel)])
        assert "E4" not in _by_residue(without)

        sel.write_text("K3\nE4\nT5\n")
        with_it, _ = _fit(dnaja1_tidy, tmp_path / "kd_with", extra=["--residues", str(sel)])
        assert "E4" in _by_residue(with_it)

    def test_an_unknown_residue_is_a_clean_error_naming_it(self, dnaja1_tidy, tmp_path, capsys):
        from lunaNMR.cli import main
        sel = _selection_file(tmp_path, ["K3", "Z999"])
        code = main(["kd", "--input", str(dnaja1_tidy), "--out", str(tmp_path / "kd"),
                     "--prefix", "t", "--p0", P0, "--conc", CONC,
                     "--residues", str(sel)])
        assert code != 0
        assert "Z999" in (capsys.readouterr().err or "")

    def test_an_empty_selection_is_an_error_not_a_success(self, dnaja1_tidy, tmp_path):
        from lunaNMR.cli import main
        sel = _selection_file(tmp_path, ["# K3", "# E4"])
        code = main(["kd", "--input", str(dnaja1_tidy), "--out", str(tmp_path / "kd"),
                     "--prefix", "t", "--p0", P0, "--conc", CONC,
                     "--residues", str(sel)])
        assert code != 0

    def test_without_the_flag_everything_real_is_still_fitted(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        assert len(data["fits"]) == DNAJA1_REAL


class TestSurvey:
    """--survey produces the plot that tells the human what not to fit. It must not
    produce the quotable number — that only comes from the selected-residue fit."""

    def _survey(self, tidy, out):
        from lunaNMR.cli import main
        code = main(["kd", "--survey", "--input", str(tidy), "--out", str(out),
                     "--prefix", "P", "--p0", P0, "--conc", CONC])
        assert code == 0
        return out

    def test_it_writes_the_figures_and_evidence_into_out(self, own_tidy, tmp_path):
        out = self._survey(own_tidy, tmp_path / "survey")
        for name in ("P_csp_vs_sequence.pdf", "P_intensity_vs_sequence.pdf"):
            assert (out / name).exists(), f"missing {name}"
        assert (out / "data" / "P_survey.csv").exists()

    def test_only_figures_and_the_report_sit_at_the_top_level(self, own_tidy, tmp_path):
        """The layout property itself, asserted directly rather than as twenty paths."""
        out = self._survey(own_tidy, tmp_path / "survey")
        stray = [f.name for f in out.iterdir()
                 if f.name != "data" and f.suffix not in (".pdf", ".txt")]
        assert stray == []

    def test_the_selection_lands_beside_the_input(self, own_tidy, tmp_path):
        """--out is a scratch path; the selection has to outlive it to be merged into."""
        out = self._survey(own_tidy, tmp_path / "survey")
        assert (own_tidy.parent / "P_residues.txt").exists()
        assert not (out / "P_residues.txt").exists()

    def test_it_writes_no_fit_json(self, own_tidy, tmp_path):
        out = self._survey(own_tidy, tmp_path / "survey")
        assert list(out.glob("**/*_kd_fit_data.json")) == []

    def test_it_surveys_every_row_in_the_input(self, own_tidy, tmp_path):
        """The survey is evidence about the input, so a placeholder appears here with a
        verdict rather than being hidden — the reader sees it was recognised."""
        out = self._survey(own_tidy, tmp_path / "survey")
        table = pd.read_csv(out / "data" / "P_survey.csv")
        assert len(table) == DNAJA1_ROWS
        dummy = table[table["residue"].str.lower().str.startswith("dummy")]
        assert len(dummy) == 1
        assert dummy.iloc[0]["verdict"] == "unusable"
        assert str(dummy.iloc[0]["reasons"]).strip()

    def test_a_selection_does_not_restrict_the_survey(self, own_tidy, tmp_path):
        """--survey is what produces a selection, so it cannot honour one."""
        from lunaNMR.cli import main
        sel = _selection_file(tmp_path, ["K3"])
        out = tmp_path / "survey"
        assert main(["kd", "--survey", "--input", str(own_tidy), "--out", str(out),
                     "--prefix", "P", "--p0", P0, "--conc", CONC,
                     "--residues", str(sel)]) == 0
        assert len(pd.read_csv(out / "data" / "P_survey.csv")) == DNAJA1_ROWS

    def test_the_suggestion_file_round_trips_into_a_fit(self, own_tidy, tmp_path):
        """The suggested selection must be directly usable as --residues input."""
        self._survey(own_tidy, tmp_path / "survey")
        data, _ = _fit(own_tidy, tmp_path / "kd",
                       extra=["--residues", str(own_tidy.parent / "P_residues.txt")])
        assert len(data["fits"]) > 0

    def test_a_suggested_exclusion_can_be_uncommented_by_hand(self, own_tidy, tmp_path):
        """The header tells the user to delete the '#'. Doing exactly that must yield a
        parseable name — not the name welded to its own reason string."""
        import sys
        sys.path.insert(0, str(REPO_ROOT / "modules" / "dynamiXs_v2o0" / "dynamiXs_Kd"))
        from kd_survey import parse_residues_file

        self._survey(own_tidy, tmp_path / "survey")
        path = own_tidy.parent / "P_residues.txt"
        before = parse_residues_file(path.read_text())

        excluded = [ln for ln in path.read_text().splitlines()
                    if ln.startswith("# ") and "excluded:" in ln
                    and not ln.lower().startswith("# dummy")]
        assert excluded, "survey suggested no exclusions to uncomment"
        name = excluded[0][2:].split()[0]

        path.write_text(path.read_text().replace(excluded[0], excluded[0][2:]))
        after = parse_residues_file(path.read_text())

        assert name in after, f"{name} did not come back after uncommenting"
        assert set(after) - set(before) == {name}
        for parsed in after:
            assert parsed.replace("_", "").isalnum(), f"garbage residue name {parsed!r}"

    def test_export_on_a_survey_output_fails_cleanly(self, own_tidy, tmp_path):
        from lunaNMR.cli import main
        out = self._survey(own_tidy, tmp_path / "survey")
        code = main(["export", "kd", "--json", str(out / "P_kd_fit_data.json"),
                     "--out", str(tmp_path / "figs")])
        assert code != 0


class TestDeterminism:
    """Prerequisite for every reproduction claim: the same input must give the same Kd."""

    def test_the_same_input_gives_the_same_kd(self, dnaja1_tidy, tmp_path):
        a, _ = _fit(dnaja1_tidy, tmp_path / "a")
        b, _ = _fit(dnaja1_tidy, tmp_path / "b")
        for obs in ("csp", "intensity"):
            assert a["global"][obs]["Kd"] == b["global"][obs]["Kd"]


class TestUnitInvariants:
    """Properties of the model, independent of how well the curve is sampled."""

    def test_scaling_the_titration_scales_the_csp_kd(self, dnaja1_tidy, tmp_path):
        """The CSP isotherm carries [P]0, so it is scale-invariant only when the
        protein concentration scales with the ligand.

        Measured with the Kd outlier gate OFF. That gate selects pool membership from
        fitted values, and an unconstrained per-residue fit does not reproduce under
        rescaling — so with it on, the two runs pool different residues and the ratio
        measures the membership change rather than the model."""
        off = ["--kd-outlier-z", "0"]
        base, _ = _fit(dnaja1_tidy, tmp_path / "base", extra=off)
        scaled, _ = _fit(dnaja1_tidy, tmp_path / "scaled", conc=CONC_X10, p0="500",
                         extra=off)
        assert (set(base["global"]["csp"]["residues"])
                == set(scaled["global"]["csp"]["residues"])), "pools differ; ratio is not the model"
        ratio = scaled["global"]["csp"]["Kd"] / base["global"]["csp"]["Kd"]
        assert ratio == pytest.approx(10.0, rel=1e-4)

    def test_the_csp_kd_is_not_scale_free_in_ligand_alone(self, dnaja1_tidy, tmp_path):
        """Scaling ligand without protein must NOT give 10x — [P]0 is a real term."""
        base, _ = _fit(dnaja1_tidy, tmp_path / "base")
        scaled, _ = _fit(dnaja1_tidy, tmp_path / "scaled", conc=CONC_X10)
        ratio = scaled["global"]["csp"]["Kd"] / base["global"]["csp"]["Kd"]
        assert ratio != pytest.approx(10.0, rel=1e-2)

    def test_scaling_the_titration_scales_the_intensity_kd(self, dnaja1_tidy, tmp_path):
        """The intensity decay is [P]0-independent, so ligand alone is enough."""
        base, _ = _fit(dnaja1_tidy, tmp_path / "base")
        scaled, _ = _fit(dnaja1_tidy, tmp_path / "scaled", conc=CONC_X10)
        ratio = scaled["global"]["intensity"]["Kd"] / base["global"]["intensity"]["Kd"]
        assert ratio == pytest.approx(10.0, rel=1e-4)

    def test_a_uniform_intensity_scale_cancels(self, dnaja1_tidy, tmp_path):
        """Scale factors correct for scan count; a uniform one cancels in I/I0."""
        base, _ = _fit(dnaja1_tidy, tmp_path / "base")
        scaled, _ = _fit(dnaja1_tidy, tmp_path / "scaled",
                         extra=["--intensity-scale", "2,2,2,2,2,2,2"])
        assert scaled["global"]["intensity"]["Kd"] == base["global"]["intensity"]["Kd"]

    def test_a_uniform_intensity_scale_leaves_csp_alone(self, dnaja1_tidy, tmp_path):
        """Scales apply to height/volume, never to positions."""
        base, _ = _fit(dnaja1_tidy, tmp_path / "base")
        scaled, _ = _fit(dnaja1_tidy, tmp_path / "scaled",
                         extra=["--intensity-scale", "2,2,2,2,2,2,2"])
        assert scaled["global"]["csp"]["Kd"] == base["global"]["csp"]["Kd"]

    def test_intensity_from_volume_is_a_working_alternative(self, dnaja1_tidy, tmp_path):
        base, _ = _fit(dnaja1_tidy, tmp_path / "base")
        vol, _ = _fit(dnaja1_tidy, tmp_path / "vol",
                      extra=["--intensity-from", "volume"])
        kd = vol["global"]["intensity"]["Kd"]
        assert math.isfinite(kd) and kd > 0
        assert vol["global"]["csp"]["Kd"] == base["global"]["csp"]["Kd"]


class TestTheSurveyToFitSeam:
    """Every piece of this workflow is tested in isolation: the survey writes a file, the
    file parses, a parsed selection restricts the pool. The seam is what those miss — that
    the file the survey ACTUALLY writes, edited the way a human ACTUALLY edits it, fits
    exactly the residues left uncommented. A '#'-handling bug lives here with every
    isolated piece still green."""

    def _survey(self, tidy, out):
        from lunaNMR.cli import main
        assert main(["kd", "--survey", "--input", str(tidy), "--out", str(out),
                     "--prefix", "P", "--p0", P0, "--conc", CONC]) == 0
        return tidy.parent / "P_residues.txt"      # sibling of the input, not --out

    @staticmethod
    def _parse(path):
        import sys
        sys.path.insert(0, str(REPO_ROOT / "modules" / "dynamiXs_v2o0" / "dynamiXs_Kd"))
        from kd_survey import parse_residues_file
        return parse_residues_file(path.read_text())

    @staticmethod
    def _edit_like_a_human(path):
        """Put one suggested exclusion back, and take one suggestion out."""
        lines = path.read_text().splitlines()
        restored = dropped = None
        for i, line in enumerate(lines):
            if (restored is None and line.startswith("# ") and "excluded:" in line
                    and not line.lower().startswith("# dummy")):
                lines[i] = line[2:]
                restored = line[2:].split()[0]
            elif dropped is None and line and not line.startswith("#"):
                lines[i] = "# " + line
                dropped = line.split()[0]
        assert restored and dropped, "survey file had nothing to edit both ways"
        path.write_text("\n".join(lines) + "\n")
        return restored, dropped

    def test_the_fit_covers_exactly_the_uncommented_residues(self, own_tidy, tmp_path):
        path = self._survey(own_tidy, tmp_path / "survey")
        self._edit_like_a_human(path)
        selected = set(self._parse(path))

        data, _ = _fit(own_tidy, tmp_path / "kd", extra=["--residues", str(path)])
        fitted = set(_by_residue(data))

        # Mechanical exclusions still win, so selecting a placeholder cannot resurrect it.
        assert fitted == {r for r in selected if not r.lower().startswith("dummy")}

    def test_a_restored_residue_is_actually_fitted(self, own_tidy, tmp_path):
        path = self._survey(own_tidy, tmp_path / "survey")
        restored, _dropped = self._edit_like_a_human(path)
        data, _ = _fit(own_tidy, tmp_path / "kd", extra=["--residues", str(path)])
        assert restored in _by_residue(data)

    def test_a_removed_residue_is_actually_absent(self, own_tidy, tmp_path):
        path = self._survey(own_tidy, tmp_path / "survey")
        _restored, dropped = self._edit_like_a_human(path)
        data, _ = _fit(own_tidy, tmp_path / "kd", extra=["--residues", str(path)])
        assert dropped not in _by_residue(data)

    def test_the_global_pool_never_exceeds_the_selection(self, own_tidy, tmp_path):
        """The quotable number must be traceable to residues the human chose."""
        path = self._survey(own_tidy, tmp_path / "survey")
        self._edit_like_a_human(path)
        selected = set(self._parse(path))
        data, _ = _fit(own_tidy, tmp_path / "kd", extra=["--residues", str(path)])
        for observable in ("csp", "intensity"):
            pooled = set(data["global"][observable]["residues"])
            assert pooled <= selected, f"{observable} pooled {pooled - selected}"

    def test_the_unedited_suggestion_fits_what_it_suggested(self, own_tidy, tmp_path):
        """Accepting the suggestion verbatim is the common path and must also hold."""
        path = self._survey(own_tidy, tmp_path / "survey")
        selected = set(self._parse(path))
        data, _ = _fit(own_tidy, tmp_path / "kd", extra=["--residues", str(path)])
        assert set(_by_residue(data)) == {r for r in selected
                                          if not r.lower().startswith("dummy")}


class TestSurveyMergesRatherThanClobbers:
    """The workflow's premise is that the human picks. A survey that rewrites the
    selection file wholesale erases those picks and re-suggests exclusions the human
    already rejected, so the loop cannot iterate."""

    def _survey(self, tidy, out):
        from lunaNMR.cli import main
        assert main(["kd", "--survey", "--input", str(tidy), "--out", str(out),
                     "--prefix", "P", "--p0", P0, "--conc", CONC]) == 0
        return tidy.parent / "P_residues.txt"      # sibling of the input, not --out

    @staticmethod
    def _parse(path):
        import sys
        sys.path.insert(0, str(REPO_ROOT / "modules" / "dynamiXs_v2o0" / "dynamiXs_Kd"))
        from kd_survey import parse_residues_file
        return parse_residues_file(path.read_text())

    @staticmethod
    def _edit(path):
        """Uncomment one suggested exclusion, comment out one suggested inclusion."""
        lines = path.read_text().splitlines()
        restored = dropped = None
        for i, line in enumerate(lines):
            if (restored is None and line.startswith("# ") and "excluded:" in line
                    and not line.lower().startswith("# dummy")):
                lines[i] = line[2:]
                restored = line[2:].split()[0]
            elif dropped is None and line and not line.startswith("#"):
                lines[i] = "# " + line
                dropped = line.split()[0]
        assert restored and dropped
        path.write_text("\n".join(lines) + "\n")
        return restored, dropped

    def test_both_human_edits_survive_a_second_survey(self, own_tidy, tmp_path):
        out = tmp_path / "survey"
        path = self._survey(own_tidy, out)
        restored, dropped = self._edit(path)

        self._survey(own_tidy, out)          # same directory, existing file present
        after = set(self._parse(path))

        assert restored in after, f"{restored} was re-excluded, overwriting the human"
        assert dropped not in after, f"{dropped} was re-included, overwriting the human"

    def test_the_whole_selection_is_preserved_not_just_the_edits(self, own_tidy, tmp_path):
        out = tmp_path / "survey"
        path = self._survey(own_tidy, out)
        self._edit(path)
        expected = set(self._parse(path))

        self._survey(own_tidy, out)
        assert set(self._parse(path)) == expected

    def test_the_merged_file_still_parses_to_clean_names(self, own_tidy, tmp_path):
        """Merge is a fresh chance to reintroduce the '#'-welding failure."""
        out = tmp_path / "survey"
        path = self._survey(own_tidy, out)
        self._edit(path)
        self._survey(own_tidy, out)
        for name in self._parse(path):
            assert name.replace("_", "").isalnum(), f"corrupt residue name {name!r}"

    def test_an_override_is_still_honoured_on_a_third_survey(self, own_tidy, tmp_path):
        """Durability across more than one cycle — a merge that only survives once is
        a merge that loses the decision on the run after."""
        out = tmp_path / "survey"
        path = self._survey(own_tidy, out)
        restored, dropped = self._edit(path)

        self._survey(own_tidy, out)
        self._survey(own_tidy, out)

        after = set(self._parse(path))
        assert restored in after
        assert dropped not in after

    def test_a_kept_override_is_marked_in_the_file(self, own_tidy, tmp_path):
        out = tmp_path / "survey"
        path = self._survey(own_tidy, out)
        restored, _dropped = self._edit(path)
        self._survey(own_tidy, out)

        line = next(ln for ln in path.read_text().splitlines()
                    if ln.strip().startswith(restored))
        assert "[kept]" in line

    def test_a_dropped_override_is_marked_in_the_file(self, own_tidy, tmp_path):
        out = tmp_path / "survey"
        path = self._survey(own_tidy, out)
        _restored, dropped = self._edit(path)
        self._survey(own_tidy, out)

        line = next(ln for ln in path.read_text().splitlines()
                    if ln.lstrip("# ").startswith(dropped))
        assert "[dropped]" in line

    def test_the_marker_does_not_break_parsing(self, own_tidy, tmp_path):
        """The marker is display-only; nothing downstream reads it back."""
        out = tmp_path / "survey"
        path = self._survey(own_tidy, out)
        restored, _dropped = self._edit(path)
        self._survey(own_tidy, out)
        assert restored in self._parse(path)

    def test_a_first_run_is_structurally_unchanged(self, own_tidy, tmp_path):
        """Merging must not alter behaviour when there is nothing to merge into.
        Structure, not bytes — reason wording is free to change."""
        a = self._parse(self._survey(own_tidy, tmp_path / "one"))
        b = self._parse(self._survey(own_tidy, tmp_path / "two"))
        assert a == b

    def test_selecting_an_unusable_residue_is_flagged_in_the_file(self, own_tidy, tmp_path):
        """Human selection is honoured in the file, data validity at the fit — but the
        file must say the residue will not fit, or the two disagree silently."""
        out = tmp_path / "survey"
        path = self._survey(own_tidy, out)
        lines = path.read_text().splitlines()
        target = None
        for i, line in enumerate(lines):
            if line.startswith("# ") and "reference point unusable" in line:
                lines[i] = line[2:]
                target = line[2:].split()[0]
                break
        assert target, "survey flagged no unusable reference to select"
        path.write_text("\n".join(lines) + "\n")

        self._survey(own_tidy, out)
        line = next(ln for ln in path.read_text().splitlines()
                    if ln.strip().startswith(target))
        assert "will not fit" in line

    def test_the_fit_names_selected_residues_it_could_not_fit(self, own_tidy, tmp_path):
        """The disagreement must be visible from the fit end too, not only discoverable
        by diffing a selection against a pool — and named, since a count tells you
        something went wrong without telling you what to reconcile."""
        out = tmp_path / "survey"
        path = self._survey(own_tidy, out)
        lines = path.read_text().splitlines()
        selected = None
        for i, line in enumerate(lines):
            if line.startswith("# ") and "reference point unusable" in line:
                lines[i] = line[2:]
                selected = line[2:].split()[0]
                break
        if selected is None:
            pytest.skip("survey flagged no unusable reference to select")
        path.write_text("\n".join(lines) + "\n")

        data, _ = _fit(own_tidy, tmp_path / "kd", extra=["--residues", str(path)])
        unfitted = data["metadata"]["unfitted"]
        assert selected in unfitted, f"{selected} was selected, not fitted, and not reported"
        assert unfitted[selected], "reported with no reason given"

    def test_a_half_fitted_residue_is_reported(self, own_tidy, tmp_path):
        """T5, E43 and L29 fit for CSP and fail for intensity. Without them the CSP and
        intensity pools differ in size with nothing to explain the difference."""
        data, _ = _fit(own_tidy, tmp_path / "kd")
        unfitted = data["metadata"]["unfitted"]
        for res, fit in _by_residue(data).items():
            failed = [obs for obs in ("csp", "intensity")
                      if not (fit.get(obs) or {}).get("success")]
            if failed and len(failed) < 2:
                assert res in unfitted, f"{res} failed only {failed} and was not reported"
                assert len(unfitted[res]) == 1, f"{res} listed {unfitted[res]}"

    def test_a_report_names_only_the_observables_that_failed(self, own_tidy, tmp_path):
        data, _ = _fit(own_tidy, tmp_path / "kd")
        for res, reasons in data["metadata"]["unfitted"].items():
            fit = _by_residue(data)[res]
            for obs in ("csp", "intensity"):
                named = any(r.startswith(f"{obs}:") for r in reasons)
                failed = not (fit.get(obs) or {}).get("success")
                assert named == failed, f"{res}: {obs} named={named} but failed={failed}"

    def test_nothing_is_reported_unfitted_when_everything_selected_fits(
            self, own_tidy, tmp_path):
        """The dict must stay empty on a clean selection, or it reads as noise."""
        sel = _selection_file(tmp_path, ["K24", "K44", "E4"])
        data, _ = _fit(own_tidy, tmp_path / "kd", extra=["--residues", str(sel)])
        assert data["metadata"].get("unfitted", {}) == {}


class TestTheSelectionFileLivesBesideTheData:
    """Merge only protects a human's decisions if both runs point at the same file, and
    --out is a scratch path that changes between runs. The selection belongs with the
    dataset, following the sibling convention kd_params already sets."""

    @staticmethod
    def _own_copy(tidy, tmp_path):
        """The default selection path is a sibling of the input, so each test needs its
        own input directory — writing beside the session fixture would leak between
        tests and make them order dependent."""
        import shutil
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        local = data_dir / tidy.name
        shutil.copy(tidy, local)
        return local

    def _survey(self, tidy, out, extra=()):
        from lunaNMR.cli import main
        assert main(["kd", "--survey", "--input", str(tidy), "--out", str(out),
                     "--prefix", "P", "--p0", P0, "--conc", CONC, *extra]) == 0

    def test_it_is_written_beside_the_input_not_into_out(self, dnaja1_tidy, tmp_path):
        tidy = self._own_copy(dnaja1_tidy, tmp_path)
        out = tmp_path / "scratch"
        self._survey(tidy, out)
        assert (tidy.parent / "P_residues.txt").exists()
        assert not (out / "P_residues.txt").exists()

    def test_a_second_survey_to_a_different_out_still_merges(self, dnaja1_tidy, tmp_path):
        """The whole reason for the move: two runs with different --out must share one
        selection, or the merge fix protects nothing in practice."""
        import sys
        sys.path.insert(0, str(REPO_ROOT / "modules" / "dynamiXs_v2o0" / "dynamiXs_Kd"))
        from kd_survey import parse_residues_file

        tidy = self._own_copy(dnaja1_tidy, tmp_path)
        self._survey(tidy, tmp_path / "run_one")
        path = tidy.parent / "P_residues.txt"

        lines = path.read_text().splitlines()
        dropped = next(ln.split()[0] for ln in lines if ln and not ln.startswith("#"))
        path.write_text("\n".join("# " + ln if ln.startswith(dropped) else ln
                                   for ln in lines) + "\n")

        self._survey(tidy, tmp_path / "run_two")
        assert dropped not in parse_residues_file(path.read_text())

    def test_selection_overrides_the_default_location(self, dnaja1_tidy, tmp_path):
        tidy = self._own_copy(dnaja1_tidy, tmp_path)
        chosen = tmp_path / "picked" / "mine.txt"
        chosen.parent.mkdir(parents=True)
        self._survey(tidy, tmp_path / "scratch", extra=["--selection", str(chosen)])
        assert chosen.exists()
        assert not (tidy.parent / "P_residues.txt").exists()

    def test_a_selection_is_never_picked_up_implicitly(self, dnaja1_tidy, tmp_path):
        """Silently applying a selection someone forgot about is worse than making them
        type a path, so a fit with no --residues must still fit everything."""
        tidy = self._own_copy(dnaja1_tidy, tmp_path)
        self._survey(tidy, tmp_path / "scratch")
        assert (tidy.parent / "P_residues.txt").exists()
        data, _ = _fit(tidy, tmp_path / "kd")
        assert len(data["fits"]) == DNAJA1_REAL


class TestThresholdsArePersistedAndOverridable:
    """A dataset should remember its own thresholds and a flag should still win.
    Asserted behaviourally — via what the survey actually flags — because a threshold
    that resolves correctly and is never applied is a stored preference, not a setting.
    (The survey does not yet WRITE a params JSON, so persistence is exercised with a
    hand-placed sibling, which is the documented discovery path.)"""

    @staticmethod
    def _own_copy(tidy, tmp_path):
        import shutil
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        local = data_dir / tidy.name
        shutil.copy(tidy, local)
        return local

    @staticmethod
    def _persist(tidy, **thresholds):
        (tidy.parent / "kd_params.json").write_text(json.dumps(thresholds))

    def _survey(self, tidy, out, extra=()):
        from lunaNMR.cli import main
        assert main(["kd", "--survey", "--input", str(tidy), "--out", str(out),
                     "--prefix", "P", "--p0", P0, "--conc", CONC, *extra]) == 0
        return out

    @staticmethod
    def _flagged(out, needle):
        reasons = pd.read_csv(out / "data" / "P_survey.csv")["reasons"].fillna("")
        return int(reasons.str.contains(needle).sum())

    def test_the_default_reference_threshold_is_applied(self, dnaja1_tidy, tmp_path):
        """At 10.0 the four broken references on DNAJA1 are caught."""
        tidy = self._own_copy(dnaja1_tidy, tmp_path)
        out = self._survey(tidy, tmp_path / "survey")
        assert self._flagged(out, "reference") == 4

    def test_a_flag_loosens_the_reference_threshold(self, dnaja1_tidy, tmp_path):
        """L29 sits at 49.87, so a threshold above it lets exactly L29 through."""
        tidy = self._own_copy(dnaja1_tidy, tmp_path)
        out = self._survey(tidy, tmp_path / "survey", extra=["--ref-max-ratio", "100"])
        assert self._flagged(out, "reference") == 3

    def test_a_persisted_value_beats_the_default(self, dnaja1_tidy, tmp_path):
        tidy = self._own_copy(dnaja1_tidy, tmp_path)
        self._persist(tidy, ref_max_ratio=100.0)
        out = self._survey(tidy, tmp_path / "survey")
        assert self._flagged(out, "reference") == 3

    def test_a_flag_beats_a_persisted_value(self, dnaja1_tidy, tmp_path):
        """The precedence chain's weak link: flag > persisted > default."""
        tidy = self._own_copy(dnaja1_tidy, tmp_path)
        self._persist(tidy, ref_max_ratio=100.0)
        out = self._survey(tidy, tmp_path / "survey", extra=["--ref-max-ratio", "10"])
        assert self._flagged(out, "reference") == 4

    def test_the_plateau_threshold_is_applied(self, dnaja1_tidy, tmp_path):
        """12 of 69 residues exceed a dd_max ratio of 10 on DNAJA1; raising it flags fewer."""
        tidy = self._own_copy(dnaja1_tidy, tmp_path)
        tight = self._survey(tidy, tmp_path / "tight", extra=["--dd-runaway-ratio", "10"])
        loose = self._survey(tidy, tmp_path / "loose", extra=["--dd-runaway-ratio", "1000"])
        assert self._flagged(tight, "plateau") > self._flagged(loose, "plateau")

    def test_the_noise_quantile_is_applied(self, dnaja1_tidy, tmp_path):
        """The floor is a quantile of the dataset's own CSP spread, so raising it must
        call more residues non-movers."""
        tidy = self._own_copy(dnaja1_tidy, tmp_path)
        low = self._survey(tidy, tmp_path / "low", extra=["--noise-quantile", "0.10"])
        high = self._survey(tidy, tmp_path / "high", extra=["--noise-quantile", "0.50"])
        assert self._flagged(high, "noise floor") > self._flagged(low, "noise floor")

    def test_the_fit_uses_the_same_reference_threshold_as_the_survey(
            self, dnaja1_tidy, tmp_path):
        """A dataset judged at one threshold and fitted at another is the asymmetry that
        would let the survey clear a residue the fit then silently drops."""
        loose, _ = _fit(dnaja1_tidy, tmp_path / "loose",
                        extra=["--ref-max-ratio", "1e12"])
        tight, _ = _fit(dnaja1_tidy, tmp_path / "tight",
                        extra=["--ref-max-ratio", "10"])
        assert len(tight["metadata"]["unfitted"]) > len(loose["metadata"]["unfitted"])


# The spectra are named in equivalents; [L] = eq * P0.
EQUIVALENTS = "0,0.125,0.25,0.5,1,2,3"


class TestConcentrationUnits:
    """Without --conc, equivalents parsed from spectrum names were fed in as molar
    concentrations against P0=50 — a confident Kd wrong by the factor P0, success=True,
    no warning. --conc-units names what the numbers are."""

    def _fit_raw(self, tidy, out, extra):
        from lunaNMR.cli import main
        assert main(["kd", "--input", str(tidy), "--out", str(out), "--prefix", "t",
                     "--p0", P0, "--observable", "csp,intensity", *extra]) == 0
        return json.loads((out / "data" / "t_kd_fit_data.json").read_text())

    def test_equivalents_without_conc_matches_the_explicit_absolute_list(
            self, dnaja1_tidy, tmp_path):
        """The acceptance test, designed before the feature existed."""
        eq = self._fit_raw(dnaja1_tidy, tmp_path / "eq", ["--conc-units", "equivalents"])
        ab = self._fit_raw(dnaja1_tidy, tmp_path / "ab", ["--conc", CONC])
        assert eq["metadata"]["concentrations"] == ab["metadata"]["concentrations"]
        for observable in ("csp", "intensity"):
            assert (eq["global"][observable]["Kd"]
                    == pytest.approx(ab["global"][observable]["Kd"], rel=1e-9))

    def test_the_per_residue_fits_match_too(self, dnaja1_tidy, tmp_path):
        """A matching global Kd could hide compensating per-residue differences."""
        eq = self._fit_raw(dnaja1_tidy, tmp_path / "eq", ["--conc-units", "equivalents"])
        ab = self._fit_raw(dnaja1_tidy, tmp_path / "ab", ["--conc", CONC])
        assert _by_residue(eq).keys() == _by_residue(ab).keys()
        for name, fit in _by_residue(eq).items():
            other = _by_residue(ab)[name]
            for observable in ("csp", "intensity"):
                a = (fit.get(observable) or {}).get("Kd")
                b = (other.get(observable) or {}).get("Kd")
                if a is None or b is None:
                    assert a is b, f"{name}/{observable} succeeded in only one run"
                else:
                    assert a == pytest.approx(b, rel=1e-9), f"{name}/{observable}"

    def test_absolute_is_the_default_and_is_unchanged(self, dnaja1_tidy, tmp_path):
        """Getting this wrong moves the headline number by the factor P0, so the
        default must be provably untouched by the feature."""
        plain = self._fit_raw(dnaja1_tidy, tmp_path / "plain", ["--conc", CONC])
        named = self._fit_raw(dnaja1_tidy, tmp_path / "named",
                              ["--conc", CONC, "--conc-units", "absolute"])
        assert plain["metadata"]["concentrations"] == named["metadata"]["concentrations"]
        for observable in ("csp", "intensity"):
            assert (plain["global"][observable]["Kd"]
                    == named["global"][observable]["Kd"])

    def test_equivalents_scales_an_explicit_conc_list(self, dnaja1_tidy, tmp_path):
        """The flag describes the numbers, not where they came from."""
        given = self._fit_raw(dnaja1_tidy, tmp_path / "given",
                              ["--conc", EQUIVALENTS, "--conc-units", "equivalents"])
        assert given["metadata"]["concentrations"] == pytest.approx(
            [0.0, 6.25, 12.5, 25.0, 50.0, 100.0, 150.0])

    def test_the_silent_wrong_answer_is_still_reachable_by_default(
            self, dnaja1_tidy, tmp_path):
        """Documents the hazard the flag exists for: absolute is the default, so a
        titration named in equivalents still needs one of --conc or --conc-units."""
        wrong = self._fit_raw(dnaja1_tidy, tmp_path / "wrong", [])
        right = self._fit_raw(dnaja1_tidy, tmp_path / "right",
                              ["--conc-units", "equivalents"])
        assert wrong["metadata"]["concentrations"] != right["metadata"]["concentrations"]

    def test_a_persisted_conc_units_is_applied(self, own_tidy, tmp_path):
        """A stored setting nobody reads is a preference, not a setting."""
        (own_tidy.parent / "kd_params.json").write_text(
            json.dumps({"conc_units": "equivalents"}))
        persisted = self._fit_raw(own_tidy, tmp_path / "persisted", [])
        assert persisted["metadata"]["concentrations"] == pytest.approx(
            [0.0, 6.25, 12.5, 25.0, 50.0, 100.0, 150.0])

    def test_a_flag_beats_a_persisted_conc_units(self, own_tidy, tmp_path):
        (own_tidy.parent / "kd_params.json").write_text(
            json.dumps({"conc_units": "equivalents"}))
        forced = self._fit_raw(own_tidy, tmp_path / "forced",
                               ["--conc-units", "absolute"])
        assert forced["metadata"]["concentrations"] == pytest.approx(
            [0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 3.0])


class TestTheSurveyRecordsWhatTheFitNeeds:
    """AFFINITY_PLAYBOOK tells an agent the survey persists its settings and not to
    re-enter them on the fit. Whatever the survey fails to record, the fit re-derives in
    silence: concentrations fall back to the CSV's own point labels, and conc_units falls
    back to 'absolute', which reads equivalents as molar — a Kd wrong by the factor P0
    with success=True. This is the round-trip that sentence promises.
    """

    def _survey(self, tidy, out, extra):
        from lunaNMR.cli import main
        assert main(["kd", "--survey", "--input", str(tidy), "--out", str(out),
                     "--prefix", "t", "--p0", P0, *extra]) == 0

    def _fit(self, tidy, out, extra):
        from lunaNMR.cli import main
        assert main(["kd", "--input", str(tidy), "--out", str(out), "--prefix", "t",
                     "--p0", P0, "--observable", "csp,intensity", *extra]) == 0
        return json.loads((out / "data" / "t_kd_fit_data.json").read_text())

    def test_a_fit_after_a_survey_inherits_its_units(self, own_tidy, tmp_path):
        """The documented workflow: survey names the units, the fit does not repeat them."""
        self._survey(own_tidy, tmp_path / "s", ["--conc-units", "equivalents"])
        inherited = self._fit(own_tidy, tmp_path / "f", [])
        assert inherited["metadata"]["concentrations"] == pytest.approx(
            [0.0, 6.25, 12.5, 25.0, 50.0, 100.0, 150.0])

    def test_a_fit_after_a_survey_inherits_its_concentrations(self, own_tidy, tmp_path):
        """An explicit --conc list is the other half of the same promise."""
        self._survey(own_tidy, tmp_path / "s", ["--conc", CONC_X10])
        inherited = self._fit(own_tidy, tmp_path / "f", [])
        assert inherited["metadata"]["concentrations"] == pytest.approx(
            [0.0, 62.5, 125.0, 250.0, 500.0, 1000.0, 1500.0])

    def test_re_stating_the_units_does_not_convert_twice(self, own_tidy, tmp_path):
        """The trap in persisting concentrations: if the stored list were already
        converted, an agent that also passes --conc-units would multiply by P0 again.
        Restating what the survey recorded must be a no-op, not an error."""
        self._survey(own_tidy, tmp_path / "s",
                     ["--conc", CONC_X10, "--conc-units", "equivalents"])
        restated = self._fit(own_tidy, tmp_path / "f", ["--conc-units", "equivalents"])
        silent = self._fit(own_tidy, tmp_path / "g", [])
        assert restated["metadata"]["concentrations"] == pytest.approx(
            silent["metadata"]["concentrations"])

    def test_the_survey_records_the_thresholds_it_judged_with(self, own_tidy, tmp_path):
        """A threshold the survey used to reject residues, but does not record, makes
        the fit's pool disagree with the selection file the survey just wrote."""
        from kd_params import find_params_source, load_params
        self._survey(own_tidy, tmp_path / "s", ["--csp-sigma-multiple", "2.0"])
        stored = load_params(find_params_source(str(own_tidy)))
        assert stored["csp_sigma_multiple"] == 2.0


# Measured at the last titration point, top 10% trimmed so the binders do not set the
# threshold meant to identify them. Untrimmed sigma is materially larger.
DNAJA1_SIGMA_TRIMMED = 0.0202
DNAJA1_SIGMA_UNTRIMMED = 0.0288
DNAJA1_SIGNIFICANT = 25
DNAJA1_MEASURED = 53


class TestCspSignificanceCutoff:
    """Only residues whose CSP exceeds sigma at the highest titrant point contribute to
    the shared CSP Kd. A residue with no CSP at that point is unmeasured, which is a
    third state — calling it 'not significant' would assert something never observed."""

    def _sig(self, data):
        return data["metadata"]["csp_significance"]

    def test_sigma_is_the_trimmed_population(self, dnaja1_tidy, tmp_path):
        """The specific way decision 2 reverts quietly: using sigma over all residues
        lets the binders inflate the threshold meant to find them."""
        block = self._sig(_fit(dnaja1_tidy, tmp_path / "kd")[0])
        assert block["sigma"] == pytest.approx(DNAJA1_SIGMA_TRIMMED, rel=0.05)
        assert block["sigma"] != pytest.approx(DNAJA1_SIGMA_UNTRIMMED, rel=0.05)

    def test_the_trim_fraction_is_recorded(self, dnaja1_tidy, tmp_path):
        assert self._sig(_fit(dnaja1_tidy, tmp_path / "kd")[0])["trim_fraction"] \
            == pytest.approx(0.1)

    def test_the_threshold_is_the_multiple_times_sigma(self, dnaja1_tidy, tmp_path):
        block = self._sig(_fit(dnaja1_tidy, tmp_path / "kd")[0])
        assert block["threshold"] == pytest.approx(block["multiple"] * block["sigma"])

    def test_the_expected_number_of_residues_is_significant(self, dnaja1_tidy, tmp_path):
        block = self._sig(_fit(dnaja1_tidy, tmp_path / "kd")[0])
        assert len(block["significant"]) == DNAJA1_SIGNIFICANT

    def test_the_three_states_are_disjoint(self, dnaja1_tidy, tmp_path):
        block = self._sig(_fit(dnaja1_tidy, tmp_path / "kd")[0])
        sig, notsig, unmeasured = (set(block["significant"]),
                                   set(block["not_significant"]),
                                   set(block["unmeasured"]))
        assert sig & notsig == set()
        assert sig & unmeasured == set()
        assert notsig & unmeasured == set()

    def test_the_three_states_cover_every_residue(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        block = self._sig(data)
        covered = (set(block["significant"]) | set(block["not_significant"])
                   | set(block["unmeasured"]))
        assert covered == set(_by_residue(data))

    def test_unmeasured_is_not_reported_as_insignificant(self, dnaja1_tidy, tmp_path):
        """A residue with no CSP at the last point was never observed there; calling it
        insignificant asserts a measurement nobody made."""
        block = self._sig(_fit(dnaja1_tidy, tmp_path / "kd")[0])
        assert block["unmeasured"], "no unmeasured residues — the distinction is untested"
        assert block["n_measured"] == (len(block["significant"])
                                       + len(block["not_significant"]))

    def test_the_per_residue_flag_agrees_with_the_lists(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        block = self._sig(data)
        expected = {name: True for name in block["significant"]}
        expected.update({name: False for name in block["not_significant"]})
        expected.update({name: None for name in block["unmeasured"]})
        for name, fit in _by_residue(data).items():
            assert (fit.get("csp") or {}).get("significant") is expected[name], name

    def test_the_global_csp_pool_is_drawn_only_from_significant_residues(
            self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        pooled = set(data["global"]["csp"]["residues"])
        assert pooled <= set(self._sig(data)["significant"])

    def test_insignificant_residues_keep_their_own_kd(self, dnaja1_tidy, tmp_path):
        """Gating the pool must not suppress per-residue reporting — the decision was
        to gate the shared Kd only."""
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        by_residue = _by_residue(data)
        reported = [r for r in self._sig(data)["not_significant"]
                    if (by_residue[r].get("csp") or {}).get("success")]
        assert reported, "every insignificant residue lost its per-residue CSP fit"

    def test_intensity_is_untouched_by_the_csp_cutoff(self, dnaja1_tidy, tmp_path):
        """CSP only. If the intensity pool were also gated it would be a subset."""
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        intensity_pool = set(data["global"]["intensity"]["residues"])
        assert intensity_pool - set(self._sig(data)["significant"]), \
            "the intensity pool carries no CSP-insignificant residue — CSP gating leaked"

    def test_a_larger_multiple_admits_fewer_residues(self, dnaja1_tidy, tmp_path):
        """2 sigma is a common variant, so the knob has to actually move the cutoff."""
        one = self._sig(_fit(dnaja1_tidy, tmp_path / "one")[0])
        two = self._sig(_fit(dnaja1_tidy, tmp_path / "two",
                             extra=["--csp-sigma-multiple", "2.0"])[0])
        assert two["multiple"] == pytest.approx(2.0)
        assert two["threshold"] == pytest.approx(2.0 * one["sigma"])
        assert len(two["significant"]) < len(one["significant"])

    def test_a_persisted_multiple_is_applied(self, own_tidy, tmp_path):
        (own_tidy.parent / "kd_params.json").write_text(
            json.dumps({"csp_sigma_multiple": 2.0}))
        assert self._sig(_fit(own_tidy, tmp_path / "kd")[0])["multiple"] \
            == pytest.approx(2.0)

    def test_a_flag_beats_a_persisted_multiple(self, own_tidy, tmp_path):
        (own_tidy.parent / "kd_params.json").write_text(
            json.dumps({"csp_sigma_multiple": 2.0}))
        block = self._sig(_fit(own_tidy, tmp_path / "kd",
                               extra=["--csp-sigma-multiple", "1.0"])[0])
        assert block["multiple"] == pytest.approx(1.0)

    def test_the_expected_number_of_residues_is_measured(self, dnaja1_tidy, tmp_path):
        block = self._sig(_fit(dnaja1_tidy, tmp_path / "kd")[0])
        assert block["n_measured"] == DNAJA1_MEASURED

    def test_the_second_dataset_is_judged_on_its_own_spread(self, erdj6_tidy, tmp_path):
        """Sigma is a property of the dataset, so ERDJ6 must not inherit DNAJA1's."""
        block = self._sig(_fit(erdj6_tidy, tmp_path / "kd")[0])
        assert block["sigma"] != pytest.approx(DNAJA1_SIGMA_TRIMMED, rel=0.02)


class TestTheCspPoolReconcilesInOnePlace:
    """A residue can now be absent from the shared CSP Kd for four different reasons.
    Without one map naming which gate removed each, reconciling a small pool means
    reconstructing it from four mechanisms after the fact."""

    def test_every_fitted_residue_is_either_pooled_or_has_a_stated_reason(
            self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        pooled = set(data["global"]["csp"]["residues"])
        excluded = data["metadata"]["csp_pool_excluded"]
        assert pooled | set(excluded) == set(_by_residue(data))

    def test_no_residue_is_both_pooled_and_excluded(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        assert set(data["global"]["csp"]["residues"]) & set(
            data["metadata"]["csp_pool_excluded"]) == set()

    def test_every_exclusion_carries_a_reason(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        for residue, why in data["metadata"]["csp_pool_excluded"].items():
            assert isinstance(why, str) and why.strip(), residue

    def test_the_reasons_distinguish_the_gates(self, dnaja1_tidy, tmp_path):
        """One reason string for everything would reconcile the count while hiding
        which mechanism did the removing."""
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        assert len(set(data["metadata"]["csp_pool_excluded"].values())) > 1


# A real fit JSON produced before csp_significance and global['csp']['residues'] existed.
LEGACY_KD_JSON = REPO_ROOT / "data_example" / "tmp" / "kd" / "kd_fit" / "kd_kd_fit_data.json"


class TestOlderFitJsonsStillExport:
    """Two independent places an older JSON hits the new figure code: the significance
    line reads metadata['csp_significance'], and the kd-bars filter reads
    global[obs]['residues']. Neither existed before today, and a saved analysis must
    still render rather than raise."""

    @pytest.fixture
    def legacy(self):
        if not LEGACY_KD_JSON.exists():
            pytest.skip("no pre-significance fit JSON available")
        return json.loads(LEGACY_KD_JSON.read_text())

    def test_the_fixture_really_predates_both_keys(self, legacy):
        """Guards the premise: if a regenerated fixture gains these keys, the tests
        below silently stop testing degradation."""
        assert "csp_significance" not in legacy.get("metadata", {})
        assert not any("residues" in (legacy.get("global", {}).get(obs) or {})
                       for obs in ("csp", "intensity"))

    def test_it_exports_without_raising(self, tmp_path):
        from lunaNMR.cli import main
        if not LEGACY_KD_JSON.exists():
            pytest.skip("no pre-significance fit JSON available")
        out = tmp_path / "legacy_figs"
        assert main(["export", "kd", "--json", str(LEGACY_KD_JSON),
                     "--out", str(out), "--prefix", "old"]) == 0
        assert list(out.glob("*.pdf")), "no figures produced from an older fit"

    def test_the_significance_line_is_simply_absent(self, tmp_path):
        """No csp_significance means no threshold to draw — the page renders without it."""
        from lunaNMR.cli import main
        if not LEGACY_KD_JSON.exists():
            pytest.skip("no pre-significance fit JSON available")
        out = tmp_path / "legacy_refbars"
        assert main(["export", "kd", "--json", str(LEGACY_KD_JSON), "--out", str(out),
                     "--prefix", "old", "--kind", "ref-bars"]) == 0
        assert (out / "old_csp_ref_vs_point.pdf").exists()

    def test_the_kd_bars_fall_back_to_every_residue(self, tmp_path):
        """With no recorded pool, filtering to it would draw nothing; an older fit
        reported every residue's Kd and the figure must still say so."""
        from lunaNMR.cli import main
        if not LEGACY_KD_JSON.exists():
            pytest.skip("no pre-significance fit JSON available")
        out = tmp_path / "legacy_kdbars"
        assert main(["export", "kd", "--json", str(LEGACY_KD_JSON), "--out", str(out),
                     "--prefix", "old", "--kind", "kd-bars"]) == 0
        pdf = out / "old_csp_kd_vs_residue.pdf"
        assert pdf.exists() and pdf.stat().st_size > 0


def _seq_key(name):
    """Residue number, so K3 precedes K10. A plain string sort puts K10 first."""
    match = re.match(r"([A-Za-z]+)(\d+)", str(name))
    return (int(match.group(2)), match.group(1)) if match else (10**9, str(name))


def _in_sequence_order(names):
    names = list(names)
    return names == sorted(names, key=_seq_key)


class TestResiduesAreOrderedBySequenceNumber:
    """Every sorted() over residue names is a string sort unless it says otherwise, so
    K10 lands before K3 and the sequence silently scrambles at two digits. The survey
    already used residue_sort_key; the fit did not, so the two disagreed."""

    def test_the_fits_are_in_sequence_order(self, dnaja1_tidy, tmp_path):
        """This one drives everything else — results.txt, every figure's x-axis, and
        the order of csp_pool_excluded all follow the fits order."""
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        assert _in_sequence_order(f["residue"] for f in data["fits"])

    @pytest.mark.parametrize("key", ["significant", "not_significant", "unmeasured"])
    def test_each_significance_list_is_in_sequence_order(self, dnaja1_tidy, tmp_path, key):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        assert _in_sequence_order(data["metadata"]["csp_significance"][key])

    @pytest.mark.parametrize("observable", ["csp", "intensity"])
    def test_the_global_pool_is_in_sequence_order(self, dnaja1_tidy, tmp_path, observable):
        """Derived from fits order, so this should follow — asserted rather than assumed."""
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        assert _in_sequence_order(data["global"][observable]["residues"])

    def test_the_exclusion_map_is_in_sequence_order(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        assert _in_sequence_order(data["metadata"]["csp_pool_excluded"])

    def test_the_results_table_is_in_sequence_order(self, dnaja1_tidy, tmp_path):
        _data, results = _fit(dnaja1_tidy, tmp_path / "kd")
        names = [ln.split("\t")[0] for ln in results.splitlines()
                 if ln and not ln.startswith("#") and not ln.startswith("Residue")]
        assert _in_sequence_order(names)

    def test_the_exported_tables_follow_the_json(self, dnaja1_tidy, tmp_path):
        """The figures' x-axis comes from the same order as their CSV mirror, so the
        CSV is the machine-readable check on the figure."""
        from lunaNMR.cli import main
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        figs = tmp_path / "figs"
        assert main(["export", "kd", "--json",
                     str(tmp_path / "kd" / "data" / "t_kd_fit_data.json"),
                     "--out", str(figs), "--prefix", "t"]) == 0
        table = pd.read_csv(figs / "data" / "t_csp_ref_vs_point.csv")
        column = table.columns[0]
        assert _in_sequence_order(table[column])


class TestOrderingSurvivesADecadeBoundary:
    """A dataset whose residues are all single-digit passes an alphabetical sort by
    luck, so the property has to be pinned where a string sort actually differs."""

    @staticmethod
    def _tidy(tmp_path, residues=("K3", "K10", "K11", "K2")):
        rows = ["spectrum_name,assignment,peak_number,ppm_x,ppm_y,height,volume,snr,quality,r_squared"]
        for point_index, point in enumerate((0.0, 6.25, 12.5, 25.0, 50.0, 100.0, 150.0)):
            for peak, name in enumerate(residues):
                shift = 8.0 + 0.01 * peak + 0.002 * point_index
                height = 1e9 * (1.0 - 0.1 * point_index)
                rows.append(f"{point},{name},{peak},{shift:.4f},120.0,{height},{height/10},50,Excellent,0.99")
        path = tmp_path / "decade_tidy.csv"
        path.write_text("\n".join(rows) + "\n")
        return path

    def test_a_two_digit_residue_does_not_sort_before_a_one_digit_one(self, tmp_path):
        data, _ = _fit(self._tidy(tmp_path), tmp_path / "kd")
        names = [f["residue"] for f in data["fits"]]
        assert names == ["K2", "K3", "K10", "K11"]

    def test_a_string_sort_would_have_failed_this(self):
        """Guards the fixture: if these ever sort the same way, the test proves nothing."""
        names = ["K3", "K10", "K11", "K2"]
        assert sorted(names) != sorted(names, key=_seq_key)


class TestKdOutlierExclusion:
    """Kd spans 25 decades on this data, so mean/stdev cannot find outliers: the largest
    value is what sets sigma, and it then sits inside its own 2-sigma band. The gate uses
    median + MAD on log10(Kd)."""

    def _block(self, data):
        return data["metadata"]["csp_kd_outliers"]

    def test_the_gate_records_its_statistics(self, dnaja1_tidy, tmp_path):
        block = self._block(_fit(dnaja1_tidy, tmp_path / "kd")[0])
        for key in ("median_log10", "mad_log10", "z_max", "median_Kd", "excluded"):
            assert key in block, key

    def test_the_gate_statistics_are_scale_invariant(self, dnaja1_tidy, tmp_path):
        """Scaling every concentration shifts log10(Kd) by a constant, so the median
        moves by exactly that and the MAD does not move at all. If the MAD changed, the
        gate would be measuring the units rather than the spread."""
        base = self._block(_fit(dnaja1_tidy, tmp_path / "base")[0])
        scaled = self._block(_fit(dnaja1_tidy, tmp_path / "scaled",
                                  conc=CONC_X10, p0="500")[0])
        assert scaled["median_log10"] - base["median_log10"] == pytest.approx(1.0, abs=1e-6)
        # Not exact: per-residue fits that are poorly determined do not reproduce
        # bit-for-bit under rescaling, so the spread wobbles in the fifth figure.
        assert scaled["mad_log10"] == pytest.approx(base["mad_log10"], rel=1e-4)

    def test_extreme_residues_are_excluded(self, dnaja1_tidy, tmp_path):
        assert self._block(_fit(dnaja1_tidy, tmp_path / "kd")[0])["excluded"]

    def test_a_spread_based_test_would_have_sheltered_them(self, dnaja1_tidy, tmp_path):
        """Why physics runs before statistics. The implausible values inflate the spread
        enough to sit inside it: on this population mean/2-sigma shelters E20 and K23,
        and only the credible window catches them."""
        import statistics
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        excluded = set(self._block(data)["excluded"])
        assert excluded, "no outliers to test the comparison against"
        by_residue = _by_residue(data)

        # The population the gate sees is the pool PLUS whatever it removed — using the
        # survivors alone would exclude the very values under test.
        population = set(data["global"]["csp"]["residues"]) | excluded
        kds = [k for k in ((by_residue[r].get("csp") or {}).get("Kd") for r in population)
               if k is not None and math.isfinite(k) and k > 0]
        mean, sigma = statistics.fmean(kds), statistics.stdev(kds)
        sheltered = [r for r in excluded
                     if abs((by_residue[r]["csp"]["Kd"]) - mean) <= 2 * sigma]
        assert sheltered, ("mean/2-sigma would have caught every outlier here, so this "
                           "dataset does not demonstrate why MAD is needed")

    def test_the_excluded_residues_are_named_with_their_z_scores(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        excluded = self._block(data)["excluded"]
        reasons = data["metadata"]["csp_pool_excluded"]
        for residue in excluded:   # excluded maps residue -> robust z-score
            assert residue in reasons, f"{residue} excluded but absent from the map"
            assert re.search(r"\d", str(reasons[residue])), \
                f"{residue} reason carries no z-score: {reasons[residue]!r}"

    def test_the_window_exclusions_do_not_depend_on_the_z_threshold(
            self, dnaja1_tidy, tmp_path):
        """Residues outside what the titration can resolve are excluded on physics, so
        moving the statistical threshold must not change them."""
        sets = {z: frozenset(r for r, score in
                    self._block(_fit(dnaja1_tidy, tmp_path / f"z{z}",
                                     extra=["--kd-outlier-z", str(z)])[0])["excluded"].items()
                    if score is None)
                for z in (2.5, 3.0, 3.5)}
        assert len(set(sets.values())) == 1, sets

    def test_disabling_z_leaves_the_physics_gate_running(self, dnaja1_tidy, tmp_path):
        """z_max governs only the statistical test. A residue outside what the titration
        can resolve was not measured, and no flag makes that untrue."""
        off = self._block(_fit(dnaja1_tidy, tmp_path / "off",
                               extra=["--kd-outlier-z", "0"])[0])
        assert off["excluded"], "the window gate stopped running when z was disabled"
        assert all(score is None for score in off["excluded"].values()), \
            "a residue was scored by z after z was disabled"

    def test_disabling_z_still_reports_the_centre(self, dnaja1_tidy, tmp_path):
        """Credibility is a property of the population, not of the gate."""
        off = self._block(_fit(dnaja1_tidy, tmp_path / "off",
                               extra=["--kd-outlier-z", "0"])[0])
        assert off["median_credible"] is True
        assert "credible_window" in off

    def test_disabling_it_does_not_disturb_the_other_gates(self, dnaja1_tidy, tmp_path):
        """Turning one gate off must not quietly re-admit residues the others removed."""
        off = _fit(dnaja1_tidy, tmp_path / "off", extra=["--kd-outlier-z", "0"])[0]
        significant = set(off["metadata"]["csp_significance"]["significant"])
        assert set(off["global"]["csp"]["residues"]) <= significant

    def test_intensity_is_untouched(self, dnaja1_tidy, tmp_path):
        """CSP only — asserted against the gate being off, so a leak shows as a change."""
        on = _fit(dnaja1_tidy, tmp_path / "on")[0]
        off = _fit(dnaja1_tidy, tmp_path / "off", extra=["--kd-outlier-z", "0"])[0]
        assert (on["global"]["intensity"]["residues"]
                == off["global"]["intensity"]["residues"])
        assert on["global"]["intensity"]["Kd"] == off["global"]["intensity"]["Kd"]

    def test_five_gates_still_reconcile_in_one_place(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        pooled = set(data["global"]["csp"]["residues"])
        excluded = set(data["metadata"]["csp_pool_excluded"])
        assert pooled | excluded == set(_by_residue(data))
        assert pooled & excluded == set()

    def test_the_reasons_now_distinguish_five_gates(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        kinds = {str(v).split(":")[0].split(" below")[0]
                 for v in data["metadata"]["csp_pool_excluded"].values()}
        assert len(kinds) >= 4, kinds

    def test_a_persisted_threshold_is_applied(self, own_tidy, tmp_path):
        (own_tidy.parent / "kd_params.json").write_text(
            json.dumps({"kd_outlier_z": 0.0}))
        block = self._block(_fit(own_tidy, tmp_path / "kd")[0])
        assert all(score is None for score in block["excluded"].values())


class TestAnImplausibleCentreIsDeclaredNotSkipped:
    """If the typical residue fits a Kd outside what the titration can resolve, the
    finding is that most fits are meaningless — a whole-dataset statement. Skipping the
    outlier gate quietly would leave the output looking normal."""

    def _block(self, data):
        return data["metadata"]["csp_kd_outliers"]

    def test_erdj6_trips_it(self, erdj6_tidy, tmp_path):
        block = self._block(_fit(erdj6_tidy, tmp_path / "kd")[0])
        assert block["median_credible"] is False

    def test_dnaja1_does_not(self, dnaja1_tidy, tmp_path):
        """It has to discriminate, or it is just a banner on every run."""
        block = self._block(_fit(dnaja1_tidy, tmp_path / "kd")[0])
        assert block["median_credible"] is True

    def test_a_tripped_dataset_excludes_nothing(self, erdj6_tidy, tmp_path):
        """Warned-and-gated-anyway is the worst outcome: on this dataset four of five
        'outliers' had negative z, i.e. they were nearer a plausible Kd than the bulk."""
        assert self._block(_fit(erdj6_tidy, tmp_path / "kd")[0])["excluded"] == {}

    def test_the_verdict_reaches_the_metadata(self, erdj6_tidy, tmp_path):
        warnings = _fit(erdj6_tidy, tmp_path / "kd")[0]["metadata"]["quality_warnings"]
        assert warnings and any("MAKES NO SENSE" in w.upper() for w in warnings)

    def test_quality_warnings_is_always_a_list(self, dnaja1_tidy, tmp_path):
        """A consumer should never have to handle a missing key."""
        warnings = _fit(dnaja1_tidy, tmp_path / "kd")[0]["metadata"]["quality_warnings"]
        assert isinstance(warnings, list)

    def test_a_clean_dataset_carries_no_verdict(self, dnaja1_tidy, tmp_path):
        warnings = _fit(dnaja1_tidy, tmp_path / "kd")[0]["metadata"]["quality_warnings"]
        assert not any("MAKES NO SENSE" in w.upper() for w in warnings)

    def test_the_credible_window_is_one_decade_around_the_titration(
            self, erdj6_tidy, tmp_path):
        """Recorded only on the dataset that trips — see the note to the architect."""
        data, _ = _fit(erdj6_tidy, tmp_path / "kd")
        nonzero = [c for c in data["metadata"]["concentrations"] if c]
        low, high = self._block(data)["credible_window"]
        assert low == pytest.approx(min(nonzero) / 10.0)
        assert high == pytest.approx(max(nonzero) * 10.0)

    def test_the_window_matches_what_the_global_fit_bounds_to(self, erdj6_tidy, tmp_path):
        """The shared Kd is bounded to the same window. If these drift apart the tool
        contradicts itself about what the experiment can resolve — ERDJ6 pins its
        shared Kd at that bound, so the two are directly comparable here."""
        data, _ = _fit(erdj6_tidy, tmp_path / "kd")
        _low, high = self._block(data)["credible_window"]
        assert data["global"]["csp"]["Kd"] == pytest.approx(high, rel=1e-3)
        assert data["global"]["csp"]["reliable"] is False

    def test_disabling_the_gate_does_not_suppress_the_verdict(self, erdj6_tidy, tmp_path):
        """The warning is about the data, not about the gate."""
        data, _ = _fit(erdj6_tidy, tmp_path / "kd", extra=["--kd-outlier-z", "0"])
        assert any("MAKES NO SENSE" in w.upper()
                   for w in data["metadata"]["quality_warnings"])

    def test_the_verdict_names_the_number_and_the_window(self, erdj6_tidy, tmp_path):
        """A verdict without its measurement is an opinion."""
        warnings = _fit(erdj6_tidy, tmp_path / "kd")[0]["metadata"]["quality_warnings"]
        verdict = next(w for w in warnings if "MAKES NO SENSE" in w.upper())
        assert re.search(r"\d", verdict)


class TestPhysicsGatesBeforeStatistics:
    """MAD adapts to the population, so a spread-based test cannot exclude values that
    are themselves inflating the spread — |z|>3 on a 0.65-decade MAD permits 0.11 to
    71,000 on a titration topping out at 150. A residue whose own Kd lies outside what
    the titration can resolve was not measured, whatever the population looks like."""

    def _block(self, data):
        return data["metadata"]["csp_kd_outliers"]

    def test_residues_outside_the_window_are_excluded(self, dnaja1_tidy, tmp_path):
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        assert {"E20", "E41", "K23"} <= set(self._block(data)["excluded"])

    def test_the_reason_names_the_window_not_a_z_score(self, dnaja1_tidy, tmp_path):
        """Window and z are different claims and must not read the same."""
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        reasons = data["metadata"]["csp_pool_excluded"]
        for residue in ("E20", "E41", "K23"):
            assert "resolve" in str(reasons[residue]), reasons[residue]

    def test_a_window_exclusion_carries_no_z_score(self, dnaja1_tidy, tmp_path):
        """They were never scored, so reporting a z would invent a measurement."""
        block = self._block(_fit(dnaja1_tidy, tmp_path / "kd")[0])
        for residue in ("E20", "E41", "K23"):
            assert block["excluded"][residue] is None

    def test_the_centre_is_measured_on_survivors(self, dnaja1_tidy, tmp_path):
        """Taken before the window gate the median sits at 89.4; on the survivors it is
        lower, because the implausible values are gone."""
        block = self._block(_fit(dnaja1_tidy, tmp_path / "kd")[0])
        low, high = block["credible_window"]
        assert block["median_Kd"] < 89.0
        assert low < block["median_Kd"] < high

    def test_erdj6_raises_the_verdict_and_gates_nothing(self, erdj6_tidy, tmp_path):
        """Untrustworthy centre and implausible individuals are different failures and
        get different responses: refuse to gate and shout, versus gate the individuals."""
        data, _ = _fit(erdj6_tidy, tmp_path / "kd")
        block = self._block(data)
        assert block["median_credible"] is False
        assert block["excluded"] == {}
        assert any("MAKES NO SENSE" in w.upper()
                   for w in data["metadata"]["quality_warnings"])

    @pytest.mark.parametrize("tidy_fixture", ["dnaja1_tidy", "erdj6_tidy"])
    def test_one_window_is_used_everywhere(self, request, tidy_fixture, tmp_path):
        """The window bounds the shared Kd, decides whether the median is credible, and
        gates individual residues. If those three ever diverge the tool contradicts
        itself about what this experiment can resolve."""
        tidy = request.getfixturevalue(tidy_fixture)
        data, _ = _fit(tidy, tmp_path / "kd")
        nonzero = [c for c in data["metadata"]["concentrations"] if c]
        low, high = self._block(data)["credible_window"]
        assert low == pytest.approx(min(nonzero) / 10.0)
        assert high == pytest.approx(max(nonzero) * 10.0)
        assert low <= data["global"]["csp"]["Kd"] <= high


CLI_AGENT_DOC = REPO_ROOT / "docs" / "CLI_AGENT.md"


class TestTheDocumentedWorkflowIsReal:
    """CLI_AGENT.md is read by agents that will run what it says verbatim. A flag that
    does not exist, or a path that moved, is a bug in the deliverable."""

    @pytest.fixture(scope="class")
    def doc_text(self):
        if not CLI_AGENT_DOC.exists():
            pytest.skip("docs/CLI_AGENT.md not present")
        return CLI_AGENT_DOC.read_text()

    @pytest.fixture(scope="class")
    def documented_kd_flags(self, doc_text):
        """Only the runnable command lines — those are what an agent copies."""
        flags = set()
        for line in doc_text.splitlines():
            if "lunaNMR kd" in line or "lunaNMR export kd" in line:
                flags.update(re.findall(r"(?<![\w-])--[a-z][a-z0-9-]+", line))
        assert flags, "no kd command lines found in the document"
        return flags

    def _accepted(self, subcommand):
        import subprocess
        out = subprocess.run([sys.executable, "-m", "lunaNMR", *subcommand, "--help"],
                             cwd=str(REPO_ROOT), capture_output=True, text=True)
        assert out.returncode == 0, out.stderr
        return set(re.findall(r"(?<![\w-])--[a-z][a-z0-9-]+", out.stdout))

    def test_every_documented_kd_flag_exists(self, documented_kd_flags):
        """The flags added today are the ones most likely to be documented wrong."""
        real = self._accepted(["kd"]) | self._accepted(["export", "kd"]) \
            | self._accepted(["series"]) | {"--help", "--format", "--dry-run"}
        unknown = {f for f in documented_kd_flags if f not in real}
        assert unknown == set(), f"documented but not accepted: {sorted(unknown)}"

    @pytest.mark.parametrize("flag", ["--survey", "--residues", "--selection",
                                      "--conc-units", "--csp-sigma-multiple",
                                      "--kd-outlier-z", "--ref-max-ratio",
                                      "--dd-runaway-ratio", "--noise-quantile"])
    def test_todays_flags_are_documented(self, flag, doc_text):
        """A flag nobody can discover may as well not exist."""
        assert flag in doc_text

    def test_the_documented_json_path_matches_the_layout(self, doc_text):
        """export kd reads <out>/data/<prefix>_kd_fit_data.json since the layout moved.
        Scoped to runnable command lines — prose and table rows name the flag without a
        path, and an agent copies commands."""
        commands = [ln for ln in doc_text.splitlines()
                    if "python -m lunaNMR export kd" in ln and "--json" in ln]
        assert commands, "no runnable export kd command documented"
        assert all("/data/" in ln for ln in commands), commands

    def test_the_equivalents_trap_is_stated_as_open(self):
        """It is avoidable, not closed — absolute is still the default."""
        text = CLI_AGENT_DOC.read_text().lower()
        assert "equivalent" in text
        assert "--conc-units" in text


class TestTheFigureAndTheFitAgreeOnAUsableReference:
    """A broken reference was rejected by the fit path and accepted by the figure path,
    so residues the fitter had already dropped were plotted at ratios up to 3.5e+09 and
    one bar flattened all sixty-eight real ones. Two paths, one convention."""

    def _models(self):
        sys.path.insert(0, str(REPO_ROOT / "modules" / "dynamiXs_v2o0" / "dynamiXs_Kd"))
        import kd_models
        return kd_models

    def test_one_definition_of_a_usable_reference_exists(self):
        """kd_input must import it rather than keep a second copy that can drift."""
        models = self._models()
        import kd_input
        assert kd_input.REF_MAX_RATIO is models.REF_MAX_RATIO
        assert kd_input.reference_usable is models.reference_usable

    def test_the_two_paths_agree_on_the_same_residues(self, dnaja1_tidy, tmp_path):
        """The property that matters is agreement, not that each behaves in isolation."""
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        dropped_by_fit = {r for r, f in _by_residue(data).items()
                          if not (f.get("intensity") or {}).get("success")}
        assert dropped_by_fit, "no residue was dropped, so agreement is untested here"

        import numpy as np
        models = self._models()
        checked = 0
        for residue in dropped_by_fit:
            series = _by_residue(data)[residue]["series"]
            heights = np.asarray([h if h is not None else np.nan
                                  for h in series["height"]], float)
            if heights.size < 2 or not np.isfinite(heights[0]):
                continue
            if models.reference_usable(heights, heights[0]):
                continue                      # the fit dropped it for another reason
            checked += 1
            ratio = models.pair_observable(series, 0, heights.size - 1, "intensity")
            assert not np.isfinite(ratio), \
                f"{residue}: reference rejected by the fit, still plotted at {ratio:.4g}"
        assert checked, "no broken-reference residue reached the comparison"

    def test_no_exported_ratio_exceeds_the_shared_limit(self, dnaja1_tidy, tmp_path):
        """A property of the deliverable, so it survives however the figure path is
        later refactored — the CSV is what the figure draws."""
        from lunaNMR.cli import main
        data, _ = _fit(dnaja1_tidy, tmp_path / "kd")
        figs = tmp_path / "figs"
        assert main(["export", "kd", "--json",
                     str(tmp_path / "kd" / "data" / "t_kd_fit_data.json"),
                     "--out", str(figs), "--prefix", "t", "--kind", "ref-bars"]) == 0

        table = pd.read_csv(figs / "data" / "t_intensity_ref_vs_point.csv")
        numeric = table.select_dtypes("number").stack().dropna()
        assert not numeric.empty
        assert numeric.max() <= self._models().REF_MAX_RATIO, \
            f"max exported I/I0 = {numeric.max():.4g}"

    def test_a_real_ratio_is_still_plotted(self, dnaja1_tidy, tmp_path):
        """Guards over-correction: rejecting everything would also satisfy the cap."""
        from lunaNMR.cli import main
        _fit(dnaja1_tidy, tmp_path / "kd")
        figs = tmp_path / "figs"
        assert main(["export", "kd", "--json",
                     str(tmp_path / "kd" / "data" / "t_kd_fit_data.json"),
                     "--out", str(figs), "--prefix", "t", "--kind", "ref-bars"]) == 0
        table = pd.read_csv(figs / "data" / "t_intensity_ref_vs_point.csv")
        numeric = table.select_dtypes("number").stack().dropna()
        assert (numeric > 0).sum() > 40, "almost nothing survived — the cap is too tight"
