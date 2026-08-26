# ABOUTME: Unit tests for the Kd survey's exclusion evidence, verdicts and selection file.
# ABOUTME: Drives kd_survey directly on constructed rows so reasons are testable without argparse.

import sys
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_Kd"
sys.path.insert(0, str(_DIR))

from kd_survey import (_add_verdicts, format_residues_file, parse_residues_file,
                       reason_text, survey_table, DD_RUNAWAY_RATIO)

# The noise floor is the 25th percentile of the rows' own max_csp, so a row is only
# "moving" relative to its companions. Tests supply a spread and target one member.
_COMPANIONS = [0.040, 0.050, 0.060, 0.070]


def _row(residue, max_csp=0.06, **kw):
    """A residue that trips nothing, overridden per test."""
    row = {'residue': residue, 'max_csp': max_csp, 'intensity_final': 0.3,
           'dd_ratio': 1.2, 'n_failed_points': 0, 'n_points_csp': 7,
           'n_points_intensity': 7, 'csp_r_squared': 0.95,
           'intensity_r_squared': 0.95}
    row.update(kw)
    return row


def _verdicts(subject, **kw):
    """Score `subject` against a normal-looking cohort and return its row."""
    rows = [_row(subject, **kw)] + [_row(f"X{i}", max_csp=v)
                                    for i, v in enumerate(_COMPANIONS)]
    return _add_verdicts(rows)[0]


def _keys(subject, **kw):
    return set(_verdicts(subject, **kw)['reason_keys'])


class TestEvidenceIsSeparateFromWording:
    """reason_keys must be assertable without matching prose, so severity ordering and
    wording can change without rewriting every test."""

    def test_a_clean_residue_trips_nothing(self):
        assert _keys('K3') == set()

    def test_an_unconstrained_plateau_is_its_own_key(self):
        assert _keys('E41', dd_ratio=DD_RUNAWAY_RATIO * 50) == {'plateau_unconstrained'}

    def test_a_residue_below_the_noise_floor_is_its_own_key(self):
        assert _keys('A17', max_csp=0.001) == {'no_motion'}

    def test_a_rejected_reference_is_its_own_key(self):
        """intensity_ratio_series blanks every point when the reference is unusable."""
        assert _keys('T5', n_points_intensity=0) == {'reference_unusable'}

    def test_missing_positions_are_their_own_key(self):
        assert _keys('D36', n_points_csp=1) == {'no_position_data'}

    def test_lost_points_are_their_own_key(self):
        assert _keys('G68', n_failed_points=2) == {'failed_points'}

    def test_a_placeholder_is_recognised_by_name(self):
        assert 'dummy' in _keys('dummy_001')

    def test_the_wording_renders_from_the_keys(self):
        text = reason_text(_verdicts('E41', dd_ratio=400.0))
        assert text == ['plateau 400x its own largest CSP — not constrained by this titration']

    def test_the_wording_quotes_the_number_behind_the_key(self):
        """A reason without its measurement is an opinion."""
        assert '0.0010' in reason_text(_verdicts('A17', max_csp=0.001))[0]


class TestTwoReasonsAreBothReported:
    """A residue failing on two counts must say both — citing one hides the other."""

    def test_both_keys_fire(self):
        assert _keys('K3', max_csp=0.001, dd_ratio=400.0) == {'no_motion',
                                                              'plateau_unconstrained'}

    def test_both_are_rendered(self):
        assert len(reason_text(_verdicts('K3', max_csp=0.001, dd_ratio=400.0))) == 2

    def test_unusable_reasons_are_ordered_before_debatable_ones(self):
        keys = _verdicts('T5', n_points_intensity=0, dd_ratio=400.0)['reason_keys']
        assert keys.index('reference_unusable') < keys.index('plateau_unconstrained')


class TestVerdictSeverity:
    """Unusable residues are commented out; debatable ones are annotated and kept.
    A suggestion that excludes everything it can criticise is not a suggestion."""

    def test_a_rejected_reference_is_unusable(self):
        assert _verdicts('T5', n_points_intensity=0)['verdict'] == 'unusable'

    def test_a_placeholder_is_unusable(self):
        assert _verdicts('dummy_001')['verdict'] == 'unusable'

    def test_an_unconstrained_plateau_is_only_worth_checking(self):
        assert _verdicts('Y8', dd_ratio=400.0)['verdict'] == 'check'

    def test_a_non_mover_is_only_worth_checking(self):
        assert _verdicts('A17', max_csp=0.001)['verdict'] == 'check'

    def test_a_clean_residue_is_ok(self):
        assert _verdicts('E4')['verdict'] == 'ok'

    def test_a_debatable_residue_survives_the_round_trip(self):
        rows = _add_verdicts([_row('E4'), _row('Y8', dd_ratio=400.0),
                              _row('T5', n_points_intensity=0)])
        selected = parse_residues_file(format_residues_file(rows))
        assert 'Y8' in selected, "a debatable residue was excluded, not annotated"
        assert 'T5' not in selected


class TestSelectionFileFormat:
    def test_a_selected_residue_round_trips(self):
        rows = _add_verdicts([_row('K3'), _row('E4')])
        assert parse_residues_file(format_residues_file(rows)) == ['K3', 'E4']

    def test_uncommenting_an_exclusion_yields_a_clean_name(self):
        """The reason must sit behind its own '#', or deleting the leading one welds
        the residue name to its reason and the parser blames the user's input."""
        rows = _add_verdicts([_row('K3'), _row('T5', n_points_intensity=0)])
        text = format_residues_file(rows)
        line = next(ln for ln in text.splitlines()
                    if ln.startswith("# ") and ln[2:].startswith("T5"))
        assert parse_residues_file(text.replace(line, line[2:])) == ['K3', 'T5']

    def test_an_annotated_residue_keeps_its_evidence_visible(self):
        rows = _add_verdicts([_row('Y8', dd_ratio=400.0)])
        line = next(ln for ln in format_residues_file(rows).splitlines()
                    if ln.startswith('Y8'))
        assert '#' in line and 'plateau' in line

    def test_a_comment_only_file_selects_nothing(self):
        assert parse_residues_file("# K3\n# E4\n") == []

    def test_blank_lines_are_ignored(self):
        assert parse_residues_file("K3\n\n\nE4\n") == ['K3', 'E4']

    def test_surrounding_whitespace_is_stripped(self):
        assert parse_residues_file("  K3  \n\tE4\n") == ['K3', 'E4']


class TestSurveyTable:
    """survey_table returns (column_names, rows) ready for CSV writing."""

    def test_every_surveyed_residue_gets_a_row(self):
        rows = _add_verdicts([_row('K3'), _row('E4'), _row('T5', n_points_intensity=0)])
        _columns, table = survey_table(rows)
        assert len(table) == 3

    def test_the_numbers_behind_a_reason_are_carried_as_fields(self):
        """Evidence must be inspectable as data, not only readable as prose."""
        columns, table = survey_table(_add_verdicts([_row('E41', dd_ratio=400.0)]))
        assert float(dict(zip(columns, table[0]))['dd_ratio']) == pytest.approx(400.0)

    def test_the_verdict_and_reasons_are_columns(self):
        columns, _table = survey_table(_add_verdicts([_row('K3')]))
        assert 'verdict' in columns and 'reasons' in columns

    def test_an_empty_dataset_does_not_raise(self):
        columns, table = survey_table(_add_verdicts([]))
        assert table == []
        assert columns
