# ABOUTME: Tests that peak-list assignments are stripped of surrounding whitespace on load.
# ABOUTME: A list written "3LysH, 8.2, 126.3" carried the space into every output CSV.

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from lunaNMR.utils.file_manager import NMRFileManager


class TestAssignmentWhitespace:
    def test_space_after_comma_is_stripped(self, tmp_path):
        """The conventional "name, x, y" spacing puts a leading space in every
        assignment. Residue matching across datasets is by exact string, so
        ' 3LysH' and '3LysH' silently share no residues at all."""
        peaks = tmp_path / "spaced.txt"
        peaks.write_text(
            "Assignment, Position_X, Position_Y\n"
            " 3LysH, 8.20507, 126.30044\n"
            " 10ValH, 8.22083, 120.02435\n"
        )

        df = NMRFileManager().load_peak_list(str(peaks))

        assert df['Assignment'].tolist() == ['3LysH', '10ValH']

    def test_unspaced_list_is_unchanged(self, tmp_path):
        peaks = tmp_path / "tight.txt"
        peaks.write_text(
            "Assignment,Position_X,Position_Y\n"
            "3LysH,8.179464,126.462847\n"
        )

        df = NMRFileManager().load_peak_list(str(peaks))

        assert df['Assignment'].tolist() == ['3LysH']

    def test_trailing_whitespace_is_stripped(self, tmp_path):
        peaks = tmp_path / "trailing.txt"
        peaks.write_text(
            "Assignment,Position_X,Position_Y\n"
            "3LysH\t,8.179464,126.462847\n"
        )

        df = NMRFileManager().load_peak_list(str(peaks))

        assert df['Assignment'].tolist() == ['3LysH']

    def test_positions_survive_the_strip(self, tmp_path):
        peaks = tmp_path / "spaced.txt"
        peaks.write_text(
            "Assignment, Position_X, Position_Y\n"
            " 3LysH, 8.20507, 126.30044\n"
        )

        df = NMRFileManager().load_peak_list(str(peaks))

        assert df['Position_X'].iloc[0] == 8.20507
        assert df['Position_Y'].iloc[0] == 126.30044
