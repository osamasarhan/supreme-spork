"""Tests for casp16_eval.alignment - US-align output parser."""

import pytest

from casp16_eval.alignment import _parse_us_align_output


# ---------------------------------------------------------------------------
# Sample US-align stdout (abridged, realistic format)
# ---------------------------------------------------------------------------

SAMPLE_OUTPUT = """\
 **************************************************************************
 *                        US-align (Version 20220511)                     *
 * Universal Structure Alignment of Proteins and Nucleic Acids            *
 **************************************************************************

Name of Structure_1: target.pdb
Name of Structure_2: template.pdb
Length of Structure_1: 73 residues
Length of Structure_2: 76 residues

Aligned length= 71, RMSD=   1.93, Seq_ID=n_identical/n_aligned= 0.831

TM-score= 0.9234 (if normalized by length of Structure_1, i.e., L=73, d0=2.77)
TM-score= 0.8917 (if normalized by length of Structure_2, i.e., L=76, d0=2.84)

Number of residues in common=  71
"""

MISSING_FIELDS_OUTPUT = """\
 **************************************************************************
 *                        US-align (Version 20220511)                     *
 **************************************************************************

Name of Structure_1: target.pdb
Name of Structure_2: template.pdb
"""


class TestParseUsAlignOutput:
    def test_tm_score_1_parsed(self):
        result = _parse_us_align_output(SAMPLE_OUTPUT)
        assert result["tm_score_1"] == pytest.approx(0.9234)

    def test_tm_score_2_parsed(self):
        result = _parse_us_align_output(SAMPLE_OUTPUT)
        assert result["tm_score_2"] == pytest.approx(0.8917)

    def test_rmsd_parsed(self):
        result = _parse_us_align_output(SAMPLE_OUTPUT)
        assert result["rmsd"] == pytest.approx(1.93)

    def test_seq_id_parsed(self):
        result = _parse_us_align_output(SAMPLE_OUTPUT)
        assert result["seq_id"] == pytest.approx(0.831)

    def test_aligned_len_parsed(self):
        result = _parse_us_align_output(SAMPLE_OUTPUT)
        assert result["aligned_len"] == 71

    def test_stdout_preserved(self):
        result = _parse_us_align_output(SAMPLE_OUTPUT)
        assert result["stdout"] == SAMPLE_OUTPUT

    def test_missing_fields_return_none(self):
        result = _parse_us_align_output(MISSING_FIELDS_OUTPUT)
        assert result["tm_score_1"] is None
        assert result["tm_score_2"] is None
        assert result["rmsd"] is None
        assert result["seq_id"] is None
        assert result["aligned_len"] is None
