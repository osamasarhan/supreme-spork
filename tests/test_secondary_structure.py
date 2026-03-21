"""Tests for casp16_eval.secondary_structure - base-pair F1 and DSSR parser."""

import pytest

from casp16_eval.secondary_structure import (
    _parse_dssr_base_pairs,
    compute_f1_score,
)


# ---------------------------------------------------------------------------
# Sample DSSR tabular output (abridged)
# ---------------------------------------------------------------------------

DSSR_TABULAR = """\
# DSSR (v1.9.9-2020feb06) output
# base pairs
A.G1   A.C20
A.C2   A.G19
A.U3   A.A18
# stacking
A.G1   A.C2
"""


class TestParseDssrBasePairs:
    def test_parses_three_pairs(self):
        pairs = _parse_dssr_base_pairs(DSSR_TABULAR)
        assert len(pairs) == 3

    def test_pair_content(self):
        pairs = _parse_dssr_base_pairs(DSSR_TABULAR)
        assert frozenset({"A.G1", "A.C20"}) in pairs

    def test_stacking_not_included(self):
        pairs = _parse_dssr_base_pairs(DSSR_TABULAR)
        # After "# stacking" the section changes; no pair from that block
        assert frozenset({"A.G1", "A.C2"}) not in pairs

    def test_empty_output_returns_empty_list(self):
        assert _parse_dssr_base_pairs("") == []


class TestComputeF1Score:
    def test_perfect_match(self):
        pairs = [frozenset({"A.G1", "A.C10"}), frozenset({"A.U2", "A.A9"})]
        assert compute_f1_score(pairs, pairs) == pytest.approx(1.0)

    def test_no_match(self):
        pred = [frozenset({"A.G1", "A.C10"})]
        ref = [frozenset({"A.U2", "A.A9"})]
        assert compute_f1_score(pred, ref) == pytest.approx(0.0)

    def test_partial_match(self):
        pred = [frozenset({"A.G1", "A.C10"}), frozenset({"A.U2", "A.A9"})]
        ref = [frozenset({"A.G1", "A.C10"}), frozenset({"A.C3", "A.G8"})]
        # TP=1, precision=0.5, recall=0.5 → F1=0.5
        assert compute_f1_score(pred, ref) == pytest.approx(0.5)

    def test_empty_predicted_returns_zero(self):
        ref = [frozenset({"A.G1", "A.C10"})]
        assert compute_f1_score([], ref) == pytest.approx(0.0)

    def test_empty_reference_returns_zero(self):
        pred = [frozenset({"A.G1", "A.C10"})]
        assert compute_f1_score(pred, []) == pytest.approx(0.0)

    def test_both_empty_returns_zero(self):
        assert compute_f1_score([], []) == pytest.approx(0.0)

    def test_recall_dominated(self):
        # pred has only 1 of 3 ref pairs → recall=1/3, precision=1/1
        pred = [frozenset({"A.G1", "A.C10"})]
        ref = [
            frozenset({"A.G1", "A.C10"}),
            frozenset({"A.U2", "A.A9"}),
            frozenset({"A.C3", "A.G8"}),
        ]
        expected = 2 * (1 / 1) * (1 / 3) / (1 / 1 + 1 / 3)
        assert compute_f1_score(pred, ref) == pytest.approx(expected)
