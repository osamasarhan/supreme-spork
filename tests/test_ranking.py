"""Tests for casp16_eval.ranking - Z-score computation and ranking."""

import math
import pytest

from casp16_eval.ranking import (
    calculate_z_scores,
    compute_na_monomer_ranking,
    compute_multimer_ranking,
    rank_groups,
)


class TestCalculateZScores:
    def test_single_element_returns_zero(self):
        assert calculate_z_scores([0.8]) == [0.0]

    def test_two_equal_elements_return_zeros(self):
        z = calculate_z_scores([0.5, 0.5])
        assert all(v == pytest.approx(0.0) for v in z)

    def test_symmetric_two_elements(self):
        z = calculate_z_scores([1.0, 0.0])
        # Mean=0.5, pstdev=0.5 → z = (x-0.5)/0.5
        assert z[0] == pytest.approx(1.0)
        assert z[1] == pytest.approx(-1.0)

    def test_outlier_excluded_from_reference_stats(self):
        # 9 tightly-clustered scores (0.5–0.9) plus one outlier (0.0).
        # initial mean ≈ 0.63, pstdev ≈ 0.24, threshold ≈ 0.14
        # 0.0 is below threshold → filtered out of reference stats
        # final mean ≈ 0.70, final std ≈ 0.13 → z(0.0) ≈ -5.4
        scores = [0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.0]
        z = calculate_z_scores(scores)
        # The outlier (0.0) should receive a very negative z-score
        assert z[-1] < -2.0

    def test_length_preserved(self):
        scores = [0.9, 0.8, 0.7, 0.6, 0.5]
        z = calculate_z_scores(scores)
        assert len(z) == len(scores)

    def test_empty_returns_empty(self):
        assert calculate_z_scores([]) == []

    def test_higher_score_higher_z(self):
        scores = [0.9, 0.5, 0.1]
        z = calculate_z_scores(scores)
        assert z[0] > z[1] > z[2]


class TestComputeNaMonomerRanking:
    def test_positive_z_scores(self):
        result = compute_na_monomer_ranking(1.0, 1.0, 1.0)
        assert result == pytest.approx(1.0)

    def test_formula_weights(self):
        # 0.3*2 + 0.3*1 + 0.4*0 = 0.9
        result = compute_na_monomer_ranking(z_tm=2.0, z_gdt=1.0, z_lddt=0.0)
        assert result == pytest.approx(0.9)

    def test_negative_clamped_to_zero(self):
        result = compute_na_monomer_ranking(-1.0, -1.0, -1.0)
        assert result == 0.0

    def test_mixed_clamps_when_negative(self):
        # 0.3*(-2) + 0.3*0 + 0.4*0 = -0.6 → clamp → 0
        result = compute_na_monomer_ranking(-2.0, 0.0, 0.0)
        assert result == 0.0


class TestComputeMultimerRanking:
    def test_all_ones(self):
        # global = 0.3*1+0.3*1+0.4*1 = 1.0
        # interface = (1+1+1)/3 = 1.0
        # final = 0.3*1 + 0.7*1 = 1.0
        result = compute_multimer_ranking(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        assert result == pytest.approx(1.0)

    def test_negative_clamped(self):
        result = compute_multimer_ranking(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0)
        assert result == 0.0

    def test_global_interface_weighting(self):
        # global part = 0.3*(1+1+1) = 0.9 (wait, formula = 0.3*z_tm+0.3*z_gdt+0.4*z_lddt)
        # with all z=1: global=1, interface=1 → final=0.3+0.7=1.0
        # Test with zero interface metrics
        # global = 0.3*1+0.3*1+0.4*1 = 1; interface = 0
        # final = 0.3*1 + 0.7*0 = 0.3
        result = compute_multimer_ranking(1.0, 1.0, 1.0, 0.0, 0.0, 0.0)
        assert result == pytest.approx(0.3)


class TestRankGroups:
    def test_descending_order(self):
        scores = {"group_A": 2.5, "group_B": 1.0, "group_C": 3.0}
        ranking = rank_groups(scores)
        assert ranking[0][1] == "group_C"
        assert ranking[1][1] == "group_A"
        assert ranking[2][1] == "group_B"

    def test_rank_starts_at_one(self):
        ranking = rank_groups({"A": 1.0, "B": 0.5})
        assert ranking[0][0] == 1
        assert ranking[1][0] == 2

    def test_empty_input(self):
        assert rank_groups({}) == []
