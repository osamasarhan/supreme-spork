"""
ranking.py - Z-score based ranking system used in the CASP16 assessment.

The methodology:
1. Compute raw metric scores for all groups on a given target.
2. Remove outliers whose score is ≥ 2 standard deviations *below* the mean.
3. Re-compute mean and standard deviation on the filtered set.
4. Convert every score (including the removed outliers) to a Z-score.
5. Negative Z-scores are clamped to 0 (poor predictions are not penalised).
6. Per-target Z-scores are combined with category-specific weights to produce
   an overall ranking.

Category weights
----------------
NA Monomers
    Z = 0.3·Z_TM + 0.3·Z_GDT + 0.4·Z_lDDT

RNA Multimers / NA-Protein targets
    Z = 0.3·(0.3·Z_TM + 0.3·Z_GDT + 0.4·Z_lDDT)
      + 0.7·(⅓·Z_ICS  + ⅓·Z_IPS  + ⅓·Z_i-lDDT)
"""

import numpy as np
from typing import List, Sequence


def calculate_z_scores(scores: Sequence[float]) -> List[float]:
    """Compute Z-scores with outlier removal.

    Outliers are defined as values **more than 2 standard deviations below
    the mean** of the input *scores*.  They are excluded when computing the
    reference mean and standard deviation but are still assigned a Z-score
    using those reference statistics (which will consequently be very
    negative).

    Parameters
    ----------
    scores:
        Raw metric values for all groups on a single target.

    Returns
    -------
    List of Z-scores, one per element in *scores*, in the same order.
    If *scores* contains fewer than 2 elements, or if the filtered standard
    deviation is zero, all Z-scores are returned as 0.0.

    Examples
    --------
    >>> calculate_z_scores([0.9, 0.8, 0.7, 0.1])
    [0.707..., 0.0, -0.707..., -2.828...]
    """
    arr = np.array(scores, dtype=float)
    if arr.size < 2:
        return [0.0] * arr.size

    mean_val = np.mean(arr)
    std_val = np.std(arr)  # population std (ddof=0) for initial outlier filter

    # Remove outliers (> 2 std dev below the mean)
    threshold = mean_val - 2.0 * std_val
    filtered = arr[arr > threshold]

    if filtered.size < 2:
        return [0.0] * arr.size

    final_mean = np.mean(filtered)
    final_std = np.std(filtered)

    if final_std == 0.0:
        return [0.0] * arr.size

    return list((arr - final_mean) / final_std)


def compute_na_monomer_ranking(
    z_tm: float, z_gdt: float, z_lddt: float
) -> float:
    """Weighted Z-score for Nucleic Acid Monomer targets.

    Parameters
    ----------
    z_tm:   Z-score for TM-score.
    z_gdt:  Z-score for GDT_TS.
    z_lddt: Z-score for lDDT.

    Returns
    -------
    Combined Z-score clamped at 0 (negative values become 0).

    Notes
    -----
    Formula: Z = 0.3·Z_TM + 0.3·Z_GDT + 0.4·Z_lDDT
    """
    z_score = 0.3 * z_tm + 0.3 * z_gdt + 0.4 * z_lddt
    return max(0.0, z_score)


def compute_multimer_ranking(
    z_tm: float,
    z_gdt: float,
    z_lddt: float,
    z_ics: float,
    z_ips: float,
    z_i_lddt: float,
) -> float:
    """Weighted Z-score for RNA-multimer and NA-protein targets.

    Integrates global fold metrics (TM-score, GDT_TS, lDDT) with interface
    quality metrics (ICS, IPS, i-lDDT).

    Parameters
    ----------
    z_tm:     Z-score for TM-score.
    z_gdt:    Z-score for GDT_TS.
    z_lddt:   Z-score for lDDT.
    z_ics:    Z-score for Interface Contact Score (F1).
    z_ips:    Z-score for Interface Patch Score (Jaccard).
    z_i_lddt: Z-score for interface lDDT.

    Returns
    -------
    Combined Z-score clamped at 0 (negative values become 0).

    Notes
    -----
    Formula:
        Z_global    = 0.3·Z_TM + 0.3·Z_GDT + 0.4·Z_lDDT
        Z_interface = ⅓·Z_ICS + ⅓·Z_IPS + ⅓·Z_i-lDDT
        Z           = 0.3·Z_global + 0.7·Z_interface
    """
    global_score = 0.3 * z_tm + 0.3 * z_gdt + 0.4 * z_lddt
    interface_score = (z_ics + z_ips + z_i_lddt) / 3.0
    z_score = 0.3 * global_score + 0.7 * interface_score
    return max(0.0, z_score)


def rank_groups(group_z_scores: dict) -> List[tuple]:
    """Sort groups by their combined Z-score in descending order.

    Parameters
    ----------
    group_z_scores:
        Mapping of group identifier → combined Z-score.

    Returns
    -------
    List of ``(rank, group_id, z_score)`` tuples, starting at rank 1.
    """
    sorted_groups = sorted(
        group_z_scores.items(), key=lambda kv: kv[1], reverse=True
    )
    return [(rank + 1, gid, score) for rank, (gid, score) in enumerate(sorted_groups)]
