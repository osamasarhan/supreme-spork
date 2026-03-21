"""
casp16_eval - Python toolkit for replicating the CASP16 nucleic acid
structure prediction assessment.

Modules
-------
alignment           : US-align wrapper for template identification
metrics             : 3D-structure quality metrics (TM-score, GDT_TS, lDDT, INF)
ranking             : Z-score based ranking for NA monomers and multimers
secondary_structure : Secondary / tertiary structure analysis (DSSR, Rosetta rna_motif)
interface           : Multimer interface evaluation (ICS, IPS, DockQ) and symmetry (AnAnaS)
"""

from .alignment import run_us_align
from .metrics import compute_global_metrics
from .ranking import (
    calculate_z_scores,
    compute_na_monomer_ranking,
    compute_multimer_ranking,
)
from .secondary_structure import (
    evaluate_secondary_structure,
    evaluate_tertiary_motifs,
    compute_f1_score,
)
from .interface import check_symmetry, compute_interface_metrics

__all__ = [
    "run_us_align",
    "compute_global_metrics",
    "calculate_z_scores",
    "compute_na_monomer_ranking",
    "compute_multimer_ranking",
    "evaluate_secondary_structure",
    "evaluate_tertiary_motifs",
    "compute_f1_score",
    "check_symmetry",
    "compute_interface_metrics",
]
