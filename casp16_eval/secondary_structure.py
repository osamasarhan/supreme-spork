"""
secondary_structure.py - Secondary and tertiary RNA structure analysis.

Tools used in the CASP16 assessment
-------------------------------------
* **DSSR** (v1.9.9-2020feb06) – base-pair / secondary structure extraction.
* **Rosetta rna_motif**        – tertiary motif identification (A-Minor,
  T-Loop, GNRA Tetraloop, …).

Base-pair matching
------------------
A predicted base pair is counted as a *true positive* only when **both**
nucleotide identifiers match exactly (chain ID, residue number, and
nucleotide type) against the experimental reference.

The overall quality is summarised by the **F1-score** (harmonic mean of
precision and recall).
"""

import re
import subprocess
from typing import FrozenSet, List, Optional, Set, Tuple


# Type alias: a base pair is represented as a frozenset of two residue ids
BasePair = FrozenSet[str]


def evaluate_secondary_structure(
    model_pdb: str,
    dssr_executable: str = "x3dna-dssr",
) -> Optional[List[BasePair]]:
    """Extract base pairs from *model_pdb* using DSSR.

    The CASP16 authors manually extracted the base-pair list from the DSSR
    table output (``--tabular`` mode).

    Parameters
    ----------
    model_pdb:
        Path to the PDB model to analyse.
    dssr_executable:
        Name / full path of the ``x3dna-dssr`` binary.

    Returns
    -------
    List of base pairs as frozensets of two residue identifier strings, or
    ``None`` when DSSR is unavailable.

    Raises
    ------
    subprocess.CalledProcessError
        When DSSR exits with a non-zero status.
    """
    try:
        result = subprocess.run(
            [dssr_executable, "-i", model_pdb, "--tabular"],
            capture_output=True,
            text=True,
            check=True,
        )
    except FileNotFoundError:
        return None

    return _parse_dssr_base_pairs(result.stdout)


def evaluate_tertiary_motifs(
    model_pdb: str,
    rna_motif_executable: str = "rna_motif",
) -> Optional[dict]:
    """Identify tertiary RNA motifs using Rosetta *rna_motif*.

    Motif types searched for by the CASP16 assessment include A-Minor
    interactions, T-Loops, and GNRA Tetraloops.

    Parameters
    ----------
    model_pdb:
        Path to the PDB model to analyse.
    rna_motif_executable:
        Name / full path of the Rosetta ``rna_motif`` binary.

    Returns
    -------
    Dictionary mapping motif type → list of residue-level descriptions, or
    ``None`` when Rosetta is unavailable.

    Raises
    ------
    subprocess.CalledProcessError
        When ``rna_motif`` exits with a non-zero status.
    """
    try:
        result = subprocess.run(
            [rna_motif_executable, "-in:file:s", model_pdb],
            capture_output=True,
            text=True,
            check=True,
        )
    except FileNotFoundError:
        return None

    return _parse_rna_motif_output(result.stdout)


def compute_f1_score(
    predicted_pairs: List[BasePair],
    target_pairs: List[BasePair],
) -> float:
    """Compute the F1-score for base-pair prediction.

    A predicted pair counts as a true positive only when **both** residue
    identifiers match exactly a pair in *target_pairs*.

    Parameters
    ----------
    predicted_pairs:
        Base pairs extracted from the predicted model.
    target_pairs:
        Base pairs from the experimental reference structure.

    Returns
    -------
    F1-score in the range [0, 1].  Returns 0.0 when either list is empty.

    Examples
    --------
    >>> pred = [frozenset({"A.G1", "A.C10"}), frozenset({"A.U2", "A.A9"})]
    >>> ref  = [frozenset({"A.G1", "A.C10"}), frozenset({"A.C3", "A.G8"})]
    >>> compute_f1_score(pred, ref)
    0.5
    """
    if not predicted_pairs or not target_pairs:
        return 0.0

    predicted_set: Set[BasePair] = set(predicted_pairs)
    target_set: Set[BasePair] = set(target_pairs)

    true_positives = len(predicted_set & target_set)
    precision = true_positives / len(predicted_set)
    recall = true_positives / len(target_set)

    if precision + recall == 0.0:
        return 0.0
    return 2.0 * precision * recall / (precision + recall)


# ---------------------------------------------------------------------------
# Internal parsers
# ---------------------------------------------------------------------------

def _parse_dssr_base_pairs(stdout: str) -> List[BasePair]:
    """Parse the tabular DSSR output and return a list of base pairs.

    Each non-comment, non-header line in the base-pair section contains the
    two residue IDs in the first two whitespace-separated columns.
    """
    pairs: List[BasePair] = []
    in_bp_section = False

    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue

        # Check for "base pair(s)" section header BEFORE skipping comments so
        # that headers like "# base pairs" are still recognised.
        if re.search(r"base.pair", line, re.IGNORECASE):
            in_bp_section = True
            continue

        # Skip remaining comment / decoration lines that are not data rows.
        if line.startswith("#") or line.startswith("*") or line.startswith("-"):
            if in_bp_section:
                # A comment line inside the BP section ends it.
                in_bp_section = False
            continue

        if in_bp_section:
            parts = line.split()
            if len(parts) >= 2:
                pairs.append(frozenset({parts[0], parts[1]}))

    return pairs


def _parse_rna_motif_output(stdout: str) -> dict:
    """Parse Rosetta rna_motif stdout into a motif-type → descriptions dict."""
    motifs: dict = {}
    current_motif: Optional[str] = None

    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        # Rosetta prints motif type headers like "A_MINOR_MOTIF ..."
        motif_match = re.match(
            r"(A.MINOR|T.LOOP|GNRA.TETRALOOP|LOOP_E|KINK.TURN)", line, re.IGNORECASE
        )
        if motif_match:
            current_motif = motif_match.group(1).upper()
            motifs.setdefault(current_motif, [])
        elif current_motif and line:
            motifs[current_motif].append(line)

    return motifs
