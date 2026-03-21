"""
alignment.py - Template identification using US-align.

The CASP16 assessment used US-align (Version 20220511) with the ``-mol RNA``
flag to perform sequence-independent structural alignment against PDB
templates released before November 11, 2024.

Reference
---------
https://github.com/DasLab/CASP16_NA
https://github.com/daslab/daslab_tools/
"""

import re
import subprocess
from typing import Optional


def run_us_align(
    target_pdb: str,
    template_pdb: str,
    mol: str = "RNA",
    extra_args: Optional[list] = None,
) -> dict:
    """Run US-align between *target_pdb* and *template_pdb*.

    Parameters
    ----------
    target_pdb:
        Path to the query / predicted structure (PDB or mmCIF format).
    template_pdb:
        Path to the reference / template structure (PDB or mmCIF format).
    mol:
        Molecule type passed to ``-mol``.  Use ``"RNA"`` (default) for RNA
        targets or ``"DNA"`` for DNA targets.
    extra_args:
        Any additional command-line arguments forwarded verbatim to US-align.

    Returns
    -------
    dict with keys:
        ``tm_score_1``   – TM-score normalised by length of *target_pdb*
        ``tm_score_2``   – TM-score normalised by length of *template_pdb*
        ``rmsd``         – RMSD over aligned residues (Å)
        ``seq_id``       – Sequence identity fraction
        ``aligned_len``  – Number of aligned residue pairs
        ``stdout``       – Raw US-align standard output

    Raises
    ------
    FileNotFoundError
        When the ``USalign`` executable cannot be found on PATH.
    subprocess.CalledProcessError
        When US-align exits with a non-zero status code.
    """
    cmd = ["USalign", target_pdb, template_pdb, "-mol", mol]
    if extra_args:
        cmd.extend(extra_args)

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return _parse_us_align_output(result.stdout)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_TM1_RE = re.compile(r"TM-score=\s*([\d.]+)\s*\(if normalized by length of Structure_1")
_TM2_RE = re.compile(r"TM-score=\s*([\d.]+)\s*\(if normalized by length of Structure_2")
_RMSD_RE = re.compile(r"RMSD=\s*([\d.]+),\s*Seq_ID=n_identical/n_aligned=\s*([\d.]+)")
_ALIGNED_RE = re.compile(r"Number of residues in common=\s*(\d+)")


def _parse_us_align_output(stdout: str) -> dict:
    """Parse raw US-align stdout into a structured dictionary."""
    result: dict = {"stdout": stdout}

    m = _TM1_RE.search(stdout)
    result["tm_score_1"] = float(m.group(1)) if m else None

    m = _TM2_RE.search(stdout)
    result["tm_score_2"] = float(m.group(1)) if m else None

    m = _RMSD_RE.search(stdout)
    if m:
        result["rmsd"] = float(m.group(1))
        result["seq_id"] = float(m.group(2))
    else:
        result["rmsd"] = None
        result["seq_id"] = None

    m = _ALIGNED_RE.search(stdout)
    result["aligned_len"] = int(m.group(1)) if m else None

    return result
