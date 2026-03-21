"""
metrics.py - 3D structural quality metrics used in the CASP16 assessment.

The assessment combined four global metrics:

* **TM-score / TM-align** – calculated via US-align
* **GDT_TS**              – calculated via Local-Global Alignment (LGA) using C4′ atoms
* **lDDT** (with steric penalty) – calculated via OpenStructure
* **INF** (Interaction Network Fidelity) – calculated via Tarna-tools / ClaRNA
"""

import subprocess
from typing import Optional

from .alignment import run_us_align


def compute_global_metrics(
    prediction: str,
    reference: str,
    mol: str = "RNA",
    lga_executable: str = "lga",
    openstructure_executable: str = "ost",
    tarna_executable: str = "tarna",
) -> dict:
    """Compute all global 3D structure metrics for a single prediction.

    Parameters
    ----------
    prediction:
        Path to the predicted model (PDB format).
    reference:
        Path to the experimental reference structure (PDB format).
    mol:
        Molecule type for US-align (``"RNA"`` or ``"DNA"``).
    lga_executable:
        Name / full path of the LGA binary.
    openstructure_executable:
        Name / full path of the OpenStructure ``ost`` binary.
    tarna_executable:
        Name / full path of the Tarna-tools binary.

    Returns
    -------
    dict with keys ``tm_score``, ``gdt_ts``, ``lddt``, ``inf``.
    Any metric whose external tool is unavailable will be ``None``.
    """
    metrics: dict = {}

    # TM-score via US-align
    us_result = run_us_align(prediction, reference, mol=mol)
    metrics["tm_score"] = us_result.get("tm_score_1")

    # GDT_TS via LGA (C4' atoms)
    metrics["gdt_ts"] = _run_lga(prediction, reference, lga_executable)

    # lDDT via OpenStructure
    metrics["lddt"] = _run_openstructure_lddt(
        prediction, reference, openstructure_executable
    )

    # INF via Tarna-tools / ClaRNA
    metrics["inf"] = _run_tarna_inf(prediction, reference, tarna_executable)

    return metrics


# ---------------------------------------------------------------------------
# Thin wrappers around external tools
# ---------------------------------------------------------------------------

def _run_lga(prediction: str, reference: str, lga_executable: str) -> Optional[float]:
    """Run LGA and return GDT_TS computed over C4′ atoms."""
    try:
        result = subprocess.run(
            [lga_executable, "-3", "-sda", "-ca", prediction, reference],
            capture_output=True,
            text=True,
            check=True,
        )
        return _parse_lga_gdt_ts(result.stdout)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def _parse_lga_gdt_ts(stdout: str) -> Optional[float]:
    import re
    m = re.search(r"GDT_TS=\s*([\d.]+)", stdout)
    return float(m.group(1)) if m else None


def _run_openstructure_lddt(
    prediction: str, reference: str, ost_executable: str
) -> Optional[float]:
    """Invoke OpenStructure to calculate lDDT with steric penalty."""
    script = (
        "from ost.mol.alg import lDDTScorer; "
        "import ost.io as io; "
        f"pred = io.LoadPDB('{prediction}'); "
        f"ref = io.LoadPDB('{reference}'); "
        "scorer = lDDTScorer(ref); "
        "print('lDDT:', scorer.lDDT(pred, penalize_extra_chains=True))"
    )
    try:
        result = subprocess.run(
            [ost_executable, "--script", "/dev/stdin"],
            input=script,
            capture_output=True,
            text=True,
            check=True,
        )
        return _parse_openstructure_lddt(result.stdout)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def _parse_openstructure_lddt(stdout: str) -> Optional[float]:
    import re
    m = re.search(r"lDDT:\s*([\d.]+)", stdout)
    return float(m.group(1)) if m else None


def _run_tarna_inf(
    prediction: str, reference: str, tarna_executable: str
) -> Optional[float]:
    """Run Tarna-tools / ClaRNA to compute Interaction Network Fidelity (INF)."""
    try:
        result = subprocess.run(
            [tarna_executable, "--pred", prediction, "--ref", reference, "--inf"],
            capture_output=True,
            text=True,
            check=True,
        )
        return _parse_tarna_inf(result.stdout)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def _parse_tarna_inf(stdout: str) -> Optional[float]:
    import re
    m = re.search(r"INF\s*[=:]\s*([\d.]+)", stdout)
    return float(m.group(1)) if m else None
