"""
interface.py - Multimer interface quality and symmetry evaluation.

Interface metrics (CASP16)
--------------------------
* **ICS** – Interface Contact Score; F1 of interface contacts.
* **IPS** – Interface Patch Score; Jaccard coefficient of interface patches.
* **i-lDDT** – lDDT restricted to interface residues.
* **DockQ**  – Combined docking quality score.

All interface metrics are computed via **OpenStructure**.

Symmetry detection
------------------
Symmetry is detected via **AnAnaS** (Applied Numerical Algorithms for
Nuclear and Atomic Structures) with a 10 Å RMSD cutoff.  Because AnAnaS
operates on protein C-alpha atoms, RNA models are pre-processed by:

1. Renaming every nucleotide residue to ``ALA``.
2. Renaming the phosphorus atom (``P``) to ``CA``.

The highest-order symmetry group found within the RMSD cutoff is reported.
"""

import re
import subprocess
from typing import Optional


def compute_interface_metrics(
    prediction: str,
    reference: str,
    ost_executable: str = "ost",
) -> dict:
    """Compute ICS, IPS, i-lDDT, and DockQ for a multimer prediction.

    Parameters
    ----------
    prediction:
        Path to the predicted multimer model (PDB format).
    reference:
        Path to the experimental reference structure (PDB format).
    ost_executable:
        Name / full path of the OpenStructure ``ost`` binary.

    Returns
    -------
    dict with keys ``ics``, ``ips``, ``i_lddt``, ``dockq``.
    Values are ``None`` when OpenStructure is unavailable.
    """
    script = _build_ost_interface_script(prediction, reference)
    try:
        result = subprocess.run(
            [ost_executable, "--script", "/dev/stdin"],
            input=script,
            capture_output=True,
            text=True,
            check=True,
        )
        return _parse_ost_interface_output(result.stdout)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return {"ics": None, "ips": None, "i_lddt": None, "dockq": None}


def check_symmetry(
    model_pdb: str,
    rmsd_cutoff: float = 10.0,
    ananas_executable: str = "AnAnaS",
    tmp_dir: str = "/tmp",
) -> Optional[dict]:
    """Detect the highest-order symmetry group in *model_pdb* using AnAnaS.

    Before running AnAnaS the RNA model is converted so that:
    - every nucleotide residue type is replaced with ``ALA``
    - the phosphorus atom name ``P`` is replaced with ``CA``

    Parameters
    ----------
    model_pdb:
        Path to the RNA model (PDB format).
    rmsd_cutoff:
        Maximum RMSD (Å) for accepting a symmetry solution.  The CASP16
        authors used 10 Å.
    ananas_executable:
        Name / full path of the AnAnaS binary.
    tmp_dir:
        Directory for the temporary pre-processed PDB file.

    Returns
    -------
    dict with keys:
        ``symmetry`` – symmetry group string (e.g. ``"C3"``)
        ``order``    – integer symmetry order
        ``rmsd``     – RMSD of the best symmetry solution (Å)
    or ``None`` when AnAnaS is unavailable or no symmetry is found within the
    RMSD cutoff.
    """
    import os
    import tempfile

    converted_pdb = _convert_rna_to_protein_format(model_pdb, tmp_dir)
    if converted_pdb is None:
        return None

    try:
        result = subprocess.run(
            [ananas_executable, converted_pdb, "-t", str(rmsd_cutoff)],
            capture_output=True,
            text=True,
            check=True,
        )
        symmetry = _parse_ananas_output(result.stdout, rmsd_cutoff)
        return symmetry
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    finally:
        if converted_pdb and os.path.exists(converted_pdb):
            os.unlink(converted_pdb)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_ost_interface_script(prediction: str, reference: str) -> str:
    """Return a Python script (for ``ost --script``) that prints ICS, IPS, i-lDDT, DockQ."""
    return (
        "from ost.mol.alg import dockq as dq, scoring\n"
        "import ost.io as io\n"
        f"pred = io.LoadPDB('{prediction}')\n"
        f"ref  = io.LoadPDB('{reference}')\n"
        "result = dq.DockQ(pred, ref)\n"
        "print('ICS:', result.ics)\n"
        "print('IPS:', result.ips)\n"
        "print('i_lDDT:', result.i_lddt)\n"
        "print('DockQ:', result.dockq)\n"
    )


def _parse_ost_interface_output(stdout: str) -> dict:
    """Parse the metric values printed by the OpenStructure script."""
    def extract(key: str) -> Optional[float]:
        m = re.search(rf"{key}:\s*([\d.]+)", stdout)
        return float(m.group(1)) if m else None

    return {
        "ics": extract("ICS"),
        "ips": extract("IPS"),
        "i_lddt": extract("i_lDDT"),
        "dockq": extract("DockQ"),
    }


def _convert_rna_to_protein_format(model_pdb: str, tmp_dir: str) -> Optional[str]:
    """Write a temporary PDB where nucleotides → ALA and P atoms → CA.

    Returns the path of the temporary file, or ``None`` on I/O error.
    """
    import tempfile

    _RNA_RESIDUES = frozenset({"A", "C", "G", "U", "DA", "DC", "DG", "DT"})

    try:
        with open(model_pdb) as fh:
            lines = fh.readlines()
    except OSError:
        return None

    converted_lines = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            residue_name = line[17:20].strip()
            atom_name = line[12:16].strip()
            if residue_name in _RNA_RESIDUES:
                # Replace residue name with ALA (columns 17-19, 3 chars)
                line = line[:17] + "ALA" + line[20:]
                # Replace phosphorus with CA (columns 12-15, 4 chars padded)
                if atom_name == "P":
                    line = line[:12] + " CA " + line[16:]
        converted_lines.append(line)

    import os
    fd, tmp_path = tempfile.mkstemp(suffix=".pdb", dir=tmp_dir)
    try:
        with os.fdopen(fd, "w") as fh:
            fh.writelines(converted_lines)
    except OSError:
        return None

    return tmp_path


def _parse_ananas_output(stdout: str, rmsd_cutoff: float) -> Optional[dict]:
    """Parse AnAnaS stdout and return the highest-order symmetry within cutoff."""
    best: Optional[dict] = None

    for line in stdout.splitlines():
        # Example line: "C3   order=3   RMSD=1.23 Å"
        m = re.search(
            r"(C\d+|D\d+|T|O|I)\s+.*?order[=\s]+(\d+).*?RMSD[=\s]+([\d.]+)",
            line,
            re.IGNORECASE,
        )
        if not m:
            continue
        sym = m.group(1).upper()
        order = int(m.group(2))
        rmsd = float(m.group(3))
        if rmsd <= rmsd_cutoff:
            if best is None or order > best["order"]:
                best = {"symmetry": sym, "order": order, "rmsd": rmsd}

    return best
