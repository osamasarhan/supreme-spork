# CASP16 Nucleic Acid Evaluation Pipeline

A Python toolkit for replicating the CASP16 assessment of nucleic acid
structure prediction, as described in the accompanying paper.

Source code, supporting information tables, and data are openly available in
the dedicated DasLab repository:
[https://github.com/DasLab/CASP16_NA](https://github.com/DasLab/CASP16_NA)

Predicted models, reference structures, and abstracts can be downloaded from
the [CASP16 prediction centre](https://predictioncenter.org/download_area/CASP16/predictions/).

---

## Package structure

```
casp16_eval/
├── __init__.py             Public API
├── alignment.py            US-align wrapper (template identification)
├── metrics.py              3D structure metrics (TM-score, GDT_TS, lDDT, INF)
├── ranking.py              Z-score ranking for NA monomers and multimers
├── secondary_structure.py  Base-pair / motif analysis (DSSR, Rosetta rna_motif)
└── interface.py            Multimer interface evaluation (ICS, IPS, DockQ, AnAnaS)
```

---

## Installation

```bash
pip install .          # production
pip install ".[dev]"   # with test dependencies
```

### External tools

The following third-party programs must be available on `PATH` for their
respective modules to function:

| Tool | Used for | Version used in CASP16 |
|------|----------|----------------------|
| [US-align](https://zhanggroup.org/US-align/) | TM-score, template search | 20220511 |
| LGA | GDT_TS (C4′ atoms) | — |
| [OpenStructure](https://openstructure.org/) (`ost`) | lDDT, ICS, IPS, i-lDDT, DockQ | — |
| [Tarna-tools / ClaRNA](https://github.com/mmagnus/RNA-tools) | INF | — |
| [DSSR](https://x3dna.org/) (`x3dna-dssr`) | Base-pair extraction | v1.9.9-2020feb06 |
| Rosetta `rna_motif` | Tertiary motif identification | — |
| [AnAnaS](https://team.inria.fr/nano-d/software/ananas/) | Symmetry detection | — |

---

## Quick start

```python
from casp16_eval import (
    run_us_align,
    compute_global_metrics,
    calculate_z_scores,
    compute_na_monomer_ranking,
    compute_multimer_ranking,
    evaluate_secondary_structure,
    compute_f1_score,
    compute_interface_metrics,
    check_symmetry,
)

# 1. Template identification
alignment = run_us_align("prediction.pdb", "template.pdb", mol="RNA")
print(f"TM-score: {alignment['tm_score_1']:.4f}")

# 2. Global 3D metrics
metrics = compute_global_metrics("prediction.pdb", "reference.pdb")
print(metrics)

# 3. Z-score ranking (NA monomers)
tm_scores   = [0.95, 0.80, 0.70, 0.60, 0.10]
z_tm        = calculate_z_scores(tm_scores)
gdt_scores  = [0.90, 0.75, 0.65, 0.55, 0.08]
z_gdt       = calculate_z_scores(gdt_scores)
lddt_scores = [0.92, 0.78, 0.68, 0.58, 0.09]
z_lddt      = calculate_z_scores(lddt_scores)

for i, (ztm, zgdt, zl) in enumerate(zip(z_tm, z_gdt, z_lddt)):
    rank_score = compute_na_monomer_ranking(ztm, zgdt, zl)
    print(f"Group {i+1}: {rank_score:.4f}")

# 4. Secondary structure F1-score
pred_pairs = evaluate_secondary_structure("prediction.pdb")
ref_pairs  = evaluate_secondary_structure("reference.pdb")
if pred_pairs and ref_pairs:
    f1 = compute_f1_score(pred_pairs, ref_pairs)
    print(f"Base-pair F1: {f1:.4f}")

# 5. Multimer interface metrics
interface = compute_interface_metrics("multimer_pred.pdb", "multimer_ref.pdb")
print(interface)

# 6. Symmetry detection
sym = check_symmetry("multimer_pred.pdb", rmsd_cutoff=10.0)
if sym:
    print(f"Symmetry: {sym['symmetry']} (order {sym['order']})")
```

---

## Ranking formulae

### NA Monomers

```
Z = 0.3·Z_TM + 0.3·Z_GDT + 0.4·Z_lDDT      (negative Z clamped to 0)
```

### RNA Multimers and NA-Protein targets

```
Z_global    = 0.3·Z_TM + 0.3·Z_GDT + 0.4·Z_lDDT
Z_interface = ⅓·Z_ICS + ⅓·Z_IPS + ⅓·Z_i-lDDT
Z           = 0.3·Z_global + 0.7·Z_interface  (negative Z clamped to 0)
```

Outliers (scores ≥ 2 standard deviations below the mean) are excluded when
computing the reference mean and standard deviation used for Z-score
normalisation.

---

## Running the tests

```bash
pytest
```
