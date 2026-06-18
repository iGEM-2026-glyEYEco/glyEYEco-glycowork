# gal3: Galectin-3 Computational Validation Package

Wraps the **glEYEco × glycowork** analysis notebook into an importable
Python package so computations can be triggered programmatically and
eventually exposed through a web interface.

## Installation

```bash
# 1. Clone / copy this folder, then install in editable mode:
pip install -e .

# Optional extras
pip install -e ".[ml]"    # LectinOracle_flex binding predictions
pip install -e ".[chem]"  # IUPAC → SMILES → 3D PDB (requires RDKit)
```

## Quick start

```python
from gal3 import pipeline

# Run every stage in one call
results = pipeline.run_full()

# Or run individual stages
tissue_data  = pipeline.run_tissue()
binding_data = pipeline.run_binding()
cancer_data  = pipeline.run_cancer()
cross_data   = pipeline.run_cross_tissue()
```

## Package layout

```
gal3_package/
├── pyproject.toml          ← installs with pip
├── README.md               ← full API reference
├── example_usage.py        ← starting point
└── gal3/
    ├── __init__.py         ← public surface
    ├── pipeline.py         ← high-level orchestrator  ← start here
    ├── utils.py            ← col_contains_any, cached data loader
    ├── tissue.py           ← Section 2 (ocular glycan validation)
    ├── binding.py          ← Section 3 (Gal-3 specificity)
    ├── cancer.py           ← Section 6 (cancer glycan mining)
    ├── cross_tissue.py     ← Sections 7+9 (cardiac/muscle + signatures)
    └── docking.py          ← Section 5 (IUPAC → SMILES → PDB)
```

## Module reference

### `gal3.tissue`
| Function | Returns | Notes |
|---|---|---|
| `get_tissue_glycans(tissue_terms)` | `pd.DataFrame` | Filter df_glycan by tissue |
| `get_glycan_type_counts(df)` | `pd.Series` | N-glycan vs O-glycan counts |
| `get_core1_glycans(df)` | `pd.DataFrame` | Core 1 T-antigen-containing glycans |
| `get_beta_galactosides(df)` | `pd.DataFrame` | β-galactoside-containing glycans |
| `annotate_motifs(glycan_list)` | `pd.DataFrame` | Motif presence matrix |

### `gal3.binding`
| Function | Returns |
|---|---|
| `get_gal3_scores()` | Sorted binding scores from glycan_binding |
| `get_gal3_sequence()` | Gal-3 amino-acid sequence string |
| `motif_affinity_summary(gal3_scores)` | Mean/max scores per motif class |
| `offtarget_risk_table(ocular_df)` | On-/off-target prevalence + score lookup |

### `gal3.cancer`
| Function | Returns |
|---|---|
| `get_cancer_glycans()` | `(all_disease_df, cancer_df)` |
| `differential_expression()` | Top enriched motifs (Log2FC + p-value) |

### `gal3.cross_tissue`
| Function | Returns |
|---|---|
| `get_cardiac_glycans()` | Cardiac subset of df_glycan |
| `get_muscle_glycans()` | Muscle subset of df_glycan |
| `glycan_set_overlap(e, c, m)` | Set-algebra dict (eye_only, shared_all, …) |
| `gal3_motif_prevalence(tissue_map)` | Motif count matrix across tissues |
| `signature_comparison(...)` | Per-motif stats for a disease signature |

### `gal3.docking`
| Function | Returns | Extra dependency |
|---|---|---|
| `get_glytoucan_ids(targets)` | GlyTouCan ID lookup table | — |
| `iupac_to_smiles(targets)` | SMILES DataFrame | `glycowork[chem]` |
| `generate_pdbs(smiles_df, out_dir)` | `(generated, failed)` lists | `rdkit` |