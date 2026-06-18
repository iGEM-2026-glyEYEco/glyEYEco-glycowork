"""
example_usage.py
----------------
Demonstrates how to call the gal3 package from a plain Python script
(e.g. the backend of a web interface).

Run with:
    python example_usage.py
"""

from gal3 import pipeline
from gal3 import tissue, binding, cancer, cross_tissue, docking

# ─────────────────────────────────────────────────────────
# Option A — run everything at once
# ─────────────────────────────────────────────────────────
print("=== Running full pipeline ===")
results = pipeline.run_full()

print(f"\nOcular glycans found      : {len(results.ocular_df)}")
print(f"Core 1 (T-antigen) hits   : {len(results.core1_df)}")
print(f"β-galactoside hits        : {len(results.beta_gal_df)}")
print(f"\nGal-3 sequence length     : {len(results.gal3_sequence)} aa")
print(f"\nOff-target risk table:\n{results.offtarget_risk}\n")
print(f"Cancer glycans identified : {len(results.cancer_df)}")
print(f"\nTop differential motifs (cancer vs baseline):\n{results.differential_expr}\n")
print(f"Glycan set overlap (eye vs cardiac+muscle):")
for k, v in results.set_overlap.items():
    print(f"  {k}: {len(v)}")

# ─────────────────────────────────────────────────────────
# Option B — call individual modules directly
# ─────────────────────────────────────────────────────────

# Custom tissue terms
cornea_df = tissue.get_tissue_glycans(["cornea", "corneal epithelium"])
print(f"\nCornea-specific glycans   : {len(cornea_df)}")

# Gal-3 binding profile
scores = binding.get_gal3_scores()
print(f"Top 5 Gal-3 binders:\n{scores.head()}\n")

# Disease signature (dry eye)
dry_eye_sig = cross_tissue.signature_comparison(
    system="Dry eye",
    tissue_terms=["eye", "cornea", "tears", "ocular", "lacrimal", "conjunctiva"],
    signature_motifs=["Oglycan_core1"],
    regex_signature={
        "Tn": r"GalNAc",
        "Core 1": r"Gal\(b1-3\)GalNAc",
    },
)
print(f"Dry-eye significant motifs:\n{dry_eye_sig[['Motif','Log2FC','q']].head()}\n")

# Docking targets — SMILES conversion (requires glycowork[chem])
# Uncomment when glycowork[chem] is installed:
# smiles_df = docking.iupac_to_smiles()
# generated, failed = docking.generate_pdbs(smiles_df, out_dir="my_pdbs")
# print(f"PDB files generated: {len(generated)}, failed: {len(failed)}")