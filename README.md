# glyEYEco-glycowork

Computational validation pipeline for the **glEYEco** project using the [glycowork](https://github.com/BojarLab/glycowork) package.

glEYEco uses the glycan-binding domain (GBD) of **Galectin-3** to anchor a therapeutic hydrogel protein to the hyaluronic-acid-rich mucosal environment of the ocular surface, targeting mucins and β-galactoside-containing glycans on corneal epithelial cells.

## What this repository provides

A self-contained Jupyter notebook (`glEYEco_glycowork_analysis.ipynb`) that:

1. **Validates the target tissue** — mines the `df_glycan` dataset (~50,500 glycan sequences with ~20,000 tissue associations) to map the exact glycosylation patterns present at the ocular surface, including disease-associated changes.
2. **Confirms Galectin-3 specificity** — filters the `glycan_binding` dataset (>790,000 protein–glycan interactions) for Galectin-3 binding partners, cross-references them with ocular tissue glycans, and confirms that β-galactoside motifs (LacNAc, T-antigen) are abundant targets.
3. **Predicts fusion-protein binding with LectinOracle** — runs the engineered modular fusion protein through `LectinOracle_flex` to verify it retains high-affinity specificity for corneal epithelial glycans.
4. **Generates 3D structures for docking** — converts IUPAC-condensed glycan sequences to SMILES and then to energy-minimised 3D `.pdb` files ready for AlphaFold / AutoDock docking simulations.

## Quick start

```bash
pip install glycowork
# Optional: deep-learning models
pip install "glycowork[ml]"
# Optional: SMILES / 3D structure generation
pip install "glycowork[chem]" rdkit
```

Open the notebook:

```bash
jupyter notebook glEYEco_glycowork_analysis.ipynb
```

## References

- glycowork: [Thomes et al., 2021](https://academic.oup.com/glycob/advance-article/doi/10.1093/glycob/cwab067/6311240) | [GitHub](https://github.com/BojarLab/glycowork)
- LectinOracle: [Lundstrom et al., 2021](https://onlinelibrary.wiley.com/doi/10.1002/advs.202103807)
- GLYCAM-Web 3D builder: https://glycam.org/cb
- GlyTouCan glycan registry: https://glytoucan.org
