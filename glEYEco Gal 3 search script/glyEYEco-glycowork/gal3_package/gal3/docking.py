"""
gal3.docking
------------
Section 5 of the notebook: convert IUPAC glycan strings to SMILES
and optionally to 3D PDB conformers for AutoDock/AlphaFold pipelines.

Public API
~~~~~~~~~~
get_glytoucan_ids(targets)          -> pd.DataFrame
iupac_to_smiles(targets)            -> pd.DataFrame  (requires glycowork[chem])
generate_pdbs(smiles_df, out_dir)   -> tuple[list, list]  (requires rdkit)
"""

from __future__ import annotations
from pathlib import Path
import pandas as pd
from .utils import load_data

# Default set of glycans used for docking in the notebook
DEFAULT_DOCKING_TARGETS: list[str] = [
    "Gal(b1-3)GalNAc",                                      # ON-TARGET : Core 1
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc",                          # ON-TARGET : Sialylated Core 1
    "Gal(b1-4)GlcNAc",                                      # OFF-TARGET: LacNAc
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc",                          # OFF-TARGET: Sialyl-LacNAc
    (
        "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)"
        "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]"
        "Man(b1-4)GlcNAc(b1-4)GlcNAc"
    ),                                                        # OFF-TARGET: bi-ant. N-glycan
    "GalNAc",                                                # NEG CTRL  : Tn Antigen
]


def get_glytoucan_ids(
    targets: list[str] | None = None,
) -> pd.DataFrame:
    """Look up GlyTouCan IDs for a list of glycan IUPAC strings.

    Parameters
    ----------
    targets : list[str] | None
        Defaults to DEFAULT_DOCKING_TARGETS.

    Returns
    -------
    pd.DataFrame  — columns: glycan, glytoucan_id.
    """
    if targets is None:
        targets = DEFAULT_DOCKING_TARGETS
    df_glycan, _ = load_data()
    result = (
        df_glycan[df_glycan["glycan"].isin(targets)][["glycan", "glytoucan_id"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    return result


def iupac_to_smiles(
    targets: list[str] | None = None,
) -> pd.DataFrame:
    """Convert IUPAC glycan strings to SMILES using glycowork[chem].

    Requires:  pip install "glycowork[chem]"

    Parameters
    ----------
    targets : list[str] | None
        Defaults to DEFAULT_DOCKING_TARGETS.

    Returns
    -------
    pd.DataFrame  — columns: glycan, SMILES.

    Raises
    ------
    ImportError if glycowork[chem] is not installed.
    """
    from glycowork.motif.annotate import IUPAC_to_SMILES  # type: ignore

    if targets is None:
        targets = DEFAULT_DOCKING_TARGETS
    smiles_list = IUPAC_to_SMILES(targets)
    return pd.DataFrame({"glycan": targets, "SMILES": smiles_list})


def generate_pdbs(
    smiles_df: pd.DataFrame,
    out_dir: str | Path = "glycan_pdbs",
) -> tuple[list[tuple[str, str]], list[str]]:
    """Generate energy-minimised 3D PDB files from a SMILES DataFrame.

    Requires:  pip install rdkit

    Parameters
    ----------
    smiles_df : pd.DataFrame
        Output of iupac_to_smiles() — must have columns 'glycan' and 'SMILES'.
    out_dir : str | Path
        Directory to write PDB files into. Created if it does not exist.

    Returns
    -------
    (generated, failed)
        generated : list of (glycan_iupac, pdb_path) tuples.
        failed    : list of glycan_iupac strings that could not be processed.
    """
    from rdkit import Chem  # type: ignore
    from rdkit.Chem import AllChem, rdmolfiles  # type: ignore

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    generated: list[tuple[str, str]] = []
    failed: list[str] = []

    for _, row in smiles_df.iterrows():
        glycan: str = row["glycan"]
        smiles: str | None = row["SMILES"]

        if not smiles or str(smiles).strip() == "":
            failed.append(glycan)
            continue
        try:
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is None:
                failed.append(glycan)
                continue
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            AllChem.MMFFOptimizeMolecule(mol)
            safe_name = (
                glycan[:40]
                .replace("(", "").replace(")", "")
                .replace("[", "").replace("]", "")
                .replace(" ", "_")
            )
            pdb_path = str(out / f"{safe_name}.pdb")
            rdmolfiles.MolToPDBFile(mol, pdb_path)
            generated.append((glycan, pdb_path))
        except Exception as exc:
            print(f"Warning: could not generate 3D structure for {glycan[:50]}: {exc}")
            failed.append(glycan)

    return generated, failed