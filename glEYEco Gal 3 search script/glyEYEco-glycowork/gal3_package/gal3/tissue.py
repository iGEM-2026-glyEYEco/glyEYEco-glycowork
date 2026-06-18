"""
gal3.tissue
-----------
Section 2 of the notebook: validate the target tissue glycan landscape.

Public API
~~~~~~~~~~
get_ocular_glycans(tissue_terms)  -> pd.DataFrame
get_glycan_type_counts(df)        -> pd.Series
get_core1_glycans(df)             -> pd.DataFrame
get_beta_galactosides(df)         -> pd.DataFrame
annotate_motifs(glycan_list)      -> pd.DataFrame
"""

from __future__ import annotations
import pandas as pd
from .utils import col_contains_any, load_data

# ---------------------------------------------------------------------------
# Default tissue search terms (can be overridden by the caller)
# ---------------------------------------------------------------------------
OCULAR_TISSUE_TERMS = [
    "eye", "cornea", "tears", "ocular",
    "lacrimal", "conjunctiva", "lens",
]

CARDIAC_TISSUE_TERMS = [
    "heart", "cardiac", "myocardium", "cardiomyocyte",
    "atrium", "ventricle", "pericardium", "endocardium", "myocardial",
]

MUSCLE_TISSUE_TERMS = [
    "muscle", "skeletal muscle", "myofiber", "myoblast", "myotube",
    "smooth muscle", "striated", "sarcomere", "myositis",
]

# Motif patterns used in the notebook
CORE1_PATTERN = r"Gal\(b1-3\)GalNAc"
BETA_GAL_PATTERN = r"Gal\(b"


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

def get_tissue_glycans(
    tissue_terms: list[str] | None = None,
) -> pd.DataFrame:
    """Filter df_glycan to rows whose tissue_sample contains any of *tissue_terms*.

    Parameters
    ----------
    tissue_terms : list[str] | None
        Defaults to OCULAR_TISSUE_TERMS.

    Returns
    -------
    pd.DataFrame  — filtered subset of df_glycan, reset index.
    """
    if tissue_terms is None:
        tissue_terms = OCULAR_TISSUE_TERMS
    df_glycan, _ = load_data()
    mask = df_glycan["tissue_sample"].apply(lambda x: col_contains_any(x, tissue_terms))
    return df_glycan[mask].copy().reset_index(drop=True)


# Keep the old name as an alias for backwards compatibility
get_ocular_glycans = get_tissue_glycans


def get_glycan_type_counts(df: pd.DataFrame) -> pd.Series:
    """Count occurrences of each glycan_type in *df*.

    Returns
    -------
    pd.Series  indexed by glycan_type, values = count.
    """
    return df["glycan_type"].value_counts(dropna=False)


def get_core1_glycans(df: pd.DataFrame) -> pd.DataFrame:
    """Return rows of *df* whose glycan string contains the Core 1 O-glycan motif.

    Core 1 pattern: Gal(b1-3)GalNAc
    """
    mask = df["glycan"].str.contains(CORE1_PATTERN, regex=True, na=False)
    return df[mask].reset_index(drop=True)


def get_beta_galactosides(df: pd.DataFrame) -> pd.DataFrame:
    """Return rows of *df* whose glycan string contains any β-galactoside linkage."""
    mask = df["glycan"].str.contains(BETA_GAL_PATTERN, regex=True, na=False)
    return df[mask].reset_index(drop=True)


def annotate_motifs(
    glycan_list: list[str],
    feature_set: list[str] | None = None,
) -> pd.DataFrame:
    """Annotate a list of IUPAC glycan strings with structural motifs.

    Parameters
    ----------
    glycan_list : list[str]
        IUPAC-condensed glycan strings.
    feature_set : list[str] | None
        Passed directly to glycowork.motif.annotate.annotate_dataset.
        Defaults to ['known', 'terminal'].

    Returns
    -------
    pd.DataFrame  — motif presence matrix (rows = glycans, cols = motifs).
    """
    from glycowork.motif.annotate import annotate_dataset  # type: ignore
    if feature_set is None:
        feature_set = ["known", "terminal"]
    return annotate_dataset(glycan_list, feature_set=feature_set, condense=True)