"""
gal3.cancer
-----------
Section 6 of the notebook: extract disease/cancer glycan signatures and
run differential glycomics.

Public API
~~~~~~~~~~
get_cancer_glycans()                    -> tuple[pd.DataFrame, pd.DataFrame]
differential_expression()               -> pd.DataFrame
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from .utils import load_data

CANCER_KEYWORDS = [
    "cancer", "carcinoma", "tumor", "melanoma",
    "adenocarcinoma", "malignan",
]

FAKE_NULLS = {"nan", "none", "", "[]", "null"}


def get_cancer_glycans() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mine df_glycan for disease-associated and cancer-specific glycans.

    Returns
    -------
    (all_disease_df, cancer_df)
        all_disease_df : all rows with a non-null disease_association.
        cancer_df      : subset whose disease_association matches cancer terms.
    """
    df_glycan, _ = load_data()

    col = "disease_association" if "disease_association" in df_glycan.columns else "disease"
    df = df_glycan.copy()
    df[col] = df[col].astype(str)

    diseased = df[~df[col].str.lower().isin(FAKE_NULLS)].copy()

    pattern = "|".join(CANCER_KEYWORDS)
    cancer_mask = diseased[col].str.contains(pattern, case=False, na=False)
    cancer_df = diseased[cancer_mask].copy()

    return diseased, cancer_df


def differential_expression(
    p_threshold: float = 0.05,
    top_n: int = 10,
) -> pd.DataFrame:
    """Compute Log2FC + Welch's t-test comparing cancer vs non-cancer glycans.

    Parameters
    ----------
    p_threshold : float
        p-value cutoff for significance filtering (default 0.05).
    top_n : int
        Number of top enriched motifs to return (default 10).

    Returns
    -------
    pd.DataFrame  — columns: Motif, Log2FC, pval.
                    Filtered to p < p_threshold, sorted by Log2FC descending.
    """
    from scipy.stats import ttest_ind  # type: ignore
    from glycowork.motif.annotate import annotate_dataset  # type: ignore

    df_glycan, _ = load_data()
    col = "disease_association" if "disease_association" in df_glycan.columns else "disease"
    df = df_glycan.copy()
    df[col] = df[col].astype(str)
    clean = df[~df[col].str.lower().isin(FAKE_NULLS)].copy()

    pattern = "|".join(CANCER_KEYWORDS)
    cancer_mask = clean[col].str.contains(pattern, case=False, na=False)
    cancer_glycans = clean[cancer_mask]["glycan"].tolist()
    baseline_glycans = clean[~cancer_mask]["glycan"].tolist()

    combined = cancer_glycans + baseline_glycans
    motif_matrix = annotate_dataset(combined, feature_set="known")

    cancer_m = motif_matrix.iloc[: len(cancer_glycans)]
    baseline_m = motif_matrix.iloc[len(cancer_glycans) :]

    results = []
    for motif in motif_matrix.columns:
        d, h = cancer_m[motif], baseline_m[motif]
        lfc = np.log2((d.mean() + 1e-5) / (h.mean() + 1e-5))
        pval = (
            ttest_ind(d, h, equal_var=False)[1]
            if (d.nunique() > 1 or h.nunique() > 1)
            else 1.0
        )
        results.append({"Motif": motif, "Log2FC": lfc, "pval": pval})

    df_diff = pd.DataFrame(results)
    significant = df_diff[(df_diff["pval"] < p_threshold) & (df_diff["Log2FC"] > 0)]
    return significant.sort_values("Log2FC", ascending=False).head(top_n).reset_index(drop=True)