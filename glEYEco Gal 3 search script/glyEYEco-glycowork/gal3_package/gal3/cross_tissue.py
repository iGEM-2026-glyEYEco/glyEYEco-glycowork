"""
gal3.cross_tissue
-----------------
Sections 7 & 9 of the notebook: compare glycan repertoires across
ocular, cardiac, and skeletal-muscle tissues and run
disease-signature comparisons.

Public API
~~~~~~~~~~
get_cardiac_glycans()       -> pd.DataFrame
get_muscle_glycans()        -> pd.DataFrame
glycan_set_overlap(...)     -> dict[str, set[str]]
gal3_motif_prevalence(...)  -> pd.DataFrame
signature_comparison(...)   -> pd.DataFrame
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from .utils import col_contains_any, load_data
from .tissue import CARDIAC_TISSUE_TERMS, MUSCLE_TISSUE_TERMS, OCULAR_TISSUE_TERMS

GAL3_MOTIFS: dict[str, str] = {
    "LacNAc [Gal(b1-4)GlcNAc]":           r"Gal\(b1-4\)GlcNAc",
    "Core-1 / T-antigen [Gal(b1-3)GalNAc]": r"Gal\(b1-3\)GalNAc",
    "Poly-LacNAc repeat":                  r"Gal\(b1-4\)GlcNAc.*Gal\(b1-4\)GlcNAc",
    "Sialyl-LacNAc [Neu5Ac-Gal(b1-4)]":   r"Neu5Ac.*Gal\(b1-4\)GlcNAc",
}


def _filter_tissue(terms: list[str]) -> pd.DataFrame:
    df_glycan, _ = load_data()
    mask = df_glycan["tissue_sample"].apply(lambda x: col_contains_any(x, terms))
    return df_glycan[mask].copy().reset_index(drop=True)


def get_cardiac_glycans() -> pd.DataFrame:
    """Return df_glycan rows associated with cardiac tissue."""
    return _filter_tissue(CARDIAC_TISSUE_TERMS)


def get_muscle_glycans() -> pd.DataFrame:
    """Return df_glycan rows associated with skeletal muscle tissue."""
    return _filter_tissue(MUSCLE_TISSUE_TERMS)


def glycan_set_overlap(
    ocular_df: pd.DataFrame,
    cardiac_df: pd.DataFrame,
    muscle_df: pd.DataFrame,
) -> dict[str, set[str]]:
    """Compute set-algebra overlaps across the three tissue glycomes.

    Returns
    -------
    dict with keys:
        eye_only, cardiac_only, muscle_only, shared_all,
        shared_eye_cardiac, shared_eye_muscle, shared_cardiac_muscle.
    """
    e = set(ocular_df["glycan"])
    c = set(cardiac_df["glycan"])
    m = set(muscle_df["glycan"])
    cm = c | m

    return {
        "eye_only":            e - cm,
        "cardiac_only":        c - e - m,
        "muscle_only":         m - e - c,
        "shared_all":          e & cm,
        "shared_eye_cardiac":  e & c,
        "shared_eye_muscle":   e & m,
        "shared_cardiac_muscle": c & m,
    }


def gal3_motif_prevalence(
    tissue_map: dict[str, pd.DataFrame],
    motifs: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Count how many glycans in each tissue carry each Galectin-3 binding motif.

    Parameters
    ----------
    tissue_map : dict[str, pd.DataFrame]
        e.g. {"Eye": ocular_df, "Cardiac": cardiac_df, "Muscle": muscle_df}
    motifs : dict[str, str] | None
        {label: regex_pattern}. Defaults to GAL3_MOTIFS.

    Returns
    -------
    pd.DataFrame — rows = motifs, columns = tissue names, values = counts.
    """
    if motifs is None:
        motifs = GAL3_MOTIFS

    records: dict[str, dict] = {}
    for motif_name, pattern in motifs.items():
        row: dict[str, int] = {}
        for tissue_name, tdf in tissue_map.items():
            row[tissue_name] = int(
                tdf["glycan"].str.contains(pattern, regex=True, na=False).sum()
            )
        records[motif_name] = row

    return pd.DataFrame(records).T


def signature_comparison(
    system: str,
    tissue_terms: list[str],
    signature_motifs: list[str],
    regex_signature: dict[str, str] | None = None,
    p_threshold: float = 0.05,
) -> pd.DataFrame:
    """Split a tissue glycome by a literature glycan signature and contrast the
    two cohorts across the full known-motif vocabulary (Sections 8–11).

    Parameters
    ----------
    system : str
        Human-readable name, e.g. "Dry eye".
    tissue_terms : list[str]
        Used to filter df_glycan to the relevant tissue.
    signature_motifs : list[str]
        glycowork 'known' motif names that define the diseased cohort.
    regex_signature : dict[str, str] | None
        Extra {label: IUPAC_regex} patterns for epitopes not in glycowork
        known motifs (e.g. Tn, sialyl-Tn).
    p_threshold : float
        Significance threshold for the Fisher's exact q-value (default 0.05).

    Returns
    -------
    pd.DataFrame  — per-motif statistics (Log2FC, p-value, FDR q, etc.).
                    Filtered to q < p_threshold.
    """
    from scipy.stats import ttest_ind, fisher_exact  # type: ignore
    from statsmodels.stats.multitest import multipletests  # type: ignore
    from glycowork.motif.annotate import annotate_dataset  # type: ignore

    df_glycan, _ = load_data()
    sub = df_glycan[
        df_glycan["tissue_sample"].apply(lambda x: col_contains_any(x, tissue_terms))
    ].copy()
    glycans = sub["glycan"].tolist()

    M = (annotate_dataset(glycans, feature_set=["known"]) > 0).astype(int)

    mask = pd.Series(False, index=M.index)
    for m in [s for s in signature_motifs if s in M.columns]:
        mask |= M[m] > 0
    if regex_signature:
        gser = pd.Series(glycans, index=M.index)
        for pat in regex_signature.values():
            mask |= gser.str.contains(pat, regex=True, na=False)
    mask = mask.values
    Md, Mb = M[mask], M[~mask]

    records = []
    for motif in M.columns:
        d, h = Md[motif], Mb[motif]
        lfc = np.log2((d.mean() + 1e-3) / (h.mean() + 1e-3))
        t_p = (
            ttest_ind(d, h, equal_var=False)[1]
            if (d.nunique() > 1 or h.nunique() > 1)
            else 1.0
        )
        a, c = int(d.sum()), int(h.sum())
        try:
            orr, f_p = fisher_exact([[a, len(d) - a], [c, len(h) - c]])
        except Exception:
            orr, f_p = np.nan, 1.0
        records.append(
            {
                "Motif": motif,
                "System": system,
                "Diseased_%": round(100 * d.mean(), 1),
                "Baseline_%": round(100 * h.mean(), 1),
                "Log2FC": round(lfc, 2),
                "t_p": t_p,
                "OR": orr,
                "fisher_p": f_p,
            }
        )

    S = pd.DataFrame(records)
    S["q"] = multipletests(S["fisher_p"], method="fdr_bh")[1]
    return (
        S[S["q"] < p_threshold]
        .sort_values("fisher_p")
        .reset_index(drop=True)
    )