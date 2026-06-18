"""
gal3.binding
------------
Section 3 of the notebook: confirm Galectin-3 binding specificity.

Public API
~~~~~~~~~~
get_gal3_scores()                        -> pd.DataFrame  (glycan → binding_score)
get_gal3_sequence()                      -> str
motif_affinity_summary(gal3_scores)      -> pd.DataFrame
offtarget_risk_table(ocular_df,
                     gal3_scores)        -> pd.DataFrame
"""

from __future__ import annotations
import pandas as pd
from .utils import load_data

# ---------------------------------------------------------------------------
# Motif definitions — the therapeutic targets and controls described in §3.2
# ---------------------------------------------------------------------------
TARGET_GLYCANS: dict[str, dict] = {
    "Core 1 O-Glycan (on-target)": {
        "iupac":   "Gal(b1-3)GalNAc",
        "pattern": r"Gal\(b1-3\)GalNAc",
        "class":   "on-target",
    },
    "LacNAc (off-target)": {
        "iupac":   "Gal(b1-4)GlcNAc",
        "pattern": r"Gal\(b1-4\)GlcNAc",
        "class":   "off-target",
    },
    "Galact. bi-ant. N-glycan (off-target)": {
        "iupac":   "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
        "pattern": r"Gal\(b1-4\)GlcNAc\(b1-2\)Man",
        "class":   "off-target",
    },
    "Tn Antigen (neg ctrl)": {
        "iupac":   "GalNAc",
        "pattern": r"^GalNAc$",
        "class":   "neg-ctrl",
    },
}


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

def _find_gal3_row(glycan_binding: pd.DataFrame) -> pd.DataFrame:
    """Locate the Galectin-3 row; try exact match first, then regex."""
    rows = glycan_binding[glycan_binding["protein"] == "Gal_3"]
    if rows.empty:
        rows = glycan_binding[
            glycan_binding["protein"].str.contains(
                r"galectin.?3", case=False, regex=True, na=False
            )
        ]
    if rows.empty:
        raise ValueError(
            "Galectin-3 not found in glycan_binding['protein']. "
            "Check the exact protein name in the dataset."
        )
    return rows.iloc[[0]]


def get_gal3_scores() -> pd.DataFrame:
    """Extract and sort the Galectin-3 glycan binding scores from glycan_binding.

    Returns
    -------
    pd.DataFrame with columns ['binding_score'], index = glycan IUPAC string.
    Sorted descending by binding_score.
    """
    _, glycan_binding = load_data()
    glycan_cols = [c for c in glycan_binding.columns if c not in ("protein", "target")]
    gal3_row = _find_gal3_row(glycan_binding)

    scores = (
        gal3_row[glycan_cols]
        .T.dropna()
        .rename(columns={gal3_row.index[0]: "binding_score"})
        .sort_values("binding_score", ascending=False)
    )
    return scores


def get_gal3_sequence() -> str:
    """Return the amino-acid sequence string for Galectin-3 from glycan_binding."""
    _, glycan_binding = load_data()
    gal3_row = _find_gal3_row(glycan_binding)
    return gal3_row["target"].values[0]


def motif_affinity_summary(gal3_scores: pd.DataFrame) -> pd.DataFrame:
    """For each target motif in TARGET_GLYCANS, search glycan_binding column names
    and return mean/max binding scores (§3.3).

    Parameters
    ----------
    gal3_scores : pd.DataFrame
        Output of get_gal3_scores().

    Returns
    -------
    pd.DataFrame  — columns: Glycan, Class, Gal-3 binders with motif,
                    Mean binding score, Max binding score.
    """
    rows = []
    for name, info in TARGET_GLYCANS.items():
        matching = gal3_scores[
            gal3_scores.index.str.contains(info["pattern"], regex=True, na=False)
        ]
        rows.append(
            {
                "Glycan": name,
                "Class": info["class"],
                "Gal-3 binders with motif": len(matching),
                "Mean binding score": (
                    round(matching["binding_score"].mean(), 4) if len(matching) else None
                ),
                "Max binding score": (
                    round(matching["binding_score"].max(), 4) if len(matching) else None
                ),
            }
        )
    return pd.DataFrame(rows)


def offtarget_risk_table(
    ocular_df: pd.DataFrame,
) -> pd.DataFrame:
    """Build the off-target risk table (§3.2): binding score lookup + motif prevalence.

    Parameters
    ----------
    ocular_df : pd.DataFrame
        Tissue-filtered subset of df_glycan (output of tissue.get_ocular_glycans).

    Returns
    -------
    pd.DataFrame — one row per target/off-target.
    """
    df_glycan, _ = load_data()
    gal3_scores = get_gal3_scores()

    rows = []
    for name, info in TARGET_GLYCANS.items():
        exact_score = (
            round(gal3_scores.loc[info["iupac"], "binding_score"], 4)
            if info["iupac"] in gal3_scores.index
            else None
        )
        ocular_n = int(
            ocular_df["glycan"].str.contains(info["pattern"], regex=True, na=False).sum()
        )
        total_n = int(
            df_glycan["glycan"].str.contains(info["pattern"], regex=True, na=False).sum()
        )
        rows.append(
            {
                "Glycan": name,
                "Class": info["class"],
                "Exact DB match": exact_score is not None,
                "Binding score (if found)": exact_score,
                "Ocular glycans w/ motif": ocular_n,
                "Total DB glycans w/ motif": total_n,
            }
        )
    return pd.DataFrame(rows)