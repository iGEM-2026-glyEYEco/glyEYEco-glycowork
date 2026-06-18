"""
gal3.pipeline
-------------
High-level orchestrator: runs the full Galectin-3 analysis in one call
or lets you run individual stages.

Usage
-----
>>> from gal3 import pipeline
>>> results = pipeline.run_full()          # everything
>>> results = pipeline.run_tissue()        # just tissue validation
>>> results = pipeline.run_binding()       # just binding specificity
"""

from __future__ import annotations
from dataclasses import dataclass, field
import pandas as pd

from . import tissue, binding, cancer, cross_tissue


@dataclass
class AnalysisResults:
    """Container returned by run_full() — one attribute per pipeline stage."""

    # Stage 1 — tissue
    ocular_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    glycan_type_counts: pd.Series = field(default_factory=pd.Series)
    core1_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    beta_gal_df: pd.DataFrame = field(default_factory=pd.DataFrame)

    # Stage 2 — binding
    gal3_scores: pd.DataFrame = field(default_factory=pd.DataFrame)
    gal3_sequence: str = ""
    motif_affinity: pd.DataFrame = field(default_factory=pd.DataFrame)
    offtarget_risk: pd.DataFrame = field(default_factory=pd.DataFrame)

    # Stage 3 — cancer
    all_disease_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    cancer_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    differential_expr: pd.DataFrame = field(default_factory=pd.DataFrame)

    # Stage 4 — cross-tissue
    cardiac_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    muscle_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    set_overlap: dict = field(default_factory=dict)
    motif_prevalence: pd.DataFrame = field(default_factory=pd.DataFrame)


# ──────────────────────────────────────────────────────────────────────────────
# Individual stage runners
# ──────────────────────────────────────────────────────────────────────────────

def run_tissue(
    tissue_terms: list[str] | None = None,
) -> dict[str, object]:
    """Run Section 2 (tissue validation) and return a results dict.

    Parameters
    ----------
    tissue_terms : list[str] | None
        Custom tissue search terms. Defaults to OCULAR_TISSUE_TERMS.

    Returns
    -------
    dict with keys: ocular_df, glycan_type_counts, core1_df, beta_gal_df.
    """
    ocular_df = tissue.get_tissue_glycans(tissue_terms)
    return {
        "ocular_df":         ocular_df,
        "glycan_type_counts": tissue.get_glycan_type_counts(ocular_df),
        "core1_df":          tissue.get_core1_glycans(ocular_df),
        "beta_gal_df":       tissue.get_beta_galactosides(ocular_df),
    }


def run_binding() -> dict[str, object]:
    """Run Section 3 (binding specificity) and return a results dict.

    Returns
    -------
    dict with keys: gal3_scores, gal3_sequence, motif_affinity, offtarget_risk.
    """
    ocular_df = tissue.get_tissue_glycans()
    gal3_scores = binding.get_gal3_scores()
    return {
        "gal3_scores":    gal3_scores,
        "gal3_sequence":  binding.get_gal3_sequence(),
        "motif_affinity": binding.motif_affinity_summary(gal3_scores),
        "offtarget_risk": binding.offtarget_risk_table(ocular_df),
    }


def run_cancer() -> dict[str, object]:
    """Run Section 6 (cancer glycan mining) and return a results dict.

    Returns
    -------
    dict with keys: all_disease_df, cancer_df, differential_expr.
    """
    all_disease, cancer_df = cancer.get_cancer_glycans()
    return {
        "all_disease_df":   all_disease,
        "cancer_df":        cancer_df,
        "differential_expr": cancer.differential_expression(),
    }


def run_cross_tissue() -> dict[str, object]:
    """Run Sections 7 & 9 (cross-tissue comparison) and return a results dict.

    Returns
    -------
    dict with keys: ocular_df, cardiac_df, muscle_df, set_overlap, motif_prevalence.
    """
    ocular_df  = tissue.get_tissue_glycans()
    cardiac_df = cross_tissue.get_cardiac_glycans()
    muscle_df  = cross_tissue.get_muscle_glycans()

    tissue_map = {
        "Eye (ocular)":  ocular_df,
        "Heart/Cardiac": cardiac_df,
        "Skeletal Muscle": muscle_df,
    }

    return {
        "ocular_df":        ocular_df,
        "cardiac_df":       cardiac_df,
        "muscle_df":        muscle_df,
        "set_overlap":      cross_tissue.glycan_set_overlap(ocular_df, cardiac_df, muscle_df),
        "motif_prevalence": cross_tissue.gal3_motif_prevalence(tissue_map),
    }


# ──────────────────────────────────────────────────────────────────────────────
# Full pipeline
# ──────────────────────────────────────────────────────────────────────────────

def run_full(
    tissue_terms: list[str] | None = None,
) -> AnalysisResults:
    """Run all four pipeline stages and return a single AnalysisResults object.

    Parameters
    ----------
    tissue_terms : list[str] | None
        Custom tissue search terms for stage 1.

    Returns
    -------
    AnalysisResults
    """
    results = AnalysisResults()

    # Stage 1
    t = run_tissue(tissue_terms)
    results.ocular_df        = t["ocular_df"]
    results.glycan_type_counts = t["glycan_type_counts"]
    results.core1_df         = t["core1_df"]
    results.beta_gal_df      = t["beta_gal_df"]

    # Stage 2
    b = run_binding()
    results.gal3_scores    = b["gal3_scores"]
    results.gal3_sequence  = b["gal3_sequence"]
    results.motif_affinity = b["motif_affinity"]
    results.offtarget_risk = b["offtarget_risk"]

    # Stage 3
    c = run_cancer()
    results.all_disease_df   = c["all_disease_df"]
    results.cancer_df        = c["cancer_df"]
    results.differential_expr = c["differential_expr"]

    # Stage 4
    ct = run_cross_tissue()
    results.cardiac_df      = ct["cardiac_df"]
    results.muscle_df       = ct["muscle_df"]
    results.set_overlap     = ct["set_overlap"]
    results.motif_prevalence = ct["motif_prevalence"]

    return results