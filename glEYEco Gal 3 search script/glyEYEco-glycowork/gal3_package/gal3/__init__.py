"""
gal3 — Galectin-3 Computational Validation Package
====================================================
Wraps the glEYEco × glycowork analysis pipeline so it can be called
programmatically from any Python script or web interface.

Quick start
-----------
>>> from gal3 import pipeline
>>> results = pipeline.run_full(tissue="ocular")
"""

from .pipeline import run_full
from . import binding, cancer, cross_tissue, docking, tissue 
from . import utils

__version__ = "0.1.0"
__all__ = ["run_full", "tissue", "binding", "cancer", "cross_tissue", "docking", "utils"]