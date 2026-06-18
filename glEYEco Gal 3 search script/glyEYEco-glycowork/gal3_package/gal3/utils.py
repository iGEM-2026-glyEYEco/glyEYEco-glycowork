"""
gal3.utils
----------
Shared helpers used across all modules:
  - col_contains_any  : flexible list/string term matcher for df_glycan columns
  - load_data         : cached loader for df_glycan and glycan_binding
"""

from __future__ import annotations
import functools
from typing import Any


def col_contains_any(col_value: Any, terms: list[str]) -> bool:
    """Return True if any element of *col_value* (str or list) matches
    any term in *terms* (case-insensitive substring match).

    Parameters
    ----------
    col_value : str | list | None | float
        A single cell value from a df_glycan column (tissue_sample,
        disease_association, species …).
    terms : list[str]
        Search terms, e.g. ['eye', 'cornea', 'tears'].

    Returns
    -------
    bool
    """
    if col_value is None or isinstance(col_value, float):
        return False
    items: list[str] = col_value if isinstance(col_value, list) else [str(col_value)]
    return any(term.lower() in item.lower() for item in items for term in terms)


@functools.lru_cache(maxsize=1)
def load_data():
    """Load and cache df_glycan and glycan_binding from glycowork.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        (df_glycan, glycan_binding)
    """
    from glycowork.glycan_data.loader import df_glycan, glycan_binding  # type: ignore
    return df_glycan, glycan_binding