"""Small internal helpers shared across modules."""

from __future__ import annotations

import logging

import anndata as ad
import pandas as pd

logger = logging.getLogger(__name__)

# Warn before densifying more than this many float32 values (~0.8 GB).
DENSIFY_WARN_ELEMENTS = 2e8


def warn_if_densifying(n_obs: int, n_vars: int, what: str) -> None:
    """Log a warning when an operation will create a large dense matrix."""
    n = n_obs * n_vars
    if n > DENSIFY_WARN_ELEMENTS:
        logger.warning(
            "%s densifies a %d x %d matrix (~%.1f GB as float32)", what, n_obs, n_vars, n * 4e-9
        )


def unique_gene_names(names: pd.Index) -> pd.Index:
    """Make duplicated gene names unique (``ACTB``, ``ACTB-1``), logging how many changed."""
    if names.is_unique:
        return names
    n_dup = int(names.duplicated().sum())
    logger.warning("%d duplicated gene names made unique (gene_ids are preserved)", n_dup)
    return ad.utils.make_index_unique(names)
