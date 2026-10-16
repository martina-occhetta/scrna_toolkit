"""Small internal helpers shared across modules."""

from __future__ import annotations

import logging

import numpy as np
import scipy.sparse as sp

logger = logging.getLogger(__name__)

# Warn before densifying more than this many float32 values (~0.8 GB).
DENSIFY_WARN_ELEMENTS = 2e8


def looks_like_counts(X, n_rows: int = 1000) -> bool:
    """Cheap check on the first ``n_rows`` rows: non-negative integers (no densification)."""
    block = X[:n_rows]
    values = block.data if sp.issparse(block) else np.asarray(block).ravel()
    if values.size == 0:
        return True
    return bool((values >= 0).all() and np.allclose(values, np.round(values)))


def warn_if_densifying(n_obs: int, n_vars: int, what: str) -> None:
    """Log a warning when an operation will create a large dense matrix."""
    n = n_obs * n_vars
    if n > DENSIFY_WARN_ELEMENTS:
        logger.warning(
            "%s densifies a %d x %d matrix (~%.1f GB as float32)", what, n_obs, n_vars, n * 4e-9
        )
