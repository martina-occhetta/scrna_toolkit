"""Small internal helpers shared across modules."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def looks_like_counts(X, n_rows: int = 1000) -> bool:
    """Cheap check on the first ``n_rows`` rows: non-negative integers (no densification)."""
    block = X[:n_rows]
    values = block.data if sp.issparse(block) else np.asarray(block).ravel()
    if values.size == 0:
        return True
    return bool((values >= 0).all() and np.allclose(values, np.round(values)))
