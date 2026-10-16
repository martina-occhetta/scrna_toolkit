"""Small internal helpers shared across modules."""

from __future__ import annotations

import logging

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
