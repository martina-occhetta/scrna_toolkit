"""Guard for functions that assume raw counts."""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import scipy.sparse as sp
from anndata import AnnData

from scrna_toolkit.provenance import get_representation

logger = logging.getLogger(__name__)


def ensure_raw_counts(
    adata: AnnData,
    layer: str | None = None,
    *,
    allow_non_integer: bool = False,
    chunk_size: int = 10_000,
) -> None:
    """Raise ``ValueError`` unless ``adata.X`` (or ``layers[layer]``) holds raw counts.

    Raw counts are non-negative integers. The whole matrix is scanned in row chunks of
    ``chunk_size`` (sparse matrices are never densified). For ``X``, a tracked
    ``X_representation`` other than ``"counts"`` fails immediately.

    ``allow_non_integer=True`` accepts non-negative fractional values, e.g. ambient-corrected
    counts from SoupX. Note that this also stops normalized data from being detected, so use it
    only when you know the values are counts. Negative values always fail.
    """
    if layer is None:
        where = "adata.X"
        representation = get_representation(adata)
        if representation not in (None, "counts"):
            raise ValueError(
                f"adata.X holds {representation!r} values, not raw counts; "
                "use layer='counts' if raw counts were stored there"
            )
        matrix: Any = adata.X
    else:
        where = f"adata.layers[{layer!r}]"
        if layer not in adata.layers:
            raise KeyError(f"{where} does not exist")
        matrix = adata.layers[layer]

    n_negative = n_fractional = 0
    example: list[float] = []
    for start in range(0, adata.n_obs, chunk_size):
        block = matrix[start : start + chunk_size]
        values = block.data if sp.issparse(block) else np.asarray(block).ravel()
        n_negative += int((values < 0).sum())
        fractional = values != np.round(values)
        n_fractional += int(fractional.sum())
        if not example and fractional.any():
            example = [float(v) for v in values[fractional][:3]]

    if n_negative:
        raise ValueError(
            f"{where} has {n_negative} negative values; expected raw counts (was it scaled?)"
        )
    if n_fractional and not allow_non_integer:
        raise ValueError(
            f"{where} has {n_fractional} non-integer values (e.g. {example}); expected raw "
            "counts. If these are corrected counts (e.g. from SoupX), pass "
            "allow_non_integer=True."
        )
    if n_fractional:
        logger.warning("%s has %d non-integer values; accepted as counts", where, n_fractional)
