"""Count preservation and normalization.

Representations used by this package:

- ``adata.layers["counts"]``: the untouched raw count matrix (preferred place for counts).
- ``adata.X``: starts as counts; after :func:`normalize_total` holds library-size normalized
  values, after :func:`log1p` log-normalized expression (all genes, still sparse).
- ``adata.raw``: not used. It freezes a full AnnData snapshot and is easy to confuse with raw
  counts; ``layers["counts"]`` is explicit.

The current state of ``X`` is tracked in ``adata.uns["scrna_toolkit"]["X_representation"]``.
"""

from __future__ import annotations

import logging
from typing import Any

import scanpy as sc
from anndata import AnnData

from scrna_toolkit._utils import looks_like_counts
from scrna_toolkit.provenance import get_representation, record_step, set_representation

logger = logging.getLogger(__name__)


def store_counts(adata: AnnData, layer: str = "counts", *, overwrite: bool = False) -> AnnData:
    """Copy raw counts from ``X`` into ``adata.layers[layer]`` (once) before normalization.

    One copy of the matrix is unavoidable because normalization rewrites ``X``. If the layer
    already exists it is kept as is (no second copy) unless ``overwrite=True``.
    """
    if layer in adata.layers and not overwrite:
        logger.info("layers[%r] already exists; keeping it", layer)
        return adata
    representation = get_representation(adata)
    if representation not in (None, "counts"):
        raise ValueError(f"adata.X holds {representation!r} values, not counts")
    X: Any = adata.X
    if not looks_like_counts(X):
        logger.warning("adata.X is not integer-valued; storing it as %r anyway", layer)
    adata.layers[layer] = X.copy()
    set_representation(adata, "counts")
    record_step(adata, "store_counts", {"layer": layer})
    return adata


def normalize_total(adata: AnnData, target_sum: float | None = 1e4) -> AnnData:
    """Scale each cell to ``target_sum`` total counts (in place on ``X``)."""
    representation = get_representation(adata)
    if representation not in (None, "counts"):
        raise ValueError(f"adata.X already holds {representation!r} values; refusing to normalize")
    if representation is None and not looks_like_counts(adata.X):
        logger.warning("normalize_total expects counts but adata.X is not integer-valued")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    set_representation(adata, "normalized")
    logger.info("normalize_total: X counts -> normalized (target_sum=%s)", target_sum)
    record_step(adata, "normalize_total", {"target_sum": target_sum})
    return adata


def log1p(adata: AnnData) -> AnnData:
    """Apply ``log(1 + x)`` to ``X`` in place."""
    representation = get_representation(adata)
    if representation in ("log1p", "scaled"):
        raise ValueError(f"adata.X already holds {representation!r} values")
    if representation == "counts":
        logger.warning("log1p applied to un-normalized counts")
    sc.pp.log1p(adata)
    set_representation(adata, "log1p")
    logger.info("log1p: X -> log1p")
    record_step(adata, "log1p")
    return adata
