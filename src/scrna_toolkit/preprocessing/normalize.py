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

import numpy as np
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData

from scrna_toolkit.counts import ensure_raw_counts
from scrna_toolkit.provenance import get_representation, record_step, set_representation

logger = logging.getLogger(__name__)


def _shares_memory(a: Any, b: Any) -> bool:
    """Whether two matrices share their value buffer (in-place edits to one change the other)."""
    if sp.issparse(a) and sp.issparse(b):
        return bool(np.shares_memory(a.data, b.data))
    if sp.issparse(a) or sp.issparse(b):
        return False
    return bool(np.shares_memory(np.asarray(a), np.asarray(b)))


def store_counts(
    adata: AnnData,
    layer: str = "counts",
    *,
    overwrite: bool = False,
    allow_non_integer: bool = False,
) -> AnnData:
    """Copy raw counts from ``X`` into ``adata.layers[layer]`` (once) before normalization.

    One copy of the matrix is unavoidable because normalization rewrites ``X``. If the layer
    already exists it is kept (no second copy) unless ``overwrite=True``. Either way the stored
    layer is verified with :func:`ensure_raw_counts`.

    An existing layer that shares memory with ``X`` (e.g. after ``adata.layers["counts"] =
    adata.X``) would be rewritten by the in-place normalization of ``X``. In that case ``X`` is
    given its own copy so the stored counts stay untouched.
    """
    if layer in adata.layers and not overwrite:
        ensure_raw_counts(adata, layer, allow_non_integer=allow_non_integer)
        X_existing: Any = adata.X
        if _shares_memory(adata.layers[layer], X_existing):
            logger.warning(
                "layers[%r] shares memory with X; copying X so normalization cannot modify it",
                layer,
            )
            adata.X = X_existing.copy()
        else:
            logger.info("layers[%r] already exists; keeping it", layer)
        return adata
    ensure_raw_counts(adata, allow_non_integer=allow_non_integer)
    X: Any = adata.X
    adata.layers[layer] = X.copy()
    set_representation(adata, "counts")
    record_step(adata, "store_counts", {"layer": layer})
    return adata


def normalize_total(
    adata: AnnData, target_sum: float | None = 1e4, *, allow_non_integer: bool = False
) -> AnnData:
    """Scale each cell to ``target_sum`` total counts (in place on ``X``).

    ``X`` must hold raw counts. It is scanned with :func:`ensure_raw_counts` unless it is
    already tracked as counts (e.g. by :func:`store_counts`).
    """
    representation = get_representation(adata)
    if representation not in (None, "counts"):
        raise ValueError(f"adata.X already holds {representation!r} values; refusing to normalize")
    if representation is None:
        ensure_raw_counts(adata, allow_non_integer=allow_non_integer)
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
