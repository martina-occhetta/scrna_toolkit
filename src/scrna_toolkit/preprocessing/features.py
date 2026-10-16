"""Highly variable gene selection."""

from __future__ import annotations

import logging

import scanpy as sc
from anndata import AnnData

from scrna_toolkit.counts import ensure_raw_counts
from scrna_toolkit.provenance import get_representation, record_step

logger = logging.getLogger(__name__)

_COUNT_FLAVORS = ("seurat_v3", "seurat_v3_paper")


def select_hvgs(
    adata: AnnData,
    n_top_genes: int = 3000,
    *,
    flavor: str = "seurat",
    batch_key: str | None = None,
    layer: str | None = None,
    subset: bool = False,
    allow_non_integer: bool = False,
) -> AnnData:
    """Flag highly variable genes in ``var["highly_variable"]`` (genes are kept by default).

    ``seurat``/``cell_ranger`` flavors use log-normalized ``X``; ``seurat_v3`` flavors use raw
    counts and default to ``layers["counts"]`` (and need the optional ``scikit-misc`` package).
    """
    if n_top_genes > adata.n_vars:
        raise ValueError(f"n_top_genes={n_top_genes} exceeds the {adata.n_vars} genes available")
    if flavor in _COUNT_FLAVORS and layer is None:
        if "counts" not in adata.layers:
            raise ValueError(f"flavor={flavor!r} needs raw counts; run store_counts() first")
        layer = "counts"
    if flavor in _COUNT_FLAVORS:
        ensure_raw_counts(adata, layer, allow_non_integer=allow_non_integer)
    if flavor not in _COUNT_FLAVORS and layer is None and get_representation(adata) != "log1p":
        logger.warning("flavor=%r expects log-normalized X", flavor)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=flavor,
        batch_key=batch_key,
        layer=layer,
        subset=subset,
    )
    logger.info(
        "select_hvgs: %d HVGs (flavor=%s, input=%s)",
        int(adata.var["highly_variable"].sum()),
        flavor,
        f"layers[{layer!r}]" if layer else "X",
    )
    record_step(
        adata,
        "select_hvgs",
        {"n_top_genes": n_top_genes, "flavor": flavor, "batch_key": batch_key, "layer": layer},
    )
    return adata
