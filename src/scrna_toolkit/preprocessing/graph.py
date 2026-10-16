"""Neighbour graph, UMAP and Leiden clustering (thin, explicit Scanpy wrappers)."""

from __future__ import annotations

import logging

import scanpy as sc
from anndata import AnnData

from scrna_toolkit.provenance import record_step

logger = logging.getLogger(__name__)


def run_neighbors(
    adata: AnnData,
    n_neighbors: int = 15,
    *,
    use_rep: str = "X_pca",
    n_pcs: int | None = None,
    seed: int = 0,
) -> AnnData:
    """Build the kNN graph from ``obsm[use_rep]`` (stored in ``obsp``)."""
    if use_rep not in adata.obsm:
        raise KeyError(f"obsm[{use_rep!r}] missing; run run_pca() first")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep, n_pcs=n_pcs, random_state=seed)
    logger.info("run_neighbors: k=%d on obsm[%r]", n_neighbors, use_rep)
    record_step(
        adata,
        "run_neighbors",
        {"n_neighbors": n_neighbors, "use_rep": use_rep, "n_pcs": n_pcs, "seed": seed},
    )
    return adata


def _require_neighbors(adata: AnnData) -> None:
    if "neighbors" not in adata.uns:
        raise KeyError("Neighbour graph missing; run run_neighbors() first")


def run_umap(adata: AnnData, *, min_dist: float = 0.5, seed: int = 0) -> AnnData:
    """Compute a 2-D UMAP into ``obsm["X_umap"]`` from the neighbour graph."""
    _require_neighbors(adata)
    sc.tl.umap(adata, min_dist=min_dist, random_state=seed)
    record_step(adata, "run_umap", {"min_dist": min_dist, "seed": seed})
    return adata


def run_leiden(
    adata: AnnData, resolution: float = 1.0, *, key_added: str = "leiden", seed: int = 0
) -> AnnData:
    """Leiden clustering on the neighbour graph into ``obs[key_added]`` (igraph backend)."""
    _require_neighbors(adata)
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=key_added,
        random_state=seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    logger.info(
        "run_leiden: %d clusters (resolution=%s)", adata.obs[key_added].nunique(), resolution
    )
    record_step(
        adata, "run_leiden", {"resolution": resolution, "key_added": key_added, "seed": seed}
    )
    return adata
