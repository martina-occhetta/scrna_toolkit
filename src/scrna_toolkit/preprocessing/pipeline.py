"""The standard preprocessing pipeline, composed from the individual steps.

Stages and the representation each one uses:

1. ``store_counts``   X (counts) -> copied to ``layers["counts"]``
2. ``normalize_total`` X: counts -> library-size normalized
3. ``log1p``          X: normalized -> log1p (kept for all genes, sparse)
4. ``select_hvgs``    flags ``var["highly_variable"]`` from log1p X
5. ``run_pca``        scaled temporary copy of HVG columns -> ``obsm["X_pca"]``
6. ``run_neighbors``  ``obsm["X_pca"]`` -> kNN graph in ``obsp``
7. ``run_umap``       graph -> ``obsm["X_umap"]`` (optional)
8. ``run_leiden``     graph -> ``obs["leiden"]`` (optional)
"""

from __future__ import annotations

import dataclasses
import logging
from typing import Any

from anndata import AnnData

from scrna_toolkit.config import PreprocessingConfig
from scrna_toolkit.preprocessing.dimensionality import run_pca
from scrna_toolkit.preprocessing.features import select_hvgs
from scrna_toolkit.preprocessing.graph import run_leiden, run_neighbors, run_umap
from scrna_toolkit.preprocessing.normalize import log1p, normalize_total, store_counts
from scrna_toolkit.provenance import record_config

logger = logging.getLogger(__name__)


def preprocess(
    adata: AnnData,
    config: PreprocessingConfig | None = None,
    *,
    copy: bool = False,
    **overrides: Any,
) -> AnnData:
    """Run the standard pipeline (see module docstring) and return the processed AnnData.

    ``overrides`` replace individual :class:`PreprocessingConfig` fields, e.g.
    ``preprocess(adata, n_hvgs=2000, n_pcs=30)``. Works in place unless ``copy=True``.
    Expects QC-filtered raw counts in ``X`` (or an existing ``layers["counts"]``).
    """
    cfg = dataclasses.replace(config or PreprocessingConfig(), **overrides)
    if copy:
        adata = adata.copy()
    logger.info("preprocess: %d cells x %d genes, config=%s", adata.n_obs, adata.n_vars, cfg)

    store_counts(adata, allow_non_integer=cfg.allow_non_integer_counts)
    normalize_total(adata, target_sum=cfg.target_sum)
    log1p(adata)
    select_hvgs(
        adata,
        cfg.n_hvgs,
        flavor=cfg.hvg_flavor,
        batch_key=cfg.hvg_batch_key,
        allow_non_integer=cfg.allow_non_integer_counts,
    )
    run_pca(adata, cfg.n_pcs, scale=cfg.scale, max_value=cfg.scale_max_value, seed=cfg.seed)
    run_neighbors(adata, cfg.n_neighbors, seed=cfg.seed)
    if cfg.run_umap:
        run_umap(adata, min_dist=cfg.umap_min_dist, seed=cfg.seed)
    if cfg.run_leiden:
        run_leiden(adata, cfg.leiden_resolution, seed=cfg.seed)

    record_config(adata, "preprocessing", dataclasses.asdict(cfg))
    return adata
