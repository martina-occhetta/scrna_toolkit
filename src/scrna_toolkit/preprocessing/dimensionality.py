"""Scaling and PCA.

By default PCA is computed on a *scaled copy of the HVG columns only*, so ``adata.X`` keeps
its sparse log-normalized values for all genes and only an ``n_cells x n_hvgs`` dense block is
created temporarily.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData

from scrna_toolkit._utils import warn_if_densifying
from scrna_toolkit.provenance import record_step, set_representation

logger = logging.getLogger(__name__)


def scale(adata: AnnData, max_value: float | None = 10.0, *, zero_center: bool = True) -> AnnData:
    """Scale ``X`` to unit variance (and zero mean) in place. Densifies ``X``.

    Usually not needed directly: :func:`run_pca` with ``scale=True`` scales a temporary copy.
    """
    if zero_center:
        warn_if_densifying(adata.n_obs, adata.n_vars, "scale")
    sc.pp.scale(adata, max_value=max_value, zero_center=zero_center)
    set_representation(adata, "scaled")
    record_step(adata, "scale", {"max_value": max_value, "zero_center": zero_center})
    return adata


def run_pca(
    adata: AnnData,
    n_comps: int = 50,
    *,
    use_hvgs: bool = True,
    scale: bool = True,
    max_value: float | None = 10.0,
    seed: int = 0,
) -> AnnData:
    """Compute PCA into ``obsm["X_pca"]``, ``varm["PCs"]`` and ``uns["pca"]``.

    Parameters
    ----------
    use_hvgs
        Restrict to ``var["highly_variable"]`` genes.
    scale
        Scale (unit variance, clipped at ``max_value``) a temporary copy of the selected genes
        before PCA. ``adata.X`` is left unchanged. With ``scale=False`` PCA runs on ``X``
        directly (sparse-friendly).
    """
    if use_hvgs:
        if "highly_variable" not in adata.var:
            raise KeyError("var['highly_variable'] missing; run select_hvgs() or use_hvgs=False")
        mask = adata.var["highly_variable"].to_numpy()
    else:
        mask = np.ones(adata.n_vars, dtype=bool)
    n_features = int(mask.sum())
    if n_comps >= min(adata.n_obs, n_features):
        raise ValueError(
            f"n_comps={n_comps} must be < min(n_cells={adata.n_obs}, n_genes={n_features})"
        )

    if scale:
        warn_if_densifying(adata.n_obs, n_features, "run_pca(scale=True)")
        X: Any = adata.X[:, mask]  # type: ignore[index]
        sub = AnnData(X=X.toarray() if sp.issparse(X) else np.array(X))
        sc.pp.scale(sub, max_value=max_value)
        sc.pp.pca(sub, n_comps=n_comps, random_state=seed)
        adata.obsm["X_pca"] = sub.obsm["X_pca"]
        pcs: Any = sub.varm["PCs"]
        loadings = np.zeros((adata.n_vars, n_comps), dtype=pcs.dtype)
        loadings[mask] = pcs
        adata.varm["PCs"] = loadings
        adata.uns["pca"] = {k: sub.uns["pca"][k] for k in ("variance", "variance_ratio")}
    else:
        sc.pp.pca(adata, n_comps=n_comps, mask_var=mask if use_hvgs else None, random_state=seed)

    input_desc = f"{'scaled copy of ' if scale else ''}X[:, {n_features} genes]"
    logger.info("run_pca: %d components from %s", n_comps, input_desc)
    record_step(
        adata,
        "run_pca",
        {
            "n_comps": n_comps,
            "use_hvgs": use_hvgs,
            "n_genes": n_features,
            "scale": scale,
            "max_value": max_value,
            "seed": seed,
        },
    )
    return adata
