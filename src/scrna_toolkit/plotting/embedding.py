"""Embedding plots (thin wrapper over ``scanpy.pl.embedding``)."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import scanpy as sc
from anndata import AnnData
from matplotlib.figure import Figure

from scrna_toolkit.plotting._common import finish

logger = logging.getLogger(__name__)


def plot_umap(
    adata: AnnData,
    color: str | Sequence[str] | None = None,
    *,
    basis: str = "umap",
    ncols: int = 3,
    path: str | Path | None = None,
    **kwargs: Any,
) -> Figure | None:
    """Plot ``obsm["X_<basis>"]`` coloured by obs columns or genes; absent keys are skipped."""
    if f"X_{basis}" not in adata.obsm:
        logger.warning("plot_umap skipped: obsm['X_%s'] not found", basis)
        return None
    requested = [color] if isinstance(color, str) else list(color or [])
    available = [c for c in requested if c in adata.obs or c in adata.var_names]
    if len(available) < len(requested):
        logger.warning(
            "plot_umap: skipping unknown keys %s", sorted(set(requested) - set(available))
        )
    fig = sc.pl.embedding(
        adata,
        basis=basis,
        color=available or None,
        ncols=ncols,
        show=False,
        return_fig=True,
        **kwargs,
    )
    return finish(fig, path)
