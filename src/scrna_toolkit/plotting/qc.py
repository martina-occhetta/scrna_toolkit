"""QC metric distributions and scatter plots."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path
from typing import cast

import matplotlib.pyplot as plt
import pandas as pd
from anndata import AnnData
from matplotlib.figure import Figure

from scrna_toolkit.plotting._common import category_width, finish, label_categories, missing
from scrna_toolkit.summary.dataset import DEFAULT_QC_COLS

logger = logging.getLogger(__name__)


def plot_qc_distributions(
    adata: AnnData,
    metrics: Sequence[str] = DEFAULT_QC_COLS,
    *,
    groupby: str | None = None,
    bins: int = 50,
    path: str | Path | None = None,
) -> Figure | None:
    """Histogram per QC metric, or per-group boxplots when ``groupby`` is given."""
    present = [m for m in metrics if m in adata.obs]
    if not present:
        logger.warning("plot_qc_distributions skipped: no QC metrics; run calculate_qc_metrics()")
        return None
    if groupby is not None and missing(adata, [groupby], "grouping"):
        groupby = None

    if groupby is None:
        fig, axes = plt.subplots(1, len(present), figsize=(4 * len(present), 3.2), squeeze=False)
        for ax, metric in zip(axes[0], present, strict=True):
            ax.hist(adata.obs[metric].to_numpy(), bins=bins)
            ax.set_xlabel(metric)
            ax.set_ylabel("cells")
    else:
        groups = cast(pd.DataFrame, adata.obs).groupby(groupby, observed=True)
        labels = [str(k) for k in groups.groups]
        fig, axes = plt.subplots(
            len(present),
            1,
            figsize=(category_width(len(labels), 6), 2.8 * len(present)),
            squeeze=False,
        )
        for ax, metric in zip(axes[:, 0], present, strict=True):
            ax.boxplot(
                [g[metric].to_numpy() for _, g in groups],
                positions=range(len(labels)),
                showfliers=False,
            )
            ax.set_ylabel(metric)
            ax.set_xlabel(groupby)
            label_categories(ax, labels)
    fig.tight_layout()
    return finish(fig, path)


def plot_qc_scatter(
    adata: AnnData,
    x: str = "total_counts",
    y: str = "n_genes_by_counts",
    *,
    color: str | None = "pct_counts_mt",
    path: str | Path | None = None,
) -> Figure | None:
    """Scatter of two QC metrics, coloured by a third (skipped if absent)."""
    if missing(adata, [x, y], "plot_qc_scatter"):
        return None
    fig, ax = plt.subplots(figsize=(4.5, 4))
    c = adata.obs[color].to_numpy() if color and color in adata.obs else None
    points = ax.scatter(adata.obs[x], adata.obs[y], c=c, s=3, rasterized=True)
    if c is not None:
        fig.colorbar(points, ax=ax, label=color)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    fig.tight_layout()
    return finish(fig, path)
