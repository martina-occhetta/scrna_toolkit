"""Shared plotting helpers: saving, sizing for many categories, Scanpy-compatible colours."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from anndata import AnnData
from matplotlib.axes import Axes
from matplotlib.figure import Figure

logger = logging.getLogger(__name__)

MAX_LABELLED_CATEGORIES = 60


def finish(fig: Figure, path: str | Path | None) -> Figure:
    """Save ``fig`` if a path is given; always return it (never shows or closes it)."""
    if path is not None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, bbox_inches="tight", dpi=150)
    return fig


def missing(adata: AnnData, columns: list[str], plot: str) -> bool:
    """Log and return True if any ``obs`` column needed by ``plot`` is absent."""
    absent = [c for c in columns if c not in adata.obs]
    if absent:
        logger.warning("%s skipped: obs columns %s not found", plot, absent)
    return bool(absent)


def category_width(n: int, base: float = 4.0) -> float:
    return float(min(max(base, 0.25 * n), 30.0))


def label_categories(ax: Axes, labels: list[str]) -> None:
    """Rotated tick labels, hidden when there are too many to read."""
    ax.set_xticks(range(len(labels)))
    if len(labels) > MAX_LABELLED_CATEGORIES:
        ax.set_xticklabels([])
        ax.set_xlabel(f"{ax.get_xlabel()} ({len(labels)} categories, labels hidden)")
    else:
        ax.set_xticklabels(labels, rotation=90)


def category_colors(adata: AnnData, column: str, categories: list[str]) -> list:
    """Use Scanpy's ``uns["<column>_colors"]`` when it matches, else the matplotlib defaults."""
    stored = adata.uns.get(f"{column}_colors")
    obs_categories = getattr(adata.obs[column], "cat", None)
    if stored is not None and obs_categories is not None:
        lookup = dict(zip(obs_categories.categories.astype(str), stored, strict=False))
        if all(c in lookup for c in categories):
            return [lookup[c] for c in categories]
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if len(categories) <= len(cycle):
        return list(cycle[: len(categories)])
    cmap = plt.get_cmap(plt.rcParams["image.cmap"])
    return [cmap(x) for x in np.linspace(0, 1, len(categories))]
