"""Cell counts per group and cell-type composition plots."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData
from matplotlib.figure import Figure

from scrna_toolkit.plotting._common import (
    category_colors,
    category_width,
    finish,
    label_categories,
    missing,
)


def plot_cells_per_group(
    adata: AnnData, groupby: str, *, log: bool = False, path: str | Path | None = None
) -> Figure | None:
    """Bar chart of cell counts per category, largest first."""
    if missing(adata, [groupby], "plot_cells_per_group"):
        return None
    counts = adata.obs[groupby].value_counts()
    fig, ax = plt.subplots(figsize=(category_width(len(counts)), 3.5))
    ax.bar(range(len(counts)), counts.to_numpy())
    ax.set_xlabel(groupby)
    ax.set_ylabel("cells")
    if log:
        ax.set_yscale("log")
    label_categories(ax, counts.index.astype(str).tolist())
    fig.tight_layout()
    return finish(fig, path)


def plot_cells_per_sample(adata: AnnData, sample_col: str = "sample_id", **kwargs) -> Figure | None:
    return plot_cells_per_group(adata, sample_col, **kwargs)


def plot_cells_per_patient(
    adata: AnnData, patient_col: str = "patient_id", **kwargs
) -> Figure | None:
    return plot_cells_per_group(adata, patient_col, **kwargs)


def _composition(adata: AnnData, groupby: str, cell_type_col: str, normalize: bool) -> pd.DataFrame:
    table = pd.crosstab(adata.obs[groupby], adata.obs[cell_type_col])
    return table.div(table.sum(axis=1), axis=0) if normalize else table


def plot_cell_type_composition(
    adata: AnnData,
    groupby: str = "sample_id",
    cell_type_col: str = "cell_type",
    *,
    normalize: bool = True,
    path: str | Path | None = None,
) -> Figure | None:
    """Stacked bars of cell-type proportions (or counts) per group."""
    if missing(adata, [groupby, cell_type_col], "plot_cell_type_composition"):
        return None
    table = _composition(adata, groupby, cell_type_col, normalize)
    cell_types = table.columns.astype(str).tolist()
    colors = category_colors(adata, cell_type_col, cell_types)

    fig, ax = plt.subplots(figsize=(category_width(len(table)) + 2.5, 4))
    bottom = np.zeros(len(table))
    for cell_type, color in zip(table.columns, colors, strict=True):
        values = table[cell_type].to_numpy(dtype=float)
        ax.bar(range(len(table)), values, bottom=bottom, color=color, label=str(cell_type))
        bottom += values
    ax.set_xlabel(groupby)
    ax.set_ylabel("fraction of cells" if normalize else "cells")
    label_categories(ax, table.index.astype(str).tolist())
    ax.legend(
        title=cell_type_col,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize="small",
        ncols=max(1, len(cell_types) // 25),
    )
    fig.tight_layout()
    return finish(fig, path)


def plot_composition_heatmap(
    adata: AnnData,
    groupby: str = "sample_id",
    cell_type_col: str = "cell_type",
    *,
    normalize: bool = True,
    path: str | Path | None = None,
) -> Figure | None:
    """Heatmap of cell-type proportions (rows: groups, columns: cell types)."""
    if missing(adata, [groupby, cell_type_col], "plot_composition_heatmap"):
        return None
    table = _composition(adata, groupby, cell_type_col, normalize)
    fig, ax = plt.subplots(
        figsize=(category_width(table.shape[1]) + 1.5, category_width(table.shape[0], 3))
    )
    image = ax.imshow(table.to_numpy(dtype=float), aspect="auto")
    fig.colorbar(image, ax=ax, label="fraction of cells" if normalize else "cells")
    ax.set_xlabel(cell_type_col)
    label_categories(ax, table.columns.astype(str).tolist())
    rows = table.index.astype(str).tolist()
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows if len(rows) <= 60 else [])
    ax.set_ylabel(groupby)
    fig.tight_layout()
    return finish(fig, path)
