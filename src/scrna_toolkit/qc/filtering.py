"""Threshold-based cell and gene filtering.

Thresholds are dataset-specific; every default here is ``None`` (no filtering).
By default a filtered copy is returned and the input is left untouched.
"""

from __future__ import annotations

import logging
import operator
from collections.abc import Callable
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData

from scrna_toolkit.provenance import record_step

logger = logging.getLogger(__name__)

# parameter name -> (obs/var column, comparison that a *passing* value satisfies)
_CELL_CRITERIA: dict[str, tuple[str, Callable]] = {
    "min_genes": ("n_genes_by_counts", operator.ge),
    "max_genes": ("n_genes_by_counts", operator.le),
    "min_counts": ("total_counts", operator.ge),
    "max_counts": ("total_counts", operator.le),
    "max_pct_mt": ("pct_counts_mt", operator.le),
}


def _apply(
    adata: AnnData, mask: np.ndarray, axis: str, inplace: bool, step: str, params: dict
) -> AnnData:
    n_before = adata.n_obs if axis == "obs" else adata.n_vars
    if inplace:
        if axis == "obs":
            adata._inplace_subset_obs(mask)
        else:
            adata._inplace_subset_var(mask)
        out = adata
    else:
        out = (adata[mask] if axis == "obs" else adata[:, mask]).copy()
    n_after = int(mask.sum())
    what = "cells" if axis == "obs" else "genes"
    logger.info("%s: %d -> %d %s (%d removed)", step, n_before, n_after, what, n_before - n_after)
    record_step(out, step, {**params, "n_before": n_before, "n_after": n_after})
    return out


def _summary(rows: list[tuple[str, float, int]], n_removed: int) -> pd.DataFrame:
    summary = pd.DataFrame(rows, columns=["criterion", "threshold", "n_failing"])
    total = pd.DataFrame([("total_removed", np.nan, n_removed)], columns=summary.columns)
    return pd.concat([summary, total], ignore_index=True).set_index("criterion")


def filter_cells(
    adata: AnnData,
    *,
    min_genes: int | None = None,
    max_genes: int | None = None,
    min_counts: float | None = None,
    max_counts: float | None = None,
    max_pct_mt: float | None = None,
    inplace: bool = False,
    return_summary: bool = False,
) -> AnnData | tuple[AnnData, pd.DataFrame]:
    """Remove cells failing any of the given QC thresholds.

    Requires the columns added by :func:`calculate_qc_metrics`. ``max_pct_mt`` is a
    percentage (0-100). With ``return_summary=True`` also returns a table of how many cells
    fail each criterion (a cell can fail several) and the total removed.
    """
    thresholds = {
        "min_genes": min_genes,
        "max_genes": max_genes,
        "min_counts": min_counts,
        "max_counts": max_counts,
        "max_pct_mt": max_pct_mt,
    }
    mask = np.ones(adata.n_obs, dtype=bool)
    rows: list[tuple[str, float, int]] = []
    for name, value in thresholds.items():
        if value is None:
            continue
        column, passes = _CELL_CRITERIA[name]
        if column not in adata.obs:
            raise KeyError(f"{name} needs obs[{column!r}]; run calculate_qc_metrics() first")
        ok = passes(adata.obs[column].to_numpy(), value)
        rows.append((name, value, int((~ok).sum())))
        mask &= ok

    out = _apply(adata, mask, "obs", inplace, "filter_cells", thresholds)
    return (out, _summary(rows, int((~mask).sum()))) if return_summary else out


def filter_genes(
    adata: AnnData,
    *,
    min_cells: int | None = None,
    min_counts: float | None = None,
    layer: str | None = None,
    inplace: bool = False,
    return_summary: bool = False,
) -> AnnData | tuple[AnnData, pd.DataFrame]:
    """Remove genes detected in fewer than ``min_cells`` cells or with < ``min_counts`` counts.

    Statistics are computed from the current matrix (so they reflect any earlier cell
    filtering), sparse-friendly and without densifying.
    """
    matrix: Any = adata.layers[layer] if layer else adata.X
    mask = np.ones(adata.n_vars, dtype=bool)
    rows: list[tuple[str, float, int]] = []
    if min_cells is not None:
        ok = np.asarray((matrix > 0).sum(axis=0)).ravel() >= min_cells
        rows.append(("min_cells", min_cells, int((~ok).sum())))
        mask &= ok
    if min_counts is not None:
        ok = np.asarray(matrix.sum(axis=0)).ravel() >= min_counts
        rows.append(("min_counts", min_counts, int((~ok).sum())))
        mask &= ok

    params = {"min_cells": min_cells, "min_counts": min_counts, "layer": layer}
    out = _apply(adata, mask, "var", inplace, "filter_genes", params)
    return (out, _summary(rows, int((~mask).sum()))) if return_summary else out
