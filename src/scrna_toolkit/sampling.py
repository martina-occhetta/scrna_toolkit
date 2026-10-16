"""Donor-aware stratified downsampling, e.g. for a laptop-sized development object.

Plain random sampling of a large atlas keeps its biases: a few donors with many cells dominate
and rare cell types can vanish. :func:`stratified_sample` instead splits the cell budget
level by level (by default donor, optionally then sample/library):

1. At each level the budget is shared out by *water-filling*: small groups are taken whole and
   the rest of the budget is split equally among the larger ones.
2. Within the innermost group, cells are allocated across cell types in proportion to their
   size, with a floor so rare types are kept.

Every selected cell gets a weight: the number of cells it represents in its stratum
(innermost group x cell type). Weighted counts recover full-data composition.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from scrna_toolkit.provenance import record_step

logger = logging.getLogger(__name__)

WEIGHT_KEY = "downsample_weight"


def waterfill(sizes: pd.Series, budget: int) -> pd.Series:
    """Allocate ``budget`` as ``min(size, cap)`` with the largest cap that fits.

    Groups smaller than the cap are taken whole; any remainder goes one unit at a time to the
    remaining groups in index order. Returns integers summing to ``min(budget, sizes.sum())``.
    """
    sizes = sizes.astype(np.int64)
    if budget < 0:
        raise ValueError("budget must be non-negative")
    if sizes.sum() <= budget:
        return sizes.copy()
    lo, hi = 0, int(sizes.max())
    while lo < hi:  # largest cap with sum(min(size, cap)) <= budget
        mid = (lo + hi + 1) // 2
        if np.minimum(sizes, mid).sum() <= budget:
            lo = mid
        else:
            hi = mid - 1
    alloc = np.minimum(sizes, lo)
    leftover = budget - int(alloc.sum())
    room = alloc.index[(sizes > alloc).to_numpy()]
    alloc.loc[room[:leftover]] += 1
    return alloc


def _proportional_with_floor(sizes: pd.Series, budget: int, floor: int) -> pd.Series:
    """Allocation proportional to size, with ``min(size, floor)`` guaranteed per group."""
    base = np.minimum(sizes, floor)
    if base.sum() >= budget:
        return waterfill(sizes, budget)
    room = sizes - base
    rest = budget - int(base.sum())  # <= room.sum() because budget <= sizes.sum()
    share = room / room.sum() * rest
    extra = np.floor(share).astype(np.int64)
    remainder = rest - int(extra.sum())  # largest-remainder rounding, so the total is exact
    order = (share - extra).sort_values(ascending=False, kind="stable").index
    extra.loc[order[:remainder]] += 1
    return (base + np.minimum(extra, room)).astype(np.int64)


def _split(positions: np.ndarray, codes: np.ndarray) -> dict[int, np.ndarray]:
    """Group ``positions`` by their code, keeping file order within each group."""
    order = np.argsort(codes[positions], kind="stable")
    keys, starts = np.unique(codes[positions][order], return_index=True)
    return dict(zip(keys.tolist(), np.split(positions[order], starts[1:]), strict=True))


def _stratified_positions(
    obs: pd.DataFrame,
    n_cells: int,
    *,
    levels: Sequence[str],
    cell_type_key: str | None,
    min_cells_per_type: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    if n_cells <= 0:
        raise ValueError("n_cells must be positive")
    if min_cells_per_type < 0:
        raise ValueError("min_cells_per_type must be non-negative")
    missing = [c for c in [*levels, *([cell_type_key] if cell_type_key else [])] if c not in obs]
    if missing:
        raise KeyError(f"Columns not in obs: {missing}")
    level_codes = []
    for column in levels:
        if obs[column].isna().any():
            raise ValueError(f"obs[{column!r}] has missing values; cannot stratify by it")
        level_codes.append(pd.factorize(obs[column], sort=True)[0])
    type_codes = (
        pd.factorize(obs[cell_type_key], sort=True, use_na_sentinel=False)[0]
        if cell_type_key
        else None
    )

    rng = np.random.default_rng(seed)
    picked: list[np.ndarray] = []
    weights: list[np.ndarray] = []

    def take(pool: np.ndarray, k: int) -> None:
        if k > 0:
            picked.append(rng.choice(pool, size=k, replace=False))
            weights.append(np.full(k, len(pool) / k))

    def allocate(positions: np.ndarray, depth: int, budget: int) -> None:
        if depth < len(level_codes):
            groups = _split(positions, level_codes[depth])
            sizes = pd.Series({key: len(pos) for key, pos in groups.items()})
            for key, k in waterfill(sizes, budget).items():
                if k > 0:
                    allocate(groups[key], depth + 1, int(k))
        elif type_codes is None:
            take(positions, min(budget, len(positions)))
        else:
            strata = _split(positions, type_codes)
            sizes = pd.Series({key: len(pos) for key, pos in strata.items()})
            for key, k in _proportional_with_floor(sizes, budget, min_cells_per_type).items():
                take(strata[key], int(k))

    allocate(np.arange(len(obs)), 0, int(n_cells))
    if not picked:
        return np.array([], dtype=np.int64), np.array([], dtype=float)
    positions = np.concatenate(picked)
    order = np.argsort(positions)
    return positions[order], np.concatenate(weights)[order]


def stratified_sample(
    obs: pd.DataFrame,
    n_cells: int,
    *,
    levels: Sequence[str] = ("patient_id",),
    cell_type_key: str | None = None,
    min_cells_per_type: int = 25,
    seed: int = 0,
) -> pd.Series:
    """Choose about ``n_cells`` rows of ``obs``, balanced across donors (and other levels).

    Parameters
    ----------
    levels
        Nested grouping columns, outermost first, e.g. ``("patient_id", "sample_id")``. The
        budget is water-filled across the groups of each level in turn. Missing values raise.
    cell_type_key
        Within the innermost group, allocate across these labels in proportion to size, with
        at least ``min_cells_per_type`` cells (or all of them, if fewer) per label. Missing
        labels form their own stratum. ``None`` samples the innermost group uniformly.
    seed
        Seed for the random choice within strata; equal inputs give equal selections.

    Returns
    -------
    pandas.Series
        Weights indexed by the selected ``obs`` names, in the original row order: the number of
        cells each selected cell represents in its stratum. They sum to ``len(obs)``.

    Exactly ``min(n_cells, len(obs))`` cells are selected. When the cell-type floors alone
    exceed a group's budget, that group's budget is water-filled across cell types instead.
    """
    positions, weights = _stratified_positions(
        obs,
        n_cells,
        levels=levels,
        cell_type_key=cell_type_key,
        min_cells_per_type=min_cells_per_type,
        seed=seed,
    )
    return pd.Series(weights, index=obs.index[positions], name=WEIGHT_KEY)


def downsample_cells(
    adata: AnnData,
    n_cells: int,
    *,
    levels: Sequence[str] = ("patient_id",),
    cell_type_key: str | None = None,
    min_cells_per_type: int = 25,
    seed: int = 0,
    weight_key: str = WEIGHT_KEY,
) -> AnnData:
    """Return a donor-aware stratified subset of ``adata`` (a copy) with weights in ``obs``.

    See :func:`stratified_sample` for the parameters. The weights are stored in
    ``obs[weight_key]``; use them to recover full-data composition from the subset.
    """
    positions, weights = _stratified_positions(
        adata.obs,
        n_cells,
        levels=levels,
        cell_type_key=cell_type_key,
        min_cells_per_type=min_cells_per_type,
        seed=seed,
    )
    subset = adata[positions].copy()
    subset.obs[weight_key] = weights
    logger.info("Downsampled %d -> %d cells (levels %s)", adata.n_obs, subset.n_obs, list(levels))
    record_step(
        subset,
        "downsample_cells",
        {
            "n_cells": n_cells,
            "n_obs_before": adata.n_obs,
            "levels": list(levels),
            "cell_type_key": cell_type_key,
            "min_cells_per_type": min_cells_per_type,
            "seed": seed,
        },
    )
    return subset
