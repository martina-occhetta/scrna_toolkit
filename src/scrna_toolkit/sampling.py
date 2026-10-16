"""Donor-aware stratified downsampling, e.g. for a laptop-sized development object.

Plain random sampling of a large atlas keeps its biases: a few donors with many cells dominate
and rare cell types can vanish. :func:`stratified_sample` instead splits the cell budget
level by level (by default donor, optionally then sample/library):

1. At each level the budget is shared out by *water-filling*: small groups are taken whole and
   the rest of the budget is split equally among the larger ones.
2. Within the innermost group, cells are allocated across cell types in proportion to their
   size, with a floor so rare types are kept.

Every stratum (innermost group x cell type) gets at least one cell; a budget too small for that
raises instead of silently dropping donors. Each selected cell's weight is the number of cells
it represents in its stratum, so the weights sum to ``len(obs)`` and weighted counts recover
full-data composition.
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
    """Allocation proportional to size, with ``min(size, floor)`` guaranteed per group.

    Falls back to :func:`waterfill` when the floors alone exceed ``budget``.
    """
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
    return (base + extra).astype(np.int64)


def _split(positions: np.ndarray, codes: np.ndarray) -> dict[int, np.ndarray]:
    """Group ``positions`` by their code, keeping file order within each group."""
    order = np.argsort(codes[positions], kind="stable")
    keys, starts = np.unique(codes[positions][order], return_index=True)
    return dict(zip(keys.tolist(), np.split(positions[order], starts[1:]), strict=True))


def stratified_sample(
    obs: pd.DataFrame,
    n_cells: int,
    *,
    levels: Sequence[str] = ("patient_id",),
    cell_type_key: str | None = None,
    min_cells_per_type: int = 25,
    seed: int = 0,
) -> pd.Series:
    """Choose ``min(n_cells, len(obs))`` rows of ``obs``, balanced across donors.

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

    Raises
    ------
    ValueError
        If ``n_cells`` is too small to give every stratum (every donor, every group of the
        other levels and every cell type within them) at least one cell. Dropping a donor
        would make the sample, and its weights, unrepresentative.
    """
    if n_cells <= 0:
        raise ValueError("n_cells must be positive")
    if min_cells_per_type < 0:
        raise ValueError("min_cells_per_type must be non-negative")
    keys = [*levels, *([cell_type_key] if cell_type_key else [])]
    missing = [c for c in keys if c not in obs]
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
    n_strata = obs.groupby(keys, dropna=False, observed=True).ngroups if keys else 1

    def check_no_empty(alloc: pd.Series, what: str) -> None:
        n_empty = int((alloc == 0).sum())
        if n_empty:
            raise ValueError(
                f"n_cells={n_cells} is too small: {n_empty} of {len(alloc)} {what} groups would "
                f"get no cells, so they could not be represented. Increase n_cells (at least "
                f"{n_strata}, one per stratum, is needed) or stratify by fewer levels."
            )

    rng = np.random.default_rng(seed)
    picked: list[np.ndarray] = []
    weights: list[np.ndarray] = []

    def allocate(positions: np.ndarray, depth: int, budget: int) -> None:
        if depth < len(level_codes):
            groups = _split(positions, level_codes[depth])
            sizes = pd.Series({key: len(pos) for key, pos in groups.items()})
            alloc = waterfill(sizes, budget)
            check_no_empty(alloc, repr(levels[depth]))
            for key, k in alloc.items():
                allocate(groups[key], depth + 1, int(k))
            return
        if type_codes is None:
            strata = {0: positions}
            alloc = pd.Series({0: budget})
        else:
            strata = _split(positions, type_codes)
            sizes = pd.Series({key: len(pos) for key, pos in strata.items()})
            alloc = _proportional_with_floor(sizes, budget, min_cells_per_type)
            check_no_empty(alloc, repr(cell_type_key))
        for key, k in alloc.items():
            pool = strata[key]
            picked.append(rng.choice(pool, size=int(k), replace=False))
            weights.append(np.full(int(k), len(pool) / int(k)))

    allocate(np.arange(len(obs)), 0, min(int(n_cells), len(obs)))
    positions = np.concatenate(picked)
    order = np.argsort(positions)
    return pd.Series(
        np.concatenate(weights)[order], index=obs.index[positions[order]], name=WEIGHT_KEY
    )


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

    See :func:`stratified_sample` for the parameters and when it raises. The weights are
    stored in ``obs[weight_key]``; use them to recover full-data composition from the subset.
    """
    if not adata.obs_names.is_unique:
        raise ValueError("obs_names must be unique; call adata.obs_names_make_unique() first")
    weights = stratified_sample(
        adata.obs,
        n_cells,
        levels=levels,
        cell_type_key=cell_type_key,
        min_cells_per_type=min_cells_per_type,
        seed=seed,
    )
    subset = adata[weights.index].copy()
    subset.obs[weight_key] = weights.to_numpy()
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
