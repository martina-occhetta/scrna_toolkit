"""Storage types for ``obs`` columns."""

from __future__ import annotations

import logging
from collections.abc import Iterable
from typing import cast

import pandas as pd
from anndata import AnnData

from scrna_toolkit.provenance import record_step

logger = logging.getLogger(__name__)


def _is_string_column(values: pd.Series) -> bool:
    if pd.api.types.is_string_dtype(values.dtype) and values.dtype != object:
        return True  # pandas "string"/"str" dtypes
    return values.dtype == object and pd.api.types.infer_dtype(values, skipna=True) == "string"


def categorize_obs(
    adata: AnnData,
    *,
    max_categories: int = 5000,
    columns: Iterable[str] | None = None,
) -> AnnData:
    """Store string ``obs`` columns as categoricals (smaller files, faster ``groupby``).

    A column is converted if it holds only strings (missing values allowed) and has at most
    ``max_categories`` distinct values. Columns that are already categorical, numeric, boolean
    or of mixed types are left unchanged. ``columns`` restricts the candidates (default: all).
    Values and missing values are unchanged. Modifies ``adata.obs`` in place and returns
    ``adata``.
    """
    obs = cast(pd.DataFrame, adata.obs)
    candidates = list(obs.columns if columns is None else columns)
    absent = [c for c in candidates if c not in obs]
    if absent:
        raise KeyError(f"Columns not in adata.obs: {absent}")

    converted, too_many = [], []
    for column in candidates:
        values = obs[column]
        if not _is_string_column(values):
            continue
        if values.nunique(dropna=True) > max_categories:
            too_many.append(column)
            continue
        adata.obs[column] = values.astype("category")
        converted.append(column)
    if too_many:
        logger.info(
            "Kept %d columns with > %d values as strings: %s",
            len(too_many),
            max_categories,
            too_many,
        )

    record_step(
        adata,
        "categorize_obs",
        {"max_categories": max_categories, "converted": converted, "kept_as_strings": too_many},
    )
    return adata
