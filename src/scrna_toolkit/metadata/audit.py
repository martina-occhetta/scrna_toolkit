"""Describe ``obs`` metadata without modifying it: column levels and per-column summaries."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import cast

import pandas as pd
from anndata import AnnData

DEFAULT_LEVELS: tuple[str, ...] = ("patient_id", "sample_id")


def _obs(data: AnnData | pd.DataFrame) -> pd.DataFrame:
    return cast(pd.DataFrame, data.obs) if isinstance(data, AnnData) else data


def infer_field_levels(
    data: AnnData | pd.DataFrame,
    levels: Sequence[str] = DEFAULT_LEVELS,
    columns: Iterable[str] | None = None,
) -> pd.Series:
    """Return the coarsest level at which each column takes a single value.

    Parameters
    ----------
    data
        An ``AnnData`` (its ``obs`` is used) or a cell-level ``DataFrame``.
    levels
        Columns identifying nested units, coarsest first, e.g.
        ``("patient_id", "sample_id", "library_id")``. Levels absent from the table are skipped.
    columns
        Columns to classify (default: all). Absent columns raise a ``KeyError``.

    Each column is labelled ``"constant"`` if it has one value across all cells, else the first
    level whose every group has a single value (a level column gets its own name at the
    latest), else ``"cell"``. Missing values count as a value, and cells with a missing level
    key form a group of their own, so a column is never placed at a level by ignoring cells.
    """
    obs = _obs(data)
    cols = list(obs.columns if columns is None else columns)
    absent = [c for c in cols if c not in obs]
    if absent:
        raise KeyError(f"Columns not in obs: {absent}")

    result = pd.Series("cell", index=pd.Index(cols, name="column"), name="level", dtype=object)
    n_unique = obs[cols].nunique(dropna=False)
    result[n_unique.index[n_unique <= 1]] = "constant"
    undecided = [c for c in cols if n_unique[c] > 1]
    for key in levels:
        if not undecided:
            break
        if key not in obs:
            continue
        if key in undecided:
            result[key] = key
            undecided.remove(key)
        if not undecided:
            break
        per_group = obs.groupby(key, observed=True, dropna=False)[undecided].nunique(dropna=False)
        single = set(per_group.columns[(per_group <= 1).all()])
        result[list(single)] = key
        undecided = [c for c in undecided if c not in single]
    return result


def _examples(values: pd.Series, n_examples: int) -> str:
    present = values.dropna()
    if present.empty:
        return ""
    if pd.api.types.is_numeric_dtype(present) and not pd.api.types.is_bool_dtype(present):
        return f"range {present.min():.3g} – {present.max():.3g}"
    return ", ".join(present.astype(str).value_counts().head(n_examples).index)


def metadata_audit_table(
    data: AnnData | pd.DataFrame,
    columns: Iterable[str] | None = None,
    *,
    levels: Sequence[str] = DEFAULT_LEVELS,
    descriptions: Mapping[str, str] | None = None,
    n_examples: int = 4,
) -> pd.DataFrame:
    """One row per ``obs`` column: level, dtype, unique and missing counts, example values.

    ``level`` comes from :func:`infer_field_levels` with ``levels``. Examples are the
    ``n_examples`` most frequent values, or the range for numeric columns. With
    ``descriptions`` (``{column: text}``), a ``description`` column is added (empty where a
    column is not described). Only ``NaN``/``None`` count as missing; placeholder strings such
    as ``"NA"`` are reported by :func:`validate_metadata`.
    """
    obs = _obs(data)
    cols = list(obs.columns if columns is None else columns)
    field_levels = infer_field_levels(obs, levels=levels, columns=cols)
    table = pd.DataFrame(
        {
            "column": cols,
            "level": field_levels.to_numpy(),
            "dtype": [str(obs[c].dtype) for c in cols],
            "n_unique": [int(obs[c].nunique(dropna=True)) for c in cols],
            "n_missing": [int(obs[c].isna().sum()) for c in cols],
            "examples": [_examples(obs[c], n_examples) for c in cols],
        }
    )
    if descriptions is not None:
        table["description"] = [descriptions.get(c, "") for c in cols]
    return table
