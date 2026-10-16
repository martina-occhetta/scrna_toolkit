"""Dataset-level and grouped summaries returned as pandas objects (nothing is printed)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import cast

import numpy as np
import pandas as pd
from anndata import AnnData

from scrna_toolkit.metadata.schema import CANONICAL_FIELDS

DEFAULT_QC_COLS = ("total_counts", "n_genes_by_counts", "pct_counts_mt")


def summarize_dataset(
    adata: AnnData,
    *,
    sample_col: str = "sample_id",
    patient_col: str = "patient_id",
    cell_type_col: str = "cell_type",
    expected: Sequence[str] = CANONICAL_FIELDS,
) -> pd.Series:
    """One-line summary of a dataset. Counts for absent columns are ``NaN``.

    ``missing_metadata`` lists the ``expected`` fields absent from ``adata.obs``.
    """
    obs = cast(pd.DataFrame, adata.obs)

    def n_unique(col: str) -> float:
        return float(obs[col].nunique()) if col in obs else np.nan

    def cells_per(col: str, how: str) -> float:
        if col not in obs:
            return np.nan
        return float(obs[col].value_counts().agg(how))

    return pd.Series(
        {
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "n_samples": n_unique(sample_col),
            "n_patients": n_unique(patient_col),
            "n_cell_types": n_unique(cell_type_col),
            "median_cells_per_sample": cells_per(sample_col, "median"),
            "median_cells_per_patient": cells_per(patient_col, "median"),
            "minimum_cells_per_sample": cells_per(sample_col, "min"),
            "minimum_cells_per_patient": cells_per(patient_col, "min"),
            "missing_metadata": [f for f in expected if f not in obs],
        },
        dtype=object,
        name="value",
    )


def summarize_by(
    adata: AnnData,
    groupby: str,
    *,
    sample_col: str = "sample_id",
    patient_col: str = "patient_id",
    qc_cols: Sequence[str] = DEFAULT_QC_COLS,
) -> pd.DataFrame:
    """Per-group table: cell counts, fraction of cells, samples/patients, median QC metrics.

    Columns for absent metadata or QC metrics are omitted. Missing group labels form their own
    group so no cells silently disappear.
    """
    obs = cast(pd.DataFrame, adata.obs)
    if groupby not in obs:
        raise KeyError(f"Column {groupby!r} not in adata.obs")
    grouped = obs.groupby(groupby, observed=True, dropna=False)

    table = pd.DataFrame({"n_cells": grouped.size()})
    table["fraction_of_cells"] = table["n_cells"] / adata.n_obs
    for name, col in (("n_samples", sample_col), ("n_patients", patient_col)):
        if col in obs and col != groupby:
            table[name] = grouped[col].nunique()
    for col in qc_cols:
        if col in obs:
            table[f"median_{col}"] = grouped[col].median()
    return table.sort_values("n_cells", ascending=False)
