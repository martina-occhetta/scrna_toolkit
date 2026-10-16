"""Map dataset-specific ``obs`` columns onto canonical names without losing the originals."""

from __future__ import annotations

import logging
from collections import Counter
from collections.abc import Mapping

import pandas as pd
from anndata import AnnData

from scrna_toolkit.metadata.schema import CANONICAL_FIELDS
from scrna_toolkit.provenance import record_step

logger = logging.getLogger(__name__)


def harmonise_obs(
    adata: AnnData,
    mapping: Mapping[str, str],
    *,
    overwrite: bool = False,
    strict: bool = False,
) -> AnnData:
    """Copy ``obs`` columns to canonical names, keeping the original columns.

    Parameters
    ----------
    mapping
        ``{source_column: target_column}``, e.g. ``{"orig.ident": "sample_id"}``.
    overwrite
        Allow replacing an existing target column.
    strict
        Raise if a source column is missing (default: log a warning and skip it).

    Modifies ``adata.obs`` in place and returns ``adata``. Missing metadata is never filled in.
    """
    duplicated = [t for t, n in Counter(mapping.values()).items() if n > 1]
    if duplicated:
        raise ValueError(f"Several source columns map to the same target: {duplicated}")

    applied: dict[str, str] = {}
    for source, target in mapping.items():
        if source not in adata.obs:
            if strict:
                raise KeyError(f"Source column {source!r} not in adata.obs")
            logger.warning("Source column %r not in adata.obs; skipping", source)
            continue
        if source == target:
            continue
        if target in adata.obs and not overwrite:
            raise ValueError(
                f"Target column {target!r} already exists; pass overwrite=True to replace it"
            )
        if target not in CANONICAL_FIELDS:
            logger.info("Target column %r is not a canonical field", target)
        adata.obs[target] = adata.obs[source].copy()
        applied[source] = target

    record_step(adata, "harmonise_obs", {"mapping": applied, "overwrite": overwrite})
    return adata


def merge_sample_metadata(
    adata: AnnData,
    metadata: pd.DataFrame,
    *,
    on: str = "sample_id",
    overwrite: bool = False,
) -> AnnData:
    """Add sample-level columns (one row per sample) to every cell of that sample.

    ``metadata`` must contain ``on`` as a column with unique values. Cells whose sample is
    absent from the table get missing values (reported, never imputed). Modifies in place.
    """
    if on not in adata.obs:
        raise KeyError(f"Column {on!r} not in adata.obs; harmonise it first")
    if on not in metadata:
        raise KeyError(f"Column {on!r} not in sample metadata table")
    if metadata[on].duplicated().any():
        dups = metadata.loc[metadata[on].duplicated(), on].unique().tolist()
        raise ValueError(f"Sample metadata has duplicated {on!r} values: {dups[:10]}")

    table = metadata.set_index(on)
    clashes = [c for c in table.columns if c in adata.obs and not overwrite]
    if clashes:
        raise ValueError(f"Columns already in adata.obs: {clashes}; pass overwrite=True")

    keys = adata.obs[on].astype(str)
    table.index = table.index.astype(str)
    unmatched = sorted(set(keys) - set(table.index))
    if unmatched:
        logger.warning("%d samples have no metadata row: %s", len(unmatched), unmatched[:10])
    for column in table.columns:
        adata.obs[column] = keys.map(table[column]).to_numpy()

    record_step(
        adata,
        "merge_sample_metadata",
        {"on": on, "columns": list(table.columns), "unmatched_samples": unmatched},
    )
    return adata
