"""Report problems in ``adata.obs`` metadata without modifying or rejecting the data."""

from __future__ import annotations

import logging
import re
from collections.abc import Sequence
from typing import cast

import pandas as pd
from anndata import AnnData

from scrna_toolkit.metadata.schema import CANONICAL_FIELDS, SAMPLE_LEVEL_FIELDS

logger = logging.getLogger(__name__)

ISSUE_COLUMNS = ["check", "column", "severity", "message"]
_MISSING_LIKE = {"", "na", "nan", "none", "null"}


def _normalise_label(value: object) -> str:
    return re.sub(r"[\s_\-]+", " ", str(value).strip().lower())


def validate_metadata(
    adata: AnnData,
    expected: Sequence[str] = CANONICAL_FIELDS,
    *,
    sample_col: str = "sample_id",
) -> pd.DataFrame:
    """Return a table of metadata issues (one row per issue; empty if none).

    Checks: absent expected fields (``info``), duplicated cell IDs (``error``), missing or
    placeholder values in present fields, sample-level fields taking several values within one
    sample, and labels that differ only by case/whitespace/separators (``warning``).
    """
    obs = cast(pd.DataFrame, adata.obs)
    issues: list[tuple[str, str, str, str]] = []

    for field in expected:
        if field not in obs:
            issues.append(("missing_field", field, "info", f"{field!r} not in adata.obs"))

    n_dup = int(adata.obs_names.duplicated().sum())
    if n_dup:
        issues.append(("duplicated_cell_ids", "obs_names", "error", f"{n_dup} duplicated IDs"))

    present = [f for f in expected if f in obs]
    for field in present:
        values = obs[field]
        n_missing = int(values.isna().sum())
        n_missing += int(
            values.dropna().astype(str).str.strip().str.lower().isin(_MISSING_LIKE).sum()
        )
        if n_missing:
            issues.append(
                ("missing_values", field, "warning", f"{n_missing} cells without a value")
            )

    if sample_col in obs:
        for field in present:
            if field == sample_col or field not in SAMPLE_LEVEL_FIELDS:
                continue
            n_values = obs.groupby(sample_col, observed=True)[field].nunique()
            bad = n_values[n_values > 1]
            if len(bad):
                issues.append(
                    (
                        "one_to_many",
                        field,
                        "warning",
                        f"{len(bad)} samples map to several {field!r} values: "
                        f"{list(bad.index[:5])}",
                    )
                )

    for field in present:
        dtype = obs[field].dtype
        if not (pd.api.types.is_string_dtype(dtype) or isinstance(dtype, pd.CategoricalDtype)):
            continue
        groups: dict[str, set[str]] = {}
        for value in obs[field].dropna().unique():
            groups.setdefault(_normalise_label(value), set()).add(str(value))
        clashes = [sorted(v) for v in groups.values() if len(v) > 1]
        if clashes:
            issues.append(
                ("inconsistent_labels", field, "warning", f"Near-duplicate labels: {clashes[:5]}")
            )

    report = pd.DataFrame(issues, columns=ISSUE_COLUMNS)
    n_problems = int((report["severity"] != "info").sum())
    logger.info(
        "Metadata validation: %d problems, %d absent fields", n_problems, len(report) - n_problems
    )
    return report
