"""Per-cell QC metrics. Calculation never filters; see :mod:`scrna_toolkit.qc.filtering`."""

from __future__ import annotations

import logging
import re
from collections.abc import Sequence

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from scrna_toolkit._utils import looks_like_counts
from scrna_toolkit.provenance import record_step

logger = logging.getLogger(__name__)


def calculate_qc_metrics(
    adata: AnnData,
    mitochondrial_prefix: str = "MT-",
    *,
    mitochondrial_genes: Sequence[str] | None = None,
    ribosomal_prefixes: Sequence[str] | None = None,
    hemoglobin_regex: str | None = None,
    gene_symbol_col: str | None = None,
    layer: str | None = None,
) -> AnnData:
    """Add standard QC metrics to ``adata.obs`` (in place) and return ``adata``.

    Adds ``total_counts``, ``n_genes_by_counts``, ``total_counts_mt`` and ``pct_counts_mt``
    (a percentage, 0-100), plus ``*_ribo`` / ``*_hb`` columns when those gene sets are
    requested. Gene flags are stored in ``var["mt"]``, ``var["ribo"]``, ``var["hb"]``.

    Parameters
    ----------
    mitochondrial_prefix
        Prefix of mitochondrial gene symbols (``"MT-"`` human, ``"mt-"`` mouse).
    mitochondrial_genes
        Explicit list of mitochondrial genes; overrides the prefix.
    ribosomal_prefixes
        e.g. ``("RPS", "RPL")``. Not computed by default.
    hemoglobin_regex
        e.g. ``r"^HB[^(P)]"``. Not computed by default.
    gene_symbol_col
        ``var`` column holding gene symbols, if ``var_names`` are IDs.
    layer
        Layer holding raw counts; defaults to ``X``.
    """
    matrix = adata.layers[layer] if layer else adata.X
    if not looks_like_counts(matrix):
        logger.warning("QC metrics expect raw counts but the matrix is not integer-valued")

    symbols = pd.Index(adata.var[gene_symbol_col] if gene_symbol_col else adata.var_names).astype(
        str
    )
    if mitochondrial_genes is not None:
        adata.var["mt"] = symbols.isin(list(mitochondrial_genes))
    else:
        adata.var["mt"] = symbols.str.startswith(mitochondrial_prefix)
    qc_vars = ["mt"]
    if ribosomal_prefixes:
        adata.var["ribo"] = symbols.str.startswith(tuple(ribosomal_prefixes))
        qc_vars.append("ribo")
    if hemoglobin_regex:
        pattern = re.compile(hemoglobin_regex)
        adata.var["hb"] = np.array([bool(pattern.search(s)) for s in symbols])
        qc_vars.append("hb")

    for name in qc_vars:
        n = int(adata.var[name].sum())
        if n == 0:
            logger.warning("No %s genes found; pct_counts_%s will be 0", name, name)
        else:
            logger.info("Found %d %s genes", n, name)

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=qc_vars, percent_top=None, log1p=False, inplace=True, layer=layer
    )
    record_step(
        adata,
        "calculate_qc_metrics",
        {
            "mitochondrial_prefix": mitochondrial_prefix,
            "mitochondrial_genes": list(mitochondrial_genes) if mitochondrial_genes else None,
            "ribosomal_prefixes": list(ribosomal_prefixes) if ribosomal_prefixes else None,
            "hemoglobin_regex": hemoglobin_regex,
            "gene_symbol_col": gene_symbol_col,
            "layer": layer,
            "n_genes_flagged": {name: int(adata.var[name].sum()) for name in qc_vars},
        },
    )
    return adata
