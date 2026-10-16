"""Load common single-cell count formats into AnnData (cells x genes).

Supported formats: ``h5ad``, ``10x_h5`` (Cell Ranger .h5), ``10x_mtx`` (matrix.mtx +
barcodes.tsv + features.tsv/genes.tsv, optionally gzipped and/or with a GEO-style filename
prefix) and dense ``csv``/``tsv`` count tables.

Detection is deliberately conservative: anything ambiguous raises
:class:`FormatDetectionError` instead of guessing.
"""

from __future__ import annotations

import logging
import warnings
from collections.abc import Mapping
from pathlib import Path
from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from anndata import AnnData

from scrna_toolkit.provenance import record_step

logger = logging.getLogger(__name__)

Format = Literal["h5ad", "10x_h5", "10x_mtx", "csv", "tsv"]
Orientation = Literal["cells_x_genes", "genes_x_cells"]
FORMATS: tuple[str, ...] = ("h5ad", "10x_h5", "10x_mtx", "csv", "tsv")


class FormatDetectionError(ValueError):
    """Raised when the input format cannot be determined unambiguously."""


def _strip_gz(name: str) -> str:
    return name[:-3] if name.endswith(".gz") else name


def _find_one(directory: Path, names: list[str]) -> Path | None:
    for name in names:
        for candidate in (directory / name, directory / f"{name}.gz"):
            if candidate.is_file():
                return candidate
    return None


def _find_10x_files(directory: Path, prefix: str | None = None) -> tuple[Path, Path, Path]:
    """Return (matrix, barcodes, features) paths for a 10x-style directory."""
    if prefix is None:
        matrices = sorted(
            p for p in directory.iterdir() if _strip_gz(p.name).endswith("matrix.mtx")
        )
        if not matrices:
            contents = sorted(p.name for p in directory.iterdir())[:20]
            raise FormatDetectionError(
                f"No *matrix.mtx(.gz) file in {directory}. Contents: {contents}"
            )
        if len(matrices) > 1:
            prefixes = [_strip_gz(p.name)[: -len("matrix.mtx")] for p in matrices]
            raise FormatDetectionError(
                f"Several matrices in {directory} (prefixes {prefixes}); pass prefix=... "
                "or use load_samples() to load and concatenate them."
            )
        prefix = _strip_gz(matrices[0].name)[: -len("matrix.mtx")]

    matrix = _find_one(directory, [f"{prefix}matrix.mtx"])
    barcodes = _find_one(directory, [f"{prefix}barcodes.tsv"])
    features = _find_one(directory, [f"{prefix}features.tsv", f"{prefix}genes.tsv"])
    missing = [
        n
        for n, p in [("matrix.mtx", matrix), ("barcodes.tsv", barcodes), ("features.tsv", features)]
        if p is None
    ]
    if missing or matrix is None or barcodes is None or features is None:
        raise FormatDetectionError(
            f"Missing 10x files {missing} (prefix={prefix!r}) in {directory}"
        )
    return matrix, barcodes, features


def detect_format(path: str | Path, prefix: str | None = None) -> str:
    """Infer the format of ``path``; raise :class:`FormatDetectionError` if ambiguous."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.is_dir():
        _find_10x_files(path, prefix)
        return "10x_mtx"
    name = _strip_gz(path.name.lower())
    if name.endswith(".h5ad"):
        return "h5ad"
    if name.endswith(".h5"):
        return "10x_h5"
    if name.endswith("matrix.mtx"):
        return "10x_mtx"
    if name.endswith(".csv"):
        return "csv"
    if name.endswith(".tsv"):
        return "tsv"
    raise FormatDetectionError(
        f"Cannot infer format of {path.name!r}; pass format= one of {FORMATS}."
    )


def _unique_gene_names(names: pd.Index) -> pd.Index:
    if names.is_unique:
        return names
    n_dup = int(names.duplicated().sum())
    logger.warning("%d duplicated gene names made unique (gene_ids are preserved)", n_dup)
    return ad.utils.make_index_unique(names)


def _read_10x_mtx(path: Path, prefix: str | None, gex_only: bool) -> AnnData:
    if path.is_file():  # a matrix.mtx file was given directly
        prefix = _strip_gz(path.name)[: -len("matrix.mtx")]
        path = path.parent
    matrix_path, barcodes_path, features_path = _find_10x_files(path, prefix)

    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, dtype=str)[0].to_numpy()
    features = pd.read_csv(features_path, sep="\t", header=None, dtype=str)
    with warnings.catch_warnings():  # scipy>=1.15 warns about its sparse return type
        warnings.filterwarnings("ignore", "The default value for `spmatrix`")
        matrix = sp.csr_matrix(scipy.io.mmread(matrix_path), dtype=np.float32)

    n_bc, n_feat = len(barcodes), len(features)
    if matrix.shape == (n_feat, n_bc) and n_bc != n_feat:
        matrix = matrix.T.tocsr()  # standard 10x layout is genes x cells
    elif matrix.shape != (n_bc, n_feat) or n_bc == n_feat:
        raise FormatDetectionError(
            f"Matrix shape {matrix.shape} does not identify orientation for "
            f"{n_bc} barcodes and {n_feat} features"
        )

    var = pd.DataFrame({"gene_ids": features[0].to_numpy()})
    if features.shape[1] > 1:
        var["gene_symbols"] = features[1].to_numpy()
    if features.shape[1] > 2:
        var["feature_types"] = features[2].to_numpy()
    names = pd.Index(var.get("gene_symbols", var["gene_ids"]).to_numpy())
    var.index = _unique_gene_names(names)

    adata = AnnData(X=matrix, obs=pd.DataFrame(index=barcodes), var=var)
    if gex_only and "feature_types" in adata.var:
        keep = (adata.var["feature_types"] == "Gene Expression").to_numpy()
        if not keep.all():
            logger.info("Keeping %d Gene Expression features of %d", keep.sum(), len(keep))
            adata = adata[:, keep].copy()
    return adata


def _read_delimited(path: Path, sep: str, orientation: Orientation | None) -> AnnData:
    if orientation not in ("cells_x_genes", "genes_x_cells"):
        raise FormatDetectionError(
            "Orientation of delimited count tables cannot be determined reliably; pass "
            "orientation='cells_x_genes' or orientation='genes_x_cells'."
        )
    df = pd.read_csv(path, sep=sep, index_col=0)
    non_numeric = [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]
    if non_numeric:
        raise ValueError(f"Non-numeric columns in count table {path.name}: {non_numeric[:5]}")
    if orientation == "genes_x_cells":
        df = df.T
    return AnnData(
        X=sp.csr_matrix(df.to_numpy(dtype=np.float32)),
        obs=pd.DataFrame(index=df.index.astype(str)),
        var=pd.DataFrame(index=df.columns.astype(str)),
    )


def load_scrna(
    path: str | Path,
    format: str | None = None,
    *,
    orientation: Orientation | None = None,
    prefix: str | None = None,
    gex_only: bool = True,
) -> AnnData:
    """Load a single-cell count matrix as a cells x genes AnnData.

    Parameters
    ----------
    path
        File or (for 10x MTX) directory.
    format
        One of ``h5ad``, ``10x_h5``, ``10x_mtx``, ``csv``, ``tsv``; inferred if omitted.
    orientation
        Required for ``csv``/``tsv``: ``"cells_x_genes"`` or ``"genes_x_cells"``.
    prefix
        Filename prefix for 10x MTX files (e.g. ``"GSM123_"``) when a directory holds several.
    gex_only
        For 10x inputs with feature types, keep only Gene Expression features.

    Gene symbols become ``var_names`` (made unique) and Ensembl IDs are kept in
    ``var["gene_ids"]`` when both are available.
    """
    path = Path(path)
    fmt = format or detect_format(path, prefix)
    if fmt not in FORMATS:
        raise ValueError(f"Unknown format {fmt!r}; expected one of {FORMATS}")

    if fmt == "h5ad":
        adata = ad.read_h5ad(path)
    elif fmt == "10x_h5":
        import scanpy as sc

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Variable names are not unique")
            adata = sc.read_10x_h5(path, gex_only=gex_only)
        adata.var_names = _unique_gene_names(adata.var_names)
    elif fmt == "10x_mtx":
        adata = _read_10x_mtx(path, prefix, gex_only)
    else:
        adata = _read_delimited(path, "," if fmt == "csv" else "\t", orientation)

    if not adata.var_names.is_unique:
        adata.var_names = _unique_gene_names(adata.var_names)
    if not adata.obs_names.is_unique:
        logger.warning("Cell barcodes are not unique; use load_samples() or obs_names_make_unique")

    logger.info("Loaded %s (%s): %d cells x %d genes", path, fmt, adata.n_obs, adata.n_vars)
    record_step(adata, "load_scrna", {"path": str(path), "format": fmt})
    return adata


def load_samples(
    paths: Mapping[str, str | Path],
    *,
    sample_key: str = "sample_id",
    **load_kwargs,
) -> AnnData:
    """Load one matrix per sample and concatenate them.

    ``paths`` maps sample IDs to inputs for :func:`load_scrna`. Cell names are suffixed
    with ``-<sample_id>`` so barcodes stay unique; only genes shared by all samples are kept.
    """
    adatas = {sample: load_scrna(path, **load_kwargs) for sample, path in paths.items()}
    gene_sets = {frozenset(a.var_names) for a in adatas.values()}
    if len(gene_sets) > 1:
        logger.warning("Samples have different gene sets; keeping the intersection")
    combined = ad.concat(adatas, label=sample_key, index_unique="-", join="inner", merge="same")
    combined.obs[sample_key] = combined.obs[sample_key].astype("category")
    record_step(
        combined,
        "load_samples",
        {"paths": {k: str(v) for k, v in paths.items()}, "sample_key": sample_key},
    )
    return combined
