"""Inspect and stream AnnData-style HDF5 files, including legacy (anndata < 0.7) ones.

Files written by anndata before 0.7 (and some GEO depositions, which often use a plain ``.h5``
suffix) store the matrix as::

    X/        group, attrs h5sparse_format='csr', h5sparse_shape=(n_obs, n_vars)
      data, indices, indptr
    obs, var  compound datasets with an 'index' field; categoricals as integer codes with the
              categories in uns/<column>_categories

Current anndata cannot read them. :func:`read_h5_rows` reads them (and current CSR h5ad files)
row by row in chunks of stored entries, keeping only the selected rows, so a subset of a matrix
larger than memory can be loaded. :func:`h5_layout` tells these files apart from 10x ``.h5``
files.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from scrna_toolkit._utils import unique_gene_names
from scrna_toolkit.provenance import record_step

logger = logging.getLogger(__name__)

INT32_MAX = int(np.iinfo(np.int32).max)
DEFAULT_CHUNK_NNZ = 100_000_000  # stored entries per chunk: ~0.8-1.2 GB for data + indices


def _decode(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in ("O", "S"):
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in values])
    return values


def _attr_str(value: Any) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def h5_layout(path: str | Path) -> str:
    """Classify an HDF5 file without reading its matrix.

    Returns ``"10x_h5"`` (Cell Ranger v2/v3), ``"h5ad_legacy"`` (anndata < 0.7) or ``"h5ad"``
    (current anndata, whatever the suffix). Raises ``ValueError`` for anything else.
    """
    with h5py.File(path, "r") as h5:
        keys = set(h5.keys())
        if isinstance(h5.get("matrix"), h5py.Group) and "barcodes" in h5["matrix"]:
            return "10x_h5"  # Cell Ranger >= 3
        x = h5.get("X")
        if isinstance(x, h5py.Group):
            if "h5sparse_format" in x.attrs:
                return "h5ad_legacy"
            if "encoding-type" in x.attrs:
                return "h5ad"
        if isinstance(x, h5py.Dataset) and "obs" in keys:
            obs = h5["obs"]
            return "h5ad_legacy" if isinstance(obs, h5py.Dataset) else "h5ad"
        for key in keys:  # Cell Ranger 2: one group per genome
            group = h5[key]
            if isinstance(group, h5py.Group) and {"barcodes", "indptr"} <= set(group.keys()):
                return "10x_h5"
    raise ValueError(f"Unrecognised HDF5 layout in {path}; top-level keys {sorted(keys)[:20]}")


def h5_matrix_info(path: str | Path) -> dict[str, Any]:
    """Encoding, shape, nnz and dtypes of a sparse ``X`` without reading the matrix."""
    with h5py.File(path, "r") as h5:
        group = h5.get("X")
        if not isinstance(group, h5py.Group):
            raise TypeError(f"{path}: X is not a sparse group")
        attrs = dict(group.attrs)
        if "h5sparse_format" in attrs:
            fmt, shape = _attr_str(attrs["h5sparse_format"]), attrs["h5sparse_shape"]
            layout = "h5ad_legacy"
        else:
            fmt = _attr_str(attrs.get("encoding-type", "")).replace("_matrix", "")
            shape, layout = attrs["shape"], "h5ad"
        info: dict[str, Any] = {
            "layout": layout,
            "format": fmt,
            "n_obs": int(shape[0]),
            "n_vars": int(shape[1]),
            "nnz": int(group["data"].shape[0]),
            "data_dtype": str(group["data"].dtype),
            "indices_dtype": str(group["indices"].dtype),
            "indptr_dtype": str(group["indptr"].dtype),
        }
    info["density"] = info["nnz"] / max(1, info["n_obs"] * info["n_vars"])
    return info


def _legacy_frame(h5: h5py.File, name: str) -> pd.DataFrame:
    raw = h5[name][:]
    frame = pd.DataFrame({field: _decode(raw[field]) for field in raw.dtype.names})
    for column in list(frame.columns):
        categories_key = f"uns/{column}_categories"
        if column != "index" and categories_key in h5:
            categories = list(_decode(h5[categories_key][:]))
            frame[column] = pd.Categorical.from_codes(
                frame[column].astype(int), categories=categories
            )
    frame = frame.set_index("index")
    frame.index = frame.index.astype(str)
    frame.index.name = None
    return frame


def _read_elem(elem: Any) -> Any:
    try:
        from anndata.io import read_elem
    except ImportError:  # anndata < 0.11
        from anndata.experimental import read_elem  # type: ignore[no-redef]
    return read_elem(elem)


def read_h5_axes(path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """``obs`` and ``var`` tables of an AnnData HDF5 file (legacy or current), without ``X``."""
    with h5py.File(path, "r") as h5:
        if isinstance(h5["obs"], h5py.Dataset):
            return _legacy_frame(h5, "obs"), _legacy_frame(h5, "var")
        return _read_elem(h5["obs"]), _read_elem(h5["var"])


def _row_chunks(indptr: np.ndarray, max_nnz: int) -> Iterator[tuple[int, int]]:
    """Consecutive row ranges holding at most ``max_nnz`` entries (or a single larger row)."""
    n_rows, start = indptr.shape[0] - 1, 0
    while start < n_rows:
        stop = int(np.searchsorted(indptr, indptr[start] + max_nnz, side="right")) - 1
        stop = min(n_rows, max(stop, start + 1))
        yield start, stop
        start = stop


def _open_csr(h5: h5py.File, path: str | Path) -> tuple[h5py.Group, np.ndarray, int, int]:
    group = h5.get("X")
    if not isinstance(group, h5py.Group):
        raise TypeError(f"{path}: X is dense; read_h5_rows streams sparse CSR matrices only")
    attrs = group.attrs
    if "h5sparse_format" in attrs:
        fmt, shape = _attr_str(attrs["h5sparse_format"]), attrs["h5sparse_shape"]
    else:
        fmt, shape = (
            _attr_str(attrs.get("encoding-type", "")).replace("_matrix", ""),
            attrs["shape"],
        )
    if fmt != "csr":
        raise TypeError(f"{path}: X is stored as {fmt!r}; rows can only be streamed from CSR")
    indptr = group["indptr"][:].astype(np.int64)
    return group, indptr, int(shape[0]), int(shape[1])


def _row_mask(rows: Any, n_rows: int) -> np.ndarray:
    if rows is None:
        return np.ones(n_rows, dtype=bool)
    rows = np.asarray(rows)
    if rows.dtype == bool:
        if rows.shape != (n_rows,):
            raise ValueError(f"Boolean row mask has length {rows.shape[0]}, matrix has {n_rows}")
        return rows
    if rows.dtype.kind not in "iu":
        raise TypeError("rows must be a boolean mask or integer row positions")
    if rows.size and (rows.min() < 0 or rows.max() >= n_rows):
        raise IndexError(f"Row positions must be in [0, {n_rows})")
    if np.unique(rows).size != rows.size:
        raise ValueError("Row positions contain duplicates")
    mask = np.zeros(n_rows, dtype=bool)
    mask[rows] = True
    return mask


def read_h5_rows(
    path: str | Path,
    rows: Any = None,
    *,
    max_nnz_per_chunk: int = DEFAULT_CHUNK_NNZ,
    dtype: Any = np.float32,
) -> AnnData:
    """Load selected rows of an AnnData HDF5 file (legacy or current CSR) in bounded memory.

    Parameters
    ----------
    rows
        Boolean mask over all rows, integer row positions, or ``None`` for all rows. Rows are
        returned in file order, whatever the order of the positions given.
    max_nnz_per_chunk
        Stored entries read from disk at a time. Peak memory is the output matrix plus one
        chunk of ``data`` and ``indices``.
    dtype
        dtype of the returned values (``float32`` by default, as elsewhere in this package).

    The output arrays are pre-allocated from ``indptr``, with ``int32`` indices when they fit,
    and the file is read in one sequential pass. ``obs`` and ``var`` are subset to match.
    """
    adata, n_rows = _read_rows(path, rows, max_nnz_per_chunk=max_nnz_per_chunk, dtype=dtype)
    record_step(adata, "read_h5_rows", {"path": str(path), "n_rows_total": n_rows})
    return adata


def _read_rows(
    path: str | Path, rows: Any, *, max_nnz_per_chunk: int, dtype: Any
) -> tuple[AnnData, int]:
    obs, var = read_h5_axes(path)
    with h5py.File(path, "r") as h5:
        group, indptr, n_rows, n_cols = _open_csr(h5, path)
        if indptr.shape[0] != n_rows + 1:
            raise ValueError(f"{path}: indptr length {indptr.shape[0]} != n_obs + 1")
        keep = _row_mask(rows, n_rows)
        row_nnz = np.diff(indptr)
        out_nnz = int(row_nnz[keep].sum())
        index_dtype = np.int32 if max(out_nnz, n_cols) <= INT32_MAX else np.int64
        data = np.empty(out_nnz, dtype=dtype)
        indices = np.empty(out_nnz, dtype=index_dtype)
        out_indptr = np.zeros(int(keep.sum()) + 1, dtype=index_dtype)
        np.cumsum(row_nnz[keep], out=out_indptr[1:])

        written = 0
        for start, stop in _row_chunks(indptr, max_nnz_per_chunk):
            if not keep[start:stop].any():
                continue
            e0, e1 = int(indptr[start]), int(indptr[stop])
            entry_keep = np.repeat(keep[start:stop], row_nnz[start:stop])
            n_new = int(entry_keep.sum())
            data[written : written + n_new] = group["data"][e0:e1][entry_keep]
            indices[written : written + n_new] = group["indices"][e0:e1][entry_keep]
            written += n_new
            logger.debug("Read rows %d-%d of %d", start, stop, n_rows)
    if written != out_nnz:
        raise RuntimeError(f"Read {written} stored entries, expected {out_nnz}")

    matrix = sp.csr_matrix((data, indices, out_indptr), shape=(int(keep.sum()), n_cols))
    var.index = unique_gene_names(pd.Index(var.index.astype(str)))
    adata = AnnData(X=matrix, obs=obs.iloc[keep].copy(), var=var)
    logger.info(
        "Read %d of %d rows (%d stored entries) from %s", adata.n_obs, n_rows, out_nnz, path
    )
    return adata, n_rows
