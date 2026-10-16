"""Legacy (anndata < 0.7) HDF5 files and streamed row subsets."""

import h5py
import numpy as np
import pytest
import scipy.sparse as sp

import scrna_toolkit as sct
from scrna_toolkit.io import FormatDetectionError, h5_layout, h5_matrix_info
from scrna_toolkit.io.h5 import _row_chunks

N_OBS, N_VARS = 50, 12
SAMPLES = ["S0", "S1", "S2"]


def _matrix():
    rng = np.random.default_rng(1)
    dense = rng.integers(0, 4, size=(N_OBS, N_VARS)).astype(np.float32)
    dense[rng.random(dense.shape) < 0.5] = 0
    dense[[3, 17, 40]] = 0  # empty rows
    return sp.csr_matrix(dense)


def write_legacy_h5(path, matrix=None, fmt="csr"):
    """The layout of e.g. GSE180661_matrix.h5: h5sparse X, compound obs/var, coded categoricals."""
    matrix = (_matrix() if matrix is None else matrix).asformat(fmt)
    obs = np.zeros(N_OBS, dtype=[("index", "S12"), ("sample", "i1"), ("n_counts", "f4")])
    obs["index"] = [f"cell{i}".encode() for i in range(N_OBS)]
    obs["sample"] = np.arange(N_OBS) % 3
    obs["n_counts"] = np.asarray(_matrix().sum(1)).ravel()
    var = np.zeros(N_VARS, dtype=[("index", "S8"), ("gene_ids", "S12"), ("feature_types", "i1")])
    var["index"] = [b"MT-CO1", b"ACTB", b"ACTB"] + [f"G{i}".encode() for i in range(N_VARS - 3)]
    var["gene_ids"] = [f"ENSG{i:05d}".encode() for i in range(N_VARS)]
    with h5py.File(path, "w") as h5:
        x = h5.create_group("X")
        x.attrs["h5sparse_format"] = fmt
        x.attrs["h5sparse_shape"] = np.array(matrix.shape)
        x["data"], x["indices"] = matrix.data, matrix.indices.astype(np.int64)
        x["indptr"] = matrix.indptr.astype(np.int64)
        h5["obs"], h5["var"] = obs, var
        h5["uns/sample_categories"] = np.array(SAMPLES, dtype="S")
        h5["uns/feature_types_categories"] = np.array([b"Gene Expression"])
    return path


@pytest.fixture
def legacy(tmp_path):
    return write_legacy_h5(tmp_path / "GSE0_matrix.h5")


def test_legacy_h5_is_detected_and_loaded(legacy):
    assert sct.detect_format(legacy) == "h5ad_legacy"
    assert h5_layout(legacy) == "h5ad_legacy"

    adata = sct.load_scrna(legacy)

    assert adata.shape == (N_OBS, N_VARS)
    np.testing.assert_array_equal(adata.X.toarray(), _matrix().toarray())
    assert adata.X.dtype == np.float32 and adata.X.indices.dtype == np.int32
    assert list(adata.obs_names[:2]) == ["cell0", "cell1"]
    assert list(adata.obs["sample"].cat.categories) == SAMPLES
    assert adata.obs["sample"].iloc[4] == "S1"
    assert (adata.var["feature_types"] == "Gene Expression").all()
    assert adata.var_names.is_unique and adata.var["gene_ids"].iloc[2] == "ENSG00002"
    assert sct.get_history(adata)[-1]["params"]["format"] == "h5ad_legacy"


def test_current_anndata_file_with_h5_suffix_is_read_as_h5ad(tmp_path, adata):
    path = tmp_path / "matrix.h5"
    adata.write_h5ad(path)
    assert sct.detect_format(path) == "h5ad"
    assert sct.load_scrna(path).shape == adata.shape


def test_unrecognised_h5_files_raise(tmp_path):
    other = tmp_path / "other.h5"
    with h5py.File(other, "w") as h5:
        h5["something"] = np.arange(3)
    with pytest.raises(FormatDetectionError, match="something"):
        sct.detect_format(other)
    not_hdf5 = tmp_path / "broken.h5"
    not_hdf5.write_bytes(b"not an hdf5 file")
    with pytest.raises(FormatDetectionError, match="format="):
        sct.load_scrna(not_hdf5)


@pytest.mark.parametrize("chunk", [1, 7, 10_000])
def test_read_h5_rows_matches_in_memory_subset(legacy, chunk):
    keep = np.zeros(N_OBS, dtype=bool)
    keep[[0, 3, 4, 17, 25, 49]] = True  # includes empty rows

    subset = sct.read_h5_rows(legacy, keep, max_nnz_per_chunk=chunk)

    np.testing.assert_array_equal(subset.X.toarray(), _matrix()[keep].toarray())
    assert list(subset.obs_names) == [f"cell{i}" for i in np.flatnonzero(keep)]
    assert subset.n_vars == N_VARS
    assert sct.get_history(subset)[-1]["step"] == "read_h5_rows"


def test_read_h5_rows_positions_are_returned_in_file_order(legacy):
    subset = sct.read_h5_rows(legacy, [30, 2, 11])
    assert list(subset.obs_names) == ["cell2", "cell11", "cell30"]
    np.testing.assert_array_equal(subset.X.toarray(), _matrix()[[2, 11, 30]].toarray())
    assert sct.read_h5_rows(legacy, np.array([], dtype=int)).shape == (0, N_VARS)


def test_read_h5_rows_on_current_h5ad(tmp_path, adata):
    path = tmp_path / "x.h5ad"
    adata.write_h5ad(path)
    assert h5_matrix_info(path)["layout"] == "h5ad"

    subset = sct.read_h5_rows(
        path, adata.obs["orig.ident"].eq("S1").to_numpy(), max_nnz_per_chunk=50
    )

    expected = adata[adata.obs["orig.ident"] == "S1"]
    np.testing.assert_array_equal(subset.X.toarray(), expected.X.toarray())
    assert list(subset.obs_names) == list(expected.obs_names)
    assert list(subset.obs["patient"]) == list(expected.obs["patient"])


def test_read_h5_axes_and_matrix_info(legacy):
    obs, var = sct.read_h5_axes(legacy)
    assert obs.shape == (N_OBS, 2) and list(var.index[:3]) == ["MT-CO1", "ACTB", "ACTB"]

    info = h5_matrix_info(legacy)
    assert info["layout"] == "h5ad_legacy" and info["format"] == "csr"
    assert (info["n_obs"], info["n_vars"], info["nnz"]) == (N_OBS, N_VARS, _matrix().nnz)
    assert info["indices_dtype"] == "int64"


def test_read_h5_rows_rejects_bad_input(legacy, tmp_path):
    with pytest.raises(ValueError, match="length"):
        sct.read_h5_rows(legacy, np.ones(3, dtype=bool))
    with pytest.raises(ValueError, match="duplicates"):
        sct.read_h5_rows(legacy, [1, 1])
    with pytest.raises(IndexError):
        sct.read_h5_rows(legacy, [N_OBS])
    with pytest.raises(TypeError, match="boolean mask or integer"):
        sct.read_h5_rows(legacy, ["cell1"])
    csc = write_legacy_h5(tmp_path / "csc.h5", fmt="csc")
    with pytest.raises(TypeError, match="'csc'"):
        sct.read_h5_rows(csc, [0])


def test_row_chunks_cover_all_rows_within_the_entry_budget():
    indptr = np.concatenate([[0], np.cumsum([3, 0, 5, 20, 1, 1, 0, 4])])
    chunks = list(_row_chunks(indptr, max_nnz=6))

    assert chunks[0][0] == 0 and chunks[-1][1] == len(indptr) - 1
    assert all(a[1] == b[0] for a, b in zip(chunks[:-1], chunks[1:], strict=True))
    for start, stop in chunks:
        assert indptr[stop] - indptr[start] <= 6 or stop - start == 1  # one oversized row
