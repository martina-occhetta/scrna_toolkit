import numpy as np
import pytest
import scipy.sparse as sp

import scrna_toolkit as sct


def test_accepts_sparse_and_dense_counts(adata):
    sct.ensure_raw_counts(adata)
    adata.X = adata.X.toarray()
    sct.ensure_raw_counts(adata, chunk_size=7)  # several chunks, dense


def test_rejects_normalized_values_with_example(adata):
    adata.X = adata.X / 3
    with pytest.raises(ValueError, match=r"non-integer values \(e\.g\."):
        sct.ensure_raw_counts(adata)


def test_problem_outside_first_chunk_is_found(adata):
    X = adata.X.toarray()
    X[-1, -1] = 0.5  # last cell only
    adata.X = sp.csr_matrix(X)
    with pytest.raises(ValueError, match="1 non-integer"):
        sct.ensure_raw_counts(adata, chunk_size=10)


def test_negative_values_always_fail(adata):
    adata.X = adata.X.toarray() - 1.0  # e.g. centred/scaled data
    with pytest.raises(ValueError, match="negative"):
        sct.ensure_raw_counts(adata, allow_non_integer=True)


def test_allow_non_integer(adata, caplog):
    adata.X = adata.X * 0.5
    sct.ensure_raw_counts(adata, allow_non_integer=True)
    assert "accepted as counts" in caplog.text


def test_tracked_representation_and_layers(adata):
    sct.store_counts(adata)
    sct.normalize_total(adata)
    with pytest.raises(ValueError, match="'normalized' values"):
        sct.ensure_raw_counts(adata)
    sct.ensure_raw_counts(adata, layer="counts")
    with pytest.raises(KeyError):
        sct.ensure_raw_counts(adata, layer="spliced")


def test_preprocess_config_allows_non_integer_counts(adata):
    adata.X = adata.X * 0.5
    with pytest.raises(ValueError):
        sct.preprocess(adata, copy=True, n_hvgs=50, n_pcs=10, run_umap=False)
    out = sct.preprocess(
        adata, copy=True, n_hvgs=50, n_pcs=10, run_umap=False, allow_non_integer_counts=True
    )
    assert np.allclose(out.layers["counts"].toarray(), adata.X.toarray())
