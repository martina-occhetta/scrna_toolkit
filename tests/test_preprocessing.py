import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp

import scrna_toolkit as sct
from scrna_toolkit.config import PreprocessingConfig

SMALL = {"n_hvgs": 50, "n_pcs": 10, "n_neighbors": 10}


def test_store_counts_preserves_raw_counts_through_normalization(adata):
    raw = adata.X.toarray().copy()

    sct.store_counts(adata)
    sct.normalize_total(adata, target_sum=1e4)
    sct.log1p(adata)

    assert sp.issparse(adata.layers["counts"]) and sp.issparse(adata.X)
    np.testing.assert_array_equal(adata.layers["counts"].toarray(), raw)
    expected = np.log1p(raw / raw.sum(1, keepdims=True) * 1e4)
    np.testing.assert_allclose(adata.X.toarray(), expected, rtol=1e-4)
    assert sct.get_representation(adata) == "log1p"


def test_store_counts_does_not_overwrite_existing_layer(adata):
    adata.layers["counts"] = adata.X * 2
    sct.store_counts(adata)
    np.testing.assert_array_equal(adata.layers["counts"].toarray(), adata.X.toarray() * 2)


def test_double_normalization_is_refused(adata):
    sct.normalize_total(adata)
    sct.log1p(adata)
    with pytest.raises(ValueError, match="log1p"):
        sct.normalize_total(adata)
    with pytest.raises(ValueError):
        sct.log1p(adata)
    with pytest.raises(ValueError, match="not counts"):
        sct.store_counts(adata)


def _log_normalized(adata):
    sct.store_counts(adata)
    sct.normalize_total(adata)
    sct.log1p(adata)
    return adata


def test_select_hvgs_flags_genes_without_subsetting(adata):
    sct.select_hvgs(_log_normalized(adata), n_top_genes=50)
    assert adata.var["highly_variable"].sum() == 50
    assert adata.n_vars == 200
    with pytest.raises(ValueError, match="exceeds"):
        sct.select_hvgs(adata, n_top_genes=500)


def test_run_pca_on_scaled_hvg_copy_leaves_X_untouched(adata):
    sct.select_hvgs(_log_normalized(adata), n_top_genes=50)
    before = adata.X.copy()

    sct.run_pca(adata, n_comps=10)

    assert adata.obsm["X_pca"].shape == (120, 10)
    assert adata.varm["PCs"].shape == (200, 10)
    assert not adata.varm["PCs"][~adata.var["highly_variable"].to_numpy()].any()
    assert sp.issparse(adata.X)
    assert (adata.X != before).nnz == 0


def test_run_pca_without_scaling_and_validation(adata):
    _log_normalized(adata)
    with pytest.raises(KeyError, match="select_hvgs"):
        sct.run_pca(adata, n_comps=10)
    sct.run_pca(adata, n_comps=10, use_hvgs=False, scale=False)
    assert adata.obsm["X_pca"].shape == (120, 10)
    with pytest.raises(ValueError, match="n_comps"):
        sct.run_pca(adata, n_comps=500, use_hvgs=False)


def test_graph_steps_require_prerequisites(adata):
    with pytest.raises(KeyError, match="run_pca"):
        sct.run_neighbors(adata)
    with pytest.raises(KeyError, match="run_neighbors"):
        sct.run_umap(adata)


def test_preprocess_pipeline_outputs_and_provenance(adata, tmp_path):
    out = sct.preprocess(adata, **SMALL)

    assert out is adata
    assert {"X_pca", "X_umap"} <= set(adata.obsm)
    assert "leiden" in adata.obs and "counts" in adata.layers
    config = sct.get_config(adata, "preprocessing")
    assert config["n_hvgs"] == 50 and config["target_sum"] == 1e4
    steps = [h["step"] for h in sct.get_history(adata)]
    assert steps == [
        "store_counts",
        "normalize_total",
        "log1p",
        "select_hvgs",
        "run_pca",
        "run_neighbors",
        "run_umap",
        "run_leiden",
    ]
    adata.write_h5ad(tmp_path / "processed.h5ad")  # everything stored must be serialisable
    assert sct.get_config(ad.read_h5ad(tmp_path / "processed.h5ad"), "preprocessing") == config


def test_preprocess_is_deterministic(adata):
    a = sct.preprocess(adata, copy=True, **SMALL)
    b = sct.preprocess(adata, copy=True, **SMALL)
    np.testing.assert_allclose(a.obsm["X_umap"], b.obsm["X_umap"])
    assert (a.obs["leiden"] == b.obs["leiden"]).all()
    assert sct.get_representation(adata) is None  # copy=True left the input untouched


def test_preprocess_accepts_config_and_optional_steps(adata):
    config = PreprocessingConfig(n_hvgs=50, n_pcs=10, run_umap=False, run_leiden=False)
    sct.preprocess(adata, config)
    assert "X_umap" not in adata.obsm and "leiden" not in adata.obs
    with pytest.raises(TypeError):
        sct.preprocess(adata, n_hvg=10)
