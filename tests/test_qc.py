import logging

import numpy as np
import pytest
import scipy.sparse as sp

import scrna_toolkit as sct


def test_calculate_qc_metrics_values(adata):
    sct.calculate_qc_metrics(adata)
    X = adata.X.toarray()

    np.testing.assert_allclose(adata.obs["total_counts"], X.sum(1))
    np.testing.assert_array_equal(adata.obs["n_genes_by_counts"], (X > 0).sum(1))
    mt = X[:, :5].sum(1)
    np.testing.assert_allclose(adata.obs["total_counts_mt"], mt)
    np.testing.assert_allclose(adata.obs["pct_counts_mt"], 100 * mt / X.sum(1), rtol=1e-5)
    assert adata.var["mt"].sum() == 5
    assert adata.n_obs == 120  # no filtering during metric calculation


def test_custom_gene_sets_and_symbol_column(adata):
    adata.var["symbol"] = [f"mt-{n}" if n.startswith("MT-") else n for n in adata.var_names]
    adata.var_names = [f"ENSG{i}" for i in range(adata.n_vars)]

    sct.calculate_qc_metrics(
        adata,
        mitochondrial_prefix="mt-",
        gene_symbol_col="symbol",
        ribosomal_prefixes=["GENE1"],
        hemoglobin_regex=r"^GENE2\d$",
    )

    assert adata.var["mt"].sum() == 5
    assert {"pct_counts_ribo", "pct_counts_hb"} <= set(adata.obs.columns)
    assert adata.var["hb"].sum() == 10

    sct.calculate_qc_metrics(adata, mitochondrial_genes=["mt-MT-G0"], gene_symbol_col="symbol")
    assert adata.var["mt"].sum() == 1


def test_no_mito_genes_warns(adata, caplog):
    with caplog.at_level(logging.WARNING):
        sct.calculate_qc_metrics(adata, mitochondrial_prefix="XYZ-")
    assert "No mt genes" in caplog.text


def test_non_count_matrix_is_rejected_unless_allowed(adata, caplog):
    adata.X = adata.X * 0.5
    with pytest.raises(ValueError, match="non-integer"):
        sct.calculate_qc_metrics(adata)
    with caplog.at_level(logging.WARNING):
        sct.calculate_qc_metrics(adata, allow_non_integer=True)
    assert "accepted as counts" in caplog.text


def test_filter_cells_matches_manual_mask_and_keeps_input(adata):
    sct.calculate_qc_metrics(adata)
    obs = adata.obs
    min_genes = int(obs["n_genes_by_counts"].median())
    expected = (obs["n_genes_by_counts"] >= min_genes) & (obs["pct_counts_mt"] <= 3)

    out, summary = sct.filter_cells(adata, min_genes=min_genes, max_pct_mt=3, return_summary=True)

    assert list(out.obs_names) == list(obs.index[expected])
    assert adata.n_obs == 120  # input untouched
    assert sp.issparse(out.X)
    assert summary.loc["total_removed", "n_failing"] == (~expected).sum()
    assert summary.loc["min_genes", "n_failing"] == (obs["n_genes_by_counts"] < min_genes).sum()
    params = sct.get_history(out)[-1]["params"]
    assert params["min_genes"] == min_genes and params["n_after"] == expected.sum()


def test_filter_cells_inplace(adata):
    sct.calculate_qc_metrics(adata)
    out = sct.filter_cells(
        adata, max_counts=float(adata.obs["total_counts"].median()), inplace=True
    )
    assert out is adata
    assert adata.n_obs < 120


def test_filter_cells_without_thresholds_is_noop(adata):
    sct.calculate_qc_metrics(adata)
    assert sct.filter_cells(adata).n_obs == adata.n_obs


def test_filter_cells_requires_metrics(adata):
    with pytest.raises(KeyError, match="calculate_qc_metrics"):
        sct.filter_cells(adata, min_genes=10)


def test_filter_genes_min_cells_and_counts(adata):
    X = adata.X.tolil()
    X[:, 10] = 0  # gene never detected
    adata.X = X.tocsr()
    n_cells = np.asarray((adata.X > 0).sum(0)).ravel()

    out, summary = sct.filter_genes(adata, min_cells=50, return_summary=True)

    assert out.n_vars == (n_cells >= 50).sum()
    assert "GENE5" not in out.var_names  # column 10
    assert adata.n_vars == 200
    assert sct.filter_genes(adata, min_counts=1).n_vars == 199
    assert summary.loc["total_removed", "n_failing"] == (n_cells < 50).sum()
