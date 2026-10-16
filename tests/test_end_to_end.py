"""The whole documented workflow on synthetic 10x files: load → ... → report → h5ad."""

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

matplotlib.use("Agg")

from conftest import make_adata  # noqa: E402

import scrna_toolkit as sct  # noqa: E402


def _write_10x(adata, directory):
    directory.mkdir(parents=True)
    scipy.io.mmwrite(directory / "matrix.mtx", sp.coo_matrix(adata.X.T))
    pd.Series(adata.obs_names).to_csv(directory / "barcodes.tsv", index=False, header=False)
    features = pd.DataFrame(
        {
            "id": [f"ENSG{i:05d}" for i in range(adata.n_vars)],
            "name": adata.var_names,
            "type": "Gene Expression",
        }
    )
    features.to_csv(directory / "features.tsv", sep="\t", index=False, header=False)


def test_full_workflow(tmp_path):
    # Two samples sharing barcodes, as Cell Ranger would produce.
    raw = {s: make_adata(n_obs=80, seed=i) for i, s in enumerate(["S1", "S2"])}
    for sample, a in raw.items():
        _write_10x(a, tmp_path / sample)
    sample_table = pd.DataFrame({"sample_id": ["S1", "S2"], "patient_id": ["P1", "P2"]})

    adata = sct.load_samples({s: tmp_path / s for s in raw})
    sct.merge_sample_metadata(adata, sample_table)
    sct.harmonise_obs(adata, {"annotation": "cell_type"}, strict=False)  # absent: skipped
    issues = sct.validate_metadata(adata)
    assert not (issues["severity"] == "error").any()

    sct.calculate_qc_metrics(adata)
    adata = sct.filter_cells(adata, min_genes=20, max_pct_mt=50)
    adata = sct.filter_genes(adata, min_cells=3)
    expected_counts = adata.X.copy()

    sct.preprocess(adata, n_hvgs=50, n_pcs=10, n_neighbors=10)
    summary = sct.summarize_dataset(adata)
    report = sct.qc_report(adata, output_dir=tmp_path / "report")

    assert summary["n_samples"] == 2 and summary["n_patients"] == 2
    assert (adata.layers["counts"] != expected_counts).nnz == 0  # counts preserved exactly
    assert "umap" in report.figures and "cells_per_patient" in report.tables
    assert set(adata.var["gene_ids"]).issubset({f"ENSG{i:05d}" for i in range(200)})

    adata.write_h5ad(tmp_path / "processed.h5ad")
    reloaded = ad.read_h5ad(tmp_path / "processed.h5ad")
    steps = [h["step"] for h in sct.get_history(reloaded)]
    assert steps[:4] == [
        "load_samples",
        "merge_sample_metadata",
        "harmonise_obs",
        "calculate_qc_metrics",
    ]
    assert steps[-1] == "run_leiden"
    assert sct.get_config(reloaded, "preprocessing")["n_pcs"] == 10
    np.testing.assert_allclose(reloaded.obsm["X_umap"], adata.obsm["X_umap"])
