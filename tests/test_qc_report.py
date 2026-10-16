import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

import scrna_toolkit as sct  # noqa: E402

MAPPING = {"orig.ident": "sample_id", "patient": "patient_id", "annotation": "cell_type"}


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def test_qc_summary_quantiles(adata):
    sct.calculate_qc_metrics(adata)
    summary = sct.qc_summary(adata)
    assert list(summary.index) == ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    assert summary.loc["total_counts", "median"] == adata.obs["total_counts"].median()


def test_qc_report_full(adata, tmp_path):
    sct.harmonise_obs(adata, MAPPING)
    sct.calculate_qc_metrics(adata)
    sct.preprocess(adata, n_hvgs=50, n_pcs=10, n_neighbors=10)

    report = sct.qc_report(adata, output_dir=tmp_path)

    assert {
        "dataset_summary",
        "qc_summary",
        "cells_per_sample",
        "cells_per_patient",
        "cells_per_cell_type",
        "metadata_issues",
    } <= set(report.tables)
    assert {
        "qc_distributions",
        "qc_scatter",
        "cells_per_sample",
        "cells_per_patient",
        "cell_type_composition",
        "umap",
    } == set(report.figures)
    assert (tmp_path / "qc_summary.csv").is_file()
    assert (tmp_path / "umap.png").is_file()


def test_qc_report_without_metadata_or_embeddings(adata, tmp_path):
    report = sct.qc_report(adata)
    assert set(report.tables) == {"dataset_summary", "metadata_issues"}
    assert report.figures == {}
    assert not any(tmp_path.iterdir())  # nothing written without output_dir
