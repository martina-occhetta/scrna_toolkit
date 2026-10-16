import numpy as np
import pytest

import scrna_toolkit as sct

MAPPING = {"orig.ident": "sample_id", "patient": "patient_id", "annotation": "cell_type"}


def test_summarize_dataset_with_metadata(adata):
    sct.harmonise_obs(adata, MAPPING)

    summary = sct.summarize_dataset(adata)

    assert summary["n_cells"] == 120 and summary["n_genes"] == 200
    assert summary["n_samples"] == 4 and summary["n_patients"] == 2
    assert summary["n_cell_types"] == 3
    assert summary["median_cells_per_sample"] == 30
    assert summary["minimum_cells_per_patient"] == 60
    assert "sample_id" not in summary["missing_metadata"]
    assert "tissue" in summary["missing_metadata"]


def test_summarize_dataset_without_metadata(adata):
    summary = sct.summarize_dataset(adata)
    assert summary["n_cells"] == 120
    assert np.isnan(summary["n_samples"]) and np.isnan(summary["median_cells_per_patient"])
    assert "patient_id" in summary["missing_metadata"]


def test_summarize_by_patient_and_cell_type(adata):
    sct.harmonise_obs(adata, MAPPING)
    sct.calculate_qc_metrics(adata)

    by_patient = sct.summarize_by(adata, "patient_id")
    assert by_patient["n_cells"].to_dict() == {"P0": 60, "P1": 60}
    assert by_patient["n_samples"].tolist() == [2, 2]
    assert "n_patients" not in by_patient
    assert "median_pct_counts_mt" in by_patient

    by_type = sct.summarize_by(adata, "cell_type")
    assert by_type["n_cells"].sum() == 120
    assert by_type["fraction_of_cells"].sum() == pytest.approx(1.0)
    assert by_type["n_cells"].is_monotonic_decreasing


def test_summarize_by_keeps_missing_labels_and_requires_column(adata):
    adata.obs["cell_type"] = adata.obs["annotation"].astype(object)
    adata.obs.iloc[:5, adata.obs.columns.get_loc("cell_type")] = np.nan
    assert sct.summarize_by(adata, "cell_type")["n_cells"].sum() == 120
    with pytest.raises(KeyError):
        sct.summarize_by(adata, "tissue")
