import pandas as pd
import pytest

import scrna_toolkit as sct
from scrna_toolkit.datasets import base, registry
from scrna_toolkit.io import FormatDetectionError


def _write_counts(directory, name="GSE1_counts.csv"):
    directory.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        [[1, 0, 3], [0, 2, 1]], index=["c1", "c2"], columns=["MT-CO1", "ACTB", "CD3E"]
    )
    df.to_csv(directory / name)


@pytest.fixture
def no_network(monkeypatch):
    calls = []
    monkeypatch.setattr(
        base, "download_geo_processed_files", lambda *a, **k: calls.append((a, k)) or []
    )
    return calls


def test_get_dataset_falls_back_to_generic_geo(tmp_path):
    dataset = sct.get_dataset("gse1", data_dir=tmp_path)
    assert isinstance(dataset, sct.GEODataset)
    assert dataset.directory == tmp_path / "GSE1"
    with pytest.raises(KeyError, match="Unknown dataset"):
        sct.get_dataset("my-dataset")


def test_geo_dataset_get_loads_single_matrix_and_harmonises(tmp_path, no_network):
    _write_counts(tmp_path / "GSE1")
    dataset = sct.GEODataset(
        "GSE1",
        tmp_path,
        patterns=["*counts*"],
        load_kwargs={"orientation": "cells_x_genes"},
        obs_mapping={"dataset": "batch"},
    )

    adata = dataset.get()

    assert no_network[0][1]["patterns"] == ["*counts*"]
    assert adata.shape == (2, 3)
    assert (adata.obs["dataset"] == "GSE1").all()
    assert (adata.obs["batch"] == "GSE1").all()


def test_geo_dataset_requires_choice_between_several_matrices(tmp_path, no_network):
    _write_counts(tmp_path / "GSE1", "GSE1_a.csv")
    _write_counts(tmp_path / "GSE1", "GSE1_b.csv")
    (tmp_path / "GSE1" / "README.txt").write_text("not a matrix")

    dataset = sct.GEODataset("GSE1", tmp_path, load_kwargs={"orientation": "cells_x_genes"})
    with pytest.raises(FormatDetectionError, match="matrix_file"):
        dataset.load()
    dataset.matrix_file = "GSE1_b.csv"
    assert dataset.load().n_obs == 2


def test_geo_dataset_load_before_download(tmp_path):
    with pytest.raises(FileNotFoundError, match="download"):
        sct.GEODataset("GSE1", tmp_path).load()


def test_register_custom_adapter(monkeypatch, adata):
    monkeypatch.setattr(registry, "_REGISTRY", {})

    @sct.register_dataset("toy")
    class ToyDataset(sct.DatasetAdapter):
        def download(self):
            return []

        def load(self):
            return adata

        def harmonise_metadata(self, a):
            return sct.harmonise_obs(a, {"orig.ident": "sample_id"})

    assert sct.list_datasets() == ["toy"]
    assert "sample_id" in sct.get_dataset("toy").get().obs
    with pytest.raises(ValueError, match="already registered"):
        sct.register_dataset("toy")(ToyDataset)
