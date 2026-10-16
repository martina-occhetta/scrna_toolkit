"""Census tests with a fake ``cellxgene_census`` module (no install or network needed)."""

import contextlib
import logging
import sys
import types

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

import scrna_toolkit as sct
from scrna_toolkit.io import census


def _census_like_anndata(n_obs=5):
    """Shaped like cellxgene_census.get_anndata output: positional indexes, soma ids."""
    obs = pd.DataFrame(
        {
            "soma_joinid": np.arange(100, 100 + n_obs),
            "donor_id": ["d1"] * n_obs,
            "cell_type": ["T cell"] * n_obs,
        },
        index=[str(i) for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {"feature_id": ["ENSG1", "ENSG2", "ENSG3"], "feature_name": ["MT-CO1", "ACTB", "ACTB"]},
        index=["0", "1", "2"],
    )
    X = sp.csr_matrix(np.arange(n_obs * 3, dtype=np.float32).reshape(n_obs, 3))
    return AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def fake_census(monkeypatch):
    calls = {"open_soma": [], "get_anndata": [], "download": []}
    module = types.ModuleType("cellxgene_census")

    @contextlib.contextmanager
    def open_soma(census_version=None):
        calls["open_soma"].append(census_version)
        yield object()

    def get_anndata(census, organism, **kwargs):
        calls["get_anndata"].append({"organism": organism, **kwargs})
        return _census_like_anndata(n_obs=calls.get("n_obs", 5))

    def download_source_h5ad(dataset_id, to_path, census_version):
        calls["download"].append((dataset_id, census_version))
        _census_like_anndata().write_h5ad(to_path)

    module.open_soma = open_soma
    module.get_anndata = get_anndata
    module.download_source_h5ad = download_source_h5ad
    module.get_census_version_description = lambda v: {"release_build": "2025-11-08"}
    monkeypatch.setitem(sys.modules, "cellxgene_census", module)
    monkeypatch.setattr(census, "_count_cells", lambda c, org, f: calls.get("count", 5))
    return calls


FILTER = "tissue_general == 'blood' and is_primary_data == True"


def test_missing_optional_dependency_gives_install_hint(monkeypatch):
    monkeypatch.setitem(sys.modules, "cellxgene_census", None)
    with pytest.raises(ImportError, match=r"scrna-toolkit\[census\]"):
        sct.load_census(FILTER)


def test_load_census_follows_package_conventions(fake_census, tmp_path):
    adata = sct.load_census(
        FILTER, var_value_filter="feature_name in ['ACTB']", obs_column_names=["cell_type"]
    )

    assert list(adata.obs_names) == ["100", "101", "102", "103", "104"]  # soma_joinid
    assert adata.var_names.is_unique and adata.var_names[0] == "MT-CO1"
    assert list(adata.var["gene_ids"]) == ["ENSG1", "ENSG2", "ENSG3"]
    assert sp.issparse(adata.X)
    sct.ensure_raw_counts(adata)
    adata.write_h5ad(tmp_path / "census.h5ad")  # duplicate symbols must not break writing
    # the resolved build, not the "stable" alias, is used and recorded
    assert fake_census["open_soma"] == ["2025-11-08"]
    params = sct.get_history(adata)[-1]["params"]
    assert params["census_version"] == "2025-11-08" and params["obs_value_filter"] == FILTER
    kwargs = fake_census["get_anndata"][0]
    assert kwargs["organism"] == "Homo sapiens"
    assert kwargs["var_value_filter"] == "feature_name in ['ACTB']"
    assert kwargs["obs_column_names"] == ["cell_type"]


def test_max_cells_guard_runs_before_fetching(fake_census):
    fake_census["count"] = 2_000_000
    with pytest.raises(ValueError, match="max_cells"):
        sct.load_census(FILTER)
    assert fake_census["get_anndata"] == []
    sct.load_census(FILTER, max_cells=None)  # check disabled
    assert len(fake_census["get_anndata"]) == 1


def test_empty_query_raises(fake_census):
    fake_census["n_obs"] = 0
    with pytest.raises(ValueError, match="no cells"):
        sct.load_census(FILTER, max_cells=None)


def test_warns_when_primary_data_not_filtered(fake_census, caplog):
    with caplog.at_level(logging.WARNING):
        sct.load_census("tissue_general == 'blood'")
    assert "is_primary_data" in caplog.text


def test_download_census_h5ad_skips_existing(fake_census, tmp_path):
    path = sct.download_census_h5ad("abc-123", tmp_path / "census")
    assert path == tmp_path / "census" / "abc-123.h5ad" and path.is_file()
    assert not (tmp_path / "census" / "abc-123.h5ad.part").exists()
    sct.download_census_h5ad("abc-123", tmp_path / "census")
    assert fake_census["download"] == [("abc-123", "2025-11-08")]


def test_census_dataset_caches_query(fake_census, tmp_path):
    dataset = sct.CensusDataset(
        "blood",
        FILTER,
        tmp_path,
        census_version="2025-11-08",
        obs_mapping={"donor_id": "patient_id"},
    )

    adata = dataset.get()
    again = dataset.get()  # served from data_dir/census/blood.h5ad

    assert dataset.path == tmp_path / "census" / "blood.h5ad"
    assert len(fake_census["get_anndata"]) == 1
    assert (adata.obs["patient_id"] == "d1").all() and (adata.obs["dataset"] == "blood").all()
    assert again.shape == adata.shape
    assert [h["step"] for h in sct.get_history(again)][0] == "load_census"


def test_census_dataset_load_before_download(tmp_path):
    with pytest.raises(FileNotFoundError, match="download"):
        sct.CensusDataset("x", FILTER, tmp_path).load()
