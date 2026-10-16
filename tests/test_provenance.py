import anndata as ad

import scrna_toolkit as sct
from scrna_toolkit.provenance import UNS_KEY


def test_record_step_appends_history_in_order(adata):
    sct.record_step(adata, "first", {"a": 1, "b": None})
    sct.record_step(adata, "second", {"nested": {"x": [1, 2]}})

    history = sct.get_history(adata)
    assert [h["step"] for h in history] == ["first", "second"]
    assert history[0]["params"] == {"a": 1, "b": None}
    assert adata.uns[UNS_KEY]["package_version"] == sct.__version__


def test_history_survives_h5ad_roundtrip(adata, tmp_path):
    sct.record_step(adata, "filter", {"min_genes": 200, "max_genes": None})
    path = tmp_path / "x.h5ad"
    adata.write_h5ad(path)

    history = sct.get_history(ad.read_h5ad(path))
    assert history[0]["params"] == {"min_genes": 200, "max_genes": None}


def test_get_history_empty_without_records(adata):
    assert sct.get_history(adata) == []
