"""Live GEO tests. Excluded by default; run with ``pytest -m integration``."""

import pytest

import scrna_toolkit as sct

pytestmark = pytest.mark.integration

ACCESSION = "GSE165897"  # small, stable public series with processed matrices


@pytest.fixture(autouse=True)
def _isolated_cache(monkeypatch, tmp_path):
    monkeypatch.setenv("SCRNA_TOOLKIT_CACHE", str(tmp_path))


def test_live_geo_metadata():
    meta = sct.get_geo_metadata(ACCESSION)
    assert meta["accession"] == ACCESSION
    assert meta["title"]
    assert "Homo sapiens" in meta["organisms"]
    assert len(meta["sample_ids"]) > 0
    assert "naïve" in meta["summary"]  # UTF-8 decoding


def test_live_geo_samples_and_files():
    samples = sct.list_geo_samples(ACCESSION)
    assert samples["sample_accession"].str.startswith("GSM").all()
    files = sct.list_geo_supplementary_files(ACCESSION)
    assert (files["kind"] == "processed").any()
    assert files["url"].str.startswith("https://").all()


def test_live_unknown_accession():
    with pytest.raises(sct.io.GEOError):
        sct.get_geo_metadata("GSE999999999")


# Smallest dataset in the 2025-11-08 Census build (146 cells, all is_primary_data == False,
# i.e. duplicates of cells in another dataset); the build date is pinned so it stays valid.
CENSUS_VERSION = "2025-11-08"
CENSUS_DATASET = "0895c838-e550-48a3-a777-dbcd35d30272"


def test_live_census_query():
    pytest.importorskip("cellxgene_census")
    adata = sct.load_census(
        f"dataset_id == '{CENSUS_DATASET}'",
        census_version=CENSUS_VERSION,
        var_value_filter="feature_name in ['MT-CO1', 'ACTB', 'CD3E']",
        obs_column_names=["soma_joinid", "donor_id", "cell_type", "tissue"],
    )
    assert adata.n_obs == 146 and set(adata.var_names) == {"MT-CO1", "ACTB", "CD3E"}
    assert adata.var["gene_ids"].str.startswith("ENSG").all()
    sct.ensure_raw_counts(adata)
