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
