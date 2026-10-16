"""GEO tests with the network mocked out (see test_integration.py for live GEO tests)."""

import pytest
import requests

from scrna_toolkit.io import download, geo
from scrna_toolkit.io.download import DownloadError, download_file

SERIES = """^SERIES = GSE1
!Series_title = A test series
!Series_summary = Line one.
!Series_sample_organism = Homo sapiens
!Series_sample_id = GSM10
!Series_sample_id = GSM11
!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1/suppl/GSE1_counts.tsv.gz
!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1/suppl/GSE1_cellInfo.tsv.gz
!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1/suppl/GSE1_RAW.tar
"""

SAMPLES = """^SAMPLE = GSM10
!Sample_title = P1_pre
!Sample_source_name_ch1 = tumour
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = individual: patient P1
!Sample_characteristics_ch1 = treatment phase: Treatment-naive
!Sample_supplementary_file_1 = ftp://x/GSM10_matrix.mtx.gz
!Sample_supplementary_file_2 = ftp://x/GSM10_R1.fastq.gz
^SAMPLE = GSM11
!Sample_title = P1_post
!Sample_characteristics_ch1 = individual: patient P1
!Sample_characteristics_ch1 = treatment phase: post-NACT
!Sample_supplementary_file_1 = NONE
"""


@pytest.fixture
def fake_geo(monkeypatch, tmp_path):
    """Serve canned SOFT text and record every request; isolate the cache."""
    monkeypatch.setenv("SCRNA_TOOLKIT_CACHE", str(tmp_path / "cache"))
    calls = []

    def fake_fetch_text(url, params=None, **_):
        calls.append(params)
        if params["acc"] == "GSE404":
            return "<!DOCTYPE HTML><html>not found</html>"
        return SERIES if params["targ"] == "self" else SAMPLES

    monkeypatch.setattr(geo, "fetch_text", fake_fetch_text)
    return calls


def test_get_geo_metadata_parses_and_caches(fake_geo):
    meta = geo.get_geo_metadata("gse1")

    assert meta["accession"] == "GSE1"
    assert meta["title"] == "A test series"
    assert meta["sample_ids"] == ["GSM10", "GSM11"]
    assert meta["organisms"] == ["Homo sapiens"]
    assert meta["supplementary_files"][0].startswith("https://ftp.ncbi.nlm.nih.gov/")

    geo.get_geo_metadata("GSE1")
    assert len(fake_geo) == 1  # second call served from cache
    geo.get_geo_metadata("GSE1", refresh=True)
    assert len(fake_geo) == 2


def test_invalid_and_unknown_accessions(fake_geo):
    with pytest.raises(ValueError, match="GSE"):
        geo.get_geo_metadata("GSM10")
    with pytest.raises(geo.GEOError):
        geo.get_geo_metadata("GSE404")


def test_list_geo_samples_parses_characteristics(fake_geo):
    samples = geo.list_geo_samples("GSE1")

    assert list(samples["sample_accession"]) == ["GSM10", "GSM11"]
    assert list(samples["individual"]) == ["patient P1", "patient P1"]
    assert list(samples["treatment_phase"]) == ["Treatment-naive", "post-NACT"]


def test_list_supplementary_files_classifies_kinds(fake_geo):
    files = geo.list_geo_supplementary_files("GSE1", include_samples=True)

    kinds = dict(zip(files["filename"], files["kind"], strict=True))
    assert kinds == {
        "GSE1_counts.tsv.gz": "processed",
        "GSE1_cellInfo.tsv.gz": "processed",
        "GSE1_RAW.tar": "raw_archive",
        "GSM10_matrix.mtx.gz": "processed",
        "GSM10_R1.fastq.gz": "sequence",
    }
    assert set(files["level"]) == {"series", "sample"}


def test_download_defaults_to_processed_series_files(fake_geo, monkeypatch, tmp_path):
    downloaded = []
    monkeypatch.setattr(
        geo, "download_file", lambda url, dest, **_: downloaded.append(dest) or dest
    )

    paths = geo.download_geo_processed_files("GSE1", tmp_path / "raw")
    assert [p.name for p in paths] == ["GSE1_counts.tsv.gz", "GSE1_cellInfo.tsv.gz"]
    assert paths[0].parent == tmp_path / "raw" / "GSE1"

    downloaded.clear()
    geo.download_geo_processed_files("GSE1", tmp_path, patterns=["*cellInfo*"])
    assert [p.name for p in downloaded] == ["GSE1_cellInfo.tsv.gz"]

    downloaded.clear()
    geo.download_geo_processed_files("GSE1", tmp_path, kinds=["raw_archive"])
    assert [p.name for p in downloaded] == ["GSE1_RAW.tar"]


class FakeResponse:
    def __init__(self, status, body=b"", headers=None, encoding=None):
        self.status_code, self.body, self.headers = status, body, headers or {}
        self.content, self.encoding = body, encoding

    def iter_content(self, chunk_size):
        yield self.body


def test_download_file_retries_server_errors_and_writes_atomically(monkeypatch, tmp_path):
    responses = [FakeResponse(503), FakeResponse(200, b"data", {"Content-Length": "4"})]
    monkeypatch.setattr(download.requests, "get", lambda *a, **k: responses.pop(0))
    monkeypatch.setattr(download.time, "sleep", lambda s: None)

    path = download_file("https://example.org/f.tsv", tmp_path / "f.tsv")

    assert path.read_bytes() == b"data"
    assert not (tmp_path / "f.tsv.part").exists()
    # existing files are skipped without a request
    monkeypatch.setattr(download.requests, "get", lambda *a, **k: pytest.fail("no request"))
    download_file("https://example.org/f.tsv", tmp_path / "f.tsv")


def test_download_file_errors(monkeypatch, tmp_path):
    monkeypatch.setattr(download.time, "sleep", lambda s: None)
    monkeypatch.setattr(download.requests, "get", lambda *a, **k: FakeResponse(404))
    with pytest.raises(DownloadError, match="404"):
        download_file("https://example.org/missing", tmp_path / "x")

    def offline(*a, **k):
        raise requests.ConnectionError("offline")

    monkeypatch.setattr(download.requests, "get", offline)
    with pytest.raises(DownloadError, match="after 3 attempts"):
        download_file("https://example.org/x", tmp_path / "x")

    monkeypatch.setattr(
        download.requests, "get", lambda *a, **k: FakeResponse(200, b"ab", {"Content-Length": "9"})
    )
    with pytest.raises(DownloadError, match="Incomplete"):
        download_file("https://example.org/x", tmp_path / "x")
    assert not (tmp_path / "x").exists()


def test_fetch_text_decodes_utf8_without_charset(monkeypatch):
    body = "treatment-naïve".encode()
    monkeypatch.setattr(
        download.requests,
        "get",
        lambda *a, **k: FakeResponse(200, body, {"Content-Type": "text/plain"}, "ISO-8859-1"),
    )
    assert download.fetch_text("https://example.org") == "treatment-naïve"
    monkeypatch.setattr(
        download.requests,
        "get",
        lambda *a, **k: FakeResponse(
            200,
            "naïve".encode("latin-1"),
            {"Content-Type": "text/plain; charset=ISO-8859-1"},
            "ISO-8859-1",
        ),
    )
    assert download.fetch_text("https://example.org") == "naïve"
