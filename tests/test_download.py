"""download_file: resume, validation and manifest, with requests mocked out."""

import hashlib
import json

import pytest
import requests

from scrna_toolkit.io import download
from scrna_toolkit.io.download import DownloadError, InvalidFileError, download_file

URL = "https://example.org/data.bin"
BODY = bytes(range(256)) * 40  # 10,240 bytes


class FakeResponse:
    def __init__(self, status, body=b"", headers=None, fail_after=None):
        self.status_code, self.body, self.headers = status, body, headers or {}
        self.fail_after = fail_after

    def iter_content(self, chunk_size):
        if self.fail_after is None:
            yield self.body
            return
        yield self.body[: self.fail_after]
        raise requests.exceptions.ChunkedEncodingError("connection reset")


class FakeServer:
    """Serves BODY over GET/HEAD, honouring ``Range`` unless ``ranges=False``."""

    def __init__(self, body=BODY, ranges=True, fail_after=()):
        self.body, self.ranges = body, ranges
        self.fail_after = list(fail_after)  # truncate the next responses at these offsets
        self.requests: list[dict] = []

    def get(self, url, headers=None, **kwargs):
        headers = headers or {}
        self.requests.append(headers)
        fail = self.fail_after.pop(0) if self.fail_after else None
        total = len(self.body)
        if self.ranges and "Range" in headers:
            start = int(headers["Range"].removeprefix("bytes=").rstrip("-"))
            if start >= total:
                return FakeResponse(416, headers={"Content-Range": f"bytes */{total}"})
            chunk = self.body[start:]
            return FakeResponse(
                206,
                chunk,
                {
                    "Content-Length": str(len(chunk)),
                    "Content-Range": f"bytes {start}-{total - 1}/{total}",
                },
                fail,
            )
        return FakeResponse(200, self.body, {"Content-Length": str(total)}, fail)

    def head(self, url, headers=None, **kwargs):
        self.requests.append(headers or {})
        return FakeResponse(200, headers={"Content-Length": str(len(self.body))})


@pytest.fixture
def server(monkeypatch):
    fake = FakeServer()
    monkeypatch.setattr(download.requests, "get", fake.get)
    monkeypatch.setattr(download.requests, "head", fake.head)
    monkeypatch.setattr(download.time, "sleep", lambda s: None)
    return fake


def test_resumes_from_part_file_with_range_request(server, tmp_path):
    dest = tmp_path / "data.bin"
    (tmp_path / "data.bin.part").write_bytes(BODY[:1000])

    download_file(URL, dest, user_agent="test-agent/1.0")

    assert dest.read_bytes() == BODY
    assert server.requests == [{"User-Agent": "test-agent/1.0", "Range": "bytes=1000-"}]
    assert not (tmp_path / "data.bin.part").exists()


def test_interrupted_stream_is_resumed_within_the_same_call(server, tmp_path):
    server.fail_after = [4000, 3000]  # two connection resets, then a clean response

    download_file(URL, tmp_path / "data.bin")

    assert (tmp_path / "data.bin").read_bytes() == BODY
    assert [r.get("Range") for r in server.requests] == [None, "bytes=4000-", "bytes=7000-"]


def test_interruption_on_last_attempt_keeps_part_for_next_call(server, tmp_path):
    server.fail_after = [3000, 3000, 3000]  # retries=3: every attempt is cut short
    with pytest.raises(DownloadError, match="interrupted"):
        download_file(URL, tmp_path / "data.bin")
    assert (tmp_path / "data.bin.part").stat().st_size == 9000
    assert not (tmp_path / "data.bin").exists()

    download_file(URL, tmp_path / "data.bin")
    assert server.requests[-1]["Range"] == "bytes=9000-"
    assert (tmp_path / "data.bin").read_bytes() == BODY


def test_server_ignoring_range_restarts_from_scratch(server, tmp_path):
    server.ranges = False
    (tmp_path / "data.bin.part").write_bytes(b"stale bytes")

    download_file(URL, tmp_path / "data.bin")

    assert (tmp_path / "data.bin").read_bytes() == BODY


def test_complete_part_file_is_accepted_on_416(server, tmp_path):
    (tmp_path / "data.bin.part").write_bytes(BODY)
    download_file(URL, tmp_path / "data.bin")
    assert (tmp_path / "data.bin").read_bytes() == BODY


def test_resume_false_and_overwrite_discard_part_file(server, tmp_path):
    (tmp_path / "data.bin.part").write_bytes(b"x" * 100)
    download_file(URL, tmp_path / "data.bin", resume=False)
    assert server.requests[-1] == {}
    assert (tmp_path / "data.bin").read_bytes() == BODY

    (tmp_path / "data.bin.part").write_bytes(b"x" * 100)
    download_file(URL, tmp_path / "data.bin", overwrite=True)
    assert "Range" not in server.requests[-1]
    assert (tmp_path / "data.bin").read_bytes() == BODY


def test_validator_failure_deletes_fresh_download(server, tmp_path):
    def reject(path):
        raise ValueError("bad header")

    with pytest.raises(InvalidFileError, match="bad header"):
        download_file(URL, tmp_path / "data.bin", validator=reject)
    assert not (tmp_path / "data.bin").exists()


def test_validator_runs_on_existing_file_without_deleting_it(server, tmp_path):
    dest = tmp_path / "data.bin"
    dest.write_bytes(b"truncated")
    seen = []

    with pytest.raises(InvalidFileError, match="too small"):
        download_file(URL, dest, validator=lambda p: seen.append(p) or _min_size(p, 100))
    assert dest.exists() and seen == [dest]
    assert server.requests == []  # validated locally, nothing fetched


def _min_size(path, n):
    if path.stat().st_size < n:
        raise ValueError("too small")


def test_expected_sha256(server, tmp_path):
    good = hashlib.sha256(BODY).hexdigest()
    download_file(URL, tmp_path / "a.bin", sha256=good)
    download_file(URL, tmp_path / "a.bin", sha256=good.upper())  # existing file re-checked

    with pytest.raises(InvalidFileError, match="SHA-256"):
        download_file(URL, tmp_path / "b.bin", sha256="0" * 64)
    assert not (tmp_path / "b.bin").exists()


def test_manifest_records_and_validates_existing_files(server, tmp_path):
    dest = tmp_path / "data.bin"
    download_file(URL, dest, manifest=True)

    record = json.loads((tmp_path / "download_manifest.json").read_text())["data.bin"]
    assert record["url"] == URL and record["size"] == len(BODY)
    assert record["sha256"] == hashlib.sha256(BODY).hexdigest()
    n_requests = len(server.requests)

    download_file(URL, dest, manifest=True, verify_sha256=True)  # unchanged: passes offline
    assert len(server.requests) == n_requests

    dest.write_bytes(BODY[:-1] + b"\x00")  # same size, different content
    with pytest.raises(InvalidFileError, match="SHA-256 differs"):
        download_file(URL, dest, manifest=True, verify_sha256=True)
    download_file(URL, dest, manifest=True)  # size-only check cannot see it

    dest.write_bytes(BODY[:100])
    with pytest.raises(InvalidFileError, match="manifest recorded"):
        download_file(URL, dest, manifest=True)


def test_manifest_adopts_existing_file_and_supports_custom_path(server, tmp_path):
    dest = tmp_path / "raw" / "data.bin"
    dest.parent.mkdir()
    dest.write_bytes(BODY)
    manifest = tmp_path / "manifests" / "files.json"

    download_file(URL, dest, manifest=manifest, verify_sha256=True)

    record = json.loads(manifest.read_text())["../raw/data.bin"]
    assert record["size"] == len(BODY) and "last_validated_at" in record
    assert server.requests == []


def test_check_remote_size(server, tmp_path):
    dest = tmp_path / "data.bin"
    dest.write_bytes(BODY[:10])
    with pytest.raises(InvalidFileError, match="server reports"):
        download_file(URL, dest, check_remote_size=True, user_agent="ua")
    assert server.requests == [{"User-Agent": "ua"}]


def test_verify_sha256_needs_manifest(tmp_path):
    with pytest.raises(ValueError, match="manifest"):
        download_file(URL, tmp_path / "x", verify_sha256=True)


def test_offline_raises_for_missing_file_without_a_request(server, tmp_path):
    dest = tmp_path / "data.bin"
    with pytest.raises(FileNotFoundError, match="offline"):
        download_file(URL, dest, offline=True)
    assert server.requests == [] and not dest.exists()


def test_offline_validates_existing_file_without_a_request(server, tmp_path):
    dest = tmp_path / "data.bin"
    download_file(URL, dest, manifest=True)
    server.requests.clear()
    seen = []

    download_file(
        URL, dest, manifest=True, check_remote_size=True, validator=seen.append, offline=True
    )

    assert seen == [dest] and server.requests == []
    dest.write_bytes(BODY[:10])
    with pytest.raises(InvalidFileError, match="manifest recorded"):
        download_file(URL, dest, manifest=True, offline=True)


def test_offline_refuses_overwrite(server, tmp_path):
    dest = tmp_path / "data.bin"
    dest.write_bytes(BODY)
    with pytest.raises(FileNotFoundError, match="offline"):
        download_file(URL, dest, overwrite=True, offline=True)
    assert dest.read_bytes() == BODY and server.requests == []
