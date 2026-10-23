"""HTTP helpers with retries for transient failures and atomic file writes.

:func:`download_file` streams to ``<name>.part`` and resumes from it with an HTTP ``Range``
request after an interruption. Existing files are skipped by default. Optionally they are
validated instead: against a size/SHA-256 manifest (``download_manifest.json``), the server's
reported size, an expected checksum and/or a caller-supplied validator.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import time
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

logger = logging.getLogger(__name__)


MANIFEST_NAME = "download_manifest.json"


class DownloadError(RuntimeError):
    """Raised when a resource cannot be fetched."""


class InvalidFileError(DownloadError):
    """Raised when a downloaded or existing local file fails validation."""


def _get(
    url: str, *, retries: int, timeout: float, accept: tuple[int, ...] = (), **kwargs: Any
) -> requests.Response:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            response = requests.get(url, timeout=timeout, **kwargs)
        except requests.RequestException as err:  # connection errors, timeouts
            last_error = err
        else:
            if response.status_code < 400 or response.status_code in accept:
                return response
            if response.status_code < 500:  # client errors will not fix themselves
                raise DownloadError(f"HTTP {response.status_code} for {url}")
            last_error = DownloadError(f"HTTP {response.status_code}")
        if attempt < retries:
            logger.warning(
                "Request to %s failed (%s); retry %d/%d", url, last_error, attempt, retries - 1
            )
            time.sleep(2**attempt)
    raise DownloadError(f"Could not fetch {url} after {retries} attempts: {last_error}")


def fetch_text(
    url: str, *, params: dict[str, str] | None = None, timeout: float = 60, retries: int = 3
) -> str:
    response = _get(url, params=params, timeout=timeout, retries=retries)
    # NCBI serves text/plain without a charset, so requests would guess ISO-8859-1.
    declared = "charset=" in response.headers.get("Content-Type", "")
    encoding = response.encoding if declared and response.encoding else "utf-8"
    return response.content.decode(encoding, errors="replace")


def sha256sum(path: str | Path, chunk_size: int = 1 << 24) -> str:
    """Hex SHA-256 digest of a file, read in chunks."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: str | Path) -> dict[str, dict[str, Any]]:
    """Read a download manifest (``{}`` if it does not exist)."""
    path = Path(path)
    return json.loads(path.read_text()) if path.is_file() else {}


def _write_manifest(path: Path, manifest: dict[str, dict[str, Any]]) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    tmp.replace(path)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def remote_size(
    url: str, *, timeout: float = 60, headers: dict[str, str] | None = None
) -> int | None:
    """``Content-Length`` reported by a HEAD request, or ``None`` if unavailable."""
    try:
        response = requests.head(url, allow_redirects=True, timeout=timeout, headers=headers)
    except requests.RequestException as err:
        logger.warning("HEAD %s failed: %s", url, err)
        return None
    if response.status_code >= 400:
        logger.warning("HEAD %s returned HTTP %d", url, response.status_code)
        return None
    length = response.headers.get("Content-Length")
    return int(length) if length is not None else None


def _content_range_total(headers: Any) -> int | None:
    """Total size from a ``Content-Range: bytes a-b/total`` (or ``bytes */total``) header."""
    match = re.search(r"/(\d+)\s*$", headers.get("Content-Range", ""))
    return int(match.group(1)) if match else None


def _stream_to_part(
    url: str,
    part: Path,
    *,
    resume: bool,
    headers: dict[str, str],
    timeout: float,
    retries: int,
    chunk_size: int,
) -> int | None:
    """Fill ``part`` with the body of ``url``; return the expected total size if known."""
    expected: int | None = None
    for attempt in range(1, retries + 1):
        offset = part.stat().st_size if resume and part.exists() else 0
        request_headers = dict(headers)
        if offset:
            request_headers["Range"] = f"bytes={offset}-"
            logger.info("Resuming %s at %.1f MB", part.name, offset / 1e6)
        response = _get(
            url,
            stream=True,
            timeout=timeout,
            retries=retries,
            accept=(416,),
            headers=request_headers,
        )
        if response.status_code == 416:  # nothing left to send from ``offset``
            expected = _content_range_total(response.headers)
            if expected == offset:
                return expected
            logger.warning("Server rejected resuming %s; restarting", part.name)
            part.unlink(missing_ok=True)
            continue
        length = response.headers.get("Content-Length")
        if response.status_code == 206 and offset:
            mode = "ab"
            expected = _content_range_total(response.headers)
            if expected is None and length is not None:
                expected = offset + int(length)
        else:  # 200: the server sent the whole file (or ignored the Range header)
            mode, offset = "wb", 0
            expected = int(length) if length is not None else None
        try:
            with part.open(mode) as fh:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    fh.write(chunk)
        except requests.RequestException as err:
            if attempt == retries:
                raise DownloadError(f"Download of {url} interrupted: {err}") from err
            logger.warning("Download of %s interrupted (%s); retry %d", url, err, attempt)
            time.sleep(2**attempt)
            continue
        size = part.stat().st_size
        if expected is None or size == expected:
            return expected
        if size > expected or not resume:
            part.unlink()
            raise DownloadError(f"Incomplete download of {url}: expected {expected} bytes")
        if attempt < retries:
            logger.warning("Short read for %s (%d of %d bytes); resuming", url, size, expected)
    raise DownloadError(
        f"Incomplete download of {url}: expected {expected} bytes after {retries} attempts; "
        f"partial data kept in {part} (call again to resume)"
    )


def download_file(
    url: str,
    dest: str | Path,
    *,
    overwrite: bool = False,
    timeout: float = 60,
    retries: int = 3,
    chunk_size: int = 1 << 20,
    resume: bool = True,
    user_agent: str | None = None,
    validator: Callable[[Path], object] | None = None,
    sha256: str | None = None,
    manifest: bool | str | Path = False,
    verify_sha256: bool = False,
    check_remote_size: bool = False,
    offline: bool = False,
) -> Path:
    """Stream ``url`` to ``dest`` via a ``.part`` file; skip or validate ``dest`` if it exists.

    Parameters
    ----------
    overwrite
        Download again even if ``dest`` exists (an old ``.part`` file is discarded too).
    resume
        Continue an interrupted download from ``<dest>.part`` with an HTTP ``Range`` request.
        A server that ignores ``Range`` sends the whole file, which then replaces the part.
    user_agent
        ``User-Agent`` header for every request (some servers reject the default one).
    validator
        Called with the final path after a download and on an existing file. It should raise
        (any exception) if the file is unusable. A freshly downloaded file that fails is deleted.
    sha256
        Expected hex SHA-256 digest. Checked after a download and on an existing file.
    manifest
        ``True`` keeps ``download_manifest.json`` next to ``dest``; a path uses that file. The
        manifest records the URL, size and SHA-256 of each downloaded file. An existing file
        whose size differs from its record raises :class:`InvalidFileError`.
    verify_sha256
        Re-hash an existing file and compare it with the manifest's SHA-256 (slow for large
        files). Needs ``manifest``.
    check_remote_size
        Compare an existing file's size with the server's ``Content-Length`` (HEAD request).
    offline
        Make no network request: a missing ``dest`` raises :class:`FileNotFoundError` and
        ``check_remote_size`` is skipped. Local checks still run on an existing file.

    Without any of ``validator``/``sha256``/``manifest``/``check_remote_size`` an existing
    file is skipped without a request, as before. Validation failures raise
    :class:`InvalidFileError`; delete the file to fetch it again.
    """
    dest = Path(dest)
    if verify_sha256 and not manifest:
        raise ValueError("verify_sha256=True needs manifest")
    headers = {"User-Agent": user_agent} if user_agent else {}
    manifest_path: Path | None = None
    if manifest:
        manifest_path = dest.parent / MANIFEST_NAME if manifest is True else Path(manifest)
    records = read_manifest(manifest_path) if manifest_path else {}
    key = os.path.relpath(dest, manifest_path.parent) if manifest_path else dest.name
    record = dict(records.get(key, {}))

    if offline and (overwrite or not dest.exists()):
        raise FileNotFoundError(f"{dest} is missing (or overwrite=True) and offline=True")

    if dest.exists() and not overwrite:
        size = dest.stat().st_size
        if record.get("size") not in (None, size):
            raise InvalidFileError(
                f"{dest} is {size} bytes but the manifest recorded {record['size']}; "
                "delete it to download again"
            )
        if check_remote_size and not offline:
            remote = remote_size(url, timeout=timeout, headers=headers)
            if remote is not None and remote != size:
                raise InvalidFileError(
                    f"{dest} is {size} bytes but the server reports {remote}; "
                    "delete it to download again"
                )
        digest = sha256sum(dest) if (sha256 or verify_sha256) else None
        if sha256 and digest != sha256.lower():
            raise InvalidFileError(f"{dest}: SHA-256 {digest} does not match expected {sha256}")
        if verify_sha256:
            if record.get("sha256") not in (None, digest):
                raise InvalidFileError(f"{dest}: SHA-256 differs from the manifest")
            record["sha256"] = digest
        _run_validator(validator, dest, delete=False)
        if manifest_path is None:
            logger.info("%s exists; skipping download", dest)
            return dest
        record.setdefault("url", url)
        record["size"] = size
        record["last_validated_at"] = _now()
        logger.info("%s exists and passed validation; skipping download", dest)
    else:
        dest.parent.mkdir(parents=True, exist_ok=True)
        part = dest.with_name(dest.name + ".part")
        if overwrite or not resume:
            part.unlink(missing_ok=True)
        _stream_to_part(
            url,
            part,
            resume=resume,
            headers=headers,
            timeout=timeout,
            retries=retries,
            chunk_size=chunk_size,
        )
        part.replace(dest)
        digest = sha256sum(dest) if (sha256 or manifest_path) else None
        if sha256 and digest != sha256.lower():
            dest.unlink()
            raise InvalidFileError(f"{url}: SHA-256 {digest} does not match expected {sha256}")
        _run_validator(validator, dest, delete=True)
        logger.info("Downloaded %s (%.1f MB)", dest, dest.stat().st_size / 1e6)
        if manifest_path is None:
            return dest
        now = _now()
        record = {
            "url": url,
            "size": dest.stat().st_size,
            "sha256": digest,
            "downloaded_at": now,
            "last_validated_at": now,
        }

    assert manifest_path is not None
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    records = read_manifest(manifest_path)  # re-read: another call may have added entries
    records[key] = record
    _write_manifest(manifest_path, records)
    return dest


def _run_validator(validator: Callable[[Path], object] | None, path: Path, *, delete: bool) -> None:
    if validator is None:
        return
    try:
        validator(path)
    except InvalidFileError:
        if delete:
            path.unlink(missing_ok=True)
        raise
    except Exception as err:
        if delete:
            path.unlink(missing_ok=True)
        raise InvalidFileError(f"{path} failed validation: {err}") from err
