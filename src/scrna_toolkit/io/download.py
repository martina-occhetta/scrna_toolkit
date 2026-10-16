"""HTTP helpers with retries for transient failures and atomic file writes."""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any

import requests

logger = logging.getLogger(__name__)


class DownloadError(RuntimeError):
    """Raised when a resource cannot be fetched."""


def _get(url: str, *, retries: int, timeout: float, **kwargs: Any) -> requests.Response:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            response = requests.get(url, timeout=timeout, **kwargs)
        except requests.RequestException as err:  # connection errors, timeouts
            last_error = err
        else:
            if response.status_code < 400:
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
    return _get(url, params=params, timeout=timeout, retries=retries).text


def download_file(
    url: str,
    dest: str | Path,
    *,
    overwrite: bool = False,
    timeout: float = 60,
    retries: int = 3,
    chunk_size: int = 1 << 20,
) -> Path:
    """Stream ``url`` to ``dest`` via a ``.part`` file; skip if ``dest`` exists."""
    dest = Path(dest)
    if dest.exists() and not overwrite:
        logger.info("%s exists; skipping download", dest)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    response = _get(url, stream=True, timeout=timeout, retries=retries)
    expected = response.headers.get("Content-Length")
    with part.open("wb") as fh:
        for chunk in response.iter_content(chunk_size=chunk_size):
            fh.write(chunk)
    if expected is not None and part.stat().st_size != int(expected):
        part.unlink()
        raise DownloadError(f"Incomplete download of {url}: expected {expected} bytes")
    part.replace(dest)
    logger.info("Downloaded %s (%.1f MB)", dest, dest.stat().st_size / 1e6)
    return dest
