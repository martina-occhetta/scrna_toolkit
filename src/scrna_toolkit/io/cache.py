"""Local on-disk cache for small text resources such as GEO metadata."""

from __future__ import annotations

import os
from collections.abc import Callable
from pathlib import Path


def cache_dir() -> Path:
    """``$SCRNA_TOOLKIT_CACHE`` if set, else ``~/.cache/scrna_toolkit``."""
    return Path(os.environ.get("SCRNA_TOOLKIT_CACHE", Path.home() / ".cache" / "scrna_toolkit"))


def cached_text(key: str, fetch: Callable[[], str], *, refresh: bool = False) -> str:
    """Return the cached text for ``key``, calling ``fetch`` (and caching) on a miss."""
    path = cache_dir() / key
    if path.is_file() and not refresh:
        return path.read_text()
    text = fetch()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    tmp.replace(path)
    return text
