"""Provenance records stored in ``adata.uns["scrna_toolkit"]``.

Only h5ad-serialisable primitives are stored. Parameters are kept as JSON strings so that
``None`` values and nested structures survive ``write_h5ad`` unchanged.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Any

from anndata import AnnData

UNS_KEY = "scrna_toolkit"


def _to_json(obj: Any) -> str:
    def default(o: Any) -> Any:
        return o.item() if hasattr(o, "item") else str(o)  # numpy scalars, Paths, ...

    return json.dumps(obj, default=default)


def _store(adata: AnnData) -> dict:
    from scrna_toolkit import __version__

    store = adata.uns.setdefault(UNS_KEY, {})
    store["package_version"] = __version__
    return store


def record_step(adata: AnnData, step: str, params: dict[str, Any] | None = None) -> None:
    """Append a transformation record to ``adata.uns["scrna_toolkit"]["history"]``."""
    history = _store(adata).setdefault("history", {})
    history[f"{len(history):03d}_{step}"] = {
        "step": step,
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "params": _to_json(params or {}),
    }


def get_history(adata: AnnData) -> list[dict[str, Any]]:
    """Return recorded transformations in order, with parameters decoded from JSON."""
    history = adata.uns.get(UNS_KEY, {}).get("history", {})
    return [
        {"step": rec["step"], "timestamp": rec["timestamp"], "params": json.loads(rec["params"])}
        for _, rec in sorted(history.items())
    ]


def record_config(adata: AnnData, name: str, config: dict[str, Any]) -> None:
    """Store a configuration (e.g. ``"preprocessing"``) as ``uns[...]["<name>_config"]``."""
    _store(adata)[f"{name}_config"] = _to_json(config)


def get_config(adata: AnnData, name: str) -> dict[str, Any] | None:
    """Return a configuration stored with :func:`record_config`, or ``None``."""
    raw = adata.uns.get(UNS_KEY, {}).get(f"{name}_config")
    return None if raw is None else json.loads(raw)


def set_representation(adata: AnnData, representation: str) -> None:
    """Record what ``adata.X`` currently holds: counts, normalized, log1p or scaled."""
    _store(adata)["X_representation"] = representation


def get_representation(adata: AnnData) -> str | None:
    """Return what ``adata.X`` holds, or ``None`` if unknown (not set by this package)."""
    return adata.uns.get(UNS_KEY, {}).get("X_representation")
