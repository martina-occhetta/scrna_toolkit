"""Name -> adapter registry. Unregistered GEO accessions fall back to :class:`GEODataset`."""

from __future__ import annotations

import re
from collections.abc import Callable
from typing import Any

from scrna_toolkit.datasets.base import DatasetAdapter, GEODataset

_REGISTRY: dict[str, Callable[..., DatasetAdapter]] = {}


def register_dataset(
    name: str,
) -> Callable[[Callable[..., DatasetAdapter]], Callable[..., DatasetAdapter]]:
    """Decorator registering an adapter class (or factory) under ``name``."""

    def decorator(factory: Callable[..., DatasetAdapter]) -> Callable[..., DatasetAdapter]:
        if name in _REGISTRY:
            raise ValueError(f"Dataset {name!r} is already registered")
        _REGISTRY[name] = factory
        return factory

    return decorator


def get_dataset(name: str, **kwargs: Any) -> DatasetAdapter:
    """Return the registered adapter for ``name``, or a generic :class:`GEODataset`."""
    if name in _REGISTRY:
        return _REGISTRY[name](**kwargs)
    if re.fullmatch(r"GSE\d+", name, flags=re.IGNORECASE):
        return GEODataset(name, **kwargs)
    raise KeyError(f"Unknown dataset {name!r}. Registered: {list_datasets()}")


def list_datasets() -> list[str]:
    return sorted(_REGISTRY)
