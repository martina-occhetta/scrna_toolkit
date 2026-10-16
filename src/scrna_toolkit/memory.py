"""Opt-in guard against accidentally densifying large sparse matrices.

:func:`scrna_toolkit._utils.warn_if_densifying` only logs a warning before this package
densifies a block itself. The guard here is stricter and process-wide: while it is active,
``.toarray()`` (and therefore ``.todense()`` and ``.A``) on any SciPy sparse matrix or array
larger than the limit raises :class:`DenseConversionError`, whoever calls it. Small
densifications inside scanpy/anndata stay allowed.

The guard patches the SciPy classes, so it affects every thread in the process.
"""

from __future__ import annotations

import functools
import logging
import math
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from typing import Any

import scipy.sparse as sp

logger = logging.getLogger(__name__)

DEFAULT_MAX_ELEMENTS = 5e8  # ~2 GB as float32

_SPARSE_CLASSES = (
    "csr_matrix",
    "csc_matrix",
    "coo_matrix",
    "bsr_matrix",
    "lil_matrix",
    "dok_matrix",
    "dia_matrix",
    "csr_array",
    "csc_array",
    "coo_array",
    "bsr_array",
    "lil_array",
    "dok_array",
    "dia_array",
)

_limit: dict[str, float | None] = {"max_elements": None}
_originals: dict[type, Callable[..., Any]] = {}


class DenseConversionError(MemoryError):
    """Raised when the dense guard blocks densifying a sparse matrix above its limit."""


def _toarray_owners() -> set[type]:
    """Every class in the SciPy sparse hierarchy that defines its own ``toarray``."""
    owners: set[type] = set()
    for name in _SPARSE_CLASSES:
        cls = getattr(sp, name, None)
        if cls is None:
            continue
        owners.update(base for base in cls.__mro__ if "toarray" in base.__dict__)
    return owners


def _guarded(original: Callable[..., Any]) -> Callable[..., Any]:
    @functools.wraps(original)
    def toarray(self: Any, *args: Any, **kwargs: Any) -> Any:
        limit = _limit["max_elements"]
        n = math.prod(int(d) for d in self.shape)
        if limit is not None and n > limit:
            shape = " x ".join(f"{int(d):,}" for d in self.shape)
            raise DenseConversionError(
                f"Refusing to densify a {shape} sparse matrix (~{n * 4e-9:.1f} GB as float32; "
                f"dense guard limit {limit:.3g} elements). Work on a subset or stay sparse."
            )
        return original(self, *args, **kwargs)

    return toarray


def install_dense_guard(max_elements: float = DEFAULT_MAX_ELEMENTS) -> None:
    """Make densifying any sparse matrix with more than ``max_elements`` entries raise.

    Calling it again only changes the limit. Undo with :func:`remove_dense_guard`, or use
    :func:`dense_guard` to scope it to a ``with`` block.
    """
    if max_elements <= 0:
        raise ValueError("max_elements must be positive")
    _limit["max_elements"] = float(max_elements)
    for owner in _toarray_owners():
        if owner not in _originals:
            original = owner.__dict__["toarray"]
            _originals[owner] = original
            owner.toarray = _guarded(original)  # type: ignore[attr-defined]
    logger.info("Dense guard active (limit %.3g elements)", max_elements)


def remove_dense_guard() -> None:
    """Restore SciPy's original ``toarray`` methods. Safe to call when no guard is active."""
    for owner, original in _originals.items():
        owner.toarray = original  # type: ignore[attr-defined]
    _originals.clear()
    _limit["max_elements"] = None


def dense_guard_limit() -> float | None:
    """The active limit in elements, or ``None`` when no guard is installed."""
    return _limit["max_elements"]


@contextmanager
def dense_guard(max_elements: float = DEFAULT_MAX_ELEMENTS) -> Iterator[None]:
    """Context manager form of :func:`install_dense_guard`.

    On exit the previous state is restored: an enclosing guard keeps its own limit, and the
    guard is removed if none was active before.
    """
    previous = _limit["max_elements"]
    install_dense_guard(max_elements)
    try:
        yield
    finally:
        if previous is None:
            remove_dense_guard()
        else:
            _limit["max_elements"] = previous
