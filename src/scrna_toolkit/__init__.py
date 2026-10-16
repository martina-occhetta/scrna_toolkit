"""scrna-toolkit: reusable single-cell RNA-seq preprocessing utilities built on AnnData/Scanpy."""

import logging
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("scrna-toolkit")
except PackageNotFoundError:  # running from a source tree without installation
    __version__ = "0.0.0+unknown"

from scrna_toolkit.io import detect_format, load_samples, load_scrna
from scrna_toolkit.provenance import get_history, record_step

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = [
    "__version__",
    "detect_format",
    "get_history",
    "load_samples",
    "load_scrna",
    "record_step",
]
