"""scrna-toolkit: reusable single-cell RNA-seq preprocessing utilities built on AnnData/Scanpy."""

import logging
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("scrna-toolkit")
except PackageNotFoundError:  # running from a source tree without installation
    __version__ = "0.0.0+unknown"

from scrna_toolkit.provenance import get_history, record_step

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = ["__version__", "get_history", "record_step"]
