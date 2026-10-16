"""scrna-toolkit: reusable single-cell RNA-seq preprocessing utilities built on AnnData/Scanpy."""

import logging
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("scrna-toolkit")
except PackageNotFoundError:  # running from a source tree without installation
    __version__ = "0.0.0+unknown"

from scrna_toolkit.io import detect_format, load_samples, load_scrna
from scrna_toolkit.metadata import harmonise_obs, merge_sample_metadata, validate_metadata
from scrna_toolkit.provenance import get_history, record_step
from scrna_toolkit.qc import calculate_qc_metrics, filter_cells, filter_genes

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = [
    "__version__",
    "calculate_qc_metrics",
    "detect_format",
    "filter_cells",
    "filter_genes",
    "get_history",
    "harmonise_obs",
    "load_samples",
    "load_scrna",
    "merge_sample_metadata",
    "record_step",
    "validate_metadata",
]
