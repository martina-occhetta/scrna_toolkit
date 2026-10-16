"""Reading single-cell data into AnnData and fetching public datasets."""

from scrna_toolkit.io.download import DownloadError
from scrna_toolkit.io.geo import (
    GEOError,
    download_geo_processed_files,
    get_geo_metadata,
    list_geo_samples,
    list_geo_supplementary_files,
)
from scrna_toolkit.io.loaders import FormatDetectionError, detect_format, load_samples, load_scrna

__all__ = [
    "DownloadError",
    "FormatDetectionError",
    "GEOError",
    "detect_format",
    "download_geo_processed_files",
    "get_geo_metadata",
    "list_geo_samples",
    "list_geo_supplementary_files",
    "load_samples",
    "load_scrna",
]
