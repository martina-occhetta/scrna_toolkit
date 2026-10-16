"""Reading single-cell data into AnnData and fetching public datasets."""

from scrna_toolkit.io.census import (
    download_census_h5ad,
    list_census_datasets,
    load_census,
    resolve_census_version,
)
from scrna_toolkit.io.download import DownloadError
from scrna_toolkit.io.geo import (
    GEOError,
    download_geo_processed_files,
    get_geo_metadata,
    list_geo_samples,
    list_geo_supplementary_files,
)
from scrna_toolkit.io.h5 import h5_layout, h5_matrix_info, read_h5_axes, read_h5_rows
from scrna_toolkit.io.loaders import FormatDetectionError, detect_format, load_samples, load_scrna

__all__ = [
    "DownloadError",
    "FormatDetectionError",
    "GEOError",
    "detect_format",
    "download_census_h5ad",
    "download_geo_processed_files",
    "get_geo_metadata",
    "h5_layout",
    "h5_matrix_info",
    "list_census_datasets",
    "list_geo_samples",
    "list_geo_supplementary_files",
    "load_census",
    "load_samples",
    "load_scrna",
    "read_h5_axes",
    "read_h5_rows",
    "resolve_census_version",
]
