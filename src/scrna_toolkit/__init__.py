"""scrna-toolkit: reusable single-cell RNA-seq preprocessing utilities built on AnnData/Scanpy."""

import logging
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("scrna-toolkit")
except PackageNotFoundError:  # running from a source tree without installation
    __version__ = "0.0.0+unknown"

from scrna_toolkit.counts import ensure_raw_counts
from scrna_toolkit.datasets import (
    CensusDataset,
    DatasetAdapter,
    GEODataset,
    get_dataset,
    list_datasets,
    register_dataset,
)
from scrna_toolkit.io import (
    detect_format,
    download_census_h5ad,
    download_file,
    download_geo_processed_files,
    get_geo_metadata,
    list_census_datasets,
    list_geo_samples,
    list_geo_supplementary_files,
    load_census,
    load_samples,
    load_scrna,
    read_h5_axes,
    read_h5_rows,
)
from scrna_toolkit.memory import (
    DenseConversionError,
    dense_guard,
    dense_guard_limit,
    install_dense_guard,
    remove_dense_guard,
)
from scrna_toolkit.metadata import harmonise_obs, merge_sample_metadata, validate_metadata
from scrna_toolkit.plotting import (
    plot_cell_type_composition,
    plot_cells_per_group,
    plot_cells_per_patient,
    plot_cells_per_sample,
    plot_composition_heatmap,
    plot_qc_distributions,
    plot_qc_scatter,
    plot_umap,
)
from scrna_toolkit.preprocessing import (
    log1p,
    normalize_total,
    preprocess,
    run_leiden,
    run_neighbors,
    run_pca,
    run_umap,
    scale,
    select_hvgs,
    store_counts,
)
from scrna_toolkit.provenance import get_config, get_history, get_representation, record_step
from scrna_toolkit.qc import (
    QCReport,
    calculate_qc_metrics,
    filter_cells,
    filter_genes,
    qc_report,
    qc_summary,
)
from scrna_toolkit.sampling import downsample_cells, stratified_sample
from scrna_toolkit.summary import summarize_by, summarize_dataset

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = [
    "CensusDataset",
    "DatasetAdapter",
    "DenseConversionError",
    "GEODataset",
    "QCReport",
    "__version__",
    "calculate_qc_metrics",
    "dense_guard",
    "dense_guard_limit",
    "detect_format",
    "downsample_cells",
    "download_census_h5ad",
    "download_file",
    "download_geo_processed_files",
    "ensure_raw_counts",
    "filter_cells",
    "filter_genes",
    "get_config",
    "get_dataset",
    "get_geo_metadata",
    "get_history",
    "get_representation",
    "harmonise_obs",
    "install_dense_guard",
    "list_census_datasets",
    "list_datasets",
    "list_geo_samples",
    "list_geo_supplementary_files",
    "load_census",
    "load_samples",
    "load_scrna",
    "log1p",
    "merge_sample_metadata",
    "normalize_total",
    "plot_cell_type_composition",
    "plot_cells_per_group",
    "plot_cells_per_patient",
    "plot_cells_per_sample",
    "plot_composition_heatmap",
    "plot_qc_distributions",
    "plot_qc_scatter",
    "plot_umap",
    "preprocess",
    "qc_report",
    "qc_summary",
    "read_h5_axes",
    "read_h5_rows",
    "record_step",
    "register_dataset",
    "remove_dense_guard",
    "run_leiden",
    "run_neighbors",
    "run_pca",
    "run_umap",
    "scale",
    "select_hvgs",
    "store_counts",
    "stratified_sample",
    "summarize_by",
    "summarize_dataset",
    "validate_metadata",
]
