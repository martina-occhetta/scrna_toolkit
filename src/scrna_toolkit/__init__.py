"""scrna-toolkit: reusable single-cell RNA-seq preprocessing utilities built on AnnData/Scanpy."""

import logging
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("scrna-toolkit")
except PackageNotFoundError:  # running from a source tree without installation
    __version__ = "0.0.0+unknown"

from scrna_toolkit.io import detect_format, load_samples, load_scrna
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
from scrna_toolkit.qc import calculate_qc_metrics, filter_cells, filter_genes
from scrna_toolkit.summary import summarize_by, summarize_dataset

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = [
    "__version__",
    "calculate_qc_metrics",
    "detect_format",
    "filter_cells",
    "filter_genes",
    "get_config",
    "get_history",
    "get_representation",
    "harmonise_obs",
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
    "record_step",
    "run_leiden",
    "run_neighbors",
    "run_pca",
    "run_umap",
    "scale",
    "select_hvgs",
    "store_counts",
    "summarize_by",
    "summarize_dataset",
    "validate_metadata",
]
