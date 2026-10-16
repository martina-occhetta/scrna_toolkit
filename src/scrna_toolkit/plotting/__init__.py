"""Exploratory plots. Every function returns a matplotlib ``Figure`` (or ``None`` with a
logged warning when the required data is absent) and only saves when ``path`` is given."""

from scrna_toolkit.plotting.composition import (
    plot_cell_type_composition,
    plot_cells_per_group,
    plot_cells_per_patient,
    plot_cells_per_sample,
    plot_composition_heatmap,
)
from scrna_toolkit.plotting.embedding import plot_umap
from scrna_toolkit.plotting.qc import plot_qc_distributions, plot_qc_scatter

__all__ = [
    "plot_cell_type_composition",
    "plot_cells_per_group",
    "plot_cells_per_patient",
    "plot_cells_per_sample",
    "plot_composition_heatmap",
    "plot_qc_distributions",
    "plot_qc_scatter",
    "plot_umap",
]
