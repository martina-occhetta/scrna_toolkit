"""One-call QC overview for notebooks: tables + figures, optionally written to a directory."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import cast

import pandas as pd
from anndata import AnnData
from matplotlib.figure import Figure

from scrna_toolkit.metadata.validation import validate_metadata
from scrna_toolkit.plotting import (
    plot_cell_type_composition,
    plot_cells_per_patient,
    plot_cells_per_sample,
    plot_qc_distributions,
    plot_qc_scatter,
    plot_umap,
)
from scrna_toolkit.summary.dataset import DEFAULT_QC_COLS, summarize_by, summarize_dataset

logger = logging.getLogger(__name__)


@dataclass
class QCReport:
    tables: dict[str, pd.DataFrame] = field(default_factory=dict)
    figures: dict[str, Figure] = field(default_factory=dict)

    def save(self, output_dir: str | Path) -> list[Path]:
        """Write tables as CSV and figures as PNG; return the written paths."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        written = []
        for name, table in self.tables.items():
            written.append(output_dir / f"{name}.csv")
            table.to_csv(written[-1])
        for name, fig in self.figures.items():
            written.append(output_dir / f"{name}.png")
            fig.savefig(written[-1], bbox_inches="tight", dpi=150)
        logger.info("QC report written to %s (%d files)", output_dir, len(written))
        return written


def qc_summary(adata: AnnData, metrics: Sequence[str] = DEFAULT_QC_COLS) -> pd.DataFrame:
    """Distribution of each QC metric present (rows) as min/quantiles/max/mean."""
    obs = cast(pd.DataFrame, adata.obs)
    present = [m for m in metrics if m in obs]
    quantiles = obs[present].quantile([0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0]).T
    quantiles.columns = ["min", "q05", "q25", "median", "q75", "q95", "max"]
    quantiles["mean"] = obs[present].mean()
    return quantiles


def qc_report(
    adata: AnnData,
    *,
    sample_col: str = "sample_id",
    patient_col: str = "patient_id",
    cell_type_col: str = "cell_type",
    cluster_col: str = "leiden",
    output_dir: str | Path | None = None,
) -> QCReport:
    """Build summary tables and plots for whatever metadata/QC/embeddings exist.

    Missing pieces are skipped, never fabricated. Returns a :class:`QCReport`; when
    ``output_dir`` is given everything is also written there.
    """
    obs = adata.obs
    report = QCReport()
    report.tables["dataset_summary"] = summarize_dataset(
        adata, sample_col=sample_col, patient_col=patient_col, cell_type_col=cell_type_col
    ).to_frame()
    report.tables["metadata_issues"] = validate_metadata(adata, sample_col=sample_col)
    if any(m in obs for m in DEFAULT_QC_COLS):
        report.tables["qc_summary"] = qc_summary(adata)
    for name, col in (
        ("cells_per_sample", sample_col),
        ("cells_per_patient", patient_col),
        ("cells_per_cell_type", cell_type_col),
    ):
        if col in obs:
            report.tables[name] = summarize_by(
                adata, col, sample_col=sample_col, patient_col=patient_col
            )

    group = sample_col if sample_col in obs else None
    figures = {
        "qc_distributions": plot_qc_distributions(adata, groupby=group),
        "qc_scatter": plot_qc_scatter(adata),
        "cells_per_sample": plot_cells_per_sample(adata, sample_col) if sample_col in obs else None,
        "cells_per_patient": (
            plot_cells_per_patient(adata, patient_col) if patient_col in obs else None
        ),
        "cell_type_composition": (
            plot_cell_type_composition(adata, sample_col, cell_type_col)
            if sample_col in obs and cell_type_col in obs
            else None
        ),
    }
    if "X_umap" in adata.obsm:
        colors = [c for c in (sample_col, patient_col, cell_type_col, cluster_col) if c in obs]
        figures["umap"] = plot_umap(adata, color=colors)
    report.figures = {name: fig for name, fig in figures.items() if fig is not None}

    if output_dir is not None:
        report.save(output_dir)
    return report
