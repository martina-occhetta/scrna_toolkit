"""Query and download data from the CZ CELLxGENE Discover Census.

Requires the optional dependency: ``pip install "scrna-toolkit[census]"``.

``"stable"`` and ``"latest"`` are moving targets, so each call resolves the requested version
to a concrete build date once, uses it for every request, and records it in provenance.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd
from anndata import AnnData

from scrna_toolkit.io.loaders import _unique_gene_names
from scrna_toolkit.provenance import record_step

logger = logging.getLogger(__name__)

DEFAULT_MAX_CELLS = 1_000_000


def _import_census() -> Any:
    try:
        import cellxgene_census
    except ImportError as err:
        raise ImportError(
            'Census support needs the optional dependency: pip install "scrna-toolkit[census]"'
        ) from err
    return cellxgene_census


def resolve_census_version(census_version: str = "stable") -> str:
    """Return the concrete build date (e.g. ``"2025-11-08"``) for a version or alias."""
    description = _import_census().get_census_version_description(census_version)
    return str(description["release_build"])


def _organism_key(organism: str) -> str:
    return organism.strip().lower().replace(" ", "_")


def _read_datasets(census: Any) -> pd.DataFrame:
    return census["census_info"]["datasets"].read().concat().to_pandas()


def _count_cells(census: Any, organism: str, obs_value_filter: str) -> int:
    obs = census["census_data"][_organism_key(organism)].obs
    return len(obs.read(value_filter=obs_value_filter, column_names=["soma_joinid"]).concat())


def list_census_datasets(census_version: str = "stable") -> pd.DataFrame:
    """Table of Census source datasets (id, title, collection, cell count, ...)."""
    cc = _import_census()
    version = resolve_census_version(census_version)
    with cc.open_soma(census_version=version) as census:
        return _read_datasets(census)


def load_census(
    obs_value_filter: str,
    *,
    organism: str = "Homo sapiens",
    census_version: str = "stable",
    var_value_filter: str | None = None,
    obs_column_names: Sequence[str] | None = None,
    max_cells: int | None = DEFAULT_MAX_CELLS,
) -> AnnData:
    """Query the Census and return raw counts as a cells x genes AnnData.

    Parameters
    ----------
    obs_value_filter
        SOMA filter on cell metadata, e.g.
        ``"tissue_general == 'lung' and disease == 'normal' and is_primary_data == True"``.
        Required: the full Census holds tens of millions of cells.
    census_version
        Build date (recommended for reproducibility) or ``"stable"``/``"latest"``.
    var_value_filter
        Optional gene filter, e.g. ``"feature_name in ['CD3E', 'MS4A1']"``.
    obs_column_names
        Subset of cell metadata columns to fetch (default: all).
    max_cells
        Refuse queries matching more cells than this (checked before fetching counts).
        ``None`` disables the check.

    ``X`` holds Census raw counts. Cells are indexed by the Census-wide ``soma_joinid``,
    genes by symbol (made unique) with Ensembl IDs in ``var["gene_ids"]``.
    """
    cc = _import_census()
    version = resolve_census_version(census_version)
    if "is_primary_data" not in obs_value_filter:
        logger.warning(
            "Filter does not mention is_primary_data; cells present in several datasets may be "
            "duplicated (add 'is_primary_data == True' to keep one copy)"
        )

    with cc.open_soma(census_version=version) as census:
        if max_cells is not None:
            n_cells = _count_cells(census, organism, obs_value_filter)
            if n_cells == 0:
                raise ValueError(f"Census query matched no cells: {obs_value_filter!r}")
            if n_cells > max_cells:
                raise ValueError(
                    f"Query matches {n_cells:,} cells (> max_cells={max_cells:,}); narrow the "
                    "filter or raise max_cells"
                )
            logger.info("Census %s query matches %d cells", version, n_cells)
        adata = cc.get_anndata(
            census,
            organism=organism,
            obs_value_filter=obs_value_filter,
            var_value_filter=var_value_filter,
            obs_column_names=list(obs_column_names) if obs_column_names else None,
        )
    if adata.n_obs == 0:
        raise ValueError(f"Census query matched no cells: {obs_value_filter!r}")

    adata.obs_names = adata.obs["soma_joinid"].astype(str).to_numpy()
    adata.var["gene_ids"] = adata.var["feature_id"].to_numpy()
    adata.var_names = _unique_gene_names(pd.Index(adata.var["feature_name"].astype(str).to_numpy()))
    logger.info("Loaded %d cells x %d genes from Census %s", adata.n_obs, adata.n_vars, version)
    record_step(
        adata,
        "load_census",
        {
            "census_version": version,
            "organism": organism,
            "obs_value_filter": obs_value_filter,
            "var_value_filter": var_value_filter,
        },
    )
    return adata


def download_census_h5ad(
    dataset_id: str,
    output_dir: str | Path,
    *,
    census_version: str = "stable",
    overwrite: bool = False,
) -> Path:
    """Download the original submitted H5AD of one Census dataset to ``output_dir``.

    Source files follow the CELLxGENE schema: ``X`` may hold normalized values with raw
    counts in ``adata.raw``. Check with :func:`scrna_toolkit.ensure_raw_counts` before QC.
    """
    cc = _import_census()
    version = resolve_census_version(census_version)
    dest = Path(output_dir) / f"{dataset_id}.h5ad"
    if dest.exists() and not overwrite:
        logger.info("%s exists; skipping download", dest)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    part.unlink(missing_ok=True)  # download_source_h5ad refuses to overwrite
    cc.download_source_h5ad(dataset_id, to_path=str(part), census_version=version)
    part.replace(dest)
    logger.info("Downloaded Census %s dataset %s to %s", version, dataset_id, dest)
    return dest
