"""Typed configuration loaded from YAML.

Defaults are starting points, not recommendations: QC thresholds in particular are
dataset- and tissue-specific and should be chosen per dataset.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any

import yaml


@dataclass
class QCConfig:
    min_genes: int | None = 200
    max_genes: int | None = None
    min_counts: int | None = None
    max_counts: int | None = None
    max_pct_mt: float | None = 20.0
    min_cells_per_gene: int | None = 3
    mitochondrial_prefix: str = "MT-"
    ribosomal_prefixes: list[str] | None = None
    hemoglobin_regex: str | None = None


@dataclass
class PreprocessingConfig:
    target_sum: float = 1e4
    n_hvgs: int = 3000
    hvg_flavor: str = "seurat"
    hvg_batch_key: str | None = None
    scale: bool = True
    scale_max_value: float | None = 10.0
    n_pcs: int = 50
    n_neighbors: int = 15
    umap_min_dist: float = 0.5
    leiden_resolution: float = 1.0
    run_umap: bool = True
    run_leiden: bool = True
    seed: int = 0


# MetadataConfig field -> canonical obs column it maps onto.
_COLUMN_FIELDS = {
    "patient_col": "patient_id",
    "donor_col": "donor_id",
    "sample_col": "sample_id",
    "batch_col": "batch",
    "condition_col": "condition",
    "timepoint_col": "timepoint",
    "cell_type_col": "cell_type",
}


@dataclass
class MetadataConfig:
    """Names of the *source* columns in ``adata.obs`` for common canonical fields."""

    patient_col: str | None = None
    donor_col: str | None = None
    sample_col: str | None = None
    batch_col: str | None = None
    condition_col: str | None = None
    timepoint_col: str | None = None
    cell_type_col: str | None = None
    extra_mapping: dict[str, str] = field(default_factory=dict)

    def to_mapping(self) -> dict[str, str]:
        """Return a ``{source_column: canonical_column}`` mapping for ``harmonise_obs``."""
        mapping = {
            getattr(self, name): canonical
            for name, canonical in _COLUMN_FIELDS.items()
            if getattr(self, name) is not None
        }
        mapping.update(self.extra_mapping)
        return mapping


def _build(cls: type, values: dict[str, Any] | None) -> Any:
    values = dict(values or {})
    known = {f.name: f for f in fields(cls)}
    unknown = sorted(set(values) - set(known))
    if unknown:
        raise ValueError(f"Unknown {cls.__name__} keys: {unknown}. Valid keys: {sorted(known)}")
    defaults = cls()
    for key, value in values.items():
        # YAML parses '1e4' as a string; coerce numbers where the default is a float.
        if isinstance(getattr(defaults, key), float) and isinstance(value, int | str):
            values[key] = float(value)
    return cls(**values)


@dataclass
class Config:
    qc: QCConfig = field(default_factory=QCConfig)
    preprocessing: PreprocessingConfig = field(default_factory=PreprocessingConfig)
    metadata: MetadataConfig = field(default_factory=MetadataConfig)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Config:
        unknown = sorted(set(data) - {"qc", "preprocessing", "metadata"})
        if unknown:
            raise ValueError(f"Unknown config sections: {unknown}")
        return cls(
            qc=_build(QCConfig, data.get("qc")),
            preprocessing=_build(PreprocessingConfig, data.get("preprocessing")),
            metadata=_build(MetadataConfig, data.get("metadata")),
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def load_config(path: str | Path) -> Config:
    """Load a YAML config file. Missing sections/keys fall back to defaults."""
    with Path(path).open() as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Config file {path} must contain a mapping at the top level")
    return Config.from_dict(data)
