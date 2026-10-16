"""Dataset adapters: generic GEO datasets plus a registry for project-specific quirks."""

from scrna_toolkit.datasets.base import CensusDataset, DatasetAdapter, GEODataset
from scrna_toolkit.datasets.registry import get_dataset, list_datasets, register_dataset

__all__ = [
    "CensusDataset",
    "DatasetAdapter",
    "GEODataset",
    "get_dataset",
    "list_datasets",
    "register_dataset",
]
