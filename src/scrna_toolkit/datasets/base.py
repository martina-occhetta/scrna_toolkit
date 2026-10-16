"""Dataset adapter interface and the generic GEO implementation.

Most GEO series work with :class:`GEODataset` and a couple of arguments. Write a subclass
only for awkward cases (several matrices to combine, separate cell-metadata files, ...).
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from anndata import AnnData

from scrna_toolkit.io.geo import download_geo_processed_files
from scrna_toolkit.io.loaders import FormatDetectionError, detect_format, load_scrna
from scrna_toolkit.metadata.harmonise import harmonise_obs

logger = logging.getLogger(__name__)


class DatasetAdapter:
    """Override :meth:`download`, :meth:`load` and optionally :meth:`harmonise_metadata`."""

    def __init__(self, data_dir: str | Path = "data/raw") -> None:
        self.data_dir = Path(data_dir)

    def download(self) -> list[Path]:
        raise NotImplementedError

    def load(self) -> AnnData:
        raise NotImplementedError

    def harmonise_metadata(self, adata: AnnData) -> AnnData:
        return adata

    def get(self, *, download: bool = True) -> AnnData:
        """Download (files already present are skipped), load and harmonise."""
        if download:
            self.download()
        return self.harmonise_metadata(self.load())


class GEODataset(DatasetAdapter):
    """A GEO series whose processed files hold a single loadable count matrix.

    Parameters
    ----------
    accession
        GEO series accession, e.g. ``"GSE12345"``. Files go to ``data_dir/<accession>/``.
    patterns
        Filename globs restricting which processed files are downloaded.
    matrix_file
        File (or 10x directory) to load, relative to the dataset directory. Needed when
        several loadable files are present.
    load_kwargs
        Passed to :func:`load_scrna` (e.g. ``{"orientation": "genes_x_cells"}``).
    obs_mapping
        Passed to :func:`harmonise_obs` after loading.
    """

    def __init__(
        self,
        accession: str,
        data_dir: str | Path = "data/raw",
        *,
        patterns: Sequence[str] | None = None,
        matrix_file: str | None = None,
        load_kwargs: Mapping[str, Any] | None = None,
        obs_mapping: Mapping[str, str] | None = None,
    ) -> None:
        super().__init__(data_dir)
        self.accession = accession.upper()
        self.patterns = patterns
        self.matrix_file = matrix_file
        self.load_kwargs = dict(load_kwargs or {})
        self.obs_mapping = dict(obs_mapping or {})

    @property
    def directory(self) -> Path:
        return self.data_dir / self.accession

    def download(self) -> list[Path]:
        return download_geo_processed_files(self.accession, self.data_dir, patterns=self.patterns)

    def _matrix_path(self) -> Path:
        if self.matrix_file is not None:
            return self.directory / self.matrix_file
        if not self.directory.is_dir():
            raise FileNotFoundError(f"{self.directory} does not exist; call download() first")
        try:
            detect_format(self.directory)
            return self.directory  # a 10x MTX directory
        except FormatDetectionError:
            pass
        candidates = []
        for path in sorted(self.directory.iterdir()):
            try:
                detect_format(path)
                candidates.append(path)
            except FormatDetectionError:
                continue
        if len(candidates) != 1:
            names = [p.name for p in candidates]
            raise FormatDetectionError(
                f"Found {len(candidates)} loadable files in {self.directory} ({names}); "
                "pass matrix_file=... to choose one"
            )
        return candidates[0]

    def load(self) -> AnnData:
        adata = load_scrna(self._matrix_path(), **self.load_kwargs)
        if "dataset" not in adata.obs:
            adata.obs["dataset"] = self.accession
        return adata

    def harmonise_metadata(self, adata: AnnData) -> AnnData:
        return harmonise_obs(adata, self.obs_mapping) if self.obs_mapping else adata
