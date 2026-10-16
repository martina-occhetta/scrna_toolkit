# scrna-toolkit

A reusable Python package for common single-cell RNA-seq data ingestion, metadata harmonisation,
QC, preprocessing, and exploratory analysis, built on AnnData and Scanpy.

> Work in progress: the package is being built incrementally. See the sections below for what is
> currently available.

## Installation

```bash
pip install -e ".[dev]"
```

## Provenance

Every transformation records its parameters in `adata.uns["scrna_toolkit"]["history"]`
(parameters are stored as JSON strings so they survive `write_h5ad`). Read it back with:

```python
import scrna_toolkit as sct

sct.get_history(adata)
```

## Development

```bash
pytest
ruff check .
ruff format --check .
mypy
```
