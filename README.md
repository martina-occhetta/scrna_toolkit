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

## Configuration

Pipelines can be configured from YAML (`qc`, `preprocessing`, `metadata` sections; every key is
optional and unknown keys raise an error):

```python
from scrna_toolkit.config import load_config

config = load_config("config.yaml")
config.metadata.to_mapping()  # {"patient": "patient_id", ...}
```

## Loading data

```python
import scrna_toolkit as sct

adata = sct.load_scrna("filtered_feature_bc_matrix.h5")  # format inferred from the path
adata = sct.load_scrna("GSE123_RAW/", prefix="GSM1_")  # 10x MTX with GEO-style prefix
adata = sct.load_scrna("counts.tsv.gz", orientation="genes_x_cells")
adata = sct.load_samples({"S1": "s1/", "S2": "s2/"})  # concatenate, adds obs["sample_id"]
```

Supported: `.h5ad`, 10x `.h5`, 10x MTX directories (`matrix.mtx` + `barcodes.tsv` +
`features.tsv`/`genes.tsv`, gzipped or not), and dense CSV/TSV tables. Detection is
conservative: ambiguous inputs raise `FormatDetectionError`. Delimited tables always need an
explicit `orientation`, since it cannot be inferred reliably. Gene symbols become `var_names`
and Ensembl IDs are kept in `var["gene_ids"]`.

## Development

```bash
pytest
ruff check .
ruff format --check .
mypy
```

