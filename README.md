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

## Metadata conventions

`scrna_toolkit.metadata.CANONICAL_FIELDS` defines optional canonical `obs` columns:
`dataset`, `patient_id`, `donor_id`, `sample_id`, `batch`, `tissue`, `anatomical_site`,
`condition`, `treatment`, `treatment_status`, `timepoint`, `cell_type`, `cell_subtype`,
`response`. None are required, and missing metadata is never filled in.

```python
sct.harmonise_obs(adata, {"orig.ident": "sample_id", "patient": "patient_id"})  # originals kept
sct.merge_sample_metadata(adata, sample_table, on="sample_id")  # one row per sample
report = sct.validate_metadata(adata)  # DataFrame of issues; never raises
```

`validate_metadata` reports absent fields, duplicated cell IDs, missing/placeholder values,
samples that map to several patients (or other sample-level fields), and near-duplicate labels
such as `"T cell"` vs `"t_cell"`.

## Quality control

```python
sct.calculate_qc_metrics(adata, mitochondrial_prefix="MT-")  # adds obs columns, never filters
adata, summary = sct.filter_cells(
    adata, min_genes=200, max_pct_mt=20, return_summary=True
)  # returns a filtered copy unless inplace=True
adata = sct.filter_genes(adata, min_cells=3)
```

`pct_counts_mt` is a percentage (0-100). Use `mitochondrial_prefix="mt-"` for mouse,
`mitochondrial_genes=[...]` for an explicit list, and `gene_symbol_col=` when `var_names` are
Ensembl IDs. Ribosomal (`ribosomal_prefixes=("RPS", "RPL")`) and haemoglobin
(`hemoglobin_regex=r"^HB[^(P)]"`) fractions are opt-in.

**QC thresholds are biological and dataset-specific.** Filter functions do nothing unless you
pass thresholds; the values above are illustrations, not recommendations. Always inspect the QC
distributions for your tissue, chemistry and sequencing depth before choosing them.

## Preprocessing

```python
adata = sct.preprocess(adata, n_hvgs=3000, n_pcs=50)  # or preprocess(adata, config.preprocessing)
```

`preprocess` runs these steps, each also available on its own:

| Step | Input | Output |
| --- | --- | --- |
| `store_counts` | `X` (counts) | `layers["counts"]` (one copy) |
| `normalize_total` | `X` counts | `X` normalized to `target_sum` |
| `log1p` | `X` normalized | `X` log1p (all genes, sparse) |
| `select_hvgs` | `X` log1p (`layers["counts"]` for `seurat_v3`) | `var["highly_variable"]` |
| `run_pca` | scaled *temporary copy* of HVG columns | `obsm["X_pca"]`, `varm["PCs"]` |
| `run_neighbors` | `obsm["X_pca"]` | `obsp` graph |
| `run_umap` | graph | `obsm["X_umap"]` |
| `run_leiden` | graph | `obs["leiden"]` |

Steps work in place and return `adata` (`preprocess(..., copy=True)` works on a copy). All
stochastic steps take `seed` (default 0), so results are reproducible. The configuration is
stored in `adata.uns["scrna_toolkit"]["preprocessing_config"]` (`sct.get_config(adata, "preprocessing")`).

### How raw counts are handled

- **`adata.layers["counts"]`** holds the untouched raw count matrix. This is the canonical place
  for counts, and `store_counts` writes it once before normalization.
- **`adata.X`** starts as counts and is rewritten in place to normalized, then log1p values. It
  keeps all genes and stays sparse.
- **`adata.raw`** is not used by this package. It is a frozen snapshot of a whole AnnData and is
  often confused with raw counts.
- `adata.uns["scrna_toolkit"]["X_representation"]` records what `X` currently holds. Normalizing
  or log-transforming twice raises an error.

Scaling densifies data, so `run_pca(scale=True)` scales only a temporary dense copy of the HVG
columns and leaves `X` unchanged. A warning is logged when a dense block above ~0.8 GB is created.

## Summaries

```python
sct.summarize_dataset(adata)  # pd.Series: n_cells, n_samples, median_cells_per_sample, ...
sct.summarize_by(adata, "patient_id")  # pd.DataFrame: n_cells, fraction, n_samples, median QC
sct.summarize_by(adata, "cell_type")
```

Counts for absent metadata columns are `NaN`, and `missing_metadata` lists the absent
canonical fields.

## Development

```bash
pytest
ruff check .
ruff format --check .
mypy
```

