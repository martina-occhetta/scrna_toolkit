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

### Checking that data are raw counts

Functions that assume raw counts (`calculate_qc_metrics`, `store_counts`, `normalize_total`,
and `select_hvgs` with `seurat_v3`) call `sct.ensure_raw_counts`. It scans the whole
matrix in row chunks, without densifying, and raises if it finds negative values, non-integer
values, or a tracked `X_representation` other than counts:

```python
sct.ensure_raw_counts(adata)  # adata.X
sct.ensure_raw_counts(adata, layer="counts")
sct.ensure_raw_counts(adata, allow_non_integer=True)  # e.g. SoupX-corrected counts
```

`allow_non_integer=True` (or `allow_non_integer_counts: true` in the `qc`/`preprocessing`
config) accepts fractional values, but then normalized data can no longer be detected. Negative
values always fail.

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

## Plotting

```python
sct.plot_qc_distributions(adata, groupby="sample_id")
sct.plot_qc_scatter(adata)
sct.plot_umap(adata, color=["patient_id", "cell_type"])
sct.plot_cells_per_sample(adata)
sct.plot_cell_type_composition(adata, groupby="patient_id")
sct.plot_composition_heatmap(adata, path="figures/composition.png")
```

Every plot returns a matplotlib `Figure` and saves only when `path=` is given. If required
metadata is absent, the plot is skipped with a logged warning and returns `None`, so it is
safe to call in generic pipelines. Colours come from Scanpy's `uns["<col>_colors"]` when present,
otherwise from the matplotlib defaults. Tick labels are hidden for more than 60 categories.

## QC report

```python
report = sct.qc_report(adata, output_dir="reports/qc")  # output_dir is optional
report.tables["qc_summary"]
report.figures["qc_distributions"]
```

The report includes a dataset summary, metadata issues, QC metric quantiles, cells per
sample/patient/cell type, and QC, composition and UMAP figures. Anything whose data is
absent is skipped. With `output_dir`, tables are written as CSV and figures as PNG.

## GEO

```python
meta = sct.get_geo_metadata("GSE165897")  # title, summary, samples, files, ...
samples = sct.list_geo_samples("GSE165897")  # one row per GSM with parsed characteristics
files = sct.list_geo_supplementary_files("GSE165897", include_samples=True)
paths = sct.download_geo_processed_files("GSE165897", "data/raw", patterns=["*counts*"])
```

Metadata comes from GEO's plain-text SOFT endpoint, not HTML scraping, and is cached under
`~/.cache/scrna_toolkit` (override with `$SCRNA_TOOLKIT_CACHE`; pass `refresh=True` to
re-fetch). Downloads include only processed supplementary files by default. `*_RAW.tar`
archives (`kinds=["raw_archive"]`) and sequence files such as FASTQ/BAM (`kinds=["sequence"]`)
must be requested explicitly, and reads stored in SRA are never downloaded. Network calls
retry transient failures and raise `DownloadError`/`GEOError` with the URL involved.

## Datasets and adapters

Most GEO series need no custom code:

```python
dataset = sct.get_dataset(
    "GSE123456",
    data_dir="data/raw",
    patterns=["*filtered*"],  # which processed files to download
    obs_mapping={"orig.ident": "sample_id"},
)
adata = dataset.get()  # download (existing files skipped) -> load -> harmonise
```

`GEODataset` loads the single count matrix in `data/raw/<accession>/`. Pass
`matrix_file=` when there are several, and `load_kwargs=` for loader options. Write an adapter
only for awkward datasets, for example when the matrix and cell metadata are in separate files:

```python
import pandas as pd
import scrna_toolkit as sct


@sct.register_dataset("my-ovarian-atlas")
class MyOvarianAtlas(sct.GEODataset):
    def __init__(self, data_dir="data/raw"):
        super().__init__(
            "GSE165897",
            data_dir,
            matrix_file="GSE165897_UMIcounts_HGSOC.tsv.gz",
            load_kwargs={"orientation": "genes_x_cells"},
        )

    def harmonise_metadata(self, adata):
        # Dataset-specific quirk: cell metadata ships as a separate table keyed by barcode.
        info_path = self.directory / "GSE165897_cellInfo_HGSOC.tsv.gz"
        adata.obs = adata.obs.join(pd.read_csv(info_path, sep="\t", index_col="cell"))
        return sct.harmonise_obs(
            adata,
            {
                "sample": "sample_id",
                "treatment_phase": "treatment_status",
                "anatomical_location": "anatomical_site",
            },
        )


adata = sct.get_dataset("my-ovarian-atlas").get()
```

Adapters belong in project code. The package does not ship hard-coded biological datasets.

## Command line

```bash
scrna-toolkit geo-info GSE165897 --samples
scrna-toolkit download-geo GSE165897 --output data/raw/ --pattern "*counts*"
scrna-toolkit preprocess --input data/raw/dataset.h5 --config preprocessing.yaml \
    --output data/processed/dataset.h5ad --report-dir reports/qc
scrna-toolkit summarize data/processed/dataset.h5ad --groupby patient_id
scrna-toolkit qc-report data/processed/dataset.h5ad --output reports/qc
```

`preprocess` runs load → harmonise (`metadata` section) → QC metrics → cell/gene filtering
(`qc` section) → `preprocess` (`preprocessing` section) and writes an h5ad that records both
configs. Without `--config` the defaults are used, so pass your own thresholds. `summarize`
opens h5ad files in backed mode (only `obs` is read). Delimited inputs need `--orientation`.
Every command logs to stderr, `-v` enables debug logging with tracebacks, and failures exit
with a non-zero code. `python -m scrna_toolkit` works too.

## Development

```bash
pytest
ruff check .
ruff format --check .
mypy
```

