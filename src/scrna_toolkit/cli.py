"""Command-line interface: ``scrna-toolkit <command> ...``.

Exit code 0 on success, 1 on any error (logged), 130 on interrupt.
"""

from __future__ import annotations

import argparse
import dataclasses
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import anndata as ad
import pandas as pd
from anndata import AnnData

from scrna_toolkit import __version__
from scrna_toolkit.config import Config, load_config
from scrna_toolkit.io import (
    download_census_h5ad,
    download_geo_processed_files,
    get_geo_metadata,
    list_geo_samples,
    list_geo_supplementary_files,
    load_census,
    load_scrna,
)
from scrna_toolkit.io.loaders import FORMATS
from scrna_toolkit.metadata import harmonise_obs
from scrna_toolkit.preprocessing import preprocess
from scrna_toolkit.provenance import record_config
from scrna_toolkit.qc import calculate_qc_metrics, filter_cells, filter_genes, qc_report
from scrna_toolkit.summary import summarize_by, summarize_dataset

logger = logging.getLogger("scrna_toolkit.cli")


def _print_table(table: pd.DataFrame | pd.Series) -> None:
    with pd.option_context("display.max_rows", 200, "display.width", 160):
        print(table.to_string())


def cmd_geo_info(args: argparse.Namespace) -> None:
    meta = get_geo_metadata(args.accession, refresh=args.refresh)
    print(f"{meta['accession']}: {meta['title']}")
    print(f"Organisms: {', '.join(meta['organisms'])}")
    print(f"Platforms: {', '.join(meta['platforms'])}")
    print(f"Samples:   {len(meta['sample_ids'])}")
    print(f"\n{meta['summary']}\n")
    _print_table(list_geo_supplementary_files(args.accession, include_samples=args.include_samples))
    if args.samples:
        print()
        _print_table(list_geo_samples(args.accession))


def cmd_download_geo(args: argparse.Namespace) -> None:
    paths = download_geo_processed_files(
        args.accession,
        args.output,
        patterns=args.pattern,
        kinds=args.kind or ("processed",),
        include_samples=args.include_samples,
        overwrite=args.overwrite,
    )
    for path in paths:
        print(path)


def cmd_download_census(args: argparse.Namespace) -> None:
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(f"{args.output} exists; pass --overwrite to replace it")
    adata = load_census(
        args.filter,
        organism=args.organism,
        census_version=args.census_version,
        var_value_filter=args.var_filter,
        obs_column_names=args.obs_columns.split(",") if args.obs_columns else None,
        max_cells=args.max_cells or None,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"{args.output} ({adata.n_obs} cells x {adata.n_vars} genes)")


def cmd_download_census_dataset(args: argparse.Namespace) -> None:
    print(
        download_census_h5ad(
            args.dataset_id,
            args.output,
            census_version=args.census_version,
            overwrite=args.overwrite,
        )
    )


def _load(args: argparse.Namespace) -> AnnData:
    return load_scrna(args.input, args.format, orientation=args.orientation, prefix=args.prefix)


def cmd_preprocess(args: argparse.Namespace) -> None:
    config = load_config(args.config) if args.config else Config()
    adata = _load(args)

    mapping = config.metadata.to_mapping()
    if mapping:
        harmonise_obs(adata, mapping)

    qc = config.qc
    calculate_qc_metrics(
        adata,
        qc.mitochondrial_prefix,
        ribosomal_prefixes=qc.ribosomal_prefixes,
        hemoglobin_regex=qc.hemoglobin_regex,
    )
    filter_cells(
        adata,
        min_genes=qc.min_genes,
        max_genes=qc.max_genes,
        min_counts=qc.min_counts,
        max_counts=qc.max_counts,
        max_pct_mt=qc.max_pct_mt,
        inplace=True,
    )
    filter_genes(adata, min_cells=qc.min_cells_per_gene, inplace=True)
    record_config(adata, "qc", dataclasses.asdict(qc))
    preprocess(adata, config.preprocessing)

    if args.report_dir:
        qc_report(adata, output_dir=args.report_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    logger.info("Wrote %s (%d cells x %d genes)", args.output, adata.n_obs, adata.n_vars)


def cmd_summarize(args: argparse.Namespace) -> None:
    # Only obs is needed: open h5ad files in backed mode so X is never loaded.
    if args.input.suffix == ".h5ad":
        adata = ad.read_h5ad(args.input, backed="r")
    else:
        adata = _load(args)
    _print_table(summarize_dataset(adata))
    for column in args.groupby:
        print(f"\nBy {column}:")
        _print_table(summarize_by(adata, column))


def cmd_qc_report(args: argparse.Namespace) -> None:
    adata = _load(args)
    if "total_counts" not in adata.obs:
        calculate_qc_metrics(adata)
    report = qc_report(adata, output_dir=args.output)
    print(f"Wrote {len(report.tables)} tables and {len(report.figures)} figures to {args.output}")


def _add_input_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--format", choices=list(FORMATS))
    parser.add_argument("--orientation", choices=["cells_x_genes", "genes_x_cells"])
    parser.add_argument("--prefix", help="10x MTX filename prefix, e.g. GSM123_")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scrna-toolkit", description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("geo-info", help="show GEO series metadata and supplementary files")
    p.add_argument("accession")
    p.add_argument("--samples", action="store_true", help="also list samples")
    p.add_argument("--include-samples", action="store_true", help="list per-sample files too")
    p.add_argument("--refresh", action="store_true", help="ignore the metadata cache")
    p.set_defaults(func=cmd_geo_info)

    p = sub.add_parser("download-geo", help="download processed GEO supplementary files")
    p.add_argument("accession")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--pattern", action="append", help="filename glob (repeatable)")
    p.add_argument(
        "--kind",
        action="append",
        choices=["processed", "raw_archive", "sequence"],
        help="file kinds to include (repeatable; default: processed)",
    )
    p.add_argument("--include-samples", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(func=cmd_download_geo)

    p = sub.add_parser(
        "download-census", help="query the CELLxGENE Census (needs the [census] extra)"
    )
    p.add_argument("--filter", required=True, help="obs value filter, e.g. \"tissue == 'lung'\"")
    p.add_argument("--output", type=Path, required=True, help="output .h5ad")
    p.add_argument("--organism", default="Homo sapiens")
    p.add_argument("--census-version", default="stable", help="build date (recommended) or alias")
    p.add_argument("--var-filter", help="gene filter, e.g. \"feature_name in ['CD3E']\"")
    p.add_argument("--obs-columns", help="comma-separated obs columns to fetch (default: all)")
    p.add_argument(
        "--max-cells", type=int, default=1_000_000, help="refuse larger queries (0 disables)"
    )
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(func=cmd_download_census)

    p = sub.add_parser(
        "download-census-dataset", help="download a Census dataset's original source H5AD"
    )
    p.add_argument("dataset_id")
    p.add_argument("--output", type=Path, required=True, help="output directory")
    p.add_argument("--census-version", default="stable")
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(func=cmd_download_census_dataset)

    p = sub.add_parser("preprocess", help="load, harmonise, QC-filter and preprocess")
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True, help="output .h5ad")
    p.add_argument("--config", type=Path, help="YAML config (defaults if omitted)")
    p.add_argument("--report-dir", type=Path, help="also write a QC report here")
    _add_input_options(p)
    p.set_defaults(func=cmd_preprocess)

    p = sub.add_parser("summarize", help="print dataset summary tables")
    p.add_argument("input", type=Path)
    p.add_argument("--groupby", action="append", default=[], help="obs column (repeatable)")
    _add_input_options(p)
    p.set_defaults(func=cmd_summarize)

    p = sub.add_parser("qc-report", help="write QC tables and figures")
    p.add_argument("input", type=Path)
    p.add_argument("--output", type=Path, required=True)
    _add_input_options(p)
    p.set_defaults(func=cmd_qc_report)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    try:
        args.func(args)
    except KeyboardInterrupt:
        return 130
    except Exception as err:
        logger.error("%s: %s", type(err).__name__, err, exc_info=args.verbose)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
