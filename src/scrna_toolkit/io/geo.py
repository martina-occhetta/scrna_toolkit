"""Programmatic access to NCBI GEO series metadata and supplementary files.

Metadata comes from GEO's plain-text SOFT endpoint (``acc.cgi?form=text``), not HTML pages,
and is cached locally (see :func:`scrna_toolkit.io.cache.cache_dir`). Downloads default to
processed supplementary files only. Raw ``*_RAW.tar`` archives and sequence files
(FASTQ/SRA/BAM) are only fetched when explicitly requested via ``kinds``; reads deposited in
SRA are never downloaded by this module.
"""

from __future__ import annotations

import fnmatch
import logging
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from scrna_toolkit.io.cache import cached_text
from scrna_toolkit.io.download import download_file, fetch_text

logger = logging.getLogger(__name__)

GEO_ACC_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
_GSE = re.compile(r"^GSE\d+$")
_SEQUENCE = re.compile(r"\.(fastq|fq|sra|bam|cram)(\.gz)?$", re.IGNORECASE)


class GEOError(RuntimeError):
    """Raised when GEO returns no record for an accession."""


def _validate(accession: str) -> str:
    accession = accession.strip().upper()
    if not _GSE.match(accession):
        raise ValueError(f"Expected a GEO series accession like 'GSE12345', got {accession!r}")
    return accession


def _fetch_soft(accession: str, target: str, refresh: bool) -> str:
    def fetch() -> str:
        params = {"acc": accession, "targ": target, "form": "text", "view": "brief"}
        text = fetch_text(GEO_ACC_URL, params=params)
        if not text.lstrip().startswith("^"):  # GEO answers unknown accessions with HTML
            raise GEOError(f"GEO returned no record for {accession}")
        return text

    return cached_text(f"geo/{accession}_{target}.soft", fetch, refresh=refresh)


def parse_soft(text: str) -> list[tuple[str, str, dict[str, list[str]]]]:
    """Parse SOFT text into ``(entity_type, accession, {attribute: [values]})`` records."""
    entities: list[tuple[str, str, dict[str, list[str]]]] = []
    for line in text.splitlines():
        if line.startswith("^"):
            kind, _, acc = line[1:].partition(" = ")
            entities.append((kind.strip(), acc.strip(), {}))
        elif line.startswith("!") and entities:
            key, _, value = line[1:].partition(" = ")
            entities[-1][2].setdefault(key.strip(), []).append(value.strip())
    return entities


def _https(url: str) -> str:
    return re.sub(r"^ftp://", "https://", url)


def get_geo_metadata(accession: str, *, refresh: bool = False) -> dict[str, Any]:
    """Series-level metadata: title, summary, design, organisms, platforms, samples, files.

    ``attributes`` holds every raw SOFT field (``Series_`` prefix removed) as a list.
    """
    accession = _validate(accession)
    entities = [e for e in parse_soft(_fetch_soft(accession, "self", refresh)) if e[0] == "SERIES"]
    if not entities:
        raise GEOError(f"No SERIES record in GEO response for {accession}")
    attrs = {k.removeprefix("Series_"): v for k, v in entities[0][2].items()}

    def first(key: str) -> str:
        return attrs.get(key, [""])[0]

    return {
        "accession": accession,
        "title": first("title"),
        "summary": "\n".join(attrs.get("summary", [])),
        "overall_design": "\n".join(attrs.get("overall_design", [])),
        "organisms": sorted(set(attrs.get("sample_organism", []))),
        "platforms": attrs.get("platform_id", []),
        "sample_ids": attrs.get("sample_id", []),
        "supplementary_files": [_https(u) for u in attrs.get("supplementary_file", [])],
        "pubmed_ids": attrs.get("pubmed_id", []),
        "attributes": attrs,
    }


def _characteristics(values: list[str]) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for value in values:
        key, sep, val = value.partition(":")
        if sep:
            parsed[re.sub(r"\W+", "_", key.strip().lower()).strip("_")] = val.strip()
        else:
            parsed["characteristics"] = "; ".join(
                filter(None, [parsed.get("characteristics"), value])
            )
    return parsed


def _sample_records(accession: str, refresh: bool) -> list[tuple[str, dict[str, list[str]]]]:
    return [
        (acc, a)
        for kind, acc, a in parse_soft(_fetch_soft(accession, "gsm", refresh))
        if kind == "SAMPLE"
    ]


def list_geo_samples(accession: str, *, refresh: bool = False) -> pd.DataFrame:
    """One row per GSM sample: title, source, organism and parsed ``key: value`` characteristics.

    The table can be joined onto cells with :func:`scrna_toolkit.merge_sample_metadata`.
    """
    accession = _validate(accession)
    rows = []
    for gsm, attrs in _sample_records(accession, refresh):
        row = {
            "sample_accession": gsm,
            "title": attrs.get("Sample_title", [""])[0],
            "source_name": attrs.get("Sample_source_name_ch1", [""])[0],
            "organism": attrs.get("Sample_organism_ch1", [""])[0],
        }
        row.update(_characteristics(attrs.get("Sample_characteristics_ch1", [])))
        rows.append(row)
    return pd.DataFrame(rows)


def _kind(filename: str) -> str:
    if filename.endswith("_RAW.tar"):
        return "raw_archive"
    if _SEQUENCE.search(filename):
        return "sequence"
    return "processed"


def list_geo_supplementary_files(
    accession: str, *, include_samples: bool = False, refresh: bool = False
) -> pd.DataFrame:
    """Supplementary files with ``level`` (series/sample), ``filename``, ``url`` and ``kind``.

    ``kind`` is ``processed``, ``raw_archive`` (``*_RAW.tar`` bundles) or ``sequence``.
    """
    meta = get_geo_metadata(accession, refresh=refresh)
    rows = [("series", meta["accession"], url) for url in meta["supplementary_files"]]
    if include_samples:
        for gsm, attrs in _sample_records(meta["accession"], refresh):
            for key, values in attrs.items():
                if key.startswith("Sample_supplementary_file"):
                    rows += [("sample", gsm, _https(v)) for v in values if v.upper() != "NONE"]
    table = pd.DataFrame(rows, columns=["level", "accession", "url"])
    table.insert(2, "filename", table["url"].str.rsplit("/", n=1).str[-1])
    table["kind"] = table["filename"].map(_kind)
    return table


def download_geo_processed_files(
    accession: str,
    output_dir: str | Path,
    *,
    patterns: Sequence[str] | None = None,
    kinds: Sequence[str] = ("processed",),
    include_samples: bool = False,
    overwrite: bool = False,
    refresh: bool = False,
) -> list[Path]:
    """Download supplementary files into ``output_dir/<accession>/``.

    Parameters
    ----------
    patterns
        Optional filename globs, e.g. ``["*counts*", "*metadata*"]``.
    kinds
        File kinds to include; add ``"raw_archive"`` or ``"sequence"`` explicitly to fetch them.
    include_samples
        Also consider per-sample (GSM) supplementary files.
    """
    files = list_geo_supplementary_files(
        accession, include_samples=include_samples, refresh=refresh
    )
    selected = files[files["kind"].isin(list(kinds))]
    if patterns:
        keep = selected["filename"].map(lambda f: any(fnmatch.fnmatch(f, p) for p in patterns))
        selected = selected[keep]
    skipped = files.loc[~files.index.isin(selected.index), "filename"].tolist()
    if skipped:
        logger.info("Not downloading %d files: %s", len(skipped), skipped[:10])
    if selected.empty:
        logger.warning("No supplementary files of kinds %s matched for %s", list(kinds), accession)
        return []

    target = Path(output_dir) / _validate(accession)
    return [
        download_file(row.url, target / row.filename, overwrite=overwrite)
        for row in selected.itertuples()
    ]
