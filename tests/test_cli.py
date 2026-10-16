import anndata as ad
import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

import scrna_toolkit as sct  # noqa: E402
from scrna_toolkit import cli  # noqa: E402

CONFIG = """
qc:
  min_genes: 10
  max_pct_mt: 50
  min_cells_per_gene: 3
preprocessing:
  n_hvgs: 50
  n_pcs: 10
  n_neighbors: 10
metadata:
  sample_col: orig.ident
  patient_col: patient
  cell_type_col: annotation
"""


@pytest.fixture
def h5ad(adata, tmp_path):
    path = tmp_path / "raw.h5ad"
    adata.write_h5ad(path)
    return path


def test_preprocess_command(h5ad, tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text(CONFIG)
    output = tmp_path / "processed" / "out.h5ad"

    code = cli.main(
        [
            "preprocess",
            "--input",
            str(h5ad),
            "--config",
            str(config),
            "--output",
            str(output),
            "--report-dir",
            str(tmp_path / "report"),
        ]
    )

    assert code == 0
    result = ad.read_h5ad(output)
    assert {"sample_id", "patient_id", "cell_type", "leiden"} <= set(result.obs)
    assert {"X_pca", "X_umap"} <= set(result.obsm) and "counts" in result.layers
    assert sct.get_config(result, "qc")["min_genes"] == 10
    assert sct.get_config(result, "preprocessing")["n_hvgs"] == 50
    assert (tmp_path / "report" / "qc_summary.csv").is_file()


def test_summarize_command(h5ad, capsys):
    assert cli.main(["summarize", str(h5ad), "--groupby", "patient"]) == 0
    out = capsys.readouterr().out
    assert "n_cells" in out and "By patient" in out


def test_qc_report_command(h5ad, tmp_path):
    assert cli.main(["qc-report", str(h5ad), "--output", str(tmp_path / "qc")]) == 0
    assert (tmp_path / "qc" / "qc_distributions.png").is_file()


def test_geo_commands(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(
        cli,
        "get_geo_metadata",
        lambda acc, refresh: {
            "accession": acc,
            "title": "T",
            "organisms": ["Homo sapiens"],
            "platforms": ["GPL1"],
            "sample_ids": ["GSM1"],
            "summary": "S",
        },
    )
    monkeypatch.setattr(
        cli,
        "list_geo_supplementary_files",
        lambda acc, include_samples: pd.DataFrame({"filename": ["f.tsv.gz"]}),
    )
    monkeypatch.setattr(cli, "list_geo_samples", lambda acc: pd.DataFrame({"title": ["s1"]}))
    calls = {}
    monkeypatch.setattr(
        cli,
        "download_geo_processed_files",
        lambda acc, out, **kw: calls.update(kw) or [out / "f.tsv.gz"],
    )

    assert cli.main(["geo-info", "GSE1", "--samples"]) == 0
    out = capsys.readouterr().out
    assert "GSE1: T" in out and "f.tsv.gz" in out and "s1" in out

    assert (
        cli.main(["download-geo", "GSE1", "--output", str(tmp_path), "--pattern", "*.tsv.gz"]) == 0
    )
    assert calls["patterns"] == ["*.tsv.gz"] and calls["kinds"] == ("processed",)


def test_errors_give_nonzero_exit(tmp_path, caplog):
    assert cli.main(["summarize", str(tmp_path / "missing.h5ad")]) == 1
    assert "missing.h5ad" in caplog.text
    with pytest.raises(SystemExit):
        cli.main(["no-such-command"])


def test_census_commands(monkeypatch, tmp_path, capsys, adata):
    calls = {}

    def fake_load_census(obs_value_filter, **kwargs):
        calls["load"] = {"filter": obs_value_filter, **kwargs}
        return adata

    monkeypatch.setattr(cli, "load_census", fake_load_census)
    monkeypatch.setattr(
        cli,
        "download_census_h5ad",
        lambda dataset_id, out, **kw: calls.update(dataset=(dataset_id, kw)) or out / "x.h5ad",
    )
    output = tmp_path / "census" / "lung.h5ad"

    args = [
        "download-census",
        "--filter",
        "tissue == 'lung'",
        "--output",
        str(output),
        "--census-version",
        "2025-11-08",
        "--obs-columns",
        "cell_type,donor_id",
        "--max-cells",
        "0",
    ]
    assert cli.main(args) == 0
    assert ad.read_h5ad(output).shape == adata.shape
    assert calls["load"] == {
        "filter": "tissue == 'lung'",
        "organism": "Homo sapiens",
        "census_version": "2025-11-08",
        "var_value_filter": None,
        "obs_column_names": ["cell_type", "donor_id"],
        "max_cells": None,
    }
    assert cli.main(args) == 1  # refuses to overwrite
    assert cli.main([*args, "--overwrite"]) == 0

    assert cli.main(["download-census-dataset", "abc", "--output", str(tmp_path)]) == 0
    assert calls["dataset"] == ("abc", {"census_version": "stable", "overwrite": False})


def test_census_command_without_extra_fails_cleanly(monkeypatch, tmp_path, caplog):
    import sys

    monkeypatch.setitem(sys.modules, "cellxgene_census", None)
    code = cli.main(["download-census", "--filter", "x", "--output", str(tmp_path / "o.h5ad")])
    assert code == 1
    assert "scrna-toolkit[census]" in caplog.text
