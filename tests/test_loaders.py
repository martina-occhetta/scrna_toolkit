import gzip

import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.io
import scipy.sparse as sp

import scrna_toolkit as sct
from scrna_toolkit.io import FormatDetectionError

N_CELLS, N_GENES = 6, 4
GENE_IDS = [f"ENSG{i:05d}" for i in range(N_GENES)]
SYMBOLS = ["MT-CO1", "ACTB", "ACTB", "CD3E"]  # duplicated symbol on purpose
BARCODES = [f"AAAC{i}-1" for i in range(N_CELLS)]


def _counts():
    rng = np.random.default_rng(0)
    return rng.integers(0, 5, size=(N_CELLS, N_GENES)).astype(np.float32)


def _write(path, text, gz):
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)


def write_10x_mtx(directory, *, gz=True, legacy=False, prefix="", transpose=False):
    directory.mkdir(parents=True, exist_ok=True)
    sfx = ".gz" if gz else ""
    matrix = sp.coo_matrix(_counts() if transpose else _counts().T)  # 10x is genes x cells
    mtx = directory / f"{prefix}matrix.mtx"
    scipy.io.mmwrite(mtx, matrix)
    if gz:
        with open(mtx, "rb") as src, gzip.open(f"{mtx}.gz", "wb") as dst:
            dst.write(src.read())
        mtx.unlink()
    _write(directory / f"{prefix}barcodes.tsv{sfx}", "\n".join(BARCODES) + "\n", gz)
    if legacy:
        lines = [f"{i}\t{s}" for i, s in zip(GENE_IDS, SYMBOLS, strict=True)]
        _write(directory / f"{prefix}genes.tsv{sfx}", "\n".join(lines) + "\n", gz)
    else:
        lines = [f"{i}\t{s}\tGene Expression" for i, s in zip(GENE_IDS, SYMBOLS, strict=True)]
        _write(directory / f"{prefix}features.tsv{sfx}", "\n".join(lines) + "\n", gz)
    return directory


def write_10x_h5(path):
    matrix = sp.csc_matrix(_counts().T.astype(np.int32))  # genes x cells, CSC over cells
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g["data"], g["indices"], g["indptr"] = matrix.data, matrix.indices, matrix.indptr
        g["shape"] = np.array(matrix.shape)
        g["barcodes"] = np.array(BARCODES, dtype="S")
        feat = g.create_group("features")
        feat["id"] = np.array(GENE_IDS, dtype="S")
        feat["name"] = np.array(SYMBOLS, dtype="S")
        feat["feature_type"] = np.array(["Gene Expression"] * N_GENES, dtype="S")
        feat["genome"] = np.array(["GRCh38"] * N_GENES, dtype="S")
    return path


def _check(adata):
    assert adata.shape == (N_CELLS, N_GENES)
    assert list(adata.obs_names) == BARCODES
    assert list(adata.var["gene_ids"]) == GENE_IDS
    assert adata.var_names.is_unique
    assert sp.issparse(adata.X)
    np.testing.assert_array_equal(adata.X.toarray(), _counts())


@pytest.mark.parametrize("gz", [True, False])
@pytest.mark.parametrize("legacy", [True, False])
def test_load_10x_mtx_variants(tmp_path, gz, legacy):
    directory = write_10x_mtx(tmp_path / "mtx", gz=gz, legacy=legacy)
    assert sct.detect_format(directory) == "10x_mtx"
    _check(sct.load_scrna(directory))


def test_load_10x_mtx_geo_prefix_and_file_path(tmp_path):
    directory = write_10x_mtx(tmp_path / "geo", prefix="GSM1_")
    _check(sct.load_scrna(directory))
    _check(sct.load_scrna(directory / "GSM1_matrix.mtx.gz"))


def test_load_10x_mtx_cells_by_genes_matrix_is_detected(tmp_path):
    _check(sct.load_scrna(write_10x_mtx(tmp_path / "t", transpose=True)))


def test_directory_with_several_matrices_is_ambiguous(tmp_path):
    write_10x_mtx(tmp_path / "d", prefix="GSM1_")
    write_10x_mtx(tmp_path / "d", prefix="GSM2_")
    with pytest.raises(FormatDetectionError, match="prefix"):
        sct.load_scrna(tmp_path / "d")
    _check(sct.load_scrna(tmp_path / "d", prefix="GSM2_"))


def test_directory_without_matrix_raises(tmp_path):
    (tmp_path / "notes.txt").write_text("x")
    with pytest.raises(FormatDetectionError):
        sct.detect_format(tmp_path)


def test_load_10x_h5(tmp_path):
    path = write_10x_h5(tmp_path / "filtered_feature_bc_matrix.h5")
    assert sct.detect_format(path) == "10x_h5"
    _check(sct.load_scrna(path))


def test_h5ad_roundtrip(tmp_path, adata):
    path = tmp_path / "x.h5ad"
    adata.write_h5ad(path)
    loaded = sct.load_scrna(path)
    assert loaded.shape == adata.shape
    assert sct.get_history(loaded)[-1]["params"]["format"] == "h5ad"


@pytest.mark.parametrize("sep,suffix", [(",", "csv"), ("\t", "tsv.gz")])
def test_delimited_genes_x_cells(tmp_path, sep, suffix):
    df = pd.DataFrame(_counts().T, index=[f"G{i}" for i in range(N_GENES)], columns=BARCODES)
    path = tmp_path / f"counts.{suffix}"
    df.to_csv(path, sep=sep)

    adata = sct.load_scrna(path, orientation="genes_x_cells")

    assert adata.shape == (N_CELLS, N_GENES)
    assert list(adata.obs_names) == BARCODES
    np.testing.assert_array_equal(adata.X.toarray(), _counts())


def test_delimited_requires_orientation(tmp_path):
    path = tmp_path / "counts.csv"
    pd.DataFrame(_counts(), index=BARCODES).to_csv(path)
    with pytest.raises(FormatDetectionError, match="orientation"):
        sct.load_scrna(path)
    assert sct.load_scrna(path, orientation="cells_x_genes").shape == (N_CELLS, N_GENES)


def test_unknown_extension_is_ambiguous(tmp_path):
    path = tmp_path / "counts.txt.gz"
    path.write_bytes(b"")
    with pytest.raises(FormatDetectionError, match="format="):
        sct.load_scrna(path)


def test_load_samples_concatenates_with_sample_labels(tmp_path):
    paths = {s: write_10x_mtx(tmp_path / s) for s in ["S1", "S2"]}

    adata = sct.load_samples(paths)

    assert adata.n_obs == 2 * N_CELLS
    assert adata.obs_names.is_unique
    assert adata.obs["sample_id"].value_counts().to_dict() == {"S1": N_CELLS, "S2": N_CELLS}
