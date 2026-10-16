import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

import scrna_toolkit as sct  # noqa: E402

MAPPING = {"orig.ident": "sample_id", "patient": "patient_id", "annotation": "cell_type"}


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def annotated(adata):
    sct.harmonise_obs(adata, MAPPING)
    sct.calculate_qc_metrics(adata)
    return adata


@pytest.mark.parametrize(
    "plot",
    [
        lambda a: sct.plot_qc_distributions(a),
        lambda a: sct.plot_qc_distributions(a, groupby="sample_id"),
        lambda a: sct.plot_qc_scatter(a),
        lambda a: sct.plot_cells_per_sample(a),
        lambda a: sct.plot_cells_per_patient(a, log=True),
        lambda a: sct.plot_cell_type_composition(a),
        lambda a: sct.plot_cell_type_composition(a, "patient_id", normalize=False),
        lambda a: sct.plot_composition_heatmap(a),
    ],
)
def test_plots_return_figures(annotated, plot):
    assert isinstance(plot(annotated), Figure)


def test_plots_skip_missing_columns(adata, caplog):
    assert sct.plot_qc_distributions(adata) is None
    assert sct.plot_qc_scatter(adata) is None
    assert sct.plot_cells_per_patient(adata) is None
    assert sct.plot_cell_type_composition(adata) is None
    assert sct.plot_composition_heatmap(adata) is None
    assert sct.plot_umap(adata, color="cell_type") is None
    assert "not found" in caplog.text


def test_missing_groupby_falls_back_to_ungrouped(annotated):
    assert isinstance(sct.plot_qc_distributions(annotated, groupby="tissue"), Figure)


def test_saves_only_when_path_given(annotated, tmp_path):
    path = tmp_path / "plots" / "cells.png"
    sct.plot_cells_per_sample(annotated, path=path)
    assert path.is_file()


def test_many_categories(annotated):
    annotated.obs["sample_id"] = [f"S{i}" for i in range(annotated.n_obs)]  # 120 samples
    fig = sct.plot_cell_type_composition(annotated)
    assert fig.axes[0].get_xticklabels()[0].get_text() == ""  # labels hidden


def test_uses_scanpy_colors_when_present(annotated):
    annotated.obs["cell_type"] = annotated.obs["cell_type"].astype("category")
    annotated.uns["cell_type_colors"] = ["#111111", "#222222", "#333333"]
    fig = sct.plot_cell_type_composition(annotated)
    used = {matplotlib.colors.to_hex(p.get_facecolor()) for p in fig.axes[0].patches}
    assert used == {"#111111", "#222222", "#333333"}


def test_plot_umap_skips_unknown_keys(annotated):
    annotated.obsm["X_umap"] = np.random.default_rng(0).normal(size=(annotated.n_obs, 2))
    fig = sct.plot_umap(annotated, color=["cell_type", "not_a_column", "GENE1"])
    assert isinstance(fig, Figure)
    assert isinstance(sct.plot_umap(annotated), Figure)
