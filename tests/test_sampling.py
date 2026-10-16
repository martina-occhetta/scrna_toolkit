import numpy as np
import pandas as pd
import pytest

import scrna_toolkit as sct
from scrna_toolkit.sampling import _proportional_with_floor, waterfill


def _obs(seed=0):
    """Unbalanced cohort: one huge donor, a rare cell type, two samples per donor."""
    rng = np.random.default_rng(seed)
    sizes = {"P0": 5000, "P1": 400, "P2": 300, "P3": 40}
    rows = []
    for patient, n in sizes.items():
        types = rng.choice(["T", "B", "Myeloid", "Rare"], size=n, p=[0.6, 0.3, 0.09, 0.01])
        samples = rng.choice([f"{patient}_a", f"{patient}_b"], size=n, p=[0.8, 0.2])
        rows += list(zip([patient] * n, samples, types, strict=True))
    obs = pd.DataFrame(rows, columns=["patient_id", "sample_id", "cell_type"])
    obs.index = [f"c{i}" for i in range(len(obs))]
    return obs


def test_waterfill_takes_small_groups_whole_and_caps_large_ones():
    sizes = pd.Series({"a": 5, "b": 100, "c": 1000})
    alloc = waterfill(sizes, 105)
    assert alloc.to_dict() == {"a": 5, "b": 50, "c": 50}
    assert waterfill(sizes, 106).sum() == 106  # remainder handed out one by one
    assert waterfill(sizes, 10_000).equals(sizes)
    with pytest.raises(ValueError):
        waterfill(sizes, -1)


def test_proportional_allocation_keeps_rare_types_and_hits_budget():
    sizes = pd.Series({"T": 900, "B": 90, "Rare": 3})
    alloc = _proportional_with_floor(sizes, 100, floor=10)
    assert alloc.sum() == 100 and alloc["Rare"] == 3 and alloc["B"] >= 10
    assert alloc["T"] > alloc["B"]
    assert _proportional_with_floor(sizes, 15, floor=10).sum() == 15  # floors exceed budget


def test_stratified_sample_balances_donors_and_weights_recover_composition():
    obs = _obs()

    weights = sct.stratified_sample(obs, 600, cell_type_key="cell_type", min_cells_per_type=5)

    assert len(weights) == 600 and weights.index.is_unique
    assert list(weights.index) == [i for i in obs.index if i in set(weights.index)]  # file order
    per_patient = obs.loc[weights.index, "patient_id"].value_counts()
    assert per_patient["P3"] == 40  # small donor taken whole
    assert per_patient["P0"] <= 300  # the large donor no longer dominates
    assert weights.sum() == pytest.approx(len(obs))
    full = obs["cell_type"].value_counts()
    weighted = weights.groupby(obs.loc[weights.index, "cell_type"]).sum()
    pd.testing.assert_series_equal(
        weighted.sort_index(), full.sort_index().astype(float), check_names=False
    )
    rare = obs["cell_type"].eq("Rare").groupby(obs["patient_id"]).sum()
    rare_picked = obs.loc[weights.index, "cell_type"].eq("Rare").groupby(obs["patient_id"]).sum()
    assert (rare_picked >= np.minimum(rare, 5)).all()  # floor keeps rare types in every donor


def test_nested_levels_and_reproducibility():
    obs = _obs()
    a = sct.stratified_sample(obs, 300, levels=("patient_id", "sample_id"), seed=1)
    b = sct.stratified_sample(obs, 300, levels=("patient_id", "sample_id"), seed=1)
    c = sct.stratified_sample(obs, 300, levels=("patient_id", "sample_id"), seed=2)
    pd.testing.assert_series_equal(a, b)
    assert not a.index.equals(c.index)
    per_sample = obs.loc[a.index].groupby(["patient_id", "sample_id"]).size()
    p0 = per_sample.loc["P0"]  # donors get (300 - 40) / 3 each, split equally across samples
    assert p0.sum() == 87 and abs(p0["P0_a"] - p0["P0_b"]) <= 1


def test_budget_larger_than_data_returns_everything():
    obs = _obs()
    weights = sct.stratified_sample(obs, 10**6, cell_type_key="cell_type")
    assert len(weights) == len(obs) and (weights == 1).all()


def test_missing_values():
    obs = _obs()
    obs.loc[obs.index[:10], "cell_type"] = np.nan  # its own stratum
    weights = sct.stratified_sample(obs, 200, cell_type_key="cell_type")
    assert len(weights) == 200 and weights.sum() == pytest.approx(len(obs))

    obs.loc[obs.index[0], "patient_id"] = None
    with pytest.raises(ValueError, match="missing"):
        sct.stratified_sample(obs, 200)
    with pytest.raises(KeyError, match="donor"):
        sct.stratified_sample(obs, 200, levels=("donor",))
    with pytest.raises(ValueError, match="positive"):
        sct.stratified_sample(obs, 0, levels=("sample_id",))


def test_downsample_cells_returns_weighted_copy(adata):
    adata.obs["patient_id"] = adata.obs["patient"]

    subset = sct.downsample_cells(adata, 30, cell_type_key="annotation", min_cells_per_type=2)

    assert subset.n_obs == 30 and subset.n_vars == adata.n_vars
    assert subset.obs["downsample_weight"].sum() == pytest.approx(adata.n_obs)
    assert subset.obs["patient_id"].value_counts().to_dict() == {"P0": 15, "P1": 15}
    np.testing.assert_array_equal(subset.X.toarray(), adata[subset.obs_names].X.toarray())
    assert "downsample_weight" not in adata.obs
    assert sct.get_history(subset)[-1]["params"]["n_obs_before"] == adata.n_obs
