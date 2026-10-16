"""Shared fixtures: small synthetic AnnData objects, no downloads."""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData


def make_adata(n_obs: int = 120, n_vars: int = 200, seed: int = 0) -> AnnData:
    """Sparse negative-binomial counts with a few mitochondrial genes and sample metadata."""
    rng = np.random.default_rng(seed)
    counts = rng.negative_binomial(n=2, p=0.3, size=(n_obs, n_vars)).astype(np.float32)
    counts[rng.random(counts.shape) < 0.6] = 0
    var_names = [f"MT-G{i}" for i in range(5)] + [f"GENE{i}" for i in range(n_vars - 5)]
    obs = pd.DataFrame(
        {
            "orig.ident": [f"S{i % 4}" for i in range(n_obs)],
            "patient": [f"P{(i % 4) // 2}" for i in range(n_obs)],
            "annotation": rng.choice(["T", "B", "Myeloid"], size=n_obs),
        },
        index=[f"cell{i}" for i in range(n_obs)],
    )
    return AnnData(X=sp.csr_matrix(counts), obs=obs, var=pd.DataFrame(index=var_names))


@pytest.fixture
def adata() -> AnnData:
    return make_adata()
