import numpy as np
import pandas as pd
import pytest

import scrna_toolkit as sct

MAPPING = {"orig.ident": "sample_id", "patient": "patient_id", "annotation": "cell_type"}


def test_harmonise_obs_copies_columns_and_keeps_originals(adata):
    out = sct.harmonise_obs(adata, MAPPING)

    assert out is adata
    for source, target in MAPPING.items():
        assert source in adata.obs
        pd.testing.assert_series_equal(adata.obs[target], adata.obs[source], check_names=False)
    assert sct.get_history(adata)[-1]["params"]["mapping"] == MAPPING


def test_harmonise_obs_missing_source_skips_or_raises(adata):
    sct.harmonise_obs(adata, {"not_there": "tissue"})
    assert "tissue" not in adata.obs  # never fabricated
    with pytest.raises(KeyError):
        sct.harmonise_obs(adata, {"not_there": "tissue"}, strict=True)


def test_harmonise_obs_refuses_to_overwrite_and_duplicate_targets(adata):
    sct.harmonise_obs(adata, {"orig.ident": "sample_id"})
    with pytest.raises(ValueError, match="already exists"):
        sct.harmonise_obs(adata, {"patient": "sample_id"})
    sct.harmonise_obs(adata, {"patient": "sample_id"}, overwrite=True)
    with pytest.raises(ValueError, match="same target"):
        sct.harmonise_obs(adata, {"patient": "patient_id", "annotation": "patient_id"})


def test_merge_sample_metadata_broadcasts_and_reports_unmatched(adata):
    sct.harmonise_obs(adata, {"orig.ident": "sample_id"})
    table = pd.DataFrame({"sample_id": ["S0", "S1", "S2"], "tissue": ["tumour", "blood", "tumour"]})

    sct.merge_sample_metadata(adata, table)

    by_sample = adata.obs.groupby("sample_id", observed=True)["tissue"].first()
    assert by_sample["S0"] == "tumour" and by_sample["S1"] == "blood"
    assert adata.obs.loc[adata.obs["sample_id"] == "S3", "tissue"].isna().all()
    assert sct.get_history(adata)[-1]["params"]["unmatched_samples"] == ["S3"]


def test_merge_sample_metadata_rejects_duplicate_keys(adata):
    sct.harmonise_obs(adata, {"orig.ident": "sample_id"})
    table = pd.DataFrame({"sample_id": ["S0", "S0"], "tissue": ["a", "b"]})
    with pytest.raises(ValueError, match="duplicated"):
        sct.merge_sample_metadata(adata, table)


def _checks(report, check):
    return report[report["check"] == check]


def test_validate_clean_metadata_reports_only_absent_fields(adata):
    sct.harmonise_obs(adata, MAPPING)
    report = sct.validate_metadata(adata)
    assert set(report["severity"]) == {"info"}
    assert "tissue" in set(report["column"])
    assert "sample_id" not in set(report["column"])


def test_validate_detects_problems(adata):
    sct.harmonise_obs(adata, MAPPING)
    adata.obs["patient_id"] = adata.obs["patient_id"].astype(str)
    adata.obs.iloc[0, adata.obs.columns.get_loc("patient_id")] = "P9"  # S0 now spans 2 patients
    adata.obs.iloc[1, adata.obs.columns.get_loc("patient_id")] = np.nan
    cell_type = adata.obs["cell_type"].astype(str).to_numpy()
    cell_type[:3] = ["t", "T ", "B"]
    adata.obs["cell_type"] = cell_type
    adata.obs_names = ["dup"] * 2 + list(adata.obs_names[2:])

    report = sct.validate_metadata(adata)

    assert _checks(report, "duplicated_cell_ids")["severity"].tolist() == ["error"]
    assert _checks(report, "missing_values")["column"].tolist() == ["patient_id"]
    assert _checks(report, "one_to_many")["column"].tolist() == ["patient_id"]
    assert _checks(report, "inconsistent_labels")["column"].tolist() == ["cell_type"]


def test_validate_without_any_metadata_does_not_fail(adata):
    adata.obs = adata.obs.iloc[:, :0]
    report = sct.validate_metadata(adata)
    assert (report["check"] == "missing_field").all()
