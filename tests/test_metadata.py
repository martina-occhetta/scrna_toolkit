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


def test_merge_sample_metadata_strict_raises_on_unmatched_without_modifying(adata):
    sct.harmonise_obs(adata, {"orig.ident": "sample_id"})
    table = pd.DataFrame({"sample_id": ["S0", "S1", "S2"], "tissue": ["tumour", "blood", "tumour"]})
    n_steps = len(sct.get_history(adata))

    with pytest.raises(ValueError, match=r"1 samples .* \['S3'\]"):
        sct.merge_sample_metadata(adata, table, strict=True)

    assert "tissue" not in adata.obs
    assert len(sct.get_history(adata)) == n_steps


def test_merge_sample_metadata_strict_passes_when_all_samples_match(adata):
    sct.harmonise_obs(adata, {"orig.ident": "sample_id"})
    table = pd.DataFrame({"sample_id": ["S0", "S1", "S2", "S3", "S9"], "tissue": list("abcde")})

    sct.merge_sample_metadata(adata, table, strict=True)  # extra table rows are allowed

    assert adata.obs["tissue"].notna().all()


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


def _levelled_obs() -> pd.DataFrame:
    # Two patients, two samples each, three cells per sample.
    sample = [f"S{i // 3}" for i in range(12)]
    patient = [f"P{i // 6}" for i in range(12)]
    return pd.DataFrame(
        {
            "patient_id": patient,
            "sample_id": sample,
            "dataset": "D1",
            "subtype": ["HRD" if p == "P0" else "FBI" for p in patient],
            "site": [
                {"S0": "ovary", "S1": "omentum", "S2": "ovary", "S3": "ascites"}[s] for s in sample
            ],
            "cell_type": ["T", "B", "Myeloid"] * 4,
            "n_counts": np.arange(12.0),
        },
        index=[f"cell{i}" for i in range(12)],
    )


def test_infer_field_levels_returns_coarsest_constant_level():
    levels = sct.infer_field_levels(_levelled_obs())

    assert levels.to_dict() == {
        "patient_id": "patient_id",
        "sample_id": "sample_id",
        "dataset": "constant",
        "subtype": "patient_id",
        "site": "sample_id",
        "cell_type": "cell",
        "n_counts": "cell",
    }


def test_infer_field_levels_accepts_adata_columns_and_extra_levels(adata):
    sct.harmonise_obs(adata, {"orig.ident": "sample_id", "patient": "patient_id"})
    # Every sample (cells i % 4) holds two libraries, each from one sort fraction.
    fraction = (np.arange(adata.n_obs) // 4) % 2
    adata.obs["library_id"] = adata.obs["sample_id"].astype(str) + "_L" + fraction.astype(str)
    adata.obs["sort"] = np.where(fraction, "CD45+", "CD45-")

    levels = sct.infer_field_levels(
        adata, levels=("patient_id", "sample_id", "library_id"), columns=["sort", "orig.ident"]
    )

    assert levels.to_dict() == {"sort": "library_id", "orig.ident": "sample_id"}


def test_infer_field_levels_counts_missing_values_and_keys():
    obs = _levelled_obs()
    obs.loc["cell0", "subtype"] = np.nan  # sample S0 now has two values (NaN, HRD)
    assert sct.infer_field_levels(obs, columns=["subtype"])["subtype"] == "cell"
    obs.loc[obs["sample_id"] == "S0", "subtype"] = np.nan  # S0 missing, S1 HRD: same patient
    assert sct.infer_field_levels(obs, columns=["subtype"])["subtype"] == "sample_id"

    obs = _levelled_obs()
    obs.loc[["cell0", "cell6"], "patient_id"] = np.nan  # cells from two patients share a NaN key
    assert sct.infer_field_levels(obs, columns=["subtype"])["subtype"] == "sample_id"


def test_infer_field_levels_skips_absent_levels_and_rejects_absent_columns():
    obs = _levelled_obs().drop(columns="patient_id")
    assert sct.infer_field_levels(obs, columns=["subtype"])["subtype"] == "sample_id"
    assert sct.infer_field_levels(obs, levels=("donor",), columns=["subtype"])["subtype"] == "cell"
    with pytest.raises(KeyError, match="not_there"):
        sct.infer_field_levels(obs, columns=["not_there"])


def test_metadata_audit_table_summarises_each_column():
    obs = _levelled_obs()
    obs.loc["cell0", "cell_type"] = np.nan
    obs.loc["cell2", "cell_type"] = "B"  # B 5, T 3, Myeloid 3
    obs.loc["cell1", "n_counts"] = np.nan

    table = sct.metadata_audit_table(obs).set_index("column")

    assert list(table.columns) == ["level", "dtype", "n_unique", "n_missing", "examples"]
    assert table.loc["subtype", "level"] == "patient_id"
    assert table.loc["cell_type", "level"] == "cell"
    assert table.loc["cell_type", ["n_unique", "n_missing"]].tolist() == [3, 1]
    assert table.loc["cell_type", "examples"].split(", ")[0] == "B"  # most frequent first
    assert table.loc["n_counts", "examples"] == "range 0 – 11"
    assert table.loc["n_counts", "n_missing"] == 1
    assert table.loc["dataset", ["level", "n_unique", "examples"]].tolist() == ["constant", 1, "D1"]


def test_metadata_audit_table_columns_examples_and_descriptions(adata):
    adata.obs["flag"] = np.arange(adata.n_obs) % 3 == 0
    adata.obs["empty"] = np.nan

    table = sct.metadata_audit_table(
        adata,
        columns=["annotation", "flag", "empty"],
        levels=("patient",),
        descriptions={"annotation": "author cell type"},
        n_examples=2,
    ).set_index("column")

    assert table.index.tolist() == ["annotation", "flag", "empty"]
    assert table["description"].tolist() == ["author cell type", "", ""]
    assert len(table.loc["annotation", "examples"].split(", ")) == 2
    assert table.loc["flag", "examples"] == "False, True"  # booleans are not ranges
    assert table.loc["empty", ["n_unique", "n_missing", "examples"]].tolist() == [0, 120, ""]


SITE_MAP = {"S0": "ovary", "S1": "omentum", "S2": "ovary", "S3": "ascites"}


def test_map_values_writes_target_and_records_provenance(adata):
    out = sct.map_values(adata, "orig.ident", SITE_MAP, "anatomical_site")

    assert out is adata
    expected = adata.obs["orig.ident"].map(SITE_MAP)
    assert (adata.obs["anatomical_site"] == expected).all()
    assert "orig.ident" in adata.obs
    step = sct.get_history(adata)[-1]
    assert step["step"] == "map_values"
    assert step["params"] == {
        "column": "orig.ident",
        "target": "anatomical_site",
        "mapping": SITE_MAP,
        "unmapped": [],
    }


def test_map_values_strict_raises_on_unmapped_without_modifying(adata):
    n_steps = len(sct.get_history(adata))
    with pytest.raises(ValueError, match=r"1 values of 'orig.ident' .* \['S3'\]"):
        sct.map_values(adata, "orig.ident", {"S0": "a", "S1": "b", "S2": "c"}, "site")
    assert "site" not in adata.obs
    assert len(sct.get_history(adata)) == n_steps


def test_map_values_non_strict_leaves_unmapped_missing(adata):
    sct.map_values(adata, "orig.ident", {"S0": "a", "S1": "b", "S2": "c"}, "site", strict=False)

    s3 = adata.obs["orig.ident"] == "S3"
    assert adata.obs.loc[s3, "site"].isna().all()
    assert adata.obs.loc[~s3, "site"].notna().all()
    assert sct.get_history(adata)[-1]["params"]["unmapped"] == ["S3"]


def test_map_values_keeps_missing_and_categorical_dtype(adata):
    adata.obs["sort"] = pd.Categorical(
        np.where(np.arange(adata.n_obs) % 2, "singlet, live, CD45+", "singlet, live, CD45-")
    )
    adata.obs.loc[adata.obs_names[:3], "sort"] = np.nan
    mapping = {"singlet, live, CD45+": "CD45+", "singlet, live, CD45-": "CD45-"}

    sct.map_values(adata, "sort", mapping, "sort_fraction")

    assert isinstance(adata.obs["sort_fraction"].dtype, pd.CategoricalDtype)
    assert adata.obs["sort_fraction"].isna().sum() == 3
    assert set(adata.obs["sort_fraction"].dropna()) == {"CD45+", "CD45-"}


def test_map_values_refuses_overwrite_and_missing_column(adata):
    sct.map_values(adata, "orig.ident", SITE_MAP, "site")
    with pytest.raises(ValueError, match="already exists"):
        sct.map_values(adata, "orig.ident", SITE_MAP, "site")
    sct.map_values(adata, "orig.ident", SITE_MAP, "orig.ident", overwrite=True)  # in place
    assert set(adata.obs["orig.ident"]) == {"ovary", "omentum", "ascites"}
    with pytest.raises(KeyError, match="not_there"):
        sct.map_values(adata, "not_there", SITE_MAP, "x")


def test_categorize_obs_converts_string_columns_only(adata):
    adata.obs["site"] = np.where(np.arange(adata.n_obs) % 2, "ovary", None)
    adata.obs["mixed"] = pd.Series(
        [1, "a"] * (adata.n_obs // 2), index=adata.obs_names, dtype=object
    )
    adata.obs["n_counts"] = np.arange(adata.n_obs)
    adata.obs["flag"] = np.arange(adata.n_obs) % 2 == 0
    before = adata.obs.copy()

    out = sct.categorize_obs(adata)

    assert out is adata
    for column in ("orig.ident", "patient", "annotation", "site"):
        assert isinstance(adata.obs[column].dtype, pd.CategoricalDtype), column
        pd.testing.assert_series_equal(
            adata.obs[column].astype(object), before[column].astype(object), check_dtype=False
        )
    assert adata.obs["site"].isna().sum() == before["site"].isna().sum()
    assert adata.obs["mixed"].dtype == object
    assert pd.api.types.is_integer_dtype(adata.obs["n_counts"])
    assert pd.api.types.is_bool_dtype(adata.obs["flag"])
    assert sct.get_history(adata)[-1]["params"]["converted"] == [
        "orig.ident",
        "patient",
        "annotation",
        "site",
    ]


def test_categorize_obs_respects_limit_and_columns(adata):
    adata.obs["barcode"] = adata.obs_names.astype(str)

    sct.categorize_obs(adata, max_categories=10)
    assert not isinstance(adata.obs["barcode"].dtype, pd.CategoricalDtype)
    assert sct.get_history(adata)[-1]["params"]["kept_as_strings"] == ["barcode"]

    sct.categorize_obs(adata, columns=["barcode"])
    assert isinstance(adata.obs["barcode"].dtype, pd.CategoricalDtype)
    with pytest.raises(KeyError, match="not_there"):
        sct.categorize_obs(adata, columns=["not_there"])
