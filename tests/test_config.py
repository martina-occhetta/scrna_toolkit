import pytest

from scrna_toolkit.config import Config, MetadataConfig, load_config

YAML = """
qc:
  min_genes: 200
  max_pct_mt: 20
  min_cells_per_gene: 3

preprocessing:
  target_sum: 1e4
  n_hvgs: 3000
  n_pcs: 50
  n_neighbors: 15
  leiden_resolution: 1.0

metadata:
  patient_col: patient
  sample_col: sample
  cell_type_col: annotation
"""


def test_load_config_from_yaml(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(YAML)

    config = load_config(path)

    assert config.qc.min_genes == 200
    assert config.qc.max_pct_mt == 20.0
    assert config.preprocessing.target_sum == 10000.0
    assert config.preprocessing.n_neighbors == 15
    assert config.metadata.to_mapping() == {
        "patient": "patient_id",
        "sample": "sample_id",
        "annotation": "cell_type",
    }


def test_missing_sections_use_defaults(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text("qc:\n  min_genes: 500\n")

    config = load_config(path)

    assert config.qc.min_genes == 500
    assert config.preprocessing == Config().preprocessing


def test_unknown_keys_are_rejected():
    with pytest.raises(ValueError, match="min_gene"):
        Config.from_dict({"qc": {"min_gene": 200}})
    with pytest.raises(ValueError, match="sections"):
        Config.from_dict({"qcc": {}})


def test_extra_mapping_is_merged():
    meta = MetadataConfig(sample_col="orig.ident", extra_mapping={"site": "anatomical_site"})
    assert meta.to_mapping() == {"orig.ident": "sample_id", "site": "anatomical_site"}


def test_to_dict_roundtrip():
    config = Config()
    assert Config.from_dict(config.to_dict()) == config
