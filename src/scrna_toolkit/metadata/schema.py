"""Optional canonical column names for ``adata.obs``.

None of these fields are required. Functions across the package use these names as
defaults (e.g. ``sample_col="sample_id"``) so that harmonised datasets need no extra arguments.
"""

CANONICAL_FIELDS: tuple[str, ...] = (
    "dataset",
    "patient_id",
    "donor_id",
    "sample_id",
    "batch",
    "tissue",
    "anatomical_site",
    "condition",
    "treatment",
    "treatment_status",
    "timepoint",
    "cell_type",
    "cell_subtype",
    "response",
)

# Fields expected to take a single value per sample; used to flag one-to-many mappings.
SAMPLE_LEVEL_FIELDS: tuple[str, ...] = (
    "dataset",
    "patient_id",
    "donor_id",
    "batch",
    "tissue",
    "anatomical_site",
    "condition",
    "treatment",
    "treatment_status",
    "timepoint",
    "response",
)
