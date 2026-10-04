import pandas as pd


# ------------------------------------------------------------
# Load Step 3 and Step 4 results
# ------------------------------------------------------------

baseline = pd.read_csv(
    "ipapy2_ecmdb_context_candidates.csv"
)

biochem = pd.read_csv(
    "ipapy2_ecmdb_biochemical_candidates.csv"
)


# ------------------------------------------------------------
# Keep ECMDB-supported candidates
# ------------------------------------------------------------

baseline_ecmdb = baseline[
    baseline["ecmdb_match"] == True
].copy()


# ------------------------------------------------------------
# Add biochemical-context information
# ------------------------------------------------------------

biochem_cols = [
    "feature_id",
    "ipapy2_id",
    "name_ipapy2",
    "name_ecmdb",
    "kegg_id",
    "post",
    "biochemical_context"
]

biochem_small = biochem[
    biochem_cols
].drop_duplicates()


# ------------------------------------------------------------
# Find features where:
#
# 1. There is more than one ECMDB candidate
# 2. At least one candidate has biochemical context
# 3. Therefore biological context helps distinguish
#    among competing candidates.
# ------------------------------------------------------------

candidate_counts = (
    baseline_ecmdb
    .groupby("feature_id")
    .size()
    .reset_index(
        name="n_ecmdb_candidates"
    )
)


context_counts = (
    biochem_small
    .groupby("feature_id")
    .agg(
        n_biochemical_candidates=(
            "biochemical_context",
            "sum"
        ),
        biochemical_context=(
            "biochemical_context",
            "any"
        )
    )
    .reset_index()
)


summary = candidate_counts.merge(
    context_counts,
    on="feature_id",
    how="left"
)


# Features with ambiguity + biochemical support
examples = summary[
    (summary["n_ecmdb_candidates"] >= 2)
    &
    (summary["n_biochemical_candidates"] >= 1)
].copy()


# Prefer examples with:
# - several competing candidates
# - more than one candidate with biochemical context
examples = examples.sort_values(
    [
        "n_ecmdb_candidates",
        "n_biochemical_candidates"
    ],
    ascending=False
)


print()
print("Potential examples:")
print(
    examples.head(10).to_string(
        index=False
    )
)


# ------------------------------------------------------------
# Print candidate-level details for the best examples
# ------------------------------------------------------------

print()
print("=" * 70)
print("DETAILED EXAMPLES")
print("=" * 70)


for feature_id in examples.head(5)["feature_id"]:

    print()
    print("FEATURE:", feature_id)

    rows = biochem_small[
        biochem_small["feature_id"]
        == feature_id
    ].copy()

    rows = rows.sort_values(
        "post",
        ascending=False
    )

    print(
        rows[
            [
                "ipapy2_id",
                "name_ipapy2",
                "name_ecmdb",
                "kegg_id",
                "post",
                "biochemical_context"
            ]
        ].to_string(
            index=False
        )
    )


# ------------------------------------------------------------
# Save examples
# ------------------------------------------------------------

example_ids = examples.head(5)[
    "feature_id"
]

example_table = biochem_small[
    biochem_small["feature_id"].isin(
        example_ids
    )
].copy()

example_table.to_csv(
    "ipapy2_context_examples.csv",
    index=False
)

print()
print("Saved:")
print("  ipapy2_context_examples.csv")