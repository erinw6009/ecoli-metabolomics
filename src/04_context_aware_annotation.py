# ============================================================
# Step 3: Context-aware annotation using ECMDB
#
# Goal:
# Add E. coli metabolome knowledge to the baseline ipaPy2
# annotations.
#
# Baseline evidence:
#   ipaPy2 MS1 annotation
#
# Contextual evidence:
#   Candidate metabolite is present in ECMDB
#
# Important:
# ECMDB membership is used as biological prior/context.
# It does NOT prove that the metabolite is present in the
# experimental sample.
# ============================================================

import json
import pickle
import pandas as pd


# ------------------------------------------------------------
# 1. Load ECMDB
# ------------------------------------------------------------

ECMDB_PATH = "data/ecmdb/ecmdb.json"

with open(ECMDB_PATH, "r") as f:
    ecmdb = json.load(f)

ecmdb_df = pd.DataFrame(ecmdb)

print("ECMDB records:", len(ecmdb_df))


# ------------------------------------------------------------
# 2. Load baseline ipaPy2 annotations
# ------------------------------------------------------------

with open(
    "ipapy2_baseline_annotations.pkl",
    "rb"
) as f:
    annotations = pickle.load(f)

print("ipaPy2 annotated features:", len(annotations))


# ------------------------------------------------------------
# 3. Convert ipaPy2 dictionaries of DataFrames into
#    one long table
# ------------------------------------------------------------

rows = []

for feature_id, candidates in annotations.items():

    for _, row in candidates.iterrows():

        rows.append({
            "feature_id": feature_id,
            "ipapy2_id": row.get("id"),
            "name": row.get("name"),
            "formula": row.get("formula"),
            "adduct": row.get("adduct"),
            "mz": row.get("m/z"),
            "charge": row.get("charge"),
            "ppm": row.get("ppm"),
            "prior": row.get("prior"),
            "post": row.get("post")
        })


ipa_df = pd.DataFrame(rows)

print("Total ipaPy2 candidate rows:", len(ipa_df))


# ------------------------------------------------------------
# 4. Normalise metabolite names
#
# This is deliberately simple and transparent.
# We will record the mapping method rather than hiding it.
# ------------------------------------------------------------

def normalise_name(x):

    if pd.isna(x):
        return ""

    return (
        str(x)
        .lower()
        .strip()
    )


ipa_df["name_normalised"] = ipa_df["name"].apply(
    normalise_name
)

ecmdb_df["name_normalised"] = ecmdb_df["name"].apply(
    normalise_name
)


# ------------------------------------------------------------
# 5. Create ECMDB lookup table
# ------------------------------------------------------------

ecmdb_lookup = ecmdb_df[
    [
        "met_id",
        "name",
        "name_normalised",
        "kegg_id",
        "hmdb_id",
        "chebi_id",
        "pubchem_id",
        "moldb_inchikey",
        "moldb_formula",
        "moldb_mono_mass"
    ]
].copy()


# ------------------------------------------------------------
# 6. Match ipaPy2 candidates to ECMDB
#
# Primary mapping used here:
# exact normalised metabolite name.
#
# We keep the mapping method explicit because database
# cross-referencing can produce ambiguous matches.
# ------------------------------------------------------------

context_df = ipa_df.merge(
    ecmdb_lookup,
    on="name_normalised",
    how="left",
    suffixes=("_ipapy2", "_ecmdb")
)


# ------------------------------------------------------------
# 7. Add context flag
# ------------------------------------------------------------

context_df["ecmdb_match"] = (
    context_df["met_id"].notna()
)


context_df["context"] = context_df[
    "ecmdb_match"
].map({
    True: "ECMDB-supported",
    False: "No ECMDB exact-name match"
})


context_df["mapping_method"] = context_df[
    "ecmdb_match"
].map({
    True: "exact normalised name",
    False: "no exact normalised name match"
})


# ------------------------------------------------------------
# 8. Save complete candidate table
# ------------------------------------------------------------

context_df.to_csv(
    "ipapy2_ecmdb_context_candidates.csv",
    index=False
)


# ------------------------------------------------------------
# 9. Create feature-level summary
#
# This tells us whether each feature has at least one
# ECMDB-supported candidate.
# ------------------------------------------------------------

feature_summary = (
    context_df
    .groupby("feature_id")
    .agg(
        n_candidates=("name_ipapy2", "count"),
        n_ecmdb_candidates=("ecmdb_match", "sum"),
        ecmdb_supported=("ecmdb_match", "any")
    )
    .reset_index()
)


feature_summary.to_csv(
    "ipapy2_ecmdb_feature_summary.csv",
    index=False
)


# ------------------------------------------------------------
# 10. Print summary
# ------------------------------------------------------------

print()
print("ECMDB context summary")
print("---------------------")

print(
    "Candidate rows:",
    len(context_df)
)

print(
    "Candidates with ECMDB match:",
    int(context_df["ecmdb_match"].sum())
)

print(
    "Features with >=1 ECMDB-supported candidate:",
    int(feature_summary["ecmdb_supported"].sum())
)

print(
    "Features without ECMDB-supported candidate:",
    int(
        (~feature_summary["ecmdb_supported"]).sum()
    )
)

print()
print("Saved:")
print("  ipapy2_ecmdb_context_candidates.csv")
print("  ipapy2_ecmdb_feature_summary.csv")

print()
print("STEP 3 COMPLETE")